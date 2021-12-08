using SparseArrays
using LinearAlgebra
using Rotations, SharedArrays
using Printf
using GLMakie
inline!(true)
using Colors
using GeometryTypes
using Contour
using BSON

earliesttime = 540.; time_step = 2.; lasttime = 0.

nx = 360; ny = 180
xcoords = collect(-180. + 180/nx: 360/nx: 180. - 180/nx)
ycoords = collect(-90. + 180/nx: 360 / nx: 90. - 180/nx)
areabox = fill(0.,ny)
delta_x = fill(0.,ny); delta_y = 110.e3
for iy in 1:ny
    areabox[iy] = 110. * 110. * cos(ycoords[iy] / 180. * pi)  * 1e6 # m2
    delta_x[iy] = 110.e3 * cos(ycoords[iy] / 180. * pi) # m
end
notyetdefined = -2; notatsurface = -1; ocean_type = 1
continent_type = 2; uplifted_continent_type = 3

ntimeslices = ( earliesttime - lasttime ) /
    time_step + 1
ages = []
for i in 1:ntimeslices
    age = earliesttime - time_step * (i-1)
    push!(ages,age)   # eg [50 45 40 35 30 25 20 15 10 5 0]
end
ages = convert(Array{Float32},ages)

sediment_layer_time_bins = [635,541,250,65,0] # resolution of the sediment record on the plates
n_sediment_layer_time_bins = length(sediment_layer_time_bins)
sediment_layer_time_bins = convert(Array{Float32},sediment_layer_time_bins) # keeps netcdf library from puking

geo_interval_time_bins = [ 540, 485, 444, 419,  # these are the beginnings
    359, 299, 251, 201, 145,
    66, 56, 34, 23, 5, 2, 0 ]
geo_interval_time_bins = convert(Array{Float32},geo_interval_time_bins) # keeps netcdf library from puking
geo_interval_names = [ "Cambrian", "Ordovician", "Silurian", "Devonian",
    "Carboniferous","Permian","Triassic","Jurassic","Cretaceous",
    "Paleocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene"]
function get_geo_interval()
    age = world.age
    return get_geo_interval(age)
end
function get_geo_interval(age)
    next_i = search_sorted_below(geo_interval_time_bins,age)
    last_i = next_i - 1
    if age == 0
        return "Anthropocene"
    elseif last_i > 0
        last_time = geo_interval_time_bins[last_i]
        next_time = geo_interval_time_bins[next_i]
        if last_time - age < ( last_time - next_time ) / 3.
            qualifier = "early "
        elseif last_time - age < ( last_time - next_time ) * 2. / 3.
            qualifier = "middle "
        else
            qualifier = "late "
        end
        current_interval_name = qualifier * geo_interval_names[last_i]
    else
        current_interval_name = "preCambrian"
    end
    return current_interval_name
end

# data structures.  Restart Julia if you modify these
mutable struct rotation_struct
    id
    age
    latitudeaxis
    longitudeaxis
    angle
    parent
end
mutable struct plate_struct
    plateID         # scalar, from gplates rotation file
    crust_type       # 2d arrays
    crust_age
    crust_thickness
    crust_density
    sediment_thickness  # geologic record, dimensions of x,y,time_bin
    sediment_fractions  # x, y, sed type, time_bin
    rotationmatrix  # 3x3 heirarchy-resolved rotation matrix
    resolvetime
    parentstack     # list of parents at last resolve time
    lastparentstack
    firstappearance
    lastappearance
end
sediment_type_names = ["CaCO3","OrgC"]
n_sediment_types = length(sediment_type_names)
mutable struct orogenic_event_struct
    name
    footprint
    onset
    finish
    altitude
end

mutable struct world_struct
    # create and manipulate a global world for snapshot composite view of rotated plates
    age
    sealevel
    plateID         # pertaining to plate tectonics
    plateIDlist
    crust_type      # beginning of the 9 state variables
    crust_age
    crust_thickness # pertaining to geomorphology
    crust_density
    sediment_thickness
    sediment_fractions
    elevation_offset
    surface_elevation
    freeboard  # end of 9 state variables
    diags
end
world_diag_names = ["ocean_created_plate_area",
    "continent_2_ocean_plate_area",
    "continent_created_plate_area",
    "ocean_2_continent_plate_area",
    "ocean_subduct_plate_area",
    "ocean_subduct_age_plate_area",
    "ocean_subduct_sediment_plate_volume",
    "ocean_subduct_sediment_CaCO3_plate_volume",
    "ocean_subduct_sediment_OrgC_plate_volume",
    "continent_subduct_plate_area",
    "ocean_2_continent_world_area",
    "continent_2_ocean_world_area",
    "IDchange_plate_area",
    "continent_orogenic_uplift_rate",
    "subduction_orogenic_uplift_rate",
    "crust_weathering_rate",
    "aolean_erosion_rate",
    "seasurface_deposition_rate",  # used to couple weathering to seafloor depo
    "sediment_deposition_rate",
    "sediment_spilloff_rate",
    "sediment_redeposition_rate"]  # generates time record on the plates, sediment sink


function create_world( age )    # no plateIDs or numbers yet
    plateIDmap = fill(0.,nx,ny)
    plateIDlist = []
    crust_type = fill(notyetdefined,nx,ny)
    crust_age = fill(0.,nx,ny)
    sealevel = get_sealevel()
    crust_thickness = fill(0.,nx,ny) # pertaining to geomorphology
    crust_density = fill(0.,nx,ny)
    sediment_thickness = fill(0.,nx,ny)
    sediment_fractions = fill(0.,nx,ny,n_sediment_types)
    elevation_offset = fill(0.,nx,ny)
    surface_elevation = fill(0.,nx,ny)
    freeboard = fill(0.,nx,ny)
    n_diags = length(world_diag_names)
    diags = fill(0.,nx,ny,n_diags)
    world = world_struct(age,sealevel,plateIDmap,plateIDlist,crust_type,crust_age,
        crust_thickness,crust_density,sediment_thickness,sediment_fractions,
        elevation_offset,surface_elevation,freeboard,diags)
    world.plateID = read_plateIDs( age )
    world.plateIDlist = find_plateID_list( world.plateID )
    continentIDmap = read_continentIDs( age )
    for ix in 1:nx
        for iy in 1:ny
            if continentIDmap[ix,iy] > 0
                crust_type[ix,iy] = continent_type
                crust_thickness[ix,iy] = continent_crust_h0
                crust_density[ix,iy] = rho_continent_crust
            else
                crust_type[ix,iy] = ocean_type
                crust_thickness[ix,iy] = ocean_crust_h0
                crust_density[ix,iy] = rho_ocean_crust
            end
        end
    end
    return world
end
function create_blank_plate(plateID)
    crust_type = fill(notyetdefined,nx,ny)
    crust_age = fill(0.,nx,ny)
    crust_thickness = fill(0.,nx,ny)
    crust_density = fill(0.,nx,ny)
    sediment_thickness = fill(0.,nx,ny,n_sediment_layer_time_bins)
    sediment_fractions = fill(0.,nx,ny,n_sediment_layer_time_bins,length(sediment_type_names))
    #sediment_thickness[:,:,1] .= 10.
    #sediment_CaCO3_frac = fill(0.,nx,ny,n_sediment_layer_time_bins)
    rotationmatrix = fill(0,3,3)
    resolvetime = -1
    parentstack = []
    lastparentstack = []
    firstappearance = -1
    lastappearance = -1
    plate = plate_struct(plateID,crust_type,crust_age,
        crust_thickness,crust_density,
        sediment_thickness,sediment_fractions,
        rotationmatrix,resolvetime,
        parentstack,lastparentstack,
        firstappearance,lastappearance)
    return plate
end
function create_orogenic_event(name,onset,finish,altitude)
    worldfootprint = read_flip_csv("orogenies/" * name * ".csv")
    #plateIDmap = read_plateIDs(0.)
    #footprints = parse_world_footprint(worldfootprint,plateIDmap)
    orogenic_event = orogenic_event_struct(name,worldfootprint,onset,finish,altitude)
    return orogenic_event
end

function create_empty_plates()
    #    println(319,rotations[1])
    plates = Dict()
    for plateID in world.plateIDlist
        if haskey(plates,plateID) == false
            plate = create_blank_plate(plateID)
            plates[plateID] = plate
        end
    end
    return plates
end
function initialize_plates()
    resolve_rotation_matrices()
    for (plateID, plate) in plates
        #        println(413," ",plate.plateID)
        plate.firstappearance, plate.lastappearance = first_last_appearance(plate.plateID)
        plate.lastparentstack = deepcopy(plate.parentstack)
    end
    remask_plates()
    clear_world_process_arrays()
end
function fill_world_from_plates()
    # project plate crust* and sediment* variables to world grid.
    #Threads.@threads
    for iworld in 1:nx
        for jworld in 1:ny
            plateID = world.plateID[iworld,jworld]
            if plateID != -1 && haskey(plates,plateID)
                plate = plates[plateID]
                if plate.resolvetime != world.age && plateID != -1
                    println("updating rotation matrix in fillworld ", plateID)
                    resolve_rotation_matrix!(plateID)
                end
                iplate,jplate = nearest_plateij(plate.rotationmatrix,iworld,jworld)
                if plate.crust_type[iplate,jplate] == notatsurface # stupid edge case
                    # the closest plate point to the projection of the world point
                    # thinks it is outside of the exposed plate.
                    # dont fill it in because the projections are not reciprocal,
                    # dont map back.  the plate point would just toggle on and off.
                    world.crust_age[iworld,jworld] = 0.
                    world.crust_type[iworld,jworld] = ocean_type
                    world.crust_thickness[iworld,jworld] = ocean_crust_h0
                    world.crust_density[iworld,jworld] = rho_ocean_crust
                else
                    world.crust_age[iworld,jworld] = plate.crust_age[iplate,jplate]
                    world.crust_type[iworld,jworld] = plate.crust_type[iplate,jplate]
                    world.crust_thickness[iworld,jworld] =
                        plate.crust_thickness[iplate,jplate]
                    world.crust_density[iworld,jworld] =
                        plate.crust_density[iplate,jplate]
                    sediment_thickness = 0.
                    for ibin in 1:n_sediment_layer_time_bins
                        sediment_thickness +=
                            plate.sediment_thickness[iplate,jplate,ibin]
                    end
                    world.sediment_thickness[iworld,jworld] = sediment_thickness

                    ibin = top_filled_sediment_box(plate,iplate,jplate)
                    for itype in 1:n_sediment_types
                        world.sediment_fractions[iworld,jworld,itype] =
                            plate.sediment_fractions[iplate,jplate,ibin,itype]
                    end
                end
            end
        end
    end
    return
end
function fill_world_orogeny(footprint) # requires at least create_world(age)
    age = world.age
    plateIDmap = world.plateID
    world_orogeny = fill(0.,nx,ny)
    for ixworld in 1:nx
        for iyworld in 1:ny
            if world.crust_type[ixworld,iyworld] >= continent_type
                plateID = plateIDmap[ixworld,iyworld]
                mplate = resolve_rotation_matrix!(plateID,age)
                ixplate,iyplate = nearest_plateij(mplate,ixworld,iyworld)
                # the position on the rotated plate field
                # assume zero rotation at present-day, so look for footprint here
                if footprint[ixplate,iyplate] > 0.
                    world_orogeny[ixworld,iyworld] = footprint[ixplate,iyplate]
                end
            end
        end
    end
    return world_orogeny
end
function increment_plate_age()   # updates the age field in plate grids
    #Threads.@threads
    for plateID in world.plateIDlist
        plate = plates[plateID]
        for iplate in 1:nx
            for jplate in 1:ny
                if plate.crust_type[iplate,jplate] != notatsurface
                    plate.crust_age[iplate,jplate] += time_step
                end
            end
        end
    end
    return
end
function update_world_continents_from_file()    # continents
    world.plateID = read_plateIDs()
    continentID = read_continentIDs()
    for ix in 1:nx
        for iy in 1:ny
            new_crust_type = ocean_type
            if continentID[ix,iy] > 0
                new_crust_type = continent_type
            end
            if world.crust_type[ix,iy] != new_crust_type
                if new_crust_type >= continent_type
                    record_flux_into_world_diags("ocean_2_continent_world_area",
                        iy,ix,iy)
                    world.crust_thickness[ix,iy] = continent_crust_h0
                    world.crust_density[ix,iy] = rho_continent_crust
                else # new_crust_type == ocean_type
                    record_flux_into_world_diags("continent_2_ocean_world_area",
                        iy,ix,iy)
                    world.crust_thickness[ix,iy] = ocean_crust_h0
                    world.crust_density[ix,iy] = rho_ocean_crust
                    #record_continent_disappearing_from_world(ix,iy)
                end
                world.crust_type[ix,iy] = new_crust_type
            end
        end
    end
    return
end
function calculate_ocean_crust_depth(agefield) # from Stein and Stein 1992
    # is this done by ocean_thermal_boundary_layer?
    crustdepthfield = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if agefield[ix,iy] == continentagedefault
                crustdepthfield[ix,iy] = 0.
            else
                if agefield[ix,iy] < 20.
                    crustdepthfield[ix,iy] = -(2600. + 365. * sqrt(agefield[ix,iy]))
                else
                    crustdepthfield[ix,iy] = -(5651. - 2473. * exp( -0.0276 * agefield[ix,iy] ))
                end
            end
        end
    end
    return crustdepthfield
end
# geology routines
function get_subduction_orogeny()
    orogenic_footprint = fill(0.,nx,ny)
    subduction = get_diag("ocean_subduct_plate_area")
    for ix in 1:nx
        for iy in 1:ny
            if subduction[ix,iy] > 0.
                uplift = subduction[ix,iy] * 1e-7 # 3e-8
                coordlist = get_neighbor_coords(ix,iy,2,2)
                for coord in coordlist
                    if world.crust_type[coord[1],coord[2]] >= continent_type
                        orogenic_footprint[coord[1],coord[2]] = uplift
                    end
                end
            end
        end
    end
    return orogenic_footprint
end
function active_orogenic_events()
    orogenic_event_names = ""
    for (name,orogenic_event) in orogenic_events
        if world.age <= orogenic_event.onset && world.age >= orogenic_event.finish
            if length(orogenic_event_names) > 0
                orogenic_event_names = orogenic_event_names * "/"
            end
            orogenic_event_names = orogenic_event_names * name
        end
    end
    return orogenic_event_names
end

function smooth_orogeny(field)
    newfield = deepcopy(field)
    for ix in 1:nx
        for iy in 1:ny
            if field[ix,iy] > 0
                xys = get_neighbor_coords(ix,iy,2,2)
                for xy in xys
                    #println(xy[1]," ",xy[2])
                    if field[xy[1],xy[2]] == 0
                        newfield[xy[1],xy[2]] +=
                            (field[ix,iy] - newfield[xy[1],xy[2]]) * 0.2
                    end
                end
            end
        end
    end
    return newfield
end
function orogeny_intensity(uplift)
    total = 0.
    for ix in 1:nx
        for iy in 1:ny
            if uplift[ix,iy] > 0.
                total += uplift[ix,iy] * areabox[iy]
            end
        end
    end
    return total
end
function orogeny()
    orogeny_footprint = fill(0.,nx,ny)
    for (name,orogenic_event) in orogenic_events
        if world.age <= orogenic_event.onset && world.age >= orogenic_event.finish
            footprint = orogenic_event.footprint .* # 0 or 1
                orogenic_event.altitude*1e3 ./      # meters uplift
                ( orogenic_event.onset - orogenic_event.finish ) .*
                time_step
            rotatedfootprint = fill_world_orogeny(footprint)
            smoothed = smooth_orogeny(rotatedfootprint)
            #@printf("%s %.1e " , name,orogeny_intensity(rotatedfootprint))
            orogeny_footprint .+= rotatedfootprint
        end
    end
    #@printf("\n")
    #smooth_orogeny(orogeny_footprint)
    set_diag("continent_orogenic_uplift_rate",orogeny_footprint)
    subduction_orogeny = get_subduction_orogeny()
    set_diag("subduction_orogenic_uplift_rate",subduction_orogeny)
    world.crust_thickness .+= orogeny_footprint .+ subduction_orogeny
    isostacy()
    return
end
continent_crust_h0 = 16500.; continent_crust_h0_max = 35000.
ocean_crust_h0 = 4000.
rho_ocean_crust = 3.; rho_continent_crust = 2.5; rho_mantle = 3.5
rho_sediment = 2.; rho_seawater = 1.024
mantle_T0 = 2000. #?
ocean_T0 = 0.
function ocean_thermal_boundary_layer( )
    elevation_offset = ocean_thermal_boundary_layer(world.crust_age,
        world.crust_thickness,world.crust_density)
    return elevation_offset
end
function ocean_thermal_boundary_layer(crust_age,crust_thickness,crust_density) # sets world.elevation_offset etc for ocn
    # elevation_offset could also be used to tweak continents up and down
    # time independent, just based on crust age
    kappa = 1.e-6
    thermal_expansion = 3.e-5
    elevation_offset = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if world.crust_type[ix,iy] == ocean_type
                age = crust_age[ix,iy] * 1.e6  # years
                boundary_thickness = sqrt( kappa * age * 3.14e7) # meters, inc crust
                boundary_thickness = max(boundary_thickness,crust_thickness[ix,iy])
                mantle_boundary_thickness = boundary_thickness - crust_thickness[ix,iy]
                if mantle_boundary_thickness > 0.
                    temperature_crust_base = mantle_T0 -
                        ( mantle_T0 - ocean_T0 ) *
                        crust_thickness[ix,iy] / boundary_thickness
                else
                    temperature_crust_base = mantle_T0
                end
                mantle_bl_temp_avg = ( mantle_T0 + temperature_crust_base ) / 2.
                crust_temp_avg = ( ocean_T0 + temperature_crust_base ) / 2.
                crust_density[ix,iy] = rho_ocean_crust *
                    ( 1. - thermal_expansion * ( crust_temp_avg - mantle_T0 ) )
                mantle_bl_density = rho_mantle *
                    ( 1. - thermal_expansion * ( mantle_bl_temp_avg - mantle_T0 ) )
                elevation_offset[ix,iy] = boundary_thickness - # whole boundary layer
                    ( mantle_boundary_thickness * mantle_bl_density / rho_mantle +
                    crust_thickness[ix,iy] * crust_density[ix,iy] / rho_ocean_crust )
                #world.crust_thickness[ix,iy] = ocean_crust_h0
                #world.crust_density[ix,iy] = crust_density
            end
        end
    end
    return elevation_offset
end
function isostacy()
    (world.surface_elevation, world.freeboard) =
        isostacy( world.crust_thickness, world.crust_density,
        world.sediment_thickness, world.elevation_offset )
end
function isostacy(crust_thickness,crust_density,sediment_thickness,
    elevation_offset) # sets the elevation of the solid surface rel to mantle line
    # time independent
    surface_elevation = fill(0.,nx,ny)
    freeboard = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            surface_boundary_mass =
                crust_thickness[ix,iy] * crust_density[ix,iy] +
                sediment_thickness[ix,iy] * rho_sediment
            surface_boundary_thickness =
                crust_thickness[ix,iy] +
                sediment_thickness[ix,iy]
            equivalent_mantle_thickness =
                surface_boundary_mass /
                rho_mantle
            surface_elevation[ix,iy] = surface_boundary_thickness -
                equivalent_mantle_thickness + elevation_offset[ix,iy]
            freeboard[ix,iy] = surface_elevation[ix,iy] -
                world.sealevel
        end
    end
    return surface_elevation, freeboard
end

function get_sealevel() # will be replaced by interpolation like for rotations
    return 4700.
end
function compute_freeboard()
    for ix in 1:nx
        for iy in 1:ny
            world.freeboard[ix,iy] = world.surface_elevation[ix,iy] -
                world.sealevel
        end
    end
    return
end


# compute from 1000km and 10 Myr length and time scales
sediment_transport_coefficient = 0.3 # 1e6 * 1e6 / 1e8 # m2 / Myr
function i2dmat(ix,iy)
    return ix + (iy-1) * nx
end
function mat2ij(id)
    iy = Int(floor(id/nx))
    ix = id - iy
    return ix, iy
end
function setup_2d_transport()
    d = sediment_transport_coefficient * time_step
    n = nx * ny
    A = spzeros(n,n)
    for ix in 1:nx
        for iy in 1:ny # each grid point
            A[i2dmat(ix,iy),i2dmat(ix,iy)] = 1. + 2. * d # itself, horizontal
            # box to the left (or wrap)
            if ix > 1
                A[i2dmat(ix-1,iy),i2dmat(ix,iy)] = -d
                A[i2dmat(ix,iy),i2dmat(ix-1,iy)] = -d
            else
                A[i2dmat(nx,iy),i2dmat(ix,iy)] = -d
                A[i2dmat(ix,iy),i2dmat(nx,iy)] = -d
            end
            # box below
            if iy == 1
                A[i2dmat(ix,iy),i2dmat(ix,iy)] += d
            elseif iy == ny
                A[i2dmat(ix,iy),i2dmat(ix,iy)] += d
            else
                A[i2dmat(ix,iy),i2dmat(ix,iy)] += 2. * d
            end
            if iy > 1
                A[i2dmat(ix,iy-1),i2dmat(ix,iy)] = -d #* areabox[iy] / areabox[iy-1]
                A[i2dmat(ix,iy),i2dmat(ix,iy-1)] = -d #* areabox[iy] / areabox[iy-1]
            end
        end
    end
    return lu(A)
end
if(! @isdefined LUofAglobal )  # saves time redoing the decomp for recompiles
    LUofAglobal = setup_2d_transport()
end
function diffusion_elevation_change(elev_field)
    new_elev_field = reshape( LUofAglobal \ vec(elev_field), (nx,ny) )
    return new_elev_field .- elev_field
end

function get_continuous_areas(maskfield)
    unaccountedforfield = deepcopy(maskfield)
    blobs = []
    for ix in 1:nx
        for iy in 1:ny
            if unaccountedforfield[ix,iy] == 1 # found a blob
                blobfield = fill(0,nx,ny)
                check_point_for_blob(unaccountedforfield,blobfield,ix,iy)
                push!(blobs,blobfield)
            end
        end
    end
    return blobs
end

function check_point_for_blob(unaccountedforfield,blobfield,ix,iy)
    if unaccountedforfield[ix,iy] == 1
        #println(ix," ",iy)
        blobfield[ix,iy] = 1
        unaccountedforfield[ix,iy] = 0
        if ix == 1
            check_point_for_blob(unaccountedforfield,blobfield,nx,iy)
        else
            check_point_for_blob(unaccountedforfield,blobfield,ix-1,iy)
        end
        if ix == nx
            check_point_for_blob(unaccountedforfield,blobfield,1,iy)
        else
            check_point_for_blob(unaccountedforfield,blobfield,ix+1,iy)
        end
        if iy > 1
            check_point_for_blob(unaccountedforfield,blobfield,ix,iy-1)
        end
        if iy < ny
            check_point_for_blob(unaccountedforfield,blobfield,ix,iy+1)
        end
    end
    return
end

function get_mask_neighbors(maskfield)
    neighborfield = fill(0,nx,ny)
    for ix in 1:nx
        for iy in 2:ny-1
            if maskfield[ix,iy] == 1
                if ix == 1
                    ixleft = nx
                else
                    ixleft = ix - 1
                end
                if ix == nx
                    ixright = 1
                else
                    ixright = ix + 1
                end
                if maskfield[ixleft,iy] == 0
                    neighborfield[ixleft,iy] = 1
                end
                if maskfield[ixright,iy] == 0
                    neighborfield[ixright,iy] = 1
                end
                if maskfield[ix,iy-1] == 0
                    neighborfield[ix,iy-1] = 1
                end
                if maskfield[ix,iy+1] == 0
                    neighborfield[ix,iy+1] = 1
                end
            end
        end
    end
    return neighborfield
end

#function fill_fake_orogeny_sediment()
    #orogenic_mask = generate_mask_field(world.crust_type,uplifted_continent_type)
    #continuous_areas = get_continuous_areas(orogenic_mask)
    #for continuous_area in continuous_area
    #    orogenic_total_production_rate = volumefield_total(
    #        ( get_diag("continent_orogenic_uplift_rate") .+
    #        get_diag("subduction_orogenic_uplift_rate") ) .*
    #        continuous_area )

#        get_diag("continent_orogenic_uplift_rate") .+
#        get_diag("subduction_orogenic_uplift_rate"), 0.)
#
#        world.freeboard # fills in 0 for ocean points

function continental_erosion_transport()
    continent_elevation = greater_than_mask_field(world.crust_type,continent_type) .*
        world.freeboard # fills in 0 for ocean points

    #continent_elevation = fill_fake_orogeny_sediment(continent_elevation)
    change = diffusion_elevation_change(continent_elevation)
    for ix in 1:nx
        for iy in 1:ny
            if change[ix,iy] < 0.
                if world.sediment_thickness[ix,iy] > -change[ix,iy]
                    set_diag("sediment_deposition_rate",ix,iy,change[ix,iy])
                else
                    set_diag("crust_weathering_rate",ix,iy,change[ix,iy])
                end
            else
                if world.freeboard[ix,iy] > 0.
                    set_diag("sediment_deposition_rate",ix,iy,change[ix,iy])
                else
                    set_diag("seasurface_deposition_rate",ix,iy,change[ix,iy])
                end
            end
        end
    end
    return
end
function ocean_area()
    ocean_area = 0.
    for ix in 1:nx
        for iy in 1:ny
            if world.freeboard[ix,iy] < 0.
                ocean_area += areabox[iy] # m2
            end
        end
    end
    return ocean_area
end
aolean_erosion_time_constant = 30. # Myr
function aolean_ocean_deposition()
    for ix in 1:nx
        for iy in 1:ny
            if world.freeboard[ix,iy] > 0. &&
                world.sediment_thickness[ix,iy] > 0.1

                rate = world.freeboard[ix,iy] /
                    aolean_erosion_time_constant * time_step
                set_diag("aolean_erosion_rate",ix,iy,rate)
                accum_diag("sediment_deposition_rate",ix,iy,-rate)
                # meters
                #world.sediment_thickness[ix,iy] -= rate
            end
        end
    end

    total_aolean_erosion = volumefield_total(get_diag("aolean_erosion_rate")) # m3
    mean_aolean_deposition_rate = total_aolean_erosion / ocean_area()
    for ix in 1:nx
        for iy in 1:ny
            if world.freeboard[ix,iy] < 0.
                accum_diag("seasurface_deposition_rate",ix,iy,
                    mean_aolean_deposition_rate)
            end
        end
    end
    return
end
function seafloor_scour_fraction( freeboard )
    if freeboard > 0 # land
        return 1.
    elseif freeboard < -100.
        return 0.
    else
        scour_fraction = ( 100. + freeboard ) / 100.
    end
    return scour_fraction
end
function get_neighbor_coords(ix,iy,xbuffer,ybuffer)
    if ix > xbuffer
        ixleft = ix-xbuffer
    else
        ixleft = nx + ix - xbuffer
    end
    if ix + xbuffer <= nx
        ixright = ix + xbuffer
    else
        ixright = xbuffer
    end
    if iy > ybuffer
        iylower = iy - ybuffer
    else
        iylower = 1
    end
    if iy + ybuffer <= ny
        iyupper = iy + ybuffer
    else
        iyupper = ny
    end
    neighbor_coords = []
    for ix in ixleft:ixright
        for iy in iylower:iyupper
            push!(neighbor_coords,[ix,iy])
        end
    end
    return neighbor_coords
end
function minimum_neighbor(ix,iy,xbuffer,ybuffer)
    xys = get_neighbor_coords(ix,iy,xbuffer,ybuffer)
    min_freeboard = world.freeboard[ix,iy]
    minxy = [ix,iy]
    for xy in xys
        if world.freeboard[xy[1],xy[2]] < min_freeboard
            minxy = xy
        end
    end
    return minxy
end
function offshore_transport(mass,ixsource,iysource)
    accum_diag("sediment_spilloff_rate",ixsource,iysource,mass)
    ybuffer = 1; xbuffer = Int(floor(ybuffer / delta_x[iysource] * delta_y ))
    ixdest,iydest = minimum_neighbor(ixsource,iysource,
        xbuffer,ybuffer)
    while [ixdest,iydest] == [ixsource,iysource] && ybuffer < 4 # mid point is the minimum
        ybuffer += 1; xbuffer = Int(floor(ybuffer / delta_x[iysource] * delta_y ))
        ixdest,iydest = minimum_neighbor(ixsource,iysource,
            xbuffer,ybuffer)
    end
    new_potential_freeboard = world.freeboard[ixdest,iydest] + mass
    scour_frac = 0.
    if new_potential_freeboard > 0.
        local_deposition = world.freeboard[ixdest,iydest] / 10.
        scour_frac = ( 1. - local_deposition / mass )
        # overshoot governor
    else
        scour_frac = seafloor_scour_fraction(new_potential_freeboard)
    end
    if scour_frac > 0.
        overload = mass * scour_frac
        local_depostion = mass * ( 1. - scour_frac )
        accum_diag("sediment_deposition_rate",ixdest,iydest,local_depostion)
        accum_diag("sediment_redeposition_rate",ixdest,iydest,local_depostion)
        offshore_transport(overload,ixdest,iydest) # records the spilloff
    else
        accum_diag("sediment_deposition_rate",ixdest,iydest,mass)
    end
    return
end

function seafloor_deposition() # moves mass from seasurface_deposition_rate
    # to sediment_deposition_rate if there is space in the water column
    for ix in 1:nx
        for iy in 1:ny
            if world.freeboard[ix,iy] < 0. # ocean point or go home
                seasurface_deposition_rate =
                    get_diag("seasurface_deposition_rate",ix,iy)
                if seasurface_deposition_rate > 0.
                    new_potential_freeboard = world.freeboard[ix,iy] +
                        seasurface_deposition_rate
                    local_deposition = 0.
                    scour_frac = 0.
                    if new_potential_freeboard > 0.
                        local_deposition = world.freeboard[ix,iy] / 10.
                        scour_frac = local_deposition / seasurface_deposition_rate
                        # overshoot governor
                    else
                        scour_frac = seafloor_scour_fraction(new_potential_freeboard)
                    end

                    if scour_frac > 0.
                        overload = seasurface_deposition_rate * scour_frac
                        local_depostion = seasurface_deposition_rate *
                            ( 1. - scour_frac )
                        accum_diag("sediment_deposition_rate",ix,iy,local_depostion)
                        #println("overflowing ",ix," ",iy," ",overload)
                        #offshore_transport(overload,ix,iy)
                    else
                        accum_diag("sediment_deposition_rate",ix,iy,
                            seasurface_deposition_rate)
                    end
                end
            end
        end
    end
    return
end
function apply_fluxes_to_world()

    world.sediment_thickness .+= get_diag("sediment_deposition_rate") .-
        get_diag("aolean_erosion_rate")
    world.crust_thickness .+= get_diag("crust_weathering_rate")
    aolean = get_diag("aolean_erosion_rate")
    for ix in 1:nx
        for iy in 1:ny
            if world.sediment_thickness[ix,iy] > aolean[ix,iy]
                world.sediment_thickness[ix,iy] -= aolean[ix,iy]
            else
                world.crust_thickness[ix,iy] -= aolean[ix,iy]
            end
        end
    end
    return
end

function apply_crust_changes_to_plates()
    Threads.@threads for plateID in world.plateIDlist
        plate = plates[plateID]
        for iplate in 1:nx
            for jplate in 1:ny
                iworld, jworld = nearest_worldij(plate.rotationmatrix,iplate,jplate)
                if world.plateID[iworld,jworld] == plateID
                    #crust_weathering_rate = get_diag("crust_weathering_rate",iworld,jworld)
                    plates[plateID].crust_thickness[iplate,jplate] =
                        world.crust_thickness[iworld,jworld]
                    #if crust_weathering_rate != 0.
                    #    plates[plateID].crust_thickness[iplate,jplate] +=
                    #        crust_weathering_rate
                    #end
                end
            end
        end
    end
    return
end

function apply_sediment_fluxes_to_plates()
    ibin = search_sorted_below(sediment_layer_time_bins,world.age) # returns max value if older
    Threads.@threads for plateID in world.plateIDlist
        plate = plates[plateID]
        for iplate in 1:nx
            for jplate in 1:ny
                iworld, jworld = nearest_worldij(plate.rotationmatrix,iplate,jplate)
                if world.plateID[iworld,jworld] == plateID
                    sediment_deposition_rate =
                        get_diag("sediment_deposition_rate",iworld,jworld) -
                        get_diag("aolean_erosion_rate",iworld,jworld)
                    if plates[plateID].sediment_thickness[iplate,jplate,ibin] +
                        sediment_deposition_rate > 0. # depositional or plenty of mud to mine

                        plate.sediment_thickness[iplate,jplate,ibin] +=
                            sediment_deposition_rate  # meters, please
                    else # eroding too much for top layer
                        local_erosion_rate = - sediment_deposition_rate
                        local_erosion_rate -=
                            plates[plateID].sediment_thickness[iplate,jplate,ibin]
                        for ibinb in ibin+1:n_sediment_layer_time_bins
                            if local_erosion_rate > # depletes layer ibinb
                                plate.sediment_thickness[iplate,jplate,ibinb]

                                local_erosion_rate -=
                                    plates[plateID].sediment_thickness[iplate,jplate,ibinb]
                                plates[plateID].sediment_thickness[iplate,jplate,ibinb] = 0.
                            else # the layer is thick enough to kill it
                                plate.sediment_thickness[iplate,jplate,ibinb] -=
                                    local_erosion_rate
                                local_erosion_rate = 0. # kills further alteration downcore
                            end
                        end
                    end
                end
            end
        end
    end
    return
end

function global_plate_sediment_inventory()
    inventory = 0.
    for (plateID,plate) in plates
        for ix in 1:nx
            for iy in 1:ny
                for ibin in 1:n_sediment_layer_time_bins
                    inventory += plate.sediment_thickness[ix,iy,ibin] *
                        areabox[iy]
                end
            end
        end
    end
    return inventory
end
function global_world_sediment_inventory()
    inventory = 0.
    for ix in 1:nx
        for iy in 1:ny
            inventory += world.sediment_thickness[ix,iy] *
                areabox[iy]
        end
    end
    return inventory
end
function step_geomorph()
    world.age -= time_step
    clear_world_process_arrays()
    world.age = tectonic_time
    increment_plate_age()
    ocean_thermal_boundary_layer()
    orogeny()
    continental_erosion_transport()
    aolean_ocean_deposition()
    seafloor_deposition()
    apply_fluxes_to_world()
    isostacy( )
    return
end
function step_everything() # rebuilds the world at the new time
    world.age -= time_step
    clear_world_process_arrays() # battle stations
    #newIDs = read_plateIDs(newtime)
    resolve_rotation_matrices()
    # set rotation matrices to the new time
    increment_plate_age()
    # ocean and cont crust gets older in plate grids
    ocean_thermal_boundary_layer()
    # sets world.elevation_offset for aging ocean crust
    update_changing_plateIDs()
    # look for plates that change their identities, move stuff between them
    # requires old plateID list but new age in world.plateID
    # updates world.plateID
    fill_world_from_plates() # finish building the new world grid

    update_world_continents_from_file()
    # inflow of continent locations from gplates. initializes world crust_*
    # fields only where the world crust type has changed from last step
    # new continent points in world grid have h0 thickness for now.
    # sets the stage for remask_plates
    remask_plates()
    # finds points of crust subduction and creation within plates fields.
    # records plate subduction and creation in world diags.
    # initializes plate crust_* variables at new points in plate fields
    # new continent points on plate grid have h0 thickness

    orogeny()
    # applies orogenic events to change world.crust_thickness.
    # returns world in isostatic equilibrium

    continental_erosion_transport()
    # diffuses a sealevel-clamped elevation field to generate fields of crust
    # erosion, sediment deposition on land and seasurface deposition on ocean
    # points.

    aolean_ocean_deposition()
    # erodes continental soils at an increasing aolean_erosion_rate with elevation.
    # deposits into globally uniform seasurface_deposition_ratedeposition.
    # operates on world diag fields only
    #seafloor_deposition()
    # translate seasurface deposition into sediment_deposition_rate entries for
    # the ocean grid points.
    # recursive lateral downhill transport when accomodation space is lacking.

    apply_fluxes_to_world()
    # update world sediment_thickness and crust_thickness, just for show if looping

    apply_crust_changes_to_plates()
    apply_sediment_fluxes_to_plates()
    # updates the plate crust_thickness and sediment_thickness fields
    isostacy()
    # just for show if looping
    # clock out
    return
end
#function substep_everything()

# storage of diagnostics variables
function get_diag(varname)
    field = world.diags[:,:,findfirst(isequal(varname),world_diag_names)]
    return field
end
function get_diag(varname,ix,iy)
    value = world.diags[ix,iy,findfirst(isequal(varname),world_diag_names)]
    return value
end
function set_diag(varname,ix,iy,value)
    ibin = findfirst(isequal(varname),world_diag_names)
    world.diags[ix,iy,ibin] = value
    return world.diags[ix,iy,ibin]
end
function set_diag(varname,value)
    ibin = findfirst(isequal(varname),world_diag_names)
    world.diags[:,:,ibin] .= value
    return
end
function accum_diag(varname,ix,iy,value)
    ibin = findfirst(isequal(varname),world_diag_names)
    world.diags[ix,iy,ibin] += value
    return world.diags[ix,iy,ibin]
end
function world_diag_pos(diag_name)
    diag_pos = findfirst(isequal(diag_name),world_diag_names)
    return diag_pos
end
function record_flux_into_world_diags(fluxname,jsource,iworld,jworld,mult::Float64=1.)
    area = areabox[jsource]
    #println("accumulating ",fluxname," ", mult)
    accum_diag(fluxname,iworld,jworld,area * mult)
    return
end
function clear_world_process_arrays()
    world.diags .= 0.
    return
end

# generic utility routines
function search_sorted_nearest(a,x)
    minfound = abs(a[1]-x)
    iclosest = 1
    for i in 2:length(a)
        if abs(a[i]-x) < minfound
            iclosest = i
            minfound = abs(a[i]-x)
        end
    end
    return iclosest
end

function search_sorted_below(a,x)
    ijustbelow = 0
    for i in 1:length(a)
        if a[1] < a[end] # for example xcoords
            if a[i] <= x
                ijustbelow = i
            end
        else
            reverse_i = length(a) - i + 1
            if a[reverse_i] <= x
                ijustbelow = reverse_i
            end
        end
    end
    return ijustbelow
end

function first_last_appearance(plateID)
    filteredrotations = retrieve_rotations_for_plate(plateID)
    times = []
    for rotation in filteredrotations
        append!(times, rotation.age)
    end
    if length(times) > 0
        last = minimum(times)
        first = maximum(times)
        return first, last
    else
        return 0, 0
    end
end
# utilities operating on a 360x180 "field"
function write_field_csv(filename,field)
    f = open(filename,"w")
    for iy in 1:ny
        for ix in 1:nx-1
            @printf(f, "%0.2f,", field[ix,iy])
        end
        @printf(f,"%0.2f\n", field[nx,iy])
    end
    close(f)
    return
end
function read_my_csv(filename)
    f = open(filename)
    arr = fill(0,nx,ny)
    lines = readlines(f)
    for (iy,line) in enumerate(lines)
        iyflip = iy
        words = split(line,",")
        for (ix,word) in enumerate(words)
            value = parse(Int,word)
            arr[ix,iyflip] = value
        end
    end
    return arr
end
function read_flip_csv(filename)
    f = open(filename)
    arr = fill(0.,nx,ny)
    lines = readlines(f)
    for (iy,line) in enumerate(lines)
        iyflip = ny - iy + 1
        words = split(line,",")
        for (ix,word) in enumerate(words)
            value = parse(Int,word)
            arr[ix,iyflip] = value
        end
    end
    return arr
end
function read_plateIDs()
    return read_plateIDs( world.age )
end
function read_plateIDs( age )
    integerage = Int(ceil(age))
    filename = "platefiles/plateIDs." * string(integerage) * ".csv"
    plateIDmap = read_my_csv(filename)
    return plateIDmap
end
function read_continentIDs()
    return read_continentIDs( world.age )
end
function read_continentIDs( age )
    #    contIDmap = fill(0,nx,ny)
    #    return contIDmap
    integerage = Int(ceil(age))
    filename = "contfiles/contIDs." * string(Int(integerage)) * ".csv"
    contIDmap = read_my_csv(filename)
    return contIDmap
end
function find_plateID_list()
    return find_plateID_list( world.plateID )
end
function find_plateID_list(plateIDmap)
    plateIDlist = sort(unique(plateIDmap))
    killme = 0
    for i in 1:length(plateIDlist)
        if plateIDlist[i] == -1
            killme = i
        end
    end
    if killme > 0
        deleteat!(plateIDlist,killme)
    end
    return plateIDlist
end
function min_max_field(field)
    minval = field[1,1]; maxval = field[1,1]
    for ix in 1:nx
        for iy in 1:ny
            if field[ix,iy] < minval
                minval = field[ix,iy]
            end
            if field[ix,iy] > maxval
                maxval = field[ix,iy]
            end
        end
    end
    return minval,maxval
end
function highest_xy_field(field)
    maxval = field[1,1]
    max_ix = 1; max_iy = 1
    for ix in 1:nx
        for iy in 1:ny
            if field[ix,iy] > maxval
                maxval = field[ix,iy]
                max_ix = ix; max_iy = iy
            end
        end
    end
    return max_ix,max_iy
end
function clamp_field(field,value)
    newfield = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if field[ix,iy] > value
                newfield[ix,iy] = field[ix,iy]
            end
        end
    end
    return newfield
end
function generate_mask_field(field,value)
    newfield = fill(0.,nx,ny)
    #Threads.@threads
    for ix in 1:nx
        for iy in 1:ny
            if field[ix,iy] == value
                newfield[ix,iy] = 1
            end
        end
    end
    return newfield
end
function greater_than_mask_field(field,value)
    newfield = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if field[ix,iy] >= value
                newfield[ix,iy] = 1
            end
        end
    end
    return newfield
end
function limit_field(field,value)
    newfield = copy(field)
    for ix in 1:nx
        for iy in 1:ny
            if field[ix,iy] > value
                newfield[ix,iy] = value
            end
        end
    end
    return newfield
end
function multifield_mean_series(multifield,ages) # wtf
    series = fill(0.,length(ages))
    for i in 1:length(ages)
        series[i] = multifield[ages[i]]
    end
    return series
end
function get_continent_mask()
    continentmask = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if world.crust_type[ix,iy] >= continent_type
                continentmask[ix,iy] = 1.
            end
        end
    end
    return continentmask
end
function generate_outline_from_mask(maskfield)
    mycontour = Contour.contour(xcoords,ycoords,maskfield,0.5)
    lines = Contour.lines(mycontour)
    if length(lines) > 0
        line = [1]
        xs,ys = Contour.coordinates(line)
    else
        xs = []; ys = []
    end
    return xs,ys
end
function line_sphere2cart(xline,yline)  # x y, not lat longt
    nxl = length(xline)
    x3d = zeros(nxl);y3d = zeros(nxl); z3d = zeros(nxl)
    for i in 1:length(xline)
        azimuth = xline[i] / 180. * pi    # longitude
        elevation =  ( yline[i] + 90. ) / 180. * pi # 0 at pole, pi/2 eq, pi pole/ 180. * pi  # latitude
        x3d[i] = 1.01 * sin(elevation) * cos(azimuth) # convert(AbstractFloat,)
        y3d[i] = 1.01 * sin(elevation) * sin(azimuth)
        z3d[i] = 1.01 * - cos(elevation)
    end
    return x3d,y3d,z3d
end
function calculate_agebins(varfield,agefield,agebinvalues)
    # agebinvalues e.g. [0., 5., 10., ...]
    binnedvar = fill(0.,length(ages))
    for ix in 1:nx
        for iy in 1:ny
            for index in 2:length(ages)
                if agefield[ix,iy] < ages[index] &&
                agefield[ix,iy] > ages[index-1]
                    binnedvar[index-1] += areabox[iy] * varfield[ix,iy]
                end
            end
        end
    end
    return binnedvar
end
function calculate_area_per_age(agefield,agebinvalues)
    binnedvar = calculate_agebins(agefield,agefield,agebinvalues)
    return binnedvar
end
function areafield_total(areafield)
    areatot = 0.
    for ix in 1:nx
        for iy in 1:ny
            areatot += areafield[ix,iy]
        end
    end
    return areatot  # km2
end
function field_RMS(field)
    areatot = 0.
    valuetot = 0.
    for ix in 1:nx
        for iy in 1:ny
            areatot += areabox[iy]
            valuetot += areabox[iy] * field[ix,iy] * field[ix,iy]
        end
    end
    valueRMS = sqrt( valuetot / areatot )
    return valueRMS
end
function field_mean(field)
    areatot = 0.
    valuetot = 0.
    for ix in 1:nx
        for iy in 1:ny
            areatot += areabox[iy]
            valuetot += areabox[iy] * field[ix,iy]
        end
    end
    valueavg = valuetot / areatot
    return valueavg
end
function volumefield_total(meterfield) # used for meters of erosion to m3 globally
    total = 0.
    for ix in 1:nx
        for iy in 1:ny
            total += meterfield[ix,iy] * areabox[iy] # m3
        end
    end
    return total
end
function field_mean(field,maskfield)
    areatot = 0.
    valuetot = 0.
    for ix in 1:nx
        for iy in 1:ny
            if maskfield[ix,iy] == 1
                areatot += areabox[iy]
                valuetot += areabox[iy] * field[ix,iy]
            end
        end
    end
    valueavg = valuetot / areatot
    return valueavg
end
function field_area(field, value)
    areatot = 0.
    for ix in 1:nx
        for iy in 1:ny
            if field[ix,iy] == value
                areatot += areabox[iy]
            end
        end
    end
    return areatot/1.e6
end
function calculate_total_plate_area()
    totalarea = 0.
    for (age, plate) in plates
        for ix in 1:nx
            for iy in 1:ny
                if plate.crust_type[ix,iy] != notatsurface
                    totalarea += areabox[iy]
                end
            end
        end
    end
    return totalarea
end
# rotation routines
function sphere2cart(coords)   # lat, longt
    elevation = ( coords[1] + 90. ) / 180. * pi # 0 at pole, pi/2 eq, pi pole
    azimuth = coords[2] / 180. * pi
    x = sin(elevation) * cos(azimuth)
    y = sin(elevation) * sin(azimuth)
    z = - cos(elevation)  #
    v = [x,y,z]
    return v
end
function cart2sphere(v)   # [x,y,z]
    x = v[1]; y = v[2]; z = v[3]
    elevation = acos( -z )
    azimuth = atan(y/x)
    if x < 0
        if y < 0
            azimuth = azimuth-pi
        else
            azimuth = azimuth+pi
        end
    end
    latitude = elevation / pi * 180. - 90.
    longitude = azimuth / pi * 180.
    return [latitude, longitude]
end
function rotation2matrix(rotation)
    lataxis = rotation.latitudeaxis
    longtaxis = rotation.longitudeaxis
    angle = rotation.angle
    m = axisangle2matrix(lataxis,longtaxis,angle)
    return m
end
function axisangle2matrix(lataxis,longtaxis,angle)
    vrot = sphere2cart([lataxis,longtaxis])
    ux = vrot[1]
    uy = vrot[2]
    uz = vrot[3]
    theta = - angle / 180. * pi

    a11 = cos(theta)+ux^2*(1-cos(theta))
    a12 = ux*uy*(1-cos(theta))-uz*sin(theta)
    a13 = ux*uz*(1-cos(theta))+uy*sin(theta)

    a21 = uy*ux*(1-cos(theta))+uz*sin(theta)
    a22 = cos(theta)+uy^2*(1-cos(theta))
    a23 = uy*uz*(1-cos(theta))-ux*sin(theta)

    a31 = uz*ux*(1-cos(theta))-uy*sin(theta)
    a32 = uz*uy*(1-cos(theta))+ux*sin(theta)
    a33 = cos(theta)+uz^2*(1-cos(theta))

    m = [ a11 a12 a13; a21 a22 a23; a31 a32 a33 ]
    rotm = RotMatrix{3}(m)
    return m #   rotm
end
function matrix2axisangle(m)
    cart = []
    push!(cart, m[2,3] - m[3,2])
    push!(cart, m[3,1] - m[1,3])
    push!(cart, m[1,2] - m[2,1])
    trace = m[1,1] + m[2,2] + m[3,3]
    theta = acos((trace - 1)/2)
    magnitude = 2 * sin(theta)
    cart ./= magnitude
    spherevec = cart2sphere(cart)
    latitudeaxis = spherevec[1]
    longitudeaxis = spherevec[2]
    angle = theta / pi * 180.
    return latitudeaxis,longitudeaxis,angle
end
function matrix2rotation(m,id,age,idparent)
    latitudeaxis,longitudeaxis,angle = matrix2axisangle(m)
    rotation = rotation_struct(id,age,latitudeaxis,longitudeaxis,angle,idparent)
    return rotation
end
function read_rotation_file(filename)
    f = open(filename)
    lines = readlines(f)
    rotations = []
    for line in lines
        words = split(line)
        id = parse(Int,words[1])
        parent = parse(Int,words[6])
        v = parse.(Float64,words[2:5])
        rotation = rotation_struct(id,v[1],v[2],v[3],v[4],parent)
        push!(rotations,rotation)
    end
    sortedrotations = sort(rotations,by=x->(x.id, x.age))
    return sortedrotations# , comments
end
function retrieve_rotations_for_plate(plateID)
    rotationsforplate = []
    for rotation in rotations
        if rotation.id == plateID
            push!(rotationsforplate,rotation)
        end
    end
    if length(rotationsforplate) == 0 && plateID != -1
        #println("1327 no rotations found for ", plateID)
        rotation = rotation_struct(plateID,0,0,0,0,0)
        push!(rotationsforplate,rotation)
        rotation = rotation_struct(plateID,1000,0,0,0,0)
        push!(rotationsforplate,rotation)
    end
    return rotationsforplate
end
function slerp(qa, qb, t)
    qm = fill(0.,4)
    # Calculate angle between them.
    cosHalfTheta = qa.q.s * qb.q.s + qa.q.v1 * qb.q.v1 + qa.q.v2 * qb.q.v2 + qa.q.v3 * qb.q.v3

    if cosHalfTheta < 0 # it took the long path, so negate one end
        #println("long path")
        qa = one(QuatRotation) / qa
        cosHalfTheta = - qa.q.s * qb.q.s - qa.q.v1 * qb.q.v1 - qa.q.v2 * qb.q.v2 - qa.q.v3 * qb.q.v3
    end
    # if qa=qb or qa=-qb then theta = 0 and we can return qa
    if (abs(cosHalfTheta) >= 1.0)
        return qa
    end
    # Calculate temporary values.
    halfTheta = acos(cosHalfTheta)
    sinHalfTheta = sqrt(1.0 - cosHalfTheta*cosHalfTheta)
    # if theta = 180 degrees then result is not fully defined
    # we could rotate around any axis normal to qa or qb
    if (abs(sinHalfTheta) < 0.001)
        qm[1] = (qa.q.s * 0.5 + qb.q.s * 0.5)
        qm[2] = (qa.q.v1 * 0.5 + qb.q.v1 * 0.5)
        qm[3] = (qa.q.v2 * 0.5 + qb.q.v2 * 0.5)
        qm[4] = (qa.q.v3 * 0.5 + qb.q.v3 * 0.5)
        qret = QuatRotation(qm)
        return qret
    end
    ratioA = sin((1 - t) * halfTheta) / sinHalfTheta
    ratioB = sin(t * halfTheta) / sinHalfTheta
    # calculate Quaternion.
    qm[1] = (qa.q.s * ratioA + qb.q.s * ratioB)
    qm[2] = (qa.q.v1 * ratioA + qb.q.v1 * ratioB)
    qm[3] = (qa.q.v2 * ratioA + qb.q.v2 * ratioB)
    qm[4] = (qa.q.v3 * ratioA + qb.q.v3 * ratioB)
    qret = QuatRotation(qm)
    return qret
end
function crow_flies(lat1,lon1,lat2,lon2)
    r = 6371e3
    theta1 = lat1 * pi/180.
    theta2 = lat2 * pi/180.
    delta_theta = (lat2-lat2) * pi/180.
    delta_lambda = (lon2-lon1) * pi/180.
    a = sin(delta_theta/2.) * sin(delta_theta/2.) +
       cos(theta1) * cos(theta2) *
       sin(delta_lambda/2.) * sin(delta_lambda/2.)
    c = 2. * atan2(sqrt(a),sqrt(1-a))
    distance_meters = r * c
    return distance_meters
end
function interpolate_mplate(plateID,age)
    # computes the rotation matrix for the plate at the time,
    # and returns the parent ID
    #println(1569," ", plateID)
    timestepoffset = 0.001
    filteredrotations = retrieve_rotations_for_plate(plateID)
    #println(plateID)
    # make use of any direct hits
    mplates = Any[]; parents = Any[]
    for rotation in filteredrotations
        if rotation.age == age
            #            println(" got a hit ", plateID)
            mplate = rotation2matrix(rotation)
            parent = rotation.parent
            push!(mplates,copy(mplate))
            push!(parents,parent)
        end
    end
    if length(parents) > 1
        println(" got multiple hits ", plateID)
    end
    if length(parents) == 0
        # get here if there are no direct age hits
        if length(filteredrotations) == 0 && plateID != -1
            println("  No rotations found for ", plateID)
        end
        if age - timestepoffset < filteredrotations[1].age
            println("  Time too late for plate ", plateID)
            return [norotation],[0]
        end
        if age - timestepoffset > filteredrotations[end].age
            println("  Time too early for plate ", plateID)
            return [norotation],[0]
        end
        locallistofages = []
        for rotation in filteredrotations
            push!(locallistofages,rotation.age)
        end
        ibelow = search_sorted_below(locallistofages,age)
        fractions = fill(0.,length(locallistofages))
        if ibelow < length(locallistofages)
            if ibelow > 0       #  tbelow    time               t(below+1)
                fractions[ibelow+1] = (age - locallistofages[ibelow]) /
                    (locallistofages[ibelow+1] - locallistofages[ibelow])
                fractions[ibelow] = 1. - fractions[ibelow+1]
                if filteredrotations[ibelow].parent != filteredrotations[ibelow+1].parent
                    println("  Oops, switching parents ", plateID)
                else
                    mplate1 = rotation2matrix(filteredrotations[ibelow])
                    mplate2 = rotation2matrix(filteredrotations[ibelow+1])
                    qa = QuatRotation(mplate1)
                    qb = QuatRotation(mplate2)
                    t = fractions[ibelow+1]
                    #println(qa," ",qb)
                    qm = slerp(qa,qb,t)
                    mplate = RotMatrix(qm)

                    parent = filteredrotations[ibelow].parent
                    push!(mplates,copy(mplate))
                    push!(parents,parent)
                end
            else  # ibelow < 0
                println("how did I get here? ", ibelow, " ", locallistofages," ", plateID)
            end
        else # ibelow >= length
            if ibelow == length(locallistofages) &&
                locallistofages[ibelow] + timestepoffset >= age
                mplate = rotation2matrix(filteredrotations[ibelow])
                push!(mplates,copy(mplate))
                parent = filteredrotations[ibelow].parent
                push!(parents,parent)
                #                println(" using first rotation line", plateID)
            end
        end
    end
    return mplates, parents
end
function resolve_rotation_matrix!(plateID)
    return resolve_rotation_matrix!(plateID,world.age)
end
function resolve_rotation_matrix!(plateID,age)
    # recursively unwinds the heirarchy of parent plates
    # returns the rotation matrix, and records it in the plate structure
    # archives the parentstack into lastparentstack
    # detects rotation jumps at ancestry changes
    #    println("Resolving ", plateID," ",age)
    timestepoffset = 0.001
    hackedage = age + timestepoffset
    if plateID == -1
        return norotation
    end
    if haskey(plates,plateID) == false
        plates[plateID] = create_blank_plate(plateID)
    end
    if plates[plateID].resolvetime == age # nice round numbers for bookkeeping
        #        println(" solved ",plateID)
        return plates[plateID].rotationmatrix
    end
    mplates, idparentsbase = interpolate_mplate(plateID,hackedage)
    threadsfound = length(mplates)  # number of simultaneous rotations
    mplateoutlist = Any[]           # list of matrices, one per rotation
    if threadsfound == 1
        idparent = idparentsbase[1]
        mplate = mplates[1]
        #        println(295," starting ",idparent," ",plates[plateID].parentstack)
        plates[plateID].lastparentstack = deepcopy(plates[plateID].parentstack)
        plates[plateID].parentstack = []
        push!(plates[plateID].parentstack,idparent)
        if idparent != 0
            # gets here if it has a legit parent to deconvolve
            mparent = resolve_rotation_matrix!(idparent,age)
            mplate = mplate * mparent
            for ipar = plates[idparent].parentstack
                push!(plates[plateID].parentstack,ipar)
            end
        end
        if plates[plateID].parentstack != plates[plateID].lastparentstack &&
            plates[plateID].lastparentstack != [] &&
            plateID in world.plateID
            mplateold = plates[plateID].rotationmatrix
            mplatenew = mplate
    end
        push!(mplateoutlist,mplate)
    end
    if threadsfound > 2
        println(" more than 2 simultaneous rotations", plateID)
    end
    if threadsfound == 2
        println(" dup results ", plateID," ",idparentsbase)
        parentstackbythread = []
        for iresult in 1:threadsfound
        #        print(" -> ", idparent)
            mplate = mplates[iresult]
            idparent = idparentsbase[iresult]
            if idparent != 0
            # gets here if it has a legit parent to deconvolve
                mparent = resolve_rotation_matrix!(idparent,hackedage)
                # puts parent stack into globalparentstack
                mplate = mparent * mplate
            end
            push!(mplateoutlist, copy(mplate))
        end
    end
    if threadsfound > 2
        println(" yikes! more than two threads ", plateID)
    end
    if haskey(plates,plateID) == false
        plates[plateID] = create_blank_plate(plateID)
    end
    if length(mplateoutlist) == 0
        println(" not getting any threads ", plateID)
        push!(mplateoutlist,norotation)
    end
    plates[plateID].rotationmatrix = mplateoutlist[1]
    plates[plateID].resolvetime = age
    return mplateoutlist[1]
end
function generate_plate_velocity_field()
    platevelocityfield = fill(0.,nx,ny)
    for plateID in world.plateIDlist
        mplate1 = resolve_rotation_matrix!(plateID,world.age+1.)
        mplate2 = resolve_rotation_matrix!(plateID,world.age)
        for ix in 1:nx
            for iy in 1:ny
                if plateIDs[ix,iy] == plateID
                    lat2 = ycoords[iy]; lon2 = xcoords[ix]
                    vworldlater = [lat2,lon2]  # lat, longt at age
                    vplate = cart2sphere( mplate2 * vin )
                    vworldinitial = transpose(mplate1) * vplate
                    lat1 = vworldinitial[1]; lon1 = vworldinitial[2]
                    distance_meters = crow_flies(lat1,lon1,lat2,lon2)
                    platevelocityfield[ix,iy] = distance_meters / 1.e3 # mm / yr
                end
            end
        end
    end
    return platevelocityfield
end
function resolve_rotation_matrices()
    #    println("  Updating rotations to ", age)
    for plateID in world.plateIDlist
        if plates[plateID].resolvetime != world.age
            if plates[plateID].resolvetime != -1
                mplateold = plates[plateID].rotationmatrix
                mplatenew = resolve_rotation_matrix!(plateID)
                #distance = meanplatetraveldistance(mplateold,mplatenew)
            else
                resolve_rotation_matrix!(plateID)
            end
        end
    end
    return
end
function calculate_flow_velocity(iworld,jworld)
    plateID = world.plateID[iworld,jworld]
    latworldinitial = ycoords[jworld]; longtworldinitial = xcoords[iworld]
    vworldinitial = sphere2cart([latworldinitial,longtworldinitial])
    #    coordsworldinitial = cart2sphere(vworldinitial)
    mplate = resolve_rotation_matrix!(plateID,age)
    vplate = mplate * vworldinitial # 3d
    mplate2 = resolve_rotation_matrix!(plateID,age+0.1)
    vworldfinal = transpose(mplate2) * vplate
    coordsworldfinal = cart2sphere(vworldfinal)
    latworldfinal = vworldfinal[1]; longtworldfinal = vworldfinal[2]
    distance = crow_flies(latworldinitial,longtworldinitial,latworldfinal,longtworldfinal)
    velocity = distance * 10.  # meters per million years
    return velocity
end
function calculate_velocity_field()
    field = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            field[ix,iy] = calculate_flow_velocity(ix,iy)
        end
    end
    return field
end
# routines dealing with changing plate IDs
function copy_plate_point_plate12coord!(newplate,oldplate,ixold,iyold,ixnew,iynew)
    if oldplate.crust_type[ixold,iyold] == notatsurface #||
        #        oldplate.crust_type[ixold,iyold] == ocean_type  # moving ridge
        # edge cases, messy spreading
        newplate.crust_type[ixnew,iynew] = ocean_type
        newplate.crust_thickness[ixnew,iynew] = ocean_crust_h0
        newplate.crust_density[ixnew,iynew] = rho_ocean_crust
        newplate.crust_age[ixnew,iynew] = 0.
    else
        newplate.crust_type[ixnew,iynew] = oldplate.crust_type[ixold,iyold]
        newplate.crust_age[ixnew,iynew] = oldplate.crust_age[ixold,iyold]
        newplate.crust_thickness[ixnew,iynew] = oldplate.crust_thickness[ixold,iyold]
        newplate.crust_density[ixnew,iynew] = oldplate.crust_density[ixold,iyold]
        for ibin in 1:n_sediment_layer_time_bins
            newplate.sediment_thickness[ixnew,iynew,ibin] =
                oldplate.sediment_thickness[ixold,iyold,ibin]
            for itype in 1:n_sediment_types
                newplate.sediment_fractions[ixnew,iynew,ibin,itype] =
                    oldplate.sediment_fractions[ixold,iyold,ibin,itype]
            end
        end
    end
    return
end
function delete_plate_point!(oldplate,ixold,iyold)
    oldplate.crust_type[ixold,iyold] = notatsurface
    oldplate.crust_age[ixold,iyold] = 0.
    for ibin in 1:n_sediment_layer_time_bins
        oldplate.sediment_thickness[ixold,iyold,ibin] = 0.
        for itype in 1:n_sediment_types
            oldplate.sediment_fractions[ixold,iyold,ibin,itype] = 0.
        end
    end
    return
end
function update_two_plates!(oldplateID,newplateID)
        # called when one or more grid points change their plateIDs
    newtime = world.age; oldtime = newtime + time_step
    #println("updating ", oldplateID," ",newplateID," ",oldtime," ",newtime)
    if haskey(plates,oldplateID) == false
        println("creating plate to migrate data from ", oldplateID)
        plates[oldplateID] = create_blank_plate(newplateID)
    end
    oldplate = plates[oldplateID]
    if haskey(plates,newplateID) == false
        plates[newplateID] = create_blank_plate(newplateID)
        resolve_rotation_matrix!(newplateID,newtime)
    end
    newplate = plates[newplateID]
    if newplate.resolvetime < 0
        println("  Trying to update to an unrotated plate ", oldplate.plateID, " ", newplate.plateID)
    end
    oldIDmap = world.plateID
    newIDmap = read_plateIDs(newtime)
    oldIDmask = generate_mask_field(oldIDmap,oldplateID)
    # grid points in world grid of old plateID number
    newIDmask = generate_mask_field(newIDmap,newplateID)
    # points in world grid where the new plate ID is found
    needschangingworldmask = oldIDmask .* newIDmask
    # points in world grid where old ID changes to new ID
    firstoldtime,lastoldtime = first_last_appearance(oldplateID)
    firstnewtime,lastnewtime = first_last_appearance(newplateID)
    switchovertime = newtime # use rotations for newtime if possible
    if lastoldtime > newtime
        switchovertime = lastoldtime  # or back up to  end time of old plate
    end
    newplaterotation = resolve_rotation_matrix!(newplateID,switchovertime)
    needschangingnewplatemask = projected_plate_maskfield!(needschangingworldmask,
        newplaterotation)
    oldplaterotation = resolve_rotation_matrix!(oldplateID,switchovertime)
    needsdeletingoldplatemask = projected_plate_maskfield!(needschangingworldmask,
        oldplaterotation)
    for ixnewplate in 1:nx
        for iynewplate in 1:ny
            if needschangingnewplatemask[ixnewplate,iynewplate] > 0
                ixworld,iyworld = nearest_worldij(newplate.rotationmatrix,
                    ixnewplate,iynewplate)
                if ixworld > 0
                    ixoldplate,iyoldplate = nearest_plateij(oldplate.rotationmatrix,
                        ixworld,iyworld)
                    copy_plate_point_plate12coord!(newplate,oldplate,
                        ixoldplate,iyoldplate,ixnewplate,iynewplate)
                    record_flux_into_world_diags("IDchange_plate_area",iynewplate,
                        ixworld,iyworld)
                    #record_plateID_change_fillin(newplate,
                    #    ixworld,iyworld,
                    #    ixnewplate,iynewplate)
                end
            end
        end
    end

    for ixplate in 1:nx
        for iyplate in 1:ny
            if needsdeletingoldplatemask[ixplate,iyplate] > 0
                delete_plate_point!(oldplate,ixplate,iyplate)
                ixworld,iyworld = nearest_worldij(oldplate.rotationmatrix,
                    ixplate,iyplate)
                record_flux_into_world_diags("IDchange_plate_area",iyplate,
                    ixworld,iyworld,-1.)
                #record_plateID_change_delete(oldplate,
                #    ixworld,iyworld,
                #    ixplate,iyplate)
            end
        end
    end
    return
end
function changing_platelists(newIDmap)
    oldIDmap = world.plateID
    oldtime = world.age + time_step
    fp = 2
    changepairs = []
    for ixworld in 1+fp:nx-fp
        for iyworld in 1+fp:ny-fp
            if oldIDmap[ixworld,iyworld] !=
                newIDmap[ixworld,iyworld]
                foundone = 0
                for ifp in -fp:fp
                    if oldIDmap[ixworld+ifp,iyworld+ifp] ==
                        newIDmap[ixworld+ifp,iyworld+ifp]
                        foundone = 1
                    end
                end
                if foundone == 0 # truly a plate change not a migrating ridge
                    oldplateID = oldIDmap[ixworld,iyworld]
                    newplateID = newIDmap[ixworld,iyworld]
                    key = [ oldplateID, newplateID ]
                    # encode both IDs into single key
                    #key = oldplateID * 10000000 + newplateID
                    push!(changepairs,key)
                end
            end
        end
    end
    # list of all gridpoints in successive world grids where the ID changes
    #println("key list", changepairs)
    changepairs = unique(changepairs)
    #println("key list", changepairs)
    return changepairs
end
function update_changing_plateIDs()
    newtime = world.age
    newIDmap = read_plateIDs(world.age)
    changingplates = changing_platelists(newIDmap)
    #println("changing ", changingplates)
    Threads.@threads for platepair in changingplates
        oldplateID = platepair[1]; newplateID = platepair[2]
        update_two_plates!(oldplateID,newplateID)
        #        end
    end
    world.plateID = newIDmap
    world.plateIDlist = find_plateID_list()
    return
end
# grid navigation utilities
function nearest_plateij(mplate,iworld,jworld)  # ray trace from a world i,j to a plate
    longitudeworld = xcoords[iworld]
    latitudeworld = ycoords[jworld]
    if latitudeworld > 90.
        latitudeworld = 180. - latitudeworld
        longitudeworld += 180.
    end
    if latitudeworld < -90.
        latitudeworld = - 180. + latitudeworld
        longitudeworld += 180.
    end
    while longitudeworld < -180.; longitudeworld += 360.; end
    while longitudeworld >  180.; longitudeworld -= 360.; end
    vworld = sphere2cart([latitudeworld,longitudeworld])
    vplate = mplate * vworld  # rotate into plate frame
    mapcoords = cart2sphere(vplate)
    iplate = search_sorted_nearest(xcoords,mapcoords[2])
    jplate = search_sorted_nearest(ycoords,mapcoords[1])
    return iplate, jplate
end
function nearest_worldij(rotationmatrix,iplate,jplate)
    latitudepoint = ycoords[jplate]
    longitudepoint = xcoords[iplate]
    vpoint = sphere2cart([latitudepoint,longitudepoint])
    newvpoint = transpose(rotationmatrix) * vpoint
    # unrotating back to world grid
    mapcoords = cart2sphere(newvpoint)
    iworld = search_sorted_nearest(xcoords,mapcoords[2])
    jworld = search_sorted_nearest(ycoords,mapcoords[1])
    return iworld, jworld
end
function nearest_worldij(plate,rotation,ifield,jfield)
    m = rotation2matrix(rotation)
    iworld,jworld = nearest_worldij(m,ifield,jfield)
    return iworld,jworld
end
function nearest_other_plate_ij(iplate1,jplate1,mplate1,mplate2)
    plate1v = sphere2cart( [ ycoords[jplate1], xcoords[iplate1] ] )
    worldv = transpose(mplate1) * plate1v
    plate2v = mplate2 * worldv
    plate2coords = cart2sphere(plate2v) # [lat, longt]
    iplate2 = search_sorted_nearest(xcoords,plate2coords[2])
    jplate2 = search_sorted_nearest(ycoords,plate2coords[1])
    return iplate2,jplate2
end

function projected_plate_maskfield!(worldmaskfield,rotationmatrix)
    # find all points in plate that project to world region
    newfield = fill(0,nx,ny)
    #Threads.@threads
    for ixplate in 1:nx
        #Threads.@threads
        for iyplate in 1:ny
            ixworld,iyworld = nearest_worldij(rotationmatrix,
                ixplate,iyplate)
            if ixworld > 0
                if worldmaskfield[ixworld,iyworld] == 1
                    newfield[ixplate,iyplate] = 1
                end
            end
        end
    end
    return newfield
end
function projected_world_maskfield!(platemaskfield,rotationmatrix)
    # find all points in plate that project to world region
    worldfield = fill(0,nx,ny)
    for ixworld in 1:nx
        #Threads.@threads
        for iyworld in 1:ny

            ixplate,iyplate = nearest_plateij(rotationmatrix,
                ixworld,iyworld)

            if platemaskfield[ixplate,iyplate] == 1
                worldfield[ixworld,iyworld] = 1
            end
        end
    end
    return worldfield
end
function remask_plate!(plate)
    # fill in plate from world grid.
    # finds grid points of subduction and crust creation on plate fields.
    # imports crust_thickness as it stands in world grid.
    # when a plate point starts corresponding to a continent point in the world
    #  it gets the world crust_thickness etc, before orogeny due to the current
    #  step.
    age = world.age
    if plate.resolvetime != age
        resolve_rotation_matrix!(plate.plateID,age)
    end
    #Threads.@threads
    for iplate in 1:nx
        for jplate in 1:ny
            iworld, jworld = nearest_worldij(plate.rotationmatrix,iplate,jplate)
            if plate.crust_type[iplate,jplate] == notyetdefined # initialize
                if world.plateID[iworld,jworld] == plate.plateID
                    plate.crust_type[iplate,jplate] = world.crust_type[iworld,jworld]
                    plate.crust_age[iplate,jplate] = 0.
                    plate.crust_thickness[iplate,jplate] =
                        world.crust_thickness[iworld,jworld]
                    if plate.crust_type[iplate,jplate] == ocean_type
                        plate.crust_density[iplate,jplate] = rho_ocean_crust
                    else
                        plate.crust_density[iplate,jplate] = rho_continent_crust
                    end
                else #
                    plate.crust_type[iplate,jplate] = notatsurface
                end
            elseif world.plateID[iworld,jworld] == plate.plateID
                # our plate outcrops here in the world grid
                if world.crust_type[iworld,jworld] == ocean_type
                    if plate.crust_type[iplate,jplate] == notatsurface
                        # ridge spreading hopefully
                        record_flux_into_world_diags("ocean_created_plate_area",
                            jplate,iworld,jworld)
                        plate.crust_type[iplate,jplate] = ocean_type
                        plate.crust_age[iplate,jplate] = 0. # new ocean crust
                        plate.crust_thickness[iplate,jplate] = ocean_crust_h0
                        plate.crust_density[iplate,jplate] = rho_ocean_crust
                    elseif plate.crust_type[iplate,jplate] >= continent_type
                        record_flux_into_world_diags("continent_2_ocean_plate_area",
                            jplate,iworld,jworld)
                        #record_continent_2_ocean_from_plates(jplate,iworld,jworld)
                        #println("disappearing continent ",plate.plateID," ",
                        #    iplate," ",jplate)
                        plate.crust_type[iplate,jplate] = ocean_type
                        plate.crust_age[iplate,jplate] = 0. # new ocean crust
                        plate.crust_thickness[iplate,jplate] = ocean_crust_h0
                        plate.crust_density[iplate,jplate] = rho_ocean_crust
                    end # or else plate was already ocean, do nothing
                else # world.crust_type[iworld,jworld] == continent_type
                    if plate.crust_type[iplate,jplate] == notatsurface
                        # continent has formed where nothing was
                        record_flux_into_world_diags("continent_created_plate_area",
                            jplate,iworld,jworld)
                        #record_continent_creation_from_plates(jplate,iworld,jworld)
                    elseif plate.crust_type[iplate,jplate] == ocean_type
                        # ocean type has transformed into continental type
                        record_flux_into_world_diags("ocean_2_continent_plate_area",
                            jplate,iworld,jworld)
                        #record_ocean_2_continent_from_plates(jplate,iworld,jworld)
                    end
                    if plate.crust_type[iplate,jplate] < continent_type
                        plate.crust_type[iplate,jplate] = continent_type
                        plate.crust_age[iplate,jplate] = 0.
                        plate.crust_thickness[iplate,jplate] =
                            continent_crust_h0
                        plate.crust_density[iplate,jplate] = rho_continent_crust
                    end
                    if plate.crust_type[iplate,jplate] >= continent_type
                        plate.crust_type[iplate,jplate] =
                            min(continent_crust_h0,plate.crust_type[iplate,jplate])
                    end
                end # updated the status of the plate grid point if required
            else       # world map says we dont exist here
                if plate.crust_type[iplate,jplate] == ocean_type
                    # but we did last time step
                    record_flux_into_world_diags("ocean_subduct_plate_area",
                        jplate,iworld,jworld)
                    record_flux_into_world_diags("ocean_subduct_age_plate_area",
                        jplate,iworld,jworld,plate.crust_age[iplate,jplate])
                    for ibin in 1:n_sediment_layer_time_bins
                        record_flux_into_world_diags("ocean_subduct_sediment_plate_volume",
                            jplate,iworld,jworld,
                            plate.sediment_thickness[iplate,jplate,ibin])
                        for itype in 1:n_sediment_types
                            fluxname = "ocean_subduct_sediment_" *
                                sediment_type_names[itype] *
                                "_plate_volume"
                            flux_amt = plate.sediment_thickness[iplate,jplate,ibin] *
                                plate.sediment_fractions[iplate,jplate,ibin,itype]
                            record_flux_into_world_diags(fluxname,
                                jplate,iworld,jworld,flux_amt)
                        end
                    end
                    #record_ocean_disappearing_from_plates(plate,iplate,jplate,iworld,jworld)
                elseif plate.crust_type[iplate,jplate] >= continent_type
                    record_flux_into_world_diags("continent_subduct_plate_area",
                        jplate,iworld,jworld)
                    #record_continent_disappearing_from_plates(plate,iplate,jplate,iworld,jworld)
                end
                if plate.crust_type[iplate,jplate] != notatsurface
                    delete_plate_point!(plate,iplate,jplate)
                end
            end # initialize or update
        end # loop over plate grid
    end
    return
end
function remask_plates()   # subduction and crust creation
    Threads.@threads for plateID in world.plateIDlist     #[25:25]# 1:length(plates)
        if haskey(plates,plateID) == false
            plates[plateID] = create_blank_plate(plateID)
        end
        plate = plates[plateID]
        #        print(plate.plateID,"\n")
        #        if age <= plate.firsttime
        #            if age >= plate.lasttime
        remask_plate!(plate)
        #println(plateID," ", areafield_total(world.subducted_area))
        #                   end
        #        end
    end
    for (plateID,plate) in plates
        if !(plateID in world.plateIDlist)
            delete!(plates,plateID)
        end
    end
    return
end
function top_filled_sediment_box(plate,iplate,jplate)
    ibinmax = search_sorted_below(sediment_layer_time_bins,world.age)
    ibinmax = max(ibinmax,1)
    ibintop = 1
    for ibin in 2:ibinmax
        if plate.sediment_thickness[iplate,jplate,ibin] >= 0.
            ibintop = ibin
        end
    end
    return ibintop
end

# Makie map plots
function zoomfield(field,ix,iy,nskirt)
    ixleft = max(1,ix-nskirt); ixright = min(nx,ix+nskirt)
    iybelow = max(1,iy-nskirt); iyabove = min(ny,iy+nskirt)
    subfield = field[ixleft:ixright,iybelow:iyabove]
    subx = xcoords[ixleft:ixright]
    suby = ycoords[iybelow:iyabove]
    scene = Scene(resolution = (1000, 800))
    heatmap!(scene,subx,suby,subfield,colormap=:thermometer)
    return scene
end
nxcolorbar = 5
function setup_plotfield()
    scene = Scene(resolution = (1000, 550))  # GLMakie.lines([-180,180],[0,0],show_axis=false)
    setup_plotfield!(scene,0,0)
    return scene
end
function setup_double_plotfield()
    scene = Scene(resolution = (1000, 1000))  # GLMakie.lines([-180,180],[0,0],show_axis=false)
    setup_plotfield!(scene,0,0)
    setup_add_offset_field!(scene,0,-200)
    return scene
end

function setup_add_offset_field!(scene,xoffset,yoffset)
    setup_plotfield!(scene,xoffset,yoffset)
    return scene
end
function setup_plotfield!(scene,xoffset,yoffset)
    for i in [-180,-120,-60,0,60,120,180,180+nxcolorbar]
        x = [i,i] .+ xoffset; y = [-90,90] .+ yoffset
        lines!(scene,x,y,show_axis=false)
    end
    for j in [-90,-60,-30,0,30,60,90]
        x = [-180,180+nxcolorbar] .+ xoffset; y = [j,j] .+ yoffset
        lines!(scene,x,y)
    end
    for i in [-200,220]
        x = [i,i] .+ xoffset; y = [-110,110] .+ yoffset
        lines!(scene,x,y,color=:white)
    end
    for j in [-110,110]
        x = [-200,200] .+ xoffset; y = [j,j] .+ yoffset
        lines!(scene,x,y,color=:white)
    end
    #    println("     range in plot ",min," ",max)
    return scene
end
function clear_plotfield!(scene)
    clear_plotfield!(scene,20)
    return scene
end
function clear_plotfield!(scene,ibase)
    for i in ibase+1:length(scene)
        delete!(scene,scene[end])
    end
    return scene
end
function add_colorbar_to_field(field,minval,maxval)
    xplotcoords = fill(0.,nx+nxcolorbar)
    for i in 1:length(xcoords)
        xplotcoords[i] = xcoords[i]
    end
    for i in nx+1:nx+nxcolorbar
        xplotcoords[i] = xplotcoords[i-1] + 1.
        #println(1191,i," ",xplotcoords[i])
    end
    colorfield = fill(0.,nx+nxcolorbar,ny)
    colorinterval = ( maxval - minval ) / ny
    #println(minval," ",maxval," ",colorinterval)
    for ix in 1:nx
        for iy in 1:ny
            colorfield[ix,iy] = field[ix,iy]
        end
    end
    for ix in nx+1:nx+nxcolorbar
        for iy in 1:ny
            colorfield[ix,iy] = minval + (iy-1) * colorinterval
        end
    end
    return colorfield,xplotcoords
end
function plotfield_with_colortab!(scene,field,minval,maxval,style::String="rainbow")
    plotfield_with_colortab!(scene,field,minval,maxval,0,0,style)
    return scene
end
function plotfield_with_colortab!(scene,field,minval,maxval,xoffset,yoffset,style::String="rainbow")
    #nxcolorbar = 5
    colorfield, xplotcoords = add_colorbar_to_field(field,minval,maxval)
    colorinterval = ( maxval - minval ) / 6
    for ilabel = 1:7
        labelval = Int(floor(minval + (ilabel-1) * colorinterval))
        labelval = string(labelval)
        labelpos = Int(-95. + (ilabel-1) * 30.)
        text!(scene,labelval,position=(190 + xoffset,labelpos + yoffset),textsize=10)
        #println("printing ", labelval)
    end
    #text!(scene,string(Int(floor(maxval))),position=(190,85),textsize=10)
    cmap = Reverse(:lightrainbow)
    if style == "freeboard"
        cmap = Reverse(:gist_earth)
    end
    heatmap!(scene,xplotcoords .+ xoffset,ycoords .+ yoffset,colorfield,
        colorrange=(maxval,minval),colormap=cmap)#thermometer)#lightrainbow)
    return scene
end
function plot_add_timestamp!(scene)
    scene = plot_add_timestamp!(scene,world.age,-180,-110)
    return scene
end
function plot_add_timestamp!(scene,age,xpos,ypos)
    timestamp = string(Int(floor(world.age))) * " Myr, " * get_geo_interval(age)
    text!(scene, timestamp, position = (xpos,ypos),textsize=15)
    return scene
end
function plot_add_plate_boundaries!(scene)
    for plateID in world.plateIDlist
        maskfield = generate_mask_field(world.plateID,plateID)
        contour!(scene,xcoords,ycoords,maskfield,color=:white)
    end
    return scene
end
function generate_flowpath_line(lat, longt, plateID)
    # velocity arrow originating at lat longt in world grid
    vworldinitial = sphere2cart([lat,longt])  # x,y,z cartesian coords
    mplate = resolve_rotation_matrix!(plateID,world.age)
    vplate = mplate * vworldinitial # still xyz
    npoints = 2
    linex = [ convert(AbstractFloat,longt) ]
    liney = [ convert(AbstractFloat,lat) ]
    for i in 1:npoints-1  # 0
        ageval = world.age - i * 5. # 5 million years of travel
        mplate = resolve_rotation_matrix!(plateID,ageval)
        vworld = transpose(mplate) * vplate
        coordsworld = cart2sphere(vworld)
        latworldnow = coordsworld[1]; longtworldnow = coordsworld[2]
        push!(linex,longtworldnow); push!(liney,latworldnow)
    end
    return linex,liney
end
function plot_add_streamlines!(scene) # sprinkle flowpathlines over map
    lats = range(-80,stop=80,length=9)
    longts = range(-170,stop=170,length=18)
    for lat in lats
        for longt in longts
            iworld = search_sorted_nearest(xcoords,longt)
            jworld = search_sorted_nearest(ycoords,lat)
            originx = [longt-1,longt,longt+1,longt,longt-1]
            originy = [lat,lat+1,lat,lat-1,lat]
            lines!(scene,originx,originy,color=:black)
            plateID = world.plateID[iworld,jworld]
            linex,liney = generate_flowpath_line(lat, longt, plateID)
            lines!(scene,linex,liney,color=:black)
        end
    end
    return scene
end
function plot_add_continent_outlines!(scene)
    continentmask = get_continent_mask()
    contour!(scene,xcoords,ycoords,continentmask,color=:black)
    return scene
end
function plot_add_orogenies!(scene)
    continent_uplifting = greater_than_mask_field(get_diag("continent_orogenic_uplift_rate"),0.01)
    contour!(scene,xcoords,ycoords,continent_uplifting,color=:blue)
    subduction_uplifting = greater_than_mask_field(get_diag("subduction_orogenic_uplift_rate"),0.01)
    contour!(scene,xcoords,ycoords,subduction_uplifting,color=:green)
    orogenies = active_orogenic_events()
    text!(scene, orogenies, position = (0,-105),textsize=15)
    return scene
end
function plotfield(field)  # for a generic field, auto scaling
    scene = setup_plotfield()
    minval,maxval = min_max_field(field)
    plotfield_with_colortab!(scene,field,minval,maxval)
    return scene
end
function plotfield(field,minval,maxval,style::String="rainbow")  # intended for generic field
    scene = setup_plotfield()
    plotfield_with_colortab!(scene,field,minval,maxval,style)
    return scene
end
function plot_two_fields(field1,field2,minval,maxval,style::String="rainbow")
    scene = setup_double_plotfield()
    #setup_add_offset_field!(scene,0,-200)
    plotfield_with_colortab!(scene,field1,minval,maxval,style)
    plotfield_with_colortab!(scene,field2,minval,maxval,0,-200,style)
    return scene
end
# depth surface plot
function setup_depth_surfaceplot()
    scene = Scene()
    for i in [-180,-120,-60,0,60,120,180,180+nxcolorbar]
        x = [i,i]; y = [-90,90]; z = [0,0]
        lines!(scene,x,y,z,show_axis=false)
    end
    for j in [-90,-60,-30,0,30,60,90]
        x = [-180,180+nxcolorbar]; y = [j,j];z = [0,0]
        lines!(scene,x,y,z)
    end
    eyeposition = [0,-0.1,.1]
    #cam = AbstractPlotting.cam3d_cad!(scene, eyeposition=eyeposition, fov=30f0)
    cam = cam3d_cad!(scene, eyeposition=eyeposition, fov=30f0)
    return scene
end
function color_depth_surfaceplot!(scene,depthfield,colorfield,minval,maxval)
    scaleddepth = depthfield * 0.003
    newcolorfield, xplotcoords = add_colorbar_to_field(colorfield,minval,maxval)
    newscaleddepth, xplotcoords = add_colorbar_to_field(scaleddepth,0.,0.)

    surface!(scene,xplotcoords,ycoords,newscaleddepth,color=newcolorfield,
    colormap=:lightrainbow,shading=false)
    #lightposition = Vec3f0(100, 0, -15), ambient = Vec3f0(1.,1.,1.)
    colorinterval = ( maxval - minval ) / 6
    for ilabel = 1:7
        labelval = Int(floor(minval + (ilabel-1) * colorinterval))
        labelval = string(labelval)
        labelpos = Int(-95. + (ilabel-1) * 30.)
        text!(scene,labelval,position=(190,labelpos),textsize=10)
        #println("printing ", labelval)
    end
    timestamp = string(Int(floor(age))) * " Myr"
    text!(scene, timestamp, position = (-180.,-110.),textsize=15)
    plotfieldboundaries!(scene)
    plot_add_streamlines!(scene)
    plot_add_continent_outlines!(scene)
    return scene
end
function depthageplot()
    scene = setup_depth_surfaceplot()
    depthfield = world.freeboard
    #depthfield = fill(0.,nx,ny)
    minval,maxval = min_max_field(world.crust_age)
    if maxval > 200.
        maxval = 200.
    end
    if maxval < 10.
        maxval = 10.
    end
    minval = 0.
    color_depth_surfaceplot!(scene,world.freeboard,world.crust_age,
        minval,maxval)

    return scene
end
# spherical projection plots
function setup_sphereplot(lat,longt) # setup grid for surface plots
    scene = Scene()
    global sphereeyeposition = [lat,longt]
    u = range(-pi,stop=π,length=nx)
    v = range(0,stop=π,length=ny)
    global spheregridx = zeros(nx,ny)
    global spheregridy = zeros(nx,ny)
    global spheregridz = zeros(nx,ny)
    for i in 1:nx
        for j in 1:ny
            spheregridx[i,j] = cos.(u[i]) * sin(v[j]);
            spheregridy[i,j] = sin.(u[i]) * sin(v[j]);
            spheregridz[i,j] = cos(v[j]);
        end
    end
    #dummy_grid = fill(0.,nx,ny)
    # plot the poles
    pole_x = [0.e0,0.e0];pole_y = copy(pole_x); pole_z = [-1.5e0,0.e0]
    lines!(scene,pole_x,pole_y,pole_z,color=:green,show_axis=false) # South pole green
    pole_z = [0.e0,1.5e0];
    lines!(scene,pole_x,pole_y,pole_z,color=:red)  # North is red
    # longitude lines
    ugrid = range(-180,stop=180,length=13)
    vgrid = range(-90,stop=90,length=ny)
    for i in 1:12
        linex = [];liney = []
        for j in 1:ny
            append!(linex,ugrid[i]) # given arc all same longt
            append!(liney,vgrid[j])
        end
        arcx,arcy,arcz = line_sphere2cart(linex,liney)
        if i == 1
            lines!(scene,arcx,arcy,arcz,color=:red)  # +- 180
        elseif i == 7
            lines!(scene,arcx,arcy,arcz,color=:green)
        else
            lines!(scene,arcx,arcy,arcz)
        end
    end
    # latitude lines
    vgrid = range(-90,stop=90,length=9)
    for j in 2:8
        linex = []; liney = []
        for i in 1:nx
            push!(linex,u[i]/pi*180.)
            push!(liney,vgrid[j])
        end
        neweqx,neweqy,neweqz = line_sphere2cart(linex,liney)
        if j == 5
            lines!(scene,neweqx,neweqy,neweqz,color=:red)
        else
            lines!(scene,neweqx,neweqy,neweqz)
        end
    end
    adjust_eyeposition_sphereplot!(scene,lat,longt)
    return scene
end
function clear_sphereplot!(scene)
    for i in 22:length(scene)
        delete!(scene,scene[end])
    end
    return scene
end
function adjust_eyeposition_sphereplot!(scene,lat,longt)
    # set viewpoint
    global sphereeyeposition = [lat,longt]
    longtrad = longt / 180. * pi
    latrad = lat / 180. * pi
    scale = 1.5
    eyex = scale * cos(longtrad) * cos(latrad)
    eyey = scale * sin(longtrad) * cos(latrad)
    eyez = scale * sin(latrad)
    eyeposition = Float32[eyex,eyey,eyez]
    println("eyeposition now ", eyeposition)
    #    cam = AbstractPlotting.cam3d_cad!(scene, eyeposition=eyeposition)
    cam = cam3d_cad!(scene, eyeposition=eyeposition)
    #scene.camera.eyeposition = eyeposition
    #println(scene.camera.eyeposition)
    return scene
end
function move_eyeposition!(scene,increment) # [ dlat, dlongt ]
    adjust_eyeposition_sphereplot!(scene,sphereeyeposition[1]+increment[1],
    sphereeyeposition[2]+increment[2])
    return scene
end
function add_sphereplot_boundaries!(scene)
    plateIDfield = read_plateIDs(age)
    plateIDlist = find_plateID_list(plateIDfield)
    for plateID in world.plateIDlist
        maskfield = generate_mask_field(world.plateID,plateID)
        x2d,y2d = generate_outline_from_mask(maskfield)  # x y coords
        #x2d = x2d .+ 180.
        x3d,y3d,z3d = line_sphere2cart(x2d,y2d) # same as generate3doutline
        lines!(scene,x3d,y3d,z3d,color=:white)
    end
    return scene
end
function add_sphereplot_colors!(scene,field)
    colorfield = fill(0,nx,ny)
    for i in 1:nx
        for j in 1:ny
            colorfield[i,j] = field[i,ny-j+1]
        end
    end
    colorfield = convert(Array{Float64},colorfield)
    surface!(scene,spheregridx,spheregridy,spheregridz,
        color=colorfield,colormap=:lightrainbow )
    return scene
end
function add_sphereplot_streamlines!(scene,lat,longt)
    lats = range(lat-60,stop=lat+60,length=13)
    longts = range(longt-60,stop=longt+60,length=13)
    for lat in lats
        for longt in longts
            iworld = search_sorted_nearest(xcoords,longt)
            jworld = search_sorted_nearest(ycoords,lat)
            originx = [longt-1,longt,longt+1,longt]
            originy = [lat,lat+1,lat,lat-1]
            org3dx,org3dy,org3dz = line_sphere2cart(originx,originy)
            lines!(scene,org3dx,org3dy,org3dz,color=:white)
            plateID = world.plateID[iworld,jworld]
            x2d,y2d = generate_flowpath_line(lat, longt, plateID)
            x3d,y3d,z3d = line_sphere2cart(x2d,y2d) # same as generate3doutline
            lines!(scene,x3d,y3d,z3d,color=:white)
        end
    end
    return scene
end
function add_sphereplot_wig!(scene,field,maskfield,rotationmatrix)
    u = range(0,stop=2*π,length=nx)
    v = range(0,stop=π,length=ny)
    x = zeros(nx,ny); y = zeros(nx,ny); z = zeros(nx,ny)
    colorfield=fill(0.,nx,ny)
    for i in 1:nx
        for j in 1:ny
            if maskfield[i,ny-j+1] == 1
                elevation = 1.1
                colorfield[i,j] = field[i,ny-j+1]
            else
                elevation = 0.9
            end
            xplate = cos.(u[i]) * sin(v[j]) * elevation
            yplate = sin.(u[i]) * sin(v[j]) * elevation
            zplate = cos(v[j]) * elevation
            vpoint = [xplate,yplate,zplate]
            newvpoint = transpose(rotationmatrix) * vpoint
            # unrotating back to world grid
            x[i,j] = newvpoint[1]
            y[i,j] = newvpoint[2]
            z[i,j] = newvpoint[3]

        end
    end
    surface!(  x,y,z,color=colorfield,show_axis=false )
    return scene
end
function add_sphereplot_continent_outlines!(scene)
    continentmask = get_continent_mask()
    mycontour = Contour.contour(xcoords,ycoords,continentmask,0.5)
    lines = Contour.lines(mycontour)
    for line in lines
        xs,ys = Contour.coordinates(line)
        #xs = xs .+ 180.
        x3d,y3d,z3d = line_sphere2cart(xs,ys)
        #x3d = x3d .* 10.
        #x3d = [-1.,1.]; y3d = [-1.,1.]; z3d = [-1.,1.]
        lines!(scene,x3d,y3d,z3d,color=:black)
    end
    return scene
end
function sphere_boundary_plot(field,lat,longt)
    scene = setup_sphereplot(lat,longt)
    updatespherecolors!(scene,field)
    updatesphereboundaries!(scene)
    return scene
end
function sphere_stream_plot(field,lat,longt)
    scene = setupsphere_boundary_plot(lat,longt)
    add_sphereplot_streamlines!(scene,lat,longt)
    return scene
end

# netcdf file io
function nc_create_and_save(filename,varname,field)
    #bufferfield = deepcopy(field)
    nccreate(filename,varname,"Longitude",xcoords,"Latitude",ycoords)
    ncwrite(field,filename,varname)
    return
end
function save_world()
    runname = "SM2"
    timestamp = string(Int(ceil(world.age)))
    # in world grid, only the age field is non-trivial to calculate
    filename = "../outfiles/world/" * runname * ".world." * timestamp * ".bson"
    rm(filename, force=true)
    BSON.@save filename world
    return
end
function read_world(age)
    runname = "SM2"
    timestamp = string(Int(ceil(age)))
    # in world grid, only the age field is non-trivial to calculate
    filename = "../outfiles/world/" * runname * ".world." * timestamp * ".bson"
    BSON.@load filename world
    return
end
function save_plates()
    runname = "SM2"
    timestamp = string(Int(ceil(world.age)))
    # in world grid, only the age field is non-trivial to calculate
    filename = "../outfiles/plates/" * runname * ".plates." * timestamp * ".bson"
    rm(filename, force=true)
    BSON.@save filename plates
    return
end
function read_plates(age)
    runname = "SM2"
    timestamp = string(Int(ceil(age)))
    filename = "../outfiles/plates/" * runname * ".plates." * timestamp * ".bson"
    BSON.@load filename plates
    return
end

scotese_elevation_plot_ages = [535,530,525,520,515,510,505,500,495,490,485,480,
    475,470,465,460,455,450,445,440,435,430,425,420,415,410,405,400,395,390.5,
    385.2,380,375,370,365,355,350,345,340,335,330,325,320,315,310,305,300,295,
    290,285,280,275,270,265,260,255,245,240,235,230,225,220,215,210,205,200,
    195,190,185,180,175,170,165,160,155,150,145,140,135,130,125,120,115,110,
    100,95,90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,5,0]
function nearest_scotese_elevation()
    index = search_sorted_nearest(scotese_elevation_plot_ages,world.age)
    scotese_age = scotese_elevation_plot_ages[index]
    #println(scotese_age)
    stringage = string(Int(ceil(scotese_age)))
    #println(stringage)
    if Int(ceil(scotese_age)) != scotese_age
        stringage = string(scotese_age)
    end
    #println(stringage)
    filename = pwd() * "../scotese_elevations/scotese." * stringage * "Ma.nc"
    #println(filename)
    field = ncread(filename,"z")
    return field[1:360,1:180], scotese_age
end




if pwd() != "/Users/archer/Synched/papers/spongeball/codebase"
    cd("Synched/papers/spongeball/codebase")
end
rotations = read_rotation_file("1000_0_rotfile_Merdith_et_al.rot")
norotation = rotation2matrix(rotations[1])

orogenic_base_uplift = 20.
orogenic_events = Dict()
orogenic_events["Pan_African"] =
    create_orogenic_event("pan_african",650.,550.,orogenic_base_uplift)  # africa, south america
orogenic_events["Avalonian"] =
    create_orogenic_event("taconian",650.,500.,orogenic_base_uplift)
orogenic_events["Taconian"] =
    create_orogenic_event("taconian",490.,440.,3*orogenic_base_uplift)
orogenic_events["Calcedonian/Acadian"] =
    create_orogenic_event("calcedonian",460.,390.,orogenic_base_uplift)
orogenic_events["Hercynian/Alleghenian/Uralian"] =
    create_orogenic_event("hercynian",300.,250.,orogenic_base_uplift) # = variscan
# Amurian (Japan) should be covered by subduction_orogeny
orogenic_events["Indo Sinean"] =
    create_orogenic_event("indo_sinean",200.,180.,0.5*orogenic_base_uplift)
orogenic_events["Cimmerian"] =
    create_orogenic_event("cimmerian",180.,150.,orogenic_base_uplift)
orogenic_events["Mongol Okhotsk"] =
    create_orogenic_event("mongol_okhotsk",140.,130.,orogenic_base_uplift)
orogenic_events["Wrangellian"] =
    create_orogenic_event("wrangellian",110.,80.,orogenic_base_uplift)
orogenic_events["Verkhoyansk"] =
    create_orogenic_event("verkhoyansk",100.,70.,orogenic_base_uplift)
orogenic_events["Alpine"] =
    create_orogenic_event("alpine",50.,0.,orogenic_base_uplift)
orogenic_events["Himalayan"] =
    create_orogenic_event("himalayan",40.,0.,orogenic_base_uplift)

function create_everything( age )
    global world = create_world( age )
    global plates = create_empty_plates( )
    initialize_plates( )
    fill_world_from_plates( )
    return
end
#
create_everything( earliesttime )
print(" ")
println( "Initialization complete" )

println( Threads.nthreads() )

## big loop
#LUofAglobal = setup_2d_transport()

operation = "run"  # "run" or "movie"
image_number = 0
for age in ages[2:end]

    if operation == "run"
        step_everything()
        if true == true
            save_world()
            if floor(age/25) == age/25 && true == true
                save_plates()
            end
        end
    else
        read_world(age)
    end
    println(world.age, " ", get_geo_interval())

    if operation == "movie"  # make movie images
        if age < 535
            global image_number += 1
            scotese_elevation,scotese_age = nearest_scotese_elevation()
            scene = plot_two_fields(world.freeboard,scotese_elevation,-4000,4000)
            plot_add_plate_boundaries!(scene)
            plot_add_continent_outlines!(scene)
            plot_add_timestamp!(scene,world.age,-180,-105)
            scotesetimestamp = string(Int(floor(scotese_age))) * " Myr"
            text!(scene, scotesetimestamp, position = (-180,-305),textsize=15)
            #plot_add_timestamp!(scene,scotese_age,-180,-305)
            plot_add_orogenies!(scene)
            outfilename = pwd() * "/image_tmp/img." * lpad(image_number,3,"0") * ".png"
            #println(outfilename)
            Makie.save(outfilename,scene)
        end
        #ffmpeg -r 10 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p elevation.mp4
        #display(scene)
        #println(age," ",scotese_age, " ",get_geo_interval(age))
    end

    if true == false  # watch movie on screen
        scene = plotfield(world.freeboard)
        display(scene)
        println(highest_xy_field(world.freeboard))
    end

    if true == false # sediment conservation diagnostics
        weathering_input = volumefield_total(get_diag("crust_weathering_rate")) / time_step
        sediment_deposition =
            volumefield_total(get_diag("sediment_deposition_rate")) / time_step
        sediment_subduction =
            areafield_total(get_diag("ocean_subduct_sediment_plate_volume")) / time_step
        world_sediment_volume = global_world_sediment_inventory()
        pos_freeb = clamp_field(world.freeboard,0.)
            mean_exposed_elevation = field_mean(pos_freeb)
        ocean_seds = world.sediment_thickness .*
            generate_mask_field(world.crust_type,ocean_type)
        ocean_sediment_volume = volumefield_total(ocean_seds)
        #plate_sediment_volume = global_plate_sediment_inventory()
        @printf("%i %.2e %.2e %.2e %.2e %.2e %.1f\n", age, weathering_input,
            sediment_deposition,sediment_subduction,
            world_sediment_volume, ocean_sediment_volume,mean_exposed_elevation)
    end

end
