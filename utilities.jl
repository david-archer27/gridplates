
# Creators
function create_world( age )    # no plateIDs or numbers yet
    plateIDmap = fill(0,nx,ny)
    plateIDlist = fill(0,1)
    crust_type = fill(not_yet_defined,nx,ny)
    crust_age = fill(0.,nx,ny)
    sealevel = get_sealevel( age )
    atmCO2 = get_atmCO2( age )
    crust_thickness = fill(0.,nx,ny) # pertaining to geomorphology
    crust_density = fill(0.,nx,ny)
    geomorph_map = fill(0,nx,ny)
    tectonics_map = fill(0,nx,ny)
    sediment_thickness = fill(0.,nx,ny)
    sediment_layer_thickness = fill(0.,nx,ny,n_sediment_time_bins)
    sediment_surface_fractions = fill(0.,nx,ny,n_sediment_types)
    sediment_layer_fractions = fill(0.,nx,ny,n_sediment_types,
        n_sediment_time_bins)
    for i_sedtype in 1:n_sediment_types
        sediment_surface_fractions[:,:,i_sedtype] .= 
            initial_sediment_fractions[i_sedtype]
        sediment_layer_fractions[:,:,i_sedtype,1] .= 
            initial_sediment_fractions[i_sedtype]
    end
    elevation_offset = fill(0.,nx,ny)
    surface_elevation = fill(0.,nx,ny)
    freeboard = fill(0.,nx,ny)
    subducted_land_sediment_volumes = [0.,0.]
    subducted_ocean_sediment_volumes = [0.,0.]
    initial_ocean_sediment_inventories = fill(0.,n_sediment_types)
    initial_land_sediment_inventories = fill(0.,n_sediment_types)
    n_diags = length(world_diag_names); diags = fill(0.,nx,ny,n_diags)
    n_frac_diags = length(world_frac_diag_names)
    frac_diags = fill(0.,nx,ny,n_sediment_types,n_frac_diags)
    world = world_struct(age,sealevel,atmCO2,plateIDmap,plateIDlist,crust_type,crust_age,
        crust_thickness,crust_density,geomorph_map,tectonics_map,
        sediment_thickness,sediment_surface_fractions,
        sediment_layer_thickness,sediment_layer_fractions,
        elevation_offset,surface_elevation,freeboard,
        subducted_land_sediment_volumes,subducted_ocean_sediment_volumes,
        initial_ocean_sediment_inventories,initial_land_sediment_inventories,
        diags,frac_diags)
    world.plateID[:,:] .= read_plateIDs( age )
    world.plateIDlist = find_plateID_list( world.plateID )
    continentIDmap = read_continentIDs( age )
    for ix in 1:nx
        for iy in 1:ny
            if continentIDmap[ix,iy] > 0
                world.crust_type[ix,iy] = continent_crust
                world.crust_thickness[ix,iy] = continent_crust_h0
                world.crust_density[ix,iy] = rho_continent_crust
                world.geomorphology[ix,iy] = sedimented_land
                world.sediment_thickness[ix,iy] = initial_land_sediment_thickness
            else
                world.crust_type[ix,iy] = ocean_crust
                world.crust_thickness[ix,iy] = ocean_crust_h0
                world.crust_density[ix,iy] = rho_ocean_crust
                world.geomorphology[ix,iy] = pelagic_seafloor
            end
        end
    end
    return world
end
function create_blank_plate(plateID)
    crust_type = fill(not_yet_defined,nx,ny)
    crust_age = fill(0.,nx,ny)
    crust_thickness = fill(0.,nx,ny)
    crust_density = fill(0.,nx,ny)
    geomorphology = fill(0,nx,ny)
    tectonics_map = fill(0,nx,ny)
    #orogenic_uplift = fill(0,nx,ny)
    sediment_thickness = fill(0.,nx,ny)
    sediment_surface_fractions = fill(0.,nx,ny,n_sediment_types)
    sediment_layer_thickness = fill(0.,nx,ny,n_sediment_time_bins)
    sediment_layer_fractions = fill(0.,nx,ny,n_sediment_types,n_sediment_time_bins)
    for i_sedtype in 1:n_sediment_types
        sediment_surface_fractions[:,:,i_sedtype] .= 
            initial_sediment_fractions[i_sedtype]
        sediment_layer_fractions[:,:,i_sedtype,1] .= 
            initial_sediment_fractions[i_sedtype]
    end
    rotationmatrix = fill(0,3,3)
    resolvetime = -1
    parentstack = []
    lastparentstack = []
    firstappearance = -1
    lastappearance = -1
    plate = plate_struct(plateID,crust_type,crust_age,
        crust_thickness,crust_density,geomorphology,tectonics_map,
        #orogenic_uplift,
        sediment_thickness,sediment_surface_fractions,
        sediment_layer_thickness,sediment_layer_fractions,
        rotationmatrix,resolvetime,
        parentstack,lastparentstack,
        firstappearance,lastappearance)
    return plate
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
    #remask_plates()
    for (plateID,plate) in plates
        initial_mask_plate!(plate)
    end
    # initializes plate.crust_type = not_yet_defined points, 
    # now plate.crust_type = not_at_surface, ocean_crust, or continent_crust
    # nothing in plates.gridpoint_changes
    clear_world_process_arrays()
    # resets world.geomorphology, world.tectonics
end
function create_everything( age )
    global world = create_world( age )

    isostacy()
    new_sediment_thickness, new_sediment_surface_fractions = 
        setup_cap_carbonates()
    world.sediment_thickness[:,:] = new_sediment_thickness[:,:]
    world.sediment_surface_fractions[:,:,:] = new_sediment_surface_fractions[:,:,:]

    global plates = create_empty_plates( )
    initialize_plates( )
    #remask_plates()

    fill_world_from_plates( ) 

    #=isostacy()
    new_sediment_thickness, new_sediment_surface_fractions = 
        setup_cap_carbonates()
    world.sediment_thickness[:,:] = new_sediment_thickness
    world.sediment_surface_fractions[:,:,:] = new_sediment_surface_fractions

    apply_geomorphology_changes_to_plates()=#

    println("initial condition at ", age)
    # leaves world.geomorphology, world.geomorphology_changes alone
end
function test_regridding_conservation( d_age )
    global world = create_world( 0 )
    world.plateID .= 101
    world.plateIDlist = [ 101 ]
    global plates = create_empty_plates( )
    initialize_plates( )
    old_land_inv = world_land_sediment_inventories( )
    old_ocean_inv = world_ocean_sediment_inventories( )
    println("before rotation ", old_land_inv, old_ocean_inv )
    world.age = d_age
    fill_world_from_plates( ) 
    println("after rotation ",world_land_sediment_inventories( ), world_ocean_sediment_inventories( ))
    println("          diff ",world_land_sediment_inventories( ) - old_land_inv,
        world_ocean_sediment_inventories( ) - old_ocean_inv )
end

# Global Utilities
function apply_orogeny_fluxes_to_world( net_crust_change_rate )
    # applies continent_orogenic_uplift_rate and subduction_orogenic_uplift_rate
    # to world.crust_thickness, returns in isostatic eq
    uplift_rate = net_crust_change_rate * # meters / Myr
        main_time_step # sub_time_step  # meters / step
    for ix in 1:nx
        for iy in 1:ny
            if world.crust_type[ix,iy] == continent_crust && uplift_rate[ix,iy] > 0.
                world.crust_thickness[ix,iy] += uplift_rate[ix,iy]
            end
        end
    end
    return
end    
function apply_crust_weathering_fluxes_to_world()
    crust_erosion_rate = get_diag("crust_erosion_rate")
    for ix in 1:nx
        for iy in 1:ny
            world.crust_thickness[ix,iy] -= crust_erosion_rate[ix,iy] * main_time_step
        end
    end
end
function apply_crust_erosion_rate()
    erosion_rate = get_diag("crust_erosion_rate")
    world.crust_thickness .-= erosion_rate .* main_time_step
end
function apply_land_sediment_fluxes( land_sediment_fraction_deposition_rate_fields )
    # the sediments on land are unresolved in time, stored in 2d arrays
    # called provisionally by set_land_runoff_fluxes, then again at the
    # end of the step.  hence returning fields to be applied, rather than
    # altering world.sediment_thickness and world.sediment_surface_fractions
    new_sediment_thickness = deepcopy(world.sediment_thickness)
    for ix in 1:nx
        for iy in 1:ny
            for i_sedtype in 1:n_sediment_types
                new_sediment_thickness[ix,iy] += 
                    land_sediment_fraction_deposition_rate_fields[ix,iy,i_sedtype] * 
                    main_time_step
                if new_sediment_thickness[ix,iy] < 0.
                    #println("negative thickness ",[ix,iy,
                    #    world.sediment_thickness[ix,iy], 
                    #    new_sediment_thickness[ix,iy]])
                    new_sediment_thickness[ix,iy] = 0.
                end
            end
        end
    end
    new_sediment_fractions = deepcopy(world.sediment_surface_fractions)
    for ix in 1:nx
        for iy in 1:ny
            if new_sediment_thickness[ix,iy] > 0.
                for i_sedtype in 1:n_sediment_types
                    old_inventory = world.sediment_thickness[ix,iy] * 
                        world.sediment_surface_fractions[ix,iy,i_sedtype]
                    new_inventory = old_inventory +
                        land_sediment_fraction_deposition_rate_fields[ix,iy,i_sedtype] * main_time_step
                    new_concentration = new_inventory / 
                        new_sediment_thickness[ix,iy]
                    new_sediment_fractions[ix,iy,i_sedtype] = 
                        new_concentration
                end
            end
        end
    end
    return new_sediment_thickness, new_sediment_fractions
end
function apply_ocean_sediment_fluxes( )
    # the sediment record is stored time-resolved on ocean points in 
    # sediment_layer_thickness and sediment_layer_fractions, as well 
    # as the total thickness in world.sediment_thickness and 
    # the surface values into world.sediment_surface_fractions
    # this is only called once, at the end of the step, so there is
    # no need to externalize dealing with world.* variables
    verbose = true
    #depositing_seafloor_mask = eq_mask(world.geomorphology,pelagic_seafloor) .+
    #    eq_mask(world.geomorphology,coastal_depocenter)
    seafloor_sediment_deposition_rate = get_diag("seafloor_sediment_deposition_rate")
    sediment_fraction_deposition_rates = 
        get_frac_diag("seafloor_sediment_fraction_deposition_rate")
    new_total_sediment_thickness = deepcopy(world.sediment_thickness)
    new_sediment_layer_thickness = deepcopy(world.sediment_layer_thickness)
    new_sediment_layer_fractions = deepcopy(world.sediment_layer_fractions)
    new_total_sediment_thickness[:,:] .+=
        seafloor_sediment_deposition_rate .* main_time_step
    new_sediment_layer_thickness[:,:,current_time_bin()] .+=
        seafloor_sediment_deposition_rate .* main_time_step
    
    #if verbose == true
    #    println("apply total dep ", volume_field(seafloor_sediment_deposition_rate), 
    #        " fracs ", volume_field(sediment_fraction_deposition_rates[:,:,1]),
    #        " ", volume_field(sediment_fraction_deposition_rates[:,:,2]))
    #end

    for ix in 1:nx
        for iy in 1:ny
            if seafloor_sediment_deposition_rate[ix,iy] > 0.

                for i_sedtype in 1:n_sediment_types
                    old_inventory = world.sediment_layer_thickness[ix,iy,current_time_bin()] * 
                        world.sediment_layer_fractions[ix,iy,i_sedtype,current_time_bin()]
                    new_inventory = old_inventory +
                        sediment_fraction_deposition_rates[ix,iy,i_sedtype] * main_time_step
                    new_concentration = new_inventory / 
                        new_sediment_layer_thickness[ix,iy,current_time_bin()]
                    world.sediment_layer_fractions[ix,iy,i_sedtype,current_time_bin()] = 
                        new_concentration
                    world.sediment_surface_fractions[ix,iy,i_sedtype] = 
                        new_concentration
                end
                world.sediment_thickness[ix,iy] = new_total_sediment_thickness[ix,iy]
                world.sediment_layer_thickness[ix,iy,:] = new_sediment_layer_thickness[ix,iy,:]

            end
        end
    end
end
function update_world_continents_from_file() 
    # sets new continents
    local new_crust_type
    world.plateID = read_plateIDs() # wants new age in world.age
    continentID = read_continentIDs()
    for ix in 1:nx
        for iy in 1:ny
            if continentID[ix,iy] > 0
                new_crust_type = continent_crust
            else
                new_crust_type = ocean_crust
            end
            if world.crust_type[ix,iy] != new_crust_type
                if new_crust_type == continent_crust
                    if enable_crust_ocn2cont
                        world.geomorphology[ix,iy] = sedimented_land
                        world.tectonics[ix,iy] = ocean_to_continent
                        accum_diag("ocean_2_continent_world_area",
                            ix,iy,1.)
                        world.crust_thickness[ix,iy] = continent_crust_h0
                        world.crust_density[ix,iy] = rho_continent_crust
                        world.crust_type[ix,iy] = continent_crust 
                        sed_overflow_fractions = fill(0.,n_sediment_types)
                        for i_sedtype in 1:n_sediment_types
                            for i_bin in 1:n_sediment_time_bins
                                sed_overflow_fractions[i_sedtype] += 
                                    world.sediment_layer_thickness[ix,iy,i_bin] *
                                    world.sediment_layer_fractions[ix,iy,i_sedtype,i_bin]
                            end
                            accum_frac_diag("ocean_2_continent_sediment_fraction_displaced",
                                ix,iy,i_sedtype,sed_overflow_fractions[i_sedtype] / main_time_step )
                                # pick this up and add to ocean sedimentation fluxes
                                # so squishing ocean sediment sideways when uplift happens
                        end  
                        world.sediment_thickness[ix,iy] = initial_land_sediment_thickness
                        world.sediment_layer_thickness[ix,iy,:] .= 0.
                        world.sediment_surface_fractions[ix,iy,:] .= 0.
                        world.sediment_surface_fractions[ix,iy,1] = 1.  
                        world.sediment_layer_thickness[ix,iy,:] .= 0.
                        world.sediment_layer_thickness[ix,iy,1] = initial_land_sediment_thickness
                        world.sediment_layer_fractions[ix,iy,:,:] .= 0.
                        world.sediment_layer_fractions[ix,iy,1,:] .= 1.
                        world.elevation_offset[ix,iy] = 0.
                    end
                else # new_crust_type == ocean_crust
                    if enable_crust_cont2ocn
                        #error([ix,iy])
                        world.geomorphology[ix,iy] = pelagic_seafloor
                        world.tectonics[ix,iy] = continent_to_ocean
                        accum_diag("continent_2_ocean_world_area",
                            ix,iy,1.)
                        world.crust_thickness[ix,iy] = ocean_crust_h0
                        world.crust_density[ix,iy] = rho_ocean_crust
                        for i_sedtype in 1:n_sediment_types
                            sed_overflow_fraction = world.sediment_thickness[ix,iy] *
                                world.sediment_surface_fractions[ix,iy,i_sedtype]
                            accum_frac_diag("continent_2_ocean_sediment_fraction_displaced",
                                ix,iy,i_sedtype,sed_overflow_fraction / main_time_step )
                        end
                        world.sediment_layer_thickness[ix,iy,:] .= 0.
                        world.sediment_thickness[ix,iy] = 0.
                        world.sediment_surface_fractions[ix,iy,:] .= 0.
                        world.sediment_surface_fractions[ix,iy,1] = 1.
                        world.crust_type[ix,iy] = new_crust_type 
                    end
                end
            end
        end
    end
    return
end
function scale_global_sediment_components( scale_factors )
    equivalent_thickness = fill(0.,n_sediment_types)
    for ix in 1:nx
        for iy in 1:ny
            #if world.crust_type[ix,iy] == continent_crust
            for i_sedtype in 1:n_sediment_types
                equivalent_thickness[i_sedtype] = 
                    world.sediment_thickness[ix,iy] *
                    world.sediment_surface_fractions[ix,iy,i_sedtype] *
                    scale_factors[i_sedtype]
            end
            world.sediment_thickness[ix,iy] = 0.
            for i_sedtype in 1:n_sediment_types
                world.sediment_thickness[ix,iy] += 
                    equivalent_thickness[i_sedtype]
            end
            if world.sediment_thickness[ix,iy] > 0.
                for i_sedtype in 1:n_sediment_types
                    world.sediment_surface_fractions[ix,iy,i_sedtype] = 
                        equivalent_thickness[i_sedtype] / 
                        world.sediment_thickness[ix,iy]
                end
            end
            if world.crust_type[ix,iy] == ocean_crust
                for ibin in 1:n_sediment_time_bins
                    for i_sedtype in 1:n_sediment_types
                        equivalent_thickness[i_sedtype] = 
                            world.sediment_layer_thickness[ix,iy,ibin] *
                            world.sediment_layer_fractions[ix,iy,i_sedtype,ibin] *
                            scale_factors[i_sedtype]
                    end
                    world.sediment_layer_thickness[ix,iy,ibin] = 0.
                    for i_sedtype in 1:n_sediment_types
                        world.sediment_layer_thickness[ix,iy,ibin] += 
                            equivalent_thickness[i_sedtype]
                    end
                    if world.sediment_layer_thickness[ix,iy,ibin] > 0.
                        for i_sedtype in 1:n_sediment_types
                            world.sediment_layer_fractions[ix,iy,i_sedtype,ibin] = 
                                equivalent_thickness[i_sedtype] / 
                                world.sediment_layer_thickness[ix,iy,ibin]
                        end
                    end
                end
            end
        end
    end
end
function increment_world_age()
    world.crust_age .+= main_time_step
    return
end
function top_filled_sediment_box(plate,iplate,jplate)
    ibinmax = search_sorted_below(sediment_time_bins,world.age)
    ibinmax = max(ibinmax,1)
    ibintop = 1
    for ibin in 2:ibinmax
        if plate.sediment_thickness[iplate,jplate,ibin] > 0.
            ibintop = ibin
        end
    end
    return ibintop
end
function current_time_bin()
    ibin = 0
    for itest in 1:n_sediment_time_bins
        if sediment_time_bins[itest] >= world.age
            ibin = itest
        end
    end
    return ibin
end


function ix_left(ix)
    if ix == 1
        ix_l = nx
    else
        ix_l = ix - 1
    end
    return ix_l
end
function ix_right(ix)
    if ix == nx
        ix_r = 1
    else
        ix_r = ix + 1
    end
    return ix_r
end
function ix_eq_size_left(ix,iy)
    cell_list = []
    number_of_cells_hor = Int( floor( delta_y / delta_x[iy] ) )
    n_cells_found = 0
    ix_tmp = ix
    while n_cells_found < number_of_cells_hor
        ix_tmp -= 1
        if ix_tmp == 0
            ix_tmp += 360
        end 
        push!(cell_list,[ix_tmp,iy])
    end
    return cell_list
end
function ix_eq_size_right(ix,iy)
    cell_list = []
    number_of_cells_hor = Int( floor( delta_y / delta_x[iy] ) )
    n_cells_found = 0
    ix_tmp = ix
    while n_cells_found < number_of_cells_hor
        ix_tmp += 1
        if ix_tmp > 360
            ix_tmp -= 360
        end 
        push!(cell_list,[ix_tmp,iy])
    end
    return cell_list
end
function get_border_neighbor_eq_size_coords(ix,iy)
    #number_of_cells_hor = Int( floor( delta_y / delta_x[iy] ) )
    neighbor_coords = ix_eq_size_left(ix,iy)
    append!( neighbor_coords, ix_eq_size_right(ix,iy))
    if iy > 1 
        append!( neighbor_coords, [ix,iy-1] )
    end
    if iy < ny
        append!( neighbor_coords, [ix,iy+1] )
    end
    return neighbor_coords
end
function get_border_neighbor_coords(ix,iy)
    local neighbor_coords
    if iy == 1 
        neighbor_coords = [ [ix_left(ix), iy], [ix_right(ix),iy],
            [ix,iy+1] ]
    elseif iy == ny
        neighbor_coords = [ [ix_left(ix), iy], [ix_right(ix),iy],
            [ix,iy-1] ]
    else
        neighbor_coords = [ [ix_left(ix), iy], [ix_right(ix),iy],
            [ix,iy-1],[ix,iy+1] ]
    end
    return neighbor_coords
end
#function get_layered_neighbor_coords(ix,iy,xbuffer,ybuffer)

function get_neighbor_coords(ixcenter,iycenter,xbuffer,ybuffer) # full block includes corner points
    ixleft = ixcenter - xbuffer; ixright = ixcenter + xbuffer
    iylower = max(iycenter - ybuffer,1); iyupper = min(iycenter + ybuffer,ny)
    neighbor_coords = []
    for ix in ixleft:ixright
        ixshift = ix
        if ix < 1 
            ixshift += nx
        end
        if ix > nx
            ixshift -= nx
        end
        for iy in iylower:iyupper
            if [ixshift,iy] != [ixcenter,iycenter]
                push!(neighbor_coords,[ixshift,iy])
            end
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
function is_land( )
    is_land_field = eq_mask(world.crust_type,continent_crust)
    return is_land_field
end
function is_ocean( )
    is_ocean_field = eq_mask(world.crust_type,ocean_crust)
    return is_ocean_field
end
function eq_mask(field,value)
    newfield = fill(0,nx,ny)
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
function gt_mask(field,value)
    newfield = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if field[ix,iy] > value
                newfield[ix,iy] = 1
            end
        end
    end
    return newfield
end
function lt_mask(field,value)
    newfield = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if field[ix,iy] < value
                newfield[ix,iy] = 1
            end
        end
    end
    return newfield
end
function ge_mask(field,value)
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
function le_mask(field,value)
    newfield = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if field[ix,iy] <= value
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
function get_continent_mask()
    continentmask = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if world.crust_type[ix,iy] >= continent_crust
                continentmask[ix,iy] = 1.
            end
        end
    end
    return continentmask
end

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

# storage of diagnostics variables
function diag_index(varname)
    index = findfirst(isequal(varname),world_diag_names)
    return index
end
function frac_diag_index(varname)
    index = findfirst(isequal(varname),world_frac_diag_names)
    return index
end
function get_diag(varname)
    field = world.diags[:,:,findfirst(isequal(varname),world_diag_names)]
    return field
end
function get_diag(varname,ix,iy)
    value = world.diags[ix,iy,findfirst(isequal(varname),world_diag_names)]
    return value
end
function get_frac_diag(varname,sediment_type)
    field = world.frac_diags[:,:,sediment_type,
        findfirst(isequal(varname),world_frac_diag_names)]
    return field
end
function get_frac_diag(varname)
    fields = world.frac_diags[:,:,:,
        findfirst(isequal(varname),world_frac_diag_names)]
    return fields
end
function set_diag(varname,ix,iy,value)
    ibin = findfirst(isequal(varname),world_diag_names)
    world.diags[ix,iy,ibin] = value
end
function set_diag(varname,value)
    ibin = findfirst(isequal(varname),world_diag_names)
    world.diags[:,:,ibin] .= value
end
function set_frac_diag(varname,ix,iy,sediment_type,value)
    ibin = findfirst(isequal(varname),world_frac_diag_names)
    world.frac_diags[ix,iy,sediment_type,ibin] = value
end
function set_frac_diag(varname,sediment_type,field)
    ibin = findfirst(isequal(varname),world_frac_diag_names)
    world.frac_diags[:,:,sediment_type,ibin] = field
end
function set_frac_diag(varname,field)
    ibin = findfirst(isequal(varname),world_frac_diag_names)
    world.frac_diags[:,:,1:n_sediment_types,ibin] = 
        field[:,:,1:n_sediment_types]
end
function reset_diag(varname)
    ibin = findfirst(isequal(varname),world_diag_names)
    world.diags[:,:,ibin] = fill(0.,nx,ny)
end
function reset_frac_diag(varname)
    ibin = findfirst(isequal(varname),world_frac_diag_names)
    world.frac_diags[:,:,:,ibin] .= 0.
end
function accum_diag(varname,field)
    ibin = findfirst(isequal(varname),world_diag_names)
    world.diags[ix,iy,ibin] .+= field
end
function accum_diag(varname,ix,iy,value)
    ibin = findfirst(isequal(varname),world_diag_names)
    world.diags[ix,iy,ibin] += value
end
function accum_frac_diag(varname,ix,iy,sediment_type,value)
    ibin = findfirst(isequal(varname),world_frac_diag_names)
    world.frac_diags[ix,iy,sediment_type,ibin] += value
end
function accum_frac_diag(varname,sediment_type,field)
    ibin = findfirst(isequal(varname),world_frac_diag_names)
    world.frac_diags[:,:,sediment_type,ibin] += field
end
function accum_frac_diag(varname,field)
    ibin = findfirst(isequal(varname),world_frac_diag_names)
    world.frac_diags[:,:,1:end,ibin] += field
end
function accum_diag(varname,value)
    ibin = findfirst(isequal(varname),world_diag_names)
    world.diags[:,:,ibin] .+= value[:,:]
    return world.diags[:,:,ibin]
end
function world_diag_pos(diag_name)
    diag_pos = findfirst(isequal(diag_name),world_diag_names)
    return diag_pos
end
#=function record_flux_into_world_diags(fluxname,jsource,iworld,jworld,mult::Float64=1.)
    area_ratio = areabox[jsource] / areabox[jworld]
    accum_diag(fluxname,iworld,jworld,mult)
    return
end
function record_sediment_fraction_flux_into_world_diags(fluxname,sediment_type,
    jsource,iworld,jworld,mult::Float64=1.)
    #area = areabox[jsource]
    accum_frac_diag(fluxname,iworld,jworld,sediment_type,mult)
    return
end=#
function clear_world_process_arrays()
    world.geomorphology .= 0
    world.tectonics .= 0
    world.diags .= 0.
    world.frac_diags .= 0.
    return
end
function clear_geomorph_process_arrays()
    for idiag in [    
        "crust_erosion_rate",               # units m/Myr, calc on substep
        #"crust_clay_source_rate",
        #"aolean_clay_erosion_rate", 
        #"aolean_clay_deposition_rate", 
        #"land_orogenic_clay_flux", # in orogeny-neighboring land grid cells
        #"land_CaCO3_dissolution_rate", # subaereal erosion
        "continental_CaCO3_deposition_rate", # when flooded
        "land_sediment_deposition_rate",
        "seafloor_sediment_deposition_rate",
        "global_sediment_deposition_rate",   
        #"coastal_orogenic_clay_flux", # in coastal ocean points, boundary fluxes
        "coastal_CaCO3_flux",
        "pelagic_CaCO3_deposition_rate",
        "seafloor_delta_CO3"]
        #println(idiag)
        reset_diag(idiag)
    end
    for idiag in [    
        "denuded_land_boundary_fraction_flux",
        "denuded_coastal_boundary_fraction_flux",
        "crust_orogenic_fraction_flux",
        "land_orogenic_fraction_flux",
        "coastal_orogenic_fraction_flux",
        "land_sediment_fraction_dissolution_rate",
        "land_sediment_fraction_deposition_rate",
        "land_trapped_sediment_rate",
        "coastal_sediment_fraction_runoff_flux",
        "ocean_sediment_fraction_influx",
        "seafloor_sediment_fraction_deposition_rate", 
        "seafloor_sediment_fraction_overflow",
        "global_sediment_fraction_deposition_rate",
        "global_sediment_fraction_deposition_ratio"]

        reset_frac_diag(idiag)
    end
end
function world_land_sediment_inventories(  )
    inventories = fill(0., 0: n_sediment_types) # total, fraction1, fraction2 
    for ix in 1:nx
        for iy in 1:ny
            if world.crust_type[ix,iy] == continent_crust
                inventories[0] += world.sediment_thickness[ix,iy] *
                    areabox[iy]
                for i_sedtype in 1:n_sediment_types
                    inventories[ i_sedtype ] += 
                        world.sediment_thickness[ix,iy] *
                        world.sediment_surface_fractions[ix,iy,i_sedtype] *
                        areabox[iy]
                end 
            end
        end
    end
    return inventories
end
function world_ocean_sediment_inventories(  )
    inventories = fill(0., 0: n_sediment_types)
    sum_inventories = fill(0., 0: n_sediment_types)
    for ix in 1:nx
        for iy in 1:ny
            if world.crust_type[ix,iy] == ocean_crust
                inventories .= 0.
                for ibin in 1:n_sediment_time_bins
                    inventories[0] += world.sediment_layer_thickness[ix,iy,ibin] *
                        areabox[iy]
                    for i_sedtype in 1:n_sediment_types
                        inventories[ i_sedtype ] += 
                            world.sediment_layer_thickness[ix,iy,ibin] *
                            world.sediment_layer_fractions[ix,iy,i_sedtype,ibin] *
                            areabox[iy]
                    end
                end
                if enable_check_sediment_layer_inventory
                    unbinned_inventory = world.sediment_thickness[ix,iy] *
                        areabox[iy]
                    if abs( unbinned_inventory - inventories[0] ) > 1.e8
                        error("world ocean sed inventory out of whack ",[ix,iy],
                            unbinned_inventory," ",inventories[0])
                    end
                end
                sum_inventories .+= inventories
            end
        end
    end
    return sum_inventories
end
function plate_land_sediment_inventories( plate )
    inventories = fill(0., 0: n_sediment_types) # total, fraction1, fraction2 
    for ix in 1:nx
        for iy in 1:ny
            if plate.crust_type[ix,iy] == continent_crust
                inventories[0] += plate.sediment_thickness[ix,iy] *
                    areabox[iy]
                for i_sedtype in 1:n_sediment_types
                    inventories[ i_sedtype ] += 
                        plate.sediment_thickness[ix,iy] *
                        plate.sediment_surface_fractions[ix,iy,i_sedtype] *
                        areabox[iy]
                end 
            end
        end
    end
    return inventories
end
function plate_ocean_sediment_inventories( plate )
    inventories = fill(0., 0: n_sediment_types) # total, fraction1, fraction2 
    sum_inventories = fill(0., 0: n_sediment_types)
    for ix in 1:nx
        for iy in 1:ny
            if plate.crust_type[ix,iy] == ocean_crust
                inventories .= 0.
                for ibin in 1:n_sediment_time_bins
                    inventories[0] += plate.sediment_layer_thickness[ix,iy,ibin] *
                        areabox[iy]
                    for i_sedtype in 1:n_sediment_types
                        inventories[ i_sedtype ] += 
                            plate.sediment_layer_thickness[ix,iy,ibin] *
                            plate.sediment_layer_fractions[ix,iy,i_sedtype,ibin] *
                            areabox[iy]
                    end
                end
                if enable_check_sediment_layer_inventory
                    unbinned_inventory = plate.sediment_thickness[ix,iy] *
                        areabox[iy]
                    if abs( unbinned_inventory - inventories[0] ) > 1.e6
                        error("inventory out of whack ",plate.plateID," ",[ix,iy],unbinned_inventory," ",inventories[0])
                    end
                end
                sum_inventories .+= inventories
            end
        end
    end
    return sum_inventories
end
function world_sediment_inventories(  )
    inventories = world_land_sediment_inventories() .+ 
        world_ocean_sediment_inventories(  )
    return inventories
end
function plate_sediment_inventories( plate )
    inventories = plate_land_sediment_inventories( plate ) .+ 
        plate_ocean_sediment_inventories( plate )
    #println("plate ", plate.plateID," ", inventories)
    return inventories
end
function global_plate_sediment_inventories()
    global_inventories = fill(0., 0:n_sediment_types)
    for plateID in world.plateIDlist
        if haskey(plates,plateID) == false
            plate = create_blank_plate(plateID)
            plates[plateID] = plate
            println("creating plate in global_plate_sediment_inventories ", plateID)
        end
        plate = plates[ plateID ]
        plate_inv = plate_sediment_inventories( plate )
        global_inventories .+= plate_inv
        #println("plate ", plateID," ", plate_inv, global_inventories)
    end
    return global_inventories
end
function global_plate_land_sediment_inventories()
    global_inventories = fill(0., 0:n_sediment_types)
    for plateID in world.plateIDlist
        if haskey(plates,plateID) == false
            plate = create_blank_plate(plateID)
            plates[plateID] = plate
            println("creating plate in global_plate_sediment_inventories ", plateID)
        end
        plate = plates[ plateID ]
        plate_inv = plate_land_sediment_inventories( plate )
        global_inventories .+= plate_inv
        #println("plate ", plateID," ", plate_inv, global_inventories)
    end
    return global_inventories
end
function global_plate_ocean_sediment_inventories()
    global_inventories = fill(0., 0:n_sediment_types)
    for plateID in world.plateIDlist
        if haskey(plates,plateID) == false
            plate = create_blank_plate(plateID)
            plates[plateID] = plate
            println("creating plate in global_plate_sediment_inventories ", plateID)
        end
        plate = plates[ plateID ]
        plate_inv = plate_ocean_sediment_inventories( plate )
        global_inventories .+= plate_inv
        #println("plate ", plateID," ", plate_inv, global_inventories)
    end
    return global_inventories
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
function land_elevation()
    is_land = gt_mask(world.freeboard, 0.)
    land_elevation = world.freeboard .* is_land
    return land_elevation
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
function field_mean(field, applies_map = fill(1,nx,ny))
    areatot = 0.
    valuetot = 0.
    for ix in 1:nx
        for iy in 1:ny
            if applies_map[ix,iy] == 1
                areatot += areabox[iy]
                valuetot += areabox[iy] * field[ix,iy]
            end
        end
    end
    valueavg = valuetot / areatot
    return valueavg
end
function field_max(field)
    maxval = -1.e10
    for ix in 1:nx
        for iy in 1:ny
            if field[ix,iy] > maxval
                maxval = field[ix,iy]
            end
        end
    end
    return maxval
end
function volume_fields( meterfields ) # [ix,iy,i_sedtype] array
    results = []
    for index in eachindex(meterfields[1,1,:])
        push!(results, volume_field(meterfields[:,:,index]))
    end
    return results
end
function volume_field(meterfield) # used for meters of erosion to m3 globally
    total = 0.
    for ix in 1:nx
        for iy in 1:ny
            if meterfield[ix,iy] != meterfield[ix,iy]
                error("Nan at ",[ix,iy])
            end
            total += meterfield[ix,iy] * areabox[iy] # m3
        end
    end
    if total != total
        error("NaN detected")
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
                if plate.crust_type[ix,iy] != not_at_surface
                    totalarea += areabox[iy]
                end
            end
        end
    end
    return totalarea
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
    #println("crow ",[a,1-a])
    c = 2. * atan(sqrt(a),sqrt(1-a))
    distance_meters = r * c
    return distance_meters
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

scotese_elevation_plot_ages = [535,530,525,520,515,510,505,500,495,490,485,480,
    475,470,465,460,455,450,445,440,435,430,425,420,415,410,405,400,395, # 390.5,
    #385.2,
    380,375,370,365,355,350,345,340,335,330,325,320,315,310,305,300,295,
    290,285,280,275,270,265,260,255,245,240,235,230,225,220,215,210,205,200,
    195,190,185,180,175,170,165,160,155,150,145,140,135,130,125,120,115,110,
    100,95,90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,5,0]

function scotese_mean_elevation_timeseries()
    mean_areas = Float32[]
    cd( scotese_data_directory )
    for scotese_age in scotese_elevation_plot_ages
        stringage = string(Int(ceil(scotese_age)))
        filename = "scotese." * stringage * "Ma.bson"
        println("reading ", filename)
        BSON.@load filename field
        sum_elevation = 0.
        sum_area = 0.
        for ix in 1:nx
            for iy in 1:ny
                if field[ix,iy] > 0.
                    sum_elevation += field[ix,iy] * areabox[iy]
                    sum_area += areabox[iy]
                end
            end
        end
        mean_elevation = sum_elevation / sum_area
        println("reading ", filename, " ", mean_elevation)
        push!(mean_areas,mean_elevation)
    end
    return mean_areas
end

    
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
    filename = scotese_data_directory * "/scotese." * stringage * "Ma.bson"
    println("reading ", filename)
    BSON.@load filename field
    return field, scotese_age
end
function check_blob_completeness( blobfield, unaccountedforfield )
    neighbors = get_blob_neighbors( blobfield )
    #ixr = 0; iyr = 0
    for ix in 1:nx
        for iy in 1:ny
            if neighbors[ix,iy] == 1 && unaccountedforfield[ix,iy] == 1
                #println("found an edge case ",ix," ",iy)
                #ixr = ix; iyr = iy
                return ix,iy
            end
        end
    end
    return 0,0
end

function get_blobs(maskfield)
    unaccountedforfield = deepcopy(maskfield) # ones in the mask
    blobs = []
    n_recursions = 0; n_restarts = 0
    for ix in 1:nx
        for iy in 1:ny
            if unaccountedforfield[ix,iy] == 1 # found a blob
                #println("launching blobfinder ",ix," ",iy)
                blobfield = fill(0,nx,ny)
                #println(length(blobfield))
                recursive_check_point_for_blob(blobfield,unaccountedforfield,
                    n_recursions,ix,iy)
                ixc,iyc = check_blob_completeness( blobfield, unaccountedforfield )
                #println("launching ", ixc," ",iyc)
                while ixc > 0
                    n_restarts += 1
                    #println("blobfinder relaunching ", ixc," ",iyc)
                    recursive_check_point_for_blob(blobfield,unaccountedforfield,
                        n_recursions,ixc,iyc)
                    ixc,iyc = check_blob_completeness( blobfield, unaccountedforfield )
                end
                push!(blobs,blobfield)
            end
        end
    end
    if n_restarts > 100
        println("blobfinder restarts ", n_restarts)
    end
    return blobs # array of arrays, one per continuous area
end
function recursive_check_point_for_blob(blobfield,unaccountedforfield,n_recursions,ix,iy)
    n_recursions_max = 8000
    if unaccountedforfield[ix,iy] == 1
        blobfield[ix,iy] = 1
        unaccountedforfield[ix,iy] = 0
        n_recursions += 1
        if n_recursions <= n_recursions_max 
            coordpairs = get_border_neighbor_coords(ix,iy)
            for ic in 1:length(coordpairs)
                coordpair = coordpairs[ic]
                ixt = coordpair[1]; iyt = coordpair[2]
                if unaccountedforfield[ixt,iyt] == 1 
                    recursive_check_point_for_blob(blobfield,unaccountedforfield,
                        n_recursions,ixt,iyt)
                end
            end
        end
    end
    return
end
function get_blob_neighbors(maskfield)
    neighbor_field = fill(0,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if maskfield[ix,iy] == 1
                coordpairs = get_border_neighbor_coords(ix,iy)
                for coordpair in coordpairs
                    ixt = coordpair[1]; iyt = coordpair[2]
                    if maskfield[ixt,iyt] == 0
                        neighbor_field[ixt,iyt] = 1
                    end
                end
            end
        end
    end
    return neighbor_field
end
function get_blob_neighbor_fields( maskfields )
    blob_neighbors = []
    for i_area in 1:length(maskfields)
        neighbor_field = get_blob_neighbors( maskfields[i_area] )
        push!( blob_neighbors, neighbor_field )
    end
    return blob_neighbors
end
function get_blob_fringe_fields( maskfield )
    blob_neighbors = []
    for i_area in 1:length(maskfield)
        neighbor_field = get_blob_fringe( maskfield[i_area] )
        push!( blob_neighbors, neighbor_field )
    end
    return blob_neighbors
end
function get_maskfield_census( maskfield )
    n_points = 0
    for ix in 1:nx
        for iy in 1:ny
            if maskfield[ix,iy] == 1
                n_points += 1
            end
        end
    end
    return n_points
end
function get_blob_fringe_exterior_list( maskfield )
    fringe_list = []
    for ix in 1:nx
        for iy in 1:ny
            if maskfield[ix,iy] == 1
                coordpairs = get_border_neighbor_coords(ix,iy) # this includes corner points
                for coordpair in coordpairs
                    ixt = coordpair[1]; iyt = coordpair[2]
                    if maskfield[ixt,iyt] == 0
                        push!(fringe_list,[ixt,iyt])
                    end
                end
            end
        end
    end
    return fringe_list
end
function get_blob_fringe_interior_list( maskfield )
    fringe_list = []
    for ix in 1:nx
        for iy in 1:ny
            if maskfield[ix,iy] == 1
                coordpairs = get_border_neighbor_coords(ix,iy) # this includes corner points
                for coordpair in coordpairs
                    ixt = coordpair[1]; iyt = coordpair[2]
                    if maskfield[ixt,iyt] == 0
                        push!(fringe_list,[ix,iy])
                    end
                end
            end
        end
    end
    return fringe_list
end
function get_blob_fringe( maskfield )
    fringe_field = fill(0,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if maskfield[ix,iy] == 1
                coordpairs = get_border_neighbor_coords(ix,iy) # this includes corner points
                for coordpair in coordpairs
                    ixt = coordpair[1]; iyt = coordpair[2]
                    if maskfield[ixt,iyt] == 0
                        fringe_field[ix,iy] = 1
                    end
                end
            end
        end
    end
    return fringe_field
end
function get_blob_onion_layers( maskfield ) 
    fringe_fields = []
    choppable_maskfield = deepcopy(maskfield)
    n_layers = 0
    voltot = volume_field(choppable_maskfield)
    while voltot > 0. # && n_layers < 10
        fringe_field = get_blob_fringe( choppable_maskfield )
        push!(fringe_fields,fringe_field)
        choppable_maskfield .-= fringe_field
        voltot = volume_field(choppable_maskfield)
        n_layers += 1
    end
    return fringe_fields
end
function show_sediment_budget()
    crustprod = volume_field(get_diag("crust_clay_source_rate"))
    println("total production  ",crustprod)
    totsed = volume_field(get_diag("sediment_deposition_rate"))
    println(" total deposit    ",totsed)
    oro2land = volume_field(get_diag("land_orogenic_clay_flux"))
    println(" oro to land      ",oro2land)
    landsed = volume_field(get_diag("land_clay_deposition_rate"))
    println(" land deposition  ",landsed)
    runoff = volume_field(get_diag("seafloor_runoff_clay_flux"))
    println(" runoff           ",runoff)
    aolerod = volume_field(get_diag("aolean_erosion_rate"))
    println(" aolean erosion   ",aolerod)
    println(" land balance     ",oro2land - landsed - runoff - aolerod)
    oro2ocn = volume_field(get_diag("coastal_orogenic_clay_flux"))
    println(" to ocean         ",oro2ocn)
    ocnsed = volume_field(get_diag("seafloor_clay_deposition_rate"))
    println(" ocean deposition ",ocnsed)
    aoldep = volume_field(get_diag("aolean_deposition_rate"))
    println(" aolean deposit   ",aoldep)
    println(" ocn balance      ",oro2ocn - ocnsed + runoff + aoldep)
end
    

