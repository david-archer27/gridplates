#

using SparseArrays
using LinearAlgebra
using Rotations, SharedArrays, StaticArrays
using Printf
using GLMakie
#inline!(true)
using Colors
using GeometryTypes
using Contour
using BSON
#
include("params.jl")
include("steps.jl")
include("plates.jl")
include("orogeny.jl")
include("land_sed.jl")
include("ocean_sed.jl")
include("caco3.jl")
include("isostacy.jl")
include("utilities.jl")
include("io.jl")

cd(code_base_directory)
rotations = read_rotation_file("1000_0_rotfile_Merdith_et_al.rot")
norotation = rotation2matrix(rotations[1])

orogenic_events = create_orogenies()

create_everything( earliesttime + time_step )
##
run_timeseries()
animate_all( earliesttime )
##
inventory_timeseries, diags_timeseries, frac_diags_timeseries, time_points = 
    tabulate_diagnostics( earliesttime )

scene = GLMakie.lines(time_points,inventory_timeseries[1,:])
##
GLMakie.lines!(time_points,plate_inventory_history)
GLMakie.lines!(time_points,cumulative_deposition_history)
GLMakie.lines!(time_points,cumulative_subduction_history)
balance = cumulative_deposition_history ./ time_step .- world_inventory_history .- cumulative_subduction_history ./ time_step
GLMakie.lines!(time_points,balance)

#
##

world = read_world(150)
plates = read_plates(150)
#
step_everything()
##
world_grid_sediment_change = world.sediment_thickness
step_everything()
world_grid_sediment_change = world_grid_sediment_change .* -1. .+ world.sediment_thickness
plot_field(world_grid_sediment_change)
#step_everything()
##
log_IO = open("../outfiles/logfile.txt","w")
while world.age > 0
    step_everything()
    save_world()
    if floor(world.age/5) == world.age/5 && world.age >= 10
        save_plates()
    end
end
close(log_IO)
##

step_everything()
##
create_everything(292)
time_step = 2.; sleep_length = 2
while world.age > 280
    old_sediment_thickness = deepcopy(world.sediment_thickness)
    world.age -= time_step
    println(world.age)

    #step_everything()
    
    world = read_world(world.age)
    #=scene = plot_elevation()
    scene = plot_add_continent_outlines!(scene)
    scene = plot_add_plate_boundaries!(scene)
    scene = plot_add_orogenies!(scene)=#
    
    ix = 190; iy = 130; nskirt = 40; sleep_length = 2
    uplift_rate = get_diag("continent_orogenic_uplift_rate") .+
        get_diag("subduction_orogenic_uplift_rate")
    crust_erosion_rate = get_diag("crust_erosion_rate")
    global_sediment_deposition_rate = get_diag("global_sediment_deposition_rate")
    #=
    println("uplift rate")
    scene = zoomplot(uplift_rate,ix,iy,nskirt,0.,600.)
    scene = zoom_add_continent_outlines!(scene,ix,iy,nskirt)
    scene = zoom_add_plate_boundaries!(scene,ix,iy,nskirt)
    scene = zoom_add_orogenies!(scene,ix,iy,nskirt)
    display(scene)
    sleep(sleep_length)=#
    println("sedimentation rate")
    scene = zoomplot(global_sediment_deposition_rate,ix,iy,nskirt,0.,600.)
    scene = zoom_add_continent_outlines!(scene,ix,iy,nskirt)
    scene = zoom_add_plate_boundaries!(scene,ix,iy,nskirt)
    scene = zoom_add_orogenies!(scene,ix,iy,nskirt)
    display(scene)
    sleep(sleep_length)
        #scene = zoomplot(crust_erosion_rate,ix,iy,nskirt,0.,600.)
    #=println("crust thickness")
    scene = zoomplot(world.crust_thickness,ix,iy,nskirt)
    scene = zoom_add_continent_outlines!(scene,ix,iy,nskirt)
    scene = zoom_add_plate_boundaries!(scene,ix,iy,nskirt)
    scene = zoom_add_orogenies!(scene,ix,iy,nskirt)
    display(scene)
    sleep(sleep_length)=#
    println("freeboard")
    scene = zoomplot(world.freeboard,ix,iy,nskirt,-100,2000)
    scene = zoom_add_continent_outlines!(scene,ix,iy,nskirt)
    scene = zoom_add_plate_boundaries!(scene,ix,iy,nskirt)
    scene = zoom_add_orogenies!(scene,ix,iy,nskirt)
    display(scene)
    sleep(sleep_length)
    #=println("erosion rate")
    scene = zoomplot(crust_erosion_rate,ix,iy,nskirt)
    scene = zoom_add_continent_outlines!(scene,ix,iy,nskirt)
    scene = zoom_add_plate_boundaries!(scene,ix,iy,nskirt)
    scene = zoom_add_orogenies!(scene,ix,iy,nskirt)  
    display(scene)
    sleep(sleep_length)=#
    println("sediment thickness")
    scene = zoomplot(world.sediment_thickness,ix,iy,nskirt)
    scene = zoom_add_continent_outlines!(scene,ix,iy,nskirt)
    scene = zoom_add_plate_boundaries!(scene,ix,iy,nskirt)
    scene = zoom_add_orogenies!(scene,ix,iy,nskirt)  
    display(scene)
    sleep(sleep_length)
    println("sediment thickness change")
    scene = zoomplot(world.sediment_thickness .- old_sediment_thickness,ix,iy,nskirt)
    scene = zoom_add_continent_outlines!(scene,ix,iy,nskirt)
    scene = zoom_add_plate_boundaries!(scene,ix,iy,nskirt)
    scene = zoom_add_orogenies!(scene,ix,iy,nskirt)  
    display(scene)
    sleep(sleep_length)
    #save_world()
    #save_plates()
end
##

# partial simple step
verbose = true
local initial_plate_inventories, update_ID_plate_inventories, 
    remasked_plate_inventories, remasked_world_inventories, 
    subduction_rates, subduction_balances,
    initial_world_inventories, filled_world_inventories,
    updated_world_inventories

world.age -= time_step
println()
println("beginning ", world.age, " Myr ", get_geo_interval())
if log_IO != 0
    println(log_IO,"beginning ", world.age, " Myr ", get_geo_interval())
end

if verbose == true
    initial_plate_inventories = global_plate_sediment_inventories()
    initial_world_inventories = world_sediment_inventories(  )
end

clear_world_process_arrays() 
# resets world.geomorphology, world.tectonics
#println("interpolating from plate grids")
oldIDmap = world.plateID
newIDmap = read_plateIDs(world.age)

world.plateID = newIDmap
world.plateIDlist = find_plateID_list()

update_changing_plateIDs(oldIDmap,newIDmap)
# look for plates that change their identities, move stuff between them
# sets plate.tectonics to record plateID changes

resolve_rotation_matrices()
# set rotation matrices to the new time
increment_plate_age()
# ocean and cont crust gets older in plate grids
 
if verbose == true
    update_ID_plate_inventories = global_plate_sediment_inventories()
    #old_plate_inventories = plate_sediment_inventories( plates[801] )
    logging_println()
    logging_println("update plateIDs init ", initial_plate_inventories)
    logging_println("               final ", update_ID_plate_inventories)
    logging_println("       plate inv bal ", ( update_ID_plate_inventories .- initial_plate_inventories ) ./
        time_step )
end

fill_world_from_plates() # build the new world grid
# world.geomorphology untouched, set initially
# nothing yet in world.tectonics
if verbose == true
    filled_world_inventories = world_sediment_inventories(  )
    logging_println()
    logging_println("fill world orig ", initial_world_inventories) 
    logging_println("         filled ", filled_world_inventories)
    logging_println("  cum world chg ", ( filled_world_inventories .- initial_world_inventories ) ./
        time_step )
end
update_world_continents_from_file()
# inflow of continent locations from gplates. initializes world crust_*
# fields only where the world crust type has changed from last step
# new continent points in world grid have h0 thickness for now
# sets world.tectonics where ocean <> continent changes
if verbose == true
    updated_world_inventories = world_sediment_inventories(  )
end

remask_plates()
# finds points of crust subduction and creation within plates fields.
# records plate subduction and creation in world diags.
# initializes plate crust_* variables at new points in plate fields.
# sets plate.tectonics for subduction, appearance
# new continent points on plate grid have h0 thickness

if verbose == true
    subduction_rates = fill(0.,0:n_sediment_types)
    subduction_rates[0] = 
        volume_field(get_diag("ocean_subduct_sediment_volume")) + 
        volume_field(get_diag("continent_subduct_sediment_volume")) 
    for i_sedtype in 1:n_sediment_types
        subduction_rates[i_sedtype] =
            volume_field(get_frac_diag("ocean_subduct_sediment_fraction_volume")[:,:,i_sedtype]) +
            volume_field(get_frac_diag("continent_subduct_sediment_fraction_volume")[:,:,i_sedtype])
    end
    remasked_plate_inventories = global_plate_sediment_inventories()
    plate_cum_change = ( remasked_plate_inventories .- initial_plate_inventories ) ./
        time_step
    world_cum_change = ( updated_world_inventories .- initial_world_inventories ) ./
        time_step 
    plate_bal = plate_cum_change .+ subduction_rates
    world_bal = world_cum_change .+ subduction_rates
    logging_println()
    logging_println(" remask plates init ", initial_plate_inventories)
    logging_println("          update ID ", update_ID_plate_inventories)
    logging_println("             remask ", remasked_plate_inventories)
    logging_println("          subducted ", subduction_rates)
    logging_println("      cum plate chg ", plate_cum_change )
    logging_println("      cum plate bal ", plate_bal )
    logging_println()
    logging_println("      cum world chg ", world_cum_change )
    logging_println("      cum world bal ", world_bal )
    # looking for world_bal = 0, so cum_change + subduction = 0,
    # (new' - old)/time_step + subduction = 0
    # new' = old - subduction * time_step    
    # new' / new = (old - subduction * time_step) / new
    imbalance_ratios = fill( 1., n_sediment_types )
    for i_sedtype in 1:n_sediment_types
        if updated_world_inventories[i_sedtype] > 0.
            imbalance_ratios[i_sedtype] = ( initial_world_inventories[i_sedtype] - 
                subduction_rates[i_sedtype] * time_step ) / updated_world_inventories[i_sedtype]
        end
    end

    #tweaked_inv = updated_world_inventories .- world_cum_change .* time_step
    #imbalance_ratios = tweaked_inv ./ updated_world_inventories
    
    scale_global_sediment_components( imbalance_ratios[1:n_sediment_types] )

    fudged_world_inventories = world_sediment_inventories(  )
    logging_println("  fudge factors ", imbalance_ratios)
    logging_println("         fudged ", fudged_world_inventories)
    logging_println("  cum world bal ", ( fudged_world_inventories .- initial_world_inventories ) ./
        time_step .+ subduction_rates)
end

world.elevation_offset = ocean_thermal_boundary_layer()
verbose = true

if verbose == true
    old_land_fraction_inventories = world_sediment_inventories( )
    #println("going in ", old_land_fraction_inventories)
end

clear_geomorph_process_arrays() # only necessary in standalone mode
orogeny()
# continent and subduction uplift rates into *_orogenic_uplift_rate diags
apply_orogeny_fluxes_to_world()

smooth_continental_crust_thickness()

#=if verbose == true
    old_land_fraction_inventories = world_sediment_inventories( )
    println("apply orogeny ", old_land_fraction_inventories)
end=#
isostacy() 
check_reburial_exposed_basement() 
#=if verbose == true
    old_land_fraction_inventories = world_sediment_inventories( )
    println("reburial ", old_land_fraction_inventories)
end=#
println();print("land")

# sets world.tectonics for reburial
original_elevation_field = land_transport_elevation()
ocean_sink = eq_mask( world.crust_type, ocean_crust )
subaereal_mask = eq_mask(world.geomorphology, sedimented_land)
new_elevation_field = original_elevation_field
# land bulk sediment transport
n_denuded = -1; n_first = true

##
while n_denuded != 0 
    #elevation_field = land_transport_elevation() # start over each time
    subaereal_mask = eq_mask(world.geomorphology, sedimented_land)
    # world.geomorphology was set in update_world_continents_from_file
    #subaereal_mask[:,1] .= 0; subaereal_mask[:,ny] .= 0
    generate_orogenic_erosion_fluxes()
    # sets land_orogenic_clay_flux, coastal_orogenic_clay_flux,
    # and crust_erosion_rate.
    aolean_transport()
    # resets and fills aolean_erosion_rate and aolean_deposition_rate
    subaereal_CaCO3_dissolution()
    # sets land surface CaCO3 dissolution rate 
    bulk_sediment_source = get_diag("land_orogenic_clay_flux") .- 
        get_diag("aolean_clay_erosion_rate") .-
        get_diag("land_CaCO3_dissolution_rate") 
    new_elevation_field = land_bulk_sediment_transport( original_elevation_field, 
        subaereal_mask, bulk_sediment_source, ocean_sink )
    # sets land_sediment_deposition_rate, land_sediment_fraction_deposition_rate, 
    # land_orogenic_clay_flux, coastal_orogenic_clay_flux, crust_erosion_rate.
    n_denuded = check_subaereal_exposure()

    #println(" accumulating trapped sediment ",volume_field(get_frac_diag("land_trapped_sediment_rate",1)))
    # sets world.geomorphology to exclude for next pass, 
    # world.tectonics for land_soil_denuded.
    # accumulates sediment that disappeared into exposed crust grid points in
    #   land_trapped_sediment_rate.  
    if n_denuded > 0
        if n_first == true
            print(" denuding ",n_denuded)
            n_first = false
        else
            print(", ",n_denuded)
        end
        #=if verbose == true # watch the erosion flux grow
            crust_clay_source = volume_field(get_diag("crust_clay_source_rate"))
            clay_src_land = volume_field(get_diag("land_orogenic_clay_flux"))
            clay_src_ocn = volume_field(get_diag("coastal_orogenic_clay_flux"))
            #land_tot_deposition = volume_field(get_diag("land_sediment_deposition_rate"))
            print(" crust clay source ",crust_clay_source," land ",clay_src_land," ocn ", clay_src_ocn)
        else
            println("")
        end=#
    end
end         
println("")
sediment_sources = fill(0.,nx,ny,n_sediment_types)
sediment_sources[:,:,clay_sediment] .= 
    get_diag("land_orogenic_clay_flux") .-
    get_diag("aolean_clay_erosion_rate")
sediment_sources[:,:,CaCO3_sediment] .= 
    - get_diag("land_CaCO3_dissolution_rate")

land_sediment_fraction_transport( new_elevation_field, 
    subaereal_mask, sediment_sources, ocean_sink ) 
# updates land_fraction_deposition_rates

#= if verbose == true
    land_clay_deposition = volume_field(
        get_frac_diag("land_sediment_fraction_deposition_rate",clay_sediment ))
    land_CaCO3_deposition = volume_field(
        get_frac_diag("land_sediment_fraction_deposition_rate",CaCO3_sediment ))
    # check the fractions add up
    println("land clay dep ",land_clay_deposition,
        " CaCO3 ",land_CaCO3_deposition,
        " c+C ", land_clay_deposition+land_CaCO3_deposition)
end  =#

set_land_runoff_fluxes( new_elevation_field, subaereal_mask, ocean_sink ) 
# fills coastal_sediment_fraction_runoff_flux
# sets world.geomorphology for coastal_depocenter
#println(" original runoff ",volume_field(get_frac_diag("coastal_sediment_fraction_runoff_flux",1)))


rescue_disappeared_land_sediment( )
#println(" augmented runoff ",volume_field(get_frac_diag("coastal_sediment_fraction_runoff_flux",1)))
# captures stuff in land_trapped_sediment_rate, adds it to 
# coastal_sediment_fraction_runoff_flux

if verbose == true
    old_land_fraction_inventories = world_sediment_inventories( )
    if old_land_fraction_inventories[1] != old_land_fraction_inventories[1]
        error("NaN detected")
    end
    #println("before update ", old_land_fraction_inventories)
end

world.sediment_thickness, world.sediment_surface_fractions = 
    apply_land_sediment_fluxes()

if verbose == true
    clay_src_land = volume_field(get_diag("land_orogenic_clay_flux"))
    sources = fill(0.,0:n_sediment_types); deposition = fill(0.,0:n_sediment_types)
    runoff = fill(0.,0:n_sediment_types); aolean = fill(0.,0:n_sediment_types)
    trapped = fill(0.,0:n_sediment_types)
    sources[1] = clay_src_land
    aolean[1] = volume_field(get_diag("aolean_clay_erosion_rate"))
    for i_sedtype in 1:n_sediment_types
        deposition[i_sedtype] = volume_field(
            get_frac_diag("land_sediment_fraction_deposition_rate",i_sedtype ))
        runoff[i_sedtype] = volume_field(
            get_frac_diag("coastal_sediment_fraction_runoff_flux",i_sedtype ))
        trapped[i_sedtype] = volume_field(
                get_frac_diag("land_trapped_sediment_rate",i_sedtype ))
        sources[0] += sources[i_sedtype]
        deposition[0] += deposition[i_sedtype]
        runoff[0] += runoff[i_sedtype]
        aolean[0] += aolean[i_sedtype]
        trapped[0] += trapped[i_sedtype]
    end
    flux_balances = sources .- deposition .- runoff .- aolean .+ trapped
    new_land_fraction_inventories = world_sediment_inventories( ) 
    change_rate = ( new_land_fraction_inventories .- old_land_fraction_inventories ) / 
        time_step
    inv_balances = change_rate - deposition
    logging_println()
    logging_println("land budget init ",old_land_fraction_inventories)
    logging_println("           final ",new_land_fraction_inventories)
    logging_println("         inv bal ",inv_balances)
    logging_println("   orogen source ",sources)
    logging_println("             dep ",deposition)
    logging_println("          runoff ",runoff)
    logging_println("         trapped ",trapped)
    logging_println("          aolean ",aolean)
    #logging_println("            sink ",runoff .+ aolean)
    logging_println("     change rate ",change_rate)
    logging_println("        flux bal ",flux_balances)
end
##
distribute_CaCO3_sedimentation( )

# sets coastal_CaCO3_flux, pelagic_CaCO3_flux
distribute_ocean_sediment_fluxes(  )
# accumulates seafloor_sediment_fraction_deposition_rate, seafloor_sediment_deposition_rate
# fills world.geomorphology for pelagic_seafloor, ocean_shelf
set_diag("global_sediment_deposition_rate", 
    get_diag("land_sediment_deposition_rate") .+ 
    get_diag("seafloor_sediment_deposition_rate") )
for i_sedtype in 1:n_sediment_types
    set_frac_diag("global_sediment_fraction_deposition_rate",i_sedtype,
        get_frac_diag("land_sediment_fraction_deposition_rate",i_sedtype) .+
        get_frac_diag("seafloor_sediment_fraction_deposition_rate",i_sedtype) )
end
sed_depo = get_diag("global_sediment_deposition_rate")
sed_frac_depo = get_frac_diag("global_sediment_fraction_deposition_rate")
for ix in 1:nx
    for iy in 1:ny
        if sed_depo[ix,iy] != 0.
            for i_sedtype in 1:n_sediment_types
                set_frac_diag("global_sediment_fraction_deposition_ratio",ix,iy,i_sedtype,
                    sed_frac_depo[ix,iy,i_sedtype] /
                    sed_depo[ix,iy])
            end
        end
    end
end

#=
if verbose == true
    # ocean clay balance
    ocn_clay_influx = volume_field(
        get_frac_diag("ocean_sediment_fraction_influx",clay_sediment))
    ocn_clay_dep = volume_field(get_frac_diag(
        "seafloor_sediment_fraction_deposition_rate",clay_sediment ))
    #println("clay ocean balance, in ", ocn_clay_influx," dep ",ocn_clay_dep,
    #    " bal ", ocn_clay_influx - ocn_clay_dep )
    # CaCO3 balance
    ocn_CaCO3_influx = volume_field(
        get_frac_diag("ocean_sediment_fraction_influx",CaCO3_sediment))
    ocn_CaCO3_dep = volume_field(get_frac_diag(
        "seafloor_sediment_fraction_deposition_rate",CaCO3_sediment ))
    #println("CaCO3 ocean balance, in ", CaCO3_influx," dep ",ocn_CaCO3_dep,
    #    " bal ", CaCO3_influx - ocn_CaCO3_dep )
    ocn_influx = [ ocn_clay_influx + ocn_CaCO3_influx,ocn_clay_influx,ocn_CaCO3_influx]
    ocn_dep = [ocn_clay_dep+ocn_CaCO3_dep,ocn_clay_dep,ocn_CaCO3_dep]
    balances = ocn_influx .- ocn_dep
end
=#
if verbose == true
    old_ocean_sed_inventories = world_sediment_inventories( )
end

apply_ocean_sediment_fluxes( )
# updates world.sediment* for ocean points from seafloor_sediment*_deposition_rate

if verbose == true
    new_ocean_sed_inventories = world_sediment_inventories( )
    total_sources = fill(0.,0:n_sediment_types); coastal_sources = fill(0.,0:n_sediment_types);
    pelagic_sources = fill(0.,0:n_sediment_types); orogenic_sources = fill(0.,0:n_sediment_types);
    ocean_sed_depo_rates = fill(0.,0:n_sediment_types)
    CaCO3_coastal_source = volume_field(get_diag("coastal_CaCO3_flux")) 
    CaCO3_pelagic_source =  volume_field(get_diag("pelagic_CaCO3_deposition_rate"))
    CaCO3_source = CaCO3_coastal_source + CaCO3_pelagic_source
    clay_runoff_source = volume_field( get_frac_diag("coastal_sediment_fraction_runoff_flux",clay_sediment) )
    clay_orogenic_source = volume_field( get_diag("coastal_orogenic_clay_flux") )
    clay_aolean_source = volume_field( get_diag("aolean_clay_deposition_rate") )
    clay_source = clay_runoff_source + clay_orogenic_source + clay_aolean_source
    total_sources[clay_sediment] = clay_source; total_sources[CaCO3_sediment] = CaCO3_source
    coastal_sources[clay_sediment] = clay_runoff_source; coastal_sources[CaCO3_sediment] = CaCO3_coastal_source
    pelagic_sources[clay_sediment] = clay_aolean_source; pelagic_sources[CaCO3_sediment] = CaCO3_pelagic_source
    orogenic_sources[clay_sediment] = clay_orogenic_source
    for i_sedtype in 1:n_sediment_types
        ocean_sed_depo_rates[i_sedtype] = volume_field(
            get_frac_diag("seafloor_sediment_fraction_deposition_rate",i_sedtype))
        total_sources[0] += total_sources[i_sedtype]
        coastal_sources[0] += coastal_sources[i_sedtype]
        pelagic_sources[0] += pelagic_sources[i_sedtype]
        orogenic_sources[0] += orogenic_sources[i_sedtype]
        ocean_sed_depo_rates[0] += ocean_sed_depo_rates[i_sedtype]
    end 
    inventory_change_rate = ( new_ocean_sed_inventories .- old_ocean_sed_inventories) /
        time_step
    balances = inventory_change_rate -
        ocean_sed_depo_rates 
    logging_println()
    logging_println("ocean budget init ", old_ocean_sed_inventories)
    logging_println("              new ", new_ocean_sed_inventories)
    logging_println("           change ", inventory_change_rate)
    logging_println("      tot sources ", total_sources)
    logging_println("  coastal sources ", coastal_sources)
    logging_println("          pelagic ", pelagic_sources)
    logging_println("         orogenic ", orogenic_sources)
    logging_println("              tot ", total_sources)
    logging_println("              dep ", ocean_sed_depo_rates)
    logging_println("              bal ", balances)
end                    

isostacy()


## big loop
#log_IO = open("../outfiles/logfile.txt","w")
ages = Float32[]
for age in range(earliesttime,step=-2.,stop=0.)
    push!(ages,age)
end

operation = "movie"  # "run" "watch" or "movie"
image_number = 0
for age in ages[2:end]
    println()

    if operation == "run"
        step_everything()
        if true == true
            save_world()
            if floor(age/5) == age/5 && true == true
                println("saving output bson files. ")
                save_plates()
            end
        end
    else
        println("reading ",age) 
        world = read_world(age)
    end

    if operation == "movie"  # make movie images
        if age < 535
            global image_number += 1
            #=
            scotese_elevation,scotese_age = nearest_scotese_elevation()
            scene = plot_two_fields(world.freeboard,scotese_elevation,-4000,4000)
            scotesetimestamp = string(Int(floor(scotese_age))) * " Myr"
            text!(scene, scotesetimestamp, position = (-180,-305),textsize=15)
            plot_add_timestamp!(scene,scotese_age,-180,-305)
            =#
            #scene = plot_sediment_thickness()
            #scene = plot_field(world.crust_age,0.,600.)
            scene = plot_field(world.freeboard,-5000.,5000)
            #scene = plot_elevation()
            #scene = plot_field(world.sediment_surface_fractions[:,:,2] .* 100.,0.,100.)
            plot_add_plate_boundaries!(scene)
            plot_add_orogenies!(scene)
            plot_add_continent_outlines!(scene)
            plot_add_timestamp!(scene,world.age,-180,-105)
            #plot_add_orogenies!(scene)
            outfilename = pwd() * "./../image_tmp/img." * lpad(image_number,3,"0") * ".png"
            #println(outfilename)
            Makie.save(outfilename,scene)
        end
        #ffmpeg -r 10 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p elevation.mp4
        #display(scene)
        #println(age," ",scotese_age, " ",get_geo_interval(age))
    end

    if operation == "watch"  # watch movie on screen
        scene = plot_sediment_thickness()
        plot_add_plate_boundaries!(scene)
        plot_add_orogenies!(scene)
        plot_add_continent_outlines!(scene)
        display(scene)
        #sleep(1)
        #println(highest_xy_field(world.freeboard))
    end

    if operation == "diags" # diagnostics


        new_ocean_sed_inventories = fill(0.,n_sediment_types+1)
        for ibin in 1:n_sediment_time_bins
            new_ocean_sed_inventories[1] += volume_field(
                world.sediment_layer_thickness[:,:,ibin] )
            for i_sedtype in 1:n_sediment_types
                new_ocean_sed_inventories[i_sedtype + 1] += volume_field(
                    world.sediment_layer_thickness[:,:,ibin] .* 
                    world.sediment_layer_fractions[:,:,i_sedtype,ibin] ) 
                if isnan(new_ocean_sed_inventories[2])
                    error("stopping")
                end
            end 
        end
        ocean_sed_depo_rates = [ volume_field(get_diag("seafloor_sediment_deposition_rate")) ]
        for i_sedtype in 1:n_sediment_types
            push!(ocean_sed_depo_rates,volume_field(
                get_frac_diag("seafloor_sediment_fraction_deposition_rate",i_sedtype)))
        end 



        weathering_input = volume_field(get_diag("crust_erosion_rate")) / time_step
        land_sediment_transport =
            volume_field(get_diag("sed_transport_deposition_rate")) / time_step
        land_aolean_transport = 
            volume_field( 
                get_diag("aolean_transport_deposition_rate") .*
                get_continent_mask()
                ) / time_step
        seafloor_aolean_transport = 
                volume_field( 
                    get_diag("aolean_transport_deposition_rate") .*
                    ( 1. .- get_continent_mask() )
                    ) / time_step
        seasurface_sediment_deposition =
            volume_field(get_diag("seasurface_deposition_rate")) / time_step
        seafloor_sediment_deposition =
            volume_field(get_diag("seafloor_deposition_rate")) / time_step
        seafloor_spilloff_rate = 
            volume_field(get_diag("seafloor_spilloff_rate")) / time_step
        sediment_subduction =
            areafield_total(get_diag("ocean_subduct_sediment_plate_volume")) / time_step
        world_sediment_volume = global_world_sediment_inventory()
        pos_freeb = clamp_field(world.freeboard,0.)
        mean_exposed_elevation = field_mean(pos_freeb)
        ocean_seds = world.sediment_thickness .*
            eq_mask(world.crust_type,ocean_crust)
        ocean_sediment_volume = volume_field(ocean_seds)
        #plate_sediment_volume = global_plate_sediment_inventory()
        @printf("%.2e %.2e %.2e %.2e %.2e %.1f\n", weathering_input,
            land_sediment_transport, seasurface_sediment_deposition, 
            seafloor_sediment_deposition, seafloor_spilloff_rate)
    end

end
#close(log_IO)
##
time_points = Float64[]
plate_inventory_history = Float64[]; world_inventory_history = Float64[]
deposition_history = [ Float64[],Float64[]]
subduction_history = [ Float64[],Float64[] ]
cumulative_depositions = fill(0.,2); cumulative_subductions = fill(0.,2)
cumulative_deposition_histories = [ Float64[],  Float64[] ]
cumulative_subduction_histories = [ Float64[], Float64[] ]

for age in ages[1:end-1]
    if floor(age/5) == age/5
        println(age)
        #global plates = read_plates(age)
        global world = read_world(age)
        #plate_inventory = global_plate_sediment_inventories()[1]
        world_inventory = world_sediment_inventories( )
        subduction_rates = fill(0.,0:n_sediment_types)
        deposition_rates = fill(0.,0:n_sediment_types)
        cumulative_depositions = fill(0.,0:n_sediment_types)
        cumulative_subductions = fill(0.,0:n_sediment_types)
        #subduction_rates[0] = volume_field(get_diag("ocean_subduct_sediment_volume")) + 
        #    volume_field(get_diag("continent_subduct_sediment_volume"))
        #deposition_rates[0] = volume_field(get_diag("global_sediment_deposition_rate"))
        for i_sedtype in 1:n_sediment_types
            subduction_rates[i_sedtype] = 
                volume_field(get_frac_diag("ocean_subduct_sediment_fraction_volume")[:,:,i_sedtype]) +
                volume_field(get_frac_diag("continent_subduct_sediment_fraction_volume")[:,:,i_sedtype])
            deposition_rates[i_sedtype] =
                 volume_field(get_frac_diag("global_sediment_fraction_deposition_rate",
                i_sedtype)) )
            cumulative_depositions[i_sedtype] += deposition_rates[i_sedtype] * time_step
            cumulative_subductions[i_sedtype] += subduction_rates[i_sedtype] * time_step
            push!(deposition_history,deposition_rate_total)
            push!(subduction_history,subduction_rate)
        end
        #push!(plate_inventory_history,plate_inventory)
        push!(time_points,-age)
        push!(world_inventory_history,world_inventory)
        push!(deposition_CaCO3_history,deposition_rate_CaCO3)
        push!(subduction_CaCO3_history,subduction_rate_CaCO3)
        push!(deposition_history,deposition_rate)
        push!(cumulative_deposition_history,cumulative_deposition)
        push!(cumulative_subduction_history,cumulative_subduction)
    end
end
scene = GLMakie.lines(time_points,world_inventory_history)
##
GLMakie.lines!(time_points,plate_inventory_history)
GLMakie.lines!(time_points,cumulative_deposition_history)
GLMakie.lines!(time_points,cumulative_subduction_history)
balance = cumulative_deposition_history ./ time_step .- world_inventory_history .- cumulative_subduction_history ./ time_step
GLMakie.lines!(time_points,balance)
#plot!(scene,ages,world_inventory_history)

## debug to

