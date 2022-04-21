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

include("params.jl")
include("steps.jl")
include("plates.jl")
include("orogeny.jl")
include("land_sed.jl")
include("ocean_sed.jl")
include("seafloor.jl")
include("utilities.jl")
include("io.jl")

earliesttime = 100.; time_step = 2.

if pwd() != "/Users/archer/Synched/papers/spongeball/codebase"
    cd("Synched/papers/spongeball/codebase")
end
rotations = read_rotation_file("1000_0_rotfile_Merdith_et_al.rot")
norotation = rotation2matrix(rotations[1])

orogenic_events = create_orogenies()

create_everything( earliesttime + time_step )
#
step_everything()

## partial simple step

world.age -= time_step
clear_world_process_arrays() 
resolve_rotation_matrices()
# set rotation matrices to the new time
increment_plate_age()
# ocean and cont crust gets older in plate grids
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
# initializes plate crust_* variables at new points in plate fields.
# new continent points on plate grid have h0 thickness
ocean_thermal_boundary_layer()

clear_geomorph_process_arrays() # only necessary in standalone mode
# world. files only
orogeny()
# continent and subduction uplift rates into *_orogenic_uplift_rate diags
apply_orogeny_fluxes_to_world()
isostacy() 
update_surface_types() 

original_elevation_field = fill(0.,nx,ny)
subaereal_mask = fill(0.,nx,ny); subaereal_mask[2,2:3] .= 1.
orogenic_source = fill(0.,nx,ny);  orogenic_source[2,2] = 1.
world.sediment_thickness = fill(0.,nx,ny); world.sediment_thickness[2,2:3] .= 1.
world.sediment_fractions = fill(0.,nx,ny,2); world.sediment_fractions[2,2:4,2] .= 1.
set_diag("land_orogenic_clay_flux",orogenic_source)
ocean_sink = fill(0.,nx,ny);  ocean_sink[2,4] = 1

new_elevation_field = land_sediment_transport( original_elevation_field, 
    subaereal_mask, orogenic_source, ocean_sink )
clay_src_land = volumefield_total(get_diag("land_orogenic_clay_flux"))
land_tot_deposition = volumefield_total(get_diag("land_sediment_deposition_rate"))
#get_diag("land_sediment_deposition_rate")[1:3,1:5]

land_sediment_fraction_transport( new_elevation_field, 
    subaereal_mask, orogenic_source, ocean_sink ) 
get_frac_diag("land_sediment_fraction_deposition_rate")[1:3,1:5,1] .+
get_frac_diag("land_sediment_fraction_deposition_rate")[1:3,1:5,2]
#
land_clay_deposition = volumefield_total(
    get_frac_diag("land_sediment_fraction_deposition_rate",clay_sediment ))
land_CaCO3_deposition = volumefield_total(
    get_frac_diag("land_sediment_fraction_deposition_rate",CaCO3_sediment ))
#
set_land_runoff_fluxes( new_elevation_field, subaereal_mask, ocean_sink ) 
land_clay_runoff = volumefield_total(
    get_frac_diag("coastal_sediment_fraction_runoff_flux",clay_sediment ))
land_CaCO3_runoff = volumefield_total(
    get_frac_diag("coastal_sediment_fraction_runoff_flux",CaCO3_sediment ))
println("land tot dep ",land_tot_deposition," clay ",land_clay_deposition,
    " CaCO3 ",land_CaCO3_deposition)
println("land clay src ", clay_src_land,
    " dep ", land_clay_deposition," runoff ",land_clay_runoff,
    " d+r ",land_clay_deposition + land_clay_runoff)
println("land CaCO3 src ", 0.,
    " dep ", land_CaCO3_deposition," runoff ",land_CaCO3_runoff,
    " d+r ",land_CaCO3_deposition + land_CaCO3_runoff)
##

##
land_clay_deposition = volumefield_total(get_frac_diag("land_sediment_fraction_deposition_rate",clay_sediment ))
land_CaCO3_deposition = volumefield_total(get_frac_diag("land_sediment_fraction_deposition_rate",CaCO3_sediment ))
land_tot_deposition = volumefield_total(get_diag("land_sediment_deposition_rate"))
land_clay_runoff = volumefield_total(get_frac_diag("coastal_sediment_fraction_runoff_flux",clay_sediment ))
land_CaCO3_runoff = volumefield_total(get_frac_diag("coastal_sediment_fraction_runoff_flux",CaCO3_sediment ))
ocn_clay_deposition = volumefield_total(get_frac_diag("seafloor_sediment_fraction_deposition_rate",clay_sediment ))
ocn_CaCO3_deposition = volumefield_total(get_frac_diag("seafloor_sediment_fraction_deposition_rate",CaCO3_sediment ))
clay_balance = crust_clay_source - land_clay_deposition - ocn_clay_deposition
clay_land_balance = clay_src_land - land_clay_deposition - land_clay_runoff
clay_runoff_balance = land_clay_runoff - ocn_clay_deposition
CaCO3_balance = - land_CaCO3_deposition - ocn_CaCO3_deposition
CaCO3_runoff_balance = land_CaCO3_runoff - ocn_CaCO3_deposition




## big loop

lasttime = 720.
ntimeslices = ( earliesttime - lasttime ) /
    time_step + 1
ages = Float32[]
for i in 1:ntimeslices
    age = earliesttime - time_step * (i-1)
    push!(ages,age)   # eg [50 45 40 35 30 25 20 15 10 5 0]
end

operation = "run"  # "run" or "movie"
image_number = 0
for age in ages[2:end]
    println(world.age, " ", get_geo_interval())

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
            outfilename = pwd() * "../image_tmp/img." * lpad(image_number,3,"0") * ".png"
            #println(outfilename)
            Makie.save(outfilename,scene)
        end
        #ffmpeg -r 10 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p elevation.mp4
        #display(scene)
        #println(age," ",scotese_age, " ",get_geo_interval(age))
    end

    if true == true  # watch movie on screen
        scene = plotfield(world.freeboard,-5000,5000)
        display(scene)
        #println(highest_xy_field(world.freeboard))
    end

    if true == false # diagnostics
        weathering_input = volumefield_total(get_diag("crust_erosion_rate")) / time_step
        land_sediment_transport =
            volumefield_total(get_diag("sed_transport_deposition_rate")) / time_step
        land_aolean_transport = 
            volumefield_total( 
                get_diag("aolean_transport_deposition_rate") .*
                get_continent_mask()
                ) / time_step
        seafloor_aolean_transport = 
                volumefield_total( 
                    get_diag("aolean_transport_deposition_rate") .*
                    ( 1. .- get_continent_mask() )
                    ) / time_step
        seasurface_sediment_deposition =
            volumefield_total(get_diag("seasurface_deposition_rate")) / time_step
        seafloor_sediment_deposition =
            volumefield_total(get_diag("seafloor_deposition_rate")) / time_step
        seafloor_spilloff_rate = 
            volumefield_total(get_diag("seafloor_spilloff_rate")) / time_step
        sediment_subduction =
            areafield_total(get_diag("ocean_subduct_sediment_plate_volume")) / time_step
        world_sediment_volume = global_world_sediment_inventory()
        pos_freeb = clamp_field(world.freeboard,0.)
        mean_exposed_elevation = field_mean(pos_freeb)
        ocean_seds = world.sediment_thickness .*
            generate_mask_field(world.crust_type,ocean_crust)
        ocean_sediment_volume = volumefield_total(ocean_seds)
        #plate_sediment_volume = global_plate_sediment_inventory()
        @printf("%.2e %.2e %.2e %.2e %.2e %.1f\n", weathering_input,
            land_sediment_transport, seasurface_sediment_deposition, 
            seafloor_sediment_deposition, seafloor_spilloff_rate)
    end

end


## debug to

