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

if pwd() != "/Users/archer/Synched/papers/spongeball/codebase"
    cd("Synched/papers/spongeball/codebase")
end
rotations = read_rotation_file("1000_0_rotfile_Merdith_et_al.rot")
norotation = rotation2matrix(rotations[1])

orogenic_events = create_orogenies()
#
create_everything( earliesttime + time_step )
#
#world = read_world(132.)
#plates = read_plates(132.)
#
#old_sed_inv = global_world_sediment_inventory()
for i in 1:2
    step_everything()
    scene = plotfield(world.sediment_surface_fractions[:,:,2] .* 100.,0.,100.)
    plot_add_continent_outlines!(scene)
    plot_add_orogenies!(scene)
    display(scene)
    #total_sed_inv = global_world_sediment_inventory()
    #clay_src_land = volumefield_total(get_diag("land_orogenic_clay_flux"))
    #CaCO3_src = volumefield_total(
    #    get_frac_diag("ocean_sediment_fraction_influx",CaCO3_sediment))
    #balance = total_sed_inv - old_sed_inv - clay_src_land - CaCO3_src

    #println("total sed inv ", total_sed, " clay ",clay_src_land," CaCO3 ", CaCO3_src,
    #    " bal ", balance)
    #old_sed_inv = total_sed_inv
    #get_diag("land_sediment_deposition_rate")[333:337,145:149]
end
##
# partial simple step

# stuff in step_everything

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
print("land")
check_reburial_exposed_basement() 
original_elevation_field = land_transport_elevation()
new_elevation_field = fill(0.,nx,ny)
ocean_sink = generate_mask_field( world.crust_type, ocean_crust )

subaereal_mask = generate_mask_field(world.surface_type, sedimented_land)
#subaereal_mask[:,1] .= 0; subaereal_mask[:,ny] .= 0
# land bulk sediment transport
n_denuded = -1; n_first = true

#=
ocean_sink = fill(0.,nx,ny); ocean_sink[2,4] = 1
subaereal_mask = fill(0,nx,ny); subaereal_mask[2,2:3] .= 1
bulk_sediment_source = fill(0.,nx,ny); bulk_sediment_source[2,2] = 1.
original_elevation_field = fill(0.,nx,ny); original_elevation_field[2,2:3] .= 1.
new_elevation_field = land_bulk_sediment_transport( original_elevation_field, 
    subaereal_mask, bulk_sediment_source, ocean_sink )
#
sediment_sources = fill(0.,nx,ny,n_sediment_types)
sediment_sources[2,2,1] = 1.
land_sediment_fraction_transport( new_elevation_field, 
    subaereal_mask, sediment_sources, ocean_sink ) 
world.sediment_thickness, world.sediment_surface_fractions = 
    apply_land_sediment_fluxes()

original_elevation_field = deepcopy(new_elevation_field)
new_elevation_field = land_bulk_sediment_transport( original_elevation_field, 
    subaereal_mask, bulk_sediment_source, ocean_sink )
land_sediment_fraction_transport( new_elevation_field, 
    subaereal_mask, sediment_sources, ocean_sink ) 
world.sediment_thickness, world.sediment_surface_fractions = 
    apply_land_sediment_fluxes()

println(world.sediment_surface_fractions[2,1:4,2])
=#

while n_denuded != 0 
    #elevation_field = land_transport_elevation() # start over each time
    subaereal_mask = generate_mask_field(world.surface_type, sedimented_land)
    #subaereal_mask[:,1] .= 0; subaereal_mask[:,ny] .= 0
    generate_orogenic_erosion_fluxes()
    # sets land_orogenic_clay_flux, coastal_orogenic_clay_flux,
    # and crust_erosion_rate.
    aolean_transport()
    # resets and fills aolean_deposition_rate and aolean_deposition_rate
    subaereal_CaCO3_dissolution()
    # sets land surface CaCO3 dissolution rate 
    bulk_sediment_source = get_diag("land_orogenic_clay_flux") .- 
        get_diag("aolean_clay_erosion_rate")# .-
        #get_diag("land_CaCO3_dissolution_rate") 
    new_elevation_field = land_bulk_sediment_transport( original_elevation_field, 
        subaereal_mask, bulk_sediment_source, ocean_sink )
    # sets land_sediment_deposition_rate, land_sediment_fraction_deposition_rate, 
    # land_orogenic_clay_flux, coastal_orogenic_clay_flux, crust_erosion_rate.
    n_denuded = check_subaereal_exposure()
    # sets world.surface_type to exclude for next pass
    if n_denuded > 0
        if n_first == true
            print(" denuding ",n_denuded)
            n_first = false
        else
            print(", ",n_denuded)
        end
    end
end         
println("")
sediment_sources = fill(0.,nx,ny,n_sediment_types)
sediment_sources[:,:,clay_sediment] .= 
    get_diag("land_orogenic_clay_flux") .-
    get_diag("aolean_clay_erosion_rate")
sediment_sources[:,:,CaCO3_sediment] .= 
    - get_diag("land_CaCO3_dissolution_rate")

#subaereal_blobs = get_blobs( subaereal_mask )
#
#subaereal_blob = subaereal_blobs[1]
#world.sediment_thickness .= 1.
#world.sediment_surface_fractions[:,:,1] .= 1.
#world.sediment_surface_fractions[:,:,2] .= 0.
#subaereal_blob[1:120,:] .= 0#; subaereal_blob[300:360,:] .= 0
#subaereal_blob[:,1:4] .= 0
#land_area_sediment_fraction_transport( 
#    new_elevation_field, 
#    subaereal_blob, sediment_sources, ocean_sink ) 
#println(get_frac_diag("land_sediment_fraction_deposition_rate")[180,10:15,1])
#
#plotfield(get_diag("land_sediment_deposition_rate")[:,:])

#
land_sediment_fraction_transport( new_elevation_field, 
    #new_total_sediment_thickness,
    subaereal_mask, sediment_sources, ocean_sink ) 
world.sediment_thickness, world.sediment_surface_fractions = 
    apply_land_sediment_fluxes()
#

    # updates land_fraction_deposition_rates
#
set_land_runoff_fluxes( new_elevation_field, subaereal_mask, ocean_sink ) 
land_CaCO3_dep = get_frac_diag("land_sediment_fraction_deposition_rate",CaCO3_sediment )




# fills coastal_sediment_fraction_runoff_flux
clay_src_land = volumefield_total(get_diag("land_orogenic_clay_flux"))
land_clay_deposition = volumefield_total(
    land_CaCO3_dep)
land_CaCO3_deposition = volumefield_total(
    get_frac_diag("land_sediment_fraction_deposition_rate",CaCO3_sediment ))
land_clay_runoff = volumefield_total(
    get_frac_diag("coastal_sediment_fraction_runoff_flux",clay_sediment ))
land_CaCO3_runoff = volumefield_total(
    get_frac_diag("coastal_sediment_fraction_runoff_flux",CaCO3_sediment ))
# land clay mass balance
println("land clay balance, src ", clay_src_land,
    " dep ", land_clay_deposition," runoff ",land_clay_runoff,
    " bal ",land_clay_deposition + land_clay_runoff - clay_src_land )
# land CaCO3 mass balance
println("land CaCO3 balance, src ", 0.,
   " dep ", land_CaCO3_deposition," runoff ",land_CaCO3_runoff,
   " bal ",land_CaCO3_deposition + land_CaCO3_runoff)
#

distribute_CaCO3_sedimentation( )
#
# sets coastal_CaCO3_flux, pelagic_CaCO3_flux
distribute_ocean_sediment_fluxes(  )
# accumulates seafloor_sediment_fraction_deposition_rate, seafloor_sediment_deposition_rate
#
world.sediment_thickness, world.sediment_surface_fractions = 
    apply_land_sediment_fluxes()

apply_ocean_sediment_fluxes( )
 
set_diag("global_sediment_deposition_rate", 
    get_diag("land_sediment_deposition_rate") .+ 
    get_diag("land_sediment_deposition_rate") )
for i_sedtype in 1:n_sediment_types
    set_frac_diag("global_sediment_fraction_deposition_rate",i_sedtype,
        get_frac_diag("land_sediment_fraction_deposition_rate",i_sedtype) .+
        get_frac_diag("seafloor_sediment_fraction_deposition_rate",i_sedtype) )
end
for ix in 1:nx
    for iy in 1:ny
        if get_diag("global_sediment_deposition_rate")[ix,iy] != 0.
            for i_sedtype in 1:n_sediment_types
                set_frac_diag("global_sediment_fraction_deposition_ratio",ix,iy,i_sedtype,
                    get_frac_diag("global_sediment_fraction_deposition_rate",i_sedtype)[ix,iy] /
                    get_diag("global_sediment_deposition_rate")[ix,iy])
            end
        end
    end
end
#scene = plotfield(world.sediment_surface_fractions[:,:,2] .* 100.)
#plot_add_continent_outlines!(scene)

land_dep = get_diag("land_sediment_deposition_rate")[180:184,103:108]
land_C_dep = get_frac_diag("land_sediment_fraction_deposition_rate")[180:184,103:108,1]


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

ages = Float32[]
for age in range(earliesttime,step=-2.,stop=130.)
    push!(ages,age)
end
operation = "run"  # "run" or "movie"
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

    if true == false  # watch movie on screen
        scene = plotfield(world.sediment_fractions[:,:,current_time_bin(),1],0,1)
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

