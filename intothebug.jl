using SparseArrays
using LinearAlgebra
using Rotations, SharedArrays, StaticArrays
using Printf
using GLMakie
using Colors
using GeometryTypes
using Contour
using BSON
using StaticArraysCore
using Plots

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

#enable_crust_cont2ocn = false

cd( base_directory * "/" * code_base_directory )
rotations = read_rotation_file("1000_0_rotfile_Merdith_et_al.rot")
norotation = rotation2matrix(rotations[1])
orogenic_events = create_orogenies()
create_everything( earliesttime ) #+ time_step )
setup_working_directories( )

step_tectonics( true )
world.elevation_offset = ocean_thermal_boundary_layer()
world.crust_age .+= time_step


if enable_orogeny
    orogeny()
end

apply_orogeny_fluxes_to_world()

smooth_continental_crust_thickness()

isostacy() 
    ocean_sink = eq_mask( world.crust_type, ocean_crust )
    if enable_land_runoff == false
        ocean_sink .= 0
    end
       original_elevation_field = land_transport_elevation()
 sediment_source_fields = fill(0.,nx,ny,0:n_sediment_types)
    orogenic_sediment_source_fields = fill(0.,nx,ny,0:n_sediment_types)
    denuded_sediment_source_fields = fill(0.,nx,ny,0:n_sediment_types)

  diffusive_mask = eq_mask(world.geomorphology, sedimented_land)
        orogenic_mask = eq_mask(world.geomorphology, exposed_basement)
 
        orogenic_sediment_source_fields[:,:,clay_sediment] = generate_orogenic_erosion_fluxes() .* is_land()
        update_flux_totals!(orogenic_sediment_source_fields)
        orogenic_boundary_source_fields = 
            distribute_fluxes_uniformly_outside_boundary( orogenic_sediment_source_fields, 
                orogenic_mask )

        denuded_sediment_source_fields = get_denuded_sediment_fluxes()
        update_flux_totals!(denuded_sediment_source_fields)
        #println("denuded src ", volume_fields(denuded_sediment_source_fields))
        denuded_boundary_source_fields = 
            distribute_fluxes_uniformly_outside_boundary( denuded_sediment_source_fields, 
                orogenic_mask )
        #denuded_land_boundary_source_fields = denuded_boundary_source_fields .* is_land()
        #println(" = denuded bdy? ", volume_fields(denuded_boundary_source_fields))
 
        aolean_erosion_field, aolean_deposition_field = aolean_transport()
        land_CaCO3_dissolution_rate_field = subaereal_CaCO3_dissolution()

        sediment_source_fields = ( orogenic_boundary_source_fields .+
            denuded_boundary_source_fields ) .* is_land()
        sediment_source_fields[:,:,clay_sediment] -= aolean_erosion_field 
        sediment_source_fields[:,:,CaCO3_sediment] -= land_CaCO3_dissolution_rate_field
        update_flux_totals!( sediment_source_fields ) 
        #println(" sediment_source_fields ",sediment_source_fields)

        new_elevation_field, land_sediment_deposition_rate_field = 
            land_bulk_sediment_transport( original_elevation_field, 
                diffusive_mask, sediment_source_fields[:,:,0], ocean_sink )









