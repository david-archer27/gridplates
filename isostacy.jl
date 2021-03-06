# Isostacy
function isostacy()
    (world.surface_elevation, world.freeboard) =
        isostacy( world.crust_thickness, world.crust_density,
        world.sediment_thickness, world.elevation_offset )
end
function muddy_isostatic_freeboard( sediment_thickness_change )
    new_sediment_thickness = sediment_thickness_change .+ world.sediment_thickness
    surface_elevation, freeboard = isostacy( world.crust_thickness,
        world.crust_density, new_sediment_thickness,
        world.elevation_offset)
    return freeboard
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
function crust_level()
    surface_boundary_mass =
        continent_crust_h0 * rho_continent_crust
    surface_boundary_thickness =
        continent_crust_h0
    equivalent_mantle_thickness =
        surface_boundary_mass /
        rho_mantle
    crust_level = surface_boundary_thickness -
        equivalent_mantle_thickness
    return crust_level
end
function calculate_ocean_crust_depth_unused(agefield) # from Stein and Stein 1992
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
            if world.crust_type[ix,iy] == ocean_crust
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
function compute_freeboard()
    for ix in 1:nx
        for iy in 1:ny
            world.freeboard[ix,iy] = world.surface_elevation[ix,iy] -
                world.sealevel
        end
    end
    return
end
# Utilities
function get_sealevel() # will be replaced by interpolation like for rotations
    return 4714.7
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


