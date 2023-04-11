# Isostacy
function isostacy( )
    lithosphere_thickness = fill(0.,nx,ny)
    new_freeboard = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            new_freeboard[ix,iy], lithosphere_thickness[ix,iy], new_crust_density = 
                isostacy( 
                    world.crust_type[ix,iy], 
                    world.crust_thickness[ix,iy],
                    world.crust_density[ix,iy], 
                    world.sediment_thickness[ix,iy],
                    world.crust_age[ix,iy])
            if world.crust_type[ix,iy] == ocean_crust
                world.crust_density[ix,iy] = new_crust_density
            end
        end
    end
    world.freeboard[:,:] = new_freeboard[:,:]
    set_diag("lithosphere_thickness",lithosphere_thickness)
end
function isostacy( crust_type, crust_thickness, crust_density, 
    sediment_thickness, crust_age )

    column_thickness_ref = continent_crust_h0 + continent_lithosphere_thickness
    column_mass_ref = continent_crust_h0 * rho_continent_crust +
        continent_lithosphere_thickness * rho_lithosphere
    column_ref_mean_density = column_mass_ref / column_thickness_ref
    column_mass_ref_in_water = ( column_ref_mean_density - rho_seawater ) *
        column_thickness_ref
    elevation = 0; lithosphere_thickness = 0;
    if crust_type == continent_crust
        #tot_cap_thickness = crust_thickness +
        #    sediment_thickness
        lithosphere_thickness = continent_lithosphere_thickness
        #uncompacted_sediment_thickness = sediment_thickness / 
        #    ( 1. - porosity_continental_sediment )
        column_thickness = crust_thickness +
            sediment_thickness + lithosphere_thickness
        column_mass = crust_thickness * crust_density +
            lithosphere_thickness * rho_lithosphere +
            sediment_thickness * rho_bulk_continental_sediment
        aestheno_thickness_to_ref = 
            ( column_mass_ref - column_mass) / 
            rho_mantle
        thickness_relative_to_ref = column_thickness +
            aestheno_thickness_to_ref
        elevation = thickness_relative_to_ref - 
            column_thickness_ref - 
            world.sealevel
        if elevation < 0.
            column_mean_density = column_mass / column_thickness
            column_mass_in_water = ( column_mean_density - rho_seawater ) * 
                column_thickness
            aestheno_thickness_to_ref = 
                ( column_mass_ref_in_water - column_mass_in_water ) / 
                ( rho_mantle - rho_seawater )
            thickness_relative_to_ref = column_thickness + 
                aestheno_thickness_to_ref
            elevation = thickness_relative_to_ref - 
                column_thickness_ref - 
                world.sealevel
        end
    end
    if crust_type == ocean_crust
        thermal_boundary_thickness = 2. * sqrt( thermal_diffusivity * crust_age )
        thermal_boundary_thickness = min(thermal_boundary_thickness, max_lithosphere_thickness)
        lithosphere_thickness = thermal_boundary_thickness - crust_thickness
        lithosphere_thickness = max(0.,lithosphere_thickness)
        thickness_to_bottom_of_lithosphere = ocean_crust_h0 + lithosphere_thickness
        crust_temperature = crust_thickness / 
            thickness_to_bottom_of_lithosphere * mantle_T0 / 2.
        crust_density = rho_ocean_crust * 
            ( 1. - crust_temperature * thermal_expansion )
        #uncompacted_sediment_thickness = sediment_thickness /
        #    ( 1.0 - porosity_ocean_sediment )
        column_thickness = crust_thickness +
            sediment_thickness + lithosphere_thickness 
        column_mass = crust_thickness * crust_density +
            sediment_thickness * rho_bulk_ocean_sediment +
            lithosphere_thickness * rho_lithosphere
        column_mean_density = column_mass / column_thickness
        column_mass_in_water = 
            crust_thickness * ( crust_density - rho_seawater ) +
            sediment_thickness * ( rho_bulk_ocean_sediment - rho_seawater ) +
            lithosphere_thickness * ( rho_lithosphere - rho_seawater )
        aestheno_thickness_to_ref = 
            ( column_mass_ref_in_water - column_mass_in_water ) / 
            ( rho_mantle - rho_seawater )
        thickness_relative_to_ref = column_thickness + 
            aestheno_thickness_to_ref
        elevation = thickness_relative_to_ref - 
            column_thickness_ref - 
            world.sealevel
    end
    return elevation, lithosphere_thickness, crust_density
end




function isostacy_orig()
    (world.freeboard) =
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
#function isostacy_point(crust_thickness,crust_density,sediment_thickness,
#    crust_age,)
#=function isostacy_point_old(crust_thickness,crust_density,sediment_thickness,
    elevation_offset)
    surface_boundary_mass =
        crust_thickness * crust_density +
        sediment_thickness * rho_sediment
    surface_boundary_thickness =
        crust_thickness + 
        sediment_thickness
    surface_boundary_density = surface_boundary_mass /
        surface_boundary_thickness
    root_over_foot_subaereally = surface_boundary_density /
        ( rho_mantle - surface_boundary_density )
    elevation_rel_mantle = surface_boundary_thickness /
        ( 1. + root_over_foot_subaereally )
    elevation_relative_cont_0 = elevation_rel_mantle - 
        reference_elevation_cont_0
    freeboard = elevation_relative_cont_0 -
        world.sealevel 
    if freeboard < 0
        submerged_root_over_foot = 
            ( surface_boundary_density - rho_seawater ) / 
            ( rho_mantle - surface_boundary_density )
        freeboard *= 
            ( 1. + root_over_foot_subaereally ) / 
            ( 1. + submerged_root_over_foot )
    end
    return freeboard
end=#
function isostacy(crust_thickness,crust_density,sediment_thickness,
    elevation_offset) # sets the elevation of the solid surface rel to mantle line
    # time independent
    #surface_elevation = fill(0.,nx,ny)
    freeboard = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            #surface_elevation[ix,iy], 
            freeboard[ix,iy] = 
                isostacy_point( crust_thickness[ix,iy], crust_density[ix,iy],
                    sediment_thickness[ix,iy], elevation_offset[ix,iy] )
        end
    end
    return freeboard
end
function hack_crust_thickness_to_elevation( ix,iy,target_elevation )
    # crust_thickness, elevation, crust_density, sediment_thickness, elevation_offset )
    this_crust_freeboard_expression = 1. - world.crust_density[ix,iy] / rho_mantle
    elevation_misfit = world.freeboard[ix,iy] - target_elevation
    world.crust_thickness[ix, iy] += elevation_misfit *
        this_crust_freeboard_expression
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
#=function ocean_thermal_boundary_layer( )
    elevation_offset = ocean_thermal_boundary_layer(world.crust_age,
        world.crust_thickness,world.crust_density)
    return elevation_offset
end=#

function ocean_thermal_boundary_layer( ) # crust_age,crust_thickness,crust_density) # sets world.elevation_offset etc for ocn
    # elevation_offset could also be used to tweak continents up and down
    # time independent, just based on crust age
    kappa = 1.e-6
    thermal_expansion = 3.3e-5
    #elevation_offset = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if world.crust_type[ix,iy] == ocean_crust
                age = world.crust_age[ix,iy] * 1.e6  # years
                boundary_thickness = sqrt( kappa * age * 3.14e7) # meters, inc crust
                boundary_thickness = max(boundary_thickness,world.crust_thickness[ix,iy])
                mantle_boundary_thickness = boundary_thickness - world.crust_thickness[ix,iy]
                temperature_crust_base = mantle_T0
                if mantle_boundary_thickness > 0.
                    #temperature_crust_base = mantle_T0 -
                    #    ( mantle_T0 - ocean_T0 ) *
                    #    crust_thickness[ix,iy] / boundary_thickness
                    temperature_crust_base = mantle_T0 * world.crust_thickness[ix,iy] /
                        boundary_thickness
                end
                mantle_bl_temp_avg = ( mantle_T0 + temperature_crust_base ) / 2.
                crust_temp_avg = ( ocean_T0 + temperature_crust_base ) / 2.
                world.crust_density[ix,iy] = rho_ocean_crust *
                    ( 1. - thermal_expansion * ( crust_temp_avg - mantle_T0 / 2. ) )
                #crust_density[ix,iy] = rho_ocean_crust * 
                #    ( 1. + thermal_expansion * ( crust_temp_avg - ocean_T0 ) )
                mantle_bl_density = rho_mantle *
                    ( 1. - thermal_expansion * ( mantle_bl_temp_avg - mantle_T0 ) )
                crust_elevation_offset = world.crust_thickness[ix,iy] * 
                    ( 1. - rho_ocean_crust / world.crust_density[ix,iy] )
                mantle_elevation_offset = mantle_boundary_thickness *
                    ( 1. - rho_mantle / mantle_bl_density )
                world.elevation_offset[ix,iy] = - crust_elevation_offset - mantle_elevation_offset
            end
        end
    end
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

function get_sealevel( age ) 
    sea_level = 0.
    if enable_sealevel_change
        sea_level = get_interpolated_time_value( sealevel_values, sealevel_timepoints, age )
    end
    return sea_level
end
function get_atmCO2( age )
    atmCO2 = 0.
    if enable_atmCO2_change
        get_interpolated_time_value( atmCO2_values, atmCO2_timepoints, age )
    end
    return atmCO2
end

function get_interpolated_time_value( values, timepoints, age )
    ibelow = search_sorted_below( timepoints, age )  # eg 542, 450, 250,
    n_timepoints = length( timepoints )
    fractions = fill( 0., n_timepoints )
    if ibelow == 0
        #println("weird age ", age)
        fractions[n_timepoints] = 1.
    elseif ibelow == 1 # age is older than first 
        fractions[ibelow] = 1.
    else
        iabove = ibelow - 1
        fractions[iabove] = (age - timepoints[ibelow]) /
            (timepoints[iabove] - timepoints[ibelow])
        fractions[ibelow] = 1. - fractions[iabove]
    end
    interpolated_value = 0. 
    for i_timepoint in 1:n_timepoints
        interpolated_value += fractions[i_timepoint] * values[i_timepoint]
    end
    return interpolated_value
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
function land_area()
    land_area = 0.
    for ix in 1:nx
        for iy in 1:ny
            if world.freeboard[ix,iy] > 0.
                land_area += areabox[iy] # m2
            end
        end
    end
    return land_area
end


