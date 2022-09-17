function find_steady_state_ocean_CO3( )
    local ocean_CO3, total_deposition
    terrestrial_CaCO3_dissolution = volume_field(
        get_diag("land_CaCO3_dissolution_rate") ) # m3 / Myr
    silicate_weathering_CaCO3_rate =
        global_CO2_degassing_rate * # mol / yr
        100. * 1.e6 / rho_sediment / 1.e6 # m3 / Myr
    target_deposition = 1.e15 # today fluxes
    #silicate_weathering_CaCO3_rate # + terrestrial_CaCO3_dissolution
    
    # ground truth stuff for diagnostics
    # from Cartapanis 2018
    coastal_deposition_today = 150.e15 / # g C / kyr 
        12. * 100. * # g CaCO3 / yr 
        1.e3 / # g CaCO3 / Myr
        rho_sediment / 1.e6 # m3 / Myr
    pelagic_deposition_today = 130.e15 / # g C / kyr 
        12. * 100. * # g CaCO3 / yr 
        1.e3 / # g CaCO3 / Myr
        rho_sediment / 1.e6 # m3 / Myr
    # from Makenzie 2004 BGD
    terrestrial_CaCO3_dissolution_today = 11.e12 * # mol / yr
        100. * # g CaCO3 / yr
        1.e6 / # g CaCO3 / Myr
        rho_sediment / 1.e6 # m3 / Myr
    another_silicate_weathering_CaCO3_rate = 21.e12 * # mol / yr
        100. * # g CaCO3 / yr
        1.e6 / # g CaCO3 / Myr
        rho_sediment / 1.e6 # m3 / Myr
    # higher than balance with degassing, acc. Mackenzie 

    ocean_CO3_low = 80.
    ocean_CO3_high = 100.
    imbalance = 1.e12; n_iters = 0
    depo_low = get_total_CaCO3_deposition_rate( ocean_CO3_low ) 
    depo_high = get_total_CaCO3_deposition_rate( ocean_CO3_high ) 
    while n_iters < 100 && imbalance > 1.e6
    #for iter in 1:10
        n_iters += 1
        ocean_CO3 = ocean_CO3_low + 
            ( target_deposition - depo_low ) *
            ( ocean_CO3_high - ocean_CO3_low ) / 
            ( depo_high - depo_low )
        total_deposition = get_total_CaCO3_deposition_rate( ocean_CO3 )
        if total_deposition > target_deposition # too much, 
            ocean_CO3_high = ocean_CO3; depo_high = total_deposition
        else
            ocean_CO3_low = ocean_CO3; depo_low = total_deposition
        end
        imbalance = abs( total_deposition - target_deposition )
        #println(" CO3 ",ocean_CO3," " )
    end
    println()
    println("  CaCO3 iters ", n_iters, " tot ", total_deposition, 
    " CO3 ",ocean_CO3," bal ",imbalance)
    return ocean_CO3
end
function get_total_CaCO3_deposition_rate( ocean_CO3 )
    coastal_CaCO3_deposition = volume_field( 
        get_coastal_CaCO3_deposition_field( ocean_CO3 ) )
    pelagic_CaCO3_deposition = volume_field( 
        get_pelagic_CaCO3_deposition_field( ocean_CO3 ) )
    continental_CaCO3_deposition = volume_field( 
        get_continental_CaCO3_deposition_field( ocean_CO3 ) )
    total_deposition = coastal_CaCO3_deposition + 
        pelagic_CaCO3_deposition + continental_CaCO3_deposition
    #imbalance = abs( total_deposition - target_deposition )
    #=println("ocean CO3 ",ocean_CO3,
        " coastal ",coastal_CaCO3_deposition,
        " pelagic ", pelagic_CaCO3_deposition,
        " continental ",continental_CaCO3_deposition,
        " tot ", total_deposition)
        =#
    return total_deposition
end
function get_continental_CaCO3_deposition_field( ocean_CO3 ) 
    shallow_rates = fill(0.,nx,ny)
    #= for ix in 1:nx
        for iy in 1:ny
            if world.crust_type[ix,iy] == continent_crust
                if world.freeboard[ix,iy] < shelf_depth_CaCO3
                    accom_meters = ( shelf_depth_CaCO3 - world.freeboard[ix,iy] ) /
                        sediment_freeboard_expression
                    shallow_rates[ix,iy] = accom_meters / time_step
                end
            end
        end
    end =#
    return shallow_rates
end
function get_pelagic_CaCO3_deposition_field( ocean_CO3 )
    pelagic_rates = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if world.geomorphology[ix,iy] <= coastal_depocenter # or pelagic
                depth = - world.freeboard[ix,iy]
                delta_CO3 = ocean_CO3 - co3_saturation( depth )
                set_diag("seafloor_delta_CO3",ix,iy,delta_CO3)
                pelagic_rates[ix,iy] = get_pelagic_deposition_rate( ocean_CO3, delta_CO3, iy )
            end
        end
    end
    return pelagic_rates
end
function co3_saturation( depth )
    #    Calcite solubility, from sayles, and from millero
    t_k = 275.   
    kprime = 4.75E-7        #                mol2 / kg2
    delv = -44.              #                cm3 / mol
    rr = 83.14              #                cm3 bar / K mol
    dk = -.0133             #                cm3 / bar mol
    pres = depth/10.         #                bar
    kpres = log(kprime) -
        delv / ( rr * t_k ) * pres +
        0.5 * dk / (rr * t_k ) * pres^2
    kpres = exp(kpres)
    csat = kpres / 0.01
    return csat * 1.e6 # micromol/l 
end
function get_pelagic_deposition_rate( ocean_CO3, delta_CO3, iy )
    latitude_scale = CaCO3_latitude_scale( iy )
    #depo = 4. * ( delta_CO3 / 40. )^2. # meters / Myr
    burial_max = 1.e-2 * ocean_CO3 * 
        exp( ( ocean_CO3 - 100. ) / 500. ) * 10. 
        # gives 1 cm/kyr (10 m/Myr) when ocean_CO3 = 100, 0.22 @ 25, 2.44 @ 200
    preservation_scale = ( atan( delta_CO3 / 5. ) + pi/2. ) / pi
    deposition = burial_max * latitude_scale * preservation_scale
    return deposition
end
function get_coastal_CaCO3_deposition_field( ocean_CO3 )
    # confined to coastal_depocenter cells, scaled by latitude
    depo_mask = eq_mask( world.geomorphology, coastal_depocenter )
    depo_rate = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            depo_rate[ix,iy] = depo_mask[ix,iy] * 
                get_coastal_CaCO3_deposition_rate( ocean_CO3, iy )
        end
    end
    return depo_rate
end
function get_coastal_CaCO3_deposition_rate( ocean_CO3, iy )
    latitude_scale = CaCO3_latitude_scale( iy )
    depo_rate_max = ocean_CO3 * exp( ( ocean_CO3 - 100. ) / 100. )
    # gives 100 m/Myr at ocean_CO3 = 100., about 1000. at 400. 
    return depo_rate_max * latitude_scale
end
function CaCO3_latitude_scale( iy )
    lat = ycoords[iy]
    lat_scale = exp( - ( ( lat / CaCO3_deposition_lat_scale )^2 ) )
    return lat_scale
end
function subaereal_CaCO3_dissolution()
    land_CaCO3_dissolution_rates = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if world.freeboard[ix,iy] > land_CaCO3_dissolving_altitude
                diss_rate = world.sediment_surface_fractions[ix,iy,CaCO3_sediment] *
                    land_CaCO3_dissolution_rate # meters/Myr
                available_caco3 =  world.sediment_thickness[ix,iy] * 
                    world.sediment_surface_fractions[ix,iy,CaCO3_sediment]
                diss_rate = min( diss_rate, available_caco3 )
                land_CaCO3_dissolution_rates[ix,iy] = diss_rate
            end
        end
    end
    return land_CaCO3_dissolution_rates
end
function subaereal_CaCO3_dissolution_original()
    reset_diag("land_CaCO3_dissolution_rate")
    for ix in 1:nx
        for iy in 1:ny
            if world.freeboard[ix,iy] > land_CaCO3_dissolving_altitude
                diss_rate = world.sediment_surface_fractions[ix,iy,CaCO3_sediment] *
                    land_CaCO3_dissolution_rate # meters/Myr
                available_caco3 =  world.sediment_thickness[ix,iy] * 
                    world.sediment_surface_fractions[ix,iy,CaCO3_sediment]
                diss_rate = min( diss_rate, available_caco3 )
                set_diag("land_CaCO3_dissolution_rate",ix,iy,diss_rate)
            end
        end
    end
end
