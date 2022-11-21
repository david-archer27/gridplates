function find_steady_state_ocean_CO3( target_CaCO3_deposition_rate )
    local ocean_CO3 
    #=terrestrial_CaCO3_dissolution = volume_field(
        get_diag("land_CaCO3_dissolution_rate") ) # m3 / Myr
    silicate_weathering_CaCO3_rate =
        global_CO2_degassing_rate * # mol / yr
        100. * 1.e6 / rho_sediment / 1.e6 # m3 / Myr
    
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
    =#
    # target_CaCO3_deposition_rate = 1.e15 # today fluxes
    ocean_CO3_low = 80.
    ocean_CO3_high = 100.
    imbalance = 1.e12; n_iters = 0; total_deposition = 0.
    depo_low = get_total_CaCO3_deposition_rate( ocean_CO3_low ) 
    depo_high = get_total_CaCO3_deposition_rate( ocean_CO3_high ) 
    while n_iters < 100 && imbalance > 1.e6
    #for iter in 1:10
        n_iters += 1
        ocean_CO3 = ocean_CO3_low + 
            ( target_CaCO3_deposition_rate - depo_low ) *
            ( ocean_CO3_high - ocean_CO3_low ) / 
            ( depo_high - depo_low )
        total_deposition = get_total_CaCO3_deposition_rate( ocean_CO3 )
        if total_deposition > target_CaCO3_deposition_rate # too much, 
            ocean_CO3_high = ocean_CO3; depo_high = total_deposition
        else
            ocean_CO3_low = ocean_CO3; depo_low = total_deposition
        end
        imbalance = abs( total_deposition - target_CaCO3_deposition_rate )
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
    if enable_continental_CaCO3_deposition
        for ix in 1:nx
            for iy in 1:ny
                if world.crust_type[ix,iy] == continent_crust
                    if world.freeboard[ix,iy] < shelf_depth_CaCO3
                        accom_meters = ( shelf_depth_CaCO3 - world.freeboard[ix,iy] ) /
                            sediment_freeboard_expression
                        shallow_rates[ix,iy] = accom_meters / main_time_step * CaCO3_latitude_scale( iy )
                    end
                end
            end
        end
    end
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
function generate_runoff_transect()
    runoff_tropics_width = 15.
    runoff_tropical_bulge_base = 1. # m / yr
    runoff_background_base = 0.1
    runoff_transect = fill(0.,ny)
    for iy in 1:ny
        tropical_bulge = runoff_tropical_bulge_base * 
            exp( - ycoords[iy]^2 / ( runoff_tropics_width )^2)
        runoff_transect[iy] = tropical_bulge + runoff_background_base
    end     
    q_transect = runoff_transect ./ 3.14e7 * # m3 / s / m2
        1000. * # l / s / m2
        1.e6 # l / km2 s
    return q_transect
end

function subaereal_sediment_dissolution()
    land_sediment_dissolution_rates = fill(0.,nx,ny,0:n_sediment_types)
    land_dissolved_Ca_prod_rates = fill(0.,nx,ny)
    q_transect = generate_runoff_transect()
    bedrock_FCO2_coeff = 0.25 # Suchet 2003 combined shield,basalt,acid volc
    carbonate_FCO2_coeff = 1.6 
    shale_FCO2_coeff = 0.6
    for ix in 1:nx
        for iy in 1:ny
            if world.freeboard[ix,iy] > 0.
                #=if world.geomorphology[ix,iy] == exposed_basement
                    land_dissolved_Ca_prod_rates[ix,iy] = 
                        bedrock_FCO2_coeff *
                        q_transect[iy] * 1.E-3 * # mol / km2 s 
                        areabox[iy] / 1.e6 * # mol / s
                        3.14E7 * 1.e6 # mol / Myr
                end=#
                if world.geomorphology[ix,iy] == sedimented_land
                    land_sediment_dissolution_rates[ix,iy,CaO_sediment] = 
                        world.sediment_surface_fractions[ix,iy,CaO_sediment] *
                        shale_FCO2_coeff * 
                        q_transect[iy] * 1.E-3 * # mol / km2 s
                        56. * 3.14e7 * 1.e6 / # g / km2 Myr
                        1.e6 / # g / m2 Myr
                        rho_sediment / 1.e6 # m3 / m2 Myr
                    land_sediment_dissolution_rates[ix,iy,CaCO3_sediment] = 
                        world.sediment_surface_fractions[ix,iy,CaCO3_sediment] *
                        carbonate_FCO2_coeff * 
                        q_transect[iy] * 1.E-3 * # mol / km2 s
                        100. * 3.14e7 * 1.e6 / # g / km2 Myr
                        1.e6 / # g / m2 Myr
                        rho_sediment / 1.e6 # m3 / m2 Myr
                        #world.sediment_thickness[ix,iy] * 
                        #world.sediment_surface_fractions[ix,iy,i_sedtype] *
                        #land_sediment_dissolution_rate_constants[i_sedtype] 
                        # * 
                        #( 1. + land_sediment_dissolution_2xCO2[i_sedtype] ) * 
                        #( world.atmCO2 / atmCO2_base )
                    for i_sedtype in 1:n_sediment_types
                        land_sediment_dissolution_rates[ix,iy,i_sedtype] = 
                            min(land_sediment_dissolution_rates[ix,iy,i_sedtype],
                            world.sediment_thickness[ix,iy] *
                                world.sediment_surface_fractions[ix,iy,i_sedtype] /
                                main_time_step * 0.3)
                            land_sediment_dissolution_rates[ix,iy,i_sedtype] =
                                max(land_sediment_dissolution_rates[ix,iy,i_sedtype],0.)
                        land_sediment_dissolution_rates[ix,iy,0] += 
                            land_sediment_dissolution_rates[ix,iy,i_sedtype]

                    end
                end
            end
        end
    end
    return land_sediment_dissolution_rates, land_dissolved_Ca_prod_rates
end
