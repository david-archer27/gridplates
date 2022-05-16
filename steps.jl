
# Meta steps
function step_everything() # rebuilds the world at the new time
    verbose = true
    local old_plate_inventories
    world.age -= time_step
    println()
    println("beginning ", world.age, " Myr ", get_geo_interval())

    if verbose == true
        println("going in ",  global_plate_sediment_inventories())
    end

    clear_world_process_arrays() 
    println("interpolating from plate grids")
    resolve_rotation_matrices()
    # set rotation matrices to the new time
    increment_plate_age()
    # ocean and cont crust gets older in plate grids
    update_changing_plateIDs()
    if verbose == true
        old_plate_inventories = global_plate_sediment_inventories()
        println("updated ", old_plate_inventories)
    end
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
    if verbose == true
        subducted_fractions = []
        for i_sedtype in 1:n_sediment_types
            push!(subducted_fractions, volumefield_total(
                get_frac_diag("ocean_subduct_sediment_fraction_volume")[:,:,i_sedtype]))
        end
        new_plate_inventories = global_plate_sediment_inventories()
        balances = new_plate_inventories .- old_plate_inventories .+ 
            ( subducted_fractions .* time_step )
        println("remask ",  global_plate_sediment_inventories()," subd ",subducted_fractions,
            " bal ", balances )
    end

    ocean_thermal_boundary_layer()
    # sets world.elevation_offset for aging ocean crust
    step_geomorph() # plugs in here or operates standalone on world grid only
    # updates world.sediment_thickness and layers, returns in isostatic equilibrium

    if verbose == true
        old_plate_inventories = global_plate_sediment_inventories()
        println("stepped ", old_plate_inventories)
    end

    apply_geomorphology_changes_to_plates()
    # updates the plate crust_thickness and sediment_thickness fields

    if verbose == true
        new_plate_inventories = global_plate_sediment_inventories()
        change_rates = ( new_plate_inventories .- old_plate_inventories ) ./
            time_step
        deposition_rates = []
        for i_sedtype in 1:n_sediment_types
            push!(deposition_rates, volumefield_total(get_frac_diag("global_sediment_fraction_deposition_rate",
                i_sedtype)))
        end
        balances = change_rates .- deposition_rates
        println("plate budgets ", old_plate_inventories," to ",new_plate_inventories,
            " bal ", balances)
    end

end

function step_geomorph() # requires previous run with tectonics to fill world diag arrays

    local new_elevation_field, subaereal_mask, 
        #crust_clay_source, land_tot_deposition,
        clay_src_land,land_clay_deposition,
        land_CaCO3_deposition, land_clay_runoff,land_CaCO3_runoff,
        old_land_clay_inventory, old_ocean_sed_inventories

    verbose = true
    clear_geomorph_process_arrays() # only necessary in standalone mode
    # world. files only
    orogeny()
    # continent and subduction uplift rates into *_orogenic_uplift_rate diags
    apply_orogeny_fluxes_to_world()
    isostacy() 
    print("land")
    check_reburial_exposed_basement() 
    original_elevation_field = land_transport_elevation()
    ocean_sink = generate_mask_field( world.crust_type, ocean_crust )

    # land bulk sediment transport
    n_denuded = -1; n_first = true
    while n_denuded != 0 
        #elevation_field = land_transport_elevation() # start over each time
        subaereal_mask = generate_mask_field(world.surface_type, sedimented_land)
        subaereal_mask[:,1] .= 0; subaereal_mask[:,ny] .= 0
        generate_orogenic_erosion_fluxes()
        # sets land_orogenic_clay_flux, coastal_orogenic_clay_flux,
        # and crust_erosion_rate.
        aolean_transport()
        # resets and fills aolean_deposition_rate and aolean_deposition_rate
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
        # sets world.surface_type to exclude for next pass
        if n_denuded > 0
            if n_first == true
                print(" denuding ",n_denuded)
                n_first = false
            else
                print(", ",n_denuded)
            end
            #=if verbose == true # watch the erosion flux grow
                crust_clay_source = volumefield_total(get_diag("crust_clay_source_rate"))
                clay_src_land = volumefield_total(get_diag("land_orogenic_clay_flux"))
                clay_src_ocn = volumefield_total(get_diag("coastal_orogenic_clay_flux"))
                #land_tot_deposition = volumefield_total(get_diag("land_sediment_deposition_rate"))
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

    #=if verbose == true
        land_clay_deposition = volumefield_total(
            get_frac_diag("land_sediment_fraction_deposition_rate",clay_sediment ))
        land_CaCO3_deposition = volumefield_total(
            get_frac_diag("land_sediment_fraction_deposition_rate",CaCO3_sediment ))
        # check the fractions add up
        println("land clay dep ",land_clay_deposition,
            " CaCO3 ",land_CaCO3_deposition,
            " c+C ", land_clay_deposition+land_CaCO3_deposition)
    end=#
    
    set_land_runoff_fluxes( new_elevation_field, subaereal_mask, ocean_sink ) 
    # fills coastal_sediment_fraction_runoff_flux

    if verbose == true
        clay_src_land = volumefield_total(get_diag("land_orogenic_clay_flux"))
        land_clay_deposition = volumefield_total(
            get_frac_diag("land_sediment_fraction_deposition_rate",clay_sediment ))
        land_CaCO3_deposition = volumefield_total(
            get_frac_diag("land_sediment_fraction_deposition_rate",CaCO3_sediment ))
        land_clay_runoff = volumefield_total(
            get_frac_diag("coastal_sediment_fraction_runoff_flux",clay_sediment ))
        land_CaCO3_runoff = volumefield_total(
            get_frac_diag("coastal_sediment_fraction_runoff_flux",CaCO3_sediment ))
        # land clay mass balance
        println("land clay balance, src ", clay_src_land,
            " dep ", land_clay_deposition," runoff ",land_clay_runoff,
            " bal ", clay_src_land - land_clay_deposition - land_clay_runoff )
        # land CaCO3 mass balance
        println("land CaCO3 balance, src ", 0.,
           " dep ", land_CaCO3_deposition," runoff ",land_CaCO3_runoff,
           " bal ",land_CaCO3_deposition + land_CaCO3_runoff)
    end

    distribute_CaCO3_sedimentation( )
    # sets coastal_CaCO3_flux, pelagic_CaCO3_flux
    distribute_ocean_sediment_fluxes(  )
    # accumulates seafloor_sediment_fraction_deposition_rate, seafloor_sediment_deposition_rate

    if verbose == true
        # ocean clay balance
        ocn_clay_influx = volumefield_total(
            get_frac_diag("ocean_sediment_fraction_influx",clay_sediment))
        ocn_clay_dep = volumefield_total(get_frac_diag(
            "seafloor_sediment_fraction_deposition_rate",clay_sediment ))
        println("clay ocean balance, in ", ocn_clay_influx," dep ",ocn_clay_dep,
            " bal ", ocn_clay_influx - ocn_clay_dep )
        # CaCO3 balance
        CaCO3_influx = volumefield_total(
            get_frac_diag("ocean_sediment_fraction_influx",CaCO3_sediment))
        ocn_CaCO3_dep = volumefield_total(get_frac_diag(
            "seafloor_sediment_fraction_deposition_rate",CaCO3_sediment ))
        println("CaCO3 balance, in ", CaCO3_influx," dep ",ocn_CaCO3_dep,
            " bal ", CaCO3_influx - ocn_CaCO3_dep )
    end

    if verbose == true
        old_land_clay_inventory = volumefield_total(
            world.sediment_thickness .* 
            world.sediment_surface_fractions[:,:,clay_sediment] )
    end

    world.sediment_thickness, world.sediment_surface_fractions = 
        apply_land_sediment_fluxes()
    # leaves the actual substitution to passback so it could be applied provisionally
    # in iteration

    if verbose == true # testing land sedimentation
        new_land_clay_inventory = volumefield_total(
            world.sediment_thickness .* 
            world.sediment_surface_fractions[:,:,clay_sediment] )
        change_rate = ( new_land_clay_inventory - old_land_clay_inventory ) /
            time_step
        balance = change_rate - land_clay_deposition
        println("land clay accum ", old_land_clay_inventory, " to ", new_land_clay_inventory,
           " chg ", change_rate," bal ", balance)
    end

    if verbose == true
        old_ocean_sed_inventories = []
        for i_sedtype in 1:n_sediment_types
            push!(old_ocean_sed_inventories, volumefield_total(
                world.sediment_layer_thickness[:,:,current_time_bin()] .* 
                world.sediment_layer_fractions[:,:,i_sedtype,current_time_bin()] ) )
        end 
    end

    apply_ocean_sediment_fluxes( )
    # updates world.sediment* for ocean points from seafloor_sediment*_deposition_rate
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

    if verbose == true
        new_ocean_sed_inventories = []
        inventory_changes = 
        ocean_sed_depo_rates = []
        balances = []
        for i_sedtype in 1:n_sediment_types
            push!(new_ocean_sed_inventories, volumefield_total(
                world.sediment_layer_thickness[:,:,current_time_bin()] .* 
                world.sediment_layer_fractions[:,:,i_sedtype,current_time_bin()] ) )
            push!(ocean_sed_depo_rates,volumefield_total(
                get_frac_diag("seafloor_sediment_fraction_deposition_rate",i_sedtype)))
        end 
        inventory_changes = ( new_ocean_sed_inventories .- old_ocean_sed_inventories) ./
            time_step
        balances = inventory_changes -
            ocean_sed_depo_rates 
        println("ocean accum ", old_ocean_sed_inventories," to ",
            new_ocean_sed_inventories, 
            " bal ", balances)
    end                    

    isostacy()

    return
end
