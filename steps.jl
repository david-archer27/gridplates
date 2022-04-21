
# Meta steps
function step_everything() # rebuilds the world at the new time
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
    # sets world.elevation_offset for aging ocean crust
    step_geomorph() # plugs in here or operates standalone on world grid only
    # updates world.sediment_thickness, returns in isostatic equilibrium
    apply_crust_changes_to_plates()
    apply_sediment_fluxes_to_plates()
    # updates the plate crust_thickness and sediment_thickness fields
    return
end
function step_geomorph() # requires previous run with tectonics to fill world diag arrays
    local new_elevation_field, subaereal_mask, orogenic_source, 
    crust_clay_source, clay_src_land,land_tot_deposition,land_clay_deposition,
    land_CaCO3_deposition, land_clay_runoff,land_CaCO3_runoff
    verbose = true
    clear_geomorph_process_arrays() # only necessary in standalone mode
    # world. files only
    orogeny()
    # continent and subduction uplift rates into *_orogenic_uplift_rate diags
    apply_orogeny_fluxes_to_world()
    isostacy() 
    update_surface_types() 
    original_elevation_field = land_transport_elevation()
    #new_elevation_field = deepcopy( original_elevation_field ) # create a memory space outside loop
    #subaereal_mask = generate_mask_field(world.surface_type, sedimented_land) # ibid
    #orogenic_source = get_diag("land_orogenic_clay_flux") # metoo
    ocean_sink = generate_mask_field( world.crust_type, ocean_crust )
    # land bulk sediment transport block
    n_denuded = -1; n_first = true
    while n_denuded != 0 
        #elevation_field = land_transport_elevation() # start over each time
        subaereal_mask = generate_mask_field(world.surface_type, sedimented_land)
        generate_orogenic_erosion_fluxes()
        # sets land_orogenic_clay_flux, coastal_orogenic_clay_flux,
        # and crust_erosion_rate.
        orogenic_source = get_diag("land_orogenic_clay_flux") # m / Myr
        new_elevation_field = land_sediment_transport( original_elevation_field, 
            subaereal_mask, orogenic_source, ocean_sink )
        # sets land_sediment_deposition_rate, land_sediment_fraction_deposition_rate, 
        # land_orogenic_clay_flux, coastal_orogenic_clay_flux, crust_erosion_rate.
        aolean_transport()
        # resets and fills aolean_deposition_rate and aolean_deposition_rate
        n_denuded = check_subaereal_exposure()
        # sets world.surface_type to exclude for next pass
        if n_denuded > 0
            if n_first == true
                print("land denuding ",n_denuded)
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
    land_sediment_fraction_transport( new_elevation_field, 
        subaereal_mask, orogenic_source, ocean_sink ) 
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
            " d+r-s ",land_clay_deposition + land_clay_runoff - clay_src_land )
        # land CaCO3 mass balance
        println("land CaCO3 balance, src ", 0.,
           " dep ", land_CaCO3_deposition," runoff ",land_CaCO3_runoff,
           " d+r ",land_CaCO3_deposition + land_CaCO3_runoff)
    end

    distribute_ocean_sediment_fluxes(  )
    # accumulates seafloor_sediment_fraction_deposition_rate, seafloor_sediment_deposition_rate

    if verbose == true
        # ocean clay balance
        clay_influx = volumefield_total(
            get_frac_diag("ocean_sediment_fraction_influx",clay_sediment))
        ocn_clay_dep = volumefield_total(get_frac_diag(
            "seafloor_sediment_fraction_deposition_rate",clay_sediment ))
        println("clay ocean balance, in ", clay_influx," dep ",ocn_clay_dep,
            " bal ", clay_influx - ocn_clay_dep )
        # CaCO3 balance
        CaCO3_influx = volumefield_total(
            get_frac_diag("ocean_sediment_fraction_influx",CaCO3_sediment))
        ocn_CaCO3_dep = volumefield_total(get_frac_diag(
            "seafloor_sediment_fraction_deposition_rate",CaCO3_sediment ))
        println("CaCO3 balance, in ", CaCO3_influx," dep ",ocn_CaCO3_dep,
            " bal ", CaCO3_influx - ocn_CaCO3_dep )
    end

    world.sediment_thickness, world.sediment_fractions = apply_sediment_fluxes_to_world()
    isostacy()

    return
end
