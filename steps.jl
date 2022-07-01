function run_timeseries()
    cd( data_output_directory * "/world" )
    file_list = readdir()
    for file in file_list
        println("deleting ", file)
        run(`rm $file`)
    end
    cd( data_output_directory * "/plates" )
    file_list = readdir()
    for file in file_list
        println("deleting ", file)
        run(`rm $file`)
    end
    cd( code_base_directory )
    log_IO = open("../outfiles/logfile.txt","w")
    while world.age > 0
        step_everything()
        save_world()
        if floor(world.age/5) == world.age/5 && world.age >= 10
            save_plates()
        end
    end
    close(log_IO)
end
function tabulate_diagnostics( age )
    world = read_world( age )
    n_time_points = world.age / time_step
    n_diags = length(world_diag_names)
    n_frac_diags = length(world_frac_diag_names)
    diags_timeseries = fill(0.,n_diags,n_time_points)
    frac_diags_timeseries = fill(0.,n_frac_diags,n_sediment_types,n_time_points)
    inventory_timeseries = fill(0.,n_sediment_types,n_time_points)
    time_points = fill(0.,n_time_points)
    i_time_point = 0
    while world.age > 0
        age -= time_step
        println(age)
        push!( time_points, age )
        i_time_point += 1
        world = read_world( age )
        inventory_timeseries[:,i_time_point] = 
            world_sediment_inventories( )[1:n_sediment_types]
        for i_diag in 1:n_diags
            diag_volume = volume_field( get_diag(world_diag_names[i_diag] ))
            diags_timeseries[i_diag,i_time_point] = diag_volume
        end
        for i_frac_diag in 1:n_frac_diags
            for i_sedtype in 1:n_sediment_types
                frac_diag_volume = volume_field( get_frac_diag(world_diag_names[i_diag], i_sedtype ))
                frac_diags_timeseries[i_diag,i_sedtype,i_time_point] = frac_diag_volume
            end
        end
    end
    return inventory_timeseries, diags_timeseries, frac_diags_timeseries, time_points
end

# Meta steps
function step_everything() # rebuilds the world at the new time
    verbose = true
    local initial_plate_inventories, update_ID_plate_inventories, 
        remasked_plate_inventories, remasked_world_inventories, 
        subduction_rates, subduction_balances,
        initial_world_inventories, filled_world_inventories,
        updated_world_inventories, land_trapped # , world_grid_sediment_change

    world.age -= time_step
    println()
    println("beginning ", world.age, " Myr ", get_geo_interval())
    if log_IO != 0
        println(log_IO,"beginning ", world.age, " Myr ", get_geo_interval())
    end

    if verbose == true
        initial_plate_inventories = global_plate_sediment_inventories()
        initial_world_inventories = world_sediment_inventories(  )
        world_grid_sediment_change = world.sediment_thickness
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
        
        #scale_global_sediment_components( imbalance_ratios[1:n_sediment_types] )

        fudged_world_inventories = world_sediment_inventories(  )
        #logging_println("  fudge factors ", imbalance_ratios)
        #logging_println("         fudged ", fudged_world_inventories)
        #logging_println("  cum world bal ", ( fudged_world_inventories .- initial_world_inventories ) ./
        #    time_step .+ subduction_rates)
    end

    world.elevation_offset = ocean_thermal_boundary_layer()
    # sets world.elevation_offset for aging ocean crust
    step_geomorph() # plugs in here or operates standalone on world grid only
    # updates world.sediment_thickness and layers, returns in isostatic equilibrium

    apply_geomorphology_changes_to_plates()
    # updates the plate crust_thickness and sediment_thickness fields

    if verbose == true
        world_grid_sediment_change = world_grid_sediment_change .* -1. .+ world.sediment_thickness
        set_diag("world_grid_sediment_change", world_grid_sediment_change)
        geomorphed_plate_inventories = global_plate_sediment_inventories()
        plate_change_rates = ( geomorphed_plate_inventories .- initial_plate_inventories ) /
            time_step
        #deposition_rates = [ volume_field(get_diag("global_sediment_deposition_rate")) ]
        #clay_flux = volume_field(get_diag("crust_clay_source_rate"))
        #CaCO3_flux = volume_field(get_frac_diag("global_sediment_fraction_deposition_rate")[:,:,CaCO3_sediment])
        source_fluxes = fill(0.,0:n_sediment_types)
        deposition_rates = fill(0.,0:n_sediment_types)
        land_trapped = fill(0.,0:n_sediment_types)
        source_fluxes[clay_sediment] = volume_field(get_diag("crust_clay_source_rate"))
        source_fluxes[CaCO3_sediment] = volume_field(get_frac_diag("global_sediment_fraction_deposition_rate")[:,:,CaCO3_sediment])
        
        for i_sedtype in 1:n_sediment_types
            #source_fluxes[i_sedtype] += volume_field(get_frac_diag("land_trapped_sediment_rate",i_sedtype))
            deposition_rates[i_sedtype] = volume_field(
                get_frac_diag("global_sediment_fraction_deposition_rate",i_sedtype))
            land_trapped[i_sedtype] = volume_field(
                get_frac_diag("land_trapped_sediment_rate",i_sedtype ))
            source_fluxes[0] += source_fluxes[i_sedtype]
            deposition_rates[0] += deposition_rates[i_sedtype]
            land_trapped[0] += land_trapped[i_sedtype]
        end 
        geomorph_world_inventories = world_sediment_inventories( )
        #world_change_rates = ( geomorph_world_inventories .- initial_world_inventories) ./ time_step
        #world_balances = world_change_rates .+ subduction_rates .- deposition_rates
        plate_balances = plate_change_rates .+ subduction_rates .- deposition_rates
        logging_println()
        logging_println(" world budget init ", initial_world_inventories)
        logging_println("            filled ", filled_world_inventories)
        logging_println("               bal ", ( filled_world_inventories .- initial_world_inventories ) ./
            time_step )
        logging_println("           updated ", updated_world_inventories)
        logging_println("           cum bal ", ( updated_world_inventories .- initial_world_inventories ) ./
            time_step .+ subduction_rates )
        logging_println("    fudged cum bal ", ( fudged_world_inventories .- initial_world_inventories ) ./
            time_step .+ subduction_rates )
        logging_println("     source inputs ", source_fluxes)
        change = ( geomorph_world_inventories .- fudged_world_inventories ) ./ time_step 
        logging_println("    final geomorph ", geomorph_world_inventories)
        logging_println("       change rate ", change )
        logging_println("      geomorph bal ", change .- source_fluxes .+ land_trapped )
        change = ( geomorph_world_inventories .- initial_world_inventories ) ./ time_step 
        logging_println("        deposition ", deposition_rates)
        logging_println("        subduction ", subduction_rates)
        logging_println("  world cum change ", change )
        logging_println("     world cum bal ", change .+ subduction_rates .- source_fluxes)
        
        logging_println()
        logging_println("plate budgets orig ", initial_plate_inventories)
        logging_println("          updateID ", update_ID_plate_inventories)
        logging_println("          remasked ", remasked_plate_inventories)
        #println("               bal ", subduction_balances)
        logging_println("          geomorph ", geomorphed_plate_inventories)
        logging_println("        subduction ", subduction_rates)
        logging_println("        deposition ", deposition_rates)
        change = ( geomorphed_plate_inventories .- initial_plate_inventories ) ./ time_step
        logging_println("            change ", change)
        logging_println("         plate bal ", plate_balances)

        is_land = gt_mask(world.freeboard,0.)
        mean_land_freeboard = field_mean(world.freeboard .* is_land)
        mean_ocean_depth = field_mean(world.freeboard .* (1. .- is_land))
        println();println("mean elevations ",[mean_land_freeboard,mean_ocean_depth])
    end

end

function step_geomorph() # requires previous run with tectonics to fill world diag arrays

    local new_elevation_field, subaereal_mask, 
        #crust_clay_source, land_tot_deposition,
        clay_src_land,land_clay_deposition,
        land_CaCO3_deposition, land_clay_runoff,land_CaCO3_runoff,
        old_land_clay_inventory, old_ocean_sed_inventories, land_trapped

    verbose = true

    if verbose == true
        old_land_fraction_inventories = world_land_sediment_inventories( )
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
    println()#;print("land")

    # sets world.tectonics for reburial
    original_elevation_field = land_transport_elevation()
    ocean_sink = eq_mask( world.crust_type, ocean_crust )
    subaereal_mask = eq_mask(world.geomorphology, sedimented_land)
    new_elevation_field = original_elevation_field
    # land bulk sediment transport
    n_denuded = -1; n_first = true
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
        #if n_denuded > 0
        depo = volume_field(get_diag("land_sediment_deposition_rate"))
        land_trapped = volume_field(get_frac_diag("land_trapped_sediment_rate")[:,:,1])
        aolean = volume_field(get_diag("aolean_clay_erosion_rate"))
        println("n, depo, aolean, trapped ",[n_denuded,depo,aolean,
            volume_field(get_frac_diag("land_trapped_sediment_rate")[:,:,1])])
            #=
            if n_first == true
                print(" denuding ",n_denuded)
                n_first = false
            else
                print(", ",n_denuded)
            end
            =#
            #=if verbose == true # watch the erosion flux grow
                crust_clay_source = volume_field(get_diag("crust_clay_source_rate"))
                clay_src_land = volume_field(get_diag("land_orogenic_clay_flux"))
                clay_src_ocn = volume_field(get_diag("coastal_orogenic_clay_flux"))
                #land_tot_deposition = volume_field(get_diag("land_sediment_deposition_rate"))
                print(" crust clay source ",crust_clay_source," land ",clay_src_land," ocn ", clay_src_ocn)
            else
                println("")
            end=#
        #end
    end         
    println("")

    area_denuded = volume_field(eq_mask(world.geomorphology,exposed_basement))
    area_covered = volume_field(eq_mask(world.geomorphology,sedimented_land) )
    println("fraction denuded ", area_denuded / ( area_denuded + area_covered ))

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
    println(" original runoff ",volume_field(get_frac_diag("coastal_sediment_fraction_runoff_flux",1)))


    redistribute_denuded_land_sediment( )
    println(" augmented runoff ",volume_field(get_frac_diag("coastal_sediment_fraction_runoff_flux",1)))
    # captures stuff in land_trapped_sediment_rate, adds it to 
    # coastal_sediment_fraction_runoff_flux

    if verbose == true
        old_land_fraction_inventories = world_land_sediment_inventories( )
        if old_land_fraction_inventories[1] != old_land_fraction_inventories[1]
            error("NaN detected")
        end
        println("before update ", old_land_fraction_inventories)
    end

    world.sediment_thickness, world.sediment_surface_fractions = 
        apply_land_sediment_fluxes()

    if verbose == true
        clay_src_land = volume_field(get_diag("land_orogenic_clay_flux"))
        sources = fill(0.,0:n_sediment_types); deposition = fill(0.,0:n_sediment_types)
        runoff = fill(0.,0:n_sediment_types); aolean = fill(0.,0:n_sediment_types)
        land_trapped = fill(0.,0:n_sediment_types)
        sources[1] = clay_src_land
        aolean[1] = volume_field(get_diag("aolean_clay_erosion_rate"))
        for i_sedtype in 1:n_sediment_types
            deposition[i_sedtype] = volume_field(
                get_frac_diag("land_sediment_fraction_deposition_rate",i_sedtype ))
            runoff[i_sedtype] = volume_field(
                get_frac_diag("coastal_sediment_fraction_runoff_flux",i_sedtype ))
            land_trapped[i_sedtype] = volume_field(
                get_frac_diag("land_trapped_sediment_rate",i_sedtype ))
            sources[0] += sources[i_sedtype]
            deposition[0] += deposition[i_sedtype]
            runoff[0] += runoff[i_sedtype]
            aolean[0] += aolean[i_sedtype]
            land_trapped[0] += land_trapped[i_sedtype]
        end
        new_land_fraction_inventories = world_land_sediment_inventories( ) 
        change_rate = ( new_land_fraction_inventories .- old_land_fraction_inventories ) / 
            time_step
        flux_balances = sources .+ land_trapped .- deposition .- runoff .- aolean 
        predicted_change_rate = deposition # includes aolean already
        inv_balances = change_rate .- predicted_change_rate
        logging_println()
        logging_println("land budget init ",old_land_fraction_inventories)
        logging_println("           final ",new_land_fraction_inventories)
        logging_println("     change rate ",change_rate)
        logging_println("             dep ",deposition)
        logging_println("   orogen source ",sources)
        logging_println("          runoff ",runoff)
        logging_println("         trapped ",land_trapped)
        logging_println("          aolean ",aolean)
        logging_println("   flux to ocean ",runoff .+ aolean)
        logging_println("        mass bal ",inv_balances)
        logging_println("        flux bal ",flux_balances)
    end

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
        old_ocean_sed_inventories = world_ocean_sediment_inventories( )
    end

    apply_ocean_sediment_fluxes( )
    # updates world.sediment* for ocean points from seafloor_sediment*_deposition_rate

    if verbose == true
        new_ocean_sed_inventories = world_ocean_sediment_inventories( )
        total_sources = fill(0.,0:n_sediment_types); coastal_sources = fill(0.,0:n_sediment_types);
        land_sources = fill(0.,0:n_sediment_types)
        pelagic_sources = fill(0.,0:n_sediment_types); orogenic_sources = fill(0.,0:n_sediment_types);
        ocean_sed_depo_rates = fill(0.,0:n_sediment_types)
        CaCO3_coastal_source = volume_field(get_diag("coastal_CaCO3_flux")) 
        CaCO3_pelagic_source =  volume_field(get_diag("pelagic_CaCO3_deposition_rate"))
        CaCO3_source = CaCO3_coastal_source + CaCO3_pelagic_source
        clay_runoff_source = volume_field( get_frac_diag("coastal_sediment_fraction_runoff_flux",clay_sediment) )
        CaCO3_runoff_source = volume_field( get_frac_diag("coastal_sediment_fraction_runoff_flux",CaCO3_sediment) )
        clay_orogenic_source = volume_field( get_diag("coastal_orogenic_clay_flux") )
        clay_aolean_source = volume_field( get_diag("aolean_clay_deposition_rate") )
        clay_source = clay_runoff_source + clay_orogenic_source + clay_aolean_source
        total_sources[clay_sediment] = clay_source
        total_sources[CaCO3_sediment] = CaCO3_coastal_source + CaCO3_pelagic_source + CaCO3_runoff_source
        land_sources[clay_sediment] = clay_runoff_source + clay_aolean_source
        coastal_sources[clay_sediment] = clay_runoff_source 
        coastal_sources[CaCO3_sediment] = CaCO3_coastal_source
        pelagic_sources[clay_sediment] = clay_aolean_source; pelagic_sources[CaCO3_sediment] = CaCO3_pelagic_source
        orogenic_sources[clay_sediment] = clay_orogenic_source
        for i_sedtype in 1:n_sediment_types
            ocean_sed_depo_rates[i_sedtype] = volume_field(
                get_frac_diag("seafloor_sediment_fraction_deposition_rate",i_sedtype))
            total_sources[0] += total_sources[i_sedtype]
            land_sources[0] += land_sources[i_sedtype]
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
        logging_println("     land sources ", land_sources)
        logging_println("  coastal sources ", coastal_sources)
        logging_println("          pelagic ", pelagic_sources)
        logging_println("         orogenic ", orogenic_sources)
        logging_println("              tot ", total_sources)
        logging_println("              dep ", ocean_sed_depo_rates)
        logging_println("         mass bal ", inventory_change_rate -
            ocean_sed_depo_rates )
        logging_println("         flux bal ", total_sources .-
            ocean_sed_depo_rates .- land_trapped )
    end                    

    isostacy()

    return
end
