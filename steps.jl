function step_everything(  ) # rebuilds the world at the new time
    local world_grid_sediment_change,initial_plate_inventories,
        initial_world_inventories, after_tectonics_world_inventories

    if enable_step_everything_diagnostics
        initial_plate_inventories = global_plate_sediment_inventories()
        initial_world_inventories = world_sediment_inventories(  )
        world_grid_sediment_change = world.sediment_thickness
    end

    world.initial_land_sediment_inventories[:] = 
        world_land_sediment_inventories()[1:end]
    world.initial_ocean_sediment_inventories[:] = 
        world_ocean_sediment_inventories()[1:end]

    step_tectonics(  )

    if enable_step_everything_diagnostics
        after_tectonics_world_inventories = world_sediment_inventories(  )
    end

    world.elevation_offset = ocean_thermal_boundary_layer()
    world.sealevel = get_sealevel( world.age )
    logging_println("Sea level ", world.sealevel)

    if enable_geomorph
        step_geomorph(  ) # plugs in here or operates standalone on world grid only
        # updates world.sediment_thickness and layers, returns in isostatic equilibrium

        apply_geomorphology_changes_to_plates()
        # updates the plate crust_thickness and sediment_thickness fields
    end

    if enable_step_everything_diagnostics
        subduction_rates = fill(0.,0:n_sediment_types)
        subduction_rates[1:n_sediment_types] .= 
            world.subducted_ocean_sediment_volumes .+
            world.subducted_land_sediment_volumes
        for i_sedtype in 1:n_sediment_types
            subduction_rates[0] += subduction_rates[i_sedtype]
        end
        world_grid_sediment_change = world_grid_sediment_change .* -1. .+ world.sediment_thickness
        set_diag("world_grid_sediment_change", world_grid_sediment_change)
        geomorphed_world_inventories = world_sediment_inventories(  )
        final_plate_inventories = global_plate_sediment_inventories()
        source_fluxes = fill(0.,0:n_sediment_types)
        dissolution_fluxes = fill(0.,0:n_sediment_types)
        deposition_rates = fill(0.,0:n_sediment_types)
        denuded_land_src = fill(0.,0:n_sediment_types)
        cont2ocn_redist_sources = fill(0.,0:n_sediment_types)
        ocn2cont_redist_sources = fill(0.,0:n_sediment_types)
        source_fluxes[1:end] = 
            volume_fields(get_frac_diag("crust_orogenic_fraction_flux")) # no denuded sediment
        source_fluxes[CaCO3_sediment] += 
            volume_field(get_diag("coastal_CaCO3_flux")) +
            volume_field(get_diag("pelagic_CaCO3_deposition_rate")) +
            volume_field(get_diag("continental_CaCO3_deposition_rate"))
            #volume_field(get_frac_diag("global_sediment_fraction_deposition_rate")[:,:,CaCO3_sediment])
        dissolution_fluxes[1:end] = 
            volume_fields(get_frac_diag("land_sediment_fraction_dissolution_rate")) # no denuded sediment
        deposition_rates[1:end] = volume_fields(
            get_frac_diag("global_sediment_fraction_deposition_rate"))
        denuded_land_src[1:end] = volume_fields(
            get_frac_diag("denuded_sediment_source_fraction_flux"))
        cont2ocn_redist_sources[1:end] = 
            volume_fields( get_frac_diag(
                "continent_2_ocean_sediment_fraction_displaced") )
        ocn2cont_redist_sources[1:end] = 
            volume_fields( get_frac_diag(
                "ocean_2_continent_sediment_fraction_displaced") )
        for i_sedtype in 1:n_sediment_types
            source_fluxes[0] += source_fluxes[i_sedtype]
            dissolution_fluxes[0] += dissolution_fluxes[i_sedtype]
            deposition_rates[0] += deposition_rates[i_sedtype]
            denuded_land_src[0] += denuded_land_src[i_sedtype]
            cont2ocn_redist_sources[0] += cont2ocn_redist_sources[i_sedtype]
            ocn2cont_redist_sources[0] += ocn2cont_redist_sources[i_sedtype]
        end 

        #source_fluxes .-= denuded_land_src .- ocn2cont_redist_sources

        #world_change_rates = ( geomorph_world_inventories .- initial_world_inventories) ./ main_time_step
        #world_balances = world_change_rates .+ subduction_rates .- deposition_rates
        plate_change = ( final_plate_inventories .- initial_plate_inventories ) ./ main_time_step
        plate_balances = plate_change .+ subduction_rates .- deposition_rates
        geomorph_change = ( geomorphed_world_inventories .- after_tectonics_world_inventories ) ./ main_time_step 
        step_change = ( geomorphed_world_inventories .- initial_world_inventories ) ./ main_time_step 
 
        logging_println()
        logging_println(" world budget init ", initial_world_inventories)
        logging_println("   after tectonics ", after_tectonics_world_inventories)
        #logging_println("               bal ", ( after_tectonics_world_inventories .- initial_world_inventories ) ./
        #    main_time_step .+ subduction_rates ) #.+ crust_transition_ocean_src )
 
        logging_println("     source inputs ", source_fluxes)
        logging_println("       dissolution ", dissolution_fluxes)
        logging_println("cont2ocn trans flx ", cont2ocn_redist_sources)
        logging_println("ocn2cont trans flx ", ocn2cont_redist_sources)
        logging_println("        deposition ", deposition_rates)
        logging_println("    final geomorph ", geomorphed_world_inventories)
        logging_println("       change rate ", geomorph_change )
        logging_println(" geomorph mass bal ", geomorph_change .- 
            deposition_rates .+ dissolution_fluxes )
        logging_println(" geomorph flux bal ", source_fluxes .- deposition_rates .+ # genuine erosion
            cont2ocn_redist_sources .+ ocn2cont_redist_sources ) #.- 
            #dissolution_fluxes )
            
        logging_println("        subduction ", subduction_rates)
        logging_println("  world cum change ", step_change )
        logging_println("     world cum bal ", step_change .+ subduction_rates .- 
            source_fluxes .+ dissolution_fluxes ) # .+ ocn2cont_redist_sources ) # .+ cont2ocn_redist_sources )
        
        logging_println()
        logging_println("plate budgets orig ", initial_plate_inventories)
        logging_println("             final ", final_plate_inventories)
        logging_println("            change ", plate_change)
        logging_println("         plate bal ", plate_balances)

        subaereal_mask = gt_mask(world.freeboard,0.)
        maxix,maxiy = highest_xy_field(world.freeboard)
        max_freeboard = world.freeboard[maxix,maxiy]
        mean_land_freeboard = field_mean(world.freeboard .* subaereal_mask)
        mean_ocean_depth = field_mean(world.freeboard .* (1 .- subaereal_mask))
        logging_println()
        logging_println("max, mean elevations ",[field_max(world.freeboard), mean_land_freeboard,mean_ocean_depth])
        c_uplift = get_diag("continent_orogenic_uplift_rate") 
        s_uplift = get_diag("subduction_orogenic_uplift_rate")  # meters / Myr
        erod = get_diag("crust_erosion_rate")
        logging_println("fastest uplift cont / subd / sum ", 
            [ field_max(c_uplift), field_max(s_uplift), 
            field_max( c_uplift + s_uplift ) ] )
        logging_println("fastest erosion, net ", 
            [ field_max(erod),
            field_max( c_uplift + s_uplift - erod) ] )
        logging_println("")
    end

end
function step_tectonics(  )
    local initial_plate_inventories, update_ID_plate_inventories, 
        remasked_plate_inventories, #remasked_world_inventories, delete_old_plate_inventories,
        subduction_rates, #subduction_balances,
        initial_world_inventories, filled_world_inventories,
        updated_world_inventories, orogenic_uplift_rates
        #land_trapped # , world_grid_sediment_change

    #world.age -= main_time_step
    #world.sealevel = get_sea_level( world.age )

    #logging_println("beginning ", world.age - main_time_step, " Myr ", get_geo_interval())
    logging_println()
    logging_println("Plate tectonics ")

    subducted_land_sediment_volumes = fill(0.,n_sediment_types)
    subducted_ocean_sediment_volumes = fill(0.,n_sediment_types)
    orogenic_uplift_rates = fill(0.,nx,ny,0:2) # total, continental, subduction
    initial_world_inventories = fill(0.,0:2)
    update_ID_plate_inventories = fill(0.,0:2)
    n_plate_substeps = Int( main_time_step / sub_time_step )
    clear_world_process_arrays() 
    initial_plate_inventories = fill(0.,0:2)
    if enable_step_tectonics_plate_diagnostics
        initial_plate_inventories = global_plate_sediment_inventories()
    end
    initial_world_inventories = fill(0.,0:2)
    if enable_step_tectonics_world_diagnostics || eliminate_regrid_drift
        initial_world_inventories = world_sediment_inventories(  )
    end

    for i_plate_substep in 1:n_plate_substeps
        
        if enable_substep_tectonics_diagnostics 
            substep_initial_world_inventories = world_sediment_inventories(  )
        end

        world.age -= sub_time_step
        #logging_println(); 
        logging_println("  substepping ",world.age," Myr")
        if enable_substep_tectonics_diagnostics
            substep_initial_plate_inventories = global_plate_sediment_inventories()
            logging_println("  substep init plate ", substep_initial_plate_inventories)
        end
        oldIDmap = world.plateID
        newIDmap = read_plateIDs(world.age)
        world.plateID = newIDmap
        world.plateIDlist = find_plateID_list()

        if enable_eyeball_changing_plateIDs
            changingplates = []
            if haskey(plateID_change_log,world.age)
                changingplates = plateID_change_log[world.age]
                if length(changingplates) > 0
                    logging_println("    migrating plate info ",changingplates)
                end
            end
            update_changing_plateIDs(oldIDmap,newIDmap,changingplates)
        else
            update_changing_plateIDs(oldIDmap,newIDmap)
        end
        if enable_substep_tectonics_diagnostics
            substep_update_ID_plate_inventories = global_plate_sediment_inventories()
            logging_println("          update IDs ", substep_update_ID_plate_inventories)
        end

    end


    resolve_rotation_matrices()
    fill_world_from_plates() # build the new world grid

    world.crust_age .+= main_time_step # sub_time_step

    substep_subducted_land_sediment_volumes, substep_subducted_ocean_sediment_volumes,
        subduction_footprint = 
        remask_plates() 

    #error("stopping")

    substep_subducted_sediment_volumes = 
        substep_subducted_land_sediment_volumes .+ 
        substep_subducted_ocean_sediment_volumes
    subducted_land_sediment_volumes .+= substep_subducted_land_sediment_volumes
    subducted_ocean_sediment_volumes .+= substep_subducted_ocean_sediment_volumes
    subducted_sediment_volumes = 
        subducted_land_sediment_volumes .+ 
        subducted_ocean_sediment_volumes
    substep_remask_plate_inventories = fill(0.,0:n_sediment_types)
    if enable_substep_tectonics_diagnostics
        substep_remask_plate_inventories = global_plate_sediment_inventories()
        logging_println("              remask ", substep_remask_plate_inventories[1:end])
        logging_println("             subduct ", substep_subducted_sediment_volumes[1:end])
        logging_println("       plate inv bal ", ( substep_remask_plate_inventories[1:end] .- 
            substep_initial_plate_inventories[1:end] .+ substep_subducted_sediment_volumes ) )# ./
            #sub_time_step )
        logging_println("")
    end
    
    if enable_orogeny
        orogenic_uplift_rates = orogeny( subduction_footprint ) # uses ocean_subduct_sediment_plate_thickness 
        # returns a global field, not just continental, 
        # because plot points change sometimes
        apply_orogeny_fluxes_to_world( orogenic_uplift_rates[:,:,0] )
    end
    #=
    orogenic_erosion_rate_field = 
        generate_orogenic_erosion_fluxes( )
        # also global
    net_crust_change_rate = orogenic_uplift_rates[:,:,0] .- 
        orogenic_erosion_rate_field
    apply_orogeny_fluxes_to_world( net_crust_change_rate )
    =#
    is_continent_crust = eq_mask(world.crust_type,continent_crust)
    smooth_masked_field!(world.crust_thickness, is_continent_crust, orogeny_smooth_coeff)
    isostacy() # for erosion rate parameterization
    
    if enable_watch_orogeny_substeps
        scene = plot_field((world.crust_thickness .- 16500) .* is_continent_crust ./ 1000. )
        plot_add_plate_boundaries!(scene)
        plot_add_continent_outlines!(scene)
        display(scene)
        sleep(2.)
    end
    apply_tectonics_changes_to_plates( orogenic_uplift_rates[:,:,0] )    
    #=if enable_step_tectonics_plate_diagnostics
        substep_applied_plate_inventories = global_plate_sediment_inventories()
        logging_println("             applied ", substep_applied_plate_inventories)
    end=#

    #smooth_continental_crust_thickness()
    #tcix,tciy = highest_xy_field(world.crust_thickness)
    #thickest_crust = world.crust_thickness[tcix,tciy]
    #uplift = g
    logging_println("    thickest, highest crust ", 
        [ field_max( world.crust_thickness ), field_max( world.freeboard )])
    logging_println("    fastest cont / subd uplift ", 
        [ field_max( orogenic_uplift_rates[:,:,1] ), field_max( orogenic_uplift_rates[:,:,2] )])
    #ulix,uliy = highest_xy_field( uplift )
    #erosion = get_diag("crust_erosion_rate")
    #logging_println("uplift fastest, net ", [field_max(uplift),field_max(erosion)])

    #=if enable_step_tectonics_world_diagnostics
        ocn_subd = subducted_ocean_sediment_volumes ./ 
            ( i_plate_substep * sub_time_step )
        cont_subd = subducted_land_sediment_volumes ./ 
            ( i_plate_substep * sub_time_step )
        logging_println("    accum subduc  ", 
        [ ocn_subd, cont_subd, ocn_subd .+ cont_subd ])
    end=#



    #  end # substep loop



    # look for plates that change their identities, move stuff between them
    # sets plate.tectonics to record plateID changes
    #delete_old_plates()
    set_diag("continent_orogenic_uplift_rate",orogenic_uplift_rates[:,:,1] .* 
        eq_mask(world.crust_type,continent_crust) ) 
    set_diag("subduction_orogenic_uplift_rate",orogenic_uplift_rates[:,:,2] .* 
        eq_mask(world.crust_type,continent_crust) )  # meters / Myr, snapshot not accumulated
    world.subducted_land_sediment_volumes = subducted_land_sediment_volumes ./ 
        main_time_step
    world.subducted_ocean_sediment_volumes = subducted_ocean_sediment_volumes ./ 
        main_time_step
    #set_frac_diag("ocean_subduct_sediment_fraction_thickness",
    #    accumulated_ocean_sediment_subduction ./ main_time_step )
    #set_frac_diag("continent_subduct_sediment_fraction_thickness",
    #    accumulated_continent_sediment_subduction ./ main_time_step )

    #resolve_rotation_matrices()
    # set rotation matrices to the new time
    substepped_world_inventories = fill(0.,0:n_sediment_types)
    if enable_step_tectonics_world_diagnostics || eliminate_regrid_drift
        substepped_world_inventories = world_sediment_inventories(  )
    end

    if enable_cont_update_from_files
        update_world_continents_from_file()
        #apply_tectonics_changes_to_plates()  done later in apply_geomorph_changes_to_plates
    end

    updated_ID_world_inventories = fill(0.,0:n_sediment_types)
    if enable_step_tectonics_world_diagnostics || eliminate_regrid_drift
        updated_ID_world_inventories = world_sediment_inventories(  )
    end
    updated_ID_plate_inventories = fill(0.,0:n_sediment_types)
    if enable_step_tectonics_plate_diagnostics
        updated_ID_plate_inventories = global_plate_sediment_inventories()
    end
    subduction_rates = fill(0.,0:n_sediment_types)
    subduction_rates[1:n_sediment_types] .= 
        world.subducted_land_sediment_volumes .+ 
        world.subducted_ocean_sediment_volumes
    for i_sedtype in 1:n_sediment_types
        subduction_rates[0] += subduction_rates[i_sedtype]
    end
    uplifted_sediment = fill(0.,0:n_sediment_types)
    demoted_sediment = fill(0.,0:n_sediment_types)
    if enable_step_tectonics_world_diagnostics || eliminate_regrid_drift
        uplifted_sediment[1:n_sediment_types] = 
            volume_fields(get_frac_diag("ocean_2_continent_sediment_fraction_displaced"))
        demoted_sediment[1:n_sediment_types] = 
            volume_fields(get_frac_diag("continent_2_ocean_sediment_fraction_displaced"))
        for i_sedtype in 1:n_sediment_types
            uplifted_sediment[0] += uplifted_sediment[i_sedtype]
            demoted_sediment[0] += demoted_sediment[i_sedtype]
        end
    end
    final_world_inventories = updated_ID_world_inventories # for init if not eliminate_regrid_drift
    if eliminate_regrid_drift
        #world_cum_change = ( substepped_world_inventories .- initial_world_inventories ) ./
        #    main_time_step 
        #tweaked_inv = substepped_world_inventories .- 
        #    ( subduction_rates .+ world_cum_change .+ uplifted_sediment .+ demoted_sediment ) .* main_time_step
        tweaked_inv = initial_world_inventories .- 
            ( subduction_rates .+ uplifted_sediment .+ demoted_sediment ) .* main_time_step
        imbalance_ratios = fill( 1., n_sediment_types )
        for i_sedtype in 1:n_sediment_types
            if updated_ID_world_inventories[i_sedtype] > 0
                imbalance_ratios[i_sedtype] = tweaked_inv[i_sedtype] ./ 
                    updated_ID_world_inventories[i_sedtype]
            end
        end
        scale_global_sediment_components( imbalance_ratios[1:n_sediment_types] )
        final_world_inventories = world_sediment_inventories(  )
        logging_println("")
        #logging_println("     world init ", initial_world_inventories)
        #logging_println("      tectonics ", after_tectonics_world_inventories)
        substepped_world_inventories
        logging_println("  fudge factors ", imbalance_ratios)
        logging_println("         fudged ", final_world_inventories)
        logging_println("           slop ", final_world_inventories .- updated_ID_world_inventories)
        logging_println("  cum world bal ", 
            ( final_world_inventories .- initial_world_inventories ) ./
            main_time_step .+ subduction_rates .+ uplifted_sediment .+ demoted_sediment )
    end

    world_cum_change = ( final_world_inventories .- initial_world_inventories ) ./
        main_time_step 
    if enable_step_tectonics_plate_diagnostics
        final_plate_inventories = global_plate_sediment_inventories()
        plate_cum_change = ( final_plate_inventories .- initial_plate_inventories ) ./
            main_time_step
        plate_bal = plate_cum_change .+ 
            subduction_rates .+ # +ve because accounts for a loss
            uplifted_sediment .+ demoted_sediment # loss that will be go into ocean
        logging_println()
        logging_println("        plates init ", initial_plate_inventories)
        logging_println("              final ", final_plate_inventories)
        logging_println("                chg ", final_plate_inventories .- initial_plate_inventories )
        #logging_println("             remask ", remasked_plate_inventories)
        #logging_println("               rate ", ( remasked_plate_inventories .- update_ID_plate_inventories) / main_time_step )
        logging_println("          subducted ", subduction_rates)
        logging_println("   uplift scattered ", uplifted_sediment)
        logging_println("  demoted scattered ", demoted_sediment)
        logging_println("             change ", plate_cum_change)
        logging_println("                bal ", plate_bal )
    end
    if enable_step_tectonics_world_diagnostics
        world_bal = world_cum_change .+ subduction_rates .+ 
            uplifted_sediment .+ demoted_sediment
        logging_println()
        logging_println("         world init ", initial_world_inventories)
        logging_println("         substepped ", substepped_world_inventories)
        logging_println("  changing plateIDs ", updated_ID_world_inventories)
        logging_println("             fudged ", final_world_inventories)
        #logging_println("         land ", world_land_sediment_inventories())
        #logging_println("        ocean ", world_ocean_sediment_inventories())
        logging_println("               rate ", ( final_world_inventories .- initial_world_inventories ) ./ main_time_step )
        logging_println("          subducted ", subduction_rates)
        logging_println("   uplift scattered ", uplifted_sediment)
        logging_println("  demoted scattered ", demoted_sediment)
        logging_println("      cum world chg ", world_cum_change )
        logging_println("      cum world bal ", world_bal )
        logging_println()
    end
    #error("stopping")
   
end
function step_geomorph(  ) # requires previous run with tectonics to fill world diag arrays

    #= quickly skip tectonics for debugging
    step_tectonics(  )
    world.elevation_offset = ocean_thermal_boundary_layer()
    world.crust_age .+= main_time_step
    =#

    if enable_step_geomorph_diagnostics
        incoming_land_fraction_inventories = world_land_sediment_inventories( )
    end
    logging_println()
    logging_println("Geomorphology"); logging_println()
    logging_println("  Land Bulk Sediment Transport")

    clear_geomorph_process_arrays() # only necessary in standalone mode

    # continent and subduction uplift rates into *_orogenic_uplift_rate diags

    isostacy() 
    #check_reburial_exposed_basement() 

    # set variables so they will survive the while loop ending
    
    # world.geomorphology initially set by plate interpolation, 
    # resetting points submerged by sea level rise or fall
    for ix in 1:nx
        for iy in 1:ny
            if world.freeboard[ix,iy] > 0. 
                if world.geomorphology[ix,iy] < sedimented_land # previously ocean
                    world.geomorphology[ix,iy] = sedimented_land
                end
            elseif world.freeboard[ix,iy] > -100.
                world.geomorphology[ix,iy] = ocean_shelf # nondepositing, ocean or cont
            else
                #error([ix,iy,world.freeboard[ix,iy]])
                world.geomorphology[ix,iy] = pelagic_seafloor
                #=if world.crust_type[ix,iy] != ocean_crust
                    error([ix,iy,world.freeboard[ix,iy]])
                end=#
            end
        end
    end

    diffusive_mask = eq_mask(world.geomorphology, sedimented_land)
    orogenic_mask = eq_mask(world.geomorphology, exposed_basement)
    #nondepositing_shelf_mask = eq_mask(world.geomorphology,ocean_shelf)
    submerged_mask = lt_mask( world.freeboard, 0. )
    if enable_land_runoff == false
        submerged_mask .= 0
    end
    subaereal_mask = ( 1 .- submerged_mask )
    #depositing_ocean_mask = submerged_mask .* 
    #    eq_mask(world.crust_type,ocean_crust)# .* 
    #    eq_mask(world.geomorphology,pelagic_seafloor)
    
    orogenic_erosion_rate_field = fill(0.,nx,ny,0:n_sediment_types)
    orogenic_sediment_source_fields = fill(0.,nx,ny,0:n_sediment_types)
    denuded_sediment_source_fields = fill(0.,nx,ny,0:n_sediment_types)
    orogenic_boundary_source_fields = fill(0.,nx,ny,0:n_sediment_types)
    denuded_boundary_source_fields = fill(0.,nx,ny,0:n_sediment_types)
    aolean_erosion_fields = fill(0.,nx,ny,0:n_sediment_types)
    aolean_deposition_fields = fill(0.,nx,ny,0:n_sediment_types)
    sediment_source_fields = fill(0.,nx,ny,0:n_sediment_types)

    land_sediment_deposition_rate_field = fill(0.,nx,ny)
    land_sediment_fraction_dissolution_rate_fields = fill(0.,nx,ny,n_sediment_types)
    original_elevation_field = world.freeboard
    new_elevation_field = fill(0.,nx,ny)
    visualize_denuding_landscape_field = fill(0.,nx,ny)
    bulk_runoff_flux = 0.
    n_denuded = -1; need_another_loop = true; last_deposition_rate = 0.
    n_denuded_tot = 0; n_loops_already = 0
    while need_another_loop == true
        n_loops_already += 1
        diffusive_mask = eq_mask(world.geomorphology, sedimented_land)
        orogenic_mask = eq_mask(world.geomorphology, exposed_basement)
        orogenic_erosion_rate_field = generate_orogenic_erosion_fluxes() .* 
            gt_mask(world.freeboard,0.)
        for i_sedtype in 1:n_sediment_types
            orogenic_sediment_source_fields[:,:,i_sedtype] = orogenic_erosion_rate_field .*
                orogenic_sediment_source_fractions[i_sedtype] .*
                rho_continent_crust ./ rho_sediment 
        end
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
 
        sediment_source_fields = ( orogenic_boundary_source_fields .+
            denuded_boundary_source_fields ) .* subaereal_mask
        #sediment_source_fields -= aolean_erosion_fields
        #sediment_source_fields -= land_sediment_fraction_dissolution_rate_fields
        update_flux_totals!( sediment_source_fields ) 
        #println(" sediment_source_fields ",sediment_source_fields)

        new_elevation_field, land_sediment_deposition_rate_field = 
            land_bulk_sediment_transport( original_elevation_field, 
                diffusive_mask, sediment_source_fields[:,:,0], submerged_mask )

        n_denuded = check_for_bedrock_exposure( land_sediment_deposition_rate_field )
        n_denuded_tot = get_maskfield_census(eq_mask(world.geomorphology,exposed_basement))
        if enable_watch_denuding_landscape
            for ix in 1:nx
                for iy in 1:ny
                    if denuded_sediment_source_fields[ix,iy,0] == 0. && # pop to 1 
                        world.geomorphology[ix,iy] == exposed_basement
                        visualize_denuding_landscape_field[ix,iy] = 1
                    else
                        visualize_denuding_landscape_field[ix,iy] *= 0.9 # gradually fade away
                    end
                end
            end
            println("plotting denuding_landscape")
            scene = plot_field(visualize_denuding_landscape_field)
            plot_add_continent_outlines!(scene)
            plot_add_orogenies!(scene)
            display(scene)
            #sleep(2)
        end
            
        for ix in 1:nx
            for iy in 1:ny
                if world.geomorphology[ix,iy] == exposed_basement
                    land_sediment_deposition_rate_field[ix,iy] =
                        - world.sediment_thickness[ix,iy] / main_time_step
                end
            end
        end
        balance = 0.
        if enable_step_geomorph_diagnostics
            bulk_source_flux = volume_field( sediment_source_fields[:,:,0] .* subaereal_mask )
            bulk_deposition_flux = volume_field( land_sediment_deposition_rate_field )
            bulk_denuded_flux = volume_field( denuded_boundary_source_fields[:,:,0] )
            # note this includes the flux to the ocean, because it has to compensate the
            # negative land deposition rates
            #aolean_flux = volume_field( aolean_erosion_field )
            #crust_conv = get_
            bulk_runoff_flux = volume_field( land_bulk_runoff_fluxes( new_elevation_field, 
                diffusive_mask, submerged_mask ))
            #bulk_aolean_flux = volume_fields( aolean_erosion_fields )[1]
            #bulk_diss_flux = volume_fields( land_sediment_fraction_dissolution_rate_fields )[1]
            bulk_cont2ocn_scatter = 0.
            for i_sedtype in 1:n_sediment_types
                bulk_cont2ocn_scatter += volume_field(get_frac_diag(
                    "continent_2_ocean_sediment_fraction_displaced")[:,:,i_sedtype])
            end
            balance = bulk_source_flux - bulk_denuded_flux - # because it had to be included in source
                bulk_deposition_flux - bulk_runoff_flux # - bulk_diss_flux - bulk_aolean_flux #- bulk_cont2ocn_scatter # 
            #mass_balance = 
            #logging_println("n, src, nud, dep, runoff, bal ",[n_denuded_tot,bulk_source_flux,
            #    bulk_denuded_flux,
            #    bulk_deposition_flux,bulk_runoff_flux, balance])
            logging_println("n, bal ",[n_denuded_tot, balance])
        else
            logging_println("    denuded ",n_denuded_tot)
        end

        need_another_loop = false
        if n_denuded > 0
            need_another_loop = true
            #println("because n_denuded")
        end
        #=if abs(balance) > 1.e6
            need_another_loop = true
            #println("because deposition ", abs( deposition - last_deposition_rate ) )
    #            last_deposition_rate = bulk_deposition_flux
        end=#
        if n_loops_already > 20
            need_another_loop = false
        end
    end
    set_diag( "crust_erosion_rate", orogenic_erosion_rate_field )
    set_frac_diag( "crust_orogenic_fraction_flux", orogenic_sediment_source_fields )
    set_frac_diag( "land_orogenic_fraction_flux", orogenic_boundary_source_fields .* subaereal_mask )
    set_frac_diag( "coastal_orogenic_fraction_flux", orogenic_boundary_source_fields .* submerged_mask )
    set_frac_diag("denuded_sediment_source_fraction_flux", denuded_sediment_source_fields[:,:,1:n_sediment_types])
    set_frac_diag("denuded_land_boundary_fraction_flux", 
        denuded_boundary_source_fields[:,:,1:n_sediment_types] .* subaereal_mask)
    set_frac_diag("denuded_coastal_boundary_fraction_flux", 
        denuded_boundary_source_fields[:,:,1:n_sediment_types] .* submerged_mask)
    set_diag("land_sediment_deposition_rate",land_sediment_deposition_rate_field)

    area_denuded = volume_field(eq_mask(world.geomorphology,exposed_basement))
    area_covered = volume_field(eq_mask(world.geomorphology,sedimented_land) )
    logging_println("    fraction denuded ", area_denuded / ( area_denuded + area_covered ))

    #println(""); println("Land Sediment Composition Calculation");println("")

    land_sediment_fraction_deposition_rate_fields = 
        land_sediment_fraction_transport( new_elevation_field, 
            land_sediment_deposition_rate_field,
            diffusive_mask, sediment_source_fields, submerged_mask ) 
                
    set_frac_diag("land_sediment_fraction_deposition_rate",
        land_sediment_fraction_deposition_rate_fields)

    coastal_sediment_fraction_runoff_flux_fields = 
        land_fraction_runoff_fluxes( new_elevation_field, 
            land_sediment_fraction_deposition_rate_fields,
            diffusive_mask, submerged_mask ) 
    set_frac_diag("coastal_sediment_fraction_runoff_flux",
        coastal_sediment_fraction_runoff_flux_fields)        

    uplifted_ocean_sediment = fill(0.,nx,ny,0:n_sediment_types)
    uplifted_ocean_sediment[:,:,1:n_sediment_types] =
        get_frac_diag("ocean_2_continent_sediment_fraction_displaced")
    update_flux_totals!(uplifted_ocean_sediment)

    demoted_land_sediment = fill(0.,nx,ny,0:n_sediment_types)
    demoted_land_sediment[:,:,1:n_sediment_types] =
        get_frac_diag("continent_2_ocean_sediment_fraction_displaced")
    update_flux_totals!(demoted_land_sediment)    
    
    submerged_slop = ( uplifted_ocean_sediment .+ demoted_land_sediment ) .*
        ( submerged_mask )
    landed_slop = ( uplifted_ocean_sediment .+ demoted_land_sediment ) .*
        ( 1 .- submerged_mask )
    
    redeposited_landed = 
        distribute_fluxes_uniformly_outside_boundary( 
            landed_slop, ( 1 .- submerged_mask ) )
            # an ocean point turns into continent, the sediment is
            # distributed into the same submerged coastal grid points
            # as the coastal runoff, adjacent to the diffusive land regions

    redeposited_submerged = 
        distribute_fluxes_uniformly_inside_boundary( 
            submerged_slop, submerged_mask )     
            # a land point suddenly finds itself oceanic.  its sediment
            # is also distributed in the coastal ocean grid cells. 
            
    if enable_extraverbose_step_geomorph_diagnostics
        uplifted = volume_fields(uplifted_ocean_sediment)
        demoted = volume_fields(demoted_land_sediment)
        landed = volume_fields(landed_slop)
        subm = volume_fields(submerged_slop)
        relanded = volume_fields(redeposited_landed)
        resubm = volume_fields(redeposited_submerged)
        logging_println("uplifted ", uplifted )
        logging_println(" demoted ", demoted )
        logging_println("     tot ", uplifted + demoted )
        logging_println(" on land ", landed)
        logging_println(" submerg ", subm)
        logging_println("     tot ", landed + subm )
        logging_println("   redep ", relanded)
        logging_println("   redep ", resubm)
        logging_println("     tot ", relanded + resubm )
        logging_println("")
    end

    new_sediment_thickness, new_sediment_surface_fractions = 
        apply_land_sediment_fluxes( land_sediment_fraction_deposition_rate_fields )
    world.sediment_thickness[:,:] = new_sediment_thickness
    world.sediment_surface_fractions[:,:,:] = new_sediment_surface_fractions
    world.crust_thickness .-= orogenic_erosion_rate_field .* main_time_step

    aolean_erosion_fields, aolean_deposition_fields = aolean_transport()
    set_frac_diag("aolean_erosion_fraction_flux",aolean_erosion_fields)
    set_frac_diag("aolean_deposition_fraction_flux",aolean_deposition_fields)
    accum_frac_diag("land_sediment_fraction_deposition_rate",
        -1. .* aolean_erosion_fields[:,:,1:end])
    accum_diag("land_sediment_deposition_rate",
        -1. .* aolean_erosion_fields[:,:,0])

    runoff_map = generate_runoff_map()
    set_diag("land_Q_runoff_field", runoff_map)
    land_sediment_fraction_dissolution_rate_fields, land_orogenic_Ca_source_rates = 
        subaereal_sediment_dissolution( runoff_map )
    set_frac_diag("land_sediment_fraction_dissolution_rate",
        land_sediment_fraction_dissolution_rate_fields)
    set_diag("land_orogenic_Ca_source_rates", land_orogenic_Ca_source_rates)

    land_sediment_weathering_index = fill(0.,nx,ny)
    CaO_CaCO3_free = fill(0.,nx,ny)
    #CaO_calcite_free = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if world.sediment_thickness[ix,iy] > 0.
                if world.sediment_surface_fractions[ix,iy,CaO_sediment] + 
                    world.sediment_surface_fractions[ix,iy,clay_sediment] > 0.
                    CaO_CaCO3_free[ix,iy] = world.sediment_surface_fractions[ix,iy,CaO_sediment]  / 
                        ( world.sediment_surface_fractions[ix,iy,CaO_sediment] 
                        + world.sediment_surface_fractions[ix,iy,clay_sediment] ) 
                    land_sediment_weathering_index[ix,iy] = 1. - 
                        CaO_CaCO3_free[ix,iy] / orogenic_sediment_source_fractions[CaO_sediment] 
                end
            end
        end
    end
    set_diag("land_sediment_weathering_index", land_sediment_weathering_index)
    #set_diag("land_sediment_CaO_CaCO3_free", CaO_CaCO3_free)
   
    logging_println("")
    logging_println("Ocean Sedimentation")
    #println("CaCO3 system")

    continental_CaCO3_dissolution_rate = volume_field(
        land_sediment_fraction_dissolution_rate_fields[:,:,CaCO3_sediment] )
    continental_CaO_dissolution_rate = volume_field(
        land_sediment_fraction_dissolution_rate_fields[:,:,CaO_sediment] )
    land_orogenic_Ca_source_rate = volume_field( land_orogenic_Ca_source_rates )
    total_Ca_sources = land_orogenic_Ca_source_rate + 
        continental_CaCO3_dissolution_rate + 
        continental_CaO_dissolution_rate
    logging_println(" Ca src tot, oro, sed CaCO3, Cao ", 
        [ total_Ca_sources,
        land_orogenic_Ca_source_rate, 
        continental_CaCO3_dissolution_rate, 
        continental_CaO_dissolution_rate,
     ] )

    ocean_CaCO3_deposition_rate = total_Ca_sources 
    #global_CaCO3_net_burial_flux # m3 / Myr
    #if enable_CaCO3_land_diss_ocean_repcp
    #    ocean_CaCO3_deposition_rate += land_dissolved_Ca_prod_rate
    #end
    ocean_CO3 = find_steady_state_ocean_CO3( ocean_CaCO3_deposition_rate )

    coastal_CaCO3_deposition_field = get_coastal_CaCO3_deposition_field( ocean_CO3 ) 
    set_diag("coastal_CaCO3_flux", coastal_CaCO3_deposition_field )
    pelagic_CaCO3_deposition_field = get_pelagic_CaCO3_deposition_field( ocean_CO3 )
    set_diag("pelagic_CaCO3_deposition_rate", pelagic_CaCO3_deposition_field )
    continental_CaCO3_deposition_field = get_continental_CaCO3_deposition_field( ocean_CO3 )
    set_diag("continental_CaCO3_deposition_rate", continental_CaCO3_deposition_field )

    accum_diag("land_sediment_deposition_rate",continental_CaCO3_deposition_field)
    accum_frac_diag("land_sediment_fraction_deposition_rate",CaCO3_sediment,continental_CaCO3_deposition_field)


    if enable_step_geomorph_diagnostics
        logging_println("coast,pelag,cont ",
            [ volume_field(coastal_CaCO3_deposition_field), 
            volume_field(pelagic_CaCO3_deposition_field), 
            volume_field(continental_CaCO3_deposition_field) ])
    end

    #logging_println("  sediment distribution calculation")

    #=clay_dep = get_frac_diag("coastal_sediment_fraction_runoff_flux",clay_sediment) .+ 
        get_frac_diag("denuded_coastal_boundary_fraction_flux",clay_sediment) .+
        get_diag("coastal_orogenic_clay_flux") .+ 
        get_diag("aolean_clay_deposition_rate")
    CaCO3_dep = get_diag("coastal_CaCO3_flux") .+  # CaCO3 fluxes
        get_diag("pelagic_CaCO3_deposition_rate") .+ 
        #get_diag("continental_CaCO3_deposition_rate") .+ 
        get_frac_diag("coastal_sediment_fraction_runoff_flux",CaCO3_sediment) .+
        get_frac_diag("denuded_coastal_boundary_fraction_flux",CaCO3_sediment)
    CaO_dep =#
    incoming_fluxes = fill(0.,nx,ny,0:n_sediment_types)
    incoming_fluxes[:,:,1:end] = get_frac_diag("coastal_sediment_fraction_runoff_flux") .+
        get_frac_diag("denuded_coastal_boundary_fraction_flux") .+
        get_frac_diag("coastal_orogenic_fraction_flux") .+
        get_frac_diag("aolean_deposition_fraction_flux")
    incoming_fluxes[:,:,CaCO3_sediment] .+= 
        get_diag("coastal_CaCO3_flux") .+ 
        get_diag("pelagic_CaCO3_deposition_rate")
    for i_sedtype in 1:n_sediment_types
        incoming_fluxes[:,:,0] .+= incoming_fluxes[:,:,i_sedtype]
    end
    
    #fill(0.,nx,ny,0:n_sediment_types) # meters / Myr
    #=incoming_fluxes[:,:,1] = clay_dep
    incoming_fluxes[:,:,2] = CaCO3_dep
    incoming_fluxes[:,:,0] = clay_dep .+ CaCO3_dep=#
    incoming_fluxes += redeposited_landed .+ 
        redeposited_submerged
    #offshelf_fluxes =  distribute_fluxes_uniformly_inside_boundary( 
    #    incoming_fluxes,  )     

    set_frac_diag( "ocean_sediment_fraction_influx",incoming_fluxes)
        #=clay_sediment, clay_dep )
    set_frac_diag( "ocean_sediment_fraction_influx",
        CaCO3_sediment, CaCO3_dep )=#
    println("distributing ocean fluxes")
    accumulating_depositing_fluxes = distribute_ocean_sediment_fluxes( incoming_fluxes )

    set_diag("seafloor_sediment_deposition_rate",
        accumulating_depositing_fluxes[:,:,0])
    for i_sedtype in 1:n_sediment_types
        set_frac_diag("seafloor_sediment_fraction_deposition_rate",i_sedtype,
            accumulating_depositing_fluxes[:,:,i_sedtype])
    end
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

    if enable_step_geomorph_diagnostics
        old_ocean_sed_inventories = world_ocean_sediment_inventories( )
    end

    apply_ocean_sediment_fluxes( )
    # updates world.sediment* for ocean points from seafloor_sediment*_deposition_rate
    
    land_sediment_fraction_tweak_rate_fields = fill(0.,nx,ny,0:n_sediment_types)
    land_sediment_fraction_tweak_rate_fields[:,:,1:end] .-= 
        get_frac_diag("aolean_erosion_fraction_flux") .+
        get_frac_diag("land_sediment_fraction_dissolution_rate")
    land_sediment_fraction_tweak_rate_fields[:,:,CaCO3_sediment] .+= 
        get_diag("continental_CaCO3_deposition_rate")

    new_sediment_thickness, new_sediment_surface_fractions = 
        apply_land_sediment_fluxes( land_sediment_fraction_tweak_rate_fields )
    world.sediment_thickness[:,:] = new_sediment_thickness
    world.sediment_surface_fractions[:,:,:] = new_sediment_surface_fractions

    if enable_step_geomorph_diagnostics
        source_fluxes = fill(0.,0:n_sediment_types)
        orogen_fluxes = fill(0.,0:n_sediment_types)
        runoff_fluxes = fill(0.,0:n_sediment_types)
        aolean_erosion_fluxes = fill(0.,0:n_sediment_types)
        denuded_fluxes_to_land = fill(0.,0:n_sediment_types)
        denuded_fluxes_to_coast = fill(0.,0:n_sediment_types)
        dissolution_fluxes = fill(0.,0:n_sediment_types)
        total_deposition_fluxes = fill(0.,0:n_sediment_types)
        transport_deposition_fluxes = fill(0.,0:n_sediment_types)
        cont2ocean_fluxes = fill(0.,0:n_sediment_types)
        ocean2cont_fluxes = fill(0.,0:n_sediment_types)
        orogen_fluxes[1:end] = volume_fields(get_frac_diag("land_orogenic_fraction_flux"))
        #orogen_fluxes[CaO_sediment] = volume_field(get_diag("land_orogenic_CaO_flux"))
        dissolution_fluxes[1:end] = volume_fields(get_frac_diag("land_sediment_fraction_dissolution_rate"))
        aolean_erosion_fluxes[1:end] = volume_fields(get_frac_diag("aolean_erosion_fraction_flux"))
        runoff_fluxes[1:n_sediment_types] = 
            volume_fields(get_frac_diag("coastal_sediment_fraction_runoff_flux"))
        total_deposition_fluxes[1:n_sediment_types] = 
            volume_fields(get_frac_diag("land_sediment_fraction_deposition_rate"))
        continent_CaCO3_deposition_flux = volume_field(get_diag("continental_CaCO3_deposition_rate"))
        transport_deposition_fluxes[:] = total_deposition_fluxes[:]
        transport_deposition_fluxes[CaCO3_sediment] -= continent_CaCO3_deposition_flux
        denuded_fluxes_to_land[1:n_sediment_types] = 
            volume_fields(get_frac_diag("denuded_land_boundary_fraction_flux") )
        denuded_fluxes_to_coast[1:n_sediment_types] = 
            volume_fields(get_frac_diag("denuded_coastal_boundary_fraction_flux") )
        cont2ocean_fluxes[1:n_sediment_types] = volume_fields(demoted_land_sediment[:,:,1:n_sediment_types])
        ocean2cont_fluxes[1:n_sediment_types] = volume_fields(uplifted_ocean_sediment[:,:,1:n_sediment_types])
        # note that this includes stuff to be deposited in the ocean, because the
        # land deposition rates reflect it, part of the mass balance for now
        for i_sedtype in 1:n_sediment_types
            orogen_fluxes[0] += orogen_fluxes[i_sedtype]
            dissolution_fluxes[0] += dissolution_fluxes[i_sedtype]
            total_deposition_fluxes[0] += total_deposition_fluxes[i_sedtype]
            transport_deposition_fluxes[0] += transport_deposition_fluxes[i_sedtype]
            aolean_erosion_fluxes[0] += aolean_erosion_fluxes[i_sedtype]
            denuded_fluxes_to_land[0] += denuded_fluxes_to_land[i_sedtype]
            denuded_fluxes_to_coast[0] += denuded_fluxes_to_coast[i_sedtype]
            runoff_fluxes[0] += runoff_fluxes[i_sedtype]
            cont2ocean_fluxes[0] += cont2ocean_fluxes[i_sedtype]
            ocean2cont_fluxes[0] += ocean2cont_fluxes[i_sedtype]
        end
        new_land_fraction_inventories = world_land_sediment_inventories( ) 
        change_rate = ( new_land_fraction_inventories .- incoming_land_fraction_inventories ) /  # incoming_land_fraction_inventories ) / 
            main_time_step
        source_fluxes = orogen_fluxes .- aolean_erosion_fluxes .- dissolution_fluxes 

        mass_balances = change_rate .+ dissolution_fluxes .- total_deposition_fluxes
        flux_balances = source_fluxes .- transport_deposition_fluxes .- runoff_fluxes .- 
            denuded_fluxes_to_coast .+ dissolution_fluxes

        logging_println("")
        logging_println("land budget init ",incoming_land_fraction_inventories)
        logging_println("           final ",new_land_fraction_inventories)
        logging_println("     change rate ",change_rate)
        logging_println("   transport dep ",transport_deposition_fluxes)
        logging_println("       total dep ",total_deposition_fluxes)
        logging_println("   orogen source ",orogen_fluxes)
        logging_println("       land diss ",dissolution_fluxes)
        logging_println(" denuded to land ",denuded_fluxes_to_land)
        logging_println("denuded to coast ",denuded_fluxes_to_coast)
        logging_println("cont 2 ocn swtch ",cont2ocean_fluxes)
        logging_println("ocn 2 cont swtch ",ocean2cont_fluxes)
        logging_println("     aolean erod ",aolean_erosion_fluxes)
        logging_println("    total source ",source_fluxes)
        logging_println("          runoff ",runoff_fluxes)
        logging_println("phs flx to ocean ",runoff_fluxes .+ aolean_erosion_fluxes)
        logging_println("tot flx to ocean ",runoff_fluxes .+ aolean_erosion_fluxes .+
            cont2ocean_fluxes .+ ocean2cont_fluxes)
        logging_println("        mass bal ",mass_balances)
        logging_println("        flux bal ",flux_balances)
    end

    if enable_step_geomorph_diagnostics
        new_ocean_sed_inventories = world_ocean_sediment_inventories( )
        total_sources = fill(0.,0:n_sediment_types)
        land_sources = fill(0.,0:n_sediment_types)
        ocean_sed_depo_rates = fill(0.,0:n_sediment_types)

        orogenic_sources = fill(0.,0:n_sediment_types)
        orogenic_sources[1:end] = volume_fields( get_frac_diag("coastal_orogenic_fraction_flux") )

        aolean_deposition = fill(0.,0:n_sediment_types)
        aolean_deposition[1:end] = volume_fields( get_frac_diag("aolean_deposition_fraction_flux") )

        pelagic_sources = fill(0.,0:n_sediment_types)
        #pelagic_sources[1:end] = volume_fields( get_frac_diag("aolean_deposition_fraction_flux") )
        pelagic_sources[CaCO3_sediment] += volume_field(get_diag("pelagic_CaCO3_deposition_rate"))

        coastal_sources = fill(0.,0:n_sediment_types)
        coastal_sources[CaCO3_sediment] = volume_field(get_diag("coastal_CaCO3_flux")) 

        continental_deposition = fill(0.,0:n_sediment_types)
        continental_deposition[1:end] = volume_fields(get_frac_diag("land_sediment_fraction_deposition_rate")) 

        continental_dissolution = fill(0.,0:n_sediment_types)
        continental_dissolution[1:end] = volume_fields(get_frac_diag("land_sediment_fraction_dissolution_rate")) 

        runoff_sources = fill(0.,0:n_sediment_types)
        runoff_sources[1:end] = volume_fields( get_frac_diag("coastal_sediment_fraction_runoff_flux") ) 
        #runoff_sources[CaCO3_sediment] = volume_field( get_frac_diag("coastal_sediment_fraction_runoff_flux",CaCO3_sediment) )

        ocn2cont_redist_sources = fill(0.,0:n_sediment_types)
        cont2ocn_redist_sources = fill(0.,0:n_sediment_types)
        denude_sources = fill(0.,0:n_sediment_types)
        #total_sources[clay_sediment] = runoff_sources[clay_sediment] + 
        #    orogenic_sources[clay_sediment] + pelagic_sources[clay_sediment] 
        #total_sources[CaCO3_sediment] = coastal_sources[CaCO3_sediment] + 
        #    pelagic_sources[CaCO3_sediment] + runoff_sources[CaCO3_sediment]
        ocn2cont_redist_sources[1:end] = 
            volume_fields( get_frac_diag("ocean_2_continent_sediment_fraction_displaced") )
        cont2ocn_redist_sources[1:end] =  
            volume_fields( get_frac_diag("continent_2_ocean_sediment_fraction_displaced") )
        denude_sources[1:end] = 
            volume_fields(get_frac_diag("denuded_coastal_boundary_fraction_flux"))
        ocean_sed_depo_rates[1:end] = volume_fields(
            get_frac_diag("seafloor_sediment_fraction_deposition_rate"))

        for i_sedtype in 1:n_sediment_types
            orogenic_sources[0] += orogenic_sources[i_sedtype]
            runoff_sources[0] += runoff_sources[i_sedtype]
            aolean_deposition[0] += aolean_deposition[i_sedtype]
            coastal_sources[0] += coastal_sources[i_sedtype]
            pelagic_sources[0] += pelagic_sources[i_sedtype]
            continental_deposition[0] += continental_deposition[i_sedtype]
            continental_dissolution[0] += continental_dissolution[i_sedtype]
            ocn2cont_redist_sources[0] += ocn2cont_redist_sources[i_sedtype]
            cont2ocn_redist_sources[0] += cont2ocn_redist_sources[i_sedtype]
            ocean_sed_depo_rates[0] += ocean_sed_depo_rates[i_sedtype]
            denude_sources[0] += denude_sources[i_sedtype]
        end 

        from_land_sources = runoff_sources + aolean_deposition
        phys_total_sources = from_land_sources + orogenic_sources + 
            coastal_sources + pelagic_sources
        full_total_sources = phys_total_sources + 
            cont2ocn_redist_sources + ocn2cont_redist_sources + 
            denude_sources
        inventory_change_rate = ( new_ocean_sed_inventories .- old_ocean_sed_inventories) /
            main_time_step
        mass_balance = inventory_change_rate - ocean_sed_depo_rates #+ 
            #ocn2cont_redist_sources # makes inv go down, so added
        flux_balance = full_total_sources .- ocean_sed_depo_rates 
        logging_println("")
        logging_println("ocean budget init ", old_ocean_sed_inventories)
        logging_println("              new ", new_ocean_sed_inventories)
        logging_println("           change ", inventory_change_rate)
        logging_println(" orogenic sources ", orogenic_sources)
        logging_println("           runoff ", runoff_sources)
        logging_println("       aolean dep ", aolean_deposition)
        logging_println("    tot from land ", from_land_sources)
        logging_println("          coastal ", coastal_sources)
        logging_println("          pelagic ", pelagic_sources)
        logging_println("         phys tot ", phys_total_sources)
        logging_println(" ocn 2 cont redst ", ocn2cont_redist_sources)
        logging_println(" cont 2 ocn redst ", cont2ocn_redist_sources)
        logging_println("  denude to coast ", denude_sources)
        logging_println("     full tot src ", full_total_sources)
        logging_println("              dep ", ocean_sed_depo_rates)
        logging_println("         mass bal ", mass_balance )
        logging_println("         flux bal ", flux_balance )
    end                    



    isostacy()
end
function step_budgets()
    #land_source_fluxes = fill(0.,0:n_sediment_types)
    land_orogen_fluxes = fill(0.,0:n_sediment_types)
    land_dissolution_fluxes = fill(0.,0:n_sediment_types)
    runoff_fluxes = fill(0.,0:n_sediment_types)
    aolean_erosion_fluxes = fill(0.,0:n_sediment_types)
    denuded_fluxes_to_land = fill(0.,0:n_sediment_types)
    denuded_fluxes_to_coast = fill(0.,0:n_sediment_types)
    demoted_land_sediment = fill(0.,0:n_sediment_types)
    uplifted_ocean_sediment = fill(0.,0:n_sediment_types)
    land_transport_deposition_fluxes = fill(0.,0:n_sediment_types)
    land_total_deposition_fluxes = fill(0.,0:n_sediment_types)
    #cont2ocean_fluxes = fill(0.,0:n_sediment_types)
    #ocean2cont_fluxes = fill(0.,0:n_sediment_types)

    #ocean_phys_total_sources = fill(0.,0:n_sediment_types)
    #ocean_land_sources = fill(0.,0:n_sediment_types)
    ocean_sed_depo_rates = fill(0.,0:n_sediment_types)
    aolean_deposition = fill(0.,0:n_sediment_types)
    ocean_orogenic_sources = fill(0.,0:n_sediment_types)
    ocean_pelagic_sources = fill(0.,0:n_sediment_types)
    ocean_coastal_sources = fill(0.,0:n_sediment_types)

    world_subduction_rates = fill(0.,0:n_sediment_types)
    world_source_fluxes = fill(0.,0:n_sediment_types)
    world_dissolution_fluxes = fill(0.,0:n_sediment_types)
    world_deposition_rates = fill(0.,0:n_sediment_types)

    land_orogen_fluxes[1:end] = volume_fields(get_frac_diag("land_orogenic_fraction_flux"))
    land_dissolution_fluxes[1:end] = volume_fields(get_frac_diag("land_sediment_fraction_dissolution_rate"))
    aolean_erosion_fluxes[1:end] = volume_fields(get_frac_diag("aolean_erosion_fraction_flux"))
    runoff_fluxes[1:n_sediment_types] = 
        volume_fields(get_frac_diag("coastal_sediment_fraction_runoff_flux"))
    land_total_deposition_fluxes[1:n_sediment_types] = 
        volume_fields(get_frac_diag("land_sediment_fraction_deposition_rate"))
    continent_CaCO3_deposition_flux = volume_field(get_diag("continental_CaCO3_deposition_rate"))
    land_transport_deposition_fluxes[:] = land_total_deposition_fluxes[:]
    land_transport_deposition_fluxes[CaCO3_sediment] -= continent_CaCO3_deposition_flux
    denuded_fluxes_to_land[1:n_sediment_types] = 
        volume_fields(get_frac_diag("denuded_land_boundary_fraction_flux") )
    denuded_fluxes_to_coast[1:n_sediment_types] = 
        volume_fields(get_frac_diag("denuded_coastal_boundary_fraction_flux") )
    demoted_land_sediment[1:n_sediment_types] =
        volume_fields(get_frac_diag("continent_2_ocean_sediment_fraction_displaced"))
    uplifted_ocean_sediment[1:n_sediment_types] =
        volume_fields(get_frac_diag("ocean_2_continent_sediment_fraction_displaced"))
    #cont2ocean_fluxes[1:n_sediment_types] = volume_fields(demoted_land_sediment[:,:,1:n_sediment_types])
    #ocean2cont_fluxes[1:n_sediment_types] = volume_fields(uplifted_ocean_sediment[:,:,1:n_sediment_types])

    ocean_orogenic_sources[1:end] = volume_fields( get_frac_diag("coastal_orogenic_fraction_flux") )
    aolean_deposition[1:end] = volume_fields( get_frac_diag("aolean_deposition_fraction_flux") )
    ocean_pelagic_sources[CaCO3_sediment] += volume_field(get_diag("pelagic_CaCO3_deposition_rate"))
    ocean_coastal_sources[CaCO3_sediment] = volume_field(get_diag("coastal_CaCO3_flux")) 
    ocean_sed_depo_rates[1:end] = volume_fields(
        get_frac_diag("seafloor_sediment_fraction_deposition_rate"))

    world_subduction_rates[1:n_sediment_types] .= 
        world.subducted_ocean_sediment_volumes .+
        world.subducted_land_sediment_volumes
    world_source_fluxes[1:end] = 
        volume_fields(get_frac_diag("crust_orogenic_fraction_flux")) # no denuded sediment
    world_source_fluxes[CaCO3_sediment] += 
        volume_field(get_diag("coastal_CaCO3_flux")) +
        volume_field(get_diag("pelagic_CaCO3_deposition_rate")) +
        volume_field(get_diag("continental_CaCO3_deposition_rate"))
        #volume_field(get_frac_diag("global_sediment_fraction_deposition_rate")[:,:,CaCO3_sediment])
    world_dissolution_fluxes[1:end] = 
        volume_fields(get_frac_diag("land_sediment_fraction_dissolution_rate")) # no denuded sediment
    world_deposition_rates[1:end] = volume_fields(
        get_frac_diag("global_sediment_fraction_deposition_rate"))


    # note that this includes stuff to be deposited in the ocean, because the
    # land deposition rates reflect it, part of the mass balance for now
    for i_sedtype in 1:n_sediment_types
        land_orogen_fluxes[0] += land_orogen_fluxes[i_sedtype]
        land_dissolution_fluxes[0] += land_dissolution_fluxes[i_sedtype]
        aolean_erosion_fluxes[0] += aolean_erosion_fluxes[i_sedtype]
        runoff_fluxes[0] += runoff_fluxes[i_sedtype]
        land_total_deposition_fluxes[0] += land_total_deposition_fluxes[i_sedtype]
        land_transport_deposition_fluxes[0] += land_transport_deposition_fluxes[i_sedtype]

        denuded_fluxes_to_land[0] += denuded_fluxes_to_land[i_sedtype]
        denuded_fluxes_to_coast[0] += denuded_fluxes_to_coast[i_sedtype]
        uplifted_ocean_sediment[0] += uplifted_ocean_sediment[i_sedtype]
        demoted_land_sediment[0] += demoted_land_sediment[i_sedtype]
        #cont2ocean_fluxes[0] += cont2ocean_fluxes[i_sedtype]
        #ocean2cont_fluxes[0] += ocean2cont_fluxes[i_sedtype]

        ocean_orogenic_sources[0] += ocean_orogenic_sources[i_sedtype]
        aolean_deposition[0] += aolean_deposition[i_sedtype]
        ocean_pelagic_sources[0] += ocean_pelagic_sources[i_sedtype]
        ocean_coastal_sources[0] += ocean_coastal_sources[i_sedtype]
        ocean_sed_depo_rates[0] += ocean_sed_depo_rates[i_sedtype]

        world_subduction_rates[0] += world_subduction_rates[i_sedtype]
        world_source_fluxes[0] += world_source_fluxes[i_sedtype]
        world_dissolution_fluxes[0] += world_dissolution_fluxes[i_sedtype]
        world_deposition_rates[0] += world_deposition_rates[i_sedtype]
    end

    land_initial_inventories = fill(0.,0:n_sediment_types)
    ocean_initial_inventories = fill(0.,0:n_sediment_types)
    #world_initial_inventories = fill(0.,0:n_sediment_types)
    land_initial_inventories[1:end] = world.initial_land_sediment_inventories
    ocean_initial_inventories[1:end] = world.initial_ocean_sediment_inventories
    #world_initial_inventories[1:end] = world.initial_land_sediment_inventories +
    #    world.initial_ocean_sediment_inventories
    for i_sedtype in 1:n_sediment_types
        land_initial_inventories[0] += land_initial_inventories[i_sedtype]
        ocean_initial_inventories[0] += ocean_initial_inventories[i_sedtype]
    end
    world_initial_inventories = land_initial_inventories + ocean_initial_inventories

    new_land_inventories = world_land_sediment_inventories( ) 
    new_ocean_inventories = world_ocean_sediment_inventories( )
    new_world_inventories = new_land_inventories + new_ocean_inventories

    land_change_rate = ( new_land_inventories .- land_initial_inventories ) /  # incoming_land_fraction_inventories ) / 
        main_time_step
    ocean_change_rate = ( new_ocean_inventories .- ocean_initial_inventories) /
        main_time_step
    world_change_rate = ( new_world_inventories .- world_initial_inventories ) /
        main_time_step

    land_source_fluxes = land_orogen_fluxes .- aolean_erosion_fluxes .- land_dissolution_fluxes 
    land_phys_flux_to_ocean = runoff_fluxes .+ aolean_erosion_fluxes
    land_tot_flux_to_ocean = land_phys_flux_to_ocean +
        demoted_land_sediment + uplifted_ocean_sediment + denuded_fluxes_to_coast
    ocean_phys_total_sources = ocean_orogenic_sources + runoff_fluxes + aolean_deposition +
        ocean_coastal_sources + ocean_pelagic_sources
    ocean_full_total_sources = ocean_phys_total_sources + 
        demoted_land_sediment + uplifted_ocean_sediment + 
        denuded_fluxes_to_coast
    

    land_mass_balances = land_change_rate .- land_dissolution_fluxes .+ 
        aolean_erosion_fluxes .- land_total_deposition_fluxes
    land_flux_balances = land_source_fluxes .- land_transport_deposition_fluxes .- 
        runoff_fluxes .- denuded_fluxes_to_coast .+ 
        aolean_erosion_fluxes .+ land_dissolution_fluxes
    ocean_mass_balances = ocean_change_rate - ocean_sed_depo_rates #+ 
    ocean_flux_balances = ocean_full_total_sources - ocean_sed_depo_rates 
    world_geomorph_mass_balances = world_change_rate - 
        world_deposition_rates - world_dissolution_fluxes
    world_geomorph_flux_balances = world_source_fluxes - 
        world_deposition_rates + world_dissolution_fluxes +
        demoted_land_sediment + uplifted_ocean_sediment

    logging_println("")
    logging_println("land budget init ",land_initial_inventories)
    logging_println("           final ",new_land_fraction_inventories)
    logging_println("     change rate ",land_change_rate)
    logging_println("   transport dep ",land_transport_deposition_fluxes)
    logging_println("       total dep ",land_total_deposition_fluxes)
    logging_println("   orogen source ",land_orogen_fluxes)
    logging_println("       land diss ",land_dissolution_fluxes)
    logging_println(" denuded to land ",denuded_fluxes_to_land)
    #logging_println("denuded to coast ",denuded_fluxes_to_coast)
    logging_println("demoted land sed ",demoted_land_sediment)
    logging_println("uplifted ocn sed ",uplifted_ocean_sediment)
    logging_println("     aolean erod ",aolean_erosion_fluxes)
    logging_println("    total source ",land_source_fluxes)
    logging_println("          runoff ",runoff_fluxes)
    logging_println("phs flx to ocean ",land_phys_flux_to_ocean)
    logging_println("tot flx to ocean ",land_tot_flux_to_ocean)
    logging_println("        mass bal ",land_mass_balances)
    logging_println("        flux bal ",land_flux_balances)
    
    logging_println("")
    logging_println("ocean budget init ", ocean_initial_inventories)
    logging_println("            final ", new_ocean_sed_inventories)
    logging_println("           change ", ocean_change_rate)
    logging_println(" orogenic sources ", ocean_orogenic_sources)
    #logging_println("           runoff ", runoff_fluxes)
    logging_println("       aolean dep ", aolean_deposition)
    logging_println("          coastal ", ocean_coastal_sources)
    logging_println("          pelagic ", ocean_pelagic_sources)
    logging_println("         phys tot ", ocean_phys_total_sources)
    logging_println(" demoted land sed ", demoted_land_sediment)
    logging_println(" uplifted ocn sed ", uplifted_ocean_sediment)
    logging_println(" denuded to coast ", denuded_fluxes_to_coast)
    logging_println("     full tot src ", ocean_full_total_sources)
    logging_println("              dep ", ocean_sed_depo_rates)
    logging_println("         mass bal ", ocean_mass_balances )
    logging_println("         flux bal ", ocean_flux_balances )

    logging_println()
    logging_println(" world budget init ", world_initial_inventories)
    #logging_println("   after tectonics ", after_tectonics_world_inventories)
    #logging_println("               bal ", ( after_tectonics_world_inventories .- initial_world_inventories ) ./
    #    main_time_step .+ subduction_rates ) #.+ crust_transition_ocean_src )

    logging_println("     source inputs ", world_source_fluxes)
    logging_println("       dissolution ", world_dissolution_fluxes)
    #logging_println("cont2ocn trans flx ", cont2ocn_redist_sources)
    #logging_println("ocn2cont trans flx ", ocn2cont_redist_sources)
    logging_println("        deposition ", world_deposition_rates)
    logging_println("    final geomorph ", new_world_inventories)
    logging_println("       change rate ", world_change_rate )
    logging_println(" geomorph mass bal ", world_geomorph_mass_balances )
    logging_println(" geomorph flux bal ", world_geomorph_flux_balances )
        
    logging_println("        subduction ", world_subduction_rates)
    logging_println("  world cum change ", world_change_rate )
    logging_println("     world cum bal ", world_change_rate .+ world_subduction_rates .- 
        world_source_fluxes .+ world_dissolution_fluxes ) # .+ ocn2cont_redist_sources ) # .+ cont2ocn_redist_sources )

end
