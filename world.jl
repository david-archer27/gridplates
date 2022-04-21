
# Meta steps
function step_everything() # rebuilds the world at the new time
    world.age -= time_step
    clear_world_process_arrays() # battle stations
    #newIDs = read_plateIDs(newtime)
    resolve_rotation_matrices()
    # set rotation matrices to the new time
    increment_plate_age()
    # ocean and cont crust gets older in plate grids
    ocean_thermal_boundary_layer()
    # sets world.elevation_offset for aging ocean crust
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
    # initializes plate crust_* variables at new points in plate fields
    # new continent points on plate grid have h0 thickness

    orogeny()
    # continent and subduction uplift rates into *_orogenic_uplift_rate diags

    apply_orogeny_fluxes_to_world()
    # updates world.crust_thickness and demotes orogenic regions if buried
    # returns world in isostatic equilibrium

    for i_step in 1:n_sediment_transport_substeps
        clear_geomorph_process_arrays()
        continental_erosion_transport()
        # diffuses a sealevel-clamped elevation field to generate fields of crust
        # erosion, sediment deposition on land and seasurface deposition on ocean
        # points.
        aolean_transport()
        # erodes continental soils at an increasing aolean_erosion_rate with elevation.
        # deposits into globally uniform seasurface_deposition_ratedeposition.
        # operates on world diag fields only
        seafloor_deposition()
        # translate seasurface deposition into sediment_deposition_rate entries for
        # the ocean grid points.
        # recursive lateral downhill transport when accomodation space is lacking.

        apply_sediment_transport_fluxes_to_world()
        # update world sediment_thickness and crust_thickness, just for show if looping
    end

    apply_crust_changes_to_plates()
    apply_sediment_fluxes_to_plates()
    # updates the plate crust_thickness and sediment_thickness fields
    isostacy()
    # just for show if looping
    # clock out
    return
end
function step_geomorph() # requires previous run with tectonics to fill world diag arrays
    world.age -= time_step
    clear_geomorph_process_arrays()
    increment_world_age() # world. files only
    ocean_thermal_boundary_layer()
    apply_orogeny_fluxes_to_world() # from world diags set by previous tectonics run
    continental_erosion_transport() 
    # resets sed_transport_deposition_rate, seasurface_deposition_rate, crust_weathering_rate
    aolean_transport() # resets aolean_transport_deposition_rate
    seafloor_deposition() # resets seafloor_deposition_rate, seafloor_scour_fraction,seafloor_spilloff_rate
    apply_sediment_transport_fluxes_to_world()
    isostacy( )
    return
end

