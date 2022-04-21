function aolean_transport()
    for ix in 1:nx
        for iy in 1:ny
            if world.freeboard[ix,iy] > 0. &&
                world.crust_type[ix,iy]  == continent_crust  # excludes orogenic
                rate = - world.freeboard[ix,iy] /
                    aolean_erosion_time_constant * time_step *
                    rho_continent_crust / rho_sediment # units of volume sediment
                set_diag("aolean_transport_deposition_rate",ix,iy,rate)
            end
        end
    end
    total_aolean_erosion = - volumefield_total(get_diag("aolean_transport_deposition_rate")) # m3
    mean_aolean_deposition_rate = total_aolean_erosion / ocean_area()
    for ix in 1:nx
        for iy in 1:ny
            if world.freeboard[ix,iy] < 0.
                set_diag("aolean_transport_deposition_rate",ix,iy,
                    mean_aolean_deposition_rate)
            end
        end
    end
    return
end

# Orogeny
function orogeny()
    orogeny_footprint = fill(0.,nx,ny)
    for (name,orogenic_event) in orogenic_events
        if world.age <= orogenic_event.onset && world.age >= orogenic_event.finish
            footprint = orogenic_event.footprint .* # 0 or 1
                orogenic_event.altitude*1e3 ./      # meters uplift
                ( orogenic_event.onset - orogenic_event.finish ) .*
                time_step
            rotatedfootprint = fill_world_orogeny(footprint)
            #smoothed = smooth_orogeny(rotatedfootprint)
            #@printf("%s %.1e " , name,orogeny_intensity(smoothed))
            orogeny_footprint .+= rotatedfootprint # smoothed
        end
    end
    orogeny_footprint .*= get_continent_mask()
    set_diag("continent_orogenic_uplift_rate",orogeny_footprint)
    subduction_orogeny = get_subduction_orogeny() # smooth_orogeny(get_subduction_orogeny())
    #for ix in 1:nx
    #    for iy in 1:ny
    #        if orogeny_footprint[ix,iy] > 0.
    #            subduction_orogeny[ix,iy] = 0.
    #        end
    #    end
    #end
    subduction_orogeny .*= get_continent_mask()
    set_diag("subduction_orogenic_uplift_rate",subduction_orogeny)
    #isostacy()
    return
end
function fill_world_orogeny(footprint) # requires at least create_world(age)
    age = world.age
    plateIDmap = world.plateID
    world_orogeny = fill(0.,nx,ny)
    for ixworld in 1:nx
        for iyworld in 1:ny
            if world.crust_type[ixworld,iyworld] >= continent_crust
                plateID = plateIDmap[ixworld,iyworld]
                mplate = resolve_rotation_matrix!(plateID,age)
                ixplate,iyplate = nearest_plateij(mplate,ixworld,iyworld)
                # the position on the rotated plate field
                # assume zero rotation at present-day, so look for footprint here
                if footprint[ixplate,iyplate] > 0.
                    world_orogeny[ixworld,iyworld] = footprint[ixplate,iyplate]
                    world.crust_type[ixworld,iyworld] = uplifted_continent_crust
                end
            end
        end
    end
    return world_orogeny
end

function continental_erosion_transport()
    exposed_crust_mask = generate_mask_field(world.sediment_thickness,0.) .*
        get_continent_mask( )
    exposed_areas = get_blobs( exposed_crust_mask )
    exposed_neighbors = get_blob_neighbor_fields( exposed_areas )
    orogenic_production_rates = get_integrated_orogenic_production_rates( exposed_areas )
    average_fringe_elevations = get_average_fringe_elevations( exposed_areas, faked_elevation_field )
    faked_elevation_field = world.freeboard .* get_continent_mask( ) .* 
        (1. .- exposed_crust_mask)
    # zeros for ocean and orogenic areas
    faked_elevation_field = setup_faked_elevation_field(faked_elevation_field, 
        exposed_areas, average_fringe_elevations) # elevate orogenic areas to avg
    elevation_perturbations = 
        get_driving_elevation_perturbations( faked_elevation_field,
            average_fringe_elevations, exposed_areas, exposed_neighbors, orogenic_production_rates )
    faked_elevation_field = world.freeboard .* get_continent_mask( ) .* 
        (1. .- exposed_crust_mask)
    faked_elevation_field = 
        setup_faked_elevation_field(faked_elevation_field, exposed_areas, 
            elevation_perturbations )
    change = diffusion_elevation_change(faked_elevation_field)
    for ix in 1:nx
        for iy in 1:ny
            if world.freeboard[ix,iy] > 0. 
                accum_diag("sed_transport_deposition_rate",ix,iy,change[ix,iy])
            else
                accum_diag("seasurface_deposition_rate",ix,iy,change[ix,iy])
            end
        end
    end 
    return
end
function continental_erosion_transport_substepped()
    exposed_crust_mask = generate_mask_field(world.sediment_thickness,0.) .*
        get_continent_mask( )
    faked_elevation_field = world.freeboard .* get_continent_mask( ) .* 
        (1. .- exposed_crust_mask)
    # zeros for ocean and orogenic areas
    exposed_areas = get_blobs( exposed_crust_mask )
    exposed_neighbors = get_blob_neighbor_fields( exposed_areas )
    orogenic_production_rates = get_integrated_orogenic_production_rates( exposed_areas )
    average_fringe_elevations = get_average_fringe_elevations( exposed_areas, faked_elevation_field )
    faked_elevation_field = setup_faked_elevation_field(faked_elevation_field, 
        exposed_areas, average_fringe_elevations) # elevate orogenic areas to avg
    elevation_perturbations = 
        get_driving_elevation_perturbations( faked_elevation_field,
            average_fringe_elevations, exposed_areas, exposed_neighbors, orogenic_production_rates )
    for i_substep in 1:n_sediment_transport_substeps
        faked_elevation_field = world.freeboard .* get_continent_mask( ) .* 
            (1. .- exposed_crust_mask)
        faked_elevation_field = 
            setup_faked_elevation_field(faked_elevation_field, exposed_areas, 
                elevation_perturbations ./ n_sediment_transport_substeps )
        change = diffusion_elevation_change(faked_elevation_field)
        for ix in 1:nx
            for iy in 1:ny
                if world.freeboard[ix,iy] > 0. 
                    if change[ix,iy] > - world.sediment_thickness[ix,iy] # enuf sediment cover
                        accum_diag("sed_transport_deposition_rate",ix,iy,change[ix,iy])
                        world.sediment_thickness[ix,iy] += change[ix,iy]
                        world.crust_type[ix,iy] = continent_crust
                    else
                        world.sediment_thickness[ix,iy] = 0.
                        world.crust_type[ix,iy] = uplifted_continent_crust
                        # but I dont want to redo all the blobs until next step
                    end
                else
                    accum_diag("seasurface_deposition_rate",ix,iy,change[ix,iy])
                end
            end
        end
        isostacy()
    end 
    return
end
function get_driving_elevation_perturbations( base_faked_elevation_field,
    average_fringe_elevations, exposed_areas, exposed_neighbors, orogenic_production_rates )
    n_areas = length( exposed_areas )
    fake_high_value = 1000.; fake_low_value = 0.
    elev_perts_low = deepcopy(average_fringe_elevations)
    elev_perts_high = elev_perts_low .+ fake_high_value
    elev_perts_interp = fill(0.,n_areas)
    orogeny_fluxes_guess = fill(0.,n_areas)
    orogeny_fluxes_low = get_orogenic_diffusive_fluxes( base_faked_elevation_field, 
        exposed_areas, exposed_neighbors, elev_perts_low )
    orogeny_fluxes_high = get_orogenic_diffusive_fluxes( base_faked_elevation_field, 
        exposed_areas, exposed_neighbors, elev_perts_high )
    for i_area in 1:n_areas
        elev_perts_interp[i_area] = average_fringe_elevations[i_area] +
            ( fake_high_value - fake_low_value ) *
            ( orogenic_production_rates[i_area] - orogeny_fluxes_low[i_area] ) /
            ( orogeny_fluxes_high[i_area] - orogeny_fluxes_low[i_area] )
    end
    return elev_perts_interp
end
function get_orogenic_diffusive_fluxes(incoming_elevation_field, orogenic_areas, 
    orogenic_neighbors, elev_perturbations)
    faked_elevation_field = 
       setup_faked_elevation_field(incoming_elevation_field, orogenic_areas, elev_perturbations)
    diffused_field = diffusion_elevation_change( faked_elevation_field )
    diffused_totals = []
    for i_area = 1:length( orogenic_areas )
        push!( diffused_totals, 
            volumefield_total( ( diffused_field .- faked_elevation_field ) .* orogenic_neighbors[i_area] ) )
    end
    return diffused_totals
end

function get_subduction_orogeny()
    orogenic_footprint = fill(0.,nx,ny)
    subduction = get_diag("ocean_subduct_plate_area")
    for ix in 1:nx
        for iy in 1:ny
            if subduction[ix,iy] > 0.
                uplift = subduction[ix,iy] * 1e-7 # 3e-8
                coordlist = get_neighbor_coords(ix,iy,2,2)
                for coord in coordlist
                    if world.crust_type[coord[1],coord[2]] >= continent_crust
                        orogenic_footprint[coord[1],coord[2]] = uplift
                        world.crust_type[coord[1],coord[2]] = uplifted_continent_crust
                    end
                end
            end
        end
    end
    return orogenic_footprint
end
function active_orogenic_events()
    orogenic_event_names = ""
    for (name,orogenic_event) in orogenic_events
        if world.age <= orogenic_event.onset && world.age >= orogenic_event.finish
            if length(orogenic_event_names) > 0
                orogenic_event_names = orogenic_event_names * "/"
            end
            orogenic_event_names = orogenic_event_names * name
        end
    end
    return orogenic_event_names
end
function create_orogenic_event(name,onset,finish,altitude)
    worldfootprint = read_flip_csv("orogenies/" * name * ".csv")
    #plateIDmap = read_plateIDs(0.)
    #footprints = parse_world_footprint(worldfootprint,plateIDmap)
    orogenic_event = orogenic_event_struct(name,worldfootprint,onset,finish,altitude)
    return orogenic_event
end
function create_orogenies()
    orogenic_base_uplift = 20.
    orogenic_events = Dict()
    orogenic_events["Pan_African"] =
        create_orogenic_event("pan_african",650.,550.,orogenic_base_uplift)  # africa, south america
    orogenic_events["Avalonian"] =
        create_orogenic_event("taconian",650.,500.,orogenic_base_uplift)
    orogenic_events["Taconian"] =
        create_orogenic_event("taconian",490.,440.,3*orogenic_base_uplift)
    orogenic_events["Calcedonian/Acadian"] =
        create_orogenic_event("calcedonian",460.,390.,orogenic_base_uplift)
    orogenic_events["Hercynian/Alleghenian/Uralian"] =
        create_orogenic_event("hercynian",300.,250.,orogenic_base_uplift) # = variscan
    # Amurian (Japan) should be covered by subduction_orogeny
    orogenic_events["Indo Sinean"] =
        create_orogenic_event("indo_sinean",200.,180.,0.5*orogenic_base_uplift)
    orogenic_events["Cimmerian"] =
        create_orogenic_event("cimmerian",180.,150.,orogenic_base_uplift)
    orogenic_events["Mongol Okhotsk"] =
        create_orogenic_event("mongol_okhotsk",140.,130.,orogenic_base_uplift)
    orogenic_events["Wrangellian"] =
        create_orogenic_event("wrangellian",110.,80.,orogenic_base_uplift)
    orogenic_events["Verkhoyansk"] =
        create_orogenic_event("verkhoyansk",100.,70.,orogenic_base_uplift)
    orogenic_events["Alpine"] =
        create_orogenic_event("alpine",50.,0.,orogenic_base_uplift)
    orogenic_events["Himalayan"] =
        create_orogenic_event("himalayan",40.,0.,orogenic_base_uplift)
    return orogenic_events
end
function smooth_orogeny(field)
    newfield = deepcopy(field)
    for ix in 1:nx
        for iy in 1:ny
            if field[ix,iy] > 0
                xys = get_neighbor_coords(ix,iy,2,2)
                for xy in xys
                    #println(xy[1]," ",xy[2])
                    if field[xy[1],xy[2]] == 0
                        newfield[xy[1],xy[2]] +=
                            (field[ix,iy] - newfield[xy[1],xy[2]]) * 0.2
                    end
                end
            end
        end
    end
    return newfield
end
function orogeny_intensity(uplift)
    total = 0.
    for ix in 1:nx
        for iy in 1:ny
            if uplift[ix,iy] > 0.
                total += uplift[ix,iy] * areabox[iy]
            end
        end
    end
    return total
end
function orogenic_erosion_rate( freeboard )
    initial_tau = 250; longterm_tau = 250.; cutoff_elevation = 2000.
    if freeboard > cutoff_elevation
        erosion_rate = max(2000.,freeboard) / initial_tau * time_step
    else
        erosion_rate = max(2000.,freeboard) / longterm_tau * time_step
    end
    return erosion_rate # units of volume of continental crust
end        
function get_integrated_orogenic_production_rates( orogenic_areas )
    n_orogenies = length(orogenic_areas)
    area_rates = []
    for i_area = 1:n_orogenies
        erosion_field = fill(0.,nx,ny)
        for ix in 1:nx
            for iy in 1:ny
                if orogenic_areas[i_area][ix,iy] == 1
                    rate = orogenic_erosion_rate(world.freeboard[ix,iy]) 
                    erosion_field[ix,iy] = rate # units of volume crust
                    set_diag("crust_weathering_rate",ix,iy,rate )
                end
            end
        end
        push!( area_rates, volumefield_total( erosion_field ) * 
            rho_continent_crust / rho_sediment ) # units volume of sediment
    end
    return area_rates
end
function get_average_fringe_elevations( orogenic_areas, faked_elevation_field )
    average_fringe_elevations = []
    for orogenic_area in orogenic_areas
        neighbor_field = get_blob_neighbors( orogenic_area )
        average_fringe_elevation = 
            volumefield_total( faked_elevation_field .* neighbor_field ) / 
            volumefield_total( neighbor_field ) 
        push!( average_fringe_elevations, average_fringe_elevation )
    end
    return average_fringe_elevations
end    
function setup_faked_elevation_field(incoming_elevation_field, orogenic_areas, elev_perturbations)
    n_areas = length(orogenic_areas)
    perturbed_field = deepcopy(incoming_elevation_field)
    for i_area in 1:n_areas
        perturbed_field .+= orogenic_areas[i_area] .* elev_perturbations[i_area]
    end
    return perturbed_field
end
# Diffusive sediment transport
function diffusion_elevation_change(elev_field)
    new_elev_field = reshape( LUofAglobal \ vec(elev_field), (nx,ny) )
    return new_elev_field .- elev_field
end
function setup_2d_transport()
    d = sediment_transport_coefficient * time_step
    n = nx * ny
    A = spzeros(n,n)
    for ix in 1:nx
        for iy in 1:ny # each grid point
            A[i2dmat(ix,iy),i2dmat(ix,iy)] = 1. + 2. * d # itself, horizontal
            # box to the left (or wrap)
            if ix > 1
                A[i2dmat(ix-1,iy),i2dmat(ix,iy)] = -d
                A[i2dmat(ix,iy),i2dmat(ix-1,iy)] = -d
            else
                A[i2dmat(nx,iy),i2dmat(ix,iy)] = -d
                A[i2dmat(ix,iy),i2dmat(nx,iy)] = -d
            end
            # box below
            if iy == 1
                A[i2dmat(ix,iy),i2dmat(ix,iy)] += d
            elseif iy == ny
                A[i2dmat(ix,iy),i2dmat(ix,iy)] += d
            else
                A[i2dmat(ix,iy),i2dmat(ix,iy)] += 2. * d
            end
            if iy > 1
                A[i2dmat(ix,iy-1),i2dmat(ix,iy)] = -d #* areabox[iy] / areabox[iy-1]
                A[i2dmat(ix,iy),i2dmat(ix,iy-1)] = -d #* areabox[iy] / areabox[iy-1]
            end
        end
    end
    return lu(A)
end
function i2dmat(ix,iy)
    return ix + (iy-1) * nx
end
function mat2ij(id)
    iy = Int(floor(id/nx))
    ix = id - iy
    return ix, iy
end

