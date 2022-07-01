function orogeny()
    orogeny_footprint = fill(0.,nx,ny)
    for (name,orogenic_event) in orogenic_events
        if world.age <= orogenic_event.onset && world.age >= orogenic_event.finish
            footprint = orogenic_event.footprint .* # 0 or 1
                orogenic_event.altitude*1e3 ./      # meters uplift
                ( orogenic_event.onset - orogenic_event.finish ) .*
                time_step
            rotatedfootprint = fill_world_orogeny(footprint)
            smoothed = smooth_orogeny(rotatedfootprint)
            #@printf("%s %.1e " , name,orogeny_intensity(smoothed))
            orogeny_footprint .+= smoothed
        end
    end
    orogeny_footprint .*= get_continent_mask()
    set_diag("continent_orogenic_uplift_rate",orogeny_footprint)
    subduction_orogeny = smooth_orogeny(get_subduction_orogeny())
    for ix in 1:nx
        for iy in 1:ny
            if orogeny_footprint[ix,iy] > 0.
                subduction_orogeny[ix,iy] = 0.
            end
        end
    end
    subduction_orogeny .*= get_continent_mask()
    set_diag("subduction_orogenic_uplift_rate",subduction_orogeny)
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
                end
            end
        end
    end
    return world_orogeny
end
function get_subduction_orogeny()
    orogenic_footprint = fill(0.,nx,ny)
    subduction = get_diag("ocean_subduct_plate_area")
    for ix in 1:nx
        for iy in 1:ny
            if subduction[ix,iy] > 0.
                uplift = subduction[ix,iy] * areabox[iy] * time_step * subduction_orogeny_parameter
                coordlist = get_neighbor_coords(ix,iy,2,2)
                for coord in coordlist
                    if world.crust_type[coord[1],coord[2]] >= continent_crust
                        orogenic_footprint[coord[1],coord[2]] = uplift
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
    orogenic_base_uplift_rate = 7.5
    orogenic_events = Dict()
    orogenic_events["Pan_African"] =
        create_orogenic_event("pan_african",650.,550.,orogenic_base_uplift_rate)  # africa, south america
    orogenic_events["Avalonian"] =
        create_orogenic_event("taconian",650.,500.,orogenic_base_uplift_rate)
    orogenic_events["Taconian"] =
        create_orogenic_event("taconian",490.,440.,3*orogenic_base_uplift_rate)
    orogenic_events["Calcedonian/Acadian"] =
        create_orogenic_event("calcedonian",460.,390.,orogenic_base_uplift_rate)
    orogenic_events["Hercynian/Alleghenian/Uralian"] =
        create_orogenic_event("hercynian",300.,250.,orogenic_base_uplift_rate) # = variscan
    # Amurian (Japan) should be covered by subduction_orogeny
    orogenic_events["Indo Sinean"] =
        create_orogenic_event("indo_sinean",200.,180.,0.5*orogenic_base_uplift_rate)
    orogenic_events["Cimmerian"] =
        create_orogenic_event("cimmerian",180.,150.,orogenic_base_uplift_rate)
    orogenic_events["Mongol Okhotsk"] =
        create_orogenic_event("mongol_okhotsk",140.,130.,orogenic_base_uplift_rate)
    orogenic_events["Wrangellian"] =
        create_orogenic_event("wrangellian",110.,80.,orogenic_base_uplift_rate)
    orogenic_events["Verkhoyansk"] =
        create_orogenic_event("verkhoyansk",100.,70.,orogenic_base_uplift_rate)
    orogenic_events["Alpine"] =
        create_orogenic_event("alpine",50.,0.,orogenic_base_uplift_rate)
    orogenic_events["Himalayan"] =
        create_orogenic_event("himalayan",40.,0.,orogenic_base_uplift_rate)
    return orogenic_events
end
function smooth_orogeny(field)
    newfield = deepcopy(field)
    for ix in 1:nx
        for iy in 1:ny
            if field[ix,iy] > 0
                xys = get_neighbor_coords(ix,iy,2,2)
                for xy in xys
                    ixn = xy[1]; iyn = xy[2]
                    if field[ix,iy] > newfield[ixn,iyn]
                        difference = field[ix,iy] - newfield[ixn,iyn]
                        newfield[ixn,iyn] += difference * 0.2
                    end
                end
            end
        end
    end
    return newfield
end
function land_sediment_fraction_transport( new_elevation_field, 
    subaereal_mask, 
    sediment_sources, ocean_sink )
    # sets land_sediment_fraction_deposition_rate on land areas
    subaereal_blobs = get_blobs( subaereal_mask )
 
    #Threads.@threads 
    for i_blob in 1:length(subaereal_blobs) 
        #println("beginning ",i_blob)
        subaereal_blob = subaereal_blobs[i_blob]
        land_area_sediment_fraction_transport( 
            new_elevation_field, 
            subaereal_blob, sediment_sources, ocean_sink ) 
        # sets land_sediment_fraction_deposition_rate on land points
        #println("i_blob ",i_blob," ",
        #    volume_field(get_frac_diag("land_sediment_fraction_deposition_rate")[:,:,CaCO3_sediment,1]))
    end 
end
function smooth_continental_crust_thickness( )
    crust_mask = eq_mask(world.crust_type,continent_crust)
    crust_blobs = get_blobs( crust_mask )
    for blob in crust_blobs
        smooth_continental_crust_blob( blob )
    end
end
function smooth_continental_crust_blob( blob )
    if get_maskfield_census( blob ) > 1
        list_pos = fill(0,nx,ny)
        list_pos_filled = 0
        list_x = []; list_y = []
        for ix in 1:nx
            for iy in 1:ny
                if blob[ix,iy] == 1
                    list_pos_filled += 1
                    list_pos[ix,iy] = list_pos_filled
                    push!( list_x, ix ); push!( list_y, iy )
                end
            end
        end
        A = spzeros( list_pos_filled, list_pos_filled )
        for i_pos = 1:list_pos_filled
            ix = list_x[i_pos]; iy = list_y[i_pos]
            h_diffcoeffs = get_horizontal_diffcoeff(crust_smooth_coeff,iy)
            v_diffcoeffs = get_vertical_diffcoeff(crust_smooth_coeff,iy)
            A[i_pos,i_pos] = 1.
            if blob[ix_left(ix),iy] == 1 
                A[i_pos,i_pos] += h_diffcoeffs[2] # left = -ve 
                i_neighbor_pos = list_pos[ix_left(ix),iy]
                A[i_pos,i_neighbor_pos] = - h_diffcoeffs[2] # only fill one direction
            end
            if blob[ix_right(ix),iy] == 1 
                A[i_pos,i_pos] += h_diffcoeffs[1]
                i_neighbor_pos = list_pos[ix_right(ix),iy]
                A[i_pos,i_neighbor_pos] = - h_diffcoeffs[1]
            end
            if iy == 1
                ixn = ( ix + 179 ) % 360 + 1
                if blob[ixn,iy] == 1
                    A[i_pos,i_pos] += v_diffcoeffs[2]
                    i_neighbor_pos = list_pos[ixn,iy]
                    A[i_pos,i_neighbor_pos] = - v_diffcoeffs[2]
                end
            else # iy > 1
                if blob[ix,iy-1] == 1
                    A[i_pos,i_pos] += v_diffcoeffs[2]
                    i_neighbor_pos = list_pos[ix,iy-1]
                    A[i_pos,i_neighbor_pos] = - v_diffcoeffs[2]
                end
            end
            if iy == ny
                ixn = ( ix + 179 ) % 360 + 1
                if blob[ixn,iy] == 1
                    A[i_pos,i_pos] += v_diffcoeffs[1]
                    i_neighbor_pos = list_pos[ixn,iy]
                    A[i_pos,i_neighbor_pos] = - v_diffcoeffs[1]
                end
            else # iy < ny
                if blob[ix,iy+1] == 1
                    A[i_pos,i_pos] += v_diffcoeffs[1]
                    i_neighbor_pos = list_pos[ix,iy+1]
                    A[i_pos,i_neighbor_pos] = - v_diffcoeffs[1]
                end
            end
        end
        lu_A = lu(A)
        R = Float64[]; 
        for i_pos = 1:list_pos_filled
            ix = list_x[i_pos]; iy = list_y[i_pos]
            push!(R, world.crust_thickness[ix,iy])
        end
        new_values_list = lu_A \ R
        for i_pos = 1:list_pos_filled
            ix = list_x[i_pos]; iy = list_y[i_pos]
            world.crust_thickness[ix,iy] = new_values_list[i_pos]
        end
    end # diffuse or no
end
#function smooth_orogenic_topography()
function generate_orogenic_erosion_fluxes()
    # resets then accumulates land_orogenic_clay_flux and coastal_orogenic_clay_flux
    # also sets cruse_erosion_rate,
    # doesnt need to accumulate since there will be no overlap of orogenic areas.
    # uses exposed_basement
    verbose = false
    reset_diag("land_orogenic_clay_flux")
    reset_diag("coastal_orogenic_clay_flux")
    orogenic_crust_mask = eq_mask( world.geomorphology,exposed_basement ) .*
        eq_mask( world.crust_type,continent_crust )
    orogenic_blobs = get_blobs( orogenic_crust_mask )
    orogenic_blob_areas = get_blob_areas( orogenic_blobs )
    orogenic_production_rates = get_integrated_orogenic_production_rates( orogenic_blobs )
    # fills crust_erosion_rate and crust_clay_source_rate, 
    # returns blob-integrated rates in an array[nblobs]
    neighbor_blobs = get_blob_neighbor_fields( orogenic_blobs )
    neighbor_blob_areas = get_blob_areas( neighbor_blobs )
    for i_area in 1:length( orogenic_blobs )
        for ix in 1:nx
            for iy in 1:ny
                if neighbor_blobs[i_area][ix,iy] > 0
                    massflux = orogenic_production_rates[i_area] / # m3 / Myr
                        neighbor_blob_areas[i_area]  # meters / Myr
                    if world.geomorphology[ix,iy] == sedimented_land
                        accum_diag("land_orogenic_clay_flux",ix,iy,massflux)
                    elseif world.geomorphology[ix,iy] < sedimented_land # ocean
                        accum_diag("coastal_orogenic_clay_flux",ix,iy,massflux)
                    end
                end
            end
        end
    end
    if verbose == true
        orogenic_area = volume_field(orogenic_crust_mask)
        avg_oro_elev = field_mean( world.freeboard .* orogenic_crust_mask )
        crust_clay_source = volume_field(get_diag("land_orogenic_clay_flux")) +
            volume_field(get_diag("coastal_orogenic_clay_flux"))
        println(" area ",orogenic_area," avg elevation ",avg_oro_elev,
            " clay src ", crust_clay_source)
    end
    return
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
function orogenic_erosion_rate( freeboard ) # meters / Myr
#= 
    initial_tau = 100; longterm_tau = 250.; cutoff_elevation = 2000.
    if freeboard > cutoff_elevation
        erosion_rate = max(2000.,freeboard) / initial_tau 
    else
        erosion_rate = freeboard / longterm_tau
    end
=#
    #tau_apparent = 30.  # Myr
    tau = orogenic_erosion_tau_apparent * crust_freeboard_expression
    erosion_rate = freeboard / tau
    return erosion_rate # meters of crust / Myr
end        
function get_integrated_orogenic_production_rates( orogenic_blobs ) # m3 sediment / Myr
    # fills crust_erosion_rate and crust_clay_source_rate, 
    # returns blob-integrated rates in an array[nblobs]
    n_orogenies = length(orogenic_blobs)
    total_sediment_thickness = world.sediment_thickness
    area_rates = []
    for i_area = 1:n_orogenies
        crust_erosion_field = fill(0.,nx,ny)
        for ix in 1:nx
            for iy in 1:ny
                if orogenic_blobs[i_area][ix,iy] == 1
                    crust_erosion_rate = orogenic_erosion_rate(world.freeboard[ix,iy]) # meters / Myr
                    crust_erosion_field[ix,iy] = crust_erosion_rate
                    #=
                    if total_sediment_thickness[ix,iy] > 0. # presumably recently exposed, leftover
                        # pays off negative sediment debt if any
                        #world.tectonics[ix,iy] = land_soil_denuded_with_excess
                        sediment_source_rate = total_sediment_thickness[ix,iy] / 
                            time_step # meters / Myr
                        world.sediment_thickness[ix,iy] = 0.
                        world.sediment_layer_thickness[ix,iy,:] .= 0.
                        # turned this off because its already getting there thru dispersion of erosion
                        #= for i_sedfrac in 1:n_sediment_types
                            fraction = world.sediment_fractions[ix,iy,i_sedfrac]
                            set_frac_diag("land_sediment_fraction_deposition_rate",ix,iy,i_sedfrac,
                                - sediment_source_rate * fraction )
                        end =#
                        # get rid of leftover sediment through proper channels
                        crust_erosion_field[ix,iy] += sediment_source_rate * 
                            rho_sediment / rho_continent_crust # meters crust / Myr
                        # add it to crust_erosion which will be converted back to sediment
                        # and added to the orogenic source
                    end
                    =#
                    set_diag("crust_erosion_rate",ix,iy,crust_erosion_rate ) # meters crust / Myr
                    set_diag("crust_clay_source_rate",ix,iy,crust_erosion_rate *
                        rho_continent_crust / rho_sediment )
                end
            end
        end
        total_crust_erosion = volume_field( crust_erosion_field )
        total_sediment_source = total_crust_erosion *
           rho_continent_crust / rho_sediment
        push!( area_rates, total_sediment_source ) # m3 sediment / Myr
    end
    return area_rates
end
function get_blob_areas( blobs )
    areas = []
    for blob in blobs
        area = volume_field( blob )
        push!( areas, area )
    end
    return areas
end
