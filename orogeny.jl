function get_scotese_mountains( )
    age_string = lpad(Int(ceil(world.age)),3,"0")
    CS_file_name = code_base_directory * 
        "/data/scotese_elevation_files_today/scotese_today." *
        age_string * "Ma.bson"
    BSON.@load CS_file_name twisted_elevation
    target_elevation = fill_world_orogeny( twisted_elevation )
    
    #target_elevation = rotate_field( twisted_elevation, 
    #    eq_mask(world.crust_type,continent_crust), cont_ID, cont_start, world.age )
    #smooth_world!( target_elevation, 1000. )
    return target_elevation
end
function rotate_field_0_to_current_age(field_at_time_0, start_time_field )
    #continent_crust_past = eq_mask(world.crust_type,continent_crust)
    field_past = rotate_field_0_to_age( field_at_time_0, start_time_field,
        world.age, world.continentID )
    return field_past
end
function mask_field_from_start_time_field( start_time_field )
    mask_field_time_0 = fill(1.,nx,ny)
    field_past = rotate_field_0_to_current_age(mask_field_time_0, start_time_field)
    return field_past
end

function rotate_field_0_to_age( field_at_time_0, start_time_field, age,
    paleo_rotation_ID_field )
    # inner core unreliant on world.machinery 
    field_past = fill(0.,nx,ny)
    for ixpast in 1:nx
        for iypast in 1:ny
            if paleo_rotation_ID_field[ixpast,iypast] > 0
                rot_ID = paleo_rotation_ID_field[ixpast,iypast]
                rotation_matrix = resolve_rotation_matrix!( rot_ID, age )
                ix0, iy0 = nearest_plateij(rotation_matrix, ixpast, iypast)
                # nearest plate because that is the 0-rotation location 
                if start_time_field[ix0,iy0] > age 
                    field_past[ixpast, iypast] = field_at_time_0[ix0, iy0]
                end
            end
        end
    end
    #=rotIDlist = find_plateID_list( rot_ids_0 )
    for rotID in rotIDlist
        if rotID in rot_ids_list_past
            #error(rotID)
            
            for ixpast in 1:nx # need the world grid at time=age to be complete
                for iypast in 1:ny
                    ix0,iy0 = nearest_plateij(rotation_matrix,ixpast,iypast)
                    # the position on the rotated plate field
                    # assume zero rotation at time=0, so look for footprint here
                    if rot_ids_0[ix0,iy0] == rotID && 
                        continent_crust_past[ixpast,iypast] == 1 &&
                        start_time_0[ix0,iy0] >= age

                        #error(ixpast," ",iypast)
                        field_past[ixpast,iypast] = field_0[ix0,iy0]
                    end
                end
            end
        end
    end =#
    return field_past
end
function rotate_field_age_to_0(field_at_age, start_time_field, age)
    rotation_ID_time_0 = read_continentIDs(0.0)
    # inner core unreliant on world.machinery 
    field_at_time_0 = fill(0.0, nx, ny)
    for ix0 in 1:nx
        for iy0 in 1:ny
            if rotation_ID_time_0[ix0, iy0] > 0
                rot_ID = rotation_ID_time_0[ix0, iy0]
                rotation_matrix = resolve_rotation_matrix!(rot_ID, age)
                ixpast, iypast = nearest_worldij(rotation_matrix, ix0, iy0)
                # nearest plate because that is the 0-rotation location 
                if start_time_field[ix0, iy0] > age
                    field_at_time_0[ix0, iy0] = field_past[ixpast, iypast]
                end
            end
        end
    end
    return field_at_time_0
end
function filtered_start_time_field( start_time_field )
    filtered_start_time_field = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if start_time_field[ix,iy] - world.age > 0 &&
                start_time_field[ix,iy] - world.age <= main_time_step

                filtered_start_time_field[ix,iy] = start_time_field[ix,iy]
            end
        end
    end
    return filtered_start_time_field
end
function get_mafics()
    filtered_large_igneous_province_start_field =
        filtered_start_time_field(large_igneous_province_start_field)
    current_large_igneous_provinces =
        mask_field_from_start_time_field(filtered_large_igneous_province_start_field)
    filtered_ophiolite_start_field =
        filtered_start_time_field(large_igneous_province_start_field)
    current_ophiolites =
        mask_field_from_start_time_field(filtered_ophiolite_start_field)
    for ix in 1:nx
        for iy in 1:ny
            if current_large_igneous_provinces[ix, iy] == 1.0 ||
                current_ophiolites[ix, iy] == 1.0

                world.crust_composition[ix, iy] = mafic_crust
                world.crust_thickness[ix, iy] +=
                    (1500.0 - world.freeboard[ix, iy]) /
                    crust_freeboard_expression
                world.geomorphology[ix, iy] = exposed_basement
                # triggers cleanup of sediment in get_denuded_sediment_fluxes()
            end
        end
    end
end
#function get_ice_sheets()

function uplift_scotese_mountains_by_crust_thickening( scotese_target_elevation )
    for ix in 1:nx
        for iy in 1:ny
            if scotese_target_elevation[ix, iy] > 1500.0
                world.crust_thickness[ix,iy] += 
                    ( scotese_target_elevation[ix, iy] - world.freeboard[ix, iy] ) /
                    crust_freeboard_expression
            end
        end
    end
end
function uplift_scotese_mountains_by_magic( scotese_target_elevation )
    for ix in 1:nx
        for iy in 1:ny
            if scotese_target_elevation[ix,iy] > 1500.
                world.elevation_offset[ix,iy] += scotese_target_elevation[ix,iy] -
                    world.freeboard[ix,iy]
            end
        end
    end
end
              
                




function old_orogeny( subduction_footprint )
    orogenic_uplift_rates = fill(0.,nx,ny,0:2) # sum, cont, subd
    # unmasked for land
    world_mask = fill(1,nx,ny)
    for (name,orogenic_event) in orogenic_events
        if world.age <= orogenic_event.onset && world.age >= orogenic_event.finish
            println( orogenic_event.name  )
            footprint = orogenic_event.footprint .* # 0 or 1
                orogenic_uplift_parameter .* 
                #orogenic_base_uplift_meters .*
                orogenic_event.uplift_scale #./ # meters of crust over the event
                #( orogenic_event.onset - orogenic_event.finish ) # Myr
            rotatedfootprint = fill_world_orogeny(footprint)
            #smooth_world!( rotatedfootprint, orogeny_smooth_coeff )
            #smooth_masked_field!(rotatedfootprint, world_mask, orogeny_smooth_coeff)
            #@printf("%s %.1e " , name,orogeny_intensity(smoothed))
            orogenic_uplift_rates[:,:,1] .+= rotatedfootprint
        end
    end
    #orogeny_footprint .*= get_continent_mask()
    #set_diag("continent_orogenic_uplift_rate",orogeny_footprint)
    #subduction_orogeny = fill(0.,nx,ny)
    #if enable_subduction_orogeny
    #    subduction_orogeny = get_subduction_orogeny( subduction_footprint )
    #end
    #smooth_world!(subduction_orogeny,subduction_orogeny_smooth_coeff)
    for ix in 1:nx
        for iy in 1:ny
            if world.crust_type[ix,iy] == continent_crust 
                orogenic_uplift_rates[ix,iy,1] -= ( 
                    world.crust_thickness[ix,iy] - continent_crust_h0 ) /
                    orogenic_erosion_tau_apparent
            end
        end
    end
    #orogenic_uplift_rates[:,:,2] .= subduction_orogeny
    orogenic_uplift_rates[:,:,0] = orogenic_uplift_rates[:,:,1] .+
        orogenic_uplift_rates[:,:,2]
    smooth_world_fields!( orogenic_uplift_rates, orogeny_smooth_coeff )
    #subduction_orogeny .*= get_continent_mask()
    #set_diag("subduction_orogenic_uplift_rate",subduction_orogeny)
    return orogenic_uplift_rates
end
function fill_world_orogeny(footprint) # requires at least create_world(age)
    age = world.age
    continentIDmap = world.continentID
    world_orogeny = fill(0.,nx,ny)
    for ixworld in 1:nx
        for iyworld in 1:ny
            if world.crust_type[ixworld,iyworld] >= continent_crust
                continentID = continentIDmap[ixworld,iyworld]
                mplate = resolve_rotation_matrix!(continentID,age)
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

function get_subduction_orogeny( subduction_footprint )
    #orogenic_footprint = fill(0.,nx,ny)
    crust_convergence = deepcopy(subduction_footprint) #.* # meters (vol/area) / main_time_step
        #subduction_orogeny_parameter
    #new_crust_convergence = fill(0.,nx,ny)
    is_continent = eq_mask(world.crust_type,continent_crust)
    fat_continent = deepcopy( is_continent )
    for i_skin in 1:5 # expand fat_continent
        neighbor_cell_field = get_blob_neighbors( fat_continent )
        fat_continent .+= neighbor_cell_field
    end
    for i_skin in 1:6 # capture everything and concentrate it on continent
        #println(i_skin)
        fringe = get_blob_fringe(fat_continent)
        fat_continent .-= fringe
        #=
        println("before ", volume_field(crust_convergence .* fringe),
            " ", volume_field(crust_convergence .* fat_continent))
        =#
        for ix in 1:nx
            for iy in 1:ny
                if fringe[ix,iy] == 1 && crust_convergence[ix,iy] > 0
                    #error(ix,iy)
                    distribute_to_nearby!(crust_convergence,ix,iy,fat_continent)
                end
            end
        end
        #=
        println("after ", volume_field(crust_convergence .* fringe),
            " ", volume_field(crust_convergence .* fat_continent))
        println("cont ", volume_field(crust_convergence .* fringe),
            " ", volume_field(crust_convergence .* is_continent))
        =#
    end

    #println(volume_field(crust_convergence))       
    for plateID in world.plateIDlist
        plate_footprint = eq_mask(world.plateID,plateID)
        #if get_maskfield_census( plate_footprint ) > 
        #println(plateID)
        smooth_masked_field!( crust_convergence, plate_footprint,
            subduction_orogeny_smooth_coeff)
    end
    #smooth_world!( crust_convergence, subduction_orogeny_smooth_coeff)
    for plateID in world.plateIDlist
        #plateID = 201
        #println("plateID ",plateID)
        original_plate_blob = eq_mask(world.plateID,plateID)
        plate_blob = deepcopy( original_plate_blob )
        expanded_plate_blob = deepcopy( original_plate_blob )
        n_onion_skins = Int( floor( subduction_orogeny_hor_offset / delta_y ) )
        for i_skin in 1:n_onion_skins # major shell 
            neighbor_cell_field = get_blob_neighbors( plate_blob )
            expanded_plate_blob = plate_blob .+ neighbor_cell_field
            #println(volume_field(crust_convergence))
            for ix in 1:nx
                for iy in 1:ny
                    if neighbor_cell_field[ix,iy] == 1 && crust_convergence[ix,iy] > 0.
                        distribute_to_nearby!(crust_convergence,ix,iy,
                            1 .- expanded_plate_blob)
                    end
                end
            end
            #crust_convergence .+= new_crust_convergence
            plate_blob = deepcopy(expanded_plate_blob)
            #new_crust_convergence = fill(0.,nx,ny)
            #println("skin ", volume_field(crust_convergence))  
            #=
            n_extra_cells_max = 5
            for i_sub_cell_move in 1:n_extra_cells_max
                #i_sub_cell_move = 1            
                for iy in 1:ny # fill in horizontally in higher lat
                    if n_equiv_boxes[iy] >= i_sub_cell_move
                        for ix in 1:nx
                            if neighbor_cell_field[ix,iy] == 1 && crust_convergence[ix,iy] > 0.
                                coordpairs = []
                                if plate_blob[ix_left(ix),iy] == 0
                                    push!( coordpairs, [ix_left(ix),iy] )
                                end
                                if plate_blob[ix_right(ix),iy] == 0
                                    push!( coordpairs, [ix_right(ix),iy] )
                                end
                                
                                if length(coordpairs) > 0
                                    area_neighbors = 0.
                                    #error(ix," ",iy)
                                    for coordpair in coordpairs
                                        ixn = coordpair[1]; iyn = coordpair[2]
                                        area_neighbors += areabox[iyn]
                                    end
                                    for coordpair in coordpairs
                                        ixn = coordpair[1]; iyn = coordpair[2]
                                        new_crust_convergence[ixn,iyn] += crust_convergence[ix,iy] * 
                                            areabox[iy] / area_neighbors
                                        expanded_plate_blob[ixn,iyn] = 1
                                    end
                                    crust_convergence[ix,iy] = 0.
                                else
                                    new_crust_convergence[ix,iy] = crust_convergence[ix,iy]
                                    crust_convergence[ix,iy] = 0.
                                end
                            end
                        end
                    end
                end
                crust_convergence .+= new_crust_convergence
                plate_blob = deepcopy(expanded_plate_blob)
                new_crust_convergence = fill(0.,nx,ny)
                #println(" sub ", volume_field(crust_convergence))  
            end 
            #crust_convergence .+= new_crust_convergence
            =#
        end
        #println(volume_field(crust_convergence))
    end
    smooth_world!(crust_convergence,subduction_orogeny_smooth_coeff)
    crust_convergence .*= subduction_orogeny_parameter
    return crust_convergence
end
function distribute_to_nearby!(value_field,ix,iy,ok_field)
    coordpairs = get_border_neighbor_coords(ix,iy)
    filtered_coordpairs = []
    for coordpair in coordpairs
        ixn = coordpair[1]; iyn = coordpair[2]
        if ok_field[ixn,iyn] == 1
            push!(filtered_coordpairs, coordpair)
        end
    end
    if length(filtered_coordpairs) > 0
        area_neighbors = 0.
        for coordpair in filtered_coordpairs
            ixn = coordpair[1]; iyn = coordpair[2]
            area_neighbors += areabox[iyn]
        end
        for coordpair in filtered_coordpairs
            ixn = coordpair[1]; iyn = coordpair[2]
            value_field[ixn,iyn] += value_field[ix,iy] * 
                areabox[iy] / area_neighbors
            #new_plate_blob[ixn,iyn] = 1
        end
        value_field[ix,iy] = 0.
    #else
    #    new_crust_convergence[ix,iy] = crust_convergence[ix,iy]
    #    crust_convergence[ix,iy] = 0.
        # just for conservation; itll be ocean probably
    end
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
function create_orogenic_event(name,onset,finish,uplift_rate_scale)
    worldfootprint = read_flip_csv("orogenies/" * name * ".csv")
    #plateIDmap = read_plateIDs(0.)
    #footprints = parse_world_footprint(worldfootprint,plateIDmap)
    orogenic_event = orogenic_event_struct(name,worldfootprint,onset,finish,uplift_rate_scale)
    return orogenic_event
end

function smooth_masked_field!( values, mask, diffcoeff )
    # alters sed_fluxes
    blobs = get_blobs( mask )
    Threads.@threads for i_blob in 1:length(blobs)
        #println("starting ",i_blob)
        blob = blobs[i_blob]
        smooth_area!( values, blob, diffcoeff )
    end
end
function smooth_area!( values, blob, diffcoeff )
    # alters sed_fluxes
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
                h_diffcoeffs = get_horizontal_diffcoeff(diffcoeff,iy)
                v_diffcoeffs = get_vertical_diffcoeff(diffcoeff,iy)
            #push!(original_values_list, in_field[ix,iy])
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

        #for i_sedtype in 0:n_sediment_types
            #input_fraction_fluxes = 
            #    get_frac_diag("seafloor_sediment_fraction_deposition_rate",i_sedtype)
        R = Float64[]; 
        for i_pos = 1:list_pos_filled
            ix = list_x[i_pos]; iy = list_y[i_pos]
            push!(R, values[ix,iy])
        end
        new_values_list = lu_A \ R
        for i_pos = 1:list_pos_filled
            ix = list_x[i_pos]; iy = list_y[i_pos]
            values[ix,iy] = new_values_list[i_pos]
            #set_frac_diag("seafloor_sediment_fraction_deposition_rate",
            #    ix,iy,i_sedtype,new_values_list[i_pos])
            #accum_diag("seafloor_sediment_deposition_rate",
            #    ix,iy,new_values_list[i_pos])
        end
        #end
    end # diffuse or no
end
function smooth_world!( values, diffcoeff )
    # used for orogeny footprint
    n_columns = nx * ny
    A = spzeros( n_columns, n_columns )
    for ix in 1:nx
        for iy in 1:ny    
            i_row = ix + ( iy - 1 ) * nx
            h_diffcoeffs = get_horizontal_diffcoeff(diffcoeff,iy)
            v_diffcoeffs = get_vertical_diffcoeff(diffcoeff,iy)

            A[i_row,i_row] = 1.
            # left
            A[i_row,i_row] += h_diffcoeffs[2] # left = -ve 
            i_neighbor_row = ix_left(ix) + ( iy - 1 ) * nx
            A[i_row,i_neighbor_row] = - h_diffcoeffs[2] # only fill one direction
            # right
            A[i_row,i_row] += h_diffcoeffs[1]
            i_neighbor_row = ix_right(ix) + ( iy - 1 ) * nx
            A[i_row,i_neighbor_row] = - h_diffcoeffs[1]
            if iy > 1
                A[i_row,i_row] += v_diffcoeffs[2]
                i_neighbor_row = ix + ( iy - 2 ) * nx
                A[i_row,i_neighbor_row] = - v_diffcoeffs[2]
            end
            if iy < ny
                A[i_row,i_row] += v_diffcoeffs[1]
                i_neighbor_row = ix + ( iy ) * nx
                A[i_row,i_neighbor_row] = - v_diffcoeffs[1]
            end
        end
    end
    
    lu_A = lu(A)

    R = fill(0.,n_columns)# Float64[]; 
    for ix in 1:nx
        for iy in 1:ny
            i_row = ix + ( iy - 1 ) * nx
            #push!(R, values[ix,iy])
            R[i_row] = values[ix,iy]
        end
    end
    new_values_list = lu_A \ R
    for ix in 1:nx
        for iy in 1:ny
            i_row = ix + ( iy - 1 ) * nx
            values[ix,iy] = new_values_list[i_row]
        end
    end
end
function smooth_world_fields!( value_fields, diffcoeff )
    # used for orogeny footprint
    n_columns = nx * ny
    A = spzeros( n_columns, n_columns )
    for ix in 1:nx
        for iy in 1:ny    
            i_row = ix + ( iy - 1 ) * nx
            h_diffcoeffs = get_horizontal_diffcoeff(diffcoeff,iy)
            v_diffcoeffs = get_vertical_diffcoeff(diffcoeff,iy)

            A[i_row,i_row] = 1.
            # left
            A[i_row,i_row] += h_diffcoeffs[2] # left = -ve 
            i_neighbor_row = ix_left(ix) + ( iy - 1 ) * nx
            A[i_row,i_neighbor_row] = - h_diffcoeffs[2] # only fill one direction
            # right
            A[i_row,i_row] += h_diffcoeffs[1]
            i_neighbor_row = ix_right(ix) + ( iy - 1 ) * nx
            A[i_row,i_neighbor_row] = - h_diffcoeffs[1]
            if iy > 1
                A[i_row,i_row] += v_diffcoeffs[2]
                i_neighbor_row = ix + ( iy - 2 ) * nx
                A[i_row,i_neighbor_row] = - v_diffcoeffs[2]
            end
            if iy < ny
                A[i_row,i_row] += v_diffcoeffs[1]
                i_neighbor_row = ix + ( iy ) * nx
                A[i_row,i_neighbor_row] = - v_diffcoeffs[1]
            end
        end
    end
    
    lu_A = lu(A)

    for i_field in axes( value_fields )[3]
        R = fill(0.,n_columns)# Float64[]; 
        for ix in 1:nx
            for iy in 1:ny
                i_row = ix + ( iy - 1 ) * nx
                #push!(R, values[ix,iy])
                R[i_row] = value_fields[ix,iy,i_field]
            end
        end
        new_values_list = lu_A \ R
        for ix in 1:nx
            for iy in 1:ny
                i_row = ix + ( iy - 1 ) * nx
                value_fields[ix,iy,i_field] = new_values_list[i_row]
            end
        end
    end
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
#function predict_orogenic_erosion_fluxes()
function generate_orogenic_erosion_fluxes()
    # resets then accumulates land_orogenic_clay_flux and coastal_orogenic_clay_flux
    # also sets cruse_erosion_rate,
    # doesnt need to accumulate since there will be no overlap of orogenic areas.
    # uses exposed_basement
    verbose = false
    #reset_diag("land_orogenic_clay_flux")
    #reset_diag("coastal_orogenic_clay_flux")
    orogenic_crust_mask = eq_mask( world.geomorphology,exposed_basement )
    orogenic_clay_flux = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if orogenic_crust_mask[ix,iy] == 1
                orogenic_clay_flux[ix,iy] = 
                    max( calculate_orogenic_erosion_rate( world.freeboard[ix,iy] ), 0. )
            end
        end
    end
    return orogenic_clay_flux
end
function get_denuded_sediment_fluxes()
    orogenic_crust_mask = eq_mask( world.geomorphology,exposed_basement ) .*
        is_land() 
    denuded_sediment_fluxes = fill(0.,nx,ny,0:n_sediment_types)
    for ix in 1:nx
        for iy in 1:ny
            if orogenic_crust_mask[ix,iy] == 1
                for i_sedtype in 1:n_sediment_types
                    denuded_sediment_fluxes[ix,iy,i_sedtype] = 
                        world.sediment_thickness[ix,iy] * 
                        world.sediment_surface_fractions[ix,iy,i_sedtype] / main_time_step
                    denuded_sediment_fluxes[ix,iy,0] += denuded_sediment_fluxes[ix,iy,i_sedtype]
                end
            end
        end
    end
    return denuded_sediment_fluxes
end
function distribute_fluxes_uniformly_outside_boundary( fluxes, source_mask )
    output_fluxes = fill( 0.,nx,ny,0:n_sediment_types )
    blobs = get_blobs( source_mask )
    blob_areas = get_blob_areas( blobs )
    neighbor_blobs = get_blob_neighbor_fields( blobs )
    neighbor_blob_areas = get_blob_areas( neighbor_blobs )
    for i_area in 1:length( blobs )
        blob = blobs[ i_area ]
        neighbor_blob = neighbor_blobs[i_area] 
        if volume_field( fluxes[ :,:,0 ] .* blob ) > 0
            #println("distributing blob ",i_area)
            for i_sedtype in 1:n_sediment_types
                integrated_flux = volume_field( fluxes[ :,:,i_sedtype ] .* blob )
                deposition_rate = integrated_flux / neighbor_blob_areas[ i_area ]
                output_fluxes[:,:,i_sedtype] += neighbor_blob .* deposition_rate
                #=for ix in 1:nx
                    for iy in 1:ny
                        if neighbor_blobs[i_area][ix,iy] > 0
                            output_fluxes[ix,iy,i_sedtype] = deposition_flux
                        end
                    end
                end=#
            end
        end
    end
    update_flux_totals!( output_fluxes )
    return output_fluxes
end
#=
function expand_blob_boundary_subduction_to_fringe!( values, mask )
    for ix in 1:nx
        for iy in 1:ny
            if mask[ix,iy] == 1 && values[ix,iy] > 0.
                coordlist = get_border_neighbor_eq_size_coords(ix,iy)
                filtered_coordlist = []
                area_neighbors = 0.
                for coordpair in coordlist
                    ixn = coordpair[1]; iyn = coordpair[2]
                    if mask[ix,iy] == 0
                        push!( filtered_coordlist, coordpar )
                        area_neighbors += areabox[iyn]
                    end =#



function distribute_fluxes_uniformly_inside_boundary( fluxes, source_mask )
    output_fluxes = fill( 0.,nx,ny,0:n_sediment_types )
    blobs = get_blobs( source_mask )
    #blob_areas = get_blob_areas( blobs )
    fringe_blobs = get_blob_fringe_fields( blobs )
    fringe_blob_areas = get_blob_areas( fringe_blobs )
    for i_area in 1:length( blobs )
        blob = blobs[ i_area ]
        fringe_blob = fringe_blobs[i_area] 
        if volume_field( fluxes[ :,:,0 ] .* blob ) > 0
            #println("distributing blob ",i_area)
            for i_sedtype in 1:n_sediment_types
                integrated_flux = volume_field( fluxes[ :,:,i_sedtype ] .* blob )
                deposition_rate = integrated_flux / fringe_blob_areas[ i_area ]
                output_fluxes[:,:,i_sedtype] += fringe_blob .* deposition_rate
                #=for ix in 1:nx
                    for iy in 1:ny
                        if neighbor_blobs[i_area][ix,iy] > 0
                            output_fluxes[ix,iy,i_sedtype] = deposition_flux
                        end
                    end
                end=#
            end
        end
    end
    update_flux_totals!( output_fluxes )
    return output_fluxes
end
function distribute_fluxes_to_nearest_outside_boundary( fluxes, source_mask )
    output_fluxes = fill( 0.,nx,ny,0:n_sediment_types )
    blobs = get_blobs( source_mask )
    blob_areas = get_blob_areas( blobs )
    for i_area in 1:length( blobs )
        blob = blobs[ i_area ]
        if volume_field( fluxes[:,:,0] .* blob ) > 0.
            blob_fringe_list = get_blob_fringe_exterior_list( blob )
            ixn, iyn = find_nearest_fringe_point(ix,iy,blob_fringe_list)
            for i_sedtype in 1:n_sediment_types
                source_meter_flux = fluxes[ix,iy,i_sedtype]
                source_volume_flux = source_meter_flux  *
                    areabox[iy]
                sink_meter_flux = source_volume_flux / areabox[iyn]
                output_fluxes[ixn, iyn,i_sedtype] += 
                    sink_meter_flux
            end
        end
    end
    update_flux_totals!( output_fluxes )
    return output_fluxes
end
function update_flux_totals!(fluxfield)
    fluxfield[:,:,0] .= 0.
    for i_sedtype in 1:n_sediment_types
        fluxfield[:,:,0] += fluxfield[:,:,i_sedtype]
    end
end
#=
                denuded_sediment_flux[ix,iy,i_sedtype] -= source_meter_flux
                denuded_sediment_flux[ix,iy,0] -= source_meter_flux
        for i_sedtype in 1:n_sediment_types
            integrated_flux = volume_field( fluxes[ :,:,i_sedtype ] .* blob )
            deposition_flux = integrated_flux / neighbor_blob_areas[ i_blob ]
            for ix in 1:nx
                for iy in 1:ny
                    if neighbor_blobs[i_area][ix,iy] > 0
                        output_fluxes[ix,iy,i_sedtype] = deposition_flux
                    end
                end
            end
        end
    end
    return output_fluxes
end



        orogenic_mass_flux = orogenic_production_rates[i_area] / # m3 / Myr
            neighbor_blob_areas[i_area]  # meters / Myr
        denuded_mass_flux = denuded_sediment_production_rates[i_area] / 
            neighbor_blob_areas[i_area]
        for ix in 1:nx
            for iy in 1:ny
                if neighbor_blobs[i_area][ix,iy] > 0
                    if world.geomorphology[ix,iy] == sedimented_land
                        accum_diag("land_orogenic_clay_flux",ix,iy,orogenic_mass_flux)
                    elseif world.geomorphology[ix,iy] < sedimented_land # ocean
                        accum_diag("coastal_orogenic_clay_flux",ix,iy,orogenic_mass_flux)
                    end
                    for i_sedtype in 1:n_sediment_types
                        accum_frac_diag("denuded_land_boundary_fraction_flux",
                            ix,iy,i_sedtype,denuded_mass_flux)
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
end
=#
function generate_orogenic_erosion_fluxes_original()
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
    denuded_sediment_production_rates = get_integrated_denuded_sediment_rates( orogenic_blobs )
    # fills crust_erosion_rate and crust_clay_source_rate, 
    # returns blob-integrated rates in an array[nblobs]
    neighbor_blobs = get_blob_neighbor_fields( orogenic_blobs )
    neighbor_blob_areas = get_blob_areas( neighbor_blobs )
    for i_area in 1:length( orogenic_blobs )
        orogenic_mass_flux = orogenic_production_rates[i_area] / # m3 / Myr
            neighbor_blob_areas[i_area]  # meters / Myr
        denuded_mass_flux = denuded_sediment_production_rates[i_area] / 
            neighbor_blob_areas[i_area]
        for ix in 1:nx
            for iy in 1:ny
                if neighbor_blobs[i_area][ix,iy] > 0
                    if world.geomorphology[ix,iy] == sedimented_land
                        accum_diag("land_orogenic_clay_flux",ix,iy,orogenic_mass_flux)
                    elseif world.geomorphology[ix,iy] < sedimented_land # ocean
                        accum_diag("coastal_orogenic_clay_flux",ix,iy,orogenic_mass_flux)
                    end
                    for i_sedtype in 1:n_sediment_types
                        accum_frac_diag("denuded_land_boundary_fraction_flux",
                            ix,iy,i_sedtype,denuded_mass_flux)
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
function calculate_orogenic_erosion_rate( freeboard ) # meters / Myr
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

function get_integrated_denuded_sediment_rates( orogenic_blobs )
    denuded_sediment = get_frac_diag("denuded_sediment_source_fraction_flux")
    n_orogenies = length(orogenic_blobs)
    area_rates_by_blob = []
    for i_area = 1:n_orogenies
        blob = orogenic_blobs[i_area]
        rates_by_type = []
        for i_sedtype in 1:n_sediment_types
            total_sediment_type_source = volume_field( 
                denuded_sediment[:,:,i_sedtype] .* blob )
            push!( rates_by_type, total_sediment_type_source )
        end
        push!( area_rates_by_blob, rates_by_type )
    end
    return area_rates_by_blob
end

function get_integrated_orogenic_production_rates_unused_I_wonder( orogenic_blobs ) # m3 sediment / Myr
    # fills crust_erosion_rate and crust_clay_source_rate, 
    # returns blob-integrated rates in an array[nblobs]
    n_orogenies = length(orogenic_blobs)
    area_rates = []
    for i_area = 1:n_orogenies
        crust_erosion_field = fill(0.,nx,ny)
        for ix in 1:nx
            for iy in 1:ny
                if orogenic_blobs[i_area][ix,iy] == 1
                    crust_erosion_rate = calculate_orogenic_erosion_rate(world.freeboard[ix,iy]) # meters / Myr
                    crust_erosion_field[ix,iy] = crust_erosion_rate
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
