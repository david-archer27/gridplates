function get_scotese_mountains()
    age_string = lpad(Int(ceil(world.age)), 3, "0")
    CS_file_name = code_base_directory *
                   "/drivers/scotese_elevation_files_modern/scotese_elev_today." *
                   age_string * "Ma.bson"
    BSON.@load CS_file_name elevation
    rotated_elevation = rotate_field_0_to_current_age(elevation)

    #smooth_area!( values, blob, diffcoeff )
    smooth_world!( rotated_elevation, 1.e10 )
    set_diag("target_elevation", rotated_elevation)
    #target_elevation = fill_world_orogeny( twisted_elevation )

    return rotated_elevation
end
function rotate_field_0_to_current_age( field_at_time_0 )
    start_time_field = fill(1000.,nx,ny)
    field_past = 
        rotate_field_0_to_current_age( 
            field_at_time_0, start_time_field )
    return field_past
end

function mask_field_from_start_time_field( start_time_field )
    mask_field_time_0 = fill(1.,nx,ny)
    field_past = rotate_field_0_to_current_age(mask_field_time_0, start_time_field)
    return field_past
end

function rotate_field_0_to_current_age( field_at_time_0, start_time_field )
    field_past = fill(0.,nx,ny)
    for ixpast in 1:nx
        for iypast in 1:ny
            if world.crust_type[ixpast,iypast] == continent_crust
                #error(ixpast," ",iypast)
                rot_ID = world.continentID[ixpast,iypast]
                rotation_matrix = resolve_rotation_matrix!( rot_ID, world.age )
                ix0, iy0 = nearest_plateij(rotation_matrix, ixpast, iypast)
                # nearest plate because that is the 0-rotation location 
                if start_time_field[ix0,iy0] > world.age 
                    field_past[ixpast, iypast] = field_at_time_0[ix0, iy0]
                end
            end
        end
    end
 
    return field_past
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
                current_ophiolites[ix, iy] == 1.0 &&
                world.crust_type[ix,iy] == continent_crust

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
function get_ice_sheets()
    age_string = lpad(Int(world.age),3,"0")
    cd( "drivers/scotese_paleo_temp_files_modern")
    file_list = readdir()
    file_name = "scotese_ptemp." *
        age_string * "Ma.bson"
    ice_sheet_mask = fill(0.,nx,ny)
    if file_name in file_list
        BSON.@load file_name paleo_temp_modern
        paleo_temp = rotate_field_0_to_current_age( paleo_temp_modern )
        ice_sheet_mask = lt_mask(paleo_temp, -10) 
    end
    for ix in 1:nx
        for iy in 1:ny
            if world.geomorphology[ix,iy] == ice_sheet_covered &&
                ice_sheet_mask[ix,iy] == 0.
                    world.geomorphology[ix, iy] = exposed_basement
            end
            if ice_sheet_mask[ix, iy] > 0.0 &&
                world.crust_type[ix,iy] == continent_crust
                    world.geomorphology[ix, iy] = ice_sheet_covered
            end
        end
    end
    cd( code_base_directory )
    return ice_sheet_mask
end
function uplift_scotese_mountains_by_crust_thickening( scotese_target_elevation )
    for ix in 1:nx
        for iy in 1:ny
            if scotese_target_elevation[ix, iy] > 1500.0 && 
                world.crust_type[ix,iy] == continent_crust

                crust_thickening_rate = ( scotese_target_elevation[ix, iy] - world.freeboard[ix, iy] ) /
                    crust_freeboard_expression
                world.crust_thickness[ix,iy] += crust_thickening_rate
                set_diag("crust_thickening_rate",ix,iy,crust_thickening_rate ./ main_time_step)
                # setting diags and changing world is against style, fix it sometime
            end
        end
    end
    smooth_continental_crust_thickness()
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

function generate_crust_erosion_fluxes()
    # generated independently from chemical weathering of crust 
    # compare "river fluxes" "particulate vs dissolved" ratio
    erosional_sediment_fluxes = fill(0.0, nx, ny, 0:n_sediment_types)
    ice_sheet_sediment_fluxes = fill(0.0, nx, ny, 0:n_sediment_types)
    for ix in 1:nx
        for iy in 1:ny
            erosion_rate = 0.
            if world.geomorphology[ix, iy] == ice_sheet_covered
                erosion_rate = ice_sheet_erosion_rate # meters / Myr
                ice_sheet_sediment_fluxes[ix,iy,:] = 
                    partition_erosion_sediment_fluxes( ix,iy,erosion_rate )
            end
            if world.geomorphology[ix, iy] == exposed_basement
                tau = orogenic_erosion_tau_apparent * crust_freeboard_expression
                erosion_rate = max(world.freeboard[ix, iy] / tau, 0.0)
                erosional_sediment_fluxes[ix,iy,:] = 
                    partition_erosion_sediment_fluxes( ix,iy,erosion_rate )
            end
        end
    end
    return erosional_sediment_fluxes, ice_sheet_sediment_fluxes
end
function partition_erosion_sediment_fluxes( ix,iy,erosion_rate )
    erosion_sediment_fluxes = fill(0.,0:n_sediment_types)
    erosion_sediment_fluxes[CaO_sediment] = erosion_rate * (
        world.crust_composition[ix, iy] * felsic_CaO_fraction +
        (1 - world.crust_composition[ix, iy]) * ultramafic_CaO_fraction ) 
    erosion_sediment_fluxes[clay_sediment] = erosion_rate - 
        erosion_sediment_fluxes[CaO_sediment]
    erosion_sediment_fluxes[0] = erosion_sediment_fluxes[CaO_sediment] +
        erosion_sediment_fluxes[clay_sediment]
    return erosion_sediment_fluxes
end

function get_denuded_sediment_fluxes( denuded_mask )
    denuded_sediment_fluxes = fill(0.,nx,ny,0:n_sediment_types)
    for ix in 1:nx
        for iy in 1:ny
            if denuded_mask[ix,iy] == 1
                for i_sedtype in first_land_transported_sediment:n_sediment_types
                    denuded_sediment_fluxes[ix,iy,i_sedtype] = 
                        world.sediment_inventories[ix,iy,i_sedtype] / main_time_step
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
