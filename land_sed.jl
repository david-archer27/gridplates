function land_bulk_sediment_transport( original_elevation_field, diffusive_mask, 
    bulk_sediment_source, ocean_sink )
    land_sediment_deposition_rates = fill(0.,nx,ny)
    subaereal_blobs = get_blobs( diffusive_mask )
    altered_elevation_field = deepcopy( original_elevation_field )
    verbose = false
    Threads.@threads for i_blob in 1:length(subaereal_blobs) 
        #if verbose == true
        #    println("beginning ",i_blob)
        #end
        subaereal_blob = subaereal_blobs[i_blob]

        blob_sediment_source = volume_field( bulk_sediment_source .* subaereal_blob )
        
        land_blob_bulk_sediment_transport!( altered_elevation_field, 
            subaereal_blob, bulk_sediment_source, ocean_sink )
        for ix in 1:nx
            for iy in 1:ny
                if subaereal_blob[ix,iy] == 1
                    land_sediment_deposition_rates[ix,iy] = 
                        ( altered_elevation_field[ix,iy] - original_elevation_field[ix,iy] ) /
                        time_step / sediment_freeboard_expression
                end
            end
        end 

        blob_sediment_dep = volume_field( land_sediment_deposition_rates .* subaereal_blob )
        #println("blob sources ",[i_blob,blob_sediment_source,blob_sediment_dep,
        #    blob_sediment_source - blob_sediment_dep])

    end # blob 
    return altered_elevation_field, land_sediment_deposition_rates
end
function land_blob_bulk_sediment_transport!( altered_elevation_field, 
    subaereal_blob, bulk_sediment_source, ocean_sink )
    # resets new_elevation_field 
    list_pos = fill(0,nx,ny)
    list_pos_filled = 0
    list_x = []; list_y = []
    original_elevation_list = []
    for ix in 1:nx
        for iy in 1:ny
            if subaereal_blob[ix,iy] == 1
                list_pos_filled += 1
                list_pos[ix,iy] = list_pos_filled
                push!( list_x, ix ); push!( list_y, iy )
            end
        end
    end
    if list_pos_filled == 1
        ix = list_x[1]; iy = list_y[1]
        h_diffcoeffs_scaled = get_horizontal_diffcoeff(land_base_diffcoeff,iy) .* sediment_freeboard_expression
        v_diffcoeffs_scaled = get_vertical_diffcoeff(land_base_diffcoeff,iy) .* sediment_freeboard_expression
        A = Float64[ 1. ]
        #runoff_loss = 0.
        if ocean_sink[ix_left(ix),iy] == 1 
            A[1] += h_diffcoeffs_scaled[2] 
            #runoff_loss += elevation_field[ix,iy] * h_diffcoeffs[2] # meters / step
        end
        if ocean_sink[ix_right(ix),iy] == 1 
            A[1] += h_diffcoeffs_scaled[1]
            #runoff_loss += elevation_field[ix,iy] * h_diffcoeffs[1] # meters / step
        end
        if iy > 1
            if ocean_sink[ix,iy-1] == 1 
                A[1] += v_diffcoeffs_scaled[2]
            #runoff_loss += elevation_field[ix,iy] * v_diffcoeffs[2] # meters / step
            end
        end
        if iy < ny
            if ocean_sink[ix,iy+1] == 1 
                A[1] += v_diffcoeffs_scaled[1]
            #runoff_loss += elevation_field[ix,iy] * v_diffcoeffs[1] # meters / step
            end
        end
        R = Float64[ altered_elevation_field[ix,iy] +
            bulk_sediment_source[ix,iy] * # meters/Myr
            time_step * # meters / step
            sediment_freeboard_expression ]
        new_elevation = A \ R
        altered_elevation_field[ix,iy] = new_elevation
    else # more than one point

        # basic equation for influx point.  prime denotes end-of-step state
        # H1' = H1 + D1 ( E2' - E1' ) + Src1
        # adjacent coastal (sink) point
        # H2' = H2 + D2 ( E1' - E2') + D1 ( 0 - E2' ) 
        # expanded to deal in elevations
        # H = E / L (L = 0.4 sediment isostatic exposed iceberg fraction)
        # E1' = E1 + D1 * L ( E2' - E1' ) + Src1 * L
        # E2' = E2 + D2 * L ( E1' - E2' ) + D1 * L ( 0 - E2' ) 
        # rearrange, set D = D * L 
        # E1' * [ 1 + D1 ] + E2' * [ - D1 ] = E1 + Src1 * L
        # E1' * [  - D2 (E1'-E2')/(2E1')] + E2' * [ 1 + D2 + D1 ] = E2
        # matrix form
        #| 1 + D1, -D1      | | E1' |      | E1 + Src * L |
        #| -D2, 1 + D2 + D1 | | E2' |  =   | E2           |

        A = spzeros( list_pos_filled, list_pos_filled )
        R = Float64[]; original_elevation_list = []
        for i_pos = 1:list_pos_filled
            ix = list_x[i_pos]; iy = list_y[i_pos]
                h_diffcoeffs_scaled = get_horizontal_diffcoeff(land_base_diffcoeff,iy) .* sediment_freeboard_expression
                v_diffcoeffs_scaled = get_vertical_diffcoeff(land_base_diffcoeff,iy) .* sediment_freeboard_expression
            push!(R, altered_elevation_field[ix,iy])
            push!(original_elevation_list, altered_elevation_field[ix,iy])
            A[i_pos,i_pos] = 1.
            if subaereal_blob[ix_left(ix),iy] == 1 
                A[i_pos,i_pos] += h_diffcoeffs_scaled[2] # left = -ve 
                i_neighbor_pos = list_pos[ix_left(ix),iy]
                A[i_pos,i_neighbor_pos] = - h_diffcoeffs_scaled[2] # only fill one direction
            elseif ocean_sink[ix_left(ix),iy] == 1 
                A[i_pos,i_pos] += h_diffcoeffs_scaled[2]
            end
            if subaereal_blob[ix_right(ix),iy] == 1 
                A[i_pos,i_pos] += h_diffcoeffs_scaled[1]
                i_neighbor_pos = list_pos[ix_right(ix),iy]
                A[i_pos,i_neighbor_pos] = - h_diffcoeffs_scaled[1]
            elseif ocean_sink[ix_right(ix),iy] == 1 
                A[i_pos,i_pos] += h_diffcoeffs_scaled[1]
            end
            if iy > 1
                if subaereal_blob[ix,iy-1] == 1
                    A[i_pos,i_pos] += v_diffcoeffs_scaled[2]
                    i_neighbor_pos = list_pos[ix,iy-1]
                    A[i_pos,i_neighbor_pos] = - v_diffcoeffs_scaled[2]
                elseif ocean_sink[ix,iy-1] == 1 
                    A[i_pos,i_pos] += v_diffcoeffs_scaled[2]
                end
            end
            if iy < ny
                if subaereal_blob[ix,iy+1] == 1
                    A[i_pos,i_pos] += v_diffcoeffs_scaled[1]
                    i_neighbor_pos = list_pos[ix,iy+1]
                    A[i_pos,i_neighbor_pos] = - v_diffcoeffs_scaled[1]
                elseif ocean_sink[ix,iy+1] == 1 
                    A[i_pos,i_pos] += v_diffcoeffs_scaled[1]
                end
            end
            R[end] += bulk_sediment_source[ix,iy] * # meters / Myr
                time_step * sediment_freeboard_expression # meters
        end

        new_elevation_list = lu(A) \ R

        for i_pos = 1:list_pos_filled
            ix = list_x[i_pos]; iy = list_y[i_pos]
            #sed_depo = ( new_elevation_list[i_pos] - elevation_field[ix,iy] ) /
            #    time_step # meters / Myr
            #set_diag("sediment_deposition_rate",ix,iy,sed_depo)
            altered_elevation_field[ix,iy] = new_elevation_list[i_pos]
        end
    end # diffuse or no
end
function check_for_bedrock_exposure( land_sediment_deposition_rate_field )
    # sets world.geomorphology to exclude for next pass,
    # puts the spoil into the nearest denuded_land_boundary_fraction_flux,
    # records it in denuded_sediment_source_fraction_flux
    # called in a loop to calculate the bulk fluxes, but the fractions denuded
    # have to be preserved 
    n_denuded = 0
    sediment_deposition_meters = 
    land_sediment_deposition_rate_field .* time_step
    total_sediment_thickness = world.sediment_thickness
    new_total_sediment_thickness = total_sediment_thickness .+ 
        sediment_deposition_meters 
    for ix in 1:nx
        for iy in 1:ny
            if world.geomorphology[ix,iy] == sedimented_land &&
                new_total_sediment_thickness[ix,iy] < 0.
                #println("found ",[ix,iy])
                n_denuded += 1
                world.geomorphology[ix,iy] = exposed_basement
            end
        end
    end
    return n_denuded
end
function land_sediment_fraction_transport( new_elevation_field, 
    land_sediment_deposition_rate_field,
    diffusive_mask, 
    sediment_source_fields, ocean_sink )
    # sets land_sediment_fraction_deposition_rate on land areas
    combined_land_sediment_fraction_deposition_rates = 
        fill(0.,nx,ny,n_sediment_types)
    sediment_deposition_meters = 
        land_sediment_deposition_rate_field .* time_step
    old_total_sediment_thickness = world.sediment_thickness
    new_total_sediment_thickness = old_total_sediment_thickness .+
        sediment_deposition_meters
    covered_diffusive_mask = diffusive_mask .- le_mask( new_total_sediment_thickness, 0. )
    subaereal_blobs = get_blobs( covered_diffusive_mask )

    Threads.@threads for i_blob in 1:length(subaereal_blobs) 
        #println("beginning ",i_blob)
        subaereal_blob = subaereal_blobs[i_blob]
        land_sediment_fraction_deposition_rates = 
            land_area_sediment_fraction_transport( 
                new_elevation_field, land_sediment_deposition_rate_field,
                subaereal_blob, sediment_source_fields, ocean_sink ) 
        for ix in 1:nx
            for iy in 1:ny
                if subaereal_blob[ix,iy] == 1
                    for i_sedtype in 1:n_sediment_types
                        if land_sediment_fraction_deposition_rates[ix,iy,i_sedtype] != 0.
                            combined_land_sediment_fraction_deposition_rates[ix,iy,i_sedtype] = 
                                land_sediment_fraction_deposition_rates[ix,iy,i_sedtype]
                        end
                    end 
                end
            end
        end
    end
    for ix in 1:nx
        for iy in 1:ny
            if world.geomorphology[ix,iy] == exposed_basement
                for i_sedtype in 1:n_sediment_types
                    combined_land_sediment_fraction_deposition_rates[ix,iy,i_sedtype] =
                        - world.sediment_thickness[ix,iy] * 
                        world.sediment_surface_fractions[ix,iy,i_sedtype] / time_step
                end
            end
        end
    end 
    return combined_land_sediment_fraction_deposition_rates
end
function land_area_sediment_fraction_transport( 
    new_elevation_field, land_sediment_deposition_rate_field,
    subaereal_blob, sediment_source_fields, ocean_sink ) 
    # sets land_sediment_fraction_deposition_rate on land points
    # fractions combined in one call so they can share the 
    # transportation matrix

    # basic equation for influx point.  prime denotes end-of-step state.
    # concentrations are averaged over entire sediment column. 
    # Inv1' = Inv1 + D1 ( E2' - E1' ) * ( C1' + C2' ) / 2 + Src1
    #                -- mass flux ---   ---- avg conc ---
    # adjacent coastal (sink) point 3
    # Inv2' = Inv2 + D2 ( E1' - E2') * (C1' + C2')/2 - D1 E2' C2'

    # rearrange to deal in inventories 
    # conc for transport is averaged over entire depth 
    # C = Inv(column)/H(column)
    # Inv1' = Inv1 + D1 ( E2' - E1' ) * ( Inv1'/H1' + Inv2'/H2' ) / 2 + Src1
    # Inv2' = Inv2 + D2 ( E1' - E2' ) * ( Inv1'/H1' + Inv2'/H2' ) / 2 - D1 E2' Inv2' / H2'

    # group unknown end-step Inv terms 
    # Inv1' * [1 - D1 (E2'-E1')/(2 H1')] + Inv2' * [  - D1 (E2'-E1')/(2 H2')                ] = Inv1 + Src1
    # Inv1' * [  - D2 (E1'-E2')/(2 H1')] + Inv2' * [1 - D2 (E1'-E2')/(2 H2') + D1 E2' / H2' ] = Inv2

    # matrix form
    #| 1 - D1 (E2'-E1') / (2 H1'),   - D1 (E2'-E1') / (2 H2')                | | Inv1' |      | Inv1 + Src |
    #|   - D2 (E1'-E2') / (2 H1'), 1 - D2 (E1'-E2') / (2 H2') + D1 E2' / H2' | | Inv2' |  =   | Inv2       |

    land_sediment_fraction_deposition_rates = fill(0.,nx,ny,n_sediment_types)
    sediment_deposition_meters = 
        land_sediment_deposition_rate_field .* time_step
    old_total_sediment_thickness = world.sediment_thickness
    new_total_sediment_thickness = old_total_sediment_thickness .+
        sediment_deposition_meters
    old_column_fraction_inventory = fill(0.,nx,ny,n_sediment_types)
    for i_sedtype in 1:n_sediment_types
        old_column_fraction_inventory[:,:,i_sedtype] = 
            world.sediment_thickness .*
            world.sediment_surface_fractions[:,:,i_sedtype]
    end

    list_pos = fill(0,nx,ny)
    list_pos_filled = 0
    list_x = []; list_y = []
    #subaereal_blob -= le_mask( new_total_sediment_thickness, 0. )
    for ix in 1:nx
        for iy in 1:ny
            if subaereal_blob[ix,iy] == 1  
                #        should be prevented by check_subaereal_exposure 
                list_pos_filled += 1
                list_pos[ix,iy] = list_pos_filled
                push!( list_x, ix ); push!( list_y, iy )
            end
        end
    end
    if list_pos_filled == 1
        A = 1.
        ix = list_x[1]; iy = list_y[1]
        if new_total_sediment_thickness[ix,iy] > 0.
            h_diffcoeffs = get_horizontal_diffcoeff(land_base_diffcoeff,iy)
            v_diffcoeffs = get_vertical_diffcoeff(land_base_diffcoeff,iy)
            if ocean_sink[ix_left(ix),iy] == 1 
                A += h_diffcoeffs[2] * 
                    new_elevation_field[ix,iy] / new_total_sediment_thickness[ix,iy]
            end
            if ocean_sink[ix_right(ix),iy] == 1 
                A += h_diffcoeffs[2] * 
                    new_elevation_field[ix,iy] / new_total_sediment_thickness[ix,iy]
            end
            if iy > 1
                if ocean_sink[ix,iy-1] == 1
                    A += v_diffcoeffs[2] * 
                    new_elevation_field[ix,iy] / new_total_sediment_thickness[ix,iy]
                end
            end
            if iy < ny
                if ocean_sink[ix,iy+1] == 1
                    A += v_diffcoeffs[1] * 
                    new_elevation_field[ix,iy] / new_total_sediment_thickness[ix,iy]
                end
            end
            for i_sedtype in 1:n_sediment_types

                #=  old_inventory = world.sediment_thickness[ix,iy,current_time_bin()] * 
                        world.sediment_fractions[ix,iy,i_sedtype,current_time_bin()]
                remaining_meters = - (
                    sediment_deposition_meters[ix,iy] +
                    world.sediment_thickness[ix,iy,current_time_bin()] 
                    )
                i_bin = current_time_bin()
                while remaining_meters > 0. && i_bin > 1
                    i_bin -= 1
                    remaining_meters -= world.sediment_thickness[ix,iy,i_bin]
                    if remaining_meters < 0. # layer i_bin has enough
                        old_inventory += remaining_meters * 
                            world.sediment_fractions[ix,iy,i_sedtype,i_bin]
                    else
                        old_inventory += world.sediment_thickness[ix,iy,i_bin] * 
                            world.sediment_fractions[ix,iy,i_sedtype,i_bin]
                    end
                end =#
                old_inventory = old_column_fraction_inventory[ix,iy,i_sedtype]
                R = old_inventory + 
                    sediment_source_fields[ix,iy,i_sedtype] * time_step
                new_inventory = A \ R   
                #new_fraction = new_inventory / 
                #    new_total_sediment_thickness[ix,iy]
                fraction_deposition_rate = 
                    ( new_inventory - old_inventory ) / 
                    time_step
                land_sediment_fraction_deposition_rates[ix,iy,i_sedtype] =
                    fraction_deposition_rate
                    #set_frac_diag("land_sediment_fraction_deposition_rate",
                #    ix,iy,i_sedtype,fraction_deposition_rate)
            end
        end
    elseif list_pos_filled > 1# more than one point
        A = spzeros( list_pos_filled, list_pos_filled )
        for i_pos = 1:list_pos_filled
            ix = list_x[i_pos]; iy = list_y[i_pos]
            h_diffcoeffs = get_horizontal_diffcoeff(land_base_diffcoeff,iy)
            # no isostatic rebound factor here -- getting D dt / dx, gives thickness transport
            v_diffcoeffs = get_vertical_diffcoeff(land_base_diffcoeff,iy)
            # up (toward positive) = 1, down = 2

            A[i_pos,i_pos] = 1.
            if subaereal_blob[ix_left(ix),iy] == 1 
                A[i_pos,i_pos] -= h_diffcoeffs[2] * # diffusion to negative direction
                    ( new_elevation_field[ix_left(ix),iy] - new_elevation_field[ix,iy] ) / 
                    2 / new_total_sediment_thickness[ix,iy]
                i_neighbor_pos = list_pos[ix_left(ix),iy]
                if i_neighbor_pos == 0
                    println("ix ",ix_left(ix)," iy ",iy)
                end
                A[i_pos,i_neighbor_pos] -= h_diffcoeffs[2] * 
                    ( new_elevation_field[ix_left(ix),iy] - new_elevation_field[ix,iy] ) /  
                    2 / new_total_sediment_thickness[ix_left(ix),iy]
            elseif ocean_sink[ix_left(ix),iy] == 1 
                A[i_pos,i_pos] += h_diffcoeffs[2] * 
                    new_elevation_field[ix,iy] / new_total_sediment_thickness[ix,iy]
            end
            if subaereal_blob[ix_right(ix),iy] == 1 
                A[i_pos,i_pos] -= h_diffcoeffs[1] * 
                    ( new_elevation_field[ix_right(ix),iy] - new_elevation_field[ix,iy] ) /  
                    2 / new_total_sediment_thickness[ix,iy]
                i_neighbor_pos = list_pos[ix_right(ix),iy]
                if i_neighbor_pos == 0
                    println("ix ",ix_right(ix)," iy ",iy)
                end
                A[i_pos,i_neighbor_pos] -= h_diffcoeffs[1] * 
                    ( new_elevation_field[ix_right(ix),iy] - new_elevation_field[ix,iy] ) /  
                    2 / new_total_sediment_thickness[ix_right(ix),iy]
            elseif ocean_sink[ix_right(ix),iy] == 1 
                A[i_pos,i_pos] += h_diffcoeffs[1] * 
                    new_elevation_field[ix,iy] / new_total_sediment_thickness[ix,iy]
            end   
            if iy == 1 
                #= tie things across the pole 
                ixn = ( ix + 179 ) % 360 + 1
                if subaereal_blob[ixn,iy] == 1
                    A[i_pos,i_pos] += v_diffcoeffs[2] * 
                        ( new_elevation_field[ixn,iy] - new_elevation_field[ix,iy] ) /  
                        2 / new_sediment_thickness[ix,iy]
                    i_neighbor_pos = list_pos[ixn,iy]
                    A[i_pos,i_neighbor_pos] = - v_diffcoeffs[2] * 
                        ( new_elevation_field[ixn,iy] - new_elevation_field[ix,iy] ) /  
                        2 / new_sediment_thickness[ixn,iy]
                end =#
            else # iy > 1
                if subaereal_blob[ix,iy-1] == 1 
                    A[i_pos,i_pos] -= v_diffcoeffs[2] * # index 2 is toward -ve in iy
                        ( new_elevation_field[ix,iy-1] - new_elevation_field[ix,iy] ) /  
                        2 / new_total_sediment_thickness[ix,iy]
                    i_neighbor_pos = list_pos[ix,iy-1]
                    if i_neighbor_pos == 0
                        println("ix ",ix," iy ",iy-1)
                    end
                    A[i_pos,i_neighbor_pos] -= v_diffcoeffs[2] * 
                        ( new_elevation_field[ix,iy-1] - new_elevation_field[ix,iy] ) /  
                        2 / new_total_sediment_thickness[ix,iy-1]
                elseif ocean_sink[ix,iy-1] == 1 
                    A[i_pos,i_pos] += v_diffcoeffs[2] * 
                        new_elevation_field[ix,iy] / new_total_sediment_thickness[ix,iy]
                end    
            end        
            if iy == ny  
                #=  
                ixn = ( ix + 179 ) % 360 + 1
                if subaereal_blob[ixn,iy] == 1
                    A[i_pos,i_pos] += v_diffcoeffs[1]  * 
                        ( new_elevation_field[ixn,iy] - new_elevation_field[ix,iy] ) /  
                        2 / new_sediment_thickness[ix,iy]
                    i_neighbor_pos = list_pos[ixn,iy]
                    A[i_pos,i_neighbor_pos] = - v_diffcoeffs[1] * 
                        ( new_elevation_field[ixn,iy] - new_elevation_field[ix,iy] ) /  
                        2 / new_sediment_thickness[ixn,iy]
                end =#
            else # iy < ny
                if subaereal_blob[ix,iy+1] == 1 
                    A[i_pos,i_pos] -= v_diffcoeffs[1] * 
                        ( new_elevation_field[ix,iy+1] - new_elevation_field[ix,iy] ) /  
                        2 / new_total_sediment_thickness[ix,iy]
                    i_neighbor_pos = list_pos[ix,iy+1]
                    if i_neighbor_pos == 0
                        println("ix ",ix," iy ",iy+1)
                        println( list_pos[ix-1:ix+1,iy-1:iy+1], subaereal_blob[ix-1:ix+1,iy-1:iy+1])
                    end
                    A[i_pos,i_neighbor_pos] -= v_diffcoeffs[1] * 
                        ( new_elevation_field[ix,iy+1] - new_elevation_field[ix,iy] ) /  
                        2 / new_total_sediment_thickness[ix,iy+1]
                elseif ocean_sink[ix,iy+1] == 1 
                    A[i_pos,i_pos] += v_diffcoeffs[1] * 
                        new_elevation_field[ix,iy] / new_total_sediment_thickness[ix,iy]
                end
            end                
        end # create A

        lu_A = lu(A)

        for i_sedtype in 1:n_sediment_types 
            old_inventory_list = []
            R = Float64[]; 
            for i_pos = 1:list_pos_filled
                ix = list_x[i_pos]; iy = list_y[i_pos]

                #= old_inventory = world.sediment_thickness[ix,iy,current_time_bin()] * 
                    world.sediment_fractions[ix,iy,i_sedtype,current_time_bin()]
                remaining_meters = - (
                    sediment_deposition_meters[ix,iy] +
                    world.sediment_thickness[ix,iy,current_time_bin()] 
                    )
                i_bin = current_time_bin()
                while remaining_meters > 0. && i_bin > 1
                    i_bin -= 1
                    remaining_meters -= world.sediment_thickness[ix,iy,i_bin]
                    if remaining_meters < 0. # layer i_bin has enough
                        old_inventory += remaining_meters * 
                            world.sediment_fractions[ix,iy,i_sedtype,i_bin]
                    else
                        old_inventory += world.sediment_thickness[ix,iy,i_bin] * 
                            world.sediment_fractions[ix,iy,i_sedtype,i_bin]
                    end
                end =#
                old_inventory = old_column_fraction_inventory[ix,iy,i_sedtype]
                push!(old_inventory_list,old_inventory)
                push!(R, old_inventory + 
                    sediment_source_fields[ix,iy,i_sedtype] * time_step )
            end
            new_inventory_list = lu_A \ R
            #new_fraction_field = fill(0.,nx,ny)
            for i_pos = 1:list_pos_filled
                if new_inventory_list[i_pos] != new_inventory_list[i_pos]
                    error("NaN detected in land_area_sediment_fraction_transport")
                end
                ix = list_x[i_pos]; iy = list_y[i_pos]
                #new_fraction = new_inventory_list[i_pos] / 
                #    new_total_layer_sediment_thickness[ix,iy]
                #new_fraction_field[ix,iy] = new_fraction
                fraction_deposition_rate = 
                    ( new_inventory_list[i_pos] - old_inventory_list[i_pos] ) / 
                    time_step
                land_sediment_fraction_deposition_rates[ix,iy,i_sedtype] =
                    fraction_deposition_rate
            end
        end
    end 
    return land_sediment_fraction_deposition_rates
end
function land_transport_elevation_unused()
    water_line = min( world.sealevel, crust_level())
    crust_elevation = world.surface_elevation .- water_line
    return crust_elevation
end
function land_bulk_runoff_fluxes( new_elevation_field, subaereal_mask, ocean_sink )
        # resets and fills coastal_sediment_fraction_runoff_flux
    # sets world.geomorphology for coastal_depocenter
    #reset_frac_diag("coastal_sediment_fraction_runoff_flux")
    coastal_sediment_runoff_flux = fill(0.,nx,ny)
    # we want to get a preview without altering the real world 

    # basic equation for adjacent coastal (sink) point
    # H2' = H2 + D2 ( E1' - E2') + D1 ( 0 - E2' ) 
    # term D1 E2 is change in H. 
    # volume flux = D1 E2 area2
    # adjacent deposition flux = D1 E2 area2 / area3

    for ix in 1:nx
        for iy in 1:ny
            if subaereal_mask[ix,iy] > 0
                h_diffcoeffs = get_horizontal_diffcoeff(land_base_diffcoeff,iy)
                v_diffcoeffs = get_vertical_diffcoeff(land_base_diffcoeff,iy)
                if ocean_sink[ix_left(ix),iy] > 0.
                    volume_flux = new_elevation_field[ix,iy] * h_diffcoeffs[2] * # meters / step
                        areabox[iy] / time_step  # m3 / Myr
                    coastal_sediment_runoff_flux[ix_left(ix),iy] +=
                        volume_flux / areabox[iy]
                end
                if ocean_sink[ix_right(ix),iy] > 0.
                    volume_flux = new_elevation_field[ix,iy] * h_diffcoeffs[1] * 
                        areabox[iy] / time_step
                    coastal_sediment_runoff_flux[ix_right(ix),iy] +=
                        volume_flux / areabox[iy]
                end
                if iy > 1
                    if ocean_sink[ix,iy-1] > 0  
                        volume_flux = new_elevation_field[ix,iy] * v_diffcoeffs[2] * 
                            areabox[iy] / time_step
                        coastal_sediment_runoff_flux[ix,iy-1] +=
                            volume_flux / areabox[iy-1]
                    end
                end
                if iy < ny
                    if ocean_sink[ix,iy+1] > 0 
                        volume_flux = new_elevation_field[ix,iy] * v_diffcoeffs[1] * 
                            areabox[iy] / time_step
                        coastal_sediment_runoff_flux[ix,iy+1] +=
                            volume_flux / areabox[iy+1]
                    end
                end
            end
        end
    end 
    return coastal_sediment_runoff_flux
end
function land_fraction_runoff_fluxes( new_elevation_field, 
    land_sediment_fraction_deposition_rate_fields,
    diffusive_mask, ocean_sink )    

    new_sediment_thickness, new_sediment_fractions = 
        apply_land_sediment_fluxes( 
            land_sediment_fraction_deposition_rate_fields ) 
    coastal_sediment_fraction_runoff_flux = 
        fill(0.,nx,ny,n_sediment_types)

    # basic equation for adjacent coastal (sink) point
    # H2' = H2 + D2 ( E1' - E2') + D1 ( 0 - E2' ) 
    # term D1 E2 is change in H. 
    # volume flux = D1 E2 area2
    # adjacent deposition flux = D1 E2 area2 / area3

    for ix in 1:nx
        for iy in 1:ny
            if diffusive_mask[ix,iy] > 0
                h_diffcoeffs = get_horizontal_diffcoeff(land_base_diffcoeff,iy)
                v_diffcoeffs = get_vertical_diffcoeff(land_base_diffcoeff,iy)
                if ocean_sink[ix_left(ix),iy] > 0.
                    #println(ix," ",iy); return
                    for i_sedtype in 1:n_sediment_types
                        volume_flux = new_elevation_field[ix,iy] * h_diffcoeffs[2] * # meters / step
                            new_sediment_fractions[ix,iy,i_sedtype] * 
                            areabox[iy] / time_step  # m3 / Myr
                        coastal_sediment_fraction_runoff_flux[ix_left(ix),iy,i_sedtype] +=
                            volume_flux / areabox[iy]
                        #accum_frac_diag( "coastal_sediment_fraction_runoff_flux",
                        #    ix_left(ix),iy,i_sedtype,
                        #    volume_flux / areabox[iy] )
                    end
                    world.geomorphology[ix_left(ix),iy] = coastal_depocenter
                end
                if ocean_sink[ix_right(ix),iy] > 0.
                    for i_sedtype in 1:n_sediment_types
                        volume_flux = new_elevation_field[ix,iy] * h_diffcoeffs[1] * 
                            new_sediment_fractions[ix,iy,i_sedtype] * 
                            areabox[iy] / time_step
                        coastal_sediment_fraction_runoff_flux[ix_right(ix),iy,i_sedtype] +=
                            volume_flux / areabox[iy]
                        #accum_frac_diag( "coastal_sediment_fraction_runoff_flux",
                        #    ix_right(ix),iy,i_sedtype,
                        #    volume_flux / areabox[iy] )
                    end
                    world.geomorphology[ix_right(ix),iy] = coastal_depocenter
                end
                if iy > 1
                    if ocean_sink[ix,iy-1] > 0  
                        for i_sedtype in 1:n_sediment_types
                            volume_flux = new_elevation_field[ix,iy] * v_diffcoeffs[2] * 
                                new_sediment_fractions[ix,iy,i_sedtype] *
                                areabox[iy] / time_step
                            coastal_sediment_fraction_runoff_flux[ix,iy-1,i_sedtype] +=
                                volume_flux / areabox[iy-1]
                            #accum_frac_diag( "coastal_sediment_fraction_runoff_flux",
                            #    ix,iy-1,i_sedtype,
                            #    volume_flux / areabox[iy-1] )
                        end
                        world.geomorphology[ix,iy-1] = coastal_depocenter
                    end
                end
                if iy < ny
                    if ocean_sink[ix,iy+1] > 0 
                        for i_sedtype in 1:n_sediment_types
                            volume_flux = new_elevation_field[ix,iy] * v_diffcoeffs[1] * 
                                new_sediment_fractions[ix,iy,i_sedtype] *
                                areabox[iy] / time_step
                            coastal_sediment_fraction_runoff_flux[ix,iy+1,i_sedtype] +=
                                volume_flux / areabox[iy+1]
                            #accum_frac_diag( "coastal_sediment_fraction_runoff_flux",
                            #    ix,iy+1,i_sedtype,
                            #    volume_flux / areabox[iy+1] )
                        end
                        world.geomorphology[ix,iy+1] = coastal_depocenter
                    end
                end
            end
        end
    end 
    return coastal_sediment_fraction_runoff_flux
end
function aolean_transport()
    # resets and fills aolean_deposition_rate and aolean_deposition_rate
    aolean_erosion_rates = fill(0.,nx,ny)
    aolean_deposition_rates = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if world.freeboard[ix,iy] > 0. &&
                world.geomorphology[ix,iy] != exposed_basement
                aolean_erosion_rates[ix,iy] = world.freeboard[ix,iy] *
                    aolean_erosion_rate_constant * # meters / Myr
                    rho_continent_crust / rho_sediment # m sediment / Myr
            end
        end
    end
    total_aolean_erosion = volume_field(aolean_erosion_rates) 
    mean_aolean_deposition_rate = total_aolean_erosion / ocean_area() # meters / Myr
    # m3 / Myr
    for ix in 1:nx
        for iy in 1:ny
            if world.freeboard[ix,iy] < 0.
                aolean_deposition_rates[ix,iy] = mean_aolean_deposition_rate 
            end
        end
    end
    return aolean_erosion_rates, aolean_deposition_rates
end

function check_reburial_exposed_basement()
    n_buried = 0
    for ix in 1:nx 
        for iy in 1:ny
            if world.geomorphology[ix,iy] == exposed_basement # incoming from plate regrid
                if world.sediment_thickness[ix,iy] > 0.
                    println("why is there mud on my grid point? ",[ix,iy],world.sediment_thickness[ix,iy])
                    world.sediment_thickness[ix,iy] = 0.
                end
                neighbors = get_border_neighbor_coords(ix,iy)
                for nc in neighbors
                    ixn = nc[1]; iyn = nc[2]
                    if world.geomorphology[ixn,iyn] == sedimented_land &&
                        world.freeboard[ixn,iyn] > world.freeboard[ix,iy] 

                        world.tectonics[ix,iy] = burying_land
                        n_buried += 1
                        world.geomorphology[ix,iy] = sedimented_land
                    end
                end
            end
        end
    end
    if n_buried > 0
        println(" burying ",n_buried)
    end
    return n_buried
end

function land_cell_balance(ix,iy)
    println(diffusive_flux_cell(world.freeboard,ix,iy,"land","up"))
    println(diffusive_flux_cell(world.freeboard,ix,iy,"land","left")," ",
        diffusive_flux_cell(world.freeboard,ix,iy,"land","left"))
    println(diffusive_flux_cell(world.freeboard,ix,iy,"land","down"))
end
function diffusive_flux_cell(elevation,ix,iy,domain,direction)
    if domain == "land"
        if direction == "left"
            diffcoeff = get_horizontal_diffcoeff(land_base_diffcoeff,iy)[2]
            neighbor_elevation = elevation[ix-1,iy]
        elseif direction == "right"
            diffcoeff = get_horizontal_diffcoeff(land_base_diffcoeff,iy)[1]
            neighbor_elevation = elevation[ix+1,iy]
        elseif direction == "up"
            diffcoeff = get_vertical_diffcoeff(land_base_diffcoeff,iy)[1]
            neighbor_elevation = elevation[ix,iy+1]
        else  # presumably down
            diffcoeff = get_vertical_diffcoeff(land_base_diffcoeff,iy)[2]
            neighbor_elevation = elevation[ix,iy-1]
        end
    else # presumably seafloor
        if direction == "left"
            diffcoeff = get_horizontal_diffcoeff(land_base_diffcoeff,iy)[2]
            neighbor_elevation = elevation[ix-1,iy]
        elseif direction == "right"
            diffcoeff = get_horizontal_diffcoeff(land_base_diffcoeff,iy)[1]
            neighbor_elevation = elevation[ix+1,iy]
        elseif direction == "up"
            diffcoeff = get_vertical_diffcoeff(land_base_diffcoeff,iy)[1]
            neighbor_elevation = elevation[ix,iy+1]
        else  # presumably down
            diffcoeff = get_vertical_diffcoeff(land_base_diffcoeff,iy)[2]
            neighbor_elevation = elevation[ix,iy-1]
        end
    end
    d_elevation = elevation[ix,iy] - neighbor_elevation
    dhdt = diffcoeff * d_elevation / time_step
    flux = dhdt / sediment_freeboard_expression
    return dhdt
end
function get_vertical_diffcoeff(base_diffcoeff,iy)
    # cell value += D1 (E+ - E0) + D2 (E- - E0)
    diffusive_time_step = time_step * 1.E6 # years
    dx_gradient = delta_y; dx_selfwidth = delta_x[iy]
    if iy == ny # for mixing across the pole, tie things together
        dx_interface = delta_x[iy] / 2
        diffup = base_diffcoeff * # m2 / yr * delta_elevation would give m3 / yr
            diffusive_time_step  /    # m3 
            dx_gradient *         # to get dElev/dx
            dx_interface /        # m3 across the interface
            dx_selfwidth / dx_gradient  # m elevation change in local box
    else
        dx_interface = ( delta_x[iy] + delta_x[iy+1] ) / 2
        diffup = base_diffcoeff * # m2 / yr * delta_elevation would give m3 / yr
            diffusive_time_step  /    # m3 
            dx_gradient *         # to get dElev/dx
            dx_interface /        # m3 across the interface
            dx_selfwidth / dx_gradient  # m elevation change in local box
    end
    if iy == 1
        dx_interface = delta_x[iy] / 2
        diffdown = base_diffcoeff * # m2 / yr * delta_elevation would give m3 / yr
            diffusive_time_step  /    # m3 
            dx_gradient *         # to get dElev/dx
            dx_interface /        # m3 across the interface
            dx_selfwidth / dx_gradient  # m elevation change in local box
    else
        dx_interface = ( delta_x[iy] + delta_x[iy-1] ) / 2
        diffdown = base_diffcoeff * # m2 / yr * delta_elevation would give m3 / yr
            diffusive_time_step  /    # m3 
            dx_gradient *         # to get dElev/dx
            dx_interface /        # m3 across the interface
            dx_selfwidth / dx_gradient  # m elevation change in local box
    end
    return [ diffup, diffdown ]
end
function get_horizontal_diffcoeff(base_diffcoeff,iy)
    dx_width = delta_y; dx_interface = delta_y
    dx_gradient = max( delta_x[iy], delta_y / 2. )
    diffusive_time_step = time_step * 1.E6
    diffcoeff = base_diffcoeff * # m2 / yr * delta_elevation would give m3 / yr
        diffusive_time_step  /    # m3 
        dx_gradient *         # to get dElev/dx
        dx_interface /        # m3 across the interface
        dx_width / dx_gradient  # m elevation change in local box
    return [ diffcoeff, diffcoeff ]
end

    

   
