function land_total_sediment_transport( original_elevation_field, diffusive_mask, 
    total_sediment_source, ocean_sink )
    land_sediment_deposition_rates = fill(0.,nx,ny)
    subaereal_blobs = get_blobs( diffusive_mask )
    altered_elevation_field = deepcopy( original_elevation_field )
    verbose = false
    Threads.@threads for i_blob in 1:length(subaereal_blobs) 
        #if verbose == true
        #    println("beginning ",i_blob)
        #end
        subaereal_blob = subaereal_blobs[i_blob]

        blob_sediment_source = volume_field( total_sediment_source .* subaereal_blob )
        
        land_blob_total_sediment_transport!( altered_elevation_field, 
            subaereal_blob, total_sediment_source, ocean_sink )
        for ix in 1:nx
            for iy in 1:ny
                if subaereal_blob[ix,iy] == 1
                    land_sediment_deposition_rates[ix,iy] = 
                        ( altered_elevation_field[ix,iy] - original_elevation_field[ix,iy] ) *
                        ( 1. - world.sediment_porosity[ix,iy] ) /
                        main_time_step / subaereal_sediment_freeboard_expression
                end
            end
        end 

        blob_sediment_dep = volume_field( land_sediment_deposition_rates .* subaereal_blob )
        #println("blob sources ",[i_blob,blob_sediment_source,blob_sediment_dep,
        #    blob_sediment_source - blob_sediment_dep])

    end # blob 
    return altered_elevation_field, land_sediment_deposition_rates
end
function land_blob_total_sediment_transport!( altered_elevation_field, 
    subaereal_blob, total_sediment_source, ocean_sink )
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
        h_diffcoeffs_scaled = get_horizontal_diffcoeff(land_base_diffcoeff,iy) .* 
            subaereal_sediment_freeboard_expression # .* ( 1. - porosity_continental_sediment )
        v_diffcoeffs_scaled = get_vertical_diffcoeff(land_base_diffcoeff,iy) .* 
            subaereal_sediment_freeboard_expression # .* ( 1. - porosity_continental_sediment )
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
            total_sediment_source[ix,iy] /              # meters compacted/Myr
            ( 1. - world.sediment_porosity[ix,iy] ) *   # meters bulk
            main_time_step * # meters / step
            subaereal_sediment_freeboard_expression ]
        new_elevation = A \ R
        altered_elevation_field[ix,iy] = new_elevation
    else # more than one point

        # original recipe, constant influx at orogenic boundary
        # basic equation for influx point.  prime denotes end-of-step state
        # H1' = H1 + D1 ( E2' - E1' ) + Src1
        # adjacent coastal (sink) point
        # H2' = H2 + D2 ( E1' - E2') + D1 ( 0 - E2' )
        # expanded to deal in elevations
        # H = E / L         (L = 0.4 sediment isostatic exposed fraction)
        # E1' = E1 + D1 * L ( E2' - E1' ) + Src1 * L
        # E2' = E2 + D2 * L ( E1' - E2' ) + D1 * L ( 0 - E2' )
        # rearrange, set D = D * L 
        # E1' * [ 1 + D1 ] + E2' * [ - D1 ] = E1 + Src1 * L
        # E1' * [   - D2 ] + E2' * [ 1 + D2 + D1 ] = E2
        # matrix form
        # | 1 + D1, -D1      | | E1' |      | E1 + Src * L |
        # | -D2, 1 + D2 + D1 | | E2' |  =   | E2           |
        # where Src is bulk (including porosity), as are deposition fluxes,
        # which have to be scaled back to compacted for model bookkeeping

        # revised to regulate influx at the orogenic boundary
        # basic equation for influx point.  prime denotes end-of-step state
        # H1' = H1 + D1 ( E2' - E1' ) + ( E1' - E0 )
        # adjacent coastal (sink) point
        # H2' = H2 + D2 ( E1' - E2') + D1 ( 0 - E2' )
        # expanded to deal in elevations
        # H = E / L         (L = 0.4 sediment isostatic exposed fraction)
        # E1' = E1 + D1 * L ( E2' - E1' ) + Src1 * L
        # E2' = E2 + D2 * L ( E1' - E2' ) + D1 * L ( 0 - E2' )
        # rearrange, set D = D * L 
        # E1' * [ 1 + D1 ] + E2' * [ - D1 ] = E1 + Src1 * L
        # E1' * [   - D2 ] + E2' * [ 1 + D2 + D1 ] = E2
        # matrix form
        # | 1 + D1, -D1      | | E1' |      | E1 + Src * L |
        # | -D2, 1 + D2 + D1 | | E2' |  =   | E2           |
        # where Src is bulk (including porosity), as are deposition fluxes,
        # which have to be scaled back to compacted for model bookkeeping

        A = spzeros( list_pos_filled, list_pos_filled )
        R = Float64[]; original_elevation_list = []
        for i_pos = 1:list_pos_filled
            ix = list_x[i_pos]; iy = list_y[i_pos]
                h_diffcoeffs_scaled = get_horizontal_diffcoeff(land_base_diffcoeff,iy) .* 
                    subaereal_sediment_freeboard_expression # .* ( 1. - porosity_continental_sediment )
                v_diffcoeffs_scaled = get_vertical_diffcoeff(land_base_diffcoeff,iy) .* 
                    subaereal_sediment_freeboard_expression # .* ( 1. - porosity_continental_sediment )
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
            R[end] += total_sediment_source[ix,iy] /         # meters compacted / Myr
                ( 1. - world.sediment_porosity[ix,iy] ) *    # meters bulk / Myr
                main_time_step * subaereal_sediment_freeboard_expression # meters
        end

        new_elevation_list = lu(A) \ R

        for i_pos = 1:list_pos_filled
            ix = list_x[i_pos]; iy = list_y[i_pos]
            #sed_depo = ( new_elevation_list[i_pos] - elevation_field[ix,iy] ) /
            #    main_time_step # meters / Myr
            #set_diag("sediment_deposition_rate",ix,iy,sed_depo)
            altered_elevation_field[ix,iy] = new_elevation_list[i_pos]
        end
    end # diffuse or no
end
function check_for_excess_land_deposition_unused!( land_sediment_deposition_rate_field,
    denuded_mask )
    excess_land_deposition = fill(0.,nx,ny)
    potential_elevation = world.freeboard .+
        land_sediment_deposition_rate_field .* main_time_step .*
        subaereal_sediment_freeboard_expression ./ ( 1. .- world.sediment_porosity )
    for ix in 1:nx
        for iy in 1:ny
            for neighbor in get_border_neighbor_coords(ix,iy)
                ixn = neighbor[1]; iyn = neighbor[2]
                if land_sediment_deposition_rate_field[ix,iy] > 0 &&
                    denuded_mask[ixn,iyn] == 1 &&
                    potential_elevation[ix,iy] > potential_elevation[ixn,iyn]
                    # an unphysical tower of mud next to a denuded grid point 
                    #println(ix," ",iy," ",land_sediment_deposition_rate_field[ix,iy])
                    excess_sediment = ( potential_elevation[ix,iy] - 
                        potential_elevation[ixn,iyn] ) / 
                        subaereal_sediment_freeboard_expression *
                        ( 1. - world.sediment_porosity[ix,iy] ) / main_time_step
                    land_sediment_deposition_rate_field[ix,iy] -=
                        excess_sediment 
                    excess_land_deposition[ix,iy] = excess_sediment 
                    potential_elevation[ix,iy] = world.freeboard[ix,iy] +
                        land_sediment_deposition_rate_field[ix,iy] * main_time_step *
                        subaereal_sediment_freeboard_expression /
                        ( 1. - world.sediment_porosity[ix,iy] )
                end
            end
        end
    end
    return excess_land_deposition
end
function check_for_bedrock_or_CaCO3_exposure(
    land_sediment_deposition_rate_field, new_elevation_field, 
    original_diffusing_sediment_thickness)
    # sets world.geomorphology to exclude for next pass,
    # puts the spoil into the nearest denuded_land_boundary_fraction_flux,
    # records it in denuded_sediment_source_fraction_flux
    # called in a loop to calculate the total fluxes, but the fractions denuded
    # have to be preserved 
    n_denuded = 0
    n_restored = 0
    sediment_deposition_meters =
        land_sediment_deposition_rate_field ./
        (1.0 .- world.sediment_porosity ) .*
        main_time_step
    new_total_sediment_thickness = original_diffusing_sediment_thickness .+
        sediment_deposition_meters

    original_geomorphology = deepcopy(world.geomorphology)

    #=for ix in 1:nx # check to see if the soil should advance onto the rock
        for iy in 1:ny
            for neighbor in get_border_neighbor_coords(ix, iy)
                ixn = neighbor[1]; iyn = neighbor[2]
                if original_geomorphology[ix, iy] == sedimented_land &&
                   original_geomorphology[ixn, iyn] == exposed_basement_or_CaCO3 &&
                   new_elevation_field[ix, iy] > new_elevation_field[ixn, iyn]
                    # an unphysical pile of mud next higher than a denuded grid point 
                    world.geomorphology[ixn, iyn] = sedimented_land
                    n_restored += 1
                    #println(ix," ",iy," ",land_sediment_deposition_rate_field[ix,iy])
                    #=excess_sediment = ( potential_elevation[ix,iy] - 
                        potential_elevation[ixn,iyn] ) ./ 
                        subaereal_sediment_freeboard_expression .*
                        ( 1. - porosity_continental_sediment ) / main_time_step
                    land_sediment_deposition_rate_field[ix,iy] -=
                        excess_sediment 
                    excess_land_deposition[ix,iy] = excess_sediment 
                    potential_elevation[ix,iy] = world.freeboard[ix,iy] .+
                        land_sediment_deposition_rate_field[ix,iy] .* main_time_step .*
                        subaereal_sediment_freeboard_expression ./
                        ( 1. - porosity_continental_sediment )=#
                end
            end
        end
    end=#
    for ix in 1:nx
        for iy in 1:ny
            if original_geomorphology[ix, iy] == sedimented_land &&
                new_total_sediment_thickness[ix, iy] < 0.0
                #println("found ",[ix,iy])
                n_denuded += 1
                world.geomorphology[ix, iy] = exposed_basement
            end
        end
    end
 
    #println("denuded, restored ", [n_denuded, n_restored])
    return n_denuded + n_restored
end
function land_sediment_fraction_transport(new_elevation_field,
    land_sediment_deposition_rate_field, original_diffusing_sediment_thickness,
    diffusive_mask,
    sediment_source_fields, ocean_sink)
    # sets land_sediment_fraction_transport_deposition_rate on land areas
    combined_land_sediment_fraction_transport_deposition_rates =
        fill(0.0, nx, ny, 0:n_sediment_types)
    sediment_deposition_meters =
        land_sediment_deposition_rate_field ./
        ( 1. .- world.sediment_porosity ) .* main_time_step
    new_total_sediment_thickness = original_diffusing_sediment_thickness .+
                                   sediment_deposition_meters
    covered_diffusive_mask = diffusive_mask .- le_mask(new_total_sediment_thickness, 0.0)
    subaereal_blobs = get_blobs(covered_diffusive_mask)

    Threads.@threads for i_blob in 1:length(subaereal_blobs)
        #println("beginning ",i_blob)
        subaereal_blob = subaereal_blobs[i_blob]
        land_sediment_fraction_transport_deposition_rates =
            land_area_sediment_fraction_transport(
                new_elevation_field, land_sediment_deposition_rate_field,
                original_diffusing_sediment_thickness,
                subaereal_blob, sediment_source_fields, ocean_sink)
        for ix in 1:nx
            for iy in 1:ny
                if subaereal_blob[ix, iy] == 1
                    for i_sedtype in first_land_transported_sediment:n_sediment_types
                        if land_sediment_fraction_transport_deposition_rates[ix, iy, i_sedtype] != 0.0
                            combined_land_sediment_fraction_transport_deposition_rates[ix, iy, i_sedtype] =
                                land_sediment_fraction_transport_deposition_rates[ix, iy, i_sedtype]
                        end
                    end
                end
            end
        end
    end
    #=
    for ix in 1:nx
        for iy in 1:ny
            if world.geomorphology[ix,iy] == exposed_basement
                for i_sedtype in 1:n_sediment_types
                    combined_land_sediment_fraction_transport_deposition_rates[ix,iy,i_sedtype] =
                        - world.sediment_thickness[ix,iy] * 
                        world.sediment_fractions[ix,iy,i_sedtype] / main_time_step
                end
            end
        end
    end =#
    update_flux_totals!(combined_land_sediment_fraction_transport_deposition_rates)
    return combined_land_sediment_fraction_transport_deposition_rates
end
function land_area_sediment_fraction_transport( 
    new_elevation_field, land_sediment_deposition_rate_field,
    original_diffusing_sediment_thickness, 
    subaereal_blob, sediment_source_fields, ocean_sink ) 
    # sets land_sediment_fraction_transport_deposition_rate on land points
    # fractions combined in one call so they can share the 
    # transportation matrix.  probably not too expensive to add more components?

    # basic equation for influx point.  prime denotes end-of-step state.

    #    box 1  |  box 2  |  ocean 
    #     Inv1  |   Inv 2 |     
    #          <D1
    #           D2> 

    # 
    # Inv1' = Inv1 + D1 ( 1-por ) ( E2' - E1' ) * ( C1' + C2' ) / 2 + Src1
    #                -- compacted mass flux ---   ---- avg conc ---   compacted
    # adjacent coastal (sink) point
    # Inv2' = Inv2 + D2 ( 1-por ) ( E1' - E2') * (C1' + C2')/2 - D1 ( 1-por ) E2' C2'

    # rearrange to deal in inventories; Hi is portion of H due to inv of one component 
    # Hi = Inv / ( 1-por )   
    # Inv = Hi ( 1-por )     ;   sum Inv = H ( 1-por )
    # C = Inv / ( H ( 1-por ) )
    # Inv1' = Inv1 + D1 ( 1-por ) ( E2' - E1' ) ( Inv1'/(( 1-por )H1') + Inv2'/(( 1-por )H2') ) / 2 + Src1
    # Inv2' = Inv2 + D2 ( 1-por ) ( E1' - E2' ) ( Inv1'/(( 1-por )H1') + Inv2'/(( 1-por )H2') ) / 2 - D1 ( 1-por ) E2' Inv2' / ((1-por)H2')

    # Inv1' = Inv1 + D1 ( E2' - E1' ) ( Inv1'/H1' + Inv2'/H2' ) / 2 + Src1
    # Inv2' = Inv2 + D2 ( E1' - E2' ) ( Inv1'/H1' + Inv2'/H2' ) / 2 - D1 E2' Inv2' / H2'

    # group unknown end-step Inv terms 
    # Inv1' [1 - D1(E2'-E1')/(2 H1')] + Inv2' [  - D1(E2'-E1')/(2 H2')                ]          = Inv1 + Src1
    # Inv1' [  - D2(E1'-E2')/(2 H1')] + Inv2' [1 - D2(E1'-E2')/(2 H2') + D1 E2' / H2' ] = Inv2

    # matrix form
    #| 1 - D1 (E2'-E1') / (2 H1'),   - D1 (E2'-E1') / (2 H2')                | | Inv1' |      | Inv1 + Src |
    #|   - D2 (E1'-E2') / (2 H1'), 1 - D2 (E1'-E2') / (2 H2') + D1 E2' / H2' | | Inv2' |  =   | Inv2       |

    land_sediment_fraction_transport_deposition_rates = fill(0.,nx,ny,n_sediment_types)
    sediment_deposition_meters = 
        land_sediment_deposition_rate_field ./ 
        ( 1. .- world.sediment_porosity ) .* main_time_step

    new_total_sediment_thickness = original_diffusing_sediment_thickness .+
        sediment_deposition_meters
    #old_column_fraction_inventory = fill(0.,nx,ny,n_sediment_types)
    #for i_sedtype in first_land_transported_sediment:n_sediment_types
    #    old_column_fraction_inventory[:,:,i_sedtype] = 
    #        world.sediment_thickness .*
    #        world.sediment_fractions[:,:,i_sedtype]
    #end

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
            h_diffcoeffs = get_horizontal_diffcoeff(land_base_diffcoeff,iy) # .* ( 1. - porosity_continental_sediment )^2
            v_diffcoeffs = get_vertical_diffcoeff(land_base_diffcoeff,iy) # .* ( 1. - porosity_continental_sediment )^2
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
            for i_sedtype in first_land_transported_sediment:n_sediment_types

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
                old_inventory = world.sediment_inventories[ix,iy,i_sedtype]
                R = world.sediment_inventories[ix,iy,i_sedtype] + 
                    sediment_source_fields[ix,iy,i_sedtype] * 
                    main_time_step
                new_inventory = A \ R   
                #new_fraction = new_inventory / 
                #    new_total_sediment_thickness[ix,iy]
                fraction_deposition_rate = 
                    ( new_inventory - old_inventory ) / 
                    main_time_step
                land_sediment_fraction_transport_deposition_rates[ix,iy,i_sedtype] =
                    fraction_deposition_rate
                    #set_frac_diag("land_sediment_fraction_transport_deposition_rate",
                #    ix,iy,i_sedtype,fraction_deposition_rate)
            end
        end
    elseif list_pos_filled > 1# more than one point
        A = spzeros( list_pos_filled, list_pos_filled )
        for i_pos = 1:list_pos_filled
            ix = list_x[i_pos]; iy = list_y[i_pos]
            h_diffcoeffs = get_horizontal_diffcoeff(land_base_diffcoeff,iy) #.* ( 1. - porosity_continental_sediment)^2
            # no isostatic rebound factor here -- getting D dt / dx, gives thickness transport
            v_diffcoeffs = get_vertical_diffcoeff(land_base_diffcoeff,iy) #.* ( 1. - porosity_continental_sediment)^2
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

        for i_sedtype in first_land_transported_sediment:n_sediment_types 
            old_inventory_list = []
            R = Float64[]; 
            for i_pos = 1:list_pos_filled
                ix = list_x[i_pos]; iy = list_y[i_pos]

                #=old_inventory = world.sediment_thickness[ix,iy,current_time_bin()] * 
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
                old_inventory = world.sediment_inventories[ix,iy,i_sedtype]
                push!(old_inventory_list,old_inventory)
                push!(R, old_inventory + 
                    sediment_source_fields[ix,iy,i_sedtype] * main_time_step )
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
                    main_time_step
                land_sediment_fraction_transport_deposition_rates[ix,iy,i_sedtype] =
                    fraction_deposition_rate
            end
        end

        # opportunistically update land sediment weathering age field 
        old_age_list = []
        R = Float64[]
        for i_pos = 1:list_pos_filled
            ix = list_x[i_pos]; iy = list_y[i_pos]

            old_age = world.sediment_weathering_age[ix,iy]
            push!(old_age_list,old_age)
            push!(R, old_age + main_time_step )
        end
        new_age_list = lu_A \ R
        #new_fraction_field = fill(0.,nx,ny)
        for i_pos = 1:list_pos_filled
            if new_age_list[i_pos] != new_age_list[i_pos]
                error("NaN detected in land_area_sediment_fraction_transport")
            end
            ix = list_x[i_pos]; iy = list_y[i_pos]
            #new_fraction = new_inventory_list[i_pos] / 
            #    new_total_layer_sediment_thickness[ix,iy]
            #new_fraction_field[ix,iy] = new_fraction
            world.sediment_weathering_age[ix,iy] = new_age_list[i_pos]
        end
    end 
    return land_sediment_fraction_transport_deposition_rates
end

function land_total_runoff_fluxes(new_elevation_field, subaereal_mask, ocean_sink)
    # resets and fills coastal_sediment_fraction_runoff_flux
    # sets world.geomorphology for coastal_depocenter
    #reset_frac_diag("coastal_sediment_fraction_runoff_flux")
    coastal_sediment_runoff_flux = fill(0.,nx,ny)
    
    # we want to get a preview without altering the real world 

    # basic equation for adjacent coastal (sink) point
    # H2' = H2 + D2 ( E1' - E2') + D1 ( 0 - E2' ) 
    # term D1 E2 is change in H. 
    # compressed volume flux = D1 ( 1-por ) E2 area2
    # adjacent deposition flux = D1 ( 1-por ) E2 area2 / area3

    for ix in 1:nx
        for iy in 1:ny
            if subaereal_mask[ix, iy] > 0
                h_diffcoeffs = get_horizontal_diffcoeff(land_base_diffcoeff, iy) .* 
                    (1.0 - world.sediment_porosity[ix,iy] )
                v_diffcoeffs = get_vertical_diffcoeff(land_base_diffcoeff, iy) .* 
                    (1.0 - world.sediment_porosity[ix,iy] )
                if ocean_sink[ix_left(ix), iy] > 0.0
                    volume_flux = new_elevation_field[ix, iy] * h_diffcoeffs[2] * # meters / step
                        areabox[iy] / main_time_step
                    coastal_sediment_runoff_flux[ix_left(ix), iy] +=
                        volume_flux / areabox[iy]
                    world.geomorphology[ix_left(ix),iy] = coastal_depocenter
                end
                if ocean_sink[ix_right(ix), iy] > 0.0
                    volume_flux = new_elevation_field[ix, iy] * h_diffcoeffs[1] *
                        areabox[iy] / main_time_step
                    coastal_sediment_runoff_flux[ix_right(ix), iy] +=
                        volume_flux / areabox[iy]
                    world.geomorphology[ix_right(ix),iy] = coastal_depocenter
                end
                if iy > 1
                    if ocean_sink[ix, iy-1] > 0
                        volume_flux = new_elevation_field[ix, iy] * v_diffcoeffs[2] *
                            areabox[iy] / main_time_step
                        coastal_sediment_runoff_flux[ix, iy-1] +=
                            volume_flux / areabox[iy-1]
                        world.geomorphology[ix,iy-1] = coastal_depocenter
                    end
                end
                if iy < ny
                    if ocean_sink[ix, iy+1] > 0
                        volume_flux = new_elevation_field[ix, iy] * v_diffcoeffs[1] *
                            areabox[iy] / main_time_step
                        coastal_sediment_runoff_flux[ix, iy+1] +=
                            volume_flux / areabox[iy+1]
                        world.geomorphology[ix,iy+1] = coastal_depocenter
                    end
                end
            end
        end
    end
    return coastal_sediment_runoff_flux
end
function land_fraction_runoff_fluxes( new_elevation_field, 
    land_sediment_fraction_transport_deposition_rate_fields,
    diffusive_mask, ocean_sink )    

    new_sediment_thickness, new_sediment_inventories = 
        apply_sediment_fluxes( 
            land_sediment_fraction_transport_deposition_rate_fields ) 

    # creating a dummy set of sediment fractions that are transported
    # to unmix them from the untransported fraction (CaCO3)
    transported_sediment_fractions = fill(0.,nx,ny,first_land_transported_sediment:n_sediment_types)
    for ix in 1:nx
        for iy in 1:ny
            total_sediment_inventory = 0.
            for i_sedtype in first_land_transported_sediment:n_sediment_types
                total_sediment_inventory += new_sediment_inventories[ix,iy,i_sedtype]
            end
            if total_sediment_inventory > 0
                for i_sedtype in first_land_transported_sediment:n_sediment_types
                    transported_sediment_fractions[ix, iy, i_sedtype] =
                        new_sediment_inventories[ix, iy, i_sedtype] /
                        total_sediment_inventory
                end
            end 
        end
    end

    coastal_sediment_fraction_runoff_flux = 
        fill(0.,nx,ny,0:n_sediment_types)

    # basic equation for adjacent coastal (sink) point
    # H2' = H2 + D2 ( E1' - E2') + D1 ( 0 - E2' ) 
    # term D1 E2 is change in H
    # compressed volume flux = D1 E2 ( 1-por ) area2
    # adjacent deposition flux = D1 E2 ( 1-por ) area2 / area3

    for ix in 1:nx
        for iy in 1:ny
            if diffusive_mask[ix,iy] > 0
                h_diffcoeffs = get_horizontal_diffcoeff(land_base_diffcoeff,iy) .* 
                    ( 1. - world.sediment_porosity[ix,iy] )
                v_diffcoeffs = get_vertical_diffcoeff(land_base_diffcoeff,iy) .* 
                    ( 1. - world.sediment_porosity[ix,iy] )
                if ocean_sink[ix_left(ix),iy] > 0.
                    #println(ix," ",iy); return
                    for i_sedtype in first_land_transported_sediment:n_sediment_types
                        volume_flux = new_elevation_field[ix,iy] * h_diffcoeffs[2] * # meters / step
                            transported_sediment_fractions[ix,iy,i_sedtype] * 
                            areabox[iy] / main_time_step  # m3 / Myr
                        coastal_sediment_fraction_runoff_flux[ix_left(ix),iy,i_sedtype] +=
                            volume_flux / areabox[iy]
                        #accum_frac_diag( "coastal_sediment_fraction_runoff_flux",
                        #    ix_left(ix),iy,i_sedtype,
                        #    volume_flux / areabox[iy] )
                    end
                    world.geomorphology[ix_left(ix),iy] = coastal_depocenter
                end
                if ocean_sink[ix_right(ix),iy] > 0.
                    for i_sedtype in first_land_transported_sediment:n_sediment_types
                        volume_flux = new_elevation_field[ix,iy] * h_diffcoeffs[1] * 
                            transported_sediment_fractions[ix,iy,i_sedtype] * 
                            areabox[iy] / main_time_step
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
                        for i_sedtype in first_land_transported_sediment:n_sediment_types
                            volume_flux = new_elevation_field[ix,iy] * v_diffcoeffs[2] * 
                                transported_sediment_fractions[ix,iy,i_sedtype] *
                                areabox[iy] / main_time_step
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
                        for i_sedtype in first_land_transported_sediment:n_sediment_types
                            volume_flux = new_elevation_field[ix,iy] * v_diffcoeffs[1] * 
                                transported_sediment_fractions[ix,iy,i_sedtype] *
                                areabox[iy] / main_time_step
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
    update_flux_totals!(coastal_sediment_fraction_runoff_flux)
    return coastal_sediment_fraction_runoff_flux
end
function aolean_transport()
    # resets and fills aolean_deposition_rate and aolean_deposition_rate
    aolean_erosion_rates = fill(0.,nx,ny,0:n_sediment_types)
    aolean_deposition_rates = fill(0.,nx,ny,0:n_sediment_types)
    for ix in 1:nx
        for iy in 1:ny
            if world.freeboard[ix,iy] > 0. &&
                world.geomorphology[ix,iy] == sedimented_land
                for i_sedtype in 1:n_sediment_types
                    aolean_erosion_rates[ix,iy,i_sedtype] = world.freeboard[ix,iy] *
                        world.sediment_fractions[ix,iy,i_sedtype] *
                        aolean_erosion_rate_constant * # meters / Myr
                        rho_continent_crust / rho_compacted_sediment # m sediment / Myr
                    aolean_erosion_rates[ix,iy,i_sedtype] = min( 
                        aolean_erosion_rates[ix,iy,i_sedtype],
                        world.sediment_thickness[ix,iy] *
                        world.sediment_fractions[ix,iy,i_sedtype] /
                        main_time_step / 10. )
                    #aolean_erosion_rates[ix,iy,0] += aolean_erosion_rates[ix,iy,i_sedtype]
                end
            end
        end
    end
    total_aolean_erosion_rates = volume_fields(aolean_erosion_rates[:,:,1:end]) 
    mean_aolean_deposition_rates = total_aolean_erosion_rates ./ ocean_area() # meters / Myr
    # m3 / Myr
    for ix in 1:nx
        for iy in 1:ny
            if world.freeboard[ix,iy] < 0.
                for i_sedtype in 1:n_sediment_types
                    aolean_deposition_rates[ix,iy,i_sedtype] = 
                        mean_aolean_deposition_rates[i_sedtype] 
                    #aolean_deposition_rates[ix,iy,0] += 
                    #    aolean_deposition_rates[ix,iy,i_sedtype]
                end
            end
        end
    end
    update_flux_totals!(aolean_erosion_rates)
    update_flux_totals!(aolean_deposition_rates)
    return aolean_erosion_rates, aolean_deposition_rates
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
    dhdt = diffcoeff * d_elevation / main_time_step
    flux = dhdt / subaereal_sediment_freeboard_expression
    return dhdt
end
function get_vertical_diffcoeff(base_diffcoeff,iy)
    # cell value += D1 (E+ - E0) + D2 (E- - E0)
    #diffusive_main_time_step = main_time_step * 1.E6 # years
    dx_gradient = delta_y; dx_selfwidth = delta_x[iy]
    if iy == ny # for mixing across the pole, tie things together
        dx_interface = delta_x[iy] / 2
        diffup = base_diffcoeff * # m2 / Myr * delta_elevation would give m3 / Myr
            main_time_step  /    # m3 
            dx_gradient *         # to get dElev/dx
            dx_interface /        # m3 across the interface
            dx_selfwidth / dx_gradient  # m elevation change in local box
    else
        dx_interface = ( delta_x[iy] + delta_x[iy+1] ) / 2
        diffup = base_diffcoeff * # m2 / yr * delta_elevation would give m3 / yr
            main_time_step  /    # m3 
            dx_gradient *         # to get dElev/dx
            dx_interface /        # m3 across the interface
            dx_selfwidth / dx_gradient  # m elevation change in local box
    end
    if iy == 1
        dx_interface = delta_x[iy] / 2
        diffdown = base_diffcoeff * # m2 / yr * delta_elevation would give m3 / yr
            main_time_step  /    # m3 
            dx_gradient *         # to get dElev/dx
            dx_interface /        # m3 across the interface
            dx_selfwidth / dx_gradient  # m elevation change in local box
    else
        dx_interface = ( delta_x[iy] + delta_x[iy-1] ) / 2
        diffdown = base_diffcoeff * # m2 / yr * delta_elevation would give m3 / yr
            main_time_step  /    # m3 
            dx_gradient *         # to get dElev/dx
            dx_interface /        # m3 across the interface
            dx_selfwidth / dx_gradient  # m elevation change in local box
    end
    return [ diffup, diffdown ]
end
function get_horizontal_diffcoeff(base_diffcoeff,iy)
    dx_width = delta_y; dx_interface = delta_y
    dx_gradient = max( delta_x[iy], delta_y / 2. )
    #diffusive_main_time_step = main_time_step * 1.E6
    diffcoeff = base_diffcoeff * # m2 / Myr * delta_elevation would give m3 / yr
        main_time_step  /    # m3 
        dx_gradient *         # to get dElev/dx
        dx_interface /        # m3 across the interface
        dx_width / dx_gradient  # m elevation change in local box
    return [ diffcoeff, diffcoeff ]
end

    

   
