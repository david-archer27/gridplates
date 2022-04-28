
function distribute_ocean_sediment_fluxes(  )
    # accumulates seafloor_sediment_fraction_deposition_rate, seafloor_sediment_deposition_rate.
    # smoothes the sediment deposition fluxes, tests for sediment overflows,
    # designates the overflow points as nondepositing shelf and moves the 
    # sediment load to adjacent ocean grid points if available, removes the 
    # shelf points from consideration and repeats.  if there is no adjacent
    # ocean grid points the sediment is left in overflow array which gets
    # distributed in a skirt around the enclosing land mass.  
    # check_sediment_overlows applies the non-overlowing sediment into 
    # seafloor_sediment_deposition_rate and seafloor_sediment_fraction_deposition_rate
    # a multi-point basin can fill up by converting all of its grid points to shelf points,
    # excess will be left in overflow.  
    verbose = true
    clay_dep = get_frac_diag("coastal_sediment_fraction_runoff_flux",clay_sediment) .+ 
        get_diag("coastal_orogenic_clay_flux") .+ 
        get_diag("aolean_clay_deposition_rate")

    set_frac_diag( "ocean_sediment_fraction_influx",
        clay_sediment, clay_dep )
    CaCO3_dep = get_diag("coastal_CaCO3_flux") .+  # CaCO3 fluxes
        get_diag("pelagic_CaCO3_deposition_rate") .+ 
        get_diag("continental_CaCO3_deposition_rate") .+ 
        get_frac_diag("coastal_sediment_fraction_runoff_flux",CaCO3_sediment)
    set_frac_diag( "ocean_sediment_fraction_influx",
        CaCO3_sediment, CaCO3_dep )
    sed_fluxes = fill(0.,nx,ny,n_sediment_types+1)
    sed_fluxes[:,:,1] = clay_dep
    sed_fluxes[:,:,2] = CaCO3_dep
    sed_fluxes[:,:,3] = clay_dep .+ CaCO3_dep
    submarine_mask = greater_than_mask_field(world.freeboard .* -1, 0.)

    smooth_sediment_fraction_deposition_rate!( sed_fluxes, submarine_mask )
    # alters sed_fluxes, smoothing the array only
    incoming_flux = volumefield_total( sed_fluxes[:,:,n_sediment_types+1] )

    overflow_fluxes, n_overflows = check_sediment_overflows!( sed_fluxes, submarine_mask )
    # accumulates seafloor_sediment_fraction_deposition_rate, seafloor_sediment_deposition_rate
    # overflow_fluxes are in the same structure as sed_fluxes

    if verbose == true
        deposited_flux = volumefield_total( get_diag("seafloor_sediment_deposition_rate") )
        overflow_flux = volumefield_total( overflow_fluxes[:,:,n_sediment_types+1] )
        println("  ocean in ", incoming_flux," dep ", deposited_flux, " overflow ", 
            overflow_flux, " bal ", incoming_flux-deposited_flux+overflow_flux )
    end
    while n_overflows > 0
        smooth_sediment_fraction_deposition_rate!( overflow_fluxes, submarine_mask )
        #smoothed_overflow_flux = volumefield_total( overflow_fluxes[:,:,n_sediment_types+1] )
        overflow_fluxes, n_overflows = check_sediment_overflows!( overflow_fluxes, submarine_mask )
        if verbose == true
            deposited_flux = volumefield_total( get_diag("seafloor_sediment_deposition_rate") )
            overflow_flux = volumefield_total( overflow_fluxes[:,:,n_sediment_types+1] )
            #println("redo after overflows, now ", n_overflows)
            println("R ocean in ", incoming_flux," dep ", deposited_flux, " overflow ", 
                overflow_flux, " bal ", incoming_flux-deposited_flux+overflow_flux )
        end
    end
end
function smooth_sediment_fraction_deposition_rate!( sed_fluxes, submarine_mask )
    # alters sed_fluxes
    submarine_blobs = get_blobs( submarine_mask )
    #Threads.@threads 
    for i_blob in 1:length(submarine_blobs)
        #println("starting ",i_blob)
        blob = submarine_blobs[i_blob]
        smooth_ocean_area!( sed_fluxes, blob )
    end
end
function smooth_ocean_area!( sed_fluxes, blob )
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
                h_diffcoeffs = get_horizontal_diffcoeff(seafloor_base_diffcoeff,iy)
                v_diffcoeffs = get_vertical_diffcoeff(seafloor_base_diffcoeff,iy)
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

        for i_sedtype in 1:n_sediment_types+1
            #input_fraction_fluxes = 
            #    get_frac_diag("seafloor_sediment_fraction_deposition_rate",i_sedtype)
            R = Float64[]; 
            for i_pos = 1:list_pos_filled
                ix = list_x[i_pos]; iy = list_y[i_pos]
                push!(R, sed_fluxes[ix,iy,i_sedtype])
            end
            new_values_list = lu_A \ R
            for i_pos = 1:list_pos_filled
                ix = list_x[i_pos]; iy = list_y[i_pos]
                sed_fluxes[ix,iy,i_sedtype] = new_values_list[i_pos]
                #set_frac_diag("seafloor_sediment_fraction_deposition_rate",
                #    ix,iy,i_sedtype,new_values_list[i_pos])
                #accum_diag("seafloor_sediment_deposition_rate",
                #    ix,iy,new_values_list[i_pos])
            end
        end
    end # diffuse or no
end
function check_sediment_overflows!( sed_fluxes, submarine_mask )
    # divides the fluxes in sed_fluxes into seafloor_sediment_fraction_deposition_rate
    # where it fits, and a returned analog to sed_fluxes for another round.
    # where sediment doesnt fit it gets partitioned into coastal flux in adjacent
    # deep-water grid points, or summed over a landmass and distributed around the
    # skirt of the landmass for internally trapped ocean points.
    # accumulates seafloor_sediment_fraction_deposition_rate, seafloor_sediment_deposition_rate
    # returns overflow_fluxes and the number of overflows
    verbose = false
    n_overflows = 0; n_trapped = 0
    coastal_overflow_fluxes = fill(0.,nx,ny,n_sediment_types+1)
    trapped_overflow_fluxes = fill(0.,nx,ny,n_sediment_types+1)
    already_depositing = get_diag("seafloor_sediment_deposition_rate")
    potential_freeboard = world.freeboard .+ 
        ( ( sed_fluxes[:,:,n_sediment_types+1] .+ already_depositing ) .* 
        time_step ./ sediment_freeboard_expression ) 
    # zero first time through, incremented on subsequent passes
    incoming_flux_total = volumefield_total(sed_fluxes[:,:,n_sediment_types+1])
    for ix in 1:nx
        for iy in 1:ny
            if submarine_mask[ix,iy] > 0
                if potential_freeboard[ix,iy] < shelf_depth_clastics
                    for i_sedtype in 1:n_sediment_types
                        accum_frac_diag("seafloor_sediment_fraction_deposition_rate",
                            ix,iy,i_sedtype,sed_fluxes[ix,iy,i_sedtype])
                        accum_diag("seafloor_sediment_deposition_rate",ix,iy,
                            sed_fluxes[ix,iy,i_sedtype])        
                    end
                else # overflows
                    world.surface_type[ix,iy] = ocean_shelf
                    submarine_mask[ix,iy] = 0 # exclude from next pass
                    if sed_fluxes[ix,iy,n_sediment_types+1] > 0
                        #println("overflowing ocean point ",ix," ",iy)
                        excess_meters = ( potential_freeboard[ix,iy] - shelf_depth_clastics ) /
                            sediment_freeboard_expression
                        meters_removed = min( excess_meters, 
                            sed_fluxes[ix,iy,n_sediment_types+1] * time_step )
                        fractions = []
                        for i_sedtype in 1:n_sediment_types
                            fraction = sed_fluxes[ix,iy,i_sedtype] /
                                sed_fluxes[ix,iy,n_sediment_types+1]
                            push!(fractions, fraction)
                            sed_fluxes[ix,iy,i_sedtype] -=
                                meters_removed * fraction / time_step
                            accum_frac_diag("seafloor_sediment_fraction_deposition_rate",
                                ix,iy,i_sedtype,sed_fluxes[ix,iy,i_sedtype])
                        end
                        sed_fluxes[ix,iy,n_sediment_types+1] -= meters_removed / time_step
                        accum_diag("seafloor_sediment_deposition_rate",ix,iy,
                            sed_fluxes[ix,iy,n_sediment_types+1])
                        border_neighbor_coords = get_border_neighbor_coords(ix,iy)
                        depositing_border_neighbor_coords = []
                        for pair in border_neighbor_coords
                            ixn = pair[1]; iyn = pair[2]
                            if potential_freeboard[ixn,iyn] < shelf_depth_clastics
                                push!(depositing_border_neighbor_coords, pair)
                            end
                        end
                        if length(depositing_border_neighbor_coords) > 0 # coastal
                            n_overflows += 1 
                            # indicating that the blob will need another pass
                            volume = meters_removed * areabox[iy]
                            for pair in border_neighbor_coords
                                ixn = pair[1]; iyn = pair[2]
                                n_boxes = length(border_neighbor_coords)
                                #println("dumping into ",n_boxes," neighbors")
                                coastal_overflow_fluxes[ix,iy,n_sediment_types+1] += 
                                    volume / n_boxes / areabox[iyn]
                                for i_sedtype in 1:n_sediment_types                                
                                    coastal_overflow_fluxes[ix,iy,i_sedtype] +=
                                        volume / n_boxes / areabox[iyn] * 
                                        fractions[i_sedtype]
                                end
                            end
                        else # trapped
                            #println("trapped")
                            n_trapped += 1
                            for i_sedtype in 1:n_sediment_types 
                                trapped_overflow_fluxes[ix,iy,i_sedtype] = 
                                    meters_removed / time_step * fractions[i_sedtype]
                                trapped_overflow_fluxes[ix,iy,n_sediment_types+1] += 
                                    meters_removed / time_step * fractions[i_sedtype]
                            end
                        end
                    end
                end
            end
        end
    end

    if n_trapped > 0 # prepare for another pass 
        #println(" ocean overflows ",n_overflows," trapped ", n_trapped)
        redistribute_ocean_overflow_fluxes!( coastal_overflow_fluxes, 
            trapped_overflow_fluxes, submarine_mask )
    end

    if verbose == true
        deposited_flux_total = volumefield_total(sed_fluxes[:,:,n_sediment_types+1])
        redeposited_flux_total = volumefield_total(coastal_overflow_fluxes[:,:,n_sediment_types+1])
        trapped_flux_total = volumefield_total(trapped_overflow_fluxes[:,:,n_sediment_types+1])
        new_redeposited_flux_total = volumefield_total(coastal_overflow_fluxes[:,:,n_sediment_types+1])
        println("ocean balance ",incoming_flux_total," ", deposited_flux_total," ",redeposited_flux_total,
            " ",trapped_flux_total," ",new_redeposited_flux_total," ",
            deposited_flux_total + new_redeposited_flux_total)
    end
    return coastal_overflow_fluxes, n_overflows + n_trapped
end
function redistribute_ocean_overflow_fluxes!( coastal_overflow_fluxes, 
    trapped_overflow_fluxes, submarine_mask )
    # move stuff from trapped overflow points to coastal deposition fluxes
    land_mask = 1 .- submarine_mask
    land_blobs = get_blobs(land_mask)
    for i_blob in 1:length(land_blobs)
        blob = land_blobs[i_blob]
        total_volume = 
            volumefield_total( trapped_overflow_fluxes[:,:,n_sediment_types+1] .* blob ) 
        if total_volume > 0
            #println("moving stuff ",i_blob," ",total_volume)
            total_volumes = []
            for i_sedtype in 1:n_sediment_types+1
                push!( total_volumes, 
                    volumefield_total( trapped_overflow_fluxes[:,:,i_sedtype] .* blob ) )
            end
            blob_fringe = get_blob_neighbors( blob )
            area_fringe = volumefield_total( blob_fringe )
            for ix in 1:nx
                for iy in 1:ny
                    if blob_fringe[ix,iy] == 1
                        #println("depositing to ",ix," ",iy)
                        for i_sedtype in 1:n_sediment_types+1
                            coastal_overflow_fluxes[ix,iy,i_sedtype] += 
                                total_volumes[i_sedtype] ./ area_fringe
                        end
                    end
                end
            end
        end
    end
end

