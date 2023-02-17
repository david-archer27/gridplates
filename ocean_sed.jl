
function distribute_ocean_sediment_fluxes( incoming_fluxes, submarine_mask )
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

    accumulating_depositing_fluxes = fill(0.,nx,ny,0:n_sediment_types)
 
    #submarine_mask = lt_mask(world.freeboard, 0.) # eq_mask(world.crust_type,ocean_crust)

    #println(" influxes ", [ volume_field(incoming_fluxes[:,:,1]), 
    #    volume_field(incoming_fluxes[:,:,2]), volume_field(incoming_fluxes[:,:,3]),
    #    volume_field(incoming_fluxes[:,:,1]) + volume_field(incoming_fluxes[:,:,2]) ])

    smooth_sediment_fraction_deposition_rate!( incoming_fluxes, submarine_mask )
    # alters incoming_fluxes, smoothing the array only

    #println(" smoothed ", [ volume_field(incoming_fluxes[:,:,1]), 
    #    volume_field(incoming_fluxes[:,:,2]), volume_field(incoming_fluxes[:,:,3]),
    #    volume_field(incoming_fluxes[:,:,1]) + volume_field(incoming_fluxes[:,:,2]) ])

    overflow_fluxes, depositing_fluxes, new_submarine_mask, n_overflows = 
        check_sediment_overflows!( incoming_fluxes, accumulating_depositing_fluxes, 
            submarine_mask )

    #println(" depositing ", [ volume_field(depositing_fluxes[:,:,1]), 
    #    volume_field(depositing_fluxes[:,:,2]), volume_field(depositing_fluxes[:,:,3]),
    #    volume_field(depositing_fluxes[:,:,1]) + volume_field(depositing_fluxes[:,:,2]) ])
        
    accumulating_depositing_fluxes = deepcopy(depositing_fluxes)
    overflow_flux = #volume_field( overflow_fluxes[:,:,1] ) +
        volume_field( overflow_fluxes[:,:,0] )
    if enable_distribute_ocean_fluxes_diagnostics
        incoming_flux = #volume_field( incoming_fluxes[:,:,1] ) + 
            volume_field( incoming_fluxes[:,:,0] )
        deposited_flux = #volume_field( accumulating_depositing_fluxes[:,:,1] ) + 
            volume_field( accumulating_depositing_fluxes[:,:,0] )
        balance = incoming_flux - deposited_flux - overflow_flux
        println("    ocean in ", incoming_flux," dep ", deposited_flux, " overflow ", 
            overflow_flux, " bal ", balance )
    else
        println("iterating ocean deposition ")
    end
    while n_overflows > 0 || overflow_flux > 0.
        #if n_overflows < 1000
        #    error("stopping")
        #end
        submarine_mask = deepcopy(new_submarine_mask)
        #println("smooth in ",volume_field(overflow_fluxes[:,:,1]))
        smooth_sediment_fraction_deposition_rate!( overflow_fluxes, submarine_mask )
        #println("smooth out ",volume_field(overflow_fluxes[:,:,1]))
        #println("oveflow ", overflow_fluxes[187,41,:])
        new_overflow_fluxes, new_depositing_fluxes, new_submarine_mask, n_overflows = 
            check_sediment_overflows!( overflow_fluxes, accumulating_depositing_fluxes, 
                submarine_mask )
        #println("new over ",volume_field(new_overflow_fluxes[:,:,1]))
        #println("new dep ",volume_field(new_depositing_fluxes[:,:,1]))
        accumulating_depositing_fluxes .+= new_depositing_fluxes
        #overflow_fluxes = deepcopy(new_overflow_fluxes)
        if enable_distribute_ocean_fluxes_diagnostics
            incoming_flux = #volume_field( incoming_fluxes[:,:,1] ) + 
                volume_field( incoming_fluxes[:,:,0] )
            deposited_flux = #volume_field( accumulating_depositing_fluxes[:,:,1] ) + 
                volume_field( accumulating_depositing_fluxes[:,:,0] )
            overflow_flux = #volume_field( new_overflow_fluxes[:,:,1] ) +
                volume_field( new_overflow_fluxes[:,:,0] )
            balance = incoming_flux-deposited_flux-overflow_flux
            println("  R ocean in ", incoming_flux," dep ", deposited_flux, " overflow ", 
                overflow_flux, " bal ", balance, " ",n_overflows )
            #=println(" influxes ", [ volume_field(incoming_fluxes[:,:,1]), 
                volume_field(incoming_fluxes[:,:,2]), volume_field(incoming_fluxes[:,:,3]),
                volume_field(incoming_fluxes[:,:,1]) + volume_field(incoming_fluxes[:,:,2]) ])
                println(" accum_dep ", [ volume_field(new_depositing_fluxes[:,:,1]), 
                volume_field(new_depositing_fluxes[:,:,2]), volume_field(new_depositing_fluxes[:,:,3]),
                volume_field(new_depositing_fluxes[:,:,1]) + volume_field(new_depositing_fluxes[:,:,2]) ])
                println(" overflow ", [ volume_field(new_overflow_fluxes[:,:,1]), 
                volume_field(new_overflow_fluxes[:,:,2]), volume_field(new_overflow_fluxes[:,:,3]),
                volume_field(new_overflow_fluxes[:,:,1]) + volume_field(new_overflow_fluxes[:,:,2]) ])=#
        end
        if enable_watch_ocean_offshore_transport
            scene = plot_field( eq_mask( new_submarine_mask ))
            plot_add_continent_outlines!( scene )
            depositing = gt_mask( accumulating_depositing_fluxes[:,:,0], 0. )
            Makie.contour!(scene,xcoords,ycoords,depositing,color=:yellow)
 
            display(scene)
            sleep(2)
        end
        overflow_fluxes = deepcopy(new_overflow_fluxes)
    end 
    return accumulating_depositing_fluxes

end
function smooth_sediment_fraction_deposition_rate!( sed_fluxes, submarine_mask )
    # alters sed_fluxes
    submarine_blobs = get_blobs( submarine_mask )
    Threads.@threads for i_blob in 1:length(submarine_blobs)
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

        for i_sedtype in 0:n_sediment_types
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
            end
        end
    end # diffuse or no
end
function check_sediment_overflows!( incoming_fluxes, already_depositing_fluxes, submarine_mask )
    # returns coastal_overflow_fluxes, depositing_fluxes, new_submarine_mask, n_overflows + n_trapped

    # divides the fluxes in sed_fluxes into seafloor_sediment_fraction_deposition_rate
    # where it fits, and a returned analog to sed_fluxes for another round.
    # where sediment doesnt fit it gets partitioned into coastal flux in adjacent
    # deep-water grid points, or summed over a landmass and distributed around the
    # skirt of the landmass for internally trapped ocean points.
    # accumulates seafloor_sediment_fraction_deposition_rate, seafloor_sediment_deposition_rate
    # returns overflow_fluxes, depositing_fluxes and the number of overflows
    local firstpass_coastal_flux_total, trapped_flux_total
    verbose = true
    n_overflows = 0; n_trapped = 0
    depositing_fluxes = fill(0.,nx,ny,0:n_sediment_types)
    coastal_overflow_fluxes = fill(0.,nx,ny,0:n_sediment_types)
    trapped_overflow_fluxes = fill(0.,nx,ny,0:n_sediment_types)
    new_submarine_mask = deepcopy(submarine_mask)
    potential_freeboard = world.freeboard .+ 
        ( incoming_fluxes[:,:,0] .+ 
        already_depositing_fluxes[:,:,0] ) .* 
        main_time_step .* submarine_sediment_freeboard_expression 
    # 
    for ix in 1:nx
        for iy in 1:ny
            if submarine_mask[ix,iy] == 0
                trapped_overflow_fluxes[ix,iy,:] = incoming_fluxes[ix,iy,:]
                #n_trapped += 1
            else
                if potential_freeboard[ix,iy] < shelf_depth_clastics # deeper, accumulating
                    depositing_fluxes[ix,iy,:] = incoming_fluxes[ix,iy,:]
                else # overflows
                    n_trapped += 1
                    new_submarine_mask[ix,iy] = 0 # exclude from next pass
                    if incoming_fluxes[ix,iy,0] > 0
                        #= println("overflowing ocean point ",ix," ",iy)
                        println("before ",volume_field(incoming_fluxes[:,:,0])," ",
                            volume_field(depositing_fluxes[:,:,0]), " ",
                            volume_field(coastal_overflow_fluxes[:,:,0]), " ",
                            volume_field(trapped_overflow_fluxes[:,:,0]), " ",
                            volume_field(incoming_fluxes[:,:,0]) +
                                volume_field(depositing_fluxes[:,:,0])+volume_field(coastal_overflow_fluxes[:,:,0]) +
                                volume_field(trapped_overflow_fluxes[:,:,0]))
                                =#
                        incoming_meters = incoming_fluxes[ix,iy,0] * main_time_step
                        excess_meters = ( potential_freeboard[ix,iy] - shelf_depth_clastics ) /
                            submarine_sediment_freeboard_expression
                        meters_exported = min( excess_meters, incoming_meters )
                        meters_deposited_locally = incoming_meters - meters_exported
                        fractions_exported = []; sum_check = 0.
                        for i_sedtype in 1:n_sediment_types
                            fraction = incoming_fluxes[ix,iy,i_sedtype] /
                                incoming_fluxes[ix,iy,0]
                            sum_check += fraction
                            push!(fractions_exported, fraction)
                            depositing_fluxes[ix,iy,i_sedtype] +=
                                meters_deposited_locally * fraction / main_time_step
                        end
                        #incoming_fluxes[ix,iy,0] -= meters_exported / main_time_step
                        depositing_fluxes[ix,iy,0] +=
                            meters_deposited_locally / main_time_step

                        if sum_check < 0.99 || sum_check > 1.01
                            error("sum_check error ",[ix,iy],incoming_fluxes[ix,iy,:])
                        end

                        # find eligible neighboring cells, if any
                        border_neighbor_coords = get_border_neighbor_coords(ix,iy)
                        depositing_border_neighbor_coords = []
                        for pair in border_neighbor_coords
                            ixn = pair[1]; iyn = pair[2]
                            if potential_freeboard[ixn,iyn] < shelf_depth_clastics
                                push!(depositing_border_neighbor_coords, pair)
                            end
                        end
                        if length(depositing_border_neighbor_coords) > 0 # coastal
                            #println("  coastal ",[ix,iy])
                            world.geomorphology[ix,iy] = ocean_shelf
                            n_overflows += 1 
                            # triggering another pass
                            volume_exported = meters_exported * areabox[iy]
                            area_depositing = 0.
                            for pair in depositing_border_neighbor_coords
                                ixn = pair[1]; iyn = pair[2]
                                area_depositing += areabox[iyn]
                            end
                            mean_meters_depositing = volume_exported / area_depositing
                            for pair in depositing_border_neighbor_coords
                                ixn = pair[1]; iyn = pair[2]
                                n_boxes = length(depositing_border_neighbor_coords)
                                #println("dumping into ",n_boxes," neighbors")
                                coastal_overflow_fluxes[ixn,iyn,0] += 
                                    mean_meters_depositing / main_time_step
                                for i_sedtype in 1:n_sediment_types                                
                                    coastal_overflow_fluxes[ixn,iyn,i_sedtype] +=
                                        mean_meters_depositing / main_time_step * 
                                        fractions_exported[i_sedtype]
                                end
                                sum_check = 0.
                                for i_sedtype in 1:n_sediment_types
                                    sum_check += fractions_exported[i_sedtype]
                                end
                                if sum_check < 0.999 || sum_check > 1.001
                                    error("sum_check error 2 ",[ix,iy,fractions_exported,sum_check])
                                end
                            end
                        else # trapped, no adjacent cells
                            #println("  trapped ",[ix,iy])
                            world.geomorphology[ix,iy] = filled_basin
                            n_trapped += 1
                            for i_sedtype in 1:n_sediment_types 
                                trapped_overflow_fluxes[ix,iy,i_sedtype] += 
                                    meters_exported / main_time_step * 
                                    fractions_exported[i_sedtype]
                            end
                            trapped_overflow_fluxes[ix,iy,0] += 
                                meters_exported / main_time_step

                        end
                        #=println("after ",volume_field(incoming_fluxes[:,:,0])," ",
                        volume_field(depositing_fluxes[:,:,0]), " ",
                        volume_field(coastal_overflow_fluxes[:,:,0]), " ",
                        volume_field(trapped_overflow_fluxes[:,:,0]), " ",
                        volume_field(incoming_fluxes[:,:,0]) + volume_field(depositing_fluxes[:,:,0]) +
                            volume_field(coastal_overflow_fluxes[:,:,0]) + volume_field(trapped_overflow_fluxes[:,:,0]))
                        =#
                    end
                end
            end
        end
    end

    if verbose == true
        incoming_flux_total = volume_field(incoming_fluxes[:,:,0])
        depositing_flux_total = volume_field(depositing_fluxes[:,:,0])
        firstpass_coastal_flux_total = volume_field(coastal_overflow_fluxes[:,:,0])
        trapped_flux_total = volume_field(trapped_overflow_fluxes[:,:,0])
        flux_balance = incoming_flux_total - depositing_flux_total - 
            firstpass_coastal_flux_total - trapped_flux_total
        #=println(" chk first pass in ", incoming_flux_total,
            " dep ", depositing_flux_total, 
            " coast ", firstpass_coastal_flux_total,
            " trap ", trapped_flux_total,
            " bal ", flux_balance)=#
        #println(" flux bal? ", [ volume_field(incoming_fluxes[:,:,1]), 
        #    volume_field(incoming_fluxes[:,:,2]), volume_field(incoming_fluxes[:,:,3]),
        #    volume_field(incoming_fluxes[:,:,1]) + volume_field(incoming_fluxes[:,:,2]) ])
    end

    if n_trapped > 0 || trapped_flux_total > 0.# prepare for another pass 
        #println(" ocean overflows ",n_overflows," trapped ", n_trapped)
        #=
        println(" coast bal? ", [ volume_field(coastal_overflow_fluxes[:,:,1]), 
            volume_field(coastal_overflow_fluxes[:,:,2]), volume_field(coastal_overflow_fluxes[:,:,3]),
            volume_field(coastal_overflow_fluxes[:,:,1]) + volume_field(coastal_overflow_fluxes[:,:,2]) ])
        println(" trap bal? ", [ volume_field(trapped_overflow_fluxes[:,:,1]), 
            volume_field(trapped_overflow_fluxes[:,:,2]), volume_field(trapped_overflow_fluxes[:,:,3]),
            volume_field(trapped_overflow_fluxes[:,:,1]) + volume_field(trapped_overflow_fluxes[:,:,2]) ])
        =#

        redistribute_trapped_land_fluxes_to_coast!( trapped_overflow_fluxes,
            coastal_overflow_fluxes, new_submarine_mask )

        #=
        println(" coast bal? ", [ volume_field(coastal_overflow_fluxes[:,:,1]), 
            volume_field(coastal_overflow_fluxes[:,:,2]), volume_field(coastal_overflow_fluxes[:,:,3]),
            volume_field(coastal_overflow_fluxes[:,:,1]) + volume_field(coastal_overflow_fluxes[:,:,2]) ])
        println(" trap bal? ", [ volume_field(trapped_overflow_fluxes[:,:,1]), 
            volume_field(trapped_overflow_fluxes[:,:,2]), volume_field(trapped_overflow_fluxes[:,:,3]),
            volume_field(trapped_overflow_fluxes[:,:,1]) + volume_field(trapped_overflow_fluxes[:,:,2]) ])
        =#
        trapped_flux_total = volume_field(trapped_overflow_fluxes[:,:,0])
        if verbose == true
            incoming_flux_total = volume_field(incoming_fluxes[:,:,0])
            depositing_flux_total = volume_field(depositing_fluxes[:,:,0])
            new_coastal_flux_total = volume_field(coastal_overflow_fluxes[:,:,0])
            new_trapped_flux_total = volume_field(trapped_overflow_fluxes[:,:,0])
            mass_balance = new_coastal_flux_total - firstpass_coastal_flux_total - trapped_flux_total
            flux_balance = incoming_flux_total - depositing_flux_total - 
                new_coastal_flux_total - new_trapped_flux_total
            #=println(" chk redist in ", incoming_flux_total,
                " dep ", depositing_flux_total, 
                " coast ", new_coastal_flux_total,
                " trap ", new_trapped_flux_total,
                " mass bal ", mass_balance,
                " flux bal ", flux_balance)=#
            
        end
    end
    
    return coastal_overflow_fluxes, depositing_fluxes, new_submarine_mask, n_overflows + n_trapped
end
function redistribute_trapped_land_fluxes_to_coast!( trapped_land_fluxes,
    coastal_fluxes, submarine_mask )
    # move stuff from trapped overflow points to coastal deposition fluxes
    land_mask = 1 .- submarine_mask
    land_blobs = get_blobs(land_mask)
    trapped_volume = 
        [ volume_field( trapped_land_fluxes[:,:,0]  ), volume_field( trapped_land_fluxes[:,:,1]  ),volume_field( trapped_land_fluxes[:,:,2]  )]
    coastal_volume = 
        [ volume_field( coastal_fluxes[:,:,0]  ), volume_field( coastal_fluxes[:,:,1]  ), volume_field( coastal_fluxes[:,:,2]  )]
    for i_blob in 1:length(land_blobs)
        #println("about to blow ", i_blob)
        blob = land_blobs[i_blob]
        if volume_field( trapped_land_fluxes[:,:,0] .* blob ) > 0
            #println("moving stuff ",i_blob," ",trapped_volume," ",coastal_volume," ",trapped_volume + coastal_volume)
            blob_fringe_list = get_blob_fringe_exterior_list( blob )
            for ix in 1:nx
                for iy in 1:ny
                    if blob[ix,iy] == 1
                        for i_sedtype in 1:n_sediment_types
                            if trapped_land_fluxes[ix,iy,i_sedtype] != 0.
                                ixn, iyn = find_nearest_fringe_point(ix,iy,blob_fringe_list)
                                #println("fringe pt ",[ix,iy,ixn,iyn])
                                source_meter_flux = trapped_land_fluxes[ix,iy,i_sedtype]
                                source_volume_flux = source_meter_flux  *
                                    areabox[iy]
                                sink_meter_flux = source_volume_flux / areabox[iyn]
                                coastal_fluxes[ixn, iyn,i_sedtype] += 
                                    sink_meter_flux
                                coastal_fluxes[ixn, iyn,0] += 
                                    sink_meter_flux
                                trapped_land_fluxes[ix,iy,i_sedtype] -= source_meter_flux
                                trapped_land_fluxes[ix,iy,0] -= source_meter_flux
                            end
                        end
                    end
                end
            end
        trapped_volume = 
            [ volume_field( trapped_land_fluxes[:,:,0]  ), volume_field( trapped_land_fluxes[:,:,1]  ),volume_field( trapped_land_fluxes[:,:,2]  )]
        coastal_volume = 
            [ volume_field( coastal_fluxes[:,:,0]  ), volume_field( coastal_fluxes[:,:,1]  ), volume_field( coastal_fluxes[:,:,2]  )]
    
        #println("moved stuff ",i_blob," ",trapped_volume," ",coastal_volume," ",trapped_volume + coastal_volume)
            
        end
    end

    #trapped_land_fluxes .= 0.
end
function find_nearest_fringe_point(ix,iy,blob_fringe_list)
    if length(blob_fringe_list) == 0
        error("zero blob fringe list ",[ix,iy])
        return ix,iy
    end
    ixf = blob_fringe_list[1][1]; iyf = blob_fringe_list[1][2]
    min_distance = crow_flies(ycoords[iy],xcoords[ix],ycoords[iyf],xcoords[ixf])
    min_x = ixf; min_y = iyf
    for i_pos in 2:length(blob_fringe_list)
        ixf = blob_fringe_list[i_pos][1]; iyf = blob_fringe_list[i_pos][2]
        #println("trying ",[xcoords[ix],ycoords[iy],xcoords[ixf],ycoords[iyf]])
        test_distance = crow_flies(ycoords[iy],xcoords[ix],ycoords[iyf],xcoords[ixf])
        if test_distance < min_distance
            #println("dist ",[min_distance,test_distance])
            test_distance = min_distance
            min_x = ixf; min_y = iyf
        end
    end
    return min_x,min_y
end
function redistribute_trapped_land_fluxes_original_recipe!( trapped_overflow_fluxes,
    coastal_overflow_fluxes, new_submarine_mask )
    # move stuff from trapped overflow points to coastal deposition fluxes
    land_mask = 1 .- new_submarine_mask
    land_blobs = get_blobs(land_mask)

    for i_blob in 1:length(land_blobs)
        blob = land_blobs[i_blob]
        total_volume = 
            volume_field( trapped_overflow_fluxes[:,:,0] .* blob ) 
        if total_volume > 0
            #println("moving stuff ",i_blob," ",total_volume)
            total_volumes = fill(0.,0:n_sediment_types)
            for i_sedtype in 0:n_sediment_types
                total_volumes[i_sedtype] = 
                    volume_field( trapped_overflow_fluxes[:,:,i_sedtype] .* blob )
            end
            blob_fringe = get_blob_neighbors( blob )
            area_fringe = volume_field( blob_fringe )
            for ix in 1:nx
                for iy in 1:ny
                    if blob_fringe[ix,iy] == 1
                        #println("depositing to ",ix," ",iy)
                        for i_sedtype in 0:n_sediment_types
                            coastal_overflow_fluxes[ix,iy,i_sedtype] += 
                                total_volumes[i_sedtype] / area_fringe
                        end
                    end
                end
            end
        end
    end
    trapped_overflow_fluxes .= 0.
end

