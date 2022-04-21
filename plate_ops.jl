# routines dealing with changing plate IDs
function copy_plate_point_plate12coord!(newplate,oldplate,ixold,iyold,ixnew,iynew)
    if oldplate.crust_type[ixold,iyold] == notatsurface #||
        #        oldplate.crust_type[ixold,iyold] == ocean_type  # moving ridge
        # edge cases, messy spreading
        newplate.crust_type[ixnew,iynew] = ocean_type
        newplate.crust_thickness[ixnew,iynew] = ocean_crust_h0
        newplate.crust_density[ixnew,iynew] = rho_ocean_crust
        newplate.crust_age[ixnew,iynew] = 0.
    else
        newplate.crust_type[ixnew,iynew] = oldplate.crust_type[ixold,iyold]
        newplate.crust_age[ixnew,iynew] = oldplate.crust_age[ixold,iyold]
        newplate.crust_thickness[ixnew,iynew] = oldplate.crust_thickness[ixold,iyold]
        newplate.crust_density[ixnew,iynew] = oldplate.crust_density[ixold,iyold]
        for ibin in 1:n_sediment_layer_time_bins
            newplate.sediment_thickness[ixnew,iynew,ibin] =
                oldplate.sediment_thickness[ixold,iyold,ibin]
            for itype in 1:n_sediment_types
                newplate.sediment_fractions[ixnew,iynew,ibin,itype] =
                    oldplate.sediment_fractions[ixold,iyold,ibin,itype]
            end
        end
    end
    return
end
function delete_plate_point!(oldplate,ixold,iyold)
    oldplate.crust_type[ixold,iyold] = notatsurface
    oldplate.crust_age[ixold,iyold] = 0.
    for ibin in 1:n_sediment_layer_time_bins
        oldplate.sediment_thickness[ixold,iyold,ibin] = 0.
        for itype in 1:n_sediment_types
            oldplate.sediment_fractions[ixold,iyold,ibin,itype] = 0.
        end
    end
    return
end
function update_two_plates!(oldplateID,newplateID)
        # called when one or more grid points change their plateIDs
    newtime = world.age; oldtime = newtime + time_step
    #println("updating ", oldplateID," ",newplateID," ",oldtime," ",newtime)
    if haskey(plates,oldplateID) == false
        println("creating plate to migrate data from ", oldplateID)
        plates[oldplateID] = create_blank_plate(newplateID)
    end
    oldplate = plates[oldplateID]
    if haskey(plates,newplateID) == false
        plates[newplateID] = create_blank_plate(newplateID)
        resolve_rotation_matrix!(newplateID,newtime)
    end
    newplate = plates[newplateID]
    if newplate.resolvetime < 0
        println("  Trying to update to an unrotated plate ", oldplate.plateID, " ", newplate.plateID)
    end
    oldIDmap = world.plateID
    newIDmap = read_plateIDs(newtime)
    oldIDmask = generate_mask_field(oldIDmap,oldplateID)
    # grid points in world grid of old plateID number
    newIDmask = generate_mask_field(newIDmap,newplateID)
    # points in world grid where the new plate ID is found
    needschangingworldmask = oldIDmask .* newIDmask
    # points in world grid where old ID changes to new ID
    firstoldtime,lastoldtime = first_last_appearance(oldplateID)
    firstnewtime,lastnewtime = first_last_appearance(newplateID)
    switchovertime = newtime # use rotations for newtime if possible
    if lastoldtime > newtime
        switchovertime = lastoldtime  # or back up to  end time of old plate
    end
    newplaterotation = resolve_rotation_matrix!(newplateID,switchovertime)
    needschangingnewplatemask = projected_plate_maskfield!(needschangingworldmask,
        newplaterotation)
    oldplaterotation = resolve_rotation_matrix!(oldplateID,switchovertime)
    needsdeletingoldplatemask = projected_plate_maskfield!(needschangingworldmask,
        oldplaterotation)
    for ixnewplate in 1:nx
        for iynewplate in 1:ny
            if needschangingnewplatemask[ixnewplate,iynewplate] > 0
                ixworld,iyworld = nearest_worldij(newplate.rotationmatrix,
                    ixnewplate,iynewplate)
                if ixworld > 0
                    ixoldplate,iyoldplate = nearest_plateij(oldplate.rotationmatrix,
                        ixworld,iyworld)
                    copy_plate_point_plate12coord!(newplate,oldplate,
                        ixoldplate,iyoldplate,ixnewplate,iynewplate)
                    record_flux_into_world_diags("IDchange_plate_area",iynewplate,
                        ixworld,iyworld)
                    #record_plateID_change_fillin(newplate,
                    #    ixworld,iyworld,
                    #    ixnewplate,iynewplate)
                end
            end
        end
    end

    for ixplate in 1:nx
        for iyplate in 1:ny
            if needsdeletingoldplatemask[ixplate,iyplate] > 0
                delete_plate_point!(oldplate,ixplate,iyplate)
                ixworld,iyworld = nearest_worldij(oldplate.rotationmatrix,
                    ixplate,iyplate)
                record_flux_into_world_diags("IDchange_plate_area",iyplate,
                    ixworld,iyworld,-1.)
                #record_plateID_change_delete(oldplate,
                #    ixworld,iyworld,
                #    ixplate,iyplate)
            end
        end
    end
    return
end
function changing_platelists(newIDmap)
    oldIDmap = world.plateID
    oldtime = world.age + time_step
    fp = 2
    changepairs = []
    for ixworld in 1+fp:nx-fp
        for iyworld in 1+fp:ny-fp
            if oldIDmap[ixworld,iyworld] !=
                newIDmap[ixworld,iyworld]
                foundone = 0
                for ifp in -fp:fp
                    if oldIDmap[ixworld+ifp,iyworld+ifp] ==
                        newIDmap[ixworld+ifp,iyworld+ifp]
                        foundone = 1
                    end
                end
                if foundone == 0 # truly a plate change not a migrating ridge
                    oldplateID = oldIDmap[ixworld,iyworld]
                    newplateID = newIDmap[ixworld,iyworld]
                    key = [ oldplateID, newplateID ]
                    # encode both IDs into single key
                    #key = oldplateID * 10000000 + newplateID
                    push!(changepairs,key)
                end
            end
        end
    end
    # list of all gridpoints in successive world grids where the ID changes
    #println("key list", changepairs)
    changepairs = unique(changepairs)
    #println("key list", changepairs)
    return changepairs
end
function update_changing_plateIDs()
    newtime = world.age
    newIDmap = read_plateIDs(world.age)
    changingplates = changing_platelists(newIDmap)
    #println("changing ", changingplates)
    Threads.@threads for platepair in changingplates
        oldplateID = platepair[1]; newplateID = platepair[2]
        update_two_plates!(oldplateID,newplateID)
        #        end
    end
    world.plateID = newIDmap
    world.plateIDlist = find_plateID_list()
    return
end
function find_plateID_list()
    return find_plateID_list( world.plateID )
end
function find_plateID_list(plateIDmap)
    plateIDlist = sort(unique(plateIDmap))
    killme = 0
    for i in 1:length(plateIDlist)
        if plateIDlist[i] == -1
            killme = i
        end
    end
    if killme > 0
        deleteat!(plateIDlist,killme)
    end
    return plateIDlist
end
function increment_plate_age()   # updates the age field in plate grids
    #Threads.@threads
    for plateID in world.plateIDlist
        plate = plates[plateID]
        for iplate in 1:nx
            for jplate in 1:ny
                if plate.crust_type[iplate,jplate] != notatsurface
                    plate.crust_age[iplate,jplate] += time_step
                end
            end
        end
    end
    return
end
function first_last_appearance(plateID)
    filteredrotations = retrieve_rotations_for_plate(plateID)
    times = []
    for rotation in filteredrotations
        append!(times, rotation.age)
    end
    if length(times) > 0
        last = minimum(times)
        first = maximum(times)
        return first, last
    else
        return 0, 0
    end
end

function remask_plate!(plate)
    # fill in plate from world grid.
    # finds grid points of subduction and crust creation on plate fields.
    # imports crust_thickness as it stands in world grid.
    # when a plate point starts corresponding to a continent point in the world
    #  it gets the world crust_thickness etc, before orogeny due to the current
    #  step.
    age = world.age
    if plate.resolvetime != age
        resolve_rotation_matrix!(plate.plateID,age)
    end
    #Threads.@threads
    for iplate in 1:nx
        for jplate in 1:ny
            iworld, jworld = nearest_worldij(plate.rotationmatrix,iplate,jplate)
            if plate.crust_type[iplate,jplate] == notyetdefined # initialize
                if world.plateID[iworld,jworld] == plate.plateID
                    plate.crust_type[iplate,jplate] = world.crust_type[iworld,jworld]
                    plate.crust_age[iplate,jplate] = 0.
                    plate.crust_thickness[iplate,jplate] =
                        world.crust_thickness[iworld,jworld]
                    if plate.crust_type[iplate,jplate] == ocean_type
                        plate.crust_density[iplate,jplate] = rho_ocean_crust
                    else
                        plate.crust_density[iplate,jplate] = rho_continent_crust
                    end
                else #
                    plate.crust_type[iplate,jplate] = notatsurface
                end
            elseif world.plateID[iworld,jworld] == plate.plateID
                # our plate outcrops here in the world grid
                if world.crust_type[iworld,jworld] == ocean_type
                    if plate.crust_type[iplate,jplate] == notatsurface
                        # ridge spreading hopefully
                        record_flux_into_world_diags("ocean_created_plate_area",
                            jplate,iworld,jworld)
                        plate.crust_type[iplate,jplate] = ocean_type
                        plate.crust_age[iplate,jplate] = 0. # new ocean crust
                        plate.crust_thickness[iplate,jplate] = ocean_crust_h0
                        plate.crust_density[iplate,jplate] = rho_ocean_crust
                    elseif plate.crust_type[iplate,jplate] >= continent_type
                        record_flux_into_world_diags("continent_2_ocean_plate_area",
                            jplate,iworld,jworld)
                        #record_continent_2_ocean_from_plates(jplate,iworld,jworld)
                        #println("disappearing continent ",plate.plateID," ",
                        #    iplate," ",jplate)
                        plate.crust_type[iplate,jplate] = ocean_type
                        plate.crust_age[iplate,jplate] = 0. # new ocean crust
                        plate.crust_thickness[iplate,jplate] = ocean_crust_h0
                        plate.crust_density[iplate,jplate] = rho_ocean_crust
                    end # or else plate was already ocean, do nothing
                else # world.crust_type[iworld,jworld] == continent_type
                    if plate.crust_type[iplate,jplate] == notatsurface
                        # continent has formed where nothing was
                        record_flux_into_world_diags("continent_created_plate_area",
                            jplate,iworld,jworld)
                        #record_continent_creation_from_plates(jplate,iworld,jworld)
                    elseif plate.crust_type[iplate,jplate] == ocean_type
                        # ocean type has transformed into continental type
                        record_flux_into_world_diags("ocean_2_continent_plate_area",
                            jplate,iworld,jworld)
                        #record_ocean_2_continent_from_plates(jplate,iworld,jworld)
                    end
                    if plate.crust_type[iplate,jplate] < continent_type
                        plate.crust_type[iplate,jplate] = continent_type
                        plate.crust_age[iplate,jplate] = 0.
                        plate.crust_thickness[iplate,jplate] =
                            continent_crust_h0
                        plate.crust_density[iplate,jplate] = rho_continent_crust
                    end
                    if plate.crust_type[iplate,jplate] >= continent_type
                        plate.crust_type[iplate,jplate] =
                            min(continent_crust_h0,plate.crust_type[iplate,jplate])
                    end
                end # updated the status of the plate grid point if required
            else       # world map says we dont exist here
                if plate.crust_type[iplate,jplate] == ocean_type
                    # but we did last time step
                    record_flux_into_world_diags("ocean_subduct_plate_area",
                        jplate,iworld,jworld)
                    record_flux_into_world_diags("ocean_subduct_age_plate_area",
                        jplate,iworld,jworld,plate.crust_age[iplate,jplate])
                    for ibin in 1:n_sediment_layer_time_bins
                        record_flux_into_world_diags("ocean_subduct_sediment_plate_volume",
                            jplate,iworld,jworld,
                            plate.sediment_thickness[iplate,jplate,ibin])
                        for itype in 1:n_sediment_types
                            fluxname = "ocean_subduct_sediment_" *
                                sediment_type_names[itype] *
                                "_plate_volume"
                            flux_amt = plate.sediment_thickness[iplate,jplate,ibin] *
                                plate.sediment_fractions[iplate,jplate,ibin,itype]
                            record_flux_into_world_diags(fluxname,
                                jplate,iworld,jworld,flux_amt)
                        end
                    end
                    #record_ocean_disappearing_from_plates(plate,iplate,jplate,iworld,jworld)
                elseif plate.crust_type[iplate,jplate] >= continent_type
                    record_flux_into_world_diags("continent_subduct_plate_area",
                        jplate,iworld,jworld)
                    #record_continent_disappearing_from_plates(plate,iplate,jplate,iworld,jworld)
                end
                if plate.crust_type[iplate,jplate] != notatsurface
                    delete_plate_point!(plate,iplate,jplate)
                end
            end # initialize or update
        end # loop over plate grid
    end
    return
end
function remask_plates()   # subduction and crust creation
    Threads.@threads for plateID in world.plateIDlist     #[25:25]# 1:length(plates)
        if haskey(plates,plateID) == false
            plates[plateID] = create_blank_plate(plateID)
        end
        plate = plates[plateID]
        #        print(plate.plateID,"\n")
        #        if age <= plate.firsttime
        #            if age >= plate.lasttime
        remask_plate!(plate)
        #println(plateID," ", areafield_total(world.subducted_area))
        #                   end
        #        end
    end
    for (plateID,plate) in plates
        if !(plateID in world.plateIDlist)
            delete!(plates,plateID)
        end
    end
    return
end
function top_filled_sediment_box(plate,iplate,jplate)
    ibinmax = search_sorted_below(sediment_layer_time_bins,world.age)
    ibinmax = max(ibinmax,1)
    ibintop = 1
    for ibin in 2:ibinmax
        if plate.sediment_thickness[iplate,jplate,ibin] >= 0.
            ibintop = ibin
        end
    end
    return ibintop
end
