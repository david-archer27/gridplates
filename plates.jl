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

# Rotation
function sphere2cart(coords)   # lat, longt
    elevation = ( coords[1] + 90. ) / 180. * pi # 0 at pole, pi/2 eq, pi pole
    azimuth = coords[2] / 180. * pi
    x = sin(elevation) * cos(azimuth)
    y = sin(elevation) * sin(azimuth)
    z = - cos(elevation)  #
    v = [x,y,z]
    return v
end
function cart2sphere(v)   # [x,y,z]
    x = v[1]; y = v[2]; z = v[3]
    elevation = acos( -z )
    azimuth = atan(y/x)
    if x < 0
        if y < 0
            azimuth = azimuth-pi
        else
            azimuth = azimuth+pi
        end
    end
    latitude = elevation / pi * 180. - 90.
    longitude = azimuth / pi * 180.
    return [latitude, longitude]
end
function rotation2matrix(rotation)
    lataxis = rotation.latitudeaxis
    longtaxis = rotation.longitudeaxis
    angle = rotation.angle
    m = axisangle2matrix(lataxis,longtaxis,angle)
    return m
end
function axisangle2matrix(lataxis,longtaxis,angle)
    vrot = sphere2cart([lataxis,longtaxis])
    ux = vrot[1]
    uy = vrot[2]
    uz = vrot[3]
    theta = - angle / 180. * pi

    a11 = cos(theta)+ux^2*(1-cos(theta))
    a12 = ux*uy*(1-cos(theta))-uz*sin(theta)
    a13 = ux*uz*(1-cos(theta))+uy*sin(theta)

    a21 = uy*ux*(1-cos(theta))+uz*sin(theta)
    a22 = cos(theta)+uy^2*(1-cos(theta))
    a23 = uy*uz*(1-cos(theta))-ux*sin(theta)

    a31 = uz*ux*(1-cos(theta))-uy*sin(theta)
    a32 = uz*uy*(1-cos(theta))+ux*sin(theta)
    a33 = cos(theta)+uz^2*(1-cos(theta))

    m = [ a11 a12 a13; a21 a22 a23; a31 a32 a33 ]
    rotm = RotMatrix{3}(m)
    return m #   rotm
end
function matrix2axisangle(m)
    cart = []
    push!(cart, m[2,3] - m[3,2])
    push!(cart, m[3,1] - m[1,3])
    push!(cart, m[1,2] - m[2,1])
    trace = m[1,1] + m[2,2] + m[3,3]
    theta = acos((trace - 1)/2)
    magnitude = 2 * sin(theta)
    cart ./= magnitude
    spherevec = cart2sphere(cart)
    latitudeaxis = spherevec[1]
    longitudeaxis = spherevec[2]
    angle = theta / pi * 180.
    return latitudeaxis,longitudeaxis,angle
end
function matrix2rotation(m,id,age,idparent)
    latitudeaxis,longitudeaxis,angle = matrix2axisangle(m)
    rotation = rotation_struct(id,age,latitudeaxis,longitudeaxis,angle,idparent)
    return rotation
end
function retrieve_rotations_for_plate(plateID)
    rotationsforplate = []
    for rotation in rotations
        if rotation.id == plateID
            push!(rotationsforplate,rotation)
        end
    end
    if length(rotationsforplate) == 0 && plateID != -1
        #println("1327 no rotations found for ", plateID)
        rotation = rotation_struct(plateID,0,0,0,0,0)
        push!(rotationsforplate,rotation)
        rotation = rotation_struct(plateID,1000,0,0,0,0)
        push!(rotationsforplate,rotation)
    end
    return rotationsforplate
end
function slerp(qa, qb, t)
    qm = fill(0.,4)
    # Calculate angle between them.
    cosHalfTheta = qa.q.s * qb.q.s + qa.q.v1 * qb.q.v1 + qa.q.v2 * qb.q.v2 + qa.q.v3 * qb.q.v3

    if cosHalfTheta < 0 # it took the long path, so negate one end
        #println("long path")
        qa = one(QuatRotation) / qa
        cosHalfTheta = - qa.q.s * qb.q.s - qa.q.v1 * qb.q.v1 - qa.q.v2 * qb.q.v2 - qa.q.v3 * qb.q.v3
    end
    # if qa=qb or qa=-qb then theta = 0 and we can return qa
    if (abs(cosHalfTheta) >= 1.0)
        return qa
    end
    # Calculate temporary values.
    halfTheta = acos(cosHalfTheta)
    sinHalfTheta = sqrt(1.0 - cosHalfTheta*cosHalfTheta)
    # if theta = 180 degrees then result is not fully defined
    # we could rotate around any axis normal to qa or qb
    if (abs(sinHalfTheta) < 0.001)
        qm[1] = (qa.q.s * 0.5 + qb.q.s * 0.5)
        qm[2] = (qa.q.v1 * 0.5 + qb.q.v1 * 0.5)
        qm[3] = (qa.q.v2 * 0.5 + qb.q.v2 * 0.5)
        qm[4] = (qa.q.v3 * 0.5 + qb.q.v3 * 0.5)
        qret = QuatRotation(qm)
        return qret
    end
    ratioA = sin((1 - t) * halfTheta) / sinHalfTheta
    ratioB = sin(t * halfTheta) / sinHalfTheta
    # calculate Quaternion.
    qm[1] = (qa.q.s * ratioA + qb.q.s * ratioB)
    qm[2] = (qa.q.v1 * ratioA + qb.q.v1 * ratioB)
    qm[3] = (qa.q.v2 * ratioA + qb.q.v2 * ratioB)
    qm[4] = (qa.q.v3 * ratioA + qb.q.v3 * ratioB)
    qret = QuatRotation(qm)
    return qret
end
function interpolate_mplate(plateID,age)
    # computes the rotation matrix for the plate at the time,
    # and returns the parent ID
    #println(1569," ", plateID)
    timestepoffset = 0.001
    filteredrotations = retrieve_rotations_for_plate(plateID)
    #println(plateID)
    # make use of any direct hits
    mplates = Any[]; parents = Any[]
    for rotation in filteredrotations
        if rotation.age == age
            #            println(" got a hit ", plateID)
            mplate = rotation2matrix(rotation)
            parent = rotation.parent
            push!(mplates,copy(mplate))
            push!(parents,parent)
        end
    end
    if length(parents) > 1
        println(" got multiple hits ", plateID)
    end
    if length(parents) == 0
        # get here if there are no direct age hits
        if length(filteredrotations) == 0 && plateID != -1
            println("  No rotations found for ", plateID)
        end
        if age - timestepoffset < filteredrotations[1].age
            println("  Time too late for plate ", plateID)
            return [norotation],[0]
        end
        if age - timestepoffset > filteredrotations[end].age
            println("  Time too early for plate ", plateID)
            return [norotation],[0]
        end
        locallistofages = []
        for rotation in filteredrotations
            push!(locallistofages,rotation.age)
        end
        ibelow = search_sorted_below(locallistofages,age)
        fractions = fill(0.,length(locallistofages))
        if ibelow < length(locallistofages)
            if ibelow > 0       #  tbelow    time               t(below+1)
                fractions[ibelow+1] = (age - locallistofages[ibelow]) /
                    (locallistofages[ibelow+1] - locallistofages[ibelow])
                fractions[ibelow] = 1. - fractions[ibelow+1]
                if filteredrotations[ibelow].parent != filteredrotations[ibelow+1].parent
                    println("  Oops, switching parents ", plateID)
                else
                    mplate1 = rotation2matrix(filteredrotations[ibelow])
                    mplate2 = rotation2matrix(filteredrotations[ibelow+1])
                    qa = QuatRotation(mplate1)
                    qb = QuatRotation(mplate2)
                    t = fractions[ibelow+1]
                    #println(qa," ",qb)
                    qm = slerp(qa,qb,t)
                    mplate = RotMatrix(qm)

                    parent = filteredrotations[ibelow].parent
                    push!(mplates,copy(mplate))
                    push!(parents,parent)
                end
            else  # ibelow < 0
                println("how did I get here? ", ibelow, " ", locallistofages," ", plateID)
            end
        else # ibelow >= length
            if ibelow == length(locallistofages) &&
                locallistofages[ibelow] + timestepoffset >= age
                mplate = rotation2matrix(filteredrotations[ibelow])
                push!(mplates,copy(mplate))
                parent = filteredrotations[ibelow].parent
                push!(parents,parent)
                #                println(" using first rotation line", plateID)
            end
        end
    end
    return mplates, parents
end
function resolve_rotation_matrix!(plateID)
    return resolve_rotation_matrix!(plateID,world.age)
end
function resolve_rotation_matrix!(plateID,age)
    # recursively unwinds the heirarchy of parent plates
    # returns the rotation matrix, and records it in the plate structure
    # archives the parentstack into lastparentstack
    # detects rotation jumps at ancestry changes
    #    println("Resolving ", plateID," ",age)
    timestepoffset = 0.001
    hackedage = age + timestepoffset
    if plateID == -1
        return norotation
    end
    if haskey(plates,plateID) == false
        plates[plateID] = create_blank_plate(plateID)
    end
    if plates[plateID].resolvetime == age # nice round numbers for bookkeeping
        #        println(" solved ",plateID)
        return plates[plateID].rotationmatrix
    end
    mplates, idparentsbase = interpolate_mplate(plateID,hackedage)
    threadsfound = length(mplates)  # number of simultaneous rotations
    mplateoutlist = Any[]           # list of matrices, one per rotation
    if threadsfound == 1
        idparent = idparentsbase[1]
        mplate = mplates[1]
        #        println(295," starting ",idparent," ",plates[plateID].parentstack)
        plates[plateID].lastparentstack = deepcopy(plates[plateID].parentstack)
        plates[plateID].parentstack = []
        push!(plates[plateID].parentstack,idparent)
        if idparent != 0
            # gets here if it has a legit parent to deconvolve
            mparent = resolve_rotation_matrix!(idparent,age)
            mplate = mplate * mparent
            for ipar = plates[idparent].parentstack
                push!(plates[plateID].parentstack,ipar)
            end
        end
        if plates[plateID].parentstack != plates[plateID].lastparentstack &&
            plates[plateID].lastparentstack != [] &&
            plateID in world.plateID
            mplateold = plates[plateID].rotationmatrix
            mplatenew = mplate
    end
        push!(mplateoutlist,mplate)
    end
    if threadsfound > 2
        println(" more than 2 simultaneous rotations", plateID)
    end
    if threadsfound == 2
        println(" dup results ", plateID," ",idparentsbase)
        parentstackbythread = []
        for iresult in 1:threadsfound
        #        print(" -> ", idparent)
            mplate = mplates[iresult]
            idparent = idparentsbase[iresult]
            if idparent != 0
            # gets here if it has a legit parent to deconvolve
                mparent = resolve_rotation_matrix!(idparent,hackedage)
                # puts parent stack into globalparentstack
                mplate = mparent * mplate
            end
            push!(mplateoutlist, copy(mplate))
        end
    end
    if threadsfound > 2
        println(" yikes! more than two threads ", plateID)
    end
    if haskey(plates,plateID) == false
        plates[plateID] = create_blank_plate(plateID)
    end
    if length(mplateoutlist) == 0
        println(" not getting any threads ", plateID)
        push!(mplateoutlist,norotation)
    end
    plates[plateID].rotationmatrix = mplateoutlist[1]
    plates[plateID].resolvetime = age
    return mplateoutlist[1]
end
function resolve_rotation_matrices()
    #    println("  Updating rotations to ", age)
    for plateID in world.plateIDlist
        if plates[plateID].resolvetime != world.age
            if plates[plateID].resolvetime != -1
                mplateold = plates[plateID].rotationmatrix
                mplatenew = resolve_rotation_matrix!(plateID)
                #distance = meanplatetraveldistance(mplateold,mplatenew)
            else
                resolve_rotation_matrix!(plateID)
            end
        end
    end
    return
end

# Interpolation
function fill_world_from_plates()
    # project plate crust* and sediment* variables to world grid.
    #Threads.@threads 
    for iworld in 1:nx
        for jworld in 1:ny
            plateID = world.plateID[iworld,jworld]
            if plateID != -1 && haskey(plates,plateID)
                plate = plates[plateID]
                if plate.resolvetime != world.age && plateID != -1
                    println("updating rotation matrix in fillworld ", plateID)
                    resolve_rotation_matrix!(plateID)
                end
                iplate,jplate = nearest_plateij(plate.rotationmatrix,iworld,jworld)
                if plate.crust_type[iplate,jplate] == notatsurface # stupid edge case
                    # the closest plate point to the projection of the world point
                    # thinks it is outside of the exposed plate.
                    # dont fill it in because the projections are not reciprocal,
                    # dont map back.  the plate point would just toggle on and off.
                    world.crust_age[iworld,jworld] = 0.
                    world.crust_type[iworld,jworld] = ocean_crust
                    world.crust_thickness[iworld,jworld] = ocean_crust_h0
                    world.crust_density[iworld,jworld] = rho_ocean_crust
                    world.surface_type[iworld,jworld] = pelagic_seafloor
                else
                    world.crust_age[iworld,jworld] = plate.crust_age[iplate,jplate]
                    world.crust_type[iworld,jworld] = plate.crust_type[iplate,jplate]
                    world.crust_thickness[iworld,jworld] =
                        plate.crust_thickness[iplate,jplate]
                    world.crust_density[iworld,jworld] =
                        plate.crust_density[iplate,jplate]
                    world.surface_type[iworld,jworld] =
                        plate.surface_type[iplate,jplate]
                    world.sediment_thickness[iworld,jworld] = 
                        plate.sediment_thickness[iplate,jplate]                        
                    world.sediment_surface_fractions[iworld,jworld,:] .=
                        plate.sediment_surface_fractions[iplate,jplate,:]
                    world.sediment_layer_thickness[iworld,jworld,:] .= 
                        plate.sediment_layer_thickness[iplate,jplate,:]                        
                    world.sediment_layer_fractions[iworld,jworld,:,:] .=
                        plate.sediment_layer_fractions[iplate,jplate,:,:]

                end
            end
        end
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
    for iplate in 1:nx
        for jplate in 1:ny
            iworld, jworld = nearest_worldij(plate.rotationmatrix,iplate,jplate)
            if plate.crust_type[iplate,jplate] == notyetdefined # initialize
                if world.plateID[iworld,jworld] == plate.plateID
                    plate.crust_type[iplate,jplate] = world.crust_type[iworld,jworld]
                    plate.crust_age[iplate,jplate] = 0.
                    plate.crust_thickness[iplate,jplate] =
                        world.crust_thickness[iworld,jworld]
                    plate.surface_type[iplate,jplate] =
                        world.surface_type[iworld,jworld]
                    if plate.crust_type[iplate,jplate] == ocean_crust
                        plate.crust_density[iplate,jplate] = rho_ocean_crust
                    else
                        plate.crust_density[iplate,jplate] = rho_continent_crust
                    end
                else
                    delete_plate_point!(plate,iplate,jplate)
                end
            elseif world.plateID[iworld,jworld] == plate.plateID
                # our plate outcrops here in the world grid
                if world.crust_type[iworld,jworld] == ocean_crust
                    if plate.crust_type[iplate,jplate] == notatsurface
                        # ridge spreading hopefully
                        record_flux_into_world_diags("ocean_created_plate_area",
                            jplate,iworld,jworld)
                        plate.crust_age[iplate,jplate] = 0. # new ocean crust
                        plate.crust_type[iplate,jplate] = ocean_crust
                        plate.crust_thickness[iplate,jplate] = ocean_crust_h0
                        plate.crust_density[iplate,jplate] = rho_ocean_crust
                    elseif plate.crust_type[iplate,jplate] == continent_crust
                        record_flux_into_world_diags("continent_2_ocean_plate_area",
                            jplate,iworld,jworld)
                        plate.crust_age[iplate,jplate] = 0. # new ocean crust
                        plate.crust_type[iplate,jplate] = ocean_crust
                        plate.crust_thickness[iplate,jplate] = ocean_crust_h0
                        plate.crust_density[iplate,jplate] = rho_ocean_crust
                    end
                else # world.crust_type[iworld,jworld] == continent_crust
                    if plate.crust_type[iplate,jplate] == notatsurface
                        # continent has formed where nothing was
                        record_flux_into_world_diags("continent_created_plate_area",
                            jplate,iworld,jworld)
                        plate.crust_age[iplate,jplate] = 0.
                        plate.crust_type[iplate,jplate] = continent_crust
                        plate.crust_thickness[iplate,jplate] =
                            continent_crust_h0
                        plate.crust_density[iplate,jplate] = rho_continent_crust
                    elseif plate.crust_type[iplate,jplate] == ocean_crust
                        # ocean type has transformed into continental type
                        record_flux_into_world_diags("ocean_2_continent_plate_area",
                            jplate,iworld,jworld)
                        plate.crust_age[iplate,jplate] = 0.
                        plate.crust_type[iplate,jplate] = continent_crust
                        plate.crust_thickness[iplate,jplate] =
                            continent_crust_h0
                        plate.crust_density[iplate,jplate] = rho_continent_crust
                    end
                end # updated the status of the plate grid point if required
            else       # world map says we dont exist here
                if plate.crust_type[iplate,jplate] == ocean_crust
                    # but we did last time step -- Subduction!
                    record_flux_into_world_diags("ocean_subduct_plate_area",
                        jplate,iworld,jworld)
                    record_flux_into_world_diags("ocean_subduct_age_plate_area",
                        jplate,iworld,jworld,plate.crust_age[iplate,jplate])
                    for ibin in 1:n_sediment_time_bins
                        record_flux_into_world_diags("ocean_subduct_sediment_plate_volume",
                            jplate,iworld,jworld,
                            plate.sediment_layer_thickness[iplate,jplate,ibin])
                        for i_sedtype in 1:n_sediment_types
                            flux_amt = plate.sediment_layer_thickness[iplate,jplate,ibin] *
                                plate.sediment_layer_fractions[iplate,jplate,i_sedtype,ibin]
                            record_sediment_fraction_flux_into_world_diags(
                                "ocean_subduct_sediment_fraction_volume",i_sedtype,
                                jplate,iworld,jworld,flux_amt)
                        end
                    end
                    #record_ocean_disappearing_from_plates(plate,iplate,jplate,iworld,jworld)
                elseif plate.crust_type[iplate,jplate] >= continent_crust
                    record_flux_into_world_diags("continent_subduct_plate_area",
                        jplate,iworld,jworld)
                    #record_continent_disappearing_from_plates(plate,iplate,jplate,iworld,jworld)
                end
                delete_plate_point!(plate,iplate,jplate)
            end # initialize or update
        end # loop over plate grid
    end
end
function remask_plates()   # subduction and crust creation
    Threads.@threads for plateID in world.plateIDlist     #[25:25]# 1:length(plates)
        if haskey(plates,plateID) == false
            plates[plateID] = create_blank_plate(plateID)
        end
        remask_plate!( plates[plateID] )
    end
    for (plateID,plate) in plates
        if !(plateID in world.plateIDlist)
            delete!(plates,plateID)
        end
    end
    return
end
function nearest_plateij(mplate,iworld,jworld)  # ray trace from a world i,j to a plate
    longitudeworld = xcoords[iworld]
    latitudeworld = ycoords[jworld]
    if latitudeworld > 90.
        latitudeworld = 180. - latitudeworld
        longitudeworld += 180.
    end
    if latitudeworld < -90.
        latitudeworld = - 180. + latitudeworld
        longitudeworld += 180.
    end
    while longitudeworld < -180.; longitudeworld += 360.; end
    while longitudeworld >  180.; longitudeworld -= 360.; end
    vworld = sphere2cart([latitudeworld,longitudeworld])
    vplate = mplate * vworld  # rotate into plate frame
    mapcoords = cart2sphere(vplate)
    iplate = search_sorted_nearest(xcoords,mapcoords[2])
    jplate = search_sorted_nearest(ycoords,mapcoords[1])
    return iplate, jplate
end
function nearest_worldij(rotationmatrix,iplate,jplate)
    latitudepoint = ycoords[jplate]
    longitudepoint = xcoords[iplate]
    vpoint = sphere2cart([latitudepoint,longitudepoint])
    newvpoint = transpose(rotationmatrix) * vpoint
    # unrotating back to world grid
    mapcoords = cart2sphere(newvpoint)
    iworld = search_sorted_nearest(xcoords,mapcoords[2])
    jworld = search_sorted_nearest(ycoords,mapcoords[1])
    return iworld, jworld
end
function nearest_worldij(plate,rotation,ifield,jfield)
    m = rotation2matrix(rotation)
    iworld,jworld = nearest_worldij(m,ifield,jfield)
    return iworld,jworld
end
function nearest_other_plate_ij_unused(iplate1,jplate1,mplate1,mplate2)
    plate1v = sphere2cart( [ ycoords[jplate1], xcoords[iplate1] ] )
    worldv = transpose(mplate1) * plate1v
    plate2v = mplate2 * worldv
    plate2coords = cart2sphere(plate2v) # [lat, longt]
    iplate2 = search_sorted_nearest(xcoords,plate2coords[2])
    jplate2 = search_sorted_nearest(ycoords,plate2coords[1])
    return iplate2,jplate2
end
function projected_plate_maskfield!(worldmaskfield,rotationmatrix)
    # find all points in plate that project to world region
    newfield = fill(0,nx,ny)
    for ixplate in 1:nx
        for iyplate in 1:ny
            ixworld,iyworld = nearest_worldij(rotationmatrix,
                ixplate,iyplate)
            if ixworld > 0
                if worldmaskfield[ixworld,iyworld] == 1
                    newfield[ixplate,iyplate] = 1
                end
            end
        end
    end
    return newfield
end
function projected_world_maskfield!(platemaskfield,rotationmatrix)
    # find all points in plate that project to world region
    worldfield = fill(0,nx,ny)
    for ixworld in 1:nx
        for iyworld in 1:ny

            ixplate,iyplate = nearest_plateij(rotationmatrix,
                ixworld,iyworld)

            if platemaskfield[ixplate,iyplate] == 1
                worldfield[ixworld,iyworld] = 1
            end
        end
    end
    return worldfield
end

# Routines dealing with changing plate IDs
function update_changing_plateIDs()
    newtime = world.age
    newIDmap = read_plateIDs(world.age)
    changingplates = changing_platelists(newIDmap)
    #println("changing ", changingplates)
    #Threads.@threads 
    for platepair in changingplates
        oldplateID = platepair[1]; newplateID = platepair[2]
        update_two_plates!(oldplateID,newplateID)
        #        end
    end
    world.plateID = newIDmap
    world.plateIDlist = find_plateID_list()
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
            end
        end
    end
    return
end
function copy_plate_point_plate12coord!(newplate,oldplate,ixold,iyold,ixnew,iynew)
    if oldplate.crust_type[ixold,iyold] == notatsurface #||
        #        oldplate.crust_type[ixold,iyold] == ocean_crust  # moving ridge
        # edge cases, messy spreading
        newplate.crust_type[ixnew,iynew] = ocean_crust
        newplate.crust_thickness[ixnew,iynew] = ocean_crust_h0
        newplate.crust_density[ixnew,iynew] = rho_ocean_crust
        newplate.crust_age[ixnew,iynew] = 0.
        newplate.surface_type[ixnew,iynew] = pelagic_seafloor
    else
        newplate.crust_type[ixnew,iynew] = oldplate.crust_type[ixold,iyold]
        newplate.crust_age[ixnew,iynew] = oldplate.crust_age[ixold,iyold]
        newplate.crust_thickness[ixnew,iynew] = oldplate.crust_thickness[ixold,iyold]
        newplate.crust_density[ixnew,iynew] = oldplate.crust_density[ixold,iyold]
        newplate.surface_type[ixnew,iynew] = oldplate.surface_type[ixold,iyold]
        newplate.sediment_thickness[ixnew,iynew] = 
            oldplate.sediment_thickness[ixold,iyold]
        newplate.sediment_surface_fractions[ixnew,iynew,:] .= 
            oldplate.sediment_surface_fractions[ixold,iyold,:]
        newplate.sediment_layer_thickness[ixnew,iynew,:] .= 
            oldplate.sediment_layer_thickness[ixold,iyold,:]
        newplate.sediment_layer_fractions[ixnew,iynew,:,:] .= 
            oldplate.sediment_layer_fractions[ixold,iyold,:,:]
    end
end
function delete_plate_point!(oldplate,ixold,iyold)
    oldplate.crust_type[ixold,iyold] = notatsurface
    oldplate.crust_age[ixold,iyold] = 0.
    oldplate.crust_thickness[ixold,iyold] = 0.
    oldplate.crust_density[ixold,iyold] = 0.
    oldplate.sediment_thickness[ixold,iyold] = 0.
    oldplate.sediment_surface_fractions[ixold,iyold,:] .= 0.
    oldplate.sediment_layer_thickness[ixold,iyold,:] .= 0.
    oldplate.sediment_layer_fractions[ixold,iyold,:,:] .= 0.
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

# Utilities
function apply_geomorphology_changes_to_plates()
    Threads.@threads for plateID in world.plateIDlist
        plate = plates[plateID]
        for iplate in 1:nx
            for jplate in 1:ny
                iworld, jworld = nearest_worldij(plate.rotationmatrix,iplate,jplate)
                if world.plateID[iworld,jworld] == plateID
                    plates[plateID].crust_thickness[iplate,jplate] =
                        world.crust_thickness[iworld,jworld]
                    plates[plateID].crust_type[iplate,jplate] =
                        world.crust_type[iworld,jworld]
                    plates[plateID].surface_type[iplate,jplate] =
                        world.surface_type[iworld,jworld]
                    plates[plateID].sediment_thickness[iplate,jplate] =
                        world.sediment_thickness[iworld,jworld]
                    plates[plateID].sediment_surface_fractions[iplate,jplate,:] .=
                        world.sediment_surface_fractions[iworld,jworld,:]
                    plates[plateID].sediment_layer_thickness[iplate,jplate,:] .=
                        world.sediment_layer_thickness[iworld,jworld,:]
                    plates[plateID].sediment_layer_fractions[iplate,jplate,:,:] .=
                        world.sediment_layer_fractions[iworld,jworld,:,:]
                end
            end
        end
    end
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




