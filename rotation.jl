
# rotation routines
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
                    world.crust_type[iworld,jworld] = ocean_type
                    world.crust_thickness[iworld,jworld] = ocean_crust_h0
                    world.crust_density[iworld,jworld] = rho_ocean_crust
                else
                    world.crust_age[iworld,jworld] = plate.crust_age[iplate,jplate]
                    world.crust_type[iworld,jworld] = plate.crust_type[iplate,jplate]
                    world.crust_thickness[iworld,jworld] =
                        plate.crust_thickness[iplate,jplate]
                    world.crust_density[iworld,jworld] =
                        plate.crust_density[iplate,jplate]
                    sediment_thickness = 0.
                    for ibin in 1:n_sediment_layer_time_bins
                        sediment_thickness +=
                            plate.sediment_thickness[iplate,jplate,ibin]
                    end
                    world.sediment_thickness[iworld,jworld] = sediment_thickness

                    ibin = top_filled_sediment_box(plate,iplate,jplate)
                    for itype in 1:n_sediment_types
                        world.sediment_fractions[iworld,jworld,itype] =
                            plate.sediment_fractions[iplate,jplate,ibin,itype]
                    end
                end
            end
        end
    end
    return
end
function fill_world_orogeny(footprint) # requires at least create_world(age)
    age = world.age
    plateIDmap = world.plateID
    world_orogeny = fill(0.,nx,ny)
    for ixworld in 1:nx
        for iyworld in 1:ny
            if world.crust_type[ixworld,iyworld] >= continent_type
                plateID = plateIDmap[ixworld,iyworld]
                mplate = resolve_rotation_matrix!(plateID,age)
                ixplate,iyplate = nearest_plateij(mplate,ixworld,iyworld)
                # the position on the rotated plate field
                # assume zero rotation at present-day, so look for footprint here
                if footprint[ixplate,iyplate] > 0.
                    world_orogeny[ixworld,iyworld] = footprint[ixplate,iyplate]
                    world.crust_type[ixworld,iyworld] = uplifted_continent_type
                end
            end
        end
    end
    return world_orogeny
end

# grid navigation utilities
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
function nearest_other_plate_ij(iplate1,jplate1,mplate1,mplate2)
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
    #Threads.@threads
    for ixplate in 1:nx
        #Threads.@threads
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
        #Threads.@threads
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
