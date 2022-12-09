function increment_plate_age()   # updates the age field in plate grids
    Threads.@threads for plateID in world.plateIDlist
        plate = plates[plateID]
        for iplate in 1:nx
            for jplate in 1:ny
                if plate.crust_type[iplate,jplate] != not_at_surface
                    plate.crust_age[iplate,jplate] += main_time_step
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
    timestepoffset = 0.001
    filteredrotations = retrieve_rotations_for_plate(plateID)
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
        if age - timestepoffset > filteredrotations[end].age + main_time_step
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
    timestepoffset = 0.001
    hackedage = age + timestepoffset
    if plateID == -1
        return norotation
    end
    if haskey(plates,plateID) == false
        #println("creating plate in resolve_rotation_matrix ",plateID)
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
    for plateID in world.plateIDlist
        if haskey(plates,plateID) == false
            println("    creating plate in resolve_rotation_matrices ",plateID)
            plates[plateID] = create_blank_plate(plateID)
            resolve_rotation_matrix!(plateID,world.age)
            initial_mask_plate!(plates[plateID])
    
        end
        if plates[plateID].resolvetime != world.age
            if plates[plateID].resolvetime != -1
                mplateold = plates[plateID].rotationmatrix
                mplatenew = resolve_rotation_matrix!(plateID)
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
    # leaves world.geomorphology and surface_type_changes alone

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
                if plate.crust_type[iplate,jplate] == not_at_surface # new crust
                    # the closest plate point to the projection of the world point
                    # thinks it is outside of the exposed plate.
                    # it will get filled in the plate grid in remask_plates()
                    world.crust_age[iworld,jworld] = 0.
                    world.crust_type[iworld,jworld] = ocean_crust
                    world.crust_thickness[iworld,jworld] = ocean_crust_h0
                    world.crust_density[iworld,jworld] = rho_ocean_crust
                    world.geomorphology[iworld,jworld] = pelagic_seafloor
                    world.tectonics[iworld,jworld] = new_ocean_crust
                    world.sediment_thickness[iworld,jworld] = 0.
                    world.sediment_surface_fractions[iworld,jworld,:] .= initial_sediment_fractions
                    world.sediment_layer_thickness[iworld,jworld,:] .= 0.
                    world.sediment_layer_fractions[iworld,jworld,:,:] .= 0.
                else
                    world.crust_type[iworld,jworld] = plate.crust_type[iplate,jplate]
                    world.crust_age[iworld,jworld] = plate.crust_age[iplate,jplate]
                    world.crust_thickness[iworld,jworld] =
                        plate.crust_thickness[iplate,jplate]
                    world.crust_density[iworld,jworld] =
                        plate.crust_density[iplate,jplate]
                    world.geomorphology[iworld,jworld] = 
                        plate.geomorphology[iplate,jplate]
                    if plate.geomorphology[iplate,jplate] == exposed_basement
                        world.geomorphology[iworld,jworld] = sedimented_land
                    end # start over for the land diffusion calculation

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
function initial_mask_plate!(plate)
    for iplate in 1:nx
        for jplate in 1:ny
            iworld, jworld = nearest_worldij(plate.rotationmatrix,iplate,jplate)
            if world.plateID[iworld,jworld] == plate.plateID
                plate.crust_type[iplate,jplate] = world.crust_type[iworld,jworld]
                plate.crust_age[iplate,jplate] = 0.
                plate.crust_thickness[iplate,jplate] =
                    world.crust_thickness[iworld,jworld]
                if world.crust_type[iworld,jworld] == ocean_crust
                    plate.crust_density[iplate,jplate] = rho_ocean_crust
                    plate.geomorphology[iplate,jplate] = pelagic_seafloor
                    plate.sediment_thickness[iplate,jplate] = 
                        world.sediment_thickness[iworld,jworld]
                    plate.sediment_layer_thickness[iplate,jplate,:] = 
                        world.sediment_layer_thickness[iworld,jworld,:]
                    plate.sediment_surface_fractions[iplate,jplate,:] = 
                        world.sediment_surface_fractions[iworld,jworld,:]
                    plate.sediment_layer_fractions[iplate,jplate,:,:] = 
                        world.sediment_layer_fractions[iworld,jworld,:,:]
                else # must be a continent
                    plate.crust_density[iplate,jplate] = rho_continent_crust
                    plate.geomorphology[iplate,jplate] = sedimented_land
                    plate.sediment_thickness[iplate,jplate] = 
                        world.sediment_thickness[iworld,jworld]
                    plate.sediment_surface_fractions[iplate,jplate,:] .=  
                        world.sediment_surface_fractions[iworld,jworld,:]
                end
            else # world.plateID[iworld,jworld] != plate.plateID
                delete_plate_point!(plate,iplate,jplate)
            end
        end
    end
end

function remask_plate!( plate, subduction_footprint )
    # initialize, update, and delete plate points.
    # finds grid points of subduction and crust creation on plate fields.
    # imports crust_thickness as it stands in world grid.

    age = world.age
    if plate.resolvetime != age
        resolve_rotation_matrix!(plate.plateID,age)
    end

    iplates = fill(0,nx,ny); jplates = fill(0,nx,ny)
    rotated_crust_type = fill(0.,nx,ny)  # plate projected onto world grid
    cropped_rotated_crust_type = fill(0.,nx,ny)   #  still outcrop
    subducted_crust_type = fill(0.,nx,ny)     # not there any more

    rotated_sediment_thickness = fill(0.,nx,ny)
    cropped_rotated_sediment_thickness = fill(0.,nx,ny)
    subducted_sediment_thickness = fill(0.,nx,ny)

    rotated_sediment_surface_fractions = fill(0.,nx,ny,n_sediment_types)
    cropped_rotated_sediment_surface_fractions = fill(0.,nx,ny,n_sediment_types)
    subducted_sediment_surface_fractions = fill(0.,nx,ny,n_sediment_types)

    rotated_sediment_layer_thickness = fill(0.,nx,ny,n_sediment_time_bins)
    cropped_rotated_sediment_layer_thickness = fill(0.,nx,ny,n_sediment_time_bins)
    subducted_sediment_layer_thickness = fill(0.,nx,ny,n_sediment_time_bins)

    rotated_sediment_layer_fractions = fill(0.,nx,ny,n_sediment_types,n_sediment_time_bins)
    cropped_rotated_sediment_layer_fractions = fill(0.,nx,ny,n_sediment_types,n_sediment_time_bins)
    subducted_sediment_layer_fractions = fill(0.,nx,ny,n_sediment_types,n_sediment_time_bins)

    initial_plate_crust_areas = [0.,0.] # land, ocean
    rotated_plate_crust_areas = [0.,0.]
    subducted_plate_crust_areas = [0.,0.]
    cropped_rotated_plate_crust_areas = [0.,0.]
    initial_plate_sediment_volumes = fill(0.,n_sediment_types)
    rotated_plate_sediment_volumes = fill(0.,n_sediment_types)
    subducted_plate_sediment_volumes = fill(0.,n_sediment_types)
    cropped_rotated_plate_sediment_volumes = fill(0.,n_sediment_types)
    updated_plate_sediment_volumes = fill(0.,n_sediment_types)

    area_imbalances = [0.,0.]
    inv_imbalances = fill(0.,n_sediment_types)
    if enable_remask_plate_diagnostics 
        initial_plate_crust_areas = [ volume_field(eq_mask(plate.crust_type,continent_crust)), 
            volume_field(eq_mask(plate.crust_type,ocean_crust)) ]
        initial_plate_sediment_volumes = plate_sediment_inventories( plate )[1:end]
    end

    for iworld in 1:nx # project into world grid, carve off subducted and init new pts
        for jworld in 1:ny 
            iplates[iworld,jworld],jplates[iworld,jworld] = 
                nearest_plateij(plate.rotationmatrix,iworld,jworld)
            iplate = iplates[iworld,jworld]; jplate = jplates[iworld,jworld]
            # entire plate rotated into current world position
            rotated_crust_type[iworld, jworld] = plate.crust_type[iplate,jplate]
            rotated_sediment_thickness[iworld, jworld] = 
                plate.sediment_thickness[iplate,jplate]
            rotated_sediment_surface_fractions[iworld, jworld,:] = 
                plate.sediment_surface_fractions[iplate,jplate,:]
            rotated_sediment_layer_thickness[iworld, jworld,:] = 
                plate.sediment_layer_thickness[iplate,jplate,:]
            rotated_sediment_layer_fractions[iworld, jworld,:,:] = 
                plate.sediment_layer_fractions[iplate,jplate,:,:]
            # entire world grids dividing plate points into surface exposed vs not
            if world.plateID[iworld,jworld] == plate.plateID
                cropped_rotated_crust_type[iworld, jworld] = 
                    rotated_crust_type[iworld, jworld]
                cropped_rotated_sediment_thickness[iworld, jworld] = 
                    rotated_sediment_thickness[iworld, jworld]
                cropped_rotated_sediment_surface_fractions[iworld, jworld,:] = 
                    rotated_sediment_surface_fractions[iworld, jworld,:]
                cropped_rotated_sediment_layer_thickness[iworld, jworld,:] = 
                    rotated_sediment_layer_thickness[iworld, jworld,:]
                cropped_rotated_sediment_layer_fractions[iworld, jworld,:,:] = 
                    rotated_sediment_layer_fractions[iworld, jworld,:,:]
            else  
                subducted_crust_type[iworld, jworld] = 
                    rotated_crust_type[iworld, jworld]
                subducted_sediment_thickness[iworld, jworld] = 
                    rotated_sediment_thickness[iworld, jworld]
                subducted_sediment_surface_fractions[iworld, jworld,:] = 
                    rotated_sediment_surface_fractions[iworld, jworld,:]
                subducted_sediment_layer_thickness[iworld, jworld,:] = 
                    rotated_sediment_layer_thickness[iworld, jworld,:]
                subducted_sediment_layer_fractions[iworld, jworld,:,:] = 
                    rotated_sediment_layer_fractions[iworld, jworld,:,:]
            end
        end
    end
    # compute subduction volumes for mass balance calculations
    subducted_land_sediment_volumes = fill(0.,n_sediment_types)
    subducted_ocean_sediment_volumes = fill(0.,n_sediment_types)
    for i_sedtype in 1:n_sediment_types
        subducted_land_sediment_volumes[i_sedtype] = 
            volume_field( subducted_sediment_thickness .* 
                eq_mask( subducted_crust_type, continent_crust ) .* 
                subducted_sediment_surface_fractions[:,:,i_sedtype] )
        for i_bin in 1:n_sediment_time_bins
            subducted_ocean_sediment_volumes[i_sedtype] += 
                volume_field( subducted_sediment_layer_thickness[:,:,i_bin] .* 
                    eq_mask( subducted_crust_type, ocean_crust ) .* 
                    subducted_sediment_layer_fractions[:,:,i_sedtype,i_bin] )
        end
    end
    # adjust points in plate grids, init + deleted
    deleted_sediment_volumes = fill(0.,n_sediment_types)
    for iplate in 1:nx
        for jplate in 1:ny
            iworld, jworld = nearest_worldij(plate.rotationmatrix,iplate,jplate)
            if world.plateID[iworld,jworld] == plate.plateID
                if plate.crust_type[iplate,jplate] <= not_at_surface # new crust
                    plate.crust_type[iplate,jplate] = world.crust_type[iworld,jworld]
                    plate.crust_age[iplate,jplate] = 0.
                    plate.crust_thickness[iplate,jplate] =
                        world.crust_thickness[iworld,jworld]
                    if world.crust_type[iworld,jworld] == ocean_crust
                        plate.tectonics[iplate,jplate] = new_ocean_crust 
                        plate.crust_density[iplate,jplate] = rho_ocean_crust
                        plate.geomorphology[iplate,jplate] = pelagic_seafloor
                        plate.sediment_thickness[iplate,jplate] = 0. # initial_ocean_sediment_thickness
                        plate.sediment_layer_thickness[iplate,jplate,:] .= 0.
                        plate.sediment_layer_thickness[iplate,jplate,1] = 0. # initial_ocean_sediment_thickness
                        plate.sediment_surface_fractions[iplate,jplate,:] .= initial_sediment_fractions
                        plate.sediment_layer_fractions[iplate,jplate,:,1] .= initial_sediment_fractions
                        plate.sediment_layer_fractions[iplate,jplate,:,2:end] .= 0.
                    else # must be a continent
                        plate.tectonics[iplate,jplate] = new_continent_crust 
                        plate.crust_density[iplate,jplate] = rho_continent_crust
                        plate.geomorphology[iplate,jplate] = sedimented_land
                        plate.sediment_thickness[iplate,jplate] = initial_land_sediment_thickness
                        plate.sediment_surface_fractions[iplate,jplate,:] .= initial_sediment_fractions
                    end
                end
            else # world.plateID[iworld,jworld] != plate.plateID
                if plate.crust_type[iplate,jplate] == continent_crust
                    plate.tectonics[iplate,jplate] = subducting_continent_crust 
                    deleted_sediment_volumes += plate.sediment_thickness[iplate,jplate] .*
                        plate.sediment_surface_fractions[iplate,jplate,:] .* areabox[jplate]
                end
                if plate.crust_type[iplate,jplate] == ocean_crust
                    plate.tectonics[iplate,jplate] = subducting_ocean_crust 
                    for i_bin in 1:n_sediment_time_bins
                        deleted_sediment_volumes += 
                            plate.sediment_layer_thickness[iplate,jplate,i_bin] .*
                            plate.sediment_layer_fractions[iplate,jplate,:,i_bin ] .* areabox[jplate]
                    end
                end
                delete_plate_point!(plate,iplate,jplate)
            end
        end
    end
    for iworld in 1:nx # subduction bookkeeping on world grid 
        for jworld in 1:ny 
            iplate = iplates[iworld,jworld]; jplate = jplates[iworld,jworld]
            if plate.crust_type[iplate,jplate] > not_at_surface && # used to exist
                world.plateID[iworld,jworld] != plate.plateID # but doesnt anymore
                
                if subducted_crust_type[iworld,jworld] == ocean_crust
                    subduction_footprint[iworld,jworld] += 
                        areabox[jplate] / areabox[jworld] / sub_time_step

                    # accumulated diagnostics in world grid, although world is changing 
                    accum_diag("ocean_subduct_plate_thickness",iworld,jworld,
                        areabox[jplate] / areabox[jworld] / sub_time_step )
                    accum_diag("ocean_subduct_age_thickness",iworld,jworld,
                        plate.crust_age[iplate,jplate] * 
                        areabox[jplate] / areabox[jworld] / sub_time_step )
                    for ibin in 1:n_sediment_time_bins # sediment record in the layers
                        for i_sedtype in 1:n_sediment_types
                            flux_amt = subducted_sediment_layer_thickness[iworld,jworld,ibin] *
                                subducted_sediment_layer_fractions[iworld,jworld,i_sedtype,ibin]
                            accum_frac_diag("ocean_subduct_sediment_fraction_thickness",iworld,jworld,
                                i_sedtype, flux_amt  * 
                                areabox[jplate] / areabox[jworld] / sub_time_step )
                        end
                    end 
                elseif subducted_crust_type[iworld,jworld] == continent_crust
                    accum_diag("continent_subduct_plate_thickness",iworld,jworld, 
                        areabox[jplate] / areabox[jworld] / sub_time_step )
                    for i_sedtype in 1:n_sediment_types
                        flux_amt = subducted_sediment_thickness[iworld,jworld] *
                            subducted_sediment_surface_fractions[iworld,jworld,i_sedtype]
                        accum_frac_diag("continent_subduct_sediment_fraction_thickness",iworld,jworld,
                            i_sedtype, flux_amt  * 
                            areabox[jplate] / areabox[jworld] / sub_time_step )
                    end
                end
            end
        end
    end
    if enable_remask_plate_diagnostics
        subducted_plate_crust_areas = [ 
            volume_field(eq_mask(subducted_crust_type,continent_crust)), 
            volume_field(eq_mask(subducted_crust_type,ocean_crust)) ]
        subducted_plate_sediment_volumes = subducted_land_sediment_volumes .+ 
            subducted_ocean_sediment_volumes

        rotated_plate_crust_areas = [ 
            volume_field(eq_mask(rotated_crust_type,continent_crust)), 
            volume_field(eq_mask(rotated_crust_type,ocean_crust)) ]
        cropped_rotated_plate_crust_areas = [ 
            volume_field(eq_mask(cropped_rotated_crust_type,continent_crust)), 
            volume_field(eq_mask(cropped_rotated_crust_type,ocean_crust)) ]

        rotated_plate_sediment_volumes = 
            volume_fields( eq_mask(rotated_crust_type,continent_crust) .*
            rotated_sediment_thickness .* rotated_sediment_surface_fractions )
        cropped_rotated_plate_sediment_volumes = 
            volume_fields( eq_mask(cropped_rotated_crust_type,continent_crust) .*
            cropped_rotated_sediment_thickness .* cropped_rotated_sediment_surface_fractions )
        for i_bin in 1:n_sediment_time_bins
            rotated_plate_sediment_volumes += volume_fields( 
                eq_mask(rotated_crust_type,ocean_crust) .*
                rotated_sediment_layer_thickness[:,:,i_bin] .* 
                rotated_sediment_layer_fractions[:,:,:,i_bin] )
            cropped_rotated_plate_sediment_volumes += volume_fields( 
                eq_mask(cropped_rotated_crust_type,ocean_crust) .*
                cropped_rotated_sediment_layer_thickness[:,:,i_bin] .* 
                cropped_rotated_sediment_layer_fractions[:,:,:,i_bin] )
        end
        final_plate_crust_areas = [ volume_field(eq_mask(plate.crust_type,continent_crust)), 
            volume_field(eq_mask(plate.crust_type,ocean_crust)) ]
        final_plate_sediment_volumes = plate_sediment_inventories( plate )[1:end]

        println();println("Plate ", plate.plateID)
        println("   inital area ", initial_plate_crust_areas)
        println("  rot to world ", rotated_plate_crust_areas) 
        println("          slop ", rotated_plate_crust_areas .- initial_plate_crust_areas)
        println(" subd on world ", subducted_plate_crust_areas)
        println("final on world ", cropped_rotated_plate_crust_areas)
        world_area_bal = cropped_rotated_plate_crust_areas .- 
            rotated_plate_crust_areas .+
            subducted_plate_crust_areas
        println("world area bal ", world_area_bal )
        println("final on plate ", final_plate_crust_areas)
        println("   inital seds ", initial_plate_sediment_volumes)
        println("  rot to world ", rotated_plate_sediment_volumes) 
        println("          slop ", rotated_plate_sediment_volumes .- initial_plate_sediment_volumes)
        println(" subd on world ", subducted_plate_sediment_volumes)
        println("final on world ", cropped_rotated_plate_sediment_volumes)
        world_sed_bal = cropped_rotated_plate_sediment_volumes .- 
            rotated_plate_sediment_volumes .+ 
            subducted_plate_sediment_volumes 
        println(" world sed bal ", world_sed_bal)
        println("final on plate ", final_plate_sediment_volumes)
        println(" subd on plate ", deleted_sediment_volumes)
        println("     plate bal ", final_plate_sediment_volumes .- 
            initial_plate_sediment_volumes .+
            deleted_sediment_volumes)
    end 

    if enable_remask_plate_diagnostics
        updated_plate_crust_areas = [ 
            volume_field(eq_mask(plate.crust_type,continent_crust)), 
            volume_field(eq_mask(plate.crust_type,ocean_crust))]
        updated_plate_sediment_volumes = plate_sediment_inventories( plate )[1:end]

        area_imbalances = 0.; inv_imbalances = 0.
        for i_crusttype in 1:2
            area_imbalances += world_area_bal[i_crusttype]
        end
        for i_sedtype in 1:n_sediment_types
            inv_imbalances += world_sed_bal[i_sedtype]
        end
        return subducted_land_sediment_volumes, subducted_ocean_sediment_volumes, 
            area_imbalances,inv_imbalances
    else
        return subducted_land_sediment_volumes, subducted_ocean_sediment_volumes
    end
end
function remask_plates()   # subduction and crust creation
    n_plates = length(world.plateIDlist)
    area_imbalance_list = fill(0.,n_plates)
    inv_imbalance_list = fill(0.,n_plates)
 
    subducted_land_sediment_volumes = fill(0.,n_sediment_types)
    subducted_ocean_sediment_volumes = fill(0.,n_sediment_types)
    subducted_land_sediment_volumes_list = fill(0.,n_sediment_types,n_plates)
    subducted_ocean_sediment_volumes_list = fill(0.,n_sediment_types,n_plates)
    subduction_footprint = fill(0.,nx,ny)

    Threads.@threads for i_plate = 1:n_plates

        plateID = world.plateIDlist[i_plate]  
        if haskey(plates,plateID) == false
            plates[plateID] = create_blank_plate(plateID)
            println("emergency plate creation in remask_plates ", plateID)
        end
        if enable_remask_plate_diagnostics
            subducted_land_sediment_volumes_list[:,i_plate], 
                subducted_ocean_sediment_volumes_list[:,i_plate],
                area_imbalance_list[i_plate], inv_imbalance_list[i_plate] = 
                    remask_plate!( plates[plateID], subduction_footprint )            
        else
            subducted_land_sediment_volumes_list[:,i_plate], 
                subducted_ocean_sediment_volumes_list[:,i_plate] =
                    remask_plate!( plates[plateID], subduction_footprint ) 
        end
    end
    for (plateID,plate) in plates
        if !(plateID in world.plateIDlist)
            delete!(plates,plateID)
        end
    end
    subducted_land_sediment_volumes = fill(0.,n_sediment_types)
    subducted_ocean_sediment_volumes = fill(0.,n_sediment_types)
    for i_plate = 1:n_plates
        subducted_land_sediment_volumes .+= 
            subducted_land_sediment_volumes_list[:,i_plate]
        subducted_ocean_sediment_volumes .+= 
            subducted_ocean_sediment_volumes_list[:,i_plate]
    end
    if enable_remask_plate_diagnostics
        for i_plate in 1:n_plates
            println("plate ", world.plateIDlist[i_plate]," ",area_imbalance_list[i_plate])
        end
        println("plotting area imbalances by plate")
        scene = GLMakie.lines(area_imbalance_list)
        display(scene)
        sleep(2)
        println("plotting inv imbalances by plate")
        for i_plate in 1:n_plates
            println("plate ", world.plateIDlist[i_plate]," ",inv_imbalance_list[i_plate])
        end
        scene = GLMakie.lines(inv_imbalance_list)
        display(scene)
        sleep(2)
    end
    
    return subducted_land_sediment_volumes, subducted_ocean_sediment_volumes, subduction_footprint
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
function delete_old_plates()
    for (plateID,plate) in plates
        if plateID in world.plateIDlist == false
            delete!( plates, plateID )
            println("deleted old plate ",plateID)
        end
    end
end

function update_changing_plateIDs(oldIDmap,newIDmap)
    changingplates = 
        auto_filter_changing_platelists(oldIDmap,newIDmap)
    update_changing_plateIDs(oldIDmap,newIDmap,changingplates)
end
function update_changing_plateIDs(oldIDmap,newIDmap,changingplates)
    newtime = world.age
    
    #Threads.@threads 
    needschangingglobalmask = fill(0,nx,ny)
    for platepair in changingplates
        oldplateID = platepair[1]; newplateID = platepair[2]
        initial_inv = 0.
        if enable_changing_plateID_diagnostics
            initial_inv = volume_field(plates[oldplateID].sediment_thickness)
            if haskey(plates,newplateID)
                initial_inv += volume_field(plates[newplateID].sediment_thickness)
            end
        end

        needschangingworldmask = update_two_plates!(oldplateID,newplateID,oldIDmap,newIDmap)

        if enable_changing_plateID_diagnostics
            final_inv = volume_field(plates[oldplateID].sediment_thickness) +
                volume_field(plates[newplateID].sediment_thickness)
            println("ID change ",[oldplateID,newplateID],
                [initial_inv,final_inv,final_inv - initial_inv])
            #scene = plot_field( needschangingworldmask )
            #plot_add_continent_outlines!( scene )
            #display(scene)
            #sleep(3)
        end
    end
    return
end
function looks_like_a_real_ID_change( needschangingworldmask )
    its_a_big_blob = false
    if length( get_blob_onion_layers(needschangingworldmask ) ) > 5
        its_a_big_blob = true
    end
    return its_a_big_blob
end

function update_two_plates!(oldplateID,newplateID,oldIDmap,newIDmap)
    # called when one or more grid points change their plateIDs
    newtime = world.age; oldtime = newtime + main_time_step
    #println("updating ", oldplateID," ",newplateID," ",oldtime," ",newtime)
    if haskey(plates,oldplateID) == false
        println("    creating plate to migrate data from ", oldplateID)
        plates[oldplateID] = create_blank_plate(newplateID)
    end
    oldplate = plates[oldplateID]
    if haskey(plates,newplateID) == false
        println("    creating plate to migrate data to ", newplateID)
        plates[newplateID] = create_blank_plate(newplateID)
        resolve_rotation_matrix!(newplateID,newtime)
        initial_mask_plate!(plates[newplateID])
    end
    newplate = plates[newplateID]

    oldIDmask = eq_mask(oldIDmap,oldplateID)
    # grid points in world grid of old plateID number at old time step
    newIDmask = eq_mask(newIDmap,newplateID)
    # points in new world grid where the new plate ID is found
    needschangingworldmask = oldIDmask .* newIDmask # fill(0,nx,ny)

    # points in world grid where old ID changes to new ID
    firstoldtime,lastoldtime = first_last_appearance(oldplateID)
    #firstnewtime,lastnewtime = first_last_appearance(newplateID)
    switchovertime = newtime # use rotations for newtime if possible
    if lastoldtime > newtime
        switchovertime = lastoldtime  # or back up to end time of old plate
    end
    newplate_switchrotation = resolve_rotation_matrix!(newplateID,switchovertime)
    needschangingnewplatemask = projected_plate_maskfield!(needschangingworldmask,
        newplate_switchrotation)
    oldplate_switchrotation = resolve_rotation_matrix!(oldplateID,switchovertime)
    needsdeletingoldplatemask = projected_plate_maskfield!(needschangingworldmask,
        oldplate_switchrotation)
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

                    accum_diag("IDchange_plate_area",ixworld,iyworld,
                        areabox[iynewplate] / areabox[iyworld] )
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
                accum_diag("IDchange_plate_area",ixworld,iyworld,
                    areabox[iyplate]/areabox[iyworld])
                oldplate.tectonics[ixplate,iyplate] = deleted_plateID
            end
        end
    end
    return needschangingworldmask
end
function copy_plate_point_plate12coord!(newplate,oldplate,ixoldplate,iyoldplate,ixnewplate,iynewplate)
    if oldplate.crust_type[ixoldplate,iyoldplate] == not_at_surface #||
        #        oldplate.crust_type[ixold,iyold] == ocean_crust  # moving ridge
        # edge cases, messy spreading
        newplate.crust_type[ixnewplate,iynewplate] = ocean_crust
        newplate.crust_age[ixnewplate,iynewplate] = 0.
        newplate.crust_thickness[ixnewplate,iynewplate] = ocean_crust_h0
        newplate.crust_density[ixnewplate,iynewplate] = rho_ocean_crust
        newplate.geomorphology[ixnewplate,iynewplate] = pelagic_seafloor
        newplate.tectonics[ixnewplate,iynewplate] = new_ocean_crust

        newplate.sediment_thickness[ixnewplate,iynewplate] = initial_ocean_sediment_thickness
        newplate.sediment_layer_thickness[ixnewplate,iynewplate,:] .= 0.
        newplate.sediment_layer_thickness[ixnewplate,iynewplate,1] = initial_ocean_sediment_thickness
        newplate.sediment_surface_fractions[ixnewplate,iynewplate,:] .= initial_sediment_fractions
        newplate.sediment_layer_fractions[ixnewplate,iynewplate,:,1] .= initial_sediment_fractions
        newplate.sediment_layer_fractions[ixnewplate,iynewplate,:,2:end] .= 0.
    else
        newplate.crust_type[ixnewplate,iynewplate] = oldplate.crust_type[ixoldplate,iyoldplate]
        newplate.crust_age[ixnewplate,iynewplate] = oldplate.crust_age[ixoldplate,iyoldplate]
        newplate.crust_thickness[ixnewplate,iynewplate] = oldplate.crust_thickness[ixoldplate,iyoldplate]
        newplate.crust_density[ixnewplate,iynewplate] = oldplate.crust_density[ixoldplate,iyoldplate]
        newplate.geomorphology[ixnewplate,iynewplate] = oldplate.geomorphology[ixoldplate,iyoldplate]
        newplate.sediment_thickness[ixnewplate,iynewplate] = 
            oldplate.sediment_thickness[ixoldplate,iyoldplate]
        newplate.sediment_surface_fractions[ixnewplate,iynewplate,:] .= 
            oldplate.sediment_surface_fractions[ixoldplate,iyoldplate,:]
        newplate.sediment_layer_thickness[ixnewplate,iynewplate,:] .= 
            oldplate.sediment_layer_thickness[ixoldplate,iyoldplate,:]
        newplate.sediment_layer_fractions[ixnewplate,iynewplate,:,:] .= 
            oldplate.sediment_layer_fractions[ixoldplate,iyoldplate,:,:]
    end
end
function delete_plate_point!(oldplate,ixold,iyold)
    oldplate.crust_type[ixold,iyold] = not_at_surface
    oldplate.crust_age[ixold,iyold] = 0.
    oldplate.crust_thickness[ixold,iyold] = 0.
    oldplate.crust_density[ixold,iyold] = 0.
    oldplate.sediment_thickness[ixold,iyold] = 0.
    oldplate.sediment_surface_fractions[ixold,iyold,:] .= 0.
    oldplate.sediment_layer_thickness[ixold,iyold,:] .= 0.
    oldplate.sediment_layer_fractions[ixold,iyold,:,:] .= 0.
end
function auto_filter_changing_platelists(oldIDmap,newIDmap)
    #oldIDmap = world.plateID
    oldtime = world.age + main_time_step
    fp = 0
    changepairs = []
    for ixworld in 1:nx # 1+fp:nx-fp
        for iyworld in 1:ny # 1+fp:ny-fp
            if oldIDmap[ixworld,iyworld] != newIDmap[ixworld,iyworld]

                oldplateID = oldIDmap[ixworld,iyworld]
                newplateID = newIDmap[ixworld,iyworld]
                key = [ oldplateID, newplateID ]
                push!(changepairs,key)

            end
        end
    end
    # list of large blobs of gridpoints in successive world grids where the ID changes
    changepairs = unique(changepairs)
    filtered_changepairs = []
    needschangingglobalmask = fill(0,nx,ny)
    for changepair in changepairs
        oldplateID = changepair[1]; newplateID = changepair[2]
        oldIDmask = eq_mask(oldIDmap,oldplateID)
        # grid points in world grid of old plateID number at old time step
        newIDmask = eq_mask(newIDmap,newplateID)
        # points in new world grid where the new plate ID is found
        needschangingworldmask = oldIDmask .* newIDmask # fill(0,nx,ny)
        nskins = length( get_blob_onion_layers(needschangingworldmask ) )
        if nskins > 2 # do a plate transplant
            push!(filtered_changepairs,changepair)
            
            if enable_watch_plate_transplants || enable_save_plate_transplant_images
                println("plate transplant ",changepair," ",nskins)
                scene = plot_field( needschangingworldmask,-1.,1. )
                plot_add_continent_outlines!( scene )
                plot_add_timestamp!(scene,world.age,-180,-105) 
                #plot_add_plate_boundaries!( scene )
                transaction = string(oldplateID) * " change to " * string(newplateID)
                text!(scene,transaction, position = (0,95),textsize=15)
                if enable_watch_plate_transplants
                    display(scene)
                    sleep(2)
                end
                if enable_save_plate_transplant_images
                    global save_plate_transplants_image_number += 1
                    image_file_name = "img." * lpad(save_plate_transplants_image_number,4,"0") * ".png"
                    filename = base_directory * "/" * output_directory * output_tag * 
                        "/" * animation_directory * "/plate_transplants/" * image_file_name
                    Makie.save( filename, scene )
                end
            end
        end 
    end
    #=if enable_watch_plate_transplants
        println("plotting plate transplants, red=yes, purple=no")
        scene = plot_field( needschangingglobalmask,-1.,1. )
        plot_add_continent_outlines!( scene )
        #plot_add_plate_boundaries!( scene )
        #text!(scene,"ID change", position = (0,95),textsize=15)
        display(scene)
        sleep(2)
    end =#
    return filtered_changepairs
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
function apply_tectonics_changes_to_plates( uplift_rate )
    uplift_meters = uplift_rate * # meters / Myr
        sub_time_step

    Threads.@threads for plateID in world.plateIDlist
        plate = plates[plateID]
        for iplate in 1:nx
            for jplate in 1:ny
                iworld, jworld = nearest_worldij(plate.rotationmatrix,iplate,jplate)
                if world.plateID[iworld,jworld] == plateID
                    
                    plates[plateID].crust_age[iplate,jplate] = world.crust_age[iworld,jworld]
                    
                    if world.crust_type[iworld,jworld] == continent_crust &&
                        uplift_meters[iworld,jworld] > 0.
                        #error(iplate," ",jplate)
                        plates[plateID].crust_type[iplate,jplate] = continent_crust
                        plates[plateID].crust_thickness[iplate,jplate] =
                            world.crust_thickness[iworld,jworld]   
                        plates[plateID].crust_density[iplate,jplate] =
                            world.crust_density[iworld,jworld]   
                    end
                end
            end
        end
    end
end
function apply_geomorphology_changes_to_plates()
    Threads.@threads for plateID in world.plateIDlist
        plate = plates[plateID]
        for iplate in 1:nx
            for jplate in 1:ny
                iworld, jworld = nearest_worldij(plate.rotationmatrix,iplate,jplate)
                if world.plateID[iworld,jworld] == plateID
                    plates[plateID].crust_type[iplate,jplate] =
                        world.crust_type[iworld,jworld]
                    plates[plateID].crust_age[iplate,jplate] =
                        world.crust_age[iworld,jworld]
                    plates[plateID].crust_thickness[iplate,jplate] =
                        world.crust_thickness[iworld,jworld]                    
                    plates[plateID].crust_density[iplate,jplate] =
                        world.crust_density[iworld,jworld]
                    plates[plateID].geomorphology[iplate,jplate] =
                        world.geomorphology[iworld,jworld]
                    plates[plateID].tectonics[iplate,jplate] =
                        world.tectonics[iworld,jworld]
                    plates[plateID].sediment_thickness[iplate,jplate] =
                        world.sediment_thickness[iworld,jworld]
                    plates[plateID].sediment_surface_fractions[iplate,jplate,:] =
                        world.sediment_surface_fractions[iworld,jworld,:]
                    plates[plateID].sediment_layer_thickness[iplate,jplate,:] =
                        world.sediment_layer_thickness[iworld,jworld,:]
                    plates[plateID].sediment_layer_fractions[iplate,jplate,:,:] =
                        world.sediment_layer_fractions[iworld,jworld,:,:]
                else
                    delete_plate_point!(plates[plateID],iplate,jplate)
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




