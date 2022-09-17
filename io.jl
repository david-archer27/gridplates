function write_field_csv(filename,field)
    f = open(filename,"w")
    for iy in 1:ny
        for ix in 1:nx-1
            @printf(f, "%0.2f,", field[ix,iy])
        end
        @printf(f,"%0.2f\n", field[nx,iy])
    end
    close(f)
    return
end
function read_my_csv(filename)
    f = open(filename)
    arr = fill(0,nx,ny)
    lines = readlines(f)
    for (iy,line) in enumerate(lines)
        iyflip = iy
        words = split(line,",")
        for (ix,word) in enumerate(words)
            value = parse(Int,word)
            arr[ix,iyflip] = value
        end
    end
    return arr
end
function read_globsed_xyz()
    f = open("../GlobSed-v3.xyz")
    field = fill(0.,nx,ny)
    lines = readlines(f)
    for iline in 1:length(lines)
        line = lines[iline]
        words = split(line)
        longt = parse(Float64,words[1])
        lat = parse(Float64,words[2])
        #println(longt, " ", lat)
        if longt == floor(longt) + 0.5 && lat == floor(lat) + 0.5
            ix = Int(longt - 0.5 + 181)
            iy = Int(lat - 0.5 + 91)
            thickness = parse(Float64,words[3])
            if thickness != thickness
                thickness = 0.
            end
            field[ix,iy] = thickness
            #println(longt, " ", lat," ",ix," ",iy," ",thickness)
        end
    end
    return field
end
function read_flip_csv(filename)
    f = open(filename)
    arr = fill(0.,nx,ny)
    lines = readlines(f)
    for (iy,line) in enumerate(lines)
        iyflip = ny - iy + 1
        words = split(line,",")
        for (ix,word) in enumerate(words)
            value = parse(Int,word)
            arr[ix,iyflip] = value
        end
    end
    return arr
end
function read_plateIDs()
    return read_plateIDs( world.age )
end
function read_plateIDs( age )
    integerage = Int(ceil(age))
    filename = base_directory * "/" * code_base_directory * "/platefiles/plateIDs." * string(integerage) * ".csv"
    plateIDmap = read_my_csv(filename)
    return plateIDmap
end
function read_continentIDs()
    return read_continentIDs( world.age )
end
function read_continentIDs( age )
    #    contIDmap = fill(0,nx,ny)
    #    return contIDmap
    integerage = Int(ceil(age))
    filename = base_directory * "/" * code_base_directory * "/contfiles/contIDs." * string(Int(integerage)) * ".csv"
    contIDmap = read_my_csv(filename)
    return contIDmap
end
function read_rotation_file(filename)
    f = open(filename)
    lines = readlines(f)
    rotations = []
    for line in lines
        words = split(line)
        id = parse(Int,words[1])
        parent = parse(Int,words[6])
        v = parse.(Float64,words[2:5])
        rotation = rotation_struct(id,v[1],v[2],v[3],v[4],parent)
        push!(rotations,rotation)
    end
    sortedrotations = sort(rotations,by=x->(x.id, x.age))
    return sortedrotations# , comments
end

# file io
function setup_working_directories()
    cd( base_directory )
    if output_directory * output_tag in readdir()
    else
        mkdir( output_directory * output_tag )
        println( "creating ", output_directory * output_tag )
    end
    cd( output_directory * output_tag )
    if world_directory in readdir()
    else
        mkdir( world_directory )
        println( "creating ", world_directory )
    end
    if plate_directory in readdir()
    else
        mkdir( plate_directory )
        println( "creating ", plate_directory )
    end
    
    if charts_directory in readdir()
    else
        mkdir( charts_directory )
        println( "creating ", charts_directory )
    end
    if animation_directory in readdir()
    else
        mkdir( animation_directory )
        println( "creating ", animation_directory )
    end
    cd( animation_directory )
    for plot_type in animation_directories
        if plot_type in readdir()
        else
            mkdir( plot_type )
            println( "creating ", plot_type  )
        end
    end
    if enable_save_plate_transplant_images

        if "plate_transplants" in readdir()
        else
            mkdir( "plate_transplants" )
            println("creating plate_transplants")
        end
    end
    cd( base_directory * "/" * code_base_directory )
end
function save_field(variable_name,field)
    timestamp = string(Int(ceil(world.age)))
    filename = "../outfiles/fields/" * variable_name * "." * timestamp * ".bson"
    rm(filename, force=true)
    BSON.@save filename field
end
function save_world()
    timestamp = string(Int(ceil(world.age)))
    # in world grid, only the age field is non-trivial to calculate
    filename = base_directory * "/" * output_directory * output_tag * "/" *
        world_directory * "/" * output_tag * ".world." * timestamp * ".bson"
    rm(filename, force=true)
    println("saving ", filename)
    BSON.@save filename world
    return
end
function read_world(age)
    timestamp = string(Int(ceil(age)))
    # in world grid, only the age field is non-trivial to calculate
    #filename = base_directory * output_directory * output_tag * "/" *
    #    world_subdirectory * "world." * output_tag * "." * timestamp * ".bson"
    filename = base_directory * "/" * output_directory * output_tag * "/" *
        world_directory * "/" * output_tag * ".world." * timestamp * ".bson"
    println("reading world ",filename)
    BSON.@load filename world
    return world
end
function save_plates()
    timestamp = string(Int(ceil(world.age)))
    for plateID in world.plateIDlist
        plate = plates[plateID]
        IDstring = string(Int(plate.plateID))
        filename = base_directory * "/" * output_directory * output_tag * "/" *
            plate_directory * "/" * output_tag * ".plate." * IDstring * "." * timestamp * ".bson"
        rm(filename, force=true)
        println("saving ", filename)
        BSON.@save filename plate
    end
    return
end
function save_plates_checkpoint()
    # in world grid, only the age field is non-trivial to calculate
    filename = base_directory * output_directory * output_tag *
        plate_subdirectory * output_tag * ".plates.checkpoint.bson"
    rm(filename, force=true)
    println("saving ", filename)
    BSON.@save filename plates
    return
end
function read_plates( )
    cd( base_directory * "/" * code_base_directory )
    age = world.age
    timestamp = string(Int(ceil(age)))
    for plateID in world.plateIDlist
        println("reading plate ", plateID)
        IDstring = string(Int(plateID))
        filename = base_directory * "/" * output_directory * output_tag * "/" *
            plate_directory * "/" * output_tag * ".plate." * IDstring * "." * timestamp * ".bson"
        BSON.@load filename plate
        if haskey(plates,plateID) == false
            println("creating plate ",plateID)
            plates[plateID] = create_blank_plate(plateID)
        end
        plates[plateID] = plate
    end
    return plates
end
function read_sedthick_xyz()
    f = open( base_directory * "/sedthick.xyz")
    lines = readlines(f)
    field = fill(0.,nx,ny)
    for line in lines
        words = split(line)
        longt = parse(Float64,words[1]); lat = parse(Float64,words[2])
        if abs( longt - floor( longt ) - 0.5 ) < 0.01 && abs( lat - floor( lat ) - 0.5 ) < 0.01
            ix = Int(longt + 0.5 ) + 180; iy = 90 + Int(lat + 0.5)
            if iy == 0 
                error("[longt,lat][ix,iy]",longt," ",lat," ",ix," ",iy)
            end
            if words[3] != "NaN"
                field[ix,iy] = parse(Float64,words[3])
            end
        end
    end
    return field
end
function logging_println(text,array,text2,array2)
    println(text,array,text2,array2)
    if log_IO != 0
        println(log_IO,text,array,text2,array2)
    end
end
function logging_println(text,array,text2)
    println(text,array,text2)
    if log_IO != 0
        println(log_IO,text,array,text2)
    end
end
function logging_println(text,array)
    println(text,array)
    if log_IO != 0
        println(log_IO,text,array)
    end
end
function logging_println(text)
    println(text)
    if log_IO != 0
        println(log_IO,text)
    end
end
function logging_println()
    println()
    if log_IO != 0
        println(log_IO,"")
    end
end
function generate_outline_from_mask(maskfield)
    mycontour = Contour.contour(xcoords,ycoords,maskfield,0.5)
    lines = Contour.lines(mycontour)
    if length(lines) > 0
        line = [1]
        xs,ys = Contour.coordinates(line)
    else
        xs = []; ys = []
    end
    return xs,ys
end
function line_sphere2cart(xline,yline)  # x y, not lat longt
    nxl = length(xline)
    x3d = zeros(nxl);y3d = zeros(nxl); z3d = zeros(nxl)
    for i in 1:length(xline)
        azimuth = xline[i] / 180. * pi    # longitude
        elevation =  ( yline[i] + 90. ) / 180. * pi # 0 at pole, pi/2 eq, pi pole/ 180. * pi  # latitude
        x3d[i] = 1.01 * sin(elevation) * cos(azimuth) # convert(AbstractFloat,)
        y3d[i] = 1.01 * sin(elevation) * sin(azimuth)
        z3d[i] = 1.01 * - cos(elevation)
    end
    return x3d,y3d,z3d
end
# Animations
function cleanup_image_files( )
    for plot_type in animation_directories
        cd( base_directory * "/" * output_directory * output_tag * "/" * 
            animation_directory * "/" * plot_type)
        file_list = readdir()
        for file in file_list
            if occursin( ".png",file ) 
                println("removing ", animation_directory * "/" * 
                plot_type * "/" * file )
                rm(file)
            end
        end
        srcfile = plot_type * "." * output_tag * ".mp4" 
        destfile = base_directory * "/" * output_directory * output_tag * "/" * 
            animation_directory * "/" * plot_type * "." * output_tag * ".mp4"
        println("moving ", srcfile, " ", destfile )
        mv( srcfile, destfile )
        cd( base_directory * "/" * output_directory * output_tag * "/" * 
            animation_directory )
            println("removing directory ", plot_type )
        rm( plot_type, recursive=true )
    end
end 
function animate_elevation()
    directory = base_directory * "/" * output_directory * output_tag * "/" * 
        animation_directory * "/" * "elevation" 
    cd( directory )
    starting_file_list = readdir()
    image_number = 0
    ages = [ animation_final_age ]
    for age in animation_initial_age: - time_step * animation_n_step : animation_final_age
        push!(ages,age)
    end
    for age in ages
        image_number += 1
        image_file_name = "img." * lpad(image_number,3,"0") * ".png"
        if image_file_name in starting_file_list
            println( "already done ", age," ", image_file_name )
        else
            println( "creating ", age, " ", directory * "/" * image_file_name )
            global world = read_world( age )
            scene = plot_elevation()
            plot_add_plate_boundaries!(scene)
            plot_add_orogenies!(scene)
            plot_add_continent_outlines!(scene)
            plot_add_timestamp!(scene,world.age,-180,-105) 
            plot_add_title!(scene,"Elevation L(0:5km)/O(-6:-2km)")
            Makie.save( image_file_name, scene )
        end
    end
    mp4_file =  "elevation." * output_tag * ".mp4"
    println("compiling ", mp4_file)
    rm( mp4_file,force=true)
    run(`ffmpeg -r 5 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p $mp4_file`)
    cd( base_directory * "/" * code_base_directory )
end
function animate_pct_CaCO3()
    directory = base_directory * "/" * output_directory * output_tag * "/" * 
        animation_directory * "/" * "pct_CaCO3" 
    cd( directory )
    starting_file_list = readdir()
    image_number = 0
    ages = [ animation_final_age ]
    for age in animation_initial_age: - time_step * animation_n_step : animation_final_age
        push!(ages,age)
    end
    for age in ages
        image_number += 1
        image_file_name = "img." * lpad(image_number,3,"0") * ".png"
        if image_file_name in starting_file_list
            println( "already done ", age," ", image_file_name )
        else
            println( "creating ", age, " ", directory * "/" * image_file_name )
            global world = read_world( age )
            scene = plot_field(world.sediment_surface_fractions[:,:,2] .* 100.,0.,100.)
            plot_add_plate_boundaries!(scene)
            plot_add_orogenies!(scene)
            plot_add_continent_outlines!(scene)
            plot_add_timestamp!(scene,world.age,-180,-105)    #scene = plot_field(world.freeboard,-5000.,5000)
            plot_add_title!(scene,"%CaCO3")
            Makie.save( image_file_name, scene )
        end
    end
    mp4_file =  "pct_CaCO3." * output_tag * ".mp4"
    println("compiling ", mp4_file)
    rm( mp4_file,force=true)
    run(`ffmpeg -r 5 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p $mp4_file`)
    cd( base_directory * "/" * code_base_directory )
end
function animate_sed_thickness()
    directory = base_directory * "/" * output_directory * output_tag * "/" * 
        animation_directory * "/" * "sed_thickness" 
    cd( directory )
    starting_file_list = readdir()
    image_number = 0
    ages = [ animation_final_age ]
    for age in animation_initial_age: - time_step * animation_n_step : animation_final_age
        push!(ages,age)
    end
    for age in ages
        image_number += 1
        image_file_name = "img." * lpad(image_number,3,"0") * ".png"
        if image_file_name in starting_file_list
            println( "already done ", age," ", image_file_name )
        else
            println( "creating ", age, " ", directory * "/" * image_file_name )
            global world = read_world( age )
            scene = plot_field(world.sediment_thickness ./ 1000.,0.,10.) # plot_sediment_thickness()
            plot_add_plate_boundaries!(scene)
            plot_add_orogenies!(scene)
            plot_add_continent_outlines!(scene)
            plot_add_timestamp!(scene,world.age,-180,-105)
            plot_add_title!(scene,"Sediment Thickness") # L(0:2km)/O(0:2km)")
            Makie.save( image_file_name, scene )
        end
    end
    mp4_file =  "sed_thickness." * output_tag * ".mp4"
    println("compiling ", mp4_file)
    rm( mp4_file,force=true)
    run(`ffmpeg -r 5 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p $mp4_file`)
    cd( base_directory * "/" * code_base_directory )
end
function animate_crust_age( )
    directory = base_directory * "/" * output_directory * output_tag * "/" * 
        animation_directory * "/" * "crust_age" 
    cd( directory )
    starting_file_list = readdir()
    image_number = 0
    ages = [ animation_final_age ]
    for age in animation_initial_age: - time_step * animation_n_step : animation_final_age
        push!(ages,age)
    end
    for age in ages
        image_number += 1
        image_file_name = "img." * lpad(image_number,3,"0") * ".png"
        if image_file_name in starting_file_list
            println( "already done ", age," ", image_file_name )
        else
            println( "creating ", age, " ", directory * "/" * image_file_name )
            global world = read_world( age )
            scene = plot_field(world.crust_age,0.,600.)
            plot_add_plate_boundaries!(scene)
            plot_add_orogenies!(scene)
            plot_add_continent_outlines!(scene)
            plot_add_timestamp!(scene,world.age,-180,-105)
            plot_add_title!(scene,"Crust age, Myr")
            Makie.save( image_file_name, scene )
        end
    end
    mp4_file =  "crust_age." * output_tag * ".mp4"
    println("compiling ", mp4_file)
    rm( mp4_file,force=true)
    run(`ffmpeg -r 5 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p $mp4_file`)
    cd( base_directory * "/" * code_base_directory )
end
function animate_scotese_elevation( )
    directory = base_directory * "/" * output_directory * output_tag * "/" * 
        animation_directory * "/" * "scotese_elevation" 
    cd( directory )
    starting_file_list = readdir()
    image_number = 0
    ages = [ animation_final_age ]
    for age in animation_initial_age: - time_step * animation_n_step : animation_final_age
        push!(ages,age)
    end
    for age in ages
        image_number += 1
        image_file_name = "img." * lpad(image_number,3,"0") * ".png"
        if image_file_name in starting_file_list
            println( "already done ", age," ", image_file_name )
        else
            println( "creating ", age, " ", directory * "/" * image_file_name )
            global world = read_world( age )
            scotese_elevation,scotese_age = nearest_scotese_elevation()
            scene = plot_two_fields(world.freeboard,scotese_elevation,-4000,4000)
            gridplatestimestamp = string(Int(floor(world.age))) * " Myr, " * get_geo_interval(age)
            scotesetimestamp = string(Int(floor(scotese_age))) * " Myr"
            text!(scene,gridplatestimestamp, position = (0,95),textsize=15)
            text!(scene, scotesetimestamp, position = (0,-110),textsize=15)
            plot_add_plate_boundaries!(scene)
            plot_add_orogenies!(scene)
            plot_add_continent_outlines!(scene)
            text!(scene,"GridPlates Elevation",position=(-180,95),textsize=20)
            text!(scene,"Scotese Reconstruction",position=(-180,-110),textsize=20)
            Makie.save( image_file_name, scene )
        end
    end
    mp4_file =  "scotese_elevation." * output_tag * ".mp4"
    println("compiling ", mp4_file)
    rm( mp4_file,force=true)
    run(`ffmpeg -r 5 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p $mp4_file`)
    cd( base_directory * "/" * code_base_directory )
end
function animate_all( )
    animate_elevation()
    animate_pct_CaCO3()
    animate_sed_thickness()
    animate_crust_age( )
    animate_scotese_elevation( )
end
# Makie map plots

function create_timeseries_charts(  )
    age = earliesttime 
    global world = read_world( age )
    n_time_points = Int( age / time_step) + 1
    n_diags = length(world_diag_names)
    n_frac_diags = length(world_frac_diag_names)
    diags_timeseries = fill(0.,n_diags,n_time_points)
    frac_diags_timeseries = fill(0.,n_frac_diags,n_sediment_types,n_time_points)
    cum_diags_timeseries = fill(0.,n_diags,n_time_points)
    cum_frac_diags_timeseries = fill(0.,n_frac_diags,n_sediment_types,n_time_points)
    land_inventory_timeseries = fill(0.,n_sediment_types,n_time_points)
    ocean_inventory_timeseries = fill(0.,n_sediment_types,n_time_points)
    area_timeseries = fill(0.,2,n_time_points)
    elevation_timeseries = fill(0.,n_time_points,3)
    time_points = fill(0.,n_time_points)
    i_time_point = 0
    last_step_inventories = 
        world_land_sediment_inventories( )[1:n_sediment_types] +
        world_ocean_sediment_inventories( )[1:n_sediment_types]
    fill( 0., n_sediment_types )
    step_changes = fill( 0., n_sediment_types, n_time_points )
    while age > animation_final_age
        age -= time_step
        println(age)
        i_time_point += 1
        time_points[i_time_point] = - age
        global world = read_world( age )
        area_timeseries[1,i_time_point] = land_area()
        area_timeseries[2,i_time_point] = ocean_area()
        land_inventory_timeseries[1:n_sediment_types,i_time_point] = 
            world_land_sediment_inventories( )[1:n_sediment_types]
        ocean_inventory_timeseries[1:n_sediment_types,i_time_point] = 
            world_ocean_sediment_inventories( )[1:n_sediment_types]
        Threads.@threads for i_diag in 1:n_diags
            diag_volume = volume_field( get_diag(world_diag_names[i_diag] ))
            diags_timeseries[i_diag,i_time_point] = diag_volume
            if i_time_point > 1
                cum_diags_timeseries[i_diag,i_time_point] = 
                    cum_diags_timeseries[i_diag,i_time_point-1] + diag_volume * time_step
            else
                cum_diags_timeseries[i_diag,i_time_point] = 
                    diag_volume * time_step
            end
        end
        Threads.@threads for i_frac_diag in 1:n_frac_diags
            for i_sedtype in 1:n_sediment_types
                frac_diag_volume = volume_field( get_frac_diag(world_frac_diag_names[i_frac_diag], i_sedtype ))
                frac_diags_timeseries[i_frac_diag,i_sedtype,i_time_point] = frac_diag_volume
                if i_time_point > 1
                    cum_frac_diags_timeseries[i_frac_diag,i_sedtype,i_time_point] = 
                        cum_frac_diags_timeseries[i_frac_diag,i_sedtype,i_time_point-1] + 
                        frac_diag_volume * time_step
                else
                    cum_frac_diags_timeseries[i_frac_diag,i_sedtype,i_time_point] = 
                        frac_diag_volume * time_step
                end
            end
        end

        step_changes[1:n_sediment_types,i_time_point] = (
            land_inventory_timeseries[1:n_sediment_types,i_time_point] .+
            ocean_inventory_timeseries[1:n_sediment_types,i_time_point] .-
            last_step_inventories ) ./ time_step
        last_step_inventories = 
            land_inventory_timeseries[1:n_sediment_types,i_time_point] .+
            ocean_inventory_timeseries[1:n_sediment_types,i_time_point]
        is_land = gt_mask(world.freeboard,0.)
        elevation_timeseries[i_time_point,1] = field_mean(world.freeboard .* is_land)
        elevation_timeseries[i_time_point,2] = field_mean(world.freeboard .* (1. .- is_land))
        elevation_timeseries[i_time_point,3] = field_max(world.freeboard)
    end
    n_time_points = length(time_points)
    cd( base_directory * "/" * output_directory * output_tag )
    plot_title = "clay fluxes"; variable_labels = ["prod" "subduct" "change" "bal"]
    units_label = "m3/Myr"
    n_lines = length(variable_labels)
    variable_point_array = fill(0.,n_time_points,n_lines) 
    clay_prod = diags_timeseries[
        diag_index("crust_clay_source_rate"),:]
    variable_point_array[:,1] = clay_prod[:]
    clay_subduct = frac_diags_timeseries[frac_diag_index(
        "ocean_subduct_sediment_fraction_thickness"),clay_sediment,:] +
        frac_diags_timeseries[frac_diag_index(
        "continent_subduct_sediment_fraction_thickness"),clay_sediment,:]
    variable_point_array[:,2] = clay_subduct
    variable_point_array[:,3] = step_changes[clay_sediment,:]
    variable_point_array[:,4] = step_changes[clay_sediment,:] + clay_subduct - clay_prod
    plot_time_series( time_points,variable_point_array,plot_title,variable_labels,units_label)

    plot_title = "CaCO3 fluxes"; variable_labels = ["prod" "subduct" "change" "bal"]
    units_label = "m3/Myr"
    n_lines = length(variable_labels)
    variable_point_array = fill(0.,n_time_points,n_lines) 
    CaCO3_prod = diags_timeseries[diag_index(
        "coastal_CaCO3_flux"),:] .+
        diags_timeseries[diag_index(
        "pelagic_CaCO3_deposition_rate"),:]
    variable_point_array[:,1] = CaCO3_prod
    CaCO3_subduct = frac_diags_timeseries[frac_diag_index(
        "ocean_subduct_sediment_fraction_thickness"),CaCO3_sediment,:] +
        frac_diags_timeseries[frac_diag_index(
        "continent_subduct_sediment_fraction_thickness"),CaCO3_sediment,:]
    variable_point_array[:,2] = CaCO3_subduct
    variable_point_array[:,3] = step_changes[CaCO3_sediment,:]
    variable_point_array[:,4] = step_changes[CaCO3_sediment,:] + CaCO3_subduct - CaCO3_prod
    plot_time_series( time_points,variable_point_array,plot_title,variable_labels,units_label)

    plot_title = "sediment inventories"
    variable_labels = ["Clay land" "Clay ocean" "Clay tot" "CaCO3 land" "CaCO3 ocean" "CaCO3 tot"  ]
    units_label = "m3"
    n_lines = length(variable_labels)
    variable_point_array = fill(0.,n_time_points,n_lines) 
    variable_point_array[:,1] = land_inventory_timeseries[clay_sediment,:]
    variable_point_array[:,2] = ocean_inventory_timeseries[clay_sediment,:]
    variable_point_array[:,3] = variable_point_array[:,1] + variable_point_array[:,2]
    variable_point_array[:,4] = land_inventory_timeseries[CaCO3_sediment,:]
    variable_point_array[:,5] = ocean_inventory_timeseries[CaCO3_sediment,:]
    variable_point_array[:,6] = variable_point_array[:,4] + variable_point_array[:,5]
    plot_time_series( time_points,variable_point_array,plot_title,variable_labels,units_label)

    plot_title = "cumulative clay fluxes"; variable_labels = ["inventory" "produced" "subducted" "bal"]
    units_label = "m3"
    n_lines = length(variable_labels)
    variable_point_array = fill(0.,n_time_points,n_lines) 
    clay_inv = land_inventory_timeseries[clay_sediment,:] .+ ocean_inventory_timeseries[clay_sediment,:]
    clay_prod = cum_diags_timeseries[diag_index(
        "crust_clay_source_rate"),:]
    clay_subduct = cum_frac_diags_timeseries[frac_diag_index(
        "ocean_subduct_sediment_fraction_thickness"),clay_sediment,:] +
        cum_frac_diags_timeseries[frac_diag_index(
        "continent_subduct_sediment_fraction_thickness"),clay_sediment,:]# .+ clay_inv[1]
    variable_point_array[:,1] = clay_inv
    variable_point_array[:,2] = clay_prod #.+ clay_inv[1]
    variable_point_array[:,3] = clay_subduct# .+ clay_inv[1]
    variable_point_array[:,4] = clay_inv .+ clay_subduct .- clay_prod 
    plot_time_series( time_points,variable_point_array,plot_title,variable_labels,units_label)

    plot_title = "cumulative CaCO3 fluxes"; variable_labels = ["inventory" "produced" "subducted" "bal"]
    units_label = "m3"
    n_lines = length(variable_labels)
    variable_point_array = fill(0.,n_time_points,n_lines) 
    CaCO3_inv = land_inventory_timeseries[CaCO3_sediment,:] .+ ocean_inventory_timeseries[CaCO3_sediment,:]
    CaCO3_prod = cum_diags_timeseries[diag_index(
        "coastal_CaCO3_flux"),:] .+
        cum_diags_timeseries[diag_index(
        "pelagic_CaCO3_deposition_rate"),:]
    CaCO3_subduct = cum_frac_diags_timeseries[frac_diag_index(
        "ocean_subduct_sediment_fraction_thickness"),CaCO3_sediment,:] +
        cum_frac_diags_timeseries[frac_diag_index(
        "continent_subduct_sediment_fraction_thickness"),CaCO3_sediment,:]# .+ clay_inv[1]
    variable_point_array[:,1] = CaCO3_inv
    variable_point_array[:,2] = CaCO3_prod #.+ clay_inv[1]
    variable_point_array[:,3] = CaCO3_subduct# .+ clay_inv[1]
    variable_point_array[:,4] = CaCO3_inv .+ CaCO3_subduct .- CaCO3_prod
    plot_time_series( time_points,variable_point_array,plot_title,variable_labels,units_label)

    plot_title = "average thickness"; variable_labels = ["land" "ocean"]
    units_label = "m"
    n_lines = length(variable_labels)
    variable_point_array = fill(0.,n_time_points,n_lines) 
    land_thickness = ( land_inventory_timeseries[clay_sediment,:] .+ 
        land_inventory_timeseries[CaCO3_sediment,:] ) ./ area_timeseries[1,:]
    ocean_thickness = ( ocean_inventory_timeseries[clay_sediment,:].+ 
        ocean_inventory_timeseries[CaCO3_sediment,:] ) ./ area_timeseries[2,:]
    variable_point_array[:,1] = land_thickness
    variable_point_array[:,2] = ocean_thickness
    plot_time_series( time_points,variable_point_array,plot_title,variable_labels,units_label)

    plot_title = "elevation"; variable_labels = ["land mean" "ocean mean" "max"]
    units_label = "m"
    plot_time_series( time_points,elevation_timeseries,plot_title,variable_labels,units_label)

    plot_title = "uplift rates"; variable_labels = ["continental" "subduction zone"]
    units_label = "m/Myr"
    variable_point_array[:,1] = diags_timeseries[
        diag_index("continent_orogenic_uplift_rate"),:]
    variable_point_array[:,2] = diags_timeseries[
        diag_index("subduction_orogenic_uplift_rate"),:]
    plot_time_series( time_points,variable_point_array,plot_title,variable_labels,units_label)
    
end
function plot_time_series(time_points,variable_point_array,plot_title,variable_labels,units_label)
    Plots.plot(time_points,variable_point_array,title=plot_title,label=variable_labels)
    Plots.xlabel!("Myr")
    Plots.ylabel!(units_label)
    cd( base_directory * "/" * output_directory * output_tag )
    charts_output_directory = base_directory * "/" * output_directory * output_tag * "/" *
        charts_directory
    #if charts_output_directory in readdir()
    #else
    #    mkdir("charts_output_directory")
    #end
    #rm(charts_output_directory,force=true,recursive=true)
    #mkdir(charts_output_directory)
    #cd(charts_output_directory)
    file_label = charts_output_directory * "/" *
        replace(plot_title," " => "_") * "." * output_tag * ".png"
    println("saving ", file_label)
    Plots.savefig(file_label)#=for i_line in 2:length(variable_labels) 
        println(i_line)
        Plots.plot!(time_points,variable_point_array[i_line])
    end
    scene = Figure(resolution=(1000, 800)) 
    ax = Axis(scene[1, 1], xlabel = "Mya", ylabel = "m^3 / Myr",
        title = plot_title)
    for i_line in 1:length(variable_labels) 
        lines!(ax,time_points, variable_point_array[i_line],label = plot_title[i_line])
    end
    axislegend(ax)
    return scene=#
end
function plot_add_time_series!(scene,time_points,variable)
    ax = Axis(scene[1, 1])
    lines!(ax,time_points,variable)
    return scene
end
function zoom_plot_field(field)
    minval,maxval = min_max_field(field)
    scene = zoom_plot_field(field,minval,maxval)
    return scene
end
function zoom_plot_field(field,minval,maxval)
    subx,suby,subfield = zoom_crop_field(field)
    scene = Scene(resolution = (1000, 800))
    cmap = Reverse(:lightrainbow)
    Makie.heatmap!(scene,subx,suby,subfield,
        colorrange=(maxval,minval),colormap=cmap)
    return scene
end
function zoom_works!(scene)
    scene = zoom_add_plate_boundaries!(scene)
    scene = zoom_add_continent_outlines!(scene)
    scene = zoom_add_orogenies!(scene)
    return scene
end
function zoom_crop_field(field)
    ixleft = max(1,zoomplot_ix-zoomplot_size)
    ixright = min(nx,zoomplot_ix+zoomplot_size)
    iybelow = max(1,zoomplot_iy-zoomplot_size)
    iyabove = min(ny,zoomplot_iy+zoomplot_size)
    subfield = field[ixleft:ixright,iybelow:iyabove]
    subx = xcoords[ixleft:ixright]
    suby = ycoords[iybelow:iyabove]
    return subx,suby,subfield
end
function zoom_add_plate_boundaries!(scene)
    for plateID in world.plateIDlist
        maskfield = eq_mask(world.plateID,plateID)
        subx,suby,subfield = zoom_crop_field(maskfield)
        Makie.contour!(scene,subx,suby,subfield,color=:white)
    end
    return scene
end
function zoom_add_coast_lines!(scene)
    continentmask = ge_mask(world.freeboard,0.)
    subx,suby,subfield = zoom_crop_field(continentmask)
    Makie.contour!(scene,subx,suby,subfield,color=:blue)
    return scene
end
function zoom_add_continent_outlines!(scene)
    continentmask = get_continent_mask()
    subx,suby,subfield = zoom_crop_field(continentmask)
    Makie.contour!(scene,subx,suby,subfield,color=:black)
    return scene
end
function zoom_add_orogenies!(scene)
    uplift_rate = get_diag("continent_orogenic_uplift_rate") .+
        get_diag("subduction_orogenic_uplift_rate")
    uplift_mask = gt_mask(uplift_rate,0.1)
    subx,suby,subfield = zoom_crop_field(uplift_mask)
    Makie.contour!(scene,subx,suby,subfield,color=:red)
    
    exposed_basement_field = eq_mask(world.geomorphology,exposed_basement)
    subx,suby,subfield = zoom_crop_field(exposed_basement_field)
    Makie.contour!(scene,subx,suby,subfield,color=:orange)

    #=ocean_shelf_field = eq_mask(world.geomorphology,ocean_shelf)
    subx,suby,subfield = zoomfield(ocean_shelf_field,ix,iy,nskirt)
    Makie.contour!(scene,subx,suby,subfield,color=:brown)=#
    #subduction_uplifting = gt_mask(get_diag("subduction_orogenic_uplift_rate"),0.01)
    #orogenies = active_orogenic_events()
    #text!(scene, orogenies, position = (0,-105),textsize=15)
    return scene
end

nxcolorbar = 5
function setup_plot_field()
    scene = Scene(resolution = (plot_resolution_x,plot_resolution_y)) #(1000, 550))  # GLMakie.lines([-180,180],[0,0],show_axis=false)
    setup_plot_field!(scene,0,0)
    return scene
end
function setup_double_plot_field()
    scene = Scene(resolution = (1000, 1000))  # GLMakie.lines([-180,180],[0,0],show_axis=false)
    setup_plot_field!(scene,0,0)
    setup_add_offset_field!(scene,0,-200)
    return scene
end
function setup_add_offset_field!(scene,xoffset,yoffset)
    setup_plot_field!(scene,xoffset,yoffset)
    return scene
end
function setup_plot_field!(scene,xoffset,yoffset)
    for i in [-180,-120,-60,0,60,120,180,180+nxcolorbar]
        x = [i,i] .+ xoffset; y = [-90,90] .+ yoffset
        lines!(scene,x,y)
    end
    for j in [-90,-60,-30,0,30,60,90]
        x = [-180,180+nxcolorbar] .+ xoffset; y = [j,j] .+ yoffset
        lines!(scene,x,y)
    end
    for i in [-200,220]
        x = [i,i] .+ xoffset; y = [-110,110] .+ yoffset
        lines!(scene,x,y,color=:white)
    end
    for j in [-110,110]
        x = [-200,200] .+ xoffset; y = [j,j] .+ yoffset
        lines!(scene,x,y,color=:white)
    end
    #    println("     range in plot ",min," ",max)
    return scene
end
function clear_plot_field!(scene)
    clear_plot_field!(scene,20)
    return scene
end
function clear_plot_field!(scene,ibase)
    for i in ibase+1:length(scene)
        delete!(scene,scene[end])
    end
    return scene
end
function add_colorbar_to_field(field,minval,maxval)
    xplotcoords = fill(0.,nx+nxcolorbar)
    for i in 1:length(xcoords)
        xplotcoords[i] = xcoords[i]
    end
    for i in nx+1:nx+nxcolorbar
        xplotcoords[i] = xplotcoords[i-1] + 1.
        #println(1191,i," ",xplotcoords[i])
    end
    colorfield = fill(0.,nx+nxcolorbar,ny)
    colorinterval = ( maxval - minval ) / ny
    #println(minval," ",maxval," ",colorinterval)
    for ix in 1:nx
        for iy in 1:ny
            colorfield[ix,iy] = field[ix,iy]
        end
    end
    for ix in nx+1:nx+nxcolorbar
        for iy in 1:ny
            colorfield[ix,iy] = minval + (iy-1) * colorinterval
        end
    end
    return colorfield,xplotcoords
end
function plot_field_with_colortab!(scene,field,minval,maxval,style::String="rainbow")
    plot_field_with_colortab!(scene,field,minval,maxval,0,0,style)
    return scene
end
function plot_field_with_colortab!(scene,field,minval,maxval,xoffset,yoffset,style::String="rainbow")
    #nxcolorbar = 5
    colorfield, xplotcoords = add_colorbar_to_field(field,minval,maxval)
    colorinterval = ( maxval - minval ) / 6
    for ilabel = 1:7
        labelval = Int(floor(minval + (ilabel-1) * colorinterval))
        labelval = string(labelval)
        labelpos = Int(-95. + (ilabel-1) * 30.)
        text!(scene,labelval,position=(190 + xoffset,labelpos + yoffset),textsize=10)
        #println("printing ", labelval)
    end
    #text!(scene,string(Int(floor(maxval))),position=(190,85),textsize=10)
    cmap = Reverse(:lightrainbow)
    if style == "freeboard"
        cmap = Reverse(:gist_earth)
    end
    Makie.heatmap!(scene,xplotcoords .+ xoffset,ycoords .+ yoffset,colorfield,
        colorrange=(maxval,minval),colormap=cmap)#thermometer)#lightrainbow)
    return scene
end
function plot_add_title!(scene,title)
    text!( scene, title, position = (-180, 110), textsize=25 ) 
    return scene
end
function plot_add_timestamp!(scene)
    scene = plot_add_timestamp!(scene,world.age,-180,-110)
    return scene
end
function plot_add_timestamp!(scene,age,xpos,ypos)
    timestamp = string(Int(floor(world.age))) * " Myr, " * get_geo_interval(age)
    text!(scene, timestamp, position = (xpos,ypos),textsize=15)
    return scene
end
function plot_add_plate_boundaries!(scene)
    for plateID in world.plateIDlist
        maskfield = eq_mask(world.plateID,plateID)
        Makie.contour!(scene,xcoords,ycoords,maskfield,color=:white)
    end
    return scene
end
function plot_add_alt_plate_boundaries!(scene)
    for plateID in world.plateIDlist
        maskfield = eq_mask(world.plateID,plateID)
        Makie.contour!(scene,xcoords,ycoords,maskfield,color=:red)
    end
    return scene
end
function plot_add_transport_regimes!(scene)
    maskfield = eq_mask(world.geomorphology,3)
    Makie.contour!(scene,xcoords,ycoords,maskfield,color=:red)
    maskfield = eq_mask(world.geomorphology,2)
    Makie.contour!(scene,xcoords,ycoords,maskfield,color=:green)
    maskfield = eq_mask(world.geomorphology,1)
    Makie.contour!(scene,xcoords,ycoords,maskfield,color=:blue)
    return scene
end

function generate_flowpath_line(lat, longt, plateID)
    # velocity arrow originating at lat longt in world grid
    vworldinitial = sphere2cart([lat,longt])  # x,y,z cartesian coords
    mplate = resolve_rotation_matrix!(plateID,world.age)
    vplate = mplate * vworldinitial # still xyz
    npoints = 2
    linex = [ convert(AbstractFloat,longt) ]
    liney = [ convert(AbstractFloat,lat) ]
    for i in 1:npoints-1  # 0
        ageval = world.age - i * 5. # 5 million years of travel
        mplate = resolve_rotation_matrix!(plateID,ageval)
        vworld = transpose(mplate) * vplate
        coordsworld = cart2sphere(vworld)
        latworldnow = coordsworld[1]; longtworldnow = coordsworld[2]
        push!(linex,longtworldnow); push!(liney,latworldnow)
    end
    return linex,liney
end
function plot_add_streamlines!(scene) # sprinkle flowpathlines over map
    lats = range(-80,stop=80,length=9)
    longts = range(-170,stop=170,length=18)
    for lat in lats
        for longt in longts
            iworld = search_sorted_nearest(xcoords,longt)
            jworld = search_sorted_nearest(ycoords,lat)
            originx = [longt-1,longt,longt+1,longt,longt-1]
            originy = [lat,lat+1,lat,lat-1,lat]
            lines!(scene,originx,originy,color=:black)
            plateID = world.plateID[iworld,jworld]
            linex,liney = generate_flowpath_line(lat, longt, plateID)
            lines!(scene,linex,liney,color=:black)
        end
    end
    return scene
end
function plot_add_coast_lines!(scene)
    landmask = ge_mask(world.freeboard,0.)
    Makie.contour!(scene,xcoords,ycoords,landmask,color=:black)
    return scene
end
function plot_add_continent_outlines!(scene)
    continentmask = get_continent_mask()
    Makie.contour!(scene,xcoords,ycoords,continentmask,color=:black)
    return scene
end
function plot_add_alt_continent_outlines!(scene)
    continentmask = get_continent_mask()
    Makie.contour!(scene,xcoords,ycoords,continentmask,color=:grey)
    return scene
end
function plot_add_sediment_thickness_contours!(scene,low,spacing,high)
    plot_add_seafloor_sediment_thickness_contours!(scene)
    plot_add_land_sediment_thickness_contours!(scene)
    return scene
end
function plot_add_seafloor_sediment_thickness_contours!(scene)
    is_ocean = 1. .- gt_mask(world.freeboard,0.)
    seafloor_sediment_thickness = world.sediment_thickness .* is_ocean
    Makie.contour!(scene,xcoords,ycoords,seafloor_sediment_thickness,levels=0:20:100,color=:yellow)
    Makie.contour!(scene,xcoords,ycoords,seafloor_sediment_thickness,levels=0:1000:5000,color=:white)
    return scene
end
function plot_add_land_sediment_thickness_contours!(scene)
    is_land = gt_mask(world.freeboard,0.)
    seafloor_sediment_thickness = world.sediment_thickness .* is_land
    Makie.contour!(scene,xcoords,ycoords,seafloor_sediment_thickness,levels=0:1000:5000,color=:green)
    Makie.contour!(scene,xcoords,ycoords,seafloor_sediment_thickness,levels=0:200:1000,color=:yellow)
    return scene
end
function plot_add_sediment_thickness_contours!(scene)
    plot_add_sediment_thickness_contours!(scene,0,1000,5000)
    plot_add_sediment_thickness_contours!(scene,0,5,25)
    return scene
end
function plot_add_subduction_orogeny!(scene)
    subduction_uplifting = gt_mask(get_diag("subduction_orogenic_uplift_rate"),0.01)
    Makie.contour!(scene,xcoords,ycoords,subduction_uplifting,color=:white)
    return scene
end
function plot_add_continental_orogeny!(scene)
    uplift_rate = get_diag("continent_orogenic_uplift_rate")
    Makie.contour!(scene,xcoords,ycoords,uplift_rate,color=:black)
    return scene
end
function plot_add_orogenies!(scene)
    uplift_rate = get_diag("continent_orogenic_uplift_rate")
    Makie.contour!(scene,xcoords,ycoords,uplift_rate,color=:red)

    subduction_uplifting = gt_mask(get_diag("subduction_orogenic_uplift_rate"),0.01)
    Makie.contour!(scene,xcoords,ycoords,subduction_uplifting,color=:green)

    exposed_basement_field = eq_mask(world.geomorphology,exposed_basement)
    Makie.contour!(scene,xcoords,ycoords,exposed_basement_field,color=:blue)

    orogenies = active_orogenic_events()
    text!(scene, orogenies, position = (0,-105),textsize=15)
    return scene
end
function plate_plot_add_outcrop!(scene,plate)
    world_outcrop = eq_mask(world.plateID, plate.plateID)
    plate_projection = projected_plate_maskfield!(world_outcrop,
        plate.rotationmatrix)
    Makie.contour!(scene,xcoords,ycoords,plate_projection,color=:white)
    return scene
end
function land_ocean_scale(field,land_low,land_high,ocean_low,ocean_high)
    newfield = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if world.crust_type[ix,iy] == continent_crust
                newfield[ix,iy] = ( field[ix,iy] - land_low ) / 
                    ( land_high - land_low )
            else
                newfield[ix,iy] = ( field[ix,iy] - ocean_low ) / 
                    ( ocean_high - ocean_low )
            end
        end
    end
    return newfield
end
function plot_sediment_thickness()
    scaled_thickness = land_ocean_scale(world.sediment_thickness,0.,2000.,0.,2000.)
    scene = plot_field(scaled_thickness,0.,1.)
    return scene
end
function plot_elevation()
    scaled_thickness = land_ocean_scale(world.freeboard,0.,5000.,-6000.,-2000.)
    scene = plot_field(scaled_thickness,0.,1.)
    return scene
end

function plot_field(field)  # for a generic field, auto scaling
    scene = setup_plot_field()
    minval,maxval = min_max_field(field)
    plot_field_with_colortab!(scene,field,minval,maxval)
    return scene
end
function plot_field(field,minval,maxval,style::String="rainbow")  # intended for generic field
    scene = setup_plot_field()
    plot_field_with_colortab!(scene,field,minval,maxval,style)
    return scene
end
function plot_two_fields(field1,field2,minval,maxval,style::String="rainbow")
    scene = setup_double_plot_field()
    #setup_add_offset_field!(scene,0,-200)
    plot_field_with_colortab!(scene,field1,minval,maxval,style)
    plot_field_with_colortab!(scene,field2,minval,maxval,0,-200,style)
    return scene
end
# depth surface plot
function setup_depth_surfaceplot()
    scene = Scene()
    for i in [-180,-120,-60,0,60,120,180,180+nxcolorbar]
        x = [i,i]; y = [-90,90]; z = [0,0]
        lines!(scene,x,y,z,show_axis=false)
    end
    for j in [-90,-60,-30,0,30,60,90]
        x = [-180,180+nxcolorbar]; y = [j,j];z = [0,0]
        lines!(scene,x,y,z)
    end
    eyeposition = [0,-0.1,.1]
    #cam = AbstractPlotting.cam3d_cad!(scene, eyeposition=eyeposition, fov=30f0)
    cam = cam3d_cad!(scene, eyeposition=eyeposition, fov=30f0)
    return scene
end
function color_depth_surfaceplot!(scene,depthfield,colorfield,minval,maxval)
    scaleddepth = depthfield * 0.003
    newcolorfield, xplotcoords = add_colorbar_to_field(colorfield,minval,maxval)
    newscaleddepth, xplotcoords = add_colorbar_to_field(scaleddepth,0.,0.)

    surface!(scene,xplotcoords,ycoords,newscaleddepth,color=newcolorfield,
    colormap=:lightrainbow,shading=false)
    #lightposition = Vec3f0(100, 0, -15), ambient = Vec3f0(1.,1.,1.)
    colorinterval = ( maxval - minval ) / 6
    for ilabel = 1:7
        labelval = Int(floor(minval + (ilabel-1) * colorinterval))
        labelval = string(labelval)
        labelpos = Int(-95. + (ilabel-1) * 30.)
        text!(scene,labelval,position=(190,labelpos),textsize=10)
        #println("printing ", labelval)
    end
    timestamp = string(Int(floor(age))) * " Myr"
    text!(scene, timestamp, position = (-180.,-110.),textsize=15)
    plot_fieldboundaries!(scene)
    plot_add_streamlines!(scene)
    plot_add_continent_outlines!(scene)
    return scene
end
function depthageplot()
    scene = setup_depth_surfaceplot()
    depthfield = world.freeboard
    #depthfield = fill(0.,nx,ny)
    minval,maxval = min_max_field(world.crust_age)
    if maxval > 200.
        maxval = 200.
    end
    if maxval < 10.
        maxval = 10.
    end
    minval = 0.
    color_depth_surfaceplot!(scene,world.freeboard,world.crust_age,
        minval,maxval)

    return scene
end
# spherical projection plots
function setup_sphereplot(lat,longt) # setup grid for surface plots
    scene = Scene()
    global sphereeyeposition = [lat,longt]
    u = range(-pi,stop=π,length=nx)
    v = range(0,stop=π,length=ny)
    global spheregridx = zeros(nx,ny)
    global spheregridy = zeros(nx,ny)
    global spheregridz = zeros(nx,ny)
    for i in 1:nx
        for j in 1:ny
            spheregridx[i,j] = cos.(u[i]) * sin(v[j]);
            spheregridy[i,j] = sin.(u[i]) * sin(v[j]);
            spheregridz[i,j] = cos(v[j]);
        end
    end
    #dummy_grid = fill(0.,nx,ny)
    # plot the poles
    pole_x = [0.e0,0.e0];pole_y = copy(pole_x); pole_z = [-1.5e0,0.e0]
    lines!(scene,pole_x,pole_y,pole_z,color=:green,show_axis=false) # South pole green
    pole_z = [0.e0,1.5e0];
    lines!(scene,pole_x,pole_y,pole_z,color=:red)  # North is red
    # longitude lines
    ugrid = range(-180,stop=180,length=13)
    vgrid = range(-90,stop=90,length=ny)
    for i in 1:12
        linex = [];liney = []
        for j in 1:ny
            append!(linex,ugrid[i]) # given arc all same longt
            append!(liney,vgrid[j])
        end
        arcx,arcy,arcz = line_sphere2cart(linex,liney)
        if i == 1
            lines!(scene,arcx,arcy,arcz,color=:red)  # +- 180
        elseif i == 7
            lines!(scene,arcx,arcy,arcz,color=:green)
        else
            lines!(scene,arcx,arcy,arcz)
        end
    end
    # latitude lines
    vgrid = range(-90,stop=90,length=9)
    for j in 2:8
        linex = []; liney = []
        for i in 1:nx
            push!(linex,u[i]/pi*180.)
            push!(liney,vgrid[j])
        end
        neweqx,neweqy,neweqz = line_sphere2cart(linex,liney)
        if j == 5
            lines!(scene,neweqx,neweqy,neweqz,color=:red)
        else
            lines!(scene,neweqx,neweqy,neweqz)
        end
    end
    adjust_eyeposition_sphereplot!(scene,lat,longt)
    return scene
end
function clear_sphereplot!(scene)
    for i in 22:length(scene)
        delete!(scene,scene[end])
    end
    return scene
end
function adjust_eyeposition_sphereplot!(scene,lat,longt)
    # set viewpoint
    global sphereeyeposition = [lat,longt]
    longtrad = longt / 180. * pi
    latrad = lat / 180. * pi
    scale = 1.5
    eyex = scale * cos(longtrad) * cos(latrad)
    eyey = scale * sin(longtrad) * cos(latrad)
    eyez = scale * sin(latrad)
    eyeposition = Float32[eyex,eyey,eyez]
    println("eyeposition now ", eyeposition)
    #    cam = AbstractPlotting.cam3d_cad!(scene, eyeposition=eyeposition)
    cam = cam3d_cad!(scene, eyeposition=eyeposition)
    #scene.camera.eyeposition = eyeposition
    #println(scene.camera.eyeposition)
    return scene
end
function move_eyeposition!(scene,increment) # [ dlat, dlongt ]
    adjust_eyeposition_sphereplot!(scene,sphereeyeposition[1]+increment[1],
    sphereeyeposition[2]+increment[2])
    return scene
end
function add_sphereplot_boundaries!(scene)
    plateIDfield = read_plateIDs(age)
    plateIDlist = find_plateID_list(plateIDfield)
    for plateID in world.plateIDlist
        maskfield = eq_mask(world.plateID,plateID)
        x2d,y2d = generate_outline_from_mask(maskfield)  # x y coords
        #x2d = x2d .+ 180.
        x3d,y3d,z3d = line_sphere2cart(x2d,y2d) # same as generate3doutline
        lines!(scene,x3d,y3d,z3d,color=:white)
    end
    return scene
end
function add_sphereplot_colors!(scene,field)
    colorfield = fill(0,nx,ny)
    for i in 1:nx
        for j in 1:ny
            colorfield[i,j] = field[i,ny-j+1]
        end
    end
    colorfield = convert(Array{Float64},colorfield)
    surface!(scene,spheregridx,spheregridy,spheregridz,
        color=colorfield,colormap=:lightrainbow )
    return scene
end
function add_sphereplot_streamlines!(scene,lat,longt)
    lats = range(lat-60,stop=lat+60,length=13)
    longts = range(longt-60,stop=longt+60,length=13)
    for lat in lats
        for longt in longts
            iworld = search_sorted_nearest(xcoords,longt)
            jworld = search_sorted_nearest(ycoords,lat)
            originx = [longt-1,longt,longt+1,longt]
            originy = [lat,lat+1,lat,lat-1]
            org3dx,org3dy,org3dz = line_sphere2cart(originx,originy)
            lines!(scene,org3dx,org3dy,org3dz,color=:white)
            plateID = world.plateID[iworld,jworld]
            x2d,y2d = generate_flowpath_line(lat, longt, plateID)
            x3d,y3d,z3d = line_sphere2cart(x2d,y2d) # same as generate3doutline
            lines!(scene,x3d,y3d,z3d,color=:white)
        end
    end
    return scene
end
function add_sphereplot_wig!(scene,field,maskfield,rotationmatrix)
    u = range(0,stop=2*π,length=nx)
    v = range(0,stop=π,length=ny)
    x = zeros(nx,ny); y = zeros(nx,ny); z = zeros(nx,ny)
    colorfield=fill(0.,nx,ny)
    for i in 1:nx
        for j in 1:ny
            if maskfield[i,ny-j+1] == 1
                elevation = 1.1
                colorfield[i,j] = field[i,ny-j+1]
            else
                elevation = 0.9
            end
            xplate = cos.(u[i]) * sin(v[j]) * elevation
            yplate = sin.(u[i]) * sin(v[j]) * elevation
            zplate = cos(v[j]) * elevation
            vpoint = [xplate,yplate,zplate]
            newvpoint = transpose(rotationmatrix) * vpoint
            # unrotating back to world grid
            x[i,j] = newvpoint[1]
            y[i,j] = newvpoint[2]
            z[i,j] = newvpoint[3]

        end
    end
    surface!(  x,y,z,color=colorfield,show_axis=false )
    return scene
end
function add_sphereplot_continent_outlines!(scene)
    continentmask = get_continent_mask()
    mycontour = Contour.contour(xcoords,ycoords,continentmask,0.5)
    lines = Contour.lines(mycontour)
    for line in lines
        xs,ys = Contour.coordinates(line)
        #xs = xs .+ 180.
        x3d,y3d,z3d = line_sphere2cart(xs,ys)
        #x3d = x3d .* 10.
        #x3d = [-1.,1.]; y3d = [-1.,1.]; z3d = [-1.,1.]
        lines!(scene,x3d,y3d,z3d,color=:black)
    end
    return scene
end
function sphere_boundary_plot(field,lat,longt)
    scene = setup_sphereplot(lat,longt)
    updatespherecolors!(scene,field)
    updatesphereboundaries!(scene)
    return scene
end
function sphere_stream_plot(field,lat,longt)
    scene = setupsphere_boundary_plot(lat,longt)
    add_sphereplot_streamlines!(scene,lat,longt)
    return scene
end
