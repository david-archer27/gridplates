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
function read_real_csv(filename)
    f = open(filename)
    arr = fill(0.,nx,ny)
    lines = readlines(f)
    for (iy,line) in enumerate(lines)
        iyflip = iy
        words = split(line,",")
        for (ix,word) in enumerate(words)
            value = parse(Float16,word)
            arr[ix,iyflip] = value
        end
    end
    return arr
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
    filename = plateID_input_directory * "/plateIDs." * string(integerage) * ".csv"
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
    filename = continent_input_directory * "/contIDs." * string(Int(integerage)) * ".csv"
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
    cd( output_location )
    output_directory_name = output_directory[findlast("/", output_directory )[1]+1:end]
    if output_directory_name in readdir()
    else
        mkdir( output_directory  )
        println( "creating ", output_directory  )
    end
    cd( output_directory )
    if "world" in readdir()
    else
        mkdir( world_directory )
        println( "creating ", world_directory )
    end
    if "plates" in readdir()
    else
        mkdir( plate_directory )
        println( "creating ", plate_directory )
    end
    if "code_bak" in readdir()
    else
        mkdir( code_backup_directory )
        println( "creating ", code_backup_directory )
        cd( code_base_directory )
        for file in readdir()
            if file[end-2:end] == ".jl"
               cp( file, code_backup_directory * "/" * file )
            end
        end
    end
    if "charts" in readdir()
    else
        mkdir( charts_directory )
        println( "creating ", charts_directory )
    end
    if "animations" in readdir()
    else
        mkdir( animation_directory )
        println( "creating ", animation_directory )
    end
    if enable_save_plate_transplant_images
        if "plate_transplants" in readdir()
        else
            mkdir( "plate_transplants" )
            println("creating plate_transplants")
        end
    end
    cd( code_base_directory )
end
function save_field(variable_name,field)
    timestamp = string(Int(ceil(world.age)))
    filename = output_directory * "/" * variable_name * "." * timestamp * ".bson"
    rm(filename, force=true)
    BSON.@save filename field
end
function save_world()
    timestamp = string(Int(ceil(world.age)))
    # in world grid, only the age field is non-trivial to calculate
    filename = world_directory * "/" * output_tag * ".world." * timestamp * ".bson"
    rm(filename, force=true)
    println("saving ", filename)
    BSON.@save filename world
    return
end
function read_world(age)
    #global world
    timestamp = string(Int(ceil(age)))
    # in world grid, only the age field is non-trivial to calculate
    filename = world_directory * "/" * output_tag * ".world." * timestamp * ".bson"
    println("reading world ",filename)
    BSON.@load filename world
    return world
end
function save_plates()
    timestamp = string(Int(ceil(world.age)))
    for plateID in world.plateIDlist
        plate = plates[plateID]
        IDstring = string(Int(plate.plateID))
        filename = plate_directory * "/" * output_tag * ".plate." * IDstring * "." * timestamp * ".bson"
        rm(filename, force=true)
        println("saving ", filename)
        BSON.@save filename plate
    end
    return
end
function read_plates( )
    age = world.age
    timestamp = string(Int(ceil(age)))
    for plateID in world.plateIDlist
        println("reading plate ", plateID)
        IDstring = string(Int(plateID))
        filename = plate_directory * "/" * output_tag * ".plate." * IDstring * "." * timestamp * ".bson"
        BSON.@load filename plate
        if haskey(plates,plateID) == false
            #println("creating plate ",plateID)
            plates[plateID] = create_blank_plate(plateID)
        end
        plates[plateID] = plate
    end
    return plates
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
function create_html_directory( )
    cd( output_location )
    html_directory_name = html_directory[findlast("/", output_directory )[1]+1:end]
    if html_directory_name in readdir()
    else
        mkdir( html_directory  )
        println( "creating ", html_directory  )
    end
    in_filename = utils_directory * "/index_template.html"
    out_filename = html_directory * "/index.html"
    #in_file = open(in_filename)
    out_file = open(out_filename,"w")
    println(out_file,"<HTML><HEAD><TITLE>GridPlates</TITLE></HEAD><body>")
    println(out_file,"<h1>Gridplates global sediment circulation model</h1>")
    println(out_file, "source code and description ")
    println(out_file,"<a href=\"https://github.com/david-archer27/gridplates/\">here</a></p>")
    println(out_file,"<table><tr>")
    println(out_file,"<td><h1>Animations</h1> (click to play)</td></tr><tr>")

    cd( data_directory )
    cp( "baum_compare_runoff.mp4", html_directory * "/baum_compare_runoff.mp4", force=true )
    cp( "baum_compare_runoff.png", html_directory * "/baum_compare_runoff.png", force=true )
    println(out_file,"<td><a href=\"baum_compare_runoff.mp4\">")
    println(out_file,"<img src=\"baum_compare_runoff.png\" width=\"500\"></a></p></td>")

    cd( animation_directory )
    for file in readdir()
        if file[end-3:end] == ".mp4"
            cp( file, html_directory * "/" * file, force=true )
            cp( file[1:end-4] * ".png", html_directory * "/" * file[1:end-4], force=true )
            
            println(out_file,"<td><a href=\"" * file * "\">" )
            println(out_file,"<img src=\"" * file[1:end-4] * "\".png\" width=\"500\">")
            println(out_file, "</a></p></td>" )  
        end
    end

    println(out_file, "</tr><tr><td><h1>Charts</h1></td></tr><tr>")

    cd( charts_directory )
    for file in readdir()
        if file[end-3:end] == ".png"
           cp( file, html_directory * "/" * file, force=true )
           println(out_file, "<td><img src=\"" * file * "\" width = \"500\"><p></td>")
        end
    end
    println(out_file, "</p></tr></table>")
    close(out_file)
    cd( code_base_directory )
end
#function plot_elevation_frame()
function animate( plot_function, plot_type )
    cd( animation_directory )
    if plot_type in readdir()
    else
        mkdir( plot_type )
        println( "creating ", plot_type  )
    end
    cd( plot_type )
    starting_file_list = readdir()
    image_number = 0
    ages = [ animation_final_age ]
    for age in animation_initial_age: - main_time_step * animation_n_step : animation_final_age
        push!(ages,age)
    end
    for age in ages
        image_number += 1
        image_file_name = "img." * lpad(image_number,3,"0") * ".png"
        if image_file_name in starting_file_list
            println( "already done ", age," ", plot_type * "/" * image_file_name )
        else
            println( "creating ", age, " ", plot_type * "/" * image_file_name )
            global world = read_world( age )
            scene = plot_function()
            Makie.save( image_file_name, scene )
        end
    end
    mp4_file = "../" * plot_type * "." * output_tag * ".mp4"
    println("compiling ", mp4_file)
    rm( mp4_file,force=true)
    run(`ffmpeg -r 2 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p $mp4_file`)
    cp( "img.001.png", "../" * plot_type * "." * output_tag * ".png", force=true )
    for imgfile in readdir()
        rm( imgfile )
    end
    cd( animation_directory )
    rm( plot_type, force=true )
    cd( code_base_directory )
end
function plot_elevation()
    #scaled_thickness = land_ocean_scale(world.freeboard,0.,5000.,-6000.,0.)
    scene = setup_plot_field()
    scene = plot_field_with_colortab!(scene,world.freeboard,-5000,5000,0,0,"freeboard")
    plot_add_coast_lines!(scene)
    plot_add_sediment_thickness_contours!(scene)
    plot_add_timestamp!(scene,world.age,-180,-105) 
    plot_add_title!(scene,"Elevation")
    #plot_field(scaled_thickness,0.,1.)
    return scene
end
function plot_pct_CaCO3()
    scene = plot_field(world.sediment_surface_fractions[:,:,CaCO3_sediment] .* 
        gt_mask(world.sediment_thickness,0.) .* 100.,0.,100.)
    plot_add_plate_boundaries!(scene)
    plot_add_orogenies!(scene)
    plot_add_continent_outlines!(scene)
    plot_add_timestamp!(scene,world.age,-180,-105)    #scene = plot_field(world.freeboard,-5000.,5000)
    plot_add_title!(scene,"%CaCO3")
    cont_depo = volume_field(
        get_diag("continental_CaCO3_deposition_rate")) * # m3 / myr
        rho_sediment / # g / yr 
        100. / 1.e12   #  E12 mol / yr
    cont_diss = volume_field(
        get_frac_diag("land_sediment_fraction_dissolution_rate",
        CaCO3_sediment)) * # m3 / myr
        rho_sediment / # g / yr 
        100. / 1.e12   #  E12 mol / yr
    coast_depo = volume_field(
        get_diag("coastal_CaCO3_flux")) * # m3 / myr
        rho_sediment / # g / yr 
        100. / 1.e12  #  E12 mol / yr
    pelagic_depo = volume_field(
        get_diag("pelagic_CaCO3_deposition_rate")) * # m3 / myr
        rho_sediment / # g / yr 
        100. / 1.e12 #  E12 mol / yr
    net = cont_depo-cont_diss+coast_depo+pelagic_depo
    text_out = "CaCO3 deposition: Continental " * my_string(cont_depo) * " - " * my_string(cont_diss) *
        ", Coastal " * my_string(coast_depo) * ", Pelagic " * my_string(pelagic_depo) *
        ", Net " * my_string(net )
    text!( scene, text_out, position = (-60, 110), textsize=20 ) 
    continental_CaCO3_dissolution_rate = volume_field(
        get_frac_diag("land_sediment_fraction_dissolution_rate",CaCO3_sediment )) * # m3 / myr
        rho_sediment / # g / yr 
        100. / 1.e12 #  E12 mol / yr
    continental_CaO_dissolution_rate = volume_field(
        get_frac_diag("land_sediment_fraction_dissolution_rate",CaO_sediment )) * # m3 / myr
        rho_sediment / # g / yr 
        100. / 1.e12 #  E12 mol / yr
    land_orogenic_Ca_source_rate = volume_field( 
        get_diag("land_orogenic_Ca_source_rates") ) * # m3 / myr
        rho_sediment / # g / yr 
        100. / 1.e12 #  E12 mol / yr
    total_Ca_sources = land_orogenic_Ca_source_rate + 
        #continental_CaCO3_dissolution_rate + 
        continental_CaO_dissolution_rate
    text_out = "Ca sources: hard rocks " * my_string(land_orogenic_Ca_source_rate) *
        ", clays " * my_string(continental_CaO_dissolution_rate) * 
        #", CaCO3 " * my_string(continental_CaCO3_dissolution_rate) * 
        ", total " * my_string(total_Ca_sources)
    text!( scene, text_out, position = (-60, 100), textsize=20 ) 
    return scene
end
function my_string( value )
    string_output = lpad( value, 4, " " )[1:4]
    return string_output
end 
function plot_pct_CaO_CaCO3_free()
    CaO_CaCO3_free = fill(0.,nx,ny)
    #CaO_calcite_free = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if world.sediment_thickness[ix,iy] > 0.
                if world.sediment_surface_fractions[ix,iy,CaO_sediment] + 
                    world.sediment_surface_fractions[ix,iy,clay_sediment] > 0.
                    CaO_CaCO3_free[ix,iy] = world.sediment_surface_fractions[ix,iy,CaO_sediment]  / 
                        ( world.sediment_surface_fractions[ix,iy,CaO_sediment] 
                        + world.sediment_surface_fractions[ix,iy,clay_sediment] ) 
                end
            end
        end
    end    
    high_value = orogenic_sediment_source_fractions[CaO_sediment] * 100.
    scene = plot_field(CaO_CaCO3_free.* 100,0.,high_value)
    plot_add_plate_boundaries!(scene)
    plot_add_orogenies!(scene)
    plot_add_continent_outlines!(scene)
    plot_add_coast_lines!(scene)
    plot_add_timestamp!(scene,world.age,-180,-105)    #scene = plot_field(world.freeboard,-5000.,5000)
    plot_add_title!(scene,"% CaO (CaCO3-free)")
    continental_CaO_dissolution_rate = volume_field(
        get_frac_diag("land_sediment_fraction_dissolution_rate",CaO_sediment )) * # m3 / myr
        rho_sediment / # g / yr 
        100. / 1.e12 #  E12 mol / yr
    land_orogenic_Ca_source_rate = volume_field( 
        get_diag("land_orogenic_Ca_source_rates") ) * # m3 / myr
        rho_sediment / # g / yr 
        100. / 1.e12 #  E12 mol / yr
    land_orogenic_CaO_source_rate = volume_field(
        get_frac_diag(("land_orogenic_fraction_flux"),CaO_sediment)) * # m3 / myr
        rho_sediment / # g / yr 
        100. / 1.e12 #  E12 mol / yr
    total_orogenic_CaO_source_rate = volume_field(
        get_frac_diag(("crust_orogenic_fraction_flux"),CaO_sediment)) * # m3 / myr
        rho_sediment / # g / yr 
        100. / 1.e12 #  E12 mol / yr
        total_Ca_sources = land_orogenic_Ca_source_rate + 
        #continental_CaCO3_dissolution_rate + 
        continental_CaO_dissolution_rate
    text_out = "Ca sources: hard rocks " * my_string(land_orogenic_Ca_source_rate) *
        ", clays " * my_string(continental_CaO_dissolution_rate) * 
        #", CaCO3 " * my_string(continental_CaCO3_dissolution_rate) * 
        ", total " * my_string(total_Ca_sources)
    text!( scene, text_out, position = (-60, 110), textsize=20 ) 
    tot_Ca_diss = land_orogenic_Ca_source_rate + continental_CaO_dissolution_rate
    tot_Ca_diss_max = land_orogenic_CaO_source_rate + land_orogenic_Ca_source_rate
    text_out = "Land CaO weathering efficiency " * my_string( tot_Ca_diss / tot_Ca_diss_max )
    text!( scene, text_out, position = (-60, 100), textsize=20 ) 
    return scene
end
function plot_weathering_index()
    land_sediment_weathering_index = fill(0.,nx,ny)
    CaO_CaCO3_free = fill(0.,nx,ny)
    #CaO_calcite_free = fill(0.,nx,ny)
    for ix in 1:nx
        for iy in 1:ny
            if world.sediment_thickness[ix,iy] > 0.
                if world.sediment_surface_fractions[ix,iy,CaO_sediment] + 
                    world.sediment_surface_fractions[ix,iy,clay_sediment] > 0.
                    CaO_CaCO3_free[ix,iy] = world.sediment_surface_fractions[ix,iy,CaO_sediment]  / 
                        ( world.sediment_surface_fractions[ix,iy,CaO_sediment] 
                        + world.sediment_surface_fractions[ix,iy,clay_sediment] ) 
                    land_sediment_weathering_index[ix,iy] = 1. - 
                        CaO_CaCO3_free[ix,iy] / orogenic_sediment_source_fractions[CaO_sediment] 
                end
            end
        end
    end
    scene = plot_field(land_sediment_weathering_index .* 100.,0.,100.)
    plot_add_plate_boundaries!(scene)
    plot_add_orogenies!(scene)
    plot_add_continent_outlines!(scene)
    plot_add_timestamp!(scene,world.age,-180,-105)    #scene = plot_field(world.freeboard,-5000.,5000)
    plot_add_title!(scene,"CaO depletion relative to unreactive clay")
    continental_CaO_dissolution_rate = volume_field(
        get_frac_diag("land_sediment_fraction_dissolution_rate",CaO_sediment )) * # m3 / myr
        rho_sediment / # g / yr 
        100. / 1.e12 #  E12 mol / yr
    land_orogenic_Ca_source_rate = volume_field( 
        get_diag("land_orogenic_Ca_source_rates") ) * # m3 / myr
        rho_sediment / # g / yr 
        100. / 1.e12 #  E12 mol / yr
    land_orogenic_CaO_source_rate = volume_field(
        get_frac_diag(("land_orogenic_fraction_flux"),CaO_sediment)) * # m3 / myr
        rho_sediment / # g / yr 
        100. / 1.e12 #  E12 mol / yr
    total_orogenic_CaO_source_rate = volume_field(
        get_frac_diag(("crust_orogenic_fraction_flux"),CaO_sediment)) * # m3 / myr
        rho_sediment / # g / yr 
        100. / 1.e12 #  E12 mol / yr
        total_Ca_sources = land_orogenic_Ca_source_rate + 
        #continental_CaCO3_dissolution_rate + 
        continental_CaO_dissolution_rate
    text_out = "Ca sources: hard rocks " * my_string(land_orogenic_Ca_source_rate) *
        ", clays " * my_string(continental_CaO_dissolution_rate) * 
        #", CaCO3 " * my_string(continental_CaCO3_dissolution_rate) * 
        ", total " * my_string(total_Ca_sources)
    text!( scene, text_out, position = (0, 110), textsize=20 ) 
    tot_Ca_diss = land_orogenic_Ca_source_rate + continental_CaO_dissolution_rate
    tot_Ca_diss_max = land_orogenic_CaO_source_rate + land_orogenic_Ca_source_rate
    text_out = "Land CaO weathering efficiency " * my_string( tot_Ca_diss / tot_Ca_diss_max )
    text!( scene, text_out, position = (0, 100), textsize=20 ) 
    return scene
end
function plot_CaO_weathering_rate()
    weathering_rate = get_frac_diag("land_sediment_fraction_dissolution_rate",CaO_sediment) +
        get_diag("land_orogenic_Ca_source_rates")
    scene = plot_field(weathering_rate,0.,3.)
    plot_add_plate_boundaries!(scene)
    plot_add_orogenies!(scene)
    plot_add_continent_outlines!(scene)
    plot_add_coast_lines!(scene)
    plot_add_land_sediment_thickness_contours!(scene)
    plot_add_timestamp!(scene,world.age,-180,-105)    #scene = plot_field(world.freeboard,-5000.,5000)
    plot_add_title!(scene,"CaO Weathering Flux, m/Myr")
    return scene
end
function plot_sediment_thickness()
    scene = plot_field(world.sediment_thickness ./ 1000.,0.,1.) # plot_sediment_thickness()
    plot_add_plate_boundaries!(scene)
    plot_add_orogenies!(scene)
    plot_add_continent_outlines!(scene)
    plot_add_coast_lines!(scene)
    plot_add_sediment_thickness_contours!(scene)
    plot_add_timestamp!(scene,world.age,-180,-105)
    plot_add_title!(scene,"Sediment Thickness") # L(0:2km)/O(0:2km)")
    return scene
end
function plot_sedimentation_rate()
    sed_rate = get_diag("global_sediment_deposition_rate")
    scene = plot_field(sed_rate,0.,100.) # plot_sediment_thickness()
    plot_add_plate_boundaries!(scene)
    plot_add_orogenies!(scene)
    plot_add_continent_outlines!(scene)
    plot_add_coast_lines!(scene)
    plot_add_sediment_thickness_contours!(scene)
    plot_add_timestamp!(scene,world.age,-180,-105)
    plot_add_title!(scene,"Sedimentation Rate, m/Myr") # L(0:2km)/O(0:2km)")
    return scene
end
function plot_crust_age( )
    scene = plot_field(world.crust_age,0.,600.)
    plot_add_plate_boundaries!(scene)
    plot_add_orogenies!(scene)
    plot_add_continent_outlines!(scene)
    plot_add_timestamp!(scene,world.age,-180,-105)
    plot_add_title!(scene,"Crust age, Myr")
    return scene
end
function plot_crust_thickness( )
    scene = plot_field(world.crust_thickness ./ 1000,0.,50.)
    plot_add_plate_boundaries!(scene)
    plot_add_orogenies!(scene)
    plot_add_continent_outlines!(scene)
    plot_add_timestamp!(scene,world.age,-180,-105)
    plot_add_title!(scene,"Crust thickness, km")
    return scene
end
function plot_crust_thickening_rate( )
    oro = get_diag("continent_orogenic_uplift_rate")
    scene = plot_field(oro,0.,500.)
    plot_add_plate_boundaries!(scene)
    plot_add_orogenies!(scene)
    plot_add_continent_outlines!(scene)
    plot_add_timestamp!(scene,world.age,-180,-105)
    plot_add_title!(scene,"Crust thickening rate, m / Myr")
    return scene
end
function plot_weathering_rate( )
    oro = get_diag("land_orogenic_Ca_source_rates")
    sed = get_frac_diag("land_sediment_fraction_dissolution_rate",CaO_sediment)
    weat = oro + sed
    scene = plot_field(weat)
    plot_add_plate_boundaries!(scene)
    plot_add_orogenies!(scene)
    plot_add_continent_outlines!(scene)
    plot_add_coast_lines!(scene)
    plot_add_timestamp!(scene,world.age,-180,-105)
    plot_add_title!(scene,"CaO weathering rate, m / Myr")
    return scene
end
function plot_runoff( )
    runoff = get_diag("land_Q_runoff_field")
    scene = plot_field(runoff,0.,40.)
    plot_add_plate_boundaries!(scene)
    plot_add_orogenies!(scene)
    plot_add_continent_outlines!(scene)
    plot_add_coast_lines!(scene)
    plot_add_land_sediment_thickness_contours!(scene)
    plot_add_timestamp!(scene,world.age,-180,-105)
    plot_add_title!(scene,"Runoff")
    return scene
end
function plot_sediment_thickness_age(image_file_name::String="")
    bin_width = 5.; bin_max = 200
    n_bins = Int( bin_max / bin_width )
    ages = fill(0.,n_bins)
    variable_labels = ["gridplates" "Olson 2016 3rd-order global present-day"]
    plot_title = string(Int(world.age)) * " Myr"
    n_variables = length(variable_labels)
    variable_point_array = fill(0.,n_bins,n_variables)
    area_totals = fill(0.,n_bins)
    thickness_totals = fill(0.,n_bins)
    # observational fit from Olson 
    c0 = 7.389E1
    c1 = 9.443
    c2 = -1.358e-1
    c3 = 1.396e-3
    #olson_fit = fill(0.,n_bins)
    for i_column in 1:n_bins
        age = bin_width * ( i_column - 1 )
        ages[i_column] = age
        if age < 120
            variable_point_array[i_column,2] = c0 + c1 * age + c2 * age^2 + c3 * age^3
        else
            variable_point_array[i_column,2] = NaN
        end
    end
    for ix in 1:nx
        for iy in 1:ny
            if world.crust_type[ix,iy] == ocean_crust
                age_bin = Int(floor( world.crust_age[ix,iy] / bin_width )) + 1
                if age_bin < n_bins
                    thickness_totals[age_bin] += world.sediment_thickness[ix,iy] * areabox[iy]
                    area_totals[age_bin] += areabox[iy]
                end
            end
        end
    end
    variable_point_array[:,1] = thickness_totals ./ ( area_totals .+ 1. )
    Plots.plot(ages,variable_point_array,title=plot_title,label=variable_labels)
    Plots.xlabel!("Crust age, Myr")
    Plots.ylabel!("Mean Sediment Thickness, m")
    if image_file_name == ""
        return
    else
        if image_file_name == "chart"
            image_file_name = charts_directory * "/" *
                "sed_thickness_age." * output_tag * ".png"
        else
            anim_directory = animation_directory * "/" * "sed_thickness_age/" 
            image_file_name = anim_directory * image_file_name
        end
        println("saving ", image_file_name)
        Plots.savefig(image_file_name)
    end
end
function animate_sed_thickness_age()
    cd( animation_directory )
    plot_type = "sed_thickness_age" 
    if plot_type in readdir()
    else
        mkdir( plot_type )
        println( "creating ", plot_type  )
    end
    cd( plot_type )
    starting_file_list = readdir()
    image_number = 0
    ages = [ animation_final_age ]
    for age in animation_initial_age: - main_time_step * animation_n_step : animation_final_age
        push!(ages,age)
    end
    for age in ages
        image_number += 1
        image_file_name = "img." * lpad(image_number,3,"0") * ".png"
        if image_file_name in starting_file_list
            println( "already done ", age," ", image_file_name )
        else
            println( "creating ", age, " ", image_file_name )
            global world = read_world( age )
            plot_sediment_thickness_age( image_file_name )
        end
    end
    mp4_file = "../sed_thickness_age." * output_tag * ".mp4"
    println("compiling ", mp4_file)
    rm( mp4_file,force=true)
    run(`ffmpeg -r 2 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p $mp4_file`)
    cp( "img.001.png", "../" * plot_type * "." * output_tag * ".png", force=true )
    for imgfile in readdir()
        rm( imgfile )
    end
    cd( animation_directory )
    rm( plot_type, force=true )
    cd( code_base_directory )
end
function animate_scotese_elevation( )
    cd( animation_directory )
    plot_type = "scotese_elevation" 
    if plot_type in readdir()
    else
        mkdir( plot_type )
        println( "creating ", plot_type  )
    end
    cd( plot_type )
    starting_file_list = readdir()
    image_number = 0
    ages = [ animation_final_age ]
    for age in animation_initial_age: - main_time_step * animation_n_step : animation_final_age
        push!(ages,age)
    end
    for age in ages
        image_number += 1
        image_file_name = "img." * lpad(image_number,3,"0") * ".png"
        if image_file_name in starting_file_list
            println( "already done ", age," ", image_file_name )
        else
            println( "creating ", age, " ", image_file_name )
            global world = read_world( age )
            scotese_elevation,scotese_age = nearest_scotese_elevation()
            scene = plot_two_fields(world.freeboard,scotese_elevation,-4000,4000,"freeboard")
            gridplatestimestamp = string(Int(floor(world.age))) * 
                " Myr, " * get_geo_interval(age) * 
                ", " * string(world.sealevel) * " m sea level"
            scotesetimestamp = string(Int(floor(scotese_age))) * " Myr"
            text!(scene,gridplatestimestamp, position = (0,95),textsize=15)
            text!(scene, scotesetimestamp, position = (0,-110),textsize=15)
            #plot_add_plate_boundaries!(scene)
            plot_add_elevation_contours!(scene)
            smooth_scotese = scotese_elevation[1:nx,1:ny]
            smooth_world!(smooth_scotese,1000.)
            Makie.contour!(scene,xcoords,ycoords .- 200.,smooth_scotese,levels=0:1000:4000,color=:blue)
            Makie.contour!(scene,xcoords,ycoords .- 200.,smooth_scotese,levels=-5000:1000:-3000,color=:white)
        
            plot_add_orogenies!(scene)
            #plot_add_continent_outlines!(scene)
            plot_add_coast_lines!(scene)
            text!(scene,"GridPlates Elevation",position=(-180,95),textsize=20)
            text!(scene,"Scotese Reconstruction",position=(-180,-110),textsize=20)
            Makie.save( image_file_name, scene )
        end
    end
    mp4_file = "../scotese_elevation." * output_tag * ".mp4"
    println("compiling ", mp4_file)
    rm( mp4_file,force=true)
    run(`ffmpeg -r 2 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p $mp4_file`)
    cp( "img.001.png", "../" * plot_type * "." * output_tag * ".png", force=true )
    for imgfile in readdir()
        rm( imgfile )
    end
    cd( animation_directory )
    rm( plot_type, force=true )
    cd( code_base_directory )
end

#function animate_first_priority()

function animate_until_you_almost_puke( )
    #animate( plot_elevation, "elevation" )
    animate( plot_pct_CaCO3, "pct_CaCO3" )
    animate( plot_pct_CaO_CaCO3_free, "pct_CaO_CaCO3_free" )
    animate( plot_weathering_index, "weathering_index" )
    animate( plot_CaO_weathering_rate, "CaO_weathering_rate" )
    animate( plot_sediment_thickness, "sediment_thickness" )
    animate( plot_sedimentation_rate, "sedimentation_rate" )
end
function animate_the_rest()
    #animate( plot_crust_age, "crust_age" )
    #animate( plot_crust_thickness, "crust_thickness" )
    animate( plot_crust_thickening_rate, "crust_thickening_rate" )
    animate( plot_runoff, "runoff")

    animate_sed_thickness_age()
    animate_scotese_elevation( )
end
# Makie map plots

function create_timeseries_charts(  )
    age = earliesttime 
    global world = read_world( age )
    n_time_points = Int( age / main_time_step)
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
    subduction_rates = fill(0.,0:n_sediment_types,n_time_points)
    cum_subduction_rates = fill(0.,0:n_sediment_types,0:n_time_points)
    cum_subduction_rates[:,0] .= 0.
    time_points = fill(0.,n_time_points)
    i_time_point = 0
    last_step_inventories = 
        world_land_sediment_inventories( )[1:n_sediment_types] +
        world_ocean_sediment_inventories( )[1:n_sediment_types]
    fill( 0., n_sediment_types )
    step_changes = fill( 0., n_sediment_types, n_time_points )
    fraction_denuded = fill(0.,n_time_points)
    mean_weathering_index = fill(0.,n_time_points)
    while age > animation_final_age
        age -= main_time_step
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
        subduction_rates[1:n_sediment_types,i_time_point] .= 
            world.subducted_ocean_sediment_volumes .+
            world.subducted_land_sediment_volumes
        for i_sedtype in 1:n_sediment_types
            subduction_rates[0,i_time_point] += subduction_rates[i_sedtype,i_time_point]
        end
        #Threads.@threads 
        for i_diag in 1:n_diags
            diag_volume = volume_field( get_diag(world_diag_names[i_diag] ))
            diags_timeseries[i_diag,i_time_point] = diag_volume
            #=if i_time_point > 1
                cum_diags_timeseries[i_diag,i_time_point] = 
                    cum_diags_timeseries[i_diag,i_time_point-1] + diag_volume * main_time_step
            else
                cum_diags_timeseries[i_diag,i_time_point] = 
                    diag_volume * main_time_step
            end=#
        end
        #Threads.@threads 
        for i_frac_diag in 1:n_frac_diags
            for i_sedtype in 1:n_sediment_types
                frac_diag_volume = volume_field( get_frac_diag(world_frac_diag_names[i_frac_diag], i_sedtype ))
                frac_diags_timeseries[i_frac_diag,i_sedtype,i_time_point] = frac_diag_volume
                #=if i_time_point > 1
                    cum_frac_diags_timeseries[i_frac_diag,i_sedtype,i_time_point] = 
                        cum_frac_diags_timeseries[i_frac_diag,i_sedtype,i_time_point-1] + 
                        frac_diag_volume * main_time_step
                else
                    cum_frac_diags_timeseries[i_frac_diag,i_sedtype,i_time_point] = 
                        frac_diag_volume * main_time_step
                end
                cum_subduction_rates[:,i_time_point] .= 
                    cum_subduction_rates[:,i_time_point-1] +
                    subduction_rates[:,i_time_point]=#
            end
        end

        step_changes[1:n_sediment_types,i_time_point] = (
            land_inventory_timeseries[1:n_sediment_types,i_time_point] .+
            ocean_inventory_timeseries[1:n_sediment_types,i_time_point] .-
            last_step_inventories ) ./ main_time_step
        last_step_inventories = 
            land_inventory_timeseries[1:n_sediment_types,i_time_point] .+
            ocean_inventory_timeseries[1:n_sediment_types,i_time_point]
        #is_land = gt_mask(world.freeboard,0.)
        elevation_timeseries[i_time_point,1] = field_mean(world.freeboard, is_land())
        elevation_timeseries[i_time_point,2] = field_mean(world.freeboard, (1 .- is_land()))
        elevation_timeseries[i_time_point,3] = field_max(world.freeboard)
        area_denuded = volume_field(eq_mask(world.geomorphology,exposed_basement))
        area_covered = volume_field(eq_mask(world.geomorphology,sedimented_land) )
        fraction_denuded[i_time_point] = area_denuded / ( area_denuded + area_covered )
        weathering_index = get_diag("land_sediment_weathering_index")
        mean_weathering_index[i_time_point] = field_mean(weathering_index, is_land())
    end
    n_time_points = length(time_points)

    if enable_geomorph

        plot_title = "clay fluxes"; variable_labels = ["prod" "land deposition" "runoff" "ocean deposition" "subduction" "change" "bal"]
        linestyles = [ :solid :solid :solid :dash :dash :dash :solid ]
        linewidths = [ :auto :auto :auto :auto :auto :auto 3 ]
        units_label = "m3/Myr"
        n_lines = length(variable_labels)
        variable_point_array = fill(0.,n_time_points,n_lines) 
        variable_point_array[:,1] = frac_diags_timeseries[
            frac_diag_index("crust_orogenic_fraction_flux"),clay_sediment,:]
        variable_point_array[:,2] = frac_diags_timeseries[
            frac_diag_index("land_sediment_fraction_deposition_rate"),clay_sediment,:]
        variable_point_array[:,3] = frac_diags_timeseries[
            frac_diag_index("coastal_sediment_fraction_runoff_flux"),clay_sediment,:]
        variable_point_array[:,4] = frac_diags_timeseries[
            frac_diag_index("seafloor_sediment_fraction_deposition_rate"),clay_sediment,:]
        variable_point_array[:,5] = subduction_rates[clay_sediment,:]
        variable_point_array[:,6] = step_changes[clay_sediment,:]
        variable_point_array[:,7] = variable_point_array[:,6] + # step_changes
            variable_point_array[:,5] - # subduct
            variable_point_array[:,1] # prod
        plot_time_series( time_points,variable_point_array,plot_title,
            variable_labels,
            linestyles,linewidths,units_label )# ,

        summed_variable_point_array = accumulate_values!( variable_point_array )
        plot_title = "cumulative clay fluxes"; variable_labels = ["inventory" "produced" "subducted"]
        linestyles = [ :solid :solid :solid :solid ]
        linewidths = [ :auto :auto :auto 3 ]
        units_label = "m3"
        n_lines = length(variable_labels)
        variable_point_array = fill(0.,n_time_points,n_lines) 
        variable_point_array[:,1] = land_inventory_timeseries[clay_sediment,:]
        variable_point_array[:,2] = summed_variable_point_array[:,1]
        variable_point_array[:,3] = summed_variable_point_array[:,5]
        plot_time_series( time_points,variable_point_array,plot_title,
            variable_labels,
            linestyles,linewidths,units_label) #,
    
        plot_title = "CaCO3 fluxes"; variable_labels = ["Continental" "coastal" "pelagic" "subduct"]
        linestyles = [ :solid :solid :solid :dash :dash :dash :dash :solid ]
        linewidths = [ :auto :auto :auto :auto :auto :auto :auto 3 ]
        units_label = "m3/Myr"
        n_lines = length(variable_labels)
        variable_point_array = fill(0.,n_time_points,n_lines) 
        #=CaCO3_prod = diags_timeseries[diag_index(
            "coastal_CaCO3_flux"),:] .+
            diags_timeseries[diag_index(
            "pelagic_CaCO3_deposition_rate"),:]=#
        variable_point_array[:,1] = diags_timeseries[diag_index(
                "continental_CaCO3_deposition_rate"),:] .-
            frac_diags_timeseries[frac_diag_index(
                "land_sediment_fraction_dissolution_rate"),CaCO3_sediment,:]
        variable_point_array[:,2] = diags_timeseries[diag_index(
            "coastal_CaCO3_flux"),:]
        variable_point_array[:,3] = diags_timeseries[diag_index(
            "pelagic_CaCO3_deposition_rate"),:]
        #variable_point_array[:,4] = frac_diags_timeseries[
        #    frac_diag_index("land_sediment_fraction_dissolution_rate"),CaCO3_sediment,:]
        #variable_point_array[:,4] = frac_diags_timeseries[
        #    frac_diag_index("coastal_sediment_fraction_runoff_flux"),CaCO3_sediment,:]
        variable_point_array[:,4] = subduction_rates[CaCO3_sediment,:] .* -1.
        #=variable_point_array[:,7] = step_changes[CaCO3_sediment,:]
        variable_point_array[:,8] = variable_point_array[:,7] + # step_changes
            variable_point_array[:,6] - # subduction_rates
            variable_point_array[:,1] - # continental_CaCO3_deposition_rate
            variable_point_array[:,2] - # coastal_CaCO3_flux
            variable_point_array[:,3] + # pelagic_CaCO3_deposition_rate
            variable_point_array[:,4] # land_sediment_fraction_dissolution_rate=#
        plot_time_series( time_points,variable_point_array,plot_title,
            variable_labels,
            linestyles,linewidths,units_label ) #,

        summed_variable_point_array = accumulate_values!( variable_point_array )
        plot_title = "cumulative CaCO3 fluxes"; variable_labels = ["inventory" "produced" "subducted"]
        linestyles = [ :solid :solid :solid :solid ]
        linewidths = [ :auto :auto :auto 3 ]
        linestyles = [ :solid :solid :solid :solid ]
        linewidths = [ :auto :auto :auto 3 ]
        units_label = "m3"
        n_lines = length(variable_labels)
        variable_point_array = fill(0.,n_time_points,n_lines) 
        variable_point_array[:,1] = land_inventory_timeseries[CaCO3_sediment,:]
        variable_point_array[:,2] = summed_variable_point_array[:,1] + # net continent
        summed_variable_point_array[:,2] + # coastal
        summed_variable_point_array[:,3] # pelagic
        variable_point_array[:,3] = summed_variable_point_array[:,4]
        plot_time_series( time_points,variable_point_array,plot_title,
            variable_labels,
            linestyles,linewidths,units_label )#,

        plot_title = "CaO fluxes"; variable_labels = ["orogenic diss" "orogenic prod" "clay diss" "runoff" "subduct" ]#"change" "bal"]
        linestyles = [ :solid :solid :dash :dash :dash :solid ]
        linewidths = [ :auto :auto :auto :auto :auto 3 ]
        units_label = "m3/Myr"
        n_lines = length(variable_labels)
        variable_point_array = fill(0.,n_time_points,n_lines) 
        variable_point_array[:,1] = diags_timeseries[diag_index(
            "land_orogenic_Ca_source_rates"),:] # direct orogenic Ca dissolution
        variable_point_array[:,2] = frac_diags_timeseries[
            frac_diag_index("crust_orogenic_fraction_flux"),CaO_sediment,:]
        variable_point_array[:,3] = frac_diags_timeseries[
            frac_diag_index("land_sediment_fraction_dissolution_rate"),CaO_sediment,:]
        variable_point_array[:,4] = frac_diags_timeseries[
            frac_diag_index("coastal_sediment_fraction_runoff_flux"),CaO_sediment,:]
        variable_point_array[:,5] = -1. .* subduction_rates[CaO_sediment,:]
        #=variable_point_array[:,5] = step_changes[CaO_sediment,:]
        variable_point_array[:,6] = variable_point_array[:,5] + # step_changes
            variable_point_array[:,4] - # subduct
            variable_point_array[:,1] + # crust_orogenic_fraction_flux
            variable_point_array[:,2] # land_sediment_fraction_dissolution_rate =#
        plot_time_series( time_points,variable_point_array,plot_title,
            variable_labels,
            linestyles,linewidths,units_label ) #,

        summed_variable_point_array = accumulate_values!( variable_point_array )
        plot_title = "cumulative CaO fluxes"; variable_labels = ["inventory" "produced" "subducted"]
        linestyles = [ :solid :solid :solid :solid ]
        linewidths = [ :auto :auto :auto 3 ]
        linestyles = [ :solid :solid :solid :solid ]
        linewidths = [ :auto :auto :auto 3 ]
        units_label = "m3"
        n_lines = length(variable_labels)
        variable_point_array = fill(0.,n_time_points,n_lines) 
        variable_point_array[:,1] = land_inventory_timeseries[CaO_sediment,:]
        variable_point_array[:,2] = summed_variable_point_array[:,2]
        variable_point_array[:,3] = summed_variable_point_array[:,5]
        plot_time_series( time_points,variable_point_array,plot_title,
            variable_labels,
            linestyles,linewidths,units_label)
    
        plot_title = "Ca fluxes"
        # ground truth from BLAG
        # 63% from CaCO3, 7% sed CaO, 14% oro CaO, (16% dolomite)
        variable_labels = ["Orogenic src" "Sed CaO diss src" "Cont CaCO3 net" "Pelagic sink" "Coastal dep sink" "net"]
        linestyles = [ :solid :solid :solid :dash :dash :dash :solid ]
        linewidths = [ :auto :auto :auto :auto :auto :auto 3 ]
        units_label = "m3 CaCO3 eq / Myr"
        n_lines = length(variable_labels)
        variable_point_array = fill(0.,n_time_points,n_lines) 
        variable_point_array[:,1] = diags_timeseries[diag_index(
            "land_orogenic_Ca_source_rates"),:]      
        variable_point_array[:,2] = frac_diags_timeseries[
            frac_diag_index("land_sediment_fraction_dissolution_rate"),CaO_sediment,:]
        variable_point_array[:,3] = frac_diags_timeseries[
            frac_diag_index("land_sediment_fraction_dissolution_rate"),CaCO3_sediment,:] -
            diags_timeseries[diag_index("continental_CaCO3_deposition_rate"),:]  
        variable_point_array[:,4] = - diags_timeseries[diag_index(
                "pelagic_CaCO3_deposition_rate"),:]      
        variable_point_array[:,5] = - diags_timeseries[diag_index(
                "coastal_CaCO3_flux"),:]   
        for i_var in 1:5   
            variable_point_array[:,6] += variable_point_array[:,i_var]
        end
        plot_time_series( time_points,variable_point_array,plot_title,
            variable_labels,
            linestyles,linewidths,units_label)

        plot_title = "CaO weathering efficiency"
        variable_labels = "diss / prod" # [ "oro / tot src"  ]
        linestyles = :solid 
        linewidths = :auto
        units_label = "%"
        variable_point_array = fill(0.,n_time_points)   
        tot_CaO_diss = diags_timeseries[diag_index(
                "land_orogenic_Ca_source_rates"),:] .+ 
            frac_diags_timeseries[
                frac_diag_index("land_sediment_fraction_dissolution_rate"),CaO_sediment,:]
        tot_CaO_prod = tot_CaO_diss .+
            frac_diags_timeseries[
                frac_diag_index("crust_orogenic_fraction_flux"),CaO_sediment,:]
        variable_point_array = tot_CaO_diss ./ tot_CaO_prod .* 100.
        plot_time_series( time_points,variable_point_array,plot_title,
            variable_labels,
            linestyles,linewidths,units_label )

        plot_title = "CaO weathering sources"
        variable_labels = "oro / tot src"
        linestyles = :solid 
        linewidths = :auto
        units_label = "%"
        variable_point_array = fill(0.,n_time_points)   
        #oro_CaO_prod = frac_diags_timeseries[
        #    frac_diag_index("crust_orogenic_fraction_flux"),CaO_sediment,:]
        oro_Ca_src = diags_timeseries[diag_index(
            "land_orogenic_Ca_source_rates"),:]
        variable_point_array = tot_CaO_diss ./ tot_CaO_prod .* 100.
        scatter_time_points = [ 0. ]
        scatter_point_array = [ 66. ] # BLAG eq 11,12 or table 1: igneous=14,sed Ca-sil=7
        scatter_labels = ["BLAG table 1"]
        plot_time_series_present_day_compare( 
            time_points,variable_point_array,plot_title,
            variable_labels,
            linestyles,linewidths,units_label,
            scatter_time_points,scatter_point_array,scatter_labels )
    
        plot_title = "weathering index"
        variable_labels = "mean weathering index"
        linestyles = :auto
        linewidths = :auto
        units_label = "%"
        plot_time_series( time_points,mean_weathering_index , 
            plot_title,
            variable_labels,
            linestyles,linewidths,units_label)#,

        plot_title = "CO2 fluxes"
        variable_labels = ["Orogenic weathering" "Sed CaO weathering" "CaCO3 subduction" "deux ex machina degas"]
        linestyles = [ :dash :dash :solid :solid :dash :dash :solid ]
        linewidths = [ :auto :auto :auto :auto :auto :auto 3 ]
        units_label = "m3 CaCO3 eq / Myr"
        n_lines = length(variable_labels)
        variable_point_array = fill(0.,n_time_points,n_lines) 
        variable_point_array[:,1] = - diags_timeseries[diag_index(
            "land_orogenic_Ca_source_rates"),:]      
        variable_point_array[:,2] = - frac_diags_timeseries[
            frac_diag_index("land_sediment_fraction_dissolution_rate"),CaO_sediment,:]
        variable_point_array[:,3] = subduction_rates[CaCO3_sediment,:]
        
        #=frac_diags_timeseries[
            frac_diag_index("land_sediment_fraction_dissolution_rate"),CaCO3_sediment,:] -
            diags_timeseries[diag_index("continental_CaCO3_deposition_rate"),:]  
        variable_point_array[:,4] = - diags_timeseries[diag_index(
                "pelagic_CaCO3_deposition_rate"),:]      
        variable_point_array[:,5] = - diags_timeseries[diag_index(
                "coastal_CaCO3_flux"),:]   =#
        for i_var in 1:2   
            variable_point_array[:,4] -= variable_point_array[:,i_var]
        end
        plot_time_series( time_points,variable_point_array,plot_title,
            variable_labels,
            linestyles,linewidths,units_label)

        plot_title = "sediment inventories"
        variable_labels = ["Clay land" "Clay ocean" "Clay tot" "CaCO3 land" "CaCO3 ocean" "CaCO3 tot" "CaO land" "CaO ocean" "CaO tot" ]
        linestyles = [ :solid :solid :solid :dash :dash :dash :solid :solid :solid ]
        linewidths = [ :auto :auto 3 :auto :auto 3 :auto :auto 3 ]
        units_label = "m3"
        n_lines = length(variable_labels)
        variable_point_array = fill(0.,n_time_points,n_lines) 
        variable_point_array[:,1] = land_inventory_timeseries[clay_sediment,:]
        variable_point_array[:,2] = ocean_inventory_timeseries[clay_sediment,:]
        variable_point_array[:,3] = variable_point_array[:,1] + variable_point_array[:,2]
        variable_point_array[:,4] = land_inventory_timeseries[CaCO3_sediment,:]
        variable_point_array[:,5] = ocean_inventory_timeseries[CaCO3_sediment,:]
        variable_point_array[:,6] = variable_point_array[:,4] + variable_point_array[:,5]
        variable_point_array[:,7] = land_inventory_timeseries[CaO_sediment,:]
        variable_point_array[:,8] = ocean_inventory_timeseries[CaO_sediment,:]
        variable_point_array[:,9] = variable_point_array[:,7] + variable_point_array[:,8]
        plot_time_series( time_points,variable_point_array,plot_title,
            variable_labels,
            linestyles,linewidths,units_label )

        plot_title = "fraction denuded"
        linestyles = :solid
        linewidths = :auto
        units_label = "%"; variable_labels = "gridplates"
        variable_point_array = fill(0.,n_time_points) 
        variable_point_array = fraction_denuded .* 100.
        scatter_time_points =  [0.] 
        scatter_point_array =  [25.] # , 35.
        scatter_labels = "Igneous and metamorphic (Holland 1978)" #  "Shield, acid volcanic, and basalts (Suchet 2003)"]
        plot_time_series_present_day_compare( 
            time_points,variable_point_array, 
            plot_title,
            variable_labels,
            linestyles,linewidths,units_label,
            scatter_time_points,scatter_point_array,scatter_labels )

            #n_lines = 1
        
        plot_compare_sedthick()

    end

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
    plot_time_series( time_points,variable_point_array,plot_title,
        variable_labels,
        linestyles,linewidths,units_label)

    plot_title = "elevation"; variable_labels = ["gridplates mean"] # * 10" "ocean mean" "max"]
    units_label = "m"
    #elevation_timeseries[:,1] .*= 10
    scotese_mean_elevations = scotese_mean_elevation_timeseries()
    linestyles = [ :solid :solid :solid :dash :dash :dash :solid ]
    linewidths = [ :auto :auto :auto :auto :auto :auto 3 ]
    plot_time_series_data_series_compare( time_points,elevation_timeseries[:,1],plot_title,
        variable_labels,
        linestyles,linewidths,units_label,
        scotese_elevation_plot_ages .* -1, scotese_mean_elevations, "Scotese mean" )

    plot_title = "uplift rates"; variable_labels = ["continental" "subduction zone"]
    units_label = "m3/Myr"
    variable_point_array[:,1] = diags_timeseries[
        diag_index("continent_orogenic_uplift_rate"),:]
    variable_point_array[:,2] = diags_timeseries[
        diag_index("subduction_orogenic_uplift_rate"),:]
    plot_time_series( time_points,variable_point_array,plot_title,
        variable_labels,
        linestyles,linewidths,units_label )
    
    plot_compare_elevations()

end
function accumulate_values!( variable_point_array )
    n_values = size( variable_point_array )[1]
    n_lines = size( variable_point_array )[2]
    for i_value in 2:n_values
        for i_line in 1:n_lines
            variable_point_array[i_value,i_line,] += 
                variable_point_array[i_value-1,i_line]
        end
    end
    return variable_point_array
end

function plot_time_series(time_points,variable_point_array,
    plot_title,variable_labels,
    linestyles,linewidths,units_label)

    Plots.plot(time_points,variable_point_array,
        title=plot_title,label=variable_labels,
        linestyle=linestyles,linewidth=linewidths)

    Plots.xlabel!("Myr")
    Plots.ylabel!(units_label)
    file_label = charts_directory * "/" *
        replace(plot_title," " => "_") * "." * output_tag * ".png"
    println("saving ", file_label)
    Plots.savefig(file_label)
end
function plot_time_series_present_day_compare(
    time_points,variable_point_array,
    plot_title,variable_labels,
    linestyles,linewidths,units_label,
    scatter_time_points,scatter_point_array, scatter_labels)

    #println(scatter_point_array)

    Plots.plot(time_points,variable_point_array,title=plot_title,
        label=variable_labels,linestyle=linestyles,
        linewidth=linewidths)
    Plots.scatter!(scatter_time_points,scatter_point_array,
        label=scatter_labels)

    Plots.xlabel!("Myr")
    Plots.ylabel!(units_label)
    file_label = charts_directory * "/" *
        replace(plot_title," " => "_") * "." * output_tag * ".png"
    println("saving ", file_label)
    Plots.savefig(file_label)
end
function plot_time_series_data_series_compare(
    time_points,variable_point_array,
    plot_title,variable_labels,
    linestyles,linewidths,units_label,
    scatter_time_points,scatter_point_array, scatter_labels)

    #println(scatter_point_array)

    Plots.plot(time_points,variable_point_array,title=plot_title,
        label=variable_labels,linestyle=linestyles,
        linewidth=linewidths)
    Plots.plot!(scatter_time_points,scatter_point_array,
        label=scatter_labels)

    Plots.xlabel!("Myr")
    Plots.ylabel!(units_label)
    file_label = charts_directory * "/" *
        replace(plot_title," " => "_") * "." * output_tag * ".png"
    println("saving ", file_label)
    Plots.savefig(file_label)
end
function plot_add_time_series!(scene,time_points,variable)
    ax = Axis(scene[1, 1])
    lines!(ax,time_points,variable)
    return scene
end
function plot_compare_elevations()
    etopo = read_flip_csv( 
        code_base_directory * "/data/" * "etopo.csv")
    world = read_world(0)
    obs_elevation_grid = [-4000.]
    while obs_elevation_grid[end] < 4000.
        push!(obs_elevation_grid,obs_elevation_grid[end] + 100. )
    end
    n_points = length(obs_elevation_grid)
    variable_point_list = fill(0.,n_points,2)
    sum_elevation_x_area = fill(0.,n_points)
    sum_area = fill(0.,n_points)
    for ix in 1:nx
        for iy in 1:ny
            i_bin = search_sorted_nearest(obs_elevation_grid,etopo[ix,iy])
            sum_elevation_x_area[i_bin] += world.freeboard[ix,iy] * areabox[iy]
            sum_area[i_bin] += areabox[iy]
        end
    end
    variable_point_list[:,1] = sum_elevation_x_area ./ ( sum_area .+ 1. )
    variable_point_list[:,2] = obs_elevation_grid
    plot_title = "freeboard comparison with ETOPO"
    variable_names = ["gridplates present-day" "1:1"]
    Plots.plot(obs_elevation_grid,variable_point_list,title=plot_title,label=variable_names)
    file_label = charts_directory * "/" * 
        "freeboard_compare." * output_tag * ".png"
    println("saving ", file_label)
    Plots.savefig(file_label)
end
function plot_compare_sedthick()
    sedthick = read_real_csv(code_base_directory * "/data/" * "sedthick.csv")
    
    #ratio = fill(0.,nx,ny)
    totals = [0.,0.]; total_areas = [0.,0.]
    for ix in 1:nx
        for iy in 1:ny
            if sedthick[ix,iy] > 0.
                #ratio[ix,iy] = world.sediment_thickness[ix,iy] / sedthick[ix,iy]
                totals[1] += sedthick[ix,iy] * areabox[iy]
                total_areas[1] += areabox[iy]
            end
            if world.sediment_thickness[ix,iy] > 0.
                totals[2] += world.sediment_thickness[ix,iy] * areabox[iy]
                total_areas[2] += areabox[iy]
            end
        end
    end
    totals ./= total_areas
    #scene = plot_field(ratio,0.,2.)
    seafloor_sediment_thickness = world.sediment_thickness .* lt_mask( world.freeboard,0. )
    scene = plot_two_fields(seafloor_sediment_thickness,sedthick,0,1000)
    Makie.contour!(scene,xcoords,ycoords,seafloor_sediment_thickness,levels=0:1000:10000,color=:yellow)
    Makie.contour!(scene,xcoords,ycoords,seafloor_sediment_thickness,levels=0:100:1000,color=:white)

    smooth_world!( sedthick, 1000. )
 
    Makie.contour!(scene,xcoords,ycoords .- 200.,sedthick,levels=0:1000:10000,color=:yellow)
    Makie.contour!(scene,xcoords,ycoords .- 200.,sedthick,levels=0:100:1000,color=:white)

    total_ratio = totals[2] / totals[1]
    output_string = "Gridplates ocean sediment vol / data = " * string(total_ratio)[1:4]
    text!(scene, output_string, position = (-140,105),textsize=25)
    file_label = charts_directory * "/" * 
        "sedthick_compare." * output_tag * ".png"
    println("saving ", file_label)
    Makie.save(file_label,scene)
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
        text!(scene,labelval,position=(190 + xoffset,labelpos + yoffset),textsize=20)
        #println("printing ", labelval)
    end
    #text!(scene,string(Int(floor(maxval))),position=(190,85),textsize=10)
    cmap = Reverse(:lightrainbow)
    if style == "freeboard"
        deepocean = colorant"midnightblue"
        shallowocean = colorant"turquoise1"
        ocean_cmap = range(deepocean,stop=shallowocean,length=100)
        lowland = colorant"darkgreen"
        upland = colorant"moccasin"
        land_cmap = range(lowland,stop=upland,length=100)
        cmap = vcat(ocean_cmap,land_cmap)
        cmap = Reverse(cmap)
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
    timestamp = string(Int(floor(world.age))) * " Myr, " * get_geo_interval(age) *
        ", Sea level = " * string(Int(floor(world.sealevel))) * " m"
    text!(scene, timestamp, position = (xpos,ypos),textsize=15)
    return scene
end
function plot_add_plate_boundaries!(scene)
    for plateID in world.plateIDlist
        maskfield = eq_mask(world.plateID,plateID)
        Makie.contour!(scene,xcoords,ycoords,maskfield,color=:black)
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
    Makie.contour!(scene,xcoords,ycoords,landmask,color=:blue)
    return scene
end
function plot_add_continent_outlines!(scene)
    continentmask = get_continent_mask()
    Makie.contour!(scene,xcoords,ycoords,continentmask,color=:grey)
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
    Makie.contour!(scene,xcoords,ycoords,seafloor_sediment_thickness,levels=0:100:1000,color=:yellow)
    Makie.contour!(scene,xcoords,ycoords,seafloor_sediment_thickness,levels=0:20:100,color=:white)
    return scene
end
function plot_add_land_sediment_thickness_contours!(scene)
    is_land = gt_mask(world.freeboard,0.)
    seafloor_sediment_thickness = world.sediment_thickness .* is_land
    #Makie.contour!(scene,xcoords,ycoords,seafloor_sediment_thickness,levels=0:1000:5000,color=:green)
    Makie.contour!(scene,xcoords,ycoords,seafloor_sediment_thickness,levels=0:500:2000,color=:yellow)
    return scene
end
function plot_add_sediment_thickness_contours!(scene)
    plot_add_sediment_thickness_contours!(scene,0,1000,5000)
    plot_add_sediment_thickness_contours!(scene,0,5,25)
    return scene
end
function plot_add_elevation_contours!(scene)
    #is_land = gt_mask(world.freeboard,0.)
    Makie.contour!(scene,xcoords,ycoords,world.freeboard,levels=0:500:4000,color=:blue)
    Makie.contour!(scene,xcoords,ycoords,world.freeboard,levels=-6000:500:-1000,color=:white)
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

    #subduction_uplifting = gt_mask(get_diag("subduction_orogenic_uplift_rate"),0.01)
    #Makie.contour!(scene,xcoords,ycoords,subduction_uplifting,color=:green)

    exposed_basement_field = eq_mask(world.geomorphology,exposed_basement)
    Makie.contour!(scene,xcoords,ycoords,exposed_basement_field,color=:grey)

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
    scene = Scene(resolution = (plot_resolution_x,plot_resolution_y)) #(1000, 550))  # GLMakie.lines([-180,180],[0,0],show_axis=false)
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
    scaleddepth = depthfield * 0.005
    newcolorfield, xplotcoords = add_colorbar_to_field(colorfield,minval,maxval)
    newscaleddepth, xplotcoords = add_colorbar_to_field(scaleddepth,0.,0.)

    Makie.surface!(scene,xplotcoords,ycoords,newscaleddepth,color=newcolorfield,
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
    #timestamp = string(Int(floor(age))) * " Myr"
    #text!(scene, timestamp, position = (-180.,-110.),textsize=15)
    #plot_fieldboundaries!(scene)
    #plot_add_streamlines!(scene)
    plot_add_coast_lines!(scene)
    return scene
end
function depth_age_plot()
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
function depth_CaCO3_plot()
    scene = setup_depth_surfaceplot()
    #depthfield = world.freeboard
    #depthfield = fill(0.,nx,ny)

    color_depth_surfaceplot!(scene,world.freeboard,
        world.sediment_surface_fractions[:,:,CaCO3_sediment] .* 100.,
        0.,100.)

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

 