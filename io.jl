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
    filename = "platefiles/plateIDs." * string(integerage) * ".csv"
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
    filename = "contfiles/contIDs." * string(Int(integerage)) * ".csv"
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
function save_field(variable_name,field)
    timestamp = string(Int(ceil(world.age)))
    filename = "../outfiles/fields/" * variable_name * "." * timestamp * ".bson"
    rm(filename, force=true)
    BSON.@save filename field
end

function save_world()
    runname = "SM2"
    timestamp = string(Int(ceil(world.age)))
    # in world grid, only the age field is non-trivial to calculate
    filename = "../outfiles/world/" * runname * ".world." * timestamp * ".bson"
    rm(filename, force=true)
    println("saving ", filename)
    BSON.@save filename world
    return
end
function read_world(age)
    runname = "SM2"
    timestamp = string(Int(ceil(age)))
    # in world grid, only the age field is non-trivial to calculate
    filename = data_output_directory * "/world/" * runname * ".world." * timestamp * ".bson"
    BSON.@load filename world
    return world
end
function save_plates()
    runname = "SM2"
    timestamp = string(Int(ceil(world.age)))
    # in world grid, only the age field is non-trivial to calculate
    filename = "../outfiles/plates/" * runname * ".plates." * timestamp * ".bson"
    rm(filename, force=true)
    println("saving ", filename)
    BSON.@save filename plates
    return
end
function save_plates_checkpoint()
    runname = "SM2"
    timestamp = string(Int(ceil(world.age)))
    # in world grid, only the age field is non-trivial to calculate
    filename = "../outfiles/plates/" * runname * ".plates.checkpoint.bson"
    rm(filename, force=true)
    println("saving ", filename)
    BSON.@save filename plates
    return
end
function read_plates(age)
    runname = "SM2"
    timestamp = string(Int(ceil(age)))
    filename = "../outfiles/plates/" * runname * ".plates." * timestamp * ".bson"
    BSON.@load filename plates
    return plates
end
function logging_println(text,array)
    println(text,array)
    if log_IO != 0
        println(log_IO,text,array)
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

# Makie map plots
function animate_all( age )
    plot_types = ["elevation","pct_CaCO3","crust_age","sed_thickness","scotese_elevation"]
    image_number = 0
    global world = read_world( age )
    for plot_type in plot_types
        directory = animation_output_directory * plot_type
        cd( directory )
        file_list = readdir()
        for file in file_list
            if occursin(".png",file)
                run(`rm $file`)
            end
        end
    end
    while world.age > 0
        age -= time_step
        println("reading ", age)
        global world = read_world( age )
        image_number += 1
        # elevation comparison with scotese reconstructions
        scotese_elevation,scotese_age = nearest_scotese_elevation()
        scene = plot_two_fields(world.freeboard,scotese_elevation,-4000,4000)
        scotesetimestamp = string(Int(floor(scotese_age))) * " Myr"
        text!(scene, scotesetimestamp, position = (-180,-305),textsize=15)
        plot_add_timestamp!(scene,scotese_age,-180,-305)
        outfilename = animation_output_directory * "scotese_elevation/img." * lpad(image_number,3,"0") * ".png"
        println(outfilename)
        Makie.save(outfilename,scene)
        # sediment thickness
        scene = plot_sediment_thickness()
        plot_add_plate_boundaries!(scene)
        plot_add_orogenies!(scene)
        plot_add_continent_outlines!(scene)
        plot_add_timestamp!(scene,world.age,-180,-105)
        outfilename = animation_output_directory * "sed_thickness/img." * lpad(image_number,3,"0") * ".png"
        println(outfilename)
        Makie.save(outfilename,scene)
        # crust age
        scene = plot_field(world.crust_age,0.,600.)
        plot_add_plate_boundaries!(scene)
        plot_add_orogenies!(scene)
        plot_add_continent_outlines!(scene)
        plot_add_timestamp!(scene,world.age,-180,-105)
        outfilename = animation_output_directory * "crust_age/img." * lpad(image_number,3,"0") * ".png"
        println(outfilename)
        Makie.save(outfilename,scene)
        # percent CaCO3
        scene = plot_field(world.sediment_surface_fractions[:,:,2] .* 100.,0.,100.)
        plot_add_plate_boundaries!(scene)
        plot_add_orogenies!(scene)
        plot_add_continent_outlines!(scene)
        plot_add_timestamp!(scene,world.age,-180,-105)    #scene = plot_field(world.freeboard,-5000.,5000)
        outfilename = animation_output_directory * "pct_CaCO3/img." * lpad(image_number,3,"0") * ".png"
        println(outfilename)
        Makie.save(outfilename,scene)
        # detailed elevation
        scene = plot_elevation()
        plot_add_plate_boundaries!(scene)
        plot_add_orogenies!(scene)
        plot_add_continent_outlines!(scene)
        plot_add_timestamp!(scene,world.age,-180,-105)    #scene = plot_field(world.freeboard,-5000.,5000)
        outfilename = animation_output_directory * "elevation/img." * lpad(image_number,3,"0") * ".png"
        println(outfilename)
        Makie.save(outfilename,scene)        
    end
    for plot_type in plot_types
        println("compiling ",plot_type)
        directory = animation_output_directory * plot_type
        cd( directory )
        rm(plot_type * ".mp4",force=true)
        run(`ffmpeg -r 10 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p $plot_type.mp4`)
        #=file_list = readdir()
        for file in file_list
            if occursin(".png",file)
                run(`rm $file`)
            end
        end=#
    end
    cd( code_base_directory ) 
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
    heatmap!(scene,subx,suby,subfield,
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
        contour!(scene,subx,suby,subfield,color=:white)
    end
    return scene
end
function zoom_add_continent_outlines!(scene)
    continentmask = get_continent_mask()
    subx,suby,subfield = zoom_crop_field(continentmask)
    contour!(scene,subx,suby,subfield,color=:black)
    return scene
end
function zoom_add_orogenies!(scene)
    uplift_rate = get_diag("continent_orogenic_uplift_rate") .+
        get_diag("subduction_orogenic_uplift_rate")
    uplift_mask = gt_mask(uplift_rate,0.1)
    subx,suby,subfield = zoom_crop_field(uplift_mask)
    contour!(scene,subx,suby,subfield,color=:red)
    
    exposed_basement_field = eq_mask(world.geomorphology,exposed_basement)
    subx,suby,subfield = zoom_crop_field(exposed_basement_field)
    contour!(scene,subx,suby,subfield,color=:orange)

    #=ocean_shelf_field = eq_mask(world.geomorphology,ocean_shelf)
    subx,suby,subfield = zoomfield(ocean_shelf_field,ix,iy,nskirt)
    contour!(scene,subx,suby,subfield,color=:brown)=#
    #subduction_uplifting = gt_mask(get_diag("subduction_orogenic_uplift_rate"),0.01)
    #orogenies = active_orogenic_events()
    #text!(scene, orogenies, position = (0,-105),textsize=15)
    return scene
end

nxcolorbar = 5
function setup_plot_field()
    scene = Scene(resolution = (1400, 800)) #(1000, 550))  # GLMakie.lines([-180,180],[0,0],show_axis=false)
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
    heatmap!(scene,xplotcoords .+ xoffset,ycoords .+ yoffset,colorfield,
        colorrange=(maxval,minval),colormap=cmap)#thermometer)#lightrainbow)
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
        contour!(scene,xcoords,ycoords,maskfield,color=:white)
    end
    return scene
end
function plot_add_transport_regimes!(scene)
    maskfield = eq_mask(world.geomorphology,3)
    contour!(scene,xcoords,ycoords,maskfield,color=:red)
    maskfield = eq_mask(world.geomorphology,2)
    contour!(scene,xcoords,ycoords,maskfield,color=:green)
    maskfield = eq_mask(world.geomorphology,1)
    contour!(scene,xcoords,ycoords,maskfield,color=:blue)
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
function plot_add_continent_outlines!(scene)
    continentmask = get_continent_mask()
    contour!(scene,xcoords,ycoords,continentmask,color=:black)
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
    contour!(scene,xcoords,ycoords,seafloor_sediment_thickness,levels=0:20:100,color=:yellow)
    contour!(scene,xcoords,ycoords,seafloor_sediment_thickness,levels=0:1000:5000,color=:white)
    return scene
end
function plot_add_land_sediment_thickness_contours!(scene)
    is_land = gt_mask(world.freeboard,0.)
    seafloor_sediment_thickness = world.sediment_thickness .* is_land
    contour!(scene,xcoords,ycoords,seafloor_sediment_thickness,levels=0:1000:5000,color=:green)
    contour!(scene,xcoords,ycoords,seafloor_sediment_thickness,levels=0:200:1000,color=:yellow)
    return scene
end
function plot_add_sediment_thickness_contours!(scene)
    plot_add_sediment_thickness_contours!(scene,0,1000,5000)
    plot_add_sediment_thickness_contours!(scene,0,5,25)
    return scene
end
function plot_add_orogenies!(scene)
    uplift_rate = get_diag("continent_orogenic_uplift_rate") .+
        get_diag("subduction_orogenic_uplift_rate")
    contour!(scene,xcoords,ycoords,uplift_rate,color=:red)

    exposed_basement_field = eq_mask(world.geomorphology,exposed_basement)
    contour!(scene,xcoords,ycoords,exposed_basement_field,color=:blue)

    #subduction_uplifting = gt_mask(get_diag("subduction_orogenic_uplift_rate"),0.01)
    #contour!(scene,xcoords,ycoords,subduction_uplifting,color=:green)
    orogenies = active_orogenic_events()
    text!(scene, orogenies, position = (0,-105),textsize=15)
    return scene
end
function plate_plot_add_outcrop!(scene,plate)
    world_outcrop = eq_mask(world.plateID, plate.plateID)
    plate_projection = projected_plate_maskfield!(world_outcrop,
        plate.rotationmatrix)
    contour!(scene,xcoords,ycoords,plate_projection,color=:white)
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
    scaled_thickness = land_ocean_scale(world.sediment_thickness,0.,5000.,0.,2000.)
    scene = plot_field(scaled_thickness,0.,1.)
    return scene
end
function plot_elevation()
    scaled_thickness = land_ocean_scale(world.freeboard,0.,5000.,-4500.,-2000.)
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
