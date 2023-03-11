using NCDatasets
#using GLMakie

function lores_generate_runoff_map(landfrac)

    west_coast_maritime_boost = fill(0.0, 96, 48)
    east_coast_maritime_boost = fill(0.0, 96, 48)
    runoff_map = fill(0.0, 96, 48)

    q_transect, east_coast_rainfall_penetration, west_coast_rainfall_penetration =
        lores_generate_zonal_met_transects()

    for iy in 1:48

        west_coast_maritime_boost[1, iy] = 0.0
        if landfrac[1, iy] <= 0.0
            west_coast_maritime_boost[1, iy] = 1.0
        end
        for ix in 2:96
            if landfrac[ix, iy] <= 0.0
                west_coast_maritime_boost[ix, iy] = 1.0
                #println(ix," ",1.)
            else
                #=downwind_slope = ( world.freeboard[ix,iy] - world.freeboard[ix-1,iy] ) /
                    delta_x[iy]
                rain_shadow_pull = 0.
                if downwind_slope > 0.
                    rain_shadow_pull = downwind_slope *
                        rain_shadow_slope_parameter
                end=#
                west_coast_maritime_boost[ix, iy] =
                    west_coast_maritime_boost[ix-1, iy] *
                    exp(-3.75 / west_coast_rainfall_penetration[iy])
                #println(ix," ",maritime_boost[ix,iy])
            end
        end
        if landfrac[1, iy] > 0.0 # wrap around
            west_coast_maritime_boost[1, iy] =
                west_coast_maritime_boost[96, iy] *
                exp(-3.5 / west_coast_rainfall_penetration[iy])
            ix_next = 2
            run_another_loop = true
            while ix_next < 96 && run_another_loop
                if landfrac[ix_next, iy] > 0.0 # carry along
                    west_coast_maritime_boost[ix_next, iy] =
                        west_coast_maritime_boost[ix_next-1, iy] *
                        exp(-3.5 / west_coast_rainfall_penetration[iy])
                else # its water so were done
                    run_another_loop = false
                end
                ix_next += 1
            end
        end

        east_coast_maritime_boost[96, iy] = 0.0
        if landfrac[96, iy] <= 0.0
            east_coast_maritime_boost[96, iy] = 1.0
        end
        for ix in 95:-1:1
            if landfrac[ix, iy] <= 0.0
                east_coast_maritime_boost[ix, iy] = 1.0
                #println(ix," ",1.)
            else
                #=downwind_slope = ( world.freeboard[ix,iy] - world.freeboard[ix-1,iy] ) /
                    delta_x[iy]
                rain_shadow_pull = 0.
                if downwind_slope > 0.
                    rain_shadow_pull = downwind_slope *
                        rain_shadow_slope_parameter
                end=#
                east_coast_maritime_boost[ix, iy] =
                    east_coast_maritime_boost[ix+1, iy] *
                    exp(-3.75 / east_coast_rainfall_penetration[iy])
                #println(ix," ",maritime_boost[ix,iy])
            end
        end
        if landfrac[96, iy] > 0.0 # wrap around
            east_coast_maritime_boost[96, iy] =
                east_coast_maritime_boost[1, iy] *
                exp(-3.5 / east_coast_rainfall_penetration[iy])
            ix_next = 95
            run_another_loop = true
            while ix_next > 1 && run_another_loop
                if landfrac[ix_next, iy] > 0.0 # carry along
                    east_coast_maritime_boost[ix_next, iy] =
                        east_coast_maritime_boost[ix_next+1, iy] *
                        exp(-3.5 / east_coast_rainfall_penetration[iy])
                else # its water so were done
                    run_another_loop = false
                end
                ix_next -= 1
            end
        end
    end
    #smooth_lores_world!( maritime_boost, maritime_rainfall_smoothing )

    for iy in 1:48
        runoff_map[:, iy] = q_transect[iy] .*
                            (west_coast_maritime_boost[:, iy] +
                             east_coast_maritime_boost[:, iy])
    end
    #smooth_lores_world!( runoff_map,maritime_rainfall_smoothing  )
    for ix in 1:96
        for iy in 1:48
            if landfrac[ix, iy] < 0.5
                runoff_map[ix, iy] = 0.0
            end
        end
    end
    return runoff_map
end
function lores_generate_zonal_met_transects()
    q_transect = fill(0.0, 48)
    west_coast_rainfall_penetration = fill(0.0, 48)
    east_coast_rainfall_penetration = fill(0.0, 48)
    for iy in 1:48
        tropical_bulge = runoff_tropical_bulge_base *
                         exp(-lats[iy]^2 / (runoff_tropics_width)^2)
        subpolar_bulge = runoff_polar_bulge_base * (
            exp(-(lats[iy] - lat_polar_bulge_center)^2 / (runoff_subpolar_width)^2) +
            exp(-(lats[iy] + lat_polar_bulge_center)^2 / (runoff_subpolar_width)^2))
        q_transect[iy] = tropical_bulge + minimum_rainfall + subpolar_bulge
        east_coast_rainfall_penetration[iy] = tropical_bulge *
                                              east_coast_rainfall_penetration_scale
        west_coast_rainfall_penetration[iy] = subpolar_bulge *
                                              west_coast_rainfall_penetration_scale
    end
    q_transect = q_transect / 3.14e7 * # m3 / s / m2
                 1000.0 * # l / s / m2
                 1.e6 # l / km2 s
    return q_transect, east_coast_rainfall_penetration, west_coast_rainfall_penetration
end
function smooth_lores_world!(values, diffcoeff)
    # used for orogeny footprint
    nx = 96
    ny = 48
    n_columns = nx * ny
    A = spzeros(n_columns, n_columns)
    for ix in 1:nx
        for iy in 1:ny
            i_row = ix + (iy - 1) * nx
            iy_fake = Int(floor(iy * 180 / 96))
            h_diffcoeffs = [diffcoeff, diffcoeff]
            v_diffcoeffs = [diffcoeff, diffcoeff]

            A[i_row, i_row] = 1.0
            # left
            A[i_row, i_row] += h_diffcoeffs[2] # left = -ve 
            i_neighbor_row = ix_lores_left(ix) + (iy - 1) * nx
            A[i_row, i_neighbor_row] = -h_diffcoeffs[2] # only fill one direction
            # right
            A[i_row, i_row] += h_diffcoeffs[1]
            i_neighbor_row = ix_lores_right(ix) + (iy - 1) * nx
            A[i_row, i_neighbor_row] = -h_diffcoeffs[1]
            if iy > 1
                A[i_row, i_row] += v_diffcoeffs[2]
                i_neighbor_row = ix + (iy - 2) * nx
                A[i_row, i_neighbor_row] = -v_diffcoeffs[2]
            end
            if iy < ny
                A[i_row, i_row] += v_diffcoeffs[1]
                i_neighbor_row = ix + (iy) * nx
                A[i_row, i_neighbor_row] = -v_diffcoeffs[1]
            end
        end
    end

    lu_A = lu(A)

    R = fill(0.0, n_columns)# Float64[]; 
    for ix in 1:nx
        for iy in 1:ny
            i_row = ix + (iy - 1) * nx
            #push!(R, values[ix,iy])
            R[i_row] = values[ix, iy]
        end
    end
    new_values_list = lu_A \ R
    for ix in 1:nx
        for iy in 1:ny
            i_row = ix + (iy - 1) * nx
            values[ix, iy] = new_values_list[i_row]
        end
    end
end
function ix_lores_left(ix)
    nx = 96
    if ix == 1
        ix_l = nx
    else
        ix_l = ix - 1
    end
    return ix_l
end
function ix_lores_right(ix)
    nx = 96
    if ix == nx
        ix_r = 1
    else
        ix_r = ix + 1
    end
    return ix_r
end
function qrunoff_misfit(landfrac, qrunoff, my_qrunoff)
    sum_misfit2 = 0.0
    means = [0.0, 0.0]
    for ix in 1:96
        for iy in 1:48
            if landfrac[ix, iy] > 0
                means[1] += qrunoff[ix, iy]
                means[2] += my_qrunoff[ix, iy]
                misfit2 = (qrunoff[ix, iy] - my_qrunoff[ix, iy])^2
                sum_misfit2 += misfit2
            end
        end
    end
    rms = sqrt(sum_misfit2)
    return means, rms
end

function plot_all()
    baum_directory = output_location * "/ensemble-weathering/ensemble-results/large-ens-convex/"
    cd(baum_directory)
    results_dirs = readdir()
    image_no = 0

    #dir_no = 100
    for dir_no in 1:length(results_dirs)
        results_dir = results_dirs[dir_no]
        if results_dir[1:3] == "e.e"
            image_no += 1
            println(results_dir)
            ds = Dataset(results_dir * "/ROF_landfrac.nc")

            landfrac = fill(0.0, 96, 48)
            landfrac[:, :] = replace(ds["landfrac"][:, :], missing => 0.0)
            qrunoff = fill(0.0, 96, 48)
            qrunoff[:, :] = replace(ds["QRUNOFF"][:, :], missing => 0.0) .* # mm/sec 
                            1.e-3 .* # m / sec
                            1.e6 .* # m3 / km2 sec 
                            1.e3 # l / km2 sec
            lats = fill(0.0, 48)
            lats[:] = ds["lat"][:]
            lons = fill(0.0, 96)
            lons[:] = ds["lon"][:] .- 180.0
            #Makie.heatmap!(scene,lons,lats,landfrac)
            cmap = :lightrainbow # Reverse(:lightrainbow)
            minval = 0.0
            maxval = 8.0

            q_transect, east_coast_rainfall_penetration, west_coast_rainfall_penetration =
                lores_generate_zonal_met_transects()
            Plots.plot(Plots.plot(qrunoff[48, :], lats))
            Plots.plot!(q_transect, lats)

            image_file_name = "img." * lpad(image_no, 3, "0") * ".png"
            println(image_file_name)

            scene = setup_double_plot_field()
            Makie.heatmap!(scene, lons, lats, qrunoff, colormap=cmap, colorrange=(minval, maxval))
            Makie.contour!(scene, lons, lats, landfrac, colormap=cmap, colorrange=(minval, maxval))

            my_qrunoff = lores_generate_runoff_map(landfrac)
            totals, rms_error = qrunoff_misfit(landfrac, qrunoff, my_qrunoff)
            println("totals ", totals)
            fraction_string = string(totals[2] / totals[1])[1:4]
            #Makie.heatmap!(scene,lons,lats .- 190.,maritime_boost)
            Makie.heatmap!(scene, lons, lats .- 190.0, my_qrunoff, colormap=cmap, colorrange=(minval, maxval))
            Makie.contour!(scene, lons, lats .- 190.0, landfrac, colormap=cmap, colorrange=(minval, maxval))
            text!(scene, results_dir, position=(-180, 95), textsize=15)
            text!(scene, "Runoff comparison with Baum et al 2022", position=(-180, 110), textsize=25)
            text!(scene, "Gridplates fit / Baum = " * fraction_string, position=(-20, -320), textsize=15)
            Makie.save(image_file_name, scene)
        end
    end
    run(`ffmpeg -r 2 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p baum_compare_runoff.mp4`)
end

runoff_tropics_width = 12.0
runoff_tropical_bulge_base = 1.0 # m / yr
east_coast_rainfall_penetration_scale = 30.0
runoff_polar_bulge_base = 1.5
runoff_subpolar_width = 7.5
lat_polar_bulge_center = 45.0
west_coast_rainfall_penetration_scale = 30.0
minimum_rainfall = 0.1
maritime_rainfall_smoothing = 0.0

function read_scotese_csv(file_name)
    twisted_elevations = fill(0.0, nx, ny)
    f = open(file_name)
    lines = readlines(f)
    for line in lines[2:end]
        words = split(line, ",")
        ix = parse(Int, words[1]) + 180
        iy = parse(Int, words[2]) + 90
        if ix > 0 && iy > 0
            twisted_elevations[ix, iy] = parse.(Float64, words[3])
        end
    end
    return twisted_elevations
end
function translate_scotese_ptemp_files()
    cd("../IceMaps_CRS_v23019")
    file_list = readdir()
    for file_name in file_list
        paleo_temp_modern = read_scotese_csv(file_name)
        #paleo_temp = fill_world_orogeny( paleo_temp_modern )
        file_age = file_name[1:3]
        outfile_name = code_base_directory * 
            "/drivers/scotese_paleo_temp_files_modern/scotese_ptemp." *
            file_age * "Ma.bson"
        BSON.@save outfile_name paleo_temp_modern
    end
    cd( code_base_directory )
end


function write_scotese_csv_files()
    cd("../Scotese_paleoelevation_modern_locations")
    file_list = readdir()
    image_number = 0
    for file_name in file_list
        if file_name[end-2:end] == "csv"
            age_string = file_name[1:3]
            age = parse(Float64, age_string)
            create_everything(age)
            println("reading ", age)
            scotese_elevation = get_scotese_mountains()
            world.freeboard[:, :] = scotese_elevation
            output_file_name = file_name[1:14] * "Merdith_list.csv"
            println("writing ", output_file_name)
            rm(output_file_name, force=true)
            f = open(output_file_name, "w")
            println(f, "Merdith_longtitude,Merdith_latitude,paleoelevation")
            for iy in ny:-1:1
                for ix in 1:nx
                    println(f, xcoords[ix], ",", ycoords[iy], ",", scotese_elevation[ix, iy])
                    #if ix == 160 && iy == 170
                    #    error("chunk")
                    #end
                end
            end
            close(f)
            bson_filename = "scotese_paleo." * age_string * "Ma.bson"
            rm(bson_filename, force=true)
            BSON.@save bson_filename scotese_elevation
            world.freeboard = scotese_elevation
            scene = plot_elevation()
            image_number += 1
            image_file_name = "img." * lpad(image_number, 3, "0") * ".png"
            Makie.save(image_file_name, scene)
        end
    end
    mp4_file = "scotese_elevs_merdith_positions.mp4"
    println("compiling ", mp4_file)
    rm(mp4_file, force=true)
    run(`ffmpeg -r 2 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p $mp4_file`)
    cd(code_base_directory)
end
function read_scotese_ophiolite_database()
    filename = "../Ophiolite_Database_v23009a.csv"
    f = open(filename)
    ophiolite_start_time_field = fill(0.0, nx, ny)
    lines = readlines(f)
    for line in lines[2:end]
        words = split(line, ",")
        latitude = parse(Float16, words[3])
        iy_location = search_sorted_nearest(ycoords, latitude)
        longitude = parse(Float16, words[4])
        ix_location = search_sorted_nearest(ycoords, longitude)
        start_time = parse(Float16, words[6])
        ophiolite_start_time_field[ix_location, iy_location] =
            start_time
    end
    return ophiolite_start_time_field
end
function read_all_scotese_twisted_elevation_files()
    cd("../Scotese_paleoelevation_modern_locations")
    file_list = readdir()
    for file_name in file_list[2:end]
        age_string = file_name[1:3]
        age = parse(Float64, age_string)
        elevation = read_scotese_csv(file_name)
        out_file_name = code_base_directory *
                        "/drivers/scotese_elevation_files_modern/scotese_elev_today." *
                        age_string * "Ma.bson"
        rm(out_file_name, force=true)
        BSON.@save out_file_name elevation
        println(file_name)
    end
    cd(code_base_directory)
end
function read_GLIM_lithology_map()
    lithology_grid = fill(0, nx, ny)
    f = open("../glim_wgs84_0point5deg.txt.asc")
    lines = readlines(f)
    for my_grid_iy in 1:180
        file_iy = 5 + my_grid_iy * 2
        line = lines[file_iy]
        words = split(line, " ")
        for my_grid_ix in 1:360
            file_ix = my_grid_ix * 2 - 1
            lithology_grid[my_grid_ix, 181-my_grid_iy] = parse(Int, words[file_ix])
        end
    end
    out_file_name = code_base_directory *
                    "/data/glim_lithology.bson"
    rm(out_file_name, force=true)
    BSON.@save out_file_name lithology_grid
    #return lithology_grid
end
function read_NCEP_runoff()
    ds = Dataset("../runof.mon.ltm.nc")
    nx_ncep = 192
    ny_ncep = 94
    xc_ncep = fill(0.0, nx_ncep)
    xc_ncep[:] = ds["lon"][:]
    yc_ncep = fill(0.0, ny_ncep)
    yc_ncep[:] = ds["lat"][:]
    #runoff = fill(0.,nx_ncep,ny_ncep,12)
    runoff = ds["runof"][:, :, :]
    runoff_annual = fill(0.0, nx_ncep, ny_ncep)
    for ix in 1:nx_ncep
        for iy in 1:ny_ncep
            ntimes = 0
            for itime in 1:12
                #println(ix,iy,runoff[ix,iy,itime])
                if ismissing(runoff[ix, iy, itime])
                else
                    ntimes += 1
                    runoff_annual[ix, iy] += runoff[ix, iy, itime]
                end
            end
            if ntimes > 0
                runoff_annual[ix, iy] /= ntimes
            end
        end
    end
    runoff_highres = fill(0.0, nx, ny)
    for ix in 1:nx
        for iy in 1:ny
            x_hires = xcoords[ix]
            if x_hires < 0
                x_hires += 360.0
            end

            ix_ncep_left = search_sorted_below(xc_ncep, x_hires)
            if ix_ncep_left < nx_ncep
                ix_ncep_right = ix_ncep_left + 1
            else
                ix_ncep_right = 1
            end
            #fractions = fill(0.,nx_ncep,ny_ncep)
            frac_right = (x_hires - xc_ncep[ix_ncep_left]) /
                         (xc_ncep[ix_ncep_right] - xc_ncep[ix_ncep_left])
            frac_left = 1.0 - frac_right
            y_hires = ycoords[iy]
            iy_ncep_below = max(1, search_sorted_below(yc_ncep, y_hires))
            #println( y_hires )
            if iy_ncep_below > 1
                iy_ncep_above = iy_ncep_below - 1
            else
                iy_ncep_above = iy_ncep_below
            end
            iy_ncep_above
            frac_above = (y_hires - yc_ncep[iy_ncep_below]) /
                         (yc_ncep[iy_ncep_above] - yc_ncep[iy_ncep_below])
            frac_below = 1.0 - frac_above
            #for ixpull in ix_ncep_left:ix_ncep_right
            #    for iypull in iy_ncep_below:iy_ncep_above
            runoff_highres[ix, iy] =
                runoff_annual[ix_ncep_left, iy_ncep_below] *
                frac_left * frac_below +
                runoff_annual[ix_ncep_right, iy_ncep_below] *
                frac_right * frac_below +
                runoff_annual[ix_ncep_left, iy_ncep_above] *
                frac_left * frac_above +
                runoff_annual[ix_ncep_right, iy_ncep_above] *
                frac_right * frac_above
            if runoff_highres[ix, iy] == runoff_highres[ix, iy]
            else
                runoff_highres[ix, iy] = 0.0
            end
        end
    end
    out_file_name = code_base_directory *
                    "/data/ncep_runoff.bson"
    rm(out_file_name, force=true)
    BSON.@save out_file_name runoff_highres
end
function read_laske_sedmap()
    cd("..")
    f = open("laske.sedmap.xyz")
    lines = readlines(f)
    sedthick = fill(0.0, nx, ny)
    for line in lines
        words = split(line)
        longt = parse(Float64, words[1])
        lat = parse(Float64, words[2])
        if abs(longt - floor(longt) - 0.5) < 0.01 && abs(lat - floor(lat) - 0.5) < 0.01
            ix = Int(longt + 0.5) + 180
            iy = 90 + Int(lat + 0.5)
            if iy == 0
                error("[longt,lat][ix,iy]", longt, " ", lat, " ", ix, " ", iy)
            end
            if words[3] != "NaN"
                sedthick[ix, iy] = parse(Float64, words[3])
            end
        end
    end
    out_file_name = code_base_directory *
                    "/data/laske_sedmap.bson"
    rm(out_file_name, force=true)
    BSON.@save out_file_name sedthick
    return field
end
function build_real_world()
    global world = create_world(0.0)
    BSON.@load "data/glim_lithology.bson" lithology_grid
    BSON.@load "data/etopo.bson" etopo
    #BSON.@load "data/laske_sedthick.bson" sedthick
    #etopo = read_flip_csv(code_base_directory * "/data/" * "etopo.csv")

    for ix in 1:nx
        for iy in 1:ny
            if lithology_grid[ix, iy] in [3, 5, 6]
                world.sediment_thickness[ix, iy] = 100.0
                world.geomorphology[ix, iy] = sedimented_land
                if lithology_grid[ix, iy] in [5] # unconsolidated or mixed sed 
                    world.sediment_fractions[ix, iy, reactive_clay] = 0.5
                elseif lithology_grid[ix, iy] in [3] # siliciclasic sed
                    world.sediment_fractions[ix, iy, reactive_clay] = 0.8
                else
                    world.sediment_fractions[ix, iy, reactive_clay] = 0.2
                end
                world.sediment_fractions[ix, iy, CaCO3_sediment] =
                    1.0 - world.sediment_fractions[ix, iy, reactive_clay]
            elseif lithology_grid[ix, iy] in [1, 2, 4, 7, 8, 9, 10]
                world.sediment_thickness[ix, iy] = 0.0
                world.geomorphology[ix, iy] = exposed_basement
                if lithology_grid[ix, iy] in [2, 4] # basic volcanic, plutonic rocks
                    world.crust_composition[ix, iy] = 0.0 # mafic             # mafic# mafic
                elseif lithology_grid[ix, iy] in [1, 3, 13]
                    world.crust_composition[ix, iy] = 0.5 # mixed 
                else
                    world.crust_composition[ix, iy] = 1.0
                end
                world.crust_density[ix, iy] =
                    rho_continent_crust * (1.0 - world.crust_composition[ix, iy]) +
                    rho_ocean_crust * world.crust_composition[ix, iy]
            end
            #hack_crust_thickness_to_elevation(ix, iy, etopo[ix, iy])

        end
    end
end

function regrid_scotese_files()
    image_number = 0
    age = world.age
    while age > 0
        age -= main_time_step
        create_everything(age)
        scotese_elevation = get_scotese_mountains()
        image_number += 1
        image_file_name = "../Scotese_paleoelevation_Merdith_locations/img." * lpad(image_number, 3, "0") * ".png"
        scene = plot_field(scotese_elevation, -4000.0, 4000.0)
        plot_add_continent_outlines!(scene)
        Makie.save(image_file_name, scene)
        age_string = lpad(Int(world.age), 3, "0")
        output_file_name = "../Scotese_paleoelevation_Merdith_locations/" * age_string * "Ma_Merdith_list.csv"
        println("writing ", output_file_name)
        rm(output_file_name, force=true)
        f = open(output_file_name, "w")
        println(f, "Merdith_longtitude,Merdith_latitude,paleoelevation")
        for iy in ny:-1:1
            for ix in 1:nx
                println(f, xcoords[ix], ",", ycoords[iy], ",", scotese_elevation[ix, iy])
                #if ix == 160 && iy == 170
                #    error("chunk")
                #end
            end
        end
        close(f)
        bson_filename = "../Scotese_paleoelevation_Merdith_locations/scotese_paleo." * age_string * ".bson"
        rm(bson_filename, force=true)
        BSON.@save bson_filename scotese_elevation
    end
    cd("../Scotese_paleoelevation_Merdith_locations")
    mp4_file = "scotese_paleo_elev_merdith.mp4"
    println("compiling ", mp4_file)
    rm(mp4_file, force=true)
    run(`ffmpeg -r 2 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p $mp4_file`)
    cd(code_base_directory)
end
