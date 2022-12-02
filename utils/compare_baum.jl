using NCDatasets
#using GLMakie

function lores_generate_runoff_map( landfrac )

    west_coast_maritime_boost = fill(0.,96,48)
    east_coast_maritime_boost = fill(0.,96,48)
    runoff_map = fill(0.,96,48) 

    q_transect, east_coast_rainfall_penetration, west_coast_rainfall_penetration = 
        lores_generate_zonal_met_transects( )

    for iy in 1:48

        west_coast_maritime_boost[1,iy] = 0.
        if landfrac[1,iy] <= 0.
            west_coast_maritime_boost[1,iy] = 1.
        end
        for ix in 2:96
            if landfrac[ix,iy] <= 0.
                west_coast_maritime_boost[ix,iy] = 1.
                #println(ix," ",1.)
            else
                #=downwind_slope = ( world.freeboard[ix,iy] - world.freeboard[ix-1,iy] ) /
                    delta_x[iy]
                rain_shadow_pull = 0.
                if downwind_slope > 0.
                    rain_shadow_pull = downwind_slope *
                        rain_shadow_slope_parameter
                end=#
                west_coast_maritime_boost[ix,iy] = 
                    west_coast_maritime_boost[ix-1,iy] *
                    exp( - 3.75 / west_coast_rainfall_penetration[iy] )
                    #println(ix," ",maritime_boost[ix,iy])
            end
        end
        if landfrac[1,iy] > 0. # wrap around
            west_coast_maritime_boost[1,iy] = 
                west_coast_maritime_boost[96,iy] *
                exp( - 3.5 / west_coast_rainfall_penetration[iy] )
            ix_next = 2; run_another_loop = true
            while ix_next < 96 && run_another_loop
                if landfrac[ix_next,iy] > 0. # carry along
                    west_coast_maritime_boost[ix_next,iy] = 
                        west_coast_maritime_boost[ix_next-1,iy] *
                        exp( - 3.5 / west_coast_rainfall_penetration[iy] )
                else # its water so were done
                    run_another_loop = false
                end
                ix_next += 1
            end
        end

        east_coast_maritime_boost[96,iy] = 0.
        if landfrac[96,iy] <= 0.
            east_coast_maritime_boost[96,iy] = 1.
        end
        for ix in 95:-1:1
            if landfrac[ix,iy] <= 0.
                east_coast_maritime_boost[ix,iy] = 1.
                #println(ix," ",1.)
            else
                #=downwind_slope = ( world.freeboard[ix,iy] - world.freeboard[ix-1,iy] ) /
                    delta_x[iy]
                rain_shadow_pull = 0.
                if downwind_slope > 0.
                    rain_shadow_pull = downwind_slope *
                        rain_shadow_slope_parameter
                end=#
                east_coast_maritime_boost[ix,iy] = 
                    east_coast_maritime_boost[ix+1,iy] *
                    exp( - 3.75 / east_coast_rainfall_penetration[iy] )
                    #println(ix," ",maritime_boost[ix,iy])
            end
        end
        if landfrac[96,iy] > 0. # wrap around
            east_coast_maritime_boost[96,iy] = 
                east_coast_maritime_boost[1,iy] *
                exp( - 3.5 / east_coast_rainfall_penetration[iy] )
            ix_next = 95; run_another_loop = true
            while ix_next > 1 && run_another_loop
                if landfrac[ix_next,iy] > 0. # carry along
                    east_coast_maritime_boost[ix_next,iy] = 
                        east_coast_maritime_boost[ix_next+1,iy] *
                        exp( - 3.5 / east_coast_rainfall_penetration[iy] )
                else # its water so were done
                    run_another_loop = false
                end
                ix_next -= 1
            end
        end
    end
    #smooth_lores_world!( maritime_boost, maritime_rainfall_smoothing )

    for iy in 1:48
        runoff_map[:,iy] = q_transect[iy] .*
            ( west_coast_maritime_boost[:,iy] + 
            east_coast_maritime_boost[:,iy] )
    end 
    #smooth_lores_world!( runoff_map,maritime_rainfall_smoothing  )
    for ix in 1:96
        for iy in 1:48
            if landfrac[ix,iy] < 0.5
                runoff_map[ix,iy] = 0.
            end
        end
    end
    return runoff_map
end
function lores_generate_zonal_met_transects( )
    q_transect = fill(0.,48)
    west_coast_rainfall_penetration = fill(0.,48)
    east_coast_rainfall_penetration = fill(0.,48)
    for iy in 1:48
        tropical_bulge = runoff_tropical_bulge_base * 
            exp( - lats[iy]^2 / ( runoff_tropics_width )^2 )
        subpolar_bulge = runoff_polar_bulge_base * (
            exp( - ( lats[iy] - lat_polar_bulge_center ) ^2 / ( runoff_subpolar_width )^2 ) +
            exp( - ( lats[iy] + lat_polar_bulge_center ) ^2 / ( runoff_subpolar_width )^2 ) )
        q_transect[iy] = tropical_bulge + minimum_rainfall + subpolar_bulge
        east_coast_rainfall_penetration[iy] = tropical_bulge * 
            east_coast_rainfall_penetration_scale
        west_coast_rainfall_penetration[iy] = subpolar_bulge * 
            west_coast_rainfall_penetration_scale
    end     
    q_transect = q_transect / 3.14e7 * # m3 / s / m2
        1000. * # l / s / m2
        1.e6 # l / km2 s
    return q_transect, east_coast_rainfall_penetration, west_coast_rainfall_penetration
end
function smooth_lores_world!( values, diffcoeff )
    # used for orogeny footprint
    nx = 96; ny = 48
    n_columns = nx * ny
    A = spzeros( n_columns, n_columns )
    for ix in 1:nx
        for iy in 1:ny    
            i_row = ix + ( iy - 1 ) * nx
            iy_fake = Int(floor( iy * 180 / 96 ))
            h_diffcoeffs = [ diffcoeff, diffcoeff ] 
            v_diffcoeffs = [ diffcoeff, diffcoeff ] 

            A[i_row,i_row] = 1.
            # left
            A[i_row,i_row] += h_diffcoeffs[2] # left = -ve 
            i_neighbor_row = ix_lores_left(ix) + ( iy - 1 ) * nx
            A[i_row,i_neighbor_row] = - h_diffcoeffs[2] # only fill one direction
            # right
            A[i_row,i_row] += h_diffcoeffs[1]
            i_neighbor_row = ix_lores_right(ix) + ( iy - 1 ) * nx
            A[i_row,i_neighbor_row] = - h_diffcoeffs[1]
            if iy > 1
                A[i_row,i_row] += v_diffcoeffs[2]
                i_neighbor_row = ix + ( iy - 2 ) * nx
                A[i_row,i_neighbor_row] = - v_diffcoeffs[2]
            end
            if iy < ny
                A[i_row,i_row] += v_diffcoeffs[1]
                i_neighbor_row = ix + ( iy ) * nx
                A[i_row,i_neighbor_row] = - v_diffcoeffs[1]
            end
        end
    end
    
    lu_A = lu(A)

    R = fill(0.,n_columns)# Float64[]; 
    for ix in 1:nx
        for iy in 1:ny
            i_row = ix + ( iy - 1 ) * nx
            #push!(R, values[ix,iy])
            R[i_row] = values[ix,iy]
        end
    end
    new_values_list = lu_A \ R
    for ix in 1:nx
        for iy in 1:ny
            i_row = ix + ( iy - 1 ) * nx
            values[ix,iy] = new_values_list[i_row]
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
function qrunoff_misfit(landfrac,qrunoff,my_qrunoff)
    sum_misfit2 = 0.
    means = [0.,0.]
    for ix in 1:96
        for iy in 1:48
            if landfrac[ix,iy] > 0
                means[1] += qrunoff[ix,iy]
                means[2] += my_qrunoff[ix,iy]
                misfit2 = ( qrunoff[ix,iy] - my_qrunoff[ix,iy] )^2
                sum_misfit2 += misfit2
            end
        end
    end 
    rms = sqrt( sum_misfit2 )
    return means, rms
end

function plot_all()
    baum_directory = output_location * "/ensemble-weathering/ensemble-results/large-ens-convex/"
    cd( baum_directory )
    results_dirs = readdir()
    image_no = 0

    #dir_no = 100
    for dir_no in 1:length(results_dirs)
        results_dir = results_dirs[ dir_no ]
        if results_dir[1:3] == "e.e"
            image_no += 1
            println(results_dir)
            ds = Dataset(results_dir * "/ROF_landfrac.nc")

            landfrac = fill(0.,96,48)
            landfrac[:,:] = replace( ds["landfrac"][:,:], missing => 0. )
            qrunoff = fill(0.,96,48)
            qrunoff[:,:] = replace( ds["QRUNOFF"][:,:], missing => 0. ) .* 1.e6
            lats = fill(0.,48)
            lats[:] = ds["lat"][:]
            lons = fill(0.,96)
            lons[:] = ds["lon"][:] .- 180.
            #Makie.heatmap!(scene,lons,lats,landfrac)
            cmap = :lightrainbow # Reverse(:lightrainbow)
            minval = 0.
            maxval = 8.


            q_transect, east_coast_rainfall_penetration, west_coast_rainfall_penetration =
                lores_generate_zonal_met_transects( )
            Plots.plot(Plots.plot(qrunoff[48,:],lats))
            Plots.plot!(q_transect,lats)

            image_file_name = "img." * lpad( image_no,3,"0") * ".png"
            println(image_file_name)

            scene = setup_double_plot_field() 
            Makie.heatmap!(scene,lons,lats,qrunoff,colormap=cmap,colorrange=(minval,maxval) )
            Makie.contour!(scene,lons,lats,landfrac,colormap=cmap,colorrange=(minval,maxval) )

            my_qrunoff = lores_generate_runoff_map( landfrac ) 
            totals, rms_error = qrunoff_misfit(landfrac,qrunoff,my_qrunoff) 
            println("totals ", totals)
            fraction_string = string( totals[2] / totals[1] )[1:4]
            #Makie.heatmap!(scene,lons,lats .- 190.,maritime_boost)
            Makie.heatmap!(scene,lons,lats .- 190.,my_qrunoff,colormap=cmap,colorrange=(minval,maxval))
            Makie.contour!(scene,lons,lats .- 190.,landfrac,colormap=cmap,colorrange=(minval,maxval))
            text!(scene,results_dir, position = (-180,95),textsize=15)
            text!(scene,"Comparison with Baum et al 2022",position = ( -90, 110 ), textsize=25)
            text!(scene,"Fit / Obs = " * fraction_string, position = ( -20, -320), textsize=15)
            Makie.save( image_file_name, scene )
        end
    end
    run(`ffmpeg -r 2 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p baum_compare_runoff.mp4`)
end

runoff_tropics_width = 12.
runoff_tropical_bulge_base = 1. # m / yr
east_coast_rainfall_penetration_scale = 30.
runoff_polar_bulge_base = 1.5
runoff_subpolar_width = 7.5
lat_polar_bulge_center = 45.
west_coast_rainfall_penetration_scale = 30.
minimum_rainfall = 0.1
maritime_rainfall_smoothing = 0.
