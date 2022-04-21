function land_sediment_transport!( elevation_field )
    # zeros and sets land_sediment_deposition_rate, crust_erosion_rate,
    # land_orogenic_clay_flux, coastal_orogenic_clay_flux, crust_erosion_rate
    # 
    verbose = false
    set_diag("land_sediment_deposition_rate",fill(0.,nx,ny))
    set_diag("crust_erosion_rate",fill(0.,nx,ny))
    #set_diag("crust_clay_source_rate",fill(0.,nx,ny))
    
    generate_orogenic_erosion_fluxes()
    # zeros then sets land_orogenic_clay_flux, coastal_orogenic_clay_flux,
    # and crust_erosion_rate.  
    orogenic_sources = get_diag("land_orogenic_clay_flux") # m / Myr
    
    ocean_sink = generate_mask_field(world.crust_type, ocean_crust )
    # ocean_sink = fill(0.,nx,ny) # hack for turning off fluxes to oceans
    subaereal_mask = generate_mask_field(world.surface_type, sedimented_land)
    #= 
    orogenic_sources = fill(0.,nx,ny); orogenic_sources[2,2] = 1.
    subaereal_mask = fill(0.,nx,ny); subaereal_mask[2:4,2] .= 1.
    ocean_sink = fill(0.,nx,ny)
    elevation_field = fill(0.,nx,ny); elevation_field[2:4,2] .= 1.
    =#
    subaereal_blobs = get_blobs( subaereal_mask )
    original_elevation_field = deepcopy( elevation_field )

    Threads.@threads for i_blob in 1:length(subaereal_blobs) 
        #if verbose == true
        #    println("beginning ",i_blob)
        #end
        subaereal_blob = subaereal_blobs[i_blob]
        land_area_sediment_transport!( elevation_field, 
            subaereal_blob, orogenic_sources, ocean_sink )
        for ix in 1:nx
            for iy in 1:ny
                if subaereal_blob[ix,iy] == 1
                    sediment_deposition_rate = 
                        ( elevation_field[ix,iy] - original_elevation_field[ix,iy] ) /
                        time_step / sediment_freeboard_expression
                    set_diag("land_sediment_deposition_rate", ix,iy, sediment_deposition_rate )
                end
            end
        end 
    end # blob 
end
function land_area_sediment_transport!( elevation_field, 
    subaereal_blob, incoming_flux, ocean_sink )
    list_pos = fill(0,nx,ny)
    list_pos_filled = 0
    list_x = []; list_y = []
    original_elevation_list = []
    for ix in 1:nx
        for iy in 1:ny
            if subaereal_blob[ix,iy] == 1
                list_pos_filled += 1
                list_pos[ix,iy] = list_pos_filled
                push!( list_x, ix ); push!( list_y, iy )
            end
        end
    end
    if list_pos_filled == 1
        ix = list_x[1]; iy = list_y[1]
        h_diffcoeffs = get_horizontal_diffcoeff(land_base_diffcoeff,iy) .* sediment_freeboard_expression
        v_diffcoeffs = get_vertical_diffcoeff(land_base_diffcoeff,iy) .* sediment_freeboard_expression
        A = Float64[ 1. ]
        #runoff_loss = 0.
        if ocean_sink[ix_left(ix),iy] == 1 
            A[1] += h_diffcoeffs[2] 
            #runoff_loss += elevation_field[ix,iy] * h_diffcoeffs[2] # meters / step
        end
        if ocean_sink[ix_right(ix),iy] == 1 
            A[1] += h_diffcoeffs[1]
            #runoff_loss += elevation_field[ix,iy] * h_diffcoeffs[1] # meters / step
        end
        if iy > 1
            if ocean_sink[ix,iy-1] == 1 
                A[1] += v_diffcoeffs[2]
            #runoff_loss += elevation_field[ix,iy] * v_diffcoeffs[2] # meters / step
            end
        end
        if iy < ny
            if ocean_sink[ix,iy+1] == 1 
                A[1] += v_diffcoeffs[1]
            #runoff_loss += elevation_field[ix,iy] * v_diffcoeffs[1] # meters / step
            end
        end
        R = Float64[ elevation_field[ix,iy] +
            incoming_flux[ix,iy] * # meters/Myr
                time_step * # meters / step
                sediment_freeboard_expression ]
        new_elevation = A \ R
        elevation_field[ix,iy] = new_elevation
    else # more than one point
        A = spzeros( list_pos_filled, list_pos_filled )
        R = Float64[]; original_elevation_list = []
        for i_pos = 1:list_pos_filled
            ix = list_x[i_pos]; iy = list_y[i_pos]
                h_diffcoeffs = get_horizontal_diffcoeff(land_base_diffcoeff,iy) .* sediment_freeboard_expression
                v_diffcoeffs = get_vertical_diffcoeff(land_base_diffcoeff,iy) .* sediment_freeboard_expression
            push!(R, elevation_field[ix,iy])
            push!(original_elevation_list, elevation_field[ix,iy])
            A[i_pos,i_pos] = 1.
            if subaereal_blob[ix_left(ix),iy] == 1 
                A[i_pos,i_pos] += h_diffcoeffs[2] # left = -ve 
                i_neighbor_pos = list_pos[ix_left(ix),iy]
                A[i_pos,i_neighbor_pos] = - h_diffcoeffs[2] # only fill one direction
            elseif ocean_sink[ix_left(ix),iy] == 1 
                A[i_pos,i_pos] += h_diffcoeffs[2]
            end
            if subaereal_blob[ix_right(ix),iy] == 1 
                A[i_pos,i_pos] += h_diffcoeffs[1]
                i_neighbor_pos = list_pos[ix_right(ix),iy]
                A[i_pos,i_neighbor_pos] = - h_diffcoeffs[1]
            elseif ocean_sink[ix_right(ix),iy] == 1 
                A[i_pos,i_pos] += h_diffcoeffs[1]
            end
            if iy > 1
                if subaereal_blob[ix,iy-1] == 1
                    A[i_pos,i_pos] += v_diffcoeffs[2]
                    i_neighbor_pos = list_pos[ix,iy-1]
                    A[i_pos,i_neighbor_pos] = - v_diffcoeffs[2]
                elseif ocean_sink[ix,iy-1] == 1 
                    A[i_pos,i_pos] += v_diffcoeffs[2]
                end
            end
            if iy < ny
                if subaereal_blob[ix,iy+1] == 1
                    A[i_pos,i_pos] += v_diffcoeffs[1]
                    i_neighbor_pos = list_pos[ix,iy+1]
                    A[i_pos,i_neighbor_pos] = - v_diffcoeffs[1]
                elseif ocean_sink[ix,iy+1] == 1 
                    A[i_pos,i_pos] += v_diffcoeffs[1]
                end
            end
            if incoming_flux[ix,iy] != 0. 
                R[end] += incoming_flux[ix,iy] * # meters / Myr
                    time_step * sediment_freeboard_expression # meters
            end
        end

        new_elevation_list = lu(A) \ R

        for i_pos = 1:list_pos_filled
            ix = list_x[i_pos]; iy = list_y[i_pos]
            #sed_depo = ( new_elevation_list[i_pos] - elevation_field[ix,iy] ) /
            #    time_step # meters / Myr
            #set_diag("sediment_deposition_rate",ix,iy,sed_depo)
            elevation_field[ix,iy] = new_elevation_list[i_pos]
       end
    end # diffuse or no
end
function check_subaereal_exposure()
    n_failed = 0
    sediment_deposition_meters = ( get_diag("land_sediment_deposition_rate") .- 
        get_diag("aolean_clay_erosion_rate") ) .* time_step
    for ix in 1:nx
        for iy in 1:ny
            if world.surface_type[ix,iy] == sedimented_land && 
                world.sediment_thickness[ix,iy] + sediment_deposition_meters[ix,iy] < 0.
                world.surface_type[ix,iy] = exposed_basement # call it denuded
                set_diag("surface_type_changes",ix,iy,land_soil_denuded)
                # leftover sediment gets turned into orogenic fluxes
                #println("stripping point ",ix," ",iy)
                n_failed += 1
            end
        end
    end
    return n_failed
end

function land_sediment_fraction_transport( 
        original_elevation_field, new_elevation_field )

    orogenic_source = get_diag("land_orogenic_clay_flux") # m / Myr
    CaCO3_source = fill(0.,nx,ny) # because its on land
    ocean_sink = generate_mask_field(world.crust_type, ocean_crust )
    # ocean_sink = fill(0.,nx,ny) # hack for turning off fluxes to oceans
    subaereal_mask = generate_mask_field(world.surface_type, sedimented_land)
    #= hack simple condition
    orogenic_source = fill(0.,nx,ny); orogenic_source[2,2] = 1.
    subaereal_mask = fill(0.,nx,ny); subaereal_mask[2:4,2] .= 1.
    ocean_sink = fill(0.,nx,ny)
    original_elevation_field = fill(0.,nx,ny); original_elevation_field[2:4,2] .= 1.
    =#
    incoming_fluxes = [ orogenic_source, CaCO3_source ]
    subaereal_blobs = get_blobs( subaereal_mask )
    new_sediment_thickness = world.sediment_thickness .+
        ( new_elevation_field .- original_elevation_field ) / 
        sediment_freeboard_expression
    Threads.@threads for i_blob in 1:length(subaereal_blobs) 
        #if verbose == true
        #    println("beginning ",i_blob)
        #end
        subaereal_blob = subaereal_blobs[i_blob]
        land_area_sediment_fraction_transport( 
            new_elevation_field, new_sediment_thickness,
            subaereal_blob, incoming_fluxes, ocean_sink ) # sets world.sediment_fractions on land points
    end # blob 
end
function land_area_sediment_fraction_transport( 
    new_elevation_field, new_sediment_thickness,
    subaereal_blob, incoming_fluxes, ocean_sink ) # sets world.sediment_fractions on land points
    # fractions combined so they can share the transportation matrix
    list_pos = fill(0,nx,ny)
    list_pos_filled = 0
    list_x = []; list_y = []
    for ix in 1:nx
        for iy in 1:ny
            if subaereal_blob[ix,iy] == 1
                list_pos_filled += 1
                list_pos[ix,iy] = list_pos_filled
                push!( list_x, ix ); push!( list_y, iy )
            end
        end
    end
    if list_pos_filled > 1 # more than one point
        A = spzeros( list_pos_filled, list_pos_filled )
        for i_pos = 1:list_pos_filled
            ix = list_x[i_pos]; iy = list_y[i_pos]
            h_diffcoeffs = get_horizontal_diffcoeff(land_base_diffcoeff,iy) 
            # no isostatic rebound factor here -- getting D dt / dx, gives thickness transport
            v_diffcoeffs = get_vertical_diffcoeff(land_base_diffcoeff,iy) 
            # Inv2' = Inv2 + [ D ( H2 - H2 ) * ( C2 + C1 ) / 2 + Src ]
            # Inv2' = Inv2 + [ D ( H2 - H2 ) * ( Inv2'/H2' + Inv1'/H1' ) / 2 + Src ]
            A[i_pos,i_pos] = 1.
            if subaereal_blob[ix_left(ix),iy] == 1 
                A[i_pos,i_pos] -= h_diffcoeffs[2] * 
                    ( new_elevation_field[ix_left(ix),iy] - new_elevation_field[ix,iy] ) / 
                    2 / new_sediment_thickness[ix,iy]
                i_neighbor_pos = list_pos[ix_left(ix),iy]
                A[i_pos,i_neighbor_pos] -= h_diffcoeffs[2] * 
                    ( new_elevation_field[ix_left(ix),iy] - new_elevation_field[ix,iy] ) /  
                    2 / new_sediment_thickness[ix_left(ix),iy]
            elseif ocean_sink[ix_left(ix),iy] == 1 
                A[i_pos,i_pos] -= h_diffcoeffs[2] * 
                    ( 0. - new_elevation_field[ix,iy] ) /  
                    2 / new_sediment_thickness[ix,iy]
            end
            if subaereal_blob[ix_right(ix),iy] == 1 
                A[i_pos,i_pos] -= h_diffcoeffs[1] * 
                    ( new_elevation_field[ix_right(ix),iy] - new_elevation_field[ix,iy] ) /  
                    2 / new_sediment_thickness[ix,iy]
                i_neighbor_pos = list_pos[ix_right(ix),iy]
                A[i_pos,i_neighbor_pos] -= h_diffcoeffs[1] * 
                    ( new_elevation_field[ix_right(ix),iy] - new_elevation_field[ix,iy] ) /  
                    2 / new_sediment_thickness[ix_right(ix),iy]
            elseif ocean_sink[ix_right(ix),iy] == 1 
                A[i_pos,i_pos] -= h_diffcoeffs[1] * 
                    ( 0. - new_elevation_field[ix,iy] ) /  
                    2 / new_sediment_thickness[ix,iy]
            end   
            if iy == 1
                ixn = ( ix + 179 ) % 360 + 1
                if subaereal_blob[ixn,iy] == 1
                    A[i_pos,i_pos] += v_diffcoeffs[2] * 
                        ( new_elevation_field[ixn,iy] - new_elevation_field[ix,iy] ) /  
                        2 / new_sediment_thickness[ix,iy]
                    i_neighbor_pos = list_pos[ixn,iy]
                    A[i_pos,i_neighbor_pos] = - v_diffcoeffs[2] * 
                        ( new_elevation_field[ixn,iy] - new_elevation_field[ix,iy] ) /  
                        2 / new_sediment_thickness[ixn,iy]
                end
            else # iy > 1
                if subaereal_blob[ix,iy-1] == 1 
                    A[i_pos,i_pos] -= v_diffcoeffs[2] * 
                        ( new_elevation_field[ix,iy-1] - new_elevation_field[ix,iy] ) /  
                        2 / new_sediment_thickness[ix,iy]
                    i_neighbor_pos = list_pos[ix,iy-1]
                    A[i_pos,i_neighbor_pos] -= v_diffcoeffs[2] * 
                        ( new_elevation_field[ix,iy-1] - new_elevation_field[ix,iy] ) /  
                        2 / new_sediment_thickness[ix,iy-1]
                elseif ocean_sink[ix,iy-1] == 1 
                    A[i_pos,i_pos] -= v_diffcoeffs[2] * 
                        ( 0. - new_elevation_field[ix,iy] ) /  
                        2 / new_sediment_thickness[ix,iy]
                end    
            end        
            if iy ==  ny    
                ixn = ( ix + 179 ) % 360 + 1
                if subaereal_blob[ixn,iy] == 1
                    A[i_pos,i_pos] += v_diffcoeffs[1]  * 
                        ( new_elevation_field[ixn,iy] - new_elevation_field[ix,iy] ) /  
                        2 / new_sediment_thickness[ix,iy]
                    i_neighbor_pos = list_pos[ixn,iy]
                    A[i_pos,i_neighbor_pos] = - v_diffcoeffs[1] * 
                        ( new_elevation_field[ixn,iy] - new_elevation_field[ix,iy] ) /  
                        2 / new_sediment_thickness[ixn,iy]
                end
            else # iy < ny
                if subaereal_blob[ix,iy+1] == 1 
                    A[i_pos,i_pos] -= v_diffcoeffs[1] * 
                        ( new_elevation_field[ix,iy+1] - new_elevation_field[ix,iy] ) /  
                        2 / new_sediment_thickness[ix,iy]
                    i_neighbor_pos = list_pos[ix,iy+1]
                    A[i_pos,i_neighbor_pos] -= v_diffcoeffs[1] * 
                        ( new_elevation_field[ix,iy+1] - new_elevation_field[ix,iy] ) /  
                        2 / new_sediment_thickness[ix,iy+1]
                elseif ocean_sink[ix,iy+1] == 1 
                    A[i_pos,i_pos] -= v_diffcoeffs[1] * 
                        ( 0. - new_elevation_field[ix,iy] ) /  
                        2 / new_sediment_thickness[ix,iy]
                end
            end                
        end # create A

        lu_A = lu(A)

        for i_sedtype in 1:n_sediment_types 
            old_inventory_list = []
            R = Float64[]; 
            for i_pos = 1:list_pos_filled
                ix = list_x[i_pos]; iy = list_y[i_pos]
                old_inventory = world.sediment_thickness[ix,iy] * 
                    world.sediment_fractions[i_sedtype,ix,iy]
                push!(old_inventory_list,old_inventory)
                push!(R, old_inventory + 
                    incoming_fluxes[i_sedtype][ix,iy] * time_step )
            end
            new_inventory_list = lu_A \ R
            new_fraction_field = fill(0.,nx,ny)
            for i_pos = 1:list_pos_filled
                ix = list_x[i_pos]; iy = list_y[i_pos]
                new_fraction = new_inventory_list[i_pos] / 
                    new_sediment_thickness[ix,iy]
                new_fraction_field[ix,iy] = new_fraction
                fraction_deposition_rate = 
                    ( new_inventory_list[i_pos] - old_inventory_list[i_pos] ) / 
                    time_step
                set_frac_diag("land_sediment_fraction_deposition_rate",ix,iy,i_sedtype,fraction_deposition_rate)
            end
        end
    end 
end
function land_transport_elevation()
    water_line = max( get_sealevel(), crust_level())
    crust_elevation = world.surface_elevation .- water_line
    return crust_elevation
end
function set_land_runoff_fluxes( elevation_field )        
    reset_frac_diag("coastal_sediment_fraction_flux")
    subaereal_mask = generate_mask_field(world.surface_type, sedimented_land )
    ocean_sink = generate_mask_field(world.crust_type, ocean_crust )
    for ix in 1:nx
        for iy in 1:ny
            if subaereal_mask[ix,iy] > 0
                h_diffcoeffs = get_horizontal_diffcoeff(land_base_diffcoeff,iy)
                v_diffcoeffs = get_vertical_diffcoeff(land_base_diffcoeff,iy)
                if ocean_sink[ix_left(ix),iy] > 0.
                    for i_sedtype in 1:n_sediment_types
                        volume_flux = elevation_field[ix,iy] * h_diffcoeffs[2] * # meters / step
                            world.sediment_fractions[i_sedtype,ix,iy] * # clay type, the only type moving offshore
                            areabox[iy] / time_step  # m3 / Myr
                        accum_frac_diag( "coastal_sediment_fraction_flux",
                            ix_left(ix),iy,i_sedtype,
                            volume_flux / areabox[iy] )
                    end
                end
                if ocean_sink[ix_right(ix),iy] > 0.
                    for i_sedtype in 1:n_sediment_types
                        volume_flux = elevation_field[ix,iy] * h_diffcoeffs[1] * 
                            world.sediment_fractions[i_sedtype,ix,iy] *
                            areabox[iy] / time_step
                        accum_frac_diag( "coastal_sediment_fraction_flux",
                            ix_right(ix),iy,i_sedtype,
                            volume_flux / areabox[iy] )
                    end
                end
                if iy > 1
                    if ocean_sink[ix,iy-1] > 0  
                        for i_sedtype in 1:n_sediment_types
                            volume_flux = elevation_field[ix,iy] * v_diffcoeffs[2] * 
                                world.sediment_fractions[i_sedtype,ix,iy] *
                                areabox[iy] / time_step
                            accum_frac_diag( "coastal_sediment_fraction_flux",
                                ix,iy-1,i_sedtype,
                                volume_flux / areabox[iy-1] )
                        end
                    end
                end
                if iy < ny
                    if ocean_sink[ix,iy+1] > 0 
                        for i_sedtype in 1:n_sediment_types
                            volume_flux = elevation_field[ix,iy] * v_diffcoeffs[1] * 
                                world.sediment_fractions[i_sedtype,ix,iy] *
                                areabox[iy] / time_step
                            accum_frac_diag( "coastal_sediment_fraction_flux",
                                ix,iy+1,i_sedtype,
                                volume_flux / areabox[iy+1] )
                        end
                    end
                end
            end
        end
    end 
    return
end
function aolean_transport()
    # resets and fills aolean_deposition_rate and aolean_deposition_rate
    set_diag("aolean_clay_erosion_rate",fill(0.,nx,ny))
    set_diag("aolean_clay_deposition_rate",fill(0.,nx,ny))
    for ix in 1:nx
        for iy in 1:ny
            if world.freeboard[ix,iy] > 0. &&
                world.surface_type[ix,iy] != exposed_basement
                rate = world.freeboard[ix,iy] *
                    aolean_erosion_rate_constant * # meters / Myr
                    rho_continent_crust / rho_sediment # m sediment / Myr
                set_diag("aolean_clay_erosion_rate",ix,iy,rate) # m / Myr
            end
        end
    end
    total_aolean_erosion = volumefield_total(get_diag("aolean_clay_erosion_rate")) 
    mean_aolean_deposition_rate = total_aolean_erosion / ocean_area() # meters / Myr
    # m3 / Myr
    for ix in 1:nx
        for iy in 1:ny
            if world.freeboard[ix,iy] < 0.
                set_diag("aolean_clay_deposition_rate",ix,iy,
                    mean_aolean_deposition_rate ) 
            end
        end
    end
    return
end

function update_surface_types()
    n_buried = 0
    for ix in 1:nx 
        for iy in 1:ny
            if world.surface_type[ix,iy] == exposed_basement # incoming from plate regrid
                neighbors = get_border_neighbor_coords(ix,iy)
                for nc in neighbors
                    ixn = nc[1]; iyn = nc[2]
                    if world.surface_type[ixn,iyn] == sedimented_land &&
                        world.freeboard[ixn,iyn] > world.freeboard[ix,iy] 
                        set_diag("surface_type_changes",ix,iy,burying_land)
                        n_buried += 1
                        world.surface_type[ix,iy] = sedimented_land
                    end
                end
            end
        end
    end
    if n_buried > 0
        println("  burying ",n_buried)
    end
    return n_buried
end
function setup_surface_types() 
    # build new surface type grid as initial condition, assumes soil cover everywhere
    isostacy()
    original_elevation = land_transport_elevation()
    is_land = greater_than_mask_field( original_elevation, 0. )
    depo_seafloor = greater_than_mask_field( - original_elevation, 100. )
    ocean_nondep_shelf = 1. .- is_land .- depo_seafloor
    recast_surface_type = is_land .* sedimented_land .+ 
        depo_seafloor .* pelagic_seafloor .+
        ocean_nondep_shelf .* ocean_shelf
    world.surface_type = recast_surface_type
end

function land_cell_balance(ix,iy)
    println(diffusive_flux_cell(world.freeboard,ix,iy,"land","up"))
    println(diffusive_flux_cell(world.freeboard,ix,iy,"land","left")," ",
        diffusive_flux_cell(world.freeboard,ix,iy,"land","left"))
    println(diffusive_flux_cell(world.freeboard,ix,iy,"land","down"))
end
function diffusive_flux_cell(elevation,ix,iy,domain,direction)
    if domain == "land"
        if direction == "left"
            diffcoeff = get_horizontal_diffcoeff(land_base_diffcoeff,iy)[2]
            neighbor_elevation = elevation[ix-1,iy]
        elseif direction == "right"
            diffcoeff = get_horizontal_diffcoeff(land_base_diffcoeff,iy)[1]
            neighbor_elevation = elevation[ix+1,iy]
        elseif direction == "up"
            diffcoeff = get_vertical_diffcoeff(land_base_diffcoeff,iy)[1]
            neighbor_elevation = elevation[ix,iy+1]
        else  # presumably down
            diffcoeff = get_vertical_diffcoeff(land_base_diffcoeff,iy)[2]
            neighbor_elevation = elevation[ix,iy-1]
        end
    else # presumably seafloor
        if direction == "left"
            diffcoeff = get_horizontal_diffcoeff(land_base_diffcoeff,iy)[2]
            neighbor_elevation = elevation[ix-1,iy]
        elseif direction == "right"
            diffcoeff = get_horizontal_diffcoeff(land_base_diffcoeff,iy)[1]
            neighbor_elevation = elevation[ix+1,iy]
        elseif direction == "up"
            diffcoeff = get_vertical_diffcoeff(land_base_diffcoeff,iy)[1]
            neighbor_elevation = elevation[ix,iy+1]
        else  # presumably down
            diffcoeff = get_vertical_diffcoeff(land_base_diffcoeff,iy)[2]
            neighbor_elevation = elevation[ix,iy-1]
        end
    end
    d_elevation = elevation[ix,iy] - neighbor_elevation
    dhdt = diffcoeff * d_elevation / time_step
    flux = dhdt / sediment_freeboard_expression
    return flux
end

         

function get_vertical_diffcoeff(base_diffcoeff,iy)
    diffusive_time_step = time_step * 1.E6 # years
    dx_gradient = delta_y; dx_selfwidth = delta_x[iy]
    if iy == ny
        diffup = 0.
    else
        dx_interface = ( delta_x[iy] + delta_x[iy+1] ) / 2
        diffup = base_diffcoeff * # m2 / yr * delta_elevation would give m3 / yr
            diffusive_time_step  /    # m3 
            dx_gradient *         # to get dElev/dx
            dx_interface /        # m3 across the interface
            dx_selfwidth / dx_gradient  # m elevation change in local box
    end
    if iy == 1
        diffdown = 0.
    else
        dx_interface = ( delta_x[iy] + delta_x[iy-1] ) / 2
        diffdown = base_diffcoeff * # m2 / yr * delta_elevation would give m3 / yr
            diffusive_time_step  /    # m3 
            dx_gradient *         # to get dElev/dx
            dx_interface /        # m3 across the interface
            dx_selfwidth / dx_gradient  # m elevation change in local box
    end
    return [ diffup, diffdown ]
end
function get_horizontal_diffcoeff(base_diffcoeff,iy)
    dx_width = delta_y; dx_interface = delta_y
    dx_gradient = max( delta_x[iy], delta_y / 2. )
    diffusive_time_step = time_step * 1.E6
    diffcoeff = base_diffcoeff * # m2 / yr * delta_elevation would give m3 / yr
        diffusive_time_step  /    # m3 
        dx_gradient *         # to get dElev/dx
        dx_interface /        # m3 across the interface
        dx_width / dx_gradient  # m elevation change in local box
    return [ diffcoeff, diffcoeff ]
end

    

   
