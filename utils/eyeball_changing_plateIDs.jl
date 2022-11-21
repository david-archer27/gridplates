function eyeball_all_changing_plateIDs( )
    oldIDmap = read_plateIDs( world.age + main_time_step )
    newIDmap = read_plateIDs( world.age )
    while world.age <= 550
        println("starting ", world.age )
        changepairs = eyeball_changing_plateIDs(oldIDmap,newIDmap)
        global plateID_change_log[ world.age ] = changepairs
        #world.age -= main_time_step
        oldIDmap = newIDmap
        newIDmap = read_plateIDs( world.age )
        filename = base_directory * "/" * code_base_directory * "/" * "plate_ID_changes.bson"
        rm(filename, force=true)
        println("saving ", filename)
        BSON.@save filename plateID_change_log
    end
    
    #return change_log
end

function eyeball_changing_plateIDs( age )
    global world = create_world( age )
    changepairs = eyeball_changing_plateIDs()
    plateID_change_log[world.age] = changepairs
    return changepairs
end
function eyeball_changing_plateIDs()
    oldIDmap = read_plateIDs( world.age + main_time_step )
    newIDmap = read_plateIDs( world.age )
    changepairs = eyeball_changing_plateIDs(oldIDmap,newIDmap)
    plateID_change_log[world.age] = changepairs
    return changepairs
end
function eyeball_changing_plateIDs(oldIDmap,newIDmap)
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
    changepairs = unique(changepairs)
    filtered_changepairs = []
    for changepair in changepairs
        oldplateID = changepair[1]; newplateID = changepair[2]
        oldIDmask = eq_mask(oldIDmap,oldplateID)
        newIDmask = eq_mask(newIDmap,newplateID)
        needschangingworldmask = oldIDmask .* newIDmask # fill(0,nx,ny)
        nskins = length( get_blob_onion_layers(needschangingworldmask ) )
        if nskins > 1
            push!(filtered_changepairs,[oldplateID,newplateID])
        end
    end
    changepairs = filtered_changepairs
    num_changepairs = length(filtered_changepairs)
    println("number of open possibilities = ", num_changepairs )
    needschangingworldmasks = fill(0,nx,ny,num_changepairs)

    filtered_changepairs = []
    old_plateID_list = world.plateIDlist
    old_plateID_field = world.plateID
    old_continent_mask = eq_mask(world.crust_type,continent_crust)
    global world = create_world( world.age + main_time_step )
    old_plateID_field = world.plateID
    new_plateID_list = world.plateIDlist
    new_continent_mask = eq_mask(world.crust_type,continent_crust)
    global world = create_world( world.age - main_time_step )
    i_choose_pair = 1
    while i_choose_pair <= num_changepairs
        choosechangepair = changepairs[i_choose_pair]
        if length(filtered_changepairs) > 0
            if filtered_changepairs[end] == choosechangepair # for backing up
                pop!(filtered_changepairs)
            end
        end
        oldplateID = choosechangepair[1]; newplateID = choosechangepair[2]
        viz_changes = fill(3,nx,ny) # no changes = green
        for i_pair_viz in 1:num_changepairs
            vizchangepair = changepairs[i_pair_viz]
            oldplateID = vizchangepair[1]; newplateID = vizchangepair[2]
            oldIDmask = eq_mask(oldIDmap,oldplateID)
            newIDmask = eq_mask(newIDmap,newplateID)
            needschangingworldmasks[:,:,i_pair_viz] = oldIDmask .* newIDmask # fill(0,nx,ny)
            if i_pair_viz > i_choose_pair # undetermined yet
                viz_changes -= needschangingworldmasks[:,:,i_pair_viz] # light violet
            elseif i_pair_viz == i_choose_pair
                viz_changes += 2 * needschangingworldmasks[:,:,i_pair_viz]
            elseif vizchangepair in filtered_changepairs
                viz_changes += needschangingworldmasks[:,:,i_pair_viz]
            else
                viz_changes -= 2 * needschangingworldmasks[:,:,i_pair_viz]
            end
        end        #oldIDmask = eq_mask(oldIDmap,oldplateID)
        #newIDmask = eq_mask(newIDmap,newplateID)
        #needschangingworldmask = oldIDmask .* newIDmask # fill(0,nx,ny)
        scene = plot_field( viz_changes .+ 4 * needschangingworldmasks[:,:,i_choose_pair],0.,6. )
        Makie.contour!(scene,xcoords,ycoords,needschangingworldmasks[:,:,i_choose_pair],color=:red)
        for plateID in old_plateID_list
            maskfield = eq_mask(old_plateID_field,plateID)
            Makie.contour!(scene,xcoords,ycoords,maskfield,color=:black)
        end
        Makie.contour!(scene,xcoords,ycoords,old_continent_mask,color=:black)
        for plateID in new_plateID_list
            maskfield = eq_mask(world.plateID,plateID)
            Makie.contour!(scene,xcoords,ycoords,maskfield,color=:white)
        end
        Makie.contour!(scene,xcoords,ycoords,new_continent_mask,color=:white)
        display(scene)
        
        println(i_choose_pair, " ", oldplateID," ",newplateID," y/n?")
        input_char = readline()
        if input_char == "y"
            println("got it")
            i_choose_pair += 1
            push!(filtered_changepairs,choosechangepair)
        elseif input_char == "b" # back
            i_choose_pair -= 1
        elseif input_char == "bb" # back
            i_choose_pair = 1 
        elseif input_char == "d" # done
            i_choose_pair = num_changepairs + 1
        else
            i_choose_pair += 1
        end
    end
    return filtered_changepairs
end
function rejigger_plateID_changes()
    rejiggered_plateID_changes = Dict()
    for i in 1:550
        changelist = []
        if haskey(plateID_change_log,i)
            nchanges = length(plateID_change_log[i])
            if nchanges > 0 
                for j in 1:nchanges
                    push!(changelist,[plateID_change_log[i][j][2],plateID_change_log[i][j][1]])
                end 
            end
        end
        rejiggered_plateID_changes[i] = changelist
    end
    return rejiggered_plateID_changes
end
function plot_plateID_changes(  )
    #global world = create_world( age )
    age = world.age
    new_plateID_list = world.plateIDlist
    new_plateID_field = world.plateID
    new_continent_mask = eq_mask(world.crust_type,continent_crust)
    changepairs = []
    if haskey(plateID_change_log, age)
        changepairs = plateID_change_log[ age ]
    end
    global world = create_world( world.age + main_time_step )
    old_plateID_field = world.plateID
    old_plateID_list = world.plateIDlist
    #update_world_continents_from_file() 
    old_continent_mask = eq_mask(world.crust_type,continent_crust)
    global world = create_world( world.age - main_time_step )
    viz_field = fill(0.,nx,ny)
    for ixworld in 1:nx # 1+fp:nx-fp
        for iyworld in 1:ny # 1+fp:ny-fp
            if old_plateID_field[ixworld,iyworld] != new_plateID_field[ixworld,iyworld]
                viz_field[ixworld,iyworld] = 1
            end
        end
    end
    for changepair in changepairs
        oldplateID = changepair[1]; newplateID = changepair[2]
        oldIDmask = eq_mask(old_plateID_field,oldplateID)
        newIDmask = eq_mask(new_plateID_field,newplateID)
        viz_field[:,:] += 2 .* (oldIDmask .* newIDmask )# fill(0,nx,ny)
    end
    scene = plot_field( viz_field, 0., 2. )
    for plateID in old_plateID_list
        maskfield = eq_mask(old_plateID_field,plateID)
        Makie.contour!(scene,xcoords,ycoords,maskfield,color=:white)
    end
    Makie.contour!(scene,xcoords,ycoords,old_continent_mask,color=:white)
    for plateID in new_plateID_list
        maskfield = eq_mask(new_plateID_field,plateID)
        Makie.contour!(scene,xcoords,ycoords,maskfield,color=:black)
    end
    Makie.contour!(scene,xcoords,ycoords,new_continent_mask,color=:black)
    text!(scene,string(world.age) * " Mya",position=(-180,100),textsize=20)
    #text!(scene,string(plateID_change_log[age]),position=(-180,-10),textsize=10 )
    return scene
end
    

function animate_plateID_changes()
    directory = base_directory * "/" * output_directory * output_tag * "/" * 
        animation_directory * "/" * "plateID_changes" 
    cd( directory )
    starting_file_list = readdir()
    image_number = 0
    ages = [ 0 ]
    for age in 549:-1:0 # animation_initial_age: - main_time_step * animation_n_step : 0
        push!(ages,age)
    end
    for age in ages
        image_number += 1
        image_file_name = "img." * lpad(image_number,3,"0") * ".png"
        if image_file_name in starting_file_list
            println( "already done ", age," ", image_file_name )
        else
            println( "creating ", age, " ", directory * "/" * image_file_name )
            global world = create_world( age )
            scene = plot_plateID_changes( )
            
            Makie.save( image_file_name, scene )
        end
    end
    mp4_file =  "plateID_changes.mp4"
    println("compiling ", mp4_file)
    rm( mp4_file,force=true)
    run(`ffmpeg -r 5 -f image2 -s 1920x1080 -i img.%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p $mp4_file`)
    cd( base_directory * "/" * code_base_directory )
end
