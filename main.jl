using SparseArrays
using LinearAlgebra
using Rotations, SharedArrays, StaticArrays
using Printf
using GLMakie
using Colors
using GeometryTypes
using Contour
using BSON
using StaticArraysCore
using Plots

include("params.jl")
include("steps.jl")
include("plates.jl")
include("orogeny.jl")
include("land_sed.jl")
include("ocean_sed.jl")
include("chem.jl")
include("isostacy.jl")
include("utilities.jl")
include("io.jl")

rotations = read_rotation_file( code_base_directory * "/drivers/1000_0_rotfile_Merdith_et_al.rot")
norotation = rotation2matrix(rotations[1])
orogenic_events = create_orogenies()
create_everything( earliesttime ) #+ time_step )
setup_working_directories( )
BSON.@load code_base_directory * "/drivers/plate_ID_changes.bson" plateID_change_log

function run_and_plot_simulation()
    run_timeseries()
    create_timeseries_charts()
    animate_all()
    create_html_directory( )
end
function run_and_plot_orogeny_simulation()
    global enable_geomorph = false
    run_timeseries()
    animate_scotese_elevation()
    animate_crust_thickening_rate( )
    create_timeseries_charts()
end
function run_timeseries()
    global log_IO = open( output_directory * "/logfile." * output_tag * ".txt", "w" )
    save_world()
    while world.age > 0
        verbose = true
        if floor(world.age/100) == world.age/100 || world.age == main_time_step
            verbose = true
        end
        step_everything( )
        flush( log_IO )
        save_world()
        time_interval = 250 # main_time_step * 10
        if floor(world.age/time_interval) == world.age/time_interval && 
            world.age >= time_interval
            save_plates()
        end
        flush(log_IO)
    end
    close(log_IO)
end
