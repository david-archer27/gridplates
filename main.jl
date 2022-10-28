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
include("caco3.jl")
include("isostacy.jl")
include("utilities.jl")
include("io.jl")

cd( base_directory * "/" * code_base_directory )
rotations = read_rotation_file("1000_0_rotfile_Merdith_et_al.rot")
norotation = rotation2matrix(rotations[1])
orogenic_events = create_orogenies()

setup_working_directories()
filename = base_directory * "/" * code_base_directory * "/" * "plate_ID_changes.bson"
BSON.@load filename plateID_change_log 

create_everything( earliesttime ) 
#
#world = read_world(440.)
#plates = read_plates()
#step_tectonics()
#
step_everything()
##

#@time 

# plateID_change_log = Dict()
#eyeball_all_changing_plateIDs( )
##
if enable_eyeball_changing_plateIDs
    filename = base_directory * "/" * code_base_directory * "/" * 
        "plate_ID_changes_by_" * string( main_time_step ) * ".bson"
    BSON.@load filename change_log
    global plateID_change_log = change_log
end
##
run_timeseries()
# + main_time_step )
#world = read_world(501)
#plates = read_plates()
#run_timeseries()
#
#step_everything( )
#world.sealevel += 0.02
#step_everything( )
#step_geomorph( )
#save_world(); save_plates()
#
#step_everything( true )
#

run_timeseries()

##

animate_all(  )

#

create_timeseries_charts( )
