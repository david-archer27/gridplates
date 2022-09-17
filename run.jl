
using SparseArrays
using LinearAlgebra
using Rotations, SharedArrays, StaticArrays
using Printf
using GLMakie
#inline!(true)
using Colors
using GeometryTypes
using Contour
using BSON
using StaticArraysCore
using Plots
#
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
#
cd( base_directory * code_base_subdirectory )
rotations = read_rotation_file("1000_0_rotfile_Merdith_et_al.rot")
norotation = rotation2matrix(rotations[1])

orogenic_events = create_orogenies()

create_everything( 550 ) #earliesttime + time_step )
#step_everything()
#world = read_world(150); read_plates()
#
run_timeseries()
