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

rotations = read_rotation_file("1000_0_rotfile_Merdith_et_al.rot")
norotation = rotation2matrix(rotations[1])
orogenic_events = create_orogenies()
create_everything( earliesttime ) #+ time_step )
setup_working_directories( )
BSON.@load "plate_ID_changes.bson" plateID_change_log


