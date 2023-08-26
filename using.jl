# This is the file to use necessary packages.
# Last edited by Chaos on 2023/08/26.

using Distributions # to use functions: "truncated", "Normal"
using LinearAlgebra # to use function "norm"
using Plots # to use functions: "plot3d", "scatter", "animate", "gif"
gr() # switch to the "GR" backend
using Random # to use function "rand"
using StaticArrays # for simulation speed, to use: "SVector", "MVector"
using Statistics # to use function "elapsed"