# Head file of the N-body simulator. To use necessary packages, set constants and global variables and define necessary struct(s).
# Last edited by Cha0s_MnK on 2023/09/25.

#using Distributions # to use functions: "truncated", "Normal"
using LinearAlgebra # to use function "norm"
using Plots # to use functions: "plot3d", "scatter", "animate", "gif"
gr() # switch to the "GR" backend
using Random # to use function "rand"
using StaticArrays # for simulation speed, to use: "SVector", "MVector"
using Statistics # to use functions: "mean", "elapsed"

# constants
const day = 24.0 * 60.0 * 60.0 # (s)
const G = 6.67430e-11 # gravitational constant (m^3/kg/s^2)
const Julian_year = 365.25 * 24.0 * 60.0 * 60.0 # (s)
const pc = 30856775814913773 # parsec (m)
const solar_mass = 1.98847e30 # (kg)
const km_per_s = 1.0e3 # (m/s)

# global variables
global Mean_mass = 1.0 # mass of the bodies (solar mass)
global mean_mass = Mean_mass * solar_mass
global num_bodies = 2 # number of the bodies
global num_name = 0.1 # name number of the saved files
#global standard_deviation_mass = 0.0 # standard deviation of the masses of the bodies (kg); temporarily fixed
global theta = 1.0 # opening angle (rad)
global TotalLength = 1 # side length of the cube (pc)
global total_length = TotalLength * pc
global TotalTime = 2.0e8 # total simulation time (Julian year)
global total_time = TotalTime * Julian_year
global Timestep = 8.0e3 # time step (Julian year)
global timestep = Timestep * Julian_year
global SofteningLength = 0.0 # softening length (pc)
global softening_length = SofteningLength * pc
global velocity_dispersion = (10.0 * km_per_s)^2

# struct(s)
mutable struct Body
    mass::Float64
    position::MVector{3, Float64}
    velocity::MVector{3, Float64}
    acceleration::MVector{3, Float64}
end