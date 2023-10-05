# This is the head file of the collisional N-body simulator. 
# It can use necessary packages, set constants, set global variables and define necessary struct(s).
# Last edited by Cha0s_MnK on 2023/10/05.

#using Distributions # to use functions: "truncated", "Normal"
using LinearAlgebra # to use function "norm", "dot"
using Plots # to use functions: "plot3d", "scatter", "animate", "gif"
gr() # switch to the "GR" backend
using Random # to use function "rand"
using StaticArrays # for simulation speed, to use: "SVector", "MVector"
using Statistics # to use functions: "mean", "elapsed"

# constants
const day = 24.0 * 60.0 * 60.0 # (s)
const G = 6.67430e-11 # gravitational constant (m^3/kg/s^2)
const Julian_year = 365.25 * 24.0 * 60.0 * 60.0 # Julian year (s)
const pc = 30856775814913773 # parsec (m)
const solar_mass = 1.98847e30 # solar mass (kg)
const KMperS = 1.0e3 # kilometre per second (m/s)

# global variables
global Mean_mass = 32.0 # mass of the bodies (solar mass)
global mean_mass = Mean_mass * solar_mass # mass of the bodies (kg)
global num_bodies = 2 # number of the bodies
global num_name = 2.1 # name number of the saved files
global num_steps = Int(8e4) # number of simulation steps
#global standard_deviation_mass = 0.0 # standard deviation of the masses of the bodies (kg)
global theta = 1.0 # opening angle (rad)
global TotalLength = 1 # side length of the simulation box (pc) ; temporarily fixed
global total_length = TotalLength * pc # side length of the simulation box (m); temporarily fixed
global TotalTime = 0.0 # total simulation time (Julian year)
global MaxTimestep = 1.0e4 # maximum timestep (Julian year)
global Timestep = 0.0 # current timestep (Julian year)
global timestep = 0.0 # current timestep (s)
global SofteningLength = 0.0 # softening length (pc)
global softening_length = SofteningLength * pc # softening length (m)

# struct(s)
mutable struct Body
    mass::Float64
    position::MVector{3, Float64}
    velocity::MVector{3, Float64}
    acceleration::MVector{3, Float64}
    jerk::MVector{3, Float64}
end