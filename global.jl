# This is the file to set constants and global variables
# Last edited by Chaos on 2023/08/26.

# constant(s)
const G = 6.67430e-11 # gravitational constant (m^3/kg/s^2)

# global variables
global mean_mass = 3000.0 # average mass of the bodies (kg)
global num_bodies = 32 # number of random bodies
global num_name = 3 # name number of the saved files
global standard_deviation_mass = 0.0 # standard deviation of the masses of the bodies (kg)
global theta = 0.1 # opening angle (rad)
global total_length = 100.0 # side length of the cube (m)
global total_time = 64*60.0*60.0 # total simulation time (s)
global time_step = 1.0 # time step (s)
global softening_coefficient = 0.001 # soft length = softening coefficient * total length