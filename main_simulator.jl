# The main running file of the N-body simulator.
# Last edited by Cha0s_MnK on 2023/09/06.

# "BHoctree" is short for "Barnes & Hut oct-tree"
# "OpeningAngle" is short for the "opening angle" criterion
# original paper:
# A hierarchical O(N log N) force-calculation algorithm
# https://www.nature.com/articles/324446a0

# use necessary packages as long as set constants and global variables.
include("head.jl")
# use BHoctree
include("BHoctree.jl")
# use time integration techniques
include("time_integration.jl")
# use gravitational force calculation methods
include("force_calculation.jl")
# use self-defined plot functions
include("plot.jl")

# basic functions of the N-body simulator
function generateBody()::Body
    # generate 1 random body

    # generate mass
    mass = mean_mass
    # generate random coordinate (uniform distribution)
    x0 = rand() # generate a random floating-point number between 0 and 1
    y0 = rand()
    z0 = rand()
    x = x0 * total_length .- total_length / 2 # switch the random number from between 0 and 1 to between -total_length/2 and total_length/2
    y = y0 * total_length .- total_length / 2
    z = z0 * total_length .- total_length / 2
    return Body(mass, MVector{3, Float64}(x, y, z), MVector{3, Float64}(0.0, 0.0, 0.0)) # initialize velocity to zero
end

function periodic_boundary!(bodies::Vector{Body})
    # apply periodic boundary condition to modify coordinates

    for i in 1:num_bodies 
        for j in 1:3
            bodies[i].position[j] = mod(bodies[i].position[j] + total_length/2, total_length) - total_length/2
        end
    end
end

function simulate_step!(bodies::Vector{Body}, time_integration!::Function, force_calculation!::Function)
    # 1 simulation step using selected time integration technique and gravitational force calculation method

    time_integration!(bodies, force_calculation!)
    periodic_boundary!(bodies) # apply periodic boundary condition
    bound_pairs(bodies) # criterion to avoid boud pairs
end

# main functions
function main_Test(num_steps::Int)
    # main function that simulate N-body problem using animation

    # generate random bodies and deepcopy it
    bodies1 = [generateBody() for _ in 1:num_bodies]
    bodies2 = deepcopy(bodies1)
    # create a matrix to store distances between the same body evolving in 2 different ways
    distances = zeros(num_steps, num_bodies)
    # simulate N-body problem in 2 different ways; store the distances; create an animation
    animation = @animate for i in 1:num_steps
        # plot current bodies; use explicit Euler and direct summation to simulate
        canvas1 = plotBodies(bodies1,
        "2nd order Runge-Kutta \n\n Number of bodies = $num_bodies \n Mass = $Mean_mass solar mass \n Side length = $TotalLength pc")
        simulate_step!(bodies1, RungeKutta2nd!, direct_summation!)
        # plot current bodies; use 2nd order Runge-Kutta and direct summation to simulate
        canvas2 = plotBodies(bodies2,
        "Leapfrog \n\n Total time = $TotalTime Julian year \n Time step = $TimeStep Julian year")
        simulate_step!(bodies2, Leapfrog!, direct_summation!)
        # calculate the distances; fix them; store them
        get_distances(i, distances, bodies1, bodies2)
        # plot the canvases in a 1*2 layout
        plot(canvas1, canvas2, layout = (1, 2))
    end
    # save the animation as a MP4 file
    gif(animation, "Results/TestAnimation$num_name.mp4", fps = 24)
    # plot curves for the time evolution of each body; save the plot as a PNG file
    plotCurves(distances, "Results/TestCurves$num_name.png")
end



function main_Curves(num_steps::Int)
    # main function that compares direct summation and BHoctree using curves

    # generate random bodies and deepcopy it
    bodies1 = [generateBody() for _ in 1:num_bodies]
    bodies2 = deepcopy(bodies1)
    # create a matrix to store distances between the same body evolving in 2 different ways
    distances = zeros(num_steps, num_bodies)
    # simulate N-body problem in 2 different ways; store the distances
    for i in 1:num_steps
        # use direct summation to simulate
        direct_summation_step!(bodies1)
        # use BHoctree to simulate
        BHoctree_step!(bodies2)
        # calculate the distances; fix them; store them
        get_distances(i, distances, bodies1, bodies2)
    end
    # plot curves for the time evolution of each body; save the plot as a PNG file
    plotCurves(distances, "Results/Curves$num_name.png")
end

function main_AnimationCurves(num_steps::Int)
    # main function that compares direct summation and BHoctree using animation and curves

    # generate random bodies and deepcopy it
    bodies1 = [generateBody() for _ in 1:num_bodies]
    bodies2 = deepcopy(bodies1)
    # create a matrix to store distances between the same body evolving in 2 different ways
    distances = zeros(num_steps, num_bodies)
    # simulate N-body problem in 2 different ways; store the distances; create an animation
    animation = @animate for i in 1:num_steps
        # plot current bodies; use direct summation to simulate
        canvas1 = plotBodies!(bodies1, 
        "Direct summation \n\n Number of bodies = $num_bodies \n Mass = $Mean_mass solar mass \n Side length = $TotalLength pc")
        direct_summation_step!(bodies1)
        # plot current bodies; use BHoctree to simulate
        canvas2 = plotBodies!(bodies2,
        "BHoctree \n\n Total time = $TotalTime Julian year \n Time step = $TimeStep Julian year \n Opening angle = $theta")
        BHoctree_step!(bodies2)
        # calculate the distances; fix them; store them
        get_distances(i, distances, bodies1, bodies2)
        # plot the canvases in a 1*2 layout
        plot(canvas1, canvas2, layout = (1, 2))
    end
    # save the animation as a MP4 file
    gif(animation, "Results/Animation$num_name.mp4", fps = 60)
    # plot curves for the time evolution of each body; save the plot as a PNG file
    plotCurves(distances, "Results/AnimationCurves$num_name.png")
end

# calculate total number of simulation steps
num_steps = Int(floor(total_time/time_step))

# run the selected main function and record the running time
running_time = @elapsed main_Test(num_steps) # (s)
println("Running time: $running_time seconds")