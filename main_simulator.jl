# The main running file of the N-body simulator.
# Last edited by Cha0s_MnK on 2023/09/25.

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
    # generate random coordinates (uniform distribution)
    x0 = rand() # generate a random floating-point number between 0 and 1
    y0 = rand()
    z0 = rand()
    x = x0 * total_length .- total_length / 2 # switch the random number from between 0 and 1 to between -total_length/2 and total_length/2
    y = y0 * total_length .- total_length / 2
    z = z0 * total_length .- total_length / 2
    # generate random velocities with certain velocity dispersion
    sigma = sqrt(velocity_dispersion)
    v_x = randn() * sigma
    v_y = randn() * sigma
    v_z = randn() * sigma
    return Body(mass, MVector{3, Float64}(x, y, z), MVector{3, Float64}(v_x, v_y, v_z), MVector{3, Float64}(0.0, 0.0, 0.0)) # initialize acceleration to zero
end

function generateBodies()::Vector{Body}
    # generate random bodies

    # generate N random bodies with zero acceleration stored
    bodies = [generateBody() for _ in 1:num_bodies]
    # calculate accelerations
    direct_summation!(bodies)
    return bodies
end

function generateKeplerBodies()::Vector{Body}
    # generate 2 bodies for Kepler problem

    # set initial positions and velocities
    bodies = [Body(mean_mass, MVector{3, Float64}(0.3*total_length, 0.3*total_length, 0.3*total_length), MVector{3, Float64}(-0.2*KMperS, -0.2*KMperS, 0.0), MVector{3, Float64}(0.0, 0.0, 0.0)),
                Body(mean_mass, MVector{3, Float64}(-0.3*total_length, -0.3*total_length, -0.3*total_length), MVector{3, Float64}(0.2*KMperS, 0.2*KMperS, 0.0), MVector{3, Float64}(0.0, 0.0, 0.0))]
    direct_summation!(bodies) # calculate accelerations
    return bodies
end

function generateKeplerBodies2()::Vector{Body}
    # generate 2 bodies for Kepler problem

    # set initial positions and velocities
    bodies = [Body(1000.0*mean_mass, MVector{3, Float64}(0.0, 0.0, 0.0), MVector{3, Float64}(0.0, 0.0, 0.0), MVector{3, Float64}(0.0, 0.0, 0.0), MVector{3, Float64}(0.0, 0.0, 0.0)),
                Body(mean_mass, MVector{3, Float64}(-0.3*total_length, -0.3*total_length, -0.3*total_length), MVector{3, Float64}(8.0*KMperS, 8.0*KMperS, 0.0), MVector{3, Float64}(0.0, 0.0, 0.0), MVector{3, Float64}(0.0, 0.0, 0.0))]
    direct_summation!(bodies) # calculate accelerations
    return bodies
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

    timestep_criterion(bodies)
    time_integration!(bodies, force_calculation!)
    periodic_boundary!(bodies) # apply periodic boundary condition
    #bound_pairs(bodies) # criterion to avoid boud pairs
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
        # use 2nd order Runge-Kutta and direct summation to simulate
        simulate_step!(bodies1, RungeKutta2nd!, direct_summation!)
        # use explicit Euler and direct summation to simulate
        simulate_step!(bodies2, Leapfrog!, direct_summation!)
        # calculate the distances; fix them; store them
        get_distances(i, distances, bodies1, bodies2)
    end
    # plot curves for the time evolution of each body; save the plot as a PNG file
    plotCurves(distances, "Results/TestCurves$num_name.png")
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

function main_Animation(num_steps::Int)
    # main function that simulate N-body problem using animation

    # generate random bodies
    bodies = [generateBody() for _ in 1:num_bodies]
    # create a matrix to store forces
    forces = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    # calculate gravitational forces
    direct_summation!(forces, bodies)
    # simulate N-body problem; create an animation
    animation = @animate for i in 1:num_steps
        # plot current bodies; use 4th order Runge-Kutta and BHoctree to simulate
        canvas = plotBodies(bodies,
        "4th order Runge-Kutta \n\n Number of bodies = $num_bodies \n Mass = $Mean_mass solar mass \n Side length = $TotalLength pc \n Total time = $TotalTime Julian year \n Time step = $TimeStep Julian year")
        plotEdges(bodies)
        simulate_step!(bodies, RungeKutta4th!, BHoctree!)
        # plot the canvas in a 1*1 layout
        plot(canvas, layout = (1, 1))
    end
    # save the animation as a MP4 file
    gif(animation, "Results/Scholarship$num_name.mp4", fps = 24)
end

function main_Kepler(num_steps::Int)
    # main function that simulate Kepler problem using animation

    # generate random bodies
    bodies = generateKeplerBodies2()
    # create a beautiful 3D canvas
    plot_length = 0.55*total_length
    canvas = plot3d(
        size = (1080, 720),
        legend = false, legendfontfamily = "Times New Roman", legendfontsize = 10,
        xlims = (-plot_length, plot_length), ylims = (-plot_length, plot_length), zlims = (-plot_length, plot_length),
        title = "Kepler problem (Leapfrog, direct summation) \n\n Mass = $Mean_mass solar mass \n Side length = $TotalLength pc \n Number of simulation steps = $num_steps", titlefontfamily = "Times New Roman", titlefontsize = 14,
        xticks = true, yticks = true, zticks = true, tickfontfamily = "Times New Roman", tickfontsize = 10,
        xlabel = "X", xlabelfontfamily = "Times New Roman", xlabelfontsize = 12, 
        ylabel = "Y", ylabelfontfamily = "Times New Roman", ylabelfontsize = 12,
        zlabel = "Z", zlabelfontfamily = "Times New Roman", zlabelfontsize = 12,
        background_color = RGB(0.95, 0.95, 0.95), 
        fg_legend = :white,
        linealpha = 0.7, linewidth = 1, linestyle = :auto, 
        marker = (:circle, 4),
        grid = true,
        color = :viridis,
        camera = (30, 30), elevation = 30, azimuth = 30
    )
    # create an animation
    animation = @animate for _ in 1:num_steps
        # plot current bodies
        scatter!(canvas, [bodies[1].position[1]], [bodies[1].position[2]], [bodies[1].position[3]], markersize=1, markerstrokewidth=0, color=:blue)
        scatter!(canvas, [bodies[2].position[1]], [bodies[2].position[2]], [bodies[2].position[3]], markersize=1, markerstrokewidth=0, color=:red)
        # use Leapfrog and direct summation to simulate
        #simulate_step!(bodies, RungeKutta4th!, direct_summation!)
        simulate_step!(bodies, Leapfrog!, direct_summation!)
    end
    # save the animation as a MP4 file
    gif(animation, "Results/Kepler$num_name.mp4", fps = 24)
end

# run the selected main function and record the running time
running_time = @elapsed main_Kepler(num_steps) # (s)
println("Running time: $running_time seconds")

# calculate the relaxation time
relaxation_time = num_bodies / (8 * log(num_bodies)) * total_length * sqrt(total_length / (G * num_bodies * mean_mass)) / Julian_year
# check whether the relaxation time is much larger than the time step
println("Relaxation time = $relaxation_time Julian_year")
println("Time step = $Timestep Julian_year")