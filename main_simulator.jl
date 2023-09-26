# This is the main running file of the collisional N-body simulator.
# Last edited by Cha0s_MnK on 2023/09/26.

include("head.jl") # use necessary packages as long as set constants and global variables.
include("BHoctree.jl") # use BHoctree
include("time_integration.jl") # use time integration techniques
include("force_calculation.jl") # use gravitational force calculation methods
include("plot.jl") # use self-defined plot functions

# basic functions of the N-body simulator
function generateBody()::Body # generate 1 random body
    # generate mass
    mass = mean_mass
    # generate random coordinates (uniform distribution)
    x0 = rand() # generate a random floating-point number between 0 and 1
    y0 = rand()
    z0 = rand()
    x = x0 * total_length - total_length / 2 # switch the random number from between 0 and 1 to between -total_length/2 and total_length/2
    y = y0 * total_length - total_length / 2
    z = z0 * total_length - total_length / 2
    # initialize velocity, acceleration and jerk to zero
    return Body(mass, MVector{3, Float64}(x, y, z), MVector{3, Float64}(0.0, 0.0, 0.0), MVector{3, Float64}(0.0, 0.0, 0.0), MVector{3, Float64}(0.0, 0.0, 0.0))
end

function generateBodies()::Vector{Body} # generate N random bodies
    # generate N random bodies with initial velocity, acceleration and jerk as zero
    bodies = [generateBody() for _ in 1:num_bodies]
    # calculate accelerations
    direct_summation!(bodies)
    # calculate jerks
    direct_summation_jerk!(bodies)
    return bodies
end

function generateKeplerBodies()::Vector{Body} # generate 2 bodies for Kepler problem
    # set initial positions and velocities
    bodies = [Body(1000.0*mean_mass, MVector{3, Float64}(0.0, 0.0, 0.0), MVector{3, Float64}(0.0, 0.0, 0.0), MVector{3, Float64}(0.0, 0.0, 0.0), MVector{3, Float64}(0.0, 0.0, 0.0)),
                Body(mean_mass, MVector{3, Float64}(-0.3*total_length, -0.3*total_length, -0.3*total_length), MVector{3, Float64}(8.0*KMperS, 8.0*KMperS, 0.0), MVector{3, Float64}(0.0, 0.0, 0.0), MVector{3, Float64}(0.0, 0.0, 0.0))]
    # calculate accelerations
    direct_summation!(bodies)
    return bodies
end

function periodic_boundary!(bodies::Vector{Body}) # apply periodic boundary condition to modify coordinates
    for i in 1:num_bodies 
        for j in 1:3
            bodies[i].position[j] = mod(bodies[i].position[j] + total_length/2, total_length) - total_length/2
        end
    end
end

function simulate_step!(bodies::Vector{Body}, time_integration!::Function, force_calculation!::Function) # 1 simulation step using selected time integration technique and gravitational force calculation method
    timestep_criterion(bodies)
    time_integration!(bodies, force_calculation!) # main simulation step
    periodic_boundary!(bodies) # apply periodic boundary condition to modify coordinates
    #bound_pairs(bodies) # criterion to avoid boud pairs
end

# main functions
function main_Kepler(num_steps::Int) # main function that simulate Kepler problem
    # generate random bodies
    bodies = generateKeplerBodies()
    # create a beautiful 3D canvas
    plot_length = 0.55*total_length
    canvas = plot3d(
        size = (1080, 720),
        legend = false, legendfontfamily = "Times New Roman", legendfontsize = 10,
        xlims = (-plot_length, plot_length), ylims = (-plot_length, plot_length), zlims = (-plot_length, plot_length),
        title = "Kepler problem (Leapfrog, direct summation) \n\n Mass = $Mean_mass solar mass \n Side length = $TotalLength pc \n Number of simulation steps = $num_steps", 
        titlefontfamily = "Times New Roman", titlefontsize = 14,
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
        # plot current 2 bodies
        scatter!(canvas, [bodies[1].position[1]], [bodies[1].position[2]], [bodies[1].position[3]], markersize=1, markerstrokewidth=0, color=:blue)
        scatter!(canvas, [bodies[2].position[1]], [bodies[2].position[2]], [bodies[2].position[3]], markersize=1, markerstrokewidth=0, color=:red)
        # simulation step using selected time integration technique and gravitational force calculation method
        #simulate_step!(bodies, RungeKutta4th!, direct_summation!)
        simulate_step!(bodies, Leapfrog!, direct_summation!)
        #simulate_step!(bodies, modified4thHermite!, direct_summation!)
    end
    # save the animation as a MP4 file
    gif(animation, "Results/Kepler$num_name.mp4", fps = 24)
end

# run the selected main function and record the running time
running_time = @elapsed main_Kepler(num_steps) # (s)
println("Running time: $running_time seconds")

# calculate the relaxation time
relaxation_time = num_bodies / (8 * log(num_bodies)) * total_length * sqrt(total_length / (G * num_bodies * mean_mass)) / Julian_year
# check whether the relaxation time is much larger than the maximum timestep
println("Relaxation time = $relaxation_time Julian_year >> Maximum timestep = $Timestep Julian_year")

println("Total simulation time = $TotalTime Julian year") # print the total simulation time