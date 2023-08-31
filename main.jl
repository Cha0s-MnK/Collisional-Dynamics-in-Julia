# This is the main running file of the N-body simulator.
# Last edited by Chaos on 2023/08/30.

# original paper:
# A hierarchical O(N log N) force-calculation algorithm
# https://www.nature.com/articles/324446a0

# "BHoctree" is short for "Barnes & Hut oct-tree"
# "OpeningAngle" is short for the "opening angle" criterion

# use necessary packages
include("using.jl")
# use constants and global variables
include("global.jl")
# use the simple N-body simulator
include("simulator.jl")
# use BHoctree
include("BHoctree.jl")
# use self-defined plot functions
include("plot.jl")

# main functions
function main_Curves(num_steps::Int)
    # main function that compares direct summation and BHoctree using curves

    # generate random bodies and deepcopy it
    bodies1 = [generateBody() for _ in 1:num_bodies]
    bodies2 = deepcopy(bodies1)
    # create an array to stroe distances
    distances = zeros(num_steps, num_bodies)
    # simulate N-body problem in different ways; store the distances
    for i in 1:num_steps
        # use direct summation to simulate
        direct_summation_step!(bodies1)
        # use BHoctree to simulate
        BHoctree_step!(bodies2)
        # calculate the distances; fix them; store them
        for j in 1:num_bodies
            distances[i, j] = norm(bodies1[j].position - bodies2[j].position)
            if distances[i, j] > 0.5*total_length # effects caused by periodic boundary conditions
                temporary_position1 = deepcopy(bodies1[j].position)
                temporary_position2 = deepcopy(bodies2[j].position)
                for k in 1:3
                    if (temporary_position1[k] - temporary_position2[k]) > 0.5*total_length
                        temporary_position2[k] += total_length
                    elseif (temporary_position2[k] - temporary_position1[k]) > 0.5*total_length
                        temporary_position1[k] += total_length
                    end
                end
                distances[i, j] = norm(temporary_position1 - temporary_position2)
            end
        end
    end
    # create an empty plot for curves
    curves = plotCurves()
    # add curves for the time evolution of each body
    for j in 1:num_bodies
        # change the unit of the variables
        # Time: s --> Julian year; Distances: m --> pc
        plot!(curves, [time_step.*(1:num_steps)./Julian_year], distances[:, j]./pc, label="Body $j")
    end
    # save the plot as a PNG file
    savefig(curves, "Results/Curves$num_name.png")
end

function main_Animation(num_steps::Int)
    # main function that compares direct summation and BHoctree using animation

    # generate random bodies and deepcopy it
    bodies1 = [generateBody() for _ in 1:num_bodies]
    bodies2 = deepcopy(bodies1)
    # simulate N-body problem in different ways; create an animation
    animation = @animate for _ in 1:num_steps
        # plot current bodies; use direct summation to simulate
        canvas1 = plotBodies!(bodies1, 
        "Direct summation \n\n number of bodies = $num_bodies \n mass of 1 body = $mean_mass (kg) \n time step = $time_step (s) \n total time = $total_time (s)")
        direct_summation_step!(bodies1)
        # plot current bodies; use BHoctree to simulate
        canvas2 = plotBodies!(bodies2, 
        "BHoctree \n\n Opening angle = $theta")
        BHoctree_step!(bodies2)
        # plot the canvases in a 1*2 layout
        plot(canvas1, canvas2, layout = (1, 2))
    end
    # save the animation as a MP4 file
    gif(animation, "Results/Animation$num_name.mp4", fps = 60)
end

function main_AnimationCurves(num_steps::Int)
    # main function that compares direct summation and BHoctree using animation and curves

    # generate random bodies and deepcopy it
    bodies1 = [generateBody() for _ in 1:num_bodies]
    bodies2 = deepcopy(bodies1)
    # create an array to stroe distances
    distances = zeros(num_steps, num_bodies)
    # simulate N-body problem in different ways; store the distances; create an animation
    animation = @animate for i in 1:num_steps
        # plot current bodies; use direct summation to simulate
        canvas1 = plotBodies!(bodies1, 
        "Direct summation \n\n number of bodies = $num_bodies \n mass of 1 body = $mean_mass (kg) \n time step = $time_step (s) \n total time = $total_time (s)")
        direct_summation_step!(bodies1)
        # plot current bodies; use BHoctree to simulate
        canvas2 = plotBodies!(bodies2, 
        "BHoctree \n\n Opening angle = $theta")
        BHoctree_step!(bodies2)
        # calculate the distances; fix them; store them
        for j in 1:num_bodies
            distances[i, j] = norm(bodies1[j].position - bodies2[j].position)
            if distances[i, j] > 0.5*total_length # effects caused by periodic boundary conditions
                temporary_position1 = deepcopy(bodies1[j].position)
                temporary_position2 = deepcopy(bodies2[j].position)
                for k in 1:3
                    if (temporary_position1[k] - temporary_position2[k]) > 0.5*total_length
                        temporary_position2[k] += total_length
                    elseif (temporary_position2[k] - temporary_position1[k]) > 0.5*total_length
                        temporary_position1[k] += total_length
                    end
                end
                distances[i, j] = norm(temporary_position1 - temporary_position2)
            end
        end
        # plot the canvases in a 1*2 layout
        plot(canvas1, canvas2, layout = (1, 2))
    end
    # save the animation as a MP4 file
    gif(animation, "Results/Animation$num_name.mp4", fps = 60)
    # create an empty plot for curves
    curves = plotCurves()
    # add curves for the time evolution of each body
    distances ./= pc # change the unit of the distances: m --> pc
    for j in 1:num_bodies
        plot!(curves, [time_step*(1:num_steps)/Julian_year], distances[:, j], label="Body $j")
    end
    # save the curves as a PNG file
    savefig(curves, "Results/AnimationCurves$num_name.png")
end

# calculate total number of simulation steps
num_steps = Int(floor(total_time/time_step))

# run the selected main function and record the running time
running_time = @elapsed main_Curves(num_steps) # (s)
println("Running time: $running_time seconds")