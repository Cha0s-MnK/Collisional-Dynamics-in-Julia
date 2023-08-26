# This is the main running file of N-body simulator.
# Last edited by Chaos on 2023/08/26.

# original paper:
# A hierarchical O(N log N) force-calculation algorithm
# https://www.nature.com/articles/324446a0

# "BHoctree" is short for "Barnes & Hut oct-tree"
# "OpeningAngle" is short for the "opening angle" criterion

# use necessary packages
include("using.jl")
# use constants and global variables
include("global.jl")
# to use BHoctree
include("BHoctree.jl")

function main_Bars()
    # main function that compares direct summation and BHoctree using bars

    # generate random bodies and deepcopy it
    bodies1 = [generateBody() for _ in 1:num_bodies]
    bodies2 = deepcopy(bodies1)
    # calculate total number of simulation steps
    num_steps = Int(total_time/time_step)
    # create an array to stroe distances
    distances = zeros(num_steps, num_bodies)
    # simulate N-body problem in different ways; store the distances
    for i in 1:num_steps
        # use direct summation to simulate
        direct_summation_step!(bodies1)
        # use BHoctree to simulate
        BHoctree_step!(bodies2)
        # calculate the distances and store them
        for j in 1:num_bodies
            #distances[i, j] = log(1+norm(bodies1[j].position - bodies2[j].position))
            distances[i, j] = norm(bodies1[j].position - bodies2[j].position)
        end
    end
    # create an empty bar
    bars = plot(size = (1080, 720),
                title = "Colour Bar", titlefontfamily = "Times New Roman", titlefontsize = 14,
                tickfontfamily = "Times New Roman", tickfontsize = 10, 
                legendfontfamily = "Times New Roman",
                xlabel = "Time(s)", xlabelfontfamily = "Times New Roman", xlabelfontsize = 12, 
                ylabel = "Distances(m)", ylabelfontfamily = "Times New Roman", ylabelfontsize = 12)
    # add bars for the time evolution of each body
    for j in 1:num_bodies
        bar!(bars, [time_step*(1:num_steps)], distances[:, j], label="Body $j", stack=true)
    end
    # save the bars as a PNG file
    savefig(bars, "Results/Bars$num_name.png")
end

function main_Curves()
    # main function that compares direct summation and BHoctree using curves

    # generate random bodies and deepcopy it
    bodies1 = [generateBody() for _ in 1:num_bodies]
    bodies2 = deepcopy(bodies1)
    # calculate total number of simulation steps
    num_steps = Int(total_time/time_step)
    # create an array to stroe distances
    distances = zeros(num_steps, num_bodies)
    # simulate N-body problem in different ways; store the distances
    for i in 1:num_steps
        # use direct summation to simulate
        direct_summation_step!(bodies1)
        # use BHoctree to simulate
        BHoctree_step!(bodies2)
        # calculate the distances and store them
        for j in 1:num_bodies
            #distances[i, j] = log(1+norm(bodies1[j].position - bodies2[j].position))
            distances[i, j] = norm(bodies1[j].position - bodies2[j].position)
        end
    end
    # create an empty plot
    curves = plot(
        size = (1080, 720),
        title = "Colour Bar", titlefontfamily = "Times New Roman", titlefontsize = 14,
        tickfontfamily = "Times New Roman", tickfontsize = 10, 
        legendfontfamily = "Times New Roman",
        xlabel = "Time(s)", xlabelfontfamily = "Times New Roman", xlabelfontsize = 12, 
        ylabel = "Distances(m)", ylabelfontfamily = "Times New Roman", ylabelfontsize = 12
    )
    # add curves for the time evolution of each body
    for j in 1:num_bodies
        plot!(curves, [time_step*(1:num_steps)], distances[:, j], label="Body $j")
    end
    # save the plot as a PNG file
    savefig(curves, "Results/Curves$num_name.png")
end

function main_Animation()
    # main function that compares direct summation and BHoctree using animation

    # generate random bodies and deepcopy it
    bodies1 = [generateBody() for _ in 1:num_bodies]
    bodies2 = deepcopy(bodies1)
    # calculate total number of simulation steps
    num_steps = Int(total_time/time_step)
    # simulate N-body problem in different ways; create an animation
    animation = @animate for i in 1:num_steps
        # plot current bodies; use direct summation to simulate
        canvas1 = plotBodies!(bodies1, "Direct summation\nSoft length = $softening_coefficient total length")
        direct_summation_step!(bodies1)
        # plot current bodies; use BHoctree to simulate
        canvas2 = plotBodies!(bodies2, "BHoctree\nOpening angle = $theta")
        BHoctree_step!(bodies2)
        # plot the canvases in a 1*2 layout
        plot(canvas1, canvas2, layout = (1, 2))
    end
    # save the animation as a MP4 file
    gif(animation, "Results/Animation$name.mp4", fps = 60)
end

function main_AnimationCurves()
    # main function that compares direct summation and BHoctree using animation and curves

    # generate random bodies and deepcopy it
    bodies1 = [generateBody() for _ in 1:num_bodies]
    bodies2 = deepcopy(bodies1)
    # calculate total number of simulation steps
    num_steps = Int(total_time/time_step)
    # create an array to stroe distances
    distances = zeros(num_steps, num_bodies)
    # simulate N-body problem in different ways; store the distances; create an animation
    animation = @animate for i in 1:num_steps
        # plot current bodies; use direct summation to simulate
        canvas1 = plotBodies!(bodies1, "Direct summation\nSoft length = $softening_coefficient total length")
        direct_summation_step!(bodies1)
        # plot current bodies; use BHoctree to simulate
        canvas2 = plotBodies!(bodies2, "BHoctree\nOpening angle = $theta")
        BHoctree_step!(bodies2)
        # calculate the distances and store them
        for j in 1:num_bodies
            #distances[i, j] = log(1+norm(bodies1[j].position - bodies2[j].position))
            distances[i, j] = norm(bodies1[j].position - bodies2[j].position)
        end
        # plot the canvases in a 1*2 layout
        plot(canvas1, canvas2, layout = (1, 2))
    end
    # save the animation as a MP4 file
    gif(animation, "Results/Animation$name.mp4", fps = 60)
    # create an empty plot
    curves = plot(
        size = (1080, 720),
        title = "Colour Bar", titlefontfamily = "Times New Roman", titlefontsize = 14,
        tickfontfamily = "Times New Roman", tickfontsize = 10, 
        legendfontfamily = "Times New Roman",
        xlabel = "Time(s)", xlabelfontfamily = "Times New Roman", xlabelfontsize = 12, 
        ylabel = "Distances(m)", ylabelfontfamily = "Times New Roman", ylabelfontsize = 12
    )
    # add curves for the time evolution of each body
    for j in 1:num_bodies
        plot!(curves, [time_step*(1:num_steps)], distances[:, j], label="Body $j")
    end
    # save the curves as a PNG file
    savefig(curves, "Results/AnimationCurves$num_name.png")
end

running_time = @elapsed main_Curves() # (s)
println("Running time: $running_time seconds")