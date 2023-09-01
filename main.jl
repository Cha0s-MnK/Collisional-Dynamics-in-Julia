# The main running file of the N-body simulator.
# Last edited by Cha0s_MnK on 2023/09/01.

# "BHoctree" is short for "Barnes & Hut oct-tree"
# "OpeningAngle" is short for the "opening angle" criterion
# original paper:
# A hierarchical O(N log N) force-calculation algorithm
# https://www.nature.com/articles/324446a0

# use necessary packages as long as set constants and global variables.
include("head.jl")
# use the simple N-body simulator
include("simulator.jl")
# use BHoctree
include("BHoctree.jl")
# use time integration functions
include("time_integration.jl")
# use self-defined plot functions
include("plot.jl")

# main functions
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
running_time = @elapsed main_AnimationCurves(num_steps) # (s)
println("Running time: $running_time seconds")