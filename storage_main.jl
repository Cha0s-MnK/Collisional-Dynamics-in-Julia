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