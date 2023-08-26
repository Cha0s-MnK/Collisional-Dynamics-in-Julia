# This is the file about basic features of N-body simulator.
# Last edited by Chaos on 2023/08/26.
# still under tests

# use necessary packages
include("using.jl")
# use constants and global variables
include("global.jl")

# struct(s)
struct Body
    mass::Float64 
    normalized_mass::Float64
    position::MVector{3, Float64}
    velocity::MVector{3, Float64}
end

# functions
function generateBody()
    # generate 1 random body

    # generate random mass (Gaussion distribution) and normalize it
    if standard_deviation_mass == 0 # zero deviation
        mass = mean_mass
        normalized_mass = 1
    else
        min_mass = mean_mass - 3 * standard_deviation_mass # minimum mass generated
        max_mass = mean_mass + 3 * standard_deviation_mass # maximum mass generated
        mass_distribution = truncated(Normal(mean_mass, standard_deviation_mass), min_mass, max_mass)
        mass = rand(mass_distribution)
        normalized_mass = (mass - min_mass) / (6 * standard_deviation_mass)
    end
    # generate random coordinate (uniform distribution)
    x0 = rand() # generate a random floating-point number between 0 and 1
    y0 = rand()
    z0 = rand()
    x = x0 * total_length .- total_length / 2 # switch the random number from between 0 and 1 to between -total_length/2 and total_length/2
    y = y0 * total_length .- total_length / 2
    z = z0 * total_length .- total_length / 2
    position = MVector{3, Float64}(x, y, z)
    # initialize velocity to zero
    velocity = MVector{3, Float64}(0.0, 0.0, 0.0)

    return Body(mass, normalized_mass, position, velocity)
end

function updateBody!(body::Body, force::MVector{3, Float64})
    # update position and velocity of a body in 1 time step
    body.position .+= body.velocity * time_step
    body.velocity .+= (force ./ body.mass) * time_step
    # apply periodic boundary condition to modify coordinates
    for i in 1:3
        body.position[i] = mod(body.position[i] + total_length/2, total_length) - total_length/2
    end
end

function add_point_force!(force::MVector{3, Float64}, body::Body, mass::Float64, position::MVector{3, Float64})
    # calculate and add 1 approximate force on body using Plummer force law
    soft_length = softening_coefficient*total_length
    if body.position != position # exclude the case that two positions overlap
        displacement = position - body.position # displacement vector
        distance = norm(displacement) # distance scalar
        force_magnitude = G * body.mass * mass / (soft_length + distance)^2
        force .+= force_magnitude .* displacement / distance
    end
end

function direct_summation_step!(bodies::Vector{Body})
    # 1 simulation step using direct summation

    # initialize forces for each body to zero
    forces = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    # calculate forces between all pairs of bodies
    temporary_force = MVector{3, Float64}(0.0, 0.0, 0.0)
    for i in 1:num_bodies
        for j in (i+1):num_bodies
            add_point_force!(temporary_force, bodies[i], bodies[j].mass, bodies[j].position)
            forces[i] .+= temporary_force
            forces[j] .-= temporary_force
            fill!(temporary_force, 0.0) # reset temporary force
        end
    end
    # update positions and velocities of all bodies
    for i in 1:num_bodies
        updateBody!(bodies[i], forces[i])
    end
end

function plotBodies!(bodies::Vector{Body}, title0::AbstractString)
    # plot current bodies on a canvas; return the canvas

    # create a beautiful 3D canvas
    plot_length = 0.6*total_length # "0.6*" to make it more beautiful
    canvas = plot3d(
        size = (1080, 720),
        legend = false,
        xlabel = "X", ylabel = "Y", zlabel = "Z",
        xlims = (-plot_length, plot_length), ylims = (-plot_length, plot_length), zlims = (-plot_length, plot_length),
        title = title0, titlefontfamily = "Times New Roman", titlefontsize = 14,
        tickfontfamily = "Times New Roman", tickfontsize = 10, 
        legendfontfamily = "Times New Roman",
        xlabelfontfamily = "Times New Roman", xlabelfontsize = 12, 
        ylabelfontfamily = "Times New Roman", ylabelfontsize = 12,
        zlabelfontfamily = "Times New Roman", zlabelfontsize = 12,
        background_color = RGB(0.95, 0.95, 0.95), 
        fg_legend = :white,
        linealpha = 0.7, linewidth = 1, linestyle = :auto, 
        marker = (:circle, 4),
        grid = true,
        color = :viridis,
        camera = (30, 30), elevation = 30, azimuth = 30
    )
    # plot current bodies on the canvas
    for body in bodies
        scatter!(canvas, [body.position[1]], [body.position[2]], [body.position[3]], markersize=5*body.normalized_mass, color=:blue) # "5*" to make it more beautiful
    end
    # return the canvas
    return canvas
end