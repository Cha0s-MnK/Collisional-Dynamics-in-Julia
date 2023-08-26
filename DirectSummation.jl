# N-body simulator: use direct summation
# Last edited by Chaos on 2023/08/20.

# original paper:
# A hierarchical O(N log N) force-calculation algorithm
# https://www.nature.com/articles/324446a0
# "BHoctree" is short for "Barnes & Hut oct-tree"
# "max" is short for "maximum"
# "min" is short for "minimum"
# "num" is short for "number"
# "OpeningAngle" is short for the "opening angle" criterion

# use necessary packages
using Distributions # to use functions: "truncated", "Normal"
using LinearAlgebra # to use function "norm"
using Plots # to use functions: "plot3d", "scatter", "animate", "gif"
gr() # switch to the "GR" backend
using Random # to use function "rand"
using StaticArrays # for simulation speed, to use: "SVector", "MVector"
using Statistics # to use function "elapsed"

# constants and global variables
const G = 6.67430e-11  # gravitational constant (m^3/kg/s^2)
const line_width = 0.5 # linewidth of the edge of cells when ploting
const mean_mass = 3000.0 # average mass of the bodies (kg)
const standard_deviation_mass = 80.0 # standard deviation of the masses of the bodies (kg)
const total_length = 100.0 # side length of the cube (m)
const total_time = 24*60.0*60.0 # total simulation time (s)
const time_step = 60.0 # time step (s)

global num_bodies = 16 # number of random bodies
global theta = 1 # opening angle (rad)

# structs
struct Body
    mass::Float64 
    normalized_mass::Float64
    position::MVector{3, Float64}
    velocity::MVector{3, Float64}
end

mutable struct Cell
    depth::Int # "0" for root cell
    center::SVector{3, Float64} # center postion of this cell
    subcells::Union{Nothing, Vector{Cell}} # subcells in this cell / Nothing (only for leaf cell)
    total_mass::Float64
    mass_center::MVector{3, Float64} # center of mass of this cell
end

# functions
function newCell(depth::Int, center::SVector{3, Float64})
    # create a new Cell
    return Cell(depth, center, nothing, 0.0, MVector(center))
end

function insertBody!(body::Body, cell::Cell)
    # insert a body into a cell
    while true # original recursion has been converted into a loop
        if cell.subcells === nothing # this cell is a leaf cell
            if cell.total_mass == 0 # this leaf cell is empty
                cell.total_mass = body.mass
                cell.mass_center = body.position
                break # insertion is over
            else # this leaf cell contains a body so we divide it into 8 subcells
                cell.subcells = Vector{Cell}(undef, 8)
                center_offset = total_length / 2^(cell.depth+2) # cell_length (total_length / 2^cell.depth) / 4
                for i in 1:8
                    x_offset = ((i-1) & 1 == 0) ? center_offset : -center_offset
                    y_offset = ((i-1) & 2 == 0) ? center_offset : -center_offset
                    z_offset = ((i-1) & 4 == 0) ? center_offset : -center_offset
                    cell.subcells[i] = newCell(cell.depth+1, SVector(cell.center[1] + x_offset, cell.center[2] + y_offset, cell.center[3] + z_offset))
                end
                insertBody!(Body(cell.total_mass, 0.0, cell.mass_center, MVector{3, Float64}(0.0,0.0,0.0)), cell) # insert old body into subcells
                cell.total_mass = 0.0 # only leaf cell can have "total_mass!=0" in this step
            end
        else # this cell is a branch cell
            cell = cell.subcells[get_subcell_index(body, cell)]
        end
    end
end

function get_subcell_index(body::Body, cell::Cell)
    # calculate index of the subcell that the body will be inserted into
    index = 1
    index += (body.position[1] > cell.center[1]) ? 0 : 1
    index += (body.position[2] > cell.center[2]) ? 0 : 2
    index += (body.position[3] > cell.center[3]) ? 0 : 4
    return index
end

function updateCell!(cell::Cell)
    # update total_mass and mass_center of a branch cell (and all its subcells, sub-subcells...); discard empty leaf cells
    if cell.subcells !== nothing # this cell is a branch cell
        if cell.total_mass == 0 # this branch cell needs to be updated
            cell.subcells = filter(keep_subcell, cell.subcells) # discard empty leaf cells
            mass_center_numerator = MVector{3, Float64}(0.0,0.0,0.0)
            for subcell in cell.subcells
                updateCell!(subcell) # ensure that each subcell has been updated
                cell.total_mass += subcell.total_mass
                mass_center_numerator .+= subcell.total_mass .* subcell.mass_center
            end
            cell.mass_center = mass_center_numerator ./ cell.total_mass
        end
    end
end

function keep_subcell(subcell::Cell)
    # select empty leaf cells and they will return false
    return subcell.total_mass != 0 || subcell.subcells !== nothing
end

function generateBody()
    # generate 1 random body

    # generate random mass (Gaussion distribution) and normalize it
    min_mass = mean_mass - 3 * standard_deviation_mass # minimum mass generated
    max_mass = mean_mass + 3 * standard_deviation_mass # maximum mass generated
    mass_distribution = truncated(Normal(mean_mass, standard_deviation_mass), min_mass, max_mass)
    mass = rand(mass_distribution)
    normalized_mass = (mass - min_mass) / (6 * standard_deviation_mass)

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
    # calculate and add 1 exact force on body using Newton's law of universal gravitation
    if body.position != position # exclude the case that two positions overlap
        displacement = position - body.position # displacement vector
        distance = norm(displacement) # distance scalar
        force_magnitude = G * body.mass * mass / distance^2
        force .+= force_magnitude .* displacement / distance
    end
end

function add_cell_force!(force::MVector{3, Float64}, body::Body, cell::Cell)
    # calculate and add approximate force on a body from a cell using BHoctree
    cell_length = total_length / 2^cell.depth
    distance = norm(cell.mass_center - body.position) # distance scalar
    if cell.subcells !== nothing && (cell_length/distance > theta || inCell(body, cell)) # this cell is a branch cell; it does not satisfy OpeningAngle or the body is in this cell
        for subcell in cell.subcells # improve the resolution of force calculation
            add_cell_force!(force, body, subcell)
        end
    else # take force from the whole cell as an approximation
        add_point_force!(force, body, cell.total_mass, cell.mass_center)
    end
end

function inCell(body::Body, cell::Cell)
    # judge whether a body is in a cell
    cell_length = total_length / 2^cell.depth
    cell_min = cell.center .- cell_length/2 # 3 minimum coordinates of the cell
    cell_max = cell.center .+ cell_length/2 # 3 maximum coordinates of the cell
    for i in 1:3
        if !(cell_min[i] < body.position[i] < cell_max[i])
            return false
        end
    end
    return true
end

function BHoctree_step!(bodies::Vector{Body})
    # 1 simulation step using BHoctree

    # initialize root cell
    root = newCell(0, SVector(0.0,0.0,0.0))
    # establish BHoctree and update it
    for body in bodies
        insertBody!(body, root)
    end
    updateCell!(root)
    # initialize forces for each body to zero
    forces = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    # # calculate approximate forces; update positions and velocities of all bodies
    for i in 1:num_bodies
        add_cell_force!(forces[i], bodies[i], root)
        updateBody!(bodies[i], forces[i])
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

function plotBodies!(canvas::Plots.Plot{Plots.GRBackend}, bodies::Vector{Body})
    # plot current bodies on the canvas
    for body in bodies
        scatter!(canvas, [body.position[1]], [body.position[2]], [body.position[3]], markersize=5*body.normalized_mass, color=:blue) # "5*" to make it more beautiful
    end
end

function my_plot3D(title::AbstractString)
    # a more beautiful 3D plot by Chaos
    plot_length = 0.6*total_length # "0.6*" to make it more beautiful
    plot3d(
        legend = false,
        xlabel = "X", ylabel = "Y", zlabel = "Z",
        xlims = (-plot_length, plot_length), ylims = (-plot_length, plot_length), zlims = (-plot_length, plot_length),
        title = title, titlefontfamily = "Times New Roman", titlefontsize = 16,
        tickfontfamily = "Times New Roman", tickfontsize = 8, 
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
end

function main_animation()
    # main animation function that compares direct summation and BHoctree

    # generate random bodies
    bodies1 = [generateBody() for _ in 1:num_bodies]
    bodies2 = copy(bodies1)
    bodies3 = copy(bodies1)
    # calculate total number of simulation steps
    num_steps = Int(total_time/time_step)
    # simulate N-body problem in different ways, create a animation and save it
    animation = @animate for _ in 1:num_steps
        # create canvases
        canvas1 = my_plot3D("Direct summation")
        canvas2 = my_plot3D("BHoctree\nOpening angle = 1")
        canvas3 = my_plot3D("BHoctree\nOpening angle = 10")

        # use direct summation to simulate
        plotBodies!(canvas1, bodies1)
        direct_summation_step!(bodies1)

        # use BHoctree to simulate
        global theta = 1
        plotBodies!(canvas2, bodies2)
        BHoctree_step!(bodies2)
        global theta = 10
        plotBodies!(canvas3, bodies3)
        BHoctree_step!(bodies3)

        plot(canvas1, canvas2, canvas3, layout = (1, 3)) # horizontal layout
    end
    gif(animation, "Julia/Compare1.mp4", fps = 60)
end

function main_BHoctree(N::Int)
    # main function to simulate N-body problem using BHoctree with different "N"(num_bodies)s

    # generate random bodies
    global num_bodies = N
    bodies = [generateBody() for _ in 1:num_bodies]
    # calculate total number of simulation steps
    num_steps = Int(total_time/time_step)
    # simulate N-body problem using BHoctree
    for _ in 1:num_steps
        BHoctree_step!(bodies)
    end
end

function main_direct_summation(N::Int)
    # main function to simulate N-body problem using direct summation with different "N"(num_bodies)s

    # generate random bodies
    global num_bodies = N
    bodies = [generateBody() for _ in 1:num_bodies]
    # calculate total number of simulation steps
    num_steps = Int(total_time/time_step)
    # simulate N-body problem using direct summation
    for _ in 1:num_steps
        direct_summation_step!(bodies)
    end
end

# increase N gradually and record the running times
Ns = [10, 31, 100, 316, 1000, 3162]#, 10000]
for N in Ns
    #running_time = @elapsed main_BHoctree(N) # (s)
    running_time = @elapsed main_direct_summation(N) # (s)
    println("Number of bodies: $N, Running time: $running_time seconds")
end