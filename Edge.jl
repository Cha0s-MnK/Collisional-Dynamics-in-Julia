# codes about N-body simulator using BHoctree
# Last edited by Chaos on 2023/08/14.

# "BHoctree" is short for "Barnes & Hut oct-tree"
# "OpeningAngle" is short for the "opening angle" criterion
# original paper:
# A hierarchical O(N log N) force-calculation algorithm
# https://www.nature.com/articles/324446a0

using Distributions # to use functions: "truncated", "Normal"
using LinearAlgebra # to use function "norm"
using Plots # to use functions: "plot3d", "scatter", "animate", "gif"
gr() # switch to the "GR" backend
using Random # to use function "rand"
using StaticArrays # for simulation speed, to use: "SVector", "MVector"
using Statistics # to use function "elapsed"

# constants
const G = 6.67430e-11  # gravitational constant (m^3/kg/s^2)
const total_length = 100.0 # side length of the cube (m)
const theta = 1 # opening angle (rad)
const mean_mass = 3000.0 # average mass of the bodies (kg)
const standard_deviation_mass = 100.0 # standard deviation of the masses of the bodies (kg)
const time_step = 60.0 # time step (s)
const total_time = 24.0*60.0*time_step # total simulation time (s)
const line_width = 0.5 # linewidth of the edge of cells when ploting

global number_bodies = 10 # number of random bodies

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
    supercell::Union{Nothing, Cell} # supercell of this cell / Nothing (only for root cell)
    subcells::Union{Nothing, Vector{Cell}} # subcells in this cell / Nothing (only for leaf cell)
    total_mass::Float64
    mass_center::MVector{3, Float64} # center of mass of this cell
end

# functions
function new_cell(depth::Int, center::SVector{3, Float64}, supercell::Union{Nothing, Cell})
    # create a new Cell
    return Cell(depth, center, supercell, nothing, 0.0, center)
end

function insert_body!(body::Body, cell::Cell)
    # insert a Body into a cell of a BHoctree
    if cell.subcells == nothing # this is a leaf cell
        if cell.total_mass == 0 # this is a empty leaf cell
            cell.total_mass = body.mass
            cell.mass_center = body.position
        else # this is leaf cell containing a body so divide it into 8 subcells
            cell.subcells = Vector{Cell}(undef, 8)
            center_offset = total_length / 2^(cell.depth+2)
            for i in 1:8
                x_offset = ((i-1) & 1 == 0) ? center_offset : -center_offset
                y_offset = ((i-1) & 2 == 0) ? center_offset : -center_offset
                z_offset = ((i-1) & 4 == 0) ? center_offset : -center_offset
                subcell_center = SVector(cell.center[1] + x_offset, cell.center[2] + y_offset, cell.center[3] + z_offset)
                cell.subcells[i] = new_cell(cell.depth+1, subcell_center, cell)
            end
            insert_body!(Body(cell.total_mass, 0.0, cell.mass_center, MVector(0.0,0.0,0.0)), cell) # produce a pseudo body to insert
            cell.total_mass = 0.0 # only leaf cell can have "total_mass!=0" in this step
            insert_body!(body, cell)
        end
    else # this is a branch cell
        insert_body!(body, cell.subcells[get_subcell_index(body, cell)])
    end
end

function insert_body!!(body::Body, cell::Cell)
    # insert a Body into a cell of a BHoctree
    while true
        if cell.subcells == nothing # this is a leaf cell
            if cell.total_mass == 0 # this is a empty leaf cell
                cell.total_mass = body.mass
                cell.mass_center = body.position
                break
            else # this is leaf cell containing a body so divide it into 8 subcells
                cell.subcells = Vector{Cell}(undef, 8)
                center_offset = total_length / 2^(cell.depth+2)
                for i in 1:8
                    x_offset = ((i-1) & 1 == 0) ? center_offset : -center_offset
                    y_offset = ((i-1) & 2 == 0) ? center_offset : -center_offset
                    z_offset = ((i-1) & 4 == 0) ? center_offset : -center_offset
                    subcell_center = SVector(cell.center[1] + x_offset, cell.center[2] + y_offset, cell.center[3] + z_offset)
                    cell.subcells[i] = new_cell(cell.depth+1, subcell_center, cell)
                end
                insert_body!!(Body(cell.total_mass, 0.0, cell.mass_center, MVector(0.0,0.0,0.0)), cell)
                cell.total_mass = 0.0 # only leaf cell can have "total_mass!=0" in this step
            end
        else # this is a branch cell
            cell = cell.subcells[get_subcell_index(body, cell)]
        end
    end
end

function get_subcell_index(body::Body, cell::Cell)
    # calculate the index of the subcell which the body goes into
    index = 1
    index += ifelse(body.position[1] > cell.center[1], 0, 1)
    index += ifelse(body.position[2] > cell.center[2], 0, 2)
    index += ifelse(body.position[3] > cell.center[3], 0, 4)
    return index
end

function updateCell!(cell::Cell)
    # update total_mass and mass_center of a branch cell; discard empty leaf cells
    if cell.subcells != nothing # this is a branch cell
        cell.subcells = filter(keep_subcell, cell.subcells) # discard empty leaf cells
        if cell.total_mass == 0 # this branch cell needs to be mass-updated
            mass_center_numerator = MVector(0.0,0.0,0.0)
            for subcell in cell.subcells
                updateCell!(subcell)
                cell.total_mass += subcell.total_mass
                mass_center_numerator .+= subcell.total_mass .* subcell.mass_center
            end
            cell.mass_center = mass_center_numerator / cell.total_mass
        end
    end
end

function keep_subcell(subcell::Cell)
    # choose empty leaf cells; they return false
    return subcell.total_mass != 0 || subcell.subcells != nothing
end

function generate_body()
    # generate 1 random body

    # generate random mass (Gaussion distribution) and normalize it
    min_mass = mean_mass - 3 * standard_deviation_mass
    max_mass = mean_mass + 3 * standard_deviation_mass
    mass_distribution = truncated(Normal(mean_mass, standard_deviation_mass), min_mass, max_mass)
    mass = rand(mass_distribution)
    normalized_mass = (mass - min_mass) / (6 * standard_deviation_mass)

    # generate random coordinate (uniform distribution)
    x0 = rand()
    y0 = rand()
    z0 = rand()
    x = x0 * total_length .- total_length / 2
    y = y0 * total_length .- total_length / 2
    z = z0 * total_length .- total_length / 2
    position = [x, y, z]

    # initialize velocity to zero
    velocity = [0.0, 0.0, 0.0]

    return Body(mass, normalized_mass, position, velocity)
end

function update_position_velocity!(body::Body, force::Vector{Float64})
    # update position and velocity in "dt" time 
    body.velocity .+= (force ./ body.mass) * time_step
    body.position .+= body.velocity * time_step
    for i in 1:3
        body.position[i] = periodic_boundary(body.position[i])
    end
end

function periodic_boundary(coordinate::Float64)
    # apply periodic boundary condition to modify coordinates
    return (coordinate + total_length/2) % total_length - total_length/2
end

function get1force(body::Body, mass::Float64, position::MVector{3, Float64})
    # calculate 1 exact force on "body" using Newton's law of universal gravitation
    force = zeros(3)
    if body.position != position # exclude the case that two positions coincide
        displacement = position - body.position
        distance = norm(displacement)
        force_magnitude = G * body.mass * mass / distance^2
        force .= force_magnitude .* displacement / distance
    end
    return force
end

function get_cell_force(body::Body, cell::Cell)
    # calculate approximate force on "body" from "cell" using BHoctree
    force = zeros(3)
    cell_length = total_length / 2^cell.depth
    distance = norm(cell.mass_center - body.position)
    if cell.subcells === nothing || cell_length/distance < theta # this is a leaf cell or it satisfies OpeningAngle
        force .= get1force(body, cell.total_mass, cell.mass_center)
    else
        for subcell in cell.subcells
            force .+= get_cell_force(body, subcell)
        end
    end
    return force
end

function BHoctree_step!(bodies::Vector{Body}, root::Cell)
    # 1 simulation step using BHoctree

    # calculate approximate forces and update positions and velocities for all bodies
    forces = [zeros(3) for _ in 1:number_bodies]
    for i in 1:number_bodies
        forces[i] .= get_cell_force(bodies[i], root)
        update_position_velocity!(bodies[i], forces[i])
    end
end

function NSquare_step!(bodies::Vector{Body})
    # 1 simulation step

    # initialize forces for each body to zero
    forces = [zeros(3) for _ in 1:number_bodies]
    temporary_force = zeros(3)

    # calculate forces between all pairs of bodies
    for i in 1:number_bodies
        for j in (i+1):number_bodies
            temporary_force .= get1force(bodies[i], bodies[j].mass, bodies[j].position)
            forces[i] .+= temporary_force
            forces[j] .-= temporary_force
        end
    end

    # update positions and velocities for all bodies
    for i in 1:number_bodies
        update_position_velocity!(bodies[i], forces[i])
    end
end

function plot_bodies(bodies)
    # plot current bodies
    plot_length = 0.6*total_length # "0.6*" to make it more beautiful
    plot3d(
        legend=false,
        xlabel="X", ylabel="Y", zlabel="Z",
        xlims=(-plot_length, plot_length), ylims=(-plot_length, plot_length), zlims=(-plot_length, plot_length)
    )
    for body in bodies
        scatter!([body.position[1]], [body.position[2]], [body.position[3]], markersize=5*body.normalized_mass, color=:blue) # "5*" to make it more beautiful
    end
end

function plot_leaf_cells(cell::Cell)
    # preorder traversal BHoctree and plot non-empty leaf cells 
    if cell.subcells != nothing # this is a branch cell
        for subcell in cell.subcells
            plot_leaf_cells(subcell)
        end
    else # this is a leaf cell (after cell update, it can not be empty)
        # calculate the corners of the cell
        cell_length = total_length / 2^cell.depth
        corners = [cell.center .+ cell_length / 2 .* MVector(x, y, z) for x in [-1, 1] for y in [-1, 1] for z in [-1, 1]]

        # plot the edges of the cell
        edges = [(1, 2), (1, 3), (2, 4), (3, 4), (5, 6), (5, 7), (6, 8), (7, 8), (1, 5), (2, 6), (3, 7), (4, 8)]
        for edge in edges
            i, j = edge
            plot!([corners[i][1], corners[j][1]],[corners[i][2], corners[j][2]],[corners[i][3], corners[j][3]], color=:blue, linewidth=line_width)
        end
    end
end

function main_animation()
    # main animation function

    # generate random bodies
    bodies = [generate_body() for _ in 1:number_bodies]

    # simulate N-body problem, create a animation and save it
    total_times = Int(total_time/time_step) # total simulation total_times
    animation = @animate for _ in 1:total_times
        plot_bodies(bodies)

        # create a BHoctree and show it
        root = new_cell(0, SVector(0.0,0.0,0.0), nothing)
        for body in bodies
            insert_body!!(body, root)
        end
        updateCell!(root)
        plot_leaf_cells(root)
        
        BHoctree_step!(bodies, root)
    end
    gif(animation, "BHoctree.mp4", fps = 60)
end

function main(N::Int)
    # main function

    # generate random bodies
    global number_bodies = N
    bodies = [generate_body() for _ in 1:number_bodies]

    # simulate N-body problem
    total_times = Int(total_time/time_step) # total simulation total_times
    for _ in 1:total_times
        # create a BHoctree
        root = new_cell(0, SVector(0.0,0.0,0.0), nothing)
        for body in bodies
            insert_body!!(body, root)
        end
        updateCell!(root)
        
        BHoctree_step!(bodies, root)
    end
end

#main_animation()

# increase N gradually and record the running times
Ns = [10, 31, 100, 316, 1000, 3162, 10000]
for N in Ns
    running_time = @elapsed main(N) # (s)
    println("Number of bodies: $N, Time: $running_time seconds")
end