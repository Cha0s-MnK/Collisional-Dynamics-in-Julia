# These are some functions for gravitational force / acceleration calculation.
# Last edited by Cha0s_MnK on 2023/09/28.

# reference:
# Star Formation in Galaxy Evolution: Connecting Numerical Models to Reality (Nickolay Y. Gnedin, Simon C.O. Glover, Ralf S. Klessen, Volker Springel)
# High Performance Computing and Numerical Modelling (Volker Springel)
# 4 Gravitational Force Calculation

function get_force(body::Body, mass::Float64, position::MVector{3, Float64})::MVector{3, Float64}
    # calculate approximate force on a body using gravitational softening
 
    force = MVector{3, Float64}(0.0, 0.0, 0.0) # initialize force to 0
    if body.position != position # exclude the case that two positions overlap
        displacement = position - body.position # displacement vector
        distance = norm(displacement) # distance scalar
        force_magnitude = G * body.mass * mass / (softening_length + distance)^2
        force += force_magnitude .* displacement ./ distance
    end
    return force
end

function get2acceleration!(body1::Body, body2::Body) # get exact accelerations between 2 bodies
    displacement = body2.position - body1.position # displacement vector
    distance = norm(displacement) # distance scalar
    force = G * body1.mass * body2.mass .* displacement ./ distance^3 # force vector
    body1.acceleration += force ./ body1.mass
    body2.acceleration -= force ./ body2.mass # Newton 3rd law
end

function bound_pairs(bodies::Vector{Body})
    # criterion to avoid boud pairs

    velocity_square = [sum(body.velocity .^ 2) for body in bodies] # calculate all velocity squares
    mean_velocity_square = mean(velocity_square) # calculate mean velocity square
    varepsilon = G * mean_mass / mean_velocity_square / pc # varepsilon must be much larger than this lower limit (pc)
    println("varepsilon's lower limit: $varepsilon")
end

function direct_summation!(bodies::Vector{Body}) # direct summation method
    for body in bodies # reset accelerations
        body.acceleration .= 0.0
    end
    for i in 1:num_bodies
        for j in (i+1):num_bodies
            get2acceleration!(bodies[i], bodies[j])
        end
    end
end

function BHoctree!(bodies::Vector{Body})
    # calculate gravitational forces / accelerations for each body using BHoctree

    # initialize accelerations for each body to zero
    for i in 1:num_bodies
        bodies[i].acceleration .= 0.0
    end
    # establish BHoctree
    root = newCell(0, SVector(0.0,0.0,0.0)) # initialize root cell
    for body in bodies # insert bodies into BHoctree sequentially
        insertBody!(body, root)
    end
    updateCell!(root) # calculate monopoles and discard empty cells
    # calculate gravitational forces / accelerations for each body
    for i in 1:num_bodies
        bodies[i].acceleration = get_cell_force(bodies[i], root) ./ bodies[i].mass
    end
end

# Some struct(s) and functions to establish BHoctree.

# "BHoctree" is short for "Barnes & Hut oct-tree"
# "OpeningAngle" is short for the "opening angle" criterion
# original paper:
# A hierarchical O(N log N) force-calculation algorithm
# https://www.nature.com/articles/324446a0

# struct(s)
mutable struct Cell
    depth::Int # "0" for root cell
    center::SVector{3, Float64} # center postion of this cell
    subcells::Union{Nothing, Vector{Cell}} # subcells in this cell / Nothing (only for leaf cell)
    total_mass::Float64
    mass_center::MVector{3, Float64} # center position of mass of this cell
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
                insertBody!(Body(cell.total_mass, cell.mass_center, MVector{3, Float64}(0.0,0.0,0.0)), cell) # insert old body into subcells
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
                mass_center_numerator += subcell.total_mass .* subcell.mass_center
            end
            cell.mass_center = mass_center_numerator ./ cell.total_mass
        end
    end
end

function keep_subcell(subcell::Cell)
    # select empty leaf cells and they will return false

    return subcell.total_mass != 0 || subcell.subcells !== nothing
end

function get_cell_force(body::Body, cell::Cell)::MVector{3, Float64}
    # calculate approximate force on a body from a cell using BHoctree

    force = MVector{3, Float64}(0.0, 0.0, 0.0) # initialize force to 0
    cell_length = total_length / 2^cell.depth # calculate side length of this cell
    distance = norm(cell.mass_center - body.position) # calculate distance scalar
    if cell.subcells !== nothing && (cell_length/distance > theta || inCell(body, cell)) 
        # this cell is a branch cell; it does not satisfy OpeningAngle or the body is in this cell
        for subcell in cell.subcells # improve the resolution
            force += get_cell_force(body, subcell)
        end
    else # consider force from the whole cell as an approximation
        force += get_force(body, cell.total_mass, cell.mass_center)
    end
    return force
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