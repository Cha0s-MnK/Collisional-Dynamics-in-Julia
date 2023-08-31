# This is the file to use BHoctree to simulate.
# Last edited by Chaos on 2023/08/28.

# original paper:
# A hierarchical O(N log N) force-calculation algorithm
# https://www.nature.com/articles/324446a0

# "BHoctree" is short for "Barnes & Hut oct-tree"
# "OpeningAngle" is short for the "opening angle" criterion

# struct(s)
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
    # calculate approximate forces
    for i in 1:num_bodies
        add_cell_force!(forces[i], bodies[i], root)
    end
    # update positions and velocities of all bodies
    for i in 1:num_bodies
        updateBody!(bodies[i], forces[i])
    end
end