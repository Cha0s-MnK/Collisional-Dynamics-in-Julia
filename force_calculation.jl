# Some functions to calculate gravitational forces for each body.
# Last edited by Cha0s_MnK on 2023/09/25.

# reference:
# Star Formation in Galaxy Evolution: Connecting Numerical Models to Reality
# (Nickolay Y. Gnedin, Simon C.O. Glover, Ralf S. Klessen, Volker Springel)
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

function bound_pairs(bodies::Vector{Body})
    # criterion to avoid boud pairs

    velocity_square = [sum(body.velocity .^ 2) for body in bodies] # calculate all velocity squares
    mean_velocity_square = mean(velocity_square) # calculate mean velocity square
    varepsilon = G * mean_mass / mean_velocity_square / pc # varepsilon must be much larger than this lower limit (pc)
    println("varepsilon's lower limit: $varepsilon")
end

function direct_summation!!(forces::Vector{MVector{3, Float64}}, bodies::Vector{Body})
    # calculate gravitational forces for each body using direct summation

    # initialize forces for each body to zero
    for i in 1:num_bodies
        forces[i] .= 0.0
    end
    # calculate gravitational forces for each body using direct summation
    temporary_force = MVector{3, Float64}(0.0, 0.0, 0.0) # create a temporary vector to store force
    for i in 1:num_bodies
        for j in (i+1):num_bodies
            temporary_force = get_force(bodies[i], bodies[j].mass, bodies[j].position)
            forces[i] += temporary_force
            forces[j] -= temporary_force # Newton’s 3rd law
        end
    end
    # calculate accelerations for each body
    for i in 1:num_bodies
        bodies[i].acceleration = forces[i] / bodies[i].mass
    end
end

function direct_summation!(bodies::Vector{Body})
    # calculate gravitational forces / accelerations for each body using direct summation

    # initialize accelerations for each body to zero
    for i in 1:num_bodies
        bodies[i].acceleration .= 0.0
    end
    # calculate gravitational forces / accelerations for each body
    temporary_force = MVector{3, Float64}(0.0, 0.0, 0.0) # create a temporary vector to store force
    for i in 1:num_bodies
        for j in (i+1):num_bodies
            temporary_force = get_force(bodies[i], bodies[j].mass, bodies[j].position)
            bodies[i].acceleration += temporary_force ./ bodies[i].mass
            bodies[j].acceleration -= temporary_force ./ bodies[j].mass # Newton’s 3rd law
        end
    end
end

function BHoctree!!(forces::Vector{MVector{3, Float64}}, bodies::Vector{Body})
    # calculate gravitational forces for each body using BHoctree

    # initialize forces for each body to zero
    for i in 1:num_bodies
        forces[i] .= 0.0
    end
    # establish BHoctree
    root = newCell(0, SVector(0.0,0.0,0.0)) # initialize root cell
    for body in bodies # insert bodies into BHoctree sequentially
        insertBody!(body, root)
    end
    updateCell!(root) # calculate monopoles and discard empty cells
    # calculate gravitational forces for each body using BHoctree
    for i in 1:num_bodies
        forces[i] = get_cell_force(bodies[i], root)
    end
    # calculate accelerations for each body
    for i in 1:num_bodies
        bodies[i].acceleration = forces[i] / bodies[i].mass
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