# Some functions to calculate gravitational forces for each body.
# Last edited by Cha0s_MnK on 2023/09/26.

# reference:
# Star Formation in Galaxy Evolution: Connecting Numerical Models to Reality
# Nickolay Y. Gnedin, Simon C.O. Glover, Ralf S. Klessen, Volker Springel
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

function get_2jerks!(body1::Body, body2::Body)
    # calculate exact jerks on 2 bodies

    displacement = body2.position - body1.position # displacement vector
    distance = norm(displacement) # distance scalar
    relative_velocity = body2.velocity - body1.velocity # relative velocity vector
    radial_velocity = dot(relative_velocity, displacement) / distance # radial velocity scalar
    pre_jerk = - G * radial_velocity / distance^3
    body1.jerk += body2.mass .* pre_jerk
    body2.jerk -= body1.mass .* pre_jerk
end

function bound_pairs(bodies::Vector{Body})
    # criterion to avoid boud pairs

    velocity_square = [sum(body.velocity .^ 2) for body in bodies] # calculate all velocity squares
    mean_velocity_square = mean(velocity_square) # calculate mean velocity square
    varepsilon = G * mean_mass / mean_velocity_square / pc # varepsilon must be much larger than this lower limit (pc)
    println("varepsilon's lower limit: $varepsilon")
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