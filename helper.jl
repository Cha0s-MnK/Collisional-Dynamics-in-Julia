# These are some very useful functions.
# Last edited by Cha0s_MnK on 2023/09/26.

function get_energy(bodies::Vector{Body})::Float64 # get the total energy of the N-body system (J)
    # calculate kinetic energy
    T = sum(0.5 * body.mass * sum(body.velocity.^2) for body in bodies)
    # calculate potential energy
    V = 0.0
    for i in 1:num_bodies
        for j in (i+1):num_bodies
            distance = norm(bodies[i].position - bodies[j].position)
            V += - G * bodies[i].mass * bodies[j].mass / distance
        end
    end
    return T + V
end

function get2jerk!(body1::Body, body2::Body) # get jerks between 2 bodies
    displacement = body2.position - body1.position # displacement vector
    distance = norm(displacement) # distance scalar
    relative_velocity = body2.velocity - body1.velocity # relative velocity vector
    radial_velocity = dot(relative_velocity, displacement) / distance # radial velocity scalar
    pre_jerk = - G * radial_velocity / distance^3
    body1.jerk += body2.mass .* pre_jerk
    body2.jerk -= body1.mass .* pre_jerk
end

function get_jerks!(bodies::Vector{Body}) # get jerks for each body using direct summation
    # initialize jerks for each body to zero
    for i in 1:num_bodies
        bodies[i].jerk .= 0.0
    end
    # calculate jerks for each body
    for i in 1:num_bodies
        for j in (i+1):num_bodies
            get2jerk!(bodies[i], bodies[j])
        end
    end
end