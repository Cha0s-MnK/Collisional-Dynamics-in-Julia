# These are some time integration techniques to update the state of the bodies in 1 simulation step.
# Last edited by Cha0s_MnK on 2023/09/26.

# reference:
# Star Formation in Galaxy Evolution: Connecting Numerical Models to Reality
# Nickolay Y. Gnedin, Simon C.O. Glover, Ralf S. Klessen, Volker Springel
# High Performance Computing and Numerical Modelling (Volker Springel)
# 3 Time Integration Techniques

function timestep_criterion(bodies::Vector{Body}) # timestep criterion
    global timestep = MaxTimestep * Julian_year # any timestep should be no larger than the user-specified maximum timestep (s)
    for i in 1:num_bodies # find the minimum timestep
        individual_timestep = 1.0e2 / norm(bodies[i].acceleration) # a free parameter to control the precision of simulation
        global timestep = min(timestep, individual_timestep)
    end
    global Timestep = timestep / Julian_year # current timestep (Julian year)
    if Timestep < MaxTimestep # print the shortened timestep
        println("New timestep = $Timestep Julian year")
    end
    global TotalTime += Timestep # add current timestep to the total simulation time (Julian year)
end

function explicitEuler!(bodies::Vector{Body}, force_calculation!::Function) # explicit Euler method
    for (_, body) in enumerate(bodies)
        body.position += body.velocity * timestep
        body.velocity += body.acceleration * timestep
    end
    force_calculation!(bodies) # calculate accelerations for future use
end

function RungeKutta2nd!(bodies::Vector{Body}, force_calculation!::Function) # 2nd order Runge-Kutta method
    last_positions = [copy(body.position) for body in bodies] # store last positions
    last_velocities = [copy(body.velocity) for body in bodies] # store last velocities
    # calculate k_1 of 2nd order Runge-Kutta method
    k_12 = [copy(body.acceleration) for body in bodies]
    # calculate k_2 of 2nd order Runge-Kutta method
    for (i, body) in enumerate(bodies)
        body.position = last_positions[i] + last_velocities[i] * timestep
    end
    force_calculation!(bodies) # recalculate accelerations
    k_22 = [copy(body.acceleration) for body in bodies]
    # update positions, velocities and accelerations of all bodies
    for (i, body) in enumerate(bodies)
        body.position = last_positions[i] + last_velocities[i] * timestep + k_12[i] * timestep^2 / 2
        body.velocity = last_velocities[i] + (k_12[i] + k_22[i]) * timestep / 2
    end
    force_calculation!(bodies) # calculate accelerations for future use
end

function RungeKutta4th!(bodies::Vector{Body}, force_calculation!::Function) # 4th order Runge-Kutta method
    last_positions = [copy(body.position) for body in bodies] # store last positions
    last_velocities = [copy(body.velocity) for body in bodies] # store last velocities
    # calculate k_1 of 4th order Runge-Kutta method
    k_12 = [copy(body.acceleration) for body in bodies]
    # calculate k_2 of 4th order Runge-Kutta method
    for (i, body) in enumerate(bodies)
        body.position = last_positions[i] + last_velocities[i] * timestep / 2
    end
    force_calculation!(bodies) # recalculate accelerations
    k_22 = [copy(body.acceleration) for body in bodies]
    # calculate k_3 of 4th order Runge-Kutta method
    for (i, body) in enumerate(bodies)
        body.position = last_positions[i] + last_velocities[i] * timestep / 2 + k_22[i] * timestep^2 / 4
    end
    force_calculation!(bodies) # recalculate accelerations
    k_32 = [copy(body.acceleration) for body in bodies]
    # calculate k_4 of 4th order Runge-Kutta method
    for (i, body) in enumerate(bodies)
        body.position = last_positions[i] + last_velocities[i] * timestep + k_32[i] * timestep^2 / 2
    end
    force_calculation!(bodies) # recalculate accelerations
    k_42 = [copy(body.acceleration) for body in bodies]
    # update positions, velocities and accelerations of all bodies
    for (i, body) in enumerate(bodies)
        body.position = last_positions[i] + last_velocities[i] * timestep + (k_12[i] + k_22[i] + k_32[i]) / 6 * timestep^2
        body.velocity = last_velocities[i] + (k_12[i] + 2*k_22[i] + 2*k_32[i] + k_42[i]) / 6 * timestep
    end
    force_calculation!(bodies) # calculate accelerations for future use
end

function Leapfrog!(bodies::Vector{Body}, force_calculation!::Function) # use Leapfrog method (KDK)
    for body in bodies
        body.velocity += body.acceleration * timestep / 2 # half step velocity update (Kick)
        body.position += body.velocity * timestep # full step position update (Drift)
    end
    force_calculation!(bodies) # recalculate accelerations
    for body in bodies
        body.velocity += body.acceleration * timestep / 2 # half step velocity update (Kick)
    end
    force_calculation!(bodies) # calculate accelerations for future use
end

function modified4thHermite!(bodies::Vector{Body}, force_calculation!::Function) # modified 4th-order Hermite method
    # calculate jerks
    direct_summation_jerk!(bodies)
    # calculate the intermediate state
    temporary_bodies = deepcopy(bodies)
    for i in 1:num_bodies
        temporary_bodies[i].position += temporary_bodies[i].position + temporary_bodies[i].velocity .* timestep + temporary_bodies[i].acceleration .* timestep^2 ./ 2 + temporary_bodies[i].jerk .* timestep^3 ./ 6
        temporary_bodies[i].velocity += temporary_bodies[i].velocity + temporary_bodies[i].acceleration .* timestep + temporary_bodies[i].jerk .* timestep^2 ./ 2
    end
    # calculate accelerations and jerks again
    force_calculation!(temporary_bodies)
    direct_summation_jerk!(temporary_bodies)
    # perform the correction step
    for i in 1:num_bodies
        bodies[i].velocity += bodies[i].velocity + (bodies[i].acceleration + temporary_bodies[i].acceleration) .* timestep ./ 2 + (bodies[i].jerk - temporary_bodies[i].jerk) .* timestep^2 ./ 12
        bodies[i].position += bodies[i].position + (bodies[i].velocity + temporary_bodies[i].velocity) .* timestep ./ 2 + (bodies[i].acceleration - temporary_bodies[i].acceleration) .* timestep^2 ./ 12
    end
    # calculate accelerations again again
    force_calculation!(bodies)
end