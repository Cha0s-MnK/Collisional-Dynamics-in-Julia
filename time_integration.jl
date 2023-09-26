# Some time integration techniques to update the position and velocity of the bodies in 1 time step.
# Last edited by Cha0s_MnK on 2023/09/26.

# reference:
# Star Formation in Galaxy Evolution: Connecting Numerical Models to Reality
# Nickolay Y. Gnedin, Simon C.O. Glover, Ralf S. Klessen, Volker Springel
# High Performance Computing and Numerical Modelling (Volker Springel)
# 3 Time Integration Techniques

function timestep_criterion(bodies::Vector{Body})
    #

    timestep_min = Timestep * Julian_year
    for i in 1:num_bodies
        individual_timestep = 1.0e2 / norm(bodies[i].acceleration)
        timestep_min = min(timestep_min, individual_timestep)
    end
    global timestep = timestep_min / Julian_year # (Julian year)
    global TotalTime += timestep 
    if timestep < Timestep
        println("New timestep = $timestep Julian year")
    end
    global timestep = timestep_min
    global total_time += timestep
end

function explicitEuler!(bodies::Vector{Body}, force_calculation!::Function)
    # use explicit Euler method

    for i in 1:num_bodies
        bodies[i].position += bodies[i].velocity * timestep
        bodies[i].velocity += bodies[i].acceleration * timestep
    end
    force_calculation!(bodies) # calculate accelerations
end

function RungeKutta2nd!(bodies::Vector{Body}, force_calculation!::Function)
    # use 2nd order Runge-Kutta method

    # calculate k_1 of 2nd order Runge-Kutta method
    k_11 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    k_12 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    force_calculation!(bodies) # calculate accelerations
    for i in 1:num_bodies
        k_11[i] = bodies[i].velocity
        k_12[i] = bodies[i].acceleration
    end
    # calculate k_2 of 2nd order Runge-Kutta method
    k_21 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    k_22 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    temporary_bodies = deepcopy(bodies) # calculate the intermediate state of the bodies
    for i in 1:num_bodies
        temporary_bodies[i].position += k_11[i] * timestep
    end
    force_calculation!(temporary_bodies) # calculate accelerations
    for i in 1:num_bodies
        k_21[i] = bodies[i].velocity + k_12[i] * timestep
        k_22[i] = temporary_bodies[i].acceleration
    end
    # update positions, velocities and accelerations of all bodies
    for i in 1:num_bodies
        bodies[i].position += (k_11[i] + k_21[i]) * timestep / 2
        bodies[i].velocity += (k_12[i] + k_22[i]) * timestep / 2
    end
    force_calculation!(bodies) # calculate accelerations
end

function RungeKutta4th!(bodies::Vector{Body}, force_calculation!::Function)
    # use 4th order Runge-Kutta method

    # calculate k_1 of 4th order Runge-Kutta method
    k_11 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    k_12 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    force_calculation!(bodies) # calculate accelerations
    for i in 1:num_bodies
        k_11[i] = bodies[i].velocity
        k_12[i] = bodies[i].acceleration
    end
    # calculate k_2 of 4th order Runge-Kutta method
    k_21 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    k_22 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    temporary_bodies = deepcopy(bodies) # calculate the intermediate state of the bodies
    for i in 1:num_bodies
        temporary_bodies[i].position += k_11[i] * timestep / 2
    end
    force_calculation!(temporary_bodies) # calculate accelerations
    for i in 1:num_bodies
        k_21[i] = bodies[i].velocity + k_12[i] * timestep / 2
        k_22[i] = temporary_bodies[i].acceleration
    end
    # calculate k_3 of 4th order Runge-Kutta method
    k_31 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    k_32 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    temporary_bodies = deepcopy(bodies) # calculate the intermediate state of the bodies
    for i in 1:num_bodies
        temporary_bodies[i].position += k_21[i] * timestep / 2
    end
    force_calculation!(temporary_bodies) # calculate accelerations
    for i in 1:num_bodies
        k_31[i] = bodies[i].velocity + k_22[i] * timestep / 2
        k_32[i] = temporary_bodies[i].acceleration
    end
    # calculate k_4 of 4th order Runge-Kutta method
    k_41 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    k_42 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    temporary_bodies = deepcopy(bodies) # calculate the intermediate state of the bodies
    for i in 1:num_bodies
        temporary_bodies[i].position += k_31[i] * timestep
    end
    force_calculation!(temporary_bodies) # calculate accelerations
    for i in 1:num_bodies
        k_41[i] = bodies[i].velocity + k_32[i] * timestep
        k_42[i] = temporary_bodies[i].acceleration
    end
    # update positions, velocities and accelerations of all bodies
    for i in 1:num_bodies
        bodies[i].position += (k_11[i]/6 + k_21[i]/3 + k_31[i]/3 + k_41[i]/6) * timestep
        bodies[i].velocity += (k_12[i]/6 + k_22[i]/3 + k_32[i]/3 + k_42[i]/6) * timestep
    end
    force_calculation!(bodies) # calculate accelerations
end

function Leapfrog!(bodies::Vector{Body}, force_calculation!::Function)
    # use Leapfrog method

    for i in 1:num_bodies
        bodies[i].velocity += bodies[i].acceleration .* timestep ./ 2
        bodies[i].position += bodies[i].velocity .* timestep
    end
    force_calculation!(bodies) # calculate accelerations
    for i in 1:num_bodies
        bodies[i].velocity += bodies[i].acceleration .* timestep ./ 2
    end
    force_calculation!(bodies) # calculate accelerations again
end

function modified4thHermite!(bodies::Vector{Body}, force_calculation!::Function)
    # use modified 4th-order Hermite method

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