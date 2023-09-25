# Some time integration techniques to update the position and velocity of the bodies in 1 time step.
# Last edited by Cha0s_MnK on 2023/09/25.

# reference:
# Star Formation in Galaxy Evolution: Connecting Numerical Models to Reality
# Nickolay Y. Gnedin, Simon C.O. Glover, Ralf S. Klessen, Volker Springel
# High Performance Computing and Numerical Modelling (Volker Springel)
# 3 Time Integration Techniques

function timestep_criterion(bodies::Vector{Body})
    timestep_min = Timestep * Julian_year
    for i in 1:num_bodies
        individual_timestep = 1.0e5 / norm(bodies[i].acceleration)
        timestep_min = min(timestep_min, individual_timestep)
    end
    global timestep = timestep_min
    println("New timestep = $timestep seconds")
end

function explicitEuler!!(bodies::Vector{Body}, force_calculation!::Function)
    # use explicit Euler method

    # create a matrix to store forces
    forces = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    # calculate gravitational forces
    force_calculation!(forces, bodies)
    # update positions and velocities of all bodies
    for i in 1:num_bodies
        bodies[i].position += bodies[i].velocity * timestep
        bodies[i].velocity += forces[i] / bodies[i].mass * timestep
    end
end

function explicitEuler!(bodies::Vector{Body}, force_calculation!::Function)
    # use explicit Euler method

    for i in 1:num_bodies
        bodies[i].position += bodies[i].velocity * timestep
        bodies[i].velocity += bodies[i].acceleration * timestep
    end
    force_calculation!(bodies) # calculate accelerations
end

function RungeKutta2nd!!(bodies::Vector{Body}, force_calculation!::Function)
    # use 2nd order Runge-Kutta method

    # create a matrix to store forces
    forces = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    # calculate k_1 of 2nd order Runge-Kutta method
    k_11 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    k_12 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    force_calculation!(forces, bodies) # calculate gravitational forces
    for i in 1:num_bodies
        k_11[i] = bodies[i].velocity
        k_12[i] = forces[i] / bodies[i].mass
    end
    # calculate k_2 of 2nd order Runge-Kutta method
    k_21 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    k_22 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    temporary_bodies = deepcopy(bodies) # calculate the intermediate state of the bodies
    for i in 1:num_bodies
        temporary_bodies[i].position += k_11[i] * timestep
    end
    force_calculation!(forces, temporary_bodies) # calculate gravitational forces
    for i in 1:num_bodies
        k_21[i] = bodies[i].velocity + k_12[i] * timestep
        k_22[i] = forces[i] / bodies[i].mass
    end
    # update positions and velocities of all bodies
    for i in 1:num_bodies
        bodies[i].position += (k_11[i] + k_21[i]) * timestep / 2
        bodies[i].velocity += (k_12[i] + k_22[i]) * timestep / 2
    end
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

function RungeKutta4th!!(bodies::Vector{Body}, force_calculation!::Function)
    # use 4th order Runge-Kutta method

    # create a matrix to store forces
    forces = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    # calculate k_1 of 4th order Runge-Kutta method
    k_11 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    k_12 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    force_calculation!(forces, bodies) # calculate gravitational forces
    for i in 1:num_bodies
        k_11[i] = bodies[i].velocity
        k_12[i] = forces[i] / bodies[i].mass
    end
    # calculate k_2 of 4th order Runge-Kutta method
    k_21 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    k_22 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    temporary_bodies = deepcopy(bodies) # calculate the intermediate state of the bodies
    for i in 1:num_bodies
        temporary_bodies[i].position += k_11[i] * timestep / 2
    end
    force_calculation!(forces, temporary_bodies) # calculate gravitational forces
    for i in 1:num_bodies
        k_21[i] = bodies[i].velocity + k_12[i] * timestep / 2
        k_22[i] = forces[i] / bodies[i].mass
    end
    # calculate k_3 of 4th order Runge-Kutta method
    k_31 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    k_32 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    temporary_bodies = deepcopy(bodies) # calculate the intermediate state of the bodies
    for i in 1:num_bodies
        temporary_bodies[i].position += k_21[i] * timestep / 2
    end
    force_calculation!(forces, temporary_bodies) # calculate gravitational forces
    for i in 1:num_bodies
        k_31[i] = bodies[i].velocity + k_22[i] * timestep / 2
        k_32[i] = forces[i] / bodies[i].mass
    end
    # calculate k_4 of 4th order Runge-Kutta method
    k_41 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    k_42 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    temporary_bodies = deepcopy(bodies) # calculate the intermediate state of the bodies
    for i in 1:num_bodies
        temporary_bodies[i].position += k_31[i] * timestep
    end
    force_calculation!(forces, temporary_bodies) # calculate gravitational forces
    for i in 1:num_bodies
        k_41[i] = bodies[i].velocity + k_32[i] * timestep
        k_42[i] = forces[i] / bodies[i].mass
    end
    # update positions and velocities of all bodies
    for i in 1:num_bodies
        bodies[i].position += (k_11[i]/6 + k_21[i]/3 + k_31[i]/3 + k_41[i]/6) * timestep
        bodies[i].velocity += (k_12[i]/6 + k_22[i]/3 + k_32[i]/3 + k_42[i]/6) * timestep
    end
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

function Leapfrog!!(bodies::Vector{Body}, force_calculation!::Function)
    # use Leapfrog method

    forces = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies] # create a matrix to store forces
    force_calculation!(forces, bodies) # calculate gravitational forces
    for i in 1:num_bodies
        bodies[i].velocity += forces[i] / bodies[i].mass * timestep / 2
        bodies[i].position += bodies[i].velocity * timestep
    end
    force_calculation!(forces, bodies) # calculate gravitational forces again
    for i in 1:num_bodies
        bodies[i].velocity += forces[i] / bodies[i].mass * timestep / 2
    end
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