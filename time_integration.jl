# Some time integration techniques to update the position and velocity of the bodies in 1 time step.
# Last edited by Cha0s_MnK on 2023/09/05.

# reference:
# Star Formation in Galaxy Evolution: Connecting Numerical Models to Reality
# (Nickolay Y. Gnedin, Simon C.O. Glover, Ralf S. Klessen, Volker Springel)
# High Performance Computing and Numerical Modelling (Volker Springel)
# 3 Time Integration Techniques

function explicitEuler!(bodies::Vector{Body}, force_calculation!::Function)
    # use explicit Euler method

    # create a matrix to store forces
    forces = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    # calculate gravitational forces
    force_calculation!(forces, bodies)
    # update positions and velocities of all bodies
    for i in 1:num_bodies
        bodies[i].position += bodies[i].velocity * time_step
        bodies[i].velocity += forces[i] / bodies[i].mass * time_step
    end
end

function RungeKutta2nd!(bodies::Vector{Body}, force_calculation!::Function)
    # use 2nd order Runge-Kutta method

    # create a matrix to store forces
    forces = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    # calculate k_1 of 2nd order Runge-Kutta
    k_11 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    k_12 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    force_calculation!(forces, bodies) # calculate gravitational forces
    for i in 1:num_bodies
        k_11[i] = bodies[i].velocity
        k_12[i] = forces[i] / bodies[i].mass
    end
    # calculate k_2 of 2nd order Runge-Kutta
    k_21 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    k_22 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    temporary_bodies = deepcopy(bodies) # calculate the intermediate state of the bodies
    for i in 1:num_bodies
        temporary_bodies[i].position += temporary_bodies[i].velocity * time_step
    end
    force_calculation!(forces, temporary_bodies) # calculate gravitational forces
    for i in 1:num_bodies
        k_21[i] = bodies[i].velocity + k_11[i] * time_step
        k_22[i] = forces[i] / bodies[i].mass
    end
    # update positions and velocities of all bodies
    for i in 1:num_bodies
        bodies[i].position += (k_11[i] + k_21[i]) * time_step / 2
        bodies[i].velocity += (k_12[i] + k_22[i]) * time_step / 2
    end
end

function RungeKutta4th!(bodies::Vector{Body}, force_calculation::Function)
    # use 4th order Runge-Kutta

    # calculate k_1 of 4th order Runge-Kutta
    k_11 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    k_12 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    forces = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    for i in 1:num_bodies
        force_calculation(forces[i], bodies, bodies[i])
        k_11[i] .= bodies[i].velocity
        k_12[i] .= forces[i] ./ bodies[i].mass
    end
    # calculate k_2 of 2nd order Runge-Kutta
    k_21 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    k_22 = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    forces = [MVector{3, Float64}(0.0,0.0,0.0) for _ in 1:num_bodies]
    temporary_bodies = deepcopy(bodies) # calculate the intermediate state of the bodies
    for i in 1:num_bodies
        temporary_bodies[i].position .+= temporary_bodies[i].velocity .* time_step
    end
    for i in 1:num_bodies
        force_calculation(forces[i], temporary_bodies, temporary_bodies[i])
        k_21[i] .= bodies[i].velocity .+ k_11[i] .* time_step
        k_22[i] .= forces[i] ./ bodies[i].mass
    end
    # update positions and velocities of all bodies
    for i in 1:num_bodies
        bodies[i].position .+= (k_11[i] + k_21[i]) .* time_step ./ 2
        bodies[i].velocity .+= (k_12[i] + k_22[i]) .* time_step ./ 2
    end
end

function Leapfrog(body::Body,)
    intermediate_velocity = body.velocity .+ force ./ body.mass .* time_step ./ 2
    body.position .+= intermediate_velocity .* time_step

end