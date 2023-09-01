# Some time integration functions to update the position and velocity of a body in 1 time step
# Last edited by Cha0s_MnK on 2023/09/01.

# reference:
# Star Formation in Galaxy Evolution: Connecting Numerical Models to Reality
# (Nickolay Y. Gnedin, Simon C.O. Glover, Ralf S. Klessen, Volker Springel)
# High Performance Computing and Numerical Modelling (Volker Springel)
# 3 Time Integration Techniques

function periodic_boundary(body::Body)
    # apply periodic boundary condition to modify coordinates

    for i in 1:3
        body.position[i] = mod(body.position[i] + total_length/2, total_length) - total_length/2
    end
end

function explicitEuler(body::Body, force::MVector{3, Float64})
    # use explicit Euler method

    body.position .+= body.velocity * time_step
    body.velocity .+= (force ./ body.mass) * time_step
end