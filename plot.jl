# This is the file about plots in a N-body simulator.
# Last edited by Chaos on 2023/08/30.

function plotBodies!(bodies::Vector{Body}, title0::AbstractString)
    # plot current bodies on a canvas; return the canvas

    # create a beautiful 3D canvas
    plot_length = 0.6*total_length # "0.6*" to make it more beautiful
    canvas = plot3d(
        size = (1080, 720),
        legend = false, legendfontfamily = "Times New Roman", legendfontsize = 10,
        xlims = (-plot_length, plot_length), ylims = (-plot_length, plot_length), zlims = (-plot_length, plot_length),
        title = title0, titlefontfamily = "Times New Roman", titlefontsize = 12,
        xticks = true, yticks = true, zticks = true, tickfontfamily = "Times New Roman", tickfontsize = 8,
        xlabel = "X", xlabelfontfamily = "Times New Roman", xlabelfontsize = 10, 
        ylabel = "Y", ylabelfontfamily = "Times New Roman", ylabelfontsize = 10,
        zlabel = "Z", zlabelfontfamily = "Times New Roman", zlabelfontsize = 10,
        background_color = RGB(0.95, 0.95, 0.95), 
        fg_legend = :white,
        linealpha = 0.7, linewidth = 1, linestyle = :auto, 
        marker = (:circle, 4),
        grid = true,
        color = :viridis,
        camera = (30, 30), elevation = 30, azimuth = 30
    )
    # plot current bodies on the canvas
    for body in bodies
        scatter!(canvas, [body.position[1]], [body.position[2]], [body.position[3]], markersize=3, color=:blue)
    end
    # return the canvas
    return canvas
end

function get_distances(i::Int, distances::Matrix{Float64}, bodies1::Vector{Body}, bodies2::Vector{Body})
    # calculate the distances; fix them; store them

    for j in 1:num_bodies
        distances[i, j] = norm(bodies2[j].position - bodies1[j].position)
        if distances[i, j] > 0.5 * total_length # effects caused by periodic boundary conditions
            temporary_position1 = bodies1[j].position
            temporary_position2 = bodies2[j].position
            # calculate the difference
            temporary_difference = temporary_position2 - temporary_position1
            # adjust positions based on periodic boundary conditions
            temporary_position1 .+= (temporary_difference .> 0.5*total_length) .* total_length
            temporary_position2 .+= (temporary_difference .< -0.5*total_length) .* total_length
            distances[i, j] = norm(temporary_position2 - temporary_position1)
        end
    end
end

function plotCurves(distances::Matrix{Float64}, PNGname::AbstractString)
    # plot curves for the time evolution of each body; save the plot as a PNG file

    # create an empty plot for curves
    curves = plot(
        size = (1080, 720),
        title = "Mass = $Mean_mass solar mass \n Side length = $TotalLength pc \n Total time = $TotalTime Julian year \n Time step = $TimeStep Julian year \n Opening angle = $theta", titlefontfamily = "Times New Roman", titlefontsize = 12,
        legendfontfamily = "Times New Roman", legendfontsize = 8,
        tickfontfamily = "Times New Roman", tickfontsize = 8,
        xlabel = "Time (Julian year)", xlabelfontfamily = "Times New Roman", xlabelfontsize = 10, 
        ylabel = "Distances (pc)", ylabelfontfamily = "Times New Roman", ylabelfontsize = 10
    )
    # add curves for the time evolution of each body
    for j in 1:num_bodies
        # change the unit of the variables
        # Time: s --> Julian year; Distances: m --> pc
        plot!(curves, [time_step.*(1:num_steps)./Julian_year], distances[:, j]./pc, label="Body $j")
    end
    # save the plot as a PNG file
    savefig(curves, PNGname)
end

function plot_edges(cell::Cell)
    # preorder traversal the BHoctree and plot edges of non-empty leaf cells 
    if cell.subcells != nothing # this is a branch cell
        for subcell in cell.subcells
            plot_edges(subcell)
        end
    else # this is a leaf cell (after cell update, it can not be empty)
        # calculate the corners of the cell
        cell_length = total_length / 2^cell.depth
        corners = [cell.center .+ cell_length / 2 .* MVector{3, Float64}(x, y, z) for x in [-1.0, 1.0] for y in [-1.0, 1.0] for z in [-1.0, 1.0]]

        # plot edges of the cell
        edges = [(1, 2), (1, 3), (2, 4), (3, 4), (5, 6), (5, 7), (6, 8), (7, 8), (1, 5), (2, 6), (3, 7), (4, 8)]
        for edge in edges
            i, j = edge
            plot!([corners[i][1], corners[j][1]],[corners[i][2], corners[j][2]],[corners[i][3], corners[j][3]], color=:blue, linewidth=line_width)
        end
    end
end