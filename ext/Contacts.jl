using MolSimToolkit: ContactMap
import PDBTools
import Plots: heatmap, cgrad

# the size of the plot should be proportional to the number 
# of residues in each dimension, with a 100 residues for 
# a 500px axis, by default
function _plot_size(map::ContactMap)
    nres1 = size(map.matrix, 1)
    nres2 = size(map.matrix, 2)
    xpixels = min(600, max(300, round(Int,500 * nres1 / 100)))
    ypixels = min(600, max(300, round(Int,500 * nres2 / 100)))
    return (xpixels, ypixels)
end

function heatmap(
    map::ContactMap{<:Real}; 
    xstep=max(1, div(size(map.matrix, 1), 20)), 
    ystep=max(1, div(size(map.matrix, 2), 20)),
    xticks=PDBTools.residue_ticks(map.residues1; stride=xstep, serial=true),
    yticks=PDBTools.residue_ticks(map.residues2; stride=ystep, serial=true),
    xrotation=60,
    xlabel="residue",
    ylabel="residue",
    colorbar_title="distance (Å)",
    aspect_ratio=1,
    xlims=(1,size(map.matrix, 1)),
    ylims=(1,size(map.matrix, 2)),
    color=:grayC,
    size=_plot_size(map),
    framestyle=:box,
    grid=false,
    clims=(1,2) .* extrema(skipmissing(map.matrix)),
    kargs...
)
    return heatmap(map.matrix; 
        xlabel, ylabel, xticks, yticks, xrotation,
        colorbar_title, color, aspect_ratio, xlims, ylims,
        size, framestyle, grid, clims,
        kargs...
    )
end

function heatmap(
    map::ContactMap{<:Bool}; 
    xstep=max(1, div(size(map.matrix, 1), 20)), 
    ystep=max(1, div(size(map.matrix, 2), 20)),
    xticks=PDBTools.residue_ticks(map.residues1; stride=xstep, serial=true),
    yticks=PDBTools.residue_ticks(map.residues2; stride=ystep, serial=true),
    xrotation=60,
    xlabel="residue",
    ylabel="residue",
    aspect_ratio=1,
    xlims=(1,size(map.matrix, 1)),
    ylims=(1,size(map.matrix, 2)),
    color=:Greys_9,
    size=_plot_size(map),
    framestyle=:box,
    grid=false,
    clims=(0,1),
    colorbar=:none,
    kargs...
)
    return heatmap(map.matrix; 
        xlabel, ylabel, xticks, yticks, xrotation,
        color, colorbar, aspect_ratio, xlims, ylims,
        size, framestyle, grid, clims,
        kargs...
    )
end
