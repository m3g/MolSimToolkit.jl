import Plots: heatmap

function MolSimToolkit.ss_heatmap(ssmap::AbstractMatrix{<:Integer})
    plt = heatmap(MolSimStyle,ssmap,
        xlabel="frame",
        ylabel="residue",
        framestyle=:box,
        color=Plots.palette(:tab20c,10),
        clims=(0.5,10.5),
    )
    return plt
end