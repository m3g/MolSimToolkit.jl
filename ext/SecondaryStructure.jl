import Plots: heatmap

function ss_heatmap(
    ssmap,
    selection,
)
    selection = at -> PDBTools.isprotein(at) && Select(selection)(at)
    plt = heatmap(ssmap,
        xlabel="frame",
        ylabel="residue",
        framestyle=:box,
        color=palette(:tab20c,10)
    )
    return plt
end