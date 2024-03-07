module REMDPlotting

using MolSimToolkit: GromacsREMDlog
import Plots: plot, plot!, annotate!, heatmap, cgrad, text

"""
    heatmap(data::GromacsREMDlog;
        xlabel="replica",
        ylabel="level",
        colorbar_title="probability ÷ (1/nreplicas)",
        fontfamily="Computer Modern",
        kargs...
    )

Function to plot the probability matrix of a (H)REMD simulation. The probability matrix
contains the probability of finding each replica at each level of perturbation.

# Example

```julia-repl
julia> using MolSimToolkit, Plots

julia> data = remd_data(MolSimToolkit.gmx2019_4_log)

julia> heatmap(data)
```

"""
function heatmap(data::GromacsREMDlog;
    xlabel="replica",
    ylabel="level",
    colorbar_title="probability ÷ (1/$(size(data.probability_matrix,1)))",
    fontfamily="Computer Modern",
    kargs...
)
    m = data.probability_matrix
    nreplicas = size(m, 1)
    plt = heatmap(
        transpose(m ./ (1 / nreplicas));
        xlabel, ylabel, colorbar_title, fontfamily,
        clims=(0, maximum(m ./ (1 / nreplicas))),
        color=cgrad(:tempo),
        xticks=(1:nreplicas, 0:nreplicas-1),
        yticks=(1:nreplicas, 0:nreplicas-1),
        framestyle=:box,
    )
    for i in axes(m, 2), j in axes(m, 1)
        color = m[j, i] > 0.2 ? :white : :black
        annotate!(plt, j, i,
            text("$(round(m[j, i] / (1/nreplicas),digits=1))", :center, color, 7, "Computer Modern")
        )
    end
    return plt
end

end # module REMDPlotting