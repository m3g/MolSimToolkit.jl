module REMDPlotting

using MolSimToolkit: GromacsREMDlog
import Plots: plot, plot!, annotate!, heatmap, cgrad, text

"""
    heatmap(data::GromacsREMDlog;
        xlabel="replica",
        ylabel="level",
        fontfamily="Computer Modern",
        probability_type::Symbol=:relative,
        kargs...
    )

Function to plot the probability matrix of a (H)REMD simulation. The probability matrix
contains the probability of finding each replica at each level of perturbation.

Use `probability_type=:relative` to plot the probability matrix normalized by the number
of replicas, or `probability_type=:absolute` to plot the probability matrix without normalization.

# Example

```julia-repl
julia> using MolSimToolkit, Plots

julia> data = remd_data(MolSimToolkit.gmx2019_4_log)

julia> heatmap(data)
```

"""
function heatmap(
    data::GromacsREMDlog;
    xlabel="replica",
    ylabel="level",
    fontfamily="Computer Modern",
    probability_type::Symbol=:relative,
    colorbar_title=nothing,
    kargs...
)
    nreplicas = size(data.probability_matrix, 1)
    if probability_type == :relative
        normalization = 1 / nreplicas
        digits = 1
        isnothing(colorbar_title) && (colorbar_title = "probability ÷ (1/$nreplicas)")
    elseif probability_type == :absolute
        normalization = 1
        digits = 2
        isnothing(colorbar_title) && (colorbar_title = "probability - (ideal = $(round(1/nreplicas; digits=3)))")
    else
        throw(ArgumentError("Invalid probability type."))
    end
    m = data.probability_matrix ./ normalization
    plt = heatmap(
        transpose(m);
        xlabel, ylabel, colorbar_title, fontfamily,
        clims=(0, maximum(m)),
        color=cgrad(:tempo),
        xticks=(1:nreplicas, 0:nreplicas-1),
        yticks=(1:nreplicas, 0:nreplicas-1),
        framestyle=:box,
    )
    for i in axes(m, 2), j in axes(m, 1)
        color = m[j, i] > 0.6 * maximum(m) ? :white : :black
        annotate!(plt, j, i,
            text("$(round(m[j, i]; digits))", :center, color, 7, "Computer Modern")
        )
    end
    return plt
end

end # module REMDPlotting