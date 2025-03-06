using MolSimToolkit: GromacsREMDlog
import Plots: plot, plot!, annotate!, heatmap, cgrad, text

"""
    heatmap(data::GromacsREMDlog;
        xlabel="replica",
        ylabel="initial replica level",
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
    xlabel="initial replica level",
    ylabel="level",
    fontfamily="Computer Modern",
    probability_type::Symbol=:relative,
    colorbar_title=nothing,
    kargs...
)
    nreplicas = size(data.probability_matrix, 1)
    if probability_type == :relative
        m = data.probability_matrix ./ (1 / nreplicas)
        digits = 1
        isnothing(colorbar_title) && (colorbar_title = "probability ÷ (1/$nreplicas)")
        clims = (0.5, 1.5)
        colorscale=:bwr
    elseif probability_type == :absolute
        m = data.probability_matrix
        digits = 2
        isnothing(colorbar_title) && (colorbar_title = "probability - (ideal = $(round(1/nreplicas; digits=3)))")
        clims=(0, maximum(m))
        colorscale=:tempo
    else
        throw(ArgumentError("Invalid probability type."))
    end

    plt = heatmap(
        transpose(m);
        xlabel, ylabel, colorbar_title, fontfamily,
        clims=clims,
        color=cgrad(colorscale),
        xticks=(1:nreplicas, 0:nreplicas-1),
        yticks=(1:nreplicas, 0:nreplicas-1),
        framestyle=:box,
        levels=12,
    )
    for i in axes(m, 2), j in axes(m, 1)
        color = if probability_type == :relative
            color = 1.2*clims[1] < m[j, i] < 0.8*clims[2] ? :black : :white
            color = :black
        elseif probability_type == :absolute
            color = m[j, i] > 0.6 * clims[2] ? :white : :black
        end
        annotate!(plt, j, i,
            text("$(round(m[j, i]; digits))", :center, color, 7, "Computer Modern")
        )
    end
    return plt
end

@testitem "REMD" begin
    using MolSimToolkit
    using Plots
    data = remd_data(MolSimToolkit.gmx2019_4_log)
    @test remd_replica_path(data, 0; stride=50) == [0, 0, 4, 5, 3, 1]
    @test remd_replica_path(data, 9; stride=50) == [9, 5, 1, 6, 6, 0]
    plt = heatmap(data; probability_type=:relative)
    tmp = tempname()*".png"
    savefig(plt, tmp)
    @test isfile(tmp)
    plt = heatmap(data; probability_type=:absolute)
    savefig(plt, tmp)
    @test isfile(tmp)
    data = remd_data(MolSimToolkit.gmx5_0_4_log)
    plt = heatmap(data; probability_type=:relative)
    savefig(plt, tmp)
    @test isfile(tmp)
    plt = heatmap(data; probability_type=:absolute)
    savefig(plt, tmp)
    @test isfile(tmp)
end