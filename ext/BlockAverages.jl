using LaTeXStrings
import Plots: histogram, plot
using Plots: plot!, histogram!, hline!, annotate!, text, @layout
using Plots.Measures: cm

"""
    histogram(md::BlockDistribution; bins=:auto)

Function that creates a histogram of the value of a property in the blocks, computed
with the `block_distribution` function.  

# Example

```julia-repl
julia> using MolSimToolkit, Plots

julia> x = BlockAverages.test_data(10^6);

julia> md = block_distribution(x; block_size=10^5);

julia> histogram(md)
```
"""
function histogram(md::BlockDistribution; bins=:auto)
    p = plot(
        size=(400, 300),
        framestyle=:box,
        fontfamily="Computer Modern",
        leftmargin=0.5cm,
        rightmargin=0.5cm,
        xlabel="value",
        ylabel="density"
    )
    histogram!(p, md.block_mean, bins=bins, normalize=:pdf, color=:gray, label=:none)
    σ = md.std_of_the_mean^2
    plot!(p, x -> (1 / sqrt(π * σ)) * exp(-(x - md.mean)^2 / σ), linewidth=2, color=:black, label=:none)
    return p
end

"""
    plot(
        data::BlockAverageData; 
        xlims=nothing, 
        ylims=nothing,
        xscale=:identity,
        title="",
    )

Function that creates a plot of output data from `block_average`. The optional 
`xlims`, `ylims`, `xscale`, `title`, can be set to adjust the apparency of the plot.    

## Example

```julia-repl
julia> using MolSimToolkit, Plots

julia> x = BlockAverages.test_data(10^6);

julia> b = block_average(x, lags=0:100:10^5);

julia> plot(b)
```

"""
function plot(
    data::BlockAverageData;
    xlims=:auto,
    ylims=:auto,
    xscale=:identity,
    title=""
)
    l = @layout [ a{0.2h} ; b c ; d e]
    p = plot(layout=l)
    plot!(subplot=1, 
        data.x,
        xlabel="step",
        ylabel="value",
        label="",
        color=:black,
    )
    hline!(
        [data.xmean],
        linestyle=:dash,
        alpha=0.5,
        color=:black,
        label=:none,
        title=title,
        subplot=2,
    )
    plot!(
        data.blocksize, data.xmean_maxerr,
        ylabel="worst block value",
        xlabel=L"\textrm{block~size~}(N)",
        label=nothing,
        linewidth=2,
        marker=:circle,
        color=:black,
        xscale=xscale,
        subplot=2
    )
    annotate!(
        maximum(data.blocksize) - 0.1 * maximum(data.blocksize),
        max(data.xmean_maxerr[end], maximum(data.xmean_maxerr)) - 0.1 * (maximum(data.xmean_maxerr) - minimum(data.xmean_maxerr)),
        text("mean = $(round(data.xmean, digits=2))", "Computer Modern", 12, :right),
        subplot=2,
    )
    plot!(data.blocksize, data.xmean_stderr,
        ylabel=L"\sigma^2 / \sqrt{N}",
        xlabel=L"\textrm{block~size~}(N)",
        label=nothing,
        linewidth=2,
        marker=:circle,
        color=:black,
        xscale=xscale,
        subplot=3
    )
    # Auto correlation function
    plot!(
        data.lags,
        data.autocor,
        ylabel=L"c(\Delta t)",
        xlabel=L"\Delta t",
        label=nothing,
        linewidth=2,
        color=:black,
        subplot=4
    )
    t95 = 1.96 / sqrt(length(data.x))
    hline!([t95], subplot=4, ls=:dash, label="", color=:grey)
    exp_fit = exp.(-inv(data.tau) .* data.lags)
    plot!(
        data.lags,
        exp_fit,
        label=nothing,
        linewidth=2,
        color=:black,
        alpha=0.5,
        subplot=4,
    )
    annotate!(
        data.lags[end] - 0.2 * data.lags[end],
        0.8 * max(maximum(data.autocor), maximum(exp_fit)),
        text("τ = $(round(data.tau, digits=2))", "Computer Modern", 12, :right),
        subplot=4,
    )
    plot!(
        p,
        size=(800, 600),
        framestyle=:box,
        fontfamily="Computer Modern",
        xlims=xlims,
        ylims=ylims,
        leftmargin=0.1cm,
        rightmargin=0.1cm,
    )
    plot!(subplot=5,
        xticks=nothing,
        yticks=nothing,
        grid=false,
        framestyle=:none,
        legend=:topleft,
        legendborder=nothing,
        legendtitle="Summary",
        legendfontsize=8,
    )
    i95 = findfirst(i -> data.autocor[i] <= t95, eachindex(data.lags))
    isnothing(i95) && (i95 = length(data.lags))
    i95 -= 1
    plot!((1,1), subplot=5, lc=:white, label="\n"*latexstring("\\Delta t (0.95) = $(data.lags[i95])"))
    plot!((1,1), subplot=5, lc=:white, label=latexstring("\\textrm{Integrated-}\\tau  = $(round(data.tau_int; digits=4))"))
    plot!((1,1), subplot=5, lc=:white, label=latexstring("N = $(length(data.x))"))
    plot!((1,1), subplot=5, lc=:white, label=latexstring("N_{eff}  = $(round(Int, data.n_effective))"))
    plot!((1,1), subplot=5, lc=:white, label=latexstring("SEM(N_{eff}) = $(round(data.xmean_stderr_neff; digits=4))"))
    return p
end

@testitem "blockaverages plotting" begin
    using MolSimToolkit, Plots
    x = BlockAverages.test_data(10^6);
    md = block_distribution(x; block_size=10^5);
    tempfile = tempname()*".png"
    plt = histogram(md)
    savefig(plt, tempfile)
    @test isfile(tempfile)
    b = block_average(x, lags=0:100:10^5)
    plt = plot(b)
    savefig(plt, tempfile)
    @test isfile(tempfile)
end