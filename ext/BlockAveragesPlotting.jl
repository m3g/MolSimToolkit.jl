using LaTeXStrings
import Plots: plot, plot!, histogram!, hline!, annotate!, text, histogram
import Plots.Measures: cm

"""
    histogram(md::BlockDistribution; bins=:auto)

Function that creates a histogram of the value of a property in the blocks, computed
with the `block_distribution` function.  

# Example

```julia-repl
julia> using MolSimToolKit, Plots

julia> x = BlockAverage.test_data(10^6);

julia> md = block_distribution(x; nblocks=1000);

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
julia> using MolSimToolKit, Plots

julia> x = BlockAverage.test_data(10^6);

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
    p = plot(layout=(3, 1))
    hline!(
        [data.xmean],
        linestyle=:dash,
        alpha=0.5,
        color=:black,
        label=:none,
        title=title,
        subplot=1,
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
        subplot=1
    )
    annotate!(
        maximum(data.blocksize) - 0.1 * maximum(data.blocksize),
        max(data.xmean_maxerr[end], maximum(data.xmean_maxerr)) - 0.1 * (maximum(data.xmean_maxerr) - minimum(data.xmean_maxerr)),
        text("mean = $(round(data.xmean, digits=2))", "Computer Modern", 12, :right),
        subplot=1,
    )
    plot!(data.blocksize, data.xmean_stderr,
        ylabel=L"\sigma^2 / \sqrt{N}",
        xlabel=L"\textrm{block~size~}(N)",
        label=nothing,
        linewidth=2,
        marker=:circle,
        color=:black,
        xscale=xscale,
        subplot=2
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
        subplot=3
    )
    exp_fit = exp.(-inv(data.tau) .* data.lags)
    plot!(
        data.lags,
        exp_fit,
        label=nothing,
        linewidth=2,
        color=:black,
        alpha=0.5,
        subplot=3,
    )
    annotate!(
        data.lags[end] - 0.2 * data.lags[end],
        0.8 * max(maximum(data.autocor), maximum(exp_fit)),
        text("τ = $(round(data.tau, digits=2))", "Computer Modern", 12, :right),
        subplot=3,
    )
    plot!(
        p,
        size=(400, 600),
        framestyle=:box,
        fontfamily="Computer Modern",
        xlims=xlims,
        ylims=ylims,
        leftmargin=0.5cm,
        rightmargin=0.5cm,
    )
    return p
end