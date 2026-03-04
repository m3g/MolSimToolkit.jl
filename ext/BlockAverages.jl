using LaTeXStrings
import Plots: histogram, plot
using Plots: plot!, histogram!, hline!, @layout
using Plots.Measures: cm

_round(x::Real; kargs...) = round(x; kargs...)
_round(x; kargs...) = round(typeof(x), x; kargs...)

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
    plot!(p, x -> (1 / sqrt(π * σ)) * exp(-(oneunit(md.mean) * x - md.mean)^2 / σ), linewidth=2, color=:black, label=:none)
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
    tu = data.dt / oneunit(data.dt)
    l = @layout [ a{0.2h} ; b c ; d e]
    p = plot(layout=l)
    plot!(subplot=1, 
        collect(1:length(data.x)) * data.dt,
        data.x,
        xlabel="time",
        ylabel="value",
        label="",
        color=:black,
        title="\n$title",
        topmargin=0.3cm,
    )
    hline!(
        [data.xmean],
        linestyle=:dash,
        alpha=0.5,
        color=:black,
        label=:none,
        subplot=2,
    )
    plot!(
        data.dt * data.blocksize, data.xmean_maxerr,
        ylabel="worst block value",
        xlabel="block size",
        label=nothing,
        linewidth=2,
        marker=:circle,
        color=:black,
        xscale=xscale,
        subplot=2,
        legendtitle="mean = $(_round(data.xmean, digits=2))",
    )
    plot!(data.dt * data.blocksize, data.xmean_stderr,
        ylabel=L"SD / \sqrt{N_{blocks}}",
        xlabel="block size",
        label=nothing,
        linewidth=2,
        marker=:circle,
        color=:black,
        xscale=xscale,
        subplot=3
    )
    # Auto correlation function
    plot!(
        data.lags * data.dt,
        data.autocor,
        ylabel=L"c(\Delta t)",
        xlabel=L"\Delta t",
        label=nothing,
        linewidth=2,
        color=:black,
        subplot=4
    )
    t95 = 1.96 / sqrt(length(data.x))
    i95 = findfirst(i -> data.autocor[i] <= t95, eachindex(data.lags))
    isnothing(i95) && (i95 = length(data.lags))
    i95 -= 1
    hline!([data.autocor[i95]], subplot=4, ls=:dash, label="", color=:grey)
    exp_fit = exp.(-inv((data.tau/oneunit(data.tau))) .* tu .* data.lags )
    plot!(
        data.lags * data.dt,
        exp_fit,
        label=nothing,
        linewidth=2,
        color=:black,
        alpha=0.5,
        subplot=4,
        legendtitle="τ = $(_round(data.tau; digits=2))",
    )
    plot!(
        p,
        size=(800, 600),
        framestyle=:box,
        fontfamily="Computer Modern",
        xlims=xlims,
        ylims=ylims,
        leftmargin=0.3cm,
        rightmargin=0.3cm,
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
    plot!((1,1), subplot=5, lc=:white, label="\n"*latexstring("\\textrm{\\Delta t (0.95) = $(_round(data.lags[i95] * data.dt; digits=4))}"))
    plot!((1,1), subplot=5, lc=:white, label=latexstring("\\textrm{Integrated-\\tau = $(_round(data.tau_int; digits=4))}~~~~"))
    plot!((1,1), subplot=5, lc=:white, label=latexstring("\\textrm{N = $(length(data.x))}"))
    plot!((1,1), subplot=5, lc=:white, label=latexstring("\\textrm{N_{eff}  = $(round(Int, data.n_effective))}"))
    plot!((1,1), subplot=5, lc=:white, label=latexstring("\\textrm{SEM(N_{eff}) = $(_round(data.xmean_stderr_neff; digits=4))}"))
    return p
end

@testitem "blockaverages plotting" begin
    using MolSimToolkit, Plots, Unitful
    ENV["GKSwstype"] = "nul"

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

    x = x .* 1u"cm"
    md = block_distribution(x; block_size=10^5);
    tempfile = tempname()*".png"
    plt = histogram(md)
    savefig(plt, tempfile)
    @test isfile(tempfile)
    b = block_average(x, lags=0:100:10^5, dt=1u"s")
    plt = plot(b)
    savefig(plt, tempfile)
    @test isfile(tempfile)

end