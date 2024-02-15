import Plots: plot, plot!

function plot(
    perturb_range::Vector{Float64},
    results::Vector{ReweightResults};
    weight_cutoff::Union{Nothing,Vector{<:Real}} = nothing,
    labels=["perturbation", "frames fraction (%)", " "],
    color_palette::Symbol=:okabe_ito,
    dimensions::Vector{Int64}=[7, 10, 10, 5, 2]
    )
    if isnothing(weight_cutoff)
        weight_cutoff = zeros(length(results))
        weight_cutoff[1] = 1/(100*length(results[1]))
        for cut in 2:length(results)
            weight_cutoff[cut] = 5*(cut-1)/(100*length(results[1]))
        end
    end
    plt = plot()
    frames_fraction = zeros(length(perturb_range))
    for cut in weight_cutoff
        for i in eachindex(perturb_range, results)
            probs = results[i].probability
            frames_fraction[i] = count(number_of_frames -> (number_of_frames <= cut), probs) / length(probs)
        end
        plot!(plt, perturb_range, frames_fraction * 100, label="$(round(cut*100,sigdigits = 2))%")
    end
    plot!(plt,
        xlabel=labels[1],
        ylabel=labels[2],
        title=labels[3],
        legendfontsize=dimensions[1],
        titlefontsize=dimensions[2],
        guidefontsize=dimensions[3],
        minorticks=dimensions[4],
        linewidth=dimensions[5],
        palette=color_palette
    )
    return plt
end