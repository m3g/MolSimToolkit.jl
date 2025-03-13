using MolSimToolkit: MapFractionsResult, MDLovoFitResult
using PDBTools: residue_ticks, name, isprotein
import Plots: plot, plot!

function Plots.plot(
    mf::MapFractionsResult;
    xlabel="fraction of aligned residues",
    ylabel="RMSD",
    kargs...
) 
    plt = plot(MolSimStyle)
    plt = plot!(
        plt,
        mf.fraction,
        [mf.rmsd_low, mf.rmsd_high, mf.rmsd_all];
        label=[ "best"  "worse"  "all" ],
        xlabel, ylabel,
        lw=2,
        kargs...
    )
    hline!(plt, [1.0 2.0 3.0]; label="", color=:black, linestyle=:dash, linewidth=1, alpha=0.2)
    plot!(plt, 
        xlims=(0.2, 1.0),
        yticks=(0:1:maximum(mf.rmsd_high),0:1:maximum(mf.rmsd_high)),
    )
    return plt
end

function Plots.plot(
    md::MDLovoFitResult;
    xlabel=["frame"  "residue"],
    ylabel=["RMSD / Å" "RMSF / Å"],
    yticks=[(0:1:maximum(md.rmsd_high),0:1:maximum(md.rmsd_high)) (0:1:maximum(md.rmsf),0:1:maximum(md.rmsf))],
    ylims=[(0, maximum(md.rmsd_high)) (0, maximum(md.rmsf))],
    kargs...
) 
    plt = plot(MolSimStyle, layout=(2, 1))
    plot!(
        plt,
        md.iframes,
        [md.rmsd_low, md.rmsd_high, md.rmsd_all];
        label=[ "best"  "worse"  "all" ],
        lw=2,
        subplot=1,
        kargs...
    )
    hline!(plt, [1.0 2.0 3.0]; label="", color=:black, linestyle=:dash, linewidth=1, alpha=0.2, subplot=1)
    xticks2 = residue_ticks(
                filter(at -> isprotein(at) && name(at) == "CA", atoms(md.simulation));
                serial=true,
                stride=max(1, div(length(md.rmsf), 20)),
    )
    plot!(
        plt,
        1:length(md.rmsf),
        md.rmsf,
        lw=2,
        subplot=2,
        xticks=xticks2,
        xrotation=60,
        label="",
        kargs...
    )
    plot!(plt; 
        size=(600,600),
        margin=0.1Plots.Measures.cm,
        leftmargin=0.3Plots.Measures.cm,
        ylims, yticks,
        xlabel, ylabel,

    )
    hline!(plt, [1.0 2.0 3.0]; label="", color=:black, linestyle=:dash, linewidth=1, alpha=0.2, subplot=2)
    return plt
end
