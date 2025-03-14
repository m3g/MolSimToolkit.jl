using MolSimToolkit: MapFractionsResult, MDLovoFitResult
using PDBTools: residue_ticks, name, isprotein
import Plots: plot, plot!

"""
    plot(mf::MapFractionsResult; xlabel="fraction of aligned residues", ylabel="RMSD", kargs...)

Plot the results of a map fractions calculation.

"""
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

"""
    plot(md::MDLovoFitResult; 
        xlabel=["frame"  "residue"], 
        ylabel=["RMSD / Å" "RMSF / Å"], 
        yticks=[(0:1:maximum(md.rmsd_high),0:1:maximum(md.rmsd_high)) (0:1:maximum(md.rmsf),0:1:maximum(md.rmsf))], 
        ylims=[(0, maximum(md.rmsd_high)) (0, maximum(md.rmsf))], 
        stride=nothing, 
        kargs...
    )

Plot the results of a MDLovoFit calculation. The `stride` argument can be used to control the number of xtick residues plotted in the RMSF plot.

"""
function Plots.plot(
    md::MDLovoFitResult;
    xlabel=["frame"  "residue"],
    ylabel=["RMSD / Å" "RMSF / Å"],
    yticks=[(0:1:maximum(md.rmsd_high),0:1:maximum(md.rmsd_high)) (0:1:maximum(md.rmsf),0:1:maximum(md.rmsf))],
    ylims=[(0, maximum(md.rmsd_high)) (0, maximum(md.rmsf))],
    stride=nothing,
    kargs...
) 
    plt = plot(MolSimStyle, layout=(2, 1))
    plot!(
        plt,
        md.frame_indices,
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
                stride= isnothing(stride) ? max(1, div(length(md.rmsf), 20)) : stride,
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

@testitem "mdlovofit - plots" begin
    using MolSimToolkit, PDBTools, Plots
    sim = Simulation(
        MolSimToolkit.Testing.namd_pdb,
        MolSimToolkit.Testing.namd_traj,
    )
    mf = map_fractions(sim)
    tmpfile = tempname()*".png"
    plt = plot(mf)
    savefig(plt, tmpfile)
    @test isfile(tmpfile)
    md = mdlovofit(sim; fraction=0.7)
    plt = plot(md)
    savefig(plt, tmpfile)
    @test isfile(tmpfile)
    plt = plot(md; stride=30)
    savefig(plt, tmpfile)
    @test isfile(tmpfile)
end