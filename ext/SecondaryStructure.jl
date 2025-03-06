function MolSimToolkit.ss_heatmap(
    ssmap::AbstractMatrix{<:Integer}; 
    scalex=1.0, 
    kargs...
)
    default = Dict(
        :framestyle => :box,
        :color => Plots.palette(:tab20c,10),
        :clims => (0.5,10.5),
        :xlabel => "frame",
        :ylabel => "residue",
        :fontfamily => "Computer Modern",
        :adjust_latex_font => true,
    )
    plt = Plots.heatmap(MolSimStyle, scalex*(1:size(ssmap,2)), 1:size(ssmap,1), ssmap; _kargs(default; kargs)...)
    return plt
end

@testitem "ss_heatmap" begin
    using Plots
    using MolSimToolkit, MolSimToolkit.Testing
    using PDBTools
    simulation = Simulation(Testing.namd_pdb, Testing.namd_traj)
    ssmap = ss_map(simulation; ss_method=stride_run, show_progress=false)
    plt = ss_heatmap(ssmap; scalex=0.1, xlabel="time / ns")
    @test plt isa Plots.Plot{Plots.GRBackend}
    prot = select(atoms(simulation), isprotein)
    plt = ss_heatmap(ssmap; yticks=residue_ticks(prot; stride = 5))
    @test plt isa Plots.Plot{Plots.GRBackend}
end

