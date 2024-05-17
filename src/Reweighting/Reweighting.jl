module Reweighting

using CellListMap.PeriodicSystems
using ..MolSimToolkit: Simulation, positions, unitcell

export reweight, lennard_jones_perturbation, poly_decay_perturbation, gaussian_decay_perturbation

const testdir = "$(@__DIR__)/test"

include("./Reweight.jl")
include("./r_equations.jl")

end #Module Reweighting

@testitem "Reweighting one frame" begin
    import PDBTools
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_one_frame.xtc")

    i1 = PDBTools.selindex(atoms(simulation), "resname TFE and name O")

    i2 = PDBTools.selindex(atoms(simulation), "residue 11")

    sum_of_dist = reweight(simulation, (i,j,r) -> r, [i1[239]], i2; cutoff = 25.0)
    @test sum_of_dist.energy ≈ [7.4295543149]
end

@testitem "Reweighting small trajectory" begin
    import PDBTools
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    i1 = PDBTools.selindex(atoms(simulation), "index 97 or index 106")

    i2 = PDBTools.selindex(atoms(simulation), "residue 15 and name HB3")

    sum_of_dist = reweight(simulation, (i,j,r) -> r, i1, i2, cutoff = 25.0)
    @test sum_of_dist.energy ≈ [
        1.773896547670759, 1.5923698293115915, 1.716614676290554, 
        1.933003841107648, 1.602329229247863, 1.9639005665480983, 
        3.573986006775934, 2.188798265022823, 2.066180657974777, 
        1.6845109623700647
    ]
end