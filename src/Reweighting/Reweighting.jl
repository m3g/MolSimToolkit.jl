module Reweighting

using CellListMap.PeriodicSystems
using ..MolSimToolkit: Simulation, positions, unitcell

export reweight, L_J, poly_decay, gaussian_decay

const testdir = "$(@__DIR__)/test"

include("./reweight.jl")
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
    @test sum_of_dist.energy ≈ [74.295543]
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
        17.73896547670759, 15.923698293115915, 17.16614676290554, 
        19.33003841107648, 16.02329229247863, 19.639005665480983, 
        35.73986006775934, 21.88798265022823, 20.66180657974777, 
        16.845109623700647
    ]
end