module Reweighting

using LinearAlgebra: diag
using CellListMap: ParticleSystem, map_pairwise, map_pairwise!
using ..MolSimToolkit: Simulation, positions, unitcell
using ..MolSimToolkit.MolecularMinimumDistances
import OrderedCollections
import PDBTools

export reweight, full_reweight, lennard_jones_perturbation, Perturbation, SystemPerturbations, gauss

const testdir = "$(@__DIR__)/test"

include("./Reweight.jl")
include("./perturbation_examples.jl")

end #Module Reweighting

@testitem "Reweighting one frame" begin
    import PDBTools
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_one_frame.xtc")

    i1 = at -> at.resname == "TFE" && (at.resnum == 239 && at.name == "O")

    i2 = "protein and residue 11"

    sum_of_dist = reweight(simulation, r -> r/10, i1, 1, i2, 1; all_dist = true, cutoff = 25.0)

    @test sum_of_dist.energy ≈ [7.4295543149]
end

@testitem "Reweighting small trajectory" begin
    import PDBTools
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    i1 = "index 97 or index 106"

    i2 = "residue 15 and name HB3"

    sum_of_dist = reweight(simulation, r -> r/10, i1, 1, i2, 1, all_dist = true, cutoff = 25.0)
    @test sum_of_dist.energy ≈ [
        1.773896547670759, 
        1.5923698293115915, 
        1.716614676290554, 
        1.933003841107648, 
        1.602329229247863, 
        1.9639005665480983, 
        3.573986006775934, 
        2.188798265022823, 
        2.066180657974777, 
        1.6845109623700647
    ]
end

@testitem "Reweighting small trajectory min dist" begin
    import PDBTools
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    i1 = "resname TFE"

    i2 = "protein"

    probs_test = reweight(simulation, r -> r, i1, 9, i2, 174)
    @test probs_test.energy ≈ [
        477.47530500203885
        459.8720850333111
        471.09562927166274
        484.763923006718
        436.1279250800145
        441.9556126997521
        759.6642496372428
        492.2758540693646
        361.8844625464858
        371.2874949349868
    ]
end

@testitem "Reweighting small trajectory min dist with contributions" begin
    import PDBTools
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    i1 = "resname TFE"

    i2 = "protein"

    probs_test = reweight(simulation, r -> r, i1, 9, i2, 174; mol_1_contrib = at -> at.resname  == "TFE", mol_2_contrib = at -> PDBTools.isprotein(at))
    @test probs_test.energy ≈ [
        477.47530500203885
        459.8720850333111
        471.09562927166274
        484.763923006718
        436.1279250800145
        441.9556126997521
        759.6642496372428
        492.2758540693646
        361.8844625464858
        371.2874949349868
    ]
end

@testitem "Reweighting small trajectory min dist with contributions 2" begin
    import PDBTools
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    i1 = "resname TFE"

    i2 = "protein"

    c1 = at -> at.resname == "TFE" && at.name == "H"

    c2 = at -> PDBTools.isprotein(at) && at.name == "O"

    probs_test = reweight(simulation, r -> r, i1, 9, i2, 174; mol_1_contrib = c1, mol_2_contrib = c2)
    @test probs_test.energy ≈ [
        0.0
        0.0
        8.239321724742876
        4.881649215140045
        0.0
        0.0
       18.21051179484232
       14.826747807218615
        0.0
        0.0      
    ]
end

@testitem "Reweighting small trajectory min dist with contributions 3" begin
    import PDBTools
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    i1 = "resname TFE"

    i2 = "resname SOL"

    c1 = at -> at.name == "H"

    c2 = at -> at.name in ["HW1", "HW2"]

    probs_test = reweight(simulation, r -> r, i1, 9, i2, 4; mol_1_contrib = c1, mol_2_contrib = c2)
    @test probs_test.energy ≈ [
        111083.78832237754
        111435.30788677695
        115309.31336585627
        109960.47285835315
        115367.46522696588
        114110.87534493855
        115304.83287464424
        111044.93873149024
        111135.20921166507
        114197.60245988662                 
    ]
end

@testitem "Reweighting small trajectory min dist with contributions 4" begin
    import PDBTools
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory_2.xtc")

    i1 = "residue 4981 and name F12"

    i2 = at -> at.residue == 4981 && at.name in ["H", "H21", "H22"]

    probs_test = reweight(simulation, r -> r, i1, 1, i2, 1)
    @test probs_test.energy ≈ [
        9.79988,
        8.608059,
        8.791174,
        9.190261999999999,
        9.154655,
        9.511442,
        9.502004,
        9.15598,
        9.294438,
        9.424879,   
    ] atol = 1.e-4
end

@testitem "Reweighting small trajectory min dist with contributions 5" begin
    import PDBTools
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory_2.xtc")

    i1 = "residue 4981 and name F12"

    i2 = at -> at.residue == 4981 && at.name in ["H", "H21", "H22"]

    c2 = at -> at.name == "H"

    probs_test = reweight(simulation, r -> r, i1, 1, i2, 1; mol_2_contrib = c2)
    @test probs_test.energy ≈ [
        4.348989,
        2.590000,
        2.735364,
        3.264019,
        4.087230,
        3.694575,
        4.272202,
        3.773232,
        3.166674,
        3.212258 
    ] atol = 1.e-4
end

@testitem "Reweighting small trajectory min dist with contributions 5" begin
    import PDBTools
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory_2.xtc")

    i1 = "residue 4981 and name F12"

    i2 = at -> at.residue == 4981 && at.name in ["H", "H21", "H22"]

    c2 = at -> at.name in ["H21", "H22"]

    probs_test = reweight(simulation, r -> r, i1, 1, i2, 1; mol_2_contrib = c2)
    @test probs_test.energy ≈ [
        5.450897,
        6.018059,
        6.05581,
        5.9262429999999995,
        5.067425,
        5.816867,
        5.229801999999999,
        5.382747999999999,
        6.127764,
        6.212621      
    ] atol = 1.e-4
end

@testitem "Reweighting small trajectory min dist with contributions using full reweight" begin
    import PDBTools
    import OrderedCollections
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory_2.xtc")

    g1 = PDBTools.selindex(simulation.atoms, "residue 4981 and name F12")

    g2 = PDBTools.selindex(simulation.atoms, at -> at.residue == 4981 && at.name in ["H", "H21", "H22"])

    c1 = at -> at = true

    c2 = at -> at.residue == 4981 && at.name in ["H21", "H22"]

    dist(r) = r

    Dict = OrderedCollections.OrderedDict(1 => Perturbation(simulation.atoms, c1, c2, dist))

    input = SystemPerturbations(g1, 1, g2, 1, Dict)

    res = full_reweight(simulation,
                    input;
                    k = 1.0,
                    T = 1.0,
                    cutoff = 12.0,
            )
    @test res[1].energy ≈ [
        5.450897,
        6.018059,
        6.05581,
        5.92624299,
        5.067425,
        5.816867,
        5.22980199,
        5.38274799,
        6.127764,
        6.212621      
    ] atol = 1.e-5
end