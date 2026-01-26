module Reweighting

using ProgressMeter
using LinearAlgebra: diag
using CellListMap
using ..MolSimToolkit: Simulation, positions, unitcell
using ..MolSimToolkit.MolecularMinimumDistances
import OrderedCollections
import PDBTools

export reweight, lennard_jones_perturbation, Perturbation, SystemPerturbations, SystemPerturbationsOneGroup, gauss

const testdir = "$(@__DIR__)/test"

include("./Reweight.jl")
include("./perturbation_examples.jl")

end #Module Reweighting

@testitem "Reweight with small trajectory using minimum distances and atoms from the same molecule" begin
    import PDBTools
    import OrderedCollections
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    g1 = PDBTools.selindex(simulation.atoms, "residue 4981 and name F12")

    g2 = PDBTools.selindex(simulation.atoms, at -> at.residue == 4981 && at.name in ["H", "H21", "H22"])

    c1 = at -> at = true

    c2 = at -> at.residue == 4981 && at.name in ["H21", "H22"]

    dist(r) = r

    Dict = OrderedCollections.OrderedDict(1 => Perturbation(simulation.atoms, c1, c2, dist, [1]))

    input = SystemPerturbations(g1, 1, g2, 3, Dict)

    res = reweight(simulation,
                    input;
                    all_distances = false,
                    k = 1.0,
                    T = 1.0,
                    cutoff = 12.0,
            )

    @test res[1].energies[1] ≈ [
        2.66700,
        0.0,
        0.0,
        2.64641,
        2.42596,
        2.61084,
        2.40289,
        2.48584,
        2.86653,
        2.86769      
    ] atol = 1.e-4

    @test res[1].distances ≈ [
        1,
        0,
        0,
        1,
        1,
        1,
        1,
        1,
        1,
        1,  
    ] atol = 1.e-4
end

@testitem "Reweight with small trajectory using all distances and atoms from the same molecule" begin
    import PDBTools
    import OrderedCollections
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    g1 = PDBTools.selindex(simulation.atoms, "residue 4981 and name F12")

    g2 = PDBTools.selindex(simulation.atoms, at -> at.residue == 4981 && at.name in ["H", "H21", "H22"])

    c1 = at -> at = true

    c2 = at -> at.residue == 4981 && at.name in ["H21", "H22"]

    dist(r) = r

    Dict = OrderedCollections.OrderedDict(1 => Perturbation(simulation.atoms, c1, c2, dist, [-2, -1, 0, 1, 2]))

    input = SystemPerturbations(g1, 1, g2, 3, Dict)

    res = reweight(simulation,
                    input;
                    all_distances = true,
                    k = 1.0,
                    T = 1.0,
                    cutoff = 12.0,
            )
    @test res[1].energies ≈ [
        -2*[5.450897,6.018059,6.05581,5.926243,5.067425,5.816867,5.229802,5.382748,6.127764, 6.212621],
        -[5.450897,6.018059,6.05581,5.926243,5.067425,5.816867,5.229802,5.382748,6.127764, 6.212621],
        zeros(10),
        [5.450897,6.018059,6.05581,5.926243,5.067425,5.816867,5.229802,5.382748,6.127764, 6.212621],
        2*[5.450897,6.018059,6.05581,5.926243,5.067425,5.816867,5.229802,5.382748,6.127764, 6.212621],
    ] atol = 1.e-4
end

@testitem "Reweight with small trajectory using minimum distances and contributions from different residues" begin
    import PDBTools
    import OrderedCollections
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    g1 = PDBTools.selindex(simulation.atoms, at -> at.residue == 15 && at.name in ["H", "N", "C", "OC1", "OC2"])

    g2 = PDBTools.selindex(simulation.atoms,at -> at.residue == 11 && at.name in ["CB", "HB1", "HB2", "HB3"])

    c1 = at -> at.name in ["H"]

    c2 = at -> at.name in ["HB3"]

    c11 = at -> at.name in ["H"]

    c12 = at -> at = true

    dist(r) = r

    Dict = OrderedCollections.OrderedDict(
        "a" => Perturbation(simulation.atoms, c1, c2, dist, [-1, 0, 1]), 
        "b" => Perturbation(simulation.atoms, c11, c12, dist, [1.0])
    )

    input = SystemPerturbations(g1, 5, g2, 4, Dict)

    res = reweight(simulation,
                    input;
                    all_distances = false,
                    k = 1.0,
                    T = 1.0,
                    cutoff = 15.0,
            )

    @test res["a"].energies ≈ [
        -[6.568418, 0, 6.064767, 0, 0, 0, 9.888358, 0, 7.498538, 0],
        zeros(10),
        [6.568418, 0, 6.064767, 0, 0, 0, 9.888358, 0, 7.498538, 0],
    ] atol = 1.e-3

    @test res["a"].probabilities ≈ [
        [0.031440, 4.41428e-5, 0.0190000, 4.41428e-5, 4.41428e-5, 4.41428e-5, 0.869599, 4.41428e-5, 0.079695, 4.41428e-5],
        ones(10)/10,
        [0.000234, 0.166546, 0.000387, 0.166546, 0.166546, 0.166546, 8.4542e-6, 0.166546, 9.224899e-5, 0.166546],
    ] atol = 1.e-5

    @test res["a"].distances ≈ [
        1, 0, 1, 0, 0, 0, 1, 0, 1, 0,
    ] atol = 1.e-1

    @test res["b"].energies[1] ≈ [
        6.568418,
        5.697282,
        6.064767,
        0,
        6.03079,
        6.464387,
        9.888358,
        9.800954,
        7.498538,
        6.192204      
    ] atol = 1.e-4
end

@testitem "Reweight with small trajectory using all distances and contributions from different residues" begin
    import PDBTools
    import OrderedCollections
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    g1 = PDBTools.selindex(simulation.atoms, at -> at.residue == 15 && at.name in ["H", "N", "C", "OC1", "OC2"])

    g2 = PDBTools.selindex(simulation.atoms,at -> at.residue == 11 && at.name in ["CB", "HB1", "HB2", "HB3"])

    c1 = at -> at.name in ["H"]

    c2 = at -> at.name in ["HB3"]

    c11 = at -> at.name in ["H"]

    c12 = at -> at = true

    dist(r) = r

    Dict = OrderedCollections.OrderedDict(
        "a" => Perturbation(simulation.atoms, c1, c2, dist, [1]), 
        "b" => Perturbation(simulation.atoms, c11, c12, dist, [1])
    )

    input = SystemPerturbations(g1, 5, g2, 4, Dict)

    res = reweight(simulation,
                    input;
                    all_distances = true,
                    k = 1.0,
                    T = 1.0,
                    cutoff = 15.0,
            )

    @test res["a"].energies[1] ≈ [
        6.568418,
        6.00216,
        6.064767,
        10.650489,
        6.989886,
        7.265108,
        9.888358,
        9.993368,
        7.498538,
        6.640912     
    ] atol = 1.e-3


    @test res["b"].energies[1] ≈ [
        27.963684,
        24.046826,
        26.315003,
        41.866706,
        26.591741,
        27.698368,
        43.392474,
        40.395707,
        32.164286,
        26.216961      
    ] atol = 1.e-4
end

@testitem "Reweight with small trajectory using minimum distances and contributions from one group" begin #ADICIONAR TESTE COM CONTRIBUIÇÕES NÃO REPETIDAS (OXIGENIO E HIDROGENIOS, POR EX.)
    import PDBTools
    import OrderedCollections
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    g1 = PDBTools.selindex(simulation.atoms, at -> at.resname == "SOL" && at.residue in [100,120,140] && at.name != "MW")

    c1 = at -> at = true

    c2 = at -> at = true

    c11 = at -> at.name in ["HW1", "HW2"]

    c12 = at -> at.name in ["OW"]

    dist(r) = r

    Dict = OrderedCollections.OrderedDict(
        "a" => Perturbation(simulation.atoms, c1, c2, dist,[1]; one_gp=true), 
        "b" => Perturbation(simulation.atoms, c11, c12, dist, [-1,2]; one_gp=true)
    )

    input = SystemPerturbationsOneGroup(g1, 3, Dict)

    res = reweight(simulation,
                    input;
                    all_distances = false,
                    k = 1.0,
                    T = 1.0,
                    cutoff = 27.0,
            )

    @test res["a"].energies[1] ≈ [
        33.45344, 
        54.10159, 
        16.77879, 
        49.15888, 
        17.86500, 
        67.33274, 
        48.39516, 
        24.73005, 
        20.73088, 
        0.0
    ] atol = 1.e-3

    @test res["a"].probabilities[1] ≈ [
        2.96043532e-15
        3.19137596e-24
        5.16492547e-8
        4.47269872e-22
        1.7431271e-8
        5.7248292e-30
        9.5995092e-22
        1.8191801e-11
        9.9241468e-10
        0.9999999
    ] atol = 1.e-7


    @test res["a"].distances ≈ [
        2,
        3,
        1,
        2,
        1,
        3,
        2,
        1,
        1,
        0
    ] atol = 1.e-3

    @test res["b"].energies ≈ [
        -[12.58697, 54.10159, 16.778796108272978, 0, 0, 0, 0, 0, 0, 0],
        [25.17394, 108.20318, 33.557592216545956, 0, 0, 0, 0, 0, 0, 0]
    ] atol = 1.e-4

    @test res["b"].distances ≈ [
        1, 
        3, 
        1, 
        0, 
        0, 
        0, 
        0, 
        0, 
        0, 
        0    
    ] atol = 1.e-4
end

@testitem "Reweight with small trajectory using all distances and contributions from one group" begin
    import PDBTools
    import OrderedCollections
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    g1 = PDBTools.selindex(simulation.atoms, at -> at.resname == "SOL" && at.residue in [100,120,140] && at.name != "MW")

    c1 = at -> at = true

    c2 = at -> at = true

    c11 = at -> at.name in ["HW1", "HW2"]

    c12 = at -> at.name in ["OW"]

    dist(r) = r

    Dict = OrderedCollections.OrderedDict(
        "a" => Perturbation(simulation.atoms, c1, c2, dist, [1]; one_gp=true), 
        "b" => Perturbation(simulation.atoms, c11, c12, dist, [1]; one_gp=true)
    )

    input = SystemPerturbationsOneGroup(g1, 3, Dict)

    res = reweight(simulation,
                    input;
                    all_distances = true,
                    k = 1.0,
                    T = 1.0,
                    cutoff = 27.0,
            )

    @test res["a"].energies[1] ≈ [
        316.9570, 
        510.3088, 
        155.5428, 
        293.7871, 
        168.9915, 
        630.7433, 
        450.7814, 
        231.5664, 
        189.6638, 
        0.0
    ] atol = 1.e-3

    @test res["a"].distances == [
        18, 
        27, 
        9, 
        12, 
        9, 
        27, 
        18, 
        9, 
        9, 
        0
    ]

    @test res["b"].energies[1] ≈ [ 
        140.8590, 
        226.0973, 
        68.6984, 
        148.9806, 
        75.3751, 
        280.1446, 
        199.8053, 
        102.9400, 
        84.8140,
        0.0
    ] atol = 1.e-3

    @test res["b"].distances == [
        8, 
        12, 
        4, 
        6, 
        4, 
        12, 
        8, 
        4, 
        4, 
        0
    ]
end