module Reweighting

using ProgressMeter
using LinearAlgebra: diag
using CellListMap: ParticleSystem, map_pairwise, map_pairwise!
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

    Dict = OrderedCollections.OrderedDict(1 => Perturbation(simulation.atoms, c1, c2, dist))

    input = SystemPerturbations(g1, 1, g2, 3, Dict)

    res = reweight(simulation,
                    input;
                    all_distances = false,
                    k = 1.0,
                    T = 1.0,
                    cutoff = 12.0,
            )
    @test res[1].energy ≈ [
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

    Dict = OrderedCollections.OrderedDict(1 => Perturbation(simulation.atoms, c1, c2, dist))

    input = SystemPerturbations(g1, 1, g2, 3, Dict)

    res = reweight(simulation,
                    input;
                    all_distances = true,
                    k = 1.0,
                    T = 1.0,
                    cutoff = 12.0,
            )
    @test res[1].energy ≈ [
        5.450897,
        6.018059,
        6.05581,
        5.926243,
        5.067425,
        5.816867,
        5.229802,
        5.382748,
        6.127764,
        6.212621    
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
        "a" => Perturbation(simulation.atoms, c1, c2, dist), 
        "b" => Perturbation(simulation.atoms, c11, c12, dist)
    )

    input = SystemPerturbations(g1, 5, g2, 4, Dict)

    res = reweight(simulation,
                    input;
                    all_distances = false,
                    k = 1.0,
                    T = 1.0,
                    cutoff = 15.0,
            )

    @test res["a"].energy ≈ [
        6.568418,
        0,
        6.064767,
        0,
        0,
        0,
        9.888358,
        0,
        7.498538,
        0     
    ] atol = 1.e-3


    @test res["b"].energy ≈ [
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
        "a" => Perturbation(simulation.atoms, c1, c2, dist), 
        "b" => Perturbation(simulation.atoms, c11, c12, dist)
    )

    input = SystemPerturbations(g1, 5, g2, 4, Dict)

    res = reweight(simulation,
                    input;
                    all_distances = true,
                    k = 1.0,
                    T = 1.0,
                    cutoff = 15.0,
            )

    @test res["a"].energy ≈ [
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


    @test res["b"].energy ≈ [
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

    g1 = PDBTools.selindex(simulation.atoms, at -> at.resname == "SOL" && at.residue in [100,120,140,160,180] && at.name != "MW")

#    c1 = at -> at.name in ["HW1", "HW2"]

#    c2 = at -> at.name in ["O"]

    c1 = at -> at = true

    c2 = at -> at = true

    c11 = at -> at.name in ["HW1", "HW2"]

    c12 = at -> true

    dist(r) = r

    Dict = OrderedCollections.OrderedDict(
        "a" => Perturbation(simulation.atoms, c1, c2, dist), 
        "b" => Perturbation(simulation.atoms, c11, c12, dist)
    )

    input = SystemPerturbationsOneGroup(g1, 3, Dict)

    res = reweight(simulation,
                    input;
                    all_distances = false,
                    k = 1.0,
                    T = 1.0,
                    cutoff = 27.0,
            )

    @test res["a"].energy ≈ [
        136.05554345097482, 
        181.0981720027807, 
        152.37416783247951, 
        87.58817228480571, 
        76.6141583934713, 
        179.4770542337321, 
        48.39516109621037, 
        189.14883805230835, 
        110.43236754714565,
        63.52665469980029
    ] atol = 1.e-3

    @test res["a"].probability ≈ [
    8.503144223910598e-39
    2.332461665989855e-58
    6.958105738162862e-46
    9.521193580872743e-18
    5.554501296481423e-13
    1.1799321005648997e-57
    0.9999997317889862
    7.437974307845558e-62
    1.141782675765464e-27
    2.682104584745969e-7
    ] atol = 1.e-3


    @test res["a"].distances ≈ [
        7.0, 
        9.0, 
        7.0, 
        4.0, 
        4.0, 
        8.0, 
        2.0, 
        8.0, 
        5.0, 
        4.0
    ] atol = 1.e-3

    @test res["b"].energy ≈ [
        129.7620583116274, 
        114.46581597869775, 
        124.40288124305596, 
        87.58817228480571, 
        68.21332921937714, 
        121.37077327618356, 
        22.697185701332867, 
        176.103612512724, 
        101.97801390530358, 
        45.92474904398998    
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
        "a" => Perturbation(simulation.atoms, c1, c2, dist), 
        "b" => Perturbation(simulation.atoms, c11, c12, dist)
    )

    input = SystemPerturbationsOneGroup(g1, 3, Dict)

    res = reweight(simulation,
                    input;
                    all_distances = true,
                    k = 1.0,
                    T = 1.0,
                    cutoff = 27.0,
            )

    @test res["a"].energy ≈ [
        316.95709175011706, 
        510.30880475614055, 
        155.54282197934864, 
        293.78717961659873, 
        168.9915457063237, 
        630.7433027192693, 
        450.78141452462495, 
        231.56648096974462, 
        189.66385308799272, 
        0.0
    ] atol = 1.e-3

    @test res["a"].distances == [
        18.0, 
        27.0, 
        9.0, 
        12.0, 
        9.0, 
        27.0, 
        18.0, 
        9.0, 
        9.0, 
        0.0
    ]

    @test res["b"].energy ≈ [ #VIZINHOS MAIS PRÓXIMOS DO MolSimToolkit NÃO SERVIU PARA TESTAR OS RESULTADOS!!!!!!! FALAR COM LEANDRO
        140.8590610379914, 
        226.09738285953762, 
        68.69849119519456, 
        148.9806945785561, 
        75.37515540697854, 
        280.14468842660045, 
        199.80533848971726, 
        102.94000222371807, 
        84.81409883203321,
        0.0
    ] atol = 1.e-3

        @test res["b"].distances == [
        8.0, 
        12.0, 
        4.0, 
        6.0, 
        4.0, 
        12.0, 
        8.0, 
        4.0, 
        4.0, 
        0.0
    ]
end