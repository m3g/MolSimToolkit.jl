module Resampling

using CellListMap.PeriodicSystems
using ..MolSimToolkit: Simulation, positions, unitcell

export new_weigths, intermol_perturb, L_J, dist, poly_decay, gaussian_decay

const testdir = "$(@__DIR__)/test"

include("./newweights.jl")
include("./ij_equations.jl")
include("./r_equations.jl")
include("sampling_graph.jl")

end

@testitem "resampling one frame" begin
    import PDBTools
    using MolSimToolkit.Resampling
    using MolSimToolkit.Resampling: testdir

    simulation = Simulation("$testdir/Testing_resampling.pdb", "/$testdir/Testing_resampling_one_frame.xtc")

    i1 = PDBTools.selindex(atoms(simulation), "resname TFE and name O")

    i2 = PDBTools.selindex(atoms(simulation), "residue 11")

    sum_of_dist = new_weigths(simulation, (i,j,d2) -> intermol_perturb(i, j, d2, dist), [i1[239]], i2, 25.0)
    @test sum_of_dist.energy ≈ [74.295543]
end

@testitem "resampling small trajectory" begin
    import PDBTools
    using MolSimToolkit.Resampling
    using MolSimToolkit.Resampling: testdir

    simulation = Simulation("$testdir/Testing_resampling.pdb", "/$testdir/Testing_resampling_small_trajectory.xtc")

    i1 = PDBTools.selindex(atoms(simulation), "index 97 or index 106")

    i2 = PDBTools.selindex(atoms(simulation), "residue 15 and name HB3")

    sum_of_dist = new_weigths(simulation, (i,j,d2) -> intermol_perturb(i, j, d2, dist), i1, i2, 25.0)
    @test sum_of_dist.energy ≈ [
    17.738965476707595, 15.923698293115915, 17.16614676290554, 
    19.33003841107648, 16.02329229247863, 19.639005665480983, 
    35.73986006775934, 21.88798265022823, 20.66180657974777, 
    16.845109623700647, 20.114166329136705, 24.68937611002383, 
    18.71136654132259, 20.41427025641757, 15.815250733848112,
    12.588332736178291, 36.50414116409441, 21.58454409077756, 
    25.40955804417851, 25.765519091038563, 18.200035069696835
    ]
end