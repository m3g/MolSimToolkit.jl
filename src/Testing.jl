module Testing

src_dir = @__DIR__ # current directory
test_dir = joinpath(src_dir, "..", "test")

namd_pdb = joinpath(test_dir, "data/namd/protein_in_popc_membrane", "structure.pdb")
namd_traj = joinpath(test_dir, "data/namd/protein_in_popc_membrane", "trajectory.dcd")

namd2_pdb = joinpath(test_dir, "data/namd/protein_in_water_tmao/structure.pdb")
namd2_traj = joinpath(test_dir, "data/namd/protein_in_water_tmao/trajectory.dcd")

end # module Testing
