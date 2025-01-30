module Testing

test_dir = @__DIR__ # current directory

namd_pdb = joinpath(test_dir, "data/namd/protein_in_popc_membrane", "structure.pdb")
namd_traj = joinpath(test_dir, "data/namd/protein_in_popc_membrane", "trajectory.dcd")

namd2_pdb = joinpath(test_dir, "data/namd/protein_in_water_tmao/structure.pdb")
namd2_traj = joinpath(test_dir, "data/namd/protein_in_water_tmao/trajectory.dcd")
namd2_traj_wrapped = joinpath(test_dir, "data/namd/protein_in_water_tmao/trajectory_wrapped.dcd")

end # module Testing
