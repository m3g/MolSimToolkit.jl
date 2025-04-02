module Testing

test_dir = @__DIR__ # current directory

namd_pdb = joinpath(test_dir, "data/namd/protein_in_popc_membrane", "structure.pdb")
namd_traj = joinpath(test_dir, "data/namd/protein_in_popc_membrane", "trajectory.dcd")

namd2_pdb = joinpath(test_dir, "data/namd/protein_in_water_tmao/structure.pdb")
namd2_traj = joinpath(test_dir, "data/namd/protein_in_water_tmao/trajectory.dcd")

mdlovofit_pdb = joinpath(test_dir, "data/mdlovofit/structure.pdb")
mdlovofit_traj = joinpath(test_dir, "data/mdlovofit/trajectory.dcd")

short_nopbc_pdb = joinpath(test_dir, "data/namd/short_nopbc.pdb")
short_nopbc_traj = joinpath(test_dir, "data/namd/short_nopbc.dcd")

end # module Testing
