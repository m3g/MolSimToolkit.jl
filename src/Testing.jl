module Testing

src_dir = @__DIR__ # current directory
test_dir = joinpath(src_dir, "..", "test")

namd_pdb = joinpath(test_dir, "data/namd", "structure.pdb")
namd_traj = joinpath(test_dir, "data/namd", "structure.dcd")

end # module Testing
