module Testing

export test_file

#
# test_file(name): path of the test file `name`. The test files are lazy artifacts (see
# Artifacts.toml), downloaded on the first call. The method is defined in the TestFilesExt
# extension, loaded with `using LazyArtifacts`. The files are:
#
# - popc_membrane.pdb, popc_membrane.psf, popc_membrane.dcd: NAMD simulation of a protein in
#   a POPC membrane.
# - protein_tmao.pdb, protein_tmao.psf, protein_tmao.dcd: NAMD simulation of a protein in
#   water and TMAO (shared with ComplexMixtures.jl).
# - short_nopbc.pdb, short_nopbc.dcd: short NAMD trajectory without periodic boundary
#   conditions.
# - mdlovofit.pdb, mdlovofit.dcd: trajectory for structural alignment tests.
# - pressure_tensor.log (with its input files pressure_tensor.namd and tip3p.prm): NAMD log
#   of a water simulation, with the pressure tensor.
# - gromacs.pdb, gromacs.xtc, gromacs.tpr: Gromacs simulation.
# - remd_gmx-2019.4.log, remd_gmx-5.0.4.log: logs of Gromacs replica-exchange simulations.
#
function test_file end

end # module Testing
