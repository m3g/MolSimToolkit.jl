module Testing

import LazyArtifacts # must be loaded for lazy artifacts
using LazyArtifacts: @artifact_str

export test_file

"""
    test_file(name::AbstractString)

Returns the path of the test file `name`. The test files are lazy artifacts (see
https://github.com/m3g/MolSimTestFiles), downloaded on the first call. The files are:

- `popc_membrane.pdb`, `popc_membrane.psf`, `popc_membrane.dcd`: NAMD simulation of a protein
  in a POPC membrane.
- `protein_tmao.pdb`, `protein_tmao.psf`, `protein_tmao.dcd`: NAMD simulation of a protein
  in water and TMAO.
- `short_nopbc.pdb`, `short_nopbc.dcd`: short NAMD trajectory without periodic boundary
  conditions.
- `mdlovofit.pdb`, `mdlovofit.dcd`: trajectory for structural alignment tests.
- `pressure_tensor.log` (with its input files `pressure_tensor.namd` and `tip3p.prm`): NAMD
  log of a water simulation, with the pressure tensor.
- `gromacs.pdb`, `gromacs.xtc`, `gromacs.tpr`: Gromacs simulation.
- `remd_gmx-2019.4.log`, `remd_gmx-5.0.4.log`: logs of Gromacs replica-exchange simulations.

"""
function test_file(name::AbstractString)
    # The files of the protein-TMAO system are shared with ComplexMixtures.jl
    dir = startswith(name, "protein_tmao.") ? artifact"protein_tmao_test_files" : artifact"molsimtoolkit_test_files"
    return joinpath(dir, name)
end

end # module Testing
