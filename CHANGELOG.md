MolSimToolkit.jl Changelog
===========================
  
[badge-breaking]: https://img.shields.io/badge/BREAKING-red.svg
[badge-deprecation]: https://img.shields.io/badge/Deprecation-orange.svg
[badge-feature]: https://img.shields.io/badge/Feature-green.svg
[badge-experimental]: https://img.shields.io/badge/Experimental-yellow.svg
[badge-enhancement]: https://img.shields.io/badge/Enhancement-blue.svg
[badge-bugfix]: https://img.shields.io/badge/Bugfix-purple.svg
[badge-fix]: https://img.shields.io/badge/Fix-purple.svg
[badge-info]: https://img.shields.io/badge/Info-gray.svg

Version 1.29.9
-------------
- ![BUGFIX][badge-bugfix] Fix call to custom `gmx` executable in `gmx_sasa_per_restype`. 

Version 1.29.8
-------------
- ![FEATURE][badge-experimental] Rename `gmx_sasa` to `gmx_delta_sasa_per_restype` for clarity (and also the currently internal `gmx_sasa_single` to `gmx_sasa_per_restype`).
- ![FEATURE][badge-experimental] Support providing the `gmx` executable to `gmx_delta_sasa_per_restype` function as a keyword argument.
- ![BUGFIX][badge-bugfix] On `mvalue`, run over residue names and convert to types internally, because the `sasas` data will be based on residue names.

Version 1.29.7
-------------
- ![INFO][badge-experimental] use `threeletter` to convert from residue names to residue types when computing mvalues, thus supporting all alternate names `PDBTools.protein_residues` define.

Version 1.29.6
-------------
- ![INFO][badge-info] use `dihedral` and `dihedrals` from MolSimToolkitShared.

Version 1.29.5
-------------
- ![INFO][badge-experimental] renamed `run_gmx_delta_sasa_per_restype` to `gmx_delta_sasa_per_restype`.
- ![INFO][badge-info] documentation updated for mvalue.

Version 1.29.4
-------------
- ![FEATURE][badge-experimental] Support for the AutoBolen model of mvalue calcuation.
- ![FEATURE][badge-experimental] Code support for other cosolvents than urea in the calculation of m-values.
- ![INFO][badge-info] removed `atomic_mass` overload (with piracy) of PDBTools function.

Version 1.29.3
-------------
- ![INFO][badge-experimental] renamed `parse_server_output` to `parse_mvalue_server_sasa`, for clarity.
- ![INFO][badge-info] updated mvalue docs.

Version 1.29.2
-------------
- ![FEATURE][badge-experimental] Experimental: mvalue calculator.
- ![INFO][badge-info] document the fields of the `UnitCell` type.
- ![INFO][badge-info] increase PDBTools.jl compat bound to 3.2.
- ![INFO][badge-info] install gromacs on CI runs on Linux, to test mvalue SASA computations.

Version 1.29.1
-------------
- ![BUGFIX][badge-bugfix] fix diagonal unitcell test for when cell is not orthorhombic and contains negative vector entries.

