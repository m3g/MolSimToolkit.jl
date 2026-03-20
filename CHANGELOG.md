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

Version 2.0.0-DEV
--------------
- ![BREAKING][badge-breaking] `atoms` function was removed in favor of `get_atoms`.
- ![BREAKING][badge-breaking] `get_frame` returns a `Frame` object, not anymore a vector of atoms. Use `get_atoms` instead.
- ![INFO][badge-info] The internal representation of frames and trajectories is done by new `MolSimToolkit` structs `Frame` and `Trajectory`. `Chemfiles` objects are not exposed anymore to the user. 

Version 1.33.0
--------------
- ![FEATURE][badge-feature] Introduce the option to compute the RMSD of a structure, given the alignment of other structure (or parts of); i. e. the `rmsd_of` option of `rmsd`. 
- ![FEATURE][badge-feature] Support string selections in `rmsd` and `rmsd_matrix` functions.
- ![ENHANCEMENT][badge-enhancement] Reconstruct polymer and complex structure before computing rmsds, to avoid issues with broken molecule coordinates through PBCs.

Version 1.32.3
--------------
- ![FIX][badge-fix] fix units of autocorrelation plot in block average plotting.
- ![INFO][badge-info]: add reference to the section about computing the number of effective samples.
- ![INFO][badge-info]: add space on top of title of block average plot to improve figure, adjust the number of digits of plotted mean according to effective standard error.

Version 1.32.2
--------------
- ![FIX][badge-fix] fix plot of t95 data in block-average plotting extension.

Version 1.32.1
--------------
- ![FIX][badge-fix] all BlockAverage plots are now printed with correct time units, and the use of fractional time units (fractional delays between samples) is properly handled.

Version 1.32.0
--------------
- ![FEATURE][badge-feature] Support units in data values of block-average functions, and add the `dt` keyword parameter to define the time-step, also supporting units. 

Version 1.31.2
--------------
- ![ENHANCEMENT][badge-enhancement] Add safeguard for rounding errors in angle computation in hydrogen-bonds.

Version 1.31.1
--------------
- ![INFO][badge-info] Fix tests for compatibility with PDBTools. 3.11.0 (removed mvalue-related functions).
- ![INFO][badge-info] Explicit imports for all PDBTools functions. 

Version 1.31.0
--------------
- ![FEATURE][badge-feature] Compute integrated correlation time and effective number of samples in BlockAverage, to plot. The exponential fit is new performed for the set of data in the 95% confidence interval. 
- ![ENHANCEMENT][badge-enhancement] Reconstruct structure before computing secondary structure.

Version 1.30.1
--------------
- ![ENHANCEMENT][badge-enhancement] Computation of hydrogen bonds faster after adding function barriers.
- ![INFO][badge-info] Organize `hydrogen_bonds` code.

Version 1.30.0
--------------
- ![FEATURE][badge-feature] Make `hydrogen_bonds` a stable interface, and it accepts a multiple selections and pairs of selections, to compute many h-bond data at once.

Version 1.29.15
---------------
- ![ENHANCEMENT][badge-enhancement] Computation of hydrogen bonds is much faster when running in parallel.

Version 1.29.14
---------------
- ![FEATURE][badge-experimental] Add `hydrogen_bonds` function.

Version 1.29.13
---------------
- ![FEATURE][badge-experimental] Support for PBCs in `delta_sasa_per_residue` with the `unitcell` keyword.
- ![FEATURE][badge-experimental] Changed interface of `mvalue` functions: the funcions now receive vectors of atoms, and thus the keyword `pdbname`  was renamed to `atoms`. 

Version 1.29.12
---------------
- ![FEATURE][badge-experimental] Add backbone and sidechain selection functions, `n_dots` and `ignore_hydrogen` as optional parameters to `[gmx_]delta_sasa_per_residue` functions.

Version 1.29.11
---------------
- ![FEATURE][badge-experimental] Implement `delta_sasa_per_restype` to compute SASAs using PDBTools.
- ![INFO][badge-info] Requires PDBTools 3.5.5

Version 1.29.10
---------------
- ![INFO][badge-info] Support EasyFit 0.7.1+

Version 1.29.9
--------------
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

