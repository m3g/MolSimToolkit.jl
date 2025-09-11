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

Version 1.29.3
-------------
- ![INFO][badge-info] renamed `parse_server_output` to `parse_mvalue_server_sasa`, for clarity.
- ![INFO][badge-info] updated mvalue docs.

Version 1.29.2
-------------
- ![FEATURE][badge-feature] Experimental: mvalue calculator.
- ![INFO][badge-info] document the fields of the `UnitCell` type.
- ![INFO][badge-info] increase PDBTools.jl compat bound to 3.2.
- ![INFO][badge-info] install gromacs on CI runs on Linux, to test mvalue SASA computations.

Version 1.29.1
-------------
- ![BUGFIX][badge-bugfix] fix diagonal unitcell test for when cell is not orthorhombic and contains negative vector entries.

