[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://m3g.github.io/MolSimToolkit.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://m3g.github.io/MolSimToolkit.jl/dev)
[![Tests](https://img.shields.io/badge/build-passing-green)](https://github.com/m3g/MolSimToolkit.jl/actions)
[![codecov](https://codecov.io/gh/m3g/MolSimToolkit.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/m3g/MolSimToolkit.jl)
[![Aqua QA](https://JuliaTesting.github.io/Aqua.jl/dev/assets/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

# MyDevMolSimToolkit

[MolSimToolkit.jl](https://github.com/m3g/MolSimToolkit.jl) provides a set of tools to 
analyse molecular dynamics simulations, and a framework for the development of custom
analysis tools. This repository is a experimental version of the package for Work in Progress functionalities
# MolSimToolkit.jl

**MolSimToolkit.jl** is a Julia package for analyzing molecular dynamics (MD) simulations. It provides fast, flexible tools for trajectory analysis, structural alignment, secondary structure mapping, and more. The toolkit is designed to be lightweight and extensible, making it easy to build custom analysis workflows.

## Features

- **Robust structural alignment (MDLovoFit)**
- **Secondary structure mapping and plotting helpers**
- **Hydrogen bond analysis**
- **Pairwise and molecular distance analysis / coordination numbers**
- **Utilities for system setup (integrated with Packmol)** 
- **Block averaging and replica exchange analysis**
- **Integration with [PDBTools.jl](https://github.com/m3g/PDBTools.jl)**

## Installation

Install via Julia's package manager:

```julia
import Pkg
Pkg.add("MolSimToolkit")
```

It is recommended to also install `PDBTools` for enhanced functionality:

```julia
Pkg.add("PDBTools")
```

## Quick Start

```julia
using MolSimToolkit
```

### Example: Robust Alignment (MDLovoFit)

```julia
using MolSimToolkit, MolSimToolkit.Testing
sim = Simulation(Testing.mdlovofit_pdb, Testing.mdlovofit_traj)
fractions = map_fractions(sim)                # Scan fractions and RMSDs
aligned = mdlovofit(sim; fraction=0.8)        # Align trajectory using 80% least-mobile atoms
```

### Example: Secondary Structure Map

```julia
using MolSimToolkit, MolSimToolkit.Testing
sim = Simulation(Testing.namd_pdb, Testing.namd_traj)
ssmap = ss_map(sim; selection="residue >= 30 and residue <= 35")
hcontent = ss_mean(ssmap; class="H")
```

### Example: Hydrogen Bonds

```julia
using MolSimToolkit, PDBTools
sim = Simulation(Testing.namd_pdb, Testing.namd_traj)
hbs = hydrogen_bonds(sim, "protein" => "water")
```

### Example: Molecular Minimum Distances

```julia
using MolSimToolkit
list = MolSimToolkit.MolecularMinimumDistances.minimum_distances(
        xpositions=water, 
        ypositions=protein, 
        xn_atoms_per_molecule=3, 
        cutoff=12.0, 
        unitcell=[83.115,83.044,83.063]
)
```

## Documentation

- [Stable Docs](https://m3g.github.io/MolSimToolkit.jl/stable)
- [Dev Docs](https://m3g.github.io/MolSimToolkit.jl/dev)

## Contributing

Contributions are welcome! Please open issues or pull requests. 

## License

MIT License. See [LICENSE](LICENSE) for details.
