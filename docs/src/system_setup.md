
# System setup

## Packmol Input Creator

This module helps the setup of a `Packmol` input file, by computing the number of molecules
and sizes necessary to build different systems. Currently (as of version 1.2.0) the module
supports the construction of systems of a solute solvated by a mixture of two solvents. 

### How to use it

```julia-repl
julia> using MolSimToolkit.PackmolInputCreator
```

### Setting up the system properties

Here, we setup a system of a polymer solvated by water and ethanol. The densities as 
a function of the molar fraction of ethanol are available in a data table:

```julia
density_table = [
# x ethanol     density
   0.0000       0.9981
   0.0416       0.9820
   0.0890       0.9685
   0.1434       0.9537
   0.2066       0.9369
   0.2809       0.9151
   0.3695       0.8923
   0.4769       0.8685
   0.6098       0.8450
   0.7786       0.8195
   1.0000       0.7906
]
```

Next, we initialize the system data structure, given the PDB files of *one molecule* of the
polymer (`poly_h.pdb`), and *one molecule* of water and ethanol:

```julia
# Directory of test files
test_dir = MolSimToolkit.PackmolInputCreator.PackmolInputCreatorDirectory*"/test"
# Construction of system data structure
system = SolutionBoxUSC(
    solute_pdbfile = "$test_dir/data/poly_h.pdb",
    solvent_pdbfile = "$test_dir/data/water.pdb",
    cossolvent_pdbfile = "$test_dir/data/ethanol.pdb",
    density_table = density_table,
)
```

Here, `SolutionBoxUSC` stands for `Solute (U)`, `Solvent (S)`, and `Cossolvent (C)`. 

Finally, we can generate an input file for `Packmol` with:

```julia
write_packmol_input(
    system; 
    concentration = 0.5, cunit = "x", 
    margin = 20.0, 
    input = "box.inp",
    output = "system.pdb"
)
```

The concentration can be given in molar fraction (`x`), molarity (`mol/L`), or volume fraction (`vv`). 

The `input` parameter provides the name of the input file for `Packmol` that will be generated. 

The `margin` parameter sets the size of the box, which will take into consideration the maximum and
minimum dimensions of the solute. Alternatively, the size of the box can be provided explicitly
with the `box_sides = [ a, b, c ]` parameters, where `a`, `b`, and `c` are the lengths of the box
in each dimension. 

## Running Packmol

`Packmol` can be run directly from within Julia using the `Packmol` Julia package:
```julia
using Packmol
run_packmol("./box.inp")
```
If everything runs correctly, the output file `system.pdb` will be generated.

## Help entries

```@autodocs
Modules = [MolSimToolkit.PackmolInputCreator]
Order = [:function, :type]
```



