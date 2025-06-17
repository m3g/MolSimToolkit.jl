```@meta
CollapsedDocStrings = true
```

# System setup

## Packmol Input Creator

This module helps the setup of a `Packmol` input file, by computing the number of molecules
and sizes necessary to build different systems. Currently (as of version 1.27.0) the module
supports the construction of systems of a solute (macromolecule for example) solvated by a 
a single solvent or a mixture of two solvents. 

- [Solute-Solvent system](@ref)
- [Solute-Solvent-Cossolvent system](@ref)

### How to use it

```julia-repl
julia> using MolSimToolkit.PackmolInputCreator
```

## Running Packmol

`Packmol` can be run directly from within Julia using the `Packmol` Julia package:
```julia
using Packmol
run_packmol("./box.inp")
```
If everything runs correctly, the output PDB file of the system will be created.

## Solute-Solvent system

Here, `SolutionBoxUS` stands for `Solute (U)` and `Solvent (S)`.
The density of the solvent can be one of `"g/mL"` or `"mol/L"` (molarity).

```@docs
PackmolInputCreator.SolutionBoxUS
write_packmol_input(::SolutionBoxUS)
```

### Setting up the system properties

We initialize the system data structure, given the PDB files of *one molecule* of the
polymer (`poly_h.pdb`), and *one molecule* of water:

```julia
# Directory of test files
test_dir = MolSimToolkit.PackmolInputCreator.PackmolInputCreatorDirectory*"/test"
# Construction of system data structure
system = SolutionBoxUSC(
    solute_pdbfile = "$test_dir/data/poly_h.pdb",
    solvent_pdbfile = "$test_dir/data/water.pdb",
    density = 1.0,
    density_units = "g/mL",
    solute_molar_mass = nothing, # optional
    solvent_molar_mass = nothing, # optional
)
```

The molar masses of the components can be provided explicitly by the user. If not, they
will be computed from the atom types in the PDB files, but this may fail if the mass
of some atom type is unknown.

Finally, we can generate an input file for `Packmol` with:

```julia
write_packmol_input(
    system; 
    margin = 20.0, 
    cubic = true,
    input = "box.inp",
    output = "system.pdb"
)
```

The `input` parameter provides the name of the input file for `Packmol` that will be generated. 

The `margin` parameter sets the size of the box, which will take into consideration the maximum and
minimum dimensions of the solute.  If `cubic` is set to `true` the box will be cubic, otherwise
it will be orthorhombic but with different length in each direction, respecting the margin
provided.

Alternatively, the size of the box can be provided explicitly
with the `box_sides = [ a, b, c ]` parameters, where `a`, `b`, and `c` are the lengths of the box
in each dimension. 

## Solute-Solvent-Cossolvent system

Here, `SolutionBoxUSC` stands for `Solute (U)`, `Solvent (S)`, and `Cossolvent (C)`.  
The concentration units can be one of `"mol/L"` (molarity), `"x"` (molar fraction),
`"vv"` (volume fraction), and `"mm"` (mass fraction). The density is assumed
to be in `g/mL`. 

```@docs
SolutionBoxUSC
write_packmol_input(::SolutionBoxUSC)
convert_concentration
convert_density_table!
```

### Setting up the system properties

Here, we setup a system of a polymer solvated by water and ethanol. The densities as 
a function of the molar fraction of ethanol are available in a data table:

```julia
density_table = [
# x cossolvent (ethanol)     density (g/mL)
         0.0000                 0.9981
         0.0416                 0.9820
         0.0890                 0.9685
         0.1434                 0.9537
         0.2066                 0.9369
         0.2809                 0.9151
         0.3695                 0.8923
         0.4769                 0.8685
         0.6098                 0.8450
         0.7786                 0.8195
         1.0000                 0.7906
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
    concentration_units = "x", # molar fraction
    solute_molar_mass = nothing, # optional
    solvent_molar_mass = nothing, # optional
    cossolvent_molar_mass = nothing, # optional
)
```

The molar masses of the components can be provided explicitly by the user. If not, they
will be computed from the atom types in the PDB files, but this may fail if the mass
of some atom type is unknown.

!!! tip
    The density table can be converted among different units with the function `convert_density_table!`,
    which acts on the `SystemBox` object. For example:
    ```julia-repl
    julia> convert_density_table!(system, "mol/L")
    ```

Finally, we can generate an input file for `Packmol` with:

```julia
write_packmol_input(
    system; 
    concentration = 0.5,
    margin = 20.0, 
    cubic = true,
    input = "box.inp",
    output = "system.pdb"
)
```

The concentration can be given in molar fraction (`x`), molarity (`mol/L`), or volume fraction (`vv`). 

The `input` parameter provides the name of the input file for `Packmol` that will be generated. 

The `margin` parameter sets the size of the box, which will take into consideration the maximum and
minimum dimensions of the solute. If `cubic` is set to `true` the box will be cubic, otherwise
it will be orthorhombic but with different length in each direction, respecting the margin
provided.

Alternatively, the size of the box can be provided explicitly
with the `box_sides = [ a, b, c ]` parameters, where `a`, `b`, and `c` are the lengths of the box
in each dimension. 
