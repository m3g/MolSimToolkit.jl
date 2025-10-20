```@meta
CollapsedDocStrings = true
```

# [Hydrogen bonds](@id hbonds)

Computes the number of hydrogen bonds of a set of atoms, or between two sets 
of atoms, for each frame in a simulation.

```@docs
hydrogen_bonds
```
## Example

```@example hbonds
using MolSimToolkit, PDBTools
using MolSimToolkit.Testing # to load test files
using Plots
# Build Simulation object
sim = Simulation(Testing.namd_pdb, Testing.namd_traj) 
# Compute h-bonds of the protein with itself
hbs_prot = hydrogen_bonds(sim, "protein")
# Compute h-bonds between protein and water
hbs_prot_water = hydrogen_bonds(sim, "protein" => "water")
# Plot 
plot(MolSimStyle, 
    [hbs_prot["protein => protein"] hbs_prot_water["protein => water"]];
    xlabel="frame",
    ylabel="number of hydrogen bonds",
    label=["protein-protein" "protein-water"],
    legend=:outertopright,
)
```

Alternativelly, multiple selections, or pairs of selections can be provided, for faster computations,
```@example hbonds
hbs = hydrogen_bonds(sim, 
        "protein", 
        "protein" => "water",
        "protein" => "resname POPC",
)
``` 
The result can be converted directly to a `DataFrame`
```@example hbonds
using DataFrames, CSV
df = DataFrame(hbs)
```
and saved to CSV file with `CSV.write("hbonds.csv", df)`.

!!! note
    The order of the pairs, e. g. `"protein" => "water"` or `"water" => "protein"`, does not affect the
    result, as electronegative atoms of both groups will be considered as possible hydrogen bond donnors and/or acceptors.

    When a single selection is provided, e. g. `"protein"`, the hydrogen bonds within that selection are
    computed, with no repetitions.



