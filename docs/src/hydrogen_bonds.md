```@meta
CollapsedDocStrings = true
```

# [Hydrogen bonds](@id hbonds)

Computes the number of hydrogen bonds of a set of atoms, or between two sets 
of atoms, for each frame in a simualtion.

!!! warning
    This is an experimental feature. Breaking changes may occur without 
    a breaking package release.

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
hbs_prot_water = hydrogen_bonds(sim, "protein", "water")
# Plot 
plot(MolSimStyle, 
    [hbs_prot hbs_prot_water];
    xlabel="frame",
    ylabel="number of hydrogen bonds",
    label=["protein-protein" "protein-water"],
    legend=:outertopright,
)
```


