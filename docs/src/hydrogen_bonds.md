```@meta
CollapsedDocStrings = true
```

# [Hydrogen bonds](@id hbonds)

!!! warning
    This is an experimental feature. Breaking changes may occur without 
    a breaking package release.

```@docs
hydrogen_bonds
```
```@example hbonds
using MolSimToolkit, PDBTools
using MolSimToolkit.Testing # to load test files
using Plots
# Build Simulation object
sim = Simulation(Testing.namd_pdb, Testing.namd_traj) 
# Compute h-bonds between protein and water
hbs = hydrogen_bonds(sim, "protein", "water")
plot(MolSimStyle, 1:length(sim), hbs;
    xlabel="frame",
    ylabel="number of H-bonds",
    label=""
)
```


