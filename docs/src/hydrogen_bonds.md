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

## Hydrogen bond occupancy

The `hydrogen_bond_occupancy` function performs the same computation as
`hydrogen_bonds`, but instead of only counting the hydrogen bonds of each
frame, it keeps track of the identity (donnor, polar hydrogen, and acceptor
atoms) of every hydrogen bond found. This makes it possible to follow
individual hydrogen bonds along the trajectory, in the same way that
[`occupancy`](@ref) tracks individual solvent molecules at a binding site.

```@docs
hydrogen_bond_occupancy
HydrogenBondOccupancy
HBond
mean(::HydrogenBondOccupancy)
```

### Example

```@example hbond_occupancy
using MolSimToolkit, MolSimToolkit.Testing, Plots

sim = Simulation(Testing.namd_pdb, Testing.namd_traj)

hbo = hydrogen_bond_occupancy(sim, 
    "protein", 
    "protein" => "water"; 
    show_progress=false
)
```

The `list` field contains, for each frame, the [`HBond`](@ref)s found in that frame:

```@example hbond_occupancy
hbo["protein => water"].list
```

and the average number of hydrogen bonds per frame is obtained, as usual, with `mean`:

```@example hbond_occupancy
mean(hbo["protein => protein"])
```

## [Persistence of hydrogen bonds](@id hbond-persistence)

Because each `HBond` retains the identity of the donnor, hydrogen, and
acceptor atoms involved, [`intermittent_correlation`](@ref) can be applied
directly to a `HydrogenBondOccupancy`, to estimate the probability that a
given hydrogen bond found at some frame is still present `delta` frames
later (allowing the bond to break and reform in between):

```@example hbond_occupancy
ic_prot_prot = intermittent_correlation(
    hbo["protein => protein"]; maxdelta=4, show_progress=false
)
ic_prot_water = intermittent_correlation(
    hbo["protein => water"]; maxdelta=4, show_progress=false
)
```

```@example hbond_occupancy
plot(MolSimStyle,
    0:4,
    [ic_prot_prot ic_prot_water], 
    xlabel="Delta (frames)", ylabel="Probability",
    linewidth=2, marker=:circle,
    label=[ "protein-protein" "protein-water" ],
)
```

In this case, protein-water hydrogen bonds do not survive enough to be found on consecutive
frames saved from the trajectory, but protein-protein hydrogen bonds do. The characteristic
times of these decays can be obtained by fitting exponential decays, as shown in the
[Characteristic residence time](@ref occupancy-residence-time) section.

!!! note
    A hydrogen bond is considered to be the same along the trajectory only if
    it is formed by the same donnor, polar hydrogen, and acceptor atoms: if
    the bridging polar hydrogen changes while the donnor and acceptor remain
    the same, that counts as a different hydrogen bond.


