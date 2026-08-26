```@meta
CollapsedDocStrings = true
```
# Site occupancy

The `occupancy` function computes, for each frame of a simulation, which
solvent molecules are found within a given cutoff distance of a binding
site, defined by a set of atoms (for example, the atoms lining a cavity, or
the surface of a protein). 

The distance considered is the **minimum distance** between the site atoms
and the atoms of each solvent molecule, computed with the same
`minimum_distances!` machinery used by [`coordination_number`](@ref).

The function returns an `Occupancy` object, wrapping, for each frame, the
list of solvent molecules found at the site. This object can then be used to
compute simple summary statistics (with `mean`) or to characterize how long,
on average, individual solvent molecules remain at the site (with
[`intermittent_correlation`](@ref)).

## Example: TMAO occupancy at a protein binding site

Here we use the `MolSimToolkit.Testing` test data, a short NAMD trajectory of
a protein solvated by water and TMAO, and compute the occupancy of the
protein surface by TMAO molecules, using a 3 Å cutoff:

```@example occupancy
using MolSimToolkit, PDBTools, MolSimToolkit.Testing, Plots

sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj)
protein = select(get_atoms(sim), "protein")
tmao = select(get_atoms(sim), "resname TMAO")

occ = occupancy(
    sim, protein, tmao;
    solvent_natomspermol=14, cutoff=3.0,
    show_progress=false,
)
```

The `list` field of the `Occupancy` object contains, for each frame, the
indices of the TMAO molecules (1-based, relative to the `tmao` selection)
found within the cutoff distance of the protein:

```@example occupancy
occ.list
```

and `n_solvent_molecules` is the total number of TMAO molecules considered:

```@example occupancy
occ.n_solvent_molecules
```

## Average occupancy

The average number of TMAO molecules found at the protein surface, per
frame, is obtained with `mean`:

```@example occupancy
mean(occ)
```

## Intermittent correlation of the occupancy

The [`intermittent_correlation`](@ref) function, applied to an `Occupancy`
object, estimates the probability that a TMAO molecule found at the site at
some frame is still found there `delta` frames later. This probability
decays as molecules exchange between the site and the bulk solvent, and can
be used to characterize the residence time of the solvent at the site:

```@example occupancy
ic = intermittent_correlation(occ; maxdelta=4, show_progress=false)
```

```@example occupancy
plot(MolSimStyle,
    0:4, parent(ic), # parent(c) is required for c is an OffsetArray.
    xlabel="Delta (frames)", ylabel="Probability",
    linewidth=2, marker=:circle,
    label="TMAO at protein surface",
)
```

### [Characteristic residence time](@id occupancy-residence-time)

In this case, the decay of the intermittent correlation function can be fit to a double 
exponential, $$c(\delta) = a\exp(-\delta/\tau)$$, using `EasyFit.fitexpdecay`
(already a dependency of MolSimToolkit.jl), to extract a characteristic
residence time $$\tau$$, in units of frames:

```@example occupancy
using EasyFit 
using JuMP, Ipopt # required for constrained non-linear fitting

fit = fitexpdecay(ic; n=2, c=mean(occ)/occ.n_solvent_molecules)
tau = fit.b
```

The constant term is set to the ratio of the mean coordination number and
the total number of solvent molecules (correlation at long times).

```@example occupancy
plot(MolSimStyle,
    0:4, parent(ic),
    seriestype=:scatter,
    xlabel="Delta (frames)", ylabel="Probability",
    marker=:circle,
    label="TMAO at protein surface",
)
plot!(fit.x, fit.y, linewidth=2, label="Bi-Exponential fit")
```

Here the trajectory is very short (20 frames), so this characteristic time is
shown only for illustration: with longer, production-quality trajectories,
$$\tau$$ provides a quantitative estimate of how long a solvent molecule
typically remains bound to the site.

## Reference functions

```@autodocs
Modules = [ MolSimToolkit ]
Pages = [
    "occupancy.jl",
]
```
