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
list of solvent molecules found at the site, and the corresponding
site-solvent distances. This object can then be used to compute simple
summary statistics (with `mean`), to characterize how long, on average,
individual solvent molecules remain at the site (with
[`intermittent_correlation`](@ref)), or to resolve that residence time as a
function of the distance to the site (with
[`intermittent_correlation_profile`](@ref) and [`residence_time`](@ref)).

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

The `distances` field carries the corresponding minimum distances between
the protein and each of those molecules, such that `occ.distances[i][j]` is
the distance of the molecule `occ.list[i][j]`:

```@example occupancy
occ.distances
```

`n_solvent_molecules` is the total number of TMAO molecules considered, and
`dmax` is the cutoff used, that is, the maximum distance at which a molecule
is considered to be at the site:

```@example occupancy
occ.n_solvent_molecules, occ.dmax
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

## [Residence time as a function of the distance](@id occupancy-residence-profile)

The correlation function above mixes molecules that are tightly bound to the
protein with molecules that are only loosely associated with it, at the edge
of the cutoff. Since the `Occupancy` object also stores the site-solvent
distances, the correlation function can be resolved by distance, with
[`intermittent_correlation_profile`](@ref).

The distances from `0` to `dmax` are split into bins of constant width
`delta_r`, and the lower edges of consecutive bins are displaced by `step_r`.
Using `step_r` smaller than `delta_r` produces overlapping bins, and thus a
smooth (quasi-continuous) profile. For each bin, the correlation at `delta` is
the probability that a molecule found in that bin at some frame is found in
the *same* bin `delta` frames later.

Here we use a larger cutoff (6 Å), such that the profile covers the first
solvation shells of the protein:

```@example occupancy
occ6 = occupancy(
    sim, protein, tmao;
    solvent_natomspermol=14, cutoff=6.0,
    show_progress=false,
)

profile = intermittent_correlation_profile(occ6;
    delta_r=1.0, step_r=0.25, maxdelta=4,
    show_progress=false,
)
```

The `r` field contains the center of each bin, `correlations` the
corresponding correlation functions, and `counts` the number of observations
found in each bin (bins with no observations produce `NaN` correlations):

```@example occupancy
plot(MolSimStyle,
    0:4, [parent(c) for c in profile.correlations[5:5:end]],
    xlabel="Delta (frames)", ylabel="Probability",
    linewidth=2,
    labels=hcat(["r = $(r) Å" for r in profile.r[5:5:end]]...),
)
```

### The 50% residence time

The [`residence_time`](@ref) function returns, for each bin, the time at
which the correlation function falls below a given `threshold`. With the
default `threshold=0.5`, this is the time after which half of the molecules
initially in the bin have left it: the 50% residence time. The value is
obtained by linear interpolation between the two delta-steps that bracket the
threshold, and is `NaN` for bins in which the correlation never falls below
it (or for empty bins).

```@example occupancy
t50 = residence_time(profile; threshold=0.5)

plot(MolSimStyle,
    profile.r, t50,
    xlabel="Distance to the protein / Å",
    ylabel="50% residence time / frames",
    linewidth=2, marker=:circle,
    label=nothing,
)
```

Use the `dt` keyword to convert the result from frames to time units, for
instance `residence_time(profile; dt=0.1)` if consecutive frames are
separated by 0.1 ns.

As above, the trajectory used here is very short (20 frames), so the profile
is shown only for illustration.

## Reference functions

```@autodocs
Modules = [ MolSimToolkit ]
Pages = [
    "occupancy.jl",
]
```
