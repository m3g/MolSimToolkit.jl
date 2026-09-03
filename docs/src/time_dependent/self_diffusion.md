```@meta
CollapsedDocStrings = true
```
# Self-diffusion coefficient

The `mean_square_displacement` function computes, for a set of molecules, the
mean square displacement (MSD) of their centers of mass as a function of the
time lag between frames. The `self_diffusion_coefficient` function then
estimates the self-diffusion coefficient from the linear (diffusive) regime
of that curve, using the Einstein relation.

!!! note
    Trajectory files typically store coordinates wrapped back into the
    simulation box by periodic boundary conditions. Computing a mean square
    displacement directly from such coordinates would introduce spurious
    jumps whenever a molecule crosses a periodic boundary between two
    frames. To avoid this, `mean_square_displacement` reconstructs, for each
    molecule independently, a continuous ("unwrapped") trajectory: at each
    frame, the periodic image closest to that molecule's position at the
    *previous* frame is chosen. This assumes that consecutive frames are
    much closer in time than the time it takes a molecule to diffuse across
    half a box length, which is the usual case.

```@docs
mean_square_displacement
self_diffusion_coefficient
```

## Example: self-diffusion of TMAO

Here we use the `MolSimToolkit.Testing` test data, a short NAMD trajectory of
a protein solvated by water and TMAO, and compute the mean square
displacement of the TMAO molecules:

```@example self_diffusion
using MolSimToolkit, PDBTools, MolSimToolkit.Testing, Plots

sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj)
tmao = select(get_atoms(sim), "resname TMAO")

msd = mean_square_displacement(sim, tmao; natomspermol=14, maxdelta=6, show_progress=false)
```

```@example self_diffusion
plot(MolSimStyle,
    0:6, parent(msd), # parent(msd) is required for msd is an OffsetArray.
    xlabel="Delta (frames)", ylabel="MSD / Å²",
    linewidth=2, marker=:circle,
    label="TMAO",
)
```

## Estimating the self-diffusion coefficient

The `self_diffusion_coefficient` function fits a straight line to (a chosen
range of) the MSD curve, and returns the self-diffusion coefficient obtained
from its slope via the Einstein relation, `MSD(t) = 2 * dim * D * t`:

```@example self_diffusion
D = self_diffusion_coefficient(msd)
```

By default the whole curve (except the trivial `delta = 0` point) is used in
the fit. The `dt` keyword converts `delta` from frames into physical time
units, for instance `self_diffusion_coefficient(msd; dt=0.1)` if consecutive
frames are separated by `0.1` ns, which would return `D` in units of Å²/ns.

!!! note
    Only a limited range of `delta` typically falls in the diffusive
    (linear) regime: very short times are dominated by ballistic motion, and
    very long times become increasingly noisy, since fewer pairs of frames
    remain available to average over as `delta` approaches the length of the
    trajectory. Always inspect the `msd` curve (as plotted above) before
    trusting a fit, and use the `mindelta` and `maxdelta` keywords of
    `self_diffusion_coefficient` to restrict the fit to the region that
    actually looks linear.

As above, the trajectory used here is very short (20 frames), so this
estimate is shown only for illustration: with longer, production-quality
trajectories, a clear diffusive regime should be identifiable in the MSD
curve before fitting.
