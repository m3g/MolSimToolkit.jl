```@meta
CollapsedDocStrings = true
```

# [Protein transfer free energy calculator](@id mvalues)

!!! warning
    This is an experimental feature. Breaking changes may occur without 
    a breaking package release.

## Transfer free energies throughout a simulation

```@docs
transfer_free_energy(::Simulation, ::String)
transfer_free_energy_frames
```

The `transfer_free_energy` function computes the transfer free energy of a protein (in kcal/mol at 1M cosolvent),
averaged over all frames of a simulation, using the Tanford transfer model:

```@example tfe
using MolSimToolkit, MolSimToolkit.Testing
import PDBTools # to choose transfer model
sim = Simulation(Testing.namd_pdb, Testing.namd_traj)
tfe = transfer_free_energy(sim, "urea"; model=PDBTools.MoeserHorinek)
```

Per-residue contributions are available in `tfe.residue_contributions_bb` and `tfe.residue_contributions_sc`.

To obtain the transfer free energy for each frame individually (rather than an average),
use `transfer_free_energy_frames`. It returns a `TransferFreeEnergyFrames` object that
supports indexing and iteration:

```@example tfe
tfe_frames = transfer_free_energy_frames(sim, "urea"; model=PDBTools.MoeserHorinek)
```

```@example tfe
tfe_frames[1].tot  # total TFE of the first frame
```

```@example tfe
tot_per_frame = [ t.tot for t in tfe_frames ]
```

Note that the contributions of each residue, at each frame, are stored in `tfe_frames[i].residue_contributions_bb` (or `_sc`), for
backbone and side-chains, and thus this array might require a significant amount of memory. But, for example, to 
obtain the contributions to the transfer free energies of a subset of the protein, at each frame, do:

```@example tfe
tfe_10_15 = [ 
    sum(t.residue_contributions_bb[10:15]) + sum(t.residue_contributions_sc[10:15]) 
    for t in tfe_frames
] 
```

## Computing m-values 

The *m*-value of a frame relative to a reference structure is the difference in transfer free energy
between that frame and the reference. Using `transfer_free_energy_frames` together with a reference
TFE computed for the first frame via `PDBTools.transfer_free_energy`, this is straightforward:

```@example mvalues
using MolSimToolkit, MolSimToolkit.Testing
import PDBTools
sim = Simulation(Testing.namd_pdb, Testing.namd_traj)
```

First, we compute TFE of the reference structure (first frame):
```@example mvalues
atoms = get_atoms(sim)
protein_ref = PDBTools.select(atoms, "protein and not element H")
tfe_ref = PDBTools.transfer_free_energy(protein_ref, "urea")
```

Now, compute TFE for each frame of the trajectory:
```@example mvalues
tfe_frames = transfer_free_energy_frames(sim, "urea")
```
The m-values are difference in TFE between each frame and the reference
```@example mvalues
m = (
    tot = [ t.tot - tfe_ref.tot for t in tfe_frames ],
    bb  = [ t.bb  - tfe_ref.bb  for t in tfe_frames ],
    sc  = [ t.sc  - tfe_ref.sc  for t in tfe_frames ],
)
```

Plotting the results, we obtain the $\Delta_{\textrm{ref}\rightarrow\textrm{target}}\Delta G^{T}$, the difference in transfer free energy from water to urea at 1M, of the target structure at each frame relative to the reference structure:

```@example mvalues
using Plots, Statistics
plt = plot(MolSimStyle, layout=(1,2))
plot!(plt, m.tot; label="Total", lw=2, sp=1)
plot!(plt, m.bb; label="Backbone", lw=2, sp=1)
plot!(plt, m.sc; label="Sidechains", lw=2, sp=1)
plot!(xlabel="frame", ylabel="m-value / (kJ/mol)", sp=1)
bar!(plt, [1  2  3],
    [mean(m.tot)  mean(m.bb)  mean(m.sc)],
    yerr=[std(m.tot)/5  std(m.bb)/5  std(m.sc)/5], # just illustrative
    xticks=([1,2,3],["Total","Backbone","Sidechain"]),
    sp=2, label="", xlabel="",
    ylabel="m-value (kJ/mol)",
)
```

Thus, it seems that urea is favoring the evolution of the conformations of the 
protein in time (target structures have more negative transfer free energies than the
reference structure) with the interactions with
the backbone being more relevant than those of the sidechains. The errors are
of course only illustrative.
