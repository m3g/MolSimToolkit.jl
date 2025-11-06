```@meta
CollapsedDocStrings = true
```

# [m-values (protein transfer free energy) calculator](@id mvalues)

!!! warning
    This is an experimental feature. Breaking changes may occur without 
    a breaking package release.

## Computing m-values for single structure pairs

The core infrastructure to compute *m*-values was moved to [PDBTools.jl](https://m3g.github.io/PDBTools.jl/stable/mvalue/). 
Here we import the relevant functions of PDBTools.jl to compute *m*-values along a trajectory:

```@example mvalues
using MolSimToolkit
#using PDBTools: read_pdb, index, set_position!, 
#                sasa_particles, sasa, mvalue, index, set_position!
using PDBTools
```

Here we show an example of how to compute m-values along a trajectory. 
We define a function to iterate over the simulation frames, and compute the m-values
for each frame, considering the first frame as the reference state:

```@example mvalues
function mvalue_traj(sim::Simulation, protein_ref::AbstractVector{<:PDBTools.Atom})
    # indices of the protein atoms, to fetch frame coordinates
    inds_protein = index.(protein_ref)
    # Create a copy to update coordinates at each frame
    protein_at_frame = copy.(protein_ref)
    # Create arrays to store total and bb and sc contributions
    tot = Float32[]
    bb = Float32[]
    sc = Float32[]
    sasa_reference = sasa_particles(protein_ref)
    for frame in sim
        # fetch protein coordinates and unitcell 
        p = positions(frame)[inds_protein] 
        uc = unitcell(frame)
        # update coordinates (note the dot for broadcast)
        set_position!.(protein_at_frame, p)
        # Compute SASA of the protein in this frame
        sasa_frame = sasa_particles(protein_at_frame; unitcell=uc) 
        # Compute mvalue
        m = mvalue(sasa_reference, sasa_frame, "urea")
        # push total, backbone and sidechain mvalues to arrays
        push!(tot, m.tot)
        push!(bb, m.bb)
        push!(sc, m.sc)
    end
    return tot, bb, sc
end
```

Running the above function over a trajectory can be done with:

```@example mvalues
using MolSimToolkit.Testing # to load test files
# Build Simulation object
sim = Simulation(Testing.namd_pdb, Testing.namd_traj) 
# We are interested only in protein atoms
protein = read_pdb(Testing.namd_pdb, "protein")
# Lets get the positions of the first frame for reference
firstframe!(sim)
p_ref = positions(current_frame(sim))[index.(protein)]
set_position!.(protein, p_ref) # note the dot
# Run the mvalues-traj function
tot, bb, sc = mvalue_traj(sim, protein)
```

Plotting the results, we obtain the $\Delta_{\textrm{ref}\rightarrow\textrm{target}}\Delta G^{T}$, the difference in transfer free energy
from water to urea at 1M, of the target structure at each frame relative to the reference structure:

```@example mvalue_traj
using Plots, Statistics
plt = plot(MolSimStyle, layout=(1,2))
plot!(plt, tot; label="Total", lw=2, sp=1)
plot!(plt, bb; label="Backbone", lw=2, sp=1)
plot!(plt, sc; label="Sidechains", lw=2, sp=1)
plot!(xlabel="frame", ylabel="m-value / (kJ/mol)", sp=1)
bar!(plt, [1  2  3],
    [mean(tot)  mean(bb)  mean(sc)],
    yerr=[std(tot)/5  std(bb)/5  std(sc)/5], # just illustrative
    xticks=([1,2,3],["Total","Backbone","Sidechain"]),
    sp=2, label="", xlabel="", ylims=(-0.08,0),
    ylabel="m-value (kJ/mol)",
)
```

Thus, it seems that urea is favoring the evolution of the conformations of the 
protein in time (target structures have more negative transfer free energies than the
reference structure) with the interactions with
the backbone being more relevant than those of the sidechains. The errors are
of course only illustrative.

