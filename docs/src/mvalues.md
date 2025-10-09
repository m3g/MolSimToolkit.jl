```@meta
CollapsedDocStrings = true
```

# [m-values (protein transfer free energy) calculator](@id mvalues)

!!! warning
    This is an experimental feature. Breaking changes may occur without 
    a breaking package release.

## Computing m-values for single structure pairs

```@docs
mvalue
delta_sasa_per_restype
parse_mvalue_server_sasa
gmx_delta_sasa_per_restype
```

## Computing m-values along a trajectory

Here we show an example of how to compute m-values along a trajectory. 

```@example mvalue_traj
using MolSimToolkit, PDBTools
using MolSimToolkit.Testing
# Build Simulation object
sim = Simulation(Testing.namd_pdb, Testing.namd_traj) 
# We are interested only in protein atoms
protein = read_pdb(Testing.namd_pdb, "protein")
inds_protein = index.(protein)
# Create a copy to update coordinates at each frame
protein_at_frame = copy.(protein)
# Lets get the positions of the first frame for reference
firstframe!(sim)
p_ref = positions(current_frame(sim))[inds_protein]
set_position!.(protein, p_ref);
```
Now we will iterate over the simulation frames, and compute the m-values
for each frame, considering the first frame as the reference state:
```@example mvalue_traj
# Create arrays to store total and bb and sc contributions
mvalues_tot = Float32[]
mvalues_bb = Float32[]
mvalues_sc = Float32[]
for frame in sim
    # fetch protein coordinates and unitcell 
    p = positions(frame)[inds_protein] 
    uc = unitcell(frame)
    # update coordinates (note the dot for broadcast)
    set_position!.(protein_at_frame, p)
    # Compute variations in SASA relative to reference
    dsasa = delta_sasa_per_restype(;
        native=protein,
        desnat=protein_at_frame,
    )
    # Compute mvalue
    m = mvalue(;
        model=MoeserHorinek,
        cosolvent="urea",
        atoms=protein,
        sasas=dsasa,
    )
    # push total, backbone and sidechain mvalues to arrays
    push!(mvalues_tot, m.tot)
    push!(mvalues_bb, m.bb)
    push!(mvalues_sc, m.sc)
end
```

Now, lets plot the results:

```@example mvalue_traj
using Plots, Statistics
plt = plot(MolSimStyle, layout=(1,2))
plot!(plt, mvalues_tot; label="Total", lw=2, sp=1)
plot!(plt, mvalues_bb; label="Backbone", lw=2, sp=1)
plot!(plt, mvalues_sc; label="Sidechains", lw=2, sp=1)
plot!(xlabel="frame", ylabel="m-value / (kJ/mol)", sp=1)
bar!(plt, [1  2  3],
    [mean(mvalues_tot)  mean(mvalues_bb)  mean(mvalues_sc)],
    yerr=[std(mvalues_tot)/5  std(mvalues_bb)/5  std(mvalues_sc)/5],
    xticks=([1,2,3],["Total","Backbone","Sidechain"]),
    sp=2, label="", xlabel="", ylims=(-0.08,0),
    ylabel="m-value (kJ/mol)",
)
```

Thus, it seems that urea is favoring the evolution of the conformations of the 
protein in time (negative transfer free energies) with the interactions with
the backbone being more relevant than those of the sidechains. The errors are
of course only illustrative.

