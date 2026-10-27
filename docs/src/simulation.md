```@meta
CollapsedDocStrings = true
```

# Simulation

The `Simulation` object is the central data structure of MolSimToolkit. It holds
the atom data, the trajectory file, and the current frame. Most analysis functions
in MolSimToolkit accept a `Simulation` as their primary argument.

## Iterating over frames

The standard way to process a trajectory is to iterate over all frames with a `for` loop:

```@example simulation_iterate
using MolSimToolkit, MolSimToolkit.Testing, LazyArtifacts

sim = Simulation(test_file("popc_membrane.pdb"), test_file("popc_membrane.dcd"); first=2, step=2, last=4)

for frame in sim
    @show frame_index(sim)
    @show positions(frame)[1]
end
```

The `pairs` iterator yields `(frame_index, frame)` tuples — useful when you
need the actual frame number in the trajectory:

```@example simulation_iterate
for (i, frame) in pairs(sim)
    @show i, frame_index(sim)
end
```

`enumerate` yields `(counter, frame)` tuples — useful when you just need a
sequential count:

```@example simulation_iterate
for (i, frame) in enumerate(sim)
    @show i, frame_index(sim)
end
```

## Example: end-to-end distance

As a simple illustration, we compute the end-to-end distance of the protein (distance between
the C$_\alpha$ atoms of the first and last residues) across all frames, applying periodic
boundary conditions when wrapping the second atom position.

```@example end_to_end
using MolSimToolkit, MolSimToolkit.Testing, LazyArtifacts
using PDBTools
using Plots: plot

sim = Simulation(test_file("popc_membrane.pdb"), test_file("popc_membrane.dcd"))

# Indices of the Cα atoms of the first and last residues
ats = get_atoms(sim)
protein_ca = findall(sel"protein and name CA", ats)
i1, i2 = first(protein_ca), last(protein_ca)

# Iterate over frames and compute end-to-end distances
end_to_end = Float64[]
for frame in sim
    p = positions(frame)
    p1 = p[i1]
    p2 = wrap(p[i2], p1, unitcell(frame))  # apply PBC
    push!(end_to_end, distance(p2,  p1))
end

end_to_end
```

We can plot the result as a function of the frame number:

```@example end_to_end
plot(MolSimStyle,
    frame_range(sim), end_to_end;
    xlabel="frame",
    ylabel="end-to-end distance / Å",
    label=nothing,
)
```

## Creating a Simulation

```@docs
Simulation
Trajectory
```

## Moving around the simulation

The `Simulation` object keeps track of a *current frame*. The functions below
allow you to navigate the trajectory without necessarily iterating over every frame.

```@docs
first_frame!
current_frame
next_frame!
goto_frame!
get_frame!
restart!
set_frame_range!
```

## Querying the simulation

```@docs
frame_range
frame_index
Base.length(::Simulation)
raw_length
get_atoms
path_pdb
path_trajectory
Base.close(::Simulation)
```

## Frames and positions

Each iteration step returns a `Frame` object. The functions below give access to
atomic positions and the unit cell from a frame.

The unit cell of a frame is a `UnitCell` (the unit cell type of PDBTools), or `nothing` if
the trajectory does not contain unit cell information. A fixed unit cell, used for all frames,
can be set with the `unitcell` keyword of the `Simulation` constructor: `unitcell=true` uses
the unit cell of the PDB file (`CRYST1` field), and `unitcell=uc` uses the given `UnitCell`
(see the [PDBTools documentation](https://m3g.github.io/PDBTools.jl/stable/) for its constructors). The atoms returned by `get_atoms` are 
a `PDBTools.Structure` carrying the unit cell of the current frame, which is thus written by 
`write_pdb` and used by PDBTools functions that accept `unitcell=true`.

```@docs
Frame
Point3D
positions
UnitCell
unitcell
```

## Periodic boundary conditions

The functions for periodic boundary conditions (and other functions shared by the packages of
the ecosystem, such as `center_of_mass`, `rmsd` and `align`) are defined in `PDBToolsBase`:

```@docs
PDBToolsBase
```

```@autodocs
Modules = [ PDBToolsBase ]
Pages = [ "wrap.jl" ]
```
## Test files

The examples use the test files of the package, which are obtained with the `test_file`
function of the `MolSimToolkit.Testing` module, available after `using LazyArtifacts`. They
are not part of the package installation: they are downloaded the first time they are used.
