# Secondary structures

This package provides convenience functions to analyze the secondary structure along
molecular dynamics simulations. 

## Secondary structure map

The secondary structure map is the profile of the secondary structure computed for 
each frame of the trajectory. This computation may be costly, particularly with the 
DSSP algorithm, so it is recommended to save the result. See [Saving and loading a map](@ref)
for further information. 

```@docs
ss_map
```

A complete example for computing a secondary structure map is shown below:

```jldoctest
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> ssmap = ss_map(simulation; selection="residue >= 30 and residue <= 35", show_progress=false)
6×5 Matrix{Int64}:
 5  9  5  5  5
 5  9  5  5  5
 5  1  5  5  5
 5  1  5  5  5
 5  1  5  5  5
 9  9  9  9  9

```

Here we have computed the secondary structure map for 6 residues of the structure, along the
5 frames of the trajectory. The resulting map is a matrix, where each code represents 
a different class of secondary structure. The conversion between representations of 
the classes can be done with these three functions of the [ProteinSecondaryStructures.jl](https://github.com/m3g/ProteinSecondaryStructures.jl) package,
which are reexported here:

- `ss_code`: convert the representation to one-letter codes like `H`, `B`, `C`, etc.
- `ss_name`: convert the representation to secondary structure names like `Alpha-helix`, `Beta-bridge`, etc.
- `ss_number`: convert the representation to code numbers, like the ones used the matrix above. 

The list of classes and code associations of is available 
[here, in the ProteinSecondaryStructures.jl documentation](https://m3g.github.io/ProteinSecondaryStructures.jl/stable/overview/#Secondary-structure-classes).

For example, considering the secondary structure map matrix above, we can do:

```julia-repl
julia> ss_name.(ssmap)
6×5 Matrix{String}:
 "turn"  "coil"       "turn"  "turn"  "turn"
 "turn"  "coil"       "turn"  "turn"  "turn"
 "turn"  "310 helix"  "turn"  "turn"  "turn"
 "turn"  "310 helix"  "turn"  "turn"  "turn"
 "turn"  "310 helix"  "turn"  "turn"  "turn"
 "coil"  "coil"       "coil"  "coil"  "coil"

```

## Plotting the map

The secondary structure along the trajectory can be plotted as heatmap:

```julia-repl
julia> using Plots

julia> heatmap(ssmap, xlabel="frame", ylabel="residue",
         framestyle=:box,
         color=palette(:tab20c,10)
       )
```

The command above will produce a raw heatmap with the data, and can be customized. We provide,
nevertheless, a recipe for constructing secondary structure heatmaps, through the `ss_heatmap`
function:

```julia-repl
julia> using Plots



```

![heatmap1](./images/secondary_structure/heatmap1.svg)

## Calculation methods: STRIDE and DSSP

```julia
ssmap = ss_map(atoms, trajectory; method=stride_run)
ssmap = ss_map(atoms, trajectory; method=dssp_run)
```

## Saving and loading a map

The secondary structure map computed is just a matrix of integer codes. Thus, it can be saved or read in any preferred format.
As a suggestion, it is possible to use `writedlm` and `readdlm` function from the `DelimitedFiles` package: 

```julia
using DelimitedFiles
# save data to ssmap.dat
writedlm("ssmap.dat", ssmap)
# load data
ssmat = readdlm("ssmap.dat", Int)
```

## Trajectory secondary structure classes

From a precomputed secondary structure map, or from a trajectory, helper functions
will provide the content of a specific call of secondary structure along the simulation:

```@docs
ss_mean
```
### From the secondary structure map

Calling `ss_content` with a class identifier function and a map (as computed above), will return the content
of that class along the trajectory:

```julia-repl
julia> ss_content(is_alphahelix, ssmap)
26-element Vector{Float64}:
 0.21052631578947367
 0.15789473684210525
 ⋮
 0.13157894736842105
```

The composition of classes for a given frame can also be retrieved from the content map:

```julia
julia> ss_composition(ssmap, 6)
Dict{String, Int64} with 10 entries:
  "310 helix"   => 7
  "bend"        => 0
  "turn"        => 17
  "kappa helix" => 0
  "beta strand" => 25
  "beta bridge" => 2
  "alpha helix" => 12
  "pi helix"    => 0
  "loop"        => 0
  "coil"        => 13
```

These functions are useful, because the computation of the secondary structure along the
trajectory (the map) can be costly.

### Single class, along the trajectory

If the user wants to compute the content of a single class of secondary structures
along a trajectory, that can be done without precomputing the secondary structure map
(note, however, that the cost is similar).

For example, in the following script we compute the content of $\alpha$-helices of the
structure along the trajectory:

```julia
using ProteinSecondaryStructures 
using PDBTools: readPDB 
using Chemfiles: Trajectory

pdbfile = ProteinSecondaryStructures.Testing.data_dir*"/Gromacs/system.pdb"
trajectory_file = ProteinSecondaryStructures.Testing.data_dir*"/Gromacs/trajectory.xtc"

atoms = readPDB(pdbfile, "protein")
trajectory = Trajectory(trajectory_file)

helical_content = ss_content(is_alphahelix, atoms, trajectory)
```

The method to compute the secondary structure can be defined with the `method`
keyword: 

```julia
helical_content = ss_content(is_alphahelix, atoms, trajectory; method=stride_run)
#or
helical_content = ss_content(is_alphahelix, atoms, trajectory; method=dssp_run)
```

## Average structure per residue

Here we provide a example where we use some features of `PDBTools.jl` and `Plots`
to illustrate the average content of $\alpha$-helices for each residue
of the protein, along the simulation. 

Here, we assume that a secondary structure map, `ssmap`, was computed using the
instructions above.

The goal is to obtain a figure similar to this one, in which in the upper pannel we
show the evolution of the total $\alpha$-helical content as a function of the simulation
frames, and in the lower pannel we show the content of helices of each residue, with
appropriate indexing. 


The script to produce the figure above is a manipulation of the `ssmap` output, using
function from `PDBTools` and the plotting features of `Plots`. THe complete script is:

```julia
using Plots, PDBTools
Plots.default(fontfamily="Computer Modern",linewidth=2,framestyle=:box)
plt = plot(layout=(2,1))
ahelix = ss_content(is_alphahelix, ssmap)
plot!(plt, subplot=1, 
    ahelix, label=nothing,
    xlabel="simulation frame", 
    ylabel="α-helical content"
 )
residue_indexes=1:length(eachresidue(atoms))
one_letter_names = eachresidue(atoms) .|> resname .|> oneletter;
string_numbers = eachresidue(atoms) .|> resnum .|> string;
xlabels = one_letter_names .* string_numbers
ahelix_avg = map(mean, eachrow(is_alphahelix.(ssmap)))
xticks=(residue_indexes[begin:5:end],xlabels[begin:5:end])
plot!(plt, subplot=2,
    residue_indexes, 
    ahelix_avg, 
    label=nothing,
    xlabel="Residue",
    ylabel="α-helical content",
    xticks=xticks, xrotation=60
)
savefig("./helical_content.svg")
```

#### Step-by-step construction of the figure

First, we load the `Plots` and `PDBTools` packages, and set some default
parameters for `Plots` for prettier output.

```julia-repl
julia> using Plots, PDBTools

julia> Plots.default(fontfamily="Computer Modern",linewidth=2,framestyle=:box)
```

We then initialize a plot with two pannels. The upper supblot will contain the
$\alpha$-helical content as a function simulation frames, and the lower subplot
will contain the average content of helices for each residue.

```julia-repl
julia> plt = plot(layout=(2,1))
```

Next, we compute, from the secondary structure map, the $\alpha$-helical content,
for each frame of the trajectory, which will be printed in the first subplot of the figure:

```julia-repl
julia> ahelix = ss_content(is_alphahelix, ssmap)

julia> plot!(plt, subplot=1, 
           ahelix, label=nothing, 
           xlabel="simulation frame", 
           ylabel="α-helical content"
        )
```

For the second plot, we first define a residue range, with the number of residues of
the protein, using [`PDBTools.eachresidue`](https://m3g.github.io/PDBTools.jl/stable/selections/#Iterate-over-residues-(or-molecules)) iterator.
Here, `length(eachresidue(atoms))` is just the number of residues of the protein:

```julia-repl
julia> residue_indexes=1:length(eachresidue(atoms))
1:76
```

We the extract the names of all residues, which we will use for creating the `x`-labels of our plot. We
iterate over all residues first to extract their names, which are converted to *one-letter* codes, and these
are concateneted (with the `*` operation on strings), with the residue numbers converted to strings: 

```julia-repl
julia> one_letter_names = eachresidue(atoms) .|> resname .|> oneletter;

julia> string_numbers = eachresidue(atoms) .|> resnum .|> string;

julia> xlabels = one_letter_names .* string_numbers
76-element Vector{String}:
 "M1"
 "Q2"
 ⋮
 "G76"
```

The `y`-axis of our plot will contain the average $\alpha$-helical content for each residue.
To extract that, we will first convert the `ssmap` to matrix of `0`s and `1`s, with the
broadcast of the `is_alphahelix` function:

```julia-repl
julia> is_alphahelix.(ssmap)
76×26 BitMatrix:
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 ⋮              ⋮              ⋮              ⋮              ⋮              ⋮
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
```

The average of each row is the the average content of helices for each residue:
```julia-repl
julia> ahelix_avg = map(mean, eachrow(is_alphahelix.(ssmap)))
76-element Vector{Float64}:
 0.0
 0.0
 ⋮
 0.0
```

We can finally plot the second supblot of our figure, with the note that we have filtered
some `x`-tick labels to avoid having a crowded axis:

```julia-repl
julia> plot!(plt, subplot=2,
           residue_indexes, 
           ahelix_avg, 
           label=nothing,
           xlabel="Residue",
           ylabel="α-helical content",
           xticks=(residue_indexes[begin:5:end], xlabels[begin:5:end]),
           xrotation=60,
       )

julia> savefig("./helical_content.svg")
```

The final line saves the figure to an external file.



