# Simulation Reweighting
Computes the new weight for a frame of a simulation based on the energy difference between the perturbed and non-perturbed original sampling

This resource is based on the Free Energy Perturbation Theory (FEP) in the Molecular Dynamics context. Most of the time, each frame will contribute equally
for calculations of some thermodynamic property, however, we can apply a perturbation on one or multiple types of atomic
interactions (e.g. making water oxygen and protein carbonyl oxygen interaction more repulsive), making these frames to have different normalized statistical contributions so that we can 
possibly preview the outcome of a new simulation with these modifications.

## How to use it
```julia-repl
julia> using MolSimToolkit.Reweighting
```

## Setting initial parameters
Firstly, we define the ```simulation``` object and set the atoms that will determine which interactions will be perturbed:

```julia-repl
julia> using MolSimToolkit, PDBTools

julia> testdir = "$(@__DIR__)/test"
"/home/lucasv/.julia/dev/MolSimToolkit/src/Reweighting/test"

julia> simulation = Simulation("$testdir/Testing_reweighting.pdb", "/$testdir/Testing_reweighting_10_frames_trajectory.xtc")
Simulation 
    Atom type: PDBTools.Atom
    PDB file: /home/lucasv/.julia/dev/MolSimToolkit/src/Resampling/test/Testing_resampling.pdb
    Trajectory file: /home/lucasv/.julia/dev/MolSimToolkit/src/Resampling/test/Testing_reweighting_10_frames_trajectory.xtc
    Total number of frames: 10
    Frame range: 1:1:10
    Number of frames in range: 10
    Current frame: nothing

julia> i1 = PDBTools.selindex(atoms(simulation), "resname TFE and name O")

julia> i2 = PDBTools.selindex(atoms(simulation), "protein and name O")
```

## Setting perturbation function
In other to obtain these weights, we have to use two functions: the ```reweight``` function, which will calculate each weight and the ```perturbation``` function, responsible for taking each computated distance between atomic pairs in every frame and determine the resulted energy using theses distances in that particular frame based on the applied perturbation.

So, secondly, we define some "perturbation" function (here we call it ```gaussian decay```) and set up its parameters. Please, take a look at the interface:

```julia-repl
julia> gaussian_decay(r, α, β) = α*exp(-abs(β)*r^2)
gaussian_decay (generic function with 1 method)

julia> α = 5.e-3
0.005

julia> β = 5.e-3
0.005
```

As it can be seen, the function has to receive two parameters: `r` which corresponds to the distance between two selected atoms and some parameter to account a modification and change its magnitude, here, we inserted two of them in the same function `α` to change the maximum value of the gaussian curve and `β` to adjust its decay behaviour with a given value of `r`.

## Computing the new weights
And finally, using the ```reweight``` function, we pass both the ```simulation``` and the last function anonymously in the input. Again, watch the interface:

```julia-repl
julia> cut_off = 12.0
12.0

julia> weights = reweight(simulation, (i,j,r) -> gaussian_decay(r, α, β), i1, i2; cutoff = cut_off)
```

`i and j`: if you selected two atom types, `i` will be the index for either the first, the second, the third and so on up to the last atom of the first group and `j` will be same, but now for the second one. With these two parameters, it is possible to determine every combination of two atoms, each one coming from one group, and compute the associated dsitance `r`, so that we are taking into account all interactions between these two atom types to our perturbation. However, if we are dealing with just one group, both of them are indexes for all the atoms of the selected group. Bear this is mind because *it is possible to compute repeated combinations* (like `i,j = 1,2 or 2,1`), so your `perturbation function` ought to be able to avoid this!

`r`: the distance between the twos atoms with indexes `i` and `j` in the selected groups

`cutoff`: the maximum distance that will be computed between two atoms. The default value is `12.0` Angstrom

```julia-repl
-------------
FRAME WEIGHTS
-------------

Average probability = 0.1
standard deviation = 0.011364584999859616

-------------------------------------------------
FRAME WEIGHTS RELATIVE TO THE ORIGINAL ONES
-------------------------------------------------

Average probability = 0.6001821184861403
standard deviation = 0.06820820700931557

----------------------------------
COMPUTED ENERGY AFTER PERTURBATION
----------------------------------

Average energy = 0.5163045415662408
standard deviation = 0.11331912115883522
```

The data in ```weights``` structure is organized as it follows:

```julia
struct ReweightResults
    probability::Vector{Float64}
    relative_probability::Vector{Float64}
    energy::Vector{Float64}
end
```

As an example, if we want the absolute weights computed for our simulation:

```julia-repl
julia> weights.probability
10-element Vector{Float64}:
 0.0946703156089622
 0.08132904424357772
 0.09869474081686125
 0.10655562666294487
 0.10022245896670733
 0.09564708264572024
 0.08975220271835752
 0.11576085753347072
 0.09825654640701757
 0.11911112439638061
```

## Reference Functions
```@autodocs
Modules = [MolSimToolkit.Reweighting]
Order = [:function, :type]
```