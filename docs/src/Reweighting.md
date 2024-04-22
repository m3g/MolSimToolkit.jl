# Simulation Reweighting
Computes the new weight for a frame of a simulation based on the energy difference between the perturbed and non-perturbed original sampling

This resource is based on the Free Energy Perturbation Theory (FEP) in the Molecular Dynamics context. Most of the time, each frame will contribute equally
for calculations of some thermodynamic property, however, we can apply a perturbation on one or multiple types of atomic
interactions (e.g. making water oxygen and protein carbonyl oxygen interaction more repulsive), making these frames to have different normalized statistical contributions so that we can 
possibly preview the outcome of a new simulation with these modifications.

## How to use it
```julia-repl
julia> import Pkg; Pkg.add("MolSimToolkit")

julia> using MolSimToolkit.Reweighting
```

## Setting initial parameters
Firstly, we define the ```simulation``` object and set the atoms that will determine which interactions will be perturbed:

```julia-repl
julia> using MolSimToolkit, PDBTools

julia> testdir = "$(@__DIR__)/test"
"/home/lucasv/.julia/dev/MolSimToolkit/src/Resampling/test"

julia> simulation = Simulation("$testdir/Testing_resampling.pdb", "/$testdir/Testing_resampling_small_trajectory.xtc")
Simulation 
    Atom type: PDBTools.Atom
    PDB file: /home/lucasv/.julia/dev/MolSimToolkit/src/Resampling/test/Testing_resampling.pdb
    Trajectory file: /home/lucasv/.julia/dev/MolSimToolkit/src/Resampling/test/Testing_resampling_small_trajectory.xtc
    Total number of frames: 21
    Frame range: 1:1:21
    Number of frames in range: 21
    Current frame: nothing

julia> i1 = PDBTools.selindex(atoms(simulation), "resname TFE and name O")

julia> i2 = PDBTools.selindex(atoms(simulation), "protein and name O")
```

## Setting perturbation function
In other to obtain these weights, we have to use two functions: the ```reweight```, which will computate each weight and the ```"perturbation"``` function, responsible for taking each computed distance between atomic pairs computated in every frame and calculate the resulted energy using theses distances in that particular frame and applying the desired perturbation.

So, secondly, we define some "perturbation" function (here we call it ```gaussian decay```) and set up its parameters:

```julia-repl
julia> function gaussian_decay(r, α, β) = α*exp(-abs(β)*r^2)
gaussian_decay (generic function with 1 method)

julia> α = 2.e-5
2.e-5

julia> β = 5.e-3
5.e-3

julia> cut_off = 12.0
12.0
```

## Computing the new weights
And finally, using the ```reweight``` function, we pass both the ```simulation``` and the last function anonymously in the input:

```julia-repl
julia> weights = reweight(simulation, (i,j,r) -> gaussian_decay(r, α, β), i1, i2, cut_off)
-------------
FRAME WEIGHTS
-------------

Average probability = 0.047619047619047616
standard deviation = 0.07670453223227784

-------------------------------------------------
FRAME PROBABILITIES RELATIVE TO THE ORIGINAL ONES
-------------------------------------------------

Average probability = 0.0003644497029859351
standard deviation = 0.0005870538237843036

----------------------------------
COMPUTED ENERGY AFTER PERTURBATION
----------------------------------

Average energy = 9.439727822732474
standard deviation = 2.219086818369851
```

The data in ```weights``` structure is organized as it follows:

```julia
struct Resampling_results
    probability::Vector{Float64}
    relative_probability::Vector{Float64}
    energy::Vector{Float64}
end
```

As an example, if we want the absolute weights computed for our simulation:

```julia-repl
julia> weights.probability
21-element Vector{Float64}:
 0.0070439494907980245
 0.00010071690416721846
 0.009708829320345114
 0.03599823389430867
 0.028511157971421924
 0.013446105873515246
 0.00024718128817611135
 0.14415904378968805
 0.0008117791731669824
 0.30217604446677
 0.00046609570679062177
 0.1505884065155757
 0.04690650533601024
 0.008499321313151335
 0.0011392562627410294
 0.04458432806281825
 0.01706069848349761
 0.006018368435226036
 0.02730599234002224
 0.004453135199832921
 0.15077485017197673
```

## Reference functions
```@autodocs
Modules = [MolSimToolkit.Reweighting]
Pages = ["reweighting.jl"]
Order = [:function, :type]
```