# Simulation Reweighting

Computes the new weight for a frame of a simulation based on the energy difference between the perturbed and non-perturbed original sampling

This resource is based on the Free Energy Perturbation Theory (FEP) in the Molecular Dynamics context. Most of the time, each frame will contribute equally
for calculations of some thermodynamic property (e.g. Kirkwood-Buff Integrals). However, we can apply a perturbation on one or multiple types of atomic
interactions (e.g. making then more repulsive or attractive), making these frames to have different normalized statistical contributions so that we can 
possibly preview the outcome of a new simulation with these modifications.

## How to use it

```julia-repl
julia> import Pkg; Pkg.add("MolSimToolkit")

julia> using MolSimToolkit.Resampling
```

## Setting functions

In other to obtain these weights, we have to use three function: the ```reweight``` , which will computate each weight, an ```"ij"``` function to apply the 
```"perturbation"``` function to each distance between atomic pairs computated in the frame. The last one is responsible for calculate the resulted energy for that distance in that frame, applying the desired perturbation.

Firstly, we define the ```simulation``` object and set the atoms that will determine which interactions will be perturbed:

```julia-repl
julia> using MolSimToolkit, PDBTools, CellListMap.PeriodicSystems

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

Secondly, we define some "perturbation" function (```gaussian decay```):

```julia-repl
julia> function gaussian_decay(r::Float64, alpha::Float64, beta::Float64, cut::Float64)
            gaussian = alpha*(r^2 - cut^2)^2*exp(-abs(beta)*r^2)
            if r > cut
               gaussian = 0
            end
            return gaussian
       end
gaussian_decay (generic function with 1 method)
```

Then, we also define some "ij" function (```intermol_perturb```):

```julia-repl
julia> function intermol_perturb(i::Int64, j::Int64, d2::Float64, perturbation::Function)
           energy = perturbation(sqrt(d2))
           return energy
       end
intermol_perturb (generic function with 1 method)
```

And finally, using the ```reweight``` function, we pass both the ```simulation``` and the last function anonymously in the input:

```julia-repl
julia> alpha = 2.e-5
2.e-5

julia> beta = 5.e-3
5.e-3

julia> cut_off = 12.0
12.0

julia> weights = reweight(simulation, (i,j,d2) -> intermol_perturb(i, j, d2, r -> gaussian_decay(r, alpha, beta, cut_off)), i1, i2, cut_off)
-------------
FRAME WEIGHTS
-------------

Average probability = 0.047619047619047616
standard deviation = 0.07670453223227784
Mode = 0.0070439494907980245

-------------------------------------------------
FRAME PROBABILITIES RELATIVE TO THE ORIGINAL ONES
-------------------------------------------------

Average probability = 0.0003644497029859351
standard deviation = 0.0005870538237843036
Mode = 5.391047129515465e-5

----------------------------------
COMPUTED ENERGY AFTER PERTURBATION
----------------------------------

Average energy = 9.439727822732474
standard deviation = 2.219086818369851
```

The data in ```weights``` structure is organized as it follows:

```julia-repl
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
Pages = ["reweight.jl"]
Order = [:function, :type]
```