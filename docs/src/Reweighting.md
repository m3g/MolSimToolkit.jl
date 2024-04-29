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

So, secondly, we define some "perturbation" function (here we call it ```gaussian decay```) and set up its parameters:

```julia-repl
julia> gaussian_decay(r, α, β) = α*exp(-abs(β)*r^2)
gaussian_decay (generic function with 1 method)

julia> α = 5.e-3
0.005

julia> β = 5.e-3
0.005

julia> cut_off = 12.0
12.0
```

## Computing the new weights
And finally, using the ```reweight``` function, we pass both the ```simulation``` and the last function anonymously in the input:

```julia-repl
julia> weights = reweight(simulation, (i,j,r) -> gaussian_decay(r, α, β), i1, i2; cutoff = cut_off)
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