# Resampling

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

In other to obtain these weights, we have to use three function: the ```new_weights``` , which will computate each weight, an ```"ij"``` function to apply the 
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

And finally, using the ```new_weights``` function, we pass both the ```simulation``` and the last function anonymously in the input:

```julia-repl
julia> alpha = 2.e-5
2.e-5

julia> beta = 5.e-3
5.e-3

julia> cut_off = 12.0
12.0

julia> weights = new_weigths(simulation, (i,j,d2) -> intermol_perturb(i, j, d2, r -> gaussian_decay(r, alpha, beta, cut_off)), i1, i2, cut_off)
Results([0.1682381803849657, 0.007950747824893706, 0.034121498541249594, 0.016285477332538448, 0.04632343416825627, 0.0029827061482835906, 0.003308362899210579, 0.1698176067076115, 0.0003745981041899426, 0.015456445688650771  …  0.0008011931617518091, 0.003845100391591656, 0.023532097284969005, 0.07521218093097595, 0.10899848267492684, 0.002026112721413265, 0.09806990189406785, 0.015721263790013106, 0.0013258115999383426, 0.03860218054966671], [0.003368421475812837, 0.00015918782324477303, 0.0006831718472600203, 0.0003260636287509907, 0.0009274758567211686, 5.971897355836635e-5, 6.623918907101538e-5, 0.0034000443424682565, 7.500106670583892e-6, 0.0003094649709016859  …  1.6041282936749778e-5, 7.698560877237773e-5, 0.0004711535852576899, 0.0015058788968764902, 0.002182339520274834, 4.056630657565006e-5, 0.001963530293271981, 0.0003147670906569741, 2.6545058059324903e-5, 0.0007728829073113398], [5.693311049156978, 8.745425774635232, 7.288764123511054, 8.028418015424984, 6.983043794310789, 9.72586077301021, 9.622238290247603, 5.683966805542638, 11.800593314778633, 8.080665651884063  …  11.04034497513653, 9.471892052638426, 7.660326433751817, 6.498378566598426, 6.127357803183168, 10.112572723202948, 6.233011255864571, 8.063677587365515, 10.53666696472723, 7.165382999102666])
```

Now, with this ```weights```structure, we can extract three results: each frame's weight, using ```weights.probability```, each frame's weights relative to its original
probability using ```weights.relative_probability```, and last, but not least, the energy difference for each frame, using ```weights.energy```.

```julia-repl
julia> weights.prob
21-element Vector{Float64}:
 0.1682381803849657
 0.007950747824893706
 0.034121498541249594
 0.016285477332538448
 0.04632343416825627
 0.0029827061482835906
 0.003308362899210579
 0.1698176067076115
 0.0003745981041899426
 0.015456445688650771
 0.16700661720083548
 0.0008011931617518091
 0.003845100391591656
 0.023532097284969005
 0.07521218093097595
 0.10899848267492684
 0.002026112721413265
 0.09806990189406785
 0.015721263790013106
 0.0013258115999383426
 0.03860218054966671
```