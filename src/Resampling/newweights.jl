"""
    new_weigths(simulation, f_perturbation, group_1, cutoff; prob = true)
    new_weigths(simulation, f_perturbation, group_1, group_2, cutoff; prob = true)

Function that calculates the energy difference when a perturbation between atoms is applied.

It either returns the resulted energy for each one of the frames or their new normalized weights for calculating a system's propriety

The function needs a MolSimToolkit's simulation object, another function to compute the perturbation for each atomic pair i,j and a cut off distance

When prob is selected (it is, by default), the result is frames' new weights and if it is not, the resulted energy

Groups are a vector of atoms indexes, extracted, for example, from a pdb file. 

If you use just one group, the distance between one atom and itself will not be computed!!!

# Example

```julia-repl
julia> import PDBTools

julia> using MolSimToolkit, MolSimToolkit.Resampling

julia> simulation = Simulation("$testdir/Testing_resampling.pdb", "/$testdir/Testing_resampling_one_frame.xtc")

julia> i1 = PDBTools.selindex(atoms(simulation), "resname TFE and name O")

julia> i2 = PDBTools.selindex(atoms(simulation), "residue 11")

julia> sum_of_dist = new_weigths(simulation, (i,j,d2) -> intermol_perturb(i, j, d2, dist), [i1[239]], i2, 25.0)
results([1.0], [5.418206840773324e-33], [74.2955431492777])

julia> sum_of_dist.energy
1-element Vector{Float64}:
 74.2955431492777

This result is the energy difference between the  perturbed frame and the original one. In this case, it is the sum of distances between one oxygen atom of TFE from all other atoms of the 11th protein's residue
```
"""
struct Results
    probability::Vector{Float64}
    relative_probability::Vector{Float64}
    energy::Vector{Float64}
end

function new_weigths(simulation::Simulation, f_perturbation::Function, group_1::Vector{Int64}, cutoff::Float64, k::Float64 = 1.0, T::Float64 = 1.0)
    prob_vec = zeros(length(simulation))
    prob_rel_vec = zeros(length(simulation))
    energy_vec = zeros(length(simulation))
    for (iframe, frame) in enumerate(simulation)
        coordinates = positions(frame)
        first_coors = coordinates[group_1]
        system = PeriodicSystem(
            xpositions = first_coors,
            unitcell = unitcell(frame),
            cutoff = cutoff,
            output = 0.0,
            output_name = :total_energy
        )
        energy_vec[iframe] = map_pairwise!((x, y, i, j, d2, total_energy) -> total_energy + f_perturbation(i, j, d2), system)
    end
    @. prob_rel_vec = exp(-(energy_vec)/k*T)
    prob_vec = prob_rel_vec/sum(prob_rel_vec)
    output = results(prob_vec, prob_rel_vec, energy_vec)
    return output
end

function new_weigths(simulation::Simulation, f_perturbation::Function, group_1::Vector{Int64}, group_2::Vector{Int64}, cutoff::Float64; k::Float64 = 1.0, T::Float64 = 1.0)
    prob_vec = zeros(length(simulation))
    prob_rel_vec = zeros(length(simulation))
    energy_vec = zeros(length(simulation))
    for (iframe, frame) in enumerate(simulation)
        coordinates = positions(frame)
        first_coors = coordinates[group_1]
        second_coors = coordinates[group_2]
        system = PeriodicSystem(
            xpositions = first_coors,
            ypositions = second_coors,
            unitcell = unitcell(frame),
            cutoff = cutoff,
            output = 0.0,
            output_name = :total_energy
        )
        energy_vec[iframe] = map_pairwise!((x, y, i, j, d2, total_energy) -> total_energy + f_perturbation(i, j, d2), system)
    end
    @. prob_rel_vec = exp(-(energy_vec)/k*T)
    prob_vec = prob_rel_vec/sum(prob_rel_vec)
    output = results(prob_vec, prob_rel_vec, energy_vec)
    return output
end

