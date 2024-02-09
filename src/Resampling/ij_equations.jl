function intermol_perturb(i::Int64, j::Int64, d2::Float64, perturbation::Function; pairwise_check::Bool=false, n_atoms_per_molecule::Int64=1)
    energy = perturbation(sqrt(d2))
    return energy
end