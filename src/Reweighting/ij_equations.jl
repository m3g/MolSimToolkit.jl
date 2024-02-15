function intermol_perturb(i::Int64, j::Int64, d2::Float64, perturbation::Function)
    energy = perturbation(sqrt(d2))
    return energy
end