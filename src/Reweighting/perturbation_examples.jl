#SETTING EQUATIONS FOR PERTUBATION

#Perturbing Lennard-Jones interactions
function switching_func(r, switch, cut)
    switching = ((cut^2-r^2)^2*(cut^2-r^2-3*(switch^2-r^2)))/(cut^2-switch^2)^3
    return switching
end

function lennard_jones(r, sig, ep) 
    V = 4 * ep * ((sig/r)^12 - (sig/r)^6)
    return V
end

function lennard_jones_perturbation(r, sig1, ep1, sig2, ep2; alpha = 0.0, beta = 0.0, smooth::Bool = false, switch = 8.0, cut = 12.0)
    epsilon_original = sqrt(ep1*ep2)
    epsilon_perturbed = sqrt(ep1*(1+beta)*ep2)
    sigma_original = (sig1 + sig2)/2
    sigma_perturbed = (sig1*(1 + alpha) + sig2)/2
    V_perturbed = (lennard_jones(r, sigma_perturbed,  epsilon_perturbed) - lennard_jones(r,  sigma_original, epsilon_original))
    if smooth && switch < r <= cut
        V_perturbed = switching_func(r,switch,cut) * V_perturbed
    else
        V_perturbed =  V_perturbed - (lennard_jones(cut, sigma_perturbed,  epsilon_perturbed) - lennard_jones(cut,  sigma_original, epsilon_original))
    end
    return V_perturbed
end

#Applying a polynomial decay pertubation
poly_decay_perturbation(r, alpha, cut) = r > cut ? zero(r) : alpha*((r/cut)^2 - 1)^2

#Applying a half-gaussian pertubation
gaussian_decay_perturbation(r, alpha, beta) = beta > 0 ? alpha*exp(-beta*r^2) : error("you need to insert positive values for β")