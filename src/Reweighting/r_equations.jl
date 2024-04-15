#SETTING EQUATIONS FOR PERTUBATION

#Perturbing Lennard-Jones interactions
function L_J(r, sig1, ep1, sig2, ep2; alpha = 0.0, beta = 0.0, smooth::Bool = false, switch = 8.0, cut = 12.0)
    r = r/ 10    
    switch = switch/10
    cut = cut/10
    epsilon_original = sqrt(ep1*ep2)
    epsilon_perturbed = sqrt(ep1*(1+beta)*ep2)
    sigma_original = (sig1 + sig2)/2
    sigma_perturbed = (sig1*(1 + alpha) + sig2)/2
    function switching_func(r, switch, cut)
        switching = ((cut^2-r^2)^2*(cut^2-r^2-3*(switch^2-r^2)))/(cut^2-switch^2)^3
        return switching
    end
    function V_perturbed(r) 
        V_result = 4 * epsilon_perturbed * ((sigma_perturbed/r)^12 - (sigma_perturbed/r)^6) - 4 * epsilon_original*((sigma_original/r)^12 - (sigma_original/r)^6)
        return V_result 
    end
    if smooth
        if r > cut
            V_perturbed = 0.0
        elseif switch < r <= cut
            V_perturbed = switching_func(r,switch,cut)* V_perturbed(r)
        else
            V_perturbed = V_perturbed(r)
        end
    else
        if r > cut
            V_perturbed = 0.0
        else
            V_perturbed = V_perturbed(r) - V_perturbed(cut)
        end
    end
    return V_perturbed
end

#Applying a polynomial decay pertubation
poly_decay(r, alpha, cut) = r > cut ? zero(r) : alpha*((r/cut)^2 - 1)^2

#Applying a half-gaussian pertubation
gaussian_decay(r, alpha, beta) = beta > 0 ? alpha*exp(-beta*r^2) : error("you need to insert positive values for β")