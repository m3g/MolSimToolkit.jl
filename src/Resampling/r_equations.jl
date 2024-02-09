#SETTING EQUATIONS FOR PERTUBATION

#Perturbing Lennard-Jones interactions
function L_J(r::Float64, sig1::Float64, ep1::Float64, sig2::Float64, ep2::Float64; alpha::Float64 = 0.0, beta::Float64 = 0.0, smooth::Bool = false, switch::Float64 = 8.0, cut::Float64 = 12.0)
    r = r/ 10    
    switch = switch/10
    cut = cut/10
    epsilon_original = sqrt(ep1*ep2)
    epsilon_perturbed = sqrt(ep1*(1+beta)*ep2)
    sigma_original = (sig1 + sig2)/2
    sigma_perturbed = (sig1*(1 + alpha) + sig2)/2
    function switching_func(r::Float64, switch::Float64, cut::Float64)
        switching = ((cut^2-r^2)^2*(cut^2-r^2-3*(switch^2-r^2)))/(cut^2-switch^2)^3
        return switching
    end
    function V_perturbed(r::Float64) 
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

#Applying a half-gaussian pertubation
function gaussian_decay(r::Float64, alpha::Float64 = 0.0, beta::Float64 = 1.0, cut::Float64 = Inf)
    gaussian = alpha*(r^2 - cut^2)^2*exp(-abs(beta)*r^2)
    if r > cut
        gaussian = 0
    end
    return gaussian
end

#Applying a polynomial decay pertubation
function poly_decay(r::Float64, alpha::Float64 = 0.0, cut::Float64 = Inf)
    poly = alpha*(r^2 - cut^2)^2
    if r > cut
        poly = 0
    end
    return poly
end

#Caculating the distance between two atoms (only for testing)
function dist(r::Float64)
    return r
end