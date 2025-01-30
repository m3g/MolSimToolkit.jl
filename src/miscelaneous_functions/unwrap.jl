


function unrwap(positions, positions_last, unitcell)
    for i in eachindex(positions, positions_last)
        positions[i] = wrap(positions[i], positions_last[i], unitcell)
    end
    return positions
end

function msd(sim::Simulation, inds)
    d = zeros(length(sim))
    first_frame!(sim)
    p0 = positions(current_frame(sim))[inds]
    plast = copy(p0)
    for iframe in 2:length(sim)
        f = next_frame!(sim)
        p = positions(f)[inds]
        for i in 1:length(inds)
            pat = wrap(p[i], plast[i], unitcell(f))
            d[iframe] += sum(abs2, pat - p0[i])
        end
        d[iframe] /= length(inds)
        plast .= p  
    end
    return d
end

@testitem "displacement" begin

    dtest = zeros(length(sim1))
    for (iframe, frame) in enumerate(sim1)
        p = positions(current_frame(sim1))[22211]
        dtest[iframe] = sum(abs2, p - p0)
    end
    
end