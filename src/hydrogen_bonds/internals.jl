struct HPolarBonds
    D::Vector{Int32} # Hydrogen-bond donnor
    H::Vector{Int32} # polar hydrogen
end
CellListMap.copy_output(x::HPolarBonds) = HPolarBonds(copy(x.D), copy(x.H))
function CellListMap.reset_output!(x::HPolarBonds)
    empty!(x.D)
    empty!(x.H)
    return x
end
function CellListMap.reducer!(x::HPolarBonds, y::HPolarBonds)
    append!(x.D, y.D)
    append!(x.H, y.H)
    return x
end

function hbond_angle(D, H, A)
    v1 = H - D
    v2 = A - D
    return acosd(dot(v1, v2) / (norm(v1) * norm(v2)))
end

function find_hbond_donnors(
    ats_sel;
    positions,
    unitcell,
    cutoff,
    parallel,
    electronegative_elements,
    d_covalent_bond,
)
    sys_polar_bonds = ParticleSystem(;
        positions,
        unitcell,
        cutoff,
        output=HPolarBonds(Int32[], Int32[]),
        parallel,
        output_name=:polar_bonds,
    )
    polar_bonds = map_pairwise!(
        (x, y, i, j, d2, polar_bonds) -> begin
            at_i = ats_sel[i]
            at_j = ats_sel[j]
            el_i = element(at_i)
            el_j = element(at_j)
            if (el_i in electronegative_elements) & (el_j == "H")
                D = i
                H = j
            elseif (el_j in electronegative_elements) & (el_i == "H")
                D = j
                H = i
            else
                return polar_bonds
            end
            if d2 < (d_covalent_bond)^2
                push!(polar_bonds.D, D)
                push!(polar_bonds.H, H)
            end
            return polar_bonds
        end,
        sys_polar_bonds,
    )
    iord = sortperm(polar_bonds.D)
    polar_bonds.D .= polar_bonds.D[iord]
    polar_bonds.H .= polar_bonds.H[iord]
    return polar_bonds
end

function count_hbonds(i, x, y, polar_bonds, positions, unitcell, ang)
    number_of_hbonds = 0
    # Find if i is a donnor
    ii = searchsortedfirst(polar_bonds.D, i)
    ii > length(polar_bonds.D) && return number_of_hbonds
    # Might have more than one polar hydrogen
    while polar_bonds.D[ii] == i
        iH = polar_bonds.H[ii]
        xH = isnothing(unitcell) ? positions[iH] : wrap(positions[iH], x, unitcell)
        hbond_ang = hbond_angle(x, xH, y)
        if hbond_ang <= ang
            number_of_hbonds += 1
        end
        ii += 1
        ii > length(polar_bonds.D) && break
    end
    return number_of_hbonds
end

function count_hbonds2(i, x, y, polar_bonds, positions, unitcell, ang)
    number_of_hbonds = 0
    # Find if i is a donnor
    ii = searchsortedfirst(polar_bonds.D, i)
    ii > length(polar_bonds.D) && return number_of_hbonds
    # Might have more than one polar hydrogen
    while polar_bonds.D[ii] == i
        iH = polar_bonds.H[ii]
        xH = isnothing(unitcell) ? positions[iH] : wrap(positions[iH], x, unitcell)
        hbond_ang = hbond_angle(x, xH, y)
        if hbond_ang <= ang
            number_of_hbonds += 1
        end
        ii += 1
        ii > length(polar_bonds.D) && break
    end
    return number_of_hbonds
end

struct SelectionData{A<:PDBTools.Atom,}
    ats::Vector{A}
    inds::Vector{Int32}
    polar_bonds::HPolarBonds
end
_key_name(sel1, sel2) = "$sel1 => $sel2"