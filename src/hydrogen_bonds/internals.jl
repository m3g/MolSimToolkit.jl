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

#=

    process_selections(selections) -> Vector{Pair{String,String}}

Convert input selections into a list of pairs. If no selection is provided, returns ["all"=>"all"].

=#
function process_selections(selections)
    if isnothing(first(selections))
        return ["all" => "all"]
    end
    return [sel isa String ? sel => sel : first(sel) => last(sel) for sel in selections]
end

#=

    initialize_hbonds_data(sim, selection_pairs) -> Tuple{OrderedDict, Dict}

Initialize the results dictionary and process selection data for each unique selection.

=#
function initialize_hbonds_data(
    sim::Simulation,
    selection_pairs::Vector{Pair{String,String}};
    parallel::Bool=true,
    donnor_acceptor_distance::Real=3.5f0,
    electronegative_elements=("N", "O", "F", "S"),
    d_covalent_bond::Real=1.2f0,
)
    hbonds = OrderedDict{String,Vector{Int}}()
    selection_data = Dict{String,SelectionData}()

    p_first_frame = positions(current_frame(sim))
    uc_first_frame = unitcell(current_frame(sim))

    for selection_pair in selection_pairs
        sel1, sel2 = first(selection_pair), last(selection_pair)
        hbonds[_key_name(sel1, sel2)] = zeros(Int, length(sim))

        for sel in selection_pair
            if !haskey(selection_data, sel)
                ats_sel = select(atoms(sim), sel)
                inds_sel = index.(ats_sel)
                polar_bonds = find_hbond_donnors(
                    ats_sel;
                    positions=@view(p_first_frame[inds_sel]),
                    unitcell=uc_first_frame.matrix,
                    cutoff=donnor_acceptor_distance,
                    parallel=parallel,
                    electronegative_elements=electronegative_elements,
                    d_covalent_bond=d_covalent_bond,
                )
                selection_data[sel] = SelectionData(ats_sel, inds_sel, polar_bonds)
            end
        end
    end

    return hbonds, selection_data
end

#=

    setup_particle_systems(selection_pairs, selection_data, uc_first_frame, donnor_acceptor_distance, parallel) 
        -> Dict{String,Union{ParticleSystem1,ParticleSystem2}}


Initialize CellListMap particle systems for each selection pair.

=#
function setup_particle_systems(
    selection_pairs,
    selection_data,
    uc_first_frame,
    donnor_acceptor_distance,
    parallel
)
    systems = Dict{String,Union{CellListMap.ParticleSystem1,CellListMap.ParticleSystem2}}()

    for selection_pair in selection_pairs
        sel1, sel2 = first(selection_pair), last(selection_pair)
        key = _key_name(sel1, sel2)

        if sel1 == sel2
            s1 = selection_data[sel1]
            systems[key] = ParticleSystem(
                positions=coor.(s1.ats),
                unitcell=uc_first_frame.matrix,
                cutoff=donnor_acceptor_distance,
                output=0,
                parallel=parallel,
                output_name=:number_of_hbonds
            )
        else
            s1, s2 = selection_data[sel1], selection_data[sel2]
            systems[key] = ParticleSystem(
                xpositions=coor.(s1.ats),
                ypositions=coor.(s2.ats),
                unitcell=uc_first_frame.matrix,
                cutoff=donnor_acceptor_distance,
                output=0,
                parallel=parallel,
                output_name=:number_of_hbonds
            )
        end
    end

    return systems
end

#=
    count_hbonds_same_selection(sys, s1, angle_cutoff, electronegative_elements) -> Int

Count hydrogen bonds when donor and acceptor are in the same selection.
=#
function count_hbonds(sys::CellListMap.ParticleSystem1, s1, angle_cutoff, electronegative_elements)
    nhb = map_pairwise!(sys) do x, y, i, j, _, number_of_hbonds
        el_i = element(s1.ats[i])
        el_j = element(s1.ats[j])
        if (el_i in electronegative_elements) & (el_j in electronegative_elements)
            number_of_hbonds += count_hbonds(
                i, x, y, s1.polar_bonds,
                sys.positions, sys.unitcell,
                angle_cutoff,
            )
            number_of_hbonds += count_hbonds(
                j, y, x, s1.polar_bonds,
                sys.positions, sys.unitcell,
                angle_cutoff,
            )
        end
        return number_of_hbonds
    end
    return nhb
end

#=
    count_hbonds(sys, s1, s2, sel1, sel2, angle_cutoff, electronegative_elements) -> Int

Count hydrogen bonds between two different selections.
=#
function count_hbonds(sys::CellListMap.ParticleSystem2, s1, s2, sel1, sel2, angle_cutoff, electronegative_elements)
    nhb = map_pairwise!(sys) do x, y, i, j, _, number_of_hbonds
        at_i = s1.ats[i]
        at_j = s2.ats[j]
        if index(at_i) == index(at_j)
            throw(ArgumentError("""\n
                Different selections cannot overlap. Detected atom $(index(at_i)) in both selections \"$sel1\" and \"$sel2\".
                """))
        end
        el_i = element(at_i)
        el_j = element(at_j)
        if (el_i in electronegative_elements) & (el_j in electronegative_elements)
            number_of_hbonds += count_hbonds2(
                i, x, y, s1.polar_bonds, sys.xpositions,
                sys.unitcell, angle_cutoff
            )
            number_of_hbonds += count_hbonds2(
                j, y, x, s2.polar_bonds, sys.ypositions,
                sys.unitcell, angle_cutoff
            )
        end
        return number_of_hbonds
    end
    return nhb
end