struct ContactMap{T}
    matrix::Matrix{Union{Missing,T}}
end
Base.zero(::Type{ContactMap{T}}, n, m) where {T} = ContactMap{T}([missing for i in 1:n, j in 1:m])
Base.setindex!(map::ContactMap, value, i, j) = map.matrix[i, j] = value
Base.show(io::IO, ::MIME"text/plain", map::ContactMap) = show(io, map.matrix)
contact_matrix(map::ContactMap) = map.matrix

function _residue_contacts(
    p::AbstractVector{<:AbstractVector{T}}, 
    residues::AbstractVector{<:PDBTools.Residue}; 
    dmax::Real=7.0,
    gap=3,
    unitcell=nothing
) where {T}
    map = zero(ContactMap{T}, length(residues), length(residues))
    for ires in eachindex(residues), jres in ires+gap:length(residues)
        r1 = residues[ires]
        r2 = residues[jres]
        for (iat, jat) in Iterators.product(eachindex(r1),eachindex(r2))
            p_i = p[PDBTools.index(r1[iat])]
            p_j = p[PDBTools.index(r2[jat])]
            p_j = !isnothing(unitcell) ? wrap(p_j,p_i,unitcell) : p_j
            d = norm(p_j - p_i)
            if d < dmax 
                map[ires, jres] = d
                map[jres, ires] = d
                break
            end
        end
    end
    return map
end

function residue_contacts(sim::Simulation, atoms; dmax=7.0, ref=1)
    inds = index.(atoms)
    if ref isa Integer
        if ref > length(sim)
            throw(ArgumentError("Reference structure index is out of simulation bounds"))
        end
        for (iframe, frame) in enumerate(sim)
            if iframe == ref
                xref = coor(positions(frame))[inds]
            end
        end
    elseif ref == ":pdb"
        xref = coor(atoms(sim))[inds]
    else
        throw(ArgumentError("Reference structure must be an integer or :pdb"))
    end
    residues = collect(eachresidue(atoms))
#   ...
end