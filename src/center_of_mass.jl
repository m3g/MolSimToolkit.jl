export center_of_mass

"""
    center_of_mass(
        positions::Positions, indexes::AbstractVector{Int}; 
        masses::Union{Nothing,AbstractVector{<:Real}} = nothing, 
        unitcell::Union{Nothing,UnitCell} = nothing,
        iref::Int = indexes[1],
    )

Calculate the center of coordinates of a selection of atoms, or optionally
the center of mass if masses are provided. 

If a unitcell is provided, the positions are wrapped relative, by default, to the
first atom in the selection. This can be changed by setting `iref` to the index
of the atom to be used as the reference.

```julia-repl # to be doctest
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> using PDBTools

julia> atoms = readPDB(Testing.namd_pdb);

julia> protein_indexes = selindex(atoms, "protein");

julia> protein_masses = mass.(atoms);

julia> trajectory = Trajectory(Testing.namd_traj);

julia> positions = Positions(currentframe(trajectory));

julia> unitcell = UnitCell(currentframe(trajectory));

julia> cm = center_of_mass(positions, protein_indexes, masses=protein_masses, unitcell=unitcell)
3-element Point3D{Float64} with indices SOneTo(3):
 -3.7290442807974906
 -1.0539450244972206
 24.513632812542127

```

"""
function center_of_mass(
    positions::Positions, indexes::AbstractVector{Int}; 
    masses::Union{Nothing,AbstractVector{<:Real}} = nothing, 
    unitcell::Union{Nothing,UnitCell} = nothing,
    iref::Int = indexes[1],
)
    totmass = 0.0
    cm = MVector{3}(0.0, 0.0, 0.0)
    xref = positions[iref]
    for i in indexes
        m = isnothing(masses) ? 1.0 : masses[i]
        totmass += m
        if !isnothing(unitcell) 
            x = wrap(positions[i], xref, unitcell)
        else 
            x = positions[i]
        end
        cm += x * m
    end
    return Point3D(cm /= totmass)
end

@testitem "centerofmass" begin
    using MolSimToolkit.Testing
    using PDBTools


end
