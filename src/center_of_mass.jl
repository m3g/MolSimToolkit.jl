export center_of_mass

"""
    center_of_mass(
        indexes::AbstractVector{Int};
        simulation::Simulation,
        positions::FramePositions,
        iref::Int = max(1, div(length(indexes),2)),
    )

Calculate the center of mass of a selection of atoms in a simulation given the
positions. The selection is defined by the `indexes` vector, which is the indexes of the atoms.

The `iref` parameter is the index of the reference atom. The center of mass is calculated
by first computing the minimum-image of all atoms relative to this atom. By default,
it is the atom closest to the middle of the indexes vector.

```julia-repl # to be doctest
julia> import PDBTools 

julia> using MolSimToolkit, MolSimToolkit.Testing

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> protein_indexes = PDBTools.selindex(atoms(simulation), "protein");

julia> coor = positions(current_frame(simulation));

julia> cm = center_of_mass(protein_indexes, simulation, coor)
3-element Point3D{Float64} with indices SOneTo(3):
 -3.7290442807974906
 -1.5339226637687564
  1.960640754560446

```

"""
function center_of_mass(
    indexes::AbstractVector{Int},
    simulation::Simulation,
    p::FramePositions;
    iref::Int=max(1, div(length(indexes), 2))
)
    xref = p[iref]
    uc = unitcell(current_frame(simulation))
    totmass = 0.0
    cm = MVector{3}(0.0, 0.0, 0.0)
    for i in indexes
        m = atomic_mass(atoms(simulation)[i])
        x = wrap(p[i], xref, uc)
        cm += x * m
        totmass += m
    end
    return Point3D(cm /= totmass)
end

@testitem "centerofmass" begin
    using MolSimToolkit.Testing
    import PDBTools
    simulation = Simulation(Testing.namd_pdb, Testing.namd_traj)
    protein_indexes = PDBTools.selindex(atoms(simulation), "protein")
    coor = positions(current_frame(simulation))
    @test center_of_mass(protein_indexes, simulation, coor) ≈
          [-3.7290442807974906, -1.5339226637687564, 1.960640754560446]
end

center_of_mass(x::AbstractVector{<:AbstractVector}, ::Nothing) = sum(x) / length(x)

center_of_mass(x::AbstractVector{<:AbstractVector}, mass::AbstractVector) =
    sum(x[i] * mass[i] for i in eachindex(x, mass)) / sum(mass)

center_of_mass(atoms::AbstractVector{<:PDBTools.Atom}) =
    sum(PDBTools.coor(atom) * atomic_mass(atom) for atom in atoms) / sum(atomic_mass, atoms)
