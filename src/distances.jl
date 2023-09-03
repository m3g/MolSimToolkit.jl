"""
    distances(trajectory, indexes1, indexes2; masses=nothing)

Function that calculates the distance between the center of masses of the selections in a trajectory,
or the center of coordinates if masses are not provided.

```julia-repl # to be doctest
julia> using PDBTools

julia> using MolSimToolkit, MolSimToolkit.Testing

julia> trajectory = Trajectory(Testing.namd_traj);

julia> atoms = readPDB(Testing.namd_pdb);

julia> masses = mass.(atoms);

julia> i1 = selindex(atoms, "protein and residue 1");

julia> i2 = selindex(atoms, "protein and residue 15");

julia> distances(trajectory, i1, i2; masses = masses)
5-element Vector{Float64}:
 23.433267858947584
 30.13791365033211
 28.48617683945202
 27.92740141686934
 23.235012287435566

```

"""
function distances(
    trajectory::Trajectory, 
    indexes1::AbstractVector{Int}, 
    indexes2::AbstractVector{Int};
    masses = nothing,
)
    distances = zeros(length(trajectory))
    for (iframe, frame) in enumerate(trajectory)
        unitcell = UnitCell(frame)
        coordinates = Positions(frame)
        cm1 = center_of_mass(coordinates, indexes1; masses = masses, unitcell = unitcell)
        cm2 = center_of_mass(coordinates, indexes2; masses = masses, unitcell = unitcell)
        cm2_wrapped = wrap(cm2,cm1,unitcell)
        d = norm(cm2_wrapped - cm1)
        distances[iframe] = d
    end
    return distances
end

@testitem "distances" begin
    using PDBTools
    using MolSimToolkit.Testing
    trajectory = Trajectory(Testing.namd_traj)
    atoms = readPDB(Testing.namd_pdb)
    masses = mass.(atoms)
    i1 = selindex(atoms, "protein and residue 1")
    i2 = selindex(atoms, "protein and residue 15")
    @test distances(trajectory, i1, i2; masses = masses) ≈ 
        [23.433267858947584, 30.13791365033211, 28.48617683945202, 27.92740141686934, 23.235012287435566]

end

