"""
    distances(simulation, indexes1, indexes2)

Function that calculates the distance between the centers of mass of two selections in a simulation.

The selections are defined by the indexes1 and indexes2 vectors, which are the indexes of the atoms.

# Example

```julia-repl # to be doctest
julia> import PDBTools

julia> using MolSimToolkit, MolSimToolkit.Testing

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> i1 = PDBTools.selindex(atoms(simulation), "protein and residue 1");

julia> i2 = PDBTools.selindex(atoms(simulation), "protein and residue 15");

julia> distances(simulation, i1, i2)
5-element Vector{Float64}:
 23.433267858947584
 30.13791365033211
 28.48617683945202
 27.92740141686934
 23.235012287435566

```

"""
function distances(
    simulation::Simulation,
    indexes1::AbstractVector{Int},
    indexes2::AbstractVector{Int};
)
    distances = zeros(length(simulation))
    for (iframe, frame) in enumerate(simulation)
        coor = positions(frame)
        cm1 = center_of_mass(indexes1, simulation, coor)
        cm2 = center_of_mass(indexes2, simulation, coor)
        cm2_wrapped = wrap(cm2, cm1, unitcell(frame))
        d = norm(cm2_wrapped - cm1)
        distances[iframe] = d
    end
    return distances
end

@testitem "distances" begin
    import PDBTools
    using MolSimToolkit.Testing
    simulation = Simulation(Testing.namd_pdb, Testing.namd_traj)
    atoms = PDBTools.readPDB(Testing.namd_pdb)
    i1 = PDBTools.selindex(atoms, "protein and residue 1")
    i2 = PDBTools.selindex(atoms, "protein and residue 15")
    @test distances(simulation, i1, i2) ≈
          [23.433267858947584, 30.13791365033211, 28.48617683945202, 27.92740141686934, 23.235012287435566]

end

