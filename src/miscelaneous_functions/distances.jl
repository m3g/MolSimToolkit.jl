"""
    distances(simulation, indices1, indices2)

Function that calculates the distance between the centers of mass of two selections in a simulation.

The selections are defined by the indices1 and indices2 vectors, which are the indices of the atoms.

# Example

```jldoctest; filter = r"([0-9]+\\.[0-9]{2})[0-9]+" => s"\\1***"
julia> using PDBTools

julia> using MolSimToolkit, MolSimToolkit.Testing

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> i1 = findall(sel"protein and residue 1", atoms(simulation));

julia> i2 = findall(sel"protein and residue 15", atoms(simulation));

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
    indices1::AbstractVector{Int},
    indices2::AbstractVector{Int};
)
    distances = zeros(length(simulation))
    for (iframe, frame) in enumerate(simulation)
        coor = positions(frame)
        cm1 = center_of_mass(indices1, simulation, coor)
        cm2 = center_of_mass(indices2, simulation, coor)
        cm2_wrapped = wrap(cm2, cm1, unitcell(frame))
        d = norm(cm2_wrapped - cm1)
        distances[iframe] = d
    end
    return distances
end

@testitem "distances" begin
    using PDBTools
    using MolSimToolkit.Testing
    simulation = Simulation(Testing.namd_pdb, Testing.namd_traj)
    i1 = findall(sel"protein and residue 1", atoms(simulation))
    i2 = findall(sel"protein and residue 15", atoms(simulation))
    @test distances(simulation, i1, i2) ≈
          [23.433267858947584, 30.13791365033211, 28.48617683945202, 27.92740141686934, 23.235012287435566]
end

