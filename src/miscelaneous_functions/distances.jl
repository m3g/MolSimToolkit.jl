"""
    distances(simulation, indices1::AbstractVector{Int}, indices2::AbstractVector{Int})
    distances(simulation, selection1::AbstractVector{PDBTools.Atom}, selection2::AbstractVector{PDBTools.Atom})

Function that calculates the distance between the centers of mass of two selections in a simulation.

The selections are defined by the indices1 and indices2 vectors, which are the indices of the atoms,
or by the selection1 and selection2 vectors, which are vectors of `PDBTools.Atom` objects.

Use `silent=true` to suppress the progress bar.

# Example

```jldoctest; filter = r"([0-9]+\\.[0-9]{2})[0-9]+" => s"\\1***"
julia> using PDBTools

julia> using MolSimToolkit, MolSimToolkit.Testing

julia> sim = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> ats = atoms(sim);

julia> i1 = findall(sel"protein and residue 1", ats); # indices

julia> i2 = findall(sel"protein and residue 15", ats); # indices

julia> distances(sim, i1, i2; silent=true)
5-element Vector{Float64}:
 23.433267858947584
 30.13791365033211
 28.48617683945202
 27.92740141686934
 23.235012287435566

julia> distances(sim, 
           filter(sel"protein and residue 1", ats), # selection (PDBTools.Atom)
           filter(sel"protein and residue 15", ats); # selection (PDBTools.Atom) 
           silent=true
       )
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
    silent::Bool = false,
)
    distances = zeros(length(simulation))
    silent || (p = Progress(length(simulation)))
    for (iframe, frame) in enumerate(simulation)
        coor = positions(frame)
        cm1 = center_of_mass(indices1, simulation, coor)
        cm2 = center_of_mass(indices2, simulation, coor)
        cm2_wrapped = wrap(cm2, cm1, unitcell(frame))
        d = norm(cm2_wrapped - cm1)
        distances[iframe] = d
        silent || next!(p)
    end
    return distances
end

# Using as input vectors of `PDBTools.Atom` objects
function distances(
    simulation::Simulation,
    selection1::AbstractVector{<:PDBTools.Atom},
    selection2::AbstractVector{<:PDBTools.Atom};
    silent::Bool = false,
)
    indices1 = PDBTools.index.(selection1)
    indices2 = PDBTools.index.(selection2)
    return distances(simulation, indices1, indices2; silent=silent)
end

@testitem "distances" begin
    using PDBTools
    using MolSimToolkit.Testing
    simulation = Simulation(Testing.namd_pdb, Testing.namd_traj)
    i1 = findall(sel"protein and residue 1", atoms(simulation))
    i2 = findall(sel"protein and residue 15", atoms(simulation))
    @test distances(simulation, i1, i2; silent=true) ≈
          [23.433267858947584, 30.13791365033211, 28.48617683945202, 27.92740141686934, 23.235012287435566]
    s1 = filter(sel"protein and residue 1", atoms(simulation))
    s2 = filter(sel"protein and residue 15", atoms(simulation))
    @test distances(simulation, s1, s2; silent=true) ≈
          [23.433267858947584, 30.13791365033211, 28.48617683945202, 27.92740141686934, 23.235012287435566]
end
