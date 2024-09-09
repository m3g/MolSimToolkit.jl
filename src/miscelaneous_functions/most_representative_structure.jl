export most_representative_structure

#
# Computes the average structure given the coordinates in a vector of frames.
#
function _average_structure!(simulation, reference, indices)
    avg_struct = zero(reference)
    for frame in simulation
        p = positions(frame)[indices]
        align!(p, reference)
        avg_struct .+= p
    end
    avg_struct /= length(simulation)
end

"""
    most_representative_structure(simulation::Simulation; atoms = nothing)

Find the most representative structure in a simulation. 
The most representative structure is the one that minimizes the RMSD with respect to the average structure of the simulation. 
The average structure is calculated by aligning all frames to a reference structure and averaging the aligned structures. 
The reference structure is updated until the most representative structure found is consistent with the previous one.

# Arguments

- `simulation::Simulation`: Simulation object.
- `atoms::Union{Nothing, AbstractVector{<:PDBTools.Atom}, AbstractVector{Int}}`: Atoms to consider in the calculation.

If `atoms` is `nothing`, the function will consider all alpha carbons in proteins ("protein and name CA").

# Returns

- `imin::Int`: Index of the most representative structure.
- `rmsdmin::Float64`: RMSD of the most representative structure with respect to the average structure.

# Example

```julia-repl
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> calphas = select(atoms(simulation), "name CA");

julia> most_representative_structure(simulation; atoms = calphas)
(4, 1.1681526249035976)
```

"""
function most_representative_structure(simulation::Simulation; atoms = nothing)
    if isnothing(atoms)
        indices = PDBTools.index.(
            PDBTools.select(MolSimToolkit.atoms(simulation), "protein and name CA")
        )
    else
        if atoms isa AbstractVector{<:PDBTools.Atom}
            indices = PDBTools.index.(atoms)
        elseif atoms isa AbstractVector{Int}
            indices = atoms
        else
            throw(ArgumentError("atoms must be an AbstractVector{<:PDBTools.Atom} or an AbstractVector{<:Int}, with atom indices."))
        end
    end
    if length(indices) == 0
        throw(ArgumentError("No atoms selected."))
    end
    self_consistent = false
    firstframe!(simulation)
    reference = positions(current_frame(simulation))[indices]
    avg_structure = _average_structure!(simulation, reference, indices) 
    rmsdmin, imin = findmin(
        frame -> rmsd(positions(frame)[indices], avg_structure), 
        frame for frame in simulation
    )
    if imin == 1
        self_consistent = true 
    end
    while !self_consistent
        reference = avg_structure
        avg_structure = _average_structure!(simulation, reference, indices) 
        new_rmsdmin, new_imin = findmin(
            frame -> rmsd(positions(frame)[indices], avg_structure), 
            frame for frame in simulation
        )
        @show imin, new_imin, rmsdmin, new_rmsdmin
        if new_imin == imin
            self_consistent = true
        end
        imin = new_imin
        rmsdmin = new_rmsdmin
    end 
    return imin, rmsdmin
end

@testitem "most_representative_structure" begin
    using MolSimToolkit, MolSimToolkit.Testing, PDBTools
    simulation = Simulation(Testing.namd_pdb, Testing.namd_traj)
    calphas = select(atoms(simulation), "name CA")
    imin, rmsdmin = most_representative_structure(simulation; atoms = calphas)
    @test imin == 4
    @test isapprox(rmsdmin, 1.1681526249035976)
    imin, rmsdmin = most_representative_structure(simulation; atoms = PDBTools.index.(calphas))
    @test imin == 4
    @test isapprox(rmsdmin, 1.1681526249035976)
    @test_throws ArgumentError most_representative_structure(simulation; atoms = 1)
    @test_throws ArgumentError most_representative_structure(simulation; atoms = Int[])
end