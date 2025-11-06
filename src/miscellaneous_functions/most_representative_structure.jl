export most_representative_structure

#
# Computes the average structure given the coordinates in a vector of frames.
#
function _average_structure(simulation, reference, indices)
    avg_struct = zero(reference)
    for frame in simulation
        p = positions(frame)[indices]
        align!(p, reference)
        avg_struct .+= p
    end
    avg_struct ./= length(simulation)
end

"""
    most_representative_structure(simulation::Simulation; atoms = nothing)

Find the most representative structure in a simulation. 
The most representative structure is the one that minimizes the RMSD with respect to the average structure of the simulation. 
The average structure is defined iteratively, first by aligning all frames to the first frame, and then by averaging the aligned structures.
The structure most similar to the average is then identified and used as the reference structure for the next iteration.
The process is repeated until the structure most similar to the average is the same as the previous iteration.

# Arguments

- `simulation::Simulation`: Simulation object.
- `atoms`: Atoms to consider in the calculation:
   - `atoms` is `nothing`: the function will consider all alpha-carbons in proteins (`"protein and name CA"`).
   - `atoms` is an `AbstractVector{<:PDBTools.Atom}`: the function will consider the atoms in the vector.
   - `atoms` is an `AbstractVector{<:Int}`: the function will consider the atoms with the indices in the vector.
   - `atoms` is a `String`: the function will consider the atoms selected by the string.

# Returns

- Tuple `(Int, Float64)`, with:
  - Index of the most representative structure.
  - RMSD of the most representative structure with respect to the average structure.

# Example

```julia-repl
julia> using MolSimToolkit, MolSimToolkit.Testing, PDBTools

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> most_representative_structure(simulation) # atoms == nothing (all alpha-carbons in proteins)
(4, 1.1681526249035976)

julia> most_representative_structure(simulation; atoms = "protein and name CA") # atoms is a String
(4, 1.1681526249035976)

julia> calphas = select(atoms(simulation), "name CA");

julia> most_representative_structure(simulation; atoms = calphas) # atoms is an Vector{PDBTools.Atom}
(4, 1.1681526249035976)

julia> ica = PDBTools.index.(calphas)

julia> most_representative_structure(simulation; atoms = ica) # atoms is an vector of indices
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
        elseif atoms isa AbstractVector{<:Integer}
            indices = atoms
        elseif atoms isa AbstractString
            indices = PDBTools.index.(
                PDBTools.select(MolSimToolkit.atoms(simulation), atoms)
            )
        else
            throw(ArgumentError("atoms must be a selection string, an array of PDBTools.Atom's or vector of integers with atom indices."))
        end
    end
    if length(indices) == 0
        throw(ArgumentError("No atoms selected."))
    end
    self_consistent = false
    first_frame!(simulation)
    reference = positions(current_frame(simulation))[indices]
    avg_structure = _average_structure(simulation, reference, indices) 
    rmsdmin, imin = findmin(
        frame -> rmsd(align!(positions(frame)[indices], avg_structure), avg_structure), 
        frame for frame in simulation
    )

    if imin == 1
        self_consistent = true 
    end
    while !self_consistent
        reference .= avg_structure
        avg_structure = _average_structure(simulation, reference, indices) 
        new_rmsdmin, new_imin = findmin(
            frame -> rmsd(align!(positions(frame)[indices], avg_structure), avg_structure), 
            frame for frame in simulation
        )
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
    imin, rmsdmin = most_representative_structure(simulation)
    @test imin == 4
    @test rmsdmin ≈ 1.1681514813293536
    imin, rmsdmin = most_representative_structure(simulation; atoms = calphas)
    @test imin == 4
    @test rmsdmin ≈ 1.1681514813293536
    imin, rmsdmin = most_representative_structure(simulation; atoms = PDBTools.index.(calphas))
    @test imin == 4
    @test rmsdmin ≈ 1.1681514813293536
    imin, rmsdmin = most_representative_structure(simulation; atoms = "protein and name CA")
    @test imin == 4
    @test rmsdmin ≈ 1.1681514813293536
    @test_throws ArgumentError most_representative_structure(simulation; atoms = 1)
    @test_throws ArgumentError most_representative_structure(simulation; atoms = Int[])
end