"""
    coordination_number(
        sim::Simulation,
        solute::AbstractVector{<:PDBTools.Atom},
        solvent::AbstractVector{<:PDBTools.Atom};
        solvent_natomspermol::Int,
        cutoff::Real,
        show_progress::Bool = true
    )

Calculate the coordination number of the solute atoms with the solvent atoms.

# Positional Arguments

- `sim::Simulation`: Simulation object.
- `solute::AbstractVector{<:PDBTools.Atom}`: Vector of solute atoms.
- `solvent::AbstractVector{<:PDBTools.Atom}`: Vector of solvent atoms.

# Keyword Arguments

- `solvent_natomspermol::Int`: Number of atoms per solvent molecule.
- `cutoff::Real`: Cutoff distance.
- `show_progress::Bool`: Show progress bar. (optional, default: `true`)

# Returns

- `cn::Vector{Int}`: Vector with the coordination number of the solute atoms with the solvent atoms, at each frame.

# Example

```jldoctest
julia> using MolSimToolkit, PDBTools, MolSimToolkit.Testing

julia> sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj; frames=1:5);

julia> protein = select(atoms(sim), "protein");

julia> tmao = select(atoms(sim), "resname TMAO");

julia> coordination_number(sim, protein, tmao; solvent_natomspermol=14, cutoff=3.0, show_progress=false)
5-element Vector{Int64}:
 7
 3
 4
 5
 6

```

"""
function coordination_number(
    sim::Simulation,
    solute::AbstractVector{<:PDBTools.Atom},
    solvent::AbstractVector{<:PDBTools.Atom};
    solvent_natomspermol::Int,
    cutoff::Real,
    show_progress::Bool = true,
)
    if length(solvent) % solvent_natomspermol != 0
        throw(ArgumentError("""\n
            Number of solvent atoms is not a multiple of solvent_natomspermol

        """))
    end
    inds_solute = PDBTools.index.(solute)
    inds_solvent = PDBTools.index.(solvent)
    cn = zeros(Int, length(sim))
    first_frame!(sim)
    p = positions(current_frame(sim))
    sys = CrossPairs(
        xpositions=p[inds_solvent],
        ypositions=p[inds_solute],
        xn_atoms_per_molecule=solvent_natomspermol,
        cutoff=cutoff,
        unitcell=unitcell(current_frame(sim)),
    )
    prg = Progress(length(sim); enabled=show_progress)
    for (iframe, frame) in enumerate(sim)
        p = positions(frame)
        sys.xpositions .= @view(p[inds_solvent])
        sys.ypositions .= @view(p[inds_solute])
        sys.unitcell = unitcell(frame)
        md_list = minimum_distances!(sys)
        cn[iframe] = count(md -> md.within_cutoff, md_list)
        next!(prg)
    end
    return cn
end

@testitem "coordination_number" begin
    using MolSimToolkit, PDBTools, MolSimToolkit.Testing
    sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj)
    protein = select(atoms(sim), "protein")
    tmao = select(atoms(sim), "resname TMAO")
    cn = coordination_number(
        sim, protein, tmao; 
        solvent_natomspermol=14, cutoff=3.0,
        show_progress=false
    )
    @test cn == [7, 3, 4, 5, 6, 5, 7, 7, 10, 5, 4, 4, 5, 6, 2, 9, 4, 5, 2, 6]
end

