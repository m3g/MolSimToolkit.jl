
"""
    Occupancy

Structure that wraps the result of the `occupancy` function.

# Fields

- `list::Vector{Vector{Int}}`: for each frame of the simulation, the list
  of the indices (relative to the solvent molecules, that is, `1` is the
  first solvent molecule, `2` the second, etc.) of the solvent molecules
  found within the cutoff distance of the binding site.
- `n_solvent_molecules::Int`: the total number of solvent molecules considered.

"""
struct Occupancy
    list::Vector{Vector{Int}}
    n_solvent_molecules::Int
end

"""
    mean(occupancy::Occupancy)

Returns the average number of solvent molecules found at the binding site
per frame.

# Example

```jldoctest
julia> using MolSimToolkit, PDBTools, MolSimToolkit.Testing

julia> sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj; frames=1:5);

julia> protein = select(get_atoms(sim), "protein");

julia> tmao = select(get_atoms(sim), "resname TMAO");

julia> occ = occupancy(sim, protein, tmao; solvent_natomspermol=14, cutoff=3.0, show_progress=false);

julia> mean(occ)
5.0

```

"""
mean(occ::Occupancy) = mean(length.(occ.list))

"""
    occupancy(
        sim::Simulation,
        site::AbstractVector{<:PDBTools.Atom},
        solvent::AbstractVector{<:PDBTools.Atom};
        solvent_natomspermol::Integer,
        cutoff::Real,
        show_progress::Bool = true,
    )

Computes, for each frame of the simulation, which solvent molecules are found
within `cutoff` of the binding `site`.

!!! note
    The distance considered is the **minimum distance** between the site atoms
    and the atoms of each solvent molecule.

# Positional Arguments

- `sim::Simulation`: Simulation object.
- `site::AbstractVector{<:PDBTools.Atom}`: Vector of atoms defining the binding site.
- `solvent::AbstractVector{<:PDBTools.Atom}`: Vector of solvent atoms.

# Keyword Arguments

- `solvent_natomspermol::Integer`: Number of atoms per solvent molecule.
- `cutoff::Real`: Cutoff distance.
- `show_progress::Bool`: Show progress bar. (optional, default: `true`)

# Returns

- `Occupancy`: an object wrapping, for each frame, the list of solvent
  molecules found within `cutoff` of the site.

# Example

```jldoctest
julia> using MolSimToolkit, PDBTools, MolSimToolkit.Testing

julia> sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj; frames=1:5);

julia> protein = select(get_atoms(sim), "protein");

julia> tmao = select(get_atoms(sim), "resname TMAO");

julia> occ = occupancy(sim, protein, tmao; solvent_natomspermol=14, cutoff=3.0, show_progress=false);

julia> length.(occ.list)
5-element Vector{Int64}:
 7
 3
 4
 5
 6

```

"""
function occupancy(
    sim::Simulation,
    site::AbstractVector{<:PDBTools.Atom},
    solvent::AbstractVector{<:PDBTools.Atom};
    solvent_natomspermol::Integer,
    cutoff::Real,
    show_progress::Bool=true,
)
    if length(solvent) % solvent_natomspermol != 0
        throw(ArgumentError("""\n
            Number of solvent atoms is not a multiple of solvent_natomspermol

        """))
    end
    inds_site = PDBTools.index.(site)
    inds_solvent = PDBTools.index.(solvent)
    n_solvent_molecules = length(inds_solvent) ÷ solvent_natomspermol

    first_frame!(sim)
    p = positions(current_frame(sim))
    uc = unitcell(current_frame(sim))
    sys = CrossPairs(
        xpositions=p[inds_solvent],
        ypositions=p[inds_site],
        xn_atoms_per_molecule=solvent_natomspermol,
        cutoff=cutoff,
        unitcell=uc.orthorhombic ? diag(uc.matrix) : uc.matrix,
    )

    list = Vector{Vector{Int}}(undef, length(sim))
    prg = Progress(length(sim); enabled=show_progress)
    for (iframe, frame) in enumerate(sim)
        p = positions(frame)
        uc = unitcell(frame)
        sys.system.xpositions .= @view(p[inds_solvent])
        sys.system.ypositions .= @view(p[inds_site])
        sys.system.unitcell = uc.orthorhombic ? diag(uc.matrix) : uc.matrix
        md_list = minimum_distances!(sys)
        list[iframe] = findall(md -> md.within_cutoff, md_list)
        next!(prg)
    end
    return Occupancy(list, n_solvent_molecules)
end

@testitem "occupancy" begin
    using MolSimToolkit, PDBTools, MolSimToolkit.Testing
    sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj)
    protein = select(get_atoms(sim), "protein")
    tmao = select(get_atoms(sim), "resname TMAO")
    occ = occupancy(
        sim, protein, tmao;
        solvent_natomspermol=14, cutoff=3.0,
        show_progress=false
    )
    @test occ isa Occupancy
    @test length(occ.list) == length(sim)
    @test occ.n_solvent_molecules == length(tmao) ÷ 14
    @test length.(occ.list) == [7, 3, 4, 5, 6, 5, 7, 7, 10, 5, 4, 4, 5, 6, 2, 9, 4, 5, 2, 6]
    @test all(imol -> 1 <= imol <= occ.n_solvent_molecules, reduce(vcat, occ.list))
    @test mean(occ) ≈ 5.3
end
