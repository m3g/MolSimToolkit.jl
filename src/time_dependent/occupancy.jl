
"""
    Occupancy

Structure that wraps the result of the `occupancy` function.

# Fields

- `list::Vector{Vector{Int}}`: for each frame of the simulation, the list
  of the indices (relative to the solvent molecules, that is, `1` is the
  first solvent molecule, `2` the second, etc.) of the solvent molecules
  found within the cutoff distance of the binding site.
- `distances::Vector{Vector{T}}`: for each frame of the simulation, the
  minimum distance between the site and each of the solvent molecules
  listed in `list`. `distances[iframe][i]` is the distance associated to
  the solvent molecule `list[iframe][i]`.
- `n_solvent_molecules::Int`: the total number of solvent molecules considered.
- `dmax::T`: the cutoff distance used in the calculation, that is, the
  maximum distance for which a solvent molecule is considered to be at
  the site. All values in `distances` are smaller than `dmax`.

!!! compat
    The `distances` and `dmax` fields were added in version 2.4.0 of MolSimToolkit.

"""
struct Occupancy{T<:Real}
    list::Vector{Vector{Int}}
    distances::Vector{Vector{T}}
    n_solvent_molecules::Int
    dmax::T
end

#=
    Occupancy(list, n_solvent_molecules)

Convenience constructor for occupancy data for which the site-solvent
distances are not available. The `distances` are filled with `NaN` and
`dmax` is set to `NaN`, such that the functions that require the distances
(as `intermittent_correlation_profile`) throw an informative error.
=#
Occupancy(list::AbstractVector{<:AbstractVector{<:Integer}}, n_solvent_molecules::Integer) =
    Occupancy(list, [fill(NaN, length(l)) for l in list], n_solvent_molecules, NaN)

function Base.show(io::IO, ::MIME"text/plain", occ::Occupancy)
    print(io, chomp(
        """
        -------------------------------------------------------------------
        Occupancy data:
        -------------------------------------------------------------------
        Total number of solvent molecules: $(occ.n_solvent_molecules)
        Maximum distance to the site (dmax): $(occ.dmax)
        Mean occupancy: $(mean(occ))
        Minimum occupancy: $(minimum(length(l) for l in occ.list))
        Maximum occupancy: $(maximum(length(l) for l in occ.list))
        -------------------------------------------------------------------
        """
    ))
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
  molecules found within `cutoff` of the site, and the corresponding
  site-solvent minimum distances. The `cutoff` is stored in the `dmax`
  field of the returned object.

# Example

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2***"
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

julia> occ.dmax
3.0

julia> occ.distances[2] # distances of the molecules of occ.list[2]
3-element Vector{Float64}:
 1.9506459897129325
 2.8333266769863195
 2.899149006265755

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
    dists = Vector{Vector{Float64}}(undef, length(sim))
    prg = Progress(length(sim); enabled=show_progress)
    for (iframe, frame) in enumerate(sim)
        p = positions(frame)
        uc = unitcell(frame)
        sys.system.xpositions .= @view(p[inds_solvent])
        sys.system.ypositions .= @view(p[inds_site])
        sys.system.unitcell = uc.orthorhombic ? diag(uc.matrix) : uc.matrix
        md_list = minimum_distances!(sys)
        imols = findall(md -> md.within_cutoff, md_list)
        list[iframe] = imols
        dists[iframe] = [Float64(md_list[imol].d) for imol in imols]
        next!(prg)
    end
    return Occupancy(list, dists, n_solvent_molecules, Float64(cutoff))
end

@testitem "occupancy" begin
    using MolSimToolkit, PDBTools, MolSimToolkit.Testing
    using ShowMethodTesting
    using LinearAlgebra: diag

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

    # distances are stored alongside the molecule indices, and are all within dmax
    @test occ.dmax == 3.0
    @test length.(occ.distances) == length.(occ.list)
    @test all(d -> 0 < d < occ.dmax, reduce(vcat, occ.distances))
    # the stored distance must be the minimum distance between the site and the molecule
    d, imol = findmin(occ.distances[1])
    imol = occ.list[1][imol]
    md = minimum_distances(
        xpositions=positions(first_frame!(sim))[index.(tmao)],
        ypositions=positions(current_frame(sim))[index.(protein)],
        xn_atoms_per_molecule=14,
        cutoff=3.0,
        unitcell=diag(unitcell(current_frame(sim)).matrix),
    )
    @test md[imol].d ≈ d
    # occupancy data without distance information
    occ_nod = Occupancy(occ.list, occ.n_solvent_molecules)
    @test isnan(occ_nod.dmax)
    @test all(isnan, reduce(vcat, occ_nod.distances))

    @test parse_show(occ) ≈ """ 
    -------------------------------------------------------------------
    Occupancy data:
    -------------------------------------------------------------------
    Total number of solvent molecules: 181
    Maximum distance to the site (dmax): 3.0
    Mean occupancy: 5.3
    Minimum occupancy: 2
    Maximum occupancy: 10
    -------------------------------------------------------------------
    """

    @test_throws "not a multiple" occupancy(
        sim, protein, tmao;
        solvent_natomspermol=13, cutoff=3.0,
        show_progress=false
    )

end
