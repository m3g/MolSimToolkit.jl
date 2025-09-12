function _sum_dihedrals!(
    vinds::AbstractVector{<:AbstractVector{<:Integer}},
    positions::FramePositions,
    init::AbstractVector{<:Real}; 
    degrees=true,
)
    for (i, v4) in enumerate(vinds)
        v1 = positions[v4[1]]
        v2 = positions[v4[2]]
        v3 = positions[v4[3]]
        v4 = positions[v4[4]]
        init[i] += dihedral(v1, v2, v3, v4; degrees)
    end
    return init
end

"""
    average_dihedrals(sim::Simulation, v::AbstractVector{<:AbstractVector{<:Integer}}; degrees=true)
    average_dihedrals(sim::Simulation, v::AbstractVector{<:AbstractVector{<:PDBTools.Atom}}; degrees=true)

Computes the average dihedral angles for many sets of 4 vectors from a trajectory. The input is a vector of vectors, 
containing the indices of the atoms forming the dihedral angles, or the `PDBTools.Atom` objects. 

The function returns a vector with the average dihedral angles in radians or degrees.

## Example

```jldoctest; filter = r"([0-9]+\\.[0-9]{2})[0-9]+" => s"\\1***"
julia> using MolSimToolkit, MolSimToolkit.Testing, PDBTools

julia> atoms = read_pdb(Testing.namd2_pdb);

julia> cAs = select(atoms, "name CA and residue < 5"); # 4 atoms

julia> r1b = select(atoms, "residue 1 and backbone"); # 4 atoms

julia> inds = [ index.(cAs), index.(r1b) ]; # List of vector of indices

julia> sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj);

julia> ds = average_dihedrals(sim, inds)
2-element Vector{Float64}:
 -60.12860673875001
  -0.3398274578758668

julia> ats = [ cAs, r1b ]; # List of vectors of PDBTools.Atom

julia> ds = average_dihedrals(sim, ats)
2-element Vector{Float64}:
 -60.12860673875001
  -0.3398274578758668

```

"""
function average_dihedrals(
    sim::Simulation,
    v::AbstractVector{<:AbstractVector{<:Integer}}; 
    degrees=true
)
    first_frame!(sim)
    p1 = positions(current_frame(sim))
    T = eltype(eltype(p1))
    init = zeros(T, length(v))
    first_frame!(sim)
    for frame in sim
        p = positions(frame)
        init = _sum_dihedrals!(v, p, init; degrees)
    end
    init .= init ./ length(sim)
    return init
end

function average_dihedrals(
    sim::Simulation,
    v::AbstractVector{<:AbstractVector{<:PDBTools.Atom}};
    degrees=true
)
    inds = [ PDBTools.index.(v) for v in v ]
    return average_dihedrals(sim, inds; degrees)
end

@testitem "average_dihedrals" begin
    using MolSimToolkit
    using MolSimToolkit.Testing
    using PDBTools

    atoms = read_pdb(Testing.namd2_pdb)
    sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj)
    first_frame!(sim)
    p1 = positions(current_frame(sim))
    cAs = select(atoms, "name CA and residue < 5")
    r1b = select(atoms, "residue 1 and backbone")
    inds = [ index.(cAs), index.(r1b) ]
    v4 = [coor(cAs[1]), coor(cAs[2]), coor(cAs[3]), coor(cAs[4])]
    init = zeros(eltype(eltype(eltype(v4))), length(inds))
    MolSimToolkit._sum_dihedrals!(inds, p1, init)
    @test init ≈ [-168.41887, -55.351723] atol = 1e-3
    MolSimToolkit._sum_dihedrals!(inds, p1, init)
    @test init ≈ 2 * [-168.41887, -55.351723] atol = 1e-3

    inds = [ index.(cAs), index.(r1b) ]
    sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj)
    ds = average_dihedrals(sim, inds)
    @test ds ≈ [-60.12860673875001, -0.3398274578758668] atol = 1e-5

    ats = [ cAs, r1b ]
    ds = average_dihedrals(sim, ats)
    @test ds ≈ [-60.12860673875001, -0.3398274578758668] atol = 1e-5
end
