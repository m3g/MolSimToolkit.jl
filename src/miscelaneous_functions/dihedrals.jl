"""
    dihedral(v1, v2, v3, v4; degrees=true)
    dihedral(v::AbstractVector; degrees=true)

Computes the dihedral angle between the planes formed by the vectors `v1-v2` and `v2-v3`, and `v2-v3` and `v3-v4`.
The input vectors must have 3 elements. The function returns the dihedral angle in radians or degrees.
If the input is a vector with 4 vectors, the function computes the dihedral angle between the planes 
formed by the vectors `v[1]-v[2]` and `v[2]-v[3]`, and `v[2]-v[3]` and `v[3]-v[4]`.

The optional argument `degrees` specifies whether the output is in degrees (default) or radians.

"""
function dihedral(v1::T, v2::T, v3::T, v4::T; degrees=true) where {T<:AbstractVector}
    if !(length(v1) == 3 && length(v2) == 3 && length(v3) == 3 && length(v4) == 3)
        throw(ArgumentError("The input vectors must have 3 elements"))
    end 

    # Vectors connecting consecutive atoms
    v21 = v2 - v1 
    v32 = v3 - v2 
    v43 = v4 - v3

    # Normal vectors to the planes formed by the vectors v1 and v2, and v2 and v3
    n1 = cross(v21, v32)
    n2 = cross(v32, v43)

    n1 = n1 / norm(n1)
    n2 = n2 / norm(n2)

    # Normalized vector perpendicular to the plane formed by the vectors n1 and n2
    u1 = cross(n1, n2) / norm(n1)

    # Normalize v2
    unitary_v32 = v32 / norm(v32)

    # Computes the projection of the vector u1 in the direction of unitary_v2,
    # and of n1 in the direction of n2
    m1 = dot(u1, unitary_v32)  
    m2 = dot(n1, n2)         

    # Computes the dihedral angle in radians or degrees
    return degrees ? atand(m1, m2) : atan(m1, m2)
end

function dihedral(v::AbstractVector; degrees=true)
    length(v) == 4 || throw(ArgumentError("The input vector must have 4 elements"))
    return dihedral(v[1], v[2], v[3], v[4]; degrees)
end

@testitem "dihedral" begin
    using MolSimToolkit, PDBTools
    using MolSimToolkit.Testing
    using StaticArrays

    atoms = read_pdb(Testing.namd2_pdb)
    d = dihedral(coor(atoms[1]), coor(atoms[2]), coor(atoms[3]), coor(atoms[4]))
    @test d ≈ -34.57 atol=1e-2

    cAs = select(atoms, "name CA and residue < 5")
    d = dihedral(coor(cAs[1]), coor(cAs[2]), coor(cAs[3]), coor(cAs[4]))
    @test d ≈ 164.43 atol=1e-2

    r1b = select(atoms, "residue 1 and backbone")
    d = dihedral(coor(r1b[1]), coor(r1b[2]), coor(r1b[3]), coor(r1b[4]))
    @test d ≈ -115.83 atol=1e-2
    d2 = dihedral(coor(r1b[1]), coor(r1b[2]), coor(r1b[3]), coor(r1b[4]); degrees=false)
    @test d2 ≈ deg2rad(d) 
    @test dihedral(SVector(coor(r1b[1]), coor(r1b[2]), coor(r1b[3]), coor(r1b[4]))) ≈ d

end

#
# Internal function that sums the dihedral angles for many sets of 4 vectors
# to compute the average dihedral angle from a trajectory
#
function _sum_dihedrals!(v::AbstractVector{<:AbstractVector}, init; degrees=true)
    for (i, v4) in enumerate(v)
        init[i] += dihedral(v4; degrees)
    end
    return init
end

"""
    dihedrals(v::AbstractVector{<:AbstractVector}; degrees=true)

Computes the dihedral angles for many sets of 4 vectors. The input is a vector of vectors, where each
element is a vector with 4 vectors. The function returns a vector with the dihedral angles in radians or degrees.

## Example

```jldoctest; filter = r"([0-9]+\\.[0-9]{2})[0-9]+" => s"\\1***"
julia> using MolSimToolkit, MolSimToolkit.Testing, PDBTools

julia> atoms = read_pdb(Testing.namd2_pdb);

julia> cAs = select(atoms, "name CA and residue < 5");

julia> r1b = select(atoms, "residue 1 and backbone");

julia> ds = dihedrals([ coor(cAs), coor(r1b) ])
2-element Vector{Float32}:
  164.4348
 -115.82544
```

"""
function dihedrals(v::AbstractVector{<:AbstractVector}; degrees=true)
    init = zeros(eltype(eltype(eltype(v))), length(v))
    return _sum_dihedrals!(v, init; degrees)
end

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

@testitem "dihedrals" begin
    using MolSimToolkit
    using MolSimToolkit.Testing
    using StaticArrays
    using PDBTools

    atoms = read_pdb(Testing.namd2_pdb)
    cAs = select(atoms, "name CA and residue < 5")
    d = dihedral(coor(cAs)...)
    v4 = [coor(cAs[1]), coor(cAs[2]), coor(cAs[3]), coor(cAs[4])]
    vv4 = [v4 for _ in 1:10]
    ds = dihedrals(vv4)
    @test all(ds .≈ d)
    init = zeros(eltype(eltype(eltype(v4))), length(vv4))
    @test MolSimToolkit._sum_dihedrals!(vv4, init) ≈ ds
    @test MolSimToolkit._sum_dihedrals!(vv4, init) ≈ 2 * ds

    sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj)
    first_frame!(sim)
    p1 = positions(current_frame(sim))
    cAs = select(atoms, "name CA and residue < 5")
    r1b = select(atoms, "residue 1 and backbone")
    inds = [ index.(cAs), index.(r1b) ]
    init = zeros(eltype(eltype(eltype(v4))), length(inds))
    MolSimToolkit._sum_dihedrals!(inds, p1, init)
    @test init ≈ [-168.41887, -55.351723] atol = 1e-3
    MolSimToolkit._sum_dihedrals!(inds, p1, init)
    @test init ≈ 2 * [-168.41887, -55.351723] atol = 1e-3

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
    cAs = select(atoms, "name CA and residue < 5")
    r1b = select(atoms, "residue 1 and backbone")
    inds = [ index.(cAs), index.(r1b) ]
    sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj)

    ds = average_dihedrals(sim, inds)
    @test ds ≈ [-60.12860673875001, -0.3398274578758668] atol = 1e-5

    ats = [ cAs, r1b ]
    ds = average_dihedrals(sim, ats)
    @test ds ≈ [-60.12860673875001, -0.3398274578758668] atol = 1e-5
end
