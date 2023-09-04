export Point3D
export Positions
export positions

"""
    Point3D{T}

A point in 3D space with coordinates `x`, `y`, and `z` of type `T`.

"""
struct Point3D{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end
Point3D(x::Vector{T}) where {T} = Point3D{T}(x[1], x[2], x[3])

"""
    Positions{T<:AbstractArray}

Container for the positions of a set of atoms. The positions are stored in a matrix,
where each column corresponds to the coordinates of an atom. The container is used
such that using the coodinates from a `Chemfiles.Frame` is transparent to the user,
and the coordinates can be accessed as `p[i]` where `i` is the index of the
atom. 

The coordinates of the atom can be accessed as `p[i].x`, `p[i].y`, and `p[i].z`.

A `Positions` object can be created with the `positions` function.

# Example

```julia-repl # to be doctest
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> frame = current_frame(simulation);

julia> p = positions(frame);

julia> p[1]
3-element Point3D{Float64} with indices SOneTo(3):
  5.912472724914551
 10.768872261047363
 28.277008056640625

julia> p[1].x
5.912472724914551

julia> p[1].y
10.768872261047363

julia> p[1].z
28.277008056640625

```

"""
struct Positions{T<:AbstractArray}
    positions::T
end
positions(f::Chemfiles.Frame) = Positions(Chemfiles.positions(f))
Positions(f::Chemfiles.Frame) = Positions(Chemfiles.positions(f))
Base.getindex(x::Positions, i::Int) = Point3D(@view(x.positions[:,i]))

import Base: show
function show(io::IO, positions::Positions)
    print(io, "Positions{", eltype(positions.positions), "} with ", size(positions.positions,2), " atoms")
end

@testitem "Positions" begin
    m = rand(3,10)
    p = Positions(m)
    @test p[1] == Point3D(m[1,1], m[2,1], m[3,1])
    @test p[1].x == m[1,1]
    @test p[1].y == m[2,1]
    @test p[1].z == m[3,1]
end
