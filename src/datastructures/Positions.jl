
export Point3D, Positions

"""
    Point3D{T}

A point in 3D space with coordinates `x`, `y`, and `z` of type `T`.

"""
struct Point3D{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

"""
    Positions{T<:AbstractArray}

Container for the positions of a set of atoms. The positions are stored in a matrix,
where each column corresponds to the coordinates of an atom. The container is used
such that using the coodinates from a `Chemfiles.Frame` is transparent to the user.

# Example

```julia-repl
julia> using MolSimToolkit

julia> import Chemfiles

julia> traj = Chemfiles.Trajectory(Testing.namd_traj);

julia> frame = Chemfiles.read(traj);

julia> positions = Positions(frame);

julia> positions[1]
3-element Point3D{Float64} with indices SOneTo(3):
  5.912472724914551
 10.768872261047363
 28.277008056640625

 julia> positions[1].x
 5.912472724914551
 
 julia> positions[1].y
 10.768872261047363
 
 julia> positions[1].z
 28.277008056640625

```

"""
struct Positions{T<:AbstractArray}
    positions::T
end
Positions(f::Chemfiles.Frame) = Positions(Chemfiles.positions(f))
Base.getindex(x::Positions, i::Int) = Point3D(@view(x.positions[:,i]))

@testitem "Positions" begin
    m = rand(3,10)
    p = Positions(m)
    @test p[1] == Point3D(m[1,1], m[2,1], m[3,1])
    @test p[1].x == m[1,1]
    @test p[1].y == m[2,1]
    @test p[1].z == m[3,1]
end
