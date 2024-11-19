export Point3D
export FramePositions
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
    FramePositions{T,P<:Point3D{T},M<:AbstractArray{T}} <: AbstractVector{P}

Container for the positions of a set of atoms. The positions are stored in a matrix,
where each column corresponds to the coordinates of an atom. The container is used
such that using the coodinates from a `Chemfiles.Frame` is transparent to the user,
and the coordinates can be accessed as `p[i]` where `i` is the index of the
atom. 

The coordinates of the atom can be accessed as `p[i].x`, `p[i].y`, and `p[i].z`.

A `FramePositions` object can be created with the `positions` function. The
construction with `FramePositions` is not considered part of the public API.

# Example

```julia-repl # to be doctest
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> frame = current_frame(simulation);

julia> p = positions(frame)
20465-element FramePositions{Float64, Point3D{Float64}, Chemfiles.ChemfilesArray}:
 [5.912472724914551, 10.768872261047363, 28.277008056640625]
 [5.040304183959961, 10.810898780822754, 27.71207046508789]
 ⋮
 [11.16289234161377, -37.30374526977539, 22.80788230895996]

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

It is also possible to take slices and views of the positions:

```
julia> p[1:2]
2-element FramePositions{Float64, Point3D{Float64}, Matrix{Float64}}:
 [5.912472724914551, 10.768872261047363, 28.277008056640625]
 [5.040304183959961, 10.810898780822754, 27.71207046508789]

julia> @view(coor[1:2])
 2-element FramePositions{Float64, Point3D{Float64}, SubArray{Float64, 2, Chemfiles.ChemfilesArray, Tuple{Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}}, false}}:
  [5.912472724914551, 10.768872261047363, 28.277008056640625]
  [5.040304183959961, 10.810898780822754, 27.71207046508789]
 
```

"""
struct FramePositions{T,P<:Point3D{T},M<:AbstractArray{T}} <: AbstractVector{P}
    positions::M
end
FramePositions(m::AbstractMatrix{T}) where {T} = FramePositions{T,Point3D{T},typeof(m)}(m)
FramePositions(f::Chemfiles.Frame) = positions(Chemfiles.positions(f))
FramePositions(x::FramePositions{T,P}) where {T,P} = FramePositions{T,P,typeof(x.positions)}(x.positions)

Base.getindex(x::FramePositions, i::Integer) = Point3D(@view(x.positions[:, i]))
Base.getindex(x::FramePositions, r::AbstractUnitRange) = FramePositions(x.positions[:, r])
Base.getindex(x::FramePositions, ivec::AbstractVector{<:Integer}) = FramePositions(x.positions[:, ivec])

Base.setindex!(x::FramePositions, value::AbstractVector, i::Integer) = x.positions[:, i] .= value

Base.length(x::FramePositions) = size(x.positions, 2)
Base.size(x::FramePositions) = (length(x),)


"""
    positions(frame::Chemfiles.Frame)

Return the positions of the atoms in a `Chemfiles.Frame` as a `FramePositions` object.

This is the default way to access the positions of the atoms in a simulation.

# Example

```julia-repl # to be doctest
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> frame = current_frame(simulation);

julia> p = positions(frame);

julia> p[1].x 
5.912472724914551

```

"""
function positions(f::Chemfiles.Frame)
    p = Chemfiles.positions(f)
    return FramePositions{eltype(p),Point3D{eltype(p)},typeof(p)}(p)
end

import Base: copy
copy(positions::FramePositions) = FramePositions(copy(positions.positions))

import Base: ==, ≈
==(x::FramePositions, y::FramePositions) = ==(x.positions, y.positions)
≈(x::FramePositions, y::FramePositions; kargs...) = ≈(x.positions, y.positions; kargs...)

import Base: view
view(positions::FramePositions, r::AbstractUnitRange) = FramePositions(@view(positions.positions[:, r]))
view(positions::FramePositions, ivec::AbstractVector{<:Integer}) = FramePositions(@view(positions.positions[:, ivec]))


@testitem "FramePositions" begin
    using BenchmarkTools

    m = rand(3, 10)
    p = FramePositions(m)
    @test p[1] == Point3D(m[1, 1], m[2, 1], m[3, 1])
    @test p[1].x == m[1, 1]
    @test p[1].y == m[2, 1]
    @test p[1].z == m[3, 1]

    # test with range
    @test p[2:3] == FramePositions(m[:, 2:3])
    @test p[2:3] ≈ FramePositions(m[:, 2:3])
    @test @view(p[2:3]) == FramePositions(m[:, 2:3])
    inds = [2, 3, 5]
    @test p[inds] == FramePositions(m[:, inds])
    @test @view(p[inds]) == FramePositions(m[:, inds])
    @test iszero(@ballocated @view($p[$inds]))

end
