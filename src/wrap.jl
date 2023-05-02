export wrap 

"""
    wrap(x, xref, unit_cell_matrix::SMatrix{N,N,T}) where {N,T}
    wrap(x,xref,sides::AbstractVector)

Wraps the coordinates of point `x` such that it is the minimum image relative to `xref`. The unit cell 
may be given as a matrix of dimension `NxN` or as a vector of length `N` with the side lengths.

"""
@inline function wrap(x::AbstractVector, xref::AbstractVector, unit_cell_matrix::SMatrix{N,N,T}) where {N,T}
    invu = inv(oneunit(T))
    unit_cell_matrix = invu * unit_cell_matrix
    x_f = wrap_cell_fraction(invu*x, unit_cell_matrix)
    xref_f = wrap_cell_fraction(invu*xref, unit_cell_matrix)
    xw = wrap(x_f, xref_f, SVector{N,eltype(x_f)}(ntuple(i -> 1, N)))
    return oneunit(T) * unit_cell_matrix * (xw - xref_f) + xref
end

@inline function wrap(x::AbstractVector, xref::AbstractVector, sides::AbstractVector)
    xw = mod.(x - xref, sides)
    xw = _wrap_single_coordinate.(xw, sides)
    return xw + xref
end

"""
    fastmod1(x)

$(INTERNAL)

Computes `mod(x,1)`, quickly, using `x - floor(x)`. Maybe irrelevant.

"""
@inline fastmod1(x) = x - floor(x)

"""
    wrap_cell_fraction(x,unit_cell_matrix)

$(INTERNAL)

# Extended help

Obtaint the coordinates of `x` as a fraction of unit cell vectors, first
positive cell. `x` is a vector of dimension `N` and `cell` a matrix of 
dimension `NxN`

## Example

```julia-repl
julia> unit_cell_matrix = [ 10 0
                            0 10 ];

julia> x = [ 15, 13 ];

julia> wrap_cell_fraction(x,unit_cell_matrix)
2-element Vector{Float64}:
 0.5
 0.3
```

"""
@inline function wrap_cell_fraction(x::AbstractVector{T}, unit_cell_matrix::AbstractMatrix) where {T}
    # Division by `oneunit` is to support Unitful quantities. 
    # this workaround works here because the units cancel.
    # see: https://github.com/PainterQubits/Unitful.jl/issues/46
    x_stripped = x ./ oneunit(T)
    m_stripped = unit_cell_matrix ./ oneunit(T)
    p = fastmod1.(m_stripped \ x_stripped)
    # Boundary coordinates belong to the lower boundary
    p = ifelse.(p .== one(eltype(x_stripped)), zero(eltype(x_stripped)), p)
    return p
end

"""
    wrap_to_first(x,unit_cell_matrix)

$(INTERNAL)

# Extended help

Wraps the coordinates of point `x` such that the returning coordinates are in the
first unit cell with all-positive coordinates. 

## Example

```julia-repl
julia> unit_cell_matrix = [ 10 0
                            0 10 ];

julia> x = [ 15, 13 ];

julia> wrap_to_first(x,unit_cell_matrix)
2-element Vector{Float64}:
 5.0
 3.0000000000000004
```

"""
@inline function wrap_to_first(x, unit_cell_matrix)
    p = wrap_cell_fraction(x, unit_cell_matrix)
    return unit_cell_matrix * p
end

#
# Wrap a single coordinate
#
@inline function _wrap_single_coordinate(x, s)
    if x >= s / 2
        x = x - s
    elseif x < -s / 2
        x = x + s
    end
    return x
end


@testitem "wrap" begin
    # Test wrapping of a vector
    x = SVector{3}([0.5, 0.5, 0.5])
    x_ref = SVector{3}([0.0, 0.0, 0.0])
    unitcell = Chemfiles.UnitCell([ 10.0 0.0 0.0 
                                    0.0 10.0 0.0
                                    0.0 0.0 10.0 ])    
    @test wrap(x, x_ref, unitcell) ≈ x
    # Test if vector is not inside unit cell
    x = SVector{3}([10.5, 10.5, 10.5])
    x_ref = SVector{3}([0.0, 0.0, 0.0])
    @test wrap(x, x_ref, unitcell) ≈ SVector{3}([0.5, 0.5, 0.5])
    # Test with a triclinic cell
    unitcell = Chemfiles.UnitCell([ 10.0 5.0 5.0 
                                    0.0 10.0 5.0
                                    0.0 0.0 10.0])
    x = SVector{3}([0.5, 0.5, 0.5])
    x_ref = SVector{3}([0.0, 0.0, 0.0])
    @test wrap(x, x_ref, unitcell) ≈ SVector{3}([0.5, 0.5, 0.5])
    # Test with atoms far from each other
    x = SVector{3}([10.5, 10.5, 10.5])
    x_ref = SVector{3}([0.0, 0.0, 0.0])
    @test wrap(x, x_ref, unitcell) ≈ SVector{3}([0.5, 0.5, 0.5])
end
