module MolecularMinimumDistances

import TestItems: @testitem
import DocStringExtensions: TYPEDEF, TYPEDFIELDS

import PDBTools: distance
import StaticArrays: SVector
import CellListMap: _uround # interal function, used for showing the unitcell
import CellListMap: ParticleSystem, map_pairwise, map_pairwise!,
    update_cutoff!, update_unitcell!

export MinimumDistance
export SelfPairs, CrossPairs, AllPairs
export minimum_distances, minimum_distances!
export within_cutoff, distance, iatom, jatom

"""
    MinimumDistance{T}

The lists of minimum-distances are stored in arrays of type `Vector{MinimumDistance{T}}`. The index
of this vector corresponds to the index of the molecule in the original array.

`MinimumDistance{T}` is a simple structure that contains four fields: a boolean marker indicating
if the distance is within the cutoff, the indices `i` and `j` of the atoms of the 
molecules that are closer to each other, and the distance `d`, with type `T`, which is
the same as that of the coordinates of the input vectors of coordinates. The best way
to access the information of a `MinimumDistance` element is through the getter functions
`within_cutoff`, `distance`, `iatom`, and `jatom`.

## Example

```julia-repl
julia> md = MinimumDistance{Float32}(true, 2, 5, 1.f0)
MinimumDistance{Float32}(true, 2, 5, 1.0f0)

julia> iatom(md)
2

julia> jatom(md)
5

julia> distance(md)
1.0f0

julia> within_cutoff(md)
true
```

"""
struct MinimumDistance{T}
    within_cutoff::Bool
    i::Int
    j::Int
    d::T
end
import Base: zero, copy, convert
zero(::Type{MinimumDistance{T}}) where {T} = MinimumDistance(false, 0, 0, typemax(T))
copy(md::MinimumDistance) = MinimumDistance(md.within_cutoff, md.i, md.j, md.d)
convert(::Type{MinimumDistance{T}}, md::MinimumDistance) where {T} = MinimumDistance(md.within_cutoff, md.i, md.j, convert(T, md.d))

#
# useful getter functions
#
"""
    within_cutoff(md::MinimumDistance) = md.within_cutoff

Returns `true` if the distance is within the cutoff, and `false` otherwise.

"""
within_cutoff(md::MinimumDistance) = md.within_cutoff

# The next function adds a method to PDBTools.distance
"""
    distance(md::MinimumDistance) = md.d

Returns the distance between the atoms of the pair.

"""
distance(md::MinimumDistance) = md.d

"""
    iatom(md::MinimumDistance) = md.i

Returns the index of the atom of the first set that is closer to the atom of the second set.

"""
iatom(md::MinimumDistance) = md.i

"""
    jatom(md::MinimumDistance) = md.j

Returns the index of the atom of the second set that is closer to the atom of the first set.

"""
jatom(md::MinimumDistance) = md.j

@testitem "MinimumDistance getter functions" begin
    md = MinimumDistance(true, 2, 5, 1.0f0)
    @test within_cutoff(md) == true
    @test distance(md) == 1.0f0
    @test iatom(md) == 2
    @test jatom(md) == 5
end

import Base.show
function Base.show(io::IO, mime::MIME"text/plain", md::MinimumDistance{T}) where {T}
    print(io, chomp("""
     $(md)

     Distance within cutoff, within_cutoff = $(md.within_cutoff)
     x atom of pair, i = $(md.i)
     y atom of pair, j = $(md.j)
     Distance found, d = $(md.d)
     """))
end

#=
    _mol_indices(i_atom,n_atoms_per_molecule) = (i_atom-1) ÷ n_atoms_per_molecule + 1

Sets the index of the molecule of an atom in the simples situation, in which all 
molecules have the same number of atoms. This is the default setting, and the 
`_mol_indices` parameter of the `minimum_distance` functions must be defined manually
in other situations. 

=#
_mol_indices(i, n_atoms_per_molecule) = (i - 1) ÷ n_atoms_per_molecule + 1
_number_of_molecules(mol_indices, positions) = length(unique(mol_indices(i) for i in eachindex(positions)))

function _get_mol_indices(mol_indices, n_atoms_per_molecule; flag::String="")
    if (isnothing(n_atoms_per_molecule) && isnothing(mol_indices)) ||
       (!isnothing(n_atoms_per_molecule) && !isnothing(mol_indices))
        throw(ArgumentError("Please specify *either* $(flag)n_atoms_per_molecule *or* $(flag)mol_indices"))
    end
    if isnothing(mol_indices)
        mol_indices = i -> _mol_indices(i, n_atoms_per_molecule)
    end
    return mol_indices
end

# Simplify signature of arrays of MinimumDistance and its tuples.
List{T} = Vector{<:MinimumDistance{T}}
ListTuple{T} = Tuple{Vector{<:MinimumDistance{T}},Vector{<:MinimumDistance{T}}}

#=
    init_list(x, mol_indices::F) where F<:Function

Initializes an array of type `Vector{MinimumDistance}` with length equal to the number of 
molecules. `x` must be provided so that the type of variable of the coordinates can be 
propagated to the distances, and `mol_indices` is the function that given an atomic index `i`
returns the index of the molecule.

    init_list(::Type{T}, number_of_molecules::Integer) 

Given the type of the coordinates of the vector of atomic coordinates and the number of molecules,
returns the initialized array of minimum distances. 

# Examples

Providing the type of variable and the number of molecules:

```julia-repl
julia> x = [ rand(SVector{2,Float32}) for _ in 1:100 ]; # 50 molecules of 2 atoms

julia> init_list(Float32,50)
50-element Vector{MinimumDistance{Float32}}:
 MinimumDistance{Float32}(false, -1, -1, Inf32)
 MinimumDistance{Float32}(false, -1, -1, Inf32)
 MinimumDistance{Float32}(false, -1, -1, Inf32)
 ⋮
 MinimumDistance{Float32}(false, -1, -1, Inf32)
 MinimumDistance{Float32}(false, -1, -1, Inf32)

```

Providing the vector of coordinates and a function that returns the index of the
molecule of each element:

```julia-repl
julia> init_list(x, i -> (i-1)÷2 + 1)
50-element Vector{MinimumDistance{Float32}}:
 MinimumDistance{Float32}(false, -1, -1, Inf32)
 MinimumDistance{Float32}(false, -1, -1, Inf32)
 MinimumDistance{Float32}(false, -1, -1, Inf32)
 ⋮
 MinimumDistance{Float32}(false, -1, -1, Inf32)
 MinimumDistance{Float32}(false, -1, -1, Inf32)

```

The above annonymous function `i -> (i-1)÷2 + 1` is equivalent to `i -> mol_indices(i,2)`,
and can be generalized if the the number of atoms of each molecule is not the same.

=#
function init_list(x::AbstractVector{<:AbstractVector}, mol_indices::F) where {F<:Function}
    T = eltype(eltype(x))
    number_of_molecules = 0
    previous_molecule = 0
    for i in eachindex(x)
        imol = mol_indices(i)
        if imol != previous_molecule
            number_of_molecules += 1
            previous_molecule = imol
        end
    end
    return init_list(T, number_of_molecules)
end
init_list(::Type{T}, n::Integer) where {T} = fill(zero(MinimumDistance{T}), n)

#
# Overload of the CellListMap functions that are required for list_threaded
# computations: copy_output, reset_output!, and reducer.
#
# Note that we have two types of output variables here: List, and a tuple of List.
# The List is simply an array of `MinimumDistance{T}`, and we have defined above
# `copy` and `zero` methods for this type, such that we only need to define 
# that reseting this variable consists of returing its zero, and set up the reducer.
# The methods for abstract arrays will take care of the rest.
#
# For the Tuple of lists, we need to be more explicit, and define appropriate copy_output,
# reset_output! and reducer functions.
#
import CellListMap: copy_output, reset_output!, reducer
reset_output!(md::MinimumDistance{T}) where {T} = zero(MinimumDistance{T})
reducer(md1::MinimumDistance{T}, md2::MinimumDistance{T}) where {T} = md1.d < md2.d ? md1 : md2
copy_output(list::ListTuple) = (copy_output(list[1]), copy_output(list[2]))
copy_output(x::MinimumDistance) = copy(x)
function reset_output!(list::ListTuple)
    reset_output!(list[1])
    reset_output!(list[2])
    return list
end
function reducer(l1::ListTuple, l2::ListTuple)
    for i in eachindex(l1[1], l2[1])
        l1[1][i] = reducer(l1[1][i], l2[1][i])
    end
    for i in eachindex(l1[2], l2[2])
        l1[2][i] = reducer(l1[2][i], l2[2][i])
    end
    return l1
end

"""
    minimum_distances!(system)

Function that computes the minimum distances for an initialized system,
of `SelfPairs`, `CrossPairs`, or `AllPairs` types. 

The function returs a `Vector{MinimumDistance}` cor `SelfPairs` and `CrossPairs`
inputs, and a Tuple of two of such vectors for the `AllPairs` input types.

This function is used as an advanced alternative from preallocated system inputs. Only a few allocations 
remain on a call to `minimum_distances!`, mostly related to the launch of the multithreaded 
calculations.

# Example

```julia-repl
julia> using MolSimToolkit, StaticArrays

julia> sys = SelfPairs(
           xpositions = rand(SVector{3,Float64},1000), 
           unitcell=[1,1,1], 
           cutoff = 0.1, 
           xn_atoms_per_molecule=10,
       )
SelfPairs system with:

Number of atoms: 1000
Cutoff: 0.1
unitcell: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
Number of molecules: 100

julia> minimum_distances!(sys)
100-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(true, 8, 579, 0.03570387474690425)
 MinimumDistance{Float64}(true, 12, 534, 0.02850448652684309)
 ⋮
 MinimumDistance{Float64}(true, 996, 423, 0.03655145613454862)

julia> using BenchmarkTools

julia> @btime minimum_distances!(\$sys);
  178.468 μs (209 allocations: 22.80 KiB)
```

"""
function minimum_distances!(sys)
    map_pairwise!(
        (x, y, i, j, d2, list) -> update_list!(i, j, d2, list, sys),
        sys.system
    )
    return sys.minimum_distances
end

"""
    minimum_distances(
       xpositions::AbstractVector{<:SVector},
       # or xpositions *and* ypositions (CrossPairs or AllPairs)
       cutoff=0.1,
       unitcell=[1,1,1],
       xn_atoms_per_molecule=5
       # or xn_atoms_per_molecule (CrossPairs)
       # or xn_atoms_per_molecule *and* yn_atoms_per_molecule (AllPairs)
    )

This function computes directly the minimum distances in a set of particles. 
Depending on the number of input position arrays provided and on the number
of molecular index information provided, a different type of calculation is
performed:

- If `xpositions` and `xn_atoms_per_molecule` are provided, the minimum distances
  within the set of molecules of the set provided are computed. 

- If `xpositions` and `ypositions` are provided, and **only** `xn_atoms_per_molecule`
  is provided, the minimum distance of molecule of set `x` will be computed relative
  to set `y` (or, in other words, `ypositions` are considered a single structure) 

- If `xpositions` and `ypositions` are provided, and `xn_atoms_per_molecule` **and**
  `yn_atoms_per_molecule` are given, the minimum distances of each molecule of `x`
  to any atom of `y` are computed, and vice-versa. A tuple of vectors of minimum
  distances is returned, with lengths corresponding to the number of molecules
  of sets `x` and `y`, respectively.

As for the other functions are constructors, the `xn_atoms_per_molecule` keyword
parameters can be substituted by a general function which returns the molecular 
index of the molecule of each atom (i. e. `(i) -> (i-1)%n_atoms_per_molecule + 1`
in the simplest and default case).

# Examples

## Single set of molecules: all minimum distances within the set

Note that the output contains a vector of `MinimumDistance` elements
with a length equal to the number of molecules of the set.

```julia-repl
julia> using MolSimToolkit, StaticArrays

julia> list = minimum_distances(
           xpositions = rand(SVector{3,Float64},10^5), 
           unitcell=[1,1,1], 
           cutoff = 0.1, 
           xn_atoms_per_molecule=10)
10000-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(true, 5, 71282, 0.007669490894775502)
 MinimumDistance{Float64}(true, 19, 36374, 0.005280726329888545)
 ⋮
 MinimumDistance{Float64}(true, 99998, 44320, 0.006509632622462869)
```

## Two sets: minimum distances of one set relative to the other

Note that the output contains the number of molecules of the `x` set.
For each molecule of this set, the minimum distance to the set `y` is 
computed. This is the typical "solute-solvent" example, where
`x` contains the solvent positions, and `y` contains the solute
positions.

```julia-repl
julia> list = minimum_distances(
           xpositions = rand(SVector{3,Float64},10^5), 
           ypositions = rand(SVector{3,Float64},10^3),
           unitcell=[1,1,1], 
           cutoff = 0.1, 
           xn_atoms_per_molecule=10,
       )
10000-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(true, 5, 596, 0.025526453519907292)
 MinimumDistance{Float64}(true, 18, 391, 0.014114699969628301)
 ⋮
 MinimumDistance{Float64}(true, 99993, 289, 0.016089848937890512)
```

## Two-sets: computing all minimum distances among molecules

If the number of molecules of both sets are provided with the `xn_atoms_per_molecule`
and `yn_atoms_per_molecule` keywords, both sets are split into molecules,
and all minimum distances are computed. For each molecule of each set, 
the minimimum distance to any other molecule of the other set is
returned. The output is a tuple of lists.

```julia-repl
julia> lists = minimum_distances(
           xpositions = rand(SVector{3,Float64},10^5), 
           ypositions = rand(SVector{3,Float64},10^3),
           unitcell=[1,1,1], 
           cutoff = 0.1, 
           xn_atoms_per_molecule=10,
           yn_atoms_per_molecule=100
       );

julia> lists[1]
10000-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(true, 10, 471, 0.03211876310646438)
 MinimumDistance{Float64}(true, 13, 113, 0.0364141004391549)
 ⋮
 MinimumDistance{Float64}(true, 99992, 673, 0.0345818388567913)

julia> lists[2]
10-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(true, 81, 754, 0.002292544732548094)
 MinimumDistance{Float64}(true, 156, 17208, 0.0018147268509811352)
 ⋮
 MinimumDistance{Float64}(true, 944, 98048, 0.002902338025311851)
```

"""
function minimum_distances(;
    xpositions=nothing,
    ypositions=nothing,
    cutoff::Real,
    unitcell::AbstractVecOrMat,
    mol_indices=nothing,
    xmol_indices=nothing,
    ymol_indices=nothing,
    xn_atoms_per_molecule::Union{Nothing,Integer}=nothing,
    yn_atoms_per_molecule::Union{Nothing,Integer}=nothing,
    parallel::Bool=true,
)
    # SelfPairs
    if isnothing(ypositions)
        mol_indices = _get_mol_indices(mol_indices, xn_atoms_per_molecule)
        system = SelfPairs(;
            xpositions=xpositions,
            cutoff=cutoff,
            unitcell=unitcell,
            mol_indices=mol_indices,
            parallel=parallel
        )
        return minimum_distances!(system)
    end
    # CrossPairs
    if isnothing(ymol_indices) && isnothing(yn_atoms_per_molecule)
        xmol_indices = _get_mol_indices(xmol_indices, xn_atoms_per_molecule; flag="x")
        system = CrossPairs(;
            xpositions=xpositions,
            ypositions=ypositions,
            cutoff=cutoff,
            unitcell=unitcell,
            xmol_indices=xmol_indices,
            parallel=parallel
        )
        return minimum_distances!(system)
    end
    # AllPairs
    if !isnothing(xpositions) && (!isnothing(yn_atoms_per_molecule) || !isnothing(ymol_indices))
        xmol_indices = _get_mol_indices(xmol_indices, xn_atoms_per_molecule; flag="x")
        ymol_indices = _get_mol_indices(ymol_indices, yn_atoms_per_molecule; flag="y")
        system = AllPairs(;
            xpositions=xpositions,
            ypositions=ypositions,
            cutoff=cutoff,
            unitcell=unitcell,
            xmol_indices=xmol_indices,
            ymol_indices=ymol_indices,
            parallel=parallel
        )
        return minimum_distances!(system)
    end
    throw(ArgumentError("Incorrect set of keyword input parameters."))
end

abstract type SystemPairs end

import Base: getproperty, propertynames
getproperty(sys::SystemPairs, s::Symbol) = getproperty(sys, Val(s))
getproperty(sys::SystemPairs, ::Val{:system}) = getfield(sys, :system)
getproperty(sys::SystemPairs, ::Val{:mol_indices}) = getfield(sys, :mol_indices)
getproperty(sys::SystemPairs, ::Val{:xmol_indices}) = getfield(sys, :xmol_indices)
getproperty(sys::SystemPairs, ::Val{:ymol_indices}) = getfield(sys, :ymol_indices)
getproperty(sys::SystemPairs, ::Val{:minimum_distances}) = sys.system.output
getproperty(sys::SystemPairs, ::Val{:xpositions}) = sys.system.xpositions
getproperty(sys::SystemPairs, ::Val{:ypositions}) = sys.system.ypositions
getproperty(sys::SystemPairs, ::Val{:cutoff}) = sys.system.cutoff
getproperty(sys::SystemPairs, ::Val{:unitcell}) = sys.system.unitcell
getproperty(sys::SystemPairs, ::Val{:parallel}) = sys.system.parallel
propertynames(sys::SystemPairs, private::Bool) =
    (:system, :mol_indices, :minimum_distances, :xpositions, :ypositions, :unitcell, :cutoff)

import Base: setproperty!
setproperty!(sys::SystemPairs, s::Symbol, value) = setproperty!(sys, Val(s), value)
setproperty!(sys::SystemPairs, ::Val{:cutoff}, cutoff) = update_cutoff!(sys.system, cutoff)
setproperty!(sys::SystemPairs, ::Val{:unitcell}, unitcell) = update_unitcell!(sys.system, unitcell)
setproperty!(sys::SystemPairs, ::Val{:parallel}, parallel) = sys.system.parallel = parallel

#
# Functions for when the lists of minimum-distances is that of a single
# set of molecules (between the molecules of that set)
#
include("./SelfPairs.jl")

#
# Functions for when one wants the list of the atoms of the second set
# that are closer to each molecule of the first set (only one list is 
# returned)
#
include("./CrossPairs.jl")

# 
# Functions for when all pairs of minimum distances are desired,
# between two disjoint sets of molecules
#
include("./AllPairs.jl")

#
# Testing routines
#
include("./testing.jl")

end # MolecularMinimumDistances
