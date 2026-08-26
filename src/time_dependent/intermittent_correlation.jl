"""
    intermittent_correlation(
        data::AbstractVector; 
        maxdelta = length(data) ÷ 10, 
        types::Function = x -> true,
        show_progress::Bool = true,
    )

Calculate the intermittent correlation function of a time series. That is,
computes the probability of finding a value of the same type at a step
`i + delta` in the time series, given that it was present in step `i`.

Returns an `OffsetArray` with indices `0:maxdelta`, where the value at position
`0` is `1.0`, corresponding to the normalized count of events. 

# Arguments

- `data::AbstractVector`: The time series to be analyzed. 
- `maxdelta::Integer`: The maximum delta-step to be considered. Defaults to 
  `length(data) ÷ 10`.
- `types` (optional): A function that returns `true` for the types of data
   that should be considered. Defaults to all data, i. e. `x -> true`. For 
   example, to ignore `0` values, use `types = x -> x != 0`.  
- `show_progress::Bool`: Show progress bar. Defaults to `true`.

# Examples

Here we produce a time-series of 10,000 elements, as a sequence of 
1's and 0's (`[1, 0, 1, 0, ...]`), and calculate the intermittent correlation function.
The probability of finding the same number (0 or 1) after odd steps is 0, and
the probability of finding the same number after even steps is 1.

```jldoctest
julia> using MolSimToolkit

julia> data = [ mod(i,2) for i in 1:10^4 ];

julia> intermittent_correlation(data; maxdelta=4, show_progress=false)
5-element OffsetArray(::Vector{Float64}, 0:4) with eltype Float64 with indices 0:4:
 1.0
 0.0
 1.0
 0.0
 1.0

julia> intermittent_correlation(data; maxdelta=4, types = x -> x != 0, show_progress=false)
5-element OffsetArray(::Vector{Float64}, 0:4) with eltype Float64 with indices 0:4:
 1.0
 0.0
 1.0
 0.0
 1.0
```

In the second run, we have ignored the `0` values, and the result is the same, 
because here the correlations of the `1` values are the same as the correlations
of the `0` values.

!!! compat
    This function was added in version 1.9.0 of MolSimToolkit. The `types` argument
    was added in version 1.10.0 and the `show_progress` argument in version 1.28.0.

"""
function intermittent_correlation(
    data::AbstractVector;
    maxdelta::Integer=max(1, length(data) ÷ 10),
    types::F=x -> true,
    show_progress::Bool=true,
) where {F<:Function}
    if maxdelta > length(data) - 1
        throw(ArgumentError("maxdelta must be less than the length of the data minus 1"))
    end
    types_considered = filter!(types, unique(data))
    counts = OffsetArrays.OffsetArray(zeros(maxdelta + 1), 0:maxdelta)
    chances = copy(counts)
    np_all = 0
    for itype in types_considered
        np = count(x -> isequal(x, itype), data)
        np_all += (np * (np - 1)) ÷ 2
    end
    p = Progress(np_all; enabled=show_progress)
    for type in types_considered
        positions = findall(x -> isequal(x, type), data)
        np = length(positions)
        for i in 1:np
            delta = 0
            while positions[i] + delta <= length(data)
                chances[delta] += 1
                delta += 1
                delta > maxdelta && break
            end
            for j in i:np
                delta = positions[j] - positions[i]
                if delta <= maxdelta
                    counts[delta] += 1
                end
                next!(p)
            end
        end
    end
    # Convert counts to probabilities
    counts ./= chances
    return counts
end

@testitem "intermittent_correlation" begin
    using MolSimToolkit
    using OffsetArrays
    data = [1, 0, 1, 0, 1]
    c = intermittent_correlation(data; maxdelta=4)
    @test c == OffsetArray([1.0, 0.0, 1.0, 0.0, 1.0], 0:4)
    c = intermittent_correlation(data; types=x -> x != 0)
    @test c == OffsetArray([1.0, 0.0], 0:1)
    c = intermittent_correlation(data; maxdelta=4, types=x -> x != 0)
    @test c == OffsetArray([1.0, 0.0, 1.0, 0.0, 1.0], 0:4)
    data = [0]
    for i in 2:101
        push!(data, data[i-1] + (mod(i, 2) == 0))
    end
    c = intermittent_correlation(data; show_progress=false)
    @test c[0:1] == [1.0, 0.5]
    c = intermittent_correlation(data)
    @test c[0:1] == [1.0, 0.5]
    @test all(==(0), c[2:end])
    data = [0, 1, 1, 2, 2]
    @test intermittent_correlation(data) == OffsetArray([1.0, 0.5], 0:1)
    @test intermittent_correlation(data; types=x -> x > 0) ≈ OffsetArray([1.0, 2 / 3], 0:1)
    @test intermittent_correlation(data; types=x -> x != 1) == OffsetArray([1.0, 0.5], 0:1)
    @test intermittent_correlation(data; types=x -> x == 1) == OffsetArray([1.0, 0.5], 0:1)
    @test_throws ArgumentError intermittent_correlation(data; maxdelta=5)
end

"""
    intermittent_correlation(
        occupancy::Occupancy;
        maxdelta::Integer = length(occupancy.list) ÷ 10,
        show_progress::Bool = true,
    )

Calculate the intermittent correlation function of the occupancy of a binding
site, as computed by the [`occupancy`](@ref) function. That is, computes the
probability of finding a solvent molecule at the site at frame `i + delta`,
given that it was found at the site at frame `i`, for each solvent molecule
independently.

Returns an `OffsetArray` with indices `0:maxdelta`, where the value at position
`0` is `1.0`, corresponding to the normalized count of events.

# Arguments

- `occupancy::Occupancy`: The result of the `occupancy` function.
- `maxdelta::Integer`: The maximum delta-step to be considered. Defaults to
  `length(occupancy.list) ÷ 10`.
- `show_progress::Bool`: Show progress bar. Defaults to `true`.

# Example

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2***"
julia> using MolSimToolkit, PDBTools, MolSimToolkit.Testing

julia> sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj);

julia> protein = select(get_atoms(sim), "protein");

julia> tmao = select(get_atoms(sim), "resname TMAO");

julia> occ = occupancy(sim, protein, tmao; solvent_natomspermol=14, cutoff=3.0, show_progress=false);

julia> c = intermittent_correlation(occ; maxdelta=4, show_progress=false);

julia> c
5-element OffsetArray(::Vector{Float64}, 0:4) with eltype Float64 with indices 0:4:
 1.0
 0.39
 0.23469387755102042
 0.12903225806451613
 0.033707865168539325

```

"""
function intermittent_correlation(
    occ::Occupancy;
    maxdelta::Integer=max(1, length(occ.list) ÷ 10),
    show_progress::Bool=true,
)
    return _list_intermittent_correlation(occ.list, 1:occ.n_solvent_molecules; maxdelta, show_progress)
end

#=
    _list_intermittent_correlation(data, identities; maxdelta, show_progress) -> OffsetArray

Shared implementation behind `intermittent_correlation(occ::Occupancy)` and
`intermittent_correlation(hbo::HydrogenBondOccupancy)`. `data` is a vector
with, for each frame, the list of identities present in that frame (solvent
molecule indices, for `Occupancy`; `HBond`s, for `HydrogenBondOccupancy`).
`identities` is the collection of all distinct identities to consider. Only
`==` (and a consistent `hash`) is required from the identity type, since
membership is tested with `in`.
=#
function _list_intermittent_correlation(
    data::AbstractVector,
    identities;
    maxdelta::Integer=max(1, length(data) ÷ 10),
    show_progress::Bool=true,
)
    ndata = length(data)
    if maxdelta > ndata - 1
        throw(ArgumentError("maxdelta must be less than the length of the data minus 1"))
    end
    counts = OffsetArrays.OffsetArray(zeros(maxdelta + 1), 0:maxdelta)
    chances = copy(counts)
    np_all = 0
    for id in identities
        np = count(list -> id in list, data)
        np_all += (np * (np - 1)) ÷ 2
    end
    p = Progress(np_all; enabled=show_progress)
    for id in identities
        positions = findall(list -> id in list, data)
        np = length(positions)
        for i in 1:np
            delta = 0
            while positions[i] + delta <= ndata
                chances[delta] += 1
                delta += 1
                delta > maxdelta && break
            end
            for j in i:np
                delta = positions[j] - positions[i]
                if delta <= maxdelta
                    counts[delta] += 1
                end
                next!(p)
            end
        end
    end
    # Convert counts to probabilities
    counts ./= chances
    return counts
end

"""
    intermittent_correlation(
        hbo::HydrogenBondOccupancy;
        maxdelta::Integer = length(hbo.list) ÷ 10,
        show_progress::Bool = true,
    )

Calculate the intermittent correlation function of the hydrogen bonds found
by the [`hydrogen_bond_occupancy`](@ref) function. That is, computes the
probability of finding a given hydrogen bond (the same donnor, polar
hydrogen, and acceptor atoms) present at frame `i + delta`, given that it was
present at frame `i`, for each hydrogen bond independently. The hydrogen
bond does not need to remain present in every frame in between: it may break
and reform within the interval (that is what makes this correlation
function "intermittent", as opposed to "continuous").

Returns an `OffsetArray` with indices `0:maxdelta`, where the value at position
`0` is `1.0`, corresponding to the normalized count of events.

# Arguments

- `hbo::HydrogenBondOccupancy`: The result of the `hydrogen_bond_occupancy` function.
- `maxdelta::Integer`: The maximum delta-step to be considered. Defaults to
  `length(hbo.list) ÷ 10`.
- `show_progress::Bool`: Show progress bar. Defaults to `true`.

# Example

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2***"
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> sim = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> hbo = hydrogen_bond_occupancy(sim, "protein", show_progress=false)["protein => protein"];

julia> c = intermittent_correlation(hbo; maxdelta=4, show_progress=false);

julia> c[0]
1.0

```

"""
function intermittent_correlation(
    hbo::HydrogenBondOccupancy;
    maxdelta::Integer=max(1, length(hbo.list) ÷ 10),
    show_progress::Bool=true,
)
    return _list_intermittent_correlation(hbo.list, unique_hbonds(hbo); maxdelta, show_progress)
end

@testitem "intermittent_correlation - HydrogenBondOccupancy" begin
    using MolSimToolkit, MolSimToolkit.Testing
    using OffsetArrays

    # A hydrogen bond present in every frame must have full correlation
    b = HBond(1, 2, 3)
    hbo = HydrogenBondOccupancy([[b], [b], [b], [b], [b]])
    c = intermittent_correlation(hbo; maxdelta=4, show_progress=false)
    @test c == OffsetArray([1.0, 1.0, 1.0, 1.0, 1.0], 0:4)

    # A different bridging hydrogen counts as a different bond: alternating
    # between the two behaves like the [1, 0, 1, 0, 1] scalar case
    b1 = HBond(1, 2, 3)
    b2 = HBond(1, 4, 3)
    hbo2 = HydrogenBondOccupancy([[b1], [b2], [b1], [b2], [b1]])
    c2 = intermittent_correlation(hbo2; maxdelta=4, show_progress=false)
    @test c2 == OffsetArray([1.0, 0.0, 1.0, 0.0, 1.0], 0:4)

    sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj)
    hbo3 = hydrogen_bond_occupancy(sim, "protein"; show_progress=false)["protein => protein"]
    c3 = intermittent_correlation(hbo3; maxdelta=4, show_progress=false)
    @test c3[0] == 1.0
    @test length(c3) == 5
    @test all(x -> 0.0 <= x <= 1.0, c3)
    @test_throws ArgumentError intermittent_correlation(hbo3; maxdelta=length(hbo3.list))
end

@testitem "intermittent_correlation - Occupancy" begin
    using MolSimToolkit, PDBTools, MolSimToolkit.Testing
    using OffsetArrays

    # Equivalence with the scalar-type version, for a series with a single
    # solvent molecule occupying the site at each frame
    data = [1, 0, 1, 0, 1]
    occ = Occupancy([[v + 1] for v in data], 2)
    c = intermittent_correlation(occ; maxdelta=4, show_progress=false)
    @test c == intermittent_correlation(data; maxdelta=4, show_progress=false)
    @test c == OffsetArray([1.0, 0.0, 1.0, 0.0, 1.0], 0:4)

    sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj)
    protein = select(get_atoms(sim), "protein")
    tmao = select(get_atoms(sim), "resname TMAO")
    occ2 = occupancy(sim, protein, tmao; solvent_natomspermol=14, cutoff=3.0, show_progress=false)
    c2 = intermittent_correlation(occ2; maxdelta=4, show_progress=false)
    @test c2[0] == 1.0
    @test length(c2) == 5
    @test all(x -> 0.0 <= x <= 1.0, c2)
    @test_throws ArgumentError intermittent_correlation(occ2; maxdelta=length(occ2.list))
end