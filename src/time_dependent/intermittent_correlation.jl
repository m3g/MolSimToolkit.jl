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
"""
    CorrelationProfile

Structure that wraps the result of the [`intermittent_correlation_profile`](@ref)
function: a set of intermittent correlation functions, one for each shell
(bin) of distances to the binding site.

# Fields

- `r::Vector{T}`: the center of each distance bin.
- `correlations::Vector{OffsetVector{Float64,Vector{Float64}}}`: for each bin,
  the intermittent correlation function, an `OffsetArray` with indices
  `0:maxdelta`. Bins for which no solvent molecule was ever found contain
  `NaN` values.
- `counts::Vector{Int}`: for each bin, the number of (frame, solvent molecule)
  observations found in that bin. Useful to discard poorly sampled bins.
- `delta_r::T`: the width of each bin.
- `step_r::T`: the displacement between the lower edges of consecutive bins.
- `dmax::T`: the maximum distance considered, inherited from the input data.

!!! compat
    This structure was added in version 2.4.0 of MolSimToolkit.

"""
struct CorrelationProfile{T<:Real}
    r::Vector{T}
    correlations::Vector{OffsetArrays.OffsetVector{Float64,Vector{Float64}}}
    counts::Vector{Int}
    delta_r::T
    step_r::T
    dmax::T
end

function Base.show(io::IO, ::MIME"text/plain", p::CorrelationProfile)
    maxdelta = length(first(p.correlations)) - 1
    print(io, chomp(
        """
        -------------------------------------------------------------------
        Intermittent correlation profile:
        -------------------------------------------------------------------
        Number of distance bins: $(length(p.r))
        Bin width (delta_r): $(p.delta_r)
        Bin step (step_r): $(p.step_r)
        Maximum distance (dmax): $(p.dmax)
        Range of bin centers: $(first(p.r)) - $(last(p.r))
        Maximum delta (frames): $maxdelta
        Observations per bin: minimum = $(minimum(p.counts)), maximum = $(maximum(p.counts))
        -------------------------------------------------------------------
        """
    ))
end

"""
    intermittent_correlation_profile(
        occ::Occupancy;
        delta_r::Real,
        step_r::Real = delta_r / 10,
        maxdelta::Integer = length(occ.list) ÷ 10,
        show_progress::Bool = true,
    )

Computes the intermittent correlation function of the site occupancy resolved
by the distance between the site and the solvent molecules, as computed by the
[`occupancy`](@ref) function.

The distances from `0` to `occ.dmax` are split into bins of constant width
`delta_r`, with the lower edge of consecutive bins displaced by `step_r`. Since
`step_r` is typically smaller than `delta_r`, the bins overlap, and the profile
is a smooth (quasi-continuous) function of the distance.

For each bin, the correlation function at `delta` is the probability that a
solvent molecule found in that bin at frame `i` is found in the *same* bin at
frame `i + delta`. The molecule does not need to remain in the bin in the
frames in between: it may leave and come back (that is what makes this
correlation function "intermittent"). Thus, each bin characterizes how long
a solvent molecule remains at that distance from the site, and the decay of
the correlation function can be converted into a residence time with the
[`residence_time`](@ref) function.

# Arguments

- `occ::Occupancy`: The result of the `occupancy` function. It must carry
  distance information, that is, it must have been produced by `occupancy`
  and not built from a list of molecule indices alone.

# Keyword arguments

- `delta_r::Real`: The width of each distance bin. Required.
- `step_r::Real`: The displacement between the lower edges of consecutive
  bins. Defaults to `delta_r / 10`.
- `maxdelta::Integer`: The maximum delta-step to be considered. Defaults to
  `length(occ.list) ÷ 10`.
- `show_progress::Bool`: Show progress bar. Defaults to `true`.

# Returns

- A [`CorrelationProfile`](@ref) object.

# Example

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2***"
julia> using MolSimToolkit, PDBTools, MolSimToolkit.Testing

julia> sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj);

julia> protein = select(get_atoms(sim), "protein");

julia> tmao = select(get_atoms(sim), "resname TMAO");

julia> occ = occupancy(sim, protein, tmao; solvent_natomspermol=14, cutoff=6.0, show_progress=false);

julia> p = intermittent_correlation_profile(occ; delta_r=1.0, step_r=0.5, maxdelta=4, show_progress=false)
-------------------------------------------------------------------
Intermittent correlation profile:
-------------------------------------------------------------------
Number of distance bins: 11
Bin width (delta_r): 1.0
Bin step (step_r): 0.5
Maximum distance (dmax): 6.0
Range of bin centers: 0.5 - 5.5
Maximum delta (frames): 4
Observations per bin: minimum = 0, maximum = 89
-------------------------------------------------------------------

julia> p.r[end] # center of the outermost distance bin
5.5

julia> p.correlations[end][0:2] # correlation of the outermost bin
3-element Vector{Float64}:
 1.0
 0.05357142857142857
 0.019230769230769232

```

"""
function intermittent_correlation_profile(
    occ::Occupancy;
    delta_r::Real,
    step_r::Real=delta_r / 10,
    maxdelta::Integer=max(1, length(occ.list) ÷ 10),
    show_progress::Bool=true,
)
    if isnan(occ.dmax)
        throw(ArgumentError("""\n
            The Occupancy object does not contain distance information.

            Use the `occupancy` function to compute the occupancy, such that the
            site-solvent distances and `dmax` are stored in the result.

        """))
    end
    return _intermittent_correlation_profile(
        occ.list, occ.distances, occ.dmax;
        delta_r, step_r, maxdelta, show_progress
    )
end

#=
    _intermittent_correlation_profile(list, distances, dmax; kargs...)

Implementation of `intermittent_correlation_profile`. `list[iframe]` is the
list of identities present at frame `iframe`, and `distances[iframe][i]` is
the distance associated to `list[iframe][i]`.
=#
function _intermittent_correlation_profile(
    list::AbstractVector,
    distances::AbstractVector,
    dmax::Real;
    delta_r::Real,
    step_r::Real,
    maxdelta::Integer,
    show_progress::Bool,
)
    if !(0 < delta_r <= dmax)
        throw(ArgumentError("delta_r must be positive and not greater than dmax = $dmax"))
    end
    if !(0 < step_r)
        throw(ArgumentError("step_r must be positive"))
    end
    if maxdelta > length(list) - 1
        throw(ArgumentError("maxdelta must be less than the length of the data minus 1"))
    end
    delta_r, step_r, dmax = promote(float(delta_r), float(step_r), float(dmax))
    # Lower edges of the bins: all bins are fully contained in [0, dmax]
    lower_edges = collect(zero(step_r):step_r:(dmax-delta_r))
    nbins = length(lower_edges)
    r = [lower + delta_r / 2 for lower in lower_edges]
    correlations = Vector{OffsetArrays.OffsetVector{Float64,Vector{Float64}}}(undef, nbins)
    counts = zeros(Int, nbins)
    prg = Progress(nbins; enabled=show_progress)
    for ibin in 1:nbins
        lower = lower_edges[ibin]
        upper = lower + delta_r
        # For each frame, which identities are within this distance bin
        bin_list = [
            [id for (id, d) in zip(list[iframe], distances[iframe]) if lower <= d < upper]
            for iframe in eachindex(list)
        ]
        counts[ibin] = sum(length, bin_list)
        identities = unique(Iterators.flatten(bin_list))
        correlations[ibin] = _list_intermittent_correlation(
            bin_list, identities; maxdelta, show_progress=false
        )
        next!(prg)
    end
    return CorrelationProfile(r, correlations, counts, delta_r, step_r, dmax)
end

@testitem "intermittent_correlation_profile" begin
    using MolSimToolkit, PDBTools, MolSimToolkit.Testing
    using ShowMethodTesting

    # A single molecule sitting still at a distance of 1.5: it is always in the
    # bins that contain 1.5, and never in the others
    list = [[1] for _ in 1:11]
    distances = [[1.5] for _ in 1:11]
    occ = Occupancy(list, distances, 1, 3.0)
    p = intermittent_correlation_profile(occ; delta_r=1.0, step_r=0.5, maxdelta=4, show_progress=false)
    @test p.r == [0.5, 1.0, 1.5, 2.0, 2.5]
    # 1.5 belongs to the bins [1.0,2.0) and [1.5,2.5) only
    @test p.counts == [0, 0, 11, 11, 0]
    @test all(isnan, p.correlations[1])
    @test all(isnan, p.correlations[2])
    @test parent(p.correlations[3]) == ones(5)
    @test parent(p.correlations[4]) == ones(5)
    @test all(isnan, p.correlations[5])

    # A molecule that alternates between two bins
    distances = [[iseven(i) ? 0.5 : 2.5] for i in 1:11]
    occ = Occupancy(list, distances, 1, 3.0)
    p = intermittent_correlation_profile(occ; delta_r=1.0, step_r=1.0, maxdelta=4, show_progress=false)
    @test p.r == [0.5, 1.5, 2.5]
    @test p.counts == [5, 0, 6]
    @test parent(p.correlations[1]) == [1.0, 0.0, 1.0, 0.0, 1.0]
    @test all(isnan, p.correlations[2])
    @test parent(p.correlations[3]) == [1.0, 0.0, 1.0, 0.0, 1.0]

    sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj)
    protein = select(get_atoms(sim), "protein")
    tmao = select(get_atoms(sim), "resname TMAO")
    occ = occupancy(sim, protein, tmao; solvent_natomspermol=14, cutoff=6.0, show_progress=false)
    p = intermittent_correlation_profile(occ; delta_r=1.0, step_r=0.5, maxdelta=4, show_progress=false)
    @test length(p.r) == 11
    @test p.r ≈ 0.5:0.5:5.5
    @test p.dmax == 6.0
    @test p.delta_r == 1.0
    @test p.step_r == 0.5
    @test p.counts == [0, 0, 17, 70, 89, 60, 44, 48, 58, 66, 61]
    # bins with data start at 1.0, empty bins are all NaN
    @test all(i -> p.counts[i] == 0 ? all(isnan, p.correlations[i]) : p.correlations[i][0] == 1.0, eachindex(p.r))
    @test all(c -> all(x -> isnan(x) || 0.0 <= x <= 1.0, c), p.correlations)
    # correlations must decay with delta
    @test all(i -> p.counts[i] == 0 || p.correlations[i][4] < p.correlations[i][0], eachindex(p.r))
    @test parse_show(p) ≈ """
    -------------------------------------------------------------------
    Intermittent correlation profile:
    -------------------------------------------------------------------
    Number of distance bins: 11
    Bin width (delta_r): 1.0
    Bin step (step_r): 0.5
    Maximum distance (dmax): 6.0
    Range of bin centers: 0.5 - 5.5
    Maximum delta (frames): 4
    Observations per bin: minimum = 0, maximum = 89
    -------------------------------------------------------------------
    """

    # default step_r
    p2 = intermittent_correlation_profile(occ; delta_r=1.0, maxdelta=4, show_progress=false)
    @test p2.step_r == 0.1
    @test length(p2.r) == 51

    @test_throws ArgumentError intermittent_correlation_profile(occ; delta_r=7.0)
    @test_throws ArgumentError intermittent_correlation_profile(occ; delta_r=-1.0)
    @test_throws ArgumentError intermittent_correlation_profile(occ; delta_r=1.0, step_r=0.0)
    @test_throws ArgumentError intermittent_correlation_profile(occ; delta_r=1.0, maxdelta=length(occ.list))
    @test_throws ArgumentError intermittent_correlation_profile(
        Occupancy(occ.list, occ.n_solvent_molecules); delta_r=1.0
    )
end

"""
    residence_time(c::AbstractVector; threshold::Real = 0.5, dt::Real = 1)
    residence_time(p::CorrelationProfile; threshold::Real = 0.5, dt::Real = 1)

Returns the time at which an intermittent correlation function falls below
`threshold`. For `threshold = 0.5` (the default), this is the half-life of the
correlation, that is, the time after which half of the solvent molecules
initially present are no longer there.

When applied to a [`CorrelationProfile`](@ref), returns one residence time for
each distance bin, that is, the residence time as a function of the distance to
the site.

The time is obtained by linear interpolation between the two consecutive
delta-steps that bracket the threshold, and is thus not necessarily an integer
number of frames. `NaN` is returned if the correlation function does not fall
below `threshold` within the `maxdelta` of the input data (or if the bin
contains no data at all).

# Keyword arguments

- `threshold::Real`: The correlation value defining the residence time.
  Defaults to `0.5`.
- `dt::Real`: The time interval between consecutive frames of the input data,
  used to convert the result from frames to time units. Defaults to `1`, that
  is, the result is given in frames.

# Example

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2***"
julia> using MolSimToolkit, PDBTools, MolSimToolkit.Testing

julia> sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj);

julia> protein = select(get_atoms(sim), "protein");

julia> tmao = select(get_atoms(sim), "resname TMAO");

julia> occ = occupancy(sim, protein, tmao; solvent_natomspermol=14, cutoff=6.0, show_progress=false);

julia> ic = intermittent_correlation(occ; maxdelta=4, show_progress=false);

julia> residence_time(ic) # in frames
0.9035714285714286

julia> p = intermittent_correlation_profile(occ; delta_r=1.0, step_r=0.5, maxdelta=4, show_progress=false);

julia> residence_time(p)[end] # outermost distance bin
0.5283018867924528

```

!!! compat
    This function was added in version 2.4.0 of MolSimToolkit.

"""
function residence_time(c::AbstractVector; threshold::Real=0.5, dt::Real=1)
    ifirst = firstindex(c)
    isempty(c) && return NaN * dt
    isnan(c[ifirst]) && return NaN * dt
    c[ifirst] < threshold && return zero(float(dt))
    for i in (ifirst+1):lastindex(c)
        if c[i] < threshold
            # linear interpolation between the previous and the current step
            frac = (c[i-1] - threshold) / (c[i-1] - c[i])
            return (i - 1 - ifirst + frac) * dt
        end
    end
    return NaN * dt
end

function residence_time(p::CorrelationProfile; threshold::Real=0.5, dt::Real=1)
    return [residence_time(c; threshold, dt) for c in p.correlations]
end

@testitem "residence_time" begin
    using MolSimToolkit, PDBTools, MolSimToolkit.Testing
    using OffsetArrays

    # the crossing is exactly at delta = 1, no interpolation needed
    @test residence_time(OffsetArray([1.0, 0.5, 0.0], 0:2)) ≈ 1.0
    # interpolated crossings
    @test residence_time(OffsetArray([1.0, 1.0, 0.0], 0:2)) ≈ 1.5
    @test residence_time(OffsetArray([1.0, 0.0], 0:1)) ≈ 0.5
    # never falls below the threshold
    @test isnan(residence_time(OffsetArray([1.0, 1.0, 1.0], 0:2)))
    # empty bins
    @test isnan(residence_time(OffsetArray([NaN, NaN], 0:1)))
    @test isnan(residence_time(Float64[]))
    # other thresholds and time units
    @test residence_time(OffsetArray([1.0, 0.8, 0.2], 0:2); threshold=0.8) ≈ 1.0
    @test residence_time(OffsetArray([1.0, 0.5, 0.0], 0:2); dt=2.0) ≈ 2.0
    @test residence_time(OffsetArray([1.0, 0.5, 0.0], 0:2); threshold=0.25) ≈ 1.5
    # works on 1-based vectors as well
    @test residence_time([1.0, 0.5, 0.0]) ≈ 1.0

    sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj)
    protein = select(get_atoms(sim), "protein")
    tmao = select(get_atoms(sim), "resname TMAO")
    occ = occupancy(sim, protein, tmao; solvent_natomspermol=14, cutoff=6.0, show_progress=false)
    p = intermittent_correlation_profile(occ; delta_r=1.0, step_r=0.5, maxdelta=4, show_progress=false)
    t = residence_time(p)
    @test length(t) == length(p.r)
    @test all(x -> isnan(x) || x >= 0, t)
    # empty bins produce NaN residence times
    @test all(i -> p.counts[i] > 0 || isnan(t[i]), eachindex(t))
    # in time units
    @test all(i -> isnan(t[i]) || residence_time(p; dt=2)[i] ≈ 2 * t[i], eachindex(t))
end
