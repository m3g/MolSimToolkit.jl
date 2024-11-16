"""
    intermittent_correlation(
        data::AbstractVector; 
        maxdelta = length(data) ÷ 10, 
        types::Function = x -> true,
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

# Examples

Here we produce a time-series of 10,000 elements, as a sequence of 
1's and 0's (`[1, 0, 1, 0, ...]`), and calculate the intermittent correlation function.
The probability of finding the same number (0 or 1) after odd steps is 0, and
the probability of finding the same number after even steps is 1.

```jldoctest
julia> using MolSimToolkit

julia> data = [ mod(i,2) for i in 1:10^4 ];

julia> intermittent_correlation(data; maxdelta=4)
5-element OffsetArray(::Vector{Float64}, 0:4) with eltype Float64 with indices 0:4:
 1.0
 0.0
 1.0
 0.0
 1.0

julia> intermittent_correlation(data; maxdelta=4, types = x -> x != 0)
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
    was added in version 1.10.0.

"""
function intermittent_correlation(
    data::AbstractVector;
    maxdelta::Integer=max(1, length(data) ÷ 10),
    types::F=x -> true,
) where {F<:Function}
    if maxdelta > length(data) - 1
        throw(ArgumentError("maxdelta must be less than the length of the data minus 1"))
    end
    types_considered = filter!(types, unique(data))
    counts = OffsetArrays.OffsetArray(zeros(maxdelta + 1), 0:maxdelta)
    chances = copy(counts)
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