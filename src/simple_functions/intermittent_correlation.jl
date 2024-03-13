"""
    intermittent_correlation(data::AbstractVector; maxdelta = length(data) ÷ 10)

Calculate the intermittent correlation function of a time series. That is,
computes the probability of finding a value of the same type at a step
`i + delta` in the time series, given that it was present in step `i`.

Returns an `OffsetArray` with indices `0:maxdelta`, where the value at position
`0` is `1.0`, corresponding to the normalized count of events. 

# Arguments

- `data::AbstractVector`: The time series to be analyzed. 
- `maxdelta::Int`: The maximum delta-step to be considered. Defaults to 
  `length(data) ÷ 10`.

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
```

!!! compat
    This function was added in version 1.9.0 of MolSimToolkit.

"""
function intermittent_correlation(
    data::AbstractVector; 
    maxdelta::Integer = max(1, length(data) ÷ 10)
)
    types = unique(data)
    counts = OffsetArrays.OffsetArray(zeros(maxdelta+1), 0:maxdelta)
    for type in types
        positions = findall(x -> isequal(x, type), data)
        np = length(positions)
        for i in 1:np, j in i:np
            delta = positions[j] - positions[i]
            if delta <= maxdelta
                counts[delta] += 1
            end
        end
    end
    for i in 0:maxdelta
        counts[i] /= length(data) - i
    end
    return counts
end

@testitem "intermittent_correlation" begin
    using MolSimToolkit
    data = [ mod(i,2) for i in 1:10^3 ];
    c = intermittent_correlation(data)
    @test all(==(1), c[0:2:end])  
    @test all(==(0), c[1:2:end])  
    @test length(c) == 101
    c = intermittent_correlation(data; maxdelta=4)
    @test length(c) == 5
    data = [1]
    for i in 2:101
        push!(data, data[i-1] + (mod(i,2) == 0))
    end
    c = intermittent_correlation(data)
    @test c[0:1] == [1.0, 0.5]
    @test all(==(0), c[2:end])
end