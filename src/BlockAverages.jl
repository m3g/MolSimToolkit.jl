module BlockAverages

import LinearAlgebra: dot
import DocStringExtensions: TYPEDEF, TYPEDFIELDS
import Statistics: mean, std
import StatsBase: autocor
import EasyFit: fitexp, upper, lower
export BlockAverageData, BlockDistribution
export block_average, block_distribution

"""

$(TYPEDEF)

Structure that contains the result of the block-average analysis of the sequence. 

`x` is the original data set. 

`xmean` is the property value computed for all data (usually the mean, but not necessarily)

`blocksize` is an array of block sizes, in which the data was split. 
By default it goes from `1` to `length(x)`, with a number of points corresponding to the
number of integer divisions of `length(x)`.

`xmean_maxerr`: The property is computed for each block, and the maximum error (difference
between the property in the block and the average property) is stored in this array, for
each blocks size. 

`xmean_stderr`: The standard error of the estimates of the property, meaning the standar 
deviation of the estimates divided by the square root of the number of blocks. 

`autocor`: Is the autocorrelation function of the data, as a function of the lag. 

`lags`: Is the set of "time" lags for which the autocorrelation will be computed. Defined by
the `lags` parameter of the `block_average` function, as a range. 

`tau`: The characteristic decay time of the autocorrelation function, as obtained by fitting
of a single exponential, of the form `exp(-t/tau)` to the data. 

"""
struct BlockAverageData{T}
    x::Vector{T}
    xmean::Float64
    blocksize::Vector{Int}
    xmean_maxerr::Vector{Float64}
    xmean_stderr::Vector{Float64}
    lags::AbstractVector{Int}
    autocor::Vector{Float64}
    tau::Float64
end

function _print_block_sizes(blocksize)
    if length(blocksize) < 5
        "[$(join(blocksize, ", "))]"
    else
        "[$(blocksize[1]), $(blocksize[2]), ..., $(blocksize[end-1]), $(blocksize[end])]"
    end
end

function Base.show(io::IO, ::MIME"text/plain", b::BlockAverageData)
    merr = findmax(b.xmean_stderr)
    izerolag = findfirst(x -> x <= 0, b.autocor)
    izerolag = isnothing(izerolag) ? 1 : izerolag
    print(io, chomp(
        """
        -------------------------------------------------------------------
        $(typeof(b))
        -------------------------------------------------------------------
        Estimated value (mean by default) = $(b.xmean)
        Length of data series: $(length(b.x))

        Block sizes: $(_print_block_sizes(b.blocksize))

        Maximum standard error (error, block size): $((merr[1], b.blocksize[merr[2]]))

        Deviations in last 3 blocks:
                 percentual: $((100/b.xmean)*(b.xmean_maxerr[max(1,lastindex(b.xmean_maxerr)-2):end] .- b.xmean))  
                   absolute: $((b.xmean_maxerr[max(1,lastindex(b.xmean_maxerr)-2):end] .- b.xmean))  

        Autocorrelation is first zero with lag: $(b.lags[izerolag])
        Characteristic time of autocorrelation decay: 
                as fraction of series length: $(b.tau / length(b.x))
                                    absolute: $(b.tau)
        -------------------------------------------------------------------
        """
    ))
end

#
# Checks if the length of the data is a multiple of block_size, and if not adjust
# the length and print a warning
#
function adjust_xinput(x_input, block_size, var="max_block_size")
    if length(x_input) < block_size
        throw(ArgumentError("number of data points must be greater than $var"))
    end
    if length(x_input) % block_size != 0
        x = @view(x_input[firstindex(x_input):lastindex(x_input)-length(x_input)%block_size])
        println("""

        WARNING: number of data points is not a multiple of $var.

                This may cause poor block sampling, because the analysis
                is performed only for sets of blocks that encompass to complete
                data set. 

                >> Only the first $(length(x)) data points will be considered.

        """)
    else
        x = x_input
    end
    return x
end

"""
    block_average(
        x::AbstractVector{T};
        by = mean,
        min_block_size::Integer = 1,
        max_block_size::Integer = length(x),
        lags::Union{Nothing,AbstractVector{<:Integer}} = nothing,
    ) where {T<:Real}

This function peforms some convergence analysis for a property computed from a series of data, typically a time-series. 
The data is given in vector `x`, and `by` defines the property to be estimated, typically, and by default, the mean value.

Two analyses are performed: a block averaging, in which the data is split in to blocks, and the mean value (or `by` value)
in each block is computed idependently. The output will contain the worst estimate obtained for all blocks, and the
standard error of the estimates, as a function of the block size. 

Finally, the autocorrelation function of the data is computed, and a single exponential is fitted, to obtain the
characteristic time of the decay of the correlation. 

The output will be a structure of type `BlockAverageData{T}`. See the corresponding help entry for more information.

All results can be plot with a convenience function `BlockAverage.plot`

The `lags` keyword can be tuned to define the range of intervals and length of the autocorrelation calculation, with
important implications to the exponential fit and correlation curve shape. See the `StatsBase.autocor` help for 
further information.

## Example

```julia-repl
julia> using MolSimToolkit

julia> x = BlockAverages.test_data(10^6); # example data generator

julia> b = block_average(x, lags=0:100:10^5)
-------------------------------------------------------------------
BlockAverageData{Float64}
-------------------------------------------------------------------
Estimated value (mean by default) = -0.13673023261452855
Length of data series: 1000000

Block sizes: [1, 2, ..., 500000, 1000000]

Maximum standard error (error, block size): (0.23264202379194165, 500000)

Deviations in last 3 blocks:
         percentual: [-349.5348293165444, -170.1467329817311, -0.0]  
           absolute: [0.47791978519330647, 0.23264202379194168, 0.0]  

Autocorrelation is first zero with lag: 16400
Characteristic time of autocorrelation decay: 
        as fraction of series length: 0.0037856443348888848
                            absolute: 3785.6443348888847
-------------------------------------------------------------------

julia> using Plots

julia> plot(b) # creates a plot with the results

```
"""
function block_average(
    x_input::AbstractVector{T};
    by=mean,
    min_block_size::Integer=1,
    max_block_size::Integer=length(x_input),
    lags::Union{Nothing,AbstractVector{<:Integer}}=nothing
) where {T<:Real}

    x = adjust_xinput(x_input, max_block_size)

    n = length(x)
    xmean = by(x)
    xmean_maxerr = Float64[]
    xmean_stderr = Float64[]
    blocksize = Int[]

    for block_size in min_block_size:max_block_size
        nblocks = n ÷ block_size
        remaining = n % block_size
        if remaining != 0
            continue
        end

        # Add new point to vectors
        push!(blocksize, block_size)
        push!(xmean_maxerr, zero(T))
        push!(xmean_stderr, zero(T))

        # Compute the property in each block, and keep the maximum error
        diff_max = -Inf
        for i in 1:nblocks
            xblock = @view x[brange(i, block_size)]
            this_block_mean = by(xblock)
            diff = abs(this_block_mean - xmean)
            if diff > diff_max
                diff_max = diff
                xmean_maxerr[end] = this_block_mean
            end
            # Compute the standard deviation of the property estimate (σ²/√N)
            xmean_stderr[end] += (this_block_mean - xmean)^2
        end
        if nblocks > 1
            xmean_stderr[end] = sqrt(xmean_stderr[end] / (nblocks - 1))
            # We want the standard error, so
            xmean_stderr[end] /= sqrt(nblocks)
        end
    end

    # Compute auto-correlation function of the data
    if isnothing(lags)
        lags = 0:round(Int, min(size(x, 1) - 1, 10 * log10(size(x, 1))))
        auto_cor = autocor(x)
    else
        auto_cor = autocor(x, lags)
    end

    tau = fitexp(lags, auto_cor, c=0.0, u=upper(a=1.1), l=lower(a=0.9)).b

    return BlockAverageData{T}(
        x,
        xmean,
        blocksize,
        xmean_maxerr,
        xmean_stderr,
        lags,
        auto_cor,
        tau
    )

end

"""

$(TYPEDEF)

$(TYPEDFIELDS)

"""
struct BlockDistribution{N}
    mean::Float64
    std_of_the_mean::Float64
    block_mean::Vector{Float64}
    std_err_of_the_mean::Float64
end
function Base.show(io::IO, ::MIME"text/plain", m::BlockDistribution)
    print(
        io,
        """
        -------------------------------------------------------------------
        $(typeof(m))
        -------------------------------------------------------------------
        Number of blocks: $(length(m.block_mean))
        Estimated mean: = $(m.mean)
        Standard error of the mean: $(m.std_err_of_the_mean)
        Standard deviation of the mean: $(m.std_of_the_mean)
        > block_mean contains the mean computed for each block.
        -------------------------------------------------------------------
        """
    )
end

"""
    block_distribution(x_input::AbstractVector; block_size::Integer) = 
        block_distribution(mean, x_input, block_size::Integer)
    block_distribution(by::Function, x_input::AbstractVector, block_size::Integer)

Given the data and the block size, computes the distribution of estimates of the 
properties for each block. Returns a `BlockDistribution{NBLOCKS}` object. The
block size must be an integer.

# Example

```julia-repl
julia> using MolSimToolkit

julia> x = BlockAverages.test_data(10^7);

julia> d = block_distribution(x; block_size = 25*10^3)
-------------------------------------------------------------------
BlockDistribution{400}
-------------------------------------------------------------------
Number of blocks: 400
Estimated value: = 0.025151622077551537
Standard error of the mean: 0.05596145099711976
Standard deviation of the mean: 1.119229019942395
> block_mean contains the mean computed for each block.
-------------------------------------------------------------------
```

The distribution is stored in the `d.block_mean` vector, and can be plotted with:
```julia-repl
julia> using Plots

julia> histogram(d)
```

"""
function block_distribution(by::Function, x_input::AbstractVector, block_size::Integer)
    if block_size < 1
        throw(ArgumentError((
            "block_size not properly defined. Use, for example: block_distribution(x; block_size=10^5)"
        )))
    end
    x = adjust_xinput(x_input, block_size, "block_size")
    n = length(x)
    nblocks = n ÷ block_size
    block_mean = fill(0.0, nblocks)
    # Compute the property in each block
    for i in 1:nblocks
        xblock = @view x[brange(i, block_size)]
        block_mean[i] = by(xblock)
    end
    return BlockDistribution{nblocks}(mean(x), std(block_mean), block_mean, std(block_mean) / sqrt(nblocks))
end
block_distribution(x_input::AbstractVector; block_size::Integer=0) = block_distribution(mean, x_input, block_size)

# Range of indices for block i of size block_size
function brange(i, block_size)
    ifirst = (i - 1) * block_size + 1
    ilast = ifirst + block_size - 1
    return ifirst:ilast
end

# Generate correlated data to test
function test_data(n)
    temperature = 1.0
    x = Vector{Float64}(undef, n)
    x[1] = 0.0
    u = 0.0
    i = 1
    while i < n
        x_trial = x[i] - 0.1 + 0.2 * rand()
        u_trial = x_trial^2
        if u_trial < u || exp((u - u_trial) / temperature) > 0.5
            i += 1
            x[i] = x_trial
            u = u_trial
        end
    end
    x
end

end # module BlockAverage