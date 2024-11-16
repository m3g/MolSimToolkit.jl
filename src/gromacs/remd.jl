export remd_data
export remd_replica_path

"""
    GromacsRMEDlog

Structure to store the log from a REMD simulation performed with Gromacs. 
This structure contains three fields:

- `steps`: Vector of steps at which the exchange was performed.
- `exchange_matrix`: Matrix of exchanges performed. 
  Each row corresponds to a step and each column to a replica. 
- `probability_matrix`: Matrix of probabilities of finding each replica at level of 
  perturbation. Each column corresponds to a replica and each row to a level of
  perturbation.

"""
struct GromacsREMDlog
    steps::Vector{Int}
    exchange_matrix::Matrix{Int}
    probability_matrix::Matrix{Float64}
end
const gmx2019_4_log = "$(@__DIR__)" * "/../../test/data/remd/gmx-2019.4.log"
const gmx5_0_4_log = "$(@__DIR__)" * "/../../test/data/remd/gmx-5.0.4.log"

"""
    remd_data(log::String)

Function to read the log file from a (H)REMD simulation performed with Gromacs.

Returns a `GromacsREMDlog` structure, containing the steps at which the exchange
was tried, the exchange matrix and the probability matrix. The exchange matrix
contains the replica number at each level of perturbation for each step. The
probability matrix contains the probability of finding each replica at each level.

Tested with log files of Gromacs versions: 
    - 2019.4
    - 5.0.4

# Example

First obtaina the REMD data from the log file:

```julia-repl
julia> using MolSimToolkit

julia> data = remd_data(MolSimToolkit.gmx2019_9_log)
```

Then plot the exchange matrix, which will provide a visual inspection of the exchange process:

```julia-repl
julia> using Plots

julia> heatmap(data) 
```
"""
function remd_data(log::String)
    if !isfile(log)
        throw(ArgumentError("File $log does not exist."))
    end
    step = 0
    nreplicas = 0
    swaps = Int[]
    steps = Int[]
    exchanges = Vector{Int}[]
    open(log, "r") do io
        for line in eachline(io)
            length(line) < 7 && continue
            if occursin("Replica exchange at step", line)
                data = split(line)
                step = parse(Int, data[5])
            end
            if line[1:7] == "Repl ex"
                data = split(line[8:end])
                if length(exchanges) == 0
                    nreplicas = parse(Int, data[end]) + 1
                    push!(exchanges, [i for i in 0:nreplicas-1])
                    push!(steps, 0)
                    swaps = zeros(Int, nreplicas)
                end
                swaps .= last(exchanges)
                iswap = 1
                for i in eachindex(data)
                    if data[i] == "x"
                        swaps[iswap-1], swaps[iswap] = swaps[iswap], swaps[iswap-1]
                    else
                        iswap += 1
                    end
                end
                push!(steps, step)
                push!(exchanges, copy(swaps))
            end
        end
    end
    exchange_matrix = zeros(Int, length(exchanges), nreplicas)
    for iframe in eachindex(exchanges)
        exchange_matrix[iframe, :] .= exchanges[iframe]
    end
    probability_matrix = reduce(hcat,
        [[count(==(i), exchange_matrix[:, j]) for i in 0:nreplicas-1] for j in 1:nreplicas]
    ) ./ length(steps)
    return GromacsREMDlog(steps, exchange_matrix, probability_matrix)
end

"""
    remd_replica_path(data::GromacsREMDlog, replica::Integer; stride::Integer = 1)

Function to obtain the path of a replica in the exchange matrix.

"""
remd_replica_path(data::GromacsREMDlog, replica::Integer; stride::Integer=1) =
    [findfirst(==(replica), data.exchange_matrix[i, :]) - 1 for i in 1:stride:size(data.exchange_matrix, 1)];

@testitem "REMD" begin
    using MolSimToolkit
    using MolSimToolkit: gmx2019_4_log, gmx5_0_4_log
    data = remd_data(gmx2019_4_log)
    @test size(data.exchange_matrix) == (251, 10)
    @test size(data.probability_matrix) == (10, 10)
    @test all(≈(1), sum(data.probability_matrix, dims=1))
    @test all(≈(1), sum(data.probability_matrix, dims=2))
    @test length(data.steps) == 251
    path1 = remd_replica_path(data, 1)
    @test all(d -> d in (-1, 0, 1), path1[i+1] - path1[i] for i in 1:length(path1)-1)
    @test length(path1) == 251
    path1 = remd_replica_path(data, 1; stride=10)
    @test length(path1) == 26
    data = remd_data(gmx5_0_4_log)
    @test size(data.exchange_matrix) == (201, 16)
    @test size(data.probability_matrix) == (16, 16)
    @test all(≈(1), sum(data.probability_matrix, dims=1))
    @test all(≈(1), sum(data.probability_matrix, dims=2))
    @test length(data.steps) == 201
end