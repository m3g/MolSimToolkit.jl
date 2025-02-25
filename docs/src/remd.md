```@meta
CollapsedDocStrings = true
```

# Replica exchange analysis

The `remd_data` function reads the output of a Gromacs-generated 
replica-exchange simulation file, and provides some tools for visualization
of the quality of the exchange process.

!!! compat
    This function was introduced in `MolSimToolkit` version 1.1.0. 

    The function was tested to read log files produced by Gromacs versions: 
      - 2019.4
      - 5.0.4
    Compatibility with other versions is not guaranteed (issue reporting and contributions are welcome). 

    The `heatmap` and the support for the `stride` argument in `remd_replica_path` where introduced in version 1.6.0

## Reading REMD data

First, read the data from the Gromacs simulation log file:

```julia-repl
julia> using MolSimToolkit

julia> data = remd_data("gromacs.log")
```
where "gromacs.log" would be a `log` file produced by Gromacs.

This will result in a data structure with three fields:

- `steps`: Vector of steps at which the exchange was performed.
- `exchange_matrix`: Matrix of exchanges performed. 
  Each row corresponds to a step and each column to a replica. 
- `probability_matrix`: Matrix of probabilities of finding each replica at level of 
  perturbation. Each column corresponds to a replica and each row to a level of
  perturbation.

## Probability heatmap

One way to visualize the exchange it to produce a heatmap of expected probabilities. This
can be done with the auxiliary `heatmap` function that is provided for the output
of `remd_data`: 

```@example
using MolSimToolkit
using Plots
data = remd_data(MolSimToolkit.gmx5_0_4_log)
heatmap(data)
```

The number of replicas here is 16 (0-15), thus the expected ideal probability of finding each replica
in each level is $1/16$. The probabilities are divided by $1/16$, such that $1.0$ implies 
an optimal exchange at that replica and level. 

To produce a similar heatmap, but with the absolute (not normalized) probabilities of 
observing each replica at each level, use `heatmap(data; probability_type=:absolute)`. 

In the above example, replica 1 got trapped in the first 5 levels, 
while a block is noticeable for replicas 10-12, which display high 
probabilities to be found at levels 13 to 15. Thus, this is an example
of a inadequate exchange between the levels.

!!! compat
    The `probability_type` option of `heatmap` was introduced in version 1.7.0.

## Replica path

A heatmap as the one above suggests checking the path of the replicas along the exchange. 
This can be obtained with the `remd_replica_path` function. For example, to obtain the path
of the replicas of number 0 and 1. Replica 0 appeares to have visited reasonably well 
all levels from 0 to 12, and replica 1 appears to be trapped in leves 13 to 15.

```@example
using MolSimToolkit
using Plots
data = remd_data(MolSimToolkit.gmx5_0_4_log)

# Obtain the paths
path0 = remd_replica_path(data, 0; stride = 1)
path10 = remd_replica_path(data, 10; stride = 1)

# Plot the path
plt = plot(MolSimStyle, xlabel="step", ylabel="replica level")
plot!(plt, path0, label="replica 0", linewidth=2)
plot!(plt, path10, label="replica 10", linewidth=2)
```

The plot confirms that the replica starting at position 0 sampled preferentially the 
first 5 states, while replica 10 was trapped temporarily in the high level states,
in this sample run.

## Probability data

An alternative visualization of the exchange process is
given by the probability matrix:

```@example
using MolSimToolkit
using Plots
data = remd_data(MolSimToolkit.gmx5_0_4_log)
scatter(MolSimStyle,
    data.probability_matrix,
    labels= Ref("Replica ") .* string.((0:15)'),
    linewidth=2,
    ylims=(0,0.12), xlims=(0.7, 16.3),
    xlabel="Level", xticks=(1:16, 0:15),
    ylabel="Probability",
    alpha=0.5,
    margin=0.5Plots.Measures.cm,
    legend=:outertopright,
)
```

Ideally, the probability of each replica populaing each level should be the inverse of the number of replicas (here $$1/16 == 0.0625$$). In this case, the simulation does not provide a proper sampling
of exchanges, as it is a short extract of a longer simulation. 

## Reference functions

```@autodocs
Modules = [MolSimToolkit]
Pages = ["gromacs/remd.jl"]
Order = [:function, :type]
```




