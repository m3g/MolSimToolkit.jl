```@meta
CollapsedDocStrings = true
```
# Block averages

Performs an analysis of the convergence of some property (usually the mean) in a time-series. 

Computes the block average of time-dependent data, to compute the standard error of the mean and, to detect sampling problems. A didactical explanation of block averaging is available [here](http://sachinashanbhag.blogspot.com/2013/08/block-averaging-estimating-uncertainty.html).  

The package also outputs the autocorrelation function of the property, and the characteristic time of the correlation
decay. 

## Data that is not time-correlated

### Compute the average of a random variable `x`:

```julia-repl
julia> using MolSimToolkit 

julia> x = rand(10_000);

julia> b = block_average(x)
-------------------------------------------------------------------
BlockAverageData{Float64}
-------------------------------------------------------------------
Estimated value (mean by default) = 0.4977014924716461
Length of data series: 10000

Block size ranges: (1, 10000)

Maximum standard error (error, block size): (0.005790454921861948, 5000)

Deviations in last 3 blocks:
         percentual: [2.8893559885598195, -1.1634393325014705, 0.0]  
           absolute: [0.014380367877881106, -0.005790454921861976, 0.0]  

Autocorrelation is first zero with lag: 2
Characteristic time of autocorrelation decay: 
        as fraction of series length: 2.0182708552030272e-5
                            absolute: 0.2018270855203027
-------------------------------------------------------------------

julia> using Plots

julia> plot(b, title="Uncorrelated data")
```

Results in:

![random.png](./images/block_averages/random.png)

Thus, the worst block estimate converges very rapidly to the true mean, the standard error of the mean is very small, and the autocorrelation decays very quickly. Since the data is not correlated along the series, the characteristic time may not be meaningful. 

## Data that is time-correlated

### Poorly-sampled data

The data above is not correlated in the input `x` vector. If the data is correlated, one can observe that in the dependence of the estimates of the average and error from the data. One can generate a test data (sort of a monte-carlo simulation of a particle in an harmonic well) using:

```julia-repl
julia> x = BlockAverages.test_data(1_000);
```
Which in this run produced:

```julia-repl
julia> plot(x; xlabel="time", ylabel="value", label=nothing)
```

![bad_sampling.png](./images/block_averages/bad_sampling.png)

The error of the estimate of the mean is, now, dependent on the block size, and we cannot see any convergence of the error, indicating that the sampling is not enough to obtain a reliable estimate of the mean value of `x`:  

```julia-repl
julia> b = block_average(x, lags=1:500)
-------------------------------------------------------------------
BlockAverageData{Float64}
-------------------------------------------------------------------
Estimated value (mean by default) = -0.5616539552467762
Length of data series: 1000

Block size ranges: (1, 1000)

Maximum standard error (error, block size): (0.24081057091463817, 500)

Deviations in last 3 blocks:
         percentual: [80.94415878065803, 42.87525595877488, -0.0]  
           absolute: [-0.4546260693327965, -0.2408105709146382, 0.0]  

Autocorrelation is first zero with lag: 194
Characteristic time of autocorrelation decay: 
        as fraction of series length: 0.06981859583876429
                            absolute: 69.8185958387643
-------------------------------------------------------------------
```

Note that we have adjusted the range of `lags` of the autocorrelation function in this case.

Several characteristics of the output indicate the poor convergence of the series: 1) The mean should be `0.` for this property; 2) the maximum standard error occurs with a block size which is half the length of the series (there is no plateau); 3) the standard error of the mean is of the order of the mean value; 4) The autocorrelation is first zero at ~20% of the length of the data set. 

The corresponding plot is obtained with:
```julia-repl
julia> using Plots

julia> plot(b, title="Bad sampling")
```

![bad_sampling_result.png](./images/block_averages/bad_sampling_result.png)

### Properly sampled data

If we increase the sampling by generating longer simulation:
```julia-repl
julia> x = BlockAverages.test_data(10^6);
```

The obtained set is now much better sampled,

![good_sampling.png](./images/block_averages/good_sampling.png)

The convergence analysis of the series produces:
```julia-repl
julia> b = block_average(x, lags=1:100:10^5, max_block_size=10^5)
-------------------------------------------------------------------
BlockAverageData{Float64}
-------------------------------------------------------------------
Estimated value (mean by default) = -0.05498853009368246
Length of data series: 1000000

Block sizes: [1, 2, ..., 62500, 100000]

Maximum standard error (error, block size): (0.18706372724807982, 31250)

Deviations in last 3 blocks:
         percentual: [-2805.4693758538297, -2600.14341058853, -1393.4253407524507]  
           absolute: [1.5426863720104287, 1.4297806418103753, 0.7662241128326587]  

Autocorrelation is first zero with lag: 14701
Characteristic time of autocorrelation decay: 
        as fraction of series length: 0.0036203287638847167
                            absolute: 3620.328763884717
-------------------------------------------------------------------
```

Note that the average value of `x` here is closer to zero, and that the maximum standard error of the mean is consistent the true value being around zero. The correlation decays fast relative to the length of the data series.

The corresponding plots are:

```julia-repl
julia> using Plots

julia> plot(b, title="Good sampling")
```

![good_sampling_result.png](./images/block_averages/good_sampling_result.png)

The plateau of standard errors (second plot) in intermediate values of block sizes is typical of a properly sampled data set, and can be used as an the uncertainty in the property estimate. 

For example, for an ever better sampled data, there is a very clear plateau of standard errors, which are smaller than those of the above example:

```julia-repl
julia> x = BlockAverages.test_data(10^7)

julia> b = block_average(x, max_block_size=10^5, lags=1:100:10^5);

julia> plot(b)
```

![best_sampling.png](./images/block_averages/best_sampling.png)

(here we have computed the statistics only up to blocks of size `10^5`)

## Visualizing the distribution of the mean

Once the overall correlation is understood from the time-series block analysis, one can 
visualize the distribution of the computed value (the mean in general) for a given
block size. A relatively good fit to a gaussian distribution is expected. For instance,
let us choose a block size of `25_000`, from a set similar to the one above:

```julia-repl
julia> using MolSimToolkit, Plots

julia> x = BlockAverages.test_data(10^7) 

julia> d = block_distribution(x; block_size = 25_000)
-------------------------------------------------------------------
MeanDistribution{400}
-------------------------------------------------------------------
Number of blocks: 400
Estimated value: = 0.06462623778329132
Standard error of the mean: 0.05641314321229929
Standard deviation of the mean: 1.1282628642459858
> block_mean contains the mean computed for each block.
-------------------------------------------------------------------

julia> histogram(d)
```

The last command will produce a plot similar to the following, in which the histogram
of values is shown side by side with the gaussian function that corresponds to the 
observed standard deviation and mean.

![best_sampling.png](./images/block_averages/mean_distribution.svg)

The optimal block size should be that for which the distribution is closer to a gaussian.

## Reference functions

```@autodocs
Modules = [MolSimToolkit.BlockAverages]
Pages = ["BlockAverages.jl"]
Order = [:function, :type]
```

```@docs
plot(::BlockAverageData)
histogram(::BlockDistribution)
```

