```@meta
CollapsedDocStrings = true
```
# Block averages

Performs an analysis of the convergence of some property (usually the mean) in a time-series. 

Computes the block average of time-dependent data, to compute the standard error of the mean and, to detect sampling problems. A didactical explanation of block averaging is available [here](http://sachinashanbhag.blogspot.com/2013/08/block-averaging-estimating-uncertainty.html).  

The package also outputs the autocorrelation function of the property, and the characteristic time of the correlation decay. 

## Effective samples

The number of effective samples if obtained by computing the integrated autocorrelation time $$\tau_{INT}$$,
```math
\tau_{INT} = 1 + 2\sum_{i=1}^M c(t)
```
where $$c(t)$$ is the autocorrelation function of the data and $$M$$ is the time lag (in data points) where $$c(t) < 1.96/\sqrt{N}$$ for the first time, which is the 95% confidence interval for the autocorrelation. The effective number of samples is then
```math
N_{\textrm{eff}} = \frac{N}{\tau_{INT}}
```
This allows the estimation of the effective standard error of the mean, as $$SD(data)/\sqrt{N_{\textrm{eff}}}$$.

## Data that is not time-correlated

### Compute the average of a random variable `x`:

```@example block_averages
using MolSimToolkit, Plots 
x = rand(10_000)
b = block_average(x)
```

```@example block_averages
plot(b, title="Uncorrelated data")
```

Thus, the worst block estimate converges very rapidly to the true mean, the standard error of the mean is very small, and the autocorrelation decays very quickly. Since the data is not correlated along the series, the characteristic time may not be meaningful. 

## Data that is time-correlated

### Poorly-sampled data

The data above is not correlated in the input `x` vector. If the data is correlated, one can observe that in the dependence of the estimates of the average and error from the data. One can generate a test data (sort of a monte-carlo simulation of a particle in an harmonic well) using:

```@example block_averages
x = BlockAverages.test_data(10_000; variance=0.01)
b = block_average(x)
```

The error of the estimate of the mean is, now, dependent on the block size, and we cannot see any convergence of the error, indicating that the sampling is not enough to obtain a reliable estimate of the mean value of `x`.

```@example block_averages
plot(b, title="Poorly sampled data")
```

Several characteristics of the output indicate the poor convergence of the series: 1) the maximum standard error occurs with a block size which is half the length of the series (there is no plateau); 2) the standard error of the mean is of the order of the mean value; 3) The autocorrelation is first zero at ~20% of the length of the data set. 4) Finally, note that the number of effective samples is very small.

### Properly sampled data

If we increase the sampling by generating longer simulation with a greater variance between sucessive data points, the convergence analysis produces:
```@example block_averages
x = BlockAverages.test_data(10^5; variance=0.1);
b = block_average(x)
```

Note that the average value of `x` here is closer to zero, and that the maximum standard error of the mean is consistent the true value being around zero. The correlation decays fast relative to the length of the data series.

The corresponding plots are:

```@example block_averages
plot(b, title="Good sampling")
```

The plateau of standard errors (second plot) in intermediate values of block sizes is typical of a properly sampled data set, and can be used as an the uncertainty in the property estimate. 

For example, for an ever better sampled data, there is a very clear plateau of standard errors, which are smaller than those of the above example:

## Visualizing the distribution of the mean

Once the overall correlation is understood from the time-series block analysis, one can 
visualize the distribution of the computed value (the mean in general) for a given
block size. A relatively good fit to a gaussian distribution is expected. For instance,
let us choose a block size of `25_000`, from a set similar to the one above:

```@example block_averages
x = BlockAverages.test_data(10^5) 
d = block_distribution(x; block_size = 2_500)
```

The following command will produce a plot similar to the following, in which the histogram
of values is shown side by side with the gaussian function that corresponds to the 
observed standard deviation and mean.

```@example block_averages
histogram(d)
```

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

