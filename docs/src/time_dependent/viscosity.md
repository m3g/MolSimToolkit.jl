```@meta
CollapsedDocStrings = true
```
# Viscosity (Green-Kubo)

The shear viscosity of a fluid can be computed from the autocorrelation of
the off-diagonal components of its pressure tensor, using the Green-Kubo
relation:

```math
\eta = \frac{V}{k_B T}\int_0^\infty \left< P_{\alpha\beta}(0) P_{\alpha\beta}(t) \right> dt
```

The `read_namd_pressure_tensor` function reads the pressure tensors printed in
the log file of a NAMD simulation, and `green_kubo_viscosity` computes the
autocorrelation function and its running integral.

## Running the simulation

NAMD prints the pressure tensor in the log file only if the `outputPressure`
option is set. The autocorrelation of the pressure tensor decays within a few
picoseconds in liquid water, thus the tensor must be printed at short
intervals, for instance at every step:

```
outputPressure     1
```

The Green-Kubo integral converges slowly: simulations of several nanoseconds
are required for precise estimates. Stochastic thermostats (as Langevin
dynamics) affect the dynamics of the system, and thus the viscosity, so the
production simulation is preferably run in the NVE ensemble, after
equilibration at the desired temperature and pressure.

```@docs
read_namd_pressure_tensor
green_kubo_viscosity
MolSimToolkit.PressureTensor
MolSimToolkit.GreenKuboViscosity
```

## Example: viscosity of TIP3P water

Here we use a short NAMD simulation of 900 TIP3P water molecules, 
included in the `MolSimToolkit.Testing` test data (2000 steps of 2 fs, with 
the pressure tensor printed at every step). The simulation is much too short 
for a converged estimate of the viscosity, but illustrates the use of the 
functions:

```@example viscosity
using MolSimToolkit, MolSimToolkit.Testing, Plots

p = read_namd_pressure_tensor(Testing.namd_pressure_log)
```

```@example viscosity
gk = green_kubo_viscosity(p; tmax=1.0)
```

```@example viscosity
plot(MolSimStyle,
    gk.time, gk.viscosity,
    xlabel="time / ps", ylabel="η / mPa⋅s",
    linewidth=2, label=nothing,
)
```

The viscosity is estimated from the plateau of the running integral, for example
by averaging it over a range of time lags in which it is approximately constant:

```@example viscosity
using Statistics: mean
mean(gk.viscosity[gk.time .> 0.5])
```
