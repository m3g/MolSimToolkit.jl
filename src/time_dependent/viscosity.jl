"""
    PressureTensor

Pressure tensors of a simulation, as read by [`read_namd_pressure_tensor`](@ref).

# Fields

- `steps::Vector{Int}`: The simulation steps at which the pressure tensor was printed.
- `tensors::Vector{SMatrix{3,3,Float64,9}}`: The pressure tensor at each of these steps, in bar.
- `timestep::Float64`: The integration time step of the simulation, in fs.
- `volume::Float64`: The average volume of the system, in Å³.
- `temperature::Float64`: The average temperature of the system, in K.

!!! compat
    This structure was added in version 2.7.0 of MolSimToolkit.

"""
struct PressureTensor
    steps::Vector{Int}
    tensors::Vector{SMatrix{3,3,Float64,9}}
    timestep::Float64
    volume::Float64
    temperature::Float64
end

Base.length(p::PressureTensor) = length(p.tensors)

function Base.show(io::IO, ::MIME"text/plain", p::PressureTensor)
    print(io, chomp(
        """
        -------------------------------------------------------------------
        PressureTensor:
        -------------------------------------------------------------------
        Number of pressure tensors: $(length(p))
        Steps: $(isempty(p.steps) ? "none" : "$(first(p.steps)) to $(last(p.steps))")
        Time step: $(p.timestep) fs
        Average volume: $(@sprintf("%.4f", p.volume)) Å³
        Average temperature: $(@sprintf("%.4f", p.temperature)) K
        -------------------------------------------------------------------
        """
    ))
end

"""
    read_namd_pressure_tensor(logfile::AbstractString; group_pressure::Bool=false)

Reads the pressure tensors printed in the log file of a NAMD simulation.
Returns a [`PressureTensor`](@ref) object, which contains the steps,
the tensors (in bar), the time step (in fs), and the average volume (in Å³) and
temperature (in K) of the simulation, as reported in the `ENERGY:` lines of the log.

NAMD prints the pressure tensor only if the `outputPressure` option is set in the
input file, for instance with `outputPressure 1` to print it at every step. The tensors
are then printed in lines of the form:

```
PRESSURE: step Pxx Pxy Pxz Pyx Pyy Pyz Pzx Pzy Pzz
GPRESSURE: step Pxx Pxy Pxz Pyx Pyy Pyz Pzx Pzy Pzz
```

where `PRESSURE` is the atomic pressure tensor, and `GPRESSURE` is the group (molecular)
pressure tensor, computed from the centers of mass of hydrogen groups.

To compute the viscosity of the system from these tensors with the Green-Kubo relation,
use [`green_kubo_viscosity`](@ref). For that, the tensors must be printed at short intervals
(typically every one to a few steps) relative to the decay time of their autocorrelations,
of the order of a picosecond for liquid water.

If a step is printed more than once (as happens when the input file contains multiple
`run` commands), only its first occurrence is kept.

# Arguments

- `logfile::AbstractString`: The path to the NAMD log file.

# Optional keyword arguments

- `group_pressure::Bool`: If `true`, the group pressure tensors (`GPRESSURE` lines) are read
  instead of the atomic ones (`PRESSURE` lines). Defaults to `false`.

# Example

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2***"
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> p = read_namd_pressure_tensor(Testing.namd_pressure_log)
-------------------------------------------------------------------
PressureTensor:
-------------------------------------------------------------------
Number of pressure tensors: 2001
Steps: 0 to 2000
Time step: 2.0 fs
Average volume: 26721.5138 Å³
Average temperature: 296.3006 K
-------------------------------------------------------------------

julia> p.tensors[1]
3×3 StaticArraysCore.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 1693.76   -483.637  -258.105
 -483.64    464.369  -144.419
 -258.106  -144.43    195.043
```

!!! compat
    This function was added in version 2.7.0 of MolSimToolkit.

"""
function read_namd_pressure_tensor(logfile::AbstractString; group_pressure::Bool=false)
    key = group_pressure ? "GPRESSURE:" : "PRESSURE:"
    steps = Int[]
    tensors = SMatrix{3,3,Float64,9}[]
    timestep = NaN
    itemp = 0
    ivolume = 0
    temperature = 0.0
    volume = 0.0
    nenergy = 0
    for line in eachline(logfile)
        if startswith(line, key)
            data = split(line)
            if length(data) != 11
                throw(ArgumentError("""\n
                    Could not read pressure tensor from line:

                    $line

                """))
            end
            step = parse(Int, data[2])
            # Repeated steps occur at the beginning of every `run` command
            !isempty(steps) && step <= last(steps) && continue
            push!(steps, step)
            # NAMD prints the tensor row by row: xx xy xz yx yy yz zx zy zz
            push!(tensors, SMatrix{3,3,Float64,9}(ntuple(i -> parse(Float64, data[2+i]), 9))')
        elseif startswith(line, "Info: TIMESTEP")
            timestep = parse(Float64, last(split(line)))
        elseif startswith(line, "ETITLE:")
            titles = split(line)
            itemp = something(findfirst(==("TEMP"), titles), 0)
            ivolume = something(findfirst(==("VOLUME"), titles), 0)
        elseif startswith(line, "ENERGY:") && itemp > 0 && ivolume > 0
            data = split(line)
            temperature += parse(Float64, data[itemp])
            volume += parse(Float64, data[ivolume])
            nenergy += 1
        end
    end
    if isempty(tensors)
        throw(ArgumentError("""\n
            No $key lines found in $logfile.
            Set, for example, `outputPressure 1` in the NAMD input file to print the pressure tensor.

        """))
    end
    if nenergy == 0
        throw(ArgumentError("""\n
            Could not read the temperature and volume from the ENERGY: lines of $logfile.

        """))
    end
    return PressureTensor(steps, tensors, timestep, volume / nenergy, temperature / nenergy)
end

@testitem "read_namd_pressure_tensor" begin
    using MolSimToolkit, MolSimToolkit.Testing
    using StaticArrays: SMatrix
    p = read_namd_pressure_tensor(Testing.namd_pressure_log)
    @test length(p) == 2001
    @test p.steps == 0:2000
    @test p.timestep == 2.0
    @test p.volume ≈ 26721.5138
    @test p.temperature ≈ 296.30055714285714
    @test p.tensors[1] ≈ SMatrix{3,3}(1693.76, -483.64, -258.106, -483.637, 464.369, -144.43, -258.105, -144.419, 195.043)
    @test p.tensors[1][1, 2] == -483.637
    @test p.tensors[1][2, 1] == -483.64
    pg = read_namd_pressure_tensor(Testing.namd_pressure_log; group_pressure=true)
    @test length(pg) == 2001
    @test pg.tensors[1][1, 2] == -693.43
    @test pg.tensors[1][2, 1] == -329.522
    # Repeated steps are skipped
    tmp = tempname()
    lines = readlines(Testing.namd_pressure_log)
    write(tmp, join(vcat(lines, filter(startswith("PRESSURE:"), lines)), "\n"))
    @test length(read_namd_pressure_tensor(tmp)) == 2001
    # Errors
    write(tmp, join(filter(!startswith("PRESSURE:"), lines), "\n"))
    @test_throws ArgumentError read_namd_pressure_tensor(tmp)
    write(tmp, join(filter(!startswith("ENERGY:"), lines), "\n"))
    @test_throws ArgumentError read_namd_pressure_tensor(tmp)
    rm(tmp)
end

"""
    GreenKuboViscosity

Result of the Green-Kubo computation of the shear viscosity, as returned by
[`green_kubo_viscosity`](@ref).

# Fields

- `time::Vector{Float64}`: The time lags, in ps.
- `acf::Vector{Float64}`: The autocorrelation function of the (shear) pressure tensor
  components, at each time lag, in bar².
- `viscosity::Vector{Float64}`: The running integral of the Green-Kubo relation,
  that is, the viscosity estimated by integrating the autocorrelation function up to each
  time lag, in mPa·s (cP).

!!! compat
    This structure was added in version 2.7.0 of MolSimToolkit.

"""
struct GreenKuboViscosity
    time::Vector{Float64}
    acf::Vector{Float64}
    viscosity::Vector{Float64}
end

function Base.show(io::IO, ::MIME"text/plain", gk::GreenKuboViscosity)
    print(io, chomp(
        """
        -------------------------------------------------------------------
        GreenKuboViscosity:
        -------------------------------------------------------------------
        Number of time lags: $(length(gk.time))
        Maximum time lag: $(@sprintf("%.4f", last(gk.time))) ps
        Viscosity at maximum time lag: $(@sprintf("%.4f", last(gk.viscosity))) mPa·s
        -------------------------------------------------------------------
        """
    ))
end

"""
    green_kubo_viscosity(
        p::PressureTensor;
        temperature::Real = p.temperature,
        volume::Real = p.volume,
        tmax::Real = 10.0,
        components::Symbol = :traceless,
        parallel::Bool = true,
    )

Computes the shear viscosity of a system from the autocorrelation of its pressure tensor,
using the Green-Kubo relation:

```math
\\eta = \\frac{V}{k_B T}\\int_0^\\infty \\left< P_{\\alpha\\beta}(0) P_{\\alpha\\beta}(t) \\right> dt
```

where ``P_{\\alpha\\beta}`` are the off-diagonal (shear) components of the pressure tensor,
``V`` is the volume of the system and ``T`` its temperature. The pressure tensors are typically
obtained from [`read_namd_pressure_tensor`](@ref).

Returns a [`GreenKuboViscosity`](@ref) object, with the time lags (in ps), the autocorrelation
function (in bar²) and the running integral of the Green-Kubo relation (the viscosity, in mPa·s),
as a function of the upper limit of the integral. The viscosity of the system is estimated from the
plateau of this running integral, which must be identified by inspecting the curve: at long times
the integral becomes increasingly noisy, since the autocorrelation function is then dominated by
statistical noise.

Two estimates of the autocorrelation function can be used, controlled by `components`:

- `:offdiagonal`: the average of the autocorrelations of the three independent off-diagonal
  components (``xy``, ``xz`` and ``yz``) of the symmetrized pressure tensor, ``(P + P^T)/2``.
- `:traceless` (default): the autocorrelations of all components of the symmetrized, traceless,
  pressure tensor, ``P' = (P + P^T)/2 - \\mathrm{tr}(P)/3 \\, I``, as
  ``\\sum_{\\alpha\\beta} \\left< P'_{\\alpha\\beta}(0) P'_{\\alpha\\beta}(t)\\right>/10``.
  For an isotropic fluid this gives the same result as the off-diagonal components,
  with better statistics (Daivis and Evans, J. Chem. Phys. 100, 541 (1994)).

The time average of each component is subtracted before computing the autocorrelations,
and the autocorrelation is averaged over all time origins. The integral is computed with
the trapezoidal rule.

!!! note "Sampling"
    The Green-Kubo integral converges slowly with the length of the simulation. For liquid
    water, simulations of several nanoseconds, with the pressure tensor printed every few
    femtoseconds, are required for estimates with a precision of a few percent. Viscosities
    obtained from simulations with stochastic thermostats (e.g. Langevin dynamics) are
    affected by the friction of the thermostat; NVE simulations, or simulations with weakly
    coupled thermostats, are preferred.

# Arguments

- `p::PressureTensor`: The pressure tensors, as read by [`read_namd_pressure_tensor`](@ref).
  The tensors must be printed at regular step intervals.

# Optional keyword arguments

- `temperature::Real`: The temperature of the system, in K. Defaults to the average
  temperature read from the log file.
- `volume::Real`: The volume of the system, in Å³. Defaults to the average volume read
  from the log file.
- `tmax::Real`: The maximum time lag, in ps, up to which the autocorrelation function
  and its integral are computed. Defaults to `10.0` ps, and is limited to half of the
  length of the simulation.
- `components::Symbol`: Either `:traceless` (default) or `:offdiagonal`; see above.
- `parallel::Bool`: Defines if the computation of the autocorrelation function is run in
  parallel. Defaults to `true`. Requires starting Julia with multi-threading.

# Example

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2***"
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> p = read_namd_pressure_tensor(Testing.namd_pressure_log);

julia> gk = green_kubo_viscosity(p; tmax=1.0)
-------------------------------------------------------------------
GreenKuboViscosity:
-------------------------------------------------------------------
Number of time lags: 501
Maximum time lag: 1.0000 ps
Viscosity at maximum time lag: 0.1799 mPa·s
-------------------------------------------------------------------

julia> gk.viscosity[end]
0.1799246341696098
```

!!! compat
    This function was added in version 2.7.0 of MolSimToolkit.

"""
function green_kubo_viscosity(
    p::PressureTensor;
    temperature::Real=p.temperature,
    volume::Real=p.volume,
    tmax::Real=10.0,
    components::Symbol=:traceless,
    parallel::Bool=true,
)
    n = length(p)
    if n < 2
        throw(ArgumentError("At least two pressure tensors are required."))
    end
    stride = p.steps[2] - p.steps[1]
    if any(!=(stride), diff(p.steps))
        throw(ArgumentError("""\n
            The pressure tensors must be printed at regular step intervals.

        """))
    end
    if !(temperature > 0 && volume > 0)
        throw(ArgumentError("temperature and volume must be positive."))
    end
    dt = stride * p.timestep / 1000 # ps
    if !(dt > 0)
        throw(ArgumentError("Invalid time step: $(p.timestep) fs."))
    end
    maxlag = min(floor(Int, tmax / dt + sqrt(eps(Float64))), n ÷ 2)
    if maxlag < 1
        throw(ArgumentError("tmax must be at least the time between two pressure tensors ($dt ps)."))
    end

    # Each column of x is a component of the tensor, and w its weight in the sum
    # of autocorrelations.
    sym = [(P + P') / 2 for P in p.tensors]
    x, w = if components == :offdiagonal
        ([P[i, j] for P in sym, (i, j) in ((1, 2), (1, 3), (2, 3))], fill(1 / 3, 3))
    elseif components == :traceless
        # The off-diagonal components appear twice in the sum over all components.
        tl = [P - (sum(diag(P)) / 3) * one(P) for P in sym]
        ([P[i, j] for P in tl, (i, j) in ((1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3))],
            [1, 1, 1, 2, 2, 2] ./ 10)
    else
        throw(ArgumentError("components must be :traceless or :offdiagonal, got :$components."))
    end
    x .-= mean(x; dims=1)

    # Each lag is independent and writes to a distinct entry of `acf`. The cost
    # per lag decreases with the lag, hence the round-robin split.
    acf = zeros(maxlag + 1)
    @threads for lags in chunks(0:maxlag; n=parallel ? Threads.nthreads() : 1, split=RoundRobin())
        for lag in lags
            s = 0.0
            for ic in axes(x, 2)
                sc = 0.0
                @simd for t in 1:(n-lag)
                    @inbounds sc += x[t, ic] * x[t+lag, ic]
                end
                s += w[ic] * sc
            end
            acf[lag+1] = s / (n - lag)
        end
    end

    # η = V/(kB T) ∫ acf dt, with V in Å³, acf in bar², t in ps. Converted to mPa·s:
    # Å³ → 1e-30 m³, bar² → 1e10 Pa², ps → 1e-12 s, Pa·s → 1e3 mPa·s.
    kB = 1.380649e-23 # J/K
    factor = volume * 1e-30 * 1e10 * 1e-12 * 1e3 / (kB * temperature)
    viscosity = zeros(maxlag + 1)
    for i in 2:maxlag+1
        viscosity[i] = viscosity[i-1] + factor * dt * (acf[i-1] + acf[i]) / 2
    end
    time = [dt * i for i in 0:maxlag]
    return GreenKuboViscosity(time, acf, viscosity)
end

@testitem "green_kubo_viscosity" begin
    using MolSimToolkit, MolSimToolkit.Testing
    using MolSimToolkit: PressureTensor
    using StaticArrays: SMatrix
    using Statistics: mean

    p = read_namd_pressure_tensor(Testing.namd_pressure_log)
    gk = green_kubo_viscosity(p; tmax=1.0)
    @test length(gk.time) == 501
    @test gk.time[1] == 0.0
    @test gk.time[end] ≈ 1.0
    @test gk.viscosity[1] == 0.0
    @test gk.viscosity[end] ≈ 0.1799246341696098
    @test gk.acf ≈ green_kubo_viscosity(p; tmax=1.0, parallel=false).acf
    gko = green_kubo_viscosity(p; tmax=1.0, components=:offdiagonal)
    @test gko.viscosity[end] ≈ 0.23209654612215955
    # Linear scaling with the volume and inverse scaling with temperature
    @test green_kubo_viscosity(p; tmax=1.0, volume=2 * p.volume).viscosity ≈ 2 * gk.viscosity
    @test green_kubo_viscosity(p; tmax=1.0, temperature=2 * p.temperature).viscosity ≈ gk.viscosity / 2
    # tmax limited to half of the trajectory
    @test length(green_kubo_viscosity(p; tmax=1000.0).time) == length(p) ÷ 2 + 1

    # Check the autocorrelation function against a direct computation for a
    # synthetic, isotropic case: with only the xy component nonzero, the off-diagonal
    # estimate is 1/3 of its autocorrelation, and the traceless one is 2/10 of it.
    n = 100
    xy = sin.(0.3 .* (1:n))
    xy .-= mean(xy)
    tensors = [SMatrix{3,3,Float64,9}(0, v, 0, v, 0, 0, 0, 0, 0) for v in xy]
    pt = PressureTensor(collect(0:n-1), tensors, 1.0, 1.0, 1.0)
    direct = [sum(xy[t] * xy[t+lag] for t in 1:n-lag) / (n - lag) for lag in 0:10]
    @test green_kubo_viscosity(pt; tmax=0.010, components=:offdiagonal).acf ≈ direct ./ 3
    @test green_kubo_viscosity(pt; tmax=0.010, components=:traceless).acf ≈ direct .* (2 / 10)

    # Errors
    @test_throws ArgumentError green_kubo_viscosity(p; components=:diagonal)
    @test_throws ArgumentError green_kubo_viscosity(p; tmax=0.001)
    @test_throws ArgumentError green_kubo_viscosity(p; temperature=-1.0)
    irregular = PressureTensor([0, 1, 3], tensors[1:3], 1.0, 1.0, 1.0)
    @test_throws ArgumentError green_kubo_viscosity(irregular)
end
