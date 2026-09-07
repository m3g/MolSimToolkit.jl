using Base.Threads: @threads
using ChunkSplitters: chunks, RoundRobin

#=
    _unwrap(raw::AbstractVector, ucs::AbstractVector) -> Vector

Given the (possibly PBC-wrapped) trajectory `raw` of a single particle (or
the center of mass of a single molecule), and the unit cell `ucs[i]`
associated to each frame `i`, returns the "unwrapped" trajectory: each
position is replaced by the periodic image closest to the (already
unwrapped) position of the previous frame. This way, a particle that crosses
a periodic boundary between two frames produces a continuous displacement
instead of a spurious jump of about one box length.
=#
function _unwrap(raw::AbstractVector{T}, ucs::AbstractVector) where {T}
    n = length(raw)
    unwrapped = Vector{T}(undef, n)
    isempty(raw) && return unwrapped
    unwrapped[1] = raw[1]
    for i in 2:n
        unwrapped[i] = wrap(raw[i], unwrapped[i-1], ucs[i])
    end
    return unwrapped
end

@testitem "_warn_on_large_displacements" begin
    using MolSimToolkit: _warn_on_large_displacements
    using StaticArrays: SMatrix
    uc = UnitCell(SMatrix{3,3,Float64,9}(10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0), true, true)
    ucs = [uc, uc, uc]

    # Small displacements: nothing to report.
    coms = [Point3D(0.0, 0.0, 0.0) Point3D(1.0, 0.0, 0.0) Point3D(2.0, 0.0, 0.0)]
    @test _warn_on_large_displacements(coms, ucs) == 0

    # A displacement of 4.9 out of a cell vector of 10.0 is essentially at the
    # limit above which the minimum image is not necessarily the correct one.
    coms = [Point3D(0.0, 0.0, 0.0) Point3D(4.9, 0.0, 0.0) Point3D(9.8, 0.0, 0.0)]
    @test (@test_logs (:warn,) _warn_on_large_displacements(coms, ucs)) == 2

    # Without unit cell information there is nothing to check.
    invalid = UnitCell(SMatrix{3,3,Float64,9}(ntuple(_ -> 0.0, 9)...), false, false)
    @test _warn_on_large_displacements(coms, [invalid, invalid, invalid]) == 0
end

@testitem "_unwrap" begin
    using MolSimToolkit: _unwrap
    using StaticArrays: SMatrix
    mat = SMatrix{3,3,Float64,9}(10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0)
    uc = UnitCell(mat, true, true)
    ucs = [uc, uc, uc, uc]

    # A particle drifting steadily in +x, wrapped back into [0, 10) by the
    # trajectory file, must be reconstructed as a continuous trajectory.
    raw = [Point3D(9.0, 0.0, 0.0), Point3D(1.0, 0.0, 0.0), Point3D(3.0, 0.0, 0.0), Point3D(5.0, 0.0, 0.0)]
    unwrapped = _unwrap(raw, ucs)
    @test unwrapped ≈ [Point3D(9.0, 0.0, 0.0), Point3D(11.0, 0.0, 0.0), Point3D(13.0, 0.0, 0.0), Point3D(15.0, 0.0, 0.0)]

    # A particle that does not cross any boundary is left untouched.
    raw2 = [Point3D(1.0, 2.0, 3.0), Point3D(1.5, 2.1, 2.9), Point3D(2.0, 2.2, 2.8)]
    @test _unwrap(raw2, ucs) ≈ raw2

    # Drifting in the opposite direction must also be reconstructed correctly.
    raw3 = [Point3D(1.0, 0.0, 0.0), Point3D(9.0, 0.0, 0.0), Point3D(7.0, 0.0, 0.0)]
    unwrapped3 = _unwrap(raw3, ucs)
    @test unwrapped3 ≈ [Point3D(1.0, 0.0, 0.0), Point3D(-1.0, 0.0, 0.0), Point3D(-3.0, 0.0, 0.0)]
end

#=
    _warn_on_large_displacements(coms::AbstractMatrix, ucs::AbstractVector)

Checks whether any of the unwrapped trajectories in `coms` (one row per
molecule, one column per frame) displaces, between two consecutive frames, by
more than `_max_relative_displacement` times a unit cell vector.

The minimum-image reconstruction performed by `_unwrap` always returns
displacements whose components, in units of the cell vectors, are not greater
than `0.5` in absolute value; it is only *correct* if the true displacement is
well below that limit, since a molecule moving almost half a cell vector could
equally well be moving in the opposite direction. A single warning is emitted
reporting the worst case found.
=#
const _max_relative_displacement = 0.45

function _warn_on_large_displacements(coms::AbstractMatrix, ucs::AbstractVector)
    n_molecules, n_frames = size(coms)
    worst = 0.0
    worst_imol = 0
    worst_iframe = 0
    count = 0
    for iframe in 2:n_frames
        uc = ucs[iframe]
        !uc.valid && continue
        invmatrix = inv(uc.matrix)
        for imol in 1:n_molecules
            d = maximum(abs, invmatrix * (coms[imol, iframe] - coms[imol, iframe-1]))
            d <= _max_relative_displacement && continue
            count += 1
            if d > worst
                worst, worst_imol, worst_iframe = d, imol, iframe
            end
        end
    end
    if count > 0
        @warn """\n
            Displacements of about half a unit cell vector were found between consecutive frames.

            The unwrapping of the coordinates is ambiguous in this case, and the mean square
            displacement is likely wrong. This usually means that the frames are too far apart
            in time for the molecules considered, or that the coordinates of each molecule are
            not contiguous in the trajectory file.

            Number of occurrences: $count (out of $(n_molecules * (n_frames - 1)) displacements)
            Worst case: molecule $worst_imol, frame $worst_iframe, displacement of $(round(worst; digits=3)) cell vectors

        """ _module = nothing _file = nothing
    end
    return count
end

"""
    mean_square_displacement(
        sim::Simulation,
        selection::AbstractVector{<:PDBTools.Atom};
        natomspermol::Integer,
        maxdelta::Integer = length(sim) ÷ 10,
        parallel::Bool = true,
        show_progress::Bool = true,
    )

Computes the mean square displacement (MSD), as a function of the time lag
`delta` (in frames), of the centers of mass of the molecules defined by
`selection`.

`selection` is expected to be the concatenation of `n` identical molecules,
each with `natomspermol` atoms (the same convention used, for instance, by
[`occupancy`](@ref)).

Trajectory files typically store coordinates wrapped back into the
simulation box by periodic boundary conditions, which would introduce
spurious jumps in a molecule's displacement whenever it crosses a periodic
boundary between two frames. To avoid this, each molecule's center of mass
is reconstructed into a continuous ("unwrapped") trajectory: at each frame,
the periodic image closest to that same molecule's (already unwrapped)
position at the previous frame is chosen. Since consecutive frames are
assumed to be much closer in time than the time it takes a molecule to
diffuse across half a box length, this reconstruction is unambiguous. If that
assumption is violated, that is, if some molecule is reconstructed as moving
by about half a unit cell vector between two consecutive frames, a warning is
emitted (once per call), since the resulting MSD is then likely meaningless.

Returns an `OffsetArray` with indices `0:maxdelta`, in squared length units
(typically Å², matching the units of the input coordinates). The value at
`delta` is the average, over all molecules and all pairs of frames separated
by `delta` steps, of the squared displacement of the molecule's center of
mass.

Use [`self_diffusion_coefficient`](@ref) to estimate the self-diffusion
coefficient from the linear (diffusive) regime of the resulting curve.

# Arguments

- `sim::Simulation`: The `Simulation` object.
- `selection::AbstractVector{<:PDBTools.Atom}`: The atoms of the `n`
  molecules considered, concatenated in sequence.

# Optional keyword arguments

- `natomspermol::Integer`: Number of atoms of each molecule. Required.
- `maxdelta::Integer`: The maximum delta-step to be considered. Defaults to
  `length(sim) ÷ 10`.
- `parallel::Bool`: Defines if the averaging over time lags is run in
  parallel. Defaults to `true`. Requires starting Julia with multi-threading.
- `show_progress::Bool`: Show progress bar. Defaults to `true`.

# Example

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2***"
julia> using MolSimToolkit, PDBTools, MolSimToolkit.Testing

julia> sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj);

julia> tmao = select(get_atoms(sim), "resname TMAO");

julia> msd = mean_square_displacement(sim, tmao; natomspermol=14, maxdelta=4, show_progress=false);

julia> msd[0]
0.0

julia> msd[4]
878.4761526314077

```

!!! compat
    This function was added in version 2.4.0 of MolSimToolkit.

"""
function mean_square_displacement(
    sim::Simulation,
    selection::AbstractVector{<:PDBTools.Atom};
    natomspermol::Integer,
    maxdelta::Integer=max(1, length(sim) ÷ 10),
    parallel::Bool=true,
    show_progress::Bool=true,
)
    if length(selection) % natomspermol != 0
        throw(ArgumentError("""\n
            Number of atoms in the selection is not a multiple of natomspermol

        """))
    end
    if maxdelta > length(sim) - 1
        throw(ArgumentError("maxdelta must be less than the number of frames minus 1"))
    end
    inds = PDBTools.index.(selection)
    n_molecules = length(inds) ÷ natomspermol
    mol_indices = [inds[(i-1)*natomspermol+1:i*natomspermol] for i in 1:n_molecules]

    n_frames = length(sim)
    f = first_frame!(sim)
    com_type = typeof(center_of_mass(mol_indices[1], sim, positions(f)))
    raw_coms = Matrix{com_type}(undef, n_molecules, n_frames)
    ucs = Vector{UnitCell}(undef, n_frames)

    prg = Progress(n_frames; enabled=show_progress, desc="Computing displacements:")
    for (iframe, frame) in enumerate(sim)
        p = positions(frame)
        ucs[iframe] = unitcell(frame)
        for imol in 1:n_molecules
            raw_coms[imol, iframe] = center_of_mass(mol_indices[imol], sim, p)
        end
        next!(prg)
    end

    coms = similar(raw_coms)
    for imol in 1:n_molecules
        coms[imol, :] .= _unwrap(@view(raw_coms[imol, :]), ucs)
    end

    # The unwrapping above assumes that a molecule moves much less than half a
    # unit cell vector between two consecutive frames; otherwise the closest
    # periodic image is not necessarily the physically correct one. Check that
    # assumption and warn once if it is violated.
    _warn_on_large_displacements(coms, ucs)

    msd = OffsetArrays.OffsetArray(zeros(maxdelta + 1), 0:maxdelta)
    # The work per `delta` is proportional to `(n_frames - delta) * n_molecules`, so
    # count the progress in units of inner iterations for a bar that advances evenly.
    ninner = n_molecules * ((maxdelta + 1) * n_frames - (maxdelta * (maxdelta + 1)) ÷ 2)
    prg = Progress(ninner; enabled=show_progress, desc="Averaging per delta:")
    # A single `delta` can take a long time, so the counter is advanced from within
    # the loop over `t` rather than once per `delta`; otherwise the bar sits blank
    # and then jumps. Flushing in blocks of about 1/500 of the total keeps the
    # display smooth while leaving `next!` (which takes a lock) out of the hot loop.
    flush_every = max(n_molecules, ninner ÷ 500)
    next!(prg; step=0, force=true) # paint the bar at 0% before any work is done
    # Each `delta` is independent and writes to a single, distinct entry of `msd`,
    # so the loop can be split among threads without any reduction. The cost per
    # `delta` decreases with `delta`, hence the round-robin split, which gives every
    # chunk a similar mix of cheap and expensive time lags.
    nchunks = parallel ? Threads.nthreads() : 1
    @threads for deltas in chunks(0:maxdelta; n=nchunks, split=RoundRobin())
        for delta in deltas
            s = 0.0
            n = 0
            nsince = 0
            for t in 1:(n_frames-delta)
                for imol in 1:n_molecules
                    d = coms[imol, t+delta] - coms[imol, t]
                    s += sum(abs2, d)
                end
                n += n_molecules
                nsince += n_molecules
                if nsince >= flush_every
                    next!(prg; step=nsince)
                    nsince = 0
                end
            end
            msd[delta] = s / n
            nsince > 0 && next!(prg; step=nsince)
        end
    end
    return msd
end

"""
    self_diffusion_coefficient(
        msd::AbstractVector;
        dt::Real = 1.0,
        dim::Integer = 3,
        mindelta::Integer = 1,
        maxdelta::Integer = lastindex(msd),
    )

Estimates the self-diffusion coefficient from a mean-square-displacement
curve `msd`, as computed by [`mean_square_displacement`](@ref), using the
Einstein relation `MSD(t) = 2 * dim * D * t`.

The coefficient is obtained as the slope of an ordinary least-squares fit of
`msd[delta]` against `t = delta * dt`, for `delta` in `mindelta:maxdelta`,
divided by `2 * dim`.

!!! note "Units"
    No unit conversion is performed: `2 * dim` is a dimensionless factor, so
    `D` comes out in `[msd units] / [dt units]`. In particular, if `msd` is
    in Å² (as returned by `mean_square_displacement`, for coordinates in Å)
    and `dt` is given in ps, the returned coefficient is in Å²/ps. Convert
    it to other units (e.g. cm²/s) yourself if needed.

!!! note
    Only a limited range of `delta` typically falls in the diffusive
    (linear) regime: very short times are dominated by ballistic motion, and
    very long times become noisy, since fewer pairs of frames are available
    to average over as `delta` approaches the length of the trajectory.
    Inspect the `msd` curve (for instance, by plotting it) to choose
    `mindelta` and `maxdelta` accordingly.

# Arguments

- `msd::AbstractVector`: The mean-square-displacement curve, as computed by
  [`mean_square_displacement`](@ref).

# Optional keyword arguments

- `dt::Real`: The time interval between consecutive frames of the trajectory
  used to compute `msd`, in whatever time unit is desired for the output
  (e.g. ps). Used to convert `delta` (in frames) into that time unit.
  Defaults to `1`, that is, `delta` is used directly, and the resulting
  units of `D` are squared-length per frame.
- `dim::Integer`: The dimensionality of the diffusion process. Defaults to
  `3` (three-dimensional diffusion).
- `mindelta::Integer`: The smallest `delta` (in frames) included in the fit.
  Defaults to `1`, excluding the trivial `delta = 0` point.
- `maxdelta::Integer`: The largest `delta` (in frames) included in the fit.
  Defaults to `lastindex(msd)`, that is, the whole curve.

# Example

```jldoctest
julia> using MolSimToolkit, OffsetArrays

julia> msd = OffsetArray([0.0, 2.0, 4.0, 6.0, 8.0], 0:4); # a purely diffusive MSD, dt=1, dim=3

julia> self_diffusion_coefficient(msd)
0.3333333333333334

```

!!! compat
    This function was added in version 2.4.0 of MolSimToolkit.

"""
function self_diffusion_coefficient(
    msd::AbstractVector;
    dt::Real=1.0,
    dim::Integer=3,
    mindelta::Integer=1,
    maxdelta::Integer=lastindex(msd),
)
    if !(firstindex(msd) <= mindelta < maxdelta <= lastindex(msd))
        throw(ArgumentError("""\n
            mindelta and maxdelta must satisfy firstindex(msd) <= mindelta < maxdelta <= lastindex(msd)

        """))
    end
    deltas = mindelta:maxdelta
    t = collect(deltas) .* dt
    y = [msd[delta] for delta in deltas]
    X = [ones(length(t)) t]
    slope = (X \ y)[2]
    return slope / (2 * dim)
end

@testitem "mean_square_displacement" begin
    using MolSimToolkit, PDBTools, MolSimToolkit.Testing
    using OffsetArrays

    sim = Simulation(Testing.namd_pdb, Testing.namd_traj)
    protein = select(get_atoms(sim), "protein")

    # A single "molecule" (the whole protein): delta=0 must be exactly zero,
    # and the result must match a naive, non-unwrapped computation, since a
    # protein's center of mass does not cross the periodic boundary within
    # this short test trajectory.
    msd = mean_square_displacement(sim, protein; natomspermol=length(protein), maxdelta=3, show_progress=false)
    @test msd isa OffsetVector
    @test msd[0] == 0.0
    @test all(>=(0), msd)

    inds = PDBTools.index.(protein)
    coms = [center_of_mass(inds, sim, positions(frame)) for frame in sim]
    for delta in 1:3
        expected = sum(sum(abs2, coms[t+delta] - coms[t]) for t in 1:(length(coms)-delta)) / (length(coms) - delta)
        @test msd[delta] ≈ expected
    end

    # Multiple molecules: TMAO, using the namd2 trajectory
    sim2 = Simulation(Testing.namd2_pdb, Testing.namd2_traj)
    tmao = select(get_atoms(sim2), "resname TMAO")
    msd2 = mean_square_displacement(sim2, tmao; natomspermol=14, maxdelta=4, show_progress=false)
    @test msd2[0] == 0.0
    @test length(msd2) == 5
    @test all(>=(0), msd2)

    @test_throws "not a multiple" mean_square_displacement(sim2, tmao; natomspermol=13, show_progress=false)
    @test_throws ArgumentError mean_square_displacement(sim2, tmao; natomspermol=14, maxdelta=length(sim2), show_progress=false)

    # The center of mass of each molecule must be computed relative to an atom of
    # that same molecule; otherwise molecules straddling the boundary of the box
    # centered at the reference are split into different periodic images, and the
    # resulting MSD is much larger than the true one. Computing the MSD of a single
    # atom of each molecule bypasses that step entirely, so both must agree.
    single = select(get_atoms(sim2), "resname TMAO and name N")
    msd_single = mean_square_displacement(sim2, single; natomspermol=1, maxdelta=4, show_progress=false)
    for delta in 1:4
        @test msd2[delta] ≈ msd_single[delta] rtol = 0.2
    end

    # Frames too far apart for the unwrapping to be unambiguous must be reported
    water = select(get_atoms(sim2), "water")
    @test_logs (:warn,) match_mode = :any mean_square_displacement(sim2, water; natomspermol=3, maxdelta=2, show_progress=false)
end

@testitem "self_diffusion_coefficient" begin
    using MolSimToolkit
    using OffsetArrays

    # Purely diffusive, noise-free MSD: D = slope / (2*dim)
    msd = OffsetArray(collect(0.0:2.0:8.0), 0:4)
    @test self_diffusion_coefficient(msd) ≈ 2.0 / 6.0
    @test self_diffusion_coefficient(msd; dim=1) ≈ 2.0 / 2.0
    @test self_diffusion_coefficient(msd; dt=2.0) ≈ (2.0 / 2.0) / 6.0

    # Restricting the fit range must not change the result for a perfectly linear curve
    @test self_diffusion_coefficient(msd; mindelta=2, maxdelta=4) ≈ self_diffusion_coefficient(msd)

    # delta=0 is a legitimate (trivial) data point, and may be included explicitly
    @test self_diffusion_coefficient(msd; mindelta=0) ≈ self_diffusion_coefficient(msd)

    @test_throws ArgumentError self_diffusion_coefficient(msd; mindelta=3, maxdelta=2)
    @test_throws ArgumentError self_diffusion_coefficient(msd; maxdelta=5)
    @test_throws ArgumentError self_diffusion_coefficient(msd; mindelta=-1)

    # Integration with mean_square_displacement, using real trajectory data
    using PDBTools, MolSimToolkit.Testing
    sim = Simulation(Testing.namd2_pdb, Testing.namd2_traj)
    tmao = select(get_atoms(sim), "resname TMAO")
    msd2 = mean_square_displacement(sim, tmao; natomspermol=14, maxdelta=4, show_progress=false)
    D = self_diffusion_coefficient(msd2)
    @test D isa Float64
end
