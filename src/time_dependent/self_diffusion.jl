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

"""
    mean_square_displacement(
        sim::Simulation,
        selection::AbstractVector{<:PDBTools.Atom};
        natomspermol::Integer,
        maxdelta::Integer = length(sim) ÷ 10,
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
diffuse across half a box length, this reconstruction is unambiguous.

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
1074.6758528050218

```

!!! compat
    This function was added in version 2.4.0 of MolSimToolkit.

"""
function mean_square_displacement(
    sim::Simulation,
    selection::AbstractVector{<:PDBTools.Atom};
    natomspermol::Integer,
    maxdelta::Integer=max(1, length(sim) ÷ 10),
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

    prg = Progress(n_frames; enabled=show_progress)
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

    msd = OffsetArrays.OffsetArray(zeros(maxdelta + 1), 0:maxdelta)
    for delta in 0:maxdelta
        s = 0.0
        n = 0
        for t in 1:(n_frames-delta), imol in 1:n_molecules
            d = coms[imol, t+delta] - coms[imol, t]
            s += sum(abs2, d)
            n += 1
        end
        msd[delta] = s / n
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
  used to compute `msd`, used to convert `delta` (in frames) into time
  units. Defaults to `1`, that is, `delta` is used directly, and the
  resulting units of `D` are squared-length per frame.
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
