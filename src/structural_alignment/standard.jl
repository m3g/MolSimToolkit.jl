"""
    reconstruct_structure!(
        x::AbstractVector{<:AbstractVector},
        indices::AbstractVector{<:Integer},
        unitcell::Union{UnitCell,AbstractMatrix};
        dmax::Real = minimum(norm.(eachcol(unitcell isa UnitCell ? unitcell.matrix : unitcell))) / 3,
        aux_inds::AbstractVector{<:Integer} = similar(indices),
    )

Reconstructs, in place, the structure formed by the atoms of `x` with indices given by
`indices`, undoing the breaks caused by periodic boundary conditions. The atoms are assumed
to be ordered such that atoms that are close in the sequence of `indices` are also close in
space, which is typically the case for the atoms of a protein or polymer, taken in sequence.

The reconstruction walks `indices` in order, wrapping each atom around the previous one. If
the next atom in the sequence is farther than `dmax` from the previous one (e.g. a chain
break, or the transition from one molecule to another when `indices` was built by
concatenating more than one structure), the remaining, not yet reconstructed atoms (both those
skipped over and those still ahead in the sequence) are searched for the one closest to the
last reconstructed position, and the walk branches to it. The atoms skipped over are then
reconstructed backwards from that branch point, before resuming towards the atoms that were
still ahead. This search-and-branch step is applied recursively to whatever remains
unreconstructed at each point, so nested or repeated breaks (in either direction) are all
resolved the same way, and an atom skipped over by one branch can still be picked up by a
later one. This makes the function robust to structures formed by more than one chain (e.g. a
dimer, or a set of indices obtained by concatenating a structure and a separate complex), for
which a single branch is typically enough to reconnect the chains.

# Arguments

- `x`: the coordinates of the atoms, modified in place.
- `indices`: the indices, within `x`, of the atoms to be reconstructed, in sequence.
- `unitcell`: the unit cell of the simulation box, as an `UnitCell` object or a matrix of
  unit cell vectors.

# Optional keyword arguments

- `dmax`: maximum distance, between consecutive atoms of `indices`, below which no branch
  search is triggered. Defaults to one third of the smallest unit cell vector length.
- `aux_inds`: a preallocated buffer, with the same axes and element type as `indices`, used
  internally to avoid repeated allocations when this function is called repeatedly (e.g.
  once per trajectory frame).

# Example

```jldoctest
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> sim = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> frame = first_frame!(sim);

julia> x = positions(frame);

julia> reconstruct_structure!(x, 1:length(x), unitcell(frame));

```

"""
function reconstruct_structure!(
    x::AbstractVector{<:AbstractVector},
    indices::AbstractVector{<:Integer},
    unitcell::Union{UnitCell,AbstractMatrix};
    dmax::Real=minimum(norm.(eachcol(unitcell isa UnitCell ? unitcell.matrix : unitcell))) / 3,
    aux_inds::AbstractVector{<:Integer}=similar(indices),
)
    length(indices) <= 1 && return x
    aux_inds .= indices
    xanchor = x[aux_inds[firstindex(aux_inds)]]
    seq = @view aux_inds[firstindex(aux_inds)+1:lastindex(aux_inds)]
    _reconstruct_walk!(x, seq, xanchor, unitcell, dmax)
    return x
end

#
# Reconstructs, in place, the atoms of `x` whose indices are given by `seq`, anchoring the
# first one to `xprev`. `seq` is entirely the pool of atoms not yet reconstructed. At each
# step, if the atom at the front of the pool is within `dmax` of the last reconstructed
# position, it is accepted and the walk simply moves on to the next one. Otherwise, the
# whole of the pool is searched for the atom closest to that position, which becomes the
# branch point; the atoms that were skipped over are merged, in reverse (so that they are
# visited backwards from the branch point), with the atoms that were still ahead, and this
# combined list becomes the new pool.
#
# Because the pool considered at every step always contains every not-yet-reconstructed
# atom -- both the ones skipped by a branch and the ones still ahead in the original
# sequence -- an atom left over from one branch remains a valid candidate for a later one,
# however many breaks are nested or chained together. The pool is only ever rebuilt (by
# `vcat`) when an actual branch occurs, which is expected to be rare, so this stays an
# in-place loop rather than a recursion over one call per atom (which would, for large
# structures, overflow the call stack).
#
function _reconstruct_walk!(
    x::AbstractVector{<:AbstractVector},
    seq::AbstractVector{<:Integer},
    xprev,
    unitcell::Union{UnitCell,AbstractMatrix},
    dmax::Real,
)
    k = firstindex(seq)
    while k <= lastindex(seq)
        xk = wrap(x[seq[k]], xprev, unitcell)
        if norm(xk - xprev) < dmax
            x[seq[k]] = xk
            xprev = xk
            k += 1
            continue
        end
        # Gap: search the remaining atoms of the pool for the one closest to xprev
        kmin, dmin, xmin = k, norm(xk - xprev), xk
        for kk in k+1:lastindex(seq)
            xkk = wrap(x[seq[kk]], xprev, unitcell)
            dkk = norm(xkk - xprev)
            if dkk < dmin
                kmin, dmin, xmin = kk, dkk, xkk
            end
        end
        x[seq[kmin]] = xmin
        # Keep the skipped atoms (to be visited backwards from the branch point) and the
        # atoms still ahead in a single pool, so that a gap in either can still be bridged
        # by an atom from the other.
        seq = vcat(@view(seq[kmin-1:-1:k]), @view(seq[kmin+1:lastindex(seq)]))
        xprev = xmin
        k = firstindex(seq)
    end
    return nothing
end

#
# Builds the sequence of indices used to jointly reconstruct `indices` and
# `rmsd_indices`: `indices` followed by the atoms of `rmsd_indices` that are
# not already part of it. Reconstructing this combined sequence with
# `reconstruct_structure!` reconnects `rmsd_indices` to `indices` (branching
# from the closest atom, as it does for chain breaks within a single
# structure), avoiding the need for a separate complex-reconstruction step.
#
function _reconstruction_indices(indices::AbstractVector{<:Integer}, rmsd_indices::AbstractVector{<:Integer})
    rmsd_indices === indices && return indices
    indices_set = Set(indices)
    extra = filter(!in(indices_set), rmsd_indices)
    isempty(extra) && return indices
    return vcat(indices, extra)
end

#
# Function to read the reference coordinates from the trajectory.
# returns the coordinates, reconstructed by reconstruct_structure!.
#
function _reference_coordinates(
    simulation,
    reference_frame,
    indices,
    rmsd_indices,
    mass;
    show_progress,
)
    combined_indices = _reconstruction_indices(indices, rmsd_indices)
    aux_inds = similar(combined_indices)
    xalign, xrmsd = if isnothing(reference_frame)
        first_frame!(simulation)
        frame = current_frame(simulation)
        p, uc = positions(frame), unitcell(frame)
        reconstruct_structure!(p, combined_indices, uc.matrix; aux_inds)
        p[indices], p[rmsd_indices]
    elseif reference_frame isa Integer
        if !(reference_frame in frame_range(simulation))
            throw(ArgumentError("""\n

                reference_frame index $reference_frame is not in frame range of the simulation: $(frame_range(simulation))

            """))
        end
        restart!(simulation)
        for _ in 1:reference_frame-1
            next_frame!(simulation)
        end
        frame = current_frame(simulation)
        p, uc = positions(frame), unitcell(frame)
        reconstruct_structure!(p, combined_indices, uc.matrix; aux_inds)
        p[indices], p[rmsd_indices]
    elseif reference_frame == :average
        xm = zeros(3, length(indices))
        xp = zeros(3, length(indices))
        xalign_ref = zeros(Point3D, length(indices))
        xrmsd_ref = zeros(Point3D, length(rmsd_indices))
        prg = Progress(length(simulation); enabled=show_progress, desc="Computing average structure:")
        for (iframe, frame) in enumerate(simulation)
            next!(prg)
            p, uc = positions(frame), unitcell(frame)
            reconstruct_structure!(p, combined_indices, uc.matrix; aux_inds)
            if iframe == 1
                xalign_ref .= @view(p[indices])
                xrmsd_ref .= @view(p[rmsd_indices])
            else
                # Align atoms to reference
                xalign = @view(p[indices])
                xcm, xref_cm, u = alignment_movements(xalign, xalign_ref; mass, xm, xp)
                apply_alignment_transformation!(xalign, xcm, xref_cm, u)
                @. xalign_ref = (xalign_ref * (iframe - 1) + xalign) / iframe
                if !(indices === rmsd_indices)
                    # Move rmsd atoms with same transformation
                    xrmsd = @view(p[rmsd_indices])
                    apply_alignment_transformation!(xrmsd, xcm, xref_cm, u)
                    @. xrmsd_ref = (xrmsd_ref * (iframe - 1) + xrmsd) / iframe
                else
                    xrmsd_ref = xalign_ref
                end
            end
        end
        xalign_ref, xrmsd_ref
    else
        throw(ArgumentError("""\n

            reference_frame must be an integer, nothing or :average
            
        """))
    end
    restart!(simulation)
    return xalign, xrmsd
end

"""
    rmsd(
        simulation::Simulation, 
        indices::AbstractVector{<:Integer}; # or sel::AbstractString 
        rmsd_of::AbstractVector{<:Integer} = indices, # or sel::AbstractString
        mass = nothing, 
        reference_frame = nothing, 
        show_progress = true
    )

Computes the root mean square deviation (RMSD) between two sets of points in along a trajectory.

# Positional arguments

- `indices`: vector with indices of all the atoms of the structure to be aligned (i. e. the protein or polymer). 
  By default, they are also the atoms for which the rmsd will be computed.
or
- `sel`: Selection string, following PDBTools.jl syntax, to define the atoms to be considered.

# Optional keyword arguments

- `rmsd_of`: indices of the atoms for which the rmsd will be computed, considering the alignment of the
  atoms defined in `indices`.
- `mass` argument can be used to provide the mass of the atoms if they are not the same.
- `reference_frame` argument can be used to provide a reference frame to align the trajectory to:
    - If `reference_frame == nothing`, the first frame will be used (default behavior).
    - If `reference_frame == :average`, the average structure will be used.
    - If `reference_frame` is an integer, the frame at that index will be used as reference.

# Examples

## Computing the rmsd along a trajectory

```jldoctest; filter = r"([0-9]+\\.[0-9]{2})[0-9]+" => s"\\1***"
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> using PDBTools

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> rmsd(simulation, "protein and name CA"; show_progress=false)
5-element Vector{Float64}:
 0.0
 2.8388710154609034
 2.9776998440690385
 2.4621444212469483
 3.8035683196100796

julia> rmsd(simulation, "protein and name CA"; reference_frame=:average, show_progress=false)
5-element Vector{Float64}:
 1.8995986972454748
 2.1512244220536973
 1.5081703191869376
 1.1651111324544219
 2.757039151265317
```

"""
function rmsd(
    simulation::Simulation, indices::AbstractVector{<:Integer};
    rmsd_of::AbstractVector{<:Integer}=indices,
    mass=nothing,
    reference_frame=nothing,
    show_progress=true,
)

    # Auxiliary arrays for the alignment
    xm = zeros(3, length(indices))
    xp = zeros(3, length(indices))

    # Define reference of the alignment
    xalign_ref, xrmsd_ref =
        _reference_coordinates(simulation, reference_frame, indices, rmsd_of, mass; show_progress)

    # Indices and auxiliary buffer used to jointly reconstruct `indices` and `rmsd_of`
    combined_indices = _reconstruction_indices(indices, rmsd_of)
    aux_inds = similar(combined_indices)

    rmsds = Float64[]
    prg = Progress(length(simulation); enabled=show_progress, desc="Computing RMSDs for each frame:")
    for frame in simulation
        next!(prg)
        p, uc = positions(frame), unitcell(frame)
        reconstruct_structure!(p, combined_indices, uc.matrix; aux_inds)
        # Obtain the transformation that aligns atoms of `indices`
        xalign = @view(p[indices])
        xcm, xref_cm, u = alignment_movements(xalign, xalign_ref; mass, xm, xp)
        # Apply transformation to atoms in rmsd_of
        xrmsd = @view(p[rmsd_of])
        apply_alignment_transformation!(xrmsd, xcm, xref_cm, u)
        push!(rmsds, rmsd(xrmsd, xrmsd_ref))
    end
    return rmsds
end

function rmsd(
    simulation::Simulation, sel::AbstractString;
    rmsd_of::AbstractString=sel,
    mass=nothing,
    reference_frame=nothing,
    show_progress=true,
)
    indices = findall(PDBTools.Select(sel), get_atoms(simulation))
    rmsd_indices = if rmsd_of === sel
        indices
    else
        findall(PDBTools.Select(rmsd_of), get_atoms(simulation))
    end
    rmsd(simulation, indices; rmsd_of=rmsd_indices, mass, reference_frame, show_progress)
end

@testitem "rmsd" begin
    using MolSimToolkit
    using MolSimToolkit.Testing: namd_pdb, namd_traj
    using StaticArrays: SVector
    using PDBTools
    using Rotations: RotMatrix3

    # Load two structures
    atoms = read_pdb(namd_pdb)
    x = position.(select(atoms, "name CA"))

    # test RMSD function
    y = x .+ Ref(SVector{3}(1, 1, 1))
    @test rmsd(x, y) ≈ sqrt(length(x) * 3 / length(x))

    # apply a random rotation and translation to x
    y = x .+ Ref(SVector{3}(45.0, -15.0, 31.5))
    y .= Ref(rand(RotMatrix3)) .* y
    @test rmsd(x, y) > 0.0
    z = align(x, y)
    @test rmsd(z, y) ≈ 0.0 atol = 1e-5

    # same test, but with coordinates obtained from a simulation frame
    simulation = Simulation(namd_pdb, namd_traj)
    first_frame!(simulation)
    cas = findall(Select("name CA"), atoms)
    x = positions(current_frame(simulation))[cas]
    y = x .+ Ref(SVector{3}(45.0, -15.0, 31.5))
    y .= Ref(rand(RotMatrix3)) .* y
    @test rmsd(x, y) > 0.0
    z = align(x, y)
    @test rmsd(z, y) ≈ 0.0 atol = 1e-10

    rmsd_aligned = zeros(length(simulation))
    rmsd_notaligned = zeros(length(simulation))
    first_frame!(simulation)
    x = positions(current_frame(simulation))[cas]
    xref = copy(x)
    for (iframe, frame) in enumerate(simulation)
        local x, z
        x = positions(frame)[cas]
        rmsd_notaligned[iframe] = rmsd(x, xref)
        z = align(x, xref)
        rmsd_aligned[iframe] = rmsd(z, xref)
    end

    @test rmsd_notaligned ≈ [0.0, 3.0627719174308323, 3.1891492625876556, 2.9879924980792314, 3.9749958688486617]
    @test rmsd_aligned ≈ [0.0, 2.8388710154609034, 2.9776998440690385, 2.4621444212469483, 3.8035683196100796]
    @test all(rmsd_aligned .<= rmsd_notaligned)

    cas = findall(Select("name CA"), atoms)
    @test rmsd(simulation, cas) ≈ rmsd_aligned
    @test rmsd(simulation, "name CA") ≈ rmsd_aligned
    @test rmsd(simulation, cas; mass=mass.(atoms[cas])) ≈ rmsd_aligned
    @test rmsd(simulation, "name CA"; mass=mass.(atoms[cas])) ≈ rmsd_aligned
    @test rmsd(simulation, cas; reference_frame=5) ≈ [3.8035683196100787, 4.680280207599843, 3.4614944346303917, 2.97835421429809, 0.0]

    # Average structure
    @test rmsd(simulation, cas; reference_frame=:average, show_progress=false) ≈ [1.8995986972454748, 2.1512244220536973, 1.5081703191869376, 1.1651111324544219, 2.757039151265317]

    # Input errors
    @test_throws "index 6 is not in frame range" rmsd(simulation, cas; reference_frame=6)
    simulation = Simulation(namd_pdb, namd_traj; step=2)
    @test_throws "index 2 is not in frame range" rmsd(simulation, cas; reference_frame=2)
    @test_throws "must be an integer" rmsd(simulation, cas; reference_frame='a')

    # Test rmsd using rmsd_of (compared to output of VMD)
    simulation = Simulation(MolSimToolkit.Testing.namd2_pdb, MolSimToolkit.Testing.namd2_traj)
    r_expected = [4.782562070489359e-17, 2.383157150961262, 2.016995295094623, 1.2936259143528597, 1.6782044794091227, 1.4561259303589944, 1.8725268913794244, 2.1419325487817344, 3.055806849195156, 2.7914291413904566, 2.3027060849001626, 4.099998830641137, 3.558595497478691, 3.412962337347324, 4.685408093153975, 3.6144080835882413, 3.807043696380331, 6.217853448234399, 4.677065083960004, 4.460039973437796]
    r = rmsd(simulation, "protein and name CA"; rmsd_of="residue 47 to 53")
    @test r ≈ r_expected

    # Test with :average reference (this test is weak, not compared to a reference)
    r = rmsd(simulation, "protein and name CA"; rmsd_of="residue 47 to 53", reference_frame=:average)
    r_expected = [3.9399443319320473, 4.033526124094138, 4.14241791220504, 3.706840579105693, 3.4744894549615526, 3.752341297472445, 3.171958741785707, 3.3091658448472425, 3.6970582271206522, 3.458804341856055, 3.355020639111263, 3.2681251061187027, 3.3584150605657452, 3.250167010371905, 3.7688706682877764, 3.0829037964839388, 3.2793987265786897, 4.930089024112391, 3.6277822605047407, 3.714044966130123]
    @test r ≈ r_expected

    # Test trying that closest atom is not first (compared to VMD)
    r = rmsd(simulation, "protein and resnum 2 to 21"; rmsd_of="protein and resnum 51 to 61")
    r_expected = [7.076072742448697e-16, 2.213844930795668, 2.277059122546429, 1.991206478988441, 2.038849169791055, 1.6148629582035468, 2.4500843431718997, 1.7215090898320655, 1.9269152419232003, 1.9522321208462048, 2.6279449804206774, 3.1972399610209674, 3.193977658189353, 2.644193744743235, 3.4094503793329403, 3.172588996098775, 3.722390992410695, 3.820539452015499, 3.3358362828365897, 3.625876855303902]
    @test r ≈ r_expected

end

"""
    rmsd_matrix(
        simulation::Simulation, 
        indices::AbstractVector{<:Integer} or sel::AbstractString; 
        mass::Union{AbstractVector{<:Integer}, Nothing} = nothing,
        align::Bool = true,
        show_progress = true,
    )

Computes the RMSD matrix for a set of atoms along a trajectory.

# Positional arguments

- `simulation`: The Simulation object. 
- `indices`: vector contains the indices of the atoms to be considered. 
or
- `sel`: A selection string, e. g. "name CA", defining the atoms to be considered.

# Optional keyword arguments

- `mass`: optional mass of the atoms if they are not the same.
- `align`: align the frames before computing the RMSD.
- `show_progress`: show or not a progress bar.

# Returns 

A symmetric matrix with the RMSD values between each pair of frames. For example, in 
a trajectory with 5 frames, the matrix will be a 5x5 matrix with the RMSD values
between the structures of each pair of frames.

# Example

```jldoctest; filter = r"([0-9]+\\.[0-9]{2})[0-9]+" => s"\\1***"
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> using PDBTools

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> rmsd_matrix(simulation, "protein and name CA"; show_progress=false)
5×5 Matrix{Float64}:
 0.0      2.83887  2.9777   2.46214  3.80357
 2.83887  0.0      2.35492  2.64463  4.68028
 2.9777   2.35492  0.0      2.08246  3.46149
 2.46214  2.64463  2.08246  0.0      2.97835
 3.80357  4.68028  3.46149  2.97835  0.0
```

"""
function rmsd_matrix(
    simulation::Simulation,
    indices::AbstractVector{<:Integer};
    mass::Union{AbstractVector{<:Integer},Nothing}=nothing,
    align::Bool=true,
    show_progress::Bool=true,
)
    if !isnothing(mass) && (length(indices) != length(mass))
        throw(ArgumentError("indices and mass vectors must have the same length"))
    end
    # Auxiliary arrays for the alignment
    xm = zeros(3, length(indices))
    xp = zeros(3, length(indices))
    # This is very memory inefficient, but it is a simple way to compute the RMSD matrix
    coordinates = [positions(frame)[indices] for frame in simulation]
    n = length(simulation)
    rmsd_matrix = zeros(n, n)
    p = Progress((n * (n - 1)) ÷ 2; enabled=show_progress)
    for iframe in 1:n
        rmsd_matrix[iframe, iframe] = 0.0
        for jframe in iframe+1:n
            next!(p)
            if align
                align!(coordinates[iframe], coordinates[jframe]; mass, xm, xp)
            end
            rmsd_matrix[iframe, jframe] = rmsd(coordinates[iframe], coordinates[jframe])
            rmsd_matrix[jframe, iframe] = rmsd_matrix[iframe, jframe]
        end
    end
    return rmsd_matrix
end

function rmsd_matrix(
    simulation::Simulation,
    sel::AbstractString;
    mass::Union{AbstractVector{<:Integer},Nothing}=nothing,
    align::Bool=true,
    show_progress::Bool=true,
)
    ats = get_atoms(simulation)
    indices = findall(PDBTools.Select(sel), ats)
    rmsd_matrix(simulation, indices; mass, align, show_progress)
end

@testitem "rmsd_matrix" begin
    using MolSimToolkit, MolSimToolkit.Testing
    using PDBTools
    atoms = readPDB(Testing.namd_pdb)
    cas = findall(Select("name CA"), atoms) # CA indices
    simulation = Simulation(Testing.namd_pdb, Testing.namd_traj)
    m = rmsd_matrix(simulation, cas)
    @test all(m .- [
        0.0 2.83887 2.9777 2.46214 3.80357
        2.83887 0.0 2.35492 2.64463 4.68028
        2.9777 2.35492 0.0 2.08246 3.46149
        2.46214 2.64463 2.08246 0.0 2.97835
        3.80357 4.68028 3.46149 2.97835 0.0
    ] .< 1e-3)

    msel = rmsd_matrix(simulation, "protein and name CA")
    @test msel ≈ m

    @test_throws ArgumentError rmsd_matrix(simulation, cas; mass=[1, 2, 3, 4, 5])
end

@testitem "reconstruct_structure!" begin
    using MolSimToolkit
    using MolSimToolkit.Testing
    using PDBTools
    using StaticArrays: SVector
    using LinearAlgebra: norm

    box = 30.0
    mat = hcat(SVector(box, 0.0, 0.0), SVector(0.0, box, 0.0), SVector(0.0, 0.0, box))

    # A single, already-contiguous chain: reconstruction must not move any atom
    # relative to its neighbors.
    chain = [SVector(1.0 * i, 0.0, 0.0) for i in 0:4]
    x = deepcopy(chain)
    reconstruct_structure!(x, collect(1:5), mat)
    @test x ≈ chain

    # A chain broken by a single periodic image jump (no real branch needed,
    # since consecutive atoms are still close once wrapped).
    x = [p + SVector(0.0, 0.0, 0.0) for p in chain]
    x[end] += SVector(3 * box, -2 * box, box)
    reconstruct_structure!(x, collect(1:5), mat)
    @test all(norm(x[i+1] - x[i]) ≈ 1.0 for i in 1:4)

    # Single branch: two chains concatenated (e.g. a dimer), stored with an
    # arbitrary periodic offset between them.
    chainA = [SVector(1.0 * i, 0.0, 0.0) for i in 0:4]
    chainB = [SVector(5.0 + 1.0 * i + 3 * box, 0.0, 0.0) for i in 0:4]
    x = Vector{SVector{3,Float64}}(vcat(chainA, chainB))
    reconstruct_structure!(x, collect(1:10), mat)
    @test all(norm(x[i+1] - x[i]) < 2.0 for i in 1:9)

    # Multiple branches: three chains concatenated, each with a different offset.
    cA = [SVector(1.0 * i, 0.0, 0.0) for i in 0:3]
    cB = [SVector(5.0 + 1.0 * i + 3 * box, 0.0, 0.0) for i in 0:3]
    cC = [SVector(10.0 + 1.0 * i - 5 * box, 2 * box, -2 * box) for i in 0:3]
    x = Vector{SVector{3,Float64}}(vcat(cA, cB, cC))
    reconstruct_structure!(x, collect(1:12), mat)
    @test all(norm(x[i+1] - x[i]) < 2.0 for i in 1:11)

    # Nested branch: within the atoms skipped over by a branch, there is a
    # further break that must, itself, be resolved by branching.
    A = [SVector(1.0 * i, 0.0, 0.0) for i in 0:3]
    B1_true, B2_true, B3_true = SVector(9.0, 10.0, 0.0), SVector(10.0, 10.0, 0.0), SVector(11.0, 10.0, 0.0)
    B4_true = SVector(4.0, 0.0, 0.0) # close to A4 = (3, 0, 0)
    B1 = B1_true + SVector(2 * box, -3 * box, 5 * box)
    B2 = B2_true + SVector(-4 * box, 6 * box, -2 * box)
    B3 = B3_true + SVector(1 * box, 1 * box, 1 * box)
    B4 = B4_true + SVector(-2 * box, 4 * box, 3 * box)
    x = Vector{SVector{3,Float64}}(vcat(A, [B1, B2, B3, B4]))
    reconstruct_structure!(x, collect(1:8), mat)
    @test norm(x[6] - x[5]) ≈ 1.0 atol = 1e-6 # B1-B2
    @test norm(x[7] - x[6]) ≈ 1.0 atol = 1e-6 # B2-B3
    @test norm(x[7] - x[5]) ≈ 2.0 atol = 1e-6 # B1-B3

    # Cross-boundary branch: an atom skipped over by one branch (Z) truly belongs,
    # spatially, with a chain that only appears later in the index sequence (C),
    # not with the chain it was jumped to (B). The search pool for filling a gap
    # must therefore include atoms from both sides of a branch point.
    A = [SVector(1.0 * i, 0.0, 0.0) for i in 0:3]
    Z_true = SVector(20.0, 0.0, 0.0)
    B_true = [SVector(4.0 + 1.0 * i, 0.0, 0.0) for i in 0:3]
    C_true = [SVector(19.0 + 1.0 * i, 0.0, 0.0) for i in 0:3]
    Z = Z_true + SVector(3 * box, -2 * box, 4 * box)
    B = [p + SVector(-3 * box, 5 * box, -box) for p in B_true]
    C = [p + SVector(2 * box, box, -4 * box) for p in C_true]
    x = Vector{SVector{3,Float64}}(vcat(A, [Z], B, C))
    reconstruct_structure!(x, collect(1:13), mat)
    @test norm(x[10] - x[5]) ≈ 1.0 atol = 1e-6 # Z-C1

    # The `dmax` keyword controls the branch search. With an overly large
    # `dmax`, the search is never triggered, so the cross-boundary case above
    # (Z-C1) degenerates to the naive, sequence-only wrap and gets it wrong.
    x = Vector{SVector{3,Float64}}(vcat(A, [Z], B, C))
    reconstruct_structure!(x, collect(1:13), mat; dmax=1000.0)
    @test !isapprox(norm(x[10] - x[5]), 1.0; atol=1e-6)

    # Conversely, an artificially tiny `dmax` forces a branch search at every
    # single step; the result must still be correct.
    x = Vector{SVector{3,Float64}}(vcat(chainA, chainB))
    reconstruct_structure!(x, collect(1:10), mat; dmax=1e-3)
    @test all(norm(x[i+1] - x[i]) ≈ 1.0 for i in 1:9)

    # A single index (or none) is a no-op.
    x = [SVector(0.0, 0.0, 0.0)]
    reconstruct_structure!(x, [1], mat)
    @test x == [SVector(0.0, 0.0, 0.0)]

    # Reusing a preallocated `aux_inds` buffer across repeated calls gives the
    # same result as the default, freshly-allocated buffer.
    aux_inds = similar(collect(1:10))
    x1 = Vector{SVector{3,Float64}}(vcat(chainA, chainB))
    x2 = deepcopy(x1)
    reconstruct_structure!(x1, collect(1:10), mat)
    reconstruct_structure!(x2, collect(1:10), mat; aux_inds)
    @test x1 == x2

    # Integration test on a real trajectory: reconstructing the whole system
    # (thousands of unrelated water/lipid atoms, so many branches are expected)
    # must not error, and does not touch atoms outside of `indices`.
    simulation = Simulation(Testing.namd2_pdb, Testing.namd2_traj)
    first_frame!(simulation)
    frame = current_frame(simulation)
    p = positions(frame)
    p_before = deepcopy(p)
    reconstruct_structure!(p, 1:length(p), unitcell(frame))
    @test length(p) == length(p_before)
    protein = findall(Select("protein and name CA"), get_atoms(simulation))
    @test all(norm(p[protein[i+1]] - p[protein[i]]) < 10.0 for i in 1:length(protein)-1)
end
