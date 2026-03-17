#
# Function to reconstruct the structure of the protein/polymer, to avoid
# periodic boundary conditions issues with molecule breaking through the
# boundaries. This function assumes that atoms that are close in the 
# sequence of the structure are also close in space.
#
function _reconstruct_structure!(
    x::AbstractVector{<:AbstractVector},
    indices::AbstractVector{<:Integer},
    unitcell::AbstractMatrix
)
    xlast = x[first(indices)]
    for i in indices
        x[i] = wrap(x[i], xlast, unitcell)
        xlast = x[i]
    end
    return x
end

#
# This function reconstructs the structure of the complex formed by the two sets of indices,
# to avoid periodic boundary conditions issues with molecule breaking through the
# boundaries. This function assumes that atoms that are close in the sequence of the structure 
# are also close in space.
# 
# The function reconstruct a complex, by finding first the closest atom between the two structures,
# and wrapping the coordinates of the second structure around that closest atom.
#
function _reconstruct_complex!(
    x::AbstractVector{<:AbstractVector},
    indices1::AbstractVector{<:Integer},
    indices2::AbstractVector{<:Integer},
    unitcell::AbstractMatrix
)
    # Find atom in indices2 that is closer to indices1, considering PBCs
    # If an very close atom is found, use it and avoid the full double loop
    id, jd, d = 0, 0, +Inf
    for i in eachindex(indices1)
        x1 = x[indices1[i]]
        for j in eachindex(indices2)
            x2 = wrap(x[indices2[j]], x1, unitcell)
            dsq = sum(abs2, x2 - x1)
            if dsq < d
                id = i
                jd = j
                d = dsq
                d < 4.0 && break
            end
        end
        d < 4.0 && break
    end
    # Reconstruct structure defined by indices2, starting from
    # atom jd, and moving backwards and forward
    x[indices2[jd]] = wrap(x[indices2[jd]], x[indices1[id]], unitcell)
    xlast = x[indices2[jd]]
    for j in jd-1:-1:firstindex(indices2)
        x[indices2[j]] = wrap(x[indices2[j]], xlast, unitcell)
        xlast = x[indices2[j]]
    end
    xlast = x[indices2[jd]]
    for j in jd+1:lastindex(indices2)
        x[indices2[j]] = wrap(x[indices2[j]], xlast, unitcell)
        xlast = x[indices2[j]]
    end
    return x
end

#
# Function to read the reference coordinates from the trajectory.
# returns the coordinates, reconstructured by _reconstruct_structure!,
# and _reconstruct_complex!
#
function _reference_coordinates(
    simulation,
    reference_frame,
    indices,
    rmsd_indices,
    mass;
    show_progress,
)
    xalign, xrmsd = if isnothing(reference_frame)
        first_frame!(simulation)
        frame = current_frame(simulation)
        p, uc = positions(frame), unitcell(frame)
        _reconstruct_structure!(p, indices, uc.matrix)
        _reconstruct_complex!(p, indices, rmsd_indices, uc.matrix)
        p[indices], p[rmsd_indices]
    elseif reference_frame isa Integer
        if !(reference_frame in frame_range(simulation))
            throw(ArgumentError("""\n

                reference_frame index $reference_frame is not in frame range of the simulation: $(frame_range(simulation))

            """))
        end
        restart!(simulation)
        for _ in 1:reference_frame
            next_frame!(simulation)
        end
        frame = current_frame(simulation)
        p, uc = positions(frame), unitcell(frame)
        _reconstruct_structure!(p, indices, uc.matrix)
        _reconstruct_complex!(p, indices, rmsd_indices, uc.matrix)
        p[indices], p[rmsd_indices]
    elseif reference_frame == :average
        xm = zeros(3, length(indices))
        xp = zeros(3, length(indices))
        xalign_ref = zeros(Point3D, length(indices))
        xrmsd_ref = zeros(Point3D, length(indices))
        prg = Progress(length(simulation); enabled=show_progress, desc="Computing average structure:")
        for (iframe, frame) in enumerate(simulation)
            next!(prg)
            p, uc = positions(frame), unitcell(frame)
            _reconstruct_structure!(p, indices, uc.matrix)
            _reconstruct_complex!(p, indices, rmsd_indices, uc.matrix)
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
                    @. xrmsd_ref = (xrmsd_ref * (iframe - 1) + xalign) / iframe
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

    rmsds = Float64[]
    prg = Progress(length(simulation); enabled=show_progress, desc="Computing RMSDs for each frame:")
    for frame in simulation
        next!(prg)
        p, uc = positions(frame), unitcell(frame)
        _reconstruct_structure!(p, indices, uc.matrix)
        _reconstruct_complex!(p, indices, rmsd_of, uc.matrix)
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
    indices = findall(PDBTools.Select(sel), atoms(simulation))
    rmsd_indices = if rmsd_of === sel
        indices
    else
        findall(PDBTools.Select(rmsd_of), atoms(simulation))
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
    x = coor(atoms, "name CA")

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

    # Test rmsd using rmsd_of
    simulation = Simulation(MolSimToolkit.Testing.namd2_pdb, MolSimToolkit.Testing.namd2_traj)
    r_expected = [4.782562070489359e-17, 2.383157150961262, 2.016995295094623, 1.2936259143528597, 1.6782044794091227, 1.4561259303589944, 1.8725268913794244, 2.1419325487817344, 3.055806849195156, 2.7914291413904566, 2.3027060849001626, 4.099998830641137, 3.558595497478691, 3.412962337347324, 4.685408093153975, 3.6144080835882413, 3.807043696380331, 6.217853448234399, 4.677065083960004, 4.460039973437796]
    r = rmsd(simulation, "protein and name CA"; rmsd_of="residue 47 to 53")
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
    ats = atoms(simulation)
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
