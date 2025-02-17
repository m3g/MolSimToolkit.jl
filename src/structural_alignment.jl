"""
    rmsd(simulation::Simulation, indices::AbstractVector{Int}; mass = nothing, reference_frame = nothing, show_progress = true)

Computes the root mean square deviation (RMSD) between two sets of points in along a trajectory.

# Arguments

- `indices` vector contains the indices of the atoms to be considered. 
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

julia> atoms = readPDB(Testing.namd_pdb);

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> cas = findall(sel"name CA", atoms); # CA indices

julia> rmsd(simulation, cas; show_progress=false)
5-element Vector{Float64}:
 0.0
 2.8388710154609034
 2.9776998440690385
 2.4621444212469483
 3.8035683196100796

julia> rmsd(simulation, cas; reference_frame=:average, show_progress=false)
5-element Vector{Float64}:
 1.8995986972454748
 2.1512244220536973
 1.5081703191869376
 1.1651111324544219
 2.757039151265317
```

"""
function rmsd(
    simulation::Simulation, indices::AbstractVector{Int};
    mass=nothing,
    reference_frame=nothing,
    show_progress=true,
)

    # Auxiliary arrays for the alignment
    xm = zeros(3, length(indices))
    xp = zeros(3, length(indices))

    # Define reference of the alignment
    xref = if isnothing(reference_frame)
        first_frame!(simulation)
        positions(current_frame(simulation))[indices]
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
        positions(current_frame(simulation))[indices]
    elseif reference_frame == :average
        xref = fill(zero(Point3D), length(indices))
        p = Progress(length(simulation); enabled=show_progress, desc="Computing average structure:")
        for (iframe, frame) in enumerate(simulation)
            next!(p)
            x = @view(positions(frame)[indices])
            align!(x, xref; mass, xm, xp)
            @. xref = (xref * (iframe - 1) + x) / iframe
        end
        xref
    else
        throw(ArgumentError("""\n

            reference_frame must be an integer, nothing or :average
            
        """))
    end
    restart!(simulation)
    rmsds = Float64[]
    p = Progress(length(simulation); enabled=show_progress, desc="Computing RMSDs for each frame:")
    for frame in simulation
        next!(p)
        x = @view(positions(frame)[indices])
        align!(x, xref; mass, xm, xp)
        push!(rmsds, rmsd(x, xref))
    end
    return rmsds
end

@testitem "rmsd" begin
    using MolSimToolkit
    using MolSimToolkit.Testing: namd_pdb, namd_traj
    using StaticArrays: SVector
    using PDBTools
    using Rotations: RotMatrix3

    # Load two structures
    atoms = readPDB(namd_pdb)
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
    @test rmsd(simulation, cas; mass=mass.(atoms[cas])) ≈ rmsd_aligned
    @test rmsd(simulation, cas; reference_frame=5) ≈ [3.8035683196100787, 4.680280207599843, 3.4614944346303917, 2.97835421429809, 0.0]

    # Average structure
    @test rmsd(simulation, cas; reference_frame=:average, show_progress=false) ≈ [1.8995986972454748, 2.1512244220536973, 1.5081703191869376, 1.1651111324544219, 2.757039151265317]

    # Input errors
    @test_throws ArgumentError rmsd(simulation, cas; reference_frame=6)
    simulation = Simulation(namd_pdb, namd_traj; step=2)
    @test_throws ArgumentError rmsd(simulation, cas; reference_frame=2)

end

"""
    rmsd_matrix(
        simulation::Simulation, 
        indices::AbstractVector{Int}; 
        mass::Union{AbstractVector{Int}, Nothing} = nothing,
        align::Bool = true,
        show_progress = true,
    )

Computes the RMSD matrix for a set of atoms along a trajectory.

The `indices` vector contains the indices of the atoms to be considered. 
The `mass` argument can be used to provide the mass of the atoms if they are not the same.
The `align` argument can be used to align the frames before computing the RMSD.

The `show_progress` argument can be used to show a progress bar.

# Returns 

A symetric matrix with the RMSD values between each pair of frames. For example, in 
a trajectory with 5 frames, the matrix will be a 5x5 matrix with the RMSD values
between the structures of each pair of frames.

# Example

```jldoctest; filter = r"([0-9]+\\.[0-9]{2})[0-9]+" => s"\\1***"
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> using PDBTools

julia> atoms = readPDB(Testing.namd_pdb);

julia> cas = findall(Select("name CA"), atoms); # CA indices

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> rmsd_matrix(simulation, cas; show_progress=false)
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
    indices::AbstractVector{Int};
    mass::Union{AbstractVector{Int},Nothing}=nothing,
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

    @test_throws ArgumentError rmsd_matrix(simulation, cas; mass=[1, 2, 3, 4, 5])
end
