using LinearAlgebra: eigvecs
using StaticArrays: SVector, MMatrix, SMatrix

_center_of_mass(x::AbstractVector{<:AbstractVector}, ::Nothing) = sum(x) / length(x)
_center_of_mass(x::AbstractVector{<:AbstractVector}, mass::AbstractVector) =
    sum(x[i] * mass[i] for i in eachindex(x, mass)) / sum(mass)

"""
    align(x, y; mass = nothing)
    align!(x, y; mass = nothing)

Aligns two structures (sets of points in 3D space). Solves
the "Procrustes" problem, which is to find the best
translation, and rotation, that aligns the two
structures, minimizing the RMSD between them.

Structures are expected to be of the same size, and the 
correspondence is assumed from the vector indices. 

`align` returns a new vector containing the coordinates of x aligned to y. 
`align!` modifies the input vector `x` in place.

"""
function align end
@doc (@doc align) align!

function align(
    x::AbstractVector{<:AbstractVector},
    y::AbstractVector{<:AbstractVector};
    mass=nothing
)
    xnew = copy(x)
    return align!(xnew, y; mass)
end

function align!(
    x::AbstractVector{<:AbstractVector},
    y::AbstractVector{<:AbstractVector};
    mass=nothing,
    # Auxiliary arrays that might be preallocated
    xm=zeros(3, length(x)),
    xp=zeros(3, length(x))
)
    length(x) == length(y) || throw(DimensionMismatch("x and y must have the same length"))
    (length(x[1]) != 3 || length(x[2]) != 3) && throw(DimensionMismatch("x and y must be 3D vectors"))

    cmx = _center_of_mass(x, mass)
    cmy = _center_of_mass(y, mass)
    for i in eachindex(x, y)
        x[i] -= cmx
        y[i] -= cmy
    end

    for i in eachindex(x, y)
        xm[1:3, i] .= y[i] .- x[i]
        xp[1:3, i] .= y[i] .+ x[i]
    end

    q = zeros(MMatrix{4,4,eltype(xm),16})
    for i in eachindex(x)
        q[1, 1] = q[1, 1] + sum(abs2, @view(xm[1:3, i]))
        q[1, 2] = q[1, 2] + xp[2, i] * xm[3, i] - xm[2, i] * xp[3, i]
        q[1, 3] = q[1, 3] + xm[1, i] * xp[3, i] - xp[1, i] * xm[3, i]
        q[1, 4] = q[1, 4] + xp[1, i] * xm[2, i] - xm[1, i] * xp[2, i]
        q[2, 2] = q[2, 2] + xp[2, i]^2 + xp[3, i]^2 + xm[1, i]^2
        q[2, 3] = q[2, 3] + xm[1, i] * xm[2, i] - xp[1, i] * xp[2, i]
        q[2, 4] = q[2, 4] + xm[1, i] * xm[3, i] - xp[1, i] * xp[3, i]
        q[3, 3] = q[3, 3] + xp[1, i]^2 + xp[3, i]^2 + xm[2, i]^2
        q[3, 4] = q[3, 4] + xm[2, i] * xm[3, i] - xp[2, i] * xp[3, i]
        q[4, 4] = q[4, 4] + xp[1, i]^2 + xp[2, i]^2 + xm[3, i]^2
    end
    q[2, 1] = q[1, 2]
    q[3, 1] = q[1, 3]
    q[3, 2] = q[2, 3]
    q[4, 1] = q[1, 4]
    q[4, 2] = q[2, 4]
    q[4, 3] = q[3, 4]
    q = SMatrix(q)

    # Computing the eigenvectors 'v' of the q matrix
    v = eigvecs(q)

    # Compute rotation matrix
    u = zeros(MMatrix{3,3,Float64,9})
    u[1, 1] = v[1, 1]^2 + v[2, 1]^2 - v[3, 1]^2 - v[4, 1]^2
    u[1, 2] = 2.0 * (v[2, 1] * v[3, 1] + v[1, 1] * v[4, 1])
    u[1, 3] = 2.0 * (v[2, 1] * v[4, 1] - v[1, 1] * v[3, 1])
    u[2, 1] = 2.0 * (v[2, 1] * v[3, 1] - v[1, 1] * v[4, 1])
    u[2, 2] = v[1, 1]^2 + v[3, 1]^2 - v[2, 1]^2 - v[4, 1]^2
    u[2, 3] = 2.0 * (v[3, 1] * v[4, 1] + v[1, 1] * v[2, 1])
    u[3, 1] = 2.0 * (v[2, 1] * v[4, 1] + v[1, 1] * v[3, 1])
    u[3, 2] = 2.0 * (v[3, 1] * v[4, 1] - v[1, 1] * v[2, 1])
    u[3, 3] = v[1, 1]^2 + v[4, 1]^2 - v[2, 1]^2 - v[3, 1]^2
    u = SMatrix(u)

    # Rotate to align x to y 
    for i in eachindex(x)
        x[i] = u * x[i]
    end

    # Move aligned x to the original center of mass of y
    for i in eachindex(x, y)
        x[i] += cmy
        y[i] += cmy
    end

    return x
end

"""
    rmsd(x::AbstractVector,y::AbstractVector)
    rmsd(simulation::Simulation, indices::AbstractVector{Int}; mass = nothing, reference_frame = nothing, show_progress = true)

Computes the root mean square deviation (RMSD) between two sets of points in 
space, or along a trajectory.

The `rmsd(x,y)` method computes the RMSD between two sets of points `x` and `y`. 
The sets are expected to be of the same size, and the correspondence is assumed from the vector indices.

The `rmsd(simulation::Simulation, indices)` method computes the RMSD along a trajectory. With the 
following options: 
- The `indices` vector contains the indices of the atoms to be considered. 
- The `mass` argument can be used to provide the mass of the atoms if they are not the same.
- The `reference_frame` argument can be used to provide a reference frame to align the trajectory to:
    - If `reference_frame == nothing`, the first frame will be used (default behavior).
    - If `reference_frame == :average`, the average structure will be used.
    - If `reference_frame` is an integer, the frame at that index will be used as reference.

# Examples

## If the set is compared toi tself, the RMSD should be zero:

```jldoctest
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> using PDBTools

julia> ca = coor(readPDB(Testing.namd_pdb), "name CA");

julia> rmsd(ca, ca)
0.0
```

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
function rmsd(x::AbstractVector, y::AbstractVector)
    rmsd = 0.0
    for i in eachindex(x, y)
        rmsd += sum(abs2, x[i] .- y[i])
    end
    return sqrt(rmsd / length(x))
end

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

@testitem "procrustes" begin
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



