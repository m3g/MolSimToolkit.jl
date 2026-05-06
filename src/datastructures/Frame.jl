export Frame
export Point3D
export positions

"""
    Point3D{T}

A point in 3D space with coordinates `x`, `y`, and `z` of type `T`.

"""
struct Point3D{T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end

"""
    Frame

Structure that contains the data of a trajectory frame.

Public interface:

- `positions(::Frame)`: returns an array of positions of atoms in the frame.
- `unitcell(::Frame)`: returns the `UnitCell` object of the frame.

The fields of `Frame` are considered internal and should not be accessed directly. Use the provided functions instead.

"""
struct Frame{T<:Chemfiles.Frame}
    frame::T
    positions::Vector{Point3D{Float64}}
    unitcell_matrix::MMatrix{3,3,Float64,9}
end
function Frame(cf::Chemfiles.Frame)
    n = Int(Chemfiles.length(cf))
    positions = Vector{Point3D{Float64}}(undef, n)
    copyto!(positions, reinterpret(Point3D{Float64}, Chemfiles.positions(cf)))
    unitcell_matrix = MMatrix{3,3,Float64,9}(transpose(Chemfiles.matrix(Chemfiles.UnitCell(cf))))
    return Frame(cf, positions, unitcell_matrix)
end
Base.show(io::IO, f::Frame) = print(io, "Frame{Chemfiles.Frame} - first atom position: $(positions(f)[1])")

"""
    positions(frame::Frame)

Return the positions of the atoms in a `Frame` as a `Vector{Point3D{Float64}}`.

This is the default way to access the positions of the atoms in a simulation.

# Example

```jldoctest
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> positions(first_frame!(simulation))[1]
3-element Point3D{Float64} with indices SOneTo(3):
  5.912472724914551
 10.768872261047363
 28.277008056640625

```

"""
positions(f::Frame) = f.positions

"""
    unitcell(frame::Frame; tol=1e-10)

Returns the unit cell of the current frame in the trajectory. The tolerance defines the relative 
size of off-diagonal elements that can be ignored to consider the cell as orthorhombic, relative
the minimum diagonal element.

"""
function unitcell(f::Frame; tol=1e-10)
    mat = SMatrix(f.unitcell_matrix)
    if all(==(0), mat)
        @warn """\n
            Unit cell vectors are zero. The trajectory file may not contain proper unit cell information.
            Wrapping of coordinates will be disabled for current frame.
            
        """ _file = nothing _line = nothing
        valid = false
    else
        valid = true
    end
    s = abs(minimum(diag(mat))) # minimum diagonal element
    orthorhombic = all(abs(mat[i, j]) < tol * s for i in 1:3, j in 1:3 if i != j)
    return UnitCell(mat, valid, orthorhombic)
end

@testitem "unitcell" begin
    using MolSimToolkit
    using MolSimToolkit.Testing
    using StaticArrays
    import Chemfiles
    traj = Chemfiles.Trajectory(Testing.namd_traj)
    frame = Chemfiles.read(traj)
    uc = unitcell(Frame(frame))
    @test uc.valid == true
    @test uc.orthorhombic == true
    @test uc.matrix ≈ transpose(Chemfiles.matrix(Chemfiles.UnitCell(frame)))
    sim = Simulation(Testing.short_nopbc_pdb, Testing.short_nopbc_traj)
    first_frame!(sim)
    f = current_frame(sim)
    uc = unitcell(f)
    @test all(==(0), uc.matrix)
    @test uc.valid == false
    @test uc.orthorhombic == false
end

@testitem "Frame" begin
    using MolSimToolkit
    using MolSimToolkit.Testing
    using ShowMethodTesting

    sim = Simulation(Testing.short_nopbc_pdb, Testing.short_nopbc_traj)
    f = first_frame!(sim)
    @test parse_show(f) ≈ "Frame{Chemfiles.Frame} - first atom position: [-5.189000129699707, -4.394999980926514, 13.395999908447266]"
    f = next_frame!(sim)
    @test parse_show(f) ≈ "Frame{Chemfiles.Frame} - first atom position: [-6.46804666519165, -10.306520462036133, 14.48043441772461]"
    f = current_frame(sim)
    @test parse_show(f) ≈ "Frame{Chemfiles.Frame} - first atom position: [-6.46804666519165, -10.306520462036133, 14.48043441772461]"
end