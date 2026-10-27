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
- `unitcell(::Frame)`: returns the `UnitCell` object of the frame (or `nothing`, if the system is not periodic).

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
julia> using MolSimToolkit, MolSimToolkit.Testing, LazyArtifacts

julia> simulation = Simulation(test_file("popc_membrane.pdb"), test_file("popc_membrane.dcd"));

julia> positions(first_frame!(simulation))[1]
3-element Point3D{Float64} with indices SOneTo(3):
  5.912472724914551
 10.768872261047363
 28.277008056640625

```

"""
positions(f::Frame) = f.positions

"""
    unitcell(frame::Frame)

Returns the unit cell of the frame, as a `UnitCell` object (the unit cell type of PDBTools),
or `nothing` if the frame does not have a unit cell (the system is not periodic).

If the `Simulation` was created with a fixed unit cell (with the `unitcell` keyword), that unit
cell is returned for every frame.

"""
function unitcell(f::Frame)
    if all(==(0), f.unitcell_matrix)
        @warn """\n
            Unit cell vectors are zero. The trajectory file may not contain proper unit cell information.
            The system will be considered non-periodic. A fixed unit cell can be set with the
            `unitcell` keyword of the `Simulation` constructor.

        """ _file = nothing _line = nothing maxlog = 1
        return nothing
    end
    return UnitCell(f.unitcell_matrix)
end

# Wrapping of coordinates, ignored if the frame does not have a unit cell
_wrap(x, _, ::Nothing) = x
_wrap(x, xref, uc::UnitCell) = wrap(x, xref, uc)

@testitem "unitcell" begin
    using MolSimToolkit
    using MolSimToolkit.Testing, LazyArtifacts
    using StaticArrays
    import Chemfiles
    traj = Chemfiles.Trajectory(test_file("popc_membrane.dcd"))
    frame = Chemfiles.read(traj)
    uc = unitcell(Frame(frame))
    @test uc isa UnitCell
    @test uc ≈ transpose(Chemfiles.matrix(Chemfiles.UnitCell(frame)))
    sim = Simulation(test_file("short_nopbc.pdb"), test_file("short_nopbc.dcd"))
    first_frame!(sim)
    f = current_frame(sim)
    @test isnothing(unitcell(f))
end

@testitem "Frame" begin
    using MolSimToolkit
    using MolSimToolkit.Testing, LazyArtifacts
    using ShowMethodTesting

    sim = Simulation(test_file("short_nopbc.pdb"), test_file("short_nopbc.dcd"))
    f = first_frame!(sim)
    @test parse_show(f) ≈ "Frame{Chemfiles.Frame} - first atom position: [-5.189000129699707, -4.394999980926514, 13.395999908447266]"
    f = next_frame!(sim)
    @test parse_show(f) ≈ "Frame{Chemfiles.Frame} - first atom position: [-6.46804666519165, -10.306520462036133, 14.48043441772461]"
    f = current_frame(sim)
    @test parse_show(f) ≈ "Frame{Chemfiles.Frame} - first atom position: [-6.46804666519165, -10.306520462036133, 14.48043441772461]"
end