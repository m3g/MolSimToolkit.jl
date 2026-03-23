"""
    Frame

Structure that contains the data of a trajectory frame.

Current methods:

- `positions(::Frame)`: returns an array of positions of atoms in the frame.
- `unitcell(::Frame)`: returns the `UnitCell` object of the frame.

"""
struct Frame{T<:Chemfiles.Frame}
    frame::T
end
Base.show(io::IO, f::Frame) = print(io, "Frame{Chemfiles.Frame} - first atom position: $(positions(f)[1])")

"""
    unitcell(frame::Frame; tol=1e-10)

Returns the unit cell of the current frame in the trajectory. The tolerance defines the relative 
size of off-diagonal elements that can be ignored to consider the cell as orthorhombic, relative
the minimum diagonal element.

"""
function unitcell(f::Frame; tol=1e-10) 
    mat = SMatrix{3,3,Float64,9}(transpose((Chemfiles.matrix(Chemfiles.UnitCell(f.frame)))))
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