"""
    Trajectory

Structure that contains the data stream of trajectory file. 

All methods dispatching on `Trajectory` are considered internal.

"""
struct Trajectory{T<:Chemfiles.Trajectory}
    trajectory::T
end
Trajectory(filename::String) = Trajectory(Chemfiles.Trajectory(filename))
Base.length(t::Trajectory) = Int(Chemfiles.length(t.trajectory))
Base.read(t::Trajectory) = Frame(Chemfiles.read(t.trajectory))
Base.close(t::Trajectory) = Chemfiles.close(t.trajectory)
path_trajectory(t::Trajectory) = normpath(Chemfiles.path(t.trajectory))
raw_length(t::Trajectory) = Int(Chemfiles.length(t.trajectory))
function read_step!(t::Trajectory, i_next_frame::Int, f::Frame)
    Chemfiles.read_step!(t.trajectory, i_next_frame - 1, f.frame)
    copyto!(f.positions, reinterpret(Point3D{Float64}, Chemfiles.positions(f.frame)))
    f.unitcell_matrix .= transpose(Chemfiles.matrix(Chemfiles.UnitCell(f.frame)))
    return f
end
