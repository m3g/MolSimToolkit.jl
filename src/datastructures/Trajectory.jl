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
    return f.frame
end
