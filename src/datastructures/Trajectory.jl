export Trajectory
export length
export close

mutable struct Trajectory{F<:Chemfiles.Frame, T<:Chemfiles.Trajectory}
    frame::F
    trajectory::T
end
Trajectory(trajectory::Chemfiles.Trajectory) = Trajectory(Chemfiles.read(trajectory), trajectory)
Trajectory(file::String) = Trajectory(Chemfiles.Trajectory(file))

import Base: close, length
close(trajectory::Trajectory) = Chemfiles.close(trajectory.trajectory)
length(trajectory::Trajectory) = Int(Chemfiles.length(trajectory.trajectory))

import Base: iterate
function iterate(trajectory::Trajectory, state=nothing)
    frame = trajectory.frame
    if isnothing(state)
        return (frame, 1)
    elseif state < length(trajectory)
        Chemfiles.read!(trajectory.trajectory, frame)
        return (frame, state + 1)
    else
        #voltar
        Chemfiles.
        return nothing
    end
end

