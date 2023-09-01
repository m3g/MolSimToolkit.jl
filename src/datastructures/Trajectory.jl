export Trajectory
export length, nframes
export close
export setrange!
export nextframe!

"""
    Trajectory

A Trajectory object represents a trajectory file. It can be iterated over to
obtain the frames in the trajectory. The Trajectory object is a mutable struct
that contains the following data, that can be retrived by the corresponding
functions:

- `range(::Trajectory)`: the range of frames to be iterated over
- `index_frame(::Trajectory)`: the index of the current frame in the trajectory
- `frame(::Trajectory)`: the current frame in the trajectory
- `length(::Trajectory)`: the number of frames in the trajectory file
- `nframes(::Trajectory)`: the number of frames to be iterated over in the trajectory file, considering the current range

The Trajectory object can also be manipulated by the following functions:

- `close(::Trajectory)`: closes the trajectory file
- `restart!(::Trajectory)`: restarts the iteration over the trajectory file
- `nextframe!(::Trajectory)`: reads the next frame in the trajectory file and returns it. Moves the current frame to the next one.
- `setrange!(::Trajectory; first, last, step)`: resets the range of frames to be iterated over. 

"""
mutable struct Trajectory{R<:AbstractRange, F<:Chemfiles.Frame, T<:Chemfiles.Trajectory}
    range::R
    frame::F
    index_frame::Int
    trajectory::T
end

#=
    Trajectory(trajectory::Chemfiles.Trajectory, range::AbstractRange)

Creates a new Trajectory object from a Chemfiles.Trajectory. If `range` is not
specified, the Trajectory will iterate over all frames in the file. If `range`
is specified, the Trajectory will iterate over the frames in the range.

This function is not supposed to be called directly. Use the Trajectory(file) function.

=#
function Trajectory(
    trajectory::Chemfiles.Trajectory, 
    range::AbstractRange
)
    frame = Chemfiles.read(trajectory)
    for _ in 2:first(range)
        Chemfiles.read!(trajectory, frame)
    end
    index_frame = first(range)
    Trajectory(range, frame, index_frame, trajectory)
end

"""
    Trajectory(file::String; first=1, last=nothing, step=1)

Creates a new Trajectory object from a file. If `first`, `last`, and `step` are not specified, the
Trajectory will iterate over all frames in the file. 

# Examples

```julia
traj = Trajectory("traj.dcd"; first=1, step=5)
```

```julia
traj = Trajectory("traj.dcd"; first=1, step=5, last=100)
```

```julia
traj = Trajectory("traj.dcd"; first=1, last=100)
```

"""
function Trajectory(file::String; first=1, last=nothing, step=1) 
    if isnothing(last)
        range = first:step:Int(Chemfiles.length(Chemfiles.Trajectory(file)))
    else
        range = first:step:last
    end
    Trajectory(Chemfiles.Trajectory(file), range=range)
end

# Overloading Base.range here is sort of odd. Probably this interface will change.
"""
    range(trajectory::Trajectory)

Returns the range of frames to be iterated over.

"""
Base.range(trajectory::Trajectory) = trajectory.range

"""
    index_frame(trajectory::Trajectory)

Returns the index of the current frame in the trajectory.

"""
index_frame(trajectory::Trajectory) = trajectory.index_frame

"""
    frame(trajectory::Trajectory)

Returns the current frame in the trajectory.

"""
frame(trajectory::Trajectory) = trajectory.frame

"""
    close(trajectory::Trajectory)

Closes the trajectory file.

"""
Base.close(trajectory::Trajectory) = Chemfiles.close(trajectory.trajectory)

"""
    length(trajectory::Trajectory)

Returns the number of frames in the trajectory file.

"""
Base.length(trajectory::Trajectory) = Int(Chemfiles.length(trajectory.trajectory))

"""
    nframes(trajectory::Trajectory)

Returns the number of frames to be iterated over in the trajectory file,
considering the current range.

"""
nframes(trajectory::Trajectory) = length(trajectory.range)

"""
    path(trajectory::Trajectory)

Returns the path to the trajectory file.

"""
path(trajectory::Trajectory) = trajectory.trajectory.path

"""
    restart!(trajectory::Trajectory)

Restarts the iteration over the trajectory file.

"""
function restart!(trajectory::Trajectory)
    close(trajectory)
    trajectory.trajectory = Chemfiles.Trajectory(path(trajectory))
    return trajectory
end

"""
    nextframe!(trajectory::Trajectory)

Reads the next frame in the trajectory file and returns it. Moves the current
frame to the next one.

"""
function nextframe!(trajectory::Trajectory) 
    if index_frame(trajectory) < last(range(trajectory))
        trajectory.index_frame += 1
        Chemfiles.read!(trajectory.trajectory, trajectory.frame)
    else
        error("End of trajectory")
    end
    return frame(trajectory)
end

import Base: iterate
function iterate(trajectory::Trajectory, state=nothing)
    if isnothing(state)
        restart!(trajectory)
        frame = frame(trajectory)
        for _ in 1:first(trajectory.range)-1
           frame = nextframe!(trajectory) 
        end
        return (frame, index_frame(trajectory))
    elseif state <= last(trajectory.range)

        Chemfiles.read!(trajectory.trajectory, frame)
        return (frame, state + 1)
    else
        #voltar
#        Chemfiles.
#        return nothing
    end
end

"""
    setrange!(trajectory::Trajectory; first=1, last=nothing, step=1)

Resets the range of frames to be iterated over. This function will restart the
iteration from the first frame of the new range.

"""
function setrange!(trajectory::Trajectory; first=1, last=nothing, step=1)
    if isnothing(last)
        range = first:step:length(trajectory)
    else
        range = first:step:last
    end
    trajectory.range = range
    trajectory.index_frame = first(range)
    restart!(trajectory)
end

@testitem "Trajectory" begin
    import Chemfiles
    using MolSimToolkit.Testing

    # Read the first atom coordinate of each frame, to test iterations
    traj_cm = Chemfiles.Trajectory(Testing.namd_traj)
    cm_first_coordinates = zeros(Chemfiles.length(traj))
    frame_cm = Chemfiles.read(traj_cm)
    cm_first_coordinates[1] = Chemfiles.positions(frame_cm)[1,1]
    for i in 2:Chemfiles.length(traj)
        Chemfiles.read!(traj_cm, frame_cm)
        cm_first_coordinates[i] = Chemfiles.positions(frame_cm)[1,1]
    end
    Chemfiles.close(traj_cm)

    # Test the Trajectory iteration
    traj = Trajectory(Testing.namd_traj)






end

