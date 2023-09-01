export Trajectory
export frame_range
export index_frame
export length
export nframes
export close
export restart!
export currentframe
export nextframe!
export set_frame_range!

"""
    Trajectory(file::String; first=1, last=nothing, step=1)

Creates a new `Trajectory` object from a file. If `first`, `last`, and `step` are not specified, the
`Trajectory` will iterate over all frames in the file. 

A `Trajectory` object represents a trajectory file. It can be iterated over to
obtain the frames in the trajectory. The `Trajectory` object is a mutable struct
that contains the following data, that can be retrived by the corresponding
functions:

- `frame_range(::Trajectory)`: the range of frames to be iterated over
- `index_frame(::Trajectory)`: the index of the current frame in the trajectory
- `frame(::Trajectory)`: the current frame in the trajectory
- `length(::Trajectory)`: the number of frames in the trajectory file
- `nframes(::Trajectory)`: the number of frames to be iterated over in the trajectory file, considering the current range

The Trajectory object can also be manipulated by the following functions:

- `close(::Trajectory)`: closes the trajectory file
- `restart!(::Trajectory)`: restarts the iteration over the trajectory file
- `currentframe(::Trajectory)`: returns the current frame in the trajectory
- `nextframe!(::Trajectory)`: reads the next frame in the trajectory file and returns it. Moves the current frame to the next one.
- `set_frame_range!(::Trajectory; first, last, step)`: resets the range of frames to be iterated over. 

One important feature of the `Trajectory` object is that i can be iterated over, frame by frame.

# Examples

```jldoctest
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> trajectory = Trajectory(Testing.namd_traj, first = 2, step = 2, last = 4);

julia> for frame in trajectory 
           @show index_frame(trajectory)
           # show x coordinate of first atom 
           @show Positions(frame)[1].x
       end


```
"""
mutable struct Trajectory{R<:AbstractRange, F<:Chemfiles.Frame, T<:Chemfiles.Trajectory}
    frame_range::R
    frame::F
    index_frame::Int
    trajectory::T
end

import Base: show
function show(io::IO, trajectory::Trajectory)
    print(io, chomp("""
    Trajectory
        Trajectory file: $(path(trajectory))
        Total number of frames: $(length(trajectory))
        Frame range: $(frame_range(trajectory))
        Number of frames in range: $(nframes(trajectory))
        Current frame: $(index_frame(trajectory))
    """))
end

#=
    Trajectory(trajectory::Chemfiles.Trajectory, frame_range::AbstractRange)

Creates a new Trajectory object from a Chemfiles.Trajectory. If `frame_range` is not
specified, the Trajectory will iterate over all frames in the file. If `frame_range`
is specified, the Trajectory will iterate over the frames in the range.

This function is not supposed to be called directly. Use the Trajectory(file) function.

=#
function Trajectory(
    trajectory::Chemfiles.Trajectory, 
    frame_range::AbstractRange
)
    frame = Chemfiles.read(trajectory)
    for _ in 2:first(frame_range)
        Chemfiles.read!(trajectory, frame)
    end
    index_frame = first(frame_range)
    return Trajectory(frame_range, frame, index_frame, trajectory)
end

#= 
    Trajectory(file::String; first=1, last=nothing, step=1)

Creates a new Trajectory object from a file. If `first`, `last`, and `step` are not specified, the
Trajectory will iterate over all frames in the file. This is the default constructor expected
to be used.

=#
function Trajectory(file::String; first=1, last=nothing, step=1) 
    if isnothing(last)
        frame_range = first:step:Int(Chemfiles.length(Chemfiles.Trajectory(file)))
    else
        frame_range = first:step:last
    end
    return Trajectory(Chemfiles.Trajectory(file), frame_range)
end

"""
    frame_range(trajectory::Trajectory)

Returns the range of frames to be iterated over.

"""
frame_range(trajectory::Trajectory) = trajectory.frame_range

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
considering the current frame range.

"""
nframes(trajectory::Trajectory) = length(trajectory.frame_range)

"""
    path(trajectory::Trajectory)

Returns the path to the trajectory file.

"""
path(trajectory::Trajectory) = Chemfiles.path(trajectory.trajectory)

"""
    restart!(trajectory::Trajectory)

Restarts the iteration over the trajectory file.

"""
function restart!(trajectory::Trajectory)
    trajectory_file = path(trajectory)
    close(trajectory)
    trajectory.trajectory = Chemfiles.Trajectory(trajectory_file)
    trajectory.index_frame = first(frame_range(trajectory))
    trajectory.frame = Chemfiles.read(trajectory.trajectory)
    return trajectory
end

"""
    currentframe(trajectory::Trajectory)

Returns the current frame in the trajectory.

"""
currentframe(trajectory::Trajectory) = trajectory.frame

"""
    nextframe!(trajectory::Trajectory)

Reads the next frame in the trajectory file and returns it. Moves the current
frame to the next one in the range to be considered (given by `frame_range(trajectory)`).

"""
function nextframe!(trajectory::Trajectory) 
    if index_frame(trajectory) == last(frame_range(trajectory))
        error("End of trajectory.")
    end
    Chemfiles.read!(trajectory.trajectory, trajectory.frame)
    iframe = index_frame(trajectory) + 1
    while iframe ∉ frame_range(trajectory) && iframe < last(frame_range(trajectory)) 
        Chemfiles.read!(trajectory.trajectory, trajectory.frame)
        iframe += 1
    end
    # If the last frame was reached, check if it is in the frame range
    if iframe ∉ frame_range(trajectory)
        error("End of trajectory.")
    end
    trajectory.index_frame = iframe
    return currentframe(trajectory)
end

import Base: iterate
function iterate(trajectory::Trajectory, iframe=nothing)
    if isnothing(iframe)
        restart!(trajectory)
        return (currentframe(trajectory), index_frame(trajectory))
    elseif iframe < last(frame_range(trajectory))
        nextframe!(trajectory)
        return (currentframe(trajectory), index_frame(trajectory))
    else
        return nothing
    end
end

"""
    set_frame_range!(trajectory::Trajectory; first=1, last=nothing, step=1)

Resets the frame range to be iterated over. This function will restart the
iteration from the first frame of the new range.

"""
function set_frame_range!(trajectory::Trajectory; first=1, last=nothing, step=1)
    if isnothing(last)
        frame_range = first:step:length(trajectory)
    else
        frame_range = first:step:last
    end
    trajectory.frame_range = frame_range
    trajectory.index_frame = first(frame_range)
    restart!(trajectory)
end

@testitem "Trajectory" begin
    import Chemfiles
    using MolSimToolkit.Testing


    #
    # Test iterator by reading coordinates
    #

    # Read coordinates with Chemfiles interface
    t = Chemfiles.Trajectory(Testing.namd_traj)
    c = zeros(Chemfiles.length(t))
    f = Chemfiles.read(t)
    c[1] = Chemfiles.positions(f)[1,1]
    for i in 2:Chemfiles.length(t)
        Chemfiles.read!(t, f)
        c[i] = Chemfiles.positions(f)[1,1]
    end
    Chemfiles.close(t)

    # Read coordinates with the Trajectory interface
    trajectory = Trajectory(Testing.namd_traj)
    c2 = zeros(length(trajectory))
    i = 0
    for frame in trajectory
        i += 1
        c2[i] = Positions(frame)[1].x
    end
    close(trajectory)

    @test c == c2

end

