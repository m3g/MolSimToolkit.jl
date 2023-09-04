export Simulation
export frame_range
export frame_index
export length
export raw_length
export close
export restart!
export current_frame
export nextframe!
export set_frame_range!
export atoms
export unitcell
export path_pdb
export path_trajectory

"""
    Simulation(pdb_file::String, trajectory_file::String; first=1, last=nothing, step=1)

Creates a new `Simulation` object from a file. If `first`, `last`, and `step` are not specified, the
`Simulation` will iterate over all frames in the file. 

A `Simulation` object contains a trajectory file and a PDB data of the atoms. It can be iterated over to
obtain the frames in the trajectory. The `Simulation` object is a mutable struct
that contains the following data, that can be retrived by the corresponding
functions:

- `frame_range(::Simulation)`: the range of frames to be iterated over
- `frame_index(::Simulation)`: the index of the current frame in the trajectory
- `length(::Simulation)`: the number of frames to be iterated over in the trajectory file, considering the current range
- `raw_length(::Simulation)`: the number of frames in the trajectory file
- `atoms(::Simulation)`: the atoms in the simulation

The Simulation object can also be manipulated by the following functions:

- `close(::Simulation)`: closes the trajectory file
- `restart!(::Simulation)`: restarts the iteration over the trajectory file
- `current_frame(::Simulation)`: returns the current frame in the trajectory
- `nextframe!(::Simulation)`: reads the next frame in the trajectory file and returns it. Moves the current frame to the next one.
- `set_frame_range!(::Simulation; first, last, step)`: resets the range of frames to be iterated over. 

One important feature of the `Simulation` object is that it can be iterated over, frame by frame. 

The `pairs` iterator can also be used to iterate over the frames, returning a tuple with the frame index
and the frame itself. 

The `enumerate` iterator can also be used to iterate over the frames, returning
a tuple with the frame counter and the frame itself.

# Examples

```julia-repl # to be doctest
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj, first = 2, step = 2, last = 4);

julia> for frame in simulation 
           @show frame_index(simulation)
           # show x coordinate of first atom 
           @show Positions(frame)[1].x
       end
frame_index(simulation) = 2
((Positions(frame))[1]).x = 5.912472724914551
frame_index(simulation) = 4
((Positions(frame))[1]).x = 7.346549034118652

julia> for (i, frame) in pairs(simulation)
           @show i, frame_index(simulation)
       end
(i, frame_index(simulation)) = (2, 2)
(i, frame_index(simulation)) = (4, 4)  

julia> for (i, frame) in enumerate(simulation)
           @show i, frame_index(simulation)
       end
(i, frame_index(simulation)) = (1, 2)
(i, frame_index(simulation)) = (2, 4)

```
"""
mutable struct Simulation{
    V<:Vector{PDBTools.Atom}, 
    R<:AbstractRange, 
    F<:Chemfiles.Frame, 
    T<:Chemfiles.Trajectory,
    L<:ReentrantLock
}
    pdb_file::String
    atoms::V 
    frame_range::R
    frame::F
    frame_index::Union{Nothing,Int}
    trajectory::T
    read_lock::L
end

import Base: show
function show(io::IO, simulation::Simulation)
    print(io, chomp("""
    Simulation
        PDB file: $(path_pdb(simulation))
        Simulation file: $(path_trajectory(simulation))
        Total number of frames: $(length(simulation))
        Frame range: $(frame_range(simulation))
        Number of frames in range: $(length(simulation))
        Current frame: $(frame_index(simulation))
    """))
end

import Base: lock
lock(f::F, simulation::Simulation) where {F<:Function} = lock(f, simulation.read_lock)

#=
    Simulation(
        pdb_file::String,
        atoms::AbstractVector{PDBTools.Atom},
        trajectory::Chemfiles.Trajectory, frame_range::AbstractRange)
    )

Creates a new Simulation object from a Chemfiles.Trajectory. If `frame_range` is not
specified, the Simulation will iterate over all frames in the file. If `frame_range`
is specified, the Simulation will iterate over the frames in the range.

This function is not supposed to be called directly. Use the Simulation(file) function.

=#
function Simulation(
    pdb_file::String,
    atoms::AbstractVector{PDBTools.Atom},
    trajectory::Chemfiles.Trajectory, 
    frame_range::AbstractRange
)
    frame = Chemfiles.read(trajectory)
    read_lock = ReentrantLock()
    frame_index = nothing
    simulation = Simulation(pdb_file, atoms, frame_range, frame, frame_index, trajectory, read_lock)
    restart!(simulation)
    return simulation
end

#= 
    Simulation(pdb_file::String, trajectory_file::String; first=1, last=nothing, step=1)

Creates a new `Simulation` object. If `first`, `last`, and `step` are not specified, the
Simulation will iterate over all frames in the file. This is the default constructor expected
to be used.

=#
function Simulation(pdb_file::String, trajectory_file::String; first=1, last=nothing, step=1) 
    if isnothing(last)
        frame_range = first:step:Int(Chemfiles.length(Chemfiles.Trajectory(trajectory_file)))
    else
        frame_range = first:step:last
    end
    atoms = PDBTools.readPDB(pdb_file)
    return Simulation(pdb_file, atoms, Chemfiles.Trajectory(trajectory_file), frame_range)
end

"""
    frame_range(simulation::Simulation)

Returns the range of frames to be iterated over.

"""
frame_range(simulation::Simulation) = simulation.frame_range

"""
    frame_index(simulation::Simulation)

Returns the index of the current frame in the trajectory. Returns `nothing` 
if no frame frame from the trajectory range has been read yet.

"""
frame_index(simulation::Simulation) = simulation.frame_index

"""
    close(simulation::Simulation)

Closes the trajectory file.

"""
Base.close(simulation::Simulation) = Chemfiles.close(simulation.trajectory)

"""
    raw_length(simulation::Simulation)

Returns the number of frames in the trajectory file.

"""
raw_length(simulation::Simulation) = Int(Chemfiles.length(simulation.trajectory))

"""
    length(simulation::Simulation)

Returns the number of frames to be iterated over in the trajectory file,
considering the current frame range.

"""
Base.length(simulation::Simulation) = length(simulation.frame_range)

"""
    atoms(simulation::Simulation)

Returns the atoms in the simulation.

"""
atoms(simulation::Simulation) = simulation.atoms

"""
    path_pdb(simulation::Simulation)

Returns the path to the pdb file of the simulation.

"""
path_pdb(simulation::Simulation) = simulation.pdb_file

"""
    path_trajectory(simulation::Simulation)

Returns the path to the trajectory file of the simulation.

"""
path_trajectory(simulation::Simulation) = Chemfiles.path(simulation.trajectory)

"""
    restart!(simulation::Simulation)

Restarts the iteration over the trajectory file.

"""
function restart!(simulation::Simulation)
    lock(simulation) do 
        trajectory_file = path_trajectory(simulation)
        close(simulation)
        simulation.trajectory = Chemfiles.Trajectory(trajectory_file)
        simulation.frame_index = nothing
    end # release lock
    return simulation
end

"""
    current_frame(simulation::Simulation)

Returns the current frame in the trajectory.

"""
current_frame(simulation::Simulation) = simulation.frame

"""
    nextframe!(simulation::Simulation)

Reads the next frame in the trajectory file and returns it. Moves the current
frame to the next one in the range to be considered (given by `frame_range(simulation)`).

"""
function nextframe!(simulation::Simulation) 
    lock(simulation) do 
        if frame_index(simulation) == last(frame_range(simulation))
            error("End of trajectory.")
        end
        Chemfiles.read!(simulation.trajectory, simulation.frame)
        if isnothing(frame_index(simulation))
            iframe = first(frame_range(simulation))
        else
            iframe = frame_index(simulation) + 1
        end
        while iframe ∉ frame_range(simulation) && iframe < last(frame_range(simulation)) 
            Chemfiles.read!(simulation.trajectory, current_frame(simulation))
            iframe += 1
        end
        # If the last frame was reached, check if it is in the frame range
        if iframe ∉ frame_range(simulation)
            error("End of trajectory.")
        end
        simulation.frame_index = iframe
    end # release lock
    return current_frame(simulation)
end

"""
    unitcell(frame::Chemfiles.Frame)

Returns the unit cell of the current frame in the trajectory.

"""
unitcell(f::Chemfiles.Frame) = unitcell(Chemfiles.UnitCell(f))
unitcell(u::Chemfiles.UnitCell) = SMatrix{3,3,Float64,9}(transpose(Chemfiles.matrix(u)))

@testitem "unitcell" begin
    using MolSimToolkit
    using MolSimToolkit.Testing
    using StaticArrays
    import Chemfiles
    traj = Chemfiles.Trajectory(Testing.namd_traj)
    frame = Chemfiles.read(traj)
    u = unitcell(frame)
    @test u ≈ transpose(Chemfiles.matrix(Chemfiles.UnitCell(frame)))
end

import Base: firstindex, lastindex
firstindex(simulation::Simulation) = first(frame_range(simulation))
lastindex(simulation::Simulation) = last(frame_range(simulation))

import Base: keys
keys(simulation::Simulation) = frame_range(simulation)

import Base:eachindex 
eachindex(simulation::Simulation) = frame_range(simulation)

import Base: iterate
function iterate(simulation::Simulation, iframe=nothing)
    if iframe == last(frame_range(simulation))
        return nothing
    end
    if isnothing(iframe)
        restart!(simulation)
    end
    nextframe!(simulation)
    return (current_frame(simulation), frame_index(simulation))
end

"""
    set_frame_range!(simulation::Simulation; first=1, last=nothing, step=1)

Resets the frame range to be iterated over. This function will restart the
iteration of the simulation trajectory.

"""
function set_frame_range!(simulation::Simulation; first=1, last=nothing, step=1)
    if isnothing(last)
        frame_range = first:step:length(simulation)
    else
        frame_range = first:step:last
    end
    simulation.frame_range = frame_range
    restart!(simulation)
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
    c = [Chemfiles.positions(f)[1,1]]
    for i in 2:Chemfiles.length(t)
        Chemfiles.read!(t, f)
        push!(c, Chemfiles.positions(f)[1,1])
    end
    Chemfiles.close(t)

    # Read coordinates with the Simulation interface
    simulation = Simulation(Testing.namd_pdb, Testing.namd_traj)
    c2 = Float64[]
    for frame in simulation
        push!(c2, Positions(frame)[1].x)
    end
    close(simulation)

    @test c == c2

end

