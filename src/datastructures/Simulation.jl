export Simulation
export frame_range
export frame_index
export length
export raw_length
export close
export restart!
export first_frame! 
export current_frame
export next_frame!
export set_frame_range!
export atoms
export unitcell
export path_pdb
export path_trajectory
export get_frame


"""
    Simulation(pdb_file::String, trajectory_file::String; frames=[1,2,3,5])
    Simulation(pdb_file::String, trajectory_file::String; frames=9:2:20)
    Simulation(pdb_file::String, trajectory_file::String; first=1, last=nothing, step=1)
    Simulation(atoms::AbstractVector{<:AtomType}, trajectory_file::String; first=1, last=nothing, step=1) 

Creates a new `Simulation` object. 

The first constructor creates a `Simulation` object from a PDB or mmCIF file and a trajectory file. It will use the
`PDBTools.Atom` for the atom type, which will populate the `atoms` vector of the `Simulation` object.
Currently, other atom types are supported, if the `MolSimToolkit.atomic_mass(::AtomType)` function is defined
for the atom type.

With the second constructor, the `atoms` vector is passed as an argument. This is useful when the atoms
are provided by a different source than the PDB file. 

The `frames` or `first`, `last`, and `step` arguments can be used to specify the frames to be iterated over:

    - `frames` can be a vector of frame indices, e. g., `frames=[1,2,3,5]` or `frames=9:2:20`.
    - `first`, `last`, and `step` are Integers that specify the frames to be iterated over. 
      If `last` is not specified, the last frame in the trajectory will be used.

A `Simulation` object contains a trajectory file and a PDB data of the atoms. It can be iterated over to
obtain the frames in the trajectory. The `Simulation` object is a mutable struct
that contains the following data, that can be retrieved by the corresponding
functions:

- `frame_range(::Simulation)`: the list of frames to be iterated over
- `frame_index(::Simulation)`: the index of the current frame in the trajectory
- `length(::Simulation)`: the number of frames to be iterated over in the trajectory file, considering the current range
- `raw_length(::Simulation)`: the number of frames in the trajectory file
- `atoms(::Simulation)`: the atoms in the simulation

The Simulation object can also be manipulated by the following functions:

- `close(::Simulation)`: closes the trajectory file
- `restart!(::Simulation)`: restarts the iteration over the trajectory file
- `first_frame!(::Simulation)`: restarts the iteration over the trajectory file and places the current frame at the first frame in the trajectory
- `current_frame(::Simulation)`: returns the current frame in the trajectory
- `next_frame!(::Simulation)`: reads the next frame in the trajectory file and returns it. Moves the current frame to the next one.
- `set_frame_range!(::Simulation; first, last, step)`: resets the range of frames to be iterated over. 
- `get_frame(::Simulation, iframe)`: returns the frame at the given index in the trajectory.

One important feature of the `Simulation` object is that it can be iterated over, frame by frame. 

The `pairs` iterator can also be used to iterate over the frames, returning a tuple with the frame index
and the frame itself. 

The `enumerate` iterator can also be used to iterate over the frames, returning
a tuple with the frame counter and the frame itself.

# Examples

```julia-repl # to be doctest
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> simulation = Simulation(
           Testing.namd_pdb, Testing.namd_traj; 
           first = 2, step = 2, last = 4,
           # or frames = [2,4]
       );

julia> for frame in simulation 
           @show frame_index(simulation)
           # show x coordinate of first atom 
           @show positions(frame)[1].x
       end
frame_index(simulation) = 2
((positions(frame))[1]).x = 5.912472724914551
frame_index(simulation) = 4
((positions(frame))[1]).x = 7.346549034118652

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
    AtomType,
    V<:Vector{<:AtomType},
    F<:Chemfiles.Frame,
    T<:Chemfiles.Trajectory,
    L<:ReentrantLock
}
    pdb_file::Union{Nothing,String}
    atoms::V
    frame_range::Vector{Int}
    frame::F
    frame_index::Union{Nothing,Int}
    trajectory::T
    read_lock::L
end

@views function _print_frame_range(simulation) 
    if length(simulation) <= 6 
        join(frame_range(simulation), ", ")
    else
        join(frame_range(simulation)[1:3], ", ") * ", ... , " * join(frame_range(simulation)[end-2:end], ", ")
    end
end

import Base: show
function show(io::IO, simulation::Simulation)
    print(io, chomp("""
    Simulation 
        Atom type: $(eltype(simulation.atoms))
        PDB file: $(isnothing(simulation.pdb_file) ? "-" : path_pdb(simulation))
        Trajectory file: $(path_trajectory(simulation))
        Total number of frames: $(raw_length(simulation))
        Frames to consider: $(_print_frame_range(simulation))
        Number of frames to consider: $(length(simulation))
        Current frame: $(frame_index(simulation))
    """))
end

import Base: lock
lock(f::F, simulation::Simulation) where {F<:Function} = lock(f, simulation.read_lock)

function _set_range(trajectory, frames, first, last, step)
    if isnothing(frames)
        isnothing(first) && (first = 1)
        isnothing(step) && (step = 1)
        isnothing(last) && (last = Int(Chemfiles.length(trajectory)))
        frame_range = first:step:last
    else
        if any(!isnothing, (first, last, step))
            throw(ArgumentError("""\n
                The `frames` argument cannot be used with `first`, `last`, or `step` arguments. 

            """))
        end
        frame_range = frames
    end
    if frame_range isa AbstractRange
        frame_range = collect(frame_range)
    end
    if !issorted(frame_range)
        sort!(frame_range)
    end
    return frame_range
end

"""
    set_frame_range!(simulation::Simulation; first=1, last=nothing, step=1)

Resets the frame range to be iterated over. This function will restart the
iteration of the simulation trajectory.

"""
function set_frame_range!(simulation::Simulation; frames=nothing, first=nothing, last=nothing, step=nothing)
    simulation.frame_range = _set_range(simulation.trajectory, frames, first, last, step)
    restart!(simulation)
end

#=

Creates a new Simulation object from a Chemfiles.Trajectory. If `frame_range` is not
specified, the Simulation will iterate over all frames in the file. If `frame_range`
is specified, the Simulation will iterate over the frames in the range.

This function is not supposed to be called directly. Use the Simulation(file) function.

=#
function Simulation(
    pdb_file::Union{Nothing,String},
    atoms::AbstractVector{AtomType},
    trajectory::Chemfiles.Trajectory,
    frames, first, last, step
) where {AtomType}
    frame_range = _set_range(trajectory, frames, first, last, step)
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
function Simulation(
    pdb_file::String, trajectory_file::String; 
    frames=nothing, first=nothing, last=nothing, step=nothing
)
    atoms = PDBTools.readPDB(pdb_file)
    return Simulation(pdb_file, atoms, Chemfiles.Trajectory(trajectory_file), frames, first, last, step)
end

function Simulation(
    atoms::AbstractVector{AtomType}, trajectory_file::String; 
    frames=nothing, first=nothing, last=nothing, step=nothing
) where {AtomType}
    return Simulation(nothing, atoms, Chemfiles.Trajectory(trajectory_file), frames, first, last, step)
end

"""
    frame_range(simulation::Simulation)

Returns the list of frames to be iterated over.

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
path_pdb(simulation::Simulation) = isnothing(simulation.pdb_file) ? nothing : normpath(simulation.pdb_file)

"""
    path_trajectory(simulation::Simulation)

Returns the path to the trajectory file of the simulation.

"""
path_trajectory(simulation::Simulation) = normpath(Chemfiles.path(simulation.trajectory))

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
function current_frame(simulation::Simulation) 
    if isnothing(frame_index(simulation))
        throw(ArgumentError("From current_frame: No frame has been read yet."))
    end
    return simulation.frame
end

"""
    first_frame!(simulation::Simulation)

Restarts the trajectory buffer, and places the current frame at the first frame in the trajectory.

# Example

```julia-repl
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj; first=3);

julia> first_frame!(simulation) 
Simulation 
    Atom type: Atom{Nothing}
    PDB file: /test/data/namd/protein_in_popc_membrane/structure.pdb
    Trajectory file: /test/data/namd/protein_in_popc_membrane/trajectory.dcd
    Total number of frames: 5
    Frames to consider: 1, 2, 3, 4, 5
    Number of frames to consider: 5
    Current frame: 3

```
"""
function first_frame!(simulation::Simulation)
    restart!(simulation)
    next_frame!(simulation)
    return simulation
end

"""
    next_frame!(simulation::Simulation)

Reads the next frame in the trajectory file and returns it. Moves the current
frame to the next one in the range to be considered (given by `frame_range(simulation)`).

"""
function next_frame!(simulation::Simulation)
    lock(simulation) do
        i_frame_in_range = if isnothing(frame_index(simulation)) 
                0
        else
            searchsortedfirst(frame_range(simulation), frame_index(simulation))
        end
        if i_frame_in_range + 1 > length(frame_range(simulation))
            throw(ArgumentError("""\n
                Next frame out of range.
                Current frame: $(frame_index(simulation)) is the last in selected frames $(_print_frame_range(simulation)).
                
            """))
        end
        i_next_frame = frame_range(simulation)[i_frame_in_range + 1]
        Chemfiles.read_step!(simulation.trajectory, i_next_frame - 1, simulation.frame)
        simulation.frame_index = i_next_frame
    end
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

import Base: eachindex
eachindex(simulation::Simulation) = frame_range(simulation)

import Base: iterate
function iterate(simulation::Simulation, iframe=nothing)
    if iframe == last(frame_range(simulation))
        return nothing
    end
    if isnothing(iframe)
        restart!(simulation)
    end
    next_frame!(simulation)
    return (current_frame(simulation), frame_index(simulation))
end

"""
    get_frame(simulation::Simulation, iframe::Integer)

Returns the frame at the given index in the trajectory. 

## Example

```julia-repl
julia> using MolSimToolkit, MolSimToolkit.Testing, PDBTools

julia> sim = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> frame4 = get_frame(sim, 4)
   Array{Atoms,1} with 20465 atoms with fields:
   index name resname chain   resnum  residue        x        y        z occup  beta model segname index_pdb
       1    N     ILE     P      211        1   -0.397   12.048   37.441  1.00  0.00     1    PROT         1
       2  HT1     ILE     P      211        1   -0.779   11.123   37.726  1.00  0.00     1    PROT         2
       3  HT2     ILE     P      211        1   -0.393   12.662   38.280  1.00  0.00     1    PROT         3
                                                       ⋮ 
   20463  SOD     SOD     S       13     4374  -11.686   23.749   19.935  1.00  0.00     1     SOD     20463
   20464  SOD     SOD     S       14     4375  -34.214   38.148   55.179  1.00  0.00     1     SOD     20464
   20465  SOD     SOD     S       15     4376    7.220  -52.702   66.223  1.00  0.00     1     SOD     20465

julia> writePDB(frame4, "frame4.pdb")
```

!!! note
    The `get_frame` function will read the frames in the trajectory until the desired frame is reached. 
    This can be slow for large trajectories. If the required frame is before the current frame of the 
    simulation, the simulation will be restarted. The simulation object is returned positioned in the
    required frame. 

"""
function get_frame(simulation::Simulation, iframe::Integer)
    if !(iframe in frame_range(simulation))
        throw(ArgumentError("get_frame: Index $iframe selected frames: $(frame_range(simulation))."))
    end
    i_current_frame = frame_index(simulation)
    if isnothing(i_current_frame) || (iframe < i_current_frame)
        first_frame!(simulation)
    end
    while frame_index(simulation) != iframe
        next_frame!(simulation)
    end
    p = positions(current_frame(simulation))
    ats = atoms(simulation)
    for iat in eachindex(ats, p)
        ats[iat].x = p[iat].x
        ats[iat].y = p[iat].y
        ats[iat].z = p[iat].z
    end
    return ats
end

@testitem "Simulation" begin
    import Chemfiles
    using MolSimToolkit.Testing

    #
    # Test iterator by reading coordinates
    #

    # Read coordinates with Chemfiles interface
    t = Chemfiles.Trajectory(Testing.namd_traj)
    c = zeros(Chemfiles.length(t))
    f = Chemfiles.read(t)
    c = [Chemfiles.positions(f)[1, 1]]
    for i in 2:Chemfiles.length(t)
        Chemfiles.read!(t, f)
        push!(c, Chemfiles.positions(f)[1, 1])
    end
    Chemfiles.close(t)

    # Read coordinates with the Simulation interface
    simulation = Simulation(Testing.namd_pdb, Testing.namd_traj)
    c2 = Float64[]
    for frame in simulation
        push!(c2, positions(frame)[1].x)
    end
    close(simulation)
    @test c == c2
end

@testitem "no-pdb construct" begin
    import PDBTools
    using MolSimToolkit.Testing
    pdb = PDBTools.readPDB(Testing.namd_pdb)
    simulation = Simulation(pdb, Testing.namd_traj)
    first_frame!(simulation)
    @test positions(current_frame(simulation))[1].x == 5.912472724914551
    @test isnothing(path_pdb(simulation))
end

@testitem "get_frame" begin
    using MolSimToolkit, MolSimToolkit.Testing, PDBTools
    sim = Simulation(Testing.namd_pdb, Testing.namd_traj)
    frames = [ copy(positions(frame)) for frame in sim ]
    @test all(coor(get_frame(sim, i)) ≈ frames[i] for i in eachindex(sim))
    @test coor(get_frame(sim, 5)) ≈ frames[5]
    @test coor(get_frame(sim, 1)) ≈ frames[1]
    @test coor(get_frame(sim, 2)) ≈ frames[2]
    ats = readPDB(Testing.namd_pdb)
    sim = Simulation(ats, Testing.namd_traj)
    @test all(coor(get_frame(sim, i)) ≈ frames[i] for i in eachindex(sim))
    @test_throws ArgumentError get_frame(sim, 100)
    sim2 = Simulation(Testing.namd_pdb, Testing.namd_traj; frames=[1,2,5])
    @test coor(get_frame(sim2, 1)) ≈ frames[1]
    @test coor(get_frame(sim2, 2)) ≈ frames[2]
    @test coor(get_frame(sim2, 5)) ≈ frames[5]
    @test_throws ArgumentError get_frame(sim2, 4)
    first_frame!(sim2)
    @test frame_index(sim2) == 1
    next_frame!(sim2)
    @test frame_index(sim2) == 2
    next_frame!(sim2)
    @test frame_index(sim2) == 5
    @test_throws ArgumentError next_frame!(sim2)
    set_frame_range!(sim2; frames=1:3)
    @test frame_range(sim2) == [1,2,3]
    set_frame_range!(sim2; first=2, last=4)
    @test frame_range(sim2) == [2,3,4]
end

#
# Legacy interface (to be removed in 2.0)
#
const nextframe! = next_frame!
const firstframe! = first_frame!
export nextframe!, firstframe!
