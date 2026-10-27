export Simulation
export Trajectory
export Frame
export UnitCell
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
export unitcell
export path_pdb
export path_trajectory
export get_frame!
export goto_frame!
export get_atoms

"""
    Simulation(pdb_file::String, trajectory_file::String; frames=[1,2,3,5], unitcell=nothing)
    Simulation(pdb_file::String, trajectory_file::String; frames=9:2:20, unitcell=nothing)
    Simulation(pdb_file::String, trajectory_file::String; first=1, last=nothing, step=1, unitcell=nothing)
    Simulation(atoms::AbstractVector{<:AtomType}, trajectory_file::String; first=1, last=nothing, step=1, unitcell=nothing)
    Simulation(pdb_file::String, trajectory_file1::String, trajectory_file2::String, ...; kargs...)
    Simulation(pdb_file::String, trajectory_files::AbstractVector{<:AbstractString}; kargs...)

Creates a new `Simulation` object. 

The first constructor creates a `Simulation` object from a PDB or mmCIF file and a trajectory file. It will use the
`PDBTools.Atom` for the atom type, which will populate the `atoms` vector of the `Simulation` object.
Currently, other atom types are supported, if the `MolSimToolkit.atomic_mass(::AtomType)` function is defined
for the atom type.

With the second constructor, the `atoms` vector is passed as an argument. This is useful when the atoms
are provided by a different source than the PDB file. 

More than one trajectory file can be provided, either as additional arguments or as a vector of file names
(with either the `pdb_file` or the `atoms` vector as the first argument). The files are considered sequential
parts of the same simulation: frames are numbered continuously across files, in the order the files
were given, and all analyses can be performed on the concatenated trajectory. All files must contain
the same number of atoms.

The `frames` or `first`, `last`, and `step` arguments can be used to specify the frames to be iterated over:

- `frames` can be a vector of frame indices, e. g., `frames=[1,2,3,5]` or `frames=9:2:20`.
- `first`, `last`, and `step` are Integers that specify the frames to be iterated over. 
   If `last` is not specified, the last frame in the trajectory will be used.

By default, the unit cell of each frame is read from the trajectory file. The `unitcell` keyword
sets a fixed unit cell, which will be used for all frames, ignoring the unit cell of the trajectory
(useful for trajectories without unit cell information, of constant-volume simulations):

- `unitcell=true`: use the unit cell of the PDB file (or of the `atoms`, if they are a `PDBTools.Structure` with a unit cell).
- `unitcell=uc`: use the unit cell `uc`, a `UnitCell` object, which can be created from the sides of
   an orthorhombic box, the matrix of lattice vectors, or the lattice parameters (see the
   [`UnitCell` documentation of PDBTools](https://m3g.github.io/PDBTools.jl/stable/)). Other values
   are converted with `UnitCell(uc)`, e. g. `unitcell=[10.0, 10.0, 10.0]`.

A `Simulation` object contains a trajectory file and a PDB data of the atoms. It can be iterated over to
obtain the frames in the trajectory. The `Simulation` object is a mutable struct
that contains the following data, that can be retrieved by the corresponding
functions:

- `frame_range(::Simulation)`: the list of frames to be iterated over
- `frame_index(::Simulation)`: the index of the current frame in the trajectory
- `length(::Simulation)`: the number of frames to be iterated over in the trajectory file, considering the current range
- `raw_length(::Simulation)`: the number of frames in the trajectory file
- `get_atoms(::Simulation)`: get the atoms of the simulation, with the positions and unit cell of the current frame, as a `PDBTools.Structure`.
- `unitcell(::Frame)`: the unit cell of a frame (`unitcell(current_frame(simulation))` for the current frame).

The Simulation object can also be manipulated by the following functions:

- `close(::Simulation)`: closes the trajectory file
- `restart!(::Simulation)`: restarts the iteration over the trajectory file, current frame is the first frame.
- `first_frame!(::Simulation)`: restarts the iteration over the trajectory file and places the current frame at the first frame in the trajectory
- `current_frame(::Simulation)`: returns the current frame in the trajectory
- `next_frame!(::Simulation)`: reads the next frame in the trajectory file and returns it. Moves the current frame to the next one.
- `set_frame_range!(::Simulation; first, last, step)`: resets the range of frames to be iterated over. 
- `get_frame!(::Simulation, iframe)`: returns the frame at the given index in the trajectory.
- `goto_frame!(::Simulation, iframe)`: similar to `get_frame!`, but returns the `Simulation` object placed in the desired frame.

One important feature of the `Simulation` object is that it can be iterated over, frame by frame. 

The `pairs` iterator can also be used to iterate over the frames, returning a tuple with the frame index
and the frame itself. 

The `enumerate` iterator can also be used to iterate over the frames, returning a tuple with the frame counter and the frame itself.

# Examples

```julia-repl # to be doctest
julia> using MolSimToolkit, MolSimToolkit.Testing, LazyArtifacts

julia> simulation = Simulation(
           test_file("popc_membrane.pdb"), test_file("popc_membrane.dcd"); 
           first = 2, step = 2, last = 4,
           # or frames = [2,4]
       );

julia> for frame in simulation 
           @show frame_index(simulation)
           # show x coordinate of first atom 
           @show positions(frame)[1].x
       end
frame_index(simulation) = 2
((positions(frame))[1]).x = 2.001203775405884
frame_index(simulation) = 4
((positions(frame))[1]).x = -0.39684221148490906 

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
    V<:AbstractVector{<:AtomType},
    F<:Frame,
    T<:Trajectory,
    L<:ReentrantLock
}
    pdb_file::Union{Nothing,String}
    atoms::V
    frame_range::Vector{Int}
    frame::F
    frame_index::Int
    trajectory::T
    unitcell::Union{Nothing,UnitCell}
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
function show(io::IO, sim::Simulation)
    compact = get(io, :compact, false)::Bool
    if compact
        isnothing(path_pdb(sim)) ? pdb = "-" : pdb = basename(path_pdb(sim))
        print(io, "Simulation($pdb, $(join(basename.(_trajectory_files(sim)), ", ")))")
    else
        print(io, chomp("""
        Simulation 
            Atom type: $(eltype(sim.atoms))
            PDB file: $(isnothing(sim.pdb_file) ? "-" : path_pdb(sim))
            $(_print_trajectory_files(sim))
            Total number of frames: $(raw_length(sim))
            Frames to consider: $(_print_frame_range(sim))
            Number of frames to consider: $(length(sim))
            Current frame: $(frame_index(sim))
        """))
        isnothing(sim.unitcell) || print(io, "\n    Fixed unit cell: $(sim.unitcell)")
    end
end

_trajectory_files(sim::Simulation) = vcat(path_trajectory(sim))

function _print_trajectory_files(sim::Simulation)
    files = _trajectory_files(sim)
    length(files) == 1 && return "Trajectory file: $(only(files))"
    return "Trajectory files:\n" * join(("        " * file for file in files), "\n")
end

import Base: lock
lock(f::F, simulation::Simulation) where {F<:Function} = lock(f, simulation.read_lock)

function _set_range(trajectory::Trajectory, frames, first, last, step)
    if isnothing(frames)
        isnothing(first) && (first = 1)
        isnothing(step) && (step = 1)
        isnothing(last) && (last = length(trajectory))
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

Resets the frame range to be iterated over. This function will restart the iteration of the simulation trajectory.

"""
function set_frame_range!(simulation::Simulation; frames=nothing, first=nothing, last=nothing, step=nothing)
    lock(simulation) do
        simulation.frame_range = _set_range(simulation.trajectory, frames, first, last, step)
        restart!(simulation)
    end
    return simulation
end

#=

Creates a new Simulation object from a Trajectory. If `frame_range` is not
specified, the Simulation will iterate over all frames in the file. If `frame_range`
is specified, the Simulation will iterate over the frames in the range.

This function is not supposed to be called directly. Use the Simulation(file) function.

=#
function Simulation(
    pdb_file::Union{Nothing,String},
    atoms::AbstractVector{AtomType},
    trajectory::Trajectory,
    frames, first, last, step, unitcell
) where {AtomType}
    frame_range = _set_range(trajectory, frames, first, last, step)
    frame = read(trajectory)
    read_lock = ReentrantLock()
    frame_index = 1
    uc = _fixed_unitcell(atoms, unitcell)
    simulation = Simulation(pdb_file, atoms, frame_range, frame, frame_index, trajectory, uc, read_lock)
    restart!(simulation)
    return simulation
end

# Fixed unit cell of the simulation, from the `unitcell` keyword of the constructor
_fixed_unitcell(_, ::Nothing) = nothing
function _fixed_unitcell(atoms, unitcell::Bool)
    unitcell || return nothing
    uc = atoms isa Structure ? atoms.unitcell : nothing
    if isnothing(uc)
        throw(ArgumentError("""\n
            unitcell=true requires the atoms to have a unit cell (e. g. read from the CRYST1 field of a PDB file).
            Provide the unit cell with `unitcell=...` instead.

        """))
    end
    return uc
end
_fixed_unitcell(_, unitcell) = UnitCell(unitcell)

# Replaces the unit cell read from the trajectory by the fixed unit cell, if any
function _set_fixed_unitcell!(simulation::Simulation)
    uc = simulation.unitcell
    isnothing(uc) || (simulation.frame.unitcell_matrix .= uc)
    return simulation
end

#= 
    Simulation(pdb_file::String, trajectory_file::String; first=1, last=nothing, step=1)

Creates a new `Simulation` object. If `first`, `last`, and `step` are not specified, the
Simulation will iterate over all frames in the file. This is the default constructor expected
to be used.

=#
function Simulation(
    pdb_file::String, trajectory_files::AbstractVector{<:AbstractString}; 
    frames=nothing, first=nothing, last=nothing, step=nothing, unitcell=nothing
)
    atoms = PDBTools.read_pdb(pdb_file)
    return Simulation(pdb_file, atoms, Trajectory(trajectory_files), frames, first, last, step, unitcell)
end

function Simulation(
    atoms::AbstractVector{AtomType}, trajectory_files::AbstractVector{<:AbstractString};
    frames=nothing, first=nothing, last=nothing, step=nothing, unitcell=nothing
) where {AtomType}
    return Simulation(nothing, atoms, Trajectory(trajectory_files), frames, first, last, step, unitcell)
end

Simulation(pdb_file::String, trajectory_file::String, more_trajectory_files::String...; kargs...) = 
    Simulation(pdb_file, [trajectory_file, more_trajectory_files...]; kargs...)
Simulation(atoms::AbstractVector, trajectory_file::String, more_trajectory_files::String...; kargs...) = 
    Simulation(atoms, [trajectory_file, more_trajectory_files...]; kargs...)

"""
    frame_range(simulation::Simulation)

Returns the list of frames to be iterated over.

"""
frame_range(simulation::Simulation) = simulation.frame_range

"""
    frame_index(simulation::Simulation)

Returns the index of the current frame in the trajectory. 

"""
frame_index(simulation::Simulation) = simulation.frame_index

"""
    close(simulation::Simulation)

Closes the trajectory file.

"""
Base.close(simulation::Simulation) = close(simulation.trajectory)

"""
    raw_length(simulation::Simulation)

Returns the number of frames in the trajectory file (or the total number of frames, if
the simulation was built from multiple trajectory files).

"""
raw_length(simulation::Simulation) = raw_length(simulation.trajectory)

"""
    length(simulation::Simulation)

Returns the number of frames to be iterated over in the trajectory file,
considering the current frame range.

"""
Base.length(simulation::Simulation) = length(simulation.frame_range)

"""
    get_atoms(simulation::Simulation)

Returns the atoms of the simulation, with the positions of the current frame, as a `PDBTools.Structure`
carrying the unit cell of the current frame (`nothing` if the system is not periodic). 

The unit cell is thus written by `PDBTools.write_pdb`, and used by the PDBTools functions that 
support periodic boundary conditions with `unitcell=true`.

"""
function get_atoms(simulation::Simulation)
    ats = copy.(simulation.atoms)
    PDBTools.set_position!.(ats, positions(current_frame(simulation)))
    return Structure(ats; unitcell=unitcell(current_frame(simulation)))
end

@testitem "get_atoms" begin
    using MolSimToolkit
    import PDBTools
    using MolSimToolkit.Testing, LazyArtifacts
    simulation = Simulation(test_file("popc_membrane.pdb"), test_file("popc_membrane.dcd"); first = 2, step = 2, last = 4)
    ats = get_atoms(goto_frame!(simulation, 2))
    @test position(ats[1]) ≈ [2.0012038, 7.035769, 43.272514]
    ats = get_atoms(goto_frame!(simulation, 4))
    @test position(ats[1]) ≈ [-0.3968422, 12.047887, 37.44127] 
    ats = get_atoms(simulation)
    @test position(ats[1]) ≈ [-0.3968422, 12.047887, 37.44127] 
    # The structure carries the unit cell of the current frame
    @test ats isa PDBTools.Structure
    @test ats.unitcell ≈ unitcell(current_frame(simulation))
    sim = Simulation(test_file("short_nopbc.pdb"), test_file("short_nopbc.dcd"))
    @test isnothing(get_atoms(sim).unitcell)
end

@testitem "fixed unit cell" begin
    using MolSimToolkit
    using MolSimToolkit.Testing, LazyArtifacts
    using PDBTools
    using ShowMethodTesting: parse_show
    # Unit cell of the PDB file
    sim = Simulation(test_file("popc_membrane.pdb"), test_file("popc_membrane.dcd"); unitcell=true)
    uc_pdb = read_pdb(test_file("popc_membrane.pdb")).unitcell
    @test all(unitcell(frame) == uc_pdb for frame in sim)
    @test get_atoms(goto_frame!(sim, 3)).unitcell == uc_pdb
    @test parse_show(sim; repl=Dict(r"(PDB file:).*" => s"\1", r"(Trajectory file:).*" => s"\1", r"PDBTools(Base)?\." => "")) ≈ """
    Simulation 
    Atom type: Atom{Nothing}
    PDB file:
    Trajectory file:
    Total number of frames: 5
    Frames to consider: 1, 2, 3, 4, 5
    Number of frames to consider: 5
    Current frame: 3
    Fixed unit cell: UnitCell(a=47.467, b=47.467, c=87.903, α=90.00, β=90.00, γ=90.00)
    """
    # The unit cell of the trajectory is used by default
    sim = Simulation(test_file("popc_membrane.pdb"), test_file("popc_membrane.dcd"))
    @test !(unitcell(first_frame!(sim)) ≈ uc_pdb)
    # Unit cell given explicitly, for a trajectory without unit cell
    sim = Simulation(test_file("short_nopbc.pdb"), test_file("short_nopbc.dcd"); unitcell=[100.0, 100.0, 100.0])
    @test all(unitcell(frame) == UnitCell(100.0, 100.0, 100.0) for frame in sim)
    sim = Simulation(read_pdb(test_file("short_nopbc.pdb")), test_file("short_nopbc.dcd"); unitcell=UnitCell(100.0, 100.0, 100.0))
    @test unitcell(first_frame!(sim)) == UnitCell(100.0, 100.0, 100.0)
    # unitcell=true requires a unit cell in the atoms
    ats = read_pdb(test_file("popc_membrane.pdb"))
    ats.unitcell = nothing
    @test_throws ArgumentError Simulation(ats, test_file("popc_membrane.dcd"); unitcell=true)
    @test_throws ArgumentError Simulation(collect(ats), test_file("popc_membrane.dcd"); unitcell=true)
end

"""
    path_pdb(simulation::Simulation)

Returns the path to the pdb file of the simulation.

"""
path_pdb(simulation::Simulation) = isnothing(simulation.pdb_file) ? nothing : normpath(simulation.pdb_file)

"""
    path_trajectory(simulation::Simulation)

Returns the path to the trajectory file of the simulation. If the simulation was built from
multiple trajectory files, returns a vector with the paths to all files.

"""
path_trajectory(simulation::Simulation) = path_trajectory(simulation.trajectory)

"""
    restart!(simulation::Simulation)

Restarts the iteration over the trajectory file, placing the simulation at the first frame.
Similar to `first_frame!`, but returning the `Simulation` object.

"""
function restart!(simulation::Simulation)
    lock(simulation) do
        first_frame!(simulation)
        simulation.frame_index = first(simulation.frame_range)
    end # release lock
    return simulation
end

"""
    current_frame(simulation::Simulation)

Returns the current frame in the trajectory.

"""
current_frame(simulation::Simulation) = simulation.frame

"""
    first_frame!(simulation::Simulation)

Restarts the trajectory buffer, and places the current frame at the first frame in the simulation frame range.
Similar to `restart!`, but returning the `Frame` object.

# Example

```jldoctest
julia> using MolSimToolkit, MolSimToolkit.Testing, LazyArtifacts

julia> simulation = Simulation(test_file("popc_membrane.pdb"), test_file("popc_membrane.dcd"); first=3); # note first=3

julia> first_frame!(simulation) 
Frame{Chemfiles.Frame} - first atom position: [7.346549034118652, 2.396179437637329, 40.05739974975586]

julia> frame_index(simulation)
3
```

"""
function first_frame!(simulation::Simulation)
    lock(simulation) do
        ifirst = first(simulation.frame_range)
        simulation.frame_index = ifirst
        read_step!(simulation.trajectory, ifirst, simulation.frame)
        _set_fixed_unitcell!(simulation)
    end
    return simulation.frame
end

"""
    next_frame!(simulation::Simulation)

Reads the next frame in the trajectory file and returns it. Moves the current
frame to the next one in the range to be considered (given by `frame_range(simulation)`).

"""
function next_frame!(simulation::Simulation)
    lock(simulation) do
        i_frame_in_range = searchsortedfirst(frame_range(simulation), frame_index(simulation))
        if i_frame_in_range + 1 > length(frame_range(simulation))
            throw(ArgumentError("""\n
                Next frame out of range.
                Current frame: $(frame_index(simulation)) is the last in selected frames $(_print_frame_range(simulation)).
                
            """))
        end
        i_next_frame = frame_range(simulation)[i_frame_in_range + 1]
        read_step!(simulation.trajectory, i_next_frame, simulation.frame)
        _set_fixed_unitcell!(simulation)
        simulation.frame_index = i_next_frame
    end
    return current_frame(simulation)
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
    else
        next_frame!(simulation)
    end
    return (current_frame(simulation), frame_index(simulation))
end

"""
    get_frame!(simulation::Simulation, iframe::Integer)

Returns the frame at the given index in the trajectory, and places the trajectory at that frame.

## Example

```jldoctest
julia> using MolSimToolkit, MolSimToolkit.Testing, LazyArtifacts, PDBTools

julia> sim = Simulation(test_file("popc_membrane.pdb"), test_file("popc_membrane.dcd"));

julia> frame4 = get_frame!(sim, 4)
Frame{Chemfiles.Frame} - first atom position: [-0.39684221148490906, 12.047886848449707, 37.44126892089844]
```

!!! note
    The `get_frame!` function will read the frames in the trajectory until the desired frame is reached. 
    This can be slow for large trajectories. If the required frame is before the current frame of the 
    simulation, the simulation will be restarted. The simulation object is returned positioned in the
    required frame. 

"""
function get_frame!(simulation::Simulation, iframe::Integer)
    lock(simulation) do
        if !(iframe in frame_range(simulation))
            throw(ArgumentError("get_frame!: Index $iframe selected, but available frames are: $(frame_range(simulation))."))
        end
        i_current_frame = frame_index(simulation)
        if iframe < i_current_frame
            first_frame!(simulation)
        end
        while frame_index(simulation) < iframe
            next_frame!(simulation)
        end
    end
    return current_frame(simulation)
end

"""
    goto_frame!(simulation::Simulation, iframe::Integer)

Returns the simulation placed at the given index in the trajectory.

## Example

```jldoctest
julia> using MolSimToolkit, MolSimToolkit.Testing, LazyArtifacts, PDBTools

julia> sim = Simulation(test_file("popc_membrane.pdb"), test_file("popc_membrane.dcd"));

julia> goto_frame!(sim, 4)
Simulation 
    Atom type: Atom{Nothing}
    PDB file: /home/leandro/.julia/artifacts/86db7e8c748e1097601a090dcc7304f11dc08e7a/popc_membrane.pdb
    Trajectory file: /home/leandro/.julia/artifacts/86db7e8c748e1097601a090dcc7304f11dc08e7a/popc_membrane.dcd
    Total number of frames: 5
    Frames to consider: 1, 2, 3, 4, 5
    Number of frames to consider: 5
    Current frame: 4

```

!!! note
    The `goto_frame!` function will read the frames in the trajectory until the desired frame is reached. 
    This can be slow for large trajectories. If the required frame is before the current frame of the 
    simulation, the simulation will be restarted. The simulation object is returned positioned in the
    required frame. 

"""
function goto_frame!(simulation::Simulation, iframe::Integer)
    get_frame!(simulation, iframe)
    return simulation
end

@testitem "Simulation" begin
    import Chemfiles
    using MolSimToolkit.Testing, LazyArtifacts

    #
    # Test iterator by reading coordinates
    #

    # Read coordinates with Chemfiles interface
    t = Chemfiles.Trajectory(test_file("popc_membrane.dcd"))
    c = zeros(Chemfiles.length(t))
    f = Chemfiles.read(t)
    c = [Chemfiles.positions(f)[1, 1]]
    for i in 2:Chemfiles.length(t)
        Chemfiles.read!(t, f)
        push!(c, Chemfiles.positions(f)[1, 1])
    end
    Chemfiles.close(t)

    # Read coordinates with the Simulation interface
    simulation = Simulation(test_file("popc_membrane.pdb"), test_file("popc_membrane.dcd"))
    c2 = Float64[]
    for frame in simulation
        push!(c2, positions(frame)[1].x)
    end
    close(simulation)
    @test c == c2
end

@testitem "no-pdb construct" begin
    import PDBTools
    using MolSimToolkit.Testing, LazyArtifacts
    pdb = PDBTools.read_pdb(test_file("popc_membrane.pdb"))
    simulation = Simulation(pdb, test_file("popc_membrane.dcd"))
    first_frame!(simulation)
    @test positions(current_frame(simulation))[1].x == 5.912472724914551
    @test isnothing(path_pdb(simulation))
end

@testitem "get_frame!" begin
    using MolSimToolkit, MolSimToolkit.Testing, LazyArtifacts, PDBTools
    sim = Simulation(test_file("popc_membrane.pdb"), test_file("popc_membrane.dcd"))
    frames = [ copy(positions(frame)) for frame in sim ]
    @test all(positions(get_frame!(sim, i)) ≈ frames[i] for i in eachindex(sim))
    @test positions(get_frame!(sim, 5)) ≈ frames[5]
    @test positions(get_frame!(sim, 1)) ≈ frames[1]
    @test positions(get_frame!(sim, 2)) ≈ frames[2]
    ats = read_pdb(test_file("popc_membrane.pdb"))
    sim = Simulation(ats, test_file("popc_membrane.dcd"))
    @test all(positions(get_frame!(sim, i)) ≈ frames[i] for i in eachindex(sim))
    @test_throws ArgumentError get_frame!(sim, 100)
    sim2 = Simulation(test_file("popc_membrane.pdb"), test_file("popc_membrane.dcd"); frames=[1,2,5])
    @test positions(get_frame!(sim2, 1)) ≈ frames[1]
    @test positions(get_frame!(sim2, 2)) ≈ frames[2]
    @test positions(get_frame!(sim2, 5)) ≈ frames[5]
    @test_throws ArgumentError get_frame!(sim2, 4)
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

    simulation = Simulation(test_file("popc_membrane.pdb"), test_file("popc_membrane.dcd"); first = 2, step = 2, last = 4)
    @test first(positions(get_frame!(simulation, 2))) ≈ [2.001203775405884, 7.035768985748291, 43.27251434326172]
    @test first(positions(get_frame!(simulation, 4))) ≈ [-0.39684221148490906, 12.047886848449707, 37.44126892089844]
    @test_throws "get_frame!: Index 1 selected, but available frames are: [2, 4]" get_frame!(simulation, 1)

end

@testitem "Simulation - show" begin
    using MolSimToolkit
    using MolSimToolkit.Testing, LazyArtifacts
    using PDBTools
    using ShowMethodTesting: parse_show
    sim = Simulation(test_file("popc_membrane.pdb"), test_file("popc_membrane.dcd"))
    # replacement to avoid errors associated to type names
    repl = Dict(
        "MolSimToolkit." => "", 
        r"PDBTools(Base)?\." => "",
        r"(PDB file:).*" => s"\1",
        r"(Trajectory file:).*" => s"\1",
    )
    @test parse_show(sim; repl=repl) ≈ """
    Simulation 
    Atom type: Atom{Nothing}
    PDB file:
    Trajectory file:
    Total number of frames: 5
    Frames to consider: 1, 2, 3, 4, 5
    Number of frames to consider: 5
    Current frame: 1
    """
    @test parse_show(repr(sim; context= :compact => true), repl=repl) ≈ "Simulation(popc_membrane.pdb, popc_membrane.dcd)"
end

@testitem "Simulation - multiple trajectories" begin
    using MolSimToolkit, MolSimToolkit.Testing, LazyArtifacts
    using PDBTools
    using ShowMethodTesting: parse_show
    pdb, dcd = test_file("popc_membrane.pdb"), test_file("popc_membrane.dcd")
    sim1 = Simulation(pdb, dcd)
    frames = [copy(positions(frame)) for frame in sim1]
    sim = Simulation(pdb, dcd, dcd)
    @test raw_length(sim) == 10
    @test length(sim) == 10
    @test path_trajectory(sim) == [normpath(dcd), normpath(dcd)]
    @test all(positions(frame) ≈ frames[mod1(i, 5)] for (i, frame) in pairs(sim))
    @test positions(get_frame!(sim, 7)) ≈ frames[2]
    @test positions(get_frame!(sim, 3)) ≈ frames[3]
    sim = Simulation(pdb, [dcd, dcd, dcd]; first=4, step=3)
    @test frame_range(sim) == [4, 7, 10, 13]
    @test [copy(positions(frame)) for frame in sim] ≈ frames[[4, 2, 5, 3]]
    sim = Simulation(read_pdb(pdb), dcd, dcd)
    @test raw_length(sim) == 10
    @test positions(get_frame!(sim, 10)) ≈ frames[5]
    @test_throws ArgumentError Simulation(pdb, dcd, test_file("short_nopbc.dcd"))
    @test_throws ArgumentError Simulation(pdb, String[])
    sim = Simulation(pdb, dcd, dcd; unitcell=true)
    @test all(unitcell(frame) == read_pdb(pdb).unitcell for frame in sim)
    repl = Dict(
        "MolSimToolkit." => "", 
        r"PDBTools(Base)?\." => "",
        r"(PDB file:).*" => s"\1",
        r"/.*/popc_membrane.dcd" => "popc_membrane.dcd",
    )
    sim = Simulation(pdb, dcd, dcd)
    @test parse_show(sim; repl=repl) ≈ """
    Simulation 
    Atom type: Atom{Nothing}
    PDB file:
    Trajectory files:
        popc_membrane.dcd
        popc_membrane.dcd
    Total number of frames: 10
    Frames to consider: 1, 2, 3, ... , 8, 9, 10
    Number of frames to consider: 10
    Current frame: 1
    """
    @test parse_show(repr(sim; context= :compact => true), repl=repl) ≈ "Simulation(popc_membrane.pdb, popc_membrane.dcd, popc_membrane.dcd)"
    close(sim)
end
