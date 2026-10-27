"""
    Trajectory

Structure that contains the data stream of one or more trajectory files. If more than one
file is provided, the files are considered sequential parts of the same trajectory, and
the frames are indexed continuously across files.

All methods dispatching on `Trajectory` are considered internal.

"""
struct Trajectory{T<:Chemfiles.Trajectory}
    trajectories::Vector{T}
    # number of frames in the files preceding each file
    offsets::Vector{Int}
end
Trajectory(filename::AbstractString) = Trajectory([filename])
function Trajectory(filenames::AbstractVector{<:AbstractString})
    if isempty(filenames)
        throw(ArgumentError("At least one trajectory file must be provided."))
    end
    trajectories = [Chemfiles.Trajectory(String(file)) for file in filenames]
    lengths = [Int(Chemfiles.length(t)) for t in trajectories]
    offsets = cumsum(lengths) .- lengths
    if length(trajectories) > 1
        natoms = [Int(Chemfiles.length(Chemfiles.read_step(t, 0))) for t in trajectories if Chemfiles.length(t) > 0]
        if !allequal(natoms)
            foreach(Chemfiles.close, trajectories)
            throw(ArgumentError("""\n
                All trajectory files must have the same number of atoms.
                Got number of atoms: $(join(natoms, ", ")) for files: $(join(filenames, ", ")).

            """))
        end
    end
    return Trajectory(trajectories, offsets)
end
Base.length(t::Trajectory) = last(t.offsets) + Int(Chemfiles.length(last(t.trajectories)))
Base.read(t::Trajectory) = Frame(Chemfiles.read(first(t.trajectories)))
Base.close(t::Trajectory) = foreach(Chemfiles.close, t.trajectories)
function path_trajectory(t::Trajectory)
    paths = [normpath(Chemfiles.path(traj)) for traj in t.trajectories]
    return length(paths) == 1 ? only(paths) : paths
end
raw_length(t::Trajectory) = length(t)
function read_step!(t::Trajectory, i_next_frame::Int, f::Frame)
    ifile = searchsortedlast(t.offsets, i_next_frame - 1)
    Chemfiles.read_step!(t.trajectories[ifile], i_next_frame - 1 - t.offsets[ifile], f.frame)
    copyto!(f.positions, reinterpret(Point3D{Float64}, Chemfiles.positions(f.frame)))
    f.unitcell_matrix .= transpose(Chemfiles.matrix(Chemfiles.UnitCell(f.frame)))
    return f
end

@testitem "Trajectory" begin
    using MolSimToolkit, MolSimToolkit.Testing, LazyArtifacts
    dcd = test_file("popc_membrane.dcd")
    # single file
    t = Trajectory(dcd)
    @test t isa Trajectory
    @test length(t) == raw_length(t) == 5
    @test path_trajectory(t) == normpath(dcd)
    f = read(t)
    sim = Simulation(test_file("popc_membrane.pdb"), dcd)
    @test positions(MolSimToolkit.read_step!(t, 3, f)) ≈ positions(get_frame!(sim, 3))
    close(t)
    # multiple files
    t = Trajectory([dcd, dcd])
    @test length(t) == raw_length(t) == 10
    @test t.offsets == [0, 5]
    @test path_trajectory(t) == [normpath(dcd), normpath(dcd)]
    @test positions(MolSimToolkit.read_step!(t, 8, f)) ≈ positions(get_frame!(sim, 3))
    close(t)
    @test_throws ArgumentError Trajectory(String[])
end
