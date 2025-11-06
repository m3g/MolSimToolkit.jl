using Statistics: mean
using DelimitedFiles: readdlm
using ChunkSplitters: chunks

# Functions of the interface
export MDLovoFitResult
export mdlovofit
export MapFractionsResult
export map_fractions

import MDLovoFit_jll
const mdlovofit_runner = MDLovoFit_jll.mdlovofit()

#=

Function to write a temporary PDB file with the selected frames of a trajectory, to run the mdlovofit executable.

Input: 
- simulation: Simulation object with the trajectory of the system and the atom data.
- maxframes: number of frames to be considered. If nothing, 100 frames are considered.

Output:
- tmp_trajectory_file: name of the temporary PDB file with the selected frames.
- frame_indices: the frame index of each frame considered.

=#
function _write_temporary_trajectory(simulation::Simulation, maxframes)
    if isnothing(maxframes)
        if length(simulation) > 100
            @warn """\n
                Running map_fractions depends on writting a temporary trajectory file with the selected frames.

                Since `maxframes` was not set, the number of frames considered will be 100.

                The considered frames will be sampled equally spaced from the simulation, as possible. The indices frames
                considered will be stored in the `frame_indices` field of the output data structure.

                To avoid this warning, set the `maxframes` parameter to the desired number of frames to be considered.
                Use `maxframes=length(simulation)` to consider all frames.

            """ _file = nothing _line = nothing
            frame_indices = first.(chunks(frame_range(simulation); n=100))
        else
            frame_indices = frame_range(simulation)
        end
    else
        if maxframes > length(simulation)
            throw(ArgumentError("""\n
                maxframes=$maxframes is greater than the number of frames in the simulation: $(length(simulation)).

                To consider all frames, set `maxframes=length(simulation)`.

            """))
        end
        frame_indices = first.(chunks(frame_range(simulation); n=maxframes))
    end

    #
    # write the temporary PDB file with the frames to be considered, to run the mdlovofit executable
    #
    # Atoms to be considered (could be generalized with the -atomsfile option of mdlovofit)
    cA = filter(at -> PDBTools.isprotein(at) && PDBTools.name(at) == "CA", atoms(simulation))
    cA_inds = PDBTools.index.(cA)
    tmp_trajectory_file = tempname() * "_mdlovofit_trajectory.pdb"
    for (iframe, frame) in enumerate(simulation)
        iframe_in_range = frame_range(simulation)[iframe]
        if iframe_in_range in frame_indices
            p = positions(frame)[cA_inds]
            PDBTools.set_position!.(cA, p)
            PDBTools.write_pdb(tmp_trajectory_file, cA; append=true)
        end
    end

    return tmp_trajectory_file, frame_indices
end

"""
    MapFractionsResult

Data structure to store the output of the `map_fractions` function.

Fields:

- `simulation` is a Simulation object with the trajectory of the system and the atom data.
- `frame_indices` contains the frame index of each frame considered.
- `fraction` contains the fraction of atoms considered in the alignment.
- `rmsd_low` contains the RMSD of the fraction of the structure with the lowest RMSD.
- `rmsd_high` contains the RMSD of the fraction not considered for the alignment.
- `rmsd_all` contains the RMSD of the whole structure.

"""
struct MapFractionsResult
    simulation::Simulation
    frame_indices::Vector{Int}
    fraction::Vector{Float64}
    rmsd_low::Vector{Float64}
    rmsd_high::Vector{Float64}
    rmsd_all::Vector{Float64}
end

function Base.show(io::IO, mf::MapFractionsResult)
    print(io, chomp("""
    -------------------------------------------------------------------
    MapFractionsResult: 
    """))
    show(IOContext(io, :compact => true), mf.simulation)
    print(io, chomp("""\n
    -------------------------------------------------------------------
    Fields: 
    - fraction: fraction of atoms considered in the alignment.
    - frame_indices: the frame index of each frame considered.
    - rmsd_low: RMSD of the fraction of the structure with the lowest RMSD.
    - rmsd_high: RMSD of the fraction not considered for the alignment.
    - rmsd_all: RMSD of the whole structure.

    Greatest fraction for which the RMSD-low is smaller than 1.0: $(round(mf.fraction[findlast(<(1.0), mf.rmsd_low)],digits=2))
                                                             2.0: $(round(mf.fraction[findlast(<(2.0), mf.rmsd_low)],digits=2))
                                                             3.0: $(round(mf.fraction[findlast(<(3.0), mf.rmsd_low)],digits=2))
    -------------------------------------------------------------------
    """))
end

"""
    map_fractions(simulation::Simulation; maxframes=nothing)

Run MDLovoFit on a trajectory, to obtain the RMSD of the fractions of the structure
with the lowest RMSD, the fraction with the highest RMSD, and the whole structure.

Only CA atoms of protein residues are considered. 

Argument:

- `simulation` is a Simulation object, with the trajectory of the system and the atom data.

Keyword arguments:

- `maxframes=nothing` is the number of frames to be considered. If nothing, 100 frames are considered. 
   Set to `length(simulation)` to consider all frames. 

## Example

```jldoctest
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> sim = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> mf = map_fractions(sim)
-------------------------------------------------------------------
MapFractionsResult: Simulation(structure.pdb, trajectory.dcd)

-------------------------------------------------------------------
Fields: 
- fraction: fraction of atoms considered in the alignment.
- frame_indices: the frame index of each frame considered.
- rmsd_low: RMSD of the fraction of the structure with the lowest RMSD.
- rmsd_high: RMSD of the fraction not considered for the alignment.
- rmsd_all: RMSD of the whole structure.

Greatest fraction for which the RMSD-low is smaller than 1.0: 0.79
                                                         2.0: 0.9
                                                         3.0: 0.99
-------------------------------------------------------------------
```

"""
function map_fractions(
    simulation::Simulation;
    maxframes=nothing,
)
    # write temporary trajectory file
    tmp_trajectory_file, frame_indices = _write_temporary_trajectory(simulation, maxframes)
    # Run mdlovofit
    mapfrac_file = tempname() * "_mapfrac.dat"
    run(pipeline(`$mdlovofit_runner -mapfrac $tmp_trajectory_file`; stdout=mapfrac_file))
    data = readdlm(mapfrac_file, comments=true, comment_char='#')
    range = 1:findlast(<(1), data[:, 1])
    fraction = data[range, 1]
    rmsd_low = data[range, 2]
    rmsd_high = data[range, 3]
    rmsd_all = data[range, 4]
    return MapFractionsResult(simulation, frame_indices, fraction, rmsd_low, rmsd_high, rmsd_all)
end

"""
    MDLovoFitResult

Data structure to store the output of the `mdlovofit` function.

Fields: 

- `simulation`: Simulation object with the trajectory of the system and the atom data. 
- `frame_indices`: vector with the frame index of each frame.
- `rmsd_low`: RMSD of the fraction of the structure with the lowest RMSD.
- `rmsd_high`: RMSD of the fraction of the structure with the highest RMSD.
- `rmsd_all`: RMSD of the whole structure.
- `rmsf`: RMSF as a function of the residue or atom index. 
- `rmsf_file`: name of the file with the RMSF data.
- `rmsd_file`: name of the file with the RMSD data.
- `aligned_pdb`: name of the PDB file with the aligned structure.

"""
struct MDLovoFitResult
    simulation::Simulation
    fraction::Float64
    frame_indices::Vector{Int}
    rmsd_low::Vector{Float64}
    rmsd_high::Vector{Float64}
    rmsd_all::Vector{Float64}
    rmsf::Vector{Float64}
    rmsf_file::String
    rmsd_file::String
    aligned_pdb::String
end

function Base.show(io::IO, result::MDLovoFitResult)
    plow = round(100 * result.fraction, digits=1)
    phigh = round(100 * (1 - result.fraction), digits=1)
    av_low = round((mean(result.rmsd_low)), digits=2)
    av_high = round((mean(result.rmsd_high)), digits=2)
    av_all = round((mean(result.rmsd_all)), digits=2)
    print(io, chomp("""
    -------------------------------------------------------------------
    MDLovoFitResult: 
    """))
    show(IOContext(io, :compact => true), result.simulation)
    print(io, chomp("""\n
    -------------------------------------------------------------------

    Aligned pdb file: $(result.aligned_pdb)
    RMSF data file: $(result.rmsf_file)
    RMSD data file: $(result.rmsd_file)

    Number of frames considered: $(length(result.frame_indices))
    Average RMSD of all atoms: $av_all
    Average RMSD of the $plow% atoms of lowest RMSD: $av_low
    Average RMSD of the $phigh% atoms of highest RMSD: $av_high

    Frame indices availabe in field: frame_indices
    RMSD data availabe in fields: rmsd_low, rmsd_high, and rmsd_all

    RMSF data availabe in field: rmsf (Number of atoms: $(length(result.rmsf)))
    -------------------------------------------------------------------
    """))
end

"""
    mdlovofit(
        simulation::Simulation; 
        fraction::AbstractFloat, 
        output_name::Union{String,Nothing} = nothing, 
        reference_frame::Int = 1, 
        maxframes=nothing
    )

Run MDLovoFit on a trajectory, aligning only the atoms which are in the fraction of the structure with the lowest RMSD.
Only CA atoms of protein residues are considered.

Arguments:

- `simulation` is a Simulation object, with the trajectory of the system and the atom data.
- `fraction` is the fraction of atoms to be considered in the alignment.

Optional arguments:

- `output_name=nothing` is the base name of the output files. If nothing, default names will be used.
- `reference_frame=1` is the index of the frame to be used as reference for the alignment.
- `maxframes=nothing` is the number of frames to be considered. If nothing, 100 frames are considered. 
   Set to `length(simulation)` to consider all frames.

Output: 

- `MDLovoFitResult` data structure with the results of the alignment.

and the output files:
- `output_name_rmsf.dat` with the RMSF data.
- `output_name_rmsd.dat` with the RMSD data.
- `output_name_aligned.pdb` with the aligned structure.

## Example

```jldoctest
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> sim = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> md = mdlovofit(sim, fraction=0.7, output_name="mdlovofit_50")
-------------------------------------------------------------------
MDLovoFitResult: Simulation(structure.pdb, trajectory.dcd)

-------------------------------------------------------------------

Aligned pdb file: mdlovofit_50_aligned.pdb
RMSF data file: mdlovofit_50_rmsf.dat
RMSD data file: mdlovofit_50_rmsd.dat

Number of frames considered: 5
Average RMSD of all atoms: 1.79
Average RMSD of the 70.0% atoms of lowest RMSD: 0.53
Average RMSD of the 30.0% atoms of highest RMSD: 4.71

Frame indices availabe in field: frame_indices
RMSD data availabe in fields: rmsd_low, rmsd_high, and rmsd_all

RMSF data availabe in field: rmsf (Number of atoms: 43)
-------------------------------------------------------------------
```

"""
function mdlovofit(
    simulation::Simulation;
    fraction::AbstractFloat,
    output_name::Union{String,Nothing}=nothing,
    reference_frame::Int=1,
    maxframes=nothing,
)
    # write temporary trajectory file
    tmp_trajectory_file, frame_indices = _write_temporary_trajectory(simulation, maxframes)
    # Run MDLovoFit
    rmsf_file, rmsd_file, output_pdb = if isnothing(output_name)
        @warn """\n
            `output_name` was not provided. Default names will be used.

        """ _file = nothing _line = nothing
        tmpfile = "mdlovofit"
        tmpfile * "_rmsf.dat", tmpfile * "_rmsd.dat", tmpfile * "_aligned.pdb"
    else
        output_name * "_rmsf.dat", output_name * "_rmsd.dat", output_name * "_aligned.pdb"
    end
    try
        run(pipeline(
            `$mdlovofit_runner -f $fraction -iref $reference_frame -rmsf $rmsf_file -t $output_pdb $tmp_trajectory_file`;
            stdout=rmsd_file
        ))
    catch
        "ERROR in MDLovoFit execution"
        "Command executed: $command"
    end

    # Read RMSD file
    rmsd_data = readdlm(rmsd_file; comments=true, comment_char='#')
    rmsd_low = rmsd_data[:, 2]
    rmsd_high = rmsd_data[:, 3]
    rmsd_all = rmsd_data[:, 4]
    # Read RMSF file
    rmsf = readdlm(rmsf_file; comments=true, comment_char='#')[:, 2]
    return MDLovoFitResult(
        simulation,
        fraction,
        frame_indices,
        rmsd_low,
        rmsd_high,
        rmsd_all,
        rmsf,
        rmsf_file,
        rmsd_file,
        output_pdb
    )
end

@testitem "map_fractions" begin
    using ShowMethodTesting
    using MolSimToolkit, MolSimToolkit.Testing

    sim = Simulation(Testing.mdlovofit_pdb, Testing.mdlovofit_traj)
    mf = map_fractions(sim)
    @test length(mf.frame_indices) == 100
    @test sum(mf.frame_indices) == 7074
    @test sum(mf.rmsd_all) ≈ 65.12515215300002
    @test_throws ArgumentError map_fractions(sim; maxframes=200)
    mf = map_fractions(sim; maxframes=121)
    @test length(mf.frame_indices) == 121
    @test sum(mf.frame_indices) == 7620
    mf = map_fractions(sim; maxframes=50)
    @test length(mf.frame_indices) == 50
    @test last(mf.frame_indices) == 122
    sim = Simulation(Testing.mdlovofit_pdb, Testing.mdlovofit_traj; frames=1:70)
    mf = map_fractions(sim)
    @test length(mf.frame_indices) == 70
    @test sum(mf.rmsd_all) ≈ 66.11294932099999

    sim = Simulation(Testing.namd_pdb, Testing.namd_traj; frames=[2, 5])
    mf = map_fractions(sim)
    @test mf.frame_indices == [2, 5]
    @test sum(mf.rmsd_low) ≈ 75.73356675900001
    @test sum(mf.rmsd_high) ≈ 769.6762035190002

    sim = Simulation(Testing.mdlovofit_pdb, Testing.mdlovofit_traj)
    mf = map_fractions(sim; maxframes=50)
    @test parse_show(mf) ≈ """
    -------------------------------------------------------------------
    MapFractionsResult: Simulation(structure.pdb, trajectory.dcd)

    -------------------------------------------------------------------
    Fields: 
    - fraction: fraction of atoms considered in the alignment.
    - frame_indices: the frame index of each frame considered.
    - rmsd_low: RMSD of the fraction of the structure with the lowest RMSD.
    - rmsd_high: RMSD of the fraction not considered for the alignment.
    - rmsd_all: RMSD of the whole structure.

    Greatest fraction for which the RMSD-low is smaller than 1.0: 0.99
                                                             2.0: 0.99
                                                             3.0: 0.99
    -------------------------------------------------------------------
    """
end

@testitem "mdlovofit" begin
    using ShowMethodTesting
    using MolSimToolkit, MolSimToolkit.Testing
    using PDBTools

    sim = Simulation(Testing.mdlovofit_pdb, Testing.mdlovofit_traj)
    @test_throws UndefKeywordError mdlovofit(sim)

    fit = mdlovofit(sim; fraction=0.50)
    @test isfile("mdlovofit_aligned.pdb")
    @test isfile("mdlovofit_rmsf.dat")
    @test isfile("mdlovofit_rmsd.dat")
    @test length(fit.frame_indices) == 100
    @test sum(fit.rmsd_low) < sum(fit.rmsd_high)

    sim = Simulation(Testing.mdlovofit_pdb, Testing.mdlovofit_traj; frames=1:50)
    fit = mdlovofit(sim; fraction=0.50, output_name="mdlovofit_50")
    @test isfile("mdlovofit_50_aligned.pdb")
    @test isfile("mdlovofit_50_rmsf.dat")
    @test isfile("mdlovofit_50_rmsd.dat")
    @test length(fit.frame_indices) == 50
    @test sum(fit.rmsd_low) < sum(fit.rmsd_high)

    sim = Simulation(Testing.namd_pdb, Testing.namd_traj; frames=[2, 5])
    md = mdlovofit(sim; fraction=0.7)
    pdb = read_pdb(md.aligned_pdb)
    @test length(eachmodel(pdb)) == 2
    @test md.frame_indices == [2, 5]
    @test md.rmsd_low ≈ [0.0, 0.540949029]

    sim = Simulation(Testing.mdlovofit_pdb, Testing.mdlovofit_traj; frames=1:50)
    fit = mdlovofit(sim; fraction=0.50, output_name="mdlovofit_50")
    @test parse_show(fit) ≈ """
     -------------------------------------------------------------------
     MDLovoFitResult: Simulation(structure.pdb, trajectory.dcd)
     
     -------------------------------------------------------------------
     
     Aligned pdb file: mdlovofit_50_aligned.pdb
     RMSF data file: mdlovofit_50_rmsf.dat
     RMSD data file: mdlovofit_50_rmsd.dat
     
     Number of frames considered: 50
     Average RMSD of all atoms: 0.94
     Average RMSD of the 50.0% atoms of lowest RMSD: 0.32
     Average RMSD of the 50.0% atoms of highest RMSD: 1.49
     
     Frame indices availabe in field: frame_indices
     RMSD data availabe in fields: rmsd_low, rmsd_high, and rmsd_all
     
     RMSF data availabe in field: rmsf (Number of atoms: 19)
     -------------------------------------------------------------------
    """
end