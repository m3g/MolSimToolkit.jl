using Statistics: mean
using DelimitedFiles: readdlm
using PDBTools: Atom, read_pdb, write_pdb, write_pdb_atom
import MDLovoFit_jll
using ChunkSplitters: chunks

# Functions of the interface
export MDLovoFitResult
export mdlovofit
export MapFractionsResult
export map_fractions

#=

Function to write a temporary PDB file with the selected frames of a trajectory, to run the mdlovofit executable.

Input: 
- simulation: Simulation object with the trajectory of the system and the atom data.
- maxframes: number of frames to be considered. If nothing, 100 frames are considered.

Output:
- tmp_trajectory_file: name of the temporary PDB file with the selected frames.
- iframes: the frame index of each frame considered.

=#
function _write_temporary_trajectory(simulation::Simulation, maxframes)
    if isnothing(maxframes) 
        if length(simulation) > 100 
            @warn """\n
                Running map_fractions depends on writting a temporary trajectory file with the selected frames.
    
                Since `maxframes` was not set, the number of frames considered will be 100.
    
                The considered frames will be sampled equally spaced from the simulation, as possible. The indices frames
                considered will be stored in the `iframes` field of the output data structure.
    
                To avoid this warning, set the `maxframes` parameter to the desired number of frames to be considered.
                Use `maxframes=length(simulation)` to consider all frames.
    
            """ _file = nothing _line = nothing 
            iframes = first.(chunks(frame_range(simulation); n=100))
        else
            iframes = frame_range(simulation)
        end
    else
        if maxframes > length(simulation)
            throw(ArgumentError("""\n
                maxframes=$maxframes is greater than the number of frames in the simulation: $(length(simulation)).
    
                To consider all frames, set `maxframes=length(simulation)`.
    
            """))
        end
        iframes = first.(chunks(frame_range(simulation); n=maxframes))
    end

    #
    # write the temporary PDB file with the frames to be considered, to run the mdlovofit executable
    #
    # Atoms to be considered (could be generalized with the -atomsfile option of mdlovofit)
    cA = filter(at -> PDBTools.isprotein(at) && PDBTools.name(at) == "CA", atoms(simulation))
    cA_inds = PDBTools.index.(cA)
    tmp_trajectory_file = tempname()*"_mdlovofit_trajectory.pdb"
    for (iframe, frame) in enumerate(simulation)
        if iframe in iframes
            p = positions(frame)[cA_inds]
            for (iat, at) in enumerate(cA)
                at.x = p[iat].x
                at.y = p[iat].y
                at.z = p[iat].z
            end
            write_pdb(tmp_trajectory_file, cA; append=true)
        end
    end

    return tmp_trajectory_file, iframes 
end

"""
    MapFractionsResult

Data structure to store the output of the `map_fractions` function.

Fields:

- `fraction` contains the fraction of atoms considered in the alignment.
- `rmsd_low` contains the RMSD of the fraction of the structure with the lowest RMSD.
- `rmsd_high` contains the RMSD of the fraction not considered for the alignment.
- `rmsd_all` contains the RMSD of the whole structure.

"""
struct MapFractionsResult
    simulation::Simulation
    iframes::Vector{Int}
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
    - iframes: the frame index of each frame considered.
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
    MDLovoFitResult

Data structure to store the output of the `mdlovofit` function.

Fields: 

- `iframes` vector with the frame index of each frame.
- `rmsd_low` is the RMSD of the fraction of the structure with the lowest RMSD.
- `rmsd_high` is the RMSD of the fraction of the structure with the highest RMSD.
- `rmsd_all` is the RMSD of the whole structure.
- `rmsf` is the RMSF as a function of the residue or atom index. 
- `aligned_pdb` is the name of the PDB file with the aligned structure.

"""
struct MDLovoFitResult
    fraction::Float64
    iframes::Vector{Int}
    rmsd_low::Vector{Float64}
    rmsd_high::Vector{Float64}
    rmsd_all::Vector{Float64}
    rmsf::Vector{Float64}
    rmsf_file::String
    rmsd_file::String
    aligned_pdb::String
end

function Base.show(io::IO, result::MDLovoFitResult)
    plow = round(100*result.fraction,digits=1)
    phigh = round(100*(1-result.fraction),digits=1)
    av_low = round((mean(result.rmsd_low)), digits=2)
    av_high = round((mean(result.rmsd_high)), digits=2)
    av_all = round((mean(result.rmsd_all)), digits=2)
    print(io, chomp("""
    -------------------------------------------------------------------
    MDLovoFitResult
    -------------------------------------------------------------------

    Aligned pdb file: $(result.aligned_pdb)
    RMSF data file: $(result.rmsf_file)
    RMSD data file: $(result.rmsd_file)

    Number of frames considered: $(length(result.iframes))
    Average RMSD of all atoms: $av_all
    Average RMSD of the $plow% atoms of lowest RMSD: $av_low
    Average RMSD of the $phigh% atoms of highest RMSD: $av_high

    Frame indices availabe in result.iframe
    RMSD data availabe in rmsd_low, rmsd_high, and rmsd_all

    RMSF data availabe in result.rmsf (Number of atoms: $(length(result.rmsf)))
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

julia> sim = Simulation(Testing.namd_pdb, Testing.namd_traj)

julia> mf = map_fractions(sim)
-------------------------------------------------------------------
MapFractionsResult: Simulation(structure.pdb, trajectory.dcd)

-------------------------------------------------------------------
Fields: 
- fraction: fraction of atoms considered in the alignment.
- iframes: the frame index of each frame considered.
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
    tmp_trajectory_file, iframes = _write_temporary_trajectory(simulation, maxframes)
    # Run mdlovofit
    mapfrac_file = tempname() *"_mapfrac.dat"
    MDLovoFit_jll.mdlovofit() do exe
        run(pipeline(`$exe -mapfrac $tmp_trajectory_file`; stdout=mapfrac_file))
    end
    data = readdlm(mapfrac_file, comments=true, comment_char='#')
    range = 1:findlast(<(1), data[:,1])
    fraction = data[range,1]
    rmsd_low = data[range,2]
    rmsd_high = data[range,3]
    rmsd_all = data[range,4]
    return MapFractionsResult(simulation, iframes, fraction, rmsd_low, rmsd_high, rmsd_all)
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

julia> sim = Simulation(Testing.namd_pdb, Testing.namd_traj)

julia> result = mdlovofit(sim, fraction=0.5, output_name="mdlovofit_50")
-------------------------------------------------------------------
MDLovoFitResult
-------------------------------------------------------------------

Aligned pdb file: mdlovofit_50_aligned.pdb
RMSF data file: mdlovofit_50_rmsf.dat
RMSD data file: mdlovofit_50_rmsd.dat

Number of frames considered: 5
Average RMSD of all atoms: 1.89
Average RMSD of the 50.0% atoms of lowest RMSD: 0.0
Average RMSD of the 50.0% atoms of highest RMSD: 3.7

Frame indices availabe in result.iframe
RMSD data availabe in rmsd_low, rmsd_high, and rmsd_all

RMSF data availabe in result.rmsf (Number of atoms: 43)
-------------------------------------------------------------------
```

"""
function mdlovofit(
    simulation::Simulation;
    fraction::AbstractFloat,
    output_name::Union{String,Nothing} = nothing,
    reference_frame::Int = 1,
    maxframes=nothing,
)
    # write temporary trajectory file
    tmp_trajectory_file, iframes = _write_temporary_trajectory(simulation, maxframes)
    # Run MDLovoFit
    rmsf_file, rmsd_file, output_pdb = if isnothing(output_name) 
        @warn """\n
            `output_name` was not provided. Default names will be used.

        """ _file = nothing _line = nothing
        tmpfile = "mdlovofit"
        tmpfile*"_rmsf.dat", tmpfile*"_rmsd.dat", tmpfile*"_aligned.pdb"
    else
        output_name*"_rmsf.dat", output_name*"_rmsd.dat", output_name*"_aligned.pdb"
    end
    try
        MDLovoFit_jll.mdlovofit() do exe
            run(pipeline(`$exe -f $fraction -iref $reference_frame -rmsf $rmsf_file -t $output_pdb $tmp_trajectory_file`; stdout=rmsd_file))
        end
    catch
        "ERROR in MDLovoFit execution"
        "Command executed: $command"
    end

    # Read RMSD file
    rmsd_data = readdlm(rmsd_file; comments=true, comment_char='#')
    rmsd_low = rmsd_data[:,2]
    rmsd_high = rmsd_data[:,3]
    rmsd_all = rmsd_data[:,4]
    # Read RMSF file
    rmsf = readdlm(rmsf_file; comments=true, comment_char='#')[:,2]
    return MDLovoFitResult(fraction, iframes, rmsd_low, rmsd_high, rmsd_all, rmsf, rmsf_file, rmsd_file, output_pdb)
end

@testitem "map_fractions" begin
    using ShowMethodTesting
    using MolSimToolkit, MolSimToolkit.Testing 

    sim = Simulation(Testing.mdlovofit_pdb, Testing.mdlovofit_traj)
    mf = map_fractions(sim)
    @test length(mf.iframes) == 100
    @test sum(mf.iframes) == 7074
    @test sum(mf.rmsd_all) ≈ 65.12515215300002
    @test_throws ArgumentError map_fractions(sim; maxframes=200)
    mf = map_fractions(sim; maxframes=121)
    @test length(mf.iframes) == 121
    @test sum(mf.iframes) == 7620
    mf = map_fractions(sim; maxframes=50)
    @test length(mf.iframes) == 50
    @test last(mf.iframes) == 122
    sim = Simulation(Testing.mdlovofit_pdb, Testing.mdlovofit_traj; frames=1:70)
    mf = map_fractions(sim)
    @test length(mf.iframes) == 70
    @test sum(mf.rmsd_all) ≈ 66.11294932099999

    sim = Simulation(Testing.mdlovofit_pdb, Testing.mdlovofit_traj)
    mf = map_fractions(sim; maxframes=50)
    @test parse_show(mf) ≈ """
    -------------------------------------------------------------------
    MapFractionsResult: Simulation(structure.pdb, trajectory.dcd)
    
    -------------------------------------------------------------------
    Fields: 
    - fraction: fraction of atoms considered in the alignment.
    - iframes: the frame index of each frame considered.
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

    sim = Simulation(Testing.mdlovofit_pdb, Testing.mdlovofit_traj)
    @test_throws UndefKeywordError mdlovofit(sim)

    fit = mdlovofit(sim; fraction=0.50)
    @test isfile("mdlovofit_aligned.pdb")
    @test isfile("mdlovofit_rmsf.dat")
    @test isfile("mdlovofit_rmsd.dat")
    @test length(fit.iframes) == 100
    @test sum(fit.rmsd_low) < sum(fit.rmsd_high)

    sim = Simulation(Testing.mdlovofit_pdb, Testing.mdlovofit_traj; frames=1:50)
    fit = mdlovofit(sim; fraction=0.50, output_name="mdlovofit_50")
    @test isfile("mdlovofit_50_aligned.pdb")
    @test isfile("mdlovofit_50_rmsf.dat")
    @test isfile("mdlovofit_50_rmsd.dat")
    @test length(fit.iframes) == 50
    @test sum(fit.rmsd_low) < sum(fit.rmsd_high)
end