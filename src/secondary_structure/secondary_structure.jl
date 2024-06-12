#
# Computes the secondary structure in a single frame
# This is an internal function.
#
function _ss_frame!(
    atoms::AbstractVector{<:PDBTools.Atom},
    frame::Chemfiles.Frame;
    ss_method::F=stride_run
) where {F<:Function}
    coordinates = positions(frame)
    for at in atoms
        iatom = PDBTools.index(at)
        atoms[iatom].x = coordinates[iatom].x
        atoms[iatom].y = coordinates[iatom].y
        atoms[iatom].z = coordinates[iatom].z
    end
    return ss_method(atoms)
end

"""
    ss_map(
        simulation::Simulation; 
        selection::Union{AbstractString,Function}=PDBTools.isprotein,
        ss_method=stride_run,
        show_progress=true
    )

Calculates the secondary structure map of the trajectory. 
Returns a matrix of secondary structure codes, where each row is a residue and each column is a frame.

By default, all protein atoms are considered. The `selection` keyword argument can be used to choose a different
selection, for example, `selection="protein and chain A"`.

The `ss_method` keyword argument can be used to choose the secondary structure prediction method,
which can be either `stride_run` or `dssp_run`. The default is `stride_run`.

The `show_progress` keyword argument controls whether a progress bar is shown.

For the classes, refer to the ProteinSecondaryStructures.jl package documentation.

"""
function ss_map(
    simulation::Simulation;
    selection::Union{AbstractString,Function}=PDBTools.isprotein,
    ss_method::F=stride_run,
    show_progress=true
) where {F<:Function}
    sel = PDBTools.select(atoms(simulation), selection)
    ss_map = zeros(Int, length(PDBTools.eachresidue(sel)), length(simulation))
    p = Progress(length(simulation); enabled=show_progress)
    for (iframe, frame) in enumerate(simulation)
        ss = _ss_frame!(sel, frame; ss_method)
        for (i, ssdata) in pairs(ss)
            ss_map[i, iframe] = ss_code_to_number(ssdata.sscode)
        end
        next!(p)
    end
    return ss_map
end

@testitem "ss_map" begin
    using MolSimToolkit
    using MolSimToolkit.Testing
    import PDBTools
    simulation = Simulation(Testing.namd_pdb, Testing.namd_traj)
    # With stride
    ssmap = ss_map(simulation; ss_method=stride_run)
    @test size(ssmap) == (43, 5)
    @test sum(ssmap) == 842
    # With DSSP
    ssmap = ss_map(PDBTools.readPDB(pdbfile, "protein"), trajectory; method=dssp_run)
    @test size(ssmap) == (43, 5)
    @test sum(ssmap) == 1006
end

"""
    ss_class_content(ss_class, ssmap::AbstractMatrix{<:Integer})

Calculates the content, en each frame of a trajectory, of a specific class
of secondary structure. 

Here, `ss_class` is the secondary structure class, as a string, character or integer,
and `ssmap` is the secondary structure map of the trajectory, as returned by the `ss_map` 
function.

The classes can be found in the ProteinSecondaryStructures.jl package documentation, at:

https://m3g.github.io/ProteinSecondaryStructures.jl/stable/explanation/#Secondary-structure-classes

"""
function ss_class_content(
    ss_class::Union{String,Char,Integer},
    ssmap::AbstractMatrix{Int},
)
    ss_class = ss_class isa Integer ? ss_class : ss_code_to_number(ss_class)
    ss_class_content = zeros(Float64, size(ssmap, 2))
    for (iframe, ss_frame) in enumerate(eachcol(ssmap))
        ss_class_content[iframe] = count(==(ss_class), ss_frame) / max(1, length(ss_frame))
    end
    return ss_class_content
end

"""
    ssmap_frame(ssmap::AbstractMatrix{Int}, iframe::Int}

Calculates the secondary structure composition of a frame of the simulation,
given the secondary structure map. 

Returns a dictionary of the secondary structure types and their counts, for the chosen frame.

"""
function ss_frame(ssmap::AbstractMatrix{Int}, iframe::Int)
    ssmap_frame= Dict{String,Int}()
    sscodes = @view(ssmap[:, iframe])
    for sscode in ss_code_to_number.(keys(ss_classes))
        ssmap_frame[class(sscode)] = count(==(sscode), sscodes)
    end
    return ssmap_frame
end

@testitem "ss_content/ss_composition" begin
    using MolSimToolkit.Testing
    using PDBTools: readPDB
    using Chemfiles: Trajectory
    pdbfile = joinpath(Testing.data_dir, "Gromacs", "system.pdb")
    trajectory = Trajectory(joinpath(Testing.data_dir, "Gromacs", "trajectory.xtc"))
    # With stride
    helical_content = ss_content(is_anyhelix, readPDB(pdbfile, "protein"), trajectory; method=stride_run, show_progress=false)
    @test length(helical_content) == 26
    @test sum(helical_content) / length(helical_content) ≈ 0.2181174089068826
    # With DSSP
    helical_content = ss_content(is_anyhelix, readPDB(pdbfile, "protein"), trajectory; method=dssp_run, show_progress=false)
    @test length(helical_content) == 26
    @test sum(helical_content) / length(helical_content) ≈ 0.21204453441295545
    # With non-contiguous indexing
    atoms = readPDB(pdbfile, only=at -> (10 <= at.residue < 30) | (40 <= at.residue < 60))
    helical_content = ss_content(is_anyhelix, atoms, trajectory; method=stride_run, show_progress=false)
    @test length(helical_content) == 26
    @test sum(helical_content) / length(helical_content) ≈ 0.20288461538461539
    # From the map
    ssmap = ss_map(readPDB(pdbfile, "protein"), trajectory; method=stride_run, show_progress=false)
    helical_content = ss_content(is_anyhelix, ssmap)
    @test length(helical_content) == 26
    @test sum(helical_content) / length(helical_content) ≈ 0.2181174089068826
    # Test ss_composition function
    @test ss_composition(ssmap, 5) == Dict("310 helix" => 6, "bend" => 0, "turn" => 17, "kappa helix" => 0, "beta strand" => 25, "beta bridge" => 2, "alpha helix" => 12, "pi helix" => 0, "loop" => 0, "coil" => 14)
end