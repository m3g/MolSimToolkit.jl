#
# Computes the secondary structure in a single frame
# This is an internal function.
#
function ss_frame!(
    atoms::AbstractVector{<:PDBTools.Atom},
    frame::Simulation.Frame;
    ss_method::F=stride_run
) where {F<:Function}
    coordinates = positions(frame)
    for (icoor, pos) in enumerate(coordinates)
        iatom = findfirst(at -> index(at) == icoor, atoms)
        if !isnothing(iatom)
            atoms[iatom].x = pos.x
            atoms[iatom].y = pos.y
            atoms[iatom].z = pos.z
        else
            throw(ArgumentError("Atom of index $icoor not found in the provided atom list of atoms."))
        end
    end
    return ss_method(atoms)
end

"""
    ss_map(
        simulation::Simulation; 
        selection::Union{AbstractString,Function}=isprotein,
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
    selection::Union{AbstractString,Function}=isprotein,
    ss_method::F=stride_run,
    show_progress=true
) where {F<:Function}
    sel = select(atoms(simulation), selection)
    ss_map = zeros(Int, length(sel), length(simulation))
    p = Progress(length(simulation); enabled=show_progress)
    for (iframe, frame) in enumerate(simulation)
        ss = ss_frame!(atoms, frame; ss_method)
        for (i, ssdata) in pairs(ss)
            ss_map[i, iframe] = code_to_number[ssdata.sscode]
        end
        next!(p)
    end
    return ss_map
end

@testitem "ss_map" begin
    using MolSimToolkit
    import ProteinSecondaryStructures.Testing
    import PDBTools
    pdbfile = joinpath(Testing.data_dir, "Gromacs", "system.pdb")
    trajectory = Simulation(pdbfile, joinpath(Testing.data_dir, "Gromacs", "trajectory.xtc"))
    # With stride
    ssmap = ss_map(PDBTools.readPDB(pdbfile, "protein"), trajectory; method=stride_run)
    @test size(ssmap) == (76, 26)
    @test sum(ssmap) == 10577
    # With DSSP
    ssmap = ss_map(PDBTools.readPDB(pdbfile, "protein"), trajectory; method=dssp_run)
    @test size(ssmap) == (76, 26)
    @test sum(ssmap) == 12069
    # From the pdb file name
    ssmap = ss_map(pdbfile, trajectory; method=stride_run, selection="protein")
    @test size(ssmap) == (76, 26)
    @test sum(ssmap) == 10577
end

"""
    ss_class_content(
        f::Function, simulation::Simulation; 
        selection::Union{AbstractString,Function}=isprotein,
        ss_method::Function=stride_run,
        show_progress=true
    )

Calculates the secondary structure content of a specific class in a simulation. 
`f` is the function that returns,
for each residue, if the secondary structure is of a certain class. For example, to calculate 
the alpha helix content, use `f = is_alphahelix`.

The `selection` keyword argument can be used to choose a different selection, 
for example, `selection="protein and chain A"`.

The `ss_method` keyword argument can be used to choose the secondary structure prediction method,
either `stride_run` or `dssp_run`. The default is `stride_run`.

The `show_progress` keyword argument controls whether a progress bar is shown.

"""
function ss_class_content(
    f::F,
    simulation::Simulation;
    selection::Union{AbstractString,Function}=PDBTools.isprotein,
    ss_method::G=stride_run,
    show_progress=true
) where {F<:Function,G<:Function}
    sel = select(atoms(simulation), selection)
    ss_class_content = zeros(Float64, length(simulation))
    p = Progress(length(simulation); enabled=show_progress)
    for (iframe, frame) in enumerate(simulation)
        ss = ss_frame!(sel, frame; ss_method)
        ss_class_content[iframe] = count(f, ss) / max(1, length(ss))
        next!(p)
    end
    return ss_class_content
end

"""
    ss_class_content(f::F, ssmap::AbstractMatrix{<:Integer})

Calculates the secondary structure content of the trajectory, given the precomputed
secondary structure map. `f` is the function that returns,
for each residue, if the secondary structure is of a certain type. For example, to calculate
the alpha helix content, use `f = is_alphahelix`.

Here, `ssmap` is the secondary structure map of the trajectory, as returned by the `ss_map` function.

"""
function ss_class_content(
    f::F,
    ssmap::AbstractMatrix{Int},
) where {F<:Function}
    ss_class_content = zeros(Float64, size(ssmap, 2))
    for (iframe, ss_frame) in enumerate(eachcol(ssmap))
        ss_class_content[iframe] = count(f, ss_frame) / max(1, length(ss_frame))
    end
    return ss_class_content
end

"""
    ssmap_frame(ssmap::AbstractMatrix{Int}, iframe::Int}

Calculates the secondary structure composition of a frame of the simulation,
given the secondary structure map. 

Returns a dictionary of the secondary structure types and their counts, for the chosen frame.

"""
function ssmap_frame(ssmap::AbstractMatrix{Int}, iframe::Int)
    ssmap_frame= Dict{String,Int}()
    sscodes = @view(ssmap[:, iframe])
    for sscode in ss_code_to_number.(keys(ss_classes))
        ssmap_frame[class(sscode)] = count(==(sscode), sscodes)
    end
    return ssmap_frame
end

@testitem "ss_content/ss_composition" begin
    import ProteinSecondaryStructures.Testing
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