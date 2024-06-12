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
    for (iat, at) in pairs(atoms)
        iatom = PDBTools.index(at)
        atoms[iat].x = coordinates[iatom].x
        atoms[iat].y = coordinates[iatom].y
        atoms[iat].z = coordinates[iatom].z
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

## Example

```jldoctest
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> ssmap = ss_map(simulation; selection="residue >= 30 and residue <= 35", show_progress=false)
6×5 Matrix{Int64}:
 5  9  5  5  5
 5  9  5  5  5
 5  1  5  5  5
 5  1  5  5  5
 5  1  5  5  5
 9  9  9  9  9

```

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
    ssmap = ss_map(simulation; selection="residue >= 30 and residue <= 35")
    @test sum(ssmap) == 166
    # With DSSP
    ssmap = ss_map(simulation; method=dssp_run)
    @test size(ssmap) == (43, 5)
    @test sum(ssmap) == 1006
end

"""
    ss_class_mean(ss_class, ssmap::AbstractMatrix{<:Integer}; dims=nothing)

Calculates the mean secondary structure class content of the trajectory, given the secondary structure map.

The mean can be calculated along the residues (default) or along the frames, by setting the `dims` keyword argument.

- `dims=nothing` (default) calculates the mean occurence of `ss_class` of the whole matrix.
- `dims=1` calculates the mean occurence of `ss_class` along the frames, for each residue.
- `dims=2` calculates the mean occurence of `ss_class` along the residues, for each frame.

`ss_class` can be either a string, a character, or an integer, setting the class of secondary structure
to be consdiered. For example, for `alpha helix`, use "H".  

The classes can be found in the ProteinSecondaryStructures.jl package documentation, at:

https://m3g.github.io/ProteinSecondaryStructures.jl/stable/explanation/#Secondary-structure-classes

## Example

```jldoctest 
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> ssmap = ss_map(simulation; # 5 frames 
                   selection="residue >= 30 and residue <= 35", # 6 residues
                   show_progress=false
               );

julia> ss_class_mean("C", ssmap)
0.23333333333333334

julia> ss_class_mean("C", ssmap; dims=1) # mean coil content per residue
6-element Vector{Float64}:
 0.2
 0.2
 0.0
 0.0
 0.0
 1.0 

julia> ss_class_mean("C", ssmap; dims=2) # mean coil content per frame
5-element Vector{Float64}:
 0.16666666666666666
 0.5
 0.16666666666666666
 0.16666666666666666
 0.16666666666666666

```

"""
function ss_class_mean(
    ss_class::Union{AbstractString,AbstractChar,Integer},
    ssmap::AbstractMatrix{Int};
    dims::Union{Nothing,Int}=nothing,
)
    ss_class = ss_class isa Integer ? ss_class : ss_code_to_number(ss_class)
    if dims == 1
        ss_class_mean = zeros(Float64, size(ssmap, 1))
        for (ires, res) in enumerate(eachrow(ssmap))
            ss_class_mean[ires] = count(==(ss_class), res) / max(1, length(res))
        end
    elseif dims == 2
        ss_class_mean = zeros(Float64, size(ssmap, 2))
        for (iframe, ss_frame) in enumerate(eachcol(ssmap))
            ss_class_mean[iframe] = count(==(ss_class), ss_frame) / max(1, length(ss_frame))
        end
    elseif isnothing(dims)
        ss_class_mean = count(==(ss_class), ssmap) / max(1, length(ssmap))
    end
    return ss_class_mean
end

@testitem "secondary structure" begin
    using MolSimToolkit
    using MolSimToolkit.Testing
    simulation = Simulation(Testing.namd_pdb, Testing.namd_traj)
    # With stride
    ssmap = ss_map(simulation; ss_method=stride_run, show_progress=false)
    @test size(ssmap) == (43, 5)
    @test sum(ssmap) == 842
    helical_content = ss_class_content("H", ssmap)
    @test length(helical_content) == 5
    @test sum(helical_content) / length(helical_content) ≈ 0.6093023255813954
    @test ss_frame(ssmap, 3) == Dict("310 helix" => 0, "bend" => 0, "turn" => 10, "beta bridge" => 0, "kappa helix" => 0, "pi helix" => 0, "beta strand" => 0, "alpha helix" => 25, "loop" => 0, "coil" => 8)
    # With DSSP
    ssmap = ss_map(simulation; ss_method=dssp_run, show_progress=false)
    @test size(ssmap) == (43, 5)
    @test sum(ssmap) == 1006
    helical_content = ss_class_content("B", ssmap)
    @test length(helical_content) == 26
    @test sum(helical_content) / length(helical_content) ≈ 0.21204453441295545
    # With non-contiguous indexing
    atoms = readPDB(pdbfile, only=at -> (10 <= at.residue < 30) | (40 <= at.residue < 60))
    helical_content = ss_content(is_anyhelix, atoms, trajectory; method=stride_run, show_progress=false)
    @test length(helical_content) == 26
    @test sum(helical_content) / length(helical_content) ≈ 0.20288461538461539
end