#
# Computes the secondary structure in a single frame
# This is an internal function.
#
const replace_aa_code = Dict{String,String}(
    "HSD" => "HIS",
    "HSP" => "HIS",
    "HSE" => "HIS",
    "HID" => "HIS",
    "HIP" => "HIS",
    "GLUP" => "GLU",
    "ASPP" => "ASP",
)
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
    # This is very bad: we are writting temporary files twice for each frame
    # one here, one in the ss_method function, to adjust the PDB header
    tmpfile = tempname()*".pdb"
    PDBTools.writePDB(atoms, tmpfile)
    return ss_method(tmpfile; adjust_pdb=true)
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
selection. The `PDBTools` selection syntax can be used, for example `selection="protein and chain A"`, 
or general Julia functions, like `selection=at -> chain(at) in ('A', 'B')`.

The `ss_method` keyword argument can be used to choose the secondary structure prediction method,
which can be either `stride_run` or `dssp_run`. The default is `stride_run`. STRIDE is a faster
algorithm, while DSSP is the default one in PDB database.

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

julia> ss_name.(ssmap)
6×5 Matrix{String}:
 "turn"  "coil"       "turn"  "turn"  "turn"
 "turn"  "coil"       "turn"  "turn"  "turn"
 "turn"  "310 helix"  "turn"  "turn"  "turn"
 "turn"  "310 helix"  "turn"  "turn"  "turn"
 "turn"  "310 helix"  "turn"  "turn"  "turn"
 "coil"  "coil"       "coil"  "coil"  "coil"

```

"""
function ss_map(
    simulation::Simulation;
    selection::Union{AbstractString,Function}=PDBTools.isprotein,
    ss_method::F=stride_run,
    show_progress=true
) where {F<:Function}
    sel = PDBTools.select(atoms(simulation), selection)
    for iat in eachindex(sel)
        if haskey(replace_aa_code, PDBTools.resname(sel[iat]))
            sel[iat].resname = replace_aa_code[sel[iat].resname]
        end
    end
    ss_map = zeros(Int, length(PDBTools.eachresidue(sel)), length(simulation))
    p = Progress(length(simulation); enabled=show_progress)
    for (iframe, frame) in enumerate(simulation)
        ss = _ss_frame!(sel, frame; ss_method)
        for (i, ssdata) in pairs(ss)
            ss_map[i, iframe] = ss_number(ssdata)
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
    ssmap = ss_map(simulation; ss_method=stride_run, show_progress=false)
    @test size(ssmap) == (43, 5)
    @test sum(ssmap) == 842
    ssmap = ss_map(simulation; selection="residue >= 30 and residue <= 35")
    @test sum(ssmap) == 166
    # With DSSP
    ssmap = ss_map(simulation; ss_method=dssp_run, show_progress=false)
    @test size(ssmap) == (43, 5)
    @test sum(ssmap) == 943
end

"""
    ss_mean(ssmap::AbstractMatrix{<:Integer}; class, dims=nothing)

Calculates the mean secondary structure class content of the trajectory, given the secondary structure map.

The secondary structure class to be considered must be defined by the `class` keyword argument.

`class` can be either a string, a character, or an integer, or a set of values, setting the class(es) 
of secondary structure to be consdiered. For example, for `alpha helix`, use "H". It can also be a vector of classes, 
such as `class=["H", "E"]`.

The mean can be calculated along the residues (default) or along the frames, by setting the `dims` keyword argument.

- `dims=nothing` (default) calculates the mean occurence of `ss_class` of the whole matrix.
- `dims=1` calculates the mean occurence of `ss_class` along the frames, for each residue.
- `dims=2` calculates the mean occurence of `ss_class` along the residues, for each frame.

The classes can be found in the ProteinSecondaryStructures.jl package documentation.

## Example

```jldoctest ;filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2***"
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> simulation = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> ssmap = ss_map(simulation; # 5 frames 
                   selection="residue >= 30 and residue <= 35", # 6 residues
                   show_progress=false
               );

julia> ss_mean(ssmap; class="C")
0.23333333333333334

julia> ss_mean(ssmap; class="C", dims=1) # mean coil content per residue
5-element Vector{Float64}:
 0.16666666666666666
 0.5
 0.16666666666666666
 0.16666666666666666
 0.16666666666666666 

julia> ss_mean(ssmap; class="C", dims=2) # mean coil content per frame
6-element Vector{Float64}:
 0.2
 0.2
 0.0
 0.0
 0.0
 1.0

julia> ss_mean(ssmap; class=["C", "T"]) # mean coil or turn
0.9

```

"""
function ss_mean(
    ssmap::AbstractMatrix{Int};
    class,
    dims::Union{Nothing,Int}=nothing,
)
    if first(class) isa Union{AbstractChar,AbstractString}
        class = ss_number.(class)
    end
    ss_mean = if isnothing(dims)
        mean(in(class), ssmap)
    else
        vec(mean(in(class), ssmap; dims))
    end
    return ss_mean
end

function ss_heatmap(args...; kargs...) 
    throw(ArgumentError("""

        The ss_heatmap function requires loading the Plots package.

    """))
end

@testitem "secondary structure" begin
    using MolSimToolkit
    using MolSimToolkit.Testing
    using Statistics: mean
    using PDBTools
    simulation = Simulation(Testing.namd_pdb, Testing.namd_traj)
    # With stride
    ssmap = ss_map(simulation; ss_method=stride_run, show_progress=false)
    @test size(ssmap) == (43, 5)
    @test sum(ssmap) == 842
    helical_content = ss_mean(ssmap; class="H")
    @test helical_content ≈ 0.6093023255813954
    h_per_frame = ss_mean(ssmap; class="H", dims=1)
    @test length(h_per_frame) == 5
    @test mean(h_per_frame) ≈ helical_content
    h_per_residue = ss_mean(ssmap; class="H", dims=2)
    @test length(h_per_residue) == 43
    @test mean(h_per_residue) ≈ helical_content
    # With DSSP
    ssmap = ss_map(simulation; ss_method=dssp_run, show_progress=false)
    @test size(ssmap) == (43, 5)
    @test sum(ssmap) == 943
    helical_content = ss_mean(ssmap; class="H")
    @test helical_content ≈ 0.5813953488372093
    h_per_frame = ss_mean(ssmap; class="H", dims=1)
    @test length(h_per_frame) == 5
    @test mean(h_per_frame) ≈ helical_content
    h_per_residue = ss_mean(ssmap; class="H", dims=2)
    @test length(h_per_residue) == 43
    @test mean(h_per_residue) ≈ helical_content
    # With non-contiguous indexing
    sel(at) = isprotein(at) && (10 <= at.residue < 30) | (40 <= at.residue < 60)
    ssmap = ss_map(simulation; selection=sel, show_progress=false)
    helical_content = ss_mean(ssmap; class="H")
    @test helical_content ≈ 0.78333333333333
    h_per_frame = ss_mean(ssmap; class="H", dims=1)
    @test length(h_per_frame) == 5
    @test mean(h_per_frame) ≈ helical_content
    h_per_residue = ss_mean(ssmap; class="H", dims=2)
    @test length(h_per_residue) == 24
    @test mean(h_per_residue) ≈ helical_content
end