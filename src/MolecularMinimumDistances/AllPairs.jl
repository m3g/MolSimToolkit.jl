struct AllPairs{T,F1<:Function,F2<:Function} <: SystemPairs
    system::T
    xmol_indices::F1
    ymol_indices::F2
end

import Base.show
function Base.show(io::IO, ::MIME"text/plain", sys::AllPairs)
    print(io,chomp("""
    AllPairs system with:

    Number of atoms of first set: $(length(sys.system.xpositions))
    Number of molecules in first set: $(_number_of_molecules(sys.xmol_indices, sys.system.xpositions))

    Number of atoms of second set: $(length(sys.system.ypositions))
    Number of molecules in second set: $(_number_of_molecules(sys.ymol_indices, sys.system.ypositions))

    Cutoff: $(sys.system.cutoff)
    unitcell: [$(join(_uround.(sys.system.unitcell),", "))]
    """))
end

"""
    AllPairs(;
        xpositions::AbstractVector{<:AbstractVector{<:Real}},
        ypositions::AbstractVector{<:AbstractVector{<:Real}},
        cutoff::Real,
        unitcell::AbstractVecOrMat,
        xn_atoms_per_molecule::Integer,
        yn_atoms_per_molecule::Integer,
        parallel::Bool=true
    ) 

Initializes a particle system for the calculation of minimum distances
between one molecule and a set of other molecules. Returns a list 
minimum distances (`MinimumDistance` type), containing for each
molecule of the set the information about the closest distance to the
reference molecule.

Instead of the number of atoms per molecule, the user can also provide a 
more general `xmol_indices` and/or `ymol_indices` functions, 
which, for each atomic index, returns the corresponding
molecular index (which is `mol_indices(i) = (i-1)%n + 1` where `n` is the
number of atoms per molecule if all molecules have the same number of
atoms and are continously stored in the array of positions). 

# Examples

```julia-repl
julia> using MolSimToolkit, StaticArrays

julia> sys = AllPairs(
           xpositions=rand(SVector{3,Float64},10^5), # "solvent" (set of molecules)
           ypositions=rand(SVector{3,Float64},1000), # "solute" (target structure)
           cutoff=0.1,
           unitcell=[1,1,1],
           xn_atoms_per_molecule=5, # of the "solvent"
           yn_atoms_per_molecule=10 # of the "solvent"
       )
AllPairs system with:

Number of atoms of first set: 100000
Number of molecules in first set: 20000

Number of atoms of second set: 1000
Number of molecules in second set: 100

Cutoff: 0.1
unitcell: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]

julia> minimum_distances!(sys)[1]
20000-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(true, 1, 859, 0.037219109441123784)
 MinimumDistance{Float64}(true, 10, 117, 0.042183794688796634)
 ⋮
 MinimumDistance{Float64}(true, 99996, 168, 0.014269620784984633)
```

"""
function AllPairs(;
    xpositions::AbstractVector{<:AbstractVector{<:Real}},
    ypositions::AbstractVector{<:AbstractVector{<:Real}},
    cutoff::Real,
    unitcell::AbstractVecOrMat,
    xn_atoms_per_molecule::Union{Nothing,Integer}=nothing,
    yn_atoms_per_molecule::Union{Nothing,Integer}=nothing,
    xmol_indices::F1=nothing,
    ymol_indices::F2=nothing,
    parallel::Bool=true
) where {F1<:Union{Nothing,Function}, F2<:Union{Nothing,Function}}
    xmol_indices = _get_mol_indices(xmol_indices, xn_atoms_per_molecule; flag="x")
    ymol_indices = _get_mol_indices(ymol_indices, yn_atoms_per_molecule; flag="y")
    system = ParticleSystem(;
        xpositions=xpositions,
        ypositions=ypositions,
        cutoff=cutoff,
        unitcell=unitcell,
        output=(init_list(xpositions, xmol_indices), init_list(ypositions, ymol_indices)),
        output_name=:minimum_distances,
        parallel=parallel
    )
    return AllPairs(system, xmol_indices, ymol_indices)
end

#=
    function update_list!(
        i, j, d2,
        mol_indices_i::Fi,
        mol_indices_j::Fj,
        lists::Tuple{T,T}
    ) where {Fi<:Function, Fj<:Function, T<:AbstractVector{<:MinimumDistance}}


Function to update the minimum distance in the case where two lists are being constructed.

=#
function update_list!(i, j, d2, lists, system::AllPairs)
    x_list = lists[1]
    y_list = lists[2]
    d = sqrt(d2)
    imol = system.xmol_indices(i)
    if d < x_list[imol].d
        x_list[imol] = MinimumDistance(true, i, j, d)
    end
    jmol = system.ymol_indices(j)
    if d < y_list[jmol].d
        y_list[jmol] = MinimumDistance(true, j, i, d)
    end
    return lists
end

@testitem "MD - AllPairs" begin
    using PDBTools
    atoms = readPDB(MolSimToolkit.Testing.namd_pdb)
    popc = selindex(atoms, "resname POPC")
    protein = selindex(atoms, "protein")
    simulation = Simulation(
        MolSimToolkit.Testing.namd_pdb,
        MolSimToolkit.Testing.namd_traj,
    )
    first_frame!(simulation)
    coor = positions(current_frame(simulation))
    uc = unitcell(current_frame(simulation))
    xsolvent = zeros(eltype(coor), length(popc))
    xsolute = zeros(eltype(coor), length(protein))
    sys = AllPairs(
        xpositions = xsolvent,
        ypositions = xsolute,
        cutoff = 6.0,
        unitcell = uc, 
        xn_atoms_per_molecule = 134,
        yn_atoms_per_molecule = length(protein),
    )
    xmd_min = zeros(length(simulation))
    ymd_indices = zeros(Int, length(simulation))
    for (iframe, frame) in enumerate(simulation)
        pos = positions(frame)
        sys.xpositions .= @view(pos[popc])
        sys.ypositions .= @view(pos[protein])
        sys.unitcell = unitcell(frame)
        md = minimum_distances!(sys)
        xmd_min[iframe] = minimum(p -> p.d, md[1])
        ymd_indices[iframe] = minimum(p -> p.i, md[2])
    end
    @test xmd_min ≈ [1.9533253817007286, 1.6489800900283895, 1.6469793433779658, 1.5597538784005365, 1.5163768980822412] 
    @test ymd_indices ≈ [378, 597, 597, 577, 576]
end


