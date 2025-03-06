struct CrossPairs{T,F<:Function} <: SystemPairs
    system::T
    mol_indices::F
end
CrossPairs(system::S, mol_indices::F) where {S,F} = CrossPairs{S,F,:minimum_distances}(system,mol_indices)

import Base.show
function Base.show(io::IO, ::MIME"text/plain", sys::CrossPairs)
    print(io,chomp("""
    CrossPairs system with:

    Number of atoms of set x: $(length(sys.system.xpositions))
    Number of molecules in set x: $(_number_of_molecules(sys.mol_indices, sys.system.xpositions))
    Number of atoms of target structure y: $(length(sys.system.ypositions))
    Cutoff: $(sys.system.cutoff)
    unitcell: [$(join(_uround.(sys.system.unitcell),", "))]
    """))
end

"""
    CrossPairs(;
        xpositions::AbstractVector{<:AbstractVector{<:Real}},
        ypositions::AbstractVector{<:AbstractVector{<:Real}},
        cutoff::Real,
        unitcell::AbstractVecOrMat,
        xn_atoms_per_molecule::Integer,
        parallel::Bool=true
    )

Initializes a particle system for the calculation of minimum distances
between one molecule and a set of other molecules. Returns a list 
minimum distances (`MinimumDistance` type), containing for each
molecule of the set the information about the closest distance to the
reference molecule.

Instead of the number of atoms per molecule, the user can also provide a 
more general `mol_indices` function, which, for each atomic index, returns the 
corresponding
molecular index (which is `mol_indices(i) = (i-1)%n + 1` where `n` is the
number of atoms per molecule if all molecules have the same number of
atoms and are continously stored in the array of positions). 

# Examples

```julia-repl
julia> using MolSimToolkit, StaticArrays

julia> sys = CrossPairs(
           xpositions=rand(SVector{3,Float64},10^5), # "solvent" (set of molecules)
           ypositions=rand(SVector{3,Float64},1000), # "solute" (target structure)
           cutoff=0.1,
           unitcell=[1,1,1],
           xn_atoms_per_molecule=5 # of the "solvent"
       )
CrossPairs system with:

Number of atoms of set: 100000
Number of atoms of target structure: 1000
Cutoff: 0.1
unitcell: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
Number of molecules in set: 20000

julia> minimum_distances!(sys)
20000-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(true, 1, 859, 0.037219109441123784)
 MinimumDistance{Float64}(true, 10, 117, 0.042183794688796634)
 ⋮
 MinimumDistance{Float64}(true, 99996, 168, 0.014269620784984633)
```

"""
function CrossPairs(;
    xpositions::AbstractVector{<:AbstractVector{<:Real}},
    ypositions::AbstractVector{<:AbstractVector{<:Real}},
    cutoff::Real,
    unitcell::AbstractVecOrMat,
    xn_atoms_per_molecule::Union{Nothing,Integer}=nothing,
    xmol_indices::F1=nothing,
    parallel::Bool=true
) where {F1<:Union{Nothing,Function}}
    xmol_indices = _get_mol_indices(xmol_indices, xn_atoms_per_molecule; flag="x")
    system = ParticleSystem(;
        xpositions=xpositions,
        ypositions=ypositions,
        cutoff=cutoff,
        unitcell=unitcell,
        output=init_list(xpositions, xmol_indices),
        output_name=:minimum_distances,
        parallel=parallel
    )
    return CrossPairs(system, xmol_indices)
end

#
# Functions that return the atoms of the second set that are closer to
# each molecule of the first set (only one list is returned).
#
function update_list!(i, j, d2, list, system::CrossPairs)
    d = sqrt(d2)
    imol = system.mol_indices(i)
    if d < list[imol].d
        list[imol] = MinimumDistance(true, i, j, d)
    end
    return list
end

@testitem "MD - CrossPairs" begin
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
    sys = CrossPairs(
        xpositions = xsolvent,
        ypositions = xsolute,
        cutoff = 6.0,
        unitcell = uc, 
        xn_atoms_per_molecule = 134,
    )
    md_count = zeros(Int, length(simulation))
    for (iframe, frame) in enumerate(simulation)
        pos = positions(frame)
        sys.xpositions .= @view(pos[popc])
        sys.ypositions .= @view(pos[protein])
        sys.unitcell = unitcell(frame)
        md = minimum_distances!(sys)
        md_count[iframe] = count(p -> p.within_cutoff, md)
    end
    @test md_count == [23, 24, 24, 25, 24]
end