struct SelfPairs{T,F<:Function} <: SystemPairs
    system::T
    mol_indices::F
end

import Base.show
function Base.show(io::IO, ::MIME"text/plain", sys::SelfPairs)
    print(io,chomp("""
    SelfPairs system with:

    Number of atoms: $(length(sys.system.xpositions))
    Cutoff: $(sys.system.cutoff)
    unitcell: [$(join(_uround.(sys.system.unitcell),", "))]
    Number of molecules: $(_number_of_molecules(sys.mol_indices, sys.system.xpositions))
    """))
end

"""
    SelfPairs(;
        xpositions::AbstractVector{<:AbstractVector{<:Real}},
        cutoff::Real,
        unitcell::AbstractVecOrMat,
        xn_atoms_per_molecule::Integer,
        parallel::Bool=true
    ) where T<:Real

Initializes a particle system for the calculation of minimum distances
within a single set of molecules. The shortest distance of each molecule
to any other molecule of the same set is computed.

Instead of the number of atoms per molecule, the user can also provide a 
more general `mol_indices` function, which, for each atomic index, returns the 
corresponding
molecular index (which is `mol_indices(i) = (i-1)%n + 1` where `n` is the
number of atoms per molecule if all molecules have the same number of
atoms and are continously stored in the array of positions). 

# Examples

```julia-repl
julia> using MolSimToolkit, StaticArrays

julia> sys = SelfPairs(
           xpositions=rand(SVector{3,Float64},10^5),
           cutoff=0.1,
           unitcell=[1,1,1],
           xn_atoms_per_molecule=5
       )
SelfPairs system with:

Number of atoms: 100000
Cutoff: 0.1
unitcell: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
Number of molecules: 20000

julia> minimum_distances!(sys)
20000-element Vector{MinimumDistance{Float64}}:
 MinimumDistance{Float64}(true, 4, 24930, 0.008039482961077074)
 MinimumDistance{Float64}(true, 6, 74055, 0.0049818659155905255)
 ⋮
 MinimumDistance{Float64}(true, 99999, 75403, 0.0025051670801269433)
```

"""
function SelfPairs(;
    xpositions::AbstractVector{<:AbstractVector{T}},
    cutoff::T,
    unitcell::AbstractVecOrMat,
    xn_atoms_per_molecule::Union{Nothing,Integer}=nothing,
    mol_indices::F1=nothing,
    parallel::Bool=true,
) where {T<:Real, F1<:Union{Nothing,Function}}
    mol_indices = _get_mol_indices(mol_indices, xn_atoms_per_molecule)
    system = ParticleSystem(;
        xpositions=xpositions,
        cutoff=cutoff,
        unitcell=unitcell,
        output=init_list(xpositions, mol_indices),
        output_name=:minimum_distances,
        parallel=parallel
    )
    return SelfPairs(system, mol_indices)
end

#
# This file containst the functions for single-sets, that is for those cases where
# the list of minimum-distances is between the molecules of a single component.
#
function update_list!(i, j, d2, list::List, system::SelfPairs)
    imol = system.mol_indices(i)
    jmol = system.mol_indices(j)
    if imol != jmol
        d = sqrt(d2)
        if d < list[imol].d
            list[imol] = MinimumDistance(true, i, j, d)
        end
        if d < list[jmol].d
            list[jmol] = MinimumDistance(true, j, i, d)
        end
    end
    return list
end

@testitem "MD - SelfPairs" begin
    using PDBTools
    atoms = readPDB(MolSimToolkit.Testing.namd_pdb)
    popc = selindex(atoms, "resname POPC")
    simulation = Simulation(
        MolSimToolkit.Testing.namd_pdb,
        MolSimToolkit.Testing.namd_traj,
    )
    first_frame!(simulation)
    coor = positions(current_frame(simulation))
    uc = unitcell(current_frame(simulation))
    xsolvent = zeros(eltype(coor), length(popc))
    sys = SelfPairs(
        xpositions = xsolvent,
        cutoff = 6.0,
        unitcell = uc, 
        xn_atoms_per_molecule = 134,
    )
    md_min = zeros(length(simulation))
    for (iframe, frame) in enumerate(simulation)
        pos = positions(frame)
        sys.xpositions .= @view(pos[popc])
        sys.unitcell = unitcell(frame)
        md = minimum_distances!(sys)
        md_min[iframe] = minimum(p -> p.d, md)
    end
    @test md_min ≈ [1.7787845332699295, 1.8000887689403775, 1.7751687089934818, 1.8329921680458754, 1.6533482891665536]
end
