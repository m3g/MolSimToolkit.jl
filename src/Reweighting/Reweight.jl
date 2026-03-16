"""
Structure that contains the result of the reweighting analysis of the sequence. 

`probability` is a vector that contains the normalized weight of each frame in the simulation after applying some perturbation. 

`relative_probability` is a vector that contains the weight of each frame in the simulation relative to the original one after applying some perturbation.

`energy` is a vector that contains the energy difference for each frame in the simulation after applying some perturbation.
"""
struct ReweightResults{T<:Real}
    probabilities::Vector{Vector{T}}
    relative_probabilities::Vector{Vector{T}}
    energies::Vector{Vector{T}}
    distances::Vector{Int}
    entropy_diff::Vector{T}
    scaling::Vector{<:Real}
end

struct Perturbation{F<:Function}
    subgroup1::Vector{Int}
    subgroup2::Vector{Int}
    perturbation_function::F
    scaling::Vector{<:Real}
end

struct SystemPerturbations
    group1::Vector{Int}
    number_atoms_group1::Int
    group2::Vector{Int}
    number_atoms_group2::Int
    perturbations::OrderedCollections.OrderedDict{Any, Perturbation}
end

struct SystemPerturbationsOneGroup
    group1::Vector{Int}
    number_atoms_group1::Int
    perturbations::OrderedCollections.OrderedDict{Any, Perturbation}
end

####CellListMap configuration
mutable struct EnergyAndDistances
    energy::Float64
    distances::Int64
end

CellListMap.copy_output(x::EnergyAndDistances) = EnergyAndDistances(copy(x.energy), copy(x.distances))

function CellListMap.reset_output!(output::EnergyAndDistances)
    output.energy = 0.0
    output.distances = 0
    return output
end

function CellListMap.reducer(x::EnergyAndDistances, y::EnergyAndDistances)
    x.energy += y.energy
    x.distances += y.distances
    return EnergyAndDistances(x.energy, x.distances)
end
####

function Perturbation(
    atoms, 
    subgroup1::Union{String, Function}, 
    subgroup2::Union{String, Function}, 
    perturbation_function::Function, 
    scaling::Vector{<:Real};
    one_gp::Bool = false
)

    v1 = PDBTools.selindex(atoms, subgroup1)
    v2 = PDBTools.selindex(atoms, subgroup2)
    if one_gp == false || v1 == v2 || (v1 != v2 && isdisjoint(v1, v2))
        return Perturbation(v1, v2, perturbation_function, scaling)
    else
        return error(
            """
            Selected subgroups have overlap of atoms but they are not equal. 
            Subgroups must be either equal or completly different.
        """) 
    end
end

    """
    reweight(
        simulation::Simulation, 
        f_perturbation::Function, 
        group    for i in eachindex(perturbations)
        perturbations[i][1] = findall(i1 -> i1 in PDBTools.select(simulation.atoms, perturbations[i][1]), simulation.atoms[gp1]) 
        perturbations[i][2] = findall(i2 -> i2 in PDBTools.select(simulation.atoms, perturbations[i][2]), simulation.atoms[gp2])
    endist::Bool = false,
        cutoff::Real = 12.0, 
        k::Real = 1.0, 
        T::Real = 1.0
    )
    reweight(
        simulation::Simulation, 
        f_perturbation::Function, 
        group_1::AbstractVector{<:Integer},
        n_atoms_per_molecule_gp_1::Int,
        group_2::AbstractVector{<:Integer},
        n_atoms_per_molecule_gp_2::Int;
        all_dist::Bool = false,
        cutoff::Real = 12.0, 
        k::Real = 1.0,
        T::Real = 1.0
    )

Function that calculates the energy difference when a perturbation is applied on the system.

It returns "ReweightResults" structure that contains three results: `probability`, `relative_probability` and `energy` vectors.

The function needs a MolSimToolkit's simulation object, another function to compute the perturbation, and one or two types of atoms.

Additionally, you can also define a cutoff distance, the constant "k" (in some cases, it will be Boltzmann constant) and the temperature, T, of the system.

group_1 (and group_2) is a vector of atoms indexes, extracted, for example, from a pdb file. 

n_atoms_per_molecule is the number of atoms per molecules of each group

By default, only minimum distances are computed in order to perturb the system.
Using all_dist option allows the computation of all possible distances between the selected group of atoms

# Example

```julia-repl
julia> import PDBTools

julia> using MolSimToolkit, MolSimToolkit.Resampling

julia> simulation = Simulation("$testdir/Testing_reweighting.pdb", "/$testdir/Testing_reweighting_10_frames_trajectory.xtc")
Simulation 
    Atom type: Atom
    PDB file: /home/lucasv/.julia/dev/MolSimToolkit/src/Reweighting/test/Testing_reweighting.pdb
    Trajectory file: /home/lucasv/.julia/dev/MolSimToolkit/src/Resampling/test/Testing_reweighting_10_frames_trajectory.xtc
    Total number of frames: 10
    Frame range: 1:1:10
    Number of frames in range: 10
    Current frame: nothing

julia> i1 = PDBTools.selindex(atoms(simulation), "index 97 or index 106")
2-element AbstractVector{<:Integer}:
  97
 106

julia> i2 = PDBTools.selindex(atoms(simulation), "residue 15 and name HB3")
1-element AbstractVector{<:Integer}:
 171

julia> sum_of_dist = reweight(simulation, r -> r, i1, 2, i2, 1; all_dist = true, cutoff = 25.0)
-------------
FRAME WEIGHTS
-------------

Average probability = 0.1
standard deviation = 0.011364584999859616

-------------------------------------------
FRAME WEIGHTS RELATIVE TO THE ORIGINAL ONES
-------------------------------------------

Average probability = 0.6001821184861403
standard deviation = 0.06820820700931557

----------------------------------
COMPUTED ENERGY AFTER PERTURBATION
----------------------------------

Average energy = 0.5163045415662408
standard deviation = 0.11331912115883522

julia> sum_of_dist.energy
10-element Vector{Real}:
 17.738965476707595
 15.923698293115915
 17.16614676290554
 19.33003841107648
 16.02329229247863
 19.639005665480983
 35.73986006775934
 21.88798265022823
 20.66180657974777
 16.845109623700647

This result is the energy difference between the  perturbed frame and the original one. In this case, it is the sum of distances between the reffered atoms
```
"""

import Base.show
import Statistics
import StatsBase.mode
function Base.show(io::IO, mime::MIME"text/plain", res::ReweightResults)
    print(
    io,
    """
    -------------
    FRAME WEIGHTS
    -------------

    Average probability = $(Statistics.mean.(res.probabilities))
    standard deviation = $(Statistics.std.(res.probabilities)/sqrt(length(first(res.probabilities))))

    ----------------------------------
    AVERAGE PERTURBED ENERGIES
    ----------------------------------

    Average Energies = $(Statistics.mean.(res.energies))
    Standard Deviations = $(Statistics.std.(res.energies)/sqrt(length(first(res.energies))))

    -------------------------------------------
    ENTROPY LOSS DUE TO REWEIGHTING
    -------------------------------------------
    
    Entropy = $(res.entropy_diff)

    """)
end

#Function for two molecular entities
function reweight(
    simulation::Simulation, #Simulation object
    pert_input::SystemPerturbations; #Data structure containing molecular entities and their subgroups whose distances will be perturbed
    all_distances::Bool = false, #Flag to use CellListMap, computing all possible distances between the entities besides minimum distances
    k::Real = 1.0, #Boltzmann constant
    T::Real = 1.0, #Temperature
    cutoff::Real = 12.0, #Cutoff of distances
    tol::Real = 1.e-16 #Tolerance to consider a distance contribution
)

    #Defining results
    raw_output = OrderedCollections.OrderedDict{Any, Vector{AbstractVector}}(k =>
        [
            [[zeros(length(simulation))] for δ in pert_input.perturbations[k].scaling],
            [[zeros(length(simulation))] for δ in pert_input.perturbations[k].scaling],
            zeros(length(simulation)),
            zeros(Int, length(simulation)),        
        ] 
        for k in keys(pert_input.perturbations)
    )

    #Number of atoms per molecule
    n_molecules_gp1 = length(pert_input.group1) ÷ pert_input.number_atoms_group1
    computed_energy = 0

    #Defining function if CellListMap option is activated
    function energy_and_distances!(i,j,d, subgroup1, subgroup2, pert_func::Function, out::EnergyAndDistances)
        if is_in(subgroup1, pert_input.group1[i]) && is_in(subgroup2, pert_input.group2[j])
            out.energy += pert_func(d)
            if abs(out.energy) >= tol
                out.distances += 1
            end
        end
        return out
    end

    #Performing computation for every frame
    @showprogress for (iframe, frame) in enumerate(simulation)
        coordinates = positions(frame)
        gp_1_coord = coordinates[pert_input.group1]
        gp_2_coord = coordinates[pert_input.group2]
        uc = unitcell(frame)
        if all_distances
            system = ParticleSystem(
                xpositions = gp_1_coord,
                ypositions = gp_2_coord,
                unitcell = uc.orthorhombic ? diag(uc.matrix) : uc.matrix,
                cutoff = cutoff,
                output = EnergyAndDistances(0.0, 0),
                output_name = :energy_and_distances
            )
            for pk in keys(pert_input.perturbations)
                map_pairwise!(
                    (x, y, i, j, d2, output) -> energy_and_distances!(
                        i, 
                        j, 
                        sqrt(d2), 
                        pert_input.perturbations[pk].subgroup1, 
                        pert_input.perturbations[pk].subgroup2, 
                        pert_input.perturbations[pk].perturbation_function,
                        output
                    ),
                    system
                )
                raw_output[pk][3][iframe] = system.energy_and_distances.energy
                raw_output[pk][4][iframe] = system.energy_and_distances.distances
                system.energy_and_distances.energy = 0.0
                system.energy_and_distances.distances = 0
            end
        else
            for mol_ind in 1:n_molecules_gp1
                gp_1_coord_ref = gp_1_coord[(mol_ind - 1) * pert_input.number_atoms_group1 + 1 : mol_ind * pert_input.number_atoms_group1]
                gp_1_list, gp_2_list = minimum_distances(
                    xpositions = gp_1_coord_ref,
                    ypositions = gp_2_coord,
                    xn_atoms_per_molecule = pert_input.number_atoms_group1,
                    yn_atoms_per_molecule = pert_input.number_atoms_group2,
                    unitcell = uc.orthorhombic ? diag(uc.matrix) : uc.matrix,
                    cutoff = cutoff
                )
                for pk in keys(pert_input.perturbations)
                    for d_i in eachindex(gp_2_list)                    
                        if gp_2_list[d_i].within_cutoff && is_in(pert_input.perturbations[pk].subgroup2, pert_input.group2[gp_2_list[d_i].i]) && is_in(pert_input.perturbations[pk].subgroup1, pert_input.group1[gp_2_list[d_i].j])
                            computed_energy = pert_input.perturbations[pk].perturbation_function(gp_2_list[d_i].d)
                            raw_output[pk][3][iframe] += computed_energy
                            raw_output[pk][4][iframe] += abs(computed_energy) >= tol ? 1 : 0
                        end
                    end
                    computed_energy = 0
                end
            end
        end
    end
    for pk in keys(raw_output)
        raw_output[pk][3] = [δ * raw_output[pk][3] for δ in pert_input.perturbations[pk].scaling]
        raw_output[pk][2] = [exp.(-raw_output[pk][3][δ]/(k*T)) for δ in eachindex(raw_output[pk][3])]
        raw_output[pk][1] = [raw_output[pk][2][δ]/sum(raw_output[pk][2][δ]) for δ in eachindex(raw_output[pk][2])]
    end
    output = OrderedCollections.OrderedDict{Any, ReweightResults}(kys => 
        ReweightResults(
            raw_output[kys][1], 
            raw_output[kys][2],
            raw_output[kys][3],  
            raw_output[kys][4],
            [-k*(sum(raw_output[kys][1][δ] .* log.(raw_output[kys][1][δ])) - log(1/length(simulation))) for δ in eachindex(pert_input.perturbations[kys].scaling)],
            pert_input.perturbations[kys].scaling
        ) 
        for kys in keys(pert_input.perturbations)
    )
    return output
end

#Function for one molecular entity
function reweight(
    simulation::Simulation, #Simulation object
    pert_input::SystemPerturbationsOneGroup; #Data structure containing molecular entities and their subgroups whose distances will be perturbed
    all_distances::Bool = false, #Flag to use CellListMap, computing all possible distances between the entities besides minimum distances
    k::Real = 1.0, #Boltzmann constant
    T::Real = 1.0, #Temperature
    cutoff::Real = 12.0, #Cutoff of distances
    tol::Real = 1.e-32 #Tolerance to consider a distance contribution
)
    #Case vector for min distances
    cases = OrderedCollections.OrderedDict{Any, Vector{<:Real}}(k => 
        zeros(4)
        for k in keys(pert_input.perturbations)
    )

    for k in keys(pert_input.perturbations)
        cases[k] = check_vecs(pert_input.perturbations[k].subgroup1, pert_input.perturbations[k].subgroup2)
    end

    #Defining results
    raw_output = OrderedCollections.OrderedDict{Any, Vector{AbstractVector}}(k =>
        [
            [[zeros(length(simulation))] for δ in pert_input.perturbations[k].scaling],
            [[zeros(length(simulation))] for δ in pert_input.perturbations[k].scaling],
            zeros(length(simulation)),
            zeros(length(simulation)),        
        ] 
        for k in keys(pert_input.perturbations)
    )

    #Number of atoms per molecule
    n_molecules_gp1 = length(pert_input.group1) ÷ pert_input.number_atoms_group1
    computed_energy = 0

    #division factor for double counting
    dvf = 1

    #Defining function if CellListMap option is activated
    function energy_and_distances!(i,j,d, subgroup1, subgroup2, pert_func::Function, output::EnergyAndDistances)
        gp1 = pert_input.group1
        mol1 = (i - 1) ÷ pert_input.number_atoms_group1 + 1
        mol2 = (j - 1) ÷ pert_input.number_atoms_group1 + 1
        if (is_in(subgroup1, gp1[i]) && is_in(subgroup2, gp1[j]) && mol1 != mol2) || (is_in(subgroup1, gp1[j]) && is_in(subgroup2, gp1[i]) && mol1 != mol2)
            output.energy += pert_func(d)
            if abs(output.energy) >= tol
                output.distances += 1
            end
        end
        return output
    end

    #Performing computation for every frame
    @showprogress for (iframe, frame) in enumerate(simulation)
        coordinates = positions(frame)
        gp_coord = coordinates[pert_input.group1]
        uc = unitcell(frame)
        if all_distances
            system = ParticleSystem(
                xpositions = gp_coord,
                unitcell = uc.orthorhombic ? diag(uc.matrix) : uc.matrix,
                cutoff = cutoff,
                output = EnergyAndDistances(0.0, 0),
                output_name = :energy_and_distances
            )
            for pk in keys(pert_input.perturbations)
                map_pairwise!(
                    (x, y, i, j, d2, output) -> energy_and_distances!(
                        i, 
                        j, 
                        sqrt(d2), 
                        pert_input.perturbations[pk].subgroup1, 
                        pert_input.perturbations[pk].subgroup2, 
                        pert_input.perturbations[pk].perturbation_function,
                        output
                    ),
                    system
                )
                raw_output[pk][3][iframe] = system.energy_and_distances.energy
                raw_output[pk][4][iframe] = system.energy_and_distances.distances
                system.energy_and_distances.energy = 0.0
                system.energy_and_distances.distances = 0
            end
        else
            for mol_ind in 1:n_molecules_gp1
                i_index = (mol_ind - 1) * pert_input.number_atoms_group1 + 1
                f_index = mol_ind * pert_input.number_atoms_group1
                gp_coord_ref = gp_coord[i_index : f_index]
                gp_1_list, gp_2_list = minimum_distances(
                    xpositions = gp_coord_ref,
                    ypositions = gp_coord,
                    xn_atoms_per_molecule = pert_input.number_atoms_group1,
                    yn_atoms_per_molecule = pert_input.number_atoms_group1,
                    unitcell = uc.orthorhombic ? diag(uc.matrix) : uc.matrix,
                    cutoff = cutoff
                )
                for pk in keys(pert_input.perturbations)
                    sg1 = pert_input.perturbations[pk].subgroup1
                    sg2 = pert_input.perturbations[pk].subgroup2
                    for d_i in eachindex(gp_2_list)
                        if gp_2_list[d_i].within_cutoff && is_in(collect(i_index:1:f_index), gp_2_list[d_i].i) == false && is_in(sg2, pert_input.group1[gp_2_list[d_i].i]) && is_in(sg1, pert_input.group1[gp_2_list[d_i].j])
                            if cases[pk][1] == 1
                                dvf = 2
                            end
                            computed_energy = pert_input.perturbations[pk].perturbation_function(gp_2_list[d_i].d)
                            raw_output[pk][3][iframe] += computed_energy / dvf
                            raw_output[pk][4][iframe] += abs(computed_energy) >= tol ? 1 / dvf : 0
                        end
                    end
                    computed_energy = 0
                    dvf = 1
                end
            end
        end
    end
    for pk in keys(raw_output)
        raw_output[pk][3] = [δ * raw_output[pk][3] for δ in pert_input.perturbations[pk].scaling]
        raw_output[pk][2] = [exp.(-raw_output[pk][3][δ]/(k*T)) for δ in eachindex(raw_output[pk][3])]
        raw_output[pk][1] = [raw_output[pk][2][δ]/sum(raw_output[pk][2][δ]) for δ in eachindex(raw_output[pk][2])]
    end
    output = OrderedCollections.OrderedDict{Any, ReweightResults}(kys => 
        ReweightResults(
            raw_output[kys][1], 
            raw_output[kys][2],
            raw_output[kys][3],  
            Int.(raw_output[kys][4]),
            [-k*(sum(raw_output[kys][1][δ] .* log.(raw_output[kys][1][δ])) - log(1/length(simulation))) for δ in eachindex(pert_input.perturbations[kys].scaling)],
            pert_input.perturbations[kys].scaling
        ) 
        for kys in keys(pert_input.perturbations)
    )
    return output
end

#Checking if PDB file and input match
function check_input(simulation::Simulation, atom_group::Union{String, Function}, n_mol_per_group::Int, name::String, debug::Bool)
    check = length(unique(PDBTools.residue.(PDBTools.select(simulation.atoms, atom_group)))) #check number of residues (number of molecules) in PDB
    division = check ÷ n_mol_per_group
    quotient = mod(check, n_mol_per_group)
    if debug
        if ((division == 1 && quotient == 0) || n_mol_per_group == 1)
            println("Number of molecules in the system seems to be correct for $name ($n_mol_per_group)")
        elseif division != 1 && quotient == 0
            return @warn("""
                Number of molecules in the system ($check) seems to be a multiple of your input for $name ($n_mol_per_group).
                """)
        elseif quotient != 0
            return @warn("""
                The number of residues ($check) in the PDB file for $name is different than the number of molecules based on the number on the input ($n_mol_per_group) 
                """)
        end
    end
end

function is_in(x, i)
    j = searchsortedfirst(x, i)
    j > length(x) && return false
    if x[j] == i
        return true
    else
        return false
    end
end

function check_vecs(v1, v2)
    if v1 == v2
        return [1,0]
    end
    return [0,1]
end