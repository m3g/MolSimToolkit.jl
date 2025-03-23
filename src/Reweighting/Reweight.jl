"""
Structure that contains the result of the reweighting analysis of the sequence. 

`probability` is a vector that contains the normalized weight of each frame in the simulation after applying some perturbation. 

`relative_probability` is a vector that contains the weight of each frame in the simulation relative to the original one after applying some perturbation.

`energy` is a vector that contains the energy difference for each frame in the simulation after applying some perturbation.
"""
struct ReweightResults{T<:Real}
    probability::Vector{T}
    relative_probability::Vector{T}
    energy::Vector{T}
end

"""
    reweight(
        simulation::Simulation, 
        f_perturbation::Function, 
        group_1::AbstractVector{<:Integer},
        n_atoms_per_molecule::Int;
        all_dist::Bool = false,
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
function reweight(
    simulation::Simulation, 
    f_perturbation::Function, 
    gp::AbstractVector{<:Integer}, #Group atoms
    n_mol_gp::Int; #Number of molecules of the group
    all_dist::Bool = false, #Compute every possible distance option
    mol_contrib::Union{String, Function, Nothing} = nothing, #Which atoms of the group we are interested?
    cutoff::Real = 12.0, 
    k::Real = 1.0, 
    T::Real = 1.0
)
    #Defining results vector
    prob_vec = zeros(length(simulation))
    prob_rel_vec = zeros(length(simulation))
    energy_vec = zeros(length(simulation))
    
    #Number of atoms per molecule
    n_atoms_per_molecule_gp = length(gp) ÷ n_mol_gp

    #Checking if PDB file and input match
    check_n_mol(simulation, gp_1, n_mol_gp_1, "group")
    
    #Retrieving atoms from the contributions
    ind_contrib = retrieve_contrib_indexes(simulation, mol_contrib, gp)

    #Performing computation for every frame
    for (iframe, frame) in enumerate(simulation)
        coordinates = positions(frame)
        gp_1_coord = coordinates[group_1]
        system = nothing
        if all_dist
            system = ParticleSystem(
                xpositions = gp_1_coord,
                unitcell = unitcell(frame),
                cutoff = cutoff,
                output = 0.0,
                output_name = :total_energy
            )
            if isequal(ind_contrib, gp) == false
                if system[d_i].i in (ind_contrib) && system[d_i].j in (ind_contrib)
                    energy_vec[iframe] = map_pairwise!((x, y, i, j, d2, total_energy) -> total_energy + f_perturbation(sqrt(d2)), system)
                end
            else
                energy_vec[iframe] = map_pairwise!((x, y, i, j, d2, total_energy) -> total_energy + f_perturbation(sqrt(d2)), system)
            end
        else
            for mol_ind in 1:(length(gp_1_coord) ÷ n_atoms_per_molecule)
                gp_1_coord_ref = gp_1_coord[(mol_ind - 1) * n_atoms_per_molecule + 1 : mol_ind * n_atoms_per_molecule]
                gp_1_ref_list, gp_1_list = minimum_distances(
                    xpositions = gp_1_coord_ref,
                    ypositions = gp_1_coord,
                    xn_atoms_per_molecule = n_atoms_per_molecule,
                    yn_atoms_per_molecule = n_atoms_per_molecule,
                    unitcell = unitcell(frame),
                    cutoff = cutoff
                )
                for d_i in eachindex(gp_1_list)
                    if gp_1_list[d_i].within_cutoff && gp_1_list[d_i].d != 0
                        if isequal(ind_contrib, gp) == false
                            if gp_1_list[d_i].i in (ind_contrib) && gp_1_list[d_i].j in (ind_contrib)
                                energy_vec[iframe] += perturb_func(gp_2_list[d_i].d)
                            end
                        else
                            energy_vec[iframe] += perturb_func(gp_2_list[d_i].d)
                        end
                    end
                end
            end
        end
    end
    @. prob_rel_vec = exp(-energy_vec/(k*T))
    prob_vec = prob_rel_vec/sum(prob_rel_vec)
    output = ReweightResults(prob_vec, prob_rel_vec, energy_vec)
    return output
end
#two group of atoms method
function reweight(
    simulation::Simulation, #MolSimToolkit simulation object
    perturb_func::Function, #Perturbation function
    gp_1::AbstractVector{<:Integer}, #Group 1 atoms
    n_mol_gp_1::Int, #Number of molecules of group 1
    gp_2::AbstractVector{<:Integer},  #Group 2 atoms
    n_mol_gp_2::Int; #Number of molecules of group 2
    all_dist::Bool = false, #Compute every possible distance option
    mol_1_contrib::Union{String, Function, Nothing} = nothing, #Which atoms of group 1 we are interested?
    mol_2_contrib::Union{String, Function, Nothing} = nothing, #Which atoms of group 2 we are interested?
    cutoff::Real = 12.0,  #Cutoff distance for calculations
    k::Real = 1.0, #Boltzmann constat value
    T::Real = 1.0 #temperature of the system
)
    #Defining results vector
    prob_vec = zeros(length(simulation))
    prob_rel_vec = zeros(length(simulation))
    energy_vec = zeros(length(simulation))
    
    #Number of atoms per molecule
    n_atoms_per_molecule_gp_1 = length(gp_1) ÷ n_mol_gp_1
    n_atoms_per_molecule_gp_2 = length(gp_2) ÷ n_mol_gp_2
    
    #Checking if PDB file and input match
    check_n_mol(simulation, gp_1, n_mol_gp_1, "group 1")
    check_n_mol(simulation, gp_2, n_mol_gp_2, "group 2")
    
    #Retrieving atoms from the contributions
    ind_contrib_1 = retrieve_contrib_indexes(simulation, mol_1_contrib, gp_1)
    ind_contrib_2 = retrieve_contrib_indexes(simulation, mol_2_contrib, gp_2)
    
    #Performing computation for every frame
    for (iframe, frame) in enumerate(simulation)
        coordinates = positions(frame)
        gp_1_coord = coordinates[gp_1]
        gp_2_coord = coordinates[gp_2]
        if all_dist
            system = ParticleSystem(
                xpositions = gp_1_coord,
                ypositions = gp_2_coord,
                unitcell = unitcell(frame),
                cutoff = cutoff,
                output = 0.0,
                output_name = :total_energy
            )
            if isequal(ind_contrib_1, gp_1) == false && isequal(ind_contrib_2, gp_2) == false
                if system[d_i].i in (ind_contrib_1) && system[d_i].j in (ind_contrib_2)
                    energy_vec[iframe] = map_pairwise!((x, y, i, j, d2, total_energy) -> total_energy + f_perturbation(sqrt(d2)), system)
                end
            else
                energy_vec[iframe] = map_pairwise!((x, y, i, j, d2, total_energy) -> total_energy + perturb_func(sqrt(d2)), system)
            end
        else
            for mol_ind in 1:n_mol_gp_1
                gp_1_coord_ref = gp_1_coord[(mol_ind - 1) * n_atoms_per_molecule_gp_1 + 1 : mol_ind * n_atoms_per_molecule_gp_1]
                gp_1_list, gp_2_list = minimum_distances(
                    xpositions = gp_1_coord_ref,
                    ypositions = gp_2_coord,
                    xn_atoms_per_molecule = n_atoms_per_molecule_gp_1,
                    yn_atoms_per_molecule = n_atoms_per_molecule_gp_2,
                    unitcell = unitcell(frame),
                    cutoff = cutoff
                )
                for d_i in eachindex(gp_2_list)
                    if gp_2_list[d_i].within_cutoff && gp_2_list[d_i].d != 0
                        if isequal(ind_contrib_1, gp_1) == false && isequal(ind_contrib_2, gp_2) == false
                            if gp_2_list[d_i].i in (ind_contrib_2) && gp_2_list[d_i].j in (ind_contrib_1)
                                energy_vec[iframe] += perturb_func(gp_2_list[d_i].d)
                            end
                        else
                            energy_vec[iframe] += perturb_func(gp_2_list[d_i].d)
                        end
                    end
                end
            end
        end
    end
    @. prob_rel_vec = exp(-energy_vec/(k*T))
    prob_vec = prob_rel_vec/sum(prob_rel_vec)
    output = ReweightResults(prob_vec, prob_rel_vec, energy_vec)
    return output
end

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

    Average probability = $(Statistics.mean(res.probability))
    standard deviation = $(Statistics.std(res.probability))

    -------------------------------------------
    FRAME WEIGHTS RELATIVE TO THE ORIGINAL ONES
    -------------------------------------------

    Average probability = $(Statistics.mean(res.relative_probability))
    standard deviation = $(Statistics.std(res.relative_probability))

    ----------------------------------
    COMPUTED ENERGY AFTER PERTURBATION
    ----------------------------------

    Average energy = $(Statistics.mean(res.energy))
    standard deviation = $(Statistics.std(res.energy))
    
    """)
end

#Checking if PDB file and input match
function check_n_mol(simulation::Simulation, atom_vec_group::AbstractVector{<:Integer}, n_mol_per_group::Int, name::String)
    check = length(unique(PDBTools.residue.(simulation.atoms[atom_vec_group]))) #check number of residues (number of molecules) in PDB
    division = check ÷ n_mol_per_group
    quotient = mod(check, n_mol_per_group)
    if division == 1 && quotient == 0
        return "Number of molecules in the system seems to be correct for $name"
    elseif division != 1 && quotient == 0
        return @warn("""
            Number of molecules in the system ($check) seems to be a multiple of your input for $name.
            """)
    elseif quotient != 0
        return @warn("""
            The number of residues ($check) in the PDB file for $name is different than the number of molecules based on the number on the input   
            """)
    end
end

#Retrieving atoms from the contributions
function retrieve_contrib_indexes(simulation::Simulation, contrib_gp::Union{String, Function, Nothing}, gp::AbstractVector{<:Integer})
    if isnothing(contrib_gp)
        return gp
    else
        return findall(i -> i in PDBTools.select(simulation.atoms, contrib_gp), simulation.atoms[gp])
    end
end