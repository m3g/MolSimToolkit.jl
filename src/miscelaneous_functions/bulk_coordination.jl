#
# Function that given a vector with the indices of some continuous group of atoms
# that belong to the same type of molecule, and the number of atoms per molecule,
# and the molecule index, returns the range of atoms of the molecule
# 
# Here:
# - `imol` is the index of the molecule 
# - `n` is the number of atoms per molecule
# - `indices` is the vector with the indices of the atoms of the group
#
_imol_atoms(imol, n, indices) = first(indices) .+ (((imol-1)*n):(imol*n-1))

"""
    bulk_coordination(
        simulation::Simulation;
        reference_solute::AbstractVector{<:PDBTools.Atom}, 
        solute::AbstractVector{<:PDBTools.Atom},
        n_atoms_per_molecule_solute::Integer,
        solvent::AbstractVector{<:PDBTools.Atom}, 
        n_atoms_per_molecule_solvent::Integer,
        dmax::Real = 5.0,
        cutoff::Real = 20.0,
        bin_size::Real = 0.1,
        show_progress = true,
    )

Computes the coordination number of one type of solvent molecule relative to another 
solvent molecule, as a function of the distance to a reference solute molecule. 

For example, imagine a protein solvated in a mixture of water and TMAO. This function
allows to compute the number of water molecules that are within a given distance to the
TMAO molecules, as a function of the distance to the protein. That is, it computes the
coordination number of water relative to TMAO, as a function of the distance to the protein.

The function returns the the distances and the histogram of the coordination number as a function of 
the distance.

!!! compat
    This function was introduced in version 1.11.0 of MolSimToolkit.jl.

# Arguments

- `simulation::Simulation`: The simulation object
- `reference_solute::AbstractVector{PDBTools.Atom}`: The atoms of the reference solute molecule
   (the protein in the example above).
- `solute::AbstractVector{PDBTools.Atom}`: The atoms of the solute molecule (the TMAO in the example above).
- `n_atoms_per_molecule_solute::Integer`: The number of atoms per molecule of the solute.
- `solvent::AbstractVector{PDBTools.Atom}`: The atoms of the solvent molecule (the water in the example above).
- `n_atoms_per_molecule_solvent::Integer`: The number of atoms per molecule of the solvent.
- `dmax::Real = 5.0`: The maximum distance to the solute to consider a solvent molecule as coordinated.
- `cutoff::Real = 20.0`: The maximum distance to the reference molecule for computing the histogram.
- `bin_size::Real = 0.1`: The size of the bins for the histogram.
- `show_progress = true`: Whether to show a progress bar.

# Example

```julia-repl
julia> using MolSimToolkit, PDBTools

julia> pdb = readPDB(MolSimToolkit.Testing.namd2_pdb); # protein-tmao-water system

julia> trajectory = MolSimToolkit.Testing.namd2_traj;

julia> simulation = Simulation(pdb, trajectory)
Simulation 
    Atom type: Atom
    PDB file: -
    Trajectory file: /home/leandro/.julia/dev/MolSimToolkit/test/data/namd/protein_in_water_tmao/trajectory.dcd
    Total number of frames: 20
    Frame range: 1:1:20
    Number of frames in range: 20
    Current frame: nothing

julia> r, h = bulk_coordination(
           simulation;
           reference_solute = select(pdb, "protein"),
           solute = select(pdb, "resname TMAO"),
           n_atoms_per_molecule_solute = 14,
           solvent = select(pdb, "water"),
           n_atoms_per_molecule_solvent = 3,
       );

julia> using Plots

julia> plot(r, h, # plots `h` as a function of `r`
           xlabel = "Distance to protein (Å)",
           ylabel = "TMAO-water Coordination number",
           linewidth=2,
           label=:none, framestyle=:box, fontfamily="Computer Modern",
       )
```

"""
function bulk_coordination(
    simulation::Simulation;
    reference_solute::AbstractVector{<:PDBTools.Atom},
    solute::AbstractVector{<:PDBTools.Atom},
    n_atoms_per_molecule_solute::Integer,
    solvent::AbstractVector{<:PDBTools.Atom},
    n_atoms_per_molecule_solvent::Integer,
    dmax::Real=5.0,
    cutoff::Real=20.0,
    bin_size::Real=0.1,
    show_progress=true,
)

    # Indices of the atoms of each group
    reference_solute_indices = PDBTools.index.(reference_solute)
    solute_indices = PDBTools.index.(solute)
    solvent_indices = PDBTools.index.(solvent)

    # Initialize the systems for the minimum distances, to avoid creating
    # the systems at every frame
    first_frame!(simulation)
    uc = unitcell(current_frame(simulation))
    p = positions(current_frame(simulation))
    sys1 = CrossPairs(
        xpositions=p[solute_indices],
        ypositions=p[reference_solute_indices],
        xn_atoms_per_molecule=n_atoms_per_molecule_solute,
        cutoff=cutoff,
        unitcell=uc,
    )
    imol_atoms = _imol_atoms(1, n_atoms_per_molecule_solute, solute_indices)
    sys2 = CrossPairs(
        xpositions=p[solvent_indices],
        ypositions=p[imol_atoms],
        xn_atoms_per_molecule=n_atoms_per_molecule_solvent,
        cutoff=dmax,
        unitcell=uc,
    )

    # Initialize the histogram
    histogram = zeros(ceil(Int, cutoff / bin_size))
    nmols_in_bin = copy(histogram)

    progress = Progress(length(simulation); enabled=show_progress)
    for frame in simulation
        next!(progress)
        p = positions(frame)
        uc = unitcell(frame)
        # Compute all the minimum-distances between the reference molecule and the solute
        sys1.xpositions .= @view(p[solute_indices])
        sys1.ypositions .= @view(p[reference_solute_indices])
        sys1.unitcell = uc
        reference_list = minimum_distances!(sys1)
        # Now, we will run over the reference_list, and for each distances
        # to the reference, find the number of solvent atoms that are closer 
        # than dmax to the given solute, to add the number to the histogram
        # of the coordination number
        for (imol, md) in enumerate(reference_list)
            !(md.within_cutoff) && continue
            imol_atoms = _imol_atoms(imol, n_atoms_per_molecule_solute, solute_indices)
            sys2.xpositions .= @view(p[solvent_indices])
            sys2.ypositions .= @view(p[imol_atoms])
            sys2.unitcell = uc
            solvent_distances = minimum_distances!(sys2)
            bin = max(1, ceil(Int, md.d / bin_size))
            nmols_in_bin[bin] += 1
            histogram[bin] += count(x -> x.within_cutoff, solvent_distances)
        end
    end

    histogram ./= max.(1, nmols_in_bin)
    r = [(i-1)*bin_size + 0.5*bin_size for i in eachindex(histogram) ]
    return r, histogram
end

@testitem "bulk_coordination" begin
    using MolSimToolkit: Testing, bulk_coordination, Simulation
    using PDBTools: select, readPDB
    pdb = readPDB(Testing.namd2_pdb) # protein-tmao-water system
    trajectory = Testing.namd2_traj
    simulation = Simulation(pdb, trajectory; last=1)
    r, h = bulk_coordination(
        simulation;
        reference_solute=select(pdb, "protein"),
        solute=select(pdb, "resname TMAO"),
        n_atoms_per_molecule_solute=14,
        solvent=select(pdb, "water"),
        n_atoms_per_molecule_solvent=3,
    )
    @test r ≈ [0.05 + 0.1*(i-1) for i in eachindex(h)]
    @test length(r) == length(h)
    @test sum(h) ≈ 2678.5
end
