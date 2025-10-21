using Base.Threads: @spawn
using ChunkSplitters: index_chunks
using CellListMap: CellListMap, ParticleSystem, map_pairwise!, update_unitcell!
using OrderedCollections: OrderedDict

include("./internals.jl")

"""
    hydrogen_bonds(sim::Simulation, sel1, sel1 => sel2,...; kargs...)

Function to compute the number of hydrogen bonds per frame in a simulation.

### Arguments

- `sim::Simulation`: The `Simulation` object. 
and, optionally,
- a list of selections, or pairs of selections, given as selection strings. Examples: `"protein"`, 
  `"protein" => "water"`, etc.

If no selection is provided, the hydrogen bonds among all atoms are computed.
If two selections of a pair are different, their atoms must not overlap (an error will be thrown).

### Returns

- An ordered dictionary in which the key specifies the selections and the values are vectors with the 
  number of hydrogen bonds in each frame.

### Optional keyword arguments

- `parallel::Bool=true`: Defines if the calculation is run in parallel. Requires starting Julia with multi-threading.
- `donnor_acceptor_distance::Real=3.5f0`: Maximum distance between donnor and acceptor to consider a hydrogen bond.
- `angle_cutoff::Real=30`: Maximum angle (in degrees) between donnor-hydrogen-acceptor to consider a hydrogen bond.
- `electronegative_elements=("N", "O", "F", "S")`: Elements considered electronegative for hydrogen bonding.
- `d_covalent_bond::Real=1.2f0`: Maximum distance between donnor and hydrogen to consider a covalent bond.

# Example

```jldoctest
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> sim = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> hbs = hydrogen_bonds(sim, "protein")
OrderedCollections.OrderedDict{String, Vector{Int64}} with 1 entry:
  "protein => protein" => [32, 28, 27, 27, 26]

julia> hbs = hydrogen_bonds(sim, "protein" => "water")
OrderedCollections.OrderedDict{String, Vector{Int64}} with 1 entry:
  "protein => water" => [75, 81, 76, 68, 80]

julia> hbs = hydrogen_bonds(sim, "protein", "protein" => "water", "water" => "resname POPC")
OrderedCollections.OrderedDict{String, Vector{Int64}} with 3 entries:
  "protein => protein"    => [32, 28, 27, 27, 26]
  "protein => water"      => [75, 81, 76, 68, 80]
  "water => resname POPC" => [413, 403, 406, 392, 376]
```

!!! note
    This function does not use topology information. It identified polar hydrogens based on distance criteria only,
    where `d_covalent_bond` is the criterium for identifying covalent bonds between donnor and hydrogen atoms.

"""
function PDBTools.hydrogen_bonds(
    sim::Simulation,
    selections::Union{Nothing,String,Pair{String,String}}...=nothing;
    parallel::Bool=true,
    show_progress::Bool=true,
    donnor_acceptor_distance::Real=3.5f0,
    angle_cutoff::Real=30,
    electronegative_elements=("N", "O", "F", "S"),
    d_covalent_bond::Real=1.2f0,
)
    selection_pairs = process_selections(selections)

    # Initialize trajectory and get first frame data
    first_frame!(sim)
    uc_first_frame = unitcell(current_frame(sim))
    p_first_frame = positions(current_frame(sim))

    # Initialize results and process selections
    hbonds, selection_data = initialize_hbonds_data(
        sim, selection_pairs;
        parallel=parallel,
        donnor_acceptor_distance=donnor_acceptor_distance,
        electronegative_elements=electronegative_elements,
        d_covalent_bond=d_covalent_bond
    )

    # Process frames in parallel
    iframe = 0
    restart!(sim)
    prg = Progress(length(sim); enabled=show_progress)

    @sync for frame_inds in index_chunks(1:length(sim); n=parallel ? Threads.nthreads() : 1)
        @spawn begin
            local uc
            local index_current_frame
            local current_positions = copy(p_first_frame)
            systems = setup_particle_systems(
                selection_pairs, selection_data,
                uc_first_frame, donnor_acceptor_distance, parallel
            )
            for _ in frame_inds
                lock(sim) do
                    next_frame!(sim)
                    current_positions .= positions(current_frame(sim))
                    uc = unitcell(current_frame(sim))
                    iframe += 1
                    index_current_frame = iframe
                    next!(prg)
                end
                for selection_pair in selection_pairs
                    sel1 = first(selection_pair)
                    sel2 = last(selection_pair)
                    if sel1 == sel2
                        sys = systems[_key_name(sel1, sel2)]
                        s1 = selection_data[sel1]
                        # Update ParticleSystem for this frame
                        sys.xpositions .= @view(current_positions[s1.inds])
                        update_unitcell!(sys, uc.matrix)
                        # Compute number of hydrogen bonds
                        number_of_hbonds = map_pairwise!(sys) do x, y, i, j, _, number_of_hbonds
                            el_i = element(s1.ats[i])
                            el_j = element(s1.ats[j])
                            if (el_i in electronegative_elements) & (el_j in electronegative_elements)
                                number_of_hbonds += count_hbonds(
                                    i, x, y, s1.polar_bonds,
                                    sys.positions, sys.unitcell,
                                    angle_cutoff,
                                )
                                number_of_hbonds += count_hbonds(
                                    j, y, x, s1.polar_bonds,
                                    sys.positions, sys.unitcell,
                                    angle_cutoff,
                                )
                            end
                            return number_of_hbonds
                        end
                    elseif sel1 != sel2
                        sys = systems[_key_name(sel1, sel2)]
                        s1 = selection_data[sel1]
                        s2 = selection_data[sel2]
                        # Update ParticleSystem for this frame
                        sys.xpositions .= @view(current_positions[s1.inds])
                        sys.ypositions .= @view(current_positions[s2.inds])
                        update_unitcell!(sys, uc.matrix)
                        number_of_hbonds = map_pairwise!(sys) do x, y, i, j, _, number_of_hbonds
                            at_i = s1.ats[i]
                            at_j = s2.ats[j]
                            if index(at_i) == index(at_j)
                                throw(ArgumentError("""\n
                                    Different selections cannot overlap. Detected atom $(index(at_i)) in both selections \"$sel1\" and \"$sel2\".

                                    """))
                            end
                            el_i = element(at_i)
                            el_j = element(at_j)
                            if (el_i in electronegative_elements) & (el_j in electronegative_elements)
                                number_of_hbonds += count_hbonds2(
                                    i, x, y, s1.polar_bonds, sys.xpositions,
                                    sys.unitcell, angle_cutoff
                                )
                                number_of_hbonds += count_hbonds2(
                                    j, y, x, s2.polar_bonds, sys.ypositions,
                                    sys.unitcell, angle_cutoff
                                )
                            end
                            return number_of_hbonds
                        end
                    end
                    hbonds[_key_name(sel1, sel2)][index_current_frame] = number_of_hbonds
                end
            end
        end
    end

    return hbonds
end

@testitem "hydrogen bonds" begin
    using MolSimToolkit
    using MolSimToolkit.Testing

    # Tested vs. gmx hbond
    sim = Simulation(Testing.gmx_pdb, Testing.gmx_traj)

    hbs = hydrogen_bonds(sim, "protein")
    @test hbs["protein => protein"] == [58, 60, 54, 54, 58]

    hbs = hydrogen_bonds(sim, "protein"; parallel=false)
    @test hbs["protein => protein"] == [58, 60, 54, 54, 58]

    hbs = hydrogen_bonds(sim, "protein" => "resname SOL")
    @test hbs["protein => resname SOL"] == [152, 153, 149, 149, 157]

    hbs = hydrogen_bonds(sim, "resname SOL and residue < 7000" => "resname SOL and residue >= 7000")
    @test hbs["resname SOL and residue < 7000 => resname SOL and residue >= 7000"] == [9049, 9062, 8903, 8977, 8857]

    hbs = hydrogen_bonds(sim; electronegative_elements=("O", "N"))
    @test hbs["all => all"] == [18231, 18205, 18113, 18063, 18090]

    hbs = hydrogen_bonds(sim,
        "protein" => "protein",
        "protein" => "resname SOL",
        "resname SOL and residue < 7000" => "resname SOL and residue >= 7000",
    )
    @test hbs["protein => protein"] == [58, 60, 54, 54, 58]
    @test hbs["protein => resname SOL"] == [152, 153, 149, 149, 157]
    @test hbs["resname SOL and residue < 7000 => resname SOL and residue >= 7000"] == [9049, 9062, 8903, 8977, 8857]

    @test_throws "overlap" hydrogen_bonds(sim, "protein" => "protein and resname ARG")
    @test_throws "overlap" hydrogen_bonds(sim, "resname SOL" => "resname SOL and residue 7000")

end