using Base.Threads: @threads
using ChunkSplitters: index_chunks
using CellListMap: CellListMap, ParticleSystem, pairwise!, update!
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
    f = first_frame!(sim)
    uc_first_frame = unitcell(f)
    p_first_frame = positions(f)

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

    @threads for frame_inds in index_chunks(1:length(sim); n=parallel ? Threads.nthreads() : 1)
        local uc
        local index_current_frame
        local current_positions
        current_positions = copy(p_first_frame)
        systems = setup_particle_systems(
            selection_pairs, selection_data,
            uc_first_frame, donnor_acceptor_distance, parallel
        )
        for _ in frame_inds
            lock(sim) do
                f = current_frame(sim)
                current_positions .= positions(f)
                uc = unitcell(f)
                iframe += 1
                index_current_frame = iframe
                next!(prg)
                if frame_index(sim) < last(frame_range(sim))
                    next_frame!(sim)
                end
            end
            for selection_pair in selection_pairs
                sel1, sel2 = first(selection_pair), last(selection_pair)
                key = _key_name(sel1, sel2)
                sys = systems[key]
                if sel1 == sel2
                    s1 = selection_data[sel1]
                    sys.xpositions .= @view(current_positions[s1.inds])
                    update!(sys; unitcell=uc.matrix)
                    number_of_hbonds = count_hbonds(
                        sys, s1, angle_cutoff, electronegative_elements
                    )
                else
                    s1 = selection_data[sel1]
                    s2 = selection_data[sel2]
                    sys.xpositions .= @view(current_positions[s1.inds])
                    sys.ypositions .= @view(current_positions[s2.inds])
                    update!(sys; unitcell=uc.matrix)
                    number_of_hbonds = count_hbonds(
                        sys, s1, s2, sel1, sel2, angle_cutoff, electronegative_elements
                    )
                end
                hbonds[key][index_current_frame] = number_of_hbonds
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

"""
    HydrogenBondOccupancy

Structure that wraps the result of the `hydrogen_bond_occupancy` function.

# Fields

- `list::Vector{Vector{HBond}}`: for each frame of the simulation, the list
  of [`HBond`](@ref)s found in that frame.

!!! note
    A hydrogen bond is considered to be the same along the trajectory only if
    it is formed by the same donnor, polar hydrogen, and acceptor atoms. If
    the donnor and acceptor remain the same but the bridging polar hydrogen
    changes, that is counted as a different hydrogen bond.

"""
struct HydrogenBondOccupancy
    list::Vector{Vector{HBond}}
end

#=
    unique_hbonds(hbo::HydrogenBondOccupancy) -> Set{HBond}

The set of all distinct hydrogen bonds found over the whole trajectory.
=#
function unique_hbonds(hbo::HydrogenBondOccupancy)
    bonds = Set{HBond}()
    for l in hbo.list
        union!(bonds, l)
    end
    return bonds
end

function Base.show(io::IO, ::MIME"text/plain", hbo::HydrogenBondOccupancy)
    print(io, chomp(
        """
        -------------------------------------------------------------------
        HydrogenBondOccupancy data:
        -------------------------------------------------------------------
        Total number of distinct hydrogen bonds: $(length(unique_hbonds(hbo)))
        Mean number of hydrogen bonds per frame: $(mean(hbo))
        Minimum number of hydrogen bonds per frame: $(minimum(length(l) for l in hbo.list))
        Maximum number of hydrogen bonds per frame: $(maximum(length(l) for l in hbo.list))
        -------------------------------------------------------------------
        """
    ))
end

"""
    mean(hbo::HydrogenBondOccupancy)

Returns the average number of hydrogen bonds found per frame.

"""
mean(hbo::HydrogenBondOccupancy) = mean(length.(hbo.list))

"""
    hydrogen_bond_occupancy(sim::Simulation, sel1, sel1 => sel2,...; kargs...)

Function to compute, for each frame of the simulation, which hydrogen bonds
are present, keeping track of the identity of the donnor and acceptor atoms
involved. This is the same computation performed by [`hydrogen_bonds`](@ref),
except that instead of just counting the hydrogen bonds of each frame, their
identities are retained, so that, for instance, the persistence of a given
hydrogen bond along the trajectory can be studied with
[`intermittent_correlation`](@ref).

A hydrogen bond is considered to be the same along the trajectory only if it
is formed by the same donnor, polar hydrogen, and acceptor atoms.

### Arguments

Identical to those of [`hydrogen_bonds`](@ref).

### Returns

- An ordered dictionary in which the key specifies the selections and the
  values are [`HydrogenBondOccupancy`](@ref) objects.

### Optional keyword arguments

Identical to those of [`hydrogen_bonds`](@ref).

# Example

```jldoctest
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> sim = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> hbo = hydrogen_bond_occupancy(sim, "protein", show_progress=false);

julia> length.(hbo["protein => protein"].list)
5-element Vector{Int64}:
 32
 28
 27
 27
 26

```

"""
function hydrogen_bond_occupancy(
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
    f = first_frame!(sim)
    uc_first_frame = unitcell(f)
    p_first_frame = positions(f)

    # Selection data (polar bonds, atom indices, etc.) is reused from the
    # `hydrogen_bonds` initialization; the counts dictionary it also builds
    # is discarded here.
    _, selection_data = initialize_hbonds_data(
        sim, selection_pairs;
        parallel=parallel,
        donnor_acceptor_distance=donnor_acceptor_distance,
        electronegative_elements=electronegative_elements,
        d_covalent_bond=d_covalent_bond
    )

    raw_lists = OrderedDict{String,Vector{HBondList}}()
    for selection_pair in selection_pairs
        sel1, sel2 = first(selection_pair), last(selection_pair)
        raw_lists[_key_name(sel1, sel2)] = [HBondList(HBond[]) for _ in 1:length(sim)]
    end

    # Process frames in parallel
    iframe = 0
    restart!(sim)
    prg = Progress(length(sim); enabled=show_progress)

    @threads for frame_inds in index_chunks(1:length(sim); n=parallel ? Threads.nthreads() : 1)
        local uc
        local index_current_frame
        local current_positions
        current_positions = copy(p_first_frame)
        systems = setup_hbond_list_particle_systems(
            selection_pairs, selection_data,
            uc_first_frame, donnor_acceptor_distance, parallel
        )
        for _ in frame_inds
            lock(sim) do
                f = current_frame(sim)
                current_positions .= positions(f)
                uc = unitcell(f)
                iframe += 1
                index_current_frame = iframe
                next!(prg)
                if frame_index(sim) < last(frame_range(sim))
                    next_frame!(sim)
                end
            end
            for selection_pair in selection_pairs
                sel1, sel2 = first(selection_pair), last(selection_pair)
                key = _key_name(sel1, sel2)
                sys = systems[key]
                if sel1 == sel2
                    s1 = selection_data[sel1]
                    sys.xpositions .= @view(current_positions[s1.inds])
                    update!(sys; unitcell=uc.matrix)
                    hblist = list_hbonds(sys, s1, angle_cutoff, electronegative_elements)
                else
                    s1 = selection_data[sel1]
                    s2 = selection_data[sel2]
                    sys.xpositions .= @view(current_positions[s1.inds])
                    sys.ypositions .= @view(current_positions[s2.inds])
                    update!(sys; unitcell=uc.matrix)
                    hblist = list_hbonds(sys, s1, s2, sel1, sel2, angle_cutoff, electronegative_elements)
                end
                raw_lists[key][index_current_frame] = HBondList(copy(hblist.bonds))
            end
        end
    end

    result = OrderedDict{String,HydrogenBondOccupancy}()
    for (key, hblists) in raw_lists
        result[key] = _build_hydrogen_bond_occupancy(hblists)
    end
    return result
end

#=
    _build_hydrogen_bond_occupancy(hblists::Vector{HBondList}) -> HydrogenBondOccupancy

Builds the `HydrogenBondOccupancy` object from the per-frame `HBondList`s,
deduplicating any bond that might have been recorded more than once within
the same frame.
=#
function _build_hydrogen_bond_occupancy(hblists::Vector{HBondList})
    return HydrogenBondOccupancy([unique(hblist.bonds) for hblist in hblists])
end

@testitem "hydrogen_bond_occupancy" begin
    using MolSimToolkit
    using MolSimToolkit.Testing

    sim = Simulation(Testing.gmx_pdb, Testing.gmx_traj)

    hbo = hydrogen_bond_occupancy(sim, "protein")
    hbs = hydrogen_bonds(sim, "protein")

    # The per-frame counts must match those of `hydrogen_bonds`
    @test length.(hbo["protein => protein"].list) == hbs["protein => protein"]

    occ = hbo["protein => protein"]
    @test occ isa HydrogenBondOccupancy
    @test length(occ.list) == length(sim)
    @test all(bond -> bond isa HBond, reduce(vcat, occ.list))
    @test mean(occ) ≈ mean(hbs["protein => protein"])

    # Bonds present in every frame must compare equal (same donnor/hydrogen/acceptor)
    common_bonds = intersect(occ.list...)
    @test !isempty(common_bonds)
    for bond in common_bonds
        @test all(list -> bond in list, occ.list)
    end

    # A different bridging hydrogen must be treated as a different bond
    b1 = HBond(1, 2, 3)
    b2 = HBond(1, 4, 3)
    @test b1 != b2

    hbo2 = hydrogen_bond_occupancy(sim, "protein"; parallel=false)
    @test length.(hbo2["protein => protein"].list) == hbs["protein => protein"]

    hbo3 = hydrogen_bond_occupancy(sim, "protein" => "resname SOL")
    hbs3 = hydrogen_bonds(sim, "protein" => "resname SOL")
    @test length.(hbo3["protein => resname SOL"].list) == hbs3["protein => resname SOL"]
end