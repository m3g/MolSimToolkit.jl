using Base.Threads: @spawn
using ChunkSplitters: index_chunks
using CellListMap: CellListMap, ParticleSystem, map_pairwise, update_unitcell!

include("./internals.jl")

"""
function PDBTools.hydrogen_bonds(
    sim::Simulation,
    sel1::Union{String,Function}=at -> true
    sel2::Union{Nothing,String,Function}=nothing;
    parallel::Bool=true,
    show_progress::Bool=true,
    donnor_acceptor_distance::Real=3.5f0,
    angle_cutoff::Real=30,
    electronegative_elements=("N", "O", "F", "S"),
    d_covalent_bond::Real=1.2f0,
)

!!! warning
    Experimental feature: interface changes can occur in non-breaking releases.

Function to compute the number of hydrogen bonds per frame in a simulation.

### Arguments

- `ats::AbstractVector{<:PDBTools.Atom}`: Vector of atoms to analyze.
and
- `sel1::Union{Function, String}=at -> true`: Selection of atoms to consider. Can be a function or a selection string. This selection is optional.
or
- `sel1::Union{Function, String}`: First selection of atoms to consider. Can be a function or a selection string.
- `sel2::Union{Function, String}`: Second selection of atoms to consider. Can be a function or a selection string.

In the case of two sets, `sel1` and `sel2` must not overlap.

### Keyword Arguments

- `parallel::Bool=true`: Defines if the calculation is run in parallel. Requires starting Julia with multi-threading.
- `unitcell::Union{Nothing,AbstractVecOrMat}=nothing`: Unit cell for periodic boundary conditions.
- `donnor_acceptor_distance::Real=3.5f0`: Maximum distance between donnor and acceptor to consider a hydrogen bond.
- `angle_cutoff::Real=30`: Maximum angle (in degrees) between donnor-hydrogen-acceptor to consider a hydrogen bond.
- `electronegative_elements=("N", "O", "F", "S")`: Elements considered electronegative for hydrogen bonding.
- `d_covalent_bond::Real=1.2f0`: Maximum distance between donnor and hydrogen to consider a covalent bond.

### Returns

- A vector with the number of hydrogen bonds per frame.

# Example

```jldoctest
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> sim = Simulation(Testing.namd_pdb, Testing.namd_traj);

julia> hbs = hydrogen_bonds(sim, "protein")
5-element Vector{Int64}:
 32
 28
 27
 27
 26

julia> hbs = hydrogen_bonds(sim, "protein", "water")
5-element Vector{Int64}:
 75
 81
 76
 68
 80
```

!!! note
    This function does not use topology information. It identified polar hydrogens based on distance criteria only,
    where `d_covalent_bond` is the criterium for identifying covalent bonds between donnor and hydrogen atoms.

"""
function PDBTools.hydrogen_bonds(
    sim::Simulation,
    sel1::Union{String,Function}=at -> true;
    parallel::Bool=true,
    show_progress::Bool=true,
    donnor_acceptor_distance::Real=3.5f0,
    angle_cutoff::Real=30,
    electronegative_elements=("N", "O", "F", "S"),
    d_covalent_bond::Real=1.2f0,
)
    ats_sel = select(atoms(sim), sel1)
    inds_sel = index.(ats_sel)
    firstframe!(sim)
    uc_first_frame = unitcell(current_frame(sim))
    p = @view(positions(current_frame(sim))[inds_sel])

    polar_bonds = find_hbond_donnors(
        ats_sel;
        positions=p,
        unitcell=uc_first_frame.matrix,
        cutoff=donnor_acceptor_distance,
        parallel,
        electronegative_elements,
        d_covalent_bond,
    )

    hbonds = zeros(Int, length(sim))
    iframe = 0
    lk = ReentrantLock()
    restart!(sim)
    prg = Progress(length(sim); enabled=show_progress)
    @sync for frame_inds in index_chunks(1:length(sim); n=parallel ? Threads.nthreads() : 1)
        @spawn begin
            local index_current_frame
            local uc
            ats_sel_positions = coor.(ats_sel)
            sys = ParticleSystem(
                positions=ats_sel_positions,
                unitcell=uc_first_frame.matrix,
                cutoff=donnor_acceptor_distance,
                output=0,
                parallel=parallel,
                output_name=:number_of_hbonds
            )
            for _ in frame_inds
                lock(lk) do
                    nextframe!(sim)
                    iframe += 1
                    next!(prg)
                    index_current_frame = iframe
                    p = @view(positions(current_frame(sim))[inds_sel])
                    uc = unitcell(current_frame(sim))
                    ats_sel_positions .= p
                end
                # Update ParticleSystem for this frame
                sys.xpositions .= ats_sel_positions
                update_unitcell!(sys, uc.matrix)
                # Compute number of hydrogen bonds
                number_of_hbonds = map_pairwise(
                    (x, y, i, j, d2, number_of_hbonds) -> begin
                        el_i = element(ats_sel[i])
                        el_j = element(ats_sel[j])
                        if (el_i in electronegative_elements) & (el_j in electronegative_elements)
                            number_of_hbonds += count_hbonds(
                                i, x, y, polar_bonds,
                                sys.positions, sys.unitcell,
                                angle_cutoff,
                            )
                            number_of_hbonds += count_hbonds(
                                j, y, x, polar_bonds,
                                sys.positions, sys.unitcell,
                                angle_cutoff,
                            )
                        end
                        return number_of_hbonds
                    end,
                    sys,
                )
                hbonds[index_current_frame] = number_of_hbonds
            end
        end
    end
    return hbonds
end

#
# Function that receives the two selections
#
function PDBTools.hydrogen_bonds(
    sim::Simulation,
    sel1::Union{String,Function},
    sel2::Union{String,Function};
    parallel::Bool=true,
    show_progress::Bool=true,
    donnor_acceptor_distance::Real=3.5f0,
    angle_cutoff::Real=30,
    electronegative_elements=("N", "O", "F", "S"),
    d_covalent_bond::Real=1.2f0,
)
    ats_sel1 = select(atoms(sim), sel1)
    inds_sel1 = index.(ats_sel1)
    ats_sel2 = select(atoms(sim), sel2)
    inds_sel2 = index.(ats_sel2)
    firstframe!(sim)
    uc_first_frame = unitcell(current_frame(sim))
    p1 = @view(positions(current_frame(sim))[inds_sel1])
    p2 = @view(positions(current_frame(sim))[inds_sel2])

    polar_bonds1 = find_hbond_donnors(
        ats_sel1;
        positions=p1,
        unitcell=uc_first_frame.matrix,
        cutoff=donnor_acceptor_distance,
        parallel,
        electronegative_elements,
        d_covalent_bond,
    )
    polar_bonds2 = find_hbond_donnors(
        ats_sel2;
        positions=p2,
        unitcell=uc_first_frame.matrix,
        cutoff=donnor_acceptor_distance,
        parallel,
        electronegative_elements,
        d_covalent_bond,
    )

    hbonds = zeros(Int, length(sim))
    iframe = 0
    lk = ReentrantLock()
    restart!(sim)
    prg = Progress(length(sim); enabled=show_progress)
    @sync for frame_inds in index_chunks(1:length(sim); n=parallel ? Threads.nthreads() : 1)
        @spawn begin
            local index_current_frame
            local uc
            ats_sel_positions1 = coor.(ats_sel1)
            ats_sel_positions2 = coor.(ats_sel2)
            sys = ParticleSystem(
                xpositions=ats_sel_positions1,
                ypositions=ats_sel_positions2,
                unitcell=uc_first_frame.matrix,
                cutoff=donnor_acceptor_distance,
                output=0,
                parallel=parallel,
                output_name=:number_of_hbonds,
            )
            for _ in frame_inds
                lock(lk) do
                    nextframe!(sim)
                    iframe += 1
                    next!(prg)
                    index_current_frame = iframe
                    p1 = @view(positions(current_frame(sim))[inds_sel1])
                    p2 = @view(positions(current_frame(sim))[inds_sel2])
                    uc = unitcell(current_frame(sim))
                    ats_sel_positions1 .= p1
                    ats_sel_positions2 .= p2
                end
                # Update ParticleSystem for this frame
                sys.xpositions .= ats_sel_positions1
                sys.ypositions .= ats_sel_positions2
                update_unitcell!(sys, uc.matrix)
                # Compute number of hydrogen bonds
                number_of_hbonds = map_pairwise(
                    (x, y, i, j, d2, number_of_hbonds) -> begin
                        at_i = ats_sel1[i]
                        at_j = ats_sel2[j]
                        if index(at_i) == index(at_j)
                            throw(ArgumentError("""\n
                                Selections cannot overlap. Detected atom $(index(at_i)) in both selections.

                                """))
                        end
                        el_i = element(at_i)
                        el_j = element(at_j)
                        if (el_i in electronegative_elements) & (el_j in electronegative_elements)
                            number_of_hbonds += count_hbonds2(
                                i, x, y, polar_bonds1, sys.xpositions,
                                sys.unitcell, angle_cutoff
                            )
                            number_of_hbonds += count_hbonds2(
                                j, y, x, polar_bonds2, sys.ypositions,
                                sys.unitcell, angle_cutoff
                            )
                        end
                        return number_of_hbonds
                    end,
                    sys,
                )
                hbonds[index_current_frame] = number_of_hbonds
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
    @test hbs == [58, 60, 54, 54, 58]

    hbs = hydrogen_bonds(sim, "protein"; parallel=false)
    @test hbs == [58, 60, 54, 54, 58]

    hbs = hydrogen_bonds(sim, "protein", "resname HOH SOL")
    @test hbs == [152, 153, 149, 149, 157]

    hbs = hydrogen_bonds(sim, "resname HOH", "resname SOL")
    @test hbs == [151, 155, 152, 147, 147]

    @test_throws ArgumentError hydrogen_bonds(sim, "protein", "protein")

end