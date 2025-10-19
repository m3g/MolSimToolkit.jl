using Base.Threads: @spawn
using ChunkSplitters: index_chunks

"""
    hydrogen_bonds(
        sim::Simulation, 
        sel1::Union{String,Function}=at -> true, 
        sel2::Union{Nothing,String,Function}=nothing; 
        nthreads=Threads.nthreads(),
        show_progress::Bool=true,
        parallel::Bool=true,
        kargs...
    )

Compute the number of hydrogen bonds for each frame of a simulation. 

For further details on `kargs...` see the help entry of `PDBTools.hydrogen_bonds`.

"""
function PDBTools.hydrogen_bonds(
    sim::Simulation,
    sel1::Union{String,Function}=at->true,
    sel2::Union{String,Function,Nothing}=nothing;
    parallel::Bool=true,
    show_progress::Bool=true,
    kargs...
)
    hbonds = zeros(Int,length(sim))
    iframe = 0
    lk = ReentrantLock()
    restart!(sim)
    prg = Progress(length(sim); enabled=show_progress)
    @sync for frame_inds in index_chunks(1:length(sim); n=parallel ? Threads.nthreads() : 1)
        @spawn begin
            local index_current_frame
            local uc
            ats = copy.(atoms(sim))
            for _ in frame_inds
                lock(lk) do
                    nextframe!(sim)
                    iframe += 1
                    next!(prg)
                    index_current_frame = iframe
                    p = positions(current_frame(sim))
                    uc = unitcell(current_frame(sim))
                    set_position!.(ats, p)
                end
                if isnothing(sel2)
                    hb = hydrogen_bonds(ats, sel1; unitcell=uc.matrix, kargs...)
                else
                    hb = hydrogen_bonds(ats, sel1, sel2; unitcell=uc.matrix, kargs...)
                end
                hbonds[index_current_frame] = length(hb)
            end
        end
    end
    return hbonds
end

@testitem "hydrogen bonds" begin
    using MolSimToolkit
    using PDBTools
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

end