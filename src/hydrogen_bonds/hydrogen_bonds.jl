using Base.Threads: @spawn
using ChunkSplitters: index_chunks

"""
    hydrogen_bonds(
        sim::Simulation, 
        sel::Union{String,Function}; 
        nthreads=Threads.nthreads(),
        show_progress::Bool=true,
        kargs...
    )
    hydrogen_bonds(
        sim::Simulation, 
        sel1::Union{String,Function}, 
        sel2::Union{String,Function}; 
        nthreads=Threads.nthreads(),
        kargs...
    )

Compute the number of hydrogen bonds for each frame of a simulation. 

For further details on `kargs...` see the help entry of `PDBTools.hydrogen_bonds`.

"""
function PDBTools.hydrogen_bonds(
    sim::Simulation,
    sel1::Union{String,Function}=at->true,
    sel2::Union{String,Function,Nothing}=nothing;
    nthreads::Integer=Threads.nthreads(),
    show_progress::Bool=true,
    kargs...
)
    hbonds = zeros(Int,length(sim))
    iframe = 0
    lk = ReentrantLock()
    restart!(sim)
    prg = Progress(length(sim); enabled=show_progress)
    @sync for frame_inds in index_chunks(1:length(sim); n=nthreads)
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

    sim = Simulation(Testing.namd_pdb, Testing.namd_traj)

    hbs = hydrogen_bonds(sim, "protein")
    @test hbs == [32, 28, 27, 27, 26]

    hbs = hydrogen_bonds(sim, "protein", "water")
    @test hbs == [75, 81, 76, 68, 80]

    hbs = hydrogen_bonds(sim, "resname POPC", "water")
    @test hbs == [413, 403, 406, 392, 376]

end