import PDBTools.transfer_free_energy
export transfer_free_energy
export transfer_free_energy_frames

#
# Internal: compute the transfer free energy of a single frame.
#
function _tfe_frame(
    frame, protein_at_frame, cosolvent, model, 
    inds_protein, reconstruct_protein, backbone, sidechain
)
    # fetch protein coordinates and unitcell
    p = positions(frame)
    uc = unitcell(frame)
    # Reconstruct the protein structure, which might be broken by the PBCs
    if reconstruct_protein
        _reconstruct_structure!(p, inds_protein, uc)
    end
    # update coordinates (note the dot for broadcast)
    PDBTools.set_position!.(protein_at_frame, @view(p[inds_protein]))
    # Compute the transfer free energy
    return PDBTools.transfer_free_energy(protein_at_frame, cosolvent; model, backbone, sidechain)
end 

"""
    transfer_free_energy(
        sim::Simulation,
        cosolvent="urea";
        model=PDBTools.AutonBolen,
        protein::AbstractVector{<:PDBTools.Atom}=PDBTools.select(get_atoms(sim), "protein and not element H"),
        backbone::Function=PDBTools.isbackbone,
        sidechain::Function=PDBTools.issidechain,
        reconstruct_protein::Bool=true,
    )

Calculates the transfer free energy (in kcal/mol at 1M cosolvent) averaged over all
frames of `sim`, using the Tanford transfer model.

If `reconstruct==true`, at each frame the protein structure is reconstructed to account for periodic boundary
conditions before computing the SASA-based transfer free energy via
`PDBTools.transfer_free_energy`.

## Positional arguments

- `sim`: `Simulation` object.
- `cosolvent`: cosolvent name (case insensitive). One of: `"betaine"`, `"proline"`,
  `"sarcosine"`, `"sorbitol"`, `"sucrose"`, `"tmao"`, `"urea"`, `"urea-app"`, `"urea-mh"`.
  Default: `"urea"`.

## Keyword arguments

- `model`: thermodynamic model. Either `PDBTools.AutonBolen` (default) or
  `PDBTools.MoeserHorinek`.
- `protein`: vector of protein atoms. Defaults to all non-hydrogen protein atoms in `sim` selected by `"protein and not element H"`.
- `backbone`: function that identifies backbone atoms. Default: `PDBTools.isbackbone`.
- `sidechain`: function that identifies side-chain atoms. Default: `PDBTools.issidechain`.
- `reconstruct_protein`: boolean to define if the protein structure must be reconstructed to avoid PBC artifacts. Defaults to `true`.

## Returns

A `PDBTools.TransferFreeEnergy` object with the frame-averaged values:

- `nres::Int`: number of residues considered.
- `tot::Float32`: total transfer free energy (kcal/mol/M).
- `bb::Float32`: backbone contribution (kcal/mol/M).
- `sc::Float32`: side-chain contribution (kcal/mol/M).
- `residue_contributions_bb::Vector{Float32}`: per-residue backbone contributions.
- `residue_contributions_sc::Vector{Float32}`: per-residue side-chain contributions.
- `cosolvent::String`: the cosolvent used.

## Example

```julia-repl
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> sim = Simulation(Testing.namd_pdb, Testing.namd_traj)

julia> tfe = transfer_free_energy(sim, "urea")

julia> tfe.tot   # total transfer free energy averaged over frames
```

## References

- Auton & Bolen: https://doi.org/10.1016/s0076-6879(07)28023-1 and https://www.pnas.org/doi/10.1073/pnas.0706251104
- Moeser & Horinek: https://doi.org/10.1021/jp409934q
"""
function PDBTools.transfer_free_energy(
    sim::Simulation,
    cosolvent="urea";
    model=PDBTools.AutonBolen,
    protein::AbstractVector{<:PDBTools.Atom}=PDBTools.select(get_atoms(sim), "protein and not element H"),
    backbone::F=PDBTools.isbackbone,
    sidechain::G=PDBTools.issidechain,
    reconstruct_protein::Bool=true,
) where {F<:Function,G<:Function}
    # indices of the protein atoms, to fetch frame coordinates
    inds_protein = PDBTools.index.(protein)
    # Create a copy to update coordinates at each frame
    protein_at_frame = copy.(protein)
    # Number of residues
    nres = length(PDBTools.eachresidue(protein))
    # Totals and contributions
    tot = 0.f0
    bb = 0.f0
    sc = 0.f0
    residue_contributions_bb = zeros(Float32, nres)
    residue_contributions_sc = zeros(Float32, nres)
    for frame in sim
        tfe = _tfe_frame(
            frame, protein_at_frame, cosolvent, model, inds_protein, 
            reconstruct_protein, backbone, sidechain
        )
        # push total, backbone and sidechain mvalues to arrays
        tot += tfe.tot
        bb += tfe.bb
        sc += tfe.sc
        residue_contributions_bb .+= tfe.residue_contributions_bb
        residue_contributions_sc .+= tfe.residue_contributions_sc
    end
    tot /= length(sim)
    bb /= length(sim)
    sc /= length(sim)
    residue_contributions_bb ./= length(sim)
    residue_contributions_sc ./= length(sim)
    return PDBTools.TransferFreeEnergy{model}(
        nres,
        tot,
        bb,
        sc,
        residue_contributions_bb,
        residue_contributions_sc,
        cosolvent
    )
end

struct TransferFreeEnergyFrames{T}
    _frames::Vector{PDBTools.TransferFreeEnergy{T}}
end
Base.getindex(tfe::TransferFreeEnergyFrames, i) = tfe._frames[i]
Base.firstindex(tfe::TransferFreeEnergyFrames) = firstindex(tfe._frames)
Base.lastindex(tfe::TransferFreeEnergyFrames) = lastindex(tfe._frames)
Base.first(tfe::TransferFreeEnergyFrames) = first(tfe._frames)
Base.last(tfe::TransferFreeEnergyFrames) = last(tfe._frames)
Base.length(tfe::TransferFreeEnergyFrames) = length(tfe._frames)
Base.keys(tfe::TransferFreeEnergyFrames) = LinearIndices(tfe._frames)
Base.iterate(tfe::TransferFreeEnergyFrames, i=firstindex(tfe)) = i > length(tfe) ? nothing : (tfe._frames[i], i + 1)
function Base.show(io::IO, tfe::TransferFreeEnergyFrames{T}) where {T}
    n = length(tfe)
    cosolvent = n > 0 ? first(tfe).cosolvent : "unknown"
    print(io, "TransferFreeEnergyFrames{$T} with $n frames (cosolvent: $cosolvent)")
    if n > 0
        _avg = SVector(0.f0, 0.f0, 0.f0)
        _min = SVector(Inf32, Inf32, Inf32)
        _max = SVector(-Inf32, -Inf32, -Inf32)
        for t in tfe
            _avg += SVector(t.tot, t.bb, t.sc)
            _min = SVector(min(_min[1], t.tot), min(_min[2], t.bb), min(_min[3], t.sc))
            _max = SVector(max(_max[1], t.tot), max(_max[2], t.bb), max(_max[3], t.sc))
        end
        _avg /= length(tfe)
        println(io)
        @printf(io, "  Total TFE  (kcal/mol/M): min = %8.3f  mean = %8.3f  max = %8.3f\n",
            _min[1], _avg[1], _max[1]) 
        @printf(io, "  Backbone   (kcal/mol/M): min = %8.3f  mean = %8.3f  max = %8.3f\n",
            _min[2], _avg[2], _max[2]) 
        @printf(io, "  Side-chain (kcal/mol/M): min = %8.3f  mean = %8.3f  max = %8.3f",
            _min[3], _avg[3], _max[3]) 
    end
end

"""
    transfer_free_energy_frames(
        sim::Simulation,
        cosolvent="urea";
        model=PDBTools.AutonBolen,
        protein::AbstractVector{<:PDBTools.Atom}=PDBTools.select(get_atoms(sim), "protein and not element H"),
        backbone::Function=PDBTools.isbackbone,
        sidechain::Function=PDBTools.issidechain,
        reconstruct_protein::Bool=true,
    )

Calculates the transfer free energy (in kcal/mol at 1M cosolvent) for each frame of
`sim` individually, using the Tanford transfer model. Unlike `transfer_free_energy`,
which returns frame-averaged values, this function returns a `TransferFreeEnergyFrames`
object containing one `PDBTools.TransferFreeEnergy` result per frame.

If `reconstruct_protein==true`, at each frame the protein structure is reconstructed to
account for periodic boundary conditions before computing the SASA-based transfer free
energy via `PDBTools.transfer_free_energy`.

## Positional arguments

- `sim`: `Simulation` object.
- `cosolvent`: cosolvent name (case insensitive). One of: `"betaine"`, `"proline"`,
  `"sarcosine"`, `"sorbitol"`, `"sucrose"`, `"tmao"`, `"urea"`, `"urea-app"`, `"urea-mh"`.
  Default: `"urea"`.

## Keyword arguments

- `model`: thermodynamic model. Either `PDBTools.AutonBolen` (default) or
  `PDBTools.MoeserHorinek`.
- `protein`: vector of protein atoms. Defaults to all non-hydrogen protein atoms in `sim` selected by `"protein and not element H"`.
- `backbone`: function that identifies backbone atoms. Default: `PDBTools.isbackbone`.
- `sidechain`: function that identifies side-chain atoms. Default: `PDBTools.issidechain`.
- `reconstruct_protein`: boolean to define if the protein structure must be reconstructed to avoid PBC artifacts. Defaults to `true`.

## Returns

A `TransferFreeEnergyFrames` object, which is an indexable and iterable collection of
`PDBTools.TransferFreeEnergy` objects, one per simulation frame. Each element has the
fields:

- `nres::Int`: number of residues considered.
- `tot::Float32`: total transfer free energy (kcal/mol/M).
- `bb::Float32`: backbone contribution (kcal/mol/M).
- `sc::Float32`: side-chain contribution (kcal/mol/M).
- `residue_contributions_bb::Vector{Float32}`: per-residue backbone contributions.
- `residue_contributions_sc::Vector{Float32}`: per-residue side-chain contributions.
- `cosolvent::String`: the cosolvent used.

## Example

```julia-repl
julia> using MolSimToolkit, MolSimToolkit.Testing

julia> sim = Simulation(Testing.namd_pdb, Testing.namd_traj)

julia> tfe_frames = transfer_free_energy_frames(sim, "urea")

julia> tfe_frames[1].tot   # total TFE of the first frame

julia> [ t.tot for t in tfe_frames ]   # total TFE for each frame
```

## References

- Auton & Bolen: https://doi.org/10.1016/s0076-6879(07)28023-1 and https://www.pnas.org/doi/10.1073/pnas.0706251104
- Moeser & Horinek: https://doi.org/10.1021/jp409934q
"""
function transfer_free_energy_frames(
    sim::Simulation,
    cosolvent="urea";
    model=PDBTools.AutonBolen,
    protein::AbstractVector{<:PDBTools.Atom}=PDBTools.select(get_atoms(sim), "protein and not element H"),
    backbone::F=PDBTools.isbackbone,
    sidechain::G=PDBTools.issidechain,
    reconstruct_protein::Bool=true,
) where {F<:Function,G<:Function}
    # indices of the protein atoms, to fetch frame coordinates
    inds_protein = PDBTools.index.(protein)
    # Create a copy to update coordinates at each frame
    protein_at_frame = copy.(protein)
    # Totals and contributions
    tfe_frames = TransferFreeEnergyFrames{model}(Vector{PDBTools.TransferFreeEnergy{model}}[])
    for frame in sim
        tfe = _tfe_frame(
            frame, protein_at_frame, cosolvent, model, inds_protein, 
            reconstruct_protein, backbone, sidechain
        )
        # push total, backbone and sidechain mvalues to arrays
        push!(tfe_frames._frames, tfe)
    end
    return tfe_frames
end

@testitem "transfer_free_energy" begin
    using MolSimToolkit
    using MolSimToolkit.Testing
    using PDBTools
    using ShowMethodTesting

    # First frame only
    sim = Simulation(Testing.namd_pdb, Testing.namd_traj; first=1, last=1)
    tfe = transfer_free_energy(sim, "urea")

    ats = get_atoms(sim)
    f = first_frame!(sim)
    set_position!.(ats, positions(f))
    protein = select(ats, "protein")
    tfe_ref = transfer_free_energy(protein, "urea")

    @test tfe.tot ≈ tfe_ref.tot
    @test tfe.bb ≈ tfe_ref.bb
    @test tfe.sc ≈ tfe_ref.sc
    @test all(tfe.residue_contributions_bb .≈ tfe_ref.residue_contributions_bb)
    @test all(tfe.residue_contributions_sc .≈ tfe_ref.residue_contributions_sc)

    # Average of first and second frame
    sim = Simulation(Testing.namd_pdb, Testing.namd_traj; first=1, last=2)
    f = get_frame!(sim, 2)
    set_position!.(ats, positions(f))
    protein = select(ats, "protein")
    tfe_ref2 = transfer_free_energy(protein, "urea")

    tfe = transfer_free_energy(sim, "urea")
    @test tfe.tot ≈ (tfe_ref.tot + tfe_ref2.tot) / 2
    @test tfe.bb ≈ (tfe_ref.bb + tfe_ref2.bb) / 2
    @test tfe.sc ≈ (tfe_ref.sc + tfe_ref2.sc) / 2
    @test all(tfe.residue_contributions_bb .≈ 0.5 * (tfe_ref.residue_contributions_bb + tfe_ref2.residue_contributions_bb))
    @test all(tfe.residue_contributions_sc .≈ 0.5 * (tfe_ref.residue_contributions_sc + tfe_ref2.residue_contributions_sc))

    # test for time-dependent TFE computation
    sim = Simulation(Testing.namd_pdb, Testing.namd_traj; first=1, last=2)
    tfe_frames = transfer_free_energy_frames(sim, "urea")
    @test length(tfe_frames) == 2
    @test firstindex(tfe_frames) == 1
    @test lastindex(tfe_frames) == 2
    @test first(tfe_frames) === tfe_frames[1]
    @test last(tfe_frames) === tfe_frames[2]
    @test tfe_frames[1].tot ≈ tfe_ref.tot
    @test tfe_frames[2].tot ≈ tfe_ref2.tot
    @test keys(tfe_frames) == [1,2]
    @test eachindex(tfe_frames) == [1,2]
    @test [ t for t in tfe_frames ] == [ first(tfe_frames), last(tfe_frames) ]
    @test parse_show(tfe_frames; repl=Dict("MolSimToolkit." => "", "PDBTools." => "")) ≈ """
        TransferFreeEnergyFrames{AutonBolen} with 2 frames (cosolvent: urea)
          Total TFE  (kcal/mol/M): min =   -0.543  mean =   -0.536  max =   -0.529
          Backbone   (kcal/mol/M): min =   -0.552  mean =   -0.536  max =   -0.521
          Side-chain (kcal/mol/M): min =   -0.008  mean =    0.001  max =    0.009
    """

end