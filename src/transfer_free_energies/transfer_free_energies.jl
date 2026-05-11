import PDBTools.transfer_free_energy
export transfer_free_energy

"""
    transfer_free_energy(
        sim::Simulation,
        cosolvent="urea";
        model=PDBTools.AutonBolen,
        protein::AbstractVector{<:PDBTools.Atom}=PDBTools.select(get_atoms(sim), "protein"),
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
- `protein`: vector of protein atoms. Defaults to all atoms in `sim` selected by `"protein"`.
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
    protein::AbstractVector{<:PDBTools.Atom}=PDBTools.select(get_atoms(sim), "protein"),
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
        # fetch protein coordinates and unitcell
        p = @view(positions(frame)[inds_protein])
        uc = unitcell(frame)
        # Reconstruct the protein structure, which might be broken by the PBCs
        if reconstruct_protein
            _reconstruct_structure!(p, inds_protein, uc)
        end
        # update coordinates (note the dot for broadcast)
        PDBTools.set_position!.(protein_at_frame, p)
        # Compute the transfer free energy
        tfe = PDBTools.transfer_free_energy(protein_at_frame, cosolvent; model, backbone, sidechain)
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

@testitem "transfer_free_energy" begin
    using MolSimToolkit
    using MolSimToolkit.Testing
    using PDBTools

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

end