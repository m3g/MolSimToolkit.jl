"""
    centerofmass(selection::Vector{Atom})

Calculate the center of mass of a selection of atoms.

"""
function centerofmass(selection::Vector{PDBTools.Atom}) 
    T = Float64
    totmass = zero(T)
    cm = zeros(MVector{T,3})
    for atom in selection
        totmass += atom.mass
        cm += atom.mass * atom.position
    end
    cm /= totmass
    return SVector(cm)
end

function centerofmass(selection::Vector{PDBTools.Atom}, unitcell::Chemfiles.UnitCell; center=SVector(0.0, 0.0, 0.0))
    totmass = eltype(atom.position)
    cm = MVector{3}(0.0, 0.0, 0.0)
    for atom in selection
        totmass += atom.mass
        cm += atom.mass * wrap(atom.position, center, unitcell)
    end
    cm /= totmass
    return SVector(cm)
end

@testitem "centerofmass" begin
    using MolSimToolkit.Testing
    using PDBTools
    pdb = read("$test_dir/test/data/structure.pdb")
    protein = select(pdb, "protein")
    cm = centerofmass(protein)
    @test cm ≈ SVector{3}(0.0, 0.0, 0.0)
end
