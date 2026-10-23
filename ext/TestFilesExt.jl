module TestFilesExt

using LazyArtifacts: LazyArtifacts, @artifact_str
import MolSimToolkit.Testing: test_file

function test_file(name::AbstractString)
    # The files of the protein-TMAO system are shared with ComplexMixtures.jl
    dir = startswith(name, "protein_tmao.") ? artifact"protein_tmao_test_files" : artifact"molsimtoolkit_test_files"
    return joinpath(dir, name)
end

end
