module Ressampling

using ..MolSimToolkit 
using PDBTools
using TestItems

export ressample

""" 
    ressample(trajectory, perturbation::Function)

"""
function ressample end

# Tests
include("./test/runtests.jl")

end # module PackmolInputCreator

