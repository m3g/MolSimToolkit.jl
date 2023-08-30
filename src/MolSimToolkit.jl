module MolSimToolkit

using TestItems
using StaticArrays
import Chemfiles
import PDBTools
import LinearAlgebra: norm

export distances

# Flag for internal function doc entries
const INTERNAL = "Internal function or structure - interface may change."

# Testing module
include("./Testing.jl")

# Basic functions
include("./wrap.jl")
include("./distances.jl")
include("./centerofmass.jl")

end
