module MolSimToolkit

using TestItems
using StaticArrays
import Chemfiles
using Printf
import PDBTools
import LinearAlgebra: norm

export wrap, wrap_to_first
export distances

# Flag for internal function doc entries
const INTERNAL = "Internal function or structure - interface may change."

# Testing module
include("./Testing.jl")

# Data structures
include("./datastructures/Simulation.jl")
include("./datastructures/Positions.jl")

# Basic functions
include("./wrap.jl")
include("./distances.jl")
include("./center_of_mass.jl")

end
