module MolSimToolkit

using Printf
using TestItems
using AtomsBase
using StaticArrays
import Chemfiles
import PDBTools
import LinearAlgebra: norm
import Reexport: @reexport

import LaTeXStrings # only because Aqua complains: used in the Plotting extensions

export wrap, wrap_to_first
export distances

# Testing module
include("./Testing.jl")

# Data structures
include("./datastructures/Simulation.jl")
include("./datastructures/Positions.jl")

# Basic functions
include("./wrap.jl")
include("./distances.jl")
include("./center_of_mass.jl")

# Analysis functions and modules
include("./BlockAverages.jl")
@reexport using .BlockAverages
include("./MolecularMinimumDistances/MolecularMinimumDistances.jl")
@reexport using .MolecularMinimumDistances
include("./gromacs/remd.jl")

# Simulation setup facilities
include("./PackmolInputCreator/PackmolInputCreator.jl")

end
