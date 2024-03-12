module MolSimToolkit

using Printf
using TestItems
using AtomsBase
using StaticArrays
import Chemfiles
import PDBTools
import OffsetArrays
using LinearAlgebra: norm
using Reexport: @reexport
using ProgressMeter: Progress, next!

import LaTeXStrings # only because Aqua complains: used in the Plotting extensions

export wrap, wrap_to_first
export distances
export align, align!, rmsd, rmsd_matrix
export intermittent_correlation

# Testing module
include("./Testing.jl")

# Data structures
include("./datastructures/Simulation.jl")
include("./datastructures/Positions.jl")

# Basic functions
include("./wrap.jl")
include("./simple_functions/distances.jl")
include("./simple_functions/center_of_mass.jl")
include("./procrustes.jl")
include("./simple_functions/intermittent_correlation.jl")

# Analysis functions and modules
include("./BlockAverages.jl")
@reexport using .BlockAverages
include("./MolecularMinimumDistances/MolecularMinimumDistances.jl")
@reexport using .MolecularMinimumDistances
include("./gromacs/remd.jl")

# Simulation setup facilities
include("./PackmolInputCreator/PackmolInputCreator.jl")

end
