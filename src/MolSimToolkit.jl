module MolSimToolkit

import Chemfiles
import PDBTools
import OffsetArrays
import LaTeXStrings # only because Aqua complains: used in the Plotting extensions

using TestItems: @testitem
using AtomsBase: atomic_mass
using StaticArrays: FieldVector, SMatrix, MVector
using LinearAlgebra: norm
using Reexport: @reexport
using ProgressMeter: Progress, next!

export wrap, wrap_to_first
export distances
export align, align!, rmsd, rmsd_matrix
export intermittent_correlation
export bulk_coordination

# New method added to this function, which is reexported
import PDBTools.center_of_mass
export center_of_mass


# Testing module
include("./Testing.jl")

# Data structures
include("./datastructures/Simulation.jl")
include("./datastructures/Positions.jl")

# Coordinate PBC wrapping functions
include("./wrap.jl")

# Miscelaneous functions
include("./miscelaneous_functions/distances.jl")
include("./miscelaneous_functions/center_of_mass.jl")
include("./miscelaneous_functions/intermittent_correlation.jl")
include("./miscelaneous_functions/bulk_coordination.jl")

#  Structural alignment
include("./procrustes.jl")

# Analysis functions and modules
include("./BlockAverages.jl")
@reexport using .BlockAverages
include("./MolecularMinimumDistances/MolecularMinimumDistances.jl")
@reexport using .MolecularMinimumDistances
include("./gromacs/remd.jl")
include("./Reweighting/Reweighting.jl")

# Simulation setup facilities
include("./PackmolInputCreator/PackmolInputCreator.jl")

# Structure for plotting styles
struct MolSimStyle end
export MolSimStyle

end
