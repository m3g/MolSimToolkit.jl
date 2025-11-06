module MolSimToolkit

import Chemfiles
import PDBTools
import OffsetArrays
import LaTeXStrings # only because Aqua complains: used in the Plotting extensions

# Names, shared between different packages
import MolSimToolkitShared: center_of_mass, 
                            distances,
                            coordination_number,
                            bulk_coordination,
                            wrap,
                            wrap_to_first,
                            align, 
                            align!, 
                            rmsd,
                            dihedral,
                            dihedrals

using TestItems: @testitem
using StaticArrays: FieldVector, SMatrix, MVector
using LinearAlgebra: norm, cross, dot, diag
using Reexport: @reexport
using ProgressMeter: Progress, next!, @showprogress
using Statistics: mean
using Printf: @sprintf

export wrap, wrap_to_first
export distances
export dihedral, dihedrals, average_dihedrals
export align, align!, rmsd, rmsd_matrix
export intermittent_correlation
export bulk_coordination
export coordination_number
export center_of_mass
export dihedral, dihedrals

# PDBTools.hydrogen_bonds is overloaded and reexported here
import PDBTools.hydrogen_bonds
export hydrogen_bonds

# Reexported from ProteinSecondaryStructures for convenience
using ProteinSecondaryStructures: dssp_run, stride_run, 
    ss_code, ss_number, ss_name
export dssp_run, stride_run, ss_code, ss_number, ss_name
# SS trajectory functions
export ss_map, ss_mean, ss_heatmap

# Version of the package: used for printing in some places
const version = pkgversion(@__MODULE__)

# Testing module
include("../test/Testing.jl")

# Data structures
include("./datastructures/Simulation.jl")
include("./datastructures/Positions.jl")

# Structural properties
include("./miscellaneous_functions/distances.jl")
include("./miscellaneous_functions/dihedrals.jl")
include("./miscellaneous_functions/center_of_mass.jl")
include("./miscellaneous_functions/most_representative_structure.jl")
include("./secondary_structure/secondary_structure.jl")
include("./hydrogen_bonds/hydrogen_bonds.jl")

# Time-dependent properties
include("./miscellaneous_functions/intermittent_correlation.jl")

# Solvation and interactions
include("./miscellaneous_functions/bulk_coordination.jl")
include("./miscellaneous_functions/coordination_number.jl")

#  Structural alignment
include("./structural_alignment/standard.jl")
include("./structural_alignment/mdlovofit.jl")

# Analysis functions and modules
include("./BlockAverages.jl")
@reexport using .BlockAverages
include("./MolecularMinimumDistances/MolecularMinimumDistances.jl")
@reexport using .MolecularMinimumDistances
include("./gromacs/remd.jl")
include("./Reweighting/Reweighting.jl")
@reexport using .Reweighting

# Simulation setup facilities
include("./PackmolInputCreator/PackmolInputCreator.jl")

# Cluster management
include("./Coaraci/Coaraci.jl")

# Structure for plotting styles
struct MolSimStyle end
export MolSimStyle

end
