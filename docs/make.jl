import Pkg
Pkg.add("Documenter")
using Documenter
using Plots
push!(LOAD_PATH, "../")
push!(LOAD_PATH, "../src/")
using MolSimToolkit
using MolSimToolkitShared
makedocs(
    modules=[MolSimToolkit, MolSimToolkitShared],
    sitename="MolSimToolkit.jl",
    doctest=false,
    pages=[
        "Home" => "index.md",
        "Structural analyses" => Any[
            "Distances and misc." => "Structural_properties.md",
            "Dihedral angle analysis" => "Dihedrals.md",
            "Secondary structure" => "secondary_structures.md",
            "Structural alignment" => "structural_alignment.md",
        ],
        "Simulation statistics" => Any[ 
            " Block averages" => "block_averages.md",
            " Replica exchange" => "remd.md",
        ],
        "Interactions" => Any[
            "Coordination numbers" => "Solvation_and_interactions.md",
            "Molecular Minimum Distances" => "molecular_minimum_distances.md",
        ],
        "Time-dependent properties" => Any[
            "Intermittent correlation" => "intermittent_correlation.md",
        ],
        "System setup" => "system_setup.md",
        "Plotting style" => "plotting_style.md",
        "Developer zone" => "Developer.md",
        "Experimental" => Any[
            "Simulation Reweighting" => "Reweighting.md",
            "Cluster submission management" => "Coaraci.md",
        ], 
    ],
)
deploydocs(
    repo="github.com/m3g/MolSimToolkit.jl.git",
    target="build",
    branch="gh-pages",
    versions=["stable" => "v^", "v#.#"],
)
