using Documenter
using Plots
using MolSimToolkit
using MolSimToolkit.PackmolInputCreator
using PDBToolsBase
ENV["LINES"] = 10
ENV["COLUMNS"] = 120
ENV["GKSwstype"] = "nul"
makedocs(
    modules = [
        MolSimToolkit, PDBToolsBase,
        isdefined(Base, :get_extension) ? Base.get_extension(MolSimToolkit, :Plotting) : MolSimToolkit.Plotting
    ],
    # PDBToolsBase is listed for the docstrings of the shared functions used here; it exports
    # nothing, so checking only the exported names does not require all its docstrings
    checkdocs=:exports,
    # Source links of the docstrings of PDBToolsBase
    remotes=Dict(pkgdir(PDBToolsBase) => (Documenter.Remotes.GitHub("m3g", "PDBToolsBase.jl"), "main")),
    sitename="MolSimToolkit.jl",
    doctest=false,
    pages=[
        "Home" => "index.md",
        "Simulation" => "simulation.md",
        "Structural analyses" => Any[
            "Hydrogen bonds" => "hydrogen_bonds.md",
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
        "Dynamical properties" => Any[
            "Site occupancy" => "time_dependent/occupancy.md",
            "Intermittent correlation" => "time_dependent/intermittent_correlation.md",
            "Self-diffusion coefficient" => "time_dependent/self_diffusion.md",
            "Viscosity (Green-Kubo)" => "time_dependent/viscosity.md",
        ],
        "System setup" => "system_setup.md",
        "Plotting style" => "plotting_style.md",
#        "Developer zone" => "Developer.md",
        "Experimental" => Any[
            "Simulation Reweighting" => "Reweighting.md",
            "Cluster submission management" => "Coaraci.md",
            "Transfer free energies" => "mvalues.md",
        ], 
    ],
)
deploydocs(
    repo="github.com/m3g/MolSimToolkit.jl.git",
    target="build",
    branch="gh-pages",
    versions=["stable" => "v^", "v#.#"],
)
