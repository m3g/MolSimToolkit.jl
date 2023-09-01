import Pkg
Pkg.add("Documenter")
using Documenter
using MolSimToolkit 
push!(LOAD_PATH, "../src/")
makedocs(
    modules = [MolSimToolkit],
    sitename = "MolSimToolkit.jl",
    pages = [
        "Home" => "index.md",
        "Developer zone" => "Developer.md",
    ],
)
deploydocs(
    repo = "github.com/m3g/MolSimToolkit.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#"],
)
