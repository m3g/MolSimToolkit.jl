import Pkg
Pkg.add("Documenter")
using Documenter
push!(LOAD_PATH, "../")
push!(LOAD_PATH, "../src/")
using MolSimToolkit 
makedocs(
    modules = [MolSimToolkit],
    sitename = "MolSimToolkit.jl",
    doctest = false, 
    pages = [
        "Home" => "index.md",
        "Functions" => "functions.md",
        "Developer zone" => "Developer.md",
    ],
)
deploydocs(
    repo = "github.com/m3g/MolSimToolkit.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#"],
)
