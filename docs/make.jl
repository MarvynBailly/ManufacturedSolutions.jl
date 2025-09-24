using Documenter
using ManufacturedSolutions


DocMeta.setdocmeta!(ManufacturedSolutions, :DocTestSetup, :(using ManufacturedSolutions); recursive=true)

makedocs(
    modules  = [ManufacturedSolutions],
    sitename = "ManufacturedSolutions.jl",
    authors  = "Marvyn Bailly",
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://MarvynBailly.github.io/ManufacturedSolutions.jl;/stable/",
        assets=String[],
    ),
    
    
    
    pages = [
        "Home" => "index.md",
    ],

    # Optional quality gates once you’re ready:
    checkdocs = :exports,
)

deploydocs(
    repo      = "github.com/MarvynBailly/ManufacturedSolutions",
    devbranch = "main",
    versions = ["stable" => "v^", "v#.#.#"] 
)