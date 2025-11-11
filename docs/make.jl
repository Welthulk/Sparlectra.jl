# docs/make.jl
using Pkg
Pkg.activate(@__DIR__)  # aktiviere das docs-Project.toml
# binde dein Paket aus dem Repo-Root ins docs-Env ein
Pkg.develop(path=joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Documenter
using Sparlectra

DocMeta.setdocmeta!(Sparlectra, :DocTestSetup, :(using Sparlectra); recursive=true)

makedocs(
    sitename = "Sparlectra.jl",
    modules  = [Sparlectra],
    clean    = false,
    doctest  = true,
    checkdocs = :none,
    format   = Documenter.HTML(;
        assets     = ["assets/tablestyle.css"],
        prettyurls = get(ENV, "CI", "false") == "true",
        collapselevel = 1,
        canonical  = "https://welthulk.github.io/Sparlectra.jl",
    ),
    pages = [
        "index.md",
        "workshop.md",
        "branchmodel.md",
        "changelog.md",
        raw"API" => [
            "reference.md",
            "component_types.md",
            "networks.md",
            "import.md",
            "remove_functions.md",
        ],
    ],
)

deploydocs(;
    repo         = "github.com/welthulk/Sparlectra.jl", # ohne .git
    devbranch    = "main",    # ggf. "master"
    push_preview = true,
)
