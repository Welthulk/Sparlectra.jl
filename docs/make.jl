using Documenter, Sparlectra

makedocs(
  sitename="Sparlectra.jl",  
  clean = false,
  doctest = false,
  format = Documenter.HTML( assets=["assets/tablestyle.css"],
        prettyurls = get(ENV, "CI", nothing) == "true",
        collapselevel = 1),
  pages = [
    "index.md",
    "API" => ["networks.md"],
    "branchmodel.md",
    "changelog.md"
  ]
)

deploydocs(
    repo = "github.com/Welthulk/Sparlectra.jl.git",
    target = "build",
)
