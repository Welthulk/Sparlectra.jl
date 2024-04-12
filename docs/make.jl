using Documenter, Sparlectra

makedocs(
  sitename="Sparlectra.jl",  
  modules = [Sparlectra],
  clean = false,
  doctest = false,
  format = Documenter.HTML( assets=["assets/tablestyle.css"],
        prettyurls = get(ENV, "CI", nothing) == "true",
        collapselevel = 1),
  pages = [
    "index.md",    
    "workshop.md",
    "API" => ["reference.md", "networks.md", "import.md"],    
    "branchmodel.md",
    "changelog.md"
  ]
)
