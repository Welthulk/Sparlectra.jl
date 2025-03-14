using Documenter, Sparlectra

makedocs(
  sitename = "Sparlectra.jl",
  modules = [Sparlectra],
  clean = false,
  doctest = true,
  checkdocs = :none,  # Skip checking for missing docstrings
  format = Documenter.HTML(assets = ["assets/tablestyle.css"], prettyurls = get(ENV, "CI", nothing) == "true", collapselevel = 1),
  pages = ["index.md", "workshop.md", "branchmodel.md", "changelog.md", raw"API" => ["reference.md", "component_types.md", "networks.md", "import.md", "remove_functions.md"]],
)

# Optionally deploy docs
deploydocs(repo = "github.com/welthulk/Sparlectra.jl.git")