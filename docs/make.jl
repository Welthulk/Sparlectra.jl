using Documenter
using Sparlectra

DocMeta.setdocmeta!(Sparlectra, :DocTestSetup, :(using Sparlectra); recursive = true)

makedocs(
  sitename = "Sparlectra.jl",
  modules = [Sparlectra],
  clean = false,
  doctest = true,
  checkdocs = :none,
  format = Documenter.HTML(; assets = ["assets/tablestyle.css"], prettyurls = get(ENV, "CI", "false") == "true", collapselevel = 1, canonical = "https://welthulk.github.io/Sparlectra.jl"),
  pages = ["index.md", "feature_matrix.md", "changelog.md", "branchmodel.md", "external_solvers.md", "import.md", "links.md", "netreports.md", "powerlimits.md", "solver.md", "voltage_dependent_control.md", "state_estimation.md", "workshop.md", raw"API" => ["reference.md"]],
)

deploydocs(; repo = "github.com/welthulk/Sparlectra.jl", devbranch = "main", push_preview = true)
