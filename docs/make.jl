using Pkg
Pkg.instantiate()
using Documenter
using Sparlectra

DocMeta.setdocmeta!(Sparlectra, :DocTestSetup, :(using Sparlectra); recursive = true)

makedocs(
  sitename = "Sparlectra.jl",
  modules = [Sparlectra],
  clean = true,
  doctest = true,
  checkdocs = :none,
  format = Documenter.HTML(assets = ["assets/tablestyle.css"], prettyurls = get(ENV, "CI", "false") == "true", collapselevel = 1, canonical = "https://welthulk.github.io/Sparlectra.jl"),
  pages = [
    "Home" => "index.md",
    "Feature Matrix" => "feature_matrix.md",
    "Changelog" => "changelog.md",
    "Branch Model" => "branchmodel.md",
    "External Solvers" => "external_solvers.md",
    "Import" => "import.md",
    "Links" => "links.md",
    "Network Reports" => "netreports.md",
    "Power Limits" => "powerlimits.md",
    "Solver" => "solver.md",
    "Voltage Dependent Control" => "voltage_dependent_control.md",
    "Transformer Control" => "transformer_control.md",
    "State Estimation" => "state_estimation.md",
    "Workshop" => "workshop.md",
    "API Reference" => "reference.md",
  ],
)

deploydocs(repo = "github.com/welthulk/Sparlectra.jl", devbranch = "main", push_preview = true)
