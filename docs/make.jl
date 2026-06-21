using Pkg
Pkg.instantiate()
using Documenter
using Sparlectra
using TOML

project_toml = TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
sparlectra_version = project_toml["version"]

DocMeta.setdocmeta!(Sparlectra, :DocTestSetup, :(using Sparlectra); recursive = true)

makedocs(
  sitename = "Sparlectra.jl v$(sparlectra_version)",
  repo = "https://github.com/Welthulk/Sparlectra.jl/blob/{commit}{path}#L{line}",
  modules = [Sparlectra],
  clean = true,
  doctest = true,
  checkdocs = :none,
  format = Documenter.HTML(assets = ["assets/tablestyle.css"], prettyurls = get(ENV, "CI", "false") == "true", collapselevel = 1, canonical = "https://welthulk.github.io/Sparlectra.jl"),
  pages = [
    "Home" => "index.md",
    "Feature Matrix" => "feature_matrix.md",
    "Changelog" => "changelog.md",
    "Configuration" => "configuration.md",
    "Programmatic API" => "programmatic_api.md",
    "Local PowerFlow Service" => "powerflow_service.md",
    "Local PowerFlow Web UI" => "webui.md",
    "Power-Flow Configuration" => "powerflow_configuration.md",
    "Q-limit Switching Strategy" => "q_limit_switching_strategy.md",
    "MATPOWER Format" => "matpower_format.md",
    "MATPOWER Import Configuration" => "matpower_import.md",
    "State-Estimation Configuration" => "state_estimation_configuration.md",
    "Performance and Profiling" => "performance_profiling.md",
    "Tests" => "tests.md",
    "Branch Model" => "branchmodel.md",
    "External Solvers" => "external_solvers.md",
    "Import" => "import.md",
    "Links" => "links.md",
    "Network Reports" => "netreports.md",
    "Power Limits" => "powerlimits.md",
    "Solver" => "solver.md",
    "Synthetic Tiled Grids" => "synthetic_grids.md",
    "Voltage Dependent Control" => "voltage_dependent_control.md",
    "Control Framework" => "control_framework.md",
    "Transformer Control" => "transformer_control.md",
    "Examples Overview" => "examples_overview.md",
    "State Estimation" => "state_estimation.md",
    "Workshop" => "workshop.md",
    "Reference" => [
      "Overview" => "reference.md",
      "API and Web UI Service" => "reference_api.md",
      "ACPFlow Runner" => "reference_acpflow.md",
      "Network Model" => "reference_network.md",
      "Rectangular Power Flow" => "reference_powerflow_rectangular.md",
      "Import and Export" => "reference_import_export.md",
      "State Estimation and Measurements" => "reference_state_estimation.md",
    ],
  ],
)
