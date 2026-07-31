# Self-contained environment setup: activate the docs project regardless of
# how this script is started (julia docs/make.jl, --project=., REPL include),
# and resolve so dependency changes in the main package are picked up.
using Pkg
Pkg.activate(@__DIR__)
Pkg.resolve()
Pkg.instantiate()
using Documenter
using Sparlectra
using TOML

project_toml = TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
sparlectra_version = project_toml["version"]

# warn=false: repeated include("docs/make.jl") in one REPL session re-sets an
# identical DocTestSetup; the "already set, overwriting" warning is noise.
DocMeta.setdocmeta!(Sparlectra, :DocTestSetup, :(using Sparlectra); recursive = true, warn = false)

makedocs(
  sitename = "Sparlectra.jl v$(sparlectra_version)",
  repo = "https://github.com/Welthulk/Sparlectra.jl/blob/{commit}{path}#L{line}",
  modules = [Sparlectra],
  clean = true,
  doctest = true,
  checkdocs = :none,
  format = Documenter.HTML(
    assets = ["assets/tablestyle.css"],
    prettyurls = get(ENV, "CI", "false") == "true",
    collapselevel = 1,
    canonical = "https://welthulk.github.io/Sparlectra.jl",
    repolink = "https://github.com/Welthulk/Sparlectra.jl",
    # Raised size limits: the single-page changelog and the autodocs network
    # reference legitimately exceed Documenter's 100 KiB default.
    size_threshold_warn = 150 * 1024,
    size_threshold = 300 * 1024,
    search_size_threshold_warn = 700 * 1024,
  ),
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
    "DTF Format" => "dtf_format.md",
    "MATPOWER Import Configuration" => "matpower_import.md",
    "CGMES Import" => "cgmes_import.md",
    "CGMES Export" => "cgmes_export.md",
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
    "Remote Voltage Control" => "remote_voltage_control.md",
    "Short-Circuit Compendium" => "short_circuit.md",
    "Examples Overview" => "examples_overview.md",
    "State Estimation" => "state_estimation.md",
    "Workshop" => "workshop.md",
    "Reference" => [
      "Overview" => "reference.md",
      "API and Web UI Service" => "reference_api.md",
      "ACPFlow Runner" => "reference_acpflow.md",
      "Network Model" => "reference_network.md",
      "Rectangular Power Flow" => "reference_powerflow_rectangular.md",
      "DC Power Flow" => "reference_powerflow_dc.md",
      "Short Circuit" => "reference_short_circuit.md",
      "Import and Export" => "reference_import_export.md",
      "State Estimation and Measurements" => "reference_state_estimation.md",
    ],
  ],
)
