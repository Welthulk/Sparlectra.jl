# Copyright 2023-2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# file: docs/make.jl
# purpose: Documenter build script for the Sparlectra documentation site,
#          including doctests; the acceptance gate for docs-touching changes
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
  # Discoverability: the sitename doubles as the <title> suffix on every
  # page. The longer "Julia AC Power Flow and State Estimation" variant
  # exceeded 70 characters together with subpage-name prefixes, so the
  # short disambiguation form is used (task 2.1 fallback).
  sitename = "Sparlectra.jl v$(sparlectra_version), Julia power flow",
  repo = "https://github.com/Welthulk/Sparlectra.jl/blob/{commit}{path}#L{line}",
  modules = [Sparlectra],
  clean = true,
  doctest = true,
  checkdocs = :none,
  format = Documenter.HTML(
    assets = ["assets/tablestyle.css"],
    # Site-wide meta description (also emitted as og:description and
    # twitter:description by Documenter). 153 characters: the 160-char
    # budget forced dropping the vessel-disambiguation sentence AND the
    # word "open-source" (task 2.2). Documenter also emits og:image tags
    # automatically from assets/preview.png plus the canonical URL.
    description = "Sparlectra.jl: Julia package for AC power flow (Newton-Raphson), DC power flow, WLS state estimation, IEC 60909 short circuit, CGMES and MATPOWER import.",
    prettyurls = get(ENV, "CI", "false") == "true",
    collapselevel = 1,
    canonical = "https://welthulk.github.io/Sparlectra.jl",
    repolink = "https://github.com/Welthulk/Sparlectra.jl",
    # Raised size limits: the single-page changelog and the autodocs network
    # reference legitimately exceed Documenter's 100 KiB default.
    size_threshold_warn = 150 * 1024,
    size_threshold = 300 * 1024,
    # 700 KiB left only ~0.5 KiB headroom; the Literate-generated workshop
    # page pushed the index to ~706 KiB.
    search_size_threshold_warn = 800 * 1024,
  ),
  pages = [
    "Home" => "index.md",
    "Feature Matrix" => "feature_matrix.md",
    "Changelog" => "changelog.md",
    "Configuration" => "configuration.md",
    "Programmatic API" => "programmatic_api.md",
    "Local PowerFlow Service" => "powerflow_service.md",
    "Local PowerFlow Web UI" => "webui.md",
    "Fast Start (Sysimage)" => "fast_start.md",
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
    "Slack and External Grid Sources" => "slack_vs_source.md",
    "Synthetic Tiled Grids" => "synthetic_grids.md",
    "Voltage Dependent Control" => "voltage_dependent_control.md",
    "Control Framework" => "control_framework.md",
    "Remote Voltage Control" => "remote_voltage_control.md",
    "Series Compensation (TCSC)" => "series_compensation.md",
    "HVDC Back-to-Back" => "hvdc_back_to_back.md",
    "Short-Circuit Compendium" => "short_circuit.md",
    "N-1 Contingency Analysis" => "contingency.md",
    "Examples Overview" => "examples_overview.md",
    "State Estimation" => "state_estimation.md",
    "Workshop" => "workshop.md",
    # Literate.jl-generated pages (committed output of docs/generate_notebooks.jl,
    # not rebuilt here — edit docs/lit/*.jl and regenerate). New notebooks get
    # a nav entry inside this group.
    "Notebooks" => [
      "Workshop Tour (all in one)" => "generated/workshop_tour.md",
      "Try it in your Browser" => "generated/workshop_intro.md",
      "Slack Types and Short Circuit" => "generated/workshop_slack_short_circuit.md",
      "Distributed Slack" => "generated/workshop_distributed_slack.md",
      "Transformer Taps" => "generated/workshop_transformers.md",
      "State Estimation Basics" => "generated/workshop_state_estimation.md",
      "TCSC Flow Steering" => "generated/workshop_series_compensation.md",
    ],
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

const SITE_CANONICAL = "https://welthulk.github.io/Sparlectra.jl"

"""
    write_sitemap(build_dir, canonical)

Post-build sitemap generator (discoverability task 1.3): walk the built
site, collect every rendered page (`index.html` per directory in
prettyurls mode, every non-asset `*.html` otherwise), and write
`sitemap.xml` at the build root with `<loc>` entries under the canonical
URL and `<lastmod>` set to the build date. `search/` and `assets/` never
carry content pages and are skipped; `generated/` pages are part of the
nav and are included. No package dependency: the file is plain string
assembly, well-formedness is asserted by the acceptance check.
"""
function write_sitemap(build_dir::AbstractString, canonical::AbstractString)
  base = rstrip(canonical, '/')
  lastmod = Libc.strftime("%Y-%m-%d", time())
  locs = String[]
  for (root, _dirs, files) in walkdir(build_dir)
    rel = relpath(root, build_dir)
    parts = rel == "." ? String[] : splitpath(rel)
    if !isempty(parts) && first(parts) in ("assets", "search")
      continue
    end
    for file in files
      endswith(file, ".html") || continue
      if file == "index.html"
        push!(locs, rel == "." ? "$(base)/" : "$(base)/" * join(parts, "/") * "/")
      else
        push!(locs, "$(base)/" * (rel == "." ? file : join(parts, "/") * "/" * file))
      end
    end
  end
  sort!(unique!(locs))
  open(joinpath(build_dir, "sitemap.xml"), "w") do io
    println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
    println(io, "<urlset xmlns=\"http://www.sitemaps.org/schemas/sitemap/0.9\">")
    for loc in locs
      println(io, "  <url><loc>", loc, "</loc><lastmod>", lastmod, "</lastmod></url>")
    end
    println(io, "</urlset>")
  end
  return joinpath(build_dir, "sitemap.xml")
end

# Root placement (task 1.4): Documenter copies docs/src/assets/* to
# build/assets/; robots.txt and llms.txt must live at the SITE ROOT to be
# honored by crawlers. sitemap.xml is generated directly at the root.
let build_dir = joinpath(@__DIR__, "build")
  for rootfile in ("robots.txt", "llms.txt")
    src = joinpath(build_dir, "assets", rootfile)
    isfile(src) || error("expected $(src) from docs/src/assets/ in the build")
    mv(src, joinpath(build_dir, rootfile); force = true)
  end
  sitemap = write_sitemap(build_dir, SITE_CANONICAL)
  @info "discoverability artifacts in place" robots = joinpath(build_dir, "robots.txt") llms = joinpath(build_dir, "llms.txt") sitemap
end
