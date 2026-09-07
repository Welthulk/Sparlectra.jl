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

"""
    generate_index_page()

Write `docs/src/index.md` from `README.md` plus `docs/index_footer.md`.

The Home page used to embed the README through an `@eval` block calling
`Markdown.parse`. That is Julia's BASE Markdown parser, and it has no notion
of inline HTML: the logo, which the README places as
`<a href="..."><img align="left" ...></a>` so the intro text wraps around it,
came out on the built site as its own literal source text (reported
2026-09-07). Generating a real Markdown page instead hands the content to
Documenter's own parser, which passes raw HTML through, so the page renders
the way the README does on GitHub.

The README addresses images and pages relative to the repository ROOT, and on
the built site those targets mean two different things. A `docs/src/...`
target is a page of this very site and becomes a page-relative link;
everything else (`LICENSE`, `CONTRIBUTING.md`, `tools/install_webui.sh`,
`docs/lit/`) is a repository file that the site does not contain at all, so it
becomes an absolute GitHub URL. Those seven links used to render as dead
relative links, which nothing noticed because `Markdown.parse` output never
reached Documenter's cross-reference check. Raw HTML `img`/`a` attributes are
rewritten as well: the Markdown-link pattern does not see an `<img>` tag.
"""
function generate_index_page()
  repo_root = normpath(joinpath(@__DIR__, ".."))
  readme = read(joinpath(repo_root, "README.md"), String)
  readme = replace(readme, r"\(docs/src/([^)]+)\)" => s"(\1)")
  readme = replace(readme, r"src=\"docs/src/([^\"]+)\"" => s"src=\"\1\"")
  readme = replace(readme, r"href=\"docs/src/([^\"]+)\"" => s"href=\"\1\"")

  # A relative target that names an existing repository path, and only such a
  # target: an unknown one is left alone so a typo still fails the build
  # instead of turning into a plausible-looking 404 on github.com.
  function repository_url(target::AbstractString)
    (occursin("://", target) || startswith(target, '#') || startswith(target, "mailto:")) && return nothing
    path = joinpath(repo_root, target)
    isdir(path) && return "https://github.com/Welthulk/Sparlectra.jl/tree/main/" * target
    isfile(path) && return "https://github.com/Welthulk/Sparlectra.jl/blob/main/" * target
    return nothing
  end
  readme = replace(readme, r"\]\(([^)\s]+)\)" => function (whole)
    target = whole[3:(end-1)]
    url = repository_url(target)
    return url === nothing ? whole : "](" * url * ")"
  end)


  # Documenter escapes inline HTML; a raw tag only survives inside a
  # ```@raw html fence. The README carries such tags as whole lines of their
  # own (the logo, which floats left so the intro text wraps around it), so
  # each standalone HTML line becomes a raw block. `<https://...>` is a
  # Markdown autolink, not a tag, and must not be wrapped.
  lines = split(readme, '\n'; keepempty = true)
  wrapped = String[]
  for line in lines
    if occursin(r"^\s*<[a-zA-Z][a-zA-Z0-9]*(\s[^>]*)?>", line) && !occursin(r"^\s*<https?://", line)
      push!(wrapped, "```@raw html", String(line), "```")
    else
      push!(wrapped, String(line))
    end
  end
  readme = join(wrapped, "\n")

  footer = read(joinpath(@__DIR__, "index_footer.md"), String)
  target = joinpath(@__DIR__, "src", "index.md")
  open(target, "w") do io
    # the banner is for whoever opens the file and wonders why their edit
    # disappeared; EditURL sends the "Edit on GitHub" link to the source
    println(io, "```@meta")
    println(io, "EditURL = \"../../README.md\"")
    println(io, "```")
    println(io)
    # The banner goes through @raw html, not into the Markdown body: a bare
    # `<!-- ... -->` is not a comment to Documenter's parser, it renders as
    # visible page text (verified on the built page).
    println(io, "```@raw html")
    println(io, "<!-- GENERATED by docs/make.jl from README.md and docs/index_footer.md. Do not edit; edit those two instead. -->")
    println(io, "```")
    println(io)
    print(io, readme)
    print(io, footer)
  end
  return target
end

generate_index_page()

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
  # every exported docstring must be rendered somewhere, or the build fails;
  # :none let six dead Pages entries drop their docstrings silently
  checkdocs = :exports,
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
    "Integration Guide" => "integration.md",
    "Sparlectra Case Format" => "scf.md",
    "Shipped Demo Cases" => "demo_cases.md",
    "Programmatic API" => "programmatic_api.md",
    "Local PowerFlow Service" => "powerflow_service.md",
    "Local PowerFlow Web UI" => "webui.md",
    "Sysimage" => "sysimage.md",
    "Power-Flow Configuration" => "powerflow_configuration.md",
    "Q-limit Switching Strategy" => "q_limit_switching_strategy.md",
    "MATPOWER Format" => "matpower_format.md",
    "DTF Format" => "dtf_format.md",
    "MATPOWER Import Configuration" => "matpower_import.md",
    "CGMES Import" => "cgmes_import.md",
    "CGMES Export" => "cgmes_export.md",
    "State-Estimation Configuration" => "state_estimation_configuration.md",
    "Performance and Profiling" => "performance_profiling.md",
    "Parallel Execution" => "parallel_execution.md",
    "Tests" => "tests.md",
    "Branch Model" => "branchmodel.md",
    "Component Types" => "component_types.md",
    "Component Removal" => "remove_functions.md",
    "External Solvers" => "external_solvers.md",
    "Import" => "import.md",
    "Links" => "links.md",
    "Network Reports" => "netreports.md",
    "Power Limits" => "powerlimits.md",
    "Solver" => "solver.md",
    "Slack and External Grid Sources" => "slack_vs_source.md",
    "Synthetic Tiled Grids" => "synthetic_grids.md",
    "Voltage Dependent Control" => "voltage_dependent_control.md",
    "FACTS Devices" => "facts.md",
    "Control Framework" => "control_framework.md",
    "Remote Voltage Control" => "remote_voltage_control.md",
    "Series Compensation (TCSC)" => "series_compensation.md",
    "HVDC Back-to-Back" => "hvdc_back_to_back.md",
    "Short-Circuit Compendium" => "short_circuit.md",
    "N-1 Contingency Analysis" => "contingency.md",
    "Examples Overview" => "examples_overview.md",
    "State Estimation" => "state_estimation.md",
    "Observability" => "observability.md",
    # Literate.jl-generated pages (committed output of docs/generate_notebooks.jl,
    # not rebuilt here — edit docs/lit/*.jl and regenerate). New notebooks get
    # a nav entry inside this group.
    "Notebooks" => [
      "Workshop Tour (basic)" => "generated/workshop_tour.md",
      "Workshop Tour (advanced)" => "generated/workshop_tour_advanced.md",
      "Workshop Tour (part 3: coordinated control)" => "generated/workshop_tour_control.md",
      "Workshop Tour (CGMES)" => "generated/workshop_tour_cgmes.md",
      "Slack Types and Short Circuit" => "generated/workshop_slack_short_circuit.md",
      "Distributed Slack" => "generated/workshop_distributed_slack.md",
      "Transformer Taps" => "generated/workshop_transformers.md",
      "State Estimation Basics" => "generated/workshop_state_estimation.md",
      "State Estimation Diagnostics" => "generated/workshop_se_diagnostics.md",
      "Taps, Phasors and Topology" => "generated/workshop_se_taps.md",
      "TCSC Flow Steering" => "generated/workshop_series_compensation.md",
      "Scenarios and Screening" => "generated/workshop_scenarios.md",
    ],
    "Reference" => [
      # one page per source directory; surviving pages keep their position
      "Overview" => "reference.md",
      "API and Services" => "reference_api.md",
      "ACPFlow Runner" => "reference_acpflow.md",
      "Core Model" => "reference_core.md",
      "Rectangular Power Flow" => "reference_powerflow_rectangular.md",
      "DC Power Flow" => "reference_powerflow_dc.md",
      "Short Circuit" => "reference_shortcircuit.md",
      "Format Adapters" => "reference_adapters.md",
      "State Estimation and Measurements" => "reference_stateestimation.md",
      "Configuration Internals" => "reference_config.md",
      "Controllers" => "reference_controller.md",
      "APSLF Bridge" => "reference_apslf.md",
      "Shared Numerics" => "reference_numerics.md",
      "Contingency" => "reference_contingency.md",
      "Synthetic Grids" => "reference_synthetic.md",
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
