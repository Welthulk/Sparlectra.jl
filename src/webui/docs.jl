const _WEBUI_DOCS_ROOT = normpath(joinpath(@__DIR__, "..", "..", "docs", "src"))

const WEBUI_HELP_TOPICS = Dict(
  "webui.casefile" => (label = "MATPOWER case file", page = "webui", heading = "PowerFlow input paths", selector = "`webui.casefile`"),
  "webui.config_file" => (label = "Configuration template file", page = "webui", heading = "PowerFlow input paths", selector = "`webui.config_file`"),
  "webui.output_root" => (label = "Output root directory", page = "webui", heading = "PowerFlow input paths", selector = "`webui.output_root`"),
  "power_flow.tol" => (label = "PowerFlow tolerance", page = "powerflow_configuration", heading = "Solver core options", selector = "`power_flow.tol`"),
  "power_flow.max_iter" => (label = "Maximum iterations", page = "powerflow_configuration", heading = "Solver core options", selector = "`power_flow.max_iter`"),
  "power_flow.autodamp" => (label = "Autodamping enabled", page = "powerflow_configuration", heading = "Solver core options", selector = "`power_flow.autodamp`"),
  "power_flow.autodamp_min" => (label = "Autodamping minimum", page = "powerflow_configuration", heading = "Solver core options", selector = "`power_flow.autodamp_min`"),
  "power_flow.qlimits.enabled" => (label = "Q-limit handling enabled", page = "powerflow_configuration", heading = "Q-limit options and guard", selector = "`power_flow.qlimits.enabled`"),
  "power_flow.wrong_branch_detection" => (label = "Wrong-branch detection", page = "configuration", heading = "Wrong-branch detection semantics (rectangular PF)", selector = ""),
  "power_flow.start_mode.angle_mode" => (label = "Start angle mode", page = "powerflow_configuration", heading = "Start mode options", selector = "`power_flow.start_mode.angle_mode`"),
  "power_flow.start_mode.voltage_mode" => (label = "Start voltage mode", page = "powerflow_configuration", heading = "Start mode options", selector = "`power_flow.start_mode.voltage_mode`"),
  "output.logfile_results" => (label = "Logfile output mode", page = "performance_profiling", heading = "Output configuration", selector = "`output.logfile_results`"),
  "benchmark.enabled" => (label = "Benchmark enabled", page = "performance_profiling", heading = "Benchmark configuration", selector = "`benchmark.enabled`"),
  "benchmark.samples" => (label = "Benchmark samples", page = "performance_profiling", heading = "Benchmark configuration", selector = "`benchmark.samples`"),
  "benchmark.seconds" => (label = "Benchmark seconds", page = "performance_profiling", heading = "Benchmark configuration", selector = "`benchmark.seconds`"),
)

const WEBUI_FORM_HELP_TOPICS = Dict(
  "casefile" => "webui.casefile",
  "config_file" => "webui.config_file",
  "output_root" => "webui.output_root",
  "power_flow_tol" => "power_flow.tol",
  "power_flow_max_iter" => "power_flow.max_iter",
  "power_flow_autodamp" => "power_flow.autodamp",
  "power_flow_autodamp_min" => "power_flow.autodamp_min",
  "power_flow_qlimits_enabled" => "power_flow.qlimits.enabled",
  "power_flow_wrong_branch_detection" => "power_flow.wrong_branch_detection",
  "power_flow_start_angle_mode" => "power_flow.start_mode.angle_mode",
  "power_flow_start_voltage_mode" => "power_flow.start_mode.voltage_mode",
  "output_logfile_results" => "output.logfile_results",
  "benchmark_enabled" => "benchmark.enabled",
  "benchmark_samples" => "benchmark.samples",
  "benchmark_seconds" => "benchmark.seconds",
)

const WEBUI_DOC_PAGES = Dict(
  "configuration" => (title = "Configuration", file = "configuration.md"),
  "powerflow_configuration" => (title = "Power-Flow Configuration", file = "powerflow_configuration.md"),
  "powerflow_service" => (title = "Local PowerFlow Service", file = "powerflow_service.md"),
  "performance_profiling" => (title = "Performance and Profiling Configuration", file = "performance_profiling.md"),
  "webui" => (title = "Local PowerFlow Web UI", file = "webui.md"),
  "feature_matrix" => (title = "Feature Matrix", file = "feature_matrix.md"),
)

"""Resolve an allowlisted Web UI help topic to its documentation metadata."""
resolve_webui_help_topic(topic::AbstractString) = get(WEBUI_HELP_TOPICS, String(topic), nothing)

"""Resolve an allowlisted documentation page to its title and Markdown file."""
resolve_webui_doc_page(page::AbstractString) = get(WEBUI_DOC_PAGES, String(page), nothing)

function _webui_doc_path(page_metadata)::String
  path = normpath(joinpath(_WEBUI_DOCS_ROOT, page_metadata.file))
  dirname(path) == _WEBUI_DOCS_ROOT || throw(ArgumentError("Documentation path is outside the allowlisted documentation root."))
  return path
end

"""Load one allowlisted Markdown document used by the local Web UI."""
function load_webui_markdown_document(page::AbstractString)::Union{String,Nothing}
  metadata = resolve_webui_doc_page(page)
  metadata === nothing && return nothing
  path = _webui_doc_path(metadata)
  return isfile(path) ? read(path, String) : nothing
end

function _webui_markdown_heading(line::AbstractString)
  matched = match(r"^(#{1,6})\s+(.+?)\s*#*\s*$", line)
  matched === nothing && return nothing
  return (level = length(matched.captures[1]), text = strip(matched.captures[2]))
end

"""Extract a Markdown heading and its content through the next peer or parent heading."""
function extract_webui_markdown_section(markdown_text::AbstractString, heading::AbstractString)::Union{String,Nothing}
  lines = split(String(markdown_text), '\n'; keepempty = true)
  start_index = nothing
  heading_level = 0
  for index in eachindex(lines)
    parsed = _webui_markdown_heading(lines[index])
    if parsed !== nothing && parsed.text == String(heading)
      start_index = index
      heading_level = parsed.level
      break
    end
  end
  start_index === nothing && return nothing

  stop_index = lastindex(lines)
  for index in (start_index + 1):lastindex(lines)
    parsed = _webui_markdown_heading(lines[index])
    if parsed !== nothing && parsed.level <= heading_level
      stop_index = index - 1
      break
    end
  end
  return strip(join(lines[start_index:stop_index], "\n"))
end

function _webui_extract_markdown_table_row(section::AbstractString, selector::AbstractString)::Union{String,Nothing}
  lines = split(String(section), '\n'; keepempty = true)
  row_index = findfirst(line -> startswith(strip(line), "|") && occursin(selector, line), lines)
  row_index === nothing && return nothing
  header_indices = findall(line -> startswith(strip(line), "|"), lines[begin:(row_index - 1)])
  length(header_indices) >= 2 || return strip(lines[row_index])
  heading_index = findfirst(line -> _webui_markdown_heading(line) !== nothing, lines)
  excerpt = String[]
  heading_index !== nothing && push!(excerpt, lines[heading_index], "")
  append!(excerpt, (lines[header_indices[1]], lines[header_indices[2]], lines[row_index]))
  return join(excerpt, "\n")
end

"""Load the Markdown excerpt configured for an allowlisted Web UI help topic."""
function load_webui_help_excerpt(topic::AbstractString)::Union{String,Nothing}
  metadata = resolve_webui_help_topic(topic)
  metadata === nothing && return nothing
  markdown_text = load_webui_markdown_document(metadata.page)
  markdown_text === nothing && return nothing
  section = extract_webui_markdown_section(markdown_text, metadata.heading)
  section === nothing && return nothing
  return isempty(metadata.selector) ? section : _webui_extract_markdown_table_row(section, metadata.selector)
end

function _webui_heading_slug(heading_html::AbstractString)::String
  text = replace(String(heading_html), r"<[^>]+>" => "")
  text = lowercase(replace(text, "&amp;" => "and", "&quot;" => "", "&#39;" => ""))
  return strip(replace(text, r"[^a-z0-9]+" => "-"), '-')
end

function _webui_rewritten_doc_href(target::AbstractString; current_page::Union{Nothing,String} = nothing)::Union{String,Nothing}
  href = String(target)
  startswith(href, "https://") && return href
  startswith(href, "http://") && return href
  startswith(href, "#") && return current_page === nothing ? nothing : href
  (startswith(href, "/") || occursin('\\', href) || occursin(':', href)) && return nothing

  relative = startswith(href, "./") ? href[3:end] : href
  matched = match(r"^([A-Za-z0-9_-]+)\.md(#[A-Za-z0-9._~:%-]+)?$", relative)
  matched === nothing && return nothing
  page = matched.captures[1]
  metadata = resolve_webui_doc_page(page)
  metadata === nothing && return nothing
  metadata.file == "$(page).md" || return nothing
  fragment = something(matched.captures[2], "")
  return "/docs/$(page)$(fragment)"
end

"""Rewrite rendered Markdown links to safe, allowlisted local documentation routes."""
function rewrite_webui_doc_links(rendered_html::AbstractString; current_page::Union{Nothing,String} = nothing)::String
  heading_pattern = r"<h([1-6])>(.*?)</h[1-6]>"s
  html = replace(String(rendered_html), heading_pattern => matched_text -> begin
    matched = match(heading_pattern, String(matched_text))
    level = matched.captures[1]
    contents = matched.captures[2]
    slug = _webui_heading_slug(contents)
    isempty(slug) ? String(matched_text) : "<h$(level) id=\"$(slug)\">$(contents)</h$(level)>"
  end)
  href_pattern = Regex("href=\"([^\"]+)\"")
  return replace(html, href_pattern => matched_text -> begin
    matched = match(href_pattern, String(matched_text))
    rewritten = _webui_rewritten_doc_href(matched.captures[1]; current_page = current_page)
    rewritten === nothing ? "aria-disabled=\"true\"" : "href=\"$(rewritten)\""
  end)
end

"""Render trusted repository Markdown as HTML using Julia's Markdown standard library."""
function render_webui_markdown(markdown_text::AbstractString; current_page::Union{Nothing,String} = nothing)::String
  io = IOBuffer()
  show(io, MIME"text/html"(), Markdown.parse(String(markdown_text)))
  return rewrite_webui_doc_links(String(take!(io)); current_page = current_page)
end
