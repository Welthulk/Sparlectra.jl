const _WEBUI_DOCS_ROOT = normpath(joinpath(@__DIR__, "..", "..", "docs", "src"))

const WEBUI_HELP_TOPICS = Dict(
  "power_flow.start_mode.voltage_mode" => (
    label = "Start voltage mode",
    page = "powerflow_configuration",
    heading = "Start mode options",
    selector = "`power_flow.start_mode.voltage_mode`",
  ),
)

const WEBUI_DOC_PAGES = Dict(
  "configuration" => (title = "Configuration", file = "configuration.md"),
  "powerflow_configuration" => (title = "Power-Flow Configuration", file = "powerflow_configuration.md"),
  "powerflow_service" => (title = "Local PowerFlow Service", file = "powerflow_service.md"),
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

"""Render trusted repository Markdown as HTML using Julia's Markdown standard library."""
function render_webui_markdown(markdown_text::AbstractString)::String
  io = IOBuffer()
  show(io, MIME"text/html"(), Markdown.parse(String(markdown_text)))
  return String(take!(io))
end
