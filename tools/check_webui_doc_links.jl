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
#
# file: tools/check_webui_doc_links.jl
# purpose: verify the Web UI help registry (WEBUI_HELP_TOPICS in
#          app/src/webui/docs.jl) against the documentation sources: every
#          topic has a hover hint (non-empty, at most 160 characters, plain
#          text without backticks or dashes) and every `doc` target names a
#          page that exists in docs/src with a heading that carries the
#          anchor as a Documenter `@id`. Part of the docs gate, so a
#          documentation change that renames an anchor or drops a section
#          fails there. No server, no solve; the check is a plain function
#          so the test suite can call it too.
#          Usage: julia --project=app tools/check_webui_doc_links.jl

using SparlectraApp

"""
    check_webui_doc_links(; docs_src, io) -> Vector{String}

Return the list of failures (empty when every topic is sound). `docs_src`
is the documentation source directory (default: the repository's
`docs/src` next to the application package).
"""
function check_webui_doc_links(; docs_src::AbstractString = normpath(joinpath(dirname(pathof(SparlectraApp)), "..", "..", "docs", "src")), io::IO = stdout)
  failures = String[]
  topics = SparlectraApp.WEBUI_HELP_TOPICS
  for topic in sort!(collect(keys(topics)))
    meta = topics[topic]
    hint = String(meta.hint)
    isempty(strip(hint)) && push!(failures, "$(topic): empty hint")
    length(hint) > 160 && push!(failures, "$(topic): hint has $(length(hint)) characters (limit 160)")
    occursin('`', hint) && push!(failures, "$(topic): hint contains a backtick")
    (occursin('—', hint) || occursin('–', hint)) && push!(failures, "$(topic): hint contains a dash")
    doc = String(meta.doc)
    isempty(doc) && continue
    m = match(r"^([A-Za-z0-9_]+)/#([A-Za-z0-9_.-]+)$", doc)
    if m === nothing
      push!(failures, "$(topic): doc '$(doc)' is not of the form page/#anchor")
      continue
    end
    page, anchor = String(m.captures[1]), String(m.captures[2])
    path = joinpath(docs_src, page * ".md")
    if !isfile(path)
      push!(failures, "$(topic): doc page $(page).md does not exist")
      continue
    end
    occursin("(@id $(anchor))", read(path, String)) || push!(failures, "$(topic): no heading with (@id $(anchor)) in $(page).md")
  end
  println(io, "Web UI help registry: ", length(topics), " topic(s), ", length(failures), " failure(s)")
  for f in failures
    println(io, "  FAIL ", f)
  end
  return failures
end

if abspath(PROGRAM_FILE) == @__FILE__
  isempty(Base.invokelatest(check_webui_doc_links)) || exit(1)
end
