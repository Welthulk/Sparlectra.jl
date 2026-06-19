# Copyright 2023–2026 Udo Schmitz
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

function _artifact_kind(path::AbstractString)::Symbol
  name = lowercase(basename(path))
  name == "effective_config.yaml" && return :effective_config
  name == "run_metadata.yaml" && return :run_metadata
  name == "matpower_auto_profile.log" && return :matpower_auto_profile
  name == "result.json" && return :result_json
  name == "q_limit.log" && return :q_limit_log
  ext = lowercase(splitext(name)[2])
  ext == ".log" && return :log
  ext == ".csv" && return :csv
  ext in (".yaml", ".yml") && return :config
  ext in (".txt", ".md", ".html", ".htm", ".pdf") && return :report
  return :other
end

function _artifact_mime_type(path::AbstractString)::String
  ext = lowercase(splitext(path)[2])
  ext == ".log" && return "text/plain"
  ext == ".json" && return "application/json"
  ext == ".csv" && return "text/csv"
  ext in (".yaml", ".yml") && return "application/x-yaml"
  ext in (".txt", ".md") && return "text/plain"
  ext in (".html", ".htm") && return "text/html"
  ext == ".pdf" && return "application/pdf"
  return "application/octet-stream"
end

function _artifact_description(kind::Symbol, name::String)::String
  kind === :log && return "Power-flow execution log"
  kind === :q_limit_log && return "Q-limit diagnostic log"
  kind === :matpower_auto_profile && return "MATPOWER import auto-profile diagnostic log"
  kind === :result_json && return "Serialized Sparlectra API result"
  kind === :csv && return "Generated CSV output"
  kind === :report && return "Generated report"
  kind === :config && return "Configuration artifact"
  kind === :effective_config && return "Effective configuration used for this run"
  kind === :run_metadata && return "Runtime request and lifecycle metadata for this run"
  return "Generated Sparlectra artifact $(name)"
end

"""
    collect_sparlectra_api_artifacts(output_dir) -> Vector{SparlectraApiArtifact}

Discover files beneath an API run directory and return deterministic artifact
metadata. This is the only filename-discovery contract needed by GUI clients.
"""
function collect_sparlectra_api_artifacts(output_dir::AbstractString)::Vector{SparlectraApiArtifact}
  root = abspath(output_dir)
  isdir(root) || return SparlectraApiArtifact[]
  files = String[]
  for (directory, _, names) in walkdir(root)
    for name in names
      push!(files, joinpath(directory, name))
    end
  end
  sort!(files)
  return [begin
    kind = _artifact_kind(path)
    relative_name = relpath(path, root)
    exists = isfile(path)
    SparlectraApiArtifact(relative_name, kind, abspath(path), _artifact_mime_type(path), exists, exists ? Int(filesize(path)) : nothing, _artifact_description(kind, relative_name))
  end for path in files]
end
