"""
    SparlectraApiArtifact

Metadata for one file produced or discovered by [`run_sparlectra_api`](@ref).
Artifact paths are absolute so callers do not need to know Sparlectra's internal
file naming or working-directory conventions.
"""
struct SparlectraApiArtifact
  name::String
  kind::Symbol
  path::String
  mime_type::String
  exists::Bool
  size_bytes::Union{Int,Nothing}
  description::String
end

"""
    SparlectraApiResult

Stable, non-interactive result returned by [`run_sparlectra_api`](@ref). Each
result has a unique `run_id` and a `schema_version` for transport consumers. The
`raw_result` field contains the underlying [`SparlectraRunResult`](@ref) for a
completed solver invocation and is `nothing` for input or configuration errors.
Use the serialization helpers for transport-safe representations that omit
`raw_result` by default.
"""
struct SparlectraApiResult
  run_id::String
  schema_version::String
  status::Symbol
  success::Bool
  converged::Union{Bool,Nothing}
  solution_available::Bool
  iterations::Union{Int,Nothing}
  final_mismatch::Union{Float64,Nothing}
  reason::Union{String,Nothing}
  message::Union{String,Nothing}
  casefile::Union{String,Nothing}
  config_file::Union{String,Nothing}
  output_dir::String
  logfile::Union{String,Nothing}
  result_file::Union{String,Nothing}
  artifacts::Vector{SparlectraApiArtifact}
  service_phase_timings::Vector{Dict{String,Any}}
  raw_result::Any
end
