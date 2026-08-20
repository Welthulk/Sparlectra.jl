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

# Date: 2026-07-26
# file: src/acpflow/net_cache.jl
# purpose: opt-in binary cache of the parsed MATPOWER case and the built Net, keyed by file hash + import-option fingerprint (issue #292)

# Cache design:
# - Opt-in via matpower_import.net_cache.enabled (default false) and only
#   active when matpower_import.auto_profile == off (auto-profile derives
#   config rewrites from the parsed case; keeping the gate here makes hit and
#   miss runs behave identically).
# - v1 caches the PARSED CASE (mpc), not the built Net: measured on an
#   82k-bus case, deserializing the Net object graph (~82k Node structs plus
#   name Dicts) costs more than parsing + building from scratch (1.8 s vs
#   1.3 s), while the mpc matrices deserialize in a fraction of the parse
#   time. The network is always built fresh from the cached mpc, so import
#   options need not be part of the key. The format version is reserved for
#   a future flat Net representation.
# - Key: SHA-256 over the component list below — cache format version,
#   Sparlectra version, Julia VERSION, SHA-256 of the case-file bytes.
# - Storage: <case dir>/.sparlectra_net_cache/<case base>.<key16>.v1.jls,
#   written atomically (temp file + mv).
# - Every failure path (unreadable file, deserialize error, component
#   mismatch, read-only directory) silently falls back to a fresh parse.

using Serialization

const NET_CACHE_FORMAT_VERSION = 1
const _NET_CACHE_MARKER = "sparlectra_net_cache"

function _net_cache_components(filename::AbstractString, cfg::SparlectraConfig)::Vector{Pair{String,String}}
  return [
    "format" => string(NET_CACHE_FORMAT_VERSION),
    "sparlectra" => string(pkgversion(Sparlectra)),
    "julia" => string(VERSION),
    "case_sha256" => bytes2hex(sha256(read(filename))),
  ]
end

function _net_cache_key(components::Vector{Pair{String,String}})::String
  io = IOBuffer()
  for (k, v) in components
    print(io, k, '=', v, '\n')
  end
  return bytes2hex(sha256(take!(io)))[1:16]
end

function _net_cache_path(filename::AbstractString, key::AbstractString)::String
  base = splitext(basename(filename))[1]
  return joinpath(dirname(abspath(filename)), ".sparlectra_net_cache", string(base, ".", key, ".v", NET_CACHE_FORMAT_VERSION, ".jls"))
end

"""
    _net_cache_load(path, components) -> Union{Nothing,MatpowerIO.MatpowerCase}

Returns the cached parsed case or `nothing`. Beyond the key in the file
name, the stored component list must match exactly — a hash-prefix collision
or stale format silently misses instead of producing a wrong case.
"""
function _net_cache_load(path::AbstractString, components::Vector{Pair{String,String}})
  isfile(path) || return nothing
  try
    payload = open(deserialize, path, "r")
    payload isa Tuple{String,Vector{Pair{String,String}},MatpowerIO.MatpowerCase} || return nothing
    marker, stored_components, mpc = payload
    marker == _NET_CACHE_MARKER || return nothing
    stored_components == components || return nothing
    return mpc
  catch err
    err isa InterruptException && rethrow(err)
    @debug "net cache read failed; falling back to a fresh parse" path error = err
    return nothing
  end
end

function _net_cache_store(path::AbstractString, components::Vector{Pair{String,String}}, mpc::MatpowerIO.MatpowerCase)
  try
    mkpath(dirname(path))
    # pid AND task identity: two tasks of one process importing the same
    # case concurrently must not write the same temp file (Phase 1)
    tmp = string(path, ".", getpid(), ".", objectid(current_task()), ".tmp")
    open(tmp, "w") do io
      serialize(io, (_NET_CACHE_MARKER, components, mpc))
    end
    mv(tmp, path; force = true)
  catch err
    err isa InterruptException && rethrow(err)
    @debug "net cache write failed; continuing without cache" path error = err
  end
  return nothing
end
