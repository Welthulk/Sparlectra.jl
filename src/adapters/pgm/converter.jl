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

# file: src/adapters/pgm/converter.jl
# purpose: the power-grid-model adapter (adapter task stage 3d): a plain
#          PGM input dataset is the SCF data section without the
#          namespaced block, so the conversion IS the shared parse,
#          validate and typing pipeline. The roles are inferred at build
#          time exactly as the reader documents (an in-service `source` is
#          the reference of its island, a regulated sym_gen is PV, loads
#          are PQ); the adapter produces no start_state and no namespaced
#          measurements, because a plain dataset carries neither.

struct PGMAdapter <: FormatAdapter end

"""
    PGMAdapterOptions

The power-grid-model conversion takes no options: the dataset is already
in the case format's own vocabulary. The struct exists for the adapter
contract.
"""
Base.@kwdef struct PGMAdapterOptions end

"""
    detect(::Type{PGMAdapter}, path) -> Bool

A `.json` file with the PGM dataset markers in its first 4 KiB. An SCF
case matches too, because the canonical writer sorts the namespaced
block to the end of the file, beyond any bounded probe: the registry
resolves that tie by probing SCF FIRST (D3, SCF detection keeps
precedence), and both formats run the same parse and build pipeline, so
the tie costs correctness nothing; the format label is refined at parse
time by the presence of the `sparlectra` key.
"""
detect(::Type{PGMAdapter}, path::AbstractString)::Bool = begin
  ext = lowercase(splitext(String(path))[2])
  ext == ".json" || return false
  isfile(path) || return false
  head = open(io -> String(read(io, 4096)), String(path), "r")
  occursin("\"sparlectra\"", head) && return false
  return occursin("\"is_batch\"", head) && occursin("\"version\"", head)
end

options_type(::PGMAdapter) = PGMAdapterOptions

"""
    import_net(::PGMAdapter, path::AbstractString, opts::PGMAdapterOptions; config = active_sparlectra_config()) -> Net

The PGM importer of the adapter contract (task_import_direct): a plain
power-grid-model dataset reads through the shared SCF reader (the SCF
data section IS a PGM dataset) and builds directly through build_net.
"""
import_net(::PGMAdapter, path::AbstractString, opts::PGMAdapterOptions; config::SparlectraConfig = active_sparlectra_config())::Net = build_net(read_scf_json(String(path)); config = config)

function convert_case(::PGMAdapter, source::AbstractString, opts::PGMAdapterOptions)::SCFCase
  return read_scf_json(String(source))
end

convert_case(a::PGMAdapter, source::AbstractString) = convert_case(a, source, PGMAdapterOptions())
