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

# file: src/adapters/powsybl/converter.jl
# purpose: the PowSyBl adapter: IIDM networks arrive as a table bundle
#          that pypowsybl produced (tools/powsybl_dump.py), or through the
#          optional extension straight from an .xiidm file. Sparlectra
#          never parses IIDM itself. The importer builds the network
#          directly (build_net_from_powsybl); the converter captures that
#          network as an SCFCase on explicit request only.

"""
    PowsyblAdapterOptions

The adapter-scope options of the PowSyBl import, mirroring the
`powsybl_import` configuration scope: the system base, the HVDC model
(`:fixed_injection` or `:paired_control`), generator ids that override
the slack choice, whether every synchronous component gets its own
slack, how a remote voltage regulation is mapped (`:hold_local` keeps
the unit PV at its own bus, `:pq` makes it PQ with `target_q`, `:remote`
attaches an outer-loop machine voltage control on the regulated bus).
"""
Base.@kwdef struct PowsyblAdapterOptions
  base_mva::Float64 = 100.0
  hvdc_mode::Symbol = :fixed_injection
  slack_ids::Vector{String} = String[]
  multi_slack::Bool = true
  remote_regulation::Symbol = :hold_local
end

include("powsybl_tables.jl")
include("powsybl_csv.jl")
include("iidm_reader.jl")
include("powsybl_report.jl")
include("powsybl_mapping.jl")

struct PowsyblAdapter <: FormatAdapter end

"""
    powsybl_adapter_options(cfg::SparlectraConfig) -> PowsyblAdapterOptions

The adapter options an effective run configuration implies (the
`powsybl_import` scope).
"""
function powsybl_adapter_options(cfg::SparlectraConfig)::PowsyblAdapterOptions
  pw = cfg.powsybl
  return PowsyblAdapterOptions(base_mva = pw.base_mva, hvdc_mode = pw.hvdc_mode, slack_ids = copy(pw.slack_ids), multi_slack = pw.multi_slack, remote_regulation = pw.remote_regulation)
end

const _POWSYBL_FILE_EXTENSIONS = (".xiidm", ".xml")

# A PowSyBl source is a bundle directory (`<case>.powsybl` with a manifest
# of the bundle format) or an IIDM file (`.xiidm`, `.xiidm.bz2`, `.xml`).
_powsybl_is_bundle_dir(path::AbstractString) = isdir(path) && endswith(lowercase(String(path)), ".powsybl")

function _powsybl_is_iidm_file(path::AbstractString)::Bool
  isfile(path) || return false
  lowered = lowercase(String(path))
  (endswith(lowered, ".xiidm.bz2") || any(e -> endswith(lowered, e), _POWSYBL_FILE_EXTENSIONS)) || return false
  endswith(lowered, ".bz2") && return true   # compressed: the name is the evidence
  head = open(io -> String(read(io, 4096)), String(path), "r")
  return occursin("xmlns:iidm=", head) || occursin("<iidm:network", head)
end

"""
    detect(::Type{PowsyblAdapter}, path) -> Bool

A bundle directory whose name ends in `.powsybl` and whose `manifest.json`
carries `"format": "powsybl_tables"` in its first 4 KiB, or a regular
`.xiidm`, `.xiidm.bz2` or `.xml` file whose first 4 KiB carry the IIDM
namespace. Content based; the extension is a hint. Probed before the
CGMES detection, which claims every directory and every `.xml`.
"""
function detect(::Type{PowsyblAdapter}, path::AbstractString)::Bool
  if _powsybl_is_bundle_dir(path)
    manifest = joinpath(String(path), POWSYBL_MANIFEST_FILE)
    isfile(manifest) || return false
    head = open(io -> String(read(io, 4096)), manifest, "r")
    return occursin("\"format\"", head) && occursin("\"" * POWSYBL_BUNDLE_FORMAT * "\"", head)
  end
  return _powsybl_is_iidm_file(path)
end

options_type(::PowsyblAdapter) = PowsyblAdapterOptions

"""
    _powsybl_tables_of(path) -> PowsyblTables

The tables of a PowSyBl source: a `.powsybl` bundle directory reads
through `read_powsybl_bundle`, an IIDM file (`.xiidm`, `.xml`) through the
native reader `read_iidm_tables`, in Julia without Python.
"""
function _powsybl_tables_of(path::AbstractString)::PowsyblTables
  _powsybl_is_bundle_dir(path) && return read_powsybl_bundle(String(path))
  isdir(path) && throw(ArgumentError("PowSyBl source $(basename(String(path))) is a directory without the .powsybl suffix or a manifest.json of the bundle format"))
  return read_iidm_tables(String(path))
end

"""
    _import_powsybl(path, opts; name) -> (net, report, tables)

The one import path behind `import_net` and `import_case`: read the
tables of the source, build the network.
"""
function _import_powsybl(path::AbstractString, opts::PowsyblAdapterOptions; name::AbstractString = "")
  tables = _powsybl_tables_of(path)
  net, report = build_net_from_powsybl(tables, opts; name = name)
  return (net, report, tables)
end

"""
    import_net(::PowsyblAdapter, path, opts; config = active_sparlectra_config()) -> Net

The PowSyBl importer of the adapter contract: a `.powsybl` bundle
directory reads through `read_powsybl_bundle`, an IIDM file through the
native `read_iidm_tables`; both build through `build_net_from_powsybl`.
No Python is involved.
"""
function import_net(::PowsyblAdapter, path::AbstractString, opts::PowsyblAdapterOptions; config::SparlectraConfig = active_sparlectra_config())::Net
  net, _, _ = _import_powsybl(path, opts)
  return net
end

"""
    convert_case(::PowsyblAdapter, source, opts) -> SCFCase

Build the network and capture it as a typed case (`source_format =
"powsybl"`); the import report goes into the namespaced meta as
`powsybl_import`.
"""
function convert_case(::PowsyblAdapter, source::AbstractString, opts::PowsyblAdapterOptions)::SCFCase
  net, report, tables = _import_powsybl(source, opts)
  case = net_to_scfcase(net; source_format = "powsybl", source_reference = basename(String(source)), include_start_state = true)
  spar = case.sparlectra
  if spar !== nothing
    spar.meta["powsybl_import"] = Dict{String,Any}(
      "case" => String(get(tables.manifest, "case", "")),
      "pypowsybl_version" => String(get(tables.manifest, "pypowsybl_version", "")),
      "report" => format_powsybl_report(report),
      "skipped" => length(report.skipped),
    )
  end
  return case
end

convert_case(a::PowsyblAdapter, source::AbstractString) = convert_case(a, source, PowsyblAdapterOptions())
