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

# file: src/adapters/adapters.jl
# purpose: the format-adapter contract of the SCF pivot architecture
#          (adapter task stage 3, design decision D3): one adapter per
#          input format converts its source into the typed SCFCase, and
#          build_net is the one network constructor behind all of them.

"""
    FormatAdapter

Abstract supertype of the format adapters. The model, corrected by
task_import_direct (2026-09-04, superseding the D2/D3 pivot design): an
adapter owns BOTH output forms of its format. Its IMPORTER builds the
network directly and hands it to the solver, with no intermediate
format, no conversion, no loss; its CONVERTER turns the format into the
typed [`SCFCase`](@ref) and runs ONLY when someone asks for that file
(export to SCF, the Case page export action, the shipped-case builders,
or the scenario engine's id index). SCF is the preferred working format
but never a mandatory way station: no import path converts.

The contract per concrete adapter:

- `detect(::Type{A}, path)::Bool`: content-based detection reading at
  most the first 4 KiB of the file.
- `options_type(::A)::Type`: the adapter's option struct.
- `import_net(::A, source, opts)::Net`: the importer, used by
  `import_case`.
- `convert_case(::A, source, opts)::SCFCase`: the converter, explicit
  requests only.
- `export_case(::A, net, path; opts)`: optional, the reverse direction.
"""
abstract type FormatAdapter end

function detect end
function options_type end
function import_net end
function convert_case end
function export_case end

include("matpower/converter.jl")
include("dtf/converter.jl")
include("cgmes/converter.jl")
include("pgm/converter.jl")
