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

# file: src/hvdc_link.jl
# purpose: persistent HVDC link record carried on Net (net.hvdcLinks) so the
#          result layer knows every link regardless of how it entered the
#          model: MATPOWER dcline rows (Stage-0 or paired_control), CGMES
#          DC-topology pairs (injections or paired_control), or the
#          programmatic API. Defined before network.jl so the Net field can
#          be concretely typed; the controller (hvdc_pair_control.jl, loaded
#          later) registers itself by name via controller_name.

"""
    HvdcLink

Immutable record of one HVDC link between two converter prosumers. The
record identifies the terminals and their provenance; live electrical
values are always read from the prosumers or the attached controller, never
stored here. `controller_name` is `nothing` for Stage-0 fixed injections
and carries the `HvdcPairControl` name once a pair controller is attached
(updates replace the vector element, the record itself stays immutable).

Fields: `name`, `from_bus`/`to_bus` (bus indices), `from_prosumer`/
`to_prosumer` (indices into `net.prosumpsVec`), `status` (1 = in service),
`source` (`:matpower`, `:cgmes`, `:api`), `kind` (`:b2b` back-to-back or
`:p2p` point-to-point, CGMES: a participating `DCLineSegment` makes it
`:p2p`), `controller_name`.
"""
struct HvdcLink
  name::String
  from_bus::Int
  to_bus::Int
  from_prosumer::Int
  to_prosumer::Int
  status::Int
  source::Symbol
  kind::Symbol
  controller_name::Union{Nothing,String}
end

# replace-style update helpers: the record is immutable, so attaching or
# detaching a controller swaps the element in net.hvdcLinks
_hvdc_link_with_controller(l::HvdcLink, ctrl_name::Union{Nothing,String}) =
  HvdcLink(l.name, l.from_bus, l.to_bus, l.from_prosumer, l.to_prosumer, l.status, l.source, l.kind, ctrl_name)
