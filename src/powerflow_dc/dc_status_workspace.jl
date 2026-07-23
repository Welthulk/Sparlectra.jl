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
#
# This file is included inside module Sparlectra. Do not add a module wrapper here.
#
# DC power-flow status registry.
#
# Deliberately a separate weak-ref-keyed table from
# _RectangularPFStatusTable/rectangular_pf_status, not a shared one: a DC
# status NamedTuple has no nr_converged/wrong_branch_status/... AC-only
# fields, so keeping it out of rectangular_pf_status(net) prevents it from
# silently masquerading as an AC status to AC-only readers
# (_rectangular_status_diagnostics, _rectangular_run_status, ...).

# file: src/powerflow_dc/dc_status_workspace.jl

mutable struct _DcPFStatusTable
  entries::Vector{Tuple{UInt,WeakRef,Any}}
end

const _DC_PF_STATUS = _DcPFStatusTable(Tuple{UInt,WeakRef,Any}[])

function _prune_dc_pf_status!(table::_DcPFStatusTable)
  filter!(entry -> entry[2].value !== nothing, table.entries)
  return table
end

function _set_dc_pf_status!(net::Net, status)
  _prune_dc_pf_status!(_DC_PF_STATUS)
  key = objectid(net)
  for i in eachindex(_DC_PF_STATUS.entries)
    entry = _DC_PF_STATUS.entries[i]
    if entry[1] == key && entry[2].value === net
      _DC_PF_STATUS.entries[i] = (key, WeakRef(net), status)
      return status
    end
  end
  push!(_DC_PF_STATUS.entries, (key, WeakRef(net), status))
  return status
end

"""
    dc_pf_status(net::Net) -> Any

Retrieve the most recent DC power-flow status for a network (set by
[`rundcpf!`](@ref)), or `nothing` if `rundcpf!` has not run on this network
(or it has been garbage-collected). Kept in a registry separate from
`rectangular_pf_status` so a DC result is never mistaken for an AC one by
AC-only status readers.
"""
function dc_pf_status(net::Net)
  _prune_dc_pf_status!(_DC_PF_STATUS)
  key = objectid(net)
  for entry in _DC_PF_STATUS.entries
    entry[1] == key && entry[2].value === net && return entry[3]
  end
  return nothing
end
