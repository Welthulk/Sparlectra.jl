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
# purpose: weak-ref-keyed registry for DC power-flow status NamedTuples with
#          the public reader dc_pf_status; kept separate from
#          rectangular_pf_status so DC results never masquerade as AC statuses

# The DC solver status lives directly on the Net (net._dc_pf_status) since
# the thread-safety Phase 1 rework, same treatment as the rectangular
# registry: the former global weak-ref table raced under concurrent solves.
# Kept as a separate field so DC results never masquerade as AC statuses.

function _set_dc_pf_status!(net::Net, status)
  # Plain field write; concurrency contract "one Net, one solver task at a
  # time" (see rectangular_status_workspace.jl).
  net._dc_pf_status = status
  return status
end

"""
    dc_pf_status(net::Net) -> Any

Retrieve the most recent DC power-flow status for a network (set by
[`rundcpf!`](@ref)), or `nothing` if `rundcpf!` has not run on this network.
Stored as a field on the `Net` (separate from `rectangular_pf_status` so a
DC result is never mistaken for an AC one); `deepcopy(net)` carries it,
serializers must skip it.
"""
dc_pf_status(net::Net) = net._dc_pf_status
