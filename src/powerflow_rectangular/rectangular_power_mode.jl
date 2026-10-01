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

# file: src/powerflow_rectangular/rectangular_power_mode.jl
# purpose: power mode (0.30.1, `power_flow.power_mode`): the state a repeated
#          solve on one network keeps between solves, so the warm solve
#          neither rebuilds the Ybus nor repeats the symbolic analysis of
#          the sparse LU nor allocates its work arrays again; plus the hook
#          through which the KLU package extension supplies its factorization.

"""
    PowerModeCache

What a power-mode solve keeps on the network (`net._power_cache`) for the
next solve: the Ybus with the fingerprint of the data it was built from,
the linear-solver context (symbolic analysis, factorization object, the
Jacobian assembly buffers), the Newton iteration workspace. The Ybus is
reused only when the fingerprint matches, the context re-analyses by
itself when the Jacobian pattern differs (bus types, bus count), so a
changed topology or a changed branch parameter never reuses stale state.
"""
mutable struct PowerModeCache
  n::Int
  ybus::Any
  ybus_fingerprint::UInt64
  linear_ctx::Any
  workspace::Any
  bus_type_fingerprint::UInt64
  ybus_reuse_count::Int
  context_reuse_count::Int
end

PowerModeCache() = PowerModeCache(0, nothing, UInt64(0), nothing, nothing, UInt64(0), 0, 0)

# The constructor of the power-mode linear-solver context: nothing means the
# UMFPACK reuse context; the KLU package extension (ext/SparlectraKLUExt.jl)
# sets it to its own context constructor when `using KLU` loads it.
const _POWER_MODE_LINEAR_CONTEXT = Ref{Any}(nothing)

"""
    power_mode_linear_solver_backend() -> Symbol

`:klu` when the KLU extension is loaded (`using KLU` in the session),
`:umfpack_reuse` otherwise: the sparse LU a power-mode solve uses.
"""
power_mode_linear_solver_backend() = _POWER_MODE_LINEAR_CONTEXT[] === nothing ? :umfpack_reuse : :klu

function _power_mode_new_context()
  ctor = _POWER_MODE_LINEAR_CONTEXT[]
  return ctor === nothing ? UmfpackReuseNewtonContext() : ctor()
end

# Fingerprint of everything createYBUS reads: the stamped admittances of
# every branch with its terminals and switching state, the shunt
# admittances, the isolated buses, the base. Cheaper than the assembly by
# two orders of magnitude (arithmetic per branch against sparse insertion),
# and a tap or a switching change moves it, so the cached Ybus is never
# stale.
function _ybus_fingerprint(net::Net)::UInt64
  h = hash(length(net.nodeVec), hash(net.baseMVA, UInt64(0x5b1ec7a0)))
  h = hash(net.isoNodes, h)
  for branch in net.branchVec
    h = hash((branch.fromBus, branch.toBus, branch.status, branch.from_status, branch.to_status), h)
    h = hash(calcAdmittance(branch, branch.comp.cVN, net.baseMVA), h)
  end
  for sh in net.shuntVec
    h = hash((sh.busIdx, sh.status, sh.model, sh.y_pu_shunt), h)
  end
  return h
end

function _power_cache!(net::Net)::PowerModeCache
  cache = net._power_cache
  if !(cache isa PowerModeCache)
    cache = PowerModeCache()
    net._power_cache = cache
  end
  return cache
end

"""
    reset_power_mode!(net)

Drops the power-mode state of `net` (Ybus, factorization, workspace). The
next power-mode solve rebuilds it. Not needed after a topology or
parameter change (the fingerprint and the pattern guard cover that); for
callers that want the memory back.
"""
function reset_power_mode!(net::Net)
  net._power_cache = nothing
  return nothing
end

# A copied network starts without power-mode state: the factorization
# objects wrap native memory (UMFPACK, KLU) that must not be shared or
# double-freed, and the Ybus fingerprint is recomputed on the first solve
# anyway. deepcopy(net) therefore yields an empty cache, not a copy.
Base.deepcopy_internal(::PowerModeCache, ::IdDict) = PowerModeCache()
