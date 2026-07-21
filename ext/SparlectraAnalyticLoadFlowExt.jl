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

# file: ext/SparlectraAnalyticLoadFlowExt.jl
# Package extension: bridges Sparlectra's external-solver interface
# (`PFModel`/`PFSolution`/`AbstractExternalSolver`/`solvePf`) to AnalyticLoadFlow.jl,
# an analytic power-series (holomorphic-embedding-style) load-flow solver.
#
# This module only loads when the host session has `using AnalyticLoadFlow` in
# addition to `using Sparlectra` (weak dependency via Julia package extensions).

module SparlectraAnalyticLoadFlowExt

using Sparlectra
using AnalyticLoadFlow

"""
    ApslfSolver <: AbstractExternalSolver

Adapter that runs a `PFModel` through AnalyticLoadFlow.jl's `solve_pf_apslf`.

Fields:
- `order::Int = 40`: highest power-series coefficient to compute.
- `use_pade::Bool = true`: evaluate the voltage series via Padé `[L/M]` approximants
  instead of direct Taylor summation.
- `nr_polish::Bool = true`: run a Newton-Raphson polishing step on the series result.
- `mode::Symbol = :direct`: `:direct` (native PV handling) or `:outer` (PQ-only series
  plus an outer secant loop for PV enforcement); forwarded to `solve_pf_apslf`.

Use [`apslf_solver`](@ref) to construct one without depending on this extension
module directly.
"""
Base.@kwdef struct ApslfSolver <: Sparlectra.AbstractExternalSolver
  order::Int = 40
  use_pade::Bool = true
  nr_polish::Bool = true
  mode::Symbol = :direct
end

"""
    solvePf(solver::ApslfSolver, model::PFModel; kwargs...) -> PFSolution

Solve `model` with AnalyticLoadFlow.jl's analytic power-series solver.

Maps the canonical `PFModel` fields onto the AnalyticLoadFlow spec expected by
`solve_pf_apslf`: `Ybus → Y`, `busType → bustype`, `real/imag(Sspec) → Pspec/Qspec`,
`Vset → Vm`, `qmin_pu/qmax_pu → Qmin/Qmax` (defaulting to `-Inf`/`Inf` when `model`
carries no Q-limits), `slack_idx → slack`. The returned voltage vector is in the same
PF ordering as `model.busIdx_net`.

`model.V0` is **not** used: AnalyticLoadFlow always starts from the canonical
analytic germ `V(s=0) = 1∠0` and does not accept an external start voltage.

`meta` carries solver-specific diagnostics: the series/Padé `order`, an APSLF
stability indicator (`dmin`/`pole`/`bus`/`level`, derived from the distance of
Padé poles to the physical evaluation point `s = 1`), and NR-polish bookkeeping
(`enabled`/`success`/`improved`/`rejected`/`reject_reason`).

Extra `kwargs` (e.g. `tol` forwarded by `runpf_external!`) are accepted and ignored;
convergence/polish behavior is controlled entirely via the `ApslfSolver` fields.
"""
function Sparlectra.solvePf(solver::ApslfSolver, model::Sparlectra.PFModel; kwargs...)
  n = length(model.busIdx_net)

  qmin_pu = isempty(model.qmin_pu) ? fill(-Inf, n) : model.qmin_pu
  qmax_pu = isempty(model.qmax_pu) ? fill(Inf, n) : model.qmax_pu

  spec = (
    Y = model.Ybus,
    bustype = model.busType,
    Pspec = real.(model.Sspec),
    Qspec = imag.(model.Sspec),
    Vm = model.Vset,
    Qmin = qmin_pu,
    Qmax = qmax_pu,
    slack = model.slack_idx,
  )

  res = AnalyticLoadFlow.solve_pf_apslf(
    spec;
    mode = solver.mode,
    order = solver.order,
    use_pade = solver.use_pade,
    nr_polish = solver.nr_polish,
    return_coeffs = true,
  )

  st = AnalyticLoadFlow.stability_from_Vcoeff(res.Vcoeff; slack = model.slack_idx, order = solver.order)
  stability = (dmin = st.dmin, pole = st.pole, bus = st.bus, level = AnalyticLoadFlow.st_level(st.dmin))

  meta = (
    solver = :apslf,
    mode = res.effective_mode,
    order = solver.order,
    use_pade = solver.use_pade,
    stability = stability,
    nr_polish_enabled = res.nr_polish_enabled,
    nr_polish_success = res.nr_polish_success,
    nr_polish_improved = res.nr_polish_improved,
    nr_polish_rejected = res.nr_polish_rejected,
    nr_polish_reject_reason = res.nr_polish_reject_reason,
    outer_iters = res.outer_iters,
    bustype_final = res.bustype,
  )

  return Sparlectra.PFSolution(
    V = res.V,
    converged = res.converged,
    iters = res.outer_iters,
    residual_inf = Sparlectra.mismatchInf(model, res.V),
    meta = meta,
  )
end

end # module SparlectraAnalyticLoadFlowExt
