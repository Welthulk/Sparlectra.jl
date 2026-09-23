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
#
# file: src/acpflow/apslf_solver.jl
# purpose: bridges the external-solver interface (PFModel/PFSolution/
#          AbstractExternalSolver/solvePf) to AnalyticLoadFlow.jl, an analytic
#          power-series (holomorphic-embedding-style) load flow. This was a
#          a package extension while AnalyticLoadFlow.jl was a weak dependency;
#          it is a normal dependency now, so the solver is always available
#          and no session has to load anything extra.

"""
    ApslfSolver <: AbstractExternalSolver

Adapter that runs a `PFModel` through AnalyticLoadFlow.jl's `solve_pf_apslf`.

Fields:
- `order::Int = 24`: highest power-series coefficient to compute.
- `use_pade::Bool = true`: evaluate the voltage series via Padé `[L/M]` approximants
  instead of direct Taylor summation.
- `nr_polish::Bool = false`: run a Newton-Raphson polishing step on the series result.
  Off since 0.13.0 (AnalyticLoadFlow 0.9.15): the series alone is a load-flow
  solution, the polish is a debugging aid.
- `mode::Symbol = :direct`: `:direct` (native PV handling) or `:outer` (PQ-only series
  plus an outer secant loop for PV enforcement); forwarded to `solve_pf_apslf`.
- `convergence_radius::Bool = true`: evaluate AnalyticLoadFlow's Padé-pole
  margin (`stability_from_Vcoeff`: the distance `dmin` of the nearest Padé
  pole to the evaluation point `s = 1`, with the bus that owns it and the
  GRN/YEL/RED level). Reported as the APSLF convergence radius next to the
  Jacobian condition; costs about as much as the solve itself, so it can be
  switched off for large networks.

[`apslf_solver`](@ref) is the keyword constructor for it.
"""
Base.@kwdef struct ApslfSolver <: AbstractExternalSolver
  order::Int = 24
  use_pade::Bool = true
  nr_polish::Bool = false
  mode::Symbol = :direct
  convergence_radius::Bool = true
end

"""
    APSLF_MIN_VERSION

Lowest AnalyticLoadFlow.jl version the adapter is verified against. The
Project.toml compat carries the same bound, but Julia does not check a
Manifest against the compat at load time: an environment whose Manifest
predates the bump loads the older version without a message, and an older
version can return a wrong series result flagged as converged.
`check_apslf_version` closes that gap at load time.
"""
const APSLF_MIN_VERSION = v"0.9.15"

"""
    check_apslf_version(loaded = pkgversion(AnalyticLoadFlow)) -> VersionNumber

Return `loaded` when it satisfies `APSLF_MIN_VERSION`, otherwise
throw an `ErrorException` that names both versions and the Pkg command that
fixes the environment. Called from the module `__init__`, so an outdated
environment fails `using Sparlectra` with the cause in one line instead of
producing wrong voltages later. Tests call it with an explicit version.
"""
function check_apslf_version(loaded::VersionNumber = pkgversion(AnalyticLoadFlow))::VersionNumber
  loaded >= APSLF_MIN_VERSION && return loaded
  error(
    "AnalyticLoadFlow $(loaded) is loaded, Sparlectra $(SparlectraVersion) needs at least $(APSLF_MIN_VERSION). ",
    "The Manifest.toml of this environment predates the dependency bump. Run in the checkout directory\n",
    "    julia --project=. -e 'using Pkg; Pkg.update(\"AnalyticLoadFlow\")'\n",
    "and start Julia again (a sysimage built on the old version has to be rebuilt as well).",
  )
end

"""
    _apslf_spec_from_model(model::PFModel) -> NamedTuple

Pure mapping from the canonical `PFModel` fields onto the AnalyticLoadFlow
spec expected by `solve_pf_apslf`: `Ybus → Y`, `busType → bustype`,
`real/imag(Sspec) → Pspec/Qspec`, `Vset → Vm`, `qmin_pu/qmax_pu → Qmin/Qmax`
(defaulting to `-Inf`/`Inf` when `model` carries no Q-limits), `slack_idx →
slack`. Split out from `solvePf` so the mapping itself is independently
testable.
"""
function _apslf_spec_from_model(model::PFModel)
  n = length(model.busIdx_net)

  qmin_pu = isempty(model.qmin_pu) ? fill(-Inf, n) : model.qmin_pu
  qmax_pu = isempty(model.qmax_pu) ? fill(Inf, n) : model.qmax_pu

  return (
    Y = model.Ybus,
    bustype = model.busType,
    Pspec = real.(model.Sspec),
    Qspec = imag.(model.Sspec),
    Vm = model.Vset,
    Qmin = qmin_pu,
    Qmax = qmax_pu,
    slack = model.slack_idx,
  )
end

"""
    solvePf(solver::ApslfSolver, model::PFModel; kwargs...) -> PFSolution

Solve `model` with AnalyticLoadFlow.jl's analytic power-series solver.

Maps the canonical `PFModel` fields onto the AnalyticLoadFlow spec via
[`_apslf_spec_from_model`](@ref). The returned voltage vector is in the same
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
function solvePf(solver::ApslfSolver, model::PFModel; kwargs...)
  spec = _apslf_spec_from_model(model)

  res = AnalyticLoadFlow.solve_pf_apslf(
    spec;
    mode = solver.mode,
    order = solver.order,
    use_pade = solver.use_pade,
    nr_polish = solver.nr_polish,
    return_coeffs = true,
  )

  # Padé-pole margin (the APSLF convergence radius), optional because its
  # cost is comparable to the solve; the bus index is reported in the net's
  # numbering, never in PF ordering
  stability = if solver.convergence_radius
    st = AnalyticLoadFlow.stability_from_Vcoeff(res.Vcoeff; slack = model.slack_idx, order = solver.order)
    (dmin = st.dmin, pole = st.pole, bus = st.bus >= 1 ? Int(model.busIdx_net[st.bus]) : 0, level = AnalyticLoadFlow.st_level(st.dmin), enabled = true)
  else
    (dmin = NaN, pole = NaN + NaN * im, bus = 0, level = "off", enabled = false)
  end

  # The residual is judged against the bus types AnalyticLoadFlow ended
  # with, not the ones it started from: a PV bus that hit a reactive limit
  # is a PQ bus at that limit in the solution, and its voltage equation no
  # longer holds by design. Judging it against Vset reported 0.027 pu on
  # case118 (19 clamped machines) for a solve that met every equation of
  # the final active set to 1e-13, and the run was labelled not converged.
  bt_final = Symbol[_apslf_bus_type(s) for s in res.bustype]
  S_final = ComplexF64[complex(real(model.Sspec[i]), res.Q[i]) for i in eachindex(model.Sspec)]
  F = mismatch_rectangular(model.Ybus, res.V, S_final, bt_final, model.Vset, model.slack_idx)
  residual_final = maximum(abs.(F))
  # clamped machines, in net bus numbering, with the limit side that binds
  clamps = Tuple{Int,Symbol}[]
  for i in eachindex(bt_final)
    (model.busType[i] === :PV && bt_final[i] === :PQ) || continue
    side = abs(res.Q[i] - spec.Qmax[i]) <= abs(res.Q[i] - spec.Qmin[i]) ? :max : :min
    push!(clamps, (Int(model.busIdx_net[i]), side))
  end

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
    qlimit_clamps = clamps,
  )

  return PFSolution(
    V = res.V,
    converged = res.converged,
    iters = res.outer_iters,
    residual_inf = residual_final,
    meta = meta,
  )
end

# AnalyticLoadFlow reports bus types in lower case (:slack/:pv/:pq); the
# PFModel uses :Slack/:PV/:PQ
function _apslf_bus_type(s::Symbol)::Symbol
  t = lowercase(String(s))
  t == "pv" && return :PV
  t == "pq" && return :PQ
  return :Slack
end

