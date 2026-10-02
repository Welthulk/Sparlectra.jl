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
  GRN/YEL/RED level). Reported as the Padé margin next to the Jacobian
  condition; costs about as much as the solve itself, so it can be switched
  off for large networks. The margin is not the convergence radius of the
  series (a Froissart doublet can sit next to `s = 1` on a convergent
  series), so the result carries the coefficient-growth radius of the
  series as well; that one is a root test on the coefficients and always
  evaluated (see `_apslf_series_radius`).

[`apslf_solver`](@ref) is the keyword constructor for it.
"""
Base.@kwdef struct ApslfSolver <: AbstractExternalSolver
  order::Int = 24
  use_pade::Bool = true
  nr_polish::Bool = false
  mode::Symbol = :direct
  convergence_radius::Bool = true
end

# Limits of the outer-mode fallback (PQ series plus a secant loop on the PV
# reactive powers) that AnalyticLoadFlow.solve_pf_apslf runs after a failed
# direct solve: the defaults of its keywords outer_fallback_nbus_max and
# outer_fallback_pv_bus_max in AnalyticLoadFlow 0.9.16. The single-pass call
# in _apslf_solve switches that fallback off and runs it itself, so it needs
# the same bounds to stay equivalent. AnalyticLoadFlow 0.9.16 provides no
# binding for them (literal keyword defaults in its solver_core.jl), so they
# are mirrored here; test_apslf.jl ("AnalyticLoadFlow version guard") names
# the version and checks them against the loaded package's source. Once the
# package provides them, read them from there instead.
const _APSLF_OUTER_FALLBACK_NBUS_MAX = 2000
const _APSLF_OUTER_FALLBACK_PV_BUS_MAX = 250

"""
    APSLF_MIN_VERSION

Lowest AnalyticLoadFlow.jl version the adapter is verified against. The
Project.toml compat carries the same bound, but Julia does not check a
Manifest against the compat at load time: an environment whose Manifest
predates the bump loads the older version without a message, and an older
version can return a wrong series result flagged as converged.
`check_apslf_version` closes that gap at load time.
"""
const APSLF_MIN_VERSION = v"0.9.16"

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

`meta` carries solver-specific diagnostics: the series/Padé `order`, the
Padé margin (`stability`: `dmin`/`pole`/`bus`/`level`, derived from the
distance of Padé poles to the physical evaluation point `s = 1`), the
coefficient-growth radius of the series (`series_radius`: `radius`/`bus`, see
`_apslf_series_radius`), and NR-polish bookkeeping
(`enabled`/`success`/`improved`/`rejected`/`reject_reason`). Both bus fields
are in the numbering of the net the model was built from (`busIdx_net`).

Without a finite reactive limit in `model` the direct kernel runs one outer
pass instead of up to 30 identical ones (see `_apslf_solve`).

Extra `kwargs` (e.g. `tol` forwarded by `runpf_external!`) are accepted and ignored;
convergence/polish behavior is controlled entirely via the `ApslfSolver` fields.
"""
function solvePf(solver::ApslfSolver, model::PFModel; kwargs...)
  spec = _apslf_spec_from_model(model)

  res = _apslf_solve(solver, spec)

  # Padé-pole margin, optional because its cost is comparable to the solve;
  # the bus index is reported in the numbering of the net the model was
  # built from (busIdx_net), never in PF ordering. Mapping it to the case
  # bus number is the caller's job (_apslf_radius_status), because an
  # island net renumbers its buses.
  stability = if solver.convergence_radius
    st = AnalyticLoadFlow.stability_from_Vcoeff(res.Vcoeff; slack = model.slack_idx, order = solver.order)
    (dmin = st.dmin, pole = st.pole, bus = st.bus >= 1 ? Int(model.busIdx_net[st.bus]) : 0, level = AnalyticLoadFlow.st_level(st.dmin), enabled = true)
  else
    (dmin = NaN, pole = NaN + NaN * im, bus = 0, level = "off", enabled = false)
  end
  # coefficient-growth radius of the series, same bus numbering as above
  sr = _apslf_series_radius(res.Vcoeff, model.slack_idx)
  series_radius = (radius = sr.radius, bus = sr.idx >= 1 ? Int(model.busIdx_net[sr.idx]) : 0)

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
    series_radius = series_radius,
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

# One call of AnalyticLoadFlow for `spec`. Without a finite reactive limit
# (power_flow.qlimits off, or no machine with a limit) nothing changes
# between the outer passes of the direct kernel: no machine can switch, and
# every pass recomputes the same series from the same germ. A diverging run
# used to repeat that identical solve max_outer = 30 times (case13659pegase:
# 30 times 1.2 s, the same voltages to the last bit with max_outer = 1), so
# such a run gets one pass. The outer-mode fallback that solve_pf_apslf runs
# after a failed direct solve is different: its passes refine the PV
# reactive powers by a secant step and need them. It is therefore switched
# off in the single-pass call and, within the same size bounds, run here
# with the default number of passes, which keeps the result of the former
# single call.
# INTERIM, remove when AnalyticLoadFlow releases its own phase-shifter
# support and the compat bound moves past it (Sparlectra issue #470, "Remove the
# interim phase-shifter handling for APSLF once AnalyticLoadFlow supports
# phase shifters"). A network with phase-shifting transformers has an
# unsymmetric Y-bus (a complex tap ratio; real off-nominal taps keep Y
# symmetric). The registered AnalyticLoadFlow 0.9.16 models a fixed shift
# exactly with the :deviation and :noload germs (sp_casePST at -5.73, -20
# and +15 degrees agrees with Newton to 1e-14 pu) and not with :flat (not
# converged, 0.09 to 0.31 pu off), so such a run pins :deviation instead of
# relying on the package default staying :deviation.
_apslf_has_phase_shift(Y) = maximum(abs, Y - transpose(Y); init = 0.0) > 1.0e-12 * max(1.0, maximum(abs, Y; init = 0.0))
_apslf_germ_kwargs(spec) = _apslf_has_phase_shift(spec.Y) ? (; germ = :deviation) : (;)

function _apslf_solve(solver::ApslfSolver, spec)
  common = (mode = solver.mode, order = solver.order, use_pade = solver.use_pade, nr_polish = solver.nr_polish, return_coeffs = true, _apslf_germ_kwargs(spec)...)
  limits_free = !any(isfinite, spec.Qmin) && !any(isfinite, spec.Qmax)
  (solver.mode === :direct && limits_free) || return AnalyticLoadFlow.solve_pf_apslf(spec; common...)
  res = AnalyticLoadFlow.solve_pf_apslf(spec; common..., max_outer = 1, outer_fallback_nbus_max = 0)
  res.converged && return res
  nbus = length(spec.bustype)
  npv = count(==(:PV), spec.bustype)
  (nbus <= _APSLF_OUTER_FALLBACK_NBUS_MAX && npv <= _APSLF_OUTER_FALLBACK_PV_BUS_MAX) || return res
  res_outer = AnalyticLoadFlow.solve_pf_apslf(spec; merge(common, (mode = :outer,))...)
  # like solve_pf_apslf itself: the fallback result only when it converged
  return res_outer.converged ? res_outer : res
end

"""
    _apslf_series_radius(Vcoeff, slack) -> (radius, idx)

Coefficient-growth (root test) estimate of the radius of convergence of the
APSLF voltage series in the embedding parameter `s`: for every non-slack row
`i` of `Vcoeff` (PF ordering, column `k + 1` holds the coefficient of
`s^k`), the smallest `|c_k|^(-1/k)` over the last four orders, then the
smallest over the rows; `idx` is the row that sets it. A radius below 1 means
the plain series does not reach the physical point `s = 1` (the Padé
evaluation may still continue past it); this is the number that predicted the
outcome on the grid-bench data sets where the Padé margin did not. `Inf`
when every coefficient beyond order 0 vanishes, `NaN` (with `idx = 0`) for a
series of order 0.
"""
function _apslf_series_radius(Vcoeff::AbstractMatrix, slack::Int)
  kmax = size(Vcoeff, 2) - 1
  kmax >= 1 || return (radius = NaN, idx = 0)
  radius = Inf
  idx = 0
  for i in axes(Vcoeff, 1)
    i == slack && continue
    for k in max(1, kmax - 3):kmax
      a = abs(Vcoeff[i, k + 1])
      (a > 0.0 && isfinite(a)) || continue
      r = a^(-1 / k)
      if r < radius
        radius = r
        idx = i
      end
    end
  end
  return (radius = radius, idx = idx)
end

# AnalyticLoadFlow reports bus types in lower case (:slack/:pv/:pq); the
# PFModel uses :Slack/:PV/:PQ
function _apslf_bus_type(s::Symbol)::Symbol
  t = lowercase(String(s))
  t == "pv" && return :PV
  t == "pq" && return :PQ
  return :Slack
end

