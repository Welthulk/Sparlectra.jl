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

# file: src/acpflow/auto_powerflow.jl
# purpose: automatic power-flow mode (power_flow.mode = auto): read-only
#          network feature extraction, ordered-rule strategy selection,
#          precedence-safe option application (explicit user keys always
#          win), the bounded escalation ladder on non-convergence, and the
#          advisory Q-limit hint generator. The final result is always a
#          rectangular NR solution (APSLF only seeds start values, the DC
#          fallback is honestly labeled and never reported as AC).

## ---------------------------------------------------------------------------
## Thresholds: every numeric decision constant lives HERE, nowhere else.
## Chosen conservatively against the known MATPOWER matrix (case5..case300,
## 1354/9241pegase, SyntheticUSA experience); a calibration sweep over the
## benchmark cases is a possible follow-up (task question Q3).
## ---------------------------------------------------------------------------
const AUTO_PF_THRESHOLDS = (
  small_bus_max = 60,            # <= this many buses: small_default territory
  large_bus_min = 2000,          # >= this many buses: hard_case territory
  resistive_share = 0.5,         # share of branches with R > X above this: resistive network
  resistive_rx_median = 1.0,     # median R/X above this: resistive network
  low_rx_share = 0.15,           # share_r_gt_x below this: transmission-style, DC angles trustworthy
  many_narrow_q = 3,             # zero+narrow Q ranges at or above this: hard_case trigger
  profile_vm_lo = 0.5,           # plausible start-profile magnitude band (wide guard band, pu)
  profile_vm_hi = 1.5,
  qlimit_hard_start_iter = 5,    # hard_case: hold Q-limit switching back this long
  qlimit_hard_hysteresis = 0.02, # hard_case: raised from the 0.01 default
  qlimit_hard_cooldown = 2,      # hard_case: raised from the 1 default
  escalation_max_stages = 7,     # hard cap: base run plus at most this many escalations
)

## ---------------------------------------------------------------------------
## Task 1: network feature extraction (read-only, format independent, cheap)
## ---------------------------------------------------------------------------

"""
    AutoPfFeatures

Read-only summary of an imported network, the input of
[`select_auto_pf_strategy`](@ref). Produced by
[`collect_auto_pf_features`](@ref) in one pass over buses, branches, and
prosumers; no factorization, no solve, no mutation.
"""
struct AutoPfFeatures
  n_bus::Int
  n_branch::Int
  n_ac_islands::Int
  rx_median::Float64
  rx_p90::Float64
  share_r_gt_x::Float64
  n_phase_shifters::Int
  n_pv::Int
  pv_share::Float64
  n_q_zero_range::Int
  n_q_narrow_range::Int
  n_gen_with_q_limits::Int
  share_pv_with_limits::Float64
  has_start_profile::Bool
  start_profile_plausible::Bool
  n_voltage_levels::Int
end

"""
    collect_auto_pf_features(net::Net; min_q_range_pu = 1e-4) -> AutoPfFeatures

Extract the auto-mode decision features from an imported network. Pure
read-only and format independent: format specifics end at import, this
works the same on MATPOWER, CGMES, and DTF nets. `min_q_range_pu` is the
narrow-range threshold the Q-limit guard uses
(`power_flow.qlimits.guard_min_q_range_pu`).
"""
function collect_auto_pf_features(net::Net; min_q_range_pu::Float64 = 1e-4)
  nbus = length(net.nodeVec)
  nbranch = length(net.branchVec)
  # islands on the electrical view (closed links connect); read-only
  # No try/catch: electricalIslandComponents handles a fresh, unsolved net
  # on its own (measured 2026-09-06). A catch could only absorb a
  # programming error and would then feed "one island" into the feature
  # vector, which is how compareWithSV's island alignment stayed dead.
  n_islands = length(electricalIslandComponents(net))
  # branch R/X profile, guarded against X = 0 (skipped, not divided)
  rx = Float64[]
  n_r_gt_x = 0
  n_rx = 0
  n_shift = 0
  for br in net.branchVec
    br.status == 1 || continue
    if br.ratio != 0.0 && (br.angle != 0.0 || br.phase_shift_deg != 0.0)
      n_shift += 1
    end
    x = abs(br.x_pu)
    r = abs(br.r_pu)
    n_rx += 1
    r > x && (n_r_gt_x += 1)
    x > 1e-12 && push!(rx, r / x)
  end
  sort!(rx)
  rx_median = isempty(rx) ? 0.0 : rx[cld(length(rx), 2)]
  rx_p90 = isempty(rx) ? 0.0 : rx[max(1, ceil(Int, 0.9 * length(rx)))]
  share_r_gt_x = n_rx == 0 ? 0.0 : n_r_gt_x / n_rx
  # bus loop: PV count, start profile, voltage levels
  n_pv = 0
  n_profile = 0
  vm_min = Inf
  vm_max = -Inf
  vm_first = NaN
  all_identical = true
  levels = Set{Float64}()
  for nd in net.nodeVec
    getNodeType(nd) == PV && (n_pv += 1)
    push!(levels, nd.comp.cVN)
    vm = nd._vm_pu
    vm === nothing && continue
    isfinite(vm) || continue
    n_profile += 1
    vm < vm_min && (vm_min = vm)
    vm > vm_max && (vm_max = vm)
    isnan(vm_first) ? (vm_first = vm) : (all_identical &= vm == vm_first)
  end
  has_profile = n_profile == nbus && nbus > 0
  # plausibility: finite everywhere, inside a wide guard band, and not one
  # constant value (a flat 1.0 pu profile carries no information)
  plausible = has_profile && vm_min >= AUTO_PF_THRESHOLDS.profile_vm_lo && vm_max <= AUTO_PF_THRESHOLDS.profile_vm_hi && !all_identical
  # generator Q limits, counted per generator (the hints name generators)
  n_gen_lim = 0
  n_zero = 0
  n_narrow = 0
  n_gen_pv = 0
  for ps in net.prosumpsVec
    isGenerator(ps) || continue
    n_gen_pv += 1
    (ps.maxQ !== nothing && ps.minQ !== nothing && isfinite(ps.maxQ) && isfinite(ps.minQ)) || continue
    n_gen_lim += 1
    range = ps.maxQ - ps.minQ
    if range <= 0.0
      n_zero += 1
    elseif range < min_q_range_pu * net.baseMVA
      # prosumer limits are in MVAr terms; compare on the guard's pu scale
      n_narrow += 1
    end
  end
  return AutoPfFeatures(
    nbus,
    nbranch,
    max(n_islands, 1),
    rx_median,
    rx_p90,
    share_r_gt_x,
    n_shift,
    n_pv,
    nbus == 0 ? 0.0 : n_pv / nbus,
    n_zero,
    n_narrow,
    n_gen_lim,
    n_gen_pv == 0 ? 0.0 : n_gen_lim / n_gen_pv,
    has_profile,
    plausible,
    length(levels),
  )
end

## ---------------------------------------------------------------------------
## Task 2: decision engine (ordered rule list, first match wins)
## ---------------------------------------------------------------------------

"""
    AutoPfDecision

Result of [`select_auto_pf_strategy`](@ref): the selected `profile`, the
dotted-key `options` the profile wants, and one human-readable `reasons`
entry per triggered rule.
"""
struct AutoPfDecision
  profile::Symbol
  options::Dict{String,Any}
  reasons::Vector{String}
end

# ordered profile rules: (id, predicate, options builder, reason). The FIRST
# matching predicate selects the profile; keep the order auditable, do not
# turn this into nested ifs.
function _auto_pf_profile_rules()
  T = AUTO_PF_THRESHOLDS
  return [
    (
      :profile_assisted,
      f -> f.has_start_profile && f.start_profile_plausible && f.n_bus > T.small_bus_max,
      f -> Dict{String,Any}(
        "power_flow.start_mode.voltage_mode" => "profile_blend",
        "power_flow.start_mode.try_blend_scan" => true,
        "power_flow.start_mode.measure_candidates" => true,
        # DC angle seed stays the choice unless the network is genuinely
        # resistive: the low_rx_share gate here cost case13659pegase its
        # convergence in the calibration sweep (manual dc+blend solved in
        # 8 iterations, classic angles never recovered)
        "power_flow.start_mode.angle_mode" => f.share_r_gt_x < T.resistive_share ? "dc" : "classic",
        "power_flow.autodamp" => true,
        "power_flow.merit.enabled" => true,
        "power_flow.trust_region.enabled" => false,
      ),
      "plausible start profile present: blend it instead of starting flat",
    ),
    (
      :small_default,
      f -> f.n_bus <= T.small_bus_max && f.n_ac_islands == 1 && f.n_phase_shifters == 0,
      f -> begin
        d = Dict{String,Any}("power_flow.autodamp" => true, "power_flow.merit.enabled" => false, "power_flow.trust_region.enabled" => false)
        if f.has_start_profile && f.start_profile_plausible
          d["power_flow.start_mode.voltage_mode"] = "profile_blend"
        else
          d["power_flow.start_mode.flatstart"] = true
          d["power_flow.start_mode.voltage_mode"] = "classic"
        end
        d
      end,
      "small single-island case without phase shifters: plain autodamp",
    ),
    (
      :hard_case,
      f -> f.n_bus >= T.large_bus_min || f.n_phase_shifters > 0 || (f.n_q_zero_range + f.n_q_narrow_range) >= T.many_narrow_q,
      f -> Dict{String,Any}(
        "power_flow.start_mode.start_projection" => true,
        "power_flow.start_mode.try_blend_scan" => true,
        "power_flow.start_mode.try_dc_start" => true,
        "power_flow.start_mode.measure_candidates" => true,
        "power_flow.autodamp" => true,
        "power_flow.merit.enabled" => true,
        "power_flow.trust_region.enabled" => false,
        "power_flow.qlimits.start_iter" => T.qlimit_hard_start_iter,
        "power_flow.qlimits.hysteresis_pu" => T.qlimit_hard_hysteresis,
        "power_flow.qlimits.cooldown_iters" => T.qlimit_hard_cooldown,
      ),
      "large case, phase shifters, or many degenerate Q ranges: full projection with Armijo, trust-region held for escalation",
    ),
    (
      :resistive,
      f -> f.share_r_gt_x > T.resistive_share || f.rx_median > T.resistive_rx_median,
      f -> Dict{String,Any}(
        "power_flow.start_mode.angle_mode" => "classic",
        "power_flow.start_mode.try_dc_start" => false,
        "power_flow.start_current_iteration.enabled" => true,
        "power_flow.start_current_iteration.accept_only_if_improved" => true,
        "power_flow.autodamp" => true,
        "power_flow.merit.enabled" => false,
        "power_flow.trust_region.enabled" => false,
      ),
      "resistive R/X profile: DC angles untrustworthy, guarded current-iteration pre-solve instead",
    ),
    (
      :transmission_flat,
      f -> !(f.has_start_profile && f.start_profile_plausible) && f.share_r_gt_x < T.low_rx_share,
      f -> Dict{String,Any}(
        "power_flow.start_mode.angle_mode" => "dc",
        "power_flow.start_mode.voltage_mode" => "classic",
        "power_flow.start_mode.flatstart" => true,
        "power_flow.start_mode.try_dc_start" => true,
        "power_flow.start_mode.measure_candidates" => true,
        "power_flow.autodamp" => true,
        "power_flow.merit.enabled" => true,
        "power_flow.trust_region.enabled" => false,
      ),
      "no usable start profile on a transmission-style network: flat magnitudes with DC angle seed",
    ),
  ]
end

"""
    select_auto_pf_strategy(features::AutoPfFeatures) -> AutoPfDecision

Select the auto-mode profile and option set for a network: an ordered rule
list picks the profile (first match wins), then profile-independent
feature gates add islands handling, Q-limit strategy, and plausibility
checks. Thresholds live in `AUTO_PF_THRESHOLDS`.
"""
function select_auto_pf_strategy(features::AutoPfFeatures)
  reasons = String[]
  profile = :profile_assisted
  options = Dict{String,Any}()
  for (id, pred, build, reason) in _auto_pf_profile_rules()
    pred(features) || continue
    profile = id
    options = build(features)
    push!(reasons, string("profile ", id, ": ", reason))
    break
  end
  if isempty(options)
    # none matched (medium case, profile present but small share): the
    # profile_assisted recipe is the safest general answer
    profile = :profile_assisted
    options = _auto_pf_profile_rules()[1][3](features)
    push!(reasons, "profile profile_assisted: fallback for a medium case without a stronger signal")
  end
  # profile-independent feature gates
  if features.n_ac_islands > 1
    options["power_flow.islands.enabled"] = true
    push!(reasons, "gate islands: $(features.n_ac_islands) AC islands, island-aware solving enabled")
  end
  if features.n_gen_with_q_limits == 0
    options["power_flow.qlimits.enabled"] = false
    push!(reasons, "gate qlimits: no generator carries Q limits, the whole subsystem is skipped")
  else
    options["power_flow.qlimits.enforcement_mode"] = "active_set"
    options["power_flow.qlimits.start_mode"] = "iteration_or_auto"
    options["power_flow.qlimits.guard"] = true
    push!(reasons, "gate qlimits: active_set enforcement with iteration_or_auto timing and guard enabled")
    if features.n_q_zero_range > 0 || features.n_q_narrow_range > 0
      options["power_flow.qlimits.guard_zero_range_mode"] = "lock_pq"
      options["power_flow.qlimits.guard_narrow_range_mode"] = "lock_pq"
      push!(reasons, "gate qlimits: $(features.n_q_zero_range) zero-range and $(features.n_q_narrow_range) narrow-range generator(s), degenerate ranges locked to PQ")
    end
  end
  # wrong-branch plausibility stays at least warn; the rescue retry exists
  # but is profile-dependent, warn is the safe auto default
  options["power_flow.wrong_branch_detection"] = "warn"
  # the solver-internal rescue ladder stays ACTIVE inside every attempt
  # (calibration sweep: case13659pegase only converges through its
  # alternate-start variant); the auto escalation stages add what the
  # inner ladder does not cover (step-control switch, pre-solve, Q-limit
  # mode flip, APSLF seed, DC fallback). Nesting is bounded: at most four
  # inner variants per outer attempt, and only on failures.
  push!(reasons, "gate retries: solver-internal rescue ladder stays active per attempt; the escalation stages add step-control, pre-solve, Q-limit-mode, APSLF-seed, and DC-fallback strategies on top")
  return AutoPfDecision(profile, options, reasons)
end

## ---------------------------------------------------------------------------
## option application: explicit user keys always win
## ---------------------------------------------------------------------------

# apply one dotted auto option onto a PowerFlowConfig copy; returns the new
# config. Every key the decision engine can emit MUST have a branch here.
function _auto_pf_apply_one(pf::PowerFlowConfig, key::String, value)
  sm = pf.start_mode
  ql = pf.qlimits
  key == "power_flow.autodamp" && return _copy_powerflow_with(pf; autodamp = Bool(value))
  key == "power_flow.merit.enabled" && return _copy_powerflow_with(pf; merit = MeritLineSearchConfig(enabled = Bool(value), armijo_c1 = pf.merit.armijo_c1, scale_p = pf.merit.scale_p, scale_q = pf.merit.scale_q, scale_v = pf.merit.scale_v, fallback_max_mismatch = pf.merit.fallback_max_mismatch))
  key == "power_flow.trust_region.enabled" && return _copy_powerflow_with(pf; trust_region = Bool(value) ? TrustRegionConfig(enabled = true) : TrustRegionConfig())
  key == "power_flow.rescue" && return _copy_powerflow_with(pf; rescue = Bool(value))
  key == "power_flow.wrong_branch_detection" && return _copy_powerflow_with(pf; wrong_branch_detection = Symbol(value))
  key == "power_flow.islands.enabled" && return _copy_powerflow_with(pf; islands = IslandPowerFlowConfig(enabled = Bool(value), mode = pf.islands.mode, reference_policy = pf.islands.reference_policy, diagnostic_continue_after_failure = pf.islands.diagnostic_continue_after_failure))
  key == "power_flow.start_mode.voltage_mode" && return _copy_powerflow_with(pf; start_mode = _copy_start_mode_with(sm; voltage_mode = Symbol(value)))
  key == "power_flow.start_mode.angle_mode" && return _copy_powerflow_with(pf; start_mode = _copy_start_mode_with(sm; angle_mode = Symbol(value)))
  key == "power_flow.start_mode.flatstart" && return _copy_powerflow_with(pf; start_mode = _copy_start_mode_with(sm; flatstart = Bool(value)))
  key == "power_flow.start_mode.try_dc_start" && return _copy_powerflow_with(pf; start_mode = _copy_start_mode_with(sm; try_dc_start = Bool(value)))
  key == "power_flow.start_mode.try_blend_scan" && return _copy_powerflow_with(pf; start_mode = _copy_start_mode_with(sm; try_blend_scan = Bool(value)))
  key == "power_flow.start_mode.measure_candidates" && return _copy_powerflow_with(pf; start_mode = _copy_start_mode_with(sm; measure_candidates = Bool(value)))
  key == "power_flow.start_mode.start_projection" && return _copy_powerflow_with(pf; start_mode = _copy_start_mode_with(sm; start_projection = Bool(value)))
  key == "power_flow.start_mode.blend_lambdas" && return _copy_powerflow_with(pf; start_mode = _copy_start_mode_with(sm; blend_lambdas = Vector{Float64}(value)))
  key == "power_flow.start_mode.accept_unmeasured_dc_start" && return _copy_powerflow_with(pf; start_mode = _copy_start_mode_with(sm; accept_unmeasured_dc_start = Bool(value)))
  key == "power_flow.start_current_iteration.enabled" && return _copy_powerflow_with(pf; start_current_iteration = _copy_start_current_iteration_with(pf.start_current_iteration; enabled = Bool(value)))
  key == "power_flow.start_current_iteration.accept_only_if_improved" && return _copy_powerflow_with(pf; start_current_iteration = _copy_start_current_iteration_with(pf.start_current_iteration; accept_only_if_improved = Bool(value)))
  key == "power_flow.qlimits.enabled" && return _copy_powerflow_with(pf; qlimits = _copy_qlimits_with(ql; ignore_q_limits = !Bool(value)))
  key == "power_flow.qlimits.enforcement_mode" && return _copy_powerflow_with(pf; qlimits = _copy_qlimits_with(ql; enforcement_mode = Symbol(value)))
  key == "power_flow.qlimits.start_mode" && return _copy_powerflow_with(pf; qlimits = _copy_qlimits_with(ql; start_mode = Symbol(value)))
  key == "power_flow.qlimits.start_iter" && return _copy_powerflow_with(pf; qlimits = _copy_qlimits_with(ql; start_iter = Int(value)))
  key == "power_flow.qlimits.hysteresis_pu" && return _copy_powerflow_with(pf; qlimits = _copy_qlimits_with(ql; hysteresis_pu = Float64(value)))
  key == "power_flow.qlimits.cooldown_iters" && return _copy_powerflow_with(pf; qlimits = _copy_qlimits_with(ql; cooldown_iters = Int(value)))
  key == "power_flow.qlimits.guard" && return _copy_powerflow_with(pf; qlimits = _copy_qlimits_with(ql; guard = Bool(value)))
  key == "power_flow.qlimits.guard_zero_range_mode" && return _copy_powerflow_with(pf; qlimits = _copy_qlimits_with(ql; guard_zero_range_mode = Symbol(value)))
  key == "power_flow.qlimits.guard_narrow_range_mode" && return _copy_powerflow_with(pf; qlimits = _copy_qlimits_with(ql; guard_narrow_range_mode = Symbol(value)))
  key == "power_flow.qlimits.guard_freeze_after_repeated_switching" && return _copy_powerflow_with(pf; qlimits = _copy_qlimits_with(ql; guard_freeze_after_repeated_switching = Bool(value)))
  key == "power_flow.qlimits.lock_pv_to_pq_buses" && return _copy_powerflow_with(pf; qlimits = _copy_qlimits_with(ql; lock_pv_to_pq_buses = Vector{Int}(value)))
  key == "power_flow.apslf_start.enabled" && return _copy_powerflow_with(pf; apslf_start = ApslfStartConfig(enabled = Bool(value), order = pf.apslf_start.order))
  key == "power_flow.dc.fallback" && return pf # handled at stage level, never applied blindly
  error("auto powerflow: no application branch for option $(key); add it to _auto_pf_apply_one")
end

"""
    _apply_auto_pf_options(pf, options, user_set_keys) -> (pf2, applied, conflicts)

Apply the decision options onto a PowerFlowConfig copy. Keys the user set
explicitly (dotted keys in `user_set_keys`) are NEVER overwritten; they are
returned in `conflicts` as "key kept by user" entries for the decision log.
"""
function _apply_auto_pf_options(pf::PowerFlowConfig, options::AbstractDict, user_set_keys::Set{String})
  applied = String[]
  conflicts = String[]
  out = pf
  for key in sort!(collect(keys(options)))
    if key in user_set_keys
      push!(conflicts, string(key, " kept at the user value; auto wanted ", repr(options[key])))
      continue
    end
    out = _auto_pf_apply_one(out, key, options[key])
    push!(applied, string(key, " = ", repr(options[key])))
  end
  return out, applied, conflicts
end


## ---------------------------------------------------------------------------
## Task 3: escalation ladder (bounded, tolerance never changes)
## ---------------------------------------------------------------------------

# Q-limit evidence of the LAST attempt, read from the solver status the run
# left on the net. Never throws; missing fields read as zero evidence.
function _auto_pf_qlimit_evidence(net::Net)
  # both calls below return their "nothing to report" value for an unsolved
  # net (nothing / empty Dict, measured 2026-09-06), so neither needs a
  # catch; one here would only hide a programming error inside the evidence
  # the auto mode decides on
  status = rectangular_pf_status(net)
  g(k, d) = status === nothing ? d : (hasproperty(status, k) ? getproperty(status, k) : d)
  counts = qlimit_switch_counts(net)
  return (
    switching_events = Int(g(:pv_pq_switching_events, 0)),
    active_set_changes = Int(g(:qlimit_active_set_changes, 0)),
    reenable_events = Int(g(:qlimit_reenable_events, 0)),
    switch_counts = counts,
  )
end

# ordered escalation stages. Each entry: id, gate(evidence, pf) -> (run,
# skip_reason), transform(pf, evidence) -> (pf2, changed_keys). Stages keep
# the accumulated config of the previous stages (the transform composes).
function _auto_pf_escalation_stages(base_pf::PowerFlowConfig)
  T = AUTO_PF_THRESHOLDS
  stages = Any[]
  push!(stages, (
    id = :L1_step_control,
    gate = (ev, pf) -> (true, ""),
    transform = (pf, ev) -> begin
      if pf.trust_region.enabled
        pf2 = _copy_powerflow_with(pf; trust_region = TrustRegionConfig(), autodamp = true)
        (pf2, ["power_flow.trust_region.enabled = false", "power_flow.autodamp = true"])
      elseif pf.merit.enabled
        # Armijo already active: the next distinct mechanism is trust-region
        # (which excludes autodamp and therefore Armijo, by validation)
        pf2 = _copy_powerflow_with(pf; autodamp = false, merit = MeritLineSearchConfig(), trust_region = TrustRegionConfig(enabled = true))
        (pf2, ["power_flow.autodamp = false", "power_flow.merit.enabled = false", "power_flow.trust_region.enabled = true"])
      else
        pf2 = _copy_powerflow_with(pf; autodamp = true, merit = MeritLineSearchConfig(enabled = true))
        (pf2, ["power_flow.autodamp = true", "power_flow.merit.enabled = true"])
      end
    end,
  ))
  push!(stages, (
    id = :L2_presolve,
    gate = (ev, pf) -> pf.start_current_iteration.enabled ? (false, "skipped_already_enabled") : (true, ""),
    transform = (pf, ev) -> (_copy_powerflow_with(pf; start_current_iteration = _copy_start_current_iteration_with(pf.start_current_iteration; enabled = true, accept_only_if_improved = true)), ["power_flow.start_current_iteration.enabled = true"]),
  ))
  push!(stages, (
    id = :L3_full_projection,
    gate = (ev, pf) -> (true, ""),
    transform = (pf, ev) -> begin
      sm = _copy_start_mode_with(pf.start_mode; start_projection = true, try_blend_scan = true, try_dc_start = true, measure_candidates = true, accept_unmeasured_dc_start = true, blend_lambdas = [0.15, 0.25, 0.5, 0.75, 0.9])
      (_copy_powerflow_with(pf; start_mode = sm), ["power_flow.start_mode.start_projection = true", "power_flow.start_mode.blend_lambdas widened", "power_flow.start_mode.accept_unmeasured_dc_start = true"])
    end,
  ))
  push!(stages, (
    id = :L4_qlimit_relax,
    gate = (ev, pf) -> pf.qlimits.ignore_q_limits ? (false, "skipped_qlimits_disabled") : (true, ""),
    transform = (pf, ev) -> begin
      # buses that chattered in the previous attempt go into the lock list
      locks = sort!(collect(keys(filter(p -> p.second >= pf.qlimits.guard_max_switches, ev.switch_counts))))
      ql = _copy_qlimits_with(pf.qlimits; start_iter = max(pf.qlimits.start_iter, 8), start_mode = :iteration_or_auto, hysteresis_pu = max(pf.qlimits.hysteresis_pu, 2 * T.qlimit_hard_hysteresis), cooldown_iters = max(pf.qlimits.cooldown_iters, 3), guard = true, guard_freeze_after_repeated_switching = true, lock_pv_to_pq_buses = union(pf.qlimits.lock_pv_to_pq_buses, locks))
      changed = ["power_flow.qlimits.start_iter raised", "power_flow.qlimits.hysteresis_pu raised", "power_flow.qlimits.cooldown_iters raised", "power_flow.qlimits.guard_freeze_after_repeated_switching = true"]
      isempty(locks) || push!(changed, string("power_flow.qlimits.lock_pv_to_pq_buses += ", locks))
      (_copy_powerflow_with(pf; qlimits = ql), changed)
    end,
  ))
  push!(stages, (
    id = :L5_qlimit_mode,
    gate = (ev, pf) -> begin
      pf.qlimits.ignore_q_limits && return (false, "skipped_qlimits_disabled")
      (ev.switching_events > 0 || ev.active_set_changes > 0 || !isempty(ev.switch_counts)) || return (false, "skipped_no_qlimit_evidence")
      (true, "")
    end,
    transform = (pf, ev) -> begin
      newmode = pf.qlimits.enforcement_mode === :active_set ? :classic_simultaneous : :active_set
      (_copy_powerflow_with(pf; qlimits = _copy_qlimits_with(pf.qlimits; enforcement_mode = newmode)), [string("power_flow.qlimits.enforcement_mode = ", newmode)])
    end,
  ))
  push!(stages, (
    id = :L6_apslf_seed,
    # no gate: AnalyticLoadFlow.jl is a required dependency, so the APSLF
    # seed is always available to the ladder
    gate = (ev, pf) -> (true, ""),
    transform = (pf, ev) -> begin
      # APSLF ONLY seeds the start values; the rectangular NR that follows is
      # the solver that produces the result (FACTS are not modeled in APSLF)
      sm = _copy_start_mode_with(pf.start_mode; dc_seed_unconditional = false)
      (_copy_powerflow_with(pf; apslf_start = ApslfStartConfig(enabled = true, order = pf.apslf_start.order), start_mode = sm), ["power_flow.apslf_start.enabled = true"])
    end,
  ))
  push!(stages, (
    id = :L7_dc_fallback,
    gate = (ev, pf) -> base_pf.dc.fallback ? (true, "") : (false, "skipped_dc_fallback_disabled"),
    transform = (pf, ev) -> (_copy_powerflow_with(pf; dc = DcPowerFlowConfig(angle_reference_deg = pf.dc.angle_reference_deg, ignore_out_of_service = pf.dc.ignore_out_of_service, fallback = true)), ["power_flow.dc.fallback engaged (result labeled DC approximation, AC status stays non-converged)"]),
  ))
  return stages
end

## ---------------------------------------------------------------------------
## Task 3a: advisory Q-limit hint generator (stable ids, never changes behavior)
## ---------------------------------------------------------------------------

"""
    auto_pf_hints(features, evidence, attempts; converged) -> Vector{NamedTuple}

Generate the advisory hint list from the recorded diagnostics: each entry
is `(id = :stable_key, text = "...")`. Purely derived from existing
diagnostics; hints never mutate configuration or behavior.
"""
function auto_pf_hints(features::AutoPfFeatures, evidence, attempts::Vector; converged::Bool)
  hints = NamedTuple[]
  if features.n_gen_with_q_limits == 0
    push!(hints, (id = :q_no_limits_defined, text = "No generator Q limits found in the source data. PV buses hold voltage without reactive bounds; results may be optimistic."))
  end
  if features.n_q_zero_range > 0
    push!(hints, (id = :q_zero_range_locked, text = "$(features.n_q_zero_range) generator(s): qmin equals qmax, treated as PQ. Check generator Q limits in the source data."))
  end
  if features.n_q_narrow_range > 0
    push!(hints, (id = :q_narrow_range, text = "$(features.n_q_narrow_range) generator(s): Q range below min_q_range_pu, PV control not enforceable. Widen limits or accept PQ behavior."))
  end
  frozen = [b for (b, c) in evidence.switch_counts if c >= 10]
  if !isempty(frozen)
    sort!(frozen)
    push!(hints, (id = :q_switching_frozen, text = "Bus(es) $(join(frozen, ", ")): PV/PQ switching froze after repeated switches. Result is valid but bus type is guard-decided; consider lock_pv_to_pq_buses."))
  end
  modeflip = findfirst(a -> a.stage === :L5_qlimit_mode && a.converged, attempts)
  if modeflip !== nothing
    push!(hints, (id = :q_mode_pinned, text = "Case converged only under a different Q-limit enforcement mode. Consider pinning power_flow.qlimits.enforcement_mode for this network (see auto_mode_decision.log)."))
  end
  if !converged && evidence.switching_events > 0
    busids = sort!(collect(keys(evidence.switch_counts)))
    push!(hints, (id = :q_failure_with_switching, text = "Non-convergence coincides with $(evidence.switching_events) Q-limit switching events on buses $(join(busids, ", ")). Likely limit data problem; verify qmin/qmax of these generators."))
  end
  return hints
end

## ---------------------------------------------------------------------------
## startup latency hint (Task 5): once per process, only without sysimage/exe
## ---------------------------------------------------------------------------

const _STARTUP_LATENCY_HINT_SHOWN = Ref(false)

# one concise info line on the first run of a fresh native Julia process;
# sysimage and app sessions never see it (webui_runtime_flavor), and
# output.startup_latency_hint = false silences it entirely
function _maybe_print_startup_latency_hint(cfg::SparlectraConfig)
  _STARTUP_LATENCY_HINT_SHOWN[] && return nothing
  _STARTUP_LATENCY_HINT_SHOWN[] = true
  cfg.output.startup_latency_hint || return nothing
  flavor = try
    webui_runtime_flavor()
  catch
    (kind = :native, built = nothing)
  end
  flavor.kind === :native || return nothing
  @info "First run in a fresh Julia process includes compilation and is slow; later runs are fast. Keep the process running, or build a sysimage (buildSysimage(); usable in any session with julia --sysimage). See the integration guide. Silence with output.startup_latency_hint: false."
  return nothing
end

## ---------------------------------------------------------------------------
## orchestrator: base attempt plus bounded escalation, full recording
## ---------------------------------------------------------------------------

# auto-mode record store: keyed by the Net object so the service layer can
# read the record from result.net without changing the result struct
const _AUTO_PF_RECORDS = Base.WeakKeyDict{Any,Any}()
const _AUTO_PF_RECORDS_LOCK = ReentrantLock()

_set_auto_pf_record!(net::Net, rec) = lock(() -> (_AUTO_PF_RECORDS[net] = rec), _AUTO_PF_RECORDS_LOCK)

"""
    auto_pf_record(net) -> Union{Nothing,NamedTuple}

The auto-mode decision/escalation record of the most recent auto run on
this net (`nothing` when the run was manual).
"""
auto_pf_record(net::Net) = lock(() -> get(_AUTO_PF_RECORDS, net, nothing), _AUTO_PF_RECORDS_LOCK)

"""
    _write_auto_mode_decision_log(path, rec)

Write the human-readable auto-mode decision artifact: feature summary,
selected profile with per-rule reasons, applied options, precedence
conflicts (user key kept), every escalation attempt with its changed keys,
solver path, iteration count and Q-limit evidence, and the hint list.
"""
function _write_auto_mode_decision_log(path::AbstractString, rec)
  open(path, "w") do io
    println(io, "Auto power-flow mode decision log")
    println(io, "=================================")
    f = rec.features
    println(io, "network: ", f.n_bus, " buses, ", f.n_branch, " branches, ", f.n_ac_islands, " AC island(s), ", f.n_voltage_levels, " voltage level(s)")
    println(io, "R/X: median ", round(f.rx_median; digits = 3), ", p90 ", round(f.rx_p90; digits = 3), ", share R>X ", round(f.share_r_gt_x; digits = 3))
    println(io, "phase shifters: ", f.n_phase_shifters, "; PV buses: ", f.n_pv, " (share ", round(f.pv_share; digits = 3), ")")
    println(io, "generator Q limits: ", f.n_gen_with_q_limits, " with limits, ", f.n_q_zero_range, " zero-range, ", f.n_q_narrow_range, " narrow-range")
    println(io, "start profile: ", f.has_start_profile ? (f.start_profile_plausible ? "present and plausible" : "present but implausible") : "absent")
    println(io)
    println(io, "selected profile: ", rec.profile)
    for r in rec.reasons
      println(io, "  reason: ", r)
    end
    println(io)
    println(io, "applied options (", length(rec.applied), "):")
    for a in rec.applied
      println(io, "  ", a)
    end
    if !isempty(rec.conflicts)
      println(io, "precedence conflicts (user key kept, auto value dropped):")
      for c in rec.conflicts
        println(io, "  ", c)
      end
    end
    println(io)
    println(io, "attempts:")
    for a in rec.attempts
      if isempty(a.solver)
        println(io, "  ", a.stage, ": ", haskey(pairs(a.qlimit_evidence), :skipped) ? a.qlimit_evidence.skipped : "skipped")
      else
        println(io, "  ", a.stage, ": ", a.converged ? "converged" : "not converged", " after ", a.iterations, " iteration(s), ", a.elapsed_s, " s, solver ", a.solver)
        for c in a.changed
          println(io, "    changed: ", c)
        end
        ev = a.qlimit_evidence
        if haskey(pairs(ev), :switching_events) && (ev.switching_events > 0 || ev.active_set_changes > 0)
          println(io, "    q-limit evidence: ", ev.switching_events, " switching event(s), ", ev.active_set_changes, " active-set change(s), buses ", sort!(collect(keys(ev.switch_counts))))
        end
      end
    end
    println(io)
    println(io, "final stage: ", rec.final_stage, "; final solver: ", rec.final_solver, "; converged: ", rec.converged)
    if !isempty(rec.hints)
      println(io)
      println(io, "hints:")
      for h in rec.hints
        println(io, "  [", h.id, "] ", h.text)
      end
    end
  end
  return path
end

# island solver failures surface as ErrorException (they throw instead of
# returning a status); an auto attempt treats that as "not converged" and
# escalates instead of dying. Anything else is a programming error and is
# rethrown unchanged.
function _auto_pf_try_execute!(net::Net, run_cfg::SparlectraConfig; performance_profile = nothing)
  t0 = time()
  try
    return _execute_sparlectra_powerflow!(net, run_cfg; performance_profile = performance_profile), nothing
  catch err
    err isa ErrorException || rethrow()
    # the wall time up to the throw still belongs to the attempt record
    return (iterations = 0, erg = 2, elapsed_s = time() - t0, solver_elapsed_s = 0.0, control_status = :none), err
  end
end

function _auto_pf_solver_label(pf::PowerFlowConfig, execution)
  pf.apslf_start.enabled && return "rectangular (apslf-seeded)"
  return "rectangular"
end

"""
    _execute_auto_sparlectra_powerflow!(net, cfg; performance_profile) -> execution

Auto-mode drive of the power flow: feature extraction, strategy selection,
precedence-safe application, base attempt, and the bounded escalation
ladder. The record (features, decision, attempts, hints) is stored on the
net (`auto_pf_record`) for the service layer. Tolerance is never changed by
any stage.
"""
function _execute_auto_sparlectra_powerflow!(net::Net, cfg::SparlectraConfig; performance_profile = nothing)
  features = collect_auto_pf_features(net; min_q_range_pu = cfg.powerflow.qlimits.guard_min_q_range_pu)
  decision = select_auto_pf_strategy(features)
  pf, applied, conflicts = _apply_auto_pf_options(cfg.powerflow, decision.options, cfg.user_set_keys)
  snap = _snapshot_start_voltages(net)
  attempts = NamedTuple[]
  run_cfg = _copy_sparlectra_with_powerflow(cfg, pf)
  execution, base_err = _auto_pf_try_execute!(net, run_cfg; performance_profile = performance_profile)
  evidence = _auto_pf_qlimit_evidence(net)
  push!(attempts, (stage = :base, changed = base_err === nothing ? applied : vcat(applied, string("solver error: ", first(sprint(showerror, base_err), 160))), solver = _auto_pf_solver_label(pf, execution), iterations = execution.iterations, converged = execution.erg == 0, elapsed_s = round(execution.elapsed_s; digits = 3), qlimit_evidence = evidence))
  final_stage = :base
  if execution.erg != 0
    stages = _auto_pf_escalation_stages(cfg.powerflow)
    n_run = 0
    for st in stages
      n_run >= AUTO_PF_THRESHOLDS.escalation_max_stages && break
      runit, skip_reason = st.gate(evidence, pf)
      if !runit
        push!(attempts, (stage = st.id, changed = String[], solver = "", iterations = 0, converged = false, elapsed_s = 0.0, qlimit_evidence = (skipped = skip_reason,)))
        continue
      end
      if st.id === :L7_dc_fallback
        # last resort: one more base-strategy run with dc.fallback armed;
        # the solver labels the DC state honestly (AC stays non-converged)
        pf = _copy_powerflow_with(pf; dc = DcPowerFlowConfig(angle_reference_deg = pf.dc.angle_reference_deg, ignore_out_of_service = pf.dc.ignore_out_of_service, fallback = true))
        changed = ["power_flow.dc.fallback = true (last resort, DC approximation)"]
      else
        pf, changed = st.transform(pf, evidence)
      end
      n_run += 1
      _restore_start_voltages!(net, snap)
      run_cfg = _copy_sparlectra_with_powerflow(cfg, pf)
      execution, stage_err = _auto_pf_try_execute!(net, run_cfg; performance_profile = performance_profile)
      evidence = _auto_pf_qlimit_evidence(net)
      push!(attempts, (stage = st.id, changed = stage_err === nothing ? changed : vcat(changed, string("solver error: ", first(sprint(showerror, stage_err), 160))), solver = _auto_pf_solver_label(pf, execution), iterations = execution.iterations, converged = execution.erg == 0, elapsed_s = round(execution.elapsed_s; digits = 3), qlimit_evidence = evidence))
      final_stage = st.id
      execution.erg == 0 && break
    end
  end
  hints = auto_pf_hints(features, evidence, attempts; converged = execution.erg == 0)
  final_solver = isempty(attempts) ? "rectangular" : last(filter(a -> !isempty(a.solver), attempts)).solver
  _set_auto_pf_record!(net, (
    features = features,
    profile = decision.profile,
    reasons = decision.reasons,
    applied = applied,
    conflicts = conflicts,
    attempts = attempts,
    final_stage = final_stage,
    final_solver = final_solver,
    hints = hints,
    converged = execution.erg == 0,
  ))
  return execution
end
