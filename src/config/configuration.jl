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

# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 20.5.2026
# file: src/config/configuration.jl
# purpose: typed configuration surface: all *Config structs including
#          SparlectraConfig, YAML loading and validation, and the active
#          global configuration state

using SHA

"""
    StartModeConfig

Typed power-flow start-option configuration. These fields collect the
flat-start and rectangular start-projection controls that otherwise tend to be
forwarded as long keyword lists through example and benchmark call chains.
"""
Base.@kwdef struct StartModeConfig
  flatstart::Bool = false
  angle_mode::Symbol = :dc
  voltage_mode::Symbol = :profile_blend
  profile_source::Symbol = :matpower_reference
  start_projection::Bool = false
  try_dc_start::Bool = true
  try_blend_scan::Bool = true
  branch_guard::Bool = true
  measure_candidates::Bool = true
  accept_unmeasured_dc_start::Bool = false
  dc_seed_unconditional::Bool = false
  reuse_import_data::Bool = true
  blend_lambdas::Vector{Float64} = [0.25, 0.5, 0.75]
  dc_angle_limit_deg::Float64 = 60.0
end

"""
    StartCurrentIterationConfig

Guarded fixed-point current-injection pre-solve configuration. The stage is
disabled by default and, when enabled, only prepares the initial voltage profile
before the normal rectangular Newton-Raphson solve.
"""
Base.@kwdef struct StartCurrentIterationConfig
  enabled::Bool = false
  max_iter::Int = 10
  tol::Float64 = 1.0e-3
  damping::Float64 = 0.5
  accept_only_if_improved::Bool = true
  min_improvement_factor::Float64 = 0.98
  vm_min_pu::Float64 = 0.5
  vm_max_pu::Float64 = 1.5
  max_angle_step_deg::Float64 = 30.0
  only_for_large_cases::Bool = false
end

"""
    MeritLineSearchConfig

Typed configuration for the optional Armijo merit-function line search used
as an alternative acceptance criterion inside the rectangular Newton-Raphson
autodamp backtracking loop. Disabled by default; when disabled the solver
behaves exactly as before this feature was added. See
[Merit-Function Line Search](@ref) for the theoretical background.
"""
Base.@kwdef struct MeritLineSearchConfig
  enabled::Bool = false
  armijo_c1::Float64 = 1.0e-4
  scale_p::Float64 = 1.0
  scale_q::Float64 = 1.0
  scale_v::Float64 = 1.0
  fallback_max_mismatch::Bool = true
end

"""
    TrustRegionConfig

Typed configuration for the optional scaled-Newton trust-region step control
in the rectangular Newton-Raphson solver: an alternative to `autodamp` that
caps the Newton step norm at an adaptive radius and accepts/rejects steps by
merit decrease. Disabled by default; mutually exclusive with `autodamp`. See
[Trust-Region Step Control](@ref) for the theoretical background.

`step_mode = :scaled` (default) rescales the full Newton direction to the
radius when it exceeds it, leaving pre-existing behavior byte-for-byte
unchanged. `step_mode = :dogleg` blends the Newton direction with a
steepest-descent (Cauchy) step along the dogleg path when the radius shrinks
below the Newton step norm, trading some convergence speed for graceful
degradation when the Newton direction becomes a poor descent direction. See
[Trust-Region Step Control](@ref), "Dogleg step mode".
"""
Base.@kwdef struct TrustRegionConfig
  enabled::Bool = false
  initial_radius::Float64 = 1.0
  min_radius::Float64 = 1e-4
  max_radius::Float64 = 10.0
  eta_accept::Float64 = 0.1
  shrink_factor::Float64 = 0.5
  expand_factor::Float64 = 2.0
  expand_threshold::Float64 = 0.75
  step_mode::Symbol = :scaled
end

"""
    QLimitConfig

Typed reactive-power limit switching configuration used by power-flow runners.
"""
Base.@kwdef struct QLimitConfig
  start_iter::Int = 2
  start_mode::Symbol = :iteration
  auto_q_delta_pu::Float64 = 1e-4
  hysteresis_pu::Float64 = 0.01
  cooldown_iters::Int = 1
  guard::Bool = false
  guard_min_q_range_pu::Float64 = 1e-4
  guard_zero_range_mode::Symbol = :lock_pq
  guard_narrow_range_mode::Symbol = :prefer_pq
  guard_max_switches::Int = 10
  guard_freeze_after_repeated_switching::Bool = true
  guard_accept_bounded_violations::Bool = false
  guard_max_remaining_violations::Int = 0
  guard_violation_mode::Symbol = :delayed_switch
  guard_violation_threshold_pu::Float64 = 1e-4
  guard_log::Bool = true
  trace_buses::Vector{Int} = Int[]
  lock_pv_to_pq_buses::Vector{Int} = Int[]
  ignore_q_limits::Bool = false
  enforcement_mode::Symbol = :active_set
end

"""
    IslandPowerFlowConfig

Configuration for AC-island-aware power-flow diagnostics.
"""
Base.@kwdef struct IslandPowerFlowConfig
  enabled::Bool = true
  mode::Symbol = :solve_independent
  reference_policy::Symbol = :matpower_like
  diagnostic_continue_after_failure::Bool = true
end

"""
    DcPowerFlowConfig

Typed configuration for the standalone DC power flow ([`rundcpf!`](@ref)),
used when `power_flow.solver == :dc`.

`angle_reference_deg` is the uniform angle offset added to every bus after
the slack-referenced linear solve (the slack bus itself is fixed at this
reference) — see [`solve_dc_powerflow`](@ref) for why this post-hoc shift is
exact. `ignore_out_of_service` documents that `status == 0` branches are
always excluded from the B′ assembly (there is currently no supported way
to include them; the field exists for forward compatibility with the YAML
schema, not as a live toggle).
"""
Base.@kwdef struct DcPowerFlowConfig
  angle_reference_deg::Float64 = 0.0
  ignore_out_of_service::Bool = true
  # Run the standalone DC power flow when the AC solve (and the rescue
  # ladder, if enabled) did not converge. The AC status stays "not
  # converged" — the fallback only leaves a usable DC state in the net.
  fallback::Bool = false
end

"""
    DistributedSlackConfig

Distributed active-power slack (issue #192): the REF bus keeps the angle
reference while the island's active-power imbalance is absorbed by a set of
participating generators via one scalar `lambda_P` per island.

# Fields
- `enabled::Bool`: off by default — disabled reproduces the classical
  single-slack behavior bit-for-bit.
- `p_mode::Symbol`: how raw participation weights are derived —
  `:pg_weighted` (scheduled Pg), `:pmax_weighted` (maxP),
  `:headroom_weighted` (`max(maxP − Pg, 0)`), `:imported`
  (`ProSumer.participationFactor`, filled from MATPOWER `APF` /
  CGMES `GeneratingUnit.normalPF`), `:explicit` (the `weights` table).
- `respect_p_limits::Bool`: diagnostic only — WARN when a participant's
  corrected P leaves `[minP, maxP]`; no re-dispatch.
- `fallback::Symbol`: `:error` throws when an island has no valid participant;
  `:ref_only` falls back to the classical slack for that island with a warning.
- `weights::Dict{String,Float64}`: `:explicit` mode only — bus name (or bus
  index as string) → weight.
"""
Base.@kwdef struct DistributedSlackConfig
  enabled::Bool = false
  p_mode::Symbol = :pg_weighted
  respect_p_limits::Bool = true
  fallback::Symbol = :error
  weights::Dict{String,Float64} = Dict{String,Float64}()
end

const DISTRIBUTED_SLACK_P_MODE_VALUES = (:pg_weighted, :pmax_weighted, :headroom_weighted, :imported, :explicit)
const DISTRIBUTED_SLACK_FALLBACK_VALUES = (:error, :ref_only)

"""
    ContingencyConfig

Configuration of the N-1 contingency batch (issue #331).

# Fields
- `rescue_ladder::Vector{Symbol}`: the per-case start-value ladder, an ordered,
  duplicate-free subset of `(:warm, :apslf, :dc, :flat)`. Default `[:warm]`
  reproduces the pre-#331 single warm solve. Each stage is one bounded solve
  with a distinct start recipe, tried in order until one converges. The allowed
  stages and their recipes are documented on [`runContingencies!`](@ref); the
  set is validated (subset, no duplicates) by `_validate_contingency_ladder`.
- `screening_mode::Symbol`: contingency screening on the base factorization
  (scenario task D5): `:off` (default) gives every scenario the full solve,
  `:flag` estimates every non-islanding outage with one Woodbury-corrected
  Newton step on the base Jacobian and runs the full solve only for flagged
  scenarios, `:only` reports the estimates without full runs (islanding and
  failed-screen scenarios still get the full solve). `:off` is the default
  DELIBERATELY: the gate calibration (2026-09-03) showed classes of outages
  a residual step cannot see (Q-capability loss behind a zero residual), so
  `:flag` is an opt-in for networks where the user has checked the
  screening share and the margins once. This config key drives the SERVICE
  path; the programmatic keyword default is `:off` as well.
- `screening_margin_pct::Float64`: the flagging margin (default `10.0`): a
  scenario is flagged for the full run when an estimated loading reaches
  `100 - margin` percent, or an estimated voltage comes within `margin`
  percent of a band limit.
"""
Base.@kwdef struct ContingencyConfig
  rescue_ladder::Vector{Symbol} = [:warm]
  screening_mode::Symbol = :off
  screening_margin_pct::Float64 = 10.0
end

const CONTINGENCY_SCREENING_MODE_VALUES = (:off, :flag, :only)

"""
    ExternalGridConfig

Compute the marked slack bus as a non-ideal external-grid source (issue
#299): the reference voltage moves to a hidden internal bus behind the
feeder impedance `z = Un²/Sk''`, so the connection-bus voltage reacts to
loading instead of being held ideally stiff.

# Fields
- `enabled::Bool`: off by default — the classical ideal slack.
- `source::Symbol`: where `Sk''`/`R/X` come from — `:auto` prefers the
  values a CGMES delivery declares on the slack bus's
  `ExternalNetworkInjection` and falls back to the config numbers below
  (MATPOWER/DTF carry no such data); `:config` always uses the config
  numbers.
- `sk_MVA::Float64`: initial symmetrical short-circuit power of the feeder.
- `rx::Float64`: its R/X ratio.
"""
Base.@kwdef struct ExternalGridConfig
  enabled::Bool = false
  source::Symbol = :auto
  sk_MVA::Float64 = 2000.0
  rx::Float64 = 0.1
end

const EXTERNAL_GRID_SOURCE_VALUES = (:auto, :config)

"""
    ApslfConfig

Typed configuration for the AnalyticLoadFlow.jl-backed analytic power-series
solver (`ApslfSolver`), used when `power_flow.solver == :apslf`.
"""
Base.@kwdef struct ApslfConfig
  order::Int = 40
  use_pade::Bool = true
  nr_polish::Bool = true
end

"""
    ApslfStartConfig

Typed configuration for using the analytic power-series solver as a start-value
generator ahead of the rectangular Newton-Raphson solve. Deliberately has no
`use_pade`/`nr_polish` fields: polishing is left to the downstream NR solve, so
the generator always runs with `nr_polish=false` internally.
"""
Base.@kwdef struct ApslfStartConfig
  enabled::Bool = false
  order::Int = 40
end

"""
    PowerFlowConfig

Typed power-flow configuration. It owns solver tolerances, sparse execution
settings, automatic damping, start-mode controls, and Q-limit behavior.
"""
Base.@kwdef struct PowerFlowConfig
  method::Symbol = :rectangular
  # :manual keeps every option as configured; :auto inspects the imported
  # network and fills the start/step/Q-limit strategy for keys the user
  # did not set explicitly (see auto_powerflow.jl and integration.md)
  mode::Symbol = :manual
  solver::Symbol = :rectangular
  linear_solver::Symbol = :umfpack
  apslf::ApslfConfig = ApslfConfig()
  apslf_start::ApslfStartConfig = ApslfStartConfig()
  # convergence bound for the LARGEST SINGLE bus mismatch (infinity norm
  # over active and reactive residuals alike; PV voltage rows share the
  # same per-unit bound as a voltage quantity). Readable physically as
  # tol * baseMVA: the 1e-8 default equals 1 W at a 100 MVA base
  # (task_tol_watts; the run log and diagnostics print the equivalent)
  tol::Float64 = 1.0e-8
  # Physical spelling of the same bound (task_tol_watts part B): when set,
  # tol_MW WINS over tol and is converted with the network's own base
  # (tol = tol_MW / baseMVA) at the moment the tolerance meets the net,
  # because the base is unknown while the configuration is read. `nothing`
  # means "not set" and changes nothing.
  tol_MW::Union{Nothing,Float64} = nothing
  max_iter::Int = 30
  autodamp::Bool = false
  autodamp_min::Float64 = 0.05
  # Promote the strongest injection to slack when none is registered
  # (ensureSlack!); off by default so data errors stay visible.
  auto_slack::Bool = false
  # Retry a non-converged AC solve from the original start state with a fixed
  # strategy ladder (alternate start, autodamp, DC-seeded projection,
  # settled Q-limits). ON by default: only failed runs pay for the retries,
  # and a user who gets a result instead of a divergence message is the
  # whole point. Turn off for solver studies that must see the raw failure.
  rescue::Bool = true
  wrong_branch_detection::Symbol = :warn
  wrong_branch_rescue::Bool = false
  wrong_branch_min_vm_pu::Float64 = 0.70
  wrong_branch_max_vm_pu::Float64 = 1.30
  wrong_branch_max_angle_spread_deg::Float64 = 180.0
  wrong_branch_max_branch_angle_deg::Float64 = 90.0
  wrong_branch_min_low_vm_count::Int = 1
  wrong_branch_rescue_max_attempts::Int = 2
  rectangular_workspace_reuse::Bool = true
  rectangular_preallocate_workspace::Symbol = :auto
  rectangular_workspace_min_buses::Int = 1000
  islands_enabled::Bool = true
  islands_mode::Symbol = :solve_independent
  islands_reference_policy::Symbol = :matpower_like
  start_mode::StartModeConfig = StartModeConfig()
  start_current_iteration::StartCurrentIterationConfig = StartCurrentIterationConfig()
  merit::MeritLineSearchConfig = MeritLineSearchConfig()
  trust_region::TrustRegionConfig = TrustRegionConfig()
  qlimits::QLimitConfig = QLimitConfig()
  islands::IslandPowerFlowConfig = IslandPowerFlowConfig()
  dc::DcPowerFlowConfig = DcPowerFlowConfig()
  distributed_slack::DistributedSlackConfig = DistributedSlackConfig()
  external_grid::ExternalGridConfig = ExternalGridConfig()
end

const WRONG_BRANCH_DETECTION_VALUES = [:off, :warn, :fail, :rescue]
const POWERFLOW_SOLVER_VALUES = (:rectangular, :apslf, :dc)
const POWERFLOW_LINEAR_SOLVER_VALUES = (:umfpack, :umfpack_reuse)
# (single definition further below, next to the other island constants; a
# duplicate Vector-typed definition used to live here)

"""
    ObservabilityConfig

State-estimation observability diagnostic configuration.
"""
Base.@kwdef struct ObservabilityConfig
  enabled::Bool = true
end

"""
    StateEstimationConfig

Typed state-estimation configuration for future SE runners and diagnostics.

`imag_activation_iteration` is the first WLS iteration in which branch
current-magnitude measurements (`ImagMeas`) enter the update step; earlier
iterations run without them (flat-start protection).
`report_residual_correlation` enables the residual-correlation (K-matrix)
columns in the bad-data diagnostics report.
"""
# No `sparse` field: sparse matrices are mandatory since 0.9.x, and a field
# that can never be anything but true is not a setting. The old key is still
# REFUSED by name below, so a stored configuration carrying it says what
# happened instead of failing with "unknown key".
Base.@kwdef struct StateEstimationConfig
  enabled::Bool = true
  method::Symbol = :wls
  # tol at the finite-difference noise floor, not below it. The estimator
  # builds its Jacobian by forward differences with step `jac_eps`, so the
  # fixed point is only located to about that accuracy; a tighter tolerance
  # is unreachable and the run clamps it with a log line. The former 1.0e-8
  # therefore only LOOKED tighter than the 1.0e-6 it actually got.
  tol::Float64 = 1.0e-6
  # One cap, everywhere. The service used 30 and the Web UI form 50 while
  # this said 20, so the same run reached three different limits depending
  # on the entry point.
  #
  # The value is 50, the largest of the three, and NOT 30: a CGMES run with
  # released taps needs between 36 and 40 iterations in its first solve and
  # failed at 30 (maintainer, 2026-09-06, run a023884e). The earlier
  # reasoning here, "a run that has not converged by 30 does not converge at
  # 50 either", was measured on sets without released taps and is wrong in
  # general. What the run reports afterwards is the iteration count of the
  # LAST solve (3 in that case), so the first solve's cost is invisible in
  # the result, which is what made 30 look sufficient.
  max_iter::Int = 50
  flatstart::Bool = true
  jac_eps::Float64 = 1.0e-6
  update_net::Bool = true
  pmu_ref_offset::Symbol = :auto
  imag_activation_iteration::Int = 2
  report_residual_correlation::Bool = false
  update_shunts::Bool = false
  # tap write-back: estimation must never silently overwrite model
  # data; the FIXED mechanical positions reach the branch fields only on
  # explicit request and only from a converged, fixed run
  update_taps::Bool = false
  robust::Bool = false
  robust_start_iteration::Int = 3
  # bad-data thresholds (0.10.0, GUI-exposed). Two decisions, and since
  # task_se_bad_data_v0100 both read the SAME quantity, the normalized
  # residual rn = r_i/sqrt(Omega_ii):
  #   k_eliminate  from here a row is REMOVED from the estimate
  #   k_suppress   from here a row is DOWN-WEIGHTED (:replacement mode)
  # robust_mode selects the solve-weight family: :off plain WLS, :staged the
  # two-stage tangential/constant R modification with configurable knees
  # robust_k1/robust_k2 (3/6 reproduces the historical `robust = true`
  # bitwise; those knees deliberately keep their classic |r|/sigma scale),
  # :replacement pins rows past k_suppress to the fixed suppression_sigma
  # (in the measurement's unit) during the solve. Statistics and reports
  # always use the original sigma; virtual rows (sigma <= 1e-6) stay exempt.
  # The Bool `robust` remains a legacy alias for :staged and applies only
  # while robust_mode is :off.
  #
  # k_suppress = 4.0 is calibrated, not chosen: measured over four networks
  # (sp_case14, sp_case60, sp_case188 and case118) with gross errors seeded
  # at 5, 10 and 20 sigma. At 3.0 the false suppressions start, up to 21
  # healthy rows on the largest of them; from 5.0 upward the weak 5 sigma
  # errors are missed. 4.0 is the only value that catches the weak errors
  # without suppressing healthy rows on any of the four. The former 6.0 was a threshold on the RAW ratio |r|/sigma and
  # corresponds to rn 6 to 11 depending on redundancy, which is why the
  # suppression was weakest exactly at nearly critical rows.
  k_eliminate::Float64 = 3.0
  robust_mode::Symbol = :off
  robust_k1::Float64 = 3.0
  robust_k2::Float64 = 6.0
  k_suppress::Float64 = 4.0
  suppression_sigma::Float64 = 2000.0
  # budget of the sequential elimination: how many rows one diagnostics run
  # may remove at most. Not a threshold but a stop condition, so it is an
  # Int and 0 disables the elimination entirely. 3 is the historical
  # default; a set that needs more than three removals to become plausible
  # is usually a model or topology problem, not bad telemetry.
  max_eliminations::Int = 3
  rank_tol_factor::Float64 = 10.0
  takahashi_min_states::Int = 200
  # IaMeas activity gate fallback: threshold is 3 * sigma of the paired
  # ImagMeas at the same end, or 3 * this floor (ampere) when no magnitude
  # measurement is paired (a current angle is meaningless near zero current)
  ia_current_floor_A::Float64 = 10.0
  # topology validation (all three stages ADVISORY: findings and
  # warnings only, never a status/measurement/model mutation). The k values
  # are sigma multiples of the respective plausibility checks; open
  # (out-of-service) elements are exempt from the plausibility checks, only
  # the status-contradiction check looks at them.
  topology_precheck::Bool = true
  topology_open_flow_k::Float64 = 4.0
  topology_dead_flow_k::Float64 = 3.0
  topology_voltage_k::Float64 = 4.0
  topology_kcl_k::Float64 = 4.0
  topology_cluster_min::Int = 3
  observability::ObservabilityConfig = ObservabilityConfig()
end

"""
    CGMESImportConfig

Options of the `cgmes_import` configuration block — the ENTSO-E CGMES import
(see `importCGMES`). `path` accepts a folder, a ZIP, or several of both
(base case plus boundary set) separated by `;`.

# Fields
- `path::String`: delivery location(s); `;`-separated for multi-part deliveries.
- `base_mva::Float64`: system base, which CGMES does not define.
- `require_boundary::Bool`: fail when topology references stay unresolved.
- `tap_control::Bool`: start from the SSH tap positions and attach the
  CGMES-defined outer-loop tap controllers instead of importing the solved
  `SvTapStep` positions as fixed taps.
- `machine_control::Bool`: attach outer-loop remote voltage controllers
  (`MachineVoltageControl`) for machines whose voltage `RegulatingControl`
  points at a different bus, instead of holding those machines PV at their
  own bus.
- `ignore_connected::Bool`: diagnostic override that treats every terminal as
  connected, for snapshots whose SSH flags contradict their own SV state.
- `vset_min_pu`, `vset_max_pu::Float64`: plausibility band for a voltage
  `RegulatingControl.targetValue`, in p.u. of the regulated bus's nominal
  voltage. A target outside the band is treated as a placeholder: it is ignored
  and the unit is held PV at the bus voltage derived from the nominal data, with
  a `warning:` message. Widen the band to accept a delivery's own values, or set
  `vset_min_pu = 0` and a large `vset_max_pu` to disable the check entirely.
- `multi_slack::Bool`: give every electrical island its own SV-declared angle
  reference (at most one per island). Required for multi-island deliveries —
  without it every island beyond the primary one has no reference and the
  island-wise power flow refuses to run. Disable only to force the legacy
  single-reference behavior.
- `start_values::Symbol`: Newton-Raphson start state for CGMES runs. `:flat`
  (default) uses a synthetic flat start — the solver earns the solution
  itself; `:sv` starts from the delivery's imported `SvVoltage` state and
  force-disables the competing start-value machines for the run. On CGMES
  runs this key wins over `power_flow.flatstart` /
  `power_flow.start_mode.flatstart`; MATPOWER and DTF runs ignore it. The
  SV comparison artifact (`sv_compare.csv`) is written either way.
"""
Base.@kwdef struct CGMESImportConfig
  path::String = ""
  base_mva::Float64 = 100.0
  require_boundary::Bool = true
  tap_control::Bool = false
  machine_control::Bool = false
  ignore_connected::Bool = false
  vset_min_pu::Float64 = 0.5
  vset_max_pu::Float64 = 1.5
  multi_slack::Bool = true
  start_values::Symbol = :auto
  placeholder_guards::Symbol = :warn_skip
  # Reconstruct missing nominal voltages (SV snap + transformer ratedU +
  # propagation) when the BaseVoltage catalog is absent — see
  # _inferNominalVoltages. Pair with require_boundary = false.
  infer_base_voltages::Bool = false
  # HVDC converter handling (#297 Draft B): :injections keeps the Stage-0
  # fixed PCC injections, :paired_control additionally attaches one
  # steerable HvdcPairControl per detected converter pair.
  hvdc_mode::Symbol = :injections
end

"""
    ShortCircuitConfig

Options of the `short_circuit` configuration block — the IEC 60909 balanced
short-circuit evaluation (`runShortCircuit!`).

# Fields
- `c_factor::Float64`: scalar override for the IEC 60909-0 voltage factor
  `c`. `0.0` (default) selects the hardcoded Table-1 values by nominal
  voltage level and case (`c_max`/`c_min`); a positive value replaces the
  table for every bus — intended for expert/verification runs where a worked
  example prescribes the factor. A configurable per-voltage-level table is
  not configurable.
"""
Base.@kwdef struct ShortCircuitConfig
  c_factor::Float64 = 0.0
  # all-bus Thevenin sweep strategy: :auto picks the Takahashi selected
  # inverse for islands with at least takahashi_min_buses buses (measured
  # 34x to 264x over per-bus solves, machine-precision identical, automatic
  # per-island fallback where inapplicable) and plain solves below;
  # :solves and :takahashi force one method for every island
  sweep_method::Symbol = :auto
  takahashi_min_buses::Int = 50
end

"""
    MatpowerImportConfig

Typed MATPOWER import/example configuration for case selection and import
conventions.
"""
Base.@kwdef struct MatpowerImportConfig
  pv_voltage_source::Symbol = :gen_vg
  pv_voltage_mismatch_tol_pu::Float64 = 1.0e-4
  compare_voltage_reference::Symbol = :imported_setpoint
  shift_unit::Symbol = :deg
  shift_sign::Float64 = 1.0
  ratio::Symbol = :normal
  enable_pq_gen_controllers::Bool = true
  apply_bus_names::Bool = false
  apply_branch_names::Bool = false
  apply_branch_kind::Bool = false
  import_for001_contingencies::Bool = true
  matpower_dcline_mode::Symbol = :pf_injections
end

"""
    ModelConfig

Typed model-construction configuration, the `model:` block: how any imported
case becomes a network, independent of its format (design decision D6 of the
adapter task; the keys lived in `matpower_import` and `transformer` before,
version-0 files are rewritten through the alias table).

`bus_shunt_model` selects how bus shunts enter the admittance model.

`tap_changer_model` selects how transformer tap changers act on the
equivalent circuit:
- `:ideal` — tap steps only change the complex winding ratio; the series
  impedance keeps its neutral-position value (no impedance feedback;
  previous Sparlectra behavior).
- `:impedance_correction` — tap steps additionally re-refer the transformer
  series impedance through the tapped winding (R and X scaled with
  `|1 + f·e^(jφ)|²`, implemented centrally in `src/equicircuit.jl`).
The option applies to all transformers of an imported case and is read by
every importer.

`auto_profile` (`:off`, `:recommend`, `:apply`) drives the read-only import
analysis with convention recommendations; `auto_profile_log` writes its
artifact. `net_cache_enabled` gates the binary net cache (active only while
`auto_profile` is `:off`). `preallocate_network`/`preallocate_min_buses`
control container preallocation for large cases.
"""
Base.@kwdef struct ModelConfig
  bus_shunt_model::Symbol = :admittance
  tap_changer_model::Symbol = :ideal
  auto_profile::Symbol = :off
  auto_profile_log::Bool = true
  net_cache_enabled::Bool = false
  preallocate_network::Symbol = :auto
  preallocate_min_buses::Int = 1000
end

"""
    MatpowerExportConfig

Typed MATPOWER export configuration.

`write_solution` selects whether [`writeMatpowerCasefile`](@ref) writes the
solved AC power-flow state back into the exported case:
- `true` (default): `mpc.bus` `VM`/`VA` reflect the solved node state and
  `mpc.branch` gains the standard MATPOWER result columns 14–17
  (`PF`, `QF`, `PT`, `QT`), sourced from the existing branch-flow report path.
  A `mpc.sparlectra.solution_written = 1` marker documents this. If the network
  has not been solved, the exporter warns and falls back to the 13-column
  model-only export instead of writing empty result columns.
- `false`: the export is a pure model file. `mpc.branch` keeps its historical
  13 columns, and `VM = 1.0`/`VA = 0.0` for all non-slack/non-PV buses (slack
  and PV setpoints are preserved).
"""
Base.@kwdef struct MatpowerExportConfig
  write_solution::Bool = true
end

"""
    PerformanceConfig

Typed performance and diagnostic-volume configuration.
"""
Base.@kwdef struct PerformanceConfig
  enabled::Bool = false
  level::Symbol = :summary
  print_to_console::Bool = true
  write_to_logfile::Bool = true
  show_allocations::Bool = false
  show_iteration_table::Bool = false
  compact_logging::Bool = true
  representative_warmup_runs::Int = 0
  compare_cold_warm::Bool = false
  skip_reference_comparison::Bool = false
  skip_expensive_diagnostics::Bool = true
  skip_branch_neighborhood_report::Bool = true
  max_diagnostic_rows::Int = 25
end

"""
The `benchmark` configuration section: repeated-run measurement of a case
(`enabled`, `samples`, `seconds`).
"""
Base.@kwdef struct BenchmarkConfig
  enabled::Bool = true
  methods::Vector{Symbol} = [:rectangular]
  seconds::Float64 = 2.0
  samples::Int = 50
  show_once::Bool = false
  show_once_output::Symbol = :classic
  show_once_max_nodes::Int = 0
end

"""
    ParallelRuntimeConfig

Typed configuration of the in-process parallel execution of independent
work items (island solves, short-circuit sweeps, contingency batches).
`enabled = false` forces every parallel site onto the serial path (the
serial functions themselves, not copies). `max_tasks = "auto"` resolves to
`Threads.nthreads()`; an integer string caps the task count (applied via
chunking, so it caps `Threads.@threads` sites too). Work lists shorter than
`min_work_items` run serially to avoid task overhead on tiny cases.
"""
Base.@kwdef struct ParallelRuntimeConfig
  enabled::Bool = true
  max_tasks::String = "auto"
  min_work_items::Int = 4
end

"""
    parallel_max_tasks(cfg::ParallelRuntimeConfig) -> Int

Resolve the configured `runtime.parallel.max_tasks` to a concrete task
count: `"auto"` yields `Threads.nthreads()`, an integer string yields that
number (validation guarantees it is positive).
"""
function parallel_max_tasks(cfg::ParallelRuntimeConfig)::Int
  value = lowercase(strip(cfg.max_tasks))
  value == "auto" && return Threads.nthreads()
  return parse(Int, value)
end

"""
    RuntimeConfig

Typed runtime/threading configuration for example and benchmark entry points.
"""
Base.@kwdef struct RuntimeConfig
  julia_threads::String = "keep"
  blas_threads::String = "keep"
  print_thread_config::Bool = true
  # the configured case selection (runtime.case / runtime.cases); lived in
  # matpower_import before the version-1 layout, aliased for old files
  case::String = ""
  cases::Vector{String} = String[]
  casefile::String = ""
  case_name::String = ""
  case_source::String = ""
  configured_default_casefile::String = ""
  parallel::ParallelRuntimeConfig = ParallelRuntimeConfig()
end

"""
    DiagnosticsConfig

Typed diagnostic-output configuration shared by examples and future modules.
"""
# Console/logfile output options are owned by `output.*` alone. The
# `diagnostics` block used to carry duplicates of six of them; none of those
# was ever read (every consumer resolves `config.output.*`), so setting them
# silently did nothing. They are deprecated: the parser warns and ignores
# them (see `DiagnosticsConfig(raw)`).
"""
The `diagnostics` configuration section (deprecated keys live on in
`output`); retained so old files keep loading.
"""
Base.@kwdef struct DiagnosticsConfig
  log_effective_config::Bool = false
end

"""
    OutputConfig

Typed output and logfile-format configuration. Console and logfile result
streams are intentionally independent so example runners can disable classic
logfile tables without suppressing compact console progress.
"""
Base.@kwdef struct OutputConfig
  console_summary::Bool = true
  # Mirror the captured run output (everything that goes into run.log) live
  # to the real console while an API/service run executes. The archive stays
  # identical; default off keeps API runs quiet like before.
  console_live::Bool = false
  console_auto_profile::Symbol = :compact
  console_diagnostics::Symbol = :compact
  console_q_limit_events::Symbol = :summary
  console_max_rows::Int = 100
  logfile_results::Symbol = :off
  result_table_max_rows::Int = 200
  result_table_large_case_threshold_buses::Int = 1000
  result_table_large_case_mode::Symbol = :summary
  detailed_result_csv_write_mode::Symbol = :auto
  detailed_result_csv_exporter::Symbol = :auto
  detailed_result_csv_direct_threshold_buses::Int = 10_000
  detailed_result_csv_buffer_initial_bytes::Int = 8 * 1024 * 1024
  detailed_result_csv_buffer_max_bytes::Int = 64 * 1024 * 1024
  detailed_result_csv_streaming_threshold_rows::Int = 100_000
  logfile_diagnostics::Symbol = :compact
  logfile_performance::Symbol = :compact
  logfile_warnings::Symbol = :table
  # one info line per process when running without a sysimage or compiled
  # executable (first run pays JIT warm-up); false silences it
  startup_latency_hint::Bool = true
end

"""
    WebUIConfig

Web UI presentation preferences that are not part of the solver/API contract.
"""
Base.@kwdef struct WebUIConfig
  show_case_settings_notice::Bool = true
  # How far the operation log reaches back. Every Web UI start drops entries
  # older than this; 0 keeps only the current session. The environment
  # variable SPARLECTRA_WEBUI_OPERATION_LOG_RETENTION_DAYS still wins, for
  # headless runs that never read this file.
  operation_log_retention_days::Int = 10
end

"""
    SparlectraConfig

Central typed configuration assembled once at application or example boundaries.
Module-specific sections own their parsing and validation through constructors
such as `PowerFlowConfig(raw)` and `MatpowerImportConfig(raw)`.
"""
Base.@kwdef struct SparlectraConfig
  powerflow::PowerFlowConfig = PowerFlowConfig()
  state_estimation::StateEstimationConfig = StateEstimationConfig()
  matpower::MatpowerImportConfig = MatpowerImportConfig()
  cgmes::CGMESImportConfig = CGMESImportConfig()
  shortcircuit::ShortCircuitConfig = ShortCircuitConfig()
  matpower_export::MatpowerExportConfig = MatpowerExportConfig()
  model::ModelConfig = ModelConfig()
  performance::PerformanceConfig = PerformanceConfig()
  benchmark::BenchmarkConfig = BenchmarkConfig()
  contingency::ContingencyConfig = ContingencyConfig()
  runtime::RuntimeConfig = RuntimeConfig()
  diagnostics::DiagnosticsConfig = DiagnosticsConfig()
  output::OutputConfig = OutputConfig()
  control::ControlConfig = ControlConfig()
  webui::WebUIConfig = WebUIConfig()
  # dotted keys the USER set explicitly (user YAML file, cli_overrides, or
  # API overrides), filled by load_sparlectra_config. The auto power-flow
  # mode only fills keys NOT in this set: explicit user choices always win.
  # Empty for configs built directly from keywords (conservative: auto may
  # then fill everything).
  user_set_keys::Set{String} = Set{String}()
end

const _SPARLECTRA_CONFIG_CACHE = Ref{Union{Nothing,NamedTuple}}(nothing)
const ACTIVE_SPARLECTRA_CONFIG = Ref{SparlectraConfig}(SparlectraConfig())
const DEFAULT_SPARLECTRA_CONFIG_PATH = normpath(joinpath(@__DIR__, "configuration.yaml.example"))
const USER_SPARLECTRA_CONFIG_PATH = normpath(joinpath(@__DIR__, "..", "..", "examples", "configuration.yaml"))

const SUPPORTED_POWERFLOW_METHOD = :rectangular
const POWERFLOW_MODE_VALUES = (:manual, :auto)
const SHORT_CIRCUIT_SWEEP_METHOD_VALUES = (:auto, :solves, :takahashi)
const POWERFLOW_START_ANGLE_MODE_VALUES = (:classic, :dc, :bus_va_blend, :matpower_va)
const POWERFLOW_START_VOLTAGE_MODE_VALUES = (:classic, :pv_gen_vg, :pv_bus_vm, :all_bus_vm, :profile_blend)
const POWERFLOW_START_PROFILE_SOURCE_VALUES = (:flat, :dc, :bus_metadata, :historical_profile, :matpower_reference, :state_estimation, :se_snapshot, :scada_snapshot)
# CGMES-only start-value selection: :flat = synthetic flat start (the solver
# earns the solution itself), :sv = start from the delivery's imported
# SvVoltage state. Wins over power_flow(.start_mode).flatstart on CGMES runs.
const CGMES_START_VALUES_VALUES = (:auto, :flat, :sv)
# placeholder-guard behavior: warn_skip keeps completeness-set filler values
# out of the solve with a warning; strict aborts the import instead — for
# deliveries where silently dropping data is not acceptable
const CGMES_PLACEHOLDER_GUARDS_VALUES = (:warn_skip, :strict)
const CGMES_HVDC_MODE_VALUES = (:injections, :paired_control)
const TRUST_REGION_STEP_MODE_VALUES = (:scaled, :dogleg)
const QLIMIT_START_MODE_VALUES = (:iteration, :auto, :iteration_or_auto)
const QLIMIT_ENFORCEMENT_MODE_VALUES = (:active_set, :classic_simultaneous, :classic_one_at_a_time)
const QLIMIT_ENFORCEMENT_MODE_LEGACY_ALIASES = Dict(
  :matpower_simultaneous => :classic_simultaneous,
  :matpower_one_at_a_time => :classic_one_at_a_time,
)
const QLIMIT_GUARD_ZERO_RANGE_MODE_VALUES = (:lock_pq,)
const QLIMIT_GUARD_NARROW_RANGE_MODE_VALUES = (:prefer_pq, :lock_pq)
const QLIMIT_GUARD_VIOLATION_MODE_VALUES = (:delayed_switch, :lock_pq)
const POWERFLOW_ISLAND_MODE_VALUES = (:solve_independent, :solve_parallel)
const POWERFLOW_ISLAND_REFERENCE_POLICY_VALUES = (:matpower_like,)
const RECTANGULAR_PREALLOCATE_WORKSPACE_VALUES = (:off, :on, :auto)
const STATE_ESTIMATION_METHOD_VALUES = (:wls,)
const STATE_ESTIMATION_ROBUST_MODE_VALUES = (:off, :staged, :replacement)
const STATE_ESTIMATION_PMU_REF_OFFSET_VALUES = (:auto, :off)
const MATPOWER_PV_VOLTAGE_SOURCE_VALUES = (:gen_vg, :bus_vm, :auto, :strict_check)
const MATPOWER_COMPARE_VOLTAGE_REFERENCE_VALUES = (:bus_vm, :gen_vg, :imported_setpoint, :hybrid)
const MATPOWER_SHIFT_UNIT_VALUES = (:deg, :rad)
const MATPOWER_RATIO_VALUES = (:normal, :reciprocal)
const MATPOWER_BUS_SHUNT_MODEL_VALUES = (:admittance, :voltage_dependent_injection)
const MATPOWER_AUTO_PROFILE_VALUES = (:off, :recommend, :apply)
const MATPOWER_DCLINE_MODE_VALUES = (:reject_active, :ignore_inactive, :pf_injections, :paired_control)

# `reject_active` was the shipped template default of early releases, so
# stored user/Web UI configuration files still carry it without it ever
# having been a deliberate choice. Loading it verbatim would abort every
# DC-line case (e.g. case_SyntheticUSA) for exactly those users — the
# configuration surface therefore normalizes it to the current default with
# a warning. The strict fail-fast behavior stays available programmatically:
# createNetFromMatPowerCase(matpower_dcline_mode = :reject_active).
function _normalize_dcline_mode(requested::Symbol)::Symbol
  requested === :reject_active || return requested
  @warn "matpower_import.matpower_dcline_mode: reject_active is deprecated in configuration files and loads as pf_injections (active mpc.dcline rows import as per-terminal injection pairs). Remove the key or set pf_injections to silence this warning."
  return :pf_injections
end
const TRANSFORMER_TAP_CHANGER_MODEL_VALUES = (:ideal, :impedance_correction)
const PERFORMANCE_LEVEL_VALUES = (:off, :summary, :iteration, :full)
const OUTPUT_CONSOLE_AUTO_PROFILE_VALUES = (:off, :compact, :full)
const OUTPUT_CONSOLE_DIAGNOSTICS_VALUES = (:off, :compact, :summary, :full)
const OUTPUT_CONSOLE_Q_LIMIT_EVENTS_VALUES = (:off, :summary, :full)
const OUTPUT_LOGFILE_RESULTS_VALUES = (:off, :compact, :classic, :full)
const OUTPUT_RESULT_TABLE_LARGE_CASE_MODE_VALUES = (:summary, :classic, :full)
const OUTPUT_DETAILED_RESULT_CSV_WRITE_MODE_VALUES = (:auto, :buffered, :streaming)
const OUTPUT_DETAILED_RESULT_CSV_EXPORTER_VALUES = (:auto, :report, :direct)
const OUTPUT_LOGFILE_DIAGNOSTICS_VALUES = (:off, :compact, :full)
const OUTPUT_LOGFILE_PERFORMANCE_VALUES = (:off, :compact, :full)
const OUTPUT_LOGFILE_WARNINGS_VALUES = (:off, :summary, :table, :full)
const BENCHMARK_SHOW_ONCE_OUTPUT_VALUES = (:classic, :dataframe, :compact)

_powerflow_method_label(method)::Symbol = Symbol(method) === :polar ? :polar_full : Symbol(method)

function unsupported_powerflow_method_message(method)::String
  return "Only the rectangular AC power-flow solver is supported in this version.\nRequested method: $(_powerflow_method_label(method))."
end

_as_symbol_vector_cfg(x)::Vector{Symbol} = x isa AbstractVector ? Symbol[_as_symbol_cfg(item) for item in x] : [_as_symbol_cfg(x)]

function _validate_rectangular_powerflow_options(; method::Symbol = SUPPORTED_POWERFLOW_METHOD, sparse::Bool = true)
  requested = _powerflow_method_label(method)
  requested === SUPPORTED_POWERFLOW_METHOD || throw(ArgumentError(unsupported_powerflow_method_message(requested)))
  sparse || throw(ArgumentError("Only sparse power-flow matrices are supported in this version. Requested sparse=false."))
  return nothing
end

function _validate_powerflow_solver_config!(raw::AbstractDict)
  pf = _merged_section(raw, "powerflow")
  method = _as_symbol_cfg(_raw_get(pf, "method", SUPPORTED_POWERFLOW_METHOD))
  if haskey(pf, "sparse") || haskey(pf, :sparse) || haskey(pf, "opt_sparse") || haskey(pf, :opt_sparse)
    throw(ArgumentError("power_flow.sparse/opt_sparse is obsolete and no longer configurable. Sparse matrices are mandatory."))
  end
  _validate_rectangular_powerflow_options(method = method, sparse = true)
  benchmark = _raw_section(raw, "benchmark")
  if haskey(benchmark, "methods")
    for m in _as_symbol_vector_cfg(benchmark["methods"])
      requested = _powerflow_method_label(m)
      requested === SUPPORTED_POWERFLOW_METHOD || throw(ArgumentError(unsupported_powerflow_method_message(requested)))
    end
  end
  return nothing
end

"""
    set_sparlectra_config!(cfg) -> SparlectraConfig

Install `cfg` as the active configuration of the process, replacing the one
every `*_config()` accessor returns.
"""
function set_sparlectra_config!(cfg::SparlectraConfig)
  ACTIVE_SPARLECTRA_CONFIG[] = cfg
  return cfg
end

"""
    _resolve_parallel_runtime(enabled, max_tasks, min_work_items) -> (on, cap, min_items)

Resolve per-call parallel overrides against the ACTIVE
`runtime.parallel` configuration (`nothing` = take the configured value).
Meant to be called ONCE by the orchestrating (serial) code of a parallel
site before any task spawns; workers never touch the config globals.
"""
function _resolve_parallel_runtime(enabled::Union{Nothing,Bool}, max_tasks::Union{Nothing,Int}, min_work_items::Union{Nothing,Int})
  pcfg = runtime_config().parallel
  on = something(enabled, pcfg.enabled)
  cap = max_tasks === nothing ? parallel_max_tasks(pcfg) : max_tasks
  min_items = something(min_work_items, pcfg.min_work_items)
  return (on, cap, min_items)
end

"""
    active_sparlectra_config() -> SparlectraConfig

The configuration the process currently runs with, as set by
[`set_sparlectra_config!`](@ref) or the last [`load_sparlectra_config!`](@ref).
"""
active_sparlectra_config()::SparlectraConfig = ACTIVE_SPARLECTRA_CONFIG[]
"""
    powerflow_config() -> PowerFlowConfig

The `power_flow` section of the active configuration.
"""
powerflow_config()::PowerFlowConfig = active_sparlectra_config().powerflow
"""
    matpower_import_config() -> MatpowerImportConfig

The `matpower_import` section of the active configuration.
"""
matpower_import_config()::MatpowerImportConfig = active_sparlectra_config().matpower
"""
    matpower_export_config() -> MatpowerExportConfig

The `matpower_export` section of the active configuration.
"""
matpower_export_config()::MatpowerExportConfig = active_sparlectra_config().matpower_export
"""
    model_config() -> ModelConfig

The `model` section of the active configuration: how an imported case
becomes a network (bus shunt model, tap changer model, auto profile,
net cache, preallocation).
"""
model_config()::ModelConfig = active_sparlectra_config().model
"""
    state_estimation_config() -> StateEstimationConfig

The `state_estimation` section of the active configuration.
"""
state_estimation_config()::StateEstimationConfig = active_sparlectra_config().state_estimation
"""
    contingency_config() -> ContingencyConfig

The `contingency` section of the active configuration.
"""
contingency_config()::ContingencyConfig = active_sparlectra_config().contingency
"""
    diagnostics_config() -> DiagnosticsConfig

The `diagnostics` section of the active configuration.
"""
diagnostics_config()::DiagnosticsConfig = active_sparlectra_config().diagnostics
"""
    output_config() -> OutputConfig

The `output` section of the active configuration.
"""
output_config()::OutputConfig = active_sparlectra_config().output
"""
    performance_config() -> PerformanceConfig

The `performance` section of the active configuration.
"""
performance_config()::PerformanceConfig = active_sparlectra_config().performance
"""
    benchmark_config() -> BenchmarkConfig

The `benchmark` section of the active configuration.
"""
benchmark_config()::BenchmarkConfig = active_sparlectra_config().benchmark
"""
    runtime_config() -> RuntimeConfig

The `runtime` section of the active configuration.
"""
runtime_config()::RuntimeConfig = active_sparlectra_config().runtime
control_config()::ControlConfig = active_sparlectra_config().control

"""
    configuration_path_from_inputs(; env_var, fallback_paths) -> String

Resolve which configuration file a process should load: the environment
variable first, then the first existing fallback path, then the packaged
default template.
"""
function configuration_path_from_inputs(; env_var::AbstractString = "SPARLECTRA_CONFIGURATION_YAML", fallback_paths::AbstractVector{<:AbstractString} = String[])
  candidate = strip(get(ENV, env_var, ""))
  !isempty(candidate) && return candidate
  for path in fallback_paths
    file = strip(String(path))
    isempty(file) && continue
    isfile(file) && return file
  end
  return ""
end

_config_key(key::Symbol) = String(key)
_config_key(key) = String(key)

function _raw_get(raw::AbstractDict, key::AbstractString, default)
  haskey(raw, key) && return raw[key]
  symkey = Symbol(key)
  haskey(raw, symkey) && return raw[symkey]
  return default
end

function _raw_section(raw::AbstractDict, key::AbstractString, aliases::AbstractVector{<:AbstractString} = String[])
  value = _raw_get(raw, key, nothing)
  if value === nothing
    for alias in aliases
      value = _raw_get(raw, alias, nothing)
      value === nothing || break
    end
  end
  value === nothing && (value = Dict{String,Any}())
  value isa AbstractDict && return value
  throw(ArgumentError("Configuration section $(repr(key)) must be a dictionary."))
end

_as_string_cfg(x)::String = x isa AbstractString ? String(x) : string(x)
function _as_string_vector_cfg(x)::Vector{String}
  x isa AbstractVector || throw(ArgumentError("matpower_import.cases must be a vector of non-empty case names."))
  values = String[strip(_as_string_cfg(item)) for item in x]
  any(isempty, values) && throw(ArgumentError("matpower_import.cases must not contain empty case names."))
  return values
end
_as_symbol_cfg(x)::Symbol = x isa Symbol ? x : (x isa Bool ? (x ? :on : :off) : Symbol(lowercase(String(x))))
function _as_auto_profile_symbol_cfg(x)::Symbol
  x === false && return :off
  x === true && return :apply
  sym = _as_symbol_cfg(x)
  sym in (:false, :no, :none) && return :off
  sym in (:true, :yes) && return :apply
  return sym
end
_as_bool_cfg(x)::Bool = as_bool(x)
_as_int_cfg(x)::Int = x isa Integer ? Int(x) : parse(Int, String(x))
_as_float_cfg(x)::Float64 = x isa Real ? Float64(x) : parse(Float64, String(x))

function _as_float_vector_cfg(x)::Vector{Float64}
  x isa AbstractVector || throw(ArgumentError("Expected a vector of floating-point values; got $(repr(x))."))
  return Float64[_as_float_cfg(item) for item in x]
end

function _as_int_vector_cfg(x)::Vector{Int}
  return as_int_vector(x)
end

function _validate_nonnegative(name::AbstractString, value::Real)
  value >= 0 || throw(ArgumentError("$(name) must be non-negative; got $(value)."))
  return value
end
function _validate_allowed_symbol(name::AbstractString, value::Symbol, allowed::Tuple)
  if !(value in allowed)
    guidance = name == "power_flow.start_mode.voltage_mode" && value === :bus_vm_va_blend ?
      " The former bus_vm_va_blend alias has been removed; use voltage_mode: profile_blend with profile_source: matpower_reference." : ""
    allowed_text = join(string.(allowed), ", ")
    throw(ArgumentError("Invalid value:\n$(name) = \"$(value)\"\n\nAllowed values:\n$(allowed_text)\n\n$(name) must be one of $(collect(allowed)); got $(value).$(guidance)"))
  end
  return value
end

function _validate_allowed_symbol(name::AbstractString, value::Symbol, allowed::AbstractVector{Symbol})
  if !(value in allowed)
    allowed_text = join(string.(allowed), ", ")
    throw(ArgumentError("Invalid value:\n$(name) = \"$(value)\"\n\nAllowed values:\n$(allowed_text)\n\n$(name) must be one of $(allowed); got $(value)."))
  end
  return value
end

function _canonical_qlimit_enforcement_mode(value::Symbol)::Symbol
  canonical = get(QLIMIT_ENFORCEMENT_MODE_LEGACY_ALIASES, value, value)
  if canonical in QLIMIT_ENFORCEMENT_MODE_VALUES
    return canonical
  end
  legacy_aliases = ["$(alias) -> $(target)" for (alias, target) in sort(collect(QLIMIT_ENFORCEMENT_MODE_LEGACY_ALIASES); by = first)]
  throw(ArgumentError("power_flow.qlimits.enforcement_mode must be one of $(collect(QLIMIT_ENFORCEMENT_MODE_VALUES)); got $(value). Legacy aliases accepted: $(legacy_aliases)."))
end

## power_flow.tol_MW (task_tol_watts part B): the tolerance in physical
## units. Zero or negative is an error naming the key; an implausibly
## large bound is a WARNING, not an error, because "coarse on purpose" is
## a legitimate choice and only the user knows the network's size.
function _validate_tol_mw(raw)
  raw === nothing && return nothing
  value = _as_float_cfg(raw)
  value > 0.0 || throw(ArgumentError("power_flow.tol_MW must be > 0 (got $(value)); it is a convergence bound in MW, and zero or negative can never be reached."))
  value >= 1.0 && @warn "power_flow.tol_MW = $(value) MW is a very coarse convergence bound; a solution accepted at this mismatch may be far from solved. The per-unit equivalent depends on the case base (tol = tol_MW / baseMVA)."
  return value
end

function _validate_positive(name::AbstractString, value::Real)
  value > 0 || throw(ArgumentError("$(name) must be positive; got $(value)."))
  return value
end

function StartModeConfig(raw::AbstractDict)
  angle_mode = _validate_allowed_symbol("power_flow.start_mode.angle_mode", _as_symbol_cfg(_raw_get(raw, "angle_mode", :dc)), POWERFLOW_START_ANGLE_MODE_VALUES)
  voltage_mode = _validate_allowed_symbol("power_flow.start_mode.voltage_mode", _as_symbol_cfg(_raw_get(raw, "voltage_mode", :profile_blend)), POWERFLOW_START_VOLTAGE_MODE_VALUES)
  profile_source = _validate_allowed_symbol("power_flow.start_mode.profile_source", _as_symbol_cfg(_raw_get(raw, "profile_source", :matpower_reference)), POWERFLOW_START_PROFILE_SOURCE_VALUES)
  return StartModeConfig(
    flatstart = _as_bool_cfg(_raw_get(raw, "flatstart", _raw_get(raw, "opt_flatstart", false))),
    angle_mode = angle_mode,
    voltage_mode = voltage_mode,
    profile_source = profile_source,
    start_projection = _as_bool_cfg(_raw_get(raw, "start_projection", false)),
    try_dc_start = _as_bool_cfg(_raw_get(raw, "try_dc_start", _raw_get(raw, "start_projection_try_dc_start", true))),
    try_blend_scan = _as_bool_cfg(_raw_get(raw, "try_blend_scan", _raw_get(raw, "start_projection_try_blend_scan", true))),
    branch_guard = _as_bool_cfg(_raw_get(raw, "branch_guard", _raw_get(raw, "start_projection_branch_guard", true))),
    measure_candidates = _as_bool_cfg(_raw_get(raw, "measure_candidates", _raw_get(raw, "start_projection_measure_candidates", true))),
    accept_unmeasured_dc_start = _as_bool_cfg(_raw_get(raw, "accept_unmeasured_dc_start", _raw_get(raw, "start_projection_accept_unmeasured_dc_start", false))),
    dc_seed_unconditional = _as_bool_cfg(_raw_get(raw, "dc_seed_unconditional", false)),
    reuse_import_data = _as_bool_cfg(_raw_get(raw, "reuse_import_data", _raw_get(raw, "start_projection_reuse_import_data", true))),
    blend_lambdas = _as_float_vector_cfg(_raw_get(raw, "blend_lambdas", _raw_get(raw, "start_projection_blend_lambdas", [0.25, 0.5, 0.75]))),
    dc_angle_limit_deg = _validate_positive("start_projection_dc_angle_limit_deg", _as_float_cfg(_raw_get(raw, "dc_angle_limit_deg", _raw_get(raw, "start_projection_dc_angle_limit_deg", 60.0)))),
  )
end

function StartCurrentIterationConfig(raw::AbstractDict)
  max_iter = _as_int_cfg(_raw_get(raw, "max_iter", 10))
  max_iter >= 0 || throw(ArgumentError("power_flow.start_current_iteration.max_iter must be non-negative; got $(max_iter)."))
  damping = _as_float_cfg(_raw_get(raw, "damping", 0.5))
  isfinite(damping) && 0.0 < damping <= 1.0 || throw(ArgumentError("power_flow.start_current_iteration.damping must be finite and in (0, 1]; got $(damping)."))
  min_improvement_factor = _as_float_cfg(_raw_get(raw, "min_improvement_factor", 0.98))
  isfinite(min_improvement_factor) && min_improvement_factor >= 0.0 || throw(ArgumentError("power_flow.start_current_iteration.min_improvement_factor must be finite and non-negative; got $(min_improvement_factor)."))
  vm_min_pu = _validate_positive("power_flow.start_current_iteration.vm_min_pu", _as_float_cfg(_raw_get(raw, "vm_min_pu", 0.5)))
  vm_max_pu = _validate_positive("power_flow.start_current_iteration.vm_max_pu", _as_float_cfg(_raw_get(raw, "vm_max_pu", 1.5)))
  vm_min_pu <= vm_max_pu || throw(ArgumentError("power_flow.start_current_iteration.vm_min_pu must be <= vm_max_pu."))
  return StartCurrentIterationConfig(
    enabled = _as_bool_cfg(_raw_get(raw, "enabled", false)),
    max_iter = max_iter,
    tol = _validate_positive("power_flow.start_current_iteration.tol", _as_float_cfg(_raw_get(raw, "tol", 1.0e-3))),
    damping = damping,
    accept_only_if_improved = _as_bool_cfg(_raw_get(raw, "accept_only_if_improved", true)),
    min_improvement_factor = min_improvement_factor,
    vm_min_pu = vm_min_pu,
    vm_max_pu = vm_max_pu,
    max_angle_step_deg = _validate_positive("power_flow.start_current_iteration.max_angle_step_deg", _as_float_cfg(_raw_get(raw, "max_angle_step_deg", 30.0))),
    only_for_large_cases = _as_bool_cfg(_raw_get(raw, "only_for_large_cases", false)),
  )
end

function MeritLineSearchConfig(raw::AbstractDict)
  armijo_c1 = _as_float_cfg(_raw_get(raw, "armijo_c1", 1.0e-4))
  isfinite(armijo_c1) && 0.0 < armijo_c1 < 0.5 || throw(ArgumentError("power_flow.merit.armijo_c1 must be finite and in (0, 0.5); got $(armijo_c1)."))
  return MeritLineSearchConfig(
    enabled = _as_bool_cfg(_raw_get(raw, "enabled", false)),
    armijo_c1 = armijo_c1,
    scale_p = _validate_positive("power_flow.merit.scale_p", _as_float_cfg(_raw_get(raw, "scale_p", 1.0))),
    scale_q = _validate_positive("power_flow.merit.scale_q", _as_float_cfg(_raw_get(raw, "scale_q", 1.0))),
    scale_v = _validate_positive("power_flow.merit.scale_v", _as_float_cfg(_raw_get(raw, "scale_v", 1.0))),
    fallback_max_mismatch = _as_bool_cfg(_raw_get(raw, "fallback_max_mismatch", true)),
  )
end

function TrustRegionConfig(raw::AbstractDict)
  initial_radius = _validate_positive("power_flow.trust_region.initial_radius", _as_float_cfg(_raw_get(raw, "initial_radius", 1.0)))
  min_radius = _validate_positive("power_flow.trust_region.min_radius", _as_float_cfg(_raw_get(raw, "min_radius", 1e-4)))
  max_radius = _validate_positive("power_flow.trust_region.max_radius", _as_float_cfg(_raw_get(raw, "max_radius", 10.0)))
  min_radius < initial_radius <= max_radius || throw(ArgumentError("power_flow.trust_region must satisfy min_radius < initial_radius <= max_radius; got min_radius=$(min_radius), initial_radius=$(initial_radius), max_radius=$(max_radius)."))
  eta_accept = _validate_positive("power_flow.trust_region.eta_accept", _as_float_cfg(_raw_get(raw, "eta_accept", 0.1)))
  shrink_factor = _as_float_cfg(_raw_get(raw, "shrink_factor", 0.5))
  isfinite(shrink_factor) && 0.0 < shrink_factor < 1.0 || throw(ArgumentError("power_flow.trust_region.shrink_factor must be finite and in (0, 1); got $(shrink_factor)."))
  expand_factor = _as_float_cfg(_raw_get(raw, "expand_factor", 2.0))
  isfinite(expand_factor) && expand_factor > 1.0 || throw(ArgumentError("power_flow.trust_region.expand_factor must be finite and > 1; got $(expand_factor)."))
  expand_threshold = _as_float_cfg(_raw_get(raw, "expand_threshold", 0.75))
  isfinite(expand_threshold) && 0.0 < expand_threshold < 1.0 || throw(ArgumentError("power_flow.trust_region.expand_threshold must be finite and in (0, 1); got $(expand_threshold)."))
  step_mode = _validate_allowed_symbol("power_flow.trust_region.step_mode", _as_symbol_cfg(_raw_get(raw, "step_mode", :scaled)), TRUST_REGION_STEP_MODE_VALUES)
  return TrustRegionConfig(
    enabled = _as_bool_cfg(_raw_get(raw, "enabled", false)),
    initial_radius = initial_radius,
    min_radius = min_radius,
    max_radius = max_radius,
    eta_accept = eta_accept,
    shrink_factor = shrink_factor,
    expand_factor = expand_factor,
    expand_threshold = expand_threshold,
    step_mode = step_mode,
  )
end

function ApslfConfig(raw::AbstractDict)
  order = _as_int_cfg(_raw_get(raw, "order", 40))
  order >= 1 || throw(ArgumentError("power_flow.apslf.order must be >= 1; got $(order)."))
  return ApslfConfig(
    order = order,
    use_pade = _as_bool_cfg(_raw_get(raw, "use_pade", true)),
    nr_polish = _as_bool_cfg(_raw_get(raw, "nr_polish", true)),
  )
end

function ApslfStartConfig(raw::AbstractDict)
  order = _as_int_cfg(_raw_get(raw, "order", 40))
  order >= 1 || throw(ArgumentError("power_flow.apslf_start.order must be >= 1; got $(order)."))
  return ApslfStartConfig(
    enabled = _as_bool_cfg(_raw_get(raw, "enabled", false)),
    order = order,
  )
end

function QLimitConfig(raw::AbstractDict)
  qlimits_enabled = _raw_get(raw, "enabled", true)
  guard_raw = _raw_get(raw, "guard", Dict{String,Any}())
  guard_enabled_default = guard_raw isa AbstractDict ? _as_bool_cfg(_raw_get(guard_raw, "enabled", false)) : _as_bool_cfg(guard_raw)
  guard_cfg = guard_raw isa AbstractDict ? guard_raw : Dict{String,Any}()
  merged = merge(Dict{Any,Any}(raw), Dict{Any,Any}(guard_cfg))
  return QLimitConfig(
    start_iter = _as_int_cfg(_raw_get(merged, "start_iter", _raw_get(merged, "qlimit_start_iter", 2))),
    start_mode = _validate_allowed_symbol(
      "power_flow.qlimits.start_mode",
      _as_symbol_cfg((_raw_get(merged, "start_mode", nothing) isa AbstractDict) ? _raw_get(merged, "qlimit_start_mode", :iteration) : _raw_get(merged, "start_mode", _raw_get(merged, "qlimit_start_mode", :iteration))),
      QLIMIT_START_MODE_VALUES,
    ),
    auto_q_delta_pu = _validate_nonnegative("qlimit_auto_q_delta_pu", _as_float_cfg(_raw_get(merged, "auto_q_delta_pu", _raw_get(merged, "qlimit_auto_q_delta_pu", 1e-4)))),
    hysteresis_pu = _validate_nonnegative("power_flow.qlimits.hysteresis_pu", _as_float_cfg(_raw_get(merged, "hysteresis_pu", _raw_get(merged, "q_hyst_pu", 0.01)))),
    cooldown_iters = _as_int_cfg(_raw_get(merged, "cooldown_iters", 1)),
    guard = _as_bool_cfg(_raw_get(raw, "qlimit_guard", guard_enabled_default)),
    guard_min_q_range_pu = _validate_nonnegative("qlimit_guard_min_q_range_pu", _as_float_cfg(_raw_get(merged, "min_q_range_pu", _raw_get(merged, "guard_min_q_range_pu", _raw_get(merged, "qlimit_guard_min_q_range_pu", 1e-4))))),
    guard_zero_range_mode = _validate_allowed_symbol("power_flow.qlimits.guard.zero_range_mode", _as_symbol_cfg(_raw_get(merged, "zero_range_mode", _raw_get(merged, "guard_zero_range_mode", _raw_get(merged, "qlimit_guard_zero_range_mode", :lock_pq)))), QLIMIT_GUARD_ZERO_RANGE_MODE_VALUES),
    guard_narrow_range_mode = _validate_allowed_symbol(
      "power_flow.qlimits.guard.narrow_range_mode",
      _as_symbol_cfg(_raw_get(merged, "narrow_range_mode", _raw_get(merged, "guard_narrow_range_mode", _raw_get(merged, "qlimit_guard_narrow_range_mode", :prefer_pq)))),
      QLIMIT_GUARD_NARROW_RANGE_MODE_VALUES,
    ),
    guard_max_switches = _as_int_cfg(_raw_get(merged, "max_switches", _raw_get(merged, "guard_max_switches", _raw_get(merged, "qlimit_guard_max_switches", 10)))),
    guard_freeze_after_repeated_switching = _as_bool_cfg(_raw_get(merged, "freeze_after_repeated_switching", _raw_get(merged, "guard_freeze_after_repeated_switching", _raw_get(merged, "qlimit_guard_freeze_after_repeated_switching", true)))),
    guard_accept_bounded_violations = _as_bool_cfg(_raw_get(merged, "accept_bounded_violations", _raw_get(merged, "guard_accept_bounded_violations", _raw_get(merged, "qlimit_guard_accept_bounded_violations", false)))),
    guard_max_remaining_violations = _as_int_cfg(_raw_get(merged, "max_remaining_violations", _raw_get(merged, "guard_max_remaining_violations", _raw_get(merged, "qlimit_guard_max_remaining_violations", 0)))),
    guard_violation_mode = _validate_allowed_symbol("power_flow.qlimits.guard.violation_mode", _as_symbol_cfg(_raw_get(merged, "violation_mode", _raw_get(merged, "guard_violation_mode", _raw_get(merged, "qlimit_guard_violation_mode", :delayed_switch)))), QLIMIT_GUARD_VIOLATION_MODE_VALUES),
    guard_violation_threshold_pu = _validate_nonnegative("qlimit_guard_violation_threshold_pu", _as_float_cfg(_raw_get(merged, "violation_threshold_pu", _raw_get(merged, "guard_violation_threshold_pu", _raw_get(merged, "qlimit_guard_violation_threshold_pu", 1e-4))))),
    guard_log = _as_bool_cfg(_raw_get(merged, "log", _raw_get(merged, "guard_log", _raw_get(merged, "qlimit_guard_log", true)))),
    trace_buses = _as_int_vector_cfg(_raw_get(merged, "trace_buses", _raw_get(merged, "qlimit_trace_buses", Int[]))),
    lock_pv_to_pq_buses = _as_int_vector_cfg(_raw_get(merged, "lock_pv_to_pq_buses", Int[])),
    ignore_q_limits = _as_bool_cfg(_raw_get(raw, "ignore_q_limits", qlimits_enabled == false)),
    enforcement_mode = _canonical_qlimit_enforcement_mode(_as_symbol_cfg(_raw_get(merged, "enforcement_mode", :active_set))),
  )
end

function _merged_section(raw::AbstractDict, section_name::AbstractString)
  aliases = section_name == "powerflow" ? ["power_flow"] : section_name == "matpower" ? ["matpower_import"] : section_name == "cgmes" ? ["cgmes_import"] : String[]
  canonical = section_name == "power_flow" ? "powerflow" : section_name == "matpower_import" ? "matpower" : section_name == "cgmes_import" ? "cgmes" : section_name
  primary = Dict{Any,Any}(_raw_section(raw, canonical))
  for alias in aliases
    merge!(primary, Dict{Any,Any}(_raw_section(raw, alias)))
  end
  return primary
end

function PowerFlowConfig(raw::AbstractDict)
  merged = _merged_section(raw, "powerflow")
  method = _as_symbol_cfg(_raw_get(merged, "method", SUPPORTED_POWERFLOW_METHOD))
  if haskey(merged, "sparse") || haskey(merged, :sparse) || haskey(merged, "opt_sparse") || haskey(merged, :opt_sparse)
    throw(ArgumentError("power_flow.sparse/opt_sparse is obsolete and no longer configurable. Sparse matrices are mandatory."))
  end
  _validate_rectangular_powerflow_options(method = method, sparse = true)
  start_raw = _raw_get(merged, "start_values", _raw_section(merged, "start_mode"))
  start_current_iteration_raw = _raw_section(merged, "start_current_iteration")
  merit_raw = _raw_section(merged, "merit")
  trust_region_raw = _raw_section(merged, "trust_region")
  qlimit_raw = _raw_section(merged, "qlimits")
  islands_raw = _raw_section(merged, "islands")
  apslf_raw = _raw_section(merged, "apslf")
  apslf_start_raw = _raw_section(merged, "apslf_start")
  dc_raw = _raw_section(merged, "dc")
  distributed_slack_raw = _raw_section(merged, "distributed_slack")
  external_grid_raw = _raw_section(merged, "external_grid")
  solver = _validate_allowed_symbol("power_flow.solver", _as_symbol_cfg(_raw_get(merged, "solver", :rectangular)), POWERFLOW_SOLVER_VALUES)
  apslf_cfg = ApslfConfig(apslf_raw)
  apslf_start_cfg = ApslfStartConfig(apslf_start_raw)
  if apslf_start_cfg.enabled && solver === :apslf
    throw(ArgumentError("power_flow.apslf_start.enabled=true is incompatible with power_flow.solver=apslf. The APSLF start-value generator only makes sense ahead of the rectangular Newton-Raphson solve; set power_flow.apslf_start.enabled=false or power_flow.solver=rectangular."))
  end
  wrong_branch_min_vm_pu = _validate_nonnegative("power_flow.wrong_branch_min_vm_pu", _as_float_cfg(_raw_get(merged, "wrong_branch_min_vm_pu", 0.70)))
  wrong_branch_max_vm_pu = _validate_positive("power_flow.wrong_branch_max_vm_pu", _as_float_cfg(_raw_get(merged, "wrong_branch_max_vm_pu", 1.30)))
  wrong_branch_min_vm_pu <= wrong_branch_max_vm_pu || throw(ArgumentError("power_flow.wrong_branch_min_vm_pu must be <= power_flow.wrong_branch_max_vm_pu."))
  wrong_branch_min_low_vm_count = _as_int_cfg(_raw_get(merged, "wrong_branch_min_low_vm_count", 1))
  wrong_branch_min_low_vm_count >= 0 || throw(ArgumentError("power_flow.wrong_branch_min_low_vm_count must be >= 0."))
  wrong_branch_rescue_max_attempts = _as_int_cfg(_raw_get(merged, "wrong_branch_rescue_max_attempts", 2))
  wrong_branch_rescue_max_attempts >= 0 || throw(ArgumentError("power_flow.wrong_branch_rescue_max_attempts must be >= 0."))
  autodamp = _as_bool_cfg(_raw_get(merged, "autodamp", false))
  merit_cfg = MeritLineSearchConfig(merit_raw)
  if merit_cfg.enabled && !autodamp
    throw(ArgumentError("power_flow.merit.enabled=true requires power_flow.autodamp=true. The merit-function line search only acts as an acceptance criterion inside the autodamp backtracking loop; set power_flow.autodamp=true or power_flow.merit.enabled=false."))
  end
  trust_region_cfg = TrustRegionConfig(trust_region_raw)
  if trust_region_cfg.enabled && autodamp
    throw(ArgumentError("power_flow.trust_region.enabled=true is incompatible with power_flow.autodamp=true. Both are Newton step-control mechanisms; layering them is undefined. Set power_flow.autodamp=false or power_flow.trust_region.enabled=false."))
  end
  distributed_slack_cfg = DistributedSlackConfig(distributed_slack_raw)
  external_grid_cfg = ExternalGridConfig(external_grid_raw)
  if external_grid_cfg.enabled && distributed_slack_cfg.enabled
    throw(
      ArgumentError(
        "power_flow.external_grid.enabled=true is incompatible with power_flow.distributed_slack.enabled=true. Both decide who covers the island's power imbalance: the external-grid source imports it through the feeder, the distributed slack spreads it over the island's generators. Combined, the source's import would be forced to its participation share of zero and the source degenerates to a bare angle anchor. Enable at most one.",
      ),
    )
  end
  start_mode_cfg = StartModeConfig(merge(Dict{Any,Any}(merged), Dict{Any,Any}(start_raw)))
  if start_mode_cfg.dc_seed_unconditional
    solver === :dc && throw(ArgumentError("power_flow.start_mode.dc_seed_unconditional=true is redundant with power_flow.solver=dc: the solver already is the standalone DC power flow. Set power_flow.start_mode.dc_seed_unconditional=false or power_flow.solver=rectangular."))
    solver === :apslf && throw(ArgumentError("power_flow.start_mode.dc_seed_unconditional=true is incompatible with power_flow.solver=apslf. DC-seeding only makes sense ahead of the rectangular Newton-Raphson solve; set power_flow.start_mode.dc_seed_unconditional=false or power_flow.solver=rectangular."))
    apslf_start_cfg.enabled && throw(ArgumentError("power_flow.start_mode.dc_seed_unconditional=true is incompatible with power_flow.apslf_start.enabled=true. These are two mutually exclusive start-value sources for the rectangular Newton-Raphson solve; enable at most one."))
  end
  return PowerFlowConfig(
    method = method,
    mode = _validate_allowed_symbol("power_flow.mode", _as_symbol_cfg(_raw_get(merged, "mode", :manual)), POWERFLOW_MODE_VALUES),
    solver = solver,
    linear_solver = _validate_allowed_symbol("power_flow.linear_solver", _as_symbol_cfg(_raw_get(merged, "linear_solver", :umfpack)), POWERFLOW_LINEAR_SOLVER_VALUES),
    apslf = apslf_cfg,
    apslf_start = apslf_start_cfg,
    tol = _validate_positive("powerflow.tol", _as_float_cfg(_raw_get(merged, "tol", 1.0e-8))),
    tol_MW = _validate_tol_mw(_raw_get(merged, "tol_MW", nothing)),
    max_iter = _as_int_cfg(_raw_get(merged, "max_iter", _raw_get(merged, "max_ite", 30))),
    autodamp = autodamp,
    autodamp_min = _validate_positive("powerflow.autodamp_min", _as_float_cfg(_raw_get(merged, "autodamp_min", 0.05))),
    auto_slack = _as_bool_cfg(_raw_get(merged, "auto_slack", false)),
    rescue = _as_bool_cfg(_raw_get(merged, "rescue", false)),
    wrong_branch_detection = _validate_allowed_symbol("power_flow.wrong_branch_detection", _as_symbol_cfg(_raw_get(merged, "wrong_branch_detection", :warn)), WRONG_BRANCH_DETECTION_VALUES),
    wrong_branch_rescue = _as_bool_cfg(_raw_get(merged, "wrong_branch_rescue", false)),
    wrong_branch_min_vm_pu = wrong_branch_min_vm_pu,
    wrong_branch_max_vm_pu = wrong_branch_max_vm_pu,
    wrong_branch_max_angle_spread_deg = _validate_nonnegative("power_flow.wrong_branch_max_angle_spread_deg", _as_float_cfg(_raw_get(merged, "wrong_branch_max_angle_spread_deg", 180.0))),
    wrong_branch_max_branch_angle_deg = _validate_nonnegative("power_flow.wrong_branch_max_branch_angle_deg", _as_float_cfg(_raw_get(merged, "wrong_branch_max_branch_angle_deg", 90.0))),
    wrong_branch_min_low_vm_count = wrong_branch_min_low_vm_count,
    wrong_branch_rescue_max_attempts = wrong_branch_rescue_max_attempts,
    rectangular_workspace_reuse = _as_bool_cfg(_raw_get(merged, "rectangular_workspace_reuse", true)),
    rectangular_preallocate_workspace = _validate_allowed_symbol("power_flow.rectangular_preallocate_workspace", _as_symbol_cfg(_raw_get(merged, "rectangular_preallocate_workspace", :auto)), RECTANGULAR_PREALLOCATE_WORKSPACE_VALUES),
    rectangular_workspace_min_buses = _as_int_cfg(_raw_get(merged, "rectangular_workspace_min_buses", 1000)),
    islands_enabled = _as_bool_cfg(_raw_get(islands_raw, "enabled", true)),
    islands_mode = _validate_allowed_symbol("power_flow.islands.mode", _as_symbol_cfg(_raw_get(islands_raw, "mode", :solve_independent)), POWERFLOW_ISLAND_MODE_VALUES),
    islands_reference_policy = _validate_allowed_symbol("power_flow.islands.reference_policy", _as_symbol_cfg(_raw_get(islands_raw, "reference_policy", :matpower_like)), POWERFLOW_ISLAND_REFERENCE_POLICY_VALUES),
    start_mode = start_mode_cfg,
    start_current_iteration = StartCurrentIterationConfig(start_current_iteration_raw),
    merit = merit_cfg,
    trust_region = trust_region_cfg,
    qlimits = QLimitConfig(merge(Dict{Any,Any}(merged), Dict{Any,Any}(qlimit_raw))),
    islands = IslandPowerFlowConfig(islands_raw),
    dc = DcPowerFlowConfig(dc_raw),
    distributed_slack = distributed_slack_cfg,
    external_grid = external_grid_cfg,
  )
end

function DcPowerFlowConfig(raw::AbstractDict)
  return DcPowerFlowConfig(
    angle_reference_deg = _as_float_cfg(_raw_get(raw, "angle_reference_deg", 0.0)),
    ignore_out_of_service = _as_bool_cfg(_raw_get(raw, "ignore_out_of_service", true)),
    fallback = _as_bool_cfg(_raw_get(raw, "fallback", false)),
  )
end

function ExternalGridConfig(raw::AbstractDict)
  source = _validate_allowed_symbol("power_flow.external_grid.source", _as_symbol_cfg(_raw_get(raw, "source", :auto)), EXTERNAL_GRID_SOURCE_VALUES)
  sk = _validate_positive("power_flow.external_grid.sk_MVA", _as_float_cfg(_raw_get(raw, "sk_MVA", 2000.0)))
  rx = _validate_nonnegative("power_flow.external_grid.rx", _as_float_cfg(_raw_get(raw, "rx", 0.1)))
  return ExternalGridConfig(enabled = _as_bool_cfg(_raw_get(raw, "enabled", false)), source = source, sk_MVA = sk, rx = rx)
end

function DistributedSlackConfig(raw::AbstractDict)
  enabled = _as_bool_cfg(_raw_get(raw, "enabled", false))
  p_mode = _validate_allowed_symbol("power_flow.distributed_slack.p_mode", _as_symbol_cfg(_raw_get(raw, "p_mode", :pg_weighted)), DISTRIBUTED_SLACK_P_MODE_VALUES)
  fallback = _validate_allowed_symbol("power_flow.distributed_slack.fallback", _as_symbol_cfg(_raw_get(raw, "fallback", :error)), DISTRIBUTED_SLACK_FALLBACK_VALUES)
  weights_raw = _raw_get(raw, "weights", Dict{Any,Any}())
  # The minimal in-repo YAML reader has no flow-mapping support: an empty
  # `weights: {}` arrives as the literal string "{}", a bare `weights:` as
  # nothing/"". All of these mean "no weights"; real tables use block style.
  if weights_raw === nothing || (weights_raw isa AbstractString && strip(weights_raw) in ("", "{}"))
    weights_raw = Dict{Any,Any}()
  end
  weights_raw isa AbstractDict || throw(ArgumentError("power_flow.distributed_slack.weights must be a mapping of bus name/index to weight (block style, one `bus: weight` line per entry)."))
  weights = Dict{String,Float64}()
  for (k, v) in weights_raw
    w = _as_float_cfg(v)
    (isfinite(w) && w >= 0.0) || throw(ArgumentError("power_flow.distributed_slack.weights[$(k)] must be finite and >= 0, got $(v)."))
    weights[string(k)] = w
  end
  if enabled && p_mode === :explicit
    (!isempty(weights) && any(>(0.0), values(weights))) || throw(ArgumentError("power_flow.distributed_slack.p_mode=explicit requires non-empty weights with at least one weight > 0."))
  end
  return DistributedSlackConfig(enabled = enabled, p_mode = p_mode, respect_p_limits = _as_bool_cfg(_raw_get(raw, "respect_p_limits", true)), fallback = fallback, weights = weights)
end

function IslandPowerFlowConfig(raw::AbstractDict)
  merged = haskey(raw, "islands") || haskey(raw, :islands) ? _merged_section(raw, "islands") : raw
  return IslandPowerFlowConfig(
    enabled = _as_bool_cfg(_raw_get(merged, "enabled", true)),
    mode = _validate_allowed_symbol("power_flow.islands.mode", _as_symbol_cfg(_raw_get(merged, "mode", :solve_independent)), POWERFLOW_ISLAND_MODE_VALUES),
    reference_policy = _validate_allowed_symbol("power_flow.islands.reference_policy", _as_symbol_cfg(_raw_get(merged, "reference_policy", :matpower_like)), POWERFLOW_ISLAND_REFERENCE_POLICY_VALUES),
    diagnostic_continue_after_failure = _as_bool_cfg(_raw_get(merged, "diagnostic_continue_after_failure", true)),
  )
end

function ObservabilityConfig(raw::AbstractDict)
  merged = _merged_section(raw, "observability")
  return ObservabilityConfig(enabled = _as_bool_cfg(_raw_get(merged, "enabled", true)))
end

function StateEstimationConfig(raw::AbstractDict)
  merged = _merged_section(raw, "state_estimation")
  if haskey(merged, "sparse") || haskey(merged, :sparse)
    throw(ArgumentError("state_estimation.sparse is obsolete and no longer configurable. Sparse matrices are mandatory."))
  end
  # robust-mode precedence: an explicit robust_mode wins; otherwise the
  # legacy Bool `robust: true` selects the staged modification
  se_robust = _as_bool_cfg(_raw_get(merged, "robust", false))
  se_robust_mode_raw = _raw_get(merged, "robust_mode", nothing)
  se_robust_mode = se_robust_mode_raw === nothing ? (se_robust ? :staged : :off) : _validate_allowed_symbol("state_estimation.robust_mode", _as_symbol_cfg(se_robust_mode_raw), STATE_ESTIMATION_ROBUST_MODE_VALUES)
  se_robust_k1 = Float64(_validate_positive("state_estimation.robust_k1", _as_float_cfg(_raw_get(merged, "robust_k1", 3.0))))
  se_robust_k2 = Float64(_validate_positive("state_estimation.robust_k2", _as_float_cfg(_raw_get(merged, "robust_k2", 6.0))))
  se_robust_k2 >= se_robust_k1 || throw(ArgumentError("state_estimation.robust_k2 must be >= robust_k1 (got k1=$(se_robust_k1), k2=$(se_robust_k2))."))
  return StateEstimationConfig(
    enabled = _as_bool_cfg(_raw_get(merged, "enabled", true)),
    method = _validate_allowed_symbol("state_estimation.method", _as_symbol_cfg(_raw_get(merged, "method", :wls)), STATE_ESTIMATION_METHOD_VALUES),
    tol = _validate_positive("state_estimation.tol", _as_float_cfg(_raw_get(merged, "tol", 1.0e-6))),
    max_iter = _as_int_cfg(_raw_get(merged, "max_iter", 50)),
    flatstart = _as_bool_cfg(_raw_get(merged, "flatstart", true)),
    jac_eps = _validate_positive("state_estimation.jac_eps", _as_float_cfg(_raw_get(merged, "jac_eps", 1.0e-6))),
    update_net = _as_bool_cfg(_raw_get(merged, "update_net", true)),
    pmu_ref_offset = _validate_allowed_symbol("state_estimation.pmu_ref_offset", _as_symbol_cfg(_raw_get(merged, "pmu_ref_offset", :auto)), STATE_ESTIMATION_PMU_REF_OFFSET_VALUES),
    imag_activation_iteration = Int(_validate_positive("state_estimation.imag_activation_iteration", _as_int_cfg(_raw_get(merged, "imag_activation_iteration", 2)))),
    report_residual_correlation = _as_bool_cfg(_raw_get(merged, "report_residual_correlation", false)),
    update_shunts = _as_bool_cfg(_raw_get(merged, "update_shunts", false)),
    update_taps = _as_bool_cfg(_raw_get(merged, "update_taps", false)),
    robust = se_robust,
    robust_start_iteration = Int(_validate_positive("state_estimation.robust_start_iteration", _as_int_cfg(_raw_get(merged, "robust_start_iteration", 3)))),
    k_eliminate = Float64(_validate_positive("state_estimation.k_eliminate", _as_float_cfg(_raw_get(merged, "k_eliminate", 3.0)))),
    robust_mode = se_robust_mode,
    robust_k1 = se_robust_k1,
    robust_k2 = se_robust_k2,
    k_suppress = Float64(_validate_positive("state_estimation.k_suppress", _as_float_cfg(_raw_get(merged, "k_suppress", 4.0)))),
    suppression_sigma = Float64(_validate_positive("state_estimation.suppression_sigma", _as_float_cfg(_raw_get(merged, "suppression_sigma", 2000.0)))),
    max_eliminations = Int(_validate_nonnegative("state_estimation.max_eliminations", _as_int_cfg(_raw_get(merged, "max_eliminations", 3)))),
    rank_tol_factor = Float64(_validate_positive("state_estimation.rank_tol_factor", _as_float_cfg(_raw_get(merged, "rank_tol_factor", 10.0)))),
    takahashi_min_states = Int(_validate_positive("state_estimation.takahashi_min_states", _as_int_cfg(_raw_get(merged, "takahashi_min_states", 200)))),
    ia_current_floor_A = Float64(_validate_positive("state_estimation.ia_current_floor_A", _as_float_cfg(_raw_get(merged, "ia_current_floor_A", 10.0)))),
    topology_precheck = _as_bool_cfg(_raw_get(merged, "topology_precheck", true)),
    topology_open_flow_k = Float64(_validate_positive("state_estimation.topology_open_flow_k", _as_float_cfg(_raw_get(merged, "topology_open_flow_k", 4.0)))),
    topology_dead_flow_k = Float64(_validate_positive("state_estimation.topology_dead_flow_k", _as_float_cfg(_raw_get(merged, "topology_dead_flow_k", 3.0)))),
    topology_voltage_k = Float64(_validate_positive("state_estimation.topology_voltage_k", _as_float_cfg(_raw_get(merged, "topology_voltage_k", 4.0)))),
    topology_kcl_k = Float64(_validate_positive("state_estimation.topology_kcl_k", _as_float_cfg(_raw_get(merged, "topology_kcl_k", 4.0)))),
    topology_cluster_min = Int(_validate_positive("state_estimation.topology_cluster_min", _as_int_cfg(_raw_get(merged, "topology_cluster_min", 3)))),
    observability = ObservabilityConfig(merged),
  )
end

function CGMESImportConfig(raw::AbstractDict)
  merged = _merged_section(raw, "cgmes")
  return CGMESImportConfig(
    path = strip(_as_string_cfg(_raw_get(merged, "path", ""))),
    base_mva = _validate_positive("cgmes_import.base_mva", _as_float_cfg(_raw_get(merged, "base_mva", 100.0))),
    require_boundary = _as_bool_cfg(_raw_get(merged, "require_boundary", true)),
    tap_control = _as_bool_cfg(_raw_get(merged, "tap_control", false)),
    machine_control = _as_bool_cfg(_raw_get(merged, "machine_control", false)),
    ignore_connected = _as_bool_cfg(_raw_get(merged, "ignore_connected", false)),
    vset_min_pu = _validate_nonnegative("cgmes_import.vset_min_pu", _as_float_cfg(_raw_get(merged, "vset_min_pu", 0.5))),
    vset_max_pu = _validate_positive("cgmes_import.vset_max_pu", _as_float_cfg(_raw_get(merged, "vset_max_pu", 1.5))),
    multi_slack = _as_bool_cfg(_raw_get(merged, "multi_slack", true)),
    start_values = _validate_allowed_symbol("cgmes_import.start_values", _as_symbol_cfg(_raw_get(merged, "start_values", :auto)), CGMES_START_VALUES_VALUES),
    placeholder_guards = _validate_allowed_symbol("cgmes_import.placeholder_guards", _as_symbol_cfg(_raw_get(merged, "placeholder_guards", :warn_skip)), CGMES_PLACEHOLDER_GUARDS_VALUES),
    infer_base_voltages = _as_bool_cfg(_raw_get(merged, "infer_base_voltages", false)),
    hvdc_mode = _validate_allowed_symbol("cgmes_import.hvdc_mode", _as_symbol_cfg(_raw_get(merged, "hvdc_mode", :injections)), CGMES_HVDC_MODE_VALUES),
  )
end

function ShortCircuitConfig(raw::AbstractDict)
  merged = _merged_section(raw, "short_circuit")
  c = _as_float_cfg(_raw_get(merged, "c_factor", 0.0))
  # 0.0 is the "use the IEC 60909-0 Table 1 defaults" sentinel; a real
  # override must stay within the physically sensible band around 1.0.
  (c == 0.0 || (isfinite(c) && 0.5 <= c <= 1.2)) || throw(ArgumentError("short_circuit.c_factor must be 0 (automatic IEC 60909-0 Table 1) or within [0.5, 1.2]; got $(c)."))
  return ShortCircuitConfig(
    c_factor = c,
    sweep_method = _validate_allowed_symbol("short_circuit.sweep_method", _as_symbol_cfg(_raw_get(merged, "sweep_method", :auto)), SHORT_CIRCUIT_SWEEP_METHOD_VALUES),
    takahashi_min_buses = Int(_validate_positive("short_circuit.takahashi_min_buses", _as_int_cfg(_raw_get(merged, "takahashi_min_buses", 50)))),
  )
end

function MatpowerImportConfig(raw::AbstractDict)
  merged = _merged_section(raw, "matpower")
  return MatpowerImportConfig(
    pv_voltage_source = _validate_allowed_symbol("matpower_import.pv_voltage_source", _as_symbol_cfg(_raw_get(merged, "pv_voltage_source", _raw_get(merged, "matpower_pv_voltage_source", :gen_vg))), MATPOWER_PV_VOLTAGE_SOURCE_VALUES),
    pv_voltage_mismatch_tol_pu = _validate_nonnegative("matpower.pv_voltage_mismatch_tol_pu", _as_float_cfg(_raw_get(merged, "pv_voltage_mismatch_tol_pu", _raw_get(merged, "matpower_pv_voltage_mismatch_tol_pu", 1.0e-4)))),
    compare_voltage_reference = _validate_allowed_symbol("matpower_import.compare_voltage_reference", _as_symbol_cfg(_raw_get(merged, "compare_voltage_reference", :imported_setpoint)), MATPOWER_COMPARE_VOLTAGE_REFERENCE_VALUES),
    shift_unit = _validate_allowed_symbol("matpower_import.shift_unit", _as_symbol_cfg(_raw_get(merged, "shift_unit", _raw_get(merged, "matpower_shift_unit", :deg))), MATPOWER_SHIFT_UNIT_VALUES),
    shift_sign = _as_float_cfg(_raw_get(merged, "shift_sign", _raw_get(merged, "matpower_shift_sign", 1.0))),
    ratio = _validate_allowed_symbol("matpower_import.ratio", _as_symbol_cfg(_raw_get(merged, "ratio", _raw_get(merged, "matpower_ratio", :normal))), MATPOWER_RATIO_VALUES),
    enable_pq_gen_controllers = _as_bool_cfg(_raw_get(merged, "enable_pq_gen_controllers", true)),
    apply_bus_names = _as_bool_cfg(_raw_get(merged, "apply_bus_names", false)),
    apply_branch_names = _as_bool_cfg(_raw_get(merged, "apply_branch_names", false)),
    apply_branch_kind = _as_bool_cfg(_raw_get(merged, "apply_branch_kind", false)),
    import_for001_contingencies = _as_bool_cfg(_raw_get(merged, "import_for001_contingencies", true)),
    matpower_dcline_mode = _normalize_dcline_mode(_validate_allowed_symbol("matpower_import.matpower_dcline_mode", _as_symbol_cfg(_raw_get(merged, "matpower_dcline_mode", :pf_injections)), MATPOWER_DCLINE_MODE_VALUES)),
  )
end

function ModelConfig(raw::AbstractDict)
  merged = _merged_section(raw, "model")
  return ModelConfig(
    bus_shunt_model = _validate_allowed_symbol("model.bus_shunt_model", _as_symbol_cfg(_raw_get(merged, "bus_shunt_model", :admittance)), MATPOWER_BUS_SHUNT_MODEL_VALUES),
    tap_changer_model = _validate_allowed_symbol("model.tap_changer_model", _as_symbol_cfg(_raw_get(merged, "tap_changer_model", :ideal)), TRANSFORMER_TAP_CHANGER_MODEL_VALUES),
    auto_profile = _validate_allowed_symbol("model.auto_profile", _as_auto_profile_symbol_cfg(_raw_get(merged, "auto_profile", :recommend)), MATPOWER_AUTO_PROFILE_VALUES),
    auto_profile_log = _as_bool_cfg(_raw_get(merged, "auto_profile_log", true)),
    net_cache_enabled = _as_bool_cfg(_raw_get(merged, "net_cache_enabled", false)),
    preallocate_network = _validate_allowed_symbol("model.preallocate_network", _as_symbol_cfg(_raw_get(merged, "preallocate_network", :auto)), [:off, :on, :auto]),
    preallocate_min_buses = _as_int_cfg(_raw_get(merged, "preallocate_min_buses", 1000)),
  )
end

function MatpowerExportConfig(raw::AbstractDict)
  merged = _merged_section(raw, "matpower_export")
  return MatpowerExportConfig(
    write_solution = _as_bool_cfg(_raw_get(merged, "write_solution", true)),
  )
end

"""
    configured_matpower_cases(config) -> Vector{String}

Return configured cases in deterministic execution order. A non-empty
`runtime.cases` list takes precedence over the compatible single-case
`runtime.case` setting.
"""
function configured_matpower_cases(runtime_cfg::RuntimeConfig)::Vector{String}
  !isempty(runtime_cfg.cases) && return copy(runtime_cfg.cases)
  case_name = strip(runtime_cfg.case)
  return isempty(case_name) ? String[] : [case_name]
end

configured_matpower_cases(cfg::SparlectraConfig)::Vector{String} = configured_matpower_cases(cfg.runtime)

function PerformanceConfig(raw::AbstractDict)
  merged = _merged_section(raw, "performance")
  return PerformanceConfig(
    enabled = _as_bool_cfg(_raw_get(merged, "enabled", _raw_get(merged, "performance_enabled", false))),
    level = _validate_allowed_symbol("performance.level", _as_symbol_cfg(_raw_get(merged, "level", _raw_get(merged, "performance_level", :summary))), PERFORMANCE_LEVEL_VALUES),
    print_to_console = _as_bool_cfg(_raw_get(merged, "print_to_console", _raw_get(merged, "performance_print_to_console", true))),
    write_to_logfile = _as_bool_cfg(_raw_get(merged, "write_to_logfile", _raw_get(merged, "performance_write_to_logfile", true))),
    show_allocations = _as_bool_cfg(_raw_get(merged, "show_allocations", _raw_get(merged, "performance_show_allocations", false))),
    show_iteration_table = _as_bool_cfg(_raw_get(merged, "show_iteration_table", _raw_get(merged, "performance_show_iteration_table", false))),
    compact_logging = _as_bool_cfg(_raw_get(merged, "compact_logging", _raw_get(merged, "performance_compact_logging", true))),
    representative_warmup_runs = _as_int_cfg(_raw_get(merged, "representative_warmup_runs", 0)),
    compare_cold_warm = _as_bool_cfg(_raw_get(merged, "compare_cold_warm", false)),
    skip_reference_comparison = _as_bool_cfg(_raw_get(merged, "skip_reference_comparison", _raw_get(merged, "performance_skip_reference_comparison", false))),
    skip_expensive_diagnostics = _as_bool_cfg(_raw_get(merged, "skip_expensive_diagnostics", _raw_get(merged, "performance_skip_expensive_diagnostics", true))),
    skip_branch_neighborhood_report = _as_bool_cfg(_raw_get(merged, "skip_branch_neighborhood_report", _raw_get(merged, "performance_skip_branch_neighborhood_report", true))),
    max_diagnostic_rows = _as_int_cfg(_raw_get(merged, "max_diagnostic_rows", _raw_get(merged, "performance_max_diagnostic_rows", 25))),
  )
end

function ParallelRuntimeConfig(raw::AbstractDict)
  max_tasks = _as_string_cfg(_raw_get(raw, "max_tasks", "auto"))
  normalized = lowercase(strip(max_tasks))
  if normalized != "auto"
    parsed = tryparse(Int, normalized)
    (parsed !== nothing && parsed >= 1) || throw(ArgumentError("runtime.parallel.max_tasks must be \"auto\" or a positive integer string, got $(repr(max_tasks))."))
  end
  min_work_items = _as_int_cfg(_raw_get(raw, "min_work_items", 4))
  min_work_items >= 1 || throw(ArgumentError("runtime.parallel.min_work_items must be >= 1, got $(min_work_items)."))
  return ParallelRuntimeConfig(enabled = _as_bool_cfg(_raw_get(raw, "enabled", true)), max_tasks = normalized == "auto" ? "auto" : normalized, min_work_items = min_work_items)
end

function RuntimeConfig(raw::AbstractDict)
  merged = _merged_section(raw, "runtime")
  parallel_raw = _raw_get(merged, "parallel", Dict{String,Any}())
  parallel_raw isa AbstractDict || throw(ArgumentError("runtime.parallel must be a mapping (enabled/max_tasks/min_work_items)."))
  return RuntimeConfig(
    julia_threads = _as_string_cfg(_raw_get(merged, "julia_threads", "keep")),
    blas_threads = _as_string_cfg(_raw_get(merged, "blas_threads", "keep")),
    print_thread_config = _as_bool_cfg(_raw_get(merged, "print_thread_config", true)),
    case = strip(_as_string_cfg(_raw_get(merged, "case", ""))),
    cases = _as_string_vector_cfg(_raw_get(merged, "cases", String[])),
    casefile = _as_string_cfg(_raw_get(merged, "casefile", "")),
    case_name = _as_string_cfg(_raw_get(merged, "case_name", "")),
    case_source = _as_string_cfg(_raw_get(merged, "case_source", "")),
    configured_default_casefile = _as_string_cfg(_raw_get(merged, "configured_default_casefile", "")),
    parallel = ParallelRuntimeConfig(parallel_raw),
  )
end

function BenchmarkConfig(raw::AbstractDict)
  merged = _merged_section(raw, "benchmark")
  methods = _as_symbol_vector_cfg(_raw_get(merged, "methods", [:rectangular]))
  for method in methods
    _validate_rectangular_powerflow_options(method = method, sparse = true)
  end
  return BenchmarkConfig(
    enabled = _as_bool_cfg(_raw_get(merged, "enabled", true)),
    methods = methods,
    seconds = _validate_positive("benchmark.seconds", _as_float_cfg(_raw_get(merged, "seconds", 2.0))),
    samples = _as_int_cfg(_raw_get(merged, "samples", 50)),
    show_once = _as_bool_cfg(_raw_get(merged, "show_once", false)),
    show_once_output = _validate_allowed_symbol("benchmark.show_once_output", _as_symbol_cfg(_raw_get(merged, "show_once_output", :classic)), BENCHMARK_SHOW_ONCE_OUTPUT_VALUES),
    show_once_max_nodes = _as_int_cfg(_raw_get(merged, "show_once_max_nodes", 0)),
  )
end

function ContingencyConfig(raw::AbstractDict)
  merged = _merged_section(raw, "contingency")
  ladder = _as_symbol_vector_cfg(_raw_get(merged, "rescue_ladder", [:warm]))
  # validate the ladder the same way the runContingencies! keyword does
  # (non-empty, allowed stages only, no duplicates)
  ladder = _validate_contingency_ladder(ladder; context = "contingency.rescue_ladder")
  screening = _raw_section(merged, "screening")
  mode = _validate_allowed_symbol("contingency.screening.mode", _as_symbol_cfg(_raw_get(screening, "mode", :off)), CONTINGENCY_SCREENING_MODE_VALUES)
  margin = _as_float_cfg(_raw_get(screening, "margin_pct", 10.0))
  (isfinite(margin) && margin >= 0.0) || throw(ArgumentError("contingency.screening.margin_pct must be a finite value >= 0; got $(margin)."))
  return ContingencyConfig(rescue_ladder = ladder, screening_mode = mode, screening_margin_pct = margin)
end

# The deprecated diagnostics.* duplicates of output.* are warned about (and
# ignored) once, in `_validate_known_config_keys` — see
# `_DEPRECATED_CONFIG_KEYS`. This constructor only reads the surviving key.
function DiagnosticsConfig(raw::AbstractDict)
  merged = _merged_section(raw, "diagnostics")
  return DiagnosticsConfig(log_effective_config = _as_bool_cfg(_raw_get(merged, "log_effective_config", false)))
end

_output_nonnegative_or_default(value::Integer, default::Integer) = value < 0 ? default : value
_output_positive_or_default(value::Integer, default::Integer) = value <= 0 ? default : value

function OutputConfig(raw::AbstractDict)
  merged = _merged_section(raw, "output")
  return OutputConfig(
    console_summary = _as_bool_cfg(_raw_get(merged, "console_summary", true)),
    console_live = _as_bool_cfg(_raw_get(merged, "console_live", false)),
    console_auto_profile = _validate_allowed_symbol("output.console_auto_profile", _as_symbol_cfg(_raw_get(merged, "console_auto_profile", :compact)), OUTPUT_CONSOLE_AUTO_PROFILE_VALUES),
    console_diagnostics = _validate_allowed_symbol("output.console_diagnostics", _as_symbol_cfg(_raw_get(merged, "console_diagnostics", :compact)), OUTPUT_CONSOLE_DIAGNOSTICS_VALUES),
    console_q_limit_events = _validate_allowed_symbol("output.console_q_limit_events", _as_symbol_cfg(_raw_get(merged, "console_q_limit_events", :summary)), OUTPUT_CONSOLE_Q_LIMIT_EVENTS_VALUES),
    console_max_rows = _as_int_cfg(_raw_get(merged, "console_max_rows", 100)),
    logfile_results = _validate_allowed_symbol("output.logfile_results", _as_symbol_cfg(_raw_get(merged, "logfile_results", :off)), OUTPUT_LOGFILE_RESULTS_VALUES),
    result_table_max_rows = _as_int_cfg(_raw_get(merged, "result_table_max_rows", 200)),
    result_table_large_case_threshold_buses = _as_int_cfg(_raw_get(merged, "result_table_large_case_threshold_buses", 1000)),
    result_table_large_case_mode = _validate_allowed_symbol("output.result_table_large_case_mode", _as_symbol_cfg(_raw_get(merged, "result_table_large_case_mode", :summary)), OUTPUT_RESULT_TABLE_LARGE_CASE_MODE_VALUES),
    detailed_result_csv_write_mode = _validate_allowed_symbol("output.detailed_result_csv_write_mode", _as_symbol_cfg(_raw_get(merged, "detailed_result_csv_write_mode", :auto)), OUTPUT_DETAILED_RESULT_CSV_WRITE_MODE_VALUES),
    detailed_result_csv_exporter = _validate_allowed_symbol("output.detailed_result_csv_exporter", _as_symbol_cfg(_raw_get(merged, "detailed_result_csv_exporter", :auto)), OUTPUT_DETAILED_RESULT_CSV_EXPORTER_VALUES),
    detailed_result_csv_direct_threshold_buses = _output_positive_or_default(_as_int_cfg(_raw_get(merged, "detailed_result_csv_direct_threshold_buses", 10_000)), 10_000),
    detailed_result_csv_buffer_initial_bytes = _output_nonnegative_or_default(_as_int_cfg(_raw_get(merged, "detailed_result_csv_buffer_initial_bytes", 8 * 1024 * 1024)), 8 * 1024 * 1024),
    detailed_result_csv_buffer_max_bytes = _output_positive_or_default(_as_int_cfg(_raw_get(merged, "detailed_result_csv_buffer_max_bytes", 64 * 1024 * 1024)), 64 * 1024 * 1024),
    detailed_result_csv_streaming_threshold_rows = _output_positive_or_default(_as_int_cfg(_raw_get(merged, "detailed_result_csv_streaming_threshold_rows", 100_000)), 100_000),
    logfile_diagnostics = _validate_allowed_symbol("output.logfile_diagnostics", _as_symbol_cfg(_raw_get(merged, "logfile_diagnostics", :compact)), OUTPUT_LOGFILE_DIAGNOSTICS_VALUES),
    logfile_performance = _validate_allowed_symbol("output.logfile_performance", _as_symbol_cfg(_raw_get(merged, "logfile_performance", :compact)), OUTPUT_LOGFILE_PERFORMANCE_VALUES),
    logfile_warnings = _validate_allowed_symbol("output.logfile_warnings", _as_symbol_cfg(_raw_get(merged, "logfile_warnings", :table)), OUTPUT_LOGFILE_WARNINGS_VALUES),
    startup_latency_hint = _as_bool_cfg(_raw_get(merged, "startup_latency_hint", true)),
  )
end

function ControlConfig(raw::AbstractDict)
  merged = _merged_section(raw, "control")
  # Named-mapping schema (issue #305): controllers is a mapping of
  # name => {type, kwargs}. Normalization also accepts the empty "{}"
  # placeholder and pre-normalized vectors; structural validation (known
  # type, allowed keys, required keys) happens here at load time so a bad
  # declaration fails before any solve. Network references are validated
  # at apply time by applyConfiguredControllers!.
  controllers = _normalize_controller_entries(_raw_get(merged, "controllers", nothing))
  _validate_controller_entries(controllers)
  return ControlConfig(
    enabled = _as_bool_cfg(_raw_get(merged, "enabled", true)),
    max_outer_iterations = _as_int_cfg(_raw_get(merged, "max_outer_iterations", 20)),
    trace = _as_bool_cfg(_raw_get(merged, "trace", true)),
    log_iterations = _as_bool_cfg(_raw_get(merged, "log_iterations", true)),
    stop_on_pf_failure = _as_bool_cfg(_raw_get(merged, "stop_on_pf_failure", true)),
    controllers = Any[controllers...],
  )
end

function WebUIConfig(raw::AbstractDict)
  merged = _merged_section(raw, "webui")
  return WebUIConfig(
    show_case_settings_notice = _as_bool_cfg(_raw_get(merged, "show_case_settings_notice", true)),
    operation_log_retention_days = Int(_validate_nonnegative("webui.operation_log_retention_days", _as_int_cfg(_raw_get(merged, "operation_log_retention_days", 10)))),
  )
end

function SparlectraConfig(raw::AbstractDict)
  return SparlectraConfig(
    powerflow = PowerFlowConfig(raw),
    state_estimation = StateEstimationConfig(raw),
    matpower = MatpowerImportConfig(raw),
    cgmes = CGMESImportConfig(raw),
    shortcircuit = ShortCircuitConfig(raw),
    matpower_export = MatpowerExportConfig(raw),
    model = ModelConfig(raw),
    performance = PerformanceConfig(raw),
    benchmark = BenchmarkConfig(raw),
    contingency = ContingencyConfig(raw),
    runtime = RuntimeConfig(raw),
    diagnostics = DiagnosticsConfig(raw),
    output = OutputConfig(raw),
    control = ControlConfig(raw),
    webui = WebUIConfig(raw),
  )
end

"""
    _copy_sparlectra_with_user_keys(cfg, keys) -> SparlectraConfig

Attach the set of dotted keys the USER set explicitly. Kept out of the raw
configuration dictionary on purpose: anything in there is written to
`effective_config.yaml`, and a reload of that file would reject an unknown
key.
"""
function _copy_sparlectra_with_user_keys(cfg::SparlectraConfig, keys::AbstractSet{String})::SparlectraConfig
  fields = NamedTuple{fieldnames(SparlectraConfig)}(getfield.(Ref(cfg), fieldnames(SparlectraConfig)))
  return SparlectraConfig(; fields..., user_set_keys = Set{String}(keys))
end

_canonical_config_key(key::AbstractString)::String = key == "powerflow" ? "power_flow" : key == "matpower" ? "matpower_import" : String(key)

function _merge_config_overrides(raw::Dict{String,Any}, overrides::AbstractDict)
  merged = deepcopy(raw)
  for (key, value) in overrides
    skey = _canonical_config_key(_config_key(key))
    if value isa AbstractDict
      section = get!(merged, skey, Dict{String,Any}())
      section isa AbstractDict || throw(ArgumentError("Cannot merge override section $(repr(skey)) into scalar value."))
      merged[skey] = _merge_config_overrides(Dict{String,Any}(String(k) => v for (k, v) in section), value)
    else
      merged[skey] = value
    end
  end
  return merged
end

# Free-form mapping keys: the child keys are user data (e.g. bus names), not
# configuration keys, so they are exempt from unknown-key validation. The
# default file keeps a scalar "{}" placeholder because the minimal YAML
# reader has no flow-mapping support.
const _FREEFORM_MAPPING_CONFIG_KEYS = ("power_flow.distributed_slack.weights", "control.controllers")

# Deprecated keys stay ACCEPTED so existing user/webui configuration files
# keep loading (a hard "unknown key" error here would brick every stored
# config that still carries them). They are ignored with a warning naming the
# replacement. The diagnostics.* entries were never-read duplicates of
# output.* (logging cleanup, 2026-07-30).
# Removed keys that are ignored WITHOUT a warning: the feature became
# always-on or went away entirely, so a leftover key in an existing
# user/webui configuration file carries no intent worth nagging about
# (output.condition_number, 0.9.7: the condition line is now printed
# unconditionally; webui.warmup, 0.10.0: the startup warm-up is gone,
# the sysimage does that work now). Without these entries every stored
# configuration carrying the key would fail to load with "unknown key".
const _REMOVED_SILENT_CONFIG_KEYS = ("output.condition_number", "webui.warmup")

const _DEPRECATED_CONFIG_KEYS = Dict(
  "diagnostics.console_summary" => "output.console_summary",
  "diagnostics.console_auto_profile" => "output.console_auto_profile",
  "diagnostics.console_diagnostics" => "output.console_diagnostics",
  "diagnostics.console_q_limit_events" => "output.console_q_limit_events",
  "diagnostics.console_max_rows" => "output.console_max_rows",
  "diagnostics.logfile_diagnostics" => "output.logfile_diagnostics",
)

"""
Current version of the configuration file format. Files declare theirs with
the top-level key `config_version`; a file without it reads as version 0 and
the alias tables below are applied stepwise, so old files keep loading while
only this file knows old key names (design decision D9 of the adapter task).
"""
const CONFIG_VERSION_CURRENT = 1

# One alias list per version step, old dotted key => new dotted key. Applied
# in ascending order (0 => 1, later 1 => 2, ...) before unknown-key
# validation, each application warns once naming both keys. The 0 => 1 moves
# are the model-block and runtime-case moves of design decision D6.
const _CONFIG_ALIASES = Dict{Int,Vector{Pair{String,String}}}(
  0 => [
    "matpower_import.bus_shunt_model" => "model.bus_shunt_model",
    "matpower_import.auto_profile" => "model.auto_profile",
    "matpower_import.auto_profile_log" => "model.auto_profile_log",
    "matpower_import.net_cache.enabled" => "model.net_cache_enabled",
    "matpower_import.preallocate_network" => "model.preallocate_network",
    "matpower_import.preallocate_min_buses" => "model.preallocate_min_buses",
    "transformer.tap_changer_model" => "model.tap_changer_model",
    "matpower_import.case" => "runtime.case",
    "matpower_import.cases" => "runtime.cases",
  ],
)

_dotted_path_parts(key::AbstractString) = String.(split(String(key), '.'))

function _dotted_config_get(raw::AbstractDict, key::AbstractString)
  current = raw
  parts = _dotted_path_parts(key)
  for (i, part) in enumerate(parts)
    current isa AbstractDict || return nothing
    haskey(current, part) || return nothing
    i == length(parts) && return current[part]
    current = current[part]
  end
  return nothing
end

function _dotted_config_set!(raw::AbstractDict, key::AbstractString, value)
  parts = _dotted_path_parts(key)
  current = raw
  for part in parts[1:end-1]
    next = get(current, part, nothing)
    if !(next isa AbstractDict)
      next = Dict{String,Any}()
      current[part] = next
    end
    current = next
  end
  current[parts[end]] = value
  return raw
end

function _dotted_config_delete!(raw::AbstractDict, key::AbstractString)
  parts = _dotted_path_parts(key)
  stack = AbstractDict[raw]
  current = raw
  for part in parts[1:end-1]
    next = get(current, part, nothing)
    next isa AbstractDict || return raw
    push!(stack, next)
    current = next
  end
  haskey(current, parts[end]) && delete!(current, parts[end])
  # drop section dicts the deletion emptied so an aliased-away section does
  # not linger as an unknown empty key
  for i in length(stack):-1:2
    isempty(stack[i]) || break
    delete!(stack[i-1], parts[i-1])
  end
  return raw
end

"""
    _config_file_version(raw, context) -> Int

The declared `config_version` of a configuration dictionary. Missing reads
as 0 with one warning (pre-versioning file); a version newer than
[`CONFIG_VERSION_CURRENT`](@ref) is an error, not a guess.
"""
function _config_file_version(raw::AbstractDict, context::AbstractString)::Int
  haskey(raw, "config_version") || begin
    # maxlog: the statement is the same for every file, and a session that
    # reads dozens of version-less files (the test suite writes one per
    # case) drowned the real output in identical warnings. The first one
    # carries the information; the rest carried only noise.
    @warn "Configuration file $(context) declares no config_version; reading it as version 0 and applying the documented key aliases. Further version-less files in this session are read the same way without repeating this." maxlog = 1
    return 0
  end
  value = raw["config_version"]
  version = value isa Integer ? Int(value) : tryparse(Int, strip(string(value)))
  version === nothing && throw(ArgumentError("config_version in $(context) must be an integer, got $(repr(value))."))
  version > CONFIG_VERSION_CURRENT && throw(ArgumentError("Configuration file $(context) declares config_version $(version); this Sparlectra reads at most $(CONFIG_VERSION_CURRENT). Update Sparlectra or lower the file version."))
  version < 0 && throw(ArgumentError("config_version in $(context) must be >= 0, got $(version)."))
  return version
end

function _validate_config_scope(raw::AbstractDict, expected::AbstractString, context::AbstractString)
  haskey(raw, "scope") || return nothing
  scope = lowercase(strip(string(raw["scope"])))
  scope == expected || throw(ArgumentError("Configuration file $(context) declares scope $(repr(scope)); expected $(repr(String(expected))) here (a case configuration belongs next to its case file, not in its place)."))
  return nothing
end

"""
    _apply_config_aliases!(raw, version, context) -> raw

Rewrite a version-`version` configuration dictionary to the current layout
by applying every alias step in order. A value is moved only when the new
key is not set (an explicitly set new key wins over a stale old one).
"""
function _apply_config_aliases!(raw::AbstractDict, version::Int, context::AbstractString)
  for step in version:(CONFIG_VERSION_CURRENT-1)
    for (old_key, new_key) in get(_CONFIG_ALIASES, step, Pair{String,String}[])
      value = _dotted_config_get(raw, old_key)
      value === nothing && continue
      # maxlog per key, not per file: the interesting part is WHICH old key
      # is still in use, and that repeats across files
      @warn "Configuration key $(old_key) in $(context) is a version-$(step) name; applied as $(new_key). The file can be updated with refresh_sparlectra_config_file." maxlog = 1 _id = Symbol("cfg_alias_", old_key)
      _dotted_config_delete!(raw, old_key)
      _dotted_config_get(raw, new_key) === nothing && _dotted_config_set!(raw, new_key, value)
    end
  end
  return raw
end

# Move version-alias source keys to their current names, keeping their
# values. Unlike the loader's _apply_config_aliases! this runs over EVERY
# alias step regardless of the declared config_version: a source key is
# invalid at any version that knows the target, so finding one always
# means an unmigrated or half-migrated file. Half-migrated files exist in
# the field: the refresh before this helper existed stamped config_version
# 1 and filled the model block from template defaults while leaving the
# old matpower_import keys (and their user values) in place. Rules per
# key: the user's stale value wins over a target that is missing or still
# equals the template default; an explicitly set target wins over the
# stale source. The source key is removed either way.
function _migrate_versioned_config_aliases!(raw::AbstractDict, defaults::AbstractDict)::Vector{String}
  migrated = String[]
  for step in 0:(CONFIG_VERSION_CURRENT-1)
    for (old_key, new_key) in get(_CONFIG_ALIASES, step, Pair{String,String}[])
      value = _dotted_config_get(raw, old_key)
      value === nothing && continue
      target = _dotted_config_get(raw, new_key)
      if target === nothing || isequal(target, _dotted_config_get(defaults, new_key))
        _dotted_config_set!(raw, new_key, value)
      end
      _dotted_config_delete!(raw, old_key)
      push!(migrated, old_key)
    end
  end
  return migrated
end

function _validate_known_config_keys(user::AbstractDict, defaults::AbstractDict; path::String = "")
  for (key, value) in user
    skey = _canonical_config_key(_config_key(key))
    isempty(path) && skey in ("_config_sources", "_config_metadata") && continue
    current_path = isempty(path) ? skey : string(path, ".", skey)
    if current_path == "matpower_import.benchmark"
      throw(ArgumentError("Removed Sparlectra configuration key: matpower_import.benchmark.\nUse top-level benchmark.enabled instead, e.g.\n\nbenchmark:\n  enabled: true"))
    end
    if isempty(path) && skey == "methods"
      for m in _as_symbol_vector_cfg(value)
        requested = _powerflow_method_label(m)
        requested === SUPPORTED_POWERFLOW_METHOD || throw(ArgumentError(unsupported_powerflow_method_message(requested)))
      end
      continue
    end
    if haskey(_DEPRECATED_CONFIG_KEYS, current_path)
      @warn "Configuration key $(current_path) is deprecated and ignored — use $(_DEPRECATED_CONFIG_KEYS[current_path])."
      continue
    end
    current_path in _REMOVED_SILENT_CONFIG_KEYS && continue
    haskey(defaults, skey) || throw(ArgumentError("Unknown Sparlectra configuration key: $(current_path)"))
    if value isa AbstractDict
      current_path in _FREEFORM_MAPPING_CONFIG_KEYS && continue
      defaults[skey] isa AbstractDict || throw(ArgumentError("Configuration key $(current_path) must be a scalar value."))
      _validate_known_config_keys(value, defaults[skey]; path = current_path)
    end
  end
  return nothing
end

_config_file_hash(path::AbstractString) = isfile(path) ? bytes2hex(sha256(read(path))) : ""
_config_file_mtime(path::AbstractString) = isfile(path) ? stat(path).mtime : 0.0

function _load_and_validate_config(default_path::AbstractString, user_path::AbstractString; cli_overrides::AbstractDict, overrides::AbstractDict, user_set_out::Union{Nothing,Set{String}} = nothing, case_scope_from_defaults::Bool = false)
  isfile(default_path) || throw(ArgumentError("Default Sparlectra config file not found: $(default_path)"))
  if abspath(user_path) != abspath(USER_SPARLECTRA_CONFIG_PATH) && !isfile(user_path)
    throw(ArgumentError("Sparlectra user config file not found: $(user_path)"))
  end
  defaults = load_yaml_dict(default_path)
  user = isfile(user_path) ? load_yaml_dict(user_path) : Dict{String,Any}()
  if isfile(user_path)
    user_version = _config_file_version(user, user_path)
    _validate_config_scope(user, "general", user_path)
    user_version < CONFIG_VERSION_CURRENT && _apply_config_aliases!(user, user_version, user_path)
  end
  _validate_known_config_keys(user, defaults)
  _validate_known_config_keys(cli_overrides, defaults)
  _validate_known_config_keys(overrides, defaults)
  canon_user = Dict{String,Any}(String(_canonical_config_key(_config_key(k))) => v for (k, v) in user)
  if case_scope_from_defaults
    # Runs with a case configuration file (issue #1 point 1, decided for
    # `defaults`; format independent): the general file's CASE-scope keys
    # do not participate, so the same case-plus-config pair reads to the
    # same numbers on every installation. A case-scope key resolves
    # case configuration -> packaged defaults; only machine-scope keys
    # (output, benchmark, runtime, webui, matpower_export) still come from
    # this machine's file. Masking BEFORE the merge also keeps the masked
    # keys out of user_set below, so the auto profile stays free on them.
    for k in collect(keys(_flatten_config_values!(Dict{String,Any}(), canon_user)))
      scf_is_case_config_key(k) && _dotted_config_delete!(canon_user, k)
    end
  end
  raw = merge_yaml_dict!(deepcopy(defaults), canon_user)
  raw = _merge_config_overrides(raw, cli_overrides)
  raw = _merge_config_overrides(raw, overrides)
  _validate_powerflow_solver_config!(raw)
  if user_set_out !== nothing
    # every dotted key the USER decided; the auto power-flow mode treats
    # these as untouchable. User YAML files are usually full template
    # copies (the refresh tooling writes every key), so a file key counts
    # only when its value DIFFERS from the template default; override keys
    # always count (their presence is the decision).
    dflat = _flatten_config_values!(Dict{String,Any}(), defaults)
    uflat = _flatten_config_values!(Dict{String,Any}(), canon_user)
    for (k, v) in uflat
      (haskey(dflat, k) && isequal(dflat[k], v)) || push!(user_set_out, k)
    end
    _flatten_config_keys!(user_set_out, cli_overrides)
    _flatten_config_keys!(user_set_out, overrides)
  end
  return raw, isfile(user_path)
end

# nested config dict -> flat dotted-key => scalar value map (read-only helper
# for the user-set detection above)
function _flatten_config_values!(out::Dict{String,Any}, raw::AbstractDict, prefix::String = "")
  for (key, value) in raw
    path = isempty(prefix) ? String(key) : string(prefix, ".", key)
    value isa AbstractDict ? _flatten_config_values!(out, value, path) : (out[path] = value)
  end
  return out
end

"""
    load_sparlectra_config([user_path]; reload=false, cli_overrides=Dict(), overrides=Dict())

Load, validate, merge, and cache the central typed `SparlectraConfig`.
Configuration precedence is `src/config/configuration.yaml.example`, optional
`examples/configuration.yaml` (or an explicit `user_path`), CLI-style overrides,
then explicit Julia API overrides. Unknown user keys throw an `ArgumentError`.
"""

"""
    refresh_sparlectra_config_file(path; write=false, backup=true, normalize_deprecated=true, default_path=DEFAULT_SPARLECTRA_CONFIG_PATH)

Compare a user YAML configuration with the current Sparlectra template, add
missing keys from the template, and optionally normalize documented deprecated
aliases. Dry-run mode returns the refreshed YAML without writing. Writes are
explicit and create a timestamped backup unless `backup=false` is requested.
"""
function refresh_sparlectra_config_file(path::AbstractString; write::Bool = false, backup::Bool = true, normalize_deprecated::Bool = true, default_path::AbstractString = DEFAULT_SPARLECTRA_CONFIG_PATH)
  isfile(path) || throw(ArgumentError("Sparlectra config file not found: $(path)"))
  isfile(default_path) || throw(ArgumentError("Default Sparlectra config file not found: $(default_path)"))
  duplicate_keys = _detect_yaml_duplicate_keys(path)
  warnings = String[]
  for key in duplicate_keys
    push!(warnings, "Duplicate YAML key detected: $(key). Refresh did not write because duplicate keys require manual review.")
  end
  defaults = load_yaml_dict(default_path)
  user = load_yaml_dict(path)
  refreshed = deepcopy(user)
  # migrate version-alias keys BEFORE filling from the template: filling
  # first would stamp config_version 1 and seed the new keys with defaults
  # while the user's values are still parked under the old names
  normalized_keys = normalize_deprecated ? _migrate_versioned_config_aliases!(refreshed, defaults) : String[]
  missing_keys = String[]
  _add_missing_config_keys!(refreshed, defaults, missing_keys)
  normalize_deprecated && append!(normalized_keys, _normalize_deprecated_config_aliases!(refreshed))
  refreshed_text = _yaml_dict_text(refreshed)
  original_text = read(path, String)
  changed = !isempty(missing_keys) || !isempty(normalized_keys) || refreshed_text != original_text
  backup_path = nothing
  written = false
  success = true
  if write && changed
    if !isempty(duplicate_keys)
      success = false
    else
      if backup
        backup_path = string(path, ".bak-", Dates.format(now(), "yyyymmdd-HHMMSS"))
        cp(path, backup_path; force = false)
      end
      Base.write(path, refreshed_text)
      written = true
    end
  end
  return (;
    success,
    changed,
    written,
    backup_path,
    missing_keys = sort!(missing_keys),
    normalized_keys = sort!(normalized_keys),
    duplicate_keys = sort!(duplicate_keys),
    warnings,
    refreshed_text,
  )
end

"""
    refresh_sparlectra_config_text(text; normalize_deprecated, default_path) -> String

Rewrite a configuration YAML text against the packaged template: missing keys
are added with their defaults and comments, deprecated keys are normalized,
user values are kept. Used by the Web UI configuration refresh.
"""
function refresh_sparlectra_config_text(text::AbstractString; normalize_deprecated::Bool = true, default_path::AbstractString = DEFAULT_SPARLECTRA_CONFIG_PATH)
  mktemp() do path, io
    write(io, text)
    close(io)
    return refresh_sparlectra_config_file(path; write = false, backup = false, normalize_deprecated, default_path)
  end
end

function _detect_yaml_duplicate_keys(path::AbstractString)::Vector{String}
  duplicates = String[]
  stack_names = Dict{Int,Vector{String}}(0 => String[])
  seen = Dict{Int,Set{String}}(0 => Set{String}())
  open(path, "r") do io
    for (lineno, rawline) in enumerate(eachline(io))
      line = rstrip(_strip_yaml_comment(rawline))
      isempty(strip(line)) && continue
      indent_count = count(==(' '), first(line, something(findfirst(!=(' '), line), lastindex(line) + 1) - 1))
      indent_count % 2 == 0 || continue
      content = strip(line)
      m = match(r"^([A-Za-z0-9_\-]+)\s*:\s*(.*)$", content)
      isnothing(m) && continue
      level = indent_count ÷ 2
      key = String(m.captures[1])
      parent_path = get(stack_names, level, String[])
      level_seen = get!(seen, level, Set{String}())
      compound = isempty(parent_path) ? key : string(join(parent_path, "."), ".", key)
      if key in level_seen
        push!(duplicates, compound)
      else
        push!(level_seen, key)
      end
      if isempty(String(m.captures[2]))
        stack_names[level + 1] = [parent_path; key]
        seen[level + 1] = Set{String}()
      end
      for old_level in collect(keys(stack_names))
        old_level > level + 1 && delete!(stack_names, old_level)
      end
      for old_level in collect(keys(seen))
        old_level > level + 1 && delete!(seen, old_level)
      end
    end
  end
  return unique(duplicates)
end

function _add_missing_config_keys!(dst::Dict{String,Any}, defaults::Dict{String,Any}, missing::Vector{String}; prefix::String = "")
  for key in sort(collect(keys(defaults)))
    path = isempty(prefix) ? key : string(prefix, ".", key)
    if !haskey(dst, key)
      dst[key] = deepcopy(defaults[key])
      push!(missing, path)
    elseif dst[key] isa Dict{String,Any} && defaults[key] isa Dict{String,Any}
      _add_missing_config_keys!(dst[key], defaults[key], missing; prefix = path)
    end
  end
  return dst
end

function _ensure_config_section!(root::Dict{String,Any}, keys::AbstractVector{<:AbstractString})::Dict{String,Any}
  current = root
  for key in keys
    if !haskey(current, key) || !(current[key] isa Dict{String,Any})
      current[key] = Dict{String,Any}()
    end
    current = current[key]
  end
  return current
end

function _normalize_deprecated_config_aliases!(root::Dict{String,Any})::Vector{String}
  normalized = String[]
  pf = haskey(root, "power_flow") ? root["power_flow"] : get(root, "powerflow", nothing)
  if pf isa Dict{String,Any}
    start_mode = get(pf, "start_mode", nothing)
    if start_mode isa Dict{String,Any} && get(start_mode, "voltage_mode", nothing) in ("bus_vm_va_blend", :bus_vm_va_blend)
      start_mode["voltage_mode"] = "profile_blend"
      start_mode["profile_source"] = "matpower_reference"
      push!(normalized, "power_flow.start_mode.voltage_mode")
      push!(normalized, "power_flow.start_mode.profile_source")
    end
    qlimits = get(pf, "qlimits", nothing)
    if qlimits isa Dict{String,Any}
      mode = get(qlimits, "enforcement_mode", nothing)
      if mode in ("matpower_simultaneous", :matpower_simultaneous)
        qlimits["enforcement_mode"] = "classic_simultaneous"
        push!(normalized, "power_flow.qlimits.enforcement_mode")
      elseif mode in ("matpower_one_at_a_time", :matpower_one_at_a_time)
        qlimits["enforcement_mode"] = "classic_one_at_a_time"
        push!(normalized, "power_flow.qlimits.enforcement_mode")
      end
    end
  end
  # The deprecated diagnostics.* duplicates of output.* are DROPPED, not
  # migrated: they were never read, so copying their values into output.*
  # would silently change behavior a stored config never had.
  diag = get(root, "diagnostics", nothing)
  if diag isa Dict{String,Any}
    for path in sort!(collect(keys(_DEPRECATED_CONFIG_KEYS)))
      startswith(path, "diagnostics.") || continue
      key = split(path, '.'; limit = 2)[2]
      if haskey(diag, key)
        delete!(diag, key)
        push!(normalized, path)
      end
    end
  end
  return unique(normalized)
end

# A string is written bare only when the repository's own YAML reader
# gives it back unchanged as a string; everything else is quoted. This is
# the self-inverse criterion, not a maintained list: it covers the YAML
# 1.1 boolean family (live case: the enum value "off" was written bare
# and came back as `false`), null spellings, every number form the reader
# knows including leading-zero integers, ":symbol" spellings,
# inline-list lookalikes, and padded values, and it stays correct when
# parse_yaml_scalar learns new forms.
function _yaml_ambiguous_scalar(s::AbstractString)
  v = parse_yaml_scalar(String(s))
  return !(v isa AbstractString && String(v) == String(s))
end

function _yaml_scalar_text(value)
  value === nothing && return "null"
  value isa Bool && return value ? "true" : "false"
  value isa Symbol && return _yaml_ambiguous_scalar(String(value)) ? repr(String(value)) : String(value)
  value isa AbstractString && return (occursin(r"[:#\[\],{}]|^$", value) || _yaml_ambiguous_scalar(value)) ? repr(String(value)) : String(value)
  value isa AbstractVector && return string("[", join((_yaml_scalar_text(v) for v in value), ", "), "]")
  return string(value)
end

function _write_yaml_dict(io::IO, data::Dict{String,Any}; indent::Int = 0)
  pad = repeat(" ", indent)
  for key in sort(collect(keys(data)))
    value = data[key]
    if value isa Dict{String,Any}
      println(io, pad, key, ":")
      _write_yaml_dict(io, value; indent = indent + 2)
    else
      println(io, pad, key, ": ", _yaml_scalar_text(value))
    end
  end
end

function _yaml_dict_text(data::Dict{String,Any})::String
  io = IOBuffer()
  _write_yaml_dict(io, data)
  return String(take!(io))
end
"""
    load_sparlectra_config(user_path; default_path, reload, cli_overrides, overrides) -> SparlectraConfig

Load the configuration from the user YAML over the packaged defaults, apply
overrides, and return it WITHOUT installing it as the active configuration.
Results are cached per path unless `reload` is set.
"""
function load_sparlectra_config(user_path::AbstractString = USER_SPARLECTRA_CONFIG_PATH; default_path::AbstractString = DEFAULT_SPARLECTRA_CONFIG_PATH, reload::Bool = false, cli_overrides::AbstractDict = Dict{String,Any}(), overrides::AbstractDict = Dict{String,Any}())
  abspath_default = abspath(default_path)
  abspath_user = abspath(user_path)
  default_mtime = _config_file_mtime(abspath_default)
  user_mtime = _config_file_mtime(abspath_user)
  default_digest = _config_file_hash(abspath_default)
  user_digest = _config_file_hash(abspath_user)
  cached = _SPARLECTRA_CONFIG_CACHE[]
  if !reload &&
     cached !== nothing &&
     cached.default_path == abspath_default &&
     cached.user_path == abspath_user &&
     cached.default_mtime == default_mtime &&
     cached.user_mtime == user_mtime &&
     cached.default_hash == default_digest &&
     cached.user_hash == user_digest &&
     isempty(cli_overrides) &&
     isempty(overrides)
    return cached.config
  end

  user_set = Set{String}()
  raw, user_found = _load_and_validate_config(abspath_default, abspath_user; cli_overrides = cli_overrides, overrides = overrides, user_set_out = user_set)
  # the user-set key list never enters the raw dict: it would travel into
  # effective_config.yaml and be rejected as an unknown key on reload
  config = _copy_sparlectra_with_user_keys(SparlectraConfig(raw), user_set)
  if isempty(cli_overrides) && isempty(overrides)
    _SPARLECTRA_CONFIG_CACHE[] = (; default_path = abspath_default, user_path = abspath_user, default_mtime, user_mtime, default_hash = default_digest, user_hash = user_digest, user_found, config)
  end
  return config
end


"""
    load_sparlectra_config!(user_path; default_path, reload, cli_overrides, overrides) -> SparlectraConfig

Like [`load_sparlectra_config`](@ref), and additionally installs the result
as the active configuration of the process.
"""
function load_sparlectra_config!(user_path::AbstractString = USER_SPARLECTRA_CONFIG_PATH; default_path::AbstractString = DEFAULT_SPARLECTRA_CONFIG_PATH, reload::Bool = false, cli_overrides::AbstractDict = Dict{String,Any}(), overrides::AbstractDict = Dict{String,Any}())
  return set_sparlectra_config!(load_sparlectra_config(user_path; default_path = default_path, reload = reload, cli_overrides = cli_overrides, overrides = overrides))
end

function _print_config_section(io::IO, name::AbstractString, value; indent::AbstractString = "")
  println(io, indent, name, ":")
  child_indent = string(indent, "  ")
  for field in fieldnames(typeof(value))
    field_value = getfield(value, field)
    if field_value isa StartModeConfig || field_value isa QLimitConfig || field_value isa ObservabilityConfig
      _print_config_section(io, String(field), field_value; indent = child_indent)
    else
      println(io, child_indent, field, ": ", field_value)
    end
  end
  return nothing
end

"""
    print_effective_config([io], config::SparlectraConfig)

Print an `Effective Sparlectra Configuration` block with typed module sections.
"""
function print_effective_config(io::IO, config::SparlectraConfig)
  println(io, "==================== Effective Sparlectra Configuration ====================")
  println(io, "default_file: ", DEFAULT_SPARLECTRA_CONFIG_PATH)
  println(io, "user_file: ", USER_SPARLECTRA_CONFIG_PATH)
  println(io, "user_file_found: ", isfile(USER_SPARLECTRA_CONFIG_PATH))
  _print_config_section(io, "power_flow", config.powerflow)
  _print_config_section(io, "state_estimation", config.state_estimation)
  _print_config_section(io, "matpower_import", config.matpower)
  _print_config_section(io, "model", config.model)
  _print_config_section(io, "performance", config.performance)
  _print_config_section(io, "runtime", config.runtime)
  _print_config_section(io, "diagnostics", config.diagnostics)
  _print_config_section(io, "output", config.output)
  println(io, "==========================================================================")
  return nothing
end

print_effective_config(config::SparlectraConfig) = print_effective_config(stdout, config)
