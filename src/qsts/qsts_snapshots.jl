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

# file: src/qsts/qsts_snapshots.jl
# purpose: QSTS mode A (qsts.mode = snapshots): every step of a series is
#          solved on its own, from a fresh copy of the unsolved template
#          plus the full delta of that step, warm started from the solution
#          of the base case (SCF) or started from the state of its own
#          snapshot (CGMES). The steps run in parallel through the batch
#          loop of the scenario engine (_engine_parallel_map); the
#          ScenarioEngine itself is not used. The result is one long table
#          (step, element_id, quantity, value), written as CSV with a step
#          table, an element table and a manifest. The result types and the
#          writer serve mode B (qsts_sequential.jl) as well. The CGMES method
#          of runQsts lives in src/adapters/cgmes/glue.jl.

"""
    QSTS_BUS_QUANTITIES

Quantities per bus of a snapshot run, in the order of the second axis of
`QstsResult.bus_values`: `vm_pu`, `va_deg` (the solved state) and
`sv_dvm_pu`, `sv_dva_deg` (model minus SV of the snapshot, the angle
aligned per electrical island).
"""
const QSTS_BUS_QUANTITIES = ["vm_pu", "va_deg", "sv_dvm_pu", "sv_dva_deg"]

"""
    QSTS_BRANCH_QUANTITIES

Quantities per branch of a snapshot run, in the order of the second axis
of `QstsResult.branch_values`: the terminal flows `p_from_mw`,
`q_from_mvar`, `p_to_mw`, `q_to_mvar` (from `calcNetLosses!`),
`loading_pct` (only with a rating), `p_loss_mw`, `q_loss_mvar`,
`tap_ratio` and `shift_deg` (only for a transformer with a tap
controller), and the SV flow differences `sv_dp_from_mw`,
`sv_dq_from_mvar`, `sv_dp_to_mw`, `sv_dq_to_mvar`.
"""
const QSTS_BRANCH_QUANTITIES = ["p_from_mw", "q_from_mvar", "p_to_mw", "q_to_mvar", "loading_pct", "p_loss_mw", "q_loss_mvar", "tap_ratio", "shift_deg",
                                "sv_dp_from_mw", "sv_dq_from_mvar", "sv_dp_to_mw", "sv_dq_to_mvar"]

"""
    QSTS_STEP_QUANTITIES

Quantities per step of a QSTS run, in the order of the first axis of
`QstsResult.step_values`: `converged` (0 or 1), `iterations`,
`outer_iterations`, `n_frozen`, `total_p_loss_mw`, `total_q_loss_mvar`,
`n_low_vm` (energized buses below `power_flow.plausibility.vmin_warn` in a
converged step, no value when the threshold is 0), and the SV comparison
`sv_n_nodes`, `sv_max_dvm_pu`, `sv_rms_dvm_pu`, `sv_max_dva_deg`,
`sv_rms_dva_deg`, `sv_va_offset_deg`, `sv_tn_spread_pu`, `sv_n_terminals`,
`sv_max_dp_mw`, `sv_rms_dp_mw`, `sv_max_dq_mvar`, `sv_rms_dq_mvar`.

A sequential run (`qsts.mode = sequential`) adds, for a converged step:
`resync` (1 in a resync step), `n_moves` (the counted tap and bank steps of
the step), `n_islands` (electrical islands, one without a source included),
`n_isolated` (buses without a value because they are isolated, in an SCF or
event series also those of an island without a source) and `shed_load_mw`
(the load on them),
`max_dvm_step_pu` and `max_dp_step_mw` (largest change of a bus voltage and
of a branch flow `p_from_mw` against the step before, only when that step
converged), `n_vm_low`, `n_vm_high` (buses outside `qsts.vm_min_pu`,
`qsts.vm_max_pu`) and `n_overload` (branches above
`qsts.loading_max_pct`). A snapshot run leaves them without a value.
"""
const QSTS_STEP_QUANTITIES = ["converged", "iterations", "outer_iterations", "n_frozen", "total_p_loss_mw", "total_q_loss_mvar", "n_low_vm",
                              "sv_n_nodes", "sv_max_dvm_pu", "sv_rms_dvm_pu", "sv_max_dva_deg", "sv_rms_dva_deg", "sv_va_offset_deg", "sv_tn_spread_pu",
                              "sv_n_terminals", "sv_max_dp_mw", "sv_rms_dp_mw", "sv_max_dq_mvar", "sv_rms_dq_mvar",
                              "resync", "n_moves", "n_islands", "n_isolated", "shed_load_mw", "max_dvm_step_pu", "max_dp_step_mw", "n_vm_low", "n_vm_high", "n_overload"]

"""
    QSTS_CONTROLLER_QUANTITIES

Quantities per tap and bank controller of a sequential run, in the order of
the second axis of `QstsResult.controller_values`: `position` (model step
of a modelled transformer, else the tap position counted from the ratio 1.0
as in `buildTapControllerReportRows`; a switched bank `bs_mvar /
step_mvar`; no value for a continuous SVC), `bs_mvar` (banks), `ref_position`
(the position the series gives the step), `dpos_ref` (simulated minus
reference position in steps), `moves` (net move in the step with its sign,
0 in a resync step), `at_limit` and `locked` (0 or 1; locked by the dead
time). A
snapshot run has no controller elements.
"""
const QSTS_CONTROLLER_QUANTITIES = ["position", "bs_mvar", "ref_position", "dpos_ref", "moves", "at_limit", "locked"]

# column positions in the value arrays, resolved once from the lists above
# so the two cannot drift apart
_qsts_col(list::Vector{String}, name::String)::Int = something(findfirst(==(name), list))
const _QB_VM = _qsts_col(QSTS_BUS_QUANTITIES, "vm_pu")
const _QB_VA = _qsts_col(QSTS_BUS_QUANTITIES, "va_deg")
const _QB_SV_DVM = _qsts_col(QSTS_BUS_QUANTITIES, "sv_dvm_pu")
const _QB_SV_DVA = _qsts_col(QSTS_BUS_QUANTITIES, "sv_dva_deg")
const _QR_P_FROM = _qsts_col(QSTS_BRANCH_QUANTITIES, "p_from_mw")
const _QR_Q_FROM = _qsts_col(QSTS_BRANCH_QUANTITIES, "q_from_mvar")
const _QR_P_TO = _qsts_col(QSTS_BRANCH_QUANTITIES, "p_to_mw")
const _QR_Q_TO = _qsts_col(QSTS_BRANCH_QUANTITIES, "q_to_mvar")
const _QR_LOADING = _qsts_col(QSTS_BRANCH_QUANTITIES, "loading_pct")
const _QR_P_LOSS = _qsts_col(QSTS_BRANCH_QUANTITIES, "p_loss_mw")
const _QR_Q_LOSS = _qsts_col(QSTS_BRANCH_QUANTITIES, "q_loss_mvar")
const _QR_TAP = _qsts_col(QSTS_BRANCH_QUANTITIES, "tap_ratio")
const _QR_SHIFT = _qsts_col(QSTS_BRANCH_QUANTITIES, "shift_deg")
const _QR_SV_DP_FROM = _qsts_col(QSTS_BRANCH_QUANTITIES, "sv_dp_from_mw")
const _QR_SV_DQ_FROM = _qsts_col(QSTS_BRANCH_QUANTITIES, "sv_dq_from_mvar")
const _QR_SV_DP_TO = _qsts_col(QSTS_BRANCH_QUANTITIES, "sv_dp_to_mw")
const _QR_SV_DQ_TO = _qsts_col(QSTS_BRANCH_QUANTITIES, "sv_dq_to_mvar")
const _QS_CONVERGED = _qsts_col(QSTS_STEP_QUANTITIES, "converged")
const _QS_ITERATIONS = _qsts_col(QSTS_STEP_QUANTITIES, "iterations")
const _QS_OUTER = _qsts_col(QSTS_STEP_QUANTITIES, "outer_iterations")
const _QS_FROZEN = _qsts_col(QSTS_STEP_QUANTITIES, "n_frozen")
const _QS_P_LOSS = _qsts_col(QSTS_STEP_QUANTITIES, "total_p_loss_mw")
const _QS_Q_LOSS = _qsts_col(QSTS_STEP_QUANTITIES, "total_q_loss_mvar")
const _QS_LOW_VM = _qsts_col(QSTS_STEP_QUANTITIES, "n_low_vm")
const _QS_SV_N_NODES = _qsts_col(QSTS_STEP_QUANTITIES, "sv_n_nodes")
const _QS_SV_MAX_DVM = _qsts_col(QSTS_STEP_QUANTITIES, "sv_max_dvm_pu")
const _QS_SV_RMS_DVM = _qsts_col(QSTS_STEP_QUANTITIES, "sv_rms_dvm_pu")
const _QS_SV_MAX_DVA = _qsts_col(QSTS_STEP_QUANTITIES, "sv_max_dva_deg")
const _QS_SV_RMS_DVA = _qsts_col(QSTS_STEP_QUANTITIES, "sv_rms_dva_deg")
const _QS_SV_OFFSET = _qsts_col(QSTS_STEP_QUANTITIES, "sv_va_offset_deg")
const _QS_SV_SPREAD = _qsts_col(QSTS_STEP_QUANTITIES, "sv_tn_spread_pu")
const _QS_SV_N_TERMINALS = _qsts_col(QSTS_STEP_QUANTITIES, "sv_n_terminals")
const _QS_SV_MAX_DP = _qsts_col(QSTS_STEP_QUANTITIES, "sv_max_dp_mw")
const _QS_SV_RMS_DP = _qsts_col(QSTS_STEP_QUANTITIES, "sv_rms_dp_mw")
const _QS_SV_MAX_DQ = _qsts_col(QSTS_STEP_QUANTITIES, "sv_max_dq_mvar")
const _QS_SV_RMS_DQ = _qsts_col(QSTS_STEP_QUANTITIES, "sv_rms_dq_mvar")
const _QS_RESYNC = _qsts_col(QSTS_STEP_QUANTITIES, "resync")
const _QS_N_MOVES = _qsts_col(QSTS_STEP_QUANTITIES, "n_moves")
const _QS_N_ISLANDS = _qsts_col(QSTS_STEP_QUANTITIES, "n_islands")
const _QS_N_ISOLATED = _qsts_col(QSTS_STEP_QUANTITIES, "n_isolated")
const _QS_SHED_LOAD = _qsts_col(QSTS_STEP_QUANTITIES, "shed_load_mw")
const _QS_MAX_DVM_STEP = _qsts_col(QSTS_STEP_QUANTITIES, "max_dvm_step_pu")
const _QS_MAX_DP_STEP = _qsts_col(QSTS_STEP_QUANTITIES, "max_dp_step_mw")
const _QS_N_VM_LOW = _qsts_col(QSTS_STEP_QUANTITIES, "n_vm_low")
const _QS_N_VM_HIGH = _qsts_col(QSTS_STEP_QUANTITIES, "n_vm_high")
const _QS_N_OVERLOAD = _qsts_col(QSTS_STEP_QUANTITIES, "n_overload")
const _QC_POSITION = _qsts_col(QSTS_CONTROLLER_QUANTITIES, "position")
const _QC_BS = _qsts_col(QSTS_CONTROLLER_QUANTITIES, "bs_mvar")
const _QC_REF = _qsts_col(QSTS_CONTROLLER_QUANTITIES, "ref_position")
const _QC_DPOS_REF = _qsts_col(QSTS_CONTROLLER_QUANTITIES, "dpos_ref")
const _QC_MOVES = _qsts_col(QSTS_CONTROLLER_QUANTITIES, "moves")
const _QC_AT_LIMIT = _qsts_col(QSTS_CONTROLLER_QUANTITIES, "at_limit")
const _QC_LOCKED = _qsts_col(QSTS_CONTROLLER_QUANTITIES, "locked")

const _QSTS_RESULT_FORMAT = "sparlectra-qsts-result"
# a format break (other columns, a renamed quantity, another unit, another
# address or gap rule) raises the major version; new quantities, new
# element kinds (the controllers of mode B), new files (qsts_switching.csv,
# qsts_metrics.csv) and new manifest keys do not
const _QSTS_RESULT_FORMAT_VERSION = "1.0"
const _QSTS_VALUES_HEADER = ["step", "element_id", "quantity", "value"]
const _QSTS_STEPS_HEADER = ["step", "label", "timestamp", "is_base", "status", "control_status", "iterations", "outer_iterations", "residual", "n_frozen", "reference", "message"]
const _QSTS_ELEMENTS_HEADER = ["element_id", "kind", "name", "from", "to", "vn_kv", "rating_mva"]
const _QSTS_SWITCHING_HEADER = ["step", "element_id", "old", "new", "steps", "trigger"]
const _QSTS_METRICS_HEADER = ["element_id", "metric", "value"]

# one net move of a controller in one step (qsts_switching.csv) and one
# metric (qsts_metrics.csv); positions and steps as in
# QSTS_CONTROLLER_QUANTITIES
const _QstsSwitchRow = NamedTuple{(:step, :element_id, :old, :new, :steps, :trigger),Tuple{Int,String,Int,Int,Int,Symbol}}
const _QstsMetricRow = NamedTuple{(:element_id, :metric, :value),Tuple{String,String,Float64}}

"""
    QstsSVReference

The SV state of one snapshot, the reference of the SV comparison of a
snapshot run. Built by [`buildCGMESSeriesUnion`](@ref) from the direct
import of each snapshot, with the exclusions of [`compareWithSV`](@ref):

- `tn_buses`, `vm_pu`, `va_deg`: one entry per topological node of the
  snapshot with an `SvVoltage` (`v / vn` and `angle`), with the union buses
  that node covers;
- `branch`, `side`, `p_mw`, `q_mvar`: one entry per `SvPowerFlow` at a line
  or transformer terminal, with the branch position in the union template
  and the side (`:from`, `:to`) the snapshot import gives the terminal.
"""
struct QstsSVReference
  tn_buses::Vector{Vector{Int}}
  vm_pu::Vector{Float64}
  va_deg::Vector{Float64}
  branch::Vector{Int}
  side::Vector{Symbol}
  p_mw::Vector{Float64}
  q_mvar::Vector{Float64}
end

# The input of a snapshot run, the same for both sources: the unsolved
# template, its index and the full delta of every step against it. The
# template is never changed. `start_ops` empty: a snapshot step starts warm
# from the base solution (SCF series, event series); one list per step: a
# snapshot step starts from these bus start values instead (a CGMES union,
# the state of each snapshot), and a sequential step leaves them out.
struct _QstsInput
  name::String
  origin::Symbol                                    # :cgmes_series, :scf_series
  template::Net
  index::ScenarioIndex
  ops::Vector{Vector{PatchOp}}
  start_ops::Vector{Vector{PatchOp}}
  labels::Vector{String}
  timestamps::Vector{String}                        # ISO 8601 UTC, "" without time
  base_step::Int                                    # 0: the base case is no step (SCF)
  mark_isolated::Bool                               # false: the delta carries the isolation (CGMES)
  residual::Vector{Int}                             # differences a delta cannot express, per step
  sv::Vector{Union{Nothing,QstsSVReference}}
  step_minutes::Union{Nothing,Float64}
  event_list::Dict{String,Any}                      # source and row count of an event series, empty otherwise
end

"""
    QstsStepInfo

One solve of a snapshot run (step 0 is the base case):

- `step`, `label`, `timestamp` (ISO 8601 UTC, empty without time), `is_base`
  (the step whose snapshot is the base case);
- `status`: `:solved`, `:not_converged` (the solve, or one AC island of a
  net with several, ended without convergence) or `:solver_error` (the
  solver threw; `message` holds the first line);
- `control_status`: the status of the controller run (`:converged`,
  `:max_outer_iterations`, `:blocked`, ...) with `qsts.snapshot_controls =
  true`, `:none` without controllers; a controller status other than
  convergence does not change `status`;
- `iterations`: Newton iterations (with controllers the larger of the sum
  over all passes and the last pass; for a non-converged AC island the
  iterations of that island), `outer_iterations`;
- `residual`: differences between the snapshot and template plus delta
  that the delta cannot express (CGMES), else 0;
- `n_frozen`: buses in electrical islands without a slack and without a PV
  bus; the solver holds their start voltage, so they get no values;
- `reference`: slack buses after the solve (bus positions);
- `message`: first line of the solver error; for AC islands that did not
  converge the message of the solve (island, buses, iterations) and status,
  iterations, reason and final mismatch of each such island; empty
  otherwise.
"""
struct QstsStepInfo
  step::Int
  label::String
  timestamp::String
  is_base::Bool
  status::Symbol
  control_status::Symbol
  iterations::Int
  outer_iterations::Int
  residual::Int
  n_frozen::Int
  reference::Vector{Int}
  message::String
end

"""
    QstsElements

The elements of a QSTS run, by position in the template. Bus fields:
`bus_ids` (the address, see [`runQsts`](@ref)), `bus_names`, `bus_vn_kv`,
`bus_reported` (`false` for the star point of a three-winding transformer,
which has values in `bus_values` but no rows in the long table). Branch
fields: `branch_ids`, `branch_names`, `branch_kind` (`:line`,
`:transformer`, `:equivalent_branch` for a CGMES `EquivalentBranch`), `branch_from`, `branch_to` (bus positions) and
`branch_rating_mva` (`NaN` without a rating).

Controller fields, filled by a sequential run only, one entry per tap
controller (followers belong to their master) and per shunt bank
controller in the order of `collect_outer_controllers`:
`controller_ids` (`<actuator address>#OLTC`, `#PST`, `#SCHRAEG`, `#MSC` or
`#SVC`), `controller_names` (`control_name`), `controller_kind` (`:OLTC`,
`:PST`, `:SCHRAEG`, `:MSC`, `:SVC`), `controller_actuator` (address of the
transformer, of the bank's shunt when the series addresses it, else of the
bank's bus) and `controller_target` (address of the regulated bus, empty
without one).
"""
struct QstsElements
  bus_ids::Vector{String}
  bus_names::Vector{String}
  bus_vn_kv::Vector{Float64}
  bus_reported::BitVector
  branch_ids::Vector{String}
  branch_names::Vector{String}
  branch_kind::Vector{Symbol}
  branch_from::Vector{Int}
  branch_to::Vector{Int}
  branch_rating_mva::Vector{Float64}
  controller_ids::Vector{String}
  controller_names::Vector{String}
  controller_kind::Vector{Symbol}
  controller_actuator::Vector{String}
  controller_target::Vector{String}
end

"""
    QstsResult

The result of [`runQsts`](@ref):

- `name`, `origin` (`:cgmes_series`, `:scf_series`), `mode` (`:snapshots`,
  `:sequential`),
  `base_step` (0 for an SCF series, whose base case is no step),
  `step_minutes` (`nothing` for irregular steps);
- `base`: the base case solve as a [`QstsStepInfo`](@ref) with step 0;
- `steps`: one [`QstsStepInfo`](@ref) per step;
- `elements`: the [`QstsElements`](@ref);
- `bus_values[bus, q, step]`, `branch_values[branch, q, step]`,
  `step_values[q, step]`, `controller_values[controller, q, step]` with `q`
  the position in [`QSTS_BUS_QUANTITIES`](@ref),
  [`QSTS_BRANCH_QUANTITIES`](@ref), [`QSTS_STEP_QUANTITIES`](@ref),
  [`QSTS_CONTROLLER_QUANTITIES`](@ref); `NaN` is no value;
- `switching`: the net moves of a sequential run, see [`qstsSwitching`](@ref);
  `metrics`: its metrics, see [`qstsMetrics`](@ref); both empty in a
  snapshot run;
- `max_iter`, `tol`, `snapshot_controls`, `tap_control`,
  `map_energy_source`, `qsts` (the `qsts` section of the run),
  `event_list` (source file and row count of an event series, see
  [`qstsEventSeries`](@ref), empty otherwise): the settings the manifest
  records;
- `notes`: the changes of the step configuration against the passed one,
  a base case without convergence and similar, and one note per converged
  solve (base case or step) with an energized bus below
  `power_flow.plausibility.vmin_warn`, naming the lowest bus and its
  voltage; the same text is logged as a warning.

[`qstsLongTable`](@ref) and [`qstsValues`](@ref) read it,
[`writeQstsResults`](@ref) writes it.
"""
struct QstsResult
  name::String
  origin::Symbol
  mode::Symbol
  base_step::Int
  step_minutes::Union{Nothing,Float64}
  base::QstsStepInfo
  steps::Vector{QstsStepInfo}
  elements::QstsElements
  bus_values::Array{Float64,3}
  branch_values::Array{Float64,3}
  step_values::Matrix{Float64}
  controller_values::Array{Float64,3}
  switching::Vector{_QstsSwitchRow}
  metrics::Vector{_QstsMetricRow}
  max_iter::Int
  tol::Float64
  snapshot_controls::Bool
  tap_control::Bool
  map_energy_source::Bool
  qsts::QstsConfig
  event_list::Dict{String,Any}
  notes::Vector{String}
end

# the loosely typed keywords of runQsts behind one concrete struct, so the
# run body compiles once (compile rule 2)
struct _QstsRunOptions
  parallel_enabled::Union{Nothing,Bool}
  parallel_max_tasks::Union{Nothing,Int}
  parallel_min_work_items::Union{Nothing,Int}
  progress::Any                                     # nothing, or callable (done, total)
end

# what every step reads; built once in the orchestrator and shared read-only
struct _QstsContext
  input::_QstsInput
  step_config::SparlectraConfig
  warm::Union{Nothing,Tuple{Vector{Union{Nothing,Float64}},Vector{Union{Nothing,Float64}}}}
  controls::Bool                                    # qsts.snapshot_controls
  tap_branch::BitVector                             # transformer with a tap controller
end

struct _QstsStepOut
  info::QstsStepInfo
  bus::Matrix{Float64}
  branch::Matrix{Float64}
  step::Vector{Float64}
  plausibility_warning::String                      # "" without one; logged by the orchestrator
  isolation_note::String                            # "" without one; a snapshot step that isolated an island without a source
end
_QstsStepOut(info::QstsStepInfo, bus::Matrix{Float64}, branch::Matrix{Float64}, step::Vector{Float64}, plausibility_warning::String) = _QstsStepOut(info, bus, branch, step, plausibility_warning, "")

# The settings a snapshot run (mode :snapshots) or a sequential run (mode
# :sequential) cannot use, each with its key named.
_qsts_check_config(config::SparlectraConfig) = _qsts_check_config(config, :snapshots)
function _qsts_check_config(config::SparlectraConfig, mode::Symbol)
  run = mode === :snapshots ? "a snapshot run" : "a sequential run"
  # off (the default) reaches this check from runQsts; name both run modes
  config.qsts.mode === :off && throw(ArgumentError("runQsts: qsts.mode is off; set qsts.mode = snapshots or sequential"))
  config.qsts.mode === mode || throw(ArgumentError(string("runQsts: qsts.mode is ", config.qsts.mode, "; ", run, " needs qsts.mode = ", mode)))
  pf = config.powerflow
  # DC refuses controllers and APSLF leaves them out (execution.jl)
  pf.solver === :rectangular || throw(ArgumentError(string("runQsts: power_flow.solver is ", pf.solver, "; ", run, " needs power_flow.solver = rectangular")))
  # auto chooses the start and step strategy per net, so each step would get its own
  pf.mode === :auto && throw(ArgumentError(string("runQsts: power_flow.mode = auto chooses the strategy per step; ", run, " needs power_flow.mode = manual")))
  # the external grid adds buses inside run_sparlectra only, so the template
  # and the run would not match
  pf.external_grid.enabled && throw(ArgumentError(string("runQsts: power_flow.external_grid.enabled = true adds buses the series does not know; ", run, " needs it off")))
  return nothing
end

# one note naming the start keys a run-level flat start switched off for
# this QSTS run (_flatstart_profile_config), the counterpart of the run.log
# line of a power-flow run
function _qsts_flatstart_note!(notes::Vector{String}, overridden::Vector{String})
  isempty(overridden) && return notes
  push!(notes, string("flat start: the run starts from the flat profile, no start projection; overridden for this run: ", join(overridden, ", ")))
  return notes
end

# the configuration of the base solve: the passed one, with the controllers
# off unless qsts.snapshot_controls asks for them (taps and banks then stay
# at the values of the snapshot)
function _qsts_base_config(config::SparlectraConfig)::SparlectraConfig
  config.qsts.snapshot_controls && return config
  config.control.enabled || return config
  return _sparlectra_config_with(config; control = _struct_with(config.control; enabled = false))
end

# The configuration of every step and one note per key that differs from the
# passed configuration. One attempt per step from the warm start: the start
# generators and the rescue paths would replace or re-judge that start, power
# mode would carry a factorization across steps (the net carries links), the
# parallelism sits at the steps, and verbose output would come from the tasks.
# `has_start`: the steps start from given voltages (a converged base case, a
# CGMES snapshot state), so the flat start is off.
function _qsts_step_config(config::SparlectraConfig, has_start::Bool)::Tuple{SparlectraConfig,Vector{String}}
  base = _qsts_base_config(config)
  pf = base.powerflow
  notes = String[]
  note(key, from, to) = from == to || push!(notes, string("step configuration: ", key, " = ", to, " (configured ", from, ")"))
  sm = pf.start_mode
  flat = has_start ? false : sm.flatstart
  note("power_flow.start_mode.flatstart", sm.flatstart, flat)
  note("power_flow.start_mode.start_projection", sm.start_projection, false)
  note("power_flow.start_mode.dc_seed_unconditional", sm.dc_seed_unconditional, false)
  note("power_flow.start_current_iteration.enabled", pf.start_current_iteration.enabled, false)
  note("power_flow.rescue", pf.rescue, false)
  note("power_flow.dc.fallback", pf.dc.fallback, false)
  note("power_flow.wrong_branch_rescue", pf.wrong_branch_rescue, false)
  note("power_flow.power_mode", pf.power_mode, false)
  note("runtime.parallel.enabled", base.runtime.parallel.enabled, false)
  diag = base.output.console_diagnostics === :full ? :compact : base.output.console_diagnostics
  note("output.console_diagnostics", base.output.console_diagnostics, diag)
  config.control.enabled == base.control.enabled ||
    push!(notes, "base and step configuration: control.enabled = false (configured true; qsts.snapshot_controls = false keeps taps and banks at the snapshot values)")
  step_pf = _copy_powerflow_with(pf;
    start_mode = _copy_start_mode_with(sm; flatstart = flat, start_projection = false, dc_seed_unconditional = false),
    start_current_iteration = _copy_start_current_iteration_with(pf.start_current_iteration; enabled = false),
    rescue = false, dc = _struct_with(pf.dc; fallback = false), wrong_branch_rescue = false, power_mode = false)
  step_cfg = _sparlectra_config_with(base; powerflow = step_pf,
    runtime = _struct_with(base.runtime; parallel = _struct_with(base.runtime.parallel; enabled = false)),
    output = _struct_with(base.output; console_diagnostics = diag))
  return step_cfg, notes
end

# The AC islands that ended without convergence in an island-wise solve,
# from the island statuses the solve leaves in `profile`, in island order.
# Empty when the net was not solved island by island or when an island
# threw (that is a solver error). The island-wise solve throws an
# ErrorException where the single-island solve returns erg = 1; with
# power_flow.islands.diagnostic_continue_after_failure (the default) it
# solves every island first and names the first failure. A non-converged
# island has no exception (exception_type "") and stage :newton_iteration
# or :pre_nr_setup; an island that threw has its exception type, and an
# island skipped after a failure is no failure of its own.
function _qsts_islands_not_converged(profile::Dict{Symbol,Any})::Vector{Any}
  statuses = get(profile, :ac_island_solver_statuses, nothing)
  statuses isa AbstractDict || return Any[]
  out = Any[]
  for id in sort!(collect(keys(statuses)))
    st = statuses[id]
    all(f -> hasproperty(st, f), (:status, :stage, :exception_type, :iterations)) || continue
    (st.status === :converged || st.stage === :skipped_after_previous_failure) && continue
    isempty(st.exception_type) || return Any[]
    st.stage in (:newton_iteration, :pre_nr_setup) || return Any[]
    push!(out, st)
  end
  return out
end

# One power-flow run of a snapshot run. The named failure the catch turns into
# a :solver_error step is a failure of the SOLVE (island validation, a
# singular matrix, a participation set that cannot be built);
# _rethrow_unless_solver_failure lets interrupts, aborts and defects of the
# calling code through. An AC island that did not converge in an
# island-wise solve is no such failure: the solve throws for it, but the
# step gets the outcome a non-converged single-island solve gives (erg = 1,
# control status :pf_failed when the net has outer controllers), with the
# iterations of the first island that did not converge, and the message of
# the solve plus status, iterations, reason and final mismatch of every
# such island. `profile` is a dict without :enabled: it records no timings
# and only carries the island statuses out of the solve. Returns the run
# outcome (nothing after a solver error) and the message: the first line of
# a solver error, empty after a converged run. print_island_table = false
# keeps the island table of the run out of the console (a sequential run
# solves its steps one after the other and would print it per step).
function _qsts_solve!(net::Net, cfg::SparlectraConfig; print_island_table::Bool = true)
  profile = Dict{Symbol,Any}()
  try
    return (_execute_sparlectra_powerflow!(net, cfg; performance_profile = profile, print_island_table = print_island_table), "")
  catch err
    _rethrow_unless_solver_failure(err)
    text = sprint(showerror, err)
    failed = _qsts_islands_not_converged(profile)
    isempty(failed) && return (nothing, String(first(split(text, '\n'))))
    control_status = isempty(collect_outer_controllers(net)) ? :none : :pf_failed
    ex = (iterations = Int(first(failed).iterations), erg = 1, elapsed_s = NaN, solver_elapsed_s = NaN, control_status = control_status)
    island_text(st) = (m = hasproperty(st, :final_mismatch) ? Float64(st.final_mismatch) : NaN;
                       string("island ", st.island_id, " ", st.status, " after ", st.iterations, " iterations, reason ", hasproperty(st, :reason) ? st.reason : :unknown,
                              isfinite(m) ? string(", final mismatch ", round(m; sigdigits = 4)) : ""))
    lines = (strip(l) for l in split(text, '\n') if !isempty(strip(l)))
    return (ex, string(join(lines, " "), "; ", join((island_text(st) for st in failed), "; ")))
  end
end

# Buses in islands of the power flow (_power_flow_islands: closed-link
# groups as one node, also through a bus typed Isolated) without a slack and
# without a PV node. The solver promotes a unit there and keeps its start
# voltage, so their state says nothing about the snapshot. Split at a
# switch-only bus, the side without the slack counted as frozen although
# the solver solves it with the rest. Isolated-typed buses are not counted,
# as before: they are skipped or take their group values elsewhere.
function _qsts_frozen_buses(net::Net)::Vector{Int}
  out = Int[]
  islands, held, _ = _power_flow_islands(net)
  for (comp, has_ref) in zip(islands, held)
    has_ref && continue
    append!(out, Int[b for b in comp if !isIsolated(net.nodeVec[b])])
  end
  return sort!(out)
end

# the transformers whose tap a controller of `net` moves (followers included)
function _qsts_tap_branches(net::Net)::BitVector
  out = falses(length(net.branchVec))
  for c in _tap_controllers(net), name in vcat([c.trafo], c.followers)
    br = _find_trafo_branch(net, name)
    j = findfirst(x -> x === br, net.branchVec)
    j === nothing || (out[j] = true)
  end
  return out
end

_qsts_slack_buses(net::Net)::Vector{Int} = Int[i for i in eachindex(net.nodeVec) if getNodeType(net.nodeVec[i]) == Slack]

# iterations, outer iterations and controller status of one run; without
# controllers the solve's own count, with them the larger of the sum over
# all passes and the last pass (the sum stays 0 on two paths of run_control!)
function _qsts_iterations(net::Net, ex, controls::Bool)::Tuple{Int,Int,Symbol}
  ex === nothing && return (0, 0, :none)
  (controls && ex.control_status !== :none) || return (Int(ex.iterations), 0, :none)
  cr = net.control_result
  cr === nothing && return (Int(ex.iterations), 0, ex.control_status)
  return (max(cr.total_pf_iterations, cr.last_pf_iterations), cr.outer_iterations, ex.control_status)
end

function _qsts_status(ex)::Symbol
  ex === nothing && return :solver_error
  return ex.erg == 0 ? :solved : :not_converged
end

# The SV comparison of one solved step against its reference, written into
# the value arrays. `skip` holds the buses without a solved value (isolated,
# frozen). Voltages: every union bus of a topological node gets the
# difference to the node's SV value; the metrics count each node once, at
# its smallest bus, so they match compareWithSV on the direct import. Angles
# are aligned per electrical island with the offsets of those entries.
function _qsts_sv_compare!(bus::Matrix{Float64}, branch::Matrix{Float64}, stepv::Vector{Float64}, net::Net, ref::QstsSVReference, skip::Set{Int})
  V = ComplexF64[something(nd._vm_pu, 1.0) * cis(deg2rad(something(nd._va_deg, 0.0))) for nd in net.nodeVec]
  # per island of the power flow, as compareWithSV aligns
  island = _power_flow_island_of_bus(net)
  raw_dva = Dict{Int,Float64}()
  rep_dvm = Float64[]
  rep_dva = Float64[]
  rep_island = Int[]
  spread = 0.0
  for j in eachindex(ref.tn_buses)
    buses = Int[b for b in ref.tn_buses[j] if !(b in skip)]
    isempty(buses) && continue
    vmin = Inf
    vmax = -Inf
    for b in buses
      vm = something(net.nodeVec[b]._vm_pu, NaN)
      bus[b, _QB_SV_DVM] = vm - ref.vm_pu[j]
      raw_dva[b] = something(net.nodeVec[b]._va_deg, NaN) - ref.va_deg[j]
      vmin = min(vmin, vm)
      vmax = max(vmax, vm)
    end
    spread = max(spread, vmax - vmin)
    b0 = minimum(buses)
    push!(rep_dvm, bus[b0, _QB_SV_DVM])
    push!(rep_dva, raw_dva[b0])
    push!(rep_island, get(island, b0, 0))
  end
  stepv[_QS_SV_N_NODES] = length(rep_dvm)
  if !isempty(rep_dvm)
    alignment = _sv_align_angles(rep_dva, rep_island)
    for (b, d) in raw_dva
      bus[b, _QB_SV_DVA] = d - get(alignment.offsets, get(island, b, 0), alignment.offset_deg)
    end
    dvm = abs.(rep_dvm)
    dva = abs.(alignment.aligned)
    stepv[_QS_SV_MAX_DVM] = maximum(dvm)
    stepv[_QS_SV_RMS_DVM] = sqrt(sum(abs2, dvm) / length(dvm))
    stepv[_QS_SV_MAX_DVA] = maximum(dva)
    stepv[_QS_SV_RMS_DVA] = sqrt(sum(abs2, dva) / length(dva))
    stepv[_QS_SV_OFFSET] = alignment.offset_deg
    stepv[_QS_SV_SPREAD] = spread
  end
  dp = Float64[]
  dq = Float64[]
  for j in eachindex(ref.branch)
    pos = ref.branch[j]
    br = net.branchVec[pos]
    (br.status == 1 && !(br.fromBus in skip) && !(br.toBus in skip)) || continue
    s = _sv_branch_terminal_flow(br, V, ref.side[j], net.baseMVA)
    d_p = real(s) - ref.p_mw[j]
    d_q = imag(s) - ref.q_mvar[j]
    from = ref.side[j] === :from
    branch[pos, from ? _QR_SV_DP_FROM : _QR_SV_DP_TO] = d_p
    branch[pos, from ? _QR_SV_DQ_FROM : _QR_SV_DQ_TO] = d_q
    push!(dp, abs(d_p))
    push!(dq, abs(d_q))
  end
  stepv[_QS_SV_N_TERMINALS] = length(dp)
  if !isempty(dp)
    stepv[_QS_SV_MAX_DP] = maximum(dp)
    stepv[_QS_SV_RMS_DP] = sqrt(sum(abs2, dp) / length(dp))
    stepv[_QS_SV_MAX_DQ] = maximum(dq)
    stepv[_QS_SV_RMS_DQ] = sqrt(sum(abs2, dq) / length(dq))
  end
  return nothing
end

# The buses of `skip` that the power flow solves with their link group: a
# bus without an in-service branch of its own whose closed-link group has a
# branch (the rule of markIsolatedBuses! with links = true) and a member
# outside `skip`. The solver contracts the group to one node and writes its
# voltage to every member, so such a bus has solved values although the
# import typed it Isolated. A bus with a branch that is skipped (isolated
# explicitly, frozen) stays skipped, and so does a group with no member
# outside `skip`.
function _qsts_link_fed_buses(net::Net, skip::Set{Int})::Vector{Int}
  by_branch = Set(_isolated_buses(net))
  by_group = Set(_isolated_buses(net; links = true))
  reps = _active_link_representative_map(net)
  live = Set{Int}(reps[j] for j in eachindex(net.nodeVec) if !(j in skip))
  return Int[i for i in skip if 1 <= i <= length(reps) && i in by_branch && !(i in by_group) && reps[i] in live]
end

# The values of one step, after the solve. Rows exist only for a solved
# value: no bus values for isolated and frozen buses, no branch values for
# a branch out of service or with an end on such a bus, no loading without
# a rating, and a step without convergence keeps only its step counters.
function _qsts_capture(ctx::_QstsContext, net::Net, k::Int, ex, msg::String, frozen::Vector{Int})::_QstsStepOut
  bus = fill(NaN, length(net.nodeVec), length(QSTS_BUS_QUANTITIES))
  branch = fill(NaN, length(net.branchVec), length(QSTS_BRANCH_QUANTITIES))
  stepv = fill(NaN, length(QSTS_STEP_QUANTITIES))
  status = _qsts_status(ex)
  iterations, outer, control_status = _qsts_iterations(net, ex, ctx.controls)
  stepv[_QS_CONVERGED] = status === :solved ? 1.0 : 0.0
  stepv[_QS_ITERATIONS] = iterations
  stepv[_QS_OUTER] = outer
  stepv[_QS_FROZEN] = length(frozen)
  plausibility_warning = ""
  if status === :solved
    # the check of the direct run (_build_sparlectra_result), which a step
    # bypasses because it calls _execute_sparlectra_powerflow!
    pl = _plausibility_diagnostics(net, ctx.step_config.powerflow.plausibility, true)
    pl.plausibility_status === :not_checked || (stepv[_QS_LOW_VM] = pl.plausibility_low_vm_count)
    plausibility_warning = pl.plausibility_warning
    skip = Set{Int}(frozen)
    union!(skip, net.isoNodes)
    for (i, nd) in enumerate(net.nodeVec)
      isIsolated(nd) || continue
      # a step of an SCF series is marked by the rule of the power flow
      # (links = true, _qsts_step), and isoNodes holds the result. The
      # import types a bus Isolated by branches alone and the type sticks,
      # so a bus that the template typed Isolated is judged by isoNodes: a
      # bus behind a closed link is solved with its group and keeps its
      # values. A bus the step itself typed Isolated (a bus op) is skipped.
      (ctx.input.mark_isolated && isIsolated(ctx.input.template.nodeVec[i])) && continue
      push!(skip, i)
    end
    # a step of a CGMES union carries the isolation of its snapshot import,
    # which types a bus Isolated by branches alone; it is not marked again
    ctx.input.mark_isolated || setdiff!(skip, _qsts_link_fed_buses(net, skip))
    for (i, nd) in enumerate(net.nodeVec)
      i in skip && continue
      bus[i, _QB_VM] = something(nd._vm_pu, NaN)
      bus[i, _QB_VA] = something(nd._va_deg, NaN)
    end
    p_loss = 0.0
    q_loss = 0.0
    for (j, br) in enumerate(net.branchVec)
      (br.status == 1 && !(br.fromBus in skip) && !(br.toBus in skip)) || continue
      f = br.fBranchFlow
      t = br.tBranchFlow
      branch[j, _QR_P_FROM] = f === nothing ? NaN : something(f.pFlow, NaN)
      branch[j, _QR_Q_FROM] = f === nothing ? NaN : something(f.qFlow, NaN)
      branch[j, _QR_P_TO] = t === nothing ? NaN : something(t.pFlow, NaN)
      branch[j, _QR_Q_TO] = t === nothing ? NaN : something(t.qFlow, NaN)
      # the rule of _base_branch_loadings (contingency.jl)
      rating = br.sn_MVA
      if rating !== nothing && isfinite(rating) && rating > 0.0
        s_from = f === nothing ? 0.0 : hypot(something(f.pFlow, 0.0), something(f.qFlow, 0.0))
        s_to = t === nothing ? 0.0 : hypot(something(t.pFlow, 0.0), something(t.qFlow, 0.0))
        branch[j, _QR_LOADING] = 100.0 * max(s_from, s_to) / rating
      end
      branch[j, _QR_P_LOSS] = something(br.pLosses, NaN)
      branch[j, _QR_Q_LOSS] = something(br.qLosses, NaN)
      p_loss += something(br.pLosses, 0.0)
      q_loss += something(br.qLosses, 0.0)
      if ctx.tap_branch[j]
        # the live values the solver used
        branch[j, _QR_TAP] = br.tap_ratio
        branch[j, _QR_SHIFT] = br.phase_shift_deg
      end
    end
    stepv[_QS_P_LOSS] = p_loss
    stepv[_QS_Q_LOSS] = q_loss
    ref = ctx.input.sv[k]
    ref === nothing || _qsts_sv_compare!(bus, branch, stepv, net, ref, skip)
  end
  info = QstsStepInfo(k, ctx.input.labels[k], ctx.input.timestamps[k], k == ctx.input.base_step, status, control_status, iterations, outer,
                      ctx.input.residual[k], length(frozen), status === :solved ? _qsts_slack_buses(net) : Int[], msg)
  return _QstsStepOut(info, bus, branch, stepv, plausibility_warning)
end

# A solve below power_flow.plausibility.vmin_warn: one note in the result
# and one warning in the log, no abort. The orchestrator calls it in step
# order, so a parallel run gives the same notes and log lines as a serial one.
function _qsts_plausibility_note!(notes::Vector{String}, k::Int, label::String, warning::String)
  isempty(warning) && return nothing
  what = k == 0 ? "base case" : string("step ", k, isempty(label) ? "" : string(" (", label, ")"))
  text = string("QSTS ", what, ": ", warning)
  push!(notes, text)
  @warn text
  return nothing
end

# One step on a fresh copy of the chunk template: the warm start first, so
# the start values of the delta (buses that come back, slack-less islands)
# overwrite it; a CGMES step has no warm start and takes the start of every
# live bus from its delta (start_ops). The values are taken outside the try.
function _qsts_step(ctx::_QstsContext, tmpl::Net, k::Int)::_QstsStepOut
  # a copied factorization would be shared between tasks (threads rule)
  tmpl._power_cache === nothing || error("runQsts: the chunk template carries a power-mode cache")
  work = deepcopy(tmpl)
  ctx.warm === nothing || _restore_start_voltages!(work, ctx.warm)
  # one apply!: its op order puts the bus start values after the appliance
  # ops, so a generator vm_pu op does not overwrite the snapshot start
  ops = isempty(ctx.input.start_ops) ? ctx.input.ops[k] : vcat(ctx.input.ops[k], ctx.input.start_ops[k])
  apply!(work, ops, ctx.input.index)
  # markIsolatedBuses! starts an isolated member of a fed link group at the
  # group's voltage; a CGMES delta isolates by its own bus ops (0 V) and is
  # not marked again, so its step gets the same start here
  n_isolated = 0
  shed = Dict{Int,Float64}()
  if ctx.input.mark_isolated
    markIsolatedBuses!(net = work, log = false, links = true)
    # F2 as in a sequential step: an island the step left without a source
    # (switched off, or its only unit out) is isolated with its load
    # zeroed; a CGMES delta carries its isolation itself
    iso_before = count(isIsolated, work.nodeVec)
    shed = _qsts_isolate_sourceless!(work)
    n_isolated = count(isIsolated, work.nodeVec) - iso_before
  else
    _start_isolated_link_members!(work)
  end
  frozen = _qsts_frozen_buses(work)
  ex, msg = _qsts_solve!(work, ctx.step_config)
  ex !== nothing && ex.erg == 0 && calcNetLosses!(work)
  out = _qsts_capture(ctx, work, k, ex, msg, frozen)
  n_isolated == 0 && return out
  label = ctx.input.labels[k]
  note = string("QSTS step ", k, isempty(label) ? "" : string(" (", label, ")"), ": ", n_isolated, " bus(es) without a source isolated, ",
                round(sum(values(shed); init = 0.0); digits = 3), " MW load unsupplied")
  return _QstsStepOut(out.info, out.bus, out.branch, out.step, out.plausibility_warning, note)
end

# addresses, names and ratings of the template's buses and branches; a CGMES
# index addresses by mRID (the inverse of id_by_mrid), an SCF index by id
function _qsts_elements(T::Net, index::ScenarioIndex)::QstsElements
  nb = length(T.nodeVec)
  nbr = length(T.branchVec)
  addr_of = Dict{Int,String}(id => a for (a, id) in index.id_by_mrid)
  bus_ids = fill("", nb)
  branch_ids = fill("", nbr)
  branch_kind = Symbol[_scf_is_transformer(br) ? :transformer : :line for br in T.branchVec]
  # sorted ids: the first address of a position wins, the same in every run
  for id in sort!(collect(keys(index.kind_by_id)))
    kind = index.kind_by_id[id]
    i = get(index.internal_by_id, id, 0)
    a = get(addr_of, id, string(id))
    if kind === :bus && 1 <= i <= nb && isempty(bus_ids[i])
      bus_ids[i] = a
    elseif (kind === :branch || kind === :transformer) && 1 <= i <= nbr && isempty(branch_ids[i])
      branch_ids[i] = a
      branch_kind[i] = kind === :transformer ? :transformer : :line
    end
  end
  # a CGMES EquivalentBranch is a network equivalent, neither line nor
  # transformer (net.equivalent_branch_metadata)
  for j in keys(T.equivalent_branch_metadata)
    1 <= j <= nbr && (branch_kind[j] = :equivalent_branch)
  end
  # a position the index does not address keeps a position address
  for i in 1:nb
    isempty(bus_ids[i]) && (bus_ids[i] = string("#bus", i))
  end
  for j in 1:nbr
    isempty(branch_ids[j]) && (branch_ids[j] = string("#branch", j))
  end
  stars = _three_winding_star_points(T)
  rating(br) = (br.sn_MVA !== nothing && isfinite(br.sn_MVA) && br.sn_MVA > 0.0) ? Float64(br.sn_MVA) : NaN
  return QstsElements(bus_ids, String[nd.comp.cName for nd in T.nodeVec], Float64[getNodeVn(nd) for nd in T.nodeVec], BitVector([!(i in stars) for i in 1:nb]),
                      branch_ids, String[br.comp.cName for br in T.branchVec], branch_kind, Int[br.fromBus for br in T.branchVec], Int[br.toBus for br in T.branchVec],
                      Float64[rating(br) for br in T.branchVec], String[], String[], Symbol[], String[], String[])
end

# The run on one input: checks, validation, base solve, steps through the
# parallel map of the engine, assembly in step order.
function _run_qsts(input::_QstsInput, config::SparlectraConfig, opts::_QstsRunOptions)::QstsResult
  _qsts_check_config(config)
  # a flat start is a flat start on the base case and on every step
  config, flat_overrides = _flatstart_profile_config(config)
  # an event series has no snapshots: each step is a switching state reached
  # from the step before
  isempty(input.event_list) || throw(ArgumentError(string("runQsts: ", input.name, " is an event series (event list ", input.event_list["source"], "); it runs with qsts.mode = sequential only")))
  n = length(input.ops)
  n >= 1 || throw(ArgumentError("runQsts: the series has no step"))
  controls = config.qsts.snapshot_controls
  # the configured controllers go onto T0 once, never onto the template
  # (a shunt controller appends a shunt), so every step has the same set
  T0 = deepcopy(input.template)
  controls && config.control.enabled && !isempty(config.control.controllers) && applyConfiguredControllers!(T0, config.control)
  # the checks of the patch model on T0. A tap_pos op on a transformer with a
  # tap controller is refused only when the controllers run (the controller
  # would move the tap away from it); without them the op is the recorded
  # tap position of the step.
  scenarios = Scenario[Scenario(name = string("step ", k), ops = input.ops[k]) for k in 1:n if !isempty(input.ops[k])]
  _validate_scenarios_net(ScenarioSet(scenarios = scenarios), input.index, T0, controls && config.control.enabled)
  notes = String[]
  _qsts_flatstart_note!(notes, flat_overrides)
  # the base case: the direct run of the template, with the controllers
  # only under qsts.snapshot_controls
  B = deepcopy(T0)
  base_frozen = _qsts_frozen_buses(B)
  # one record for the base case and every step: a tap that starts outside
  # its control range is reported once per run, not once per step
  start_warnings = _TapStartWarnings()
  ex, msg = Base.ScopedValues.with(() -> _qsts_solve!(B, _qsts_base_config(config)), _TAP_START_WARNINGS => start_warnings)
  base_status = _qsts_status(ex)
  # a CGMES union: every step starts from the state of its own snapshot,
  # carried by its delta, whatever the base case gives (no warm start and no
  # fallback to one); an SCF series: warm start from the base solution
  own_start = !isempty(input.start_ops)
  warm = nothing
  if base_status === :solved
    calcNetLosses!(B)
    own_start || (warm = _snapshot_start_voltages(B))
    # the steps of an SCF series start from this solution, so the base case
    # is checked too
    _qsts_plausibility_note!(notes, 0, "", _plausibility_diagnostics(B, config.powerflow.plausibility, true).plausibility_warning)
  else
    push!(notes, string("base case ", base_status === :solver_error ? string("failed: ", msg) : string("did not converge", isempty(msg) ? "" : string(": ", msg)), own_start ? "" : "; the steps start from the start values of the template"))
  end
  own_start && push!(notes, "start voltages: each step from the state of its own snapshot (the SV, carried by the delta)")
  base_iterations, base_outer, base_control = _qsts_iterations(B, ex, controls)
  base_label = input.base_step >= 1 ? input.labels[input.base_step] : "base"
  base_time = input.base_step >= 1 ? input.timestamps[input.base_step] : ""
  base = QstsStepInfo(0, base_label, base_time, true, base_status, base_control, base_iterations, base_outer, 0, length(base_frozen),
                      base_status === :solved ? _qsts_slack_buses(B) : Int[], msg)
  # flatstart stays off whenever the steps have a start of their own
  step_config, config_notes = _qsts_step_config(config, base_status === :solved || own_start)
  append!(notes, config_notes)
  ctx = _QstsContext(input, step_config, warm, controls, _qsts_tap_branches(T0))
  # the parallel switches of the passed configuration, not of the active one
  par = config.runtime.parallel
  outs = Base.ScopedValues.with(_TAP_START_WARNINGS => start_warnings) do
    _engine_parallel_map((tmpl, k) -> _qsts_step(ctx, tmpl, k), () -> deepcopy(T0), n;
      parallel_enabled = something(opts.parallel_enabled, par.enabled), parallel_max_tasks = something(opts.parallel_max_tasks, parallel_max_tasks(par)),
      parallel_min_work_items = something(opts.parallel_min_work_items, par.min_work_items), progress = opts.progress)
  end
  nb = length(T0.nodeVec)
  nbr = length(T0.branchVec)
  bus_values = Array{Float64,3}(undef, nb, length(QSTS_BUS_QUANTITIES), n)
  branch_values = Array{Float64,3}(undef, nbr, length(QSTS_BRANCH_QUANTITIES), n)
  step_values = Matrix{Float64}(undef, length(QSTS_STEP_QUANTITIES), n)
  steps = Vector{QstsStepInfo}(undef, n)
  for k in 1:n
    o = outs[k]::_QstsStepOut
    # a patch never adds or removes buses or branches
    (size(o.bus, 1) == nb && size(o.branch, 1) == nbr) || error("runQsts: step $(k) changed the bus or branch count")
    bus_values[:, :, k] = o.bus
    branch_values[:, :, k] = o.branch
    step_values[:, k] = o.step
    steps[k] = o.info
    isempty(o.isolation_note) || push!(notes, o.isolation_note)
    _qsts_plausibility_note!(notes, k, o.info.label, o.plausibility_warning)
  end
  pf = config.powerflow
  return QstsResult(input.name, input.origin, :snapshots, input.base_step, input.step_minutes, base, steps, _qsts_elements(T0, input.index),
                    bus_values, branch_values, step_values, Array{Float64,3}(undef, 0, length(QSTS_CONTROLLER_QUANTITIES), n), _QstsSwitchRow[], _QstsMetricRow[],
                    pf.max_iter, pf.tol, controls, config.cgmes.tap_control, config.cgmes.map_energy_source, config.qsts, input.event_list, notes)
end

"""
    runQsts(series::SCFSeries; config, parallel_enabled = nothing, parallel_max_tasks = nothing, parallel_min_work_items = nothing, progress = nothing) -> QstsResult
    runQsts(union::CGMESSeriesUnion; config, ...) -> QstsResult

QSTS mode A (`qsts.mode = snapshots`): solve every step of a series on its
own. The input is not changed.

1. Checks, each an `ArgumentError` with the key named: `qsts.mode =
   snapshots`, `power_flow.solver = rectangular`, `power_flow.mode` not
   `auto`, `power_flow.external_grid.enabled = false`; for a CGMES union
   also `cgmes_import.*` equal to the configuration the union was built
   with. The deltas pass [`validate_scenarios`](@ref) against the template,
   but a `tap_pos` op on a transformer with a tap controller (of the case
   or of `control.controllers`) is refused only when the controllers run
   (`qsts.snapshot_controls = true` with `control.enabled`).
2. The base case: the template solved with `config` (for an SCF series the
   base case of the file, for a CGMES union the snapshot of its base step).
3. Every step: a fresh copy of the unsolved template, the full delta of
   the step ([`apply!`](@ref)), one solve. A step of an SCF series starts
   from the bus voltages of the base solution. A step of a CGMES union
   starts from the state of its own snapshot: its delta carries the start
   of every bus the snapshot keeps in service (the SV voltage, without one
   the neighbour start of the importer), whatever the base case gives;
   there is no warm start and no fallback to one. Steps of an SCF series are marked for
   isolation again by the rule of the power flow (`markIsolatedBuses!`
   with `links = true`: a bus behind a closed link is in service with its
   link group), and a part of the net without a source (no unit with
   P > 0 or voltage regulation) is isolated with its injections zeroed,
   as the CGMES import does; one note per such step names the buses and
   the load left unsupplied. A CGMES delta carries the isolation itself,
   and a bus it isolates only for lack of a branch of its own gets the
   values of its link group while the group is solved. An isolated bus in a closed-link
   group with a live member starts at the voltage of the first live member,
   since the solver starts the group from its smallest bus.

Taps and shunt banks stay at the values of each snapshot: base and steps
run with `control.enabled = false`. With `qsts.snapshot_controls = true`
the `control` section of `config` applies in the base case and in every
step. The step configuration switches off the start generators and the
rescue paths (`power_flow.rescue`, `dc.fallback`, `wrong_branch_rescue`),
power mode and the island parallelism, and switches the flat start off
when the steps have a start (a converged base case, a CGMES snapshot);
[`QstsResult`](@ref) `notes` list every key that differs from `config`.

The steps run in parallel over Julia threads under the
`runtime.parallel` switches of `config` (the keywords override them);
parallel and serial runs give identical results. `progress(done, total)`
is called once per finished step.

Element addresses: for a CGMES union the mRID (union topological node, line,
two-winding transformer, transformer end of a three-winding leg,
`<PowerTransformer mRID>#star`, `<mRID>#<class>` for a reused mRID); the
node addresses follow the base step. For an SCF series the component id.

With `qsts.mode = sequential` the SCF method runs QSTS mode B instead: the
steps in order, taps and banks carried from step to step; see
[`runQsts(::AbstractVector{SCFSeries})`](@ref); the CGMES method does the
same.
"""
function runQsts(series::SCFSeries; config::SparlectraConfig, parallel_enabled::Union{Nothing,Bool} = nothing, parallel_max_tasks::Union{Nothing,Integer} = nothing,
                 parallel_min_work_items::Union{Nothing,Integer} = nothing, @nospecialize(progress = nothing))::QstsResult
  config.qsts.mode === :sequential && return only(runQsts([series]; config = config, parallel_enabled = parallel_enabled, parallel_max_tasks = parallel_max_tasks,
                                                          parallel_min_work_items = parallel_min_work_items, progress = progress))
  _qsts_check_config(config)
  opts = _QstsRunOptions(parallel_enabled, parallel_max_tasks === nothing ? nothing : Int(parallel_max_tasks), parallel_min_work_items === nothing ? nothing : Int(parallel_min_work_items), progress)
  return _run_qsts(_qsts_scf_input(series, config), config, opts)
end

# the input of a run on an SCF series: the base case built under `config`,
# the deltas of the file, no SV reference
function _qsts_scf_input(series::SCFSeries, config::SparlectraConfig)::_QstsInput
  n = length(series.steps)
  return _QstsInput(series.name, :scf_series, build_net(series; config = config), series.index, series.ops, Vector{PatchOp}[], String[s.label for s in series.steps],
                    String[s.timestamp for s in series.steps], 0, true, zeros(Int, n), Union{Nothing,QstsSVReference}[nothing for _ in 1:n], series.step_minutes,
                    Dict{String,Any}(get(series.provenance, "event_list", Dict{String,Any}())))
end

# the rows of one step in the order of the long format: step quantities,
# buses by position, branches by position, controllers by position, each
# element's quantities in list order; NaN is no value and gives no row
function _qsts_step_rows(r::QstsResult, k::Int)::Vector{Tuple{Int,String,String,Float64}}
  rows = Tuple{Int,String,String,Float64}[]
  for (q, name) in enumerate(QSTS_STEP_QUANTITIES)
    v = r.step_values[q, k]
    isnan(v) || push!(rows, (k, "", name, v))
  end
  e = r.elements
  for i in eachindex(e.bus_ids)
    e.bus_reported[i] || continue
    for (q, name) in enumerate(QSTS_BUS_QUANTITIES)
      v = r.bus_values[i, q, k]
      isnan(v) || push!(rows, (k, e.bus_ids[i], name, v))
    end
  end
  for j in eachindex(e.branch_ids)
    for (q, name) in enumerate(QSTS_BRANCH_QUANTITIES)
      v = r.branch_values[j, q, k]
      isnan(v) || push!(rows, (k, e.branch_ids[j], name, v))
    end
  end
  for c in eachindex(e.controller_ids)
    for (q, name) in enumerate(QSTS_CONTROLLER_QUANTITIES)
      v = r.controller_values[c, q, k]
      isnan(v) || push!(rows, (k, e.controller_ids[c], name, v))
    end
  end
  return rows
end

"""
    qstsLongTable(r::QstsResult) -> (step, element_id, quantity, value)

The long format of a snapshot run as a named tuple of four equally long
column vectors (`Vector{Int}`, `Vector{String}`, `Vector{String}`,
`Vector{Float64}`). `step` is the position in the series (1 to n),
`element_id` the address of the bus, branch or controller (empty for a step
quantity), `quantity` a name from [`QSTS_STEP_QUANTITIES`](@ref),
[`QSTS_BUS_QUANTITIES`](@ref), [`QSTS_BRANCH_QUANTITIES`](@ref) or
[`QSTS_CONTROLLER_QUANTITIES`](@ref), the unit in its name. A row exists
only for a solved value. The rows run by step, then step quantities, buses,
branches and controllers by position, each element's quantities in list
order. A named tuple of vectors is a Tables.jl column
table, so `DataFrame(qstsLongTable(r))` works without a dependency in
Sparlectra.
"""
function qstsLongTable(r::QstsResult)
  step = Int[]
  element_id = String[]
  quantity = String[]
  value = Float64[]
  for k in eachindex(r.steps), row in _qsts_step_rows(r, k)
    push!(step, row[1])
    push!(element_id, row[2])
    push!(quantity, row[3])
    push!(value, row[4])
  end
  return (step = step, element_id = element_id, quantity = quantity, value = value)
end

"""
    qstsValues(r::QstsResult, element_id, quantity) -> Vector{Float64}

The series of one quantity of one element over all steps, `NaN` where the
long table has no row (not solved, out of service, isolated, no rating).
An empty `element_id` selects a step quantity. Throws `ArgumentError` for
an unknown element or a quantity that does not belong to its kind.
"""
function qstsValues(r::QstsResult, element_id::AbstractString, quantity::AbstractString)::Vector{Float64}
  isempty(element_id) && return r.step_values[_qsts_quantity(QSTS_STEP_QUANTITIES, quantity, "step"), :]
  e = r.elements
  i = findfirst(==(element_id), e.bus_ids)
  if i !== nothing
    q = _qsts_quantity(QSTS_BUS_QUANTITIES, quantity, "bus")
    return e.bus_reported[i] ? r.bus_values[i, q, :] : fill(NaN, length(r.steps))
  end
  j = findfirst(==(element_id), e.branch_ids)
  j === nothing || return r.branch_values[j, _qsts_quantity(QSTS_BRANCH_QUANTITIES, quantity, "branch"), :]
  c = findfirst(==(element_id), e.controller_ids)
  c === nothing && throw(ArgumentError(string("qstsValues: unknown element id ", repr(element_id))))
  return r.controller_values[c, _qsts_quantity(QSTS_CONTROLLER_QUANTITIES, quantity, "controller"), :]
end

# the position of `quantity` in `list`, or an ArgumentError naming the list
function _qsts_quantity(list::Vector{String}, quantity::AbstractString, kind::String)::Int
  q = findfirst(==(quantity), list)
  q === nothing && throw(ArgumentError(string("qstsValues: ", repr(quantity), " is no ", kind, " quantity; known: ", join(list, ", "))))
  return q
end

# the step length for the manifest: whole minutes as an integer
_qsts_json_minutes(m::Union{Nothing,Float64}) = m === nothing ? nothing : (isinteger(m) ? Int(m) : m)

"""
    writeQstsResults(dir, r::QstsResult; format = result_csv_format()) -> Vector{String}

Write a QSTS run into `dir` (created when missing) and return the paths, four
for a snapshot run:

- `qsts_values.csv`: the long format of [`qstsLongTable`](@ref), columns
  `step, element_id, quantity, value`, written one step at a time;
- `qsts_steps.csv`: `step, label, timestamp, is_base, status,
  control_status, iterations, outer_iterations, residual, n_frozen,
  reference, message`, with `reference` the addresses of the slack buses
  joined by `|`;
- `qsts_elements.csv`: `element_id, kind, name, from, to, vn_kv,
  rating_mva` for every bus (kind `bus`, or `star_point` for the star point
  of a three-winding transformer), branch (`line`, `transformer`) and, in a
  sequential run, controller (`controller`, with the actuator in `from` and
  the regulated bus in `to`);
- `qsts_manifest.json`: format `sparlectra-qsts-result`, `format_version`
  1.0, name, origin, mode, number of steps, base step, step length, the
  base case, `power_flow.max_iter` and `tol`, `qsts.snapshot_controls`,
  `cgmes_import.tap_control` and `map_energy_source`, and the notes; a
  sequential run adds the `qsts` keys of mode B and `n_failed`.

A sequential run writes two files more:

- `qsts_switching.csv`: [`qstsSwitching`](@ref), columns `step, element_id,
  old, new, steps, trigger`;
- `qsts_metrics.csv`: [`qstsMetrics`](@ref), columns `element_id, metric,
  value`.

The CSV files follow `format` (`output.csv_format`); numbers are written
with 15 significant digits. The content holds no times, so two runs of the
same input give the same bytes.
"""
writeQstsResults(dir::AbstractString, r::QstsResult; format = result_csv_format())::Vector{String} = _write_qsts_results(dir, r, format, Dict{String,Any}())

"""
    writeQstsResults(dir, rs::AbstractVector{QstsResult}; format = result_csv_format()) -> Vector{String}

Write the results of a run on several series, each into `dir/series_<k>`
as [`writeQstsResults`](@ref) writes one result; each manifest adds
`series` (k) and `n_series`. Returns all paths in series order.
"""
function writeQstsResults(dir::AbstractString, rs::AbstractVector{QstsResult}; format = result_csv_format())::Vector{String}
  n = length(rs)
  return reduce(vcat, (_write_qsts_results(joinpath(dir, string("series_", k)), rs[k], format, Dict{String,Any}("series" => k, "n_series" => n)) for k in 1:n); init = String[])
end

# the writer of one result; `extra` goes into the manifest (the series
# position of a vector run, nothing for a single result)
function _write_qsts_results(dir::AbstractString, r::QstsResult, format, extra::Dict{String,Any})::Vector{String}
  mkpath(dir)
  e = r.elements
  values_path = joinpath(dir, "qsts_values.csv")
  rows = (row for k in eachindex(r.steps) for row in _qsts_step_rows(r, k))
  write_result_csv(values_path, _QSTS_VALUES_HEADER, rows; format = format)
  steps_path = joinpath(dir, "qsts_steps.csv")
  step_rows = ((s.step, s.label, s.timestamp, s.is_base ? 1 : 0, String(s.status), String(s.control_status), s.iterations, s.outer_iterations, s.residual, s.n_frozen,
                join((e.bus_ids[b] for b in s.reference), "|"), s.message) for s in r.steps)
  write_result_csv(steps_path, _QSTS_STEPS_HEADER, step_rows; format = format)
  elements_path = joinpath(dir, "qsts_elements.csv")
  nan_free(x::Float64) = isnan(x) ? nothing : x
  bus_rows = ((e.bus_ids[i], e.bus_reported[i] ? "bus" : "star_point", e.bus_names[i], nothing, nothing, e.bus_vn_kv[i], nothing) for i in eachindex(e.bus_ids))
  branch_rows = ((e.branch_ids[j], String(e.branch_kind[j]), e.branch_names[j], e.bus_ids[e.branch_from[j]], e.bus_ids[e.branch_to[j]], nothing, nan_free(e.branch_rating_mva[j]))
                 for j in eachindex(e.branch_ids))
  controller_rows = ((e.controller_ids[c], "controller", e.controller_names[c], e.controller_actuator[c], isempty(e.controller_target[c]) ? nothing : e.controller_target[c], nothing, nothing)
                     for c in eachindex(e.controller_ids))
  write_result_csv(elements_path, _QSTS_ELEMENTS_HEADER, Iterators.flatten((bus_rows, branch_rows, controller_rows)); format = format)
  manifest = Dict{String,Any}(
    "format" => _QSTS_RESULT_FORMAT,
    "format_version" => _QSTS_RESULT_FORMAT_VERSION,
    "name" => r.name,
    "origin" => String(r.origin),
    "mode" => String(r.mode),
    "n_steps" => length(r.steps),
    "base_step" => r.base_step,
    "step_minutes" => _qsts_json_minutes(r.step_minutes),
    "base" => Dict{String,Any}("converged" => r.base.status === :solved, "status" => String(r.base.status), "iterations" => r.base.iterations,
                               "outer_iterations" => r.base.outer_iterations, "control_status" => String(r.base.control_status), "message" => r.base.message),
    "power_flow" => Dict{String,Any}("max_iter" => r.max_iter, "tol" => r.tol),
    "qsts" => Dict{String,Any}("snapshot_controls" => r.snapshot_controls),
    "cgmes_import" => Dict{String,Any}("tap_control" => r.tap_control, "map_energy_source" => r.map_energy_source),
    "notes" => r.notes,
  )
  paths = [values_path, steps_path, elements_path]
  if r.mode === :sequential
    q = r.qsts
    merge!(manifest["qsts"], Dict{String,Any}("sequential_start" => String(q.sequential_start), "resync_every" => q.resync_every, "dead_time_steps" => q.dead_time_steps,
                                              "hunting_window" => q.hunting_window, "vm_min_pu" => q.vm_min_pu, "vm_max_pu" => q.vm_max_pu, "loading_max_pct" => q.loading_max_pct))
    manifest["n_failed"] = count(st -> st.status !== :solved, r.steps)
    isempty(r.event_list) || (manifest["event_list"] = r.event_list)
    switching_path = joinpath(dir, "qsts_switching.csv")
    write_result_csv(switching_path, _QSTS_SWITCHING_HEADER, ((w.step, w.element_id, w.old, w.new, w.steps, String(w.trigger)) for w in r.switching); format = format)
    metrics_path = joinpath(dir, "qsts_metrics.csv")
    write_result_csv(metrics_path, _QSTS_METRICS_HEADER, ((m.element_id, m.metric, m.value) for m in r.metrics); format = format)
    append!(paths, [switching_path, metrics_path])
  end
  merge!(manifest, extra)
  manifest_path = joinpath(dir, "qsts_manifest.json")
  write(manifest_path, scf_json_string(manifest))
  push!(paths, manifest_path)
  return paths
end
