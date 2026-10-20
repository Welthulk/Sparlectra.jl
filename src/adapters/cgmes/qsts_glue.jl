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

# file: src/adapters/cgmes/qsts_glue.jl
# purpose: the QSTS entry points on CGMES adapter types: runQsts on a series
#          union (mode A and B), on several unions and on an event series,
#          and qstsEventSeries on a union. Adapter glue like glue.jl, kept
#          apart so its docstrings go to the QSTS reference page.

"""
    runQsts(union::CGMESImporter.CGMESSeriesUnion; config, parallel_enabled = nothing, parallel_max_tasks = nothing, parallel_min_work_items = nothing, progress = nothing) -> QstsResult

QSTS on the union network of a CGMES snapshot series; see the SCF method of
`runQsts` for the run itself. Step `k` applies the full delta of snapshot
`k` (`ops` and `control_ops`: banks take the SSH sections, taps the SSH
steps with `cgmes_import.tap_control = true` and the `SvTapStep` positions
with the default `false`), the base case is the snapshot of
`union.base_step`, and every step is compared with the SV state of its
snapshot (`union.sv`). With `qsts.mode = snapshots` each step starts from
the state of its snapshot (`start_ops` of its delta, see
[`CGMESImporter.CGMESSnapshotDelta`](@ref)), not from the base solution.
`config.cgmes` must equal `union.import_config` field by field, or
template, delta and run would not fit together; an `ArgumentError` names
the keys that differ.

With `qsts.mode = snapshots` every step is solved on its own (mode A). With
`qsts.mode = sequential` the steps run in order with the tap and bank state
carried from step to step (mode B, see
[`runQsts(::AbstractVector{SCFSeries})`](@ref)); the taps and banks of each
snapshot are then the reference positions of its step.
"""
function runQsts(u::CGMESImporter.CGMESSeriesUnion; config::SparlectraConfig, parallel_enabled::Union{Nothing,Bool} = nothing, parallel_max_tasks::Union{Nothing,Integer} = nothing,
                 parallel_min_work_items::Union{Nothing,Integer} = nothing, @nospecialize(progress = nothing))::QstsResult
  config.qsts.mode === :sequential && return only(runQsts([u]; config = config, parallel_enabled = parallel_enabled, parallel_max_tasks = parallel_max_tasks,
                                                          parallel_min_work_items = parallel_min_work_items, progress = progress))
  _qsts_check_config(config)
  _qsts_check_union_config(u, config)
  opts = _QstsRunOptions(parallel_enabled, parallel_max_tasks === nothing ? nothing : Int(parallel_max_tasks), parallel_min_work_items === nothing ? nothing : Int(parallel_min_work_items), progress)
  return _run_qsts(CGMESImporter._qsts_input(u), config, opts)
end

"""
    runQsts(unions::AbstractVector{CGMESImporter.CGMESSeriesUnion}; config, parallel_enabled = nothing, parallel_max_tasks = nothing, parallel_min_work_items = nothing, progress = nothing) -> Vector{QstsResult}

QSTS on several CGMES series unions, one result per union: one snapshot run
after the other with `qsts.mode = snapshots`, the unions in parallel and the
steps of each in order with `qsts.mode = sequential`, as for
[`runQsts(::AbstractVector{SCFSeries})`](@ref). Every union must have been
built with the `cgmes_import` settings of `config`.
"""
function runQsts(us::AbstractVector{CGMESImporter.CGMESSeriesUnion}; config::SparlectraConfig, parallel_enabled::Union{Nothing,Bool} = nothing, parallel_max_tasks::Union{Nothing,Integer} = nothing,
                 parallel_min_work_items::Union{Nothing,Integer} = nothing, @nospecialize(progress = nothing))::Vector{QstsResult}
  if config.qsts.mode !== :sequential
    return QstsResult[runQsts(u; config = config, parallel_enabled = parallel_enabled, parallel_max_tasks = parallel_max_tasks, parallel_min_work_items = parallel_min_work_items,
                              progress = progress) for u in us]
  end
  _qsts_check_config(config, :sequential)
  foreach(u -> _qsts_check_union_config(u, config), us)
  return _run_qsts_sequential(Sparlectra._QstsInput[CGMESImporter._qsts_input(u) for u in us], config, parallel_enabled, parallel_max_tasks, parallel_min_work_items, progress)
end

# the union, its delta and a run fit together only under the cgmes_import
# settings the union was built with
function _qsts_check_union_config(u::CGMESImporter.CGMESSeriesUnion, config::SparlectraConfig)
  differ = [string("cgmes_import.", f) for f in fieldnames(CGMESImportConfig) if !isequal(getfield(config.cgmes, f), getfield(u.import_config, f))]
  isempty(differ) || throw(ArgumentError(string("runQsts: ", join(differ, ", "), length(differ) == 1 ? " differs" : " differ",
                                                " from the configuration the union was built with; build the union with the same cgmes_import settings")))
  return nothing
end

"""
    CGMESEventSeries

An event series on the union network of a CGMES snapshot series, built by
[`qstsEventSeries`](@ref): `union` (its template in the state of the base
step is the base of every step), `ops` (the full delta of each step),
`labels` (`"1"` to `"n"`) and `event_list` (source file and row count).
The steps have no time reference and no SV. Run it with [`runQsts`](@ref)
and `qsts.mode = sequential`.
"""
struct CGMESEventSeries
  union::CGMESImporter.CGMESSeriesUnion
  ops::Vector{Vector{PatchOp}}
  labels::Vector{String}
  event_list::Dict{String,Any}
end

function qstsEventSeries(u::CGMESImporter.CGMESSeriesUnion, events::QstsEventList; n_steps::Integer = maximum(events.step; init = 0))::CGMESEventSeries
  n = Int(n_steps)
  n >= 1 || throw(ArgumentError("qstsEventSeries: the event list has no event and n_steps is not given"))
  return CGMESEventSeries(u, _qsts_event_ops(u.index, u.template, events, n), String[string(k) for k in 1:n], _qsts_event_record(events))
end

"""
    runQsts(series::CGMESEventSeries; config, parallel_enabled = nothing, parallel_max_tasks = nothing, parallel_min_work_items = nothing, progress = nothing) -> QstsResult

The sequential run (`qsts.mode = sequential`) of an event series on a CGMES
union; see [`runQsts(::AbstractVector{SCFSeries})`](@ref). Like a step of an
SCF series, each step is marked for isolation by the rule of the power flow
and an island without a source is isolated. `config.cgmes` must equal the
configuration the union was built with.
"""
function runQsts(e::CGMESEventSeries; config::SparlectraConfig, parallel_enabled::Union{Nothing,Bool} = nothing, parallel_max_tasks::Union{Nothing,Integer} = nothing,
                 parallel_min_work_items::Union{Nothing,Integer} = nothing, @nospecialize(progress = nothing))::QstsResult
  u = e.union
  n = length(e.ops)
  input = _QstsInput(string("events ", e.event_list["source"], " on the union of ", u.labels[u.base_step]), :cgmes_series, u.template, u.index, e.ops, Vector{PatchOp}[], copy(e.labels), fill("", n), 0, true,
                     zeros(Int, n), Union{Nothing,QstsSVReference}[nothing for _ in 1:n], nothing, e.event_list)
  if config.qsts.mode !== :sequential
    # the refusal of a snapshot run names the key
    _qsts_check_config(config)
    throw(ArgumentError(string("runQsts: ", input.name, " is an event series (event list ", e.event_list["source"], "); it runs with qsts.mode = sequential only")))
  end
  _qsts_check_config(config, :sequential)
  _qsts_check_union_config(u, config)
  return only(_run_qsts_sequential([input], config, parallel_enabled, parallel_max_tasks, parallel_min_work_items, progress))
end
