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

# file: src/qsts/qsts_events.jl
# purpose: the event list of a sequential QSTS run: a CSV
#          `step,element_id,action,value` with open/close on branches,
#          transformers and switches and in/out on units, loads and shunts,
#          read with its line in every error, and turned into a series whose
#          step k holds the full delta of the switching state after the
#          events up to k. The SCF method lives here, the CGMES method in
#          src/adapters/cgmes/glue.jl.

const _QSTS_EVENTS_HEADER = "step,element_id,action,value"
const _QSTS_EVENT_ACTIONS = [:open, :close, :in, :out]

"""
    QstsEventList

The events of a sequential QSTS run, one entry per CSV row in file order:
`step` (from 1), `element_id` (the address of the element: the component
id of an SCF case, the mRID of a CGMES union element), `action` (`:open`,
`:close`, `:in`, `:out`) and `line` (the CSV line, named in later errors).
`source` is the file name without its path. Read with
[`readQstsEvents`](@ref), turned into a series with
[`qstsEventSeries`](@ref).
"""
struct QstsEventList
  source::String
  step::Vector{Int}
  element_id::Vector{String}
  action::Vector{Symbol}
  line::Vector{Int}
end

_qsts_event_error(source::AbstractString, line::Int, msg::AbstractString) = throw(ArgumentError(string("event list ", source, " line ", line, ": ", msg)))

"""
    readQstsEvents(path) -> QstsEventList

Read an event list. The file is UTF-8 CSV with comma and LF (a byte order
mark is dropped); the first line is exactly `step,element_id,action,value`,
every other non-empty line one event:

| Column | Rule |
|---|---|
| `step` | whole number from 1 |
| `element_id` | not empty; resolved against the series by [`qstsEventSeries`](@ref) |
| `action` | `open`, `close` (branch, transformer, switch) or `in`, `out` (unit, load, shunt) |
| `value` | empty |

An element takes one event per step. A state holds until the next event of
the same element. Every error is an `ArgumentError` naming the file and the
line.
"""
function readQstsEvents(path::AbstractString)::QstsEventList
  source = basename(String(path))
  text = read(path, String)
  startswith(text, "﻿") && (text = text[nextind(text, 1):end])
  lines = split(text, '\n')
  header = isempty(lines) ? "" : rstrip(lines[1], '\r')
  header == _QSTS_EVENTS_HEADER || _qsts_event_error(source, 1, string("the header must be ", _QSTS_EVENTS_HEADER, "; got ", repr(header)))
  step = Int[]
  element_id = String[]
  action = Symbol[]
  line = Int[]
  first_at = Dict{Tuple{Int,String},Int}()
  for (i, raw) in enumerate(lines)
    i == 1 && continue
    l = rstrip(raw, '\r')
    isempty(strip(l)) && continue
    f = split(l, ',')
    length(f) == 4 || _qsts_event_error(source, i, string("expected 4 fields (", _QSTS_EVENTS_HEADER, "), got ", length(f)))
    k = tryparse(Int, strip(f[1]))
    (k !== nothing && k >= 1) || _qsts_event_error(source, i, string("step must be a whole number from 1; got ", repr(f[1])))
    id = String(strip(f[2]))
    isempty(id) && _qsts_event_error(source, i, "element_id is empty")
    a = Symbol(strip(f[3]))
    a in _QSTS_EVENT_ACTIONS || _qsts_event_error(source, i, string("action must be open, close, in or out; got ", repr(f[3])))
    isempty(strip(f[4])) || _qsts_event_error(source, i, string("value must be empty for action ", a, "; got ", repr(f[4])))
    at = get(first_at, (k, id), 0)
    at == 0 || _qsts_event_error(source, i, string("second event of ", id, " in step ", k, "; the first is at line ", at))
    first_at[(k, id)] = i
    push!(step, k)
    push!(element_id, id)
    push!(action, a)
    push!(line, i)
  end
  return QstsEventList(source, step, element_id, action, line)
end

# whether the element is in service (closed) in `net`, the base state of an
# event series; parked units count as out
function _qsts_event_base_on(net::Net, index::ScenarioIndex, id::Int, kind::Symbol, internal::Int)::Bool
  (kind === :branch || kind === :transformer) && return net.branchVec[internal].status == 1
  kind === :link && return net.linkVec[internal].status == 1
  kind === :shunt && return net.shuntVec[internal].status == 1
  return !haskey(index.parked_values, id)
end

# the ops that put one element from its base state into state `on`: branches
# of a CGMES union take the three status fields of the union vocabulary (the
# SCF status op would also drop MATPOWER bus-shunt parts), SCF branches,
# links and shunts the status op, units and loads park and unpark
function _qsts_event_element_ops(index::ScenarioIndex, id::Int, kind::Symbol, on::Bool)::Vector{PatchOp}
  v = on ? 1.0 : 0.0
  if kind === :branch || kind === :transformer
    index.origin === :cgmes_series && return [PatchOp(op = :set, target = kind, id = id, field = f, value = v) for f in (:status, :from_status, :to_status)]
    return [PatchOp(op = :status, target = kind, id = id, value = v)]
  end
  (kind === :link || kind === :shunt) && return [PatchOp(op = :status, target = kind, id = id, value = v)]
  return [PatchOp(op = :set, target = kind, id = id, field = :in_service, value = v)]
end

# The full delta of every step 1..n_steps of an event series against `net`
# (the template in its base state): the state of each element after its
# last event up to the step, an op only where that state differs from the
# base. An element is addressed by mRID (union) or by its component id
# (SCF); open/close needs a branch, transformer or link, in/out a unit,
# load or shunt.
function _qsts_event_ops(index::ScenarioIndex, net::Net, events::QstsEventList, n_steps::Int)::Vector{Vector{PatchOp}}
  src = events.source
  resolved = Tuple{Int,Symbol,Int}[]
  for i in eachindex(events.step)
    eid = events.element_id[i]
    id = get(index.id_by_mrid, eid, nothing)
    if id === nothing
      n = tryparse(Int, eid)
      (n !== nothing && haskey(index.kind_by_id, n) && index.origin !== :cgmes_series) && (id = n)
    end
    id === nothing && _qsts_event_error(src, events.line[i], string("unknown element_id ", repr(eid),
                                                                      index.origin === :cgmes_series ? " (a switch must be a link of the union, see extra_link_switches)" : ""))
    kind = index.kind_by_id[id]
    a = events.action[i]
    switchable = kind in (:branch, :transformer, :link)
    unit = kind in (:generator, :load, :external_grid, :shunt)
    ((a in (:open, :close) && switchable) || (a in (:in, :out) && unit)) ||
      _qsts_event_error(src, events.line[i], string("action ", a, " does not apply to ", eid, ", a ", kind, a in (:open, :close) ? "; open and close need a branch, transformer or switch" : "; in and out need a unit, load or shunt"))
    events.step[i] <= n_steps || _qsts_event_error(src, events.line[i], string("step ", events.step[i], " lies beyond n_steps = ", n_steps))
    push!(resolved, (id, kind, index.internal_by_id[id]))
  end
  order = sortperm(eachindex(events.step); by = i -> (events.step[i], events.line[i]))
  ids = unique(first(r) for r in resolved)
  base_on = Dict(r[1] => _qsts_event_base_on(net, index, r[1], r[2], r[3]) for r in resolved)
  kind_of = Dict(r[1] => r[2] for r in resolved)
  state = copy(base_on)
  ops = Vector{Vector{PatchOp}}(undef, n_steps)
  j = 1
  for k in 1:n_steps
    while j <= length(order) && events.step[order[j]] <= k
      i = order[j]
      state[resolved[i][1]] = events.action[i] in (:close, :in)
      j += 1
    end
    ops[k] = PatchOp[op for id in sort(ids) if state[id] != base_on[id] for op in _qsts_event_element_ops(index, id, kind_of[id], state[id])]
  end
  return ops
end

"""
    qstsEventSeries(base::SCFCase, events::QstsEventList; n_steps = largest event step, name) -> SCFSeries
    qstsEventSeries(union::CGMESSeriesUnion, events::QstsEventList; n_steps = largest event step) -> CGMESEventSeries

A series driven by an event list instead of snapshot deltas: step k holds
the full delta of the switching state after all events up to step k
(an element keeps its state until its next event; an event that restores
the base state removes the op). Nothing else changes between steps, so a
sequential run ([`runQsts`](@ref) with `qsts.mode = sequential`) shows the
voltage jumps, flow shifts, violations and islands of the switching
sequence. The SCF series has no time reference and carries the source and
the row count under `provenance["event_list"]`, which `write_scf_series`
keeps. A snapshot run refuses an event series, and so does
`qsts.sequential_start = snapshot`. Errors name the CSV line of the event.
"""
function qstsEventSeries(base::SCFCase, events::QstsEventList; n_steps::Integer = maximum(events.step; init = 0), name::AbstractString = string("events ", events.source))::SCFSeries
  n = Int(n_steps)
  n >= 1 || throw(ArgumentError("qstsEventSeries: the event list has no event and n_steps is not given"))
  steps = [SCFSeriesStep(k, "", string(k)) for k in 1:n]
  probe = SCFSeries(name = name, base = base, steps = steps, ops = [PatchOp[] for _ in 1:n])
  ops = _qsts_event_ops(probe.index, build_net(probe), events, n)
  return SCFSeries(name = name, base = base, steps = steps, ops = ops, provenance = Dict{String,Any}("event_list" => _qsts_event_record(events)))
end

# what a run records about its event list: the file name and the row count
_qsts_event_record(events::QstsEventList)::Dict{String,Any} = Dict{String,Any}("source" => events.source, "rows" => length(events.step))
