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

# file: src/adapters/scf/scf_series.jl
# purpose: the SCF time series: one zip with a manifest (series.json), an
#          unchanged SCF base case (base.scf.json) and one table of patch
#          operations for all steps (delta.csv). Each step is the full
#          delta against the base case, so a step applies on its own. The
#          series index has origin :scf_series and takes the QSTS
#          vocabulary; an appliance the base case marks with status 0 is
#          parked when the series is loaded (the plain SCF reader builds it
#          in service, and an N-1 status 0 stays a deletion). The file
#          lives with the scenario code, not with the SCF adapter, because
#          it needs the patch model, the index and apply!, which are
#          included after the adapter.

const _SCF_SERIES_FORMAT = "sparlectra-scf-series"
const _SCF_SERIES_FORMAT_VERSION = "1.0"
# the entries of a series zip in the order the writer puts them; any other
# entry is an error, as an unknown key is in a case file
const _SCF_SERIES_ENTRIES = ["series.json", "base.scf.json", "delta.csv"]
const _SCF_SERIES_MANIFEST_KEYS = ["created_by", "format", "format_version", "name", "provenance", "step_minutes", "steps"]
const _SCF_SERIES_STEP_KEYS = ["label", "step", "timestamp"]
const _SCF_SERIES_HEADER = "step,op,target,id,field,value"
# the targets a delta row may name; a tap controller has no id in an SCF
# case, so a series cannot address one
const _SCF_SERIES_TARGETS = (:load, :generator, :external_grid, :branch, :transformer, :shunt, :link, :bus)
# ISO 8601 UTC to the minute or second, as in the measurement table but
# without fractions of a second
const _SCF_SERIES_TIMESTAMP = r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}(:\d{2})?Z$"
# bytes at the end of a file the format probe reads: the end record of a
# zip (22 bytes plus a comment of at most 65535) and the directory in front
const _SCF_SERIES_PROBE_BYTES = 65_557 + 16_384

"""
    SCFSeriesStep(step, timestamp, label)

One step of an SCF time series: the step number (steps run from 1 without
gaps, the base case is no step), the time `timestamp` in ISO 8601 UTC
(`YYYY-MM-DDThh:mm[:ss]Z`, empty without time reference) and the axis
`label`.
"""
struct SCFSeriesStep
  step::Int
  timestamp::String
  label::String
end

"""
    SCFSeries

An SCF time series in memory: a base case and one list of patch operations
per step.

- `name`, `created_by`, `provenance`: from the manifest.
- `base`: the base case ([`SCFCase`](@ref)).
- `index`: the [`ScenarioIndex`](@ref) of the base case with origin
  `:scf_series`. It takes the QSTS vocabulary of [`PatchOp`](@ref), maps
  the node ids to buses and holds the values of the parked units.
- `parked`: ids of the appliances (`sym_load`, `sym_gen`, `source`) with
  status 0 in the base case. `build_net(series)` parks them with
  `set in_service 0`; a step brings one back with `set in_service 1`.
- `steps`: one [`SCFSeriesStep`](@ref) per step.
- `ops`: the patch operations of each step, the full delta against the
  base case (`scale` refers to the base value). An empty list is the base
  state.
- `step_minutes`: the step length in minutes, `nothing` for irregular steps.

The keyword constructor
`SCFSeries(; name, base, steps, ops, step_minutes = nothing, provenance = Dict(), created_by = "Sparlectra <version>")`
runs the checks of [`read_scf_series`](@ref) and names step and op of a
failure.
"""
struct SCFSeries
  name::String
  created_by::String
  base::SCFCase
  index::ScenarioIndex
  parked::Vector{Int}
  steps::Vector{SCFSeriesStep}
  ops::Vector{Vector{PatchOp}}
  step_minutes::Union{Nothing,Float64}
  provenance::Dict{String,Any}
end

_scf_series_error(at::AbstractString, msg::AbstractString) = throw(ArgumentError(string("SCF series: ", at, ": ", msg)))

# the series from its parts, with every check of the reader; `at(i)` of a
# step names op i as "step k op i"
function SCFSeries(; name::AbstractString, base::SCFCase, steps::Vector{SCFSeriesStep}, ops::Vector{Vector{PatchOp}},
                   step_minutes::Union{Nothing,Real} = nothing, provenance::AbstractDict = Dict{String,Any}(),
                   created_by::AbstractString = string("Sparlectra ", version()))
  isempty(strip(name)) && _scf_series_error("name", "a series needs a name")
  minutes = step_minutes === nothing ? nothing : Float64(step_minutes)
  _scf_series_check_steps(steps, minutes)
  length(ops) == length(steps) || _scf_series_error("ops", string(length(ops), " op lists for ", length(steps), " steps"))
  index, parked = _scf_series_index(base)
  parked_set = Set(parked)
  for k in eachindex(ops)
    _scf_series_check_step(ops[k], k, index, parked_set, i -> string("step ", k, " op ", i))
  end
  return SCFSeries(String(name), String(created_by), base, index, parked, copy(steps), [copy(o) for o in ops], minutes, Dict{String,Any}(provenance))
end

# the step list rules: numbers 1 to n in order, a timestamp on every step or
# on none, times strictly increasing, a positive step length or none
function _scf_series_check_steps(steps::Vector{SCFSeriesStep}, minutes::Union{Nothing,Float64})
  isempty(steps) && _scf_series_error("series.json", "a series needs at least one step")
  (minutes === nothing || (isfinite(minutes) && minutes > 0.0)) || _scf_series_error("series.json", string("step_minutes must be a positive number or null, got ", minutes))
  timed = !isempty(steps[1].timestamp)
  previous = nothing
  for (k, s) in enumerate(steps)
    at = string("series.json step ", k)
    s.step == k || _scf_series_error(at, string("has the number ", s.step, "; the steps run from 1 without gaps"))
    isempty(s.timestamp) == timed && _scf_series_error(at, "either every step has a timestamp or none")
    timed || continue
    ts = s.timestamp
    t = occursin(_SCF_SERIES_TIMESTAMP, ts) ? tryparse(DateTime, ts[1:(end - 1)], length(ts) == 17 ? dateformat"yyyy-mm-ddTHH:MM" : dateformat"yyyy-mm-ddTHH:MM:SS") : nothing
    t === nothing && _scf_series_error(at, string("timestamp ", repr(ts), " is not a time in ISO 8601 UTC (YYYY-MM-DDThh:mm[:ss]Z)"))
    (previous === nothing || t > previous) || _scf_series_error(at, string("timestamp ", ts, " is not after the one of step ", k - 1))
    previous = t
  end
  return nothing
end

# The series index of a base case and the ids it parks. Every component needs
# its build-order field in the extra block (element_index, branch_index,
# shunt_index, bus_index), or a patch could not find it in the net. The
# values of a parked unit come from a build of the base case: build_net takes
# only net-level parameters from the configuration, so the prosumers are the
# same under every configuration.
function _scf_series_index(base::SCFCase)::Tuple{ScenarioIndex,Vector{Int}}
  index = ScenarioIndex(base; origin = :scf_series)
  no_order = sort!([id for id in keys(index.kind_by_id) if !haskey(index.internal_by_id, id)])
  isempty(no_order) || _scf_series_error("base.scf.json", string("no build order (extra.element_index, branch_index, shunt_index or bus_index) for id ",
    join(first(no_order, 10), ", "), length(no_order) > 10 ? ", ..." : "", "; export the case again"))
  d = base.data
  parked = sort!([r.id for rows in (d.source, d.sym_gen, d.sym_load) for r in rows if something(r.status, 1) == 0])
  if !isempty(parked)
    net = build_net(base)
    for id in parked
      index.parked_values[id] = deepcopy(net.prosumpsVec[index.internal_by_id[id]])
    end
  end
  return index, parked
end

# the fields one op writes, for the one-write-per-field rule of a step; a
# scale without field writes p and q
function _scf_series_written(op::PatchOp)::Vector{Symbol}
  op.op === :status && return [:status]
  op.op === :scale && op.field === nothing && return [:p, :q]
  return [something(op.field, :none)]
end

# The checks of one step: every op through the per-op checks of
# validate_scenarios, one write per id and field (set before scale gives
# another value than scale before set), and the parking rules. apply! parks
# and unparks before every other op, so a write on a unit that stays parked
# would land on the 0 MW load (a set p turns a generator into a load of that
# size), and status 1 does not unpark. Hence: a unit the step parks with
# `set in_service 0` takes no other op in the step, and a unit the base case
# parks takes a set, scale or status 1 only when the step brings it back with
# `set in_service 1`. `at(i)` names op i.
function _scf_series_check_step(ops::Vector{PatchOp}, k::Int, index::ScenarioIndex, parked::Set{Int}, at)
  first_write = Dict{Tuple{Int,Symbol},Int}()
  unparked = Set{Int}(op.id for op in ops if op.op === :set && op.field === :in_service && op.value == 1.0)
  # the first op of the step that parks an id; a second one is a second
  # write of in_service and fails on that rule
  parking = Dict{Int,Int}()
  for (i, op) in enumerate(ops)
    (op.op === :set && op.field === :in_service && op.value == 0.0) && get!(parking, op.id, i)
  end
  for (i, op) in enumerate(ops)
    msg = _patch_op_class_error(op, index)
    msg === nothing && (msg = _patch_op_rule_error(op, index))
    msg === nothing || _scf_series_error(at(i), msg)
    for f in _scf_series_written(op)
      j = get(first_write, (op.id, f), 0)
      j == 0 || _scf_series_error(at(i), string("id ", op.id, " gets a second write of ", f, " in step ", k, "; the first is at ", at(j)))
      first_write[(op.id, f)] = i
    end
    j = get(parking, op.id, i)
    j == i || _scf_series_error(at(i), string("id ", op.id, " is parked in this step by set in_service 0 at ", at(j), "; a parked unit takes no other op in the same step"))
    (op.id in parked && !(op.id in unparked)) || continue
    if op.op === :scale || (op.op === :set && op.field !== :in_service)
      _scf_series_error(at(i), string("id ", op.id, " is parked in the base case (status 0); a ", op.op, " on it needs set in_service 1 in the same step"))
    elseif op.op === :status && op.value == 1.0
      _scf_series_error(at(i), string("id ", op.id, " is parked in the base case (status 0); status 1 does not bring it back, use set in_service 1"))
    end
  end
  return nothing
end

## --- reading ------------------------------------------------------------------

# the entries of a series zip as text, after the entry check
function _scf_series_read_entries(path::AbstractString)::Dict{String,String}
  isfile(path) || throw(ArgumentError(string("SCF series: file not found: ", path)))
  r = ZipArchives.ZipReader(read(String(path)))
  names = ZipArchives.zip_names(r)
  for nm in names
    nm in _SCF_SERIES_ENTRIES || _scf_series_error(basename(path), string("unexpected entry ", repr(nm), "; a series zip holds exactly ", join(_SCF_SERIES_ENTRIES, ", ")))
  end
  for nm in _SCF_SERIES_ENTRIES
    c = count(==(nm), names)
    c == 1 || _scf_series_error(basename(path), string(c == 0 ? "missing entry " : "duplicate entry ", repr(nm)))
  end
  return Dict{String,String}(nm => ZipArchives.zip_readentry(r, nm, String) for nm in _SCF_SERIES_ENTRIES)
end

# a parse or validation error of the SCF reader on one entry, re-thrown with
# the entry named
function _scf_series_entry_parse(f, entry::String, text::String)
  try
    return f(text)
  catch err
    err isa ArgumentError || rethrow()
    _scf_series_error(entry, err.msg)
  end
end

function _scf_series_step(raw, k::Int)::SCFSeriesStep
  at = string("series.json step ", k)
  raw isa AbstractDict || _scf_series_error(at, "a step must be an object")
  for key in keys(raw)
    String(key) in _SCF_SERIES_STEP_KEYS || _scf_series_error(at, string("unknown key ", repr(String(key)), "; known: ", join(_SCF_SERIES_STEP_KEYS, ", ")))
  end
  n = get(raw, "step", nothing)
  (n isa Integer && !(n isa Bool)) || _scf_series_error(at, "step must be an integer")
  ts = get(raw, "timestamp", "")
  ts isa AbstractString || _scf_series_error(at, "timestamp must be a string")
  label = get(raw, "label", "")
  label isa AbstractString || _scf_series_error(at, "label must be a string")
  # the axis label defaults to hh:mm of the time, else to the step number
  isempty(label) && (label = occursin(_SCF_SERIES_TIMESTAMP, ts) ? ts[12:16] : string(n))
  return SCFSeriesStep(Int(n), String(ts), String(label))
end

function _scf_series_manifest(text::String)
  at = "series.json"
  m = _scf_series_entry_parse(scf_json_parse, at, text)
  m isa AbstractDict || _scf_series_error(at, "the manifest must be a JSON object")
  for k in keys(m)
    String(k) in _SCF_SERIES_MANIFEST_KEYS || _scf_series_error(at, string("unknown key ", repr(String(k)), "; known: ", join(_SCF_SERIES_MANIFEST_KEYS, ", ")))
  end
  fmt = get(m, "format", nothing)
  fmt == _SCF_SERIES_FORMAT || _scf_series_error(at, string("format must be ", repr(_SCF_SERIES_FORMAT), ", got ", repr(fmt)))
  ver = get(m, "format_version", nothing)
  ver == _SCF_SERIES_FORMAT_VERSION || _scf_series_error(at, string("format_version ", repr(ver), " is not supported; this reader takes ", _SCF_SERIES_FORMAT_VERSION))
  name = get(m, "name", nothing)
  name isa AbstractString || _scf_series_error(at, "name must be a string")
  created_by = get(m, "created_by", "")
  created_by isa AbstractString || _scf_series_error(at, "created_by must be a string")
  minutes = get(m, "step_minutes", nothing)
  (minutes === nothing || (minutes isa Real && !(minutes isa Bool))) || _scf_series_error(at, "step_minutes must be a number or null")
  provenance = get(m, "provenance", Dict{String,Any}())
  provenance isa AbstractDict || _scf_series_error(at, "provenance must be an object")
  raw = get(m, "steps", nothing)
  raw isa AbstractVector || _scf_series_error(at, "steps must be a list")
  steps = SCFSeriesStep[_scf_series_step(s, k) for (k, s) in enumerate(raw)]
  return (name = String(name), created_by = String(created_by), step_minutes = minutes === nothing ? nothing : Float64(minutes),
          provenance = Dict{String,Any}(provenance), steps = steps)
end

# one delta row (six stripped cells) as a patch op; `at` names the CSV line
function _scf_series_row_op(cells::AbstractVector{<:AbstractString}, at::String)::PatchOp
  op = Symbol(cells[2])
  op in (:set, :scale, :status) || _scf_series_error(at, string("op must be set, scale or status, got ", repr(String(cells[2]))))
  target = Symbol(cells[3])
  target in _SCF_SERIES_TARGETS || _scf_series_error(at, string("target ", repr(String(cells[3])), " is not one of ", join(_SCF_SERIES_TARGETS, ", ")))
  id = tryparse(Int, cells[4])
  id === nothing && _scf_series_error(at, string("id ", repr(String(cells[4])), " is not an integer"))
  field = isempty(cells[5]) ? nothing : Symbol(cells[5])
  value = nothing
  if !isempty(cells[6])
    value = tryparse(Float64, cells[6])
    (value !== nothing && isfinite(value)) || _scf_series_error(at, string("value ", repr(String(cells[6])), " is not a finite number"))
  end
  if op === :status
    field === nothing || _scf_series_error(at, "status takes no field")
    value === nothing && _scf_series_error(at, "status needs the value 0 or 1")
    return PatchOp(op = :status, target = target, id = id, value = value)
  elseif op === :scale
    value === nothing && _scf_series_error(at, "scale needs a factor in the value column")
    return PatchOp(op = :scale, target = target, id = id, field = field, factor = value)
  end
  field === nothing && _scf_series_error(at, "set needs a field")
  if value === nothing
    rule = get(_QSTS_SET_RULES, field, nothing)
    (rule !== nothing && rule[2] === :optional) || _scf_series_error(at, string("set ", field, " needs a value; only ",
      join(sort!([String(f) for (f, r) in _QSTS_SET_RULES if r[2] === :optional]), ", "), " take an empty one"))
  end
  return PatchOp(op = :set, target = target, id = id, field = field, value = value)
end

# delta.csv: UTF-8 (a byte order mark is ignored), comma, decimal point, LF
# (a CR before it is dropped), `#` comments only above the header, rows
# sorted by step. Every error names the CSV line.
function _scf_series_parse_delta(text::AbstractString, n::Int, index::ScenarioIndex, parked::Set{Int})::Vector{Vector{PatchOp}}
  ops = [PatchOp[] for _ in 1:n]
  lines = [Int[] for _ in 1:n]
  header_seen = false
  last_step = 0
  for (lineno, raw) in enumerate(eachsplit(text, '\n'))
    line = rstrip(lineno == 1 ? lstrip(raw, '﻿') : raw, '\r')
    isempty(strip(line)) && continue
    at = string("delta.csv line ", lineno)
    if !header_seen
      startswith(line, "#") && continue
      strip(line) == _SCF_SERIES_HEADER || _scf_series_error(at, string("the header must be ", _SCF_SERIES_HEADER))
      header_seen = true
      continue
    end
    startswith(line, "#") && _scf_series_error(at, "a comment is allowed only above the header")
    cells = strip.(split(line, ','))
    length(cells) == 6 || _scf_series_error(at, string("expected 6 fields (", _SCF_SERIES_HEADER, "), got ", length(cells)))
    k = tryparse(Int, cells[1])
    (k !== nothing && 1 <= k <= n) || _scf_series_error(at, string("step ", repr(String(cells[1])), " is not a step of series.json (1 to ", n, ")"))
    k >= last_step || _scf_series_error(at, string("step ", k, " after step ", last_step, "; the rows are sorted by step"))
    last_step = k
    push!(ops[k], _scf_series_row_op(cells, at))
    push!(lines[k], lineno)
  end
  header_seen || _scf_series_error("delta.csv", string("no header line ", _SCF_SERIES_HEADER))
  for k in 1:n
    _scf_series_check_step(ops[k], k, index, parked, i -> string("delta.csv line ", lines[k][i]))
  end
  return ops
end

"""
    is_scf_series(path) -> Bool

`true` when `path` is a zip file whose top level holds the entry
`series.json`, the manifest of an SCF time series. The content decides,
not the file name. Only the end of the file is read unless that end names
the manifest, so the probe stays cheap on a large CGMES delivery. A file
that ends like a zip naming `series.json` but cannot be read as a zip
raises the `ArgumentError` of the zip reader.
"""
function is_scf_series(path::AbstractString)::Bool
  isfile(path) || return false
  n = filesize(path)
  n >= 22 || return false
  tail = open(String(path)) do io
    # a zip opens with a local file header
    read(io, 4) == UInt8[0x50, 0x4b, 0x03, 0x04] || return UInt8[]
    seek(io, max(0, n - _SCF_SERIES_PROBE_BYTES))
    return read(io)
  end
  findfirst(Vector{UInt8}(codeunits("series.json")), tail) === nothing && return false
  return "series.json" in ZipArchives.zip_names(ZipArchives.ZipReader(read(String(path))))
end

"""
    read_scf_series(path) -> SCFSeries

Read an SCF time series zip and check everything before the first step:

- the zip holds exactly `series.json`, `base.scf.json` and `delta.csv`;
- the manifest: format `sparlectra-scf-series`, `format_version` exactly
  1.0, steps numbered 1 to n, timestamps on every step or none and strictly
  increasing;
- the base case: the checks of [`read_scf_json`](@ref), plus the build-order
  fields every component needs;
- `delta.csv`: header `step,op,target,id,field,value`, six fields, finite
  numbers, a known step, rows sorted by step, every op through the checks
  of [`validate_scenarios`](@ref) on the series index, one write per id and
  field in a step (`scale` without field writes p and q), no other op on a
  unit the step parks with `set in_service 0`, and a set, scale or
  status 1 on a unit the base case parks only together with
  `set in_service 1`.

Every error is an `ArgumentError` that names the entry and, for
`delta.csv`, the line. Builds no net, except one build of the base case
when it parks units (their values go to `index.parked_values`).
"""
function read_scf_series(path::AbstractString)::SCFSeries
  entries = _scf_series_read_entries(path)
  m = _scf_series_manifest(entries["series.json"])
  _scf_series_check_steps(m.steps, m.step_minutes)
  base = _scf_series_entry_parse(_scfcase_from_text, "base.scf.json", entries["base.scf.json"])
  index, parked = _scf_series_index(base)
  ops = _scf_series_parse_delta(entries["delta.csv"], length(m.steps), index, Set(parked))
  return SCFSeries(m.name, m.created_by, base, index, parked, m.steps, ops, m.step_minutes, m.provenance)
end

## --- writing ------------------------------------------------------------------

# one value cell: a whole number without a trailing .0 (status 0 and 1, a
# factor of 1), everything else in the shortest round-trip form
_scf_series_number(x::Float64)::String = (isinteger(x) && abs(x) < 1.0e15) ? string(Int(x)) : string(x)

# the three entries of a series as text, in zip order
function _scf_series_texts(series::SCFSeries)::Vector{Pair{String,String}}
  minutes = series.step_minutes
  manifest = Dict{String,Any}(
    "format" => _SCF_SERIES_FORMAT,
    "format_version" => _SCF_SERIES_FORMAT_VERSION,
    "name" => series.name,
    "created_by" => series.created_by,
    "step_minutes" => minutes === nothing ? nothing : (isinteger(minutes) ? Int(minutes) : minutes),
    "steps" => Any[Dict{String,Any}("step" => s.step, "timestamp" => s.timestamp, "label" => s.label) for s in series.steps],
    "provenance" => series.provenance,
  )
  io = IOBuffer()
  println(io, _SCF_SERIES_HEADER)
  for (k, step_ops) in enumerate(series.ops), op in step_ops
    v = op.op === :scale ? op.factor : op.value
    println(io, k, ',', op.op, ',', op.target, ',', op.id, ',', op.field === nothing ? "" : op.field, ',', v === nothing ? "" : _scf_series_number(v))
  end
  return ["series.json" => scf_json_string(manifest), "base.scf.json" => scf_json_string(scf_root_dict(series.base)), "delta.csv" => String(take!(io))]
end

# write the entries as one zip, every entry deflated
function _scf_series_write_zip(path::AbstractString, entries::Vector{Pair{String,String}})::String
  p = String(path)
  mkpath(dirname(abspath(p)))
  ZipArchives.ZipWriter(p) do w
    for (name, text) in entries
      ZipArchives.zip_newfile(w, name; compress = true)
      write(w, text)
    end
  end
  return p
end

"""
    write_scf_series(path, series) -> String

Write `series` as one zip with the entries `series.json`, `base.scf.json`
and `delta.csv`, compressed, and return the path. The entry texts are
deterministic: writing what [`read_scf_series`](@ref) returns gives the
same entries again. The zip bytes depend on the zlib build, so compare
series by their entries, not by the zip file.
"""
write_scf_series(path::AbstractString, series::SCFSeries)::String = _scf_series_write_zip(path, _scf_series_texts(series))

"""
    build_net(series::SCFSeries; config = active_sparlectra_config()) -> Net

The base net of an SCF time series: the net of the base case
(`build_net(series.base; config)`), with every appliance the base case
marks with status 0 parked through `set in_service 0`. The parking is the
base state of the series and stays. A step applies with
`apply!(net, series.ops[k], series.index)` and goes back with
[`restore!`](@ref).
"""
function build_net(series::SCFSeries; config::SparlectraConfig = active_sparlectra_config())::Net
  net = build_net(series.base; config = config)
  isempty(series.parked) && return net
  park = PatchOp[PatchOp(op = :set, target = series.index.kind_by_id[id], id = id, field = :in_service, value = 0.0) for id in series.parked]
  apply!(net, park, series.index)
  return net
end
