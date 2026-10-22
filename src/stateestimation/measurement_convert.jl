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
# file: src/stateestimation/measurement_convert.jl
# purpose: measurement CSV v1 to measurement table. convertMeasurementsV1
#          reads a v1 file against its case and writes the table with the
#          same measurements (units of the table, keys of a
#          MeasurementKeyIndex, comments carried over with the row ids
#          rewritten); _measurement_table_from_measurements is the row
#          builder the Web UI generator shares.

# The table row of every measurement of `ms` against `index`, as NamedTuples
# with string numbers (writeMeasurementTable writes them as they are), plus
# the measurements that have no place in the table and the id each converted
# measurement gets on import. `lines` gives the source line of each
# measurement for the problem list (its position when omitted).
#
# Placement: Vm, Va, Pinj and Qinj at a bus go to a node key of that bus
# (U in kV as the exact text of _measurement_table_kv_text, so the import
# gives back the pu value bit for bit); flows, currents and current angles
# at a branch end to its terminal key; the bay current and the shunt Q of a
# bus to the terminal of its one shunt. Zero-injection rows (id prefix ZI)
# carry the source token zero_injection and keep their sigma. Units are
# unchanged otherwise: P and Q in MW and MVAr, I in A, angles in degrees.
# `source` is the source cell of every row (`converted;v1` for the
# converter, the SE run artifact names where its rows came from).
function _measurement_table_from_measurements(net::Net, ms::AbstractVector, index::MeasurementKeyIndex; lines::AbstractVector{<:Integer} = collect(eachindex(ms)), source::AbstractString = "converted;v1")
  Row = NamedTuple{(:timestamp, :measurement_mrid, :terminal_mrid, :node_mrid, :type, :value, :sigma, :source, :valid),Tuple{String,String,String,String,String,String,String,String,Bool}}
  Problem = NamedTuple{(:line, :id, :problem, :detail),Tuple{Int,String,Symbol,String}}
  rows = Row[]
  problems = Problem[]
  ids = Dict{String,String}()
  # reverse lookups: bus -> node key, branch end -> terminal key, bus -> the
  # terminal of its single shunt. Several keys of one bus (TopologicalNode,
  # ConnectivityNode, busbar of a CGMES case) are equivalent; the smallest
  # is taken so the output does not depend on dictionary order.
  node_key = Dict{Int,String}()
  for (k, b) in index.nodes
    b > 0 || continue
    (haskey(node_key, b) && node_key[b] <= k) || (node_key[b] = k)
  end
  end_key = Dict{Tuple{Int,Symbol},String}()
  for (k, t) in index.terminals
    t.kind === :branch || continue
    e = (t.branchIdx, t.direction)
    (haskey(end_key, e) && end_key[e] <= k) || (end_key[e] = k)
  end
  # link -> terminal key of a link end, the from end preferred (a link flow
  # measurement is positive from link.fromBus, a row at the to end carries
  # the opposite sign)
  link_key = Dict{Int,Tuple{String,Symbol}}()
  for (k, t) in index.terminals
    t.kind === :link || continue
    cur = get(link_key, t.branchIdx, nothing)
    better = cur === nothing || (t.direction === :from && cur[2] === :to) || (t.direction === cur[2] && k < cur[1])
    better && (link_key[t.branchIdx] = (k, t.direction))
  end
  shunt_members(b) = get(index.groups, (:shunt, b), String[])
  problem!(i, m, p, d) = push!(problems, (line = Int(lines[i]), id = m.id, problem = p, detail = d))
  seen = Set{Tuple{String,String,String}}()
  for (i, m) in enumerate(ms)
    zi = m.zero_injection
    row_source = zi ? string(source, ";zero_injection") : String(source)
    kind = ""
    key = ""
    typ = ""
    value = string(m.value)
    sigma = string(m.sigma)
    if m.linkIdx !== nothing
      typ = m.typ == PflowMeas ? "P" : m.typ == QflowMeas ? "Q" : ""
      lk = get(link_key, Int(m.linkIdx), nothing)
      if isempty(typ) || lk === nothing
        problem!(i, m, isempty(typ) ? :unsupported : :unknown_key, isempty(typ) ? "$(m.typ) on bus link $(m.linkIdx) is no table type" : "bus link $(m.linkIdx) has no terminal key in the index (open link or no name)")
        continue
      end
      key = lk[1]
      lk[2] === :to && (value = string(-m.value))
      kind = "terminal"
    elseif m.branchIdx !== nothing
      typ = m.typ == PflowMeas ? "P" : m.typ == QflowMeas ? "Q" : m.typ == ImagMeas ? "I" : m.typ == IaMeas ? "Ia" : ""
      if isempty(typ)
        problem!(i, m, :unsupported, "$(m.typ) at a branch end is no table type")
        continue
      end
      key = get(end_key, (Int(m.branchIdx), m.direction), "")
      if isempty(key)
        problem!(i, m, :unknown_key, "branch $(m.branchIdx) end $(m.direction) has no key in the index (open branch, no name, or a name shared by several elements)")
        continue
      end
      kind = "terminal"
    elseif m.busIdx !== nothing
      b = Int(m.busIdx)
      if m.typ in (ShuntQMeas, ImagMeas)
        members = shunt_members(b)
        if length(members) != 1 || get(index.terminals, members[1], nothing) === nothing
          problem!(i, m, :unsupported, "$(m.typ) at bus $(b): the table keys the shunt bay by the terminal of its shunt, and the bus has $(length(members)) shunt(s) with a key")
          continue
        end
        key = members[1]
        typ = m.typ == ShuntQMeas ? "Q" : "I"
        kind = "terminal"
      elseif m.typ in (VmMeas, VaMeas, PinjMeas, QinjMeas)
        key = get(node_key, b, _measurement_table_dead_node_key(net, index, b))
        if isempty(key)
          problem!(i, m, :unknown_key, "bus $(b) has no node key in the index (not energized, or a name shared by several buses)")
          continue
        end
        kind = "node"
        if m.typ == VmMeas
          typ = "U"
          vn = Float64(index.bus_kV[b])
          value = _measurement_table_kv_text(m.value, vn)
          sigma = _measurement_table_kv_text(m.sigma, vn)
        else
          typ = m.typ == VaMeas ? "Va" : m.typ == PinjMeas ? "P" : "Q"
        end
      else
        problem!(i, m, :unsupported, "$(m.typ) at a bus is no table type (a current angle at a shunt bay)")
        continue
      end
    else
      problem!(i, m, :unknown_key, "$(m.typ) without a location")
      continue
    end
    # U and I are magnitudes in the table; a v1 row with a negative value
    # (noise drawn on a near-zero current) is refused by the reader there
    if (typ == "U" || typ == "I") && m.value < 0.0
      problem!(i, m, :bad_value, "$(typ) $(m.value) is negative; the table reads $(typ) as a magnitude")
      continue
    end
    # the same key and type twice is a duplicate in the table; the v1 file
    # may hold such redundant readings, the table cannot
    if (kind, key, typ) in seen
      problem!(i, m, :duplicate, "a second $(typ) row at $(key); the table holds one row per key and type")
      continue
    end
    push!(seen, (kind, key, typ))
    push!(rows, Row(("", "", kind == "terminal" ? key : "", kind == "node" ? key : "", typ, value, sigma, row_source, m.active)))
    # the id importMeasurementTable gives the row
    new_id = if kind == "node"
      (typ == "P" || typ == "Q") ? (zi ? string(typ == "P" ? "ZI_PINJ_bus_" : "ZI_QINJ_bus_", m.busIdx) : string(typ, "bus@", key)) : string(typ, "@", key)
    elseif m.busIdx !== nothing
      string(m.typ == ShuntQMeas ? ShuntQMeas : ImagMeas, "@bus", m.busIdx)
    else
      string(typ, "@", key)
    end
    isempty(m.id) || (ids[m.id] = new_id)
  end
  return (rows = rows, problems = problems, ids = ids)
end

# The key of a bus the index knows but does not energize (an isolated bus:
# the index maps its key to 0): its mRID for a CGMES index, else its name.
# A row written with it is reported as not_modelled by the import, the way
# the SE service leaves a v1 row at an isolated bus out of the estimate;
# the row stays in the file. Empty when the index has no such key.
function _measurement_table_dead_node_key(net::Net, index::MeasurementKeyIndex, b::Int)::String
  name = get(_bus_name_by_idx(net), b, "")
  for cand in (index.canonical ? _bus_mrid(net, b) : "", name)
    isempty(cand) && continue
    k = _index_key(index, cand)
    get(index.nodes, k, -1) == 0 && return k
  end
  return ""
end

# comment lines with the v1 ids of the generator blocks rewritten to the
# ids of the table: truth_value,<id>,<value> and the id lists of the
# critical thinning. The truth line of a row that was not written
# (`dropped`) goes: its row is not in the table, the problem list names it.
function _measurement_table_comment_ids(comments::Vector{String}, ids::Dict{String,String}; dropped::AbstractSet{String} = Set{String}())::Vector{String}
  out = String[]
  for c in comments
    if startswith(c, "truth_value,")
      parts = split(chopprefix(c, "truth_value,"), ",")
      old = join(parts[1:(end-1)], ",")
      (old in dropped && !haskey(ids, old)) && continue
      c = haskey(ids, old) ? string("truth_value,", ids[old], ",", parts[end]) : c
    elseif startswith(c, "critical_rows:") || startswith(c, "critical_removed:")
      label, rest = split(c, ":"; limit = 2)
      c = string(label, ": ", join((get(ids, String(t), String(t)) for t in split(rest)), " "))
    end
    push!(out, String(c))
  end
  return out
end

# the comment lines of a measurement CSV v1 in file order, without the
# version line and without the leading "# "
function _measurement_csv_v1_comments(file::AbstractString)::Vector{String}
  out = String[]
  for (ln, raw) in enumerate(eachline(file))
    ln == 1 && continue
    s = strip(raw)
    startswith(s, "#") && push!(out, String(strip(chopprefix(s, "#"))))
  end
  return out
end

"""
    convertMeasurementsV1(v1_file, net; out, index = nameMeasurementKeys(net)) -> NamedTuple

Write the measurement table of a measurement CSV v1 file. `net` is the
case the v1 file belongs to (the v1 file addresses buses and branches by
name and number, the table by the keys of `index`); for a CGMES case pass
`index = cgmesMeasurementKeys(result)` so the table carries mRIDs, as the
SE service reads it.

Every v1 row becomes one table row with the same measurement:

- `VmMeas` at a bus: `U` at a node key of the bus in kV, value and sigma
  written so that the import divides them back to the v1 pu values bit for
  bit (up to 18 digits where a Float64 kV value cannot).
- `VaMeas`, `PinjMeas`, `QinjMeas` at a bus: `Va`, `P`, `Q` at the node
  (degrees, bus balance in MW and MVAr, injection positive). A row at an
  isolated bus is written with the key of that bus; the import reports it
  as `not_modelled`, as the SE service leaves such a v1 row out.
- `PflowMeas`, `QflowMeas`, `ImagMeas`, `IaMeas` at a branch end: `P`, `Q`,
  `I`, `Ia` at the terminal key of that end.
- `ShuntQMeas` and the bay `ImagMeas` at a bus: `Q` and `I` at the
  terminal of the one shunt of that bus.
- zero-injection rows (id prefix `ZI`): `source = converted;v1;zero_injection`,
  their sigma kept; every other row `source = converted;v1`. `valid` is the
  `active` flag of the row.

No unit is converted except `U`: P and Q stay in MW and MVAr, I in A,
angles in degrees. Rows the table cannot carry (a link flow, a shunt bay
at a bus with several shunts, a current angle at a shunt bay, a location
without a key, a negative `U` or `I`, a second row with the same key and
type) are not written;
they are listed in the result and as comment lines of the table, never
dropped silently. The comments of the v1 file are carried over in front of
the header, with the row ids in `truth_value` and critical-row lines
rewritten to the ids the import gives, so the SE service still writes
`se_deltas.csv` and reads the tap table.

Returns `(file, rows, problems, ids, report)`: the written path, the number
of table rows, the rows not converted as `(line, id, problem, detail)`, the
v1 id to table id map and a plain-text report. A file the v1 parse refuses
(version line, header, a row) throws an error naming file and line. The v1
reader `readMeasurementsCSV!` was removed in 0.50.0; this converter still
reads every v1 file.
"""
function convertMeasurementsV1(v1_file::AbstractString, net::Net; out::AbstractString, index::MeasurementKeyIndex = nameMeasurementKeys(net))
  v1 = _read_measurements_csv_v1(net, v1_file)
  conv = _measurement_table_from_measurements(net, v1.measurements, index; lines = v1.lines)
  comments = _measurement_table_comment_ids(_measurement_csv_v1_comments(v1_file), conv.ids; dropped = Set(p.id for p in conv.problems))
  push!(comments, string("converted from ", basename(v1_file), " (measurement CSV v1): ", length(conv.rows), " row(s), ", length(conv.problems), " not converted"))
  for p in conv.problems
    push!(comments, string("not converted: line ", p.line, " ", p.id, ": ", p.problem, ": ", p.detail))
  end
  writeMeasurementTable(out, conv.rows; comments = comments)
  io = IOBuffer()
  println(io, basename(v1_file), " -> ", basename(out), ": ", length(v1.measurements), " v1 row(s), ", length(conv.rows), " table row(s), ", length(conv.problems), " not converted")
  for p in conv.problems
    println(io, "  line ", p.line, " ", p.id, ": ", p.problem, ": ", p.detail)
  end
  return (file = String(out), rows = length(conv.rows), problems = conv.problems, ids = conv.ids, report = String(take!(io)))
end

# The key index a measurement table of an imported case is read with: mRIDs
# when the import kept a CGMES model (provenance "cgmes_result"), element
# names otherwise. The SE service decides the same way.
function _measurement_key_index(imported)::MeasurementKeyIndex
  cg = get(imported.provenance, "cgmes_result", nothing)
  return cg === nothing ? nameMeasurementKeys(imported.net) : cgmesMeasurementKeys(cg)
end

# A measurement file into `net.measurements`: a measurement table through
# readMeasurementTable! with `index` (no zero-injection rows added). Used
# where a case brings its measurement file along (the sidecar
# <case>.measurements.csv of the CLI and the Web UI export, save as and add
# noise). A measurement CSV v1 is refused with MEASUREMENT_CSV_V1_REMOVED,
# which names the converter; any other file with the reason of the table
# reader.
function _read_measurement_file!(net::Net, file::AbstractString; index::Union{Nothing,MeasurementKeyIndex} = nothing)
  _is_measurement_csv_v1(file) && throw(ArgumentError("$(basename(file)): " * MEASUREMENT_CSV_V1_REMOVED))
  return readMeasurementTable!(net; file = file, index = something(index, nameMeasurementKeys(net)), zero_injection = false)
end
