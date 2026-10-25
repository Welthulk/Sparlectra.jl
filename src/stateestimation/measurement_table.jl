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
# file: src/stateestimation/measurement_table.jl
# purpose: the measurement table (one CSV row per measured value, keyed by
#          measurement, terminal or node mRID or by element name, several
#          snapshots per file): reader, writer, and the import into SE measurement sets
#          against a format-independent key index. The CGMES importer builds
#          the mRID index (cgmesMeasurementKeys); nameMeasurementKeys here
#          builds the name index of any net, and
#          generateMeasurementTableFromPF writes a synthetic table from a
#          solved power flow, generatePseudoMeasurementTable the pseudo rows
#          for unobserved injections. Nothing here reads CGMES.

# the fixed header; the order is part of the format
const MEASUREMENT_TABLE_HEADER = ["timestamp", "measurement_mrid", "terminal_mrid", "node_mrid", "type", "value", "sigma", "source", "valid"]

# the version line every table carries before its header. The measurement
# CSV v1 has its own line (# sparlectra-measurements v1); the two never
# match each other.
const MEASUREMENT_TABLE_VERSION_LINE = "# sparlectra-measurement-table v1"
# any version of the line: group 1 is everything after the name (the
# version token, and text after it that the reader refuses by name)
const _MEASUREMENT_TABLE_VERSION_RE = r"^#\s*sparlectra-measurement-table(?:\s+(.*?))?\s*$"

# the types of a terminal row and of a node row. At a terminal: P, Q (MW,
# MVAr into the equipment), I (A) and the current angle Ia (degrees). At a
# node: U (kV), the bus balance P and Q (MW, MVAr, injection positive, the
# generator convention of PinjMeas) and the voltage angle Va (degrees).
const _MEASUREMENT_TABLE_TERMINAL_TYPES = ["P", "Q", "I", "Ia"]
const _MEASUREMENT_TABLE_NODE_TYPES = ["U", "P", "Q", "Va"]

# ISO 8601 in UTC: date, time to minutes or seconds (fraction allowed), "Z"
const _MEASUREMENT_TABLE_TIMESTAMP = r"^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2})(?::(\d{2})(?:\.(\d+))?)?Z$"

# The canonical form of a timestamp, `nothing` when it is not ISO 8601 UTC
# or names an impossible date or time (2025-02-30, 25:61). One instant has
# one form: minutes only when seconds and fraction are zero
# (12:00Z = 12:00:00Z = 12:00:00.000Z), otherwise seconds with the fraction
# without trailing zeros. The sets of a table are formed by this form.
function _measurement_table_timestamp(ts::AbstractString)::Union{Nothing,String}
  m = match(_MEASUREMENT_TABLE_TIMESTAMP, ts)
  m === nothing && return nothing
  y, mo, d, h, mi = (parse(Int, c) for c in m.captures[1:5])
  sec = m.captures[6] === nothing ? 0 : parse(Int, m.captures[6])
  frac = m.captures[7] === nothing ? "" : rstrip(m.captures[7], '0')
  (1 <= mo <= 12 && 1 <= d <= Dates.daysinmonth(y, mo) && h <= 23 && mi <= 59 && sec <= 59) || return nothing
  base = string(m.captures[1], "-", m.captures[2], "-", m.captures[3], "T", m.captures[4], ":", m.captures[5])
  (sec == 0 && isempty(frac)) && return base * "Z"
  return string(base, ":", lpad(sec, 2, '0'), isempty(frac) ? "" : string(".", frac), "Z")
end

"""
    MeasurementRow

One row of a measurement table as read from the file. `value` is `NaN` and
`sigma` is `nothing` when the cell was empty or not a number; `problem` is
`:none` for a well-formed row, otherwise the reason the row cannot become a
measurement (`:field_count`, `:key_count`, `:unknown_type`, `:bad_value`,
`:bad_sigma`, `:bad_valid`, `:bad_timestamp`), with `detail` as text.
`line` is the 1-based line number in the file. `value_text` and
`sigma_text` keep the cells as written (empty for a row built in code); a
`U` cell with more than 17 significant digits is divided by the nominal
voltage from this text at 256 bits, so a kV value written by
`convertMeasurementsV1` gives back the pu value of the measurement CSV v1
bit for bit. The twelve-argument constructor leaves
both empty.
"""
struct MeasurementRow
  line::Int
  timestamp::String
  measurement_mrid::String
  terminal_mrid::String
  node_mrid::String
  type::String
  value::Float64
  sigma::Union{Nothing,Float64}
  source::String
  valid::Bool
  problem::Symbol
  detail::String
  value_text::String
  sigma_text::String
end

MeasurementRow(line, timestamp, measurement_mrid, terminal_mrid, node_mrid, type, value, sigma, source, valid, problem, detail) = MeasurementRow(line, timestamp, measurement_mrid, terminal_mrid, node_mrid, type, value, sigma, source, valid, problem, detail, "", "")

# minimal RFC 4180 field split: commas separate, double quotes enclose a
# field that contains commas or quotes, "" inside quotes is one quote
function _measurement_table_fields(line::AbstractString)::Vector{String}
  fields = String[]
  buf = IOBuffer()
  inquote = false
  chars = collect(line)
  i = 1
  while i <= length(chars)
    c = chars[i]
    if inquote
      if c == '"'
        if i < length(chars) && chars[i+1] == '"'
          write(buf, '"')
          i += 1
        else
          inquote = false
        end
      else
        write(buf, c)
      end
    elseif c == '"'
      inquote = true
    elseif c == ','
      push!(fields, String(take!(buf)))
    else
      write(buf, c)
    end
    i += 1
  end
  push!(fields, String(take!(buf)))
  return fields
end

# a cell as written: without surrounding blanks (the reader strips every
# cell, so a write after a read gives the same bytes), quoted only when it
# carries a comma or a quote. The reader reads one row per line, so a line
# break inside a cell is refused here instead of writing a file the reader
# cannot read back.
function _measurement_table_cell(cell::AbstractString)::String
  any(c -> c == '\n' || c == '\r', cell) && throw(ArgumentError("writeMeasurementTable: a cell must not contain a line break, got $(repr(cell))"))
  s = strip(cell)
  any(c -> c == ',' || c == '"', s) || return String(s)
  return string('"', replace(s, "\"" => "\"\""), '"')
end

# the row-level checks of the format; resolution against a model comes later
function _measurement_table_row(lineno::Int, f::Vector{String})::MeasurementRow
  if length(f) != length(MEASUREMENT_TABLE_HEADER)
    return MeasurementRow(lineno, "", "", "", "", "", NaN, nothing, "", false, :field_count, "$(length(f)) fields, expected $(length(MEASUREMENT_TABLE_HEADER))")
  end
  raw_ts, mm, tm, nm, typ, val, sig, src, vld = strip.(f)
  value = something(tryparse(Float64, val), NaN)
  sigma = isempty(sig) ? nothing : tryparse(Float64, sig)
  valid = vld == "true" ? true : (vld == "false" ? false : nothing)
  # the row carries the canonical timestamp, so one instant is one set
  ts = isempty(raw_ts) ? "" : something(_measurement_table_timestamp(raw_ts), raw_ts)
  row(problem, detail) = MeasurementRow(lineno, ts, mm, tm, nm, typ, value, sigma, src, valid === nothing ? false : valid, problem, detail, val, sig)
  if !isempty(raw_ts) && _measurement_table_timestamp(raw_ts) === nothing
    occursin(_MEASUREMENT_TABLE_TIMESTAMP, raw_ts) && return row(:bad_timestamp, "timestamp '$(raw_ts)' is not a valid date and time")
    return row(:bad_timestamp, "timestamp '$(raw_ts)' is not ISO 8601 UTC (YYYY-MM-DDThh:mm[:ss]Z)")
  end
  nkeys = count(!isempty, (mm, tm, nm))
  nkeys == 1 || return row(:key_count, "$(nkeys) of measurement_mrid, terminal_mrid, node_mrid set, expected exactly one")
  if !isempty(mm)
    isempty(typ) || return row(:unknown_type, "type must be empty for a measurement_mrid row (it comes from the Analog), got '$(typ)'")
  elseif !isempty(tm)
    typ in _MEASUREMENT_TABLE_TERMINAL_TYPES || return row(:unknown_type, "type '$(typ)' at a terminal, expected P, Q, I or Ia")
  else
    typ in _MEASUREMENT_TABLE_NODE_TYPES || return row(:unknown_type, "type '$(typ)' at a node, expected U, P, Q or Va")
  end
  isfinite(value) || return row(:bad_value, "value '$(val)' is not a finite number")
  (typ in ("U", "I") && value < 0.0) && return row(:bad_value, "value '$(val)' is negative, $(typ) is a magnitude")
  if !isempty(sig)
    (sigma !== nothing && isfinite(sigma) && sigma > 0.0) || return row(:bad_sigma, "sigma '$(sig)' is not a positive number")
  end
  valid === nothing && return row(:bad_valid, "valid '$(vld)' is neither true nor false")
  return row(:none, "")
end

"""
    readMeasurementTable(file) -> Vector{MeasurementRow}

Read a measurement table: CSV, UTF-8, comma separated, decimal point, LF line
ends. Before the header stands the version line
`# sparlectra-measurement-table v1`, then the header
`timestamp,measurement_mrid,terminal_mrid,node_mrid,type,value,sigma,source,valid`.
Other lines starting with `#` before the header are comments, blank lines
are skipped. A file without the header, without the version line or with
another version is rejected with an `ArgumentError` that names what was
expected. Rows that break a rule of the format are returned with their
`problem` set, never dropped, so the import report can count them.
"""
function readMeasurementTable(file::AbstractString)::Vector{MeasurementRow}
  rows = MeasurementRow[]
  header_seen = false
  version = nothing
  version_line = 0
  lineno = 0
  for raw in eachline(file)
    lineno += 1
    line = lineno == 1 ? lstrip(raw, '﻿') : raw
    line = rstrip(line, '\r')
    isempty(strip(line)) && continue
    if !header_seen
      if startswith(line, "#")
        m = match(_MEASUREMENT_TABLE_VERSION_RE, strip(line))
        if m !== nothing && version === nothing
          version = something(m.captures[1], "")
          version_line = lineno
          occursin(r"\s", version) && throw(ArgumentError("$(file):$(lineno): the version line carries text after the version ('$(strip(line))'), expected '$(MEASUREMENT_TABLE_VERSION_LINE)'"))
        end
        continue
      end
      strip.(_measurement_table_fields(line)) == MEASUREMENT_TABLE_HEADER || throw(ArgumentError("$(file):$(lineno): not a measurement table, the header must be $(join(MEASUREMENT_TABLE_HEADER, ","))"))
      version === nothing && throw(ArgumentError("$(file):$(lineno): measurement table without a version line, expected '$(MEASUREMENT_TABLE_VERSION_LINE)' before the header"))
      version == "v1" || throw(ArgumentError("$(file):$(version_line): measurement table version '$(version)', this reader reads v1 ('$(MEASUREMENT_TABLE_VERSION_LINE)')"))
      header_seen = true
      continue
    end
    push!(rows, _measurement_table_row(lineno, _measurement_table_fields(line)))
  end
  header_seen || throw(ArgumentError("$(file): not a measurement table, no header line"))
  return rows
end

"""
    isMeasurementTable(file) -> Bool

`true` when `file` carries a measurement-table version line
(`# sparlectra-measurement-table <version>`) among the `#` lines before its
first other line, and that line is the measurement-table header. Decides by
content, not by file name. A table of another version counts, so the reader
can name the version it refuses; the measurement CSV v1 does not.
"""
function isMeasurementTable(file::AbstractString)::Bool
  isfile(file) || return false
  first_line = true
  versioned = false
  for raw in eachline(file)
    line = rstrip(first_line ? lstrip(raw, '﻿') : raw, '\r')
    first_line = false
    isempty(strip(line)) && continue
    if startswith(line, "#")
      versioned |= match(_MEASUREMENT_TABLE_VERSION_RE, strip(line)) !== nothing
      continue
    end
    return versioned && strip.(_measurement_table_fields(line)) == MEASUREMENT_TABLE_HEADER
  end
  return false
end

# A file whose first line after the `#` lines is the measurement-table
# header, but without any version line before it: a table that lost its
# version line. isMeasurementTable says false for it; the SE service and the
# Web UI use this to give the reader's message instead of treating the file
# as a measurement CSV v1 or as unknown.
function _measurement_table_without_version(file::AbstractString)::Bool
  isfile(file) || return false
  first_line = true
  for raw in eachline(file)
    line = rstrip(first_line ? lstrip(raw, '﻿') : raw, '\r')
    first_line = false
    isempty(strip(line)) && continue
    if startswith(line, "#")
      match(_MEASUREMENT_TABLE_VERSION_RE, strip(line)) === nothing || return false
      continue
    end
    return strip.(_measurement_table_fields(line)) == MEASUREMENT_TABLE_HEADER
  end
  return false
end

# A number cell as written: the text the row was read with when it has one
# (a MeasurementRow from the reader keeps it, and a U value in kV can carry
# more digits than a Float64 holds), a string as given (the converter writes
# the exact kV text), otherwise the shortest round-trip form of the number.
# `nothing` (no sigma) is an empty cell.
function _measurement_table_number(r, field::Symbol)::String
  text = field === :value ? :value_text : :sigma_text
  if r isa MeasurementRow
    t = getfield(r, text)
    isempty(t) || return t
  end
  x = getproperty(r, field)
  x === nothing && return ""
  if x isa AbstractString
    t = String(strip(x))
    tryparse(Float64, t) === nothing && throw(ArgumentError("writeMeasurementTable: $(field) $(repr(x)) is not a number"))
    return t
  end
  return string(Float64(x))
end

# U and its sigma in pu of the bus. A cell of up to 17 significant digits
# (every number a Float64 prints, and what people type) is divided as the
# Float64 it parses to. A Float64 kV value cannot give back every pu value
# (for about a third of them no Float64 divides back to them), so the
# converter writes a longer text where it has to (_measurement_table_kv_text);
# such a text is divided at 256 bits and rounded once, which gives back the
# pu value of the measurement CSV v1 bit for bit. Without a text (a row
# built in code) the Float64 quotient.
function _measurement_table_pu(text::AbstractString, value::Float64, vn::Float64)::Float64
  _measurement_table_digits(text) <= 17 && return value / vn
  return setprecision(BigFloat, 256) do
    Float64(BigFloat(text) / BigFloat(vn))
  end
end

# significant digits of a number cell (the mantissa without sign, point and
# leading zeros); 0 for an empty cell
function _measurement_table_digits(text::AbstractString)::Int
  m = split(lowercase(text), 'e')[1]
  d = lstrip(filter(isdigit, m), '0')
  return length(d)
end

# The decimal text of a kV value that _measurement_table_pu turns back into
# exactly `pu`: the Float64 product when that divides back, else the exact
# product to 18 significant digits (0.5e-18 relative from it, far inside the
# rounding interval of `pu`), checked like the first.
function _measurement_table_kv_text(pu::Float64, vn::Float64)::String
  u = pu * vn
  t = string(u)
  _measurement_table_pu(t, u, vn) == pu && return t
  exact = setprecision(BigFloat, 256) do
    BigFloat(pu) * BigFloat(vn)
  end
  for n in 18:40
    t = Printf.format(Printf.Format("%.$(n)g"), exact)
    _measurement_table_pu(t, parse(Float64, t), vn) == pu && return t
  end
  throw(ArgumentError("no kV text gives back $(pu) pu at Vn = $(vn) kV"))
end

"""
    writeMeasurementTable(file, rows; comments = String[])

Write `rows` (a vector of `MeasurementRow` or of NamedTuples with the header
fields) as a measurement table: the version line
`# sparlectra-measurement-table v1` first, then `comments` as `#` lines,
then the header. Numbers use the shortest round-trip form,
`sigma = nothing` an empty cell. A `MeasurementRow` read from a file keeps
the number text it was read with, and a number given as a string is
written as it is (checked to be a number, `ArgumentError` otherwise); the
converter uses that for `U` values longer than a Float64 prints. Text cells are written without surrounding
blanks, as the reader reads them; one with a comma or a double quote is
quoted (RFC 4180); a cell or comment with a line break throws an
`ArgumentError` before the file is written, since the table has one row per
line.
"""
function writeMeasurementTable(file::AbstractString, rows; comments::Vector{String} = String[])
  # every cell is formed before the file is opened, so a refused cell
  # leaves no half-written file
  for c in comments
    any(ch -> ch == '\n' || ch == '\r', c) && throw(ArgumentError("writeMeasurementTable: a comment must not contain a line break, got $(repr(c))"))
  end
  lines = String[]
  for r in rows
    cells = String[
      _measurement_table_cell(r.timestamp),
      _measurement_table_cell(r.measurement_mrid),
      _measurement_table_cell(r.terminal_mrid),
      _measurement_table_cell(r.node_mrid),
      _measurement_table_cell(r.type),
      _measurement_table_number(r, :value),
      _measurement_table_number(r, :sigma),
      _measurement_table_cell(r.source),
      r.valid ? "true" : "false",
    ]
    push!(lines, join(cells, ","))
  end
  open(file, "w") do io
    println(io, MEASUREMENT_TABLE_VERSION_LINE)
    for c in comments
      println(io, "# ", c)
    end
    println(io, join(MEASUREMENT_TABLE_HEADER, ","))
    for l in lines
      println(io, l)
    end
  end
  return file
end

# --- key index ---------------------------------------------------------------

"""
    MeasurementTarget

Where a terminal mRID points in the network: `kind` is `:branch` (a branch
end, `branchIdx` and `direction` set), `:injection` (a load, generator or
equivalent injection; summed per bus), `:shunt` (a shunt bay; summed per
bus) or `:not_modelled` (known to the model data but not part of the solved
network: a switch, skipped or disconnected equipment, a de-energized node;
`detail` names which). `busIdx` is the bus of the terminal, 0 when none.
"""
struct MeasurementTarget
  kind::Symbol
  busIdx::Int
  branchIdx::Int
  direction::Symbol
  detail::String
end

"""
    MeasurementKeyIndex

Resolution of measurement-table keys to network locations: built from the
mRIDs of a CGMES import (`cgmesMeasurementKeys`) or from the element names
of any net (`nameMeasurementKeys`). With `canonical = true` (mRID keys)
all keys are in canonical form (`_measurement_key`: no leading `_`, no
`urn:uuid:`) and so is every key of a row before the lookup; with
`canonical = false` (name keys) keys and row keys are compared as written,
without surrounding whitespace.

- `terminals`: terminal mRID to `MeasurementTarget`.
- `nodes`: TopologicalNode, ConnectivityNode or BusbarSection mRID (or the
  bus name) to bus index (0: the node exists but is not energized in the
  model; `-n`: the key names `n` buses, possible only for names).
- `analogs`: Analog mRID to `(location, type, detail)`; `location` is a
  terminal or node key, `type` one of `"P"`, `"Q"`, `"I"`, `"U"` or `""`
  for a measurement type the table does not carry (named in `detail`).
- `groups`: `(:injection | :shunt, bus)` to the terminal keys of every
  modelled injection or shunt at that bus. A bus value is formed only when
  every member is measured.
- `bus_kV`: nominal voltage per bus index (kV), the base of `U` rows.
- `canonical`: `true` for mRID keys (the default of the five-argument
  constructor), `false` for name keys.
"""
struct MeasurementKeyIndex
  terminals::Dict{String,MeasurementTarget}
  nodes::Dict{String,Int}
  analogs::Dict{String,Tuple{String,String,String}}
  groups::Dict{Tuple{Symbol,Int},Vector{String}}
  bus_kV::Vector{Float64}
  canonical::Bool
end

MeasurementKeyIndex(terminals, nodes, analogs, groups, bus_kV) = MeasurementKeyIndex(terminals, nodes, analogs, groups, bus_kV, true)

"""Canonical form of an mRID key: without `urn:uuid:` and without a leading `_`."""
function _measurement_key(s::AbstractString)::String
  k = strip(s)
  startswith(k, "urn:uuid:") && (k = k[10:end])
  startswith(k, "_") && (k = k[2:end])
  return String(k)
end

# A name key is the element name as written, only without surrounding
# whitespace (the reader strips every cell). The mRID canonicalization does
# not apply: buses "_B9" and "B9" are two buses.
_name_key(s::AbstractString)::String = String(strip(s))

# the key of a table cell in the form of the index
_index_key(index::MeasurementKeyIndex, s::AbstractString)::String = index.canonical ? _measurement_key(s) : _name_key(s)

# --- name index for cases without CGMES ------------------------------------------

# the bus name of every bus index (busDict holds name => index)
function _measurement_bus_names(net::Net)::Vector{String}
  names = fill("", length(net.nodeVec))
  for (name, b) in net.busDict
    1 <= b <= length(names) && (names[b] = name)
  end
  return names
end

# The name keys of a net with the element behind each key, as the named
# tuple (targets, unnamed). `targets` holds (key, target, element, group)
# where element is the branch, prosumer or shunt and group the bus sum the
# element belongs to ((:injection | :shunt, bus), or nothing). `unnamed`
# holds (member, group) for every prosumer or shunt without a name that is
# a member of a bus sum: it has no key a row could carry, so `member` is a
# placeholder that only stands in the group and makes that sum partial. A
# branch without a name has no key and is in no sum. A branch name that several branches share is qualified with the
# branch number, `<name>#<k>`; a prosumer or shunt name that several of them
# share with its bus, `<name>@<bus name>`, and when several share that too
# (two generators of one voltage level at one bus carry the same generated
# name) with their order at the bus, `<name>@<bus name>#<k>`. A key that is
# still not unique resolves to :not_modelled with the reason, so a row on it
# is reported, never assigned to one of them by chance; its element stays a
# member of its bus group, so the bus sum is never formed without it.
function _name_measurement_targets(net::Net)
  key = _name_key
  iso = _unfed_isolated_buses(net)
  energized(b) = b > 0 && !(b in iso)
  busname = _measurement_bus_names(net)
  cand = Tuple{String,MeasurementTarget,Any,Union{Nothing,Tuple{Symbol,Int}}}[]

  # branches: plain name when unique, otherwise qualified with the branch
  # number (the position, as branch_nr of the CSV v1); generated names of
  # parallel lines repeat (B_ACL_110_1_2 twice in sp_case14)
  nbranch = Dict{String,Int}()
  for br in net.branchVec
    n = key(br.comp.cName)
    nbranch[n] = get(nbranch, n, 0) + 1
  end
  for (k, br) in enumerate(net.branchVec)
    name = key(br.comp.cName)
    # no key for a branch without a name; it belongs to no bus sum
    isempty(name) && continue
    nbranch[name] > 1 && (name = string(name, "#", k))
    state = _branch_terminal_state(br)
    for (dir, bus) in ((:from, Int(br.fromBus)), (:to, Int(br.toBus)))
      target = if state !== :closed
        MeasurementTarget(:not_modelled, bus, 0, :none, "branch $(name) is not closed at both ends")
      elseif !energized(Int(br.fromBus)) || !energized(Int(br.toBus))
        MeasurementTarget(:not_modelled, bus, 0, :none, "branch $(name) at a bus that is not energized")
      else
        MeasurementTarget(:branch, bus, k, dir, "branch")
      end
      push!(cand, (string(name, ":", dir), target, br, nothing))
    end
  end

  # prosumers and shunts: plain name when it is unique among them
  bays = Tuple{String,Int,Symbol,Any}[]
  for ps in net.prosumpsVec
    push!(bays, (key(ps.comp.cName), getPosumerBusIndex(ps), :injection, ps))
  end
  for sh in net.shuntVec
    push!(bays, (key(sh.comp.cName), sh.busIdx, :shunt, sh))
  end
  nbay = Dict{String,Int}()
  for (name, _, _, _) in bays
    nbay[name] = get(nbay, name, 0) + 1
  end
  # the bus-qualified form, and how many elements share it
  qualified(name, bus) = nbay[name] == 1 ? name : string(name, "@", key(busname[bus]))
  nqual = Dict{String,Int}()
  for (name, bus, _, _) in bays
    isempty(name) && continue
    q = qualified(name, bus)
    nqual[q] = get(nqual, q, 0) + 1
  end
  seen_qual = Dict{String,Int}()
  unnamed = Tuple{String,Tuple{Symbol,Int}}[]
  for (i, (name, bus, kind, el)) in enumerate(bays)
    active = kind === :shunt ? el.status == 1 : true
    if isempty(name)
      # no key, but still a member of the sum at its bus: the placeholder
      # is no terminal key, so a row can never resolve to it
      (energized(bus) && active) && push!(unnamed, (string("(unnamed ", kind, " ", i, " at bus ", key(busname[bus]), ")"), (kind, bus)))
      continue
    end
    k = qualified(name, bus)
    if nqual[k] > 1
      # order at the bus (prosumers first, then shunts, each in model order)
      seen_qual[k] = get(seen_qual, k, 0) + 1
      k = string(k, "#", seen_qual[k])
    end
    target = if !energized(bus)
      MeasurementTarget(:not_modelled, bus, 0, :none, "$(kind) $(k) at a bus that is not energized")
    elseif !active
      MeasurementTarget(:not_modelled, bus, 0, :none, "$(kind) $(k) out of service")
    else
      MeasurementTarget(kind, bus, 0, :none, kind === :shunt ? "shunt" : (isGenerator(el) ? "generator" : "load"))
    end
    push!(cand, (k, target, el, target.kind === kind ? (kind, bus) : nothing))
  end

  # bus links (switches, couplers): a link is no branch of the model (it
  # joins two buses outside the Y-bus); P and Q at its ends are link flow
  # measurements, weighted WLS rows of the estimator (:link targets,
  # branchIdx = the link's position in net.linkVec). An open link resolves
  # to :not_modelled. A link key never shadows a branch or prosumer key of
  # the same text.
  taken = Set{String}(k for (k, _, _, _) in cand)
  for (li, l) in enumerate(net.linkVec)
    name = key(l.cName)
    isempty(name) && continue
    for (dir, bus) in ((:from, Int(l.fromBus)), (:to, Int(l.toBus)))
      k = string(name, ":", dir)
      k in taken && continue
      push!(taken, k)
      target = l.status == 1 ? MeasurementTarget(:link, bus, li, dir, "bus link $(name)") : MeasurementTarget(:not_modelled, bus, 0, :none, "bus link $(name) is open")
      push!(cand, (k, target, l, nothing))
    end
  end

  nkey = Dict{String,Int}()
  for (k, _, _, _) in cand
    nkey[k] = get(nkey, k, 0) + 1
  end
  out = Tuple{String,MeasurementTarget,Any,Union{Nothing,Tuple{Symbol,Int}}}[]
  for (k, t, el, g) in cand
    if nkey[k] > 1
      push!(out, (k, MeasurementTarget(:not_modelled, t.busIdx, 0, :none, "name $(k) is shared by $(nkey[k]) elements"), el, g))
    else
      push!(out, (k, t, el, g))
    end
  end
  return (targets = out, unnamed = unnamed)
end

"""
    nameMeasurementKeys(net::Net) -> MeasurementKeyIndex

The key index of a net without CGMES identities (MATPOWER, SCF, DTF,
PowSyBl, a net built in code) for `importMeasurementTable`. The table keeps
its header; its key columns hold element names:

- `node_mrid`: the bus name. A bus that is not energized resolves to bus 0;
  a bus that only closed links join to a fed bus counts as energized.
- `terminal_mrid`: `<branch name>:from` or `<branch name>:to` for a branch
  end, the prosumer name for a load or generator (member of the injection
  sum at its bus), the shunt name for a shunt (member of the shunt sum).
- A branch name that several branches share (parallel lines) is qualified
  with the branch number, `<branch name>#<k>:from`, `k` being the position
  in `net.branchVec` (the `branch_nr` of the measurement CSV v1). A
  prosumer or shunt name that several of them share is qualified with its
  bus, `<name>@<bus name>`; when several at one bus share that as well,
  with their order at the bus, `<name>@<bus name>#<k>` (prosumers first,
  then shunts, each in model order). Both `#<k>` forms depend on the
  element order of the import and are not stable across imports; unique
  names are. A key that is still not unique resolves to `:not_modelled`
  with the reason, as do branches that are not closed at both ends, shunts
  out of service and elements at a bus that is not energized.
- Every prosumer and in-service shunt at an energized bus is a member of the
  sum at its bus (`groups`), whatever its key; a member without a row
  leaves the sum partial. That includes an element without a name: it has
  no key, so its bus sum is always partial.
- `measurement_mrid` does not apply (the index has no Analogs); a row with
  it is reported as `:unknown_key`.
- A bus link has the keys `<link name>:from` and `<link name>:to`; P and Q
  rows on them are link flow measurements (positive from the link's from
  bus to its to bus), WLS rows of the estimator; `calcLinkFlowsSE!` gives
  the link flows after the estimate. An open link resolves to
  `:not_modelled`, and a branch or prosumer that already holds the key
  keeps it.

Names are compared as written, without surrounding whitespace; the mRID
canonical form (no leading `_`, no `urn:uuid:`) does not apply, so `_B9`
and `B9` are two keys. Bus names that still give one key (`B7` and ` B7`)
make that key ambiguous: a `U` row on it is reported as `:not_modelled`.
"""
function nameMeasurementKeys(net::Net)::MeasurementKeyIndex
  key = _name_key
  iso = _unfed_isolated_buses(net)
  nodes = Dict{String,Int}()
  nbus = Dict{String,Int}()
  for (name, b) in net.busDict
    k = key(name)
    nbus[k] = get(nbus, k, 0) + 1
    nodes[k] = (b > 0 && !(b in iso)) ? b : 0
  end
  # a key several buses share resolves to none of them
  for (k, n) in nbus
    n > 1 && (nodes[k] = -n)
  end
  terminals = Dict{String,MeasurementTarget}()
  groups = Dict{Tuple{Symbol,Int},Vector{String}}()
  named = _name_measurement_targets(net)
  for (k, t, _, g) in named.targets
    terminals[k] = t
    # membership comes from the element, not from the key: an element whose
    # key is still shared stays a member, so its bus sum can only be partial
    if g !== nothing
      members = get!(Vector{String}, groups, g)
      k in members || push!(members, k)
    end
  end
  # a member without a name has no key; its placeholder makes the sum
  # partial, so it is never formed without that member
  for (m, g) in named.unnamed
    push!(get!(Vector{String}, groups, g), m)
  end
  bus_kV = Float64[nd.comp.cVN for nd in net.nodeVec]
  return MeasurementKeyIndex(terminals, nodes, Dict{String,Tuple{String,String,String}}(), groups, bus_kV, false)
end

# the sigma keywords of both table generators: every written sigma has to
# pass the reader (finite and positive), so a zero floor or a zero relative
# U sigma is refused instead of writing rows the reader rejects as bad_sigma
function _check_table_sigmas(fname::AbstractString, sigma_pq_rel::Real, sigma_pq_floor::Real, sigma_u_rel::Real, bad_data_sigmas::Real)
  (isfinite(sigma_pq_rel) && sigma_pq_rel >= 0) || throw(ArgumentError("$(fname): sigma_pq_rel must be finite and >= 0, got $(sigma_pq_rel)"))
  (isfinite(sigma_pq_floor) && sigma_pq_floor > 0) || throw(ArgumentError("$(fname): sigma_pq_floor must be finite and > 0, got $(sigma_pq_floor)"))
  (isfinite(sigma_u_rel) && sigma_u_rel > 0) || throw(ArgumentError("$(fname): sigma_u_rel must be finite and > 0, got $(sigma_u_rel)"))
  # a NaN gross error would write rows the reader rejects as bad_value
  (isfinite(bad_data_sigmas) && bad_data_sigmas >= 0) || throw(ArgumentError("$(fname): bad_data_sigmas must be finite and >= 0, got $(bad_data_sigmas)"))
  return nothing
end

# --- measurement sigma by voltage level -----------------------------------------
#
# state_estimation.sigma_by_voltage_level (SESigmaLevel rows, configuration.jl)
# gives sigmas in physical units per nominal voltage. The table import uses
# them for rows without a sigma cell, the three table generators for the
# sigma they write and the noise they draw, so generated values and the
# estimate agree. Units are those of the table: P in MW, Q in MVAr, I in A,
# U in kV. The estimator keeps P and Q in MW/MVAr and I in A, so those
# sigmas enter as they are; U is turned into pu of the bus (sigma_kV / Vn).

# the rows as SESigmaLevel, from the configured vector or any form the
# configuration parser reads (a mapping, a vector of rows)
_sigma_levels(x::Vector{SESigmaLevel}) = x
_sigma_levels(x) = _parse_sigma_by_voltage_level(x)

# the row of a nominal voltage, or nothing when the table has none
function _sigma_level_row(levels::Vector{SESigmaLevel}, vn_kv::Real)::Union{Nothing,SESigmaLevel}
  for r in levels
    abs(r.vn_kv - vn_kv) < SE_SIGMA_LEVEL_VN_TOL_KV && return r
  end
  return nothing
end

# the absolute column of a table type ("P", "Q", "I", "U"); nothing for a
# type the table has no column for
function _sigma_level_column(r::SESigmaLevel, typ::AbstractString)::Union{Nothing,Float64}
  typ == "P" && return r.p_mw
  typ == "Q" && return r.q_mvar
  typ == "I" && return r.i_a
  typ == "U" && return r.u_kv
  return nothing
end

# The level sigma of one value in table units: max(absolute, rel_pct/100 *
# |value|). Nothing when the level has no row or the row leaves this type
# on the default; with an empty table always nothing, so the callers keep
# their former sigma bit for bit.
function _sigma_from_level(levels::Vector{SESigmaLevel}, vn_kv::Real, typ::AbstractString, value::Real)::Union{Nothing,Float64}
  isempty(levels) && return nothing
  r = _sigma_level_row(levels, vn_kv)
  r === nothing && return nothing
  a = _sigma_level_column(r, typ)
  a === nothing && return nothing
  return max(a, r.rel_pct / 100.0 * abs(Float64(value)))
end

# :default is weaker than :level, :level weaker than :row; a bus sum (or a
# merged member) takes the weakest of its parts
_sigma_origin_rank(o::Symbol) = o === :row ? 1 : o === :level ? 2 : 3
_weakest_sigma_origin(origins)::Symbol = reduce((a, b) -> _sigma_origin_rank(b) > _sigma_origin_rank(a) ? b : a, origins; init = :row)

# the table type of a measurement type, nothing for the angle types
function _sigma_table_type(typ::MeasurementType)::Union{Nothing,String}
  typ in (PinjMeas, PflowMeas) && return "P"
  typ in (QinjMeas, QflowMeas, ShuntQMeas) && return "Q"
  typ == ImagMeas && return "I"
  typ == VmMeas && return "U"
  return nothing
end

# The bus that decides the level of a measurement: its bus, the bus at the
# measured end of a branch, the from bus of a link; 0 without a location.
function _sigma_level_bus(net::Net, m::Measurement)::Int
  m.busIdx !== nothing && return Int(m.busIdx)
  if m.branchIdx !== nothing && 1 <= m.branchIdx <= length(net.branchVec)
    br = net.branchVec[m.branchIdx]
    return Int(m.direction === :to ? br.toBus : br.fromBus)
  end
  if m.linkIdx !== nothing && 1 <= m.linkIdx <= length(net.linkVec)
    return Int(net.linkVec[m.linkIdx].fromBus)
  end
  return 0
end

# The level sigma for measurements made by generateMeasurementsFromPF (the
# Web UI generator). `meas` are the generated rows, `truth` the same call
# without noise (same rows, same order). A row at a configured level gets
# the level sigma of its noise-free value, and its noise is the same draw
# scaled to that sigma: value = truth + (value - truth) * sigma_new /
# sigma_old. Vm is in pu, so the kV sigma is divided by the bus Vn; P and Q
# stay in MW/MVAr, I in A. Rows without a level column keep sigma and value;
# an empty table returns `meas` itself, so the generator without the key
# writes the same set bit for bit. Returns the rows and the indices of the
# rows that took a level sigma.
function _apply_sigma_levels_to_generated(net::Net, meas::Vector{Measurement}, truth::Vector{Measurement}, levels)::Tuple{Vector{Measurement},Vector{Int}}
  rows = _sigma_levels(levels)
  changed = Int[]
  isempty(rows) && return (meas, changed)
  length(truth) == length(meas) || throw(ArgumentError("_apply_sigma_levels_to_generated: $(length(truth)) truth rows for $(length(meas)) measurements"))
  out = copy(meas)
  for (k, (m, t)) in enumerate(zip(meas, truth))
    (m.typ == t.typ && m.id == t.id) || throw(ArgumentError("_apply_sigma_levels_to_generated: row $(k) is $(m.id), the truth row $(t.id)"))
    typ = _sigma_table_type(m.typ)
    typ === nothing && continue
    bus = _sigma_level_bus(net, m)
    1 <= bus <= length(net.nodeVec) || continue
    vn = net.nodeVec[bus].comp.cVN
    # the table works in kV for U; the generated Vm row is in pu
    s = _sigma_from_level(rows, vn, typ, typ == "U" ? t.value * vn : t.value)
    s === nothing && continue
    σ = typ == "U" ? s / vn : s
    z = t.value + (m.value - t.value) * (σ / m.sigma)
    out[k] = Measurement(typ = m.typ, value = z, sigma = σ, active = m.active, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id, linkIdx = m.linkIdx, zero_injection = m.zero_injection)
    push!(changed, k)
  end
  return (out, changed)
end

"""
    _sigma_level_summary(net, measurements; origins = Dict{String,Symbol}(), levels = <state_estimation.sigma_by_voltage_level>) -> NamedTuple

Per nominal voltage of `net`: how many measurements sit at that level and
where their sigma came from (`origins`, the `sigma_origin` of a table set;
a measurement without an entry carries its own sigma, as every row of a
measurement CSV v1 does), plus the absolute sigmas of a configured row in
pu: `p_mw / baseMVA`, `q_mvar / baseMVA`, `u_kv / Vn` and
`i_a / I_base` with `I_base = baseMVA / (sqrt(3) * Vn)` (kA, so
`i_a * sqrt(3) * Vn / (1000 * baseMVA)`). Zero-injection rows (id prefix
`ZI`) are constraints and not counted.

Returns `(levels, unused, unlocated)`: `levels` one NamedTuple per nominal
voltage, highest first, with `vn_kv`, `configured`, `measurements`,
`row`, `level`, `default`, `rel_pct` and `sigma_pu = (p, q, i, u)`
(`nothing` where the row leaves the type on the default); `unused` the
configured `vn_kv` that no bus of `net` has; `unlocated` the measurements
without a bus, branch or link.
"""
function _sigma_level_summary(net::Net, measurements::AbstractVector; origins::AbstractDict = Dict{String,Symbol}(), levels = state_estimation_config().sigma_by_voltage_level)
  rows = _sigma_levels(levels)
  vns = sort!(unique(Float64[nd.comp.cVN for nd in net.nodeVec]); rev = true)
  counts = Dict{Float64,Vector{Int}}(v => [0, 0, 0, 0] for v in vns)
  unlocated = 0
  for m in measurements
    m.zero_injection && continue
    bus = _sigma_level_bus(net, m)
    if !(1 <= bus <= length(net.nodeVec))
      unlocated += 1
      continue
    end
    c = counts[net.nodeVec[bus].comp.cVN]
    c[1] += 1
    c[1 + _sigma_origin_rank(get(origins, m.id, :row))] += 1
  end
  out = NamedTuple[]
  for v in vns
    r = _sigma_level_row(rows, v)
    ibase_ka = net.baseMVA / (sqrt(3.0) * v)
    pu(x, base) = x === nothing ? nothing : x / base
    sigma_pu = r === nothing ? (p = nothing, q = nothing, i = nothing, u = nothing) : (p = pu(r.p_mw, net.baseMVA), q = pu(r.q_mvar, net.baseMVA), i = pu(r.i_a, 1000.0 * ibase_ka), u = pu(r.u_kv, v))
    c = counts[v]
    push!(out, (vn_kv = v, configured = r !== nothing, measurements = c[1], row = c[2], level = c[3], default = c[4], rel_pct = r === nothing ? 0.0 : r.rel_pct, sigma_pu = sigma_pu))
  end
  unused = Float64[r.vn_kv for r in rows if !any(v -> abs(v - r.vn_kv) < SE_SIGMA_LEVEL_VN_TOL_KV, vns)]
  return (levels = out, unused = unused, unlocated = unlocated)
end

"""
    generateMeasurementTableFromPF(net::Net; kwargs...) -> Vector{NamedTuple}

Synthetic measurement table from the solved state of `net` (call `runpf!`
first), keyed by element names (`nameMeasurementKeys`): U at the buses and
P and Q at both ends of the branches, computed with the network model at the
solved voltages. Works for every net, whatever format it came from; for a
CGMES delivery keyed by mRID use `generateMeasurementTableFromSV`. Rows are
sorted by key, so the output does not depend on element order.

Keywords:
- `timestamp = ""`: the timestamp column of every row.
- `coverage = 1.0`: share of the measurement locations that get a row,
  drawn separately for the energized buses (U) and the branch ends (P and
  Q together). `1.0` measures every location.
- `include_injections = false`: also write P and Q of every prosumer (power
  into the equipment: a load positive, a generator negative; the solved
  values of `distributeBusResults!`, computed on a copy of `net`) and Q of
  every shunt (power into the shunt), so the import forms the bus sums.
- `sigma_pq_rel = 0.01`, `sigma_pq_floor = 0.1`: sigma of P and Q as
  `max(sigma_pq_rel * |value|, sigma_pq_floor)` in MW/MVAr.
- `sigma_u_rel = 0.002`: sigma of U relative to the reading (kV).
  `sigma_pq_floor` and `sigma_u_rel` must be positive, `sigma_pq_rel` and
  `bad_data_sigmas` finite and at least 0 (an `ArgumentError` otherwise),
  so every row passes the reader.
- `sigma_by_voltage_level` (default: the active
  `state_estimation.sigma_by_voltage_level`): where the level of the
  measuring bus (the bus of the node, of the branch end, of the prosumer
  or shunt) has a row with a column for the type, that sigma
  (`max(absolute, rel_pct/100 * |value|)`, MW, MVAr, kV) replaces the
  three keywords above, for the written sigma and the noise alike.
- `noise = true`: add Gaussian noise with that sigma; `false` writes the
  model values unchanged.
- `bad_data_fraction = 0.0`: share of rows that get a gross error of
  `bad_data_sigmas` sigma (random sign) instead of the noise. They carry
  `source = "synthetic;bad_data"`, all others `source = "synthetic"`.
- `bad_data_sigmas = 20.0`, `seed = 1`: size of the gross error and the
  seed of the random stream (coverage draw, noise, gross errors); equal
  arguments give equal rows.

Keys that `nameMeasurementKeys` marks `:not_modelled` (shared names, open
branches) get no row. The rows have the fields of the table header and can
be written with `writeMeasurementTable`.
"""
function generateMeasurementTableFromPF(net::Net; timestamp::AbstractString = "", coverage::Real = 1.0, include_injections::Bool = false, sigma_pq_rel::Real = 0.01, sigma_pq_floor::Real = 0.1, sigma_u_rel::Real = 0.002, sigma_by_voltage_level = state_estimation_config().sigma_by_voltage_level, noise::Bool = true, bad_data_fraction::Real = 0.0, bad_data_sigmas::Real = 20.0, seed::Integer = 1)
  0.0 <= coverage <= 1.0 || throw(ArgumentError("generateMeasurementTableFromPF: coverage must be in [0, 1], got $(coverage)"))
  0.0 <= bad_data_fraction <= 1.0 || throw(ArgumentError("generateMeasurementTableFromPF: bad_data_fraction must be in [0, 1], got $(bad_data_fraction)"))
  _check_table_sigmas("generateMeasurementTableFromPF", sigma_pq_rel, sigma_pq_floor, sigma_u_rel, bad_data_sigmas)
  levels = _sigma_levels(sigma_by_voltage_level)
  key = _name_key
  iso = _unfed_isolated_buses(net)
  V = buildVoltageVector(net)
  busname = _measurement_bus_names(net)
  targets = _name_measurement_targets(net).targets
  rng = Random.Xoshiro(seed)
  # the measured locations: energized buses (by name) and branch ends,
  # each drawn with the coverage share
  pick(xs) = xs[sort!(Random.randperm(rng, length(xs))[1:round(Int, Float64(coverage) * length(xs))])]
  buses = pick(sort!([b for b in eachindex(V) if !(b in iso) && !isempty(busname[b])]; by = b -> key(busname[b])))
  ends = pick(sort!([(k, t) for (k, t, _, _) in targets if t.kind === :branch]; by = first))

  # truth values: (key kind, key, type, value, measuring bus)
  truth = Tuple{Symbol,String,String,Float64,Int}[]
  for b in buses
    push!(truth, (:node, key(busname[b]), "U", abs(V[b]) * net.nodeVec[b].comp.cVN, b))
  end
  for (k, t) in ends
    br = net.branchVec[t.branchIdx]
    s = (t.direction === :from ? branchFlow_pu(br, Int(br.fromBus), Int(br.toBus), 1, V) : branchFlow_pu(br, Int(br.toBus), Int(br.fromBus), 2, V)) * net.baseMVA
    push!(truth, (:terminal, k, "P", real(s), t.busIdx))
    push!(truth, (:terminal, k, "Q", imag(s), t.busIdx))
  end
  if include_injections
    # the solved P and Q per prosumer from distributeBusResults!, on a copy
    # so the caller's net keeps its own result fields; per bus they add up
    # to the injection of the solved state
    solved = deepcopy(net)
    distributeBusResults!(solved)
    pos = IdDict{Any,Int}(ps => i for (i, ps) in enumerate(net.prosumpsVec))
    for (k, t, el, _) in targets
      if t.kind === :injection
        ps = solved.prosumpsVec[pos[el]]
        # injection sign (generation positive, consumption negative) turned
        # into power into the equipment; distributeBusResults! writes every
        # prosumer at an energized bus, the only ones with this target kind
        push!(truth, (:terminal, k, "P", -something(ps.pRes, 0.0), t.busIdx))
        push!(truth, (:terminal, k, "Q", -something(ps.qRes, 0.0), t.busIdx))
      elseif t.kind === :shunt
        # the analytic bay power of the ShuntQMeas prediction (q_shunt sign)
        push!(truth, (:terminal, k, "Q", imag(abs2(V[t.busIdx]) * conj(el.y_pu_shunt) * net.baseMVA), t.busIdx))
      end
    end
  end
  sort!(truth; by = x -> (x[2], x[3]))

  n = length(truth)
  nbad = round(Int, Float64(bad_data_fraction) * n)
  bad = Set(Random.randperm(rng, n)[1:nbad])
  Row = NamedTuple{(:timestamp, :measurement_mrid, :terminal_mrid, :node_mrid, :type, :value, :sigma, :source, :valid),Tuple{String,String,String,String,String,Float64,Float64,String,Bool}}
  rows = Row[]
  for (i, (kind, k, typ, v, bus)) in enumerate(truth)
    # the level of the measuring bus first, the keywords where it has none
    lv = _sigma_from_level(levels, net.nodeVec[bus].comp.cVN, typ, v)
    sigma = lv !== nothing ? lv : typ == "U" ? Float64(sigma_u_rel) * abs(v) : max(Float64(sigma_pq_rel) * abs(v), Float64(sigma_pq_floor))
    if i in bad
      value = v + (rand(rng, Bool) ? 1.0 : -1.0) * Float64(bad_data_sigmas) * sigma
      source = "synthetic;bad_data"
    else
      value = noise ? v + sigma * randn(rng) : v
      source = "synthetic"
    end
    push!(rows, Row((String(timestamp), "", kind === :terminal ? k : "", kind === :node ? k : "", typ, value, sigma, source, true)))
  end
  return rows
end

# --- pseudo measurements from the scheduled values ---------------------------------

# The scheduled P and Q of a prosumer as power into the equipment (the sign
# of the table rows: a load positive, a generator negative). pVal and qVal
# hold the setpoints of the case (for a CGMES import the SSH values, with
# the generator sign turned to injection by the importer). A missing value
# stays `nothing`: a regulated unit without a Q setpoint gets no Q row.
function _scheduled_into_equipment(ps)::NTuple{2,Union{Nothing,Float64}}
  s = isGenerator(ps) ? -1.0 : 1.0
  p = ps.pVal === nothing ? nothing : s * Float64(ps.pVal)
  q = ps.qVal === nothing ? nothing : s * Float64(ps.qVal)
  return (p, q)
end

# The (terminal key, type) pairs a table observes, per timestamp: every
# well-formed terminal row and every Analog row that resolves to a terminal
# with a type. A valid=false row counts as well: the meter exists, and a
# second row on the same terminal key would be a duplicate in the import.
function _observed_terminal_types(rows, index::MeasurementKeyIndex)::Dict{String,Set{Tuple{String,String}}}
  seen = Dict{String,Set{Tuple{String,String}}}()
  for r in rows
    hasproperty(r, :problem) && r.problem !== :none && continue
    obs = get!(Set{Tuple{String,String}}, seen, String(r.timestamp))
    if !isempty(r.measurement_mrid)
      a = get(index.analogs, _index_key(index, r.measurement_mrid), nothing)
      a === nothing && continue
      loc, typ, _ = a
      (isempty(typ) || isempty(loc)) || push!(obs, (loc, typ))
    elseif !isempty(r.terminal_mrid)
      push!(obs, (_index_key(index, r.terminal_mrid), String(r.type)))
    end
  end
  return seen
end

# the timestamps a pseudo table is written for: the given one, or every
# timestamp of the table in order of first appearance ("" for an empty table)
function _pseudo_timestamps(rows, timestamp::Union{Nothing,AbstractString})::Vector{String}
  timestamp === nothing || return [String(timestamp)]
  out = String[]
  for r in rows
    hasproperty(r, :problem) && r.problem in (:field_count, :bad_timestamp) && continue
    ts = String(r.timestamp)
    ts in out || push!(out, ts)
  end
  isempty(out) && push!(out, "")
  return out
end

# Shared core of both generatePseudoMeasurementTable methods. `schedule`
# maps the terminal key of each injection member to its scheduled
# (P, Q) into the equipment. Members of an injection group without an entry
# (a name key that is still shared) cannot be addressed by a row; they are
# counted and named in one warning, because their bus sum stays partial.
function _pseudo_measurement_rows(fname::AbstractString, rows, index::MeasurementKeyIndex, schedule::Dict{String,NTuple{2,Union{Nothing,Float64}}}, sigma::Real, sigma_q::Real, timestamp, sigma_by_voltage_level)
  (isfinite(sigma) && sigma > 0) || throw(ArgumentError("$(fname): sigma must be finite and > 0, got $(sigma)"))
  (isfinite(sigma_q) && sigma_q > 0) || throw(ArgumentError("$(fname): sigma_q must be finite and > 0, got $(sigma_q)"))
  levels = _sigma_levels(sigma_by_voltage_level)
  observed = _observed_terminal_types(rows, index)
  Row = NamedTuple{(:timestamp, :measurement_mrid, :terminal_mrid, :node_mrid, :type, :value, :sigma, :source, :valid),Tuple{String,String,String,String,String,Float64,Float64,String,Bool}}
  out = Row[]
  unaddressable = Set{String}()
  groups = sort!([(g, sort(m)) for (g, m) in index.groups if g[1] === :injection]; by = x -> x[1][2])
  for ts in _pseudo_timestamps(rows, timestamp)
    obs = get(observed, ts, Set{Tuple{String,String}}())
    for (g, members) in groups, member in members
      sched = get(schedule, member, nothing)
      if sched === nothing
        push!(unaddressable, member)
        continue
      end
      # the level of the bus of the sum, else the keywords
      vn = index.bus_kV[g[2]]
      for (typ, v, sg) in (("P", sched[1], Float64(sigma)), ("Q", sched[2], Float64(sigma_q)))
        (v === nothing || (member, typ) in obs) && continue
        lv = _sigma_from_level(levels, vn, typ, v)
        push!(out, Row((ts, "", member, "", typ, v, lv === nothing ? sg : lv, "pseudo;ssh", true)))
      end
    end
  end
  isempty(unaddressable) || @warn "$(fname): $(length(unaddressable)) injection member(s) have no unique key and get no pseudo row; their bus sums stay partial: $(join(sort!(collect(unaddressable)), ", "))"
  return out
end

"""
    generatePseudoMeasurementTable(net::Net, rows; sigma, sigma_q = sigma, timestamp = nothing, sigma_by_voltage_level = <configured>) -> Vector{NamedTuple}
    generatePseudoMeasurementTable(result::CGMESImportResult, rows; sigma, sigma_q = sigma, timestamp = nothing, sigma_by_voltage_level = <configured>) -> Vector{NamedTuple}

Pseudo measurements for the loads and generators that a measurement table
does not observe. For every member of an injection bus sum (a load,
generator, external grid, equivalent injection or SVC) without a `P` row
in `rows`, a `P` row with its scheduled value is written; the same for
`Q`. The scheduled values are the setpoints of the case (`pVal`, `qVal` of
the prosumer; for a CGMES import the SSH values), as power into the
equipment like every table row. A unit without a Q setpoint gets no `Q`
row.

- `net::Net`: keys are element names (`nameMeasurementKeys`).
- `result::CGMESImportResult`: keys are terminal mRIDs
  (`cgmesMeasurementKeys`).

`rows` are the rows of the table (`readMeasurementTable` or a generator).
A row observes a member when it is well formed and carries that terminal
key (or an Analog on that terminal) with the same type, `valid=false`
included. Members are filled one by one, so a bus whose sum was partial
becomes complete and its real rows are kept.

Keywords:
- `sigma`: sigma of the `P` rows in MW, required, finite and positive
  (otherwise an `ArgumentError`). Choose it large against the real meters
  so they win.
- `sigma_q = sigma`: sigma of the `Q` rows in MVAr, same check.
- `timestamp = nothing`: write rows for every timestamp of `rows` (an
  empty table gives one set with an empty timestamp); a string writes only
  that timestamp.
- `sigma_by_voltage_level` (default: the active
  `state_estimation.sigma_by_voltage_level`): where the level of the bus
  has a row with a `p_mw` (`q_mvar`) column, that sigma
  (`max(absolute, rel_pct/100 * |value|)`) replaces `sigma` (`sigma_q`).
  Those are the sigmas of the meters at that level; pass `SESigmaLevel[]`
  to keep the larger pseudo sigmas.

The rows carry `source = "pseudo;ssh"` and `valid = true`. Append them to
the table (`vcat(rows, pseudo)` or one `writeMeasurementTable` call); the
import treats them like any other row and adds them to the bus sums.
Shunts get no pseudo rows. A member whose name key is shared by several
elements, or that has no name, cannot be addressed by a row; such members
are named in a warning.
"""
function generatePseudoMeasurementTable(net::Net, rows::AbstractVector; sigma::Real, sigma_q::Real = sigma, timestamp::Union{Nothing,AbstractString} = nothing, sigma_by_voltage_level = state_estimation_config().sigma_by_voltage_level)
  schedule = Dict{String,NTuple{2,Union{Nothing,Float64}}}()
  for (k, t, el, _) in _name_measurement_targets(net).targets
    t.kind === :injection && (schedule[k] = _scheduled_into_equipment(el))
  end
  return _pseudo_measurement_rows("generatePseudoMeasurementTable", rows, nameMeasurementKeys(net), schedule, sigma, sigma_q, timestamp, sigma_by_voltage_level)
end

# --- import into measurement sets -----------------------------------------------

"""
    MeasurementSet

The measurements of one timestamp of a measurement table, ready for
`runse!`. `sources` maps each measurement id to the file lines it came from
(an injection or shunt sum has several). `counts` holds `:rows`, `:active` and
`:inactive` (measurements created active or inactive), `:invalid_rows`
(resolved rows with `valid=false`), `:merged` (redundant values of one
group member merged into one), `:sigma_lifted` (row sigmas read at the
floor, the smaller of `ZERO_INJECTION_SIGMA` and the zero-injection
sigma of the import), `:zero_injection` (rows `applyMeasurementSet!` added at passive
buses, present once the set was applied) and one entry per problem; `problems` lists every row that did not become a
measurement as `(line, problem, detail)`. `sigma_origin` maps each
measurement id to where its sigma came from: `:row` (the sigma cell of the
row), `:level` (`state_estimation.sigma_by_voltage_level`) or `:default`
(`default_sigmas`); a bus sum takes the weakest origin of its members
(`:default` before `:level` before `:row`).
"""
struct MeasurementSet
  timestamp::String
  measurements::Vector{Measurement}
  sources::Dict{String,Vector{Int}}
  counts::Dict{Symbol,Int}
  problems::Vector{NamedTuple{(:line, :problem, :detail),Tuple{Int,Symbol,String}}}
  sigma_origin::Dict{String,Symbol}
end

"""
    MeasurementTableImport

Result of `importMeasurementTable`: one `MeasurementSet` per timestamp in
file order of first appearance (an empty timestamp is a set of its own), plus
`rejected` for rows that cannot be assigned to any set (wrong field count,
malformed timestamp).
"""
struct MeasurementTableImport
  sets::Vector{MeasurementSet}
  rejected::Vector{NamedTuple{(:line, :problem, :detail),Tuple{Int,Symbol,String}}}
end

# measurement type of one resolved quantity at a branch end
function _measurement_table_quantity(typ::String, target::MeasurementTarget)
  if target.kind === :branch
    typ == "P" && return PflowMeas
    typ == "Q" && return QflowMeas
    typ == "I" && return ImagMeas
    typ == "Ia" && return IaMeas
  end
  return nothing
end

# measurement type of a group sum: injections give PinjMeas or QinjMeas,
# shunts ShuntQMeas or the bus-referenced bay current ImagMeas
_measurement_table_group_type(kind::Symbol, typ::AbstractString) = kind === :shunt ? (typ == "I" ? ImagMeas : ShuntQMeas) : (typ == "P" ? PinjMeas : QinjMeas)

# `source` is free text; the one token the import evaluates is
# `zero_injection` (among `;`-separated tokens) on a bus balance row
_measurement_table_zero_injection(source::AbstractString)::Bool = any(t -> strip(t) == "zero_injection", split(source, ';'))

# one pending value of a group sum (injection or shunt bay at a bus)
struct _GroupPart
  terminal::String
  value::Float64
  sigma::Float64
  line::Int
  # where the sigma came from (:row, :level, :default)
  origin::Symbol
end

"""
    importMeasurementTable(rows, index; default_sigmas = measurementStdDevs(), sigma_by_voltage_level = <state_estimation.sigma_by_voltage_level>, zero_injection_sigma = <state_estimation.zero_injection_sigma_mva>) -> MeasurementTableImport

Turn the rows of a measurement table into one `MeasurementSet` per
timestamp, resolving every key through `index` (`MeasurementKeyIndex`).

Mapping:
- terminal of a branch end: `P`, `Q`, `I`, `Ia` become `PflowMeas`,
  `QflowMeas`, `ImagMeas`, `IaMeas` (current angle in degrees) at that end. Sign: positive into the equipment at that
  terminal (the CGMES `SvPowerFlow` convention), which is the Sparlectra
  flow direction from the bus into the branch.
- terminals of injections at one bus: the `P` (or `Q`) values of ALL
  modelled injections at the bus are summed into one `PinjMeas` (`QinjMeas`)
  with the sign turned (network injection = minus the power into the
  equipment) and `sigma = sqrt(sum sigma^2)`. When every member has a
  value but not every member a valid one, the sum is created inactive.
  When a member has no row at all, the rows of the sum are counted as
  `:partial_injection`. Several values of one member (redundant Analogs)
  are first merged to their inverse-variance mean (counted as `:merged`).
- terminals of shunts at one bus: `Q` summed the same way into one
  `ShuntQMeas` (power into the shunt, the `q_shunt` convention); `I` (A)
  becomes the bus-referenced bay current `ImagMeas` when the bus has one
  shunt (with several, the rows are `:unsupported`: the bay current is not
  the sum of the magnitudes). `P` at a shunt, `I` at an injection and `Ia`
  at an injection or shunt are counted `:unsupported`.
- node (TopologicalNode, ConnectivityNode, BusbarSection) with `U` in kV:
  `VmMeas` in pu of the bus nominal voltage; sigma likewise.
- node with `P` or `Q`: the bus balance, `PinjMeas` or `QinjMeas` with the
  value as it is (MW, MVAr, injection positive: generation minus load at
  the bus). It stands next to a sum of prosumer rows at the same bus; both
  enter the set as redundant measurements. A balance row whose `source`
  carries the token `zero_injection` (`;`-separated) is a zero-injection
  constraint and gets the id `ZI_PINJ_bus_<k>` / `ZI_QINJ_bus_<k>` of
  `addZeroInjectionMeasurements!`, so it is protected like one.
- node with `Va`: `VaMeas` in degrees.
- `measurement_mrid`: type and location from the Analog in `index`.
- with a name index (`nameMeasurementKeys`) the key columns hold names:
  the bus name in `node_mrid`, `<branch name>:from`/`:to`, the prosumer
  or shunt name in `terminal_mrid`; the mapping above is the same.

Several keys at one location (two Analogs on one terminal, U at every
busbar section of one node) are redundant measurements and all enter the
set; only the same key with the same type twice is a `:duplicate`.

The sigma of a row, in this order: the `sigma` cell of the row; else the
row of `sigma_by_voltage_level` (default: the active
`state_estimation.sigma_by_voltage_level`) for the nominal voltage of the
measuring bus, the bus of the node, of the branch end or of the bus sum
(P, Q in MW/MVAr, I in A, U in kV turned into pu of that bus,
`max(absolute, rel_pct/100 * |value|)`); else `default_sigmas` of the
measurement type (pu for `VmMeas`). `set.sigma_origin` records which of
the three it was. A sigma below the floor, the smaller of
`ZERO_INJECTION_SIGMA` and `zero_injection_sigma` (default: the active
`state_estimation.zero_injection_sigma_mva`), in
the unit of the measurement (pu for `VmMeas`) is read at the floor and counted as
`:sigma_lifted`, as the measurement CSV v1 reader does. Rows with `valid=false` become inactive measurements and make a
group sum inactive when their member has no valid row; next to a valid row
of the same member they are not used and counted as `:invalid_unused`,
also in an inactive sum. Every other row
that does not become a measurement is counted by reason: `:unknown_key`,
`:not_modelled`, `:unsupported`, `:partial_injection`, `:duplicate`, plus
the format problems of `readMeasurementTable`. Nothing is dropped silently.
"""
function importMeasurementTable(rows::Vector{MeasurementRow}, index::MeasurementKeyIndex; default_sigmas::Dict = measurementStdDevs(), sigma_by_voltage_level = state_estimation_config().sigma_by_voltage_level, zero_injection_sigma::Real = state_estimation_config().zero_injection_sigma_mva)::MeasurementTableImport
  return _import_measurement_table(rows, index; default_sigmas = default_sigmas, sigma_by_voltage_level = sigma_by_voltage_level, zero_injection_sigma = zero_injection_sigma)
end

# importMeasurementTable with one addition for the run artifact of the SE
# service: `filled_sigmas`, when given, receives for every row without a
# sigma cell the sigma the import gave it, as the text of a sigma cell in
# the unit of the table (U in kV), with its origin (:level or :default).
# A row written with that text imports to the same sigma, so the artifact
# carries the sigma of every row and reads back the same under any
# configuration.
function _import_measurement_table(rows::Vector{MeasurementRow}, index::MeasurementKeyIndex; default_sigmas::Dict = measurementStdDevs(), sigma_by_voltage_level = state_estimation_config().sigma_by_voltage_level, zero_injection_sigma::Real = state_estimation_config().zero_injection_sigma_mva, filled_sigmas::Union{Nothing,Dict{Int,Tuple{String,Symbol}}} = nothing)::MeasurementTableImport
  levels = _sigma_levels(sigma_by_voltage_level)
  Problem = NamedTuple{(:line, :problem, :detail),Tuple{Int,Symbol,String}}
  rejected = Problem[]
  order = String[]
  by_ts = Dict{String,Vector{MeasurementRow}}()
  for r in rows
    if r.problem in (:field_count, :bad_timestamp)
      push!(rejected, (line = r.line, problem = r.problem, detail = r.detail))
      continue
    end
    # rows built in code may carry another form of one instant
    ts = something(_measurement_table_timestamp(r.timestamp), r.timestamp)
    haskey(by_ts, ts) || push!(order, ts)
    push!(get!(Vector{MeasurementRow}, by_ts, ts), r)
  end
  # the energized buses: every bus a node key resolves to (a node that is
  # not energized resolves to 0)
  energized = Set{Int}(b for b in values(index.nodes) if b > 0)
  # the sigma floor of every row (see _measurement_set)
  sigma_floor = min(ZERO_INJECTION_SIGMA, Float64(zero_injection_sigma))
  sets = MeasurementSet[_measurement_set(ts, by_ts[ts], index, default_sigmas, energized, levels, sigma_floor; filled_sigmas = filled_sigmas) for ts in order]
  return MeasurementTableImport(sets, rejected)
end

function _measurement_set(ts::String, rows::Vector{MeasurementRow}, index::MeasurementKeyIndex, default_sigmas::Dict, energized::Set{Int}, levels::Vector{SESigmaLevel}, sigma_floor::Float64; filled_sigmas::Union{Nothing,Dict{Int,Tuple{String,Symbol}}} = nothing)::MeasurementSet
  meas = Measurement[]
  sources = Dict{String,Vector{Int}}()
  origins = Dict{String,Symbol}()
  counts = Dict{Symbol,Int}(:rows => length(rows), :active => 0, :inactive => 0, :invalid_rows => 0)
  problems = NamedTuple{(:line, :problem, :detail),Tuple{Int,Symbol,String}}[]
  bump!(k, n = 1) = (counts[k] = get(counts, k, 0) + n)
  problem!(r, p, d) = (bump!(p); push!(problems, (line = r.line, problem = p, detail = d)))
  # a duplicate is the same KEY with the same type twice in one set; several
  # keys at one location (two Analogs on a terminal, U at every busbar of a
  # node) are redundant measurements and all stay in the set
  seen_ids = Set{String}()
  # (kind, bus, type) -> the member rows of a group sum with their terminal
  group_rows = Dict{Tuple{Symbol,Int,String},Vector{Tuple{MeasurementRow,String,Float64,Symbol}}}()

  function add!(r::MeasurementRow, typ::MeasurementType, value::Float64, sigma::Float64, origin::Symbol; busIdx = nothing, branchIdx = nothing, linkIdx = nothing, direction = :none, id::String, zero_injection::Bool = false)
    push!(meas, Measurement(typ = typ, value = value, sigma = sigma, active = r.valid, busIdx = busIdx, branchIdx = branchIdx, linkIdx = linkIdx, direction = direction, id = id, zero_injection = zero_injection))
    sources[id] = [r.line]
    origins[id] = origin
    bump!(r.valid ? :active : :inactive)
  end
  default_sigma(typ) = Float64(get(default_sigmas, typ, 1.0))
  # A sigma from the file below the floor (in the unit of the measurement,
  # pu for VmMeas) is read at the floor, as the CSV v1 reader does: a
  # weight 1/sigma^2 of 1e400 overflows to Inf and the estimate reports
  # J = Inf as converged. The floor is the tightest legitimate row, the
  # zero-injection sigma of the caller (state_estimation.zero_injection_sigma_mva
  # of its configuration) where it is set below ZERO_INJECTION_SIGMA, so a
  # ZI row written at that sigma keeps it. A larger setting does not raise
  # the floor: it would lift every real row below it.
  lift(sigma::Float64) = sigma < sigma_floor ? (bump!(:sigma_lifted); sigma_floor) : sigma
  # The sigma of one row in the unit of the table (MW, MVAr, A, kV) with its
  # origin: the row's own cell, else the row of the level of `bus`, else
  # nothing (the caller takes the default of the measurement type). U
  # leaves this in kV; the caller turns it into pu of the bus.
  function row_sigma(r::MeasurementRow, typ::AbstractString, bus::Int)
    r.sigma === nothing || return (r.sigma, :row)
    lv = (bus >= 1 && bus <= length(index.bus_kV)) ? _sigma_from_level(levels, index.bus_kV[bus], typ, r.value) : nothing
    return lv === nothing ? (nothing, :default) : (lv, :level)
  end
  # the sigma a row without a sigma cell got, as cell text (see
  # _import_measurement_table); a row with its own sigma is not recorded
  filled!(r::MeasurementRow, text::String, origin::Symbol) = (filled_sigmas === nothing || r.sigma !== nothing) || (filled_sigmas[r.line] = (text, origin))

  for r in rows
    if r.problem !== :none
      problem!(r, r.problem, r.detail)
      continue
    end
    # resolve the key to (location key, type); `at_node` tells a node key
    # (a node_mrid row, or an Analog whose location is a node) from a
    # terminal key: name keys of a bus and of a prosumer can be equal
    typ = r.type
    if !isempty(r.measurement_mrid)
      a = get(index.analogs, _index_key(index, r.measurement_mrid), nothing)
      if a === nothing
        problem!(r, :unknown_key, "measurement_mrid $(r.measurement_mrid) is no Analog of the model")
        continue
      end
      loc, typ, adetail = a
      if isempty(typ)
        problem!(r, :unknown_type, "Analog $(r.measurement_mrid): $(adetail)")
        continue
      end
      if typ in ("U", "I") && r.value < 0.0
        problem!(r, :bad_value, "Analog $(r.measurement_mrid): value $(r.value) is negative, $(typ) is a magnitude")
        continue
      end
      id = _index_key(index, r.measurement_mrid)
      at_node = haskey(index.nodes, loc) && !haskey(index.terminals, loc)
    else
      at_node = isempty(r.terminal_mrid)
      loc = _index_key(index, at_node ? r.node_mrid : r.terminal_mrid)
      # a bus balance gets its own id form, so P at a bus and P at a
      # prosumer of the same name are two keys, not a duplicate
      id = (at_node && typ in ("P", "Q")) ? string(typ, "bus@", loc) : string(typ, "@", loc)
    end
    if id in seen_ids
      problem!(r, :duplicate, "key $(id) appears twice in the set")
      continue
    end
    push!(seen_ids, id)
    if typ == "U" || (at_node && typ in ("P", "Q", "Va"))
      # a node key; for U also the node of a terminal (an Analog on a
      # terminal). A terminal carries the voltage of its node, also a switch
      # terminal that is :not_modelled itself, but only while that node is
      # energized; a terminal target keeps its real bus, so it is checked here
      bus = get(index.nodes, loc, nothing)
      via = nothing
      if bus === nothing && typ == "U"
        via = get(index.terminals, loc, nothing)
        bus = via === nothing ? nothing : (via.busIdx in energized ? via.busIdx : 0)
      end
      if bus === nothing
        problem!(r, :unknown_key, "node $(loc) is not in the model")
        continue
      end
      if bus < 0
        problem!(r, :not_modelled, "node $(loc) is ambiguous, the name key is shared by $(-bus) buses")
        continue
      end
      if bus == 0
        detail = if via === nothing
          "node $(loc) is not energized in the model"
        elseif via.busIdx == 0
          "terminal $(loc) has no node in the model"
        else
          "terminal $(loc): its node (bus $(via.busIdx)) is not energized in the model"
        end
        problem!(r, :not_modelled, detail)
        continue
      end
      r.valid || bump!(:invalid_rows)
      if typ == "U"
        vn = index.bus_kV[bus]
        s_kv, origin = row_sigma(r, "U", bus)
        # the row's own sigma from its text, a level sigma as a number
        sigma = s_kv === nothing ? default_sigma(VmMeas) : lift(origin === :row ? _measurement_table_pu(r.sigma_text, s_kv, vn) : s_kv / vn)
        filled!(r, s_kv === nothing ? _measurement_table_kv_text(sigma, vn) : string(s_kv), origin)
        add!(r, VmMeas, _measurement_table_pu(r.value_text, r.value, vn), sigma, origin; busIdx = bus, id = id)
      elseif typ == "Va"
        s_row, origin = row_sigma(r, typ, bus)
        add!(r, VaMeas, r.value, s_row === nothing ? default_sigma(VaMeas) : lift(s_row), origin; busIdx = bus, id = id)
        filled!(r, string(something(s_row, default_sigma(VaMeas))), origin)
      else
        # the bus balance: injection positive, the value of PinjMeas as it
        # is. A row marked zero_injection in `source` is the constraint the
        # generator writes at a passive bus; it gets the id of
        # addZeroInjectionMeasurements!, which keeps it out of bad-data
        # elimination and gross errors
        mt = typ == "P" ? PinjMeas : QinjMeas
        zi = _measurement_table_zero_injection(r.source)
        if zi
          zid = string(typ == "P" ? "ZI_PINJ_bus_" : "ZI_QINJ_bus_", bus)
          if zid in seen_ids
            problem!(r, :duplicate, "zero-injection row $(zid) appears twice in the set")
            continue
          end
          push!(seen_ids, zid)
          id = zid
        end
        s_row, origin = row_sigma(r, typ, bus)
        add!(r, mt, r.value, s_row === nothing ? default_sigma(mt) : lift(s_row), origin; busIdx = bus, id = id, zero_injection = zi)
        filled!(r, string(something(s_row, default_sigma(mt))), origin)
      end
      continue
    end
    t = get(index.terminals, loc, nothing)
    if t === nothing
      problem!(r, :unknown_key, "terminal $(loc) is not in the model")
      continue
    end
    if t.kind === :not_modelled
      problem!(r, :not_modelled, "terminal $(loc): $(t.detail)")
      continue
    end
    if t.kind === :branch
      r.valid || bump!(:invalid_rows)
      mt = _measurement_table_quantity(typ, t)
      # the level of a branch end is the bus at that end (on a
      # transformer the voltage of that side)
      s_row, origin = row_sigma(r, typ, t.busIdx)
      sigma = s_row === nothing ? default_sigma(mt) : lift(s_row)
      add!(r, mt, r.value, sigma, origin; branchIdx = t.branchIdx, direction = t.direction, id = id)
      filled!(r, string(something(s_row, default_sigma(mt))), origin)
      continue
    end
    # bus link end: P and Q become a link flow measurement, positive from
    # link.fromBus to link.toBus (a row at the to end flows into the link
    # from the to bus, so its sign turns). The estimator takes it as a WLS
    # row (link rows, _se_link_table)
    if t.kind === :link
      if !(typ in ("P", "Q"))
        problem!(r, :unsupported, "$(typ) at a bus link terminal ($(t.detail)) is not a modelled quantity")
        continue
      end
      r.valid || bump!(:invalid_rows)
      mt = typ == "P" ? PflowMeas : QflowMeas
      s_row, origin = row_sigma(r, typ, t.busIdx)
      sigma = s_row === nothing ? default_sigma(mt) : lift(s_row)
      add!(r, mt, t.direction === :to ? -r.value : r.value, sigma, origin; linkIdx = t.branchIdx, direction = :from, id = id)
      filled!(r, string(something(s_row, default_sigma(mt))), origin)
      continue
    end
    # injection or shunt member: collect for the group sum at the bus. A
    # shunt takes Q and the bay current I; an injection P and Q
    if typ == "Ia" || (t.kind === :injection && typ == "I") || (t.kind === :shunt && typ == "P")
      problem!(r, :unsupported, "$(typ) at a $(t.kind) terminal ($(t.detail)) is not a modelled quantity")
      continue
    end
    r.valid || bump!(:invalid_rows)
    mt = _measurement_table_group_type(t.kind, typ)
    s_row, origin = row_sigma(r, typ, t.busIdx)
    push!(get!(Vector{Tuple{MeasurementRow,String,Float64,Symbol}}, group_rows, (t.kind, t.busIdx, typ)), (r, loc, s_row === nothing ? default_sigma(mt) : lift(s_row), origin))
    filled!(r, string(something(s_row, default_sigma(mt))), origin)
  end

  # group sums. Every modelled member needs a value: all members valid gives
  # an active sum, all members present but some only invalid an inactive sum
  # (the input says "do not use", as for a single branch row), a member
  # without any row leaves the sum partial
  for gkey in sort!(collect(keys(group_rows)); by = k -> (string(k[1]), k[2], k[3]))
    kind, bus, typ = gkey
    members = Set(get(index.groups, (kind, bus), String[]))
    entries = group_rows[gkey]
    any_by_member = Dict{String,Vector{_GroupPart}}()
    valid_by_member = Dict{String,Vector{_GroupPart}}()
    for (r, loc, sigma, origin) in entries
      part = _GroupPart(loc, r.value, sigma, r.line, origin)
      push!(get!(Vector{_GroupPart}, any_by_member, loc), part)
      r.valid && push!(get!(Vector{_GroupPart}, valid_by_member, loc), part)
    end
    absent = setdiff(members, keys(any_by_member))
    # the bay current of several shunts is the magnitude of the sum of
    # their currents, not the sum of the magnitudes
    if typ == "I" && length(members) > 1
      for (r, _, _, _) in entries
        problem!(r, :unsupported, "I at a shunt terminal of bus $(bus): the bus has $(length(members)) shunts, and the bay current is not the sum of their magnitudes")
      end
      continue
    end
    if !isempty(absent)
      for (r, _, _, _) in entries
        problem!(r, :partial_injection, "$(kind) sum at bus $(bus) of type $(typ): $(length(absent)) of $(length(members)) member(s) without a value")
      end
      continue
    end
    active = issubset(members, keys(valid_by_member))
    # each member uses its valid rows, and only a member without any valid
    # row its invalid ones (the sum is then inactive); a valid=false row of a
    # member that also has a valid row is reported, not merged, in an active
    # and in an inactive sum alike
    for (r, loc, _, _) in entries
      (!r.valid && haskey(valid_by_member, loc)) && problem!(r, :invalid_unused, "valid=false next to a valid row of $(loc); the $(kind) sum at bus $(bus) of type $(typ) uses the valid rows of that member")
    end
    got = _GroupPart[]
    for m in sort!(collect(members))
      ps = haskey(valid_by_member, m) ? valid_by_member[m] : any_by_member[m]
      push!(got, _merge_group_parts(ps))
      length(ps) > 1 && bump!(:merged, length(ps) - 1)
    end
    mt = _measurement_table_group_type(kind, typ)
    # network injection is minus the power into the equipment; a shunt
    # measurement stays in the into-the-shunt (q_shunt) convention, its
    # current is a magnitude
    total = sum(p.value for p in got)
    value = kind === :shunt ? total : -total
    sigma = sqrt(sum(p.sigma^2 for p in got))
    id = string(mt, "@bus", bus)
    push!(meas, Measurement(typ = mt, value = value, sigma = sigma, active = active, busIdx = bus, id = id))
    sources[id] = sort!([r.line for (r, loc, _, _) in entries if r.valid || !haskey(valid_by_member, loc)])
    origins[id] = _weakest_sigma_origin(p.origin for p in got)
    bump!(active ? :active : :inactive)
  end
  return MeasurementSet(ts, meas, sources, counts, problems, origins)
end

# several values of one group member (redundant Analogs): inverse-variance
# mean and its sigma
function _merge_group_parts(ps::Vector{_GroupPart})::_GroupPart
  length(ps) == 1 && return ps[1]
  w = [1.0 / p.sigma^2 for p in ps]
  return _GroupPart(ps[1].terminal, sum(w .* [p.value for p in ps]) / sum(w), 1.0 / sqrt(sum(w)), ps[1].line, _weakest_sigma_origin(p.origin for p in ps))
end

"""
    applyMeasurementSet!(net, set; replace = true, zero_injection = <state_estimation.zero_injection is not off>, zero_injection_sigma = <state_estimation.zero_injection_sigma_mva>) -> Net

Put the measurements of `set` into `net.measurements` (replacing them by
default), ready for `runse!(net)`.

With `zero_injection` (default: the active `state_estimation.zero_injection`
is `on` or `constraint`) every passive bus (no generation, load
or shunt, not the slack, not isolated) that has no injection measurement
in `net.measurements` afterwards gets a zero-injection pair `ZI_PINJ_bus_<k>`
/ `ZI_QINJ_bus_<k>` with sigma `zero_injection_sigma` (MW, Mvar), the rows the measurement
generator writes for passive buses. The number of rows added is stored in
`set.counts[:zero_injection]` and printed by `measurementImportReport`.
With `zero_injection = false` no row is added. `set.measurements` itself
is not changed.
"""
function applyMeasurementSet!(net::Net, set::MeasurementSet; replace::Bool = true, zero_injection::Bool = state_estimation_config().zero_injection !== :off, zero_injection_sigma::Real = state_estimation_config().zero_injection_sigma_mva)
  target = _net_measurements(net)
  replace && empty!(target)
  append!(target, set.measurements)
  set.counts[:zero_injection] = zero_injection ? _add_missing_zero_injections!(target, net; sigma = zero_injection_sigma) : 0
  return net
end

"""
    readMeasurementTable!(net; file, index = nameMeasurementKeys(net), timestamp = nothing, replace = true, zero_injection = <state_estimation.zero_injection is not off>, zero_injection_sigma = <state_estimation.zero_injection_sigma_mva>) -> MeasurementSet

Read a measurement table into `net.measurements`: `readMeasurementTable`, then
`importMeasurementTable` with `index` (element names by default; for a
CGMES case `cgmesMeasurementKeys(result)`), then `applyMeasurementSet!` of
the set with `timestamp` (the first set when `nothing`). Returns that set;
rows that did not become a measurement are in its `problems` and in
`measurementImportReport`. Throws `ArgumentError` when the file is no
table the reader accepts, when no row can be assigned to a set, or when the
table has no set with `timestamp`.
"""
function readMeasurementTable!(net::Net; file::AbstractString, index::MeasurementKeyIndex = nameMeasurementKeys(net), timestamp::Union{Nothing,AbstractString} = nothing, replace::Bool = true, zero_injection::Bool = state_estimation_config().zero_injection !== :off, zero_injection_sigma::Real = state_estimation_config().zero_injection_sigma_mva)::MeasurementSet
  imp = importMeasurementTable(readMeasurementTable(file), index; zero_injection_sigma = zero_injection_sigma)
  isempty(imp.sets) && throw(ArgumentError("$(file): no row of the table can be assigned to a measurement set\n$(measurementImportReport(imp))"))
  set = if timestamp === nothing
    first(imp.sets)
  else
    ts = something(_measurement_table_timestamp(timestamp), String(timestamp))
    k = findfirst(s -> s.timestamp == ts, imp.sets)
    k === nothing && throw(ArgumentError("$(file): no measurement set with timestamp $(timestamp); the table has $(join((isempty(s.timestamp) ? "(none)" : s.timestamp for s in imp.sets), ", "))"))
    imp.sets[k]
  end
  applyMeasurementSet!(net, set; replace = replace, zero_injection = zero_injection, zero_injection_sigma = zero_injection_sigma)
  return set
end

"""
    measurementImportReport(imp) -> String

Plain-text report of a `MeasurementTableImport`: per set the row count, the
active and inactive measurements, the zero-injection rows added by
`applyMeasurementSet!` (for a set that was applied) and the count per
problem, followed by the rejected rows. Every problem row is listed with its line.
"""
function measurementImportReport(imp::MeasurementTableImport)::String
  io = IOBuffer()
  println(io, "measurement table import: ", length(imp.sets), " set(s), ", length(imp.rejected), " rejected row(s)")
  for s in imp.sets
    probs = sort!([(k, v) for (k, v) in s.counts if !(k in (:rows, :active, :inactive, :zero_injection)) && v > 0]; by = first)
    # the zero-injection count exists once the set was applied to a net
    zi = haskey(s.counts, :zero_injection) ? string(", zero-injection rows ", s.counts[:zero_injection]) : ""
    println(io, "set ", isempty(s.timestamp) ? "(no timestamp)" : s.timestamp, ": rows ", s.counts[:rows], ", measurements ", length(s.measurements), " (active ", s.counts[:active], ", inactive ", s.counts[:inactive], ")", zi, isempty(probs) ? "" : string(", ", join(("$(k) $(v)" for (k, v) in probs), ", ")))
    for p in s.problems
      println(io, "  line ", p.line, ": ", p.problem, ": ", p.detail)
    end
  end
  for p in imp.rejected
    println(io, "rejected line ", p.line, ": ", p.problem, ": ", p.detail)
  end
  return String(take!(io))
end
