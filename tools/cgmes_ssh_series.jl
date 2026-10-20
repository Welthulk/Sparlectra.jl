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

# file: tools/cgmes_ssh_series.jl
# purpose: synthetic SSH series from one CGMES snapshot. Reads one bundle
#          and writes one bundle per step in which only the SSH changes:
#          loads follow a daily profile, wind units a ramp, thermal units
#          keep the active power balance, switches follow an optional event
#          list. Every other file is written byte for byte as read, no object
#          is added or removed. Equal inputs and seed give identical bytes.
#
# usage: julia --project=. tools/cgmes_ssh_series.jl <base>... <outdir> --roles <roles.toml>
#            [--steps 96] [--minutes 15] [--seed 1] [--events <events.csv>] [--name <name>]
#
#        <base> is what loadCGMES reads: a zip (nested zips included), a
#        folder, an xml, or several of them (base case plus boundary). The
#        last positional argument is the output directory.
#
# Role mapping (TOML). Keys are CGMES classes, values load, wind, thermal or
# fixed:
#
#     EnergyConsumer = "load"
#     SynchronousMachine = "thermal"
#     "SynchronousMachine:WindGeneratingUnit" = "wind"
#
#     [noise]
#     load = 0.02    # relative sigma of the load factor, default 0.02
#     wind = 0.05    # absolute sigma of the wind factor, default 0.05
#
#     [range]
#     load = [0.8, 1.2]     # factor range of the load profile, default [0.6, 1.2]
#     wind = [0.0, 1.3]     # factor range of the wind profile, default [0.0, 1.3]
#     thermal = [0.5, 1.5]  # bounds of the thermal balance factor, default none
#     EnergySource = [0.8, 1.2]  # range of the objects of one role key
#
#   A range is [low, high] with 0 <= low <= high. Without a range key the
#   profile runs as described below, and the output is the same as before the
#   key existed. A key of the role mapping whose role is load or wind can
#   have its own range: its objects follow the profile of their role on that
#   range, all other objects keep the range of the profile. The thermal
#   factor is one factor for all thermal objects, so a thermal key has no
#   range of its own.
#
#   Classes: EnergyConsumer, ConformLoad, NonConformLoad, EnergySource,
#   EquivalentInjection, SynchronousMachine, AsynchronousMachine,
#   ExternalNetworkInjection, PowerElectronicsConnection (CGMES 3.0). A
#   machine can also be keyed by the class of the object its GeneratingUnit
#   reference points at, written "<machine class>:<unit class>" (quoted,
#   because of the colon); that key wins over the plain class. A
#   PowerElectronicsConnection is keyed the same way by the class of its
#   units (BatteryUnit, PhotoVoltaicUnit, PowerElectronicsWindUnit,
#   PowerElectronicsUnit), whose PowerElectronicsConnection reference points
#   at it, for example "PowerElectronicsConnection:PhotoVoltaicUnit". The
#   unit-class key applies only when all units of the connection have that
#   class; a connection with units of several classes takes the plain key and
#   gets a note. Like any key with role load or wind, such a key can have its
#   own range; with role thermal it takes part in the one thermal factor and
#   its bounds. A class without a key is fixed, and so is
#   everything else in the SSH: tap steps, shunt sections, regulation
#   targets (the controllers move them). Nothing is guessed.
#
# Profiles. The time of day is the clock time of the step's scenarioTime as
# written in the SSH; x = (k-1)/(n-1) is the position of step k in the series.
#   load     low + (high - low) * g(h) with the object's range (the range of
#            its key, else the load range), g a sum of three
#            bells (08:00 width 1.5 h weight 0.75, 13:00 width 3 h weight
#            0.5, 19:00 width 2 h weight 1) capped at 1. With the default
#            range 0.6 + 0.6 * g(h): night 0.6, morning peak 1.13, evening
#            peak 1.2. Per object and step times (1 + noise * N(0,1)),
#            clipped to the range.
#   wind     r(x) = ramp 0.2 until x = 0.2, up to 1.1 at x = 0.5, down to 0.6
#            at x = 0.9, mapped onto the object's range (the range of its
#            key, else the wind range) as
#            low + r(x) * (high - low) / 1.3 (the default range gives r(x)
#            itself); per object and step plus noise * N(0,1), clipped to the
#            range.
#   thermal  one factor per step for all thermal objects, chosen so that the
#            sum of active power over all objects of the classes above stays
#            at the base sum. An error names the step when the thermal
#            objects carry no active power. Without range.thermal an error
#            also names the step when the factor would be negative. With
#            range.thermal the factor is clipped to the bounds; in a step
#            where a bound binds the rest of the change stays with the slack,
#            and the manifest notes of that step's bundle give the needed
#            factor and the rest in MW. Without a thermal key there is no
#            balance, the slack takes it.
#   fixed    unchanged.
# The objects of every class above write p in load sign convention (positive
# when power flows out of the node), so one sum covers them all. p and q are
# scaled with the same factor: q keeps its base ratio to p, a constant power
# factor for loads; for a voltage regulating machine the SSH q is only the
# start value the power flow replaces. Values are written with at most six
# decimals; a value whose number does not change keeps its text.
#
# Reference units: CGMES referencePriority 0 means "not a reference", 1 is the
# highest priority. The objects with the smallest positive referencePriority
# (all of them on a tie) keep their SSH values whatever their role: the power
# flow makes them the slack, their output is a result.
#
# Random numbers: one Random.Xoshiro(seed), drawn per step for the load
# objects and then the wind objects, each in mRID order, one draw per object
# even when its p is zero.
#
# Events (CSV, optional, the format step,element_id,action,value of the
# QSTS event list). action is open or close on a switch (Switch and its
# subclasses); the state holds from that step until the next event on the
# same switch, value is not read. A line whose first column is "step" is the
# header and must read step,element_id,action,value; blank lines and lines
# starting with # are skipped; anything else is an error that names the line. The TP is not recomputed: a switch closed this way between
# two topological nodes is a closed, not retained switch the TP does not
# merge, and the series check reports it as SSH against TP.
#
# Output per step k: <name>_<yyyymmddTHHMM>Z.zip (name defaults to the first
# base path without extension, the stamp is the step's scenarioTime in UTC)
# with every XML of the base and the rewritten SSH. In the SSH header
# md:Model.scenarioTime is base time + (k-1) * minutes in the base format,
# md:Model.description gets step, seed, roles, the ranges the mapping sets and
# events appended to the base text (inserted after scenarioTime when the base
# has none), and the rdf:about of md:FullModel is
# a new UUID, the SHA-1 of base SSH id, step and seed in UUID version 5
# layout, behind the prefix of the base (urn:uuid:, #_ or #). Every other
# header field stays as in the base.
#
# manifest.toml in the output directory: a list "bundle", one per step, with
# name, model (the modelingAuthoritySet of the EQ, as modeling_authority),
# scenario_time, modeling_authority, eq_uuid, target (the zip name), files
# (kind, uuid, file, source), missing_dependencies (DependentOn ids no file of
# the bundle supplies) and notes, plus a table "generator" with tool, base,
# seed, steps, step_minutes, roles, noise, events and reference_units. With a
# range key the generator table also has range, and with range.thermal and at
# least one thermal object outside the reference units every bundle has
# thermal_factor and imbalance_mw (MW left to the slack, load sign, 0 where no
# bound binds). A range that reaches no object gets a note in every bundle.
# Every id is written as the importer's header reader gives it: "urn:uuid:x"
# as x, "#_x" as _x.
#
# Every input is checked and every step computed before the output directory
# is created; a failed run writes nothing.

using Sparlectra
using Dates
using Printf
using Random
using SHA
using TOML
using ZipArchives

const _SSHS_C = Sparlectra.CGMESImporter

# the p and q property of every class the role mapping knows; an inherited
# attribute is written with the class that defines it
const _SSHS_PQ = Dict{String,Tuple{String,String}}(
  "EnergyConsumer" => ("EnergyConsumer.p", "EnergyConsumer.q"),
  "ConformLoad" => ("EnergyConsumer.p", "EnergyConsumer.q"),
  "NonConformLoad" => ("EnergyConsumer.p", "EnergyConsumer.q"),
  "EnergySource" => ("EnergySource.activePower", "EnergySource.reactivePower"),
  "EquivalentInjection" => ("EquivalentInjection.p", "EquivalentInjection.q"),
  "SynchronousMachine" => ("RotatingMachine.p", "RotatingMachine.q"),
  "AsynchronousMachine" => ("RotatingMachine.p", "RotatingMachine.q"),
  "ExternalNetworkInjection" => ("ExternalNetworkInjection.p", "ExternalNetworkInjection.q"),
  "PowerElectronicsConnection" => ("PowerElectronicsConnection.p", "PowerElectronicsConnection.q"),
)
const _SSHS_MACHINES = ["SynchronousMachine", "AsynchronousMachine"]
const _SSHS_UNITS = ["GeneratingUnit", "ThermalGeneratingUnit", "HydroGeneratingUnit", "WindGeneratingUnit", "SolarGeneratingUnit", "NuclearGeneratingUnit"]
# the unit classes of a PowerElectronicsConnection, as the importer lists them
const _SSHS_PEC = "PowerElectronicsConnection"
const _SSHS_PE_UNITS = [String(c) for c in _SSHS_C._PE_UNIT_CLASSES]
const _SSHS_ROLES = ["load", "wind", "thermal", "fixed"]
const _SSHS_NOISE_DEFAULT = Dict{String,Float64}("load" => 0.02, "wind" => 0.05)
# factor ranges without a range key; thermal has none, its factor is only
# checked against 0
const _SSHS_RANGE_DEFAULT = Dict{String,Tuple{Float64,Float64}}("load" => (0.6, 1.2), "wind" => (0.0, 1.3))
# the profile names under [range]; every other key there is a role key
const _SSHS_RANGE_KEYS = ["load", "wind", "thermal"]
# the importer's switch classes plus the other subclasses of Switch
const _SSHS_SWITCHES = vcat([String(c) for c in _SSHS_C._SWITCH_CLASSES], ["Fuse", "GroundDisconnector", "Sectionaliser", "Recloser", "DisconnectingCircuitBreaker", "Cut"])
# kind of a file in the manifest: its first profile in this order
const _SSHS_KINDS = [:EQ, :EQ_BD, :TP, :TP_BD, :SSH, :SV, :DL, :GL, :DY, :EQ_OP, :EQ_SC, :EQ_BD_OP]

const _SSHS_HEAD = r"<([A-Za-z0-9_]+):FullModel\s+rdf:about=\"([^\"]*)\"\s*>"
const _SSHS_OBJ = r"<([A-Za-z0-9_]+):([A-Za-z0-9_]+)\s+rdf:(?:about|ID)=\"#?([^\"]*)\"\s*>"
const _SSHS_PROP = r"<([A-Za-z0-9_]+):([A-Za-z0-9_]+\.[A-Za-z0-9_]+)>([^<]*)</\1:\2>"
const _SSHS_TIME = r"^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2})(?::(\d{2})(?:\.(\d+))?)?(Z|[+-]\d{2}:?\d{2})?$"

# one SSH file, scanned once: the byte ranges of the header values and of
# every simple property of every object block
struct _SSHSText
  name::String
  text::String
  pfx::String                     # namespace prefix of md:FullModel
  about::UnitRange{Int}           # rdf:about value of md:FullModel
  id_prefix::String               # "urn:uuid:", "#_", "#" or ""
  base_id::String
  time::UnitRange{Int}            # md:Model.scenarioTime value
  descr::Union{Nothing,UnitRange{Int}}
  descr_at::Int                   # insertion point without a description
  descr_indent::String
  class::Dict{String,String}      # mRID -> class of its block
  values::Dict{Tuple{String,String},Vector{UnitRange{Int}}}  # (mRID, property) -> value ranges
end

# a scenarioTime as written: the instant and the format to write it back
struct _SSHSTime
  t::DateTime
  seconds::Bool
  frac::Int
  zone::String
end

# one object whose p or q a role may change
struct _SSHSUnit
  mrid::String
  class::String
  key::String                     # the role key that matched, else the class
  role::String
  reference::Bool
  text::Int                       # index into the scanned SSH files
  p::Union{Nothing,UnitRange{Int}}
  q::Union{Nothing,UnitRange{Int}}
  p0::Float64
  q0::Float64
end

_sshs_err(msg) = throw(ArgumentError(string("cgmes_ssh_series: ", msg)))

_sshs_xml_escape(s::AbstractString) = replace(String(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")

# roles by key, noise per profile, and the ranges the mapping sets (only
# those: an unset range keeps the default path and the output bytes)
function _sshs_roles(spec)::Tuple{Dict{String,String},Dict{String,Float64},Dict{String,Tuple{Float64,Float64}}}
  d = spec isa AbstractString ? TOML.parsefile(spec) : spec
  roles = Dict{String,String}()
  noise = copy(_SSHS_NOISE_DEFAULT)
  ranges = Dict{String,Tuple{Float64,Float64}}()
  for (k, v) in d
    if k == "noise"
      v isa AbstractDict || _sshs_err("role mapping: noise must be a table with load and wind")
      for (p, x) in v
        haskey(noise, p) || _sshs_err("role mapping: noise.$(p); only load and wind carry noise")
        (x isa Real && x >= 0) || _sshs_err("role mapping: noise.$(p) = $(repr(x)), expected a number >= 0")
        noise[p] = Float64(x)
      end
      continue
    end
    if k == "range"
      v isa AbstractDict || _sshs_err("role mapping: range must be a table with $(join(_SSHS_RANGE_KEYS, ", ")) or keys of the mapping")
      for (p, x) in v
        ok = x isa AbstractVector && length(x) == 2 && all(y -> y isa Real && !(y isa Bool) && isfinite(y), x) && 0 <= x[1] <= x[2]
        ok || _sshs_err("role mapping: range.$(p) = $(repr(x)), expected [low, high] with 0 <= low <= high")
        ranges[p] = (Float64(x[1]), Float64(x[2]))
      end
      continue
    end
    parts = split(k, ':'; limit = 2)
    cls = String(parts[1])
    haskey(_SSHS_PQ, cls) || _sshs_err("role mapping: unknown class $(cls) in key \"$(k)\"; known: $(join(sort!(collect(keys(_SSHS_PQ))), ", "))")
    if length(parts) == 2
      unit_classes = cls in _SSHS_MACHINES ? _SSHS_UNITS : (cls == _SSHS_PEC ? _SSHS_PE_UNITS : String[])
      isempty(unit_classes) && _sshs_err("role mapping: key \"$(k)\": only $(join(_SSHS_MACHINES, ", ")) and $(_SSHS_PEC) can be keyed by a unit class")
      parts[2] in unit_classes || _sshs_err("role mapping: key \"$(k)\": unknown unit class $(parts[2]) of $(cls); known: $(join(unit_classes, ", "))")
    end
    (v isa AbstractString && v in _SSHS_ROLES) || _sshs_err("role mapping: \"$(k)\" = $(repr(v)); roles are $(join(_SSHS_ROLES, ", "))")
    roles[k] = String(v)
  end
  # a range key names a profile or a role key whose objects follow the load
  # or the wind profile; checked after the loop because [range] may come
  # before the role keys
  for p in sort!(collect(keys(ranges)))
    p in _SSHS_RANGE_KEYS && continue
    haskey(roles, p) || _sshs_err("role mapping: range.$(p) is neither a profile ($(join(_SSHS_RANGE_KEYS, ", "))) nor a key of the mapping")
    roles[p] in ("load", "wind") || _sshs_err("role mapping: range.$(p): the key has role $(roles[p]); only keys with role load or wind have a range of their own")
  end
  return roles, noise, ranges
end

_sshs_source(f) = isdir(f.source) ? joinpath(f.source, f.name) : string(f.source, "!", f.name)

# the XML payloads of the base with their headers; an entry name that two
# containers share is kept once when the content is equal and refused when it
# differs, because the output zips are flat
function _sshs_files(paths::Vector{String})
  files = _SSHS_C.collectCGMESFiles(paths)
  keep = Int[]
  infos = Dict{Int,Any}()
  seen = Dict{String,Int}()
  notes = String[]
  for (i, f) in enumerate(files)
    nm = basename(f.name)
    j = get(seen, nm, 0)
    if j != 0
      files[j].content == f.content || _sshs_err("two different files named $(nm) in the base: $(_sshs_source(files[j])) and $(_sshs_source(f))")
      push!(notes, "identical copy of $(nm) at $(_sshs_source(f)) left out")
      continue
    end
    seen[nm] = i
    push!(keep, i)
    infos[i] = _SSHS_C.readCGMESHeader(f)
  end
  return [(file = files[i], info = infos[i], name = basename(files[i].name)) for i in keep], notes
end

function _sshs_scan(name::String, text::String)::_SSHSText
  h = match(_SSHS_HEAD, text)
  h === nothing && _sshs_err("SSH $(name) has no md:FullModel header")
  pfx = h.captures[1]
  about = h.offsets[2]:(h.offsets[2] + ncodeunits(h.captures[2]) - 1)
  about_text = h.captures[2]
  id_prefix = startswith(about_text, "urn:uuid:") ? "urn:uuid:" : (startswith(about_text, "#_") ? "#_" : (startswith(about_text, "#") ? "#" : ""))
  head_end = findnext("</$(pfx):FullModel>", text, h.offset)
  head_end === nothing && _sshs_err("SSH $(name): md:FullModel is not closed")
  hs = h.offset + ncodeunits(h.match)
  he = first(head_end) - 1
  time = nothing
  descr = nothing
  descr_at = 0
  indent = "    "
  for p in eachmatch(_SSHS_PROP, SubString(text, hs, he))
    # offsets of a match in a SubString count from its start
    vstart = hs - 1 + p.offsets[3]
    vr = vstart:(vstart + ncodeunits(p.captures[3]) - 1)
    if p.captures[2] == "Model.scenarioTime"
      time = vr
      el = hs - 1 + p.offset
      descr_at = el + ncodeunits(p.match)
      nl = findprev(==(UInt8('\n')), codeunits(text), el)
      lead = nl === nothing ? "" : text[(nl + 1):(el - 1)]
      all(isspace, lead) && (indent = lead)
    elseif p.captures[2] == "Model.description"
      descr = vr
    end
  end
  time === nothing && _sshs_err("SSH $(name) has no md:Model.scenarioTime")
  class = Dict{String,String}()
  values = Dict{Tuple{String,String},Vector{UnitRange{Int}}}()
  pos = first(head_end)
  while true
    m = match(_SSHS_OBJ, text, pos)
    m === nothing && break
    pfx_o, cls, mrid = m.captures
    bs = m.offset + ncodeunits(m.match)
    stop = findnext("</$(pfx_o):$(cls)>", text, bs)
    stop === nothing && _sshs_err("SSH $(name): block of $(cls) $(mrid) is not closed")
    class[mrid] = cls
    for p in eachmatch(_SSHS_PROP, SubString(text, bs, first(stop) - 1))
      vstart = bs - 1 + p.offsets[3]
      push!(get!(Vector{UnitRange{Int}}, values, (String(mrid), String(p.captures[2]))), vstart:(vstart + ncodeunits(p.captures[3]) - 1))
    end
    pos = last(stop) + 1
  end
  return _SSHSText(name, text, pfx, about, id_prefix, about_text[(ncodeunits(id_prefix) + 1):end], time, descr, descr_at, indent, class, values)
end

function _sshs_parse_time(s::AbstractString, where_::String)::_SSHSTime
  m = match(_SSHS_TIME, strip(s))
  m === nothing && _sshs_err("$(where_): scenarioTime \"$(s)\" is not an ISO 8601 time")
  c = m.captures
  sec = c[6] === nothing ? 0 : parse(Int, c[6])
  frac = c[7] === nothing ? "" : c[7]
  ms = isempty(frac) ? 0 : parse(Int, rpad(frac[1:min(3, end)], 3, '0'))
  t = DateTime(parse(Int, c[1]), parse(Int, c[2]), parse(Int, c[3]), parse(Int, c[4]), parse(Int, c[5]), sec, ms)
  return _SSHSTime(t, c[6] !== nothing, ncodeunits(frac), something(c[8], ""))
end

function _sshs_time_text(f::_SSHSTime, t::DateTime)::String
  s = @sprintf("%04d-%02d-%02dT%02d:%02d", year(t), month(t), day(t), hour(t), minute(t))
  f.seconds && (s *= @sprintf(":%02d", second(t)))
  if f.frac > 0
    ms = @sprintf("%03d", millisecond(t))
    s *= "." * (f.frac <= 3 ? ms[1:f.frac] : ms * "0"^(f.frac - 3))
  end
  return s * f.zone
end

# file name stamp: the step's time in UTC
function _sshs_stamp(f::_SSHSTime, t::DateTime)::String
  u = t
  if ncodeunits(f.zone) > 1
    d = replace(f.zone[2:end], ":" => "")
    off = Minute(parse(Int, d[1:2]) * 60 + parse(Int, d[3:4]))
    u = f.zone[1] == '-' ? t + off : t - off
  end
  return @sprintf("%04d%02d%02dT%02d%02dZ", year(u), month(u), day(u), hour(u), minute(u))
end

function _sshs_uuid(base_id::AbstractString, step::Integer, seed::Integer)::String
  b = sha1(string(base_id, "|ssh-series|", step, "|", seed))[1:16]
  b[7] = (b[7] & 0x0f) | 0x50
  b[9] = (b[9] & 0x3f) | 0x80
  x = bytes2hex(b)
  return string(x[1:8], "-", x[9:12], "-", x[13:16], "-", x[17:20], "-", x[21:32])
end

# the rdf:about of the step's md:FullModel, with the prefix of the base
_sshs_about(t::_SSHSText, step::Integer, seed::Integer) = string(t.id_prefix, _sshs_uuid(t.base_id, step, seed))

# an rdf:about as readCGMESHeader turns it into model_id ("#_x" gives "_x",
# "urn:uuid:x" gives "x"); the manifest and the DependentOn comparison use
# this form, so they match what a reader of the written bundle sees
_sshs_header_id(about::AbstractString) = _SSHS_C._strip_urn(_SSHS_C._strip_mrid(about))

# at most six decimals, no exponent, no trailing zeros
function _sshs_number(v::Float64)::String
  r = round(v; digits = 6)
  r == 0 && return "0"
  s = @sprintf("%.6f", r)
  return String(rstrip(rstrip(s, '0'), '.'))
end

_sshs_bell(h, c, w) = (d = abs(h - c); d = min(d, 24.0 - d); exp(-0.5 * (d / w)^2))
# the profile shapes on their range (lo, hi). The default ranges give the same
# floating point values as the fixed formulas before the range key: 1.2 - 0.6
# is exactly 0.6 and (1.3 - 0.0) / 1.3 exactly 1.0, so the bytes of an
# existing series do not change
_sshs_load_shape(h::Float64, lo::Float64, hi::Float64) = lo + (hi - lo) * min(1.0, 0.75 * _sshs_bell(h, 8.0, 1.5) + 0.5 * _sshs_bell(h, 13.0, 3.0) + _sshs_bell(h, 19.0, 2.0))
_sshs_wind_ramp(x::Float64, lo::Float64, hi::Float64) = lo + (0.2 + 0.9 * clamp((x - 0.2) / 0.3, 0.0, 1.0) - 0.5 * clamp((x - 0.7) / 0.2, 0.0, 1.0)) * ((hi - lo) / 1.3)

function _sshs_value(s::_SSHSText, r::UnitRange{Int}, what::String)::Float64
  v = tryparse(Float64, strip(s.text[r]))
  v === nothing && _sshs_err("SSH $(s.name): $(what) = \"$(s.text[r])\" is not a number")
  return v
end

# the switch events: step -> [(switch mRID, open)], checked against the SSH
function _sshs_events(path::AbstractString, steps::Int, switch_open::Dict{String,Tuple{Int,UnitRange{Int}}})
  out = Dict{Int,Vector{Tuple{String,Bool}}}()
  seen = Set{Tuple{Int,String}}()
  n = 0
  for (i, line) in enumerate(eachline(path))
    s = strip(replace(line, "\ufeff" => ""))
    (isempty(s) || startswith(s, '#')) && continue
    cols = [strip(c) for c in split(s, ',')]
    where_ = "events $(path), line $(i)"
    if lowercase(cols[1]) == "step"
      # the header is strict: the old column name element_mrid is an error
      lowercase.(cols) == ["step", "element_id", "action", "value"] || _sshs_err("$(where_): header \"$(s)\", expected step,element_id,action,value")
      continue
    end
    length(cols) in (3, 4) || _sshs_err("$(where_): expected step,element_id,action,value, got \"$(s)\"")
    step = tryparse(Int, cols[1])
    step === nothing && _sshs_err("$(where_): step \"$(cols[1])\" is not a whole number")
    1 <= step <= steps || _sshs_err("$(where_): step $(step) outside 1..$(steps)")
    mrid = String(lstrip(cols[2], '#'))
    action = lowercase(cols[3])
    action in ("open", "close") || _sshs_err("$(where_): action \"$(cols[3])\"; this generator writes open and close on switches")
    haskey(switch_open, mrid) || _sshs_err("$(where_): $(mrid) is not a switch with Switch.open in the SSH")
    (step, mrid) in seen && _sshs_err("$(where_): a second event for $(mrid) at step $(step)")
    push!(seen, (step, mrid))
    push!(get!(Vector{Tuple{String,Bool}}, out, step), (mrid, action == "open"))
    n += 1
  end
  return out, n
end

function _sshs_render(s::_SSHSText, edits::Vector{Tuple{UnitRange{Int},String}})::Vector{UInt8}
  sort!(edits; by = e -> first(e[1]))
  bytes = codeunits(s.text)
  io = IOBuffer(; sizehint = length(bytes) + 4096)
  pos = 1
  for (r, v) in edits
    write(io, view(bytes, pos:(first(r) - 1)))
    write(io, v)
    pos = last(r) + 1
  end
  write(io, view(bytes, pos:length(bytes)))
  return take!(io)
end

function _sshs_kind(info)::String
  for k in _SSHS_KINDS
    k in info.profiles && return String(k)
  end
  return "UNKNOWN"
end

"""
    generateCGMESSSHSeries(base, outdir; roles, steps = 96, minutes = 15, seed = 1, events = nothing, name = "")

Write `steps` CGMES bundles into `outdir` that differ from `base` only in
the SSH (see the file header for roles, profiles, events and the output
format). `base` is a path or a vector of paths as `loadCGMES` reads them;
`roles` is the path of the TOML role mapping or its parsed `Dict`; `events`
the path of an event CSV or `nothing`. Throws `ArgumentError` for a bad
mapping, a bad event line, a base without SSH and a thermal balance that
has no solution; every step is computed and checked before `outdir` is
created, so a failed run writes nothing.

Returns a named tuple: `bundles` (zip paths in step order), `manifest`,
`scenario_times` (as written), `role` (mRID => role of every object whose p
or q the SSH carries), `reference` (mRIDs kept as reference units),
`base_p` (the base sum of p over those objects), `p_sum` and
`thermal_factor` per step (`NaN` without thermal objects), `imbalance` per
step (MW the slack keeps because `range.thermal` clipped the thermal factor,
load sign, 0 where no bound binds), `switch_open` (per step the switches an
event set, mRID => open) and `notes`.
"""
function generateCGMESSSHSeries(base, outdir::AbstractString; roles, steps::Integer = 96, minutes::Integer = 15, seed::Integer = 1, events = nothing, name::AbstractString = "")
  steps >= 1 || _sshs_err("steps must be at least 1, got $(steps)")
  minutes >= 1 || _sshs_err("minutes must be at least 1, got $(minutes)")
  role_of_key, noise, ranges = _sshs_roles(roles)
  paths = base isa AbstractString ? [String(base)] : String[String(p) for p in base]
  isempty(paths) && _sshs_err("no base path")
  files, notes = _sshs_files(paths)
  ssh_idx = [i for (i, f) in enumerate(files) if :SSH in f.info.profiles]
  isempty(ssh_idx) && _sshs_err("the base has no SSH file")
  texts = [_sshs_scan(files[i].name, files[i].file.content) for i in ssh_idx]
  times = [_sshs_parse_time(t.text[t.time], "SSH $(t.name)") for t in texts]
  # EQ and SSH merged: class, generating unit and referencePriority per object
  store = _SSHS_C.loadCGMES(paths; profile_filter = Set([:EQ, :EQ_OP, :EQ_BD, :EQ_BD_OP, :SSH]))

  # the unit classes of every PowerElectronicsConnection: its units point at
  # it, not the other way round as for a machine and its GeneratingUnit
  pe_unit_classes = Dict{String,Set{String}}()
  for c in _SSHS_C._PE_UNIT_CLASSES, u in _SSHS_C.objectsOf(store, c)
    pec = get(u.refs, :PowerElectronicsConnection, nothing)
    pec === nothing || push!(get!(Set{String}, pe_unit_classes, pec), String(u.class))
  end

  # the objects whose p or q a role may change, and every switch's open flag
  units = _SSHSUnit[]
  mixed_pec = String[]
  used_keys = Set{String}()
  switch_open = Dict{String,Tuple{Int,UnitRange{Int}}}()
  for (ti, t) in enumerate(texts)
    for (mrid, ssh_cls) in t.class
      o = get(store.objects, mrid, nothing)
      cls = o === nothing ? ssh_cls : String(o.class)
      if cls in _SSHS_SWITCHES || ssh_cls in _SSHS_SWITCHES
        r = get(t.values, (mrid, "Switch.open"), UnitRange{Int}[])
        isempty(r) || (switch_open[mrid] = (ti, first(r)))
        continue
      end
      haskey(_SSHS_PQ, cls) || continue
      pprop, qprop = _SSHS_PQ[cls]
      pr = get(t.values, (mrid, pprop), UnitRange{Int}[])
      qr = get(t.values, (mrid, qprop), UnitRange{Int}[])
      (isempty(pr) && isempty(qr)) && continue
      key = cls
      if cls in _SSHS_MACHINES && o !== nothing
        gu = _SSHS_C.ref(store, o, :GeneratingUnit)
        gu === nothing || (haskey(role_of_key, "$(cls):$(gu.class)") && (key = "$(cls):$(gu.class)"))
      elseif cls == _SSHS_PEC
        uc = sort!(collect(get(pe_unit_classes, mrid, Set{String}())))
        if length(uc) == 1
          haskey(role_of_key, "$(cls):$(uc[1])") && (key = "$(cls):$(uc[1])")
        elseif length(uc) > 1 && any(c -> haskey(role_of_key, "$(cls):$(c)"), uc)
          # nothing is guessed: units of several classes give no one unit key
          push!(mixed_pec, "$(cls) $(mrid) has units of several classes ($(join(uc, ", "))), the plain key $(cls) applies")
        end
      end
      p = isempty(pr) ? nothing : first(pr)
      q = isempty(qr) ? nothing : first(qr)
      p0 = p === nothing ? 0.0 : _sshs_value(t, p, "$(pprop) of $(mrid)")
      q0 = q === nothing ? 0.0 : _sshs_value(t, q, "$(qprop) of $(mrid)")
      haskey(role_of_key, key) && push!(used_keys, key)
      push!(units, _SSHSUnit(mrid, cls, key, get(role_of_key, key, "fixed"), false, ti, p, q, p0, q0))
    end
  end
  append!(notes, sort!(mixed_pec))
  # a key that matches nothing gives a flat series for that class; say so
  for k in sort!(collect(keys(role_of_key)))
    k in used_keys || push!(notes, "role key $(k) matches no object with p or q in the SSH")
  end
  sort!(units; by = u -> u.mrid)
  prio(u) = (o = get(store.objects, u.mrid, nothing); o === nothing ? 0.0 : _SSHS_C.num(o, :referencePriority, 0.0))
  best = minimum((prio(u) for u in units if prio(u) > 0); init = Inf)
  # in place: units is captured by the closures below and stays one binding
  for (i, u) in enumerate(units)
    units[i] = _SSHSUnit(u.mrid, u.class, u.key, u.role, prio(u) == best, u.text, u.p, u.q, u.p0, u.q0)
  end
  reference = [u.mrid for u in units if u.reference]
  ev, n_events = events === nothing ? (Dict{Int,Vector{Tuple{String,Bool}}}(), 0) : _sshs_events(events, Int(steps), switch_open)

  loads = [i for (i, u) in enumerate(units) if u.role == "load" && !u.reference]
  winds = [i for (i, u) in enumerate(units) if u.role == "wind" && !u.reference]
  thermals = [i for (i, u) in enumerate(units) if u.role == "thermal" && !u.reference]
  # a range that reaches no object has no effect; say so, as for a role key
  # that matches nothing. A key range reaches the objects of its key, a
  # profile range the other objects of its role, range.thermal the thermal
  # objects; reference units follow no profile
  reached = Set{String}(haskey(ranges, units[i].key) ? units[i].key : units[i].role for i in vcat(loads, winds))
  isempty(thermals) || push!(reached, "thermal")
  for k in sort!(collect(keys(ranges)))
    k in reached || push!(notes, "range.$(k) reaches no object outside the reference units and has no effect")
  end
  base_p = sum(u.p0 for u in units if u.p !== nothing; init = 0.0)
  p_thermal = sum(units[i].p0 for i in thermals if units[i].p !== nothing; init = 0.0)
  role_text = isempty(role_of_key) ? "none" : join(("$(k)=$(v)" for (k, v) in sort!(collect(role_of_key))), ", ")
  # only the ranges the mapping sets: without one the description stays as it was
  range_text = isempty(ranges) ? "" : string(", ranges ", join(("$(k)=$(_sshs_number(lo))..$(_sshs_number(hi))" for (k, (lo, hi)) in sort!(collect(ranges))), ", "))
  # the base SSH ids as the header reader gives them, by position in texts
  base_ids = [files[i].info.model_id for i in ssh_idx]

  # pass 1, before anything is written: p and q of every unit and step, so a
  # balance without a solution stops the run before outdir exists
  times_k = [[tm.t + Minute((k - 1) * minutes) for tm in times] for k in 1:steps]
  P = [u.p0 for u in units, _ in 1:steps]
  Q = [u.q0 for u in units, _ in 1:steps]
  th_factor = fill(NaN, steps)
  # per load and wind object the range of its key, else that of its profile
  load_range = [get(ranges, units[i].key, get(ranges, "load", _SSHS_RANGE_DEFAULT["load"])) for i in loads]
  wind_range = [get(ranges, units[i].key, get(ranges, "wind", _SSHS_RANGE_DEFAULT["wind"])) for i in winds]
  bounds = get(ranges, "thermal", nothing)
  # per step the change of p the thermal objects do not take because a bound
  # binds (MW, load sign: positive is more consumption, the slack supplies
  # it), and the factor the balance would have needed; 0 and NaN otherwise
  imbalance = zeros(steps)
  th_needed = fill(NaN, steps)
  rng = Random.Xoshiro(seed)
  for k in 1:steps
    tk = times_k[k][1]
    clock = hour(tk) + minute(tk) / 60 + second(tk) / 3600
    x = steps == 1 ? 0.0 : (k - 1) / (steps - 1)
    scale!(i, f) = begin
      units[i].p === nothing || (P[i, k] = round(units[i].p0 * f; digits = 6))
      units[i].q === nothing || (Q[i, k] = round(units[i].q0 * f; digits = 6))
    end
    # one draw per object in mRID order whatever its range, so a key range
    # leaves the values of every other object as they were
    for (j, i) in enumerate(loads)
      lo, hi = load_range[j]
      scale!(i, clamp(_sshs_load_shape(clock, lo, hi) * (1 + noise["load"] * randn(rng)), lo, hi))
    end
    for (j, i) in enumerate(winds)
      lo, hi = wind_range[j]
      scale!(i, clamp(_sshs_wind_ramp(x, lo, hi) + noise["wind"] * randn(rng), lo, hi))
    end
    isempty(thermals) && continue
    delta = sum(P[i, k] - units[i].p0 for i in vcat(loads, winds) if units[i].p !== nothing; init = 0.0)
    p_thermal == 0.0 && delta != 0.0 && _sshs_err("step $(k): the thermal objects carry no active power in the base, no factor keeps the balance")
    fth = p_thermal == 0.0 ? 1.0 : 1.0 - delta / p_thermal
    if bounds === nothing
      fth < 0.0 && _sshs_err("step $(k): the balance needs a thermal factor of $(round(fth; digits = 3)), below 0")
    elseif !(bounds[1] <= fth <= bounds[2])
      th_needed[k] = fth
      fth = clamp(fth, bounds[1], bounds[2])
      # the thermal objects move by (fth - 1) * p_thermal, the rest of the
      # change stays with the slack
      imbalance[k] = delta + (fth - 1.0) * p_thermal
    end
    for i in thermals
      scale!(i, fth)
    end
    th_factor[k] = fth
  end

  # the files that stay as they are, compressed once; every step appends its SSH
  ssh_names = Set(files[i].name for i in ssh_idx)
  template = let io = IOBuffer()
    ZipWriter(io) do w
      for f in files
        f.name in ssh_names && continue
        zip_newfile(w, f.name; compress = true)
        write(w, f.file.content)
      end
    end
    take!(io)
  end

  eq = findfirst(f -> :EQ in f.info.profiles, files)
  authority = eq === nothing ? "" : files[eq].info.modeling_authority
  eq_uuid = eq === nothing ? "" : files[eq].info.model_id
  bname = isempty(name) ? splitext(basename(rstrip(paths[1], ['/', '\\'])))[1] : String(name)
  # pass 2: every check has passed, write the bundles
  mkpath(outdir)
  state = Dict{String,Bool}()           # switch -> open, once an event set it
  switch_base = Dict{String,Bool}(m => lowercase(strip(texts[ti].text[r])) == "true" for (m, (ti, r)) in switch_open)
  bundles = String[]
  stimes = String[]
  p_sum = Float64[]
  switch_steps = Dict{String,Bool}[]
  manifest_bundles = Dict{String,Any}[]
  for k in 1:steps
    tk = times_k[k]
    for (m, open) in get(ev, k, Tuple{String,Bool}[])
      state[m] = open
    end
    # edits per SSH file: header, values, switches
    edits = [Tuple{UnitRange{Int},String}[] for _ in texts]
    for (ti, t) in enumerate(texts)
      push!(edits[ti], (t.about, _sshs_about(t, k, seed)))
      push!(edits[ti], (t.time, _sshs_time_text(times[ti], tk[ti])))
      descr = _sshs_xml_escape("Synthetic SSH step $(k) of $(steps) from $(base_ids[ti]) (tools/cgmes_ssh_series.jl): step $(minutes) min, seed $(seed), roles $(role_text)$(range_text), $(n_events) switch event(s)")
      if t.descr === nothing
        push!(edits[ti], (t.descr_at:(t.descr_at - 1), string("\n", t.descr_indent, "<", t.pfx, ":Model.description>", descr, "</", t.pfx, ":Model.description>")))
      else
        # appended, not replaced: a base description may carry a licence or
        # ownership notice that the derived file has to keep
        push!(edits[ti], (t.descr, string(t.text[t.descr], " | ", descr)))
      end
    end
    total = 0.0
    for (i, u) in enumerate(units)
      if u.p !== nothing
        v = P[i, k]
        v == u.p0 || push!(edits[u.text], (u.p, _sshs_number(v)))
        total += v
      end
      if u.q !== nothing
        v = Q[i, k]
        v == u.q0 || push!(edits[u.text], (u.q, _sshs_number(v)))
      end
    end
    for (m, open) in state
      open == switch_base[m] && continue
      ti, r = switch_open[m]
      push!(edits[ti], (r, open ? "true" : "false"))
    end
    stamp = _sshs_stamp(times[1], tk[1])
    target = string(bname, "_", stamp, ".zip")
    path = joinpath(outdir, target)
    write(path, template)
    zip_append_archive(path) do w
      for (ti, t) in enumerate(texts)
        zip_newfile(w, t.name; compress = true)
        write(w, _sshs_render(texts[ti], edits[ti]))
      end
    end
    push!(bundles, path)
    push!(stimes, _sshs_time_text(times[1], tk[1]))
    push!(p_sum, total)
    push!(switch_steps, copy(state))

    # manifest entry: the SSH carries its new id, the others their own; all
    # ids in the form of the header reader (see _sshs_header_id)
    new_ids = Dict(i => _sshs_header_id(_sshs_about(texts[j], k, seed)) for (j, i) in enumerate(ssh_idx))
    entries = Dict{String,Any}[]
    supplied = Set{String}()
    for (i, f) in enumerate(files)
      id = get(new_ids, i, f.info.model_id)
      push!(supplied, id)
      push!(entries, Dict{String,Any}("kind" => _sshs_kind(f.info), "uuid" => id, "file" => f.name, "source" => _sshs_source(f.file)))
    end
    deps = sort!(unique(String[d for f in files for d in f.info.dependent_on if !(d in supplied)]))
    bnotes = copy(notes)
    if imbalance[k] != 0.0
      push!(bnotes, @sprintf("step %d: the thermal balance needs a factor of %.4f, clipped to %s; %.3f MW of the change stay with the slack", k, th_needed[k], _sshs_number(th_factor[k]), imbalance[k]))
    end
    for (ti, t) in enumerate(texts)
      push!(bnotes, "SSH $(t.name) generated from $(base_ids[ti]): step $(k) of $(steps), seed $(seed)")
      for f in files
        (base_ids[ti] in f.info.dependent_on && !(:SSH in f.info.profiles)) && push!(bnotes, "$(f.name) depends on the base SSH $(base_ids[ti]), which this bundle replaces; it is kept as read")
      end
    end
    entry = Dict{String,Any}("name" => splitext(target)[1], "model" => authority, "scenario_time" => stimes[end], "modeling_authority" => authority,
                             "eq_uuid" => eq_uuid, "target" => target, "files" => entries, "missing_dependencies" => deps, "notes" => bnotes)
    # with thermal bounds every step says what its thermal objects took and
    # what stays with the slack; without them the entry stays as it was
    if bounds !== nothing && !isempty(thermals)
      entry["thermal_factor"] = th_factor[k]
      entry["imbalance_mw"] = imbalance[k]
    end
    push!(manifest_bundles, entry)
  end
  generator = Dict{String,Any}("tool" => "tools/cgmes_ssh_series.jl", "base" => paths, "seed" => Int(seed), "steps" => Int(steps), "step_minutes" => Int(minutes),
                               "roles" => role_of_key, "noise" => noise, "events" => events === nothing ? "" : String(events), "reference_units" => reference)
  isempty(ranges) || (generator["range"] = Dict{String,Any}(k => [lo, hi] for (k, (lo, hi)) in ranges))
  manifest = joinpath(outdir, "manifest.toml")
  open(manifest, "w") do io
    TOML.print(io, Dict{String,Any}("generator" => generator, "bundle" => manifest_bundles); sorted = true)
  end
  return (bundles = bundles, manifest = manifest, scenario_times = stimes, role = Dict(u.mrid => u.role for u in units), reference = reference,
          base_p = base_p, p_sum = p_sum, thermal_factor = th_factor, imbalance = imbalance, switch_open = switch_steps, notes = notes)
end

function _sshs_main(args::Vector{String})
  usage = "usage: julia --project=. tools/cgmes_ssh_series.jl <base>... <outdir> --roles <roles.toml> [--steps 96] [--minutes 15] [--seed 1] [--events <events.csv>] [--name <name>]"
  pos = String[]
  opts = Dict{String,String}()
  i = 1
  while i <= length(args)
    a = args[i]
    if startswith(a, "--")
      key = a[3:end]
      key in ("roles", "steps", "minutes", "seed", "events", "name") || _sshs_err("unknown option $(a)\n$(usage)")
      i < length(args) || _sshs_err("option $(a) needs a value\n$(usage)")
      opts[key] = args[i + 1]
      i += 2
    else
      push!(pos, a)
      i += 1
    end
  end
  (length(pos) >= 2 && haskey(opts, "roles")) || _sshs_err(usage)
  int_opt(k, d) = haskey(opts, k) ? parse(Int, opts[k]) : d
  r = generateCGMESSSHSeries(pos[1:(end - 1)], pos[end]; roles = opts["roles"], steps = int_opt("steps", 96), minutes = int_opt("minutes", 15),
                             seed = int_opt("seed", 1), events = get(opts, "events", nothing), name = get(opts, "name", ""))
  println("wrote ", length(r.bundles), " bundle(s) and ", r.manifest, "; ", count(v -> v != "fixed", values(r.role)), " of ", length(r.role), " objects with a profile, reference units: ", join(r.reference, ", "))
  clipped = findall(!=(0.0), r.imbalance)
  if !isempty(clipped)
    println("thermal factor clipped to range.thermal in ", length(clipped), " step(s), at most ", @sprintf("%.3f", maximum(abs, r.imbalance)), " MW stay with the slack:")
    for k in clipped
      println(@sprintf("  step %d: factor %s, %.3f MW", k, _sshs_number(r.thermal_factor[k]), r.imbalance[k]))
    end
  end
  return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
  Base.invokelatest(_sshs_main, ARGS)
end
