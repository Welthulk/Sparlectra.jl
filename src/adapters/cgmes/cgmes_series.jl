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
# file: src/adapters/cgmes/cgmes_series.jl
# purpose: series check of CGMES snapshot deliveries: group the snapshots by
#          model, order them by time, and measure per consecutive pair how
#          much of the structure and of the injections they share. An mRID
#          that an export assigns anew is matched by a fingerprint (class,
#          name, connectivity), so a renumbering does not split a series.

# classes whose shared share decides whether two snapshots belong to one
# series (TopologicalNode is not among them: its
# mRIDs are assigned per export and switching belongs to the delta)
const _SERIES_STRUCTURAL_CLASSES = [:ACLineSegment, :PowerTransformer, :PowerTransformerEnd, :BusbarSection, :Breaker, :Disconnector, :ConnectivityNode]
# classes whose turnover is reported as injection churn; the CGMES 3.0
# PowerElectronicsConnection with its units stands next to the machine with
# its generating units (the reader keeps every subclass under its own class)
const _SERIES_INJECTION_CLASSES = [:SynchronousMachine, :EnergyConsumer, :EquivalentInjection, :EnergySource, :GeneratingUnit, :ThermalGeneratingUnit, :HydroGeneratingUnit, :WindGeneratingUnit, :SolarGeneratingUnit, :NuclearGeneratingUnit,
                                   :PowerElectronicsConnection, :PowerElectronicsUnit, :PhotoVoltaicUnit, :BatteryUnit, :PowerElectronicsWindUnit]
# reported for information only
const _SERIES_INFO_CLASSES = [:TopologicalNode, :ConformLoad, :NonConformLoad]

# attributes the parameter drift skips. The reader keeps no profile per
# attribute and a CGMESStore source holds EQ, SSH and SV values in one
# object, so the operating state is excluded by name; SV values sit in Sv*
# classes, which are not compared. Boolean and enum entries never count
# anyway (not numeric) and are listed for completeness. A key kept under its
# dotted name (suffix collision in the reader) is checked by its suffix.
const _SERIES_NON_PARAMETER_ATTRS = Set([
  # SSH and EQ switch and service states
  :open, :locked, :normalOpen, :inService, :normallyInService, :connected,
  # SSH injections and set points
  :p, :q, :activePower, :reactivePower, :voltageMagnitude, :voltageAngle, :targetValue, :targetDeadband,
  :referencePriority, :regulationStatus, :regulationTarget, :controlEnabled, :operatingMode, :normalPF,
  # EQ value that CIM defines as a stored power flow result
  # (GeneratingUnit.initialP); it follows the dispatch, SmallGrid writes
  # -SSH p of each machine there
  :initialP,
  # SSH positions of tap changers and shunts
  :step, :sections,
  # SSH state of a CGMES 3.0 BatteryUnit; the stored energy follows the
  # dispatch from step to step, ratedE stays an EQ parameter
  :storedE, :batteryState,
  # identity strings
  :mRID, :name, :shortName, :description, :energyIdentCodeEic,
])

# relative tolerance of the parameter drift: an exporter that prints in single
# precision changes the 8th digit, a real parameter change is far larger
const _SERIES_DRIFT_RTOL = 1e-6

function _seriesIsParameter(k::Symbol)::Bool
  k in _SERIES_NON_PARAMETER_ATTRS && return false
  s = String(k)
  return !(occursin('.', s) && Symbol(last(split(s, '.'))) in _SERIES_NON_PARAMETER_ATTRS)
end

"""
    CGMESSeriesInput(source; model_key = "", scenario_time = nothing, label = "")

One snapshot for [`checkCGMESSeries`](@ref): `source` is anything
`loadCGMES` reads (folder, zip, xml, or a vector of them) or a `CGMESStore`
that `loadCGMES` returned; a store is used as it is, so snapshots changed in
memory can be checked without writing them. `model_key` groups the
snapshots into models; empty takes the `modelingAuthoritySet` of the EQ
header. `scenario_time` orders them; `nothing` takes the `scenarioTime` of
the SSH header, else of the EQ header. A store keeps the headers of its
files, so these defaults apply to it as well. `label` names the snapshot in
the report (default: the file name of a path source, `"snapshot"` for a
store).
"""
struct CGMESSeriesInput
  source::Any
  model_key::String
  scenario_time::Union{Nothing,DateTime}
  label::String
end

CGMESSeriesInput(source; model_key::AbstractString = "", scenario_time::Union{Nothing,DateTime} = nothing, label::AbstractString = "") =
  CGMESSeriesInput(source, String(model_key), scenario_time, isempty(label) ? (source isa AbstractString ? basename(String(source)) : "snapshot") : String(label))

"""
    CGMESSeriesClassDiff

Comparison of one class between two snapshots: objects only in the first,
only in the second, common by mRID, `renamed` (different mRIDs, one
fingerprint on both sides, counted as common) and `ambiguous` (a
fingerprint shared by several objects on one side, left unmatched). `kind`
is `:structural`, `:injection` or `:info`.

`drifted` counts the common and renamed objects whose equipment parameters
differ between the two snapshots: a numeric attribute present on both sides
whose values differ by more than a relative tolerance of 1e-6. Operating
state (switch states, SSH injections and set points, tap and section
positions, the stored energy and state of a `BatteryUnit`, the stored
dispatch `GeneratingUnit.initialP` of the EQ) and
identity strings do not count, nor does an attribute present
on one side only. Drift is reported only; it does not change the overlap or
the churn and does not split a series.
"""
struct CGMESSeriesClassDiff
  class::Symbol
  kind::Symbol
  only_a::Int
  only_b::Int
  common::Int
  renamed::Int
  ambiguous::Int
  drifted::Int
end

"""
    CGMESSeriesPair

Two consecutive snapshots of one model (`a`, `b` index the report's
snapshots): the time step in minutes (`nothing` without times), the
structural overlap (common over union, renamed objects counted as common),
the injection churn (1 minus that share over the injection classes), the
class comparison, why a new series starts here (`:structural_overlap`,
`:max_gap` or `nothing`) and whether the churn exceeds the configured
limit.
"""
struct CGMESSeriesPair
  a::Int
  b::Int
  gap_minutes::Union{Nothing,Float64}
  structural_overlap::Float64
  injection_churn::Float64
  classes::Vector{CGMESSeriesClassDiff}
  split_reason::Union{Nothing,Symbol}
  churn_flag::Bool
end

"""
    CGMESSeriesReport

Result of [`checkCGMESSeries`](@ref): `snapshots` (label, model key, time,
object counts per class, whether a TP profile was read, and
`ssh_tp_conflicts`: switches closed in the SSH and not retained whose
terminals the TP puts on two topological nodes), `pairs` (one per
consecutive pair inside a model), `series` (per model the runs of snapshot
indices that form a series) and the `config` used.
"""
struct CGMESSeriesReport
  snapshots::Vector{NamedTuple{(:label, :model_key, :scenario_time, :counts, :has_tp, :ssh_tp_conflicts),Tuple{String,String,Union{Nothing,DateTime},Dict{Symbol,Int},Bool,Int}}}
  pairs::Vector{CGMESSeriesPair}
  series::Vector{NamedTuple{(:model_key, :snapshots),Tuple{String,Vector{Int}}}}
  config::Sparlectra.QstsSeriesConfig
end

# "2023-05-24T15:30:00Z", "2023-05-24T15:30:00.000Z", "2023-05-24T15:30Z"
function _cgmesParseTime(s::AbstractString)::Union{Nothing,DateTime}
  m = match(r"^(\d{4}-\d{2}-\d{2}T\d{2}:\d{2})(:\d{2})?", strip(s))
  m === nothing && return nothing
  return DateTime(m.captures[1] * something(m.captures[2], ":00"))
end

# what the check keeps of one snapshot: mRID sets per class (canonical, no
# boundary objects) and a fingerprint per object of the matched classes
struct _SeriesSnapshot
  label::String
  model_key::String
  scenario_time::Union{Nothing,DateTime}
  has_tp::Bool
  ids::Dict{Symbol,Set{String}}
  fingerprint::Dict{Symbol,Dict{String,String}}
  # ConnectivityNode mRID per terminal-bearing object, for the equipment
  # fingerprint after the CN matching
  cns_of::Dict{String,Vector{String}}
  # TN -> its CN mRIDs (TP), for the TN content fingerprint
  tn_cns::Dict{String,Vector{String}}
  # class -> mRID -> numeric equipment parameters, for the parameter drift;
  # an object without any is absent
  params::Dict{Symbol,Dict{String,Dict{Symbol,Float64}}}
  # switches closed in the SSH, not retained, whose terminals the TP puts on
  # two TNs: SSH and TP disagree (case 3 of the stage 3 probe)
  ssh_tp_conflicts::Int
end

# Case 3 count of one snapshot: a switch the importer would close (its
# open rule) and that is not retained, with its two terminals on two TNs.
# The importer turns such a switch into a link, so a solve stays valid; the
# count says that TP and SSH come from different switch states (an exporter
# that did not recompute the TP, or a synthetic series that writes switch
# events into the SSH only).
function _seriesSwitchConflicts(store::CGMESStore)::Int
  any(f -> :TP in f.profiles && !f.skipped, store.files) || return 0
  cn_tn = Dict{String,String}()
  for cn in objectsOf(store, :ConnectivityNode)
    tn = get(cn.refs, :TopologicalNode, nothing)
    tn === nothing || (cn_tn[cn.mrid] = tn)
  end
  # the terminal own TN first, as the importer reads it
  terms = Dict{String,Vector{Tuple{Float64,String}}}()
  for t in objectsOf(store, :Terminal)
    eq = get(t.refs, :ConductingEquipment, nothing)
    eq === nothing && continue
    tn = something(get(t.refs, :TopologicalNode, nothing), get(cn_tn, get(t.refs, :ConnectivityNode, ""), nothing), Some(nothing))
    tn === nothing && continue
    push!(get!(Vector{Tuple{Float64,String}}, terms, eq), (something(num(t, :sequenceNumber), 1.0), tn))
  end
  ctx = _MapCtx()
  n = 0
  for cls in _SWITCH_CLASSES, sw in objectsOf(store, cls)
    v = get(terms, sw.mrid, nothing)
    (v === nothing || length(v) < 2) && continue
    sort!(v)
    v[1][2] == v[2][2] && continue
    _switchIsOpen(ctx, sw) && continue
    something(boolval(sw, :retained), false) && continue
    n += 1
  end
  return n
end

function _seriesSnapshot(input::CGMESSeriesInput)::_SeriesSnapshot
  # a store is taken as it is (its file headers carry scenario time and
  # authority set like those of a path); a path is read with EQ, TP and SSH
  # (the SSH only for the switch states of the case 3 count)
  store = input.source isa CGMESStore ? input.source : redirect_stdout(() -> loadCGMES(input.source; profile_filter = Set([:EQ, :EQ_OP, :EQ_BD, :EQ_BD_OP, :TP, :TP_BD, :SSH])), devnull)
  key = cgmesCanonicalMrid
  eq = [f for f in store.files if :EQ in f.profiles]
  ssh = [f for f in store.files if :SSH in f.profiles]
  model_key = isempty(input.model_key) ? (isempty(eq) ? "" : first(eq).modeling_authority) : input.model_key
  t = input.scenario_time
  if t === nothing
    for f in vcat(ssh, eq)
      t = _cgmesParseTime(f.scenario_time)
      t === nothing || break
    end
  end
  has_tp = any(f -> :TP in f.profiles && !f.skipped, store.files)
  classes = vcat(_SERIES_STRUCTURAL_CLASSES, _SERIES_INJECTION_CLASSES, _SERIES_INFO_CLASSES)
  ids = Dict{Symbol,Set{String}}()
  for c in classes
    ids[c] = Set(key(m) for m in get(store.byclass, c, String[]) if !(m in store.boundary))
  end
  # terminals: equipment -> its CNs (sequence order), CN -> TN (TP)
  cns_of = Dict{String,Vector{Tuple{Float64,String}}}()
  tn_cns = Dict{String,Set{String}}()
  term_cn = Dict{String,String}()
  for tm in objectsOf(store, :Terminal)
    cn = get(tm.refs, :ConnectivityNode, nothing)
    eqm = get(tm.refs, :ConductingEquipment, nothing)
    tn = get(tm.refs, :TopologicalNode, nothing)
    cn === nothing && continue
    term_cn[tm.mrid] = cn
    eqm === nothing || push!(get!(Vector{Tuple{Float64,String}}, cns_of, key(eqm)), (something(num(tm, :sequenceNumber), 1.0), key(cn)))
    tn === nothing || push!(get!(Set{String}, tn_cns, key(tn)), key(cn))
  end
  for cn in objectsOf(store, :ConnectivityNode)
    tn = get(cn.refs, :TopologicalNode, nothing)
    tn === nothing || push!(get!(Set{String}, tn_cns, key(tn)), key(cn.mrid))
  end
  # a transformer end's CN through its terminal
  for e in objectsOf(store, :PowerTransformerEnd)
    tm = get(e.refs, :Terminal, nothing)
    cn = tm === nothing ? nothing : get(term_cn, tm, nothing)
    cn === nothing || push!(get!(Vector{Tuple{Float64,String}}, cns_of, key(e.mrid)), (something(num(e, :endNumber), 1.0), key(cn)))
  end
  fp = Dict{Symbol,Dict{String,String}}()
  params = Dict{Symbol,Dict{String,Dict{Symbol,Float64}}}()
  for c in classes
    d = Dict{String,String}()
    pd = Dict{String,Dict{Symbol,Float64}}()
    for m in get(store.byclass, c, String[])
      m in store.boundary && continue
      o = store.objects[m]
      name = something(str(o, :name), "")
      # CNs are matched by name and container; everything else by name
      # here, the connectivity part is added per pair after the CN matching
      d[key(m)] = c === :ConnectivityNode ? string(name, "|", key(something(get(o.refs, :ConnectivityNodeContainer, nothing), ""))) : name
      pv = Dict{Symbol,Float64}()
      for (k, v) in o.attrs
        _seriesIsParameter(k) || continue
        x = tryparse(Float64, v)
        x === nothing || (pv[k] = x)
      end
      isempty(pv) || (pd[key(m)] = pv)
    end
    fp[c] = d
    params[c] = pd
  end
  cns_sorted = Dict(k => [cn for (_, cn) in sort(v)] for (k, v) in cns_of)
  return _SeriesSnapshot(input.label, model_key, t, has_tp, ids, fp, cns_sorted, Dict(k => sort!(collect(v)) for (k, v) in tn_cns), params, _seriesSwitchConflicts(store))
end

# true when an attribute present on both sides differs beyond the tolerance;
# an attribute on one side only is a profile difference, not a drift
function _seriesDrifted(pa::Dict{String,Dict{Symbol,Float64}}, pb::Dict{String,Dict{Symbol,Float64}}, ma::String, mb::String)::Bool
  va = get(pa, ma, nothing)
  vb = get(pb, mb, nothing)
  (va === nothing || vb === nothing) && return false
  for (k, x) in va
    y = get(vb, k, nothing)
    y === nothing && continue
    (isequal(x, y) || isapprox(x, y; rtol = _SERIES_DRIFT_RTOL)) || return true
  end
  return false
end

# match the objects only in a and only in b by fingerprint; a fingerprint
# that is not unique on its side matches nothing (counted as ambiguous)
function _seriesMatch(only_a::Vector{String}, only_b::Vector{String}, fp_a::Function, fp_b::Function)
  by_a = Dict{String,Vector{String}}()
  by_b = Dict{String,Vector{String}}()
  for m in only_a
    push!(get!(Vector{String}, by_a, fp_a(m)), m)
  end
  for m in only_b
    push!(get!(Vector{String}, by_b, fp_b(m)), m)
  end
  pairs = Dict{String,String}()
  ambiguous = 0
  for (f, ms) in by_a
    other = get(by_b, f, String[])
    isempty(other) && continue
    if length(ms) == 1 && length(other) == 1
      pairs[ms[1]] = other[1]
    else
      ambiguous += length(ms)
    end
  end
  return pairs, ambiguous
end

function _seriesPair(sa::_SeriesSnapshot, sb::_SeriesSnapshot, ia::Int, ib::Int, cfg::Sparlectra.QstsSeriesConfig)::CGMESSeriesPair
  diffs = CGMESSeriesClassDiff[]
  # 1. ConnectivityNodes by name and container; the alias maps b's CN onto a's
  cn_alias = Dict{String,String}()
  function diff!(c::Symbol, kind::Symbol, fa::Function, fb::Function)
    A = sa.ids[c]
    B = sb.ids[c]
    both = intersect(A, B)
    oa = sort!(collect(setdiff(A, B)))
    ob = sort!(collect(setdiff(B, A)))
    matched, amb = _seriesMatch(oa, ob, fa, fb)
    # renamed objects count as common, so their parameters are compared too
    pa = sa.params[c]
    pb = sb.params[c]
    drifted = count(m -> _seriesDrifted(pa, pb, m, m), both) + count(((ma, mb),) -> _seriesDrifted(pa, pb, ma, mb), matched)
    push!(diffs, CGMESSeriesClassDiff(c, kind, length(oa) - length(matched), length(ob) - length(matched), length(both), length(matched), amb, drifted))
    return matched
  end
  # matched pairs come back as a-mRID => b-mRID
  for (am, bm) in diff!(:ConnectivityNode, :structural, m -> sa.fingerprint[:ConnectivityNode][m], m -> sb.fingerprint[:ConnectivityNode][m])
    cn_alias[bm] = am
  end
  # a CN of b written as the matching CN of a
  cnb(cn) = get(cn_alias, cn, cn)
  # 2. equipment: name plus the (matched) CNs of its terminals
  for c in vcat(filter(!=(:ConnectivityNode), _SERIES_STRUCTURAL_CLASSES), _SERIES_INJECTION_CLASSES, filter(!=(:TopologicalNode), _SERIES_INFO_CLASSES))
    kind = c in _SERIES_STRUCTURAL_CLASSES ? :structural : (c in _SERIES_INJECTION_CLASSES ? :injection : :info)
    fa = m -> string(sa.fingerprint[c][m], "|", join(get(sa.cns_of, m, String[]), ","))
    fb = m -> string(sb.fingerprint[c][m], "|", join(map(cnb, get(sb.cns_of, m, String[])), ","))
    diff!(c, kind, fa, fb)
  end
  # 3. TopologicalNodes by content (their CN set), only when both carry a TP
  if sa.has_tp && sb.has_tp
    diff!(:TopologicalNode, :info, m -> join(get(sa.tn_cns, m, String[]), ","), m -> join(sort!(map(cnb, get(sb.tn_cns, m, String[]))), ","))
  end
  share(kind) = begin
    ds = filter(d -> d.kind === kind, diffs)
    same = sum(d.common + d.renamed for d in ds; init = 0)
    union = sum(d.common + d.renamed + d.only_a + d.only_b for d in ds; init = 0)
    union == 0 ? 1.0 : same / union
  end
  overlap = share(:structural)
  churn = 1.0 - share(:injection)
  gap = (sa.scenario_time === nothing || sb.scenario_time === nothing) ? nothing : Dates.value(sb.scenario_time - sa.scenario_time) / 60000.0
  reason = overlap < cfg.min_structural_overlap ? :structural_overlap : ((cfg.max_gap !== nothing && gap !== nothing && gap > cfg.max_gap) ? :max_gap : nothing)
  flag = cfg.max_injection_churn !== nothing && churn > cfg.max_injection_churn
  return CGMESSeriesPair(ia, ib, gap, overlap, churn, diffs, reason, flag)
end

"""
    checkCGMESSeries(inputs; config = active_sparlectra_config().qsts.series) -> CGMESSeriesReport

Series check of CGMES snapshots (`inputs`, a vector of
[`CGMESSeriesInput`](@ref)). The snapshots are grouped by model key and
ordered by time; every consecutive pair inside a model gets two measures:

- structural overlap: common mRIDs over the union, pooled over lines,
  transformers and their ends, busbar sections, breakers, disconnectors
  and connectivity nodes; below `config.min_structural_overlap` a new
  series starts;
- injection churn: 1 minus the same share over synchronous machines,
  energy consumers, equivalent injections, energy sources, generating
  units, power electronics connections and their units (CGMES 3.0); above
  `config.max_injection_churn` the pair is flagged, the series is not
  split.

A time step above `config.max_gap` (minutes) also starts a new series.

An object whose mRID changed between two exports still counts as common
when its fingerprint matches on both sides: connectivity nodes by name and
container, the other classes by name and the (matched) connectivity nodes
of their terminals, topological nodes (reported, not part of the overlap)
by their set of connectivity nodes. A fingerprint shared by several objects
on one side matches nothing and is counted as ambiguous. Objects of the
boundary set are not counted.

Per class the check also counts the common objects whose equipment
parameters changed (`drifted` in [`CGMESSeriesClassDiff`](@ref)), for
example a line with the same mRID and a new `r`. Drift is reported only.

Each snapshot is read once (EQ and TP profiles); only the counts, the
fingerprints and the numeric equipment parameters are kept.
"""
function checkCGMESSeries(inputs::AbstractVector{CGMESSeriesInput}; config::Sparlectra.QstsSeriesConfig = Sparlectra.active_sparlectra_config().qsts.series)::CGMESSeriesReport
  snaps = [_seriesSnapshot(i) for i in inputs]
  order = sortperm(eachindex(snaps); by = i -> (snaps[i].model_key, something(snaps[i].scenario_time, DateTime(0)), i))
  pairs = CGMESSeriesPair[]
  series = NamedTuple{(:model_key, :snapshots),Tuple{String,Vector{Int}}}[]
  prev = 0
  for i in order
    s = snaps[i]
    if prev != 0 && snaps[prev].model_key == s.model_key
      p = _seriesPair(snaps[prev], s, prev, i, config)
      push!(pairs, p)
      p.split_reason === nothing ? push!(series[end].snapshots, i) : push!(series, (model_key = s.model_key, snapshots = [i]))
    else
      push!(series, (model_key = s.model_key, snapshots = [i]))
    end
    prev = i
  end
  info = [(label = s.label, model_key = s.model_key, scenario_time = s.scenario_time, counts = Dict(c => length(v) for (c, v) in s.ids), has_tp = s.has_tp, ssh_tp_conflicts = s.ssh_tp_conflicts) for s in snaps]
  return CGMESSeriesReport(info, pairs, series, config)
end

"""
    printCGMESSeriesReport(io, report)

Plain-text series report: per series its snapshots, per snapshot the
switches on which SSH and TP disagree (case 3), per pair the time
step, the structural overlap, the injection churn, the reason for a split,
and one line per class with a difference (only in a, only in b, common,
renamed, ambiguous, parameter drift). A class whose only difference is
parameter drift gets its line too.
"""
function printCGMESSeriesReport(io::IO, r::CGMESSeriesReport)
  @printf(io, "CGMES series check: %d snapshot(s), %d series (min_structural_overlap %.2f, max_injection_churn %s, max_gap %s)\n", length(r.snapshots), length(r.series), r.config.min_structural_overlap, something(r.config.max_injection_churn, "none"), something(r.config.max_gap, "none"))
  for (k, s) in enumerate(r.series)
    println(io, "series ", k, " (model ", isempty(s.model_key) ? "?" : s.model_key, "): ", join((r.snapshots[i].label for i in s.snapshots), ", "))
  end
  # SSH against TP: named per snapshot, so a synthetic series with switch
  # events (TP not recomputed) does not read as a data error
  for sn in r.snapshots
    sn.ssh_tp_conflicts > 0 || continue
    println(io, "snapshot ", sn.label, ": ", sn.ssh_tp_conflicts, " switch(es) closed in the SSH, not retained, between two topological nodes (SSH against TP, case 3); the import turns them into links")
  end
  for p in r.pairs
    a = r.snapshots[p.a]
    b = r.snapshots[p.b]
    @printf(io, "pair %s -> %s: gap %s min, structural overlap %.4f, injection churn %.4f%s%s\n", a.label, b.label, p.gap_minutes === nothing ? "?" : string(round(p.gap_minutes; digits = 1)), p.structural_overlap, p.injection_churn, p.split_reason === nothing ? "" : ", new series ($(p.split_reason))", p.churn_flag ? ", churn above limit" : "")
    for d in p.classes
      (d.only_a + d.only_b + d.renamed + d.ambiguous + d.drifted) == 0 && continue
      @printf(io, "  %-22s %-10s only a %6d  only b %6d  common %6d  renamed %5d  ambiguous %5d  parameter drift %5d\n", d.class, d.kind, d.only_a, d.only_b, d.common, d.renamed, d.ambiguous, d.drifted)
    end
  end
  return nothing
end

# a label or model key may carry a comma or a quote; such a field is quoted
# (RFC 4180), numbers stay as Julia prints them
function _series_csv_field(x)::String
  s = string(x)
  any(c -> c in (',', '"', '\n', '\r'), s) || return s
  return string('"', replace(s, "\"" => "\"\""), '"')
end

"""
    writeCGMESSeriesReport(file, report)

The series report as CSV, one row per pair and class:
`model_key,a,b,gap_minutes,structural_overlap,injection_churn,split_reason,churn_flag,class,kind,only_a,only_b,common,renamed,ambiguous,drifted`.
`drifted` is the parameter drift of the class (see
[`CGMESSeriesClassDiff`](@ref)). A field with a comma, a quote or a line
break is quoted.
"""
function writeCGMESSeriesReport(file::AbstractString, r::CGMESSeriesReport)
  open(file, "w") do io
    println(io, "model_key,a,b,gap_minutes,structural_overlap,injection_churn,split_reason,churn_flag,class,kind,only_a,only_b,common,renamed,ambiguous,drifted")
    for p in r.pairs
      a = r.snapshots[p.a]
      b = r.snapshots[p.b]
      for d in p.classes
        println(io, join((_series_csv_field(x) for x in Any[a.model_key, a.label, b.label, something(p.gap_minutes, ""), p.structural_overlap, p.injection_churn, something(p.split_reason, ""), p.churn_flag, d.class, d.kind, d.only_a, d.only_b, d.common, d.renamed, d.ambiguous, d.drifted]), ","))
      end
    end
  end
  return file
end
