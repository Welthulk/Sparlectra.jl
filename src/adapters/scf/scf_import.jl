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

# file: src/adapters/scf/scf_import.jl
# purpose: Sparlectra Case Format (SCF) reader (issue #342). Loads a case
#          file through the EXISTING network constructors (addBus!,
#          addPIModelACLine!, addPIModelTrafo!, addProsumer!, addShunt!,
#          addLink!, applyTapNameplate!) and finishes with validate!, so
#          there is no second model-building path. Validation runs in three
#          stages: schema, reference integrity, then the model plausibility
#          checks the solver already owns. Unknown keys are hard errors.

# Rev. 1 component subset; anything else in `data` is rejected by name
const _SCF_KNOWN_DATA_COMPONENTS = ("node", "line", "generic_branch", "link", "source", "sym_load", "sym_gen", "shunt", "voltage_regulator", "sym_voltage_sensor", "sym_power_sensor", "sym_current_sensor", "fault")
const _SCF_KNOWN_ROOT_KEYS = ("version", "type", "is_batch", "attributes", "data", "sparlectra")
const _SCF_KNOWN_SPARLECTRA_KEYS = ("format_version", "meta", "roles", "components", "extra", "measurements", "start_state", "contingencies", "short_circuit", "config", "transformer_types", "branch_shunt_split", "tap_changer_models", "scenarios")
# what `sparlectra.measurements` may carry besides the rows; shunt_bay_rows
# are the bus-referenced shunt-bay readings PGM has no sensor for (revision
# 1.0 only; in 1.2 they are rows of `rows`)
const _SCF_KNOWN_MEASUREMENT_KEYS = ("rows", "shunt_bay_rows", "provenance")
const _SCF_SHUNT_BAY_TYPES = ("ShuntQ", "Imag", "Ia")
const _SCF_KNOWN_COMPONENT_KINDS = ("tap_changer", "transformer3w", "sc_source", "controllers", "shunt_state", "matpower_dcline", "for001_contingency")

_scf_get(d::AbstractDict, key::AbstractString, default = nothing) = get(d, key, default)

function _scf_require(d::AbstractDict, key::AbstractString, context::AbstractString)
  haskey(d, key) || throw(ArgumentError("SCF: $(context) is missing the required key \"$(key)\"."))
  return d[key]
end

# A voltage-dependent controller in `extra.<machine>`: `qu_control` and
# `pu_control` carry the characteristic points in per unit (voltage, power on
# `meta.s_base`), the interpolation mode and the limits in MVAr / MW. The file
# is validated HERE, with the component's name: a bad characteristic that only
# surfaced as a solver error three calls later would name nothing a reader
# can find in the file.
const _SCF_CONTROL_INTERPOLATIONS = ("linear", "spline", "polynomial")

function _scf_read_characteristic(obj, key::AbstractString, name::AbstractString)
  obj isa AbstractDict || throw(ArgumentError("SCF: $(name) $(key) must be an object with `points` and optionally `interpolation` and limits"))
  raw = get(obj, "points", nothing)
  (raw isa AbstractVector && length(raw) >= 2) || throw(ArgumentError("SCF: $(name) $(key).points needs at least two [u_pu, value_pu] pairs"))
  points = Tuple{Float64,Float64}[]
  for (i, pair) in enumerate(raw)
    (pair isa AbstractVector && length(pair) == 2 && all(x -> x isa Number, pair)) || throw(ArgumentError("SCF: $(name) $(key).points[$(i)] must be a [u_pu, value_pu] pair of numbers"))
    push!(points, (Float64(pair[1]), Float64(pair[2])))
  end
  for i in Iterators.drop(eachindex(points), 1)
    points[i][1] > points[i-1][1] || throw(ArgumentError("SCF: $(name) $(key).points must have strictly increasing voltages (point $(i) has $(points[i][1]) pu after $(points[i-1][1]) pu)"))
  end
  mode = string(get(obj, "interpolation", "linear"))
  mode in _SCF_CONTROL_INTERPOLATIONS || throw(ArgumentError("SCF: $(name) $(key).interpolation \"$(mode)\" is unknown; use one of $(join(_SCF_CONTROL_INTERPOLATIONS, ", "))"))
  return make_characteristic(points; interpolation = Symbol(mode))
end

# a limit stated in the object, in MW / MVAr; absent means no limit, and the
# infinity sentinel means the same thing
function _scf_control_limits(obj, key::AbstractString, lo_key::AbstractString, hi_key::AbstractString, name::AbstractString)
  lo = haskey(obj, lo_key) ? scf_number_or_sentinel(obj[lo_key], "extra.$(key).$(lo_key)") : nothing
  hi = haskey(obj, hi_key) ? scf_number_or_sentinel(obj[hi_key], "extra.$(key).$(hi_key)") : nothing
  (lo !== nothing && hi !== nothing && lo > hi) && throw(ArgumentError("SCF: $(name) $(key) has $(lo_key) = $(lo) above $(hi_key) = $(hi)"))
  return lo, hi
end

_scf_num(v, context::AbstractString)::Float64 = v isa ScfDecimalText ? v.value : v isa Number ? Float64(v) : throw(ArgumentError("SCF: $(context) must be a number, got $(repr(v))."))

# meta.f_nom: the system frequency the file states, or nothing when absent.
# A value that is not a positive finite number is a named error: f_nom = 0
# used to set every line b to zero without a word, and JSON `true` passed
# as 1 Hz.
function _scf_meta_f_nom(meta)::Union{Nothing,Float64}
  (meta isa AbstractDict && haskey(meta, "f_nom")) || return nothing
  v = meta["f_nom"]
  (v isa Real && !(v isa Bool) && isfinite(v) && v > 0) || throw(ArgumentError("SCF: meta.f_nom must be a positive number in Hz (the system frequency), got $(repr(v))."))
  return Float64(v)
end

_scf_int(v, context::AbstractString)::Int = v isa Integer ? Int(v) : (v isa AbstractFloat && isinteger(v) ? Int(v) : throw(ArgumentError("SCF: $(context) must be an integer, got $(repr(v)).")))

# extra.<appliance>.gencost: {"p": [...], "q": [...]} with MATPOWER gencost
# rows (q optional), see GenCost
function _scf_gencost(v, name::AbstractString)::GenCost
  v isa AbstractDict || throw(ArgumentError("SCF: extra.gencost of $(name) must be an object with \"p\" (and optionally \"q\")."))
  haskey(v, "p") || throw(ArgumentError("SCF: extra.gencost of $(name) has no \"p\" row."))
  row(x, which) = (x isa AbstractVector && all(y -> y isa Real, x)) ? Float64.(x) : throw(ArgumentError("SCF: extra.gencost.$(which) of $(name) must be a list of numbers."))
  return GenCost(row(v["p"], "p"), haskey(v, "q") ? row(v["q"], "q") : nothing)
end

# extra.<machine>.dynamics: the DynamicData fields under their field names
# (#506, revision 1.3). An unknown key is refused by name: the block is
# small, and a misspelt field would otherwise vanish without a word. A value
# that breaks the rule of its field fails with the DynamicDataError text.
function _scf_dynamics(v, name::AbstractString)::DynamicData
  v isa AbstractDict || throw(ArgumentError("SCF: extra.dynamics of $(name) must be an object."))
  known = String.(collect(fieldnames(DynamicData)))
  kw = Dict{Symbol,Any}()
  for (k, x) in v
    key = String(k)
    key in known || throw(ArgumentError("SCF: extra.dynamics of $(name) has the unknown key \"$(key)\". Known keys: $(join(known, ", "))."))
    if key == "unit_type"
      x isa AbstractString || throw(ArgumentError("SCF: extra.dynamics.unit_type of $(name) must be a string, got $(repr(x))."))
      kw[:unit_type] = try
        parse_dynamic_unit_type(x)
      catch e
        # the one failure parse_dynamic_unit_type has: a value outside the
        # three spellings
        e isa ArgumentError || rethrow()
        throw(ArgumentError("SCF: extra.dynamics of $(name): $(e.msg)."))
      end
    elseif key == "black_start"
      x isa Bool || throw(ArgumentError("SCF: extra.dynamics.black_start of $(name) must be true or false, got $(repr(x))."))
      kw[:black_start] = x
    else
      (x isa Real && !(x isa Bool)) || throw(ArgumentError("SCF: extra.dynamics.$(key) of $(name) must be a number, got $(repr(x))."))
      kw[Symbol(key)] = Float64(x)
    end
  end
  return try
    DynamicData(; kw...)
  catch e
    # the constructor's own rule check, reported with the machine's name
    e isa DynamicDataError || rethrow()
    throw(ArgumentError("SCF: extra.dynamics of $(name): $(sprint(showerror, e))."))
  end
end

# extra.<appliance>.reference_priority: a non-negative integer with CGMES
# semantics (0 = no preference, 1 = strongest); anything else is a broken
# file and is named, not clamped
function _scf_reference_priority(v, name::AbstractString)::Int
  p = _scf_int(v, "extra.reference_priority of $(name)")
  p >= 0 || throw(ArgumentError("SCF: extra.reference_priority of $(name) must be >= 0 (0 = no preference, 1 = strongest), got $(p)."))
  return p
end

## --- stage 1: schema -------------------------------------------------------

"""
    _scf_require_current_revision(spar)

The revision decides whether this reader understands the file at all.
Backward compatibility across revisions is not promised while the format is
young, so a file from another revision is refused BY NAME instead of being
read with today's meaning. Every reader that feeds a RUN goes through this:
the full import, and the config/studies accessors that are consulted before
the import. A file written by power-grid-model has no namespaced block and
never reaches the check; the display-only readers (measurement provenance,
source reference) stay lenient on purpose, so a stale file still renders
enough for the user to select it and be told to re-export.
"""
function _scf_require_current_revision(spar::AbstractDict)
  found = String(_scf_get(spar, "format_version", ""))
  isempty(found) && throw(ArgumentError("SCF: the file states no sparlectra.format_version; expected $(SCF_FORMAT_VERSION). Re-export the case with this version of Sparlectra."))
  # "1.1" is the number this format carried before the release: 1.0 was never
  # published, so the two are the SAME structure under two names. Files
  # written in that window are read as 1.0 - rejecting them would cost users
  # their exported cases over a renaming they never saw. 1.0 differs from the
  # current revision only in the unit of the current sensors, which the
  # measurement reader undoes (_scf_current_divisor). 1.2 lacks only the
  # optional extra.<machine>.dynamics block of 1.3 and reads as written. Any
  # other revision is refused as before.
  if found == "1.1"
    @info "SCF: file states format revision 1.1, which is the pre-release name of 1.0 (same structure); reading it as 1.0. Re-export to update the file."
  elseif !(found in _SCF_READABLE_REVISIONS)
    throw(ArgumentError("SCF: format revision $(found) found, $(SCF_FORMAT_VERSION) expected. Backward compatibility across revisions is not promised; re-export the case with this version of Sparlectra."))
  end
  return nothing
end

# Revisions the reader accepts: the current one and the older ones whose only
# difference it knows how to undo.
const _SCF_READABLE_REVISIONS = ("1.0", "1.1", "1.2", SCF_FORMAT_VERSION)

# Divisor from a file's sym_current_sensor numbers to ampere. Revisions 1.0
# and 1.1 wrote the current as if the estimator held kA (value * 1000) while
# it holds A, so their numbers are 1000 times too large; dividing by 1000,
# as the old reader did, gives back bit for bit what it read. From 1.2 on
# the file carries A, like PGM. A file without a sparlectra block (written
# by power-grid-model) states no revision and is read as the current one.
_scf_current_divisor(revision::AbstractString)::Float64 = revision in ("1.0", "1.1") ? 1.0e3 : 1.0

function _scf_validate_schema(root::AbstractDict)
  for key in keys(root)
    String(key) in _SCF_KNOWN_ROOT_KEYS || throw(ArgumentError("SCF: unknown top-level key \"$(key)\". Known keys: $(join(_SCF_KNOWN_ROOT_KEYS, ", "))."))
  end
  _scf_require(root, "version", "the root object")
  String(_scf_get(root, "type", "input")) == "input" || throw(ArgumentError("SCF: only input datasets are supported (type must be \"input\")."))
  _scf_get(root, "is_batch", false) === false || throw(ArgumentError("SCF: batch datasets (is_batch = true) are not a case-file input form."))
  data = _scf_get(root, "data", Dict{String,Any}())
  data isa AbstractDict || throw(ArgumentError("SCF: \"data\" must be an object."))
  for key in keys(data)
    String(key) in _SCF_KNOWN_DATA_COMPONENTS || throw(ArgumentError("SCF: unknown or unsupported component type \"$(key)\" in data. Supported in this version: $(join(_SCF_KNOWN_DATA_COMPONENTS, ", "))."))
  end
  spar = _scf_get(root, "sparlectra", nothing)
  if spar !== nothing
    spar isa AbstractDict || throw(ArgumentError("SCF: \"sparlectra\" must be an object."))
    _scf_require_current_revision(spar)
    for key in keys(spar)
      String(key) in _SCF_KNOWN_SPARLECTRA_KEYS || throw(ArgumentError("SCF: unknown key \"sparlectra.$(key)\". Known keys: $(join(_SCF_KNOWN_SPARLECTRA_KEYS, ", "))."))
    end
    # a bad frequency fails on load, before any line is built from it
    _scf_meta_f_nom(_scf_get(spar, "meta", nothing))
    comps = _scf_get(spar, "components", Dict{String,Any}())
    comps isa AbstractDict || throw(ArgumentError("SCF: \"sparlectra.components\" must be an object."))
    for key in keys(comps)
      String(key) in _SCF_KNOWN_COMPONENT_KINDS || throw(ArgumentError("SCF: unknown Sparlectra component type \"$(key)\". Known types: $(join(_SCF_KNOWN_COMPONENT_KINDS, ", "))."))
    end
  end
  return nothing
end

# study definitions: validated here so a broken study block fails on load,
# not at the start of a long sweep
# the removed study block (0.40.0): still a known key so a file carrying it
# fails with the way out instead of the generic unknown-key message
const SCF_CONTINGENCIES_REMOVED_MESSAGE = "SCF: the sparlectra.contingencies block was removed in 0.40.0. Write the study as a sparlectra.scenarios block (one scenario per case, a status op with value 0 per outage); the Web UI scenario editor writes it."
const _SCF_SHORT_CIRCUIT_CASES = ("max", "min")
const _SCF_SHORT_CIRCUIT_SWEEPS = ("all_buses", "explicit")

function _scf_validate_studies(root::AbstractDict)
  spar = _scf_get(root, "sparlectra", Dict{String,Any}())
  _scf_get(spar, "contingencies", nothing) === nothing || throw(ArgumentError(SCF_CONTINGENCIES_REMOVED_MESSAGE))
  # PGM's own fault vocabulary: Sparlectra evaluates the balanced bolted
  # three-phase fault, so anything else in the file would be silently ignored
  for row in _scf_get(_scf_get(root, "data", Dict{String,Any}()), "fault", [])
    row isa AbstractDict || throw(ArgumentError("SCF: every data.fault entry must be an object."))
    _scf_require(row, "fault_object", "a fault entry")
    ft = _scf_enum_str(_scf_get(row, "fault_type", "three_phase"), "fault_type")
    ft == "three_phase" || throw(ArgumentError("SCF: fault.fault_type must be three_phase (Sparlectra evaluates the balanced fault); got $(repr(ft))."))
    for key in ("r_f", "x_f")
      v = _scf_get(row, key, 0.0)
      _scf_num(v, "fault.$(key)") == 0.0 || throw(ArgumentError("SCF: fault.$(key) must be 0 (the short-circuit calculation is bolted); got $(v)."))
    end
  end
  meas_block = _scf_get(spar, "measurements", nothing)
  if meas_block isa AbstractDict
    # shunt_bay_rows is a block of revision 1.0 only; from 1.2 on the
    # shunt-bay readings are rows of measurements.rows
    known = String(_scf_get(spar, "format_version", "")) in ("1.0", "1.1") ? _SCF_KNOWN_MEASUREMENT_KEYS : filter(!=("shunt_bay_rows"), _SCF_KNOWN_MEASUREMENT_KEYS)
    for key in keys(meas_block)
      String(key) in known || throw(ArgumentError("SCF: unknown key \"$(key)\" in sparlectra.measurements (known: $(join(known, ", ")))."))
    end
  end
  sc = _scf_get(spar, "short_circuit", nothing)
  if sc isa AbstractDict
    case = String(_scf_get(sc, "case", "max"))
    case in _SCF_SHORT_CIRCUIT_CASES || throw(ArgumentError("SCF: short_circuit.case must be max or min; got $(repr(case))."))
    sweep = String(_scf_get(sc, "sweep", "all_buses"))
    sweep in _SCF_SHORT_CIRCUIT_SWEEPS || throw(ArgumentError("SCF: short_circuit.sweep must be all_buses or explicit; got $(repr(sweep))."))
    if sweep == "explicit"
      buses = _scf_get(sc, "buses", nothing)
      (buses isa AbstractVector && !isempty(buses)) || throw(ArgumentError("SCF: short_circuit.sweep is explicit, so short_circuit.buses must list the node ids to evaluate."))
    end
    c = _scf_get(sc, "c_factor", nothing)
    if c !== nothing
      cf = _scf_num(c, "short_circuit.c_factor")
      (cf == 0.0 || (0.5 <= cf <= 1.2)) || throw(ArgumentError("SCF: short_circuit.c_factor must be 0 (IEC Table 1) or within [0.5, 1.2]; got $(cf)."))
    end
  end
  return nothing
end

## --- stage 2: reference integrity ------------------------------------------

function _scf_validate_references(root::AbstractDict)
  data = _scf_get(root, "data", Dict{String,Any}())
  ids = Dict{Int,String}()
  for (kind, rows) in data
    rows isa AbstractVector || throw(ArgumentError("SCF: data.$(kind) must be a list."))
    for row in rows
      row isa AbstractDict || throw(ArgumentError("SCF: every entry of data.$(kind) must be an object."))
      id = _scf_int(_scf_require(row, "id", "a $(kind) entry"), "$(kind).id")
      haskey(ids, id) && throw(ArgumentError("SCF: duplicate id $(id) (used by $(ids[id]) and $(kind))."))
      ids[id] = String(kind)
    end
  end
  spar = _scf_get(root, "sparlectra", Dict{String,Any}())
  for (kind, rows) in _scf_get(spar, "components", Dict{String,Any}())
    rows isa AbstractVector || continue
    for row in rows
      row isa AbstractDict && haskey(row, "id") || continue
      id = _scf_int(row["id"], "$(kind).id")
      haskey(ids, id) && throw(ArgumentError("SCF: duplicate id $(id) (used by $(ids[id]) and $(kind))."))
      ids[id] = String(kind)
    end
  end
  # every reference must resolve, and to the right kind of object
  nodes = Set(id for (id, kind) in ids if kind == "node")
  branchlike = Set(id for (id, kind) in ids if kind in ("line", "generic_branch", "link"))
  check(id, allowed, what, kind) = (id in allowed) || throw(ArgumentError("SCF: $(kind) references $(what) $(id), which does not exist (or is not a $(what))."))
  for kind in ("line", "generic_branch", "link")
    for row in _scf_get(data, kind, [])
      check(_scf_int(_scf_require(row, "from_node", "a $(kind) entry"), "from_node"), nodes, "node", kind)
      check(_scf_int(_scf_require(row, "to_node", "a $(kind) entry"), "to_node"), nodes, "node", kind)
    end
  end
  for kind in ("source", "sym_load", "sym_gen", "shunt")
    for row in _scf_get(data, kind, [])
      check(_scf_int(_scf_require(row, "node", "a $(kind) entry"), "node"), nodes, "node", kind)
    end
  end
  appliances = Set(id for (id, kind) in ids if kind in ("sym_gen", "sym_load"))
  for row in _scf_get(data, "voltage_regulator", [])
    check(_scf_int(_scf_require(row, "regulated_object", "a voltage_regulator entry"), "regulated_object"), appliances, "generator or load", "voltage_regulator")
  end
  for row in _scf_get(data, "sym_voltage_sensor", [])
    check(_scf_int(_scf_require(row, "measured_object", "a sym_voltage_sensor entry"), "measured_object"), nodes, "node", "sym_voltage_sensor")
  end
  for kind in ("sym_power_sensor", "sym_current_sensor")
    for row in _scf_get(data, kind, [])
      obj = _scf_int(_scf_require(row, "measured_object", "a $(kind) entry"), "measured_object")
      (obj in nodes || obj in branchlike) || throw(ArgumentError("SCF: $(kind) references measured_object $(obj), which is neither a node nor a branch."))
    end
  end
  for row in _scf_get(_scf_get(spar, "components", Dict{String,Any}()), "tap_changer", [])
    check(_scf_int(_scf_require(row, "branch", "a tap_changer entry"), "branch"), branchlike, "branch", "tap_changer")
  end
  types = _scf_get(spar, "transformer_types", Dict{String,Any}())
  for grp in _scf_get(_scf_get(spar, "components", Dict{String,Any}()), "transformer3w", [])
    # Two forms: the EXPLICIT grouping of three existing legs (what the
    # writer emits), and the NAMEPLATE shorthand that describes the
    # transformer in CGMES PowerTransformerEnd terms and lets the reader
    # build the star equivalent. They are exclusive: a group either points
    # at legs or describes them.
    if haskey(grp, "nameplate") || haskey(grp, "type")
      haskey(grp, "ends") && throw(ArgumentError("SCF: a transformer3w entry has both a nameplate/type and explicit ends; use one form."))
      haskey(grp, "star_node") && throw(ArgumentError("SCF: a nameplate transformer3w creates its own star node; remove star_node."))
      plate = if haskey(grp, "type")
        tname = String(grp["type"])
        haskey(types, tname) || throw(ArgumentError("SCF: transformer3w references the unknown transformer type $(repr(tname))."))
        types[tname]
      else
        grp["nameplate"]
      end
      plate isa AbstractDict || throw(ArgumentError("SCF: a transformer3w nameplate must be an object."))
      pends = _scf_require(plate, "ends", "a transformer3w nameplate")
      (pends isa AbstractVector && length(pends) == 3) || throw(ArgumentError("SCF: a transformer3w nameplate needs exactly three ends (hv, mv, lv)."))
      seen_roles = String[]
      for e in pends
        e isa AbstractDict || throw(ArgumentError("SCF: every transformer3w nameplate end must be an object."))
        role = String(_scf_require(e, "role", "a transformer3w nameplate end"))
        role in ("hv", "mv", "lv") || throw(ArgumentError("SCF: a transformer3w nameplate end role must be hv, mv, or lv; got $(repr(role))."))
        push!(seen_roles, role)
        # the node comes from the INSTANCE when a shared type supplies the data
        nid = haskey(e, "node") ? e["node"] : begin
          inst = findfirst(x -> String(_scf_get(x, "role", "")) == role, _scf_get(grp, "nodes", []))
          inst === nothing && throw(ArgumentError("SCF: the $(role) end of a transformer3w has no node (give it in the nameplate end or in `nodes`)."))
          _scf_get(grp, "nodes", [])[inst]["node"]
        end
        check(_scf_int(nid, "transformer3w node"), nodes, "node", "transformer3w nameplate end")
        for key in ("rated_u", "rated_s", "r", "x")
          _scf_require(e, key, "a transformer3w nameplate end ($(role))")
        end
      end
      sort(seen_roles) == ["hv", "lv", "mv"] || throw(ArgumentError("SCF: a transformer3w nameplate needs one hv, one mv and one lv end; got $(join(sort(seen_roles), ", "))."))
      continue
    end
    check(_scf_int(_scf_require(grp, "star_node", "a transformer3w entry"), "star_node"), nodes, "node", "transformer3w")
    ends = _scf_require(grp, "ends", "a transformer3w entry")
    length(ends) == 3 || throw(ArgumentError("SCF: a transformer3w group must have exactly three ends, got $(length(ends))."))
    for e in ends
      check(_scf_int(_scf_require(e, "branch", "a transformer3w end"), "branch"), branchlike, "branch", "transformer3w end")
      check(_scf_int(_scf_require(e, "terminal_node", "a transformer3w end"), "terminal_node"), nodes, "node", "transformer3w end")
    end
  end
  return ids
end

## --- building the net through the existing constructors ---------------------

# The name and source id of a built branch, appliance, shunt or link, taken
# back from its extra record. The constructors generate both (B_ACL_110_1_2,
# #B_ACL_110_1_2#1); the format defines extra.name as the reference name of
# every component, so a hand-written case keeps its own names in the net,
# its reports and its measurement keys. A record without the key keeps the
# generated value (a PGM dataset has no extra block). An external_id marked
# `source_id_kind: cgmes_mrid` is the delivery's mRID, which the writer put
# in place of the component id; it lives in the restored mRID registry and
# is not a component id, so the generated id stays.
function _scf_restore_comp_identity!(comp, e::AbstractDict)
  haskey(e, "name") && (comp.cName = String(e["name"]))
  if haskey(e, "external_id") && get(e, "source_id_kind", "") != "cgmes_mrid"
    comp.cID = String(e["external_id"])
  end
  return comp
end

# reference name of an id: extra[<id>].name is mandatory wherever anything
# refers to an object by name, and it is what the writer always emits
function _scf_case_names(extra::AbstractDict, ids::AbstractDict)
  out = Dict{Int,String}()
  for (id, kind) in ids
    e = _scf_get(extra, string(id), nothing)
    if e isa AbstractDict && haskey(e, "name")
      out[id] = String(e["name"])
    else
      out[id] = string(uppercase(String(kind)), "_", id)
    end
  end
  return out
end

"""
    scf_to_net(root::AbstractDict) -> Net

Build a `Net` from a parsed SCF root object. Every element is created
through the existing public constructors and the network is validated with
`validate!` at the end, so an SCF case behaves exactly like a case from any
other importer. The document is converted to its typed form
([`SCFCase`](@ref)) first; this wrapper keeps the historical dict entry
point and applies NO configuration (callers that want the run-path
behavior use [`build_net`](@ref)); only a file without `meta.f_nom` takes
the active `model.nominal_frequency_hz` as its system frequency.
"""
function scf_to_net(root::AbstractDict)::Net
  _scf_validate_schema(root)
  _scf_validate_studies(root)
  _scf_validate_references(root)
  return _scf_net_from_case(scfcase_from_root(root))
end

# reference name and kind of every data id in the typed case; the same map
# _scf_validate_references derives on the dict, rebuilt here because the
# build path works on the struct and no longer sees the document
function _scf_case_id_kinds(case::SCFCase)::Dict{Int,String}
  ids = Dict{Int,String}()
  d = case.data
  for (kind, rows) in (("node", d.node), ("line", d.line), ("generic_branch", d.generic_branch), ("link", d.link), ("source", d.source), ("sym_load", d.sym_load), ("sym_gen", d.sym_gen), ("shunt", d.shunt), ("voltage_regulator", d.voltage_regulator), ("sym_voltage_sensor", d.sym_voltage_sensor), ("sym_power_sensor", d.sym_power_sensor), ("sym_current_sensor", d.sym_current_sensor), ("fault", d.fault))
    for row in rows
      ids[row.id] = kind
    end
  end
  return ids
end

# absent and JSON null both mean "no reading" on a sensor field
_scf_has(v) = !(v === missing || v === nothing)

"""
    scf_node_ids_in_build_order(case) -> Vector{Int}

The file node ids in the order [`build_net`](@ref) creates the buses
(recorded `extra.bus_index`, id as the tie break), so a consumer can map
net node index `k` to its file id, for example to look up the raw start
values of `start_state`.
"""
function scf_node_ids_in_build_order(case::SCFCase)::Vector{Int}
  extra = case.sparlectra === nothing ? Dict{String,Any}() : case.sparlectra.extra
  rows = [(order = _scf_int(_scf_get(_scf_get(extra, string(r.id), Dict{String,Any}()), "bus_index", typemax(Int)), "extra.bus_index"), id = r.id) for r in case.data.node]
  sort!(rows; by = e -> (e.order, e.id))
  return Int[e.id for e in rows]
end

# f_nom_default: the nominal frequency of a file that does not state f_nom
# (model.nominal_frequency_hz of the configuration the caller builds with)
function _scf_net_from_case(case::SCFCase; f_nom_default::Float64 = model_config().nominal_frequency_hz)::Net
  ids = _scf_case_id_kinds(case)
  data = case.data
  spar = case.sparlectra === nothing ? SCFSparlectra() : case.sparlectra
  meta = spar.meta
  extra = spar.extra
  names = _scf_case_names(extra, ids)
  s_base = haskey(meta, "s_base") ? _scf_num(meta["s_base"], "meta.s_base") / 1.0e6 : 100.0
  # rank 1 is the file's own statement, then the caller's configuration
  stated_f_nom = _scf_meta_f_nom(meta)
  f_nom = something(stated_f_nom, f_nom_default)
  net = Net(name = String(_scf_get(meta, "case_name", "scf_case")), baseMVA = s_base,
    vmin_pu = haskey(meta, "vmin_pu") ? _scf_num(meta["vmin_pu"], "meta.vmin_pu") : 0.9,
    vmax_pu = haskey(meta, "vmax_pu") ? _scf_num(meta["vmax_pu"], "meta.vmax_pu") : 1.1,
    nominal_frequency_hz = f_nom)
  # a stated frequency survives the configuration stamp of the import
  net._frequency_from_case = stated_f_nom !== nothing

  # --- buses -----------------------------------------------------------
  aux = Set(spar.roles.aux_nodes)
  vn_by_id = Dict{Int,Float64}()
  # the impedance base is computed from the file's own numbers with the same
  # formula the writer used (u_rated in V, s_base in VA), so a written and
  # re-read impedance is bit-identical instead of drifting by an ulp
  u_rated_by_id = Dict{Int,Float64}()
  s_base_va = s_base * 1.0e6
  bus_by_id = Dict{Int,String}()
  # The bus ORDER is part of the network's identity, exactly like the branch
  # order below: it fixes the Y-bus numbering, the generated component names
  # (which embed the bus index) and therefore the ids of a re-export. Ids are
  # assigned by NAME, which on a numerically named case (case118: "1", "10",
  # "100", ...) is not the bus order, so sorting by id here silently
  # renumbered the network. `extra.bus_index` is the recorded order; without
  # one the id order wins.
  node_rows = [(order = _scf_int(_scf_get(_scf_get(extra, string(r.id), Dict{String,Any}()), "bus_index", typemax(Int)), "extra.bus_index"), id = r.id, row = r) for r in data.node]
  sort!(node_rows; by = e -> (e.order, e.id))
  for entry in node_rows
    row = entry.row
    id = entry.id
    u_rated = row.u_rated
    vn = u_rated / 1.0e3
    vn_by_id[id] = vn
    u_rated_by_id[id] = u_rated
    name = names[id]
    e = _scf_get(extra, string(id), Dict{String,Any}())
    addBus!(
      net = net,
      busName = name,
      vn_kV = vn,
      isAux = id in aux,
      zone = haskey(e, "zone") ? _scf_int(e["zone"], "extra.zone") : nothing,
      area = haskey(e, "area") ? _scf_int(e["area"], "extra.area") : nothing,
      # the source system's bus number: the generated component names embed it
      oBusIdx = haskey(e, "source_index") ? _scf_int(e["source_index"], "extra.source_index") : nothing,
      # absent means the bus takes the network's limits
      vmin_pu = haskey(e, "vmin_pu") ? _scf_num(e["vmin_pu"], "extra.vmin_pu") : nothing,
      vmax_pu = haskey(e, "vmax_pu") ? _scf_num(e["vmax_pu"], "extra.vmax_pu") : nothing,
    )
    bus_by_id[id] = name
    # the source system's original display name (MATPOWER bus_name when the
    # import did not apply it as the reference name)
    if haskey(e, "original_name")
      net.busOriginalNameDict[length(net.nodeVec)] = String(e["original_name"])
    end
    # own bus names are a first-class channel:
    # extra.name is the reference name (busName above), external_id the
    # source system's id, component_name the internal name when it differs.
    # Restoring both here keeps a re-export of the imported net byte-stable;
    # without it addBus!'s generated Bus_<idx>_<vn> id would replace the
    # file's external_id on the second export.
    node_comp = net.nodeVec[end].comp
    haskey(e, "external_id") && (node_comp.cID = String(e["external_id"]))
    node_comp.cName = haskey(e, "component_name") ? String(e["component_name"]) : name
  end

  # --- start state (start values, never results) -------------------------
  # applied RIGHT AFTER the buses, so appliance setpoints laid down below
  # override the raw start voltage of their bus exactly like a direct
  # format import does (start_state carries the source system's raw
  # values; the setpoint precedence is the build's, not the file's)
  for (id, v) in spar.start_state.nodes
    bus = get(bus_by_id, id, nothing)
    bus === nothing && continue
    nd = net.nodeVec[net.busDict[bus]]
    v.vm_pu === nothing || (nd._vm_pu = v.vm_pu)
    v.va_deg === nothing || (nd._va_deg = v.va_deg)
  end

  # --- branches --------------------------------------------------------
  # Lines and transformers share ONE branch vector, and its order belongs to
  # the network's identity: measurements reference branch indices, and the
  # generated component ids embed them. The file records each element's
  # internal index in `extra`, so the reader rebuilds the original order
  # instead of a name-sorted one; without a recorded index the id order wins.
  zbase(id) = u_rated_by_id[id]^2 / s_base_va
  # The declared unit system (sparlectra.meta.units). In a units = pu file
  # the zb-coupled electrical parameters of lines, transformers and shunts
  # are read from the pu block of their extra entry as they stand; the SI
  # data section is a derived view for PGM readers and is not read for
  # them. The only conversion of the chain stays here, in build_net, and a
  # pu file needs none.
  units = scf_case_units(case)
  pu_block = function (id)
    e = _scf_get(extra, string(id), nothing)
    e isa AbstractDict || throw(ArgumentError("SCF: units = pu, but component $(id) has no extra entry carrying its pu values."))
    pu = _scf_get(e, "pu", nothing)
    pu isa AbstractDict || throw(ArgumentError("SCF: units = pu, but component $(id) carries no pu block in its extra entry."))
    return pu
  end
  pu_field = function (pu, field, id)
    haskey(pu, field) || throw(ArgumentError("SCF: units = pu, but component $(id) carries no $(field) in its pu block."))
    return _scf_num(pu[field], string("extra.pu.", field))
  end
  # ONE precomputed angular frequency, exactly as the writer uses it: the
  # multiplication order decides the last bit and therefore byte identity
  w_nom = 2.0 * pi * f_nom
  branch_index_by_id = Dict{Int,Int}()
  branch_rows = Any[]
  for (kind, rows) in (("line", data.line), ("generic_branch", data.generic_branch))
    for row in rows
      be = _scf_get(extra, string(row.id), Dict{String,Any}())
      order = haskey(be, "branch_index") ? _scf_int(be["branch_index"], "extra.branch_index") : row.id
      push!(branch_rows, (order = order, kind = kind, row = row, id = row.id))
    end
  end
  sort!(branch_rows; by = r -> (r.order, r.id))
  trafo_by_branch_idx = Dict{Int,PowerTransformer}()
  for entry in branch_rows
    row = entry.row
    id = entry.id
    to = row.to_node
    zb = zbase(to)
    if entry.kind == "line"
      b_s = row.c1 === nothing ? 0.0 : row.c1 * w_nom
      tan1 = row.tan1 === nothing ? 0.0 : row.tan1
      line_pu = units === :pu ? pu_block(id) : nothing
      # the frequency counts as used only where a capacitance enters the
      # model: a pu file reads b_pu as it stands
      line_pu === nothing && b_s != 0.0 && (net._frequency_used = true)
      addPIModelACLine!(
        net = net,
        fromBus = bus_by_id[row.from_node],
        toBus = bus_by_id[to],
        r_pu = line_pu === nothing ? row.r1 / zb : pu_field(line_pu, "r_pu", id),
        x_pu = line_pu === nothing ? row.x1 / zb : pu_field(line_pu, "x_pu", id),
        b_pu = line_pu === nothing ? b_s * zb : pu_field(line_pu, "b_pu", id),
        # the direct field wins: tan1 * b_s is 0 whenever the writer had to
        # fall back to g1 (conductance without capacitance)
        g_pu = line_pu === nothing ? (row.g1 === nothing ? tan1 * b_s * zb : row.g1 * zb) : pu_field(line_pu, "g_pu", id),
        status = 1,
        ratedS = row.i_n === nothing ? nothing : row.i_n * sqrt(3.0) * u_rated_by_id[to] / 1.0e6,
        from_status = something(row.from_status, 1),
        to_status = something(row.to_status, 1),
      )
    else
      # k and theta describe the LIVE tap; the neutral position comes from
      # the tap_changer entry when the file carries one, and the nameplate
      # is applied afterwards through the shared path
      gb_pu = units === :pu ? pu_block(id) : nothing
      addPIModelTrafo!(
        net = net,
        fromBus = bus_by_id[row.from_node],
        toBus = bus_by_id[to],
        r_pu = gb_pu === nothing ? row.r1 / zb : pu_field(gb_pu, "r_pu", id),
        x_pu = gb_pu === nothing ? row.x1 / zb : pu_field(gb_pu, "x_pu", id),
        b_pu = gb_pu === nothing ? (row.b1 === nothing ? 0.0 : row.b1 * zb) : pu_field(gb_pu, "b_pu", id),
        status = 1,
        ratio = something(row.k, 1.0),
        shift_deg = rad2deg(something(row.theta, 0.0)),
        ratedS = row.sn === nothing ? nothing : row.sn / 1.0e6,
        from_status = something(row.from_status, 1),
        to_status = something(row.to_status, 1),
      )
      if gb_pu !== nothing
        set_branch_shunt_total!(net.branchVec[end]; g_pu = pu_field(gb_pu, "g_pu", id), b_pu = net.branchVec[end].b_pu)
      elseif row.g1 !== nothing
        set_branch_shunt_total!(net.branchVec[end]; g_pu = row.g1 * zb, b_pu = net.branchVec[end].b_pu)
      end
      trafo_by_branch_idx[length(net.branchVec)] = net.trafos[end]
    end
    branch_index_by_id[id] = length(net.branchVec)
    # nameplate data the PGM dataset cannot express: an UNLIMITED rating (no
    # i_n/sn is "no limit", and the namespaced meta says which of the two it
    # is), plus the neutral ratio/shift of a transformer
    be_all = _scf_get(extra, string(id), Dict{String,Any}())
    _scf_restore_comp_identity!(net.branchVec[end].comp, be_all)
    meta_e = _scf_get(be_all, "meta", Dict{String,Any}())
    if meta_e isa AbstractDict && !isempty(meta_e)
      br = net.branchVec[end]
      haskey(meta_e, "sn_mva") && (br.sn_MVA = scf_number_or_sentinel(meta_e["sn_mva"], "extra.meta.sn_mva"))
      haskey(meta_e, "neutral_ratio") && (br.ratio = _scf_num(meta_e["neutral_ratio"], "extra.meta.neutral_ratio"))
      haskey(meta_e, "neutral_shift_deg") && (br.angle = _scf_num(meta_e["neutral_shift_deg"], "extra.meta.neutral_shift_deg"))
    end
    # the MATPOWER branch metadata record of the direct import:
    # report and export layers read it; the loss record itself stays with
    # the electrical values (g1/b1 of the row)
    mb = _scf_get(be_all, "matpower_branch", nothing)
    if mb isa AbstractDict
      net.matpower_branch_metadata[length(net.branchVec)] = (
        orig_name = haskey(mb, "orig_name") ? String(mb["orig_name"]) : nothing,
        orig_kind = Symbol(String(_scf_get(mb, "orig_kind", "unknown"))),
        orig_index = _scf_int(_scf_get(mb, "orig_index", 0), "matpower_branch.orig_index"),
        transformer_loss = nothing,
        tap_changer_model = haskey(mb, "tap_changer_model") ? Symbol(String(mb["tap_changer_model"])) : nothing,
        tap_impedance_correction_factor = _scf_num(_scf_get(mb, "tap_impedance_correction_factor", 1.0), "matpower_branch.tap_impedance_correction_factor"),
      )
    end
    # the delivered values of a CGMES EquivalentBranch
    eb = _scf_get(be_all, "equivalent_branch", nothing)
    if eb isa AbstractDict
      ebnum(k) = _scf_num(eb[k], string("extra.equivalent_branch.", k))
      net.equivalent_branch_metadata[length(net.branchVec)] = (r = ebnum("r_ohm"), x = ebnum("x_ohm"), r21 = ebnum("r21_ohm"), x21 = ebnum("x21_ohm"), vn_kV = ebnum("vn_kv"), stamped = eb["stamped"] === true)
    end
    if haskey(be_all, "transformer_loss_g_pu") && !isempty(net.trafos)
      net.trafos[end].side1.g = _scf_num(be_all["transformer_loss_g_pu"], "extra.transformer_loss_g_pu")
    end
    # DTF branch metadata: the full importer record, tagged
    # encoding, order preserving; export and report layers read it
    if haskey(be_all, "dtf_branch")
      net.matpower_branch_metadata[length(net.branchVec)] = scf_decode_value(be_all["dtf_branch"])
    end
    # the typed phase-tap model under extra: the 0.19 spelling, read for one
    # release; 0.20.0 writes sparlectra.tap_changer_models instead
    if haskey(be_all, "phase_taps") && entry.kind == "generic_branch" && !isempty(net.trafos)
      @warn "SCF: extra.phase_taps is deprecated since 0.20.0; the typed model is read, a re-export writes sparlectra.tap_changer_models" branch = id maxlog = 1
      net.trafos[end].side1.phase_taps = scf_decode_value(be_all["phase_taps"])
      note = resolve_branch_taps!(net.branchVec[end], net.trafos[end].side1)
      note === nothing || push!(net.tapModelNotices, note)
    end
  end
  for row in sort(data.link; by = r -> r.id)
    addLink!(net = net, fromBus = bus_by_id[row.from_node], toBus = bus_by_id[row.to_node], status = something(row.from_status, 1))
    _scf_restore_comp_identity!(net.linkVec[end], _scf_get(extra, string(row.id), Dict{String,Any}()))
  end

  # --- tap-changer nameplates (shared application path) -----------------
  comps = spar.components
  # changers whose entry carries no step (no pos, no step key) or the mark
  # step_assumed: they sit at the middle position, named once per import
  # below
  assumed_steps = Tuple{String,Int}[]
  for row in comps.tap_changer
    bidx = branch_index_by_id[row.branch]
    br = net.branchVec[bidx]
    # the file's ratio_base/angle_base_deg are the NEUTRAL position; the
    # live values were written into k/theta and are restored by the
    # nameplate application from the current step
    row.ratio_base === nothing || (br.ratio = row.ratio_base)
    row.angle_base_deg === nothing || (br.angle = row.angle_base_deg)
    tstep = 0.0
    tmin = 0.0
    tmax = 0.0
    tcur = 0.0
    pstep = 0.0
    pmin = 0.0
    pmax = 0.0
    pcur = 0.0
    psi = 90.0
    dustep = 0.0
    tgrid = :fraction
    pdeclared = true
    rdeclared = true
    # the mark of a written middle position, and that position
    rmark = false
    pmark = false
    mid = 0
    for ctrl in row.controllers
      idx = ctrl.index
      if idx == 1
        # declared false: the default ratio grid of a branch no nameplate
        # row declared, kept apart from a declared ratio changer
        ctrl.declared === false && (rdeclared = false)
        tstep = something(ctrl.step, 0.0)
        tmin = Float64(something(ctrl.pos_min, 0.0))
        tmax = Float64(something(ctrl.pos_max, 0.0))
        # no pos: NaN, the nameplate application puts the changer on its
        # neutral position 0 and marks the step as assumed
        tcur = ctrl.pos === nothing ? NaN : Float64(ctrl.pos)
        # step_assumed: the writer put the middle position into pos
        if ctrl.step_assumed === true && ctrl.pos !== nothing
          rmark = true
          mid = Int(round(ctrl.pos))
        end
        # the grid the positions count on; absent is the fraction grid
        # every file before the key was written on
        g = something(ctrl.grid, "fraction")
        g in ("fraction", "linear") || throw(ArgumentError("SCF: tap_changer controller grid must be \"fraction\" or \"linear\", got $(repr(g))."))
        tgrid = Symbol(g)
      elseif idx == 2
        # a stated alpha_deg is the direction as given, 0 included (a
        # Delta-u changer in phase with the voltage); a missing key reads
        # as the 90 deg symmetric PST, which is what the shipped cases carry
        ctrl.alpha_deg === nothing || (psi = ctrl.alpha_deg)
        # declared false: the constructor grid of a branch no nameplate row
        # declared, written so the re-read stays out of the mass release
        ctrl.declared === false && (pdeclared = false)
        pmin = Float64(something(ctrl.pos_min, 0.0))
        pmax = Float64(something(ctrl.pos_max, 0.0))
        pcur = ctrl.pos === nothing ? NaN : Float64(ctrl.pos)
        if ctrl.step_assumed === true && ctrl.pos !== nothing
          pmark = true
          rmark || (mid = Int(round(ctrl.pos)))
        end
        if ctrl.step_deg !== nothing
          pstep = ctrl.step_deg
        else
          dustep = something(ctrl.step, 0.0)
        end
      else
        throw(ArgumentError("SCF: tap_changer controller index must be 1 (ratio) or 2 (phase), got $(idx)."))
      end
    end
    # A source case can carry a neutral ratio OUTSIDE its own regulation band
    # (MATPOWER case57: ratio 0.895 with a 0.9 to 1.1 band), which is a data
    # inconsistency, not a file error. The nameplate insists on a band that
    # brackets the neutral, so it is applied on a widened band and the file's
    # exact ratio band is restored afterwards: the reader stays faithful, the
    # single construction path stays intact, and the oddity is named once.
    # a stated psi of 0 stays 0 (unlike MATPOWER, where 0 reads as 90), so a
    # ratio changer next to a Delta-u changer at psi 0 or 180 is a collinear
    # pair: with a Delta-u position off 0 the writers cannot split the
    # stored branch tap, refuse it. At position 0 and psi 0 the Delta-u
    # factor is 1 and the split is exact; the SCF writer gives such a file
    # for a pair stated at 0 or whose psi a caller set to 0, and it has to
    # read back. 180 is refused at every position, as applyTapNameplate! does.
    (tstep > 0.0 && dustep > 0.0 && abs(sin(deg2rad(psi))) < 1e-9 && !(psi == 0.0 && (isnan(pcur) || pcur == 0.0))) && throw(ArgumentError("SCF: sparlectra.components.tap_changer for branch $(row.branch): a ratio changer and a Delta-u changer at position $(pcur) with alpha_deg = $(psi) act along the same real direction; their positions cannot be told apart."))
    widened = tstep > 0.0 && !(tmin <= 0.0 <= tmax)
    if widened
      @warn "SCF: the neutral tap position lies outside the regulation band of this transformer; the band is restored as stated in the file" branch = row.branch pos_min = tmin pos_max = tmax
    end
    # the internal path builds psi as given; applyTapNameplate! would read a
    # stated 0 as 90 and the model would not match the alpha_deg it stores
    _apply_tap_nameplate!(br, psi; tap_step = tstep, tap_min_step = widened ? min(tmin, 0.0) : tmin, tap_max_step = widened ? max(tmax, 0.0) : tmax, tap_current_step = tcur, phase_step_deg = pstep, phase_min_step = pmin, phase_max_step = pmax, phase_current_step = pcur, phase_du_step = dustep, tap_ratio_grid = tgrid, context = "sparlectra.components.tap_changer for branch $(row.branch)")
    pdeclared || (br.phase_tap_declared = false)
    rdeclared || (br.ratio_tap_declared = false)
    (rmark && tstep > 0.0) && (br.ratio_step_assumed = true)
    (pmark && (pstep > 0.0 || dustep > 0.0)) && (br.phase_step_assumed = true)
    (br.ratio_step_assumed || br.phase_step_assumed) && push!(assumed_steps, ("branch $(row.branch)", mid))
    if widened
      # the band ends on the grid of the file (_tap_r1_at gives the cascade
      # fraction of a position on either grid)
      br.tap_min = br.ratio / (1.0 + _tap_r1_at(br, tmax))
      br.tap_max = br.ratio / (1.0 + _tap_r1_at(br, tmin))
    end
  end

  # --- transformers the file declares without a tap changer -------------
  # The branch constructor gives every transformer a default ratio and
  # phase grid (0.9 to 1.1 in 0.00625 steps, +-30 deg), not declared. The
  # writer states every changer a net carries as a tap_changer row (the
  # default grid with declared = false), so a generic_branch with neither a
  # row nor a typed model is a transformer without a changer, as in a
  # hand-written case or a PGM dataset. Keeping the default grid gave such a
  # branch a ratio changer of -14..17 positions (0.9 to 1.1 on the fraction
  # grid) that a tap controller could move and the tap estimation listed.
  # A branch whose taps the 0.19 extra.phase_taps model derived keeps them.
  with_changer = Set{Int}(row.branch for row in comps.tap_changer)
  for sid in keys(spar.tap_changer_models)
    # a key that is no integer is refused with its message further down
    bid = tryparse(Int, String(sid))
    bid === nothing || push!(with_changer, bid)
  end
  for row in data.generic_branch
    row.id in with_changer && continue
    br = net.branchVec[branch_index_by_id[row.id]]
    br.taps_derived && continue
    br.has_ratio_tap = false
    br.has_phase_tap = false
  end

  # --- asymmetric branch shunt split (0.20.0) ---------------------------
  # a row overrides the symmetric half of that branch; the totals stay the
  # sums, which equal the data section's totals up to the writer's rounding
  for (sid, d) in spar.branch_shunt_split
    bid = tryparse(Int, String(sid))
    bid === nothing && throw(ArgumentError("SCF: branch_shunt_split keys are branch ids, got $(repr(sid))."))
    haskey(branch_index_by_id, bid) || throw(ArgumentError("SCF: branch_shunt_split references unknown branch id $(bid)."))
    d isa AbstractDict || throw(ArgumentError("SCF: branch_shunt_split entry $(bid) must be an object."))
    set_branch_shunt!(net.branchVec[branch_index_by_id[bid]]; g_from_pu = _scf_num(_scf_require(d, "g_from_pu", "branch_shunt_split $(bid)"), "g_from_pu"), b_from_pu = _scf_num(_scf_require(d, "b_from_pu", "branch_shunt_split $(bid)"), "b_from_pu"), g_to_pu = _scf_num(_scf_require(d, "g_to_pu", "branch_shunt_split $(bid)"), "g_to_pu"), b_to_pu = _scf_num(_scf_require(d, "b_to_pu", "branch_shunt_split $(bid)"), "b_to_pu"))
  end

  # --- typed tap-changer models (0.20.0) --------------------------------
  # restored onto the from-side winding after the nameplate, then the
  # resolver derives ratio, shift, impedance and grid from the model. The
  # key `estimation` (0.60.0) is the tap changer set of the state
  # estimation and goes onto the branch as it stands; an entry that has
  # only this key runs no resolver (the branch keeps its nameplate taps).
  # An entry without the key reads as before.
  for (sid, d) in spar.tap_changer_models
    bid = tryparse(Int, String(sid))
    bid === nothing && throw(ArgumentError("SCF: tap_changer_models keys are branch ids, got $(repr(sid))."))
    haskey(branch_index_by_id, bid) || throw(ArgumentError("SCF: tap_changer_models references unknown branch id $(bid)."))
    bidx = branch_index_by_id[bid]
    haskey(trafo_by_branch_idx, bidx) || throw(ArgumentError("SCF: tap_changer_models entry $(bid) is not a generic_branch (transformer)."))
    d isa AbstractDict || throw(ArgumentError("SCF: tap_changer_models entry $(bid) must be an object."))
    if haskey(d, "ratio") || haskey(d, "phase") || !haskey(d, "estimation")
      w = trafo_by_branch_idx[bidx].side1
      ratio_assumed, phase_assumed = _scf_apply_tap_models!(w, d, bid)
      correction = Symbol(String(_scf_get(d, "tap_changer_model", "ideal")))
      # calcTapImpedanceCorrectionFactor is defined for the reciprocal
      # from-side convention; under the direct convention no physical chain
      # gives |1 + f e^(j psi)|^2 (docs/src/branchmodel.md), the factor is
      # applied as it stands and the file is told so once per transformer
      pt = w.phase_taps
      if correction === :impedance_correction && pt !== nothing && pt.kind === :asymmetrical && pt.convention === :direct_regulating_vector
        @warn "SCF: tap_changer_models $(bid): the tap impedance correction is defined for the convention reciprocal_from_side; this asymmetrical changer uses direct_regulating_vector, R and X are scaled with abs(1 + f e^(j psi))^2 as for the reciprocal convention." branch = bid
      end
      note = resolve_branch_taps!(net.branchVec[bidx], w; tap_changer_model = correction)
      note === nothing || push!(net.tapModelNotices, note)
      # a model replaces the nameplate changer of its role, its flag as well
      haskey(d, "ratio") && (net.branchVec[bidx].ratio_step_assumed = ratio_assumed)
      haskey(d, "phase") && (net.branchVec[bidx].phase_step_assumed = phase_assumed)
      if ratio_assumed || phase_assumed
        filter!(x -> x[1] != "branch $(bid)", assumed_steps)
        push!(assumed_steps, ("branch $(bid)", ratio_assumed ? w.taps.step : w.phase_taps.step))
      end
    end
    if haskey(d, "estimation")
      set = _scf_read_se_tap_changers(d["estimation"], bid)
      net.branchVec[bidx].se_tap_changers = set
      for ch in set.changers
        ch.step_assumed && push!(assumed_steps, ("branch $(bid)/$(ch.name)", Int(ch.step)))
      end
    end
  end
  for line in _tap_step_assumed_lines(assumed_steps)
    @info "SCF: " * line
  end

  # --- appliances -------------------------------------------------------
  regulator_by_object = Dict{Int,SCFVoltageRegulatorRow}()
  for row in data.voltage_regulator
    regulator_by_object[row.regulated_object] = row
  end
  slack = spar.roles.slack
  slack_nodes = Set(slack === nothing ? Int[] : slack.nodes)
  # A file written by power-grid-model itself has no `sparlectra.roles`: in
  # PGM a `source` IS the reference (an ideal voltage source behind its
  # impedance), there is no separate slack flag. Without this a plain PGM
  # dataset would load with no reference bus at all and could not be solved.
  if isempty(slack_nodes)
    for row in data.source
      something(row.status, 1) == 0 && continue
      push!(slack_nodes, row.node)
    end
  end
  # does the file mark its reference appliance at all?
  reference_marked = any(v isa AbstractDict && _scf_get(v, "reference_pri", false) === true for (_, v) in extra)
  participation = Dict{Int,Float64}()
  if slack !== nothing
    for p in slack.participation
      participation[p.object] = p.factor
    end
  end
  # the prosumer vector order is part of the identity as well (the generated
  # component ids embed the position), so the recorded element index decides
  appliance_rows = Any[]
  for (kind, rows) in (("source", data.source), ("sym_gen", data.sym_gen), ("sym_load", data.sym_load))
    for row in rows
      ae = _scf_get(extra, string(row.id), Dict{String,Any}())
      order = haskey(ae, "element_index") ? _scf_int(ae["element_index"], "extra.element_index") : row.id
      push!(appliance_rows, (order = order, kind = kind, row = row, id = row.id))
    end
  end
  sort!(appliance_rows; by = r -> (r.order, r.id))
  let
    for entry in appliance_rows
      kind = entry.kind
      row = entry.row
      id = entry.id
      node = row.node
      e = _scf_get(extra, string(id), Dict{String,Any}())
      ps_type = if kind == "source"
        "ExternalNetworkInjection"
      elseif kind == "sym_gen"
        String(_scf_get(e, "prosumption_type", "Generator"))
      else
        String(_scf_get(e, "prosumption_type", "EnergyConsumer"))
      end
      reg = get(regulator_by_object, id, nothing)
      vm = if kind == "source"
        # the reference itself: its setpoint IS the source, flag or not
        something(row.u_ref, 1.0)
      elseif reg !== nothing
        something(reg.u_ref, 1.0)
      elseif haskey(e, "vm_pu")
        # A setpoint on a machine that does not regulate is a contradiction,
        # not something to resolve: acting on it silently turned unregulated
        # machines into PV generators. Name it and stop.
        _scf_get(e, "regulated", false) === true ||
          throw(ArgumentError("SCF: $(_scf_get(e, "name", string(id))) carries a voltage setpoint (extra.vm_pu) but is not marked `regulated`; a setpoint on a machine that does not regulate has no effect. Remove the setpoint or set `regulated: true` (re-export the case)."))
        _scf_num(e["vm_pu"], "extra.vm_pu")
      elseif kind == "sym_gen" && reg === nothing && row.u_ref !== nothing
        # A dataset written by power-grid-model has no namespaced block, so the
        # rule above cannot apply: there is nothing that states whether the
        # machine regulates. The setpoint is taken (PGM means it), but the
        # decision is named instead of being silent.
        @warn "SCF: a generator states a voltage setpoint and the file carries no Sparlectra block to say whether it regulates; treating it as regulating" id = id u_ref = row.u_ref
        row.u_ref
      else
        nothing
      end
      # a MATPOWER-converted PQ generator carries its limits as constant
      # P(U)/Q(U) controllers; the explicit flag keeps every
      # existing case file's behavior untouched
      pu_ctrl = nothing
      qu_ctrl = nothing
      # An explicit characteristic (`qu_control` / `pu_control`) says more
      # than the flag and wins over it; a file carrying both is told so, by
      # component name, rather than silently taking one of them.
      ctrl_name = string(_scf_get(e, "name", string(id)))
      has_characteristic = haskey(e, "qu_control") || haskey(e, "pu_control")
      if has_characteristic && _scf_get(e, "pq_gen_controller", false) === true
        @info "SCF: $(ctrl_name) carries pq_gen_controller and an explicit characteristic; the characteristic is used"
      end
      if haskey(e, "qu_control")
        qmin_c, qmax_c = _scf_control_limits(e["qu_control"], "qu_control", "qmin_mvar", "qmax_mvar", ctrl_name)
        qu_ctrl = QUController(_scf_read_characteristic(e["qu_control"], "qu_control", ctrl_name); qmin_MVAr = qmin_c, qmax_MVAr = qmax_c, sbase_MVA = s_base)
      end
      if haskey(e, "pu_control")
        pmin_c, pmax_c = _scf_control_limits(e["pu_control"], "pu_control", "pmin_mw", "pmax_mw", ctrl_name)
        pu_ctrl = PUController(_scf_read_characteristic(e["pu_control"], "pu_control", ctrl_name); pmin_MW = pmin_c, pmax_MW = pmax_c, sbase_MVA = s_base)
      end
      if !has_characteristic && _scf_get(e, "pq_gen_controller", false) === true
        p_mw = something(row.p_specified, 0.0) / 1.0e6
        q_mw = something(row.q_specified, 0.0) / 1.0e6
        p_pu = p_mw / s_base
        q_pu = q_mw / s_base
        pmin_c = haskey(e, "min_p_mw") ? scf_number_or_sentinel(e["min_p_mw"], "extra.min_p_mw") : -Inf
        pmax_c = haskey(e, "max_p_mw") ? scf_number_or_sentinel(e["max_p_mw"], "extra.max_p_mw") : Inf
        qmin_c = haskey(e, "min_q_mvar") ? scf_number_or_sentinel(e["min_q_mvar"], "extra.min_q_mvar") : -Inf
        qmax_c = haskey(e, "max_q_mvar") ? scf_number_or_sentinel(e["max_q_mvar"], "extra.max_q_mvar") : Inf
        pu_ctrl = PUController(make_characteristic([(0.0, p_pu), (2.0, p_pu)]); pmin_MW = pmin_c, pmax_MW = pmax_c, sbase_MVA = s_base)
        qu_ctrl = QUController(make_characteristic([(0.0, q_pu), (2.0, q_pu)]); qmin_MVAr = qmin_c, qmax_MVAr = qmax_c, sbase_MVA = s_base)
      end
      addProsumer!(
        net = net,
        busName = bus_by_id[node],
        type = ps_type,
        pu_controller = pu_ctrl,
        qu_controller = qu_ctrl,
        p = kind == "source" ? nothing : something(row.p_specified, 0.0) / 1.0e6,
        q = kind == "source" ? nothing : something(row.q_specified, 0.0) / 1.0e6,
        pMax = haskey(e, "max_p_mw") ? scf_number_or_sentinel(e["max_p_mw"], "extra.max_p_mw") : nothing,
        pMin = haskey(e, "min_p_mw") ? scf_number_or_sentinel(e["min_p_mw"], "extra.min_p_mw") : nothing,
        # the regulator row is the primary place; `extra` carries the band for
        # everything that has no such row (a PGM `source`, an unregulated
        # machine) and for an unlimited band, which JSON cannot express in the
        # dataset at all
        qMax = reg !== nothing && reg.q_max !== nothing ? reg.q_max / 1.0e6 : (haskey(e, "max_q_mvar") ? scf_number_or_sentinel(e["max_q_mvar"], "extra.max_q_mvar") : nothing),
        qMin = reg !== nothing && reg.q_min !== nothing ? reg.q_min / 1.0e6 : (haskey(e, "min_q_mvar") ? scf_number_or_sentinel(e["min_q_mvar"], "extra.min_q_mvar") : nothing),
        # The file states which appliance holds the reference. A file that
        # names none - a power-grid-model dataset, or one written by hand -
        # falls back to the slack node, which is the rule that applied before
        # the marker existed.
        referencePri = if reference_marked
          _scf_get(e, "reference_pri", false) === true ? bus_by_id[node] : nothing
        else
          node in slack_nodes ? bus_by_id[node] : nothing
        end,
        vm_pu = vm,
        va_deg = kind == "source" && row.u_ref_angle !== nothing ? rad2deg(row.u_ref_angle) : nothing,
        # The FLAG decides, not the presence of a regulator row: the writer
        # gives every machine such a row to carry its setpoint and Q band, so
        # deriving regulation from it turned unregulated (PQ) generators into
        # regulated ones on every read (measured on case1354: 260 regulated
        # machines became 933). A foreign PGM file has no `extra` entry at
        # all, and there the row IS the statement of regulation.
        # A `source` is the reference by definition: its u_ref is the slack
        # voltage whether or not a hand-written extra entry carries the flag
        # (a file without it solved with the slack at 1.0 instead of 1.02,
        # every bus 0.02 pu low; the exporter always writes the flag, so
        # written files are unaffected).
        isRegulated = kind == "source" ? true : (isempty(e) ? reg !== nothing : _scf_get(e, "regulated", false) === true),
        # refreshing the bus types per appliance is O(prosumers) EACH TIME, so
        # a case with n appliances pays n^2 (measured: 6.9 s for case13659,
        # 0.5 s with the deferral). The MATPOWER importer defers the same way
        # and refreshes once when the network is complete.
        defer_bus_type_refresh = true,
        participationFactor = get(participation, id, nothing),
        # optional since 0.30.2: a file without the key (every file written
        # before) states no priority, which is the default 0
        referencePriority = haskey(e, "reference_priority") ? _scf_reference_priority(e["reference_priority"], ctrl_name) : 0,
      )
      # ratings are not constructor keywords (the PV/APU flags are derived
      # from the bus type and the voltage setpoint); restore them on the
      # created prosumer so the round trip keeps the nameplate
      ps = net.prosumpsVec[end]
      _scf_restore_comp_identity!(ps.comp, e)
      haskey(e, "rated_s_mva") && (ps.ratedS = scf_number_or_sentinel(e["rated_s_mva"], "extra.rated_s_mva"))
      haskey(e, "rated_u_kv") && (ps.ratedU = scf_number_or_sentinel(e["rated_u_kv"], "extra.rated_u_kv"))
      # the APU flag is derived from the bus type AT CONSTRUCTION time, and the
      # reader sets the final bus types first, so it must come from the file
      _scf_get(e, "apu_node", false) === true && (ps.isAPUNode = true)
      # MATPOWER generator costs (#471): optional, validated like the MATPOWER
      # import (GenCost throws on a row that is not a model 1 or 2 row)
      haskey(e, "gencost") && (ps.gencost = _scf_gencost(e["gencost"], ctrl_name))
      # dynamic parameters (#506): optional, a file without the block (every
      # file before revision 1.3) leaves the field at nothing
      haskey(e, "dynamics") && (ps.dynamics = _scf_dynamics(e["dynamics"], ctrl_name))
      # addProsumer! auto-regulates anything that carries a voltage setpoint,
      # so restoring the setpoint of an UNREGULATED machine would silently
      # promote it to a PV generator. The file knows which it was.
      (isempty(e) || kind == "source") || (ps.isRegulated = _scf_get(e, "regulated", false) === true)
      # addProsumer! derives the component type from "is this a generator", so
      # the SPECIFIC type (external network injection, synchronous machine,
      # energy consumer) has to be restored here. Without it a PGM `source`
      # came back as a generator and left again as a `sym_gen`, silently
      # changing the model of a file written by power-grid-model.
      ps.comp.cTyp = toComponentTyp(kind == "source" ? "EXTERNALNETWORKINJECTION" : uppercase(ps_type))
    end
    # deferred above, once now that every appliance is in place
    refreshBusTypesFromProsumers!(net)
    _buildQLimits!(net)
  end

  # --- shunts -----------------------------------------------------------
  shunt_rows = Any[]
  shunt_index_by_id = Dict{Int,Int}()
  for row in data.shunt
    se = _scf_get(extra, string(row.id), Dict{String,Any}())
    order = haskey(se, "shunt_index") ? _scf_int(se["shunt_index"], "extra.shunt_index") : row.id
    push!(shunt_rows, (order = order, row = row, id = row.id))
  end
  sort!(shunt_rows; by = r -> (r.order, r.id))
  for entry in shunt_rows
    row = entry.row
    node = row.node
    # counterpart of the writer's division: multiply by the base itself
    zb_sh = zbase(node)
    sh_pu = units === :pu ? pu_block(entry.id) : nothing
    g_pu = sh_pu === nothing ? something(row.g1, 0.0) * zb_sh : pu_field(sh_pu, "g_pu", entry.id)
    b_pu = sh_pu === nothing ? something(row.b1, 0.0) * zb_sh : pu_field(sh_pu, "b_pu", entry.id)
    # addShunt! takes the nominal P/Q draw at 1 pu voltage in MW/MVar
    # addShunt! forms y = (P + jQ)/baseMVA, so Q carries the susceptance with
    # the SAME sign the file states; negating it turned a capacitor into a
    # reactor on every round trip
    addShunt!(net = net, busName = bus_by_id[node], pShunt = g_pu * s_base, qShunt = b_pu * s_base)
    if sh_pu !== nothing
      # addShunt! divides the MW/MVar draw by baseMVA again, and (g * s) / s
      # is not always g; a pu file states the admittance itself, so it is
      # stamped as stated
      sh = net.shuntVec[end]
      sh.y_pu_shunt = ComplexF64(g_pu, b_pu)
      sh.G_shunt = g_pu
      sh.B_shunt = b_pu
    end
    # the generated id carries the net index; a file written from a net
    # whose controller shunt sat before a plain one has a different index
    # here, so the file's external_id is restored (same rule as for nodes)
    _scf_restore_comp_identity!(net.shuntVec[end].comp, _scf_get(extra, string(entry.id), Dict{String,Any}()))
    shunt_index_by_id[entry.id] = length(net.shuntVec)
  end

  # --- three-winding transformers from a nameplate ------------------------
  # The convenience form: instead of three generic_branch legs plus a star
  # node, the file states the transformer the way CGMES does (one
  # PowerTransformerEnd per winding) and the reader builds the star
  # equivalent through add3WTPiModelTrafo!. An export always writes the
  # explicit legs again, so this is an INPUT shorthand, not a second
  # representation of an existing network.
  ttypes = spar.transformer_types
  for grp in comps.transformer3w
    (grp.nameplate !== nothing || grp.type !== nothing) || continue
    plate = grp.type !== nothing ? ttypes[grp.type] : grp.nameplate
    by_role = Dict{String,Any}()
    for e in plate["ends"]
      by_role[String(e["role"])] = e
    end
    node_of = function (role)
      e = by_role[role]
      haskey(e, "node") && return _scf_int(e["node"], "transformer3w.$(role).node")
      idx = findfirst(x -> x.role == role, grp.nodes)
      return grp.nodes[idx].node
    end
    val(role, key, default = nothing) = begin
      e = by_role[role]
      haskey(e, key) ? _scf_num(e[key], "transformer3w.$(role).$(key)") : (default === nothing ? throw(ArgumentError("SCF: the $(role) end of a transformer3w is missing $(key).")) : default)
    end
    add3WTPiModelTrafo!(
      net = net,
      HBBus = bus_by_id[node_of("hv")],
      MBBus = bus_by_id[node_of("mv")],
      LVBus = bus_by_id[node_of("lv")],
      r = (val("hv", "r"), val("mv", "r"), val("lv", "r")),
      x = (val("hv", "x"), val("mv", "x"), val("lv", "x")),
      b = (val("hv", "b", 0.0), val("mv", "b", 0.0), val("lv", "b", 0.0)),
      ratedU_kV = (val("hv", "rated_u") / 1.0e3, val("mv", "rated_u") / 1.0e3, val("lv", "rated_u") / 1.0e3),
      ratedS_MVA = (val("hv", "rated_s") / 1.0e6, val("mv", "rated_s") / 1.0e6, val("lv", "rated_s") / 1.0e6),
      status = something(grp.status, 1),
    )
  end

  # --- shunt state (what PGM's g1/b1 cannot say) --------------------------
  for row in comps.shunt_state
    _scf_apply_shunt_state!(net.shuntVec[shunt_index_by_id[row.shunt]], row)
  end

  # --- feeder data a PGM `source` declares --------------------------------
  # PGM models a source as a voltage source BEHIND its impedance: `sk` and
  # `rx_ratio` are part of the model, not decoration. They arrive here as
  # declared case data in the same shape addExternalGrid! produces, so the
  # short circuit sees the feeder and `power_flow.external_grid` (source:
  # auto) computes with the file's own numbers instead of config defaults.
  # Our own files carry the richer components.sc_source block; that one wins:
  # a bus with a network injection or EquivalentInjection record there gets
  # no second record from the source's sk (the writer states sk from that
  # record, and the short circuit would stamp the feeder twice).
  declared_buses = Set(String(f.bus) for f in net.sc_sources.external_network_injections if f.bus !== nothing)
  for row in comps.sc_source
    row.kind in ("external_network_injection", "equivalent_injection") || continue
    b = get(row.fields, "bus", nothing)
    b isa AbstractString && push!(declared_buses, String(b))
  end
  for row in data.source
    row.sk === nothing && continue
    something(row.status, 1) == 0 && continue
    node = row.node
    busName = bus_by_id[node]
    busName in declared_buses && continue
    sk_va = row.sk
    (isfinite(sk_va) && sk_va > 0.0) || throw(ArgumentError("SCF: source.sk must be finite and > 0; got $(sk_va)."))
    rx = something(row.rx_ratio, 0.1)
    vn = getNodeVn(net.nodeVec[get(net.busDict, busName, node)])
    ik_a = sk_va / (sqrt(3.0) * vn * 1.0e3)
    push!(net.sc_sources.external_network_injections, (
      mrid = "scf-source-$(row.id)",
      name = names[row.id],
      bus = busName,
      maxInitialSymShCCurrent_A = ik_a,
      minInitialSymShCCurrent_A = nothing,
      maxR1ToX1Ratio = rx,
      minR1ToX1Ratio = nothing,
      maxR0ToX0Ratio = nothing,
      maxZ0ToZ1Ratio = row.z01_ratio,
      ikSecond = nothing,
      governorSCD = nothing,
    ))
  end

  # --- IEC 60909 source data ---------------------------------------------
  # runShortCircuit!(net) reads net.sc_sources, so a case file that carries
  # source data must restore it; without this a short circuit on an SCF case
  # would look plausible and mean nothing
  for row in comps.sc_source
    kind = row.kind
    fields = Dict{Symbol,Any}()
    for (key, v) in row.fields
      # JSON null is the file's spelling of "this field is empty"; "inf"/"-inf"
      # is how a non-finite source value travels (JSON has no Inf)
      fields[Symbol(key)] = v === nothing ? nothing : (v == "inf" || v == "-inf" ? scf_number_or_sentinel(v, "sc_source.$(key)") : v)
    end
    rec = NamedTuple{Tuple(keys(fields))}(Tuple(values(fields)))
    target = if kind == "external_network_injection"
      net.sc_sources.external_network_injections
    elseif kind == "synchronous_machine"
      net.sc_sources.synchronous_machines
    elseif kind == "asynchronous_machine"
      net.sc_sources.asynchronous_machines
    elseif kind == "equivalent_injection"
      net.sc_sources.equivalent_injections
    else
      throw(ArgumentError("SCF: unknown sc_source kind $(repr(kind))."))
    end
    push!(target, rec)
  end

  # --- machines the file names no source data for ------------------------
  # A synchronous machine feeds a fault whether or not the file says how
  # strongly. A unit without a source row is therefore a source as well: it
  # enters with the engine's documented default, and every row that
  # depends on it is flagged, as in every other format. The reference unit
  # of a bus that carries a feeder record IS that feeder (the external
  # grid) and is not counted a second time. The records are marked as not
  # declared, so an export writes the file's own rows and nothing else.
  feeder_buses = Set(String(f.bus) for f in net.sc_sources.external_network_injections if f.bus !== nothing)
  declared_machines = Set((String(m.bus), String(something(m.name, ""))) for m in net.sc_sources.synchronous_machines if m.bus !== nothing)
  machine_buses = Set(first(key) for key in declared_machines)
  bus_name_of = Dict{Int,String}(idx => name for (name, idx) in net.busDict)
  for ps in net.prosumpsVec
    (ps.comp.cTyp == Sparlectra.Generator || ps.comp.cTyp == Sparlectra.SynchronousMachine) || continue
    bus = get(bus_name_of, Int(ps.comp.cFrom_bus), nothing)
    bus === nothing && continue
    (isSlack(ps) && bus in feeder_buses) && continue
    ((bus, String(ps.comp.cName)) in declared_machines || bus in machine_buses) && continue
    vn = getNodeVn(net.nodeVec[Int(ps.comp.cFrom_bus)])
    push!(net.sc_sources.synchronous_machines, (
      mrid = string("scf-machine-", ps.comp.cID),
      name = String(ps.comp.cName),
      bus = bus,
      satDirectSubtransX_pu = nothing,
      satDirectTransX_pu = nothing,
      r0_pu = nothing,
      x0_pu = nothing,
      r2_pu = nothing,
      x2_pu = nothing,
      earthing = nothing,
      ratedS_MVA = (ps.ratedS !== nothing && ps.ratedS > 0.0) ? Float64(ps.ratedS) : nothing,
      ratedU_kV = vn > 0.0 ? vn : nothing,
      pmax_MW = (ps.maxP !== nothing && isfinite(ps.maxP) && ps.maxP > 0.0) ? Float64(ps.maxP) : nothing,
      declared = false,
    ))
  end

  # --- controllers (FACTS and regulation) --------------------------------
  # the entries are the declarative control.controllers schema, so they go
  # through applyConfiguredControllers! unchanged. Transformer references
  # use the file's branch names (extra block); they are translated to
  # branch indices first, so a name two branches share still resolves to
  # the branch of that id (a self-built case says trafo = "Moorau_110_20").
  ctrl_entries = comps.controllers
  if !isempty(ctrl_entries)
    branch_by_extra_name = Dict{String,Int}()
    for (id, idx) in branch_index_by_id
      nm = get(names, id, "")
      isempty(nm) || (branch_by_extra_name[String(nm)] = idx)
    end
    _scf_apply_controllers!(net, ctrl_entries; branch_by_name = branch_by_extra_name)
  end

  # --- measurements ------------------------------------------------------
  _scf_read_measurements!(net, case, bus_by_id, branch_index_by_id, u_rated_by_id)

  # --- CGMES source identity ----------------------------------------------
  # the structural-key mRID registry travels in the namespaced meta; with
  # it restored, the CGMES exporter and the mRID-addressed measurement
  # resolution work on a case-built net exactly as on the direct import
  ids_meta = _scf_get(meta, "cgmes_ids", nothing)
  if ids_meta isa AbstractDict
    for (k, v) in ids_meta
      net.cgmes_ids[String(k)] = String(v)
    end
  end

  # --- isolated nodes (roles.isolated_nodes) ------------------------------
  # restore the source system's isolated flags exactly like the direct
  # MATPOWER import (type Isolated, zero state, isoNodes registry)
  for id in spar.roles.isolated_nodes
    bus = get(bus_by_id, id, nothing)
    bus === nothing && continue
    idx = net.busDict[bus]
    setNodeType!(net.nodeVec[idx], "Isolated")
    setVmVa!(node = net.nodeVec[idx], vm_pu = 0.0, va_deg = 0.0)
    idx in net.isoNodes || push!(net.isoNodes, idx)
  end

  # --- MATPOWER dcline records --------------------------------------------
  # the injections themselves are ordinary sym_gen rows; this block restores
  # the metadata registry and the HVDC link records the result layer and the
  # dcline toggling read, plus the paired controllers when the case was
  # imported in paired_control mode
  orig_by_idx = net.busOrigIdxDict
  idx_by_orig = Dict{Int,Int}(orig => idx for (idx, orig) in orig_by_idx)
  for entry in comps.matpower_dcline
    fbus = _scf_int(_scf_require(entry, "from_bus", "a matpower_dcline entry"), "matpower_dcline.from_bus")
    tbus = _scf_int(_scf_require(entry, "to_bus", "a matpower_dcline entry"), "matpower_dcline.to_bus")
    f_idx = get(idx_by_orig, fbus, 0)
    t_idx = get(idx_by_orig, tbus, 0)
    (f_idx == 0 || t_idx == 0) && throw(ArgumentError("SCF: matpower_dcline references unknown source bus $(f_idx == 0 ? fbus : tbus)."))
    f_name = _scf_bus_name(net, f_idx)
    t_name = _scf_bus_name(net, t_idx)
    from_prosumer = _scf_int(_scf_require(entry, "from_element_index", "a matpower_dcline entry"), "matpower_dcline.from_element_index")
    to_prosumer = _scf_int(_scf_require(entry, "to_element_index", "a matpower_dcline entry"), "matpower_dcline.to_element_index")
    optnum(k) = haskey(entry, k) && entry[k] !== nothing ? _scf_num(entry[k], string("matpower_dcline.", k)) : nothing
    rec = (
      orig_index = _scf_int(_scf_get(entry, "orig_index", 0), "matpower_dcline.orig_index"),
      from_bus = fbus,
      to_bus = tbus,
      from_bus_name = f_name,
      to_bus_name = t_name,
      status = _scf_num(_scf_get(entry, "status", 1.0), "matpower_dcline.status"),
      pf_mw = _scf_num(_scf_get(entry, "pf_mw", 0.0), "matpower_dcline.pf_mw"),
      input_pt_mw = _scf_num(_scf_get(entry, "input_pt_mw", 0.0), "matpower_dcline.input_pt_mw"),
      effective_pt_mw = _scf_num(_scf_get(entry, "effective_pt_mw", 0.0), "matpower_dcline.effective_pt_mw"),
      pt_mw = _scf_num(_scf_get(entry, "effective_pt_mw", 0.0), "matpower_dcline.effective_pt_mw"),
      loss0_mw = _scf_num(_scf_get(entry, "loss0_mw", 0.0), "matpower_dcline.loss0_mw"),
      loss1 = _scf_num(_scf_get(entry, "loss1", 0.0), "matpower_dcline.loss1"),
      loss_mw = _scf_num(_scf_get(entry, "pf_mw", 0.0), "matpower_dcline.pf_mw") - _scf_num(_scf_get(entry, "effective_pt_mw", 0.0), "matpower_dcline.effective_pt_mw"),
      qf_mvar = _scf_num(_scf_get(entry, "qf_mvar", 0.0), "matpower_dcline.qf_mvar"),
      qt_mvar = _scf_num(_scf_get(entry, "qt_mvar", 0.0), "matpower_dcline.qt_mvar"),
      vf_pu = optnum("vf_pu"),
      vt_pu = optnum("vt_pu"),
      qminf_mvar = optnum("qminf_mvar"),
      qmaxf_mvar = optnum("qmaxf_mvar"),
      qmint_mvar = optnum("qmint_mvar"),
      qmaxt_mvar = optnum("qmaxt_mvar"),
      from_prosumer = from_prosumer,
      to_prosumer = to_prosumer,
      from_voltage_controlled = _scf_get(entry, "from_voltage_controlled", false) === true,
      to_voltage_controlled = _scf_get(entry, "to_voltage_controlled", false) === true,
    )
    push!(net.matpowerDclineMetadata, rec)
    push!(net.hvdcLinks, HvdcLink(string("DCLINE_", rec.orig_index), f_idx, t_idx, from_prosumer, to_prosumer, 1, :matpower, :p2p, nothing))
  end
  if !isempty(net.matpowerDclineMetadata) && _scf_get(meta, "matpower_dcline_paired_control", false) === true
    for m in net.matpowerDclineMetadata
      addHvdcPairControl!(
        net;
        from_bus = m.from_bus_name,
        to_bus = m.to_bus_name,
        p_transfer_mw = m.pf_mw,
        loss_mw = m.loss0_mw,
        loss_fraction = m.loss1,
        from_q_mvar = m.from_voltage_controlled ? nothing : m.qf_mvar,
        to_q_mvar = m.to_voltage_controlled ? nothing : m.qt_mvar,
        from_qmin_mvar = m.qminf_mvar,
        from_qmax_mvar = m.qmaxf_mvar,
        to_qmin_mvar = m.qmint_mvar,
        to_qmax_mvar = m.qmaxt_mvar,
        name = string("DCLINE_", m.orig_index),
        from_prosumer = m.from_prosumer,
        to_prosumer = m.to_prosumer,
      )
    end
  end

  # --- FOR001 contingency labels ------------------------------------------
  # net.for001Contingencies is a Vector{String}: both the DTF importer and
  # the MATPOWER metadata carry the outage LABELS, so the rows store one
  # label each
  for rec in comps.for001_contingency
    haskey(rec, "label") && push!(net.for001Contingencies, String(rec["label"]))
  end

  # --- stage 3: the model plausibility the solver already owns ------------
  ok, msg = validate!(net = net)
  ok || @warn "SCF: the imported network did not validate cleanly" message = msg
  return net
end

# one shunt_state row onto a built shunt: a components.shunt_state row, or
# the shunt_state object of a shunt voltage controller entry
function _scf_apply_shunt_state!(sh::Shunt, row::SCFShuntStateRow)
  row.status === nothing || (sh.status = row.status)
  if row.model !== nothing
    model = Symbol(row.model)
    model in (:Y, :VoltageDependentInjection, :PQ) || throw(ArgumentError("SCF: shunt_state.model must be Y, VoltageDependentInjection, or PQ; got $(repr(row.model))."))
    sh.model = model
  end
  row.estimate === nothing || (sh.estimate = row.estimate === true)
  return sh
end

# The measurements of a case file into net.measurements. Revision 1.2 lists
# one row per measured value with the semantics of the measurement table
# (_scf_read_measurement_rows!); revisions 1.0 and 1.1 list the rows per PGM
# sensor (_scf_read_measurements_v10!). A file without a sparlectra block
# (written by power-grid-model) has sensors only and goes the 1.0 way, with
# default ids and the sensor sigmas.
function _scf_read_measurements!(net::Net, case::SCFCase, bus_by_id, branch_index_by_id, u_rated_by_id)
  spar = case.sparlectra
  if spar === nothing || spar.format_version in ("1.0", "1.1")
    _scf_read_measurements_v10!(net, case, bus_by_id, branch_index_by_id, u_rated_by_id)
  else
    _scf_read_measurement_rows!(net, case, bus_by_id, branch_index_by_id, u_rated_by_id)
  end
  return nothing
end

# Revision 1.0 (and 1.1): PGM sensors back into Sparlectra measurement rows;
# measurements.rows supplies per sensor the original row ids, sigmas, active
# flags and positions, measurements.shunt_bay_rows the rows without a sensor.
# The 1.2 reader calls this for the sensors its rows do not reference
# (`skip` holds the consumed (sensor id, reading) pairs, `meta = false`
# ignores the per-sensor rows of the 1.0 layout); those come after the rows.
function _scf_read_measurements_v10!(net::Net, case::SCFCase, bus_by_id, branch_index_by_id, u_rated_by_id; skip::AbstractSet = Set{Tuple{Int,Int}}(), meta::Bool = true)
  data = case.data
  spar = case.sparlectra === nothing ? SCFSparlectra() : case.sparlectra
  rows = meta ? _scf_get(spar.measurements, "rows", []) : []
  meta_by_sensor = Dict{Int,Any}()
  for r in rows
    r isa AbstractDict && haskey(r, "sensor") || continue
    meta_by_sensor[_scf_int(r["sensor"], "measurements.rows.sensor")] = r
  end
  sid(id, i, fallback) = begin
    m = get(meta_by_sensor, id, nothing)
    m === nothing && return fallback
    idsv = _scf_get(m, "sparlectra_ids", nothing)
    (idsv isa AbstractVector && length(idsv) >= i) ? String(idsv[i]) : fallback
  end
  # the authoritative sigma in Sparlectra units: PGM carries ONE SI sigma per
  # sensor, which cannot express two different row sigmas (P and Q) exactly,
  # and the SI conversion is not bit-reversible. The measurement weight is
  # derived from sigma, so it is never stored.
  sigma_of(id, i, fallback) = begin
    m = get(meta_by_sensor, id, nothing)
    m === nothing && return fallback
    sg = _scf_get(m, "sigmas", nothing)
    (sg isa AbstractVector && length(sg) >= i) ? _scf_num(sg[i], "measurements.rows.sigmas") : fallback
  end
  # the PGM sensor is the authority for the value; an exact value is only
  # recorded where the SI conversion would not reproduce it bit for bit
  value_of(id, i, fallback) = begin
    m = get(meta_by_sensor, id, nothing)
    m === nothing && return fallback
    vv = _scf_get(m, "values", nothing)
    (vv isa AbstractVector && length(vv) >= i) ? _scf_num(vv[i], "measurements.rows.values") : fallback
  end
  active(id, i) = begin
    m = get(meta_by_sensor, id, nothing)
    m === nothing && return true
    a = _scf_get(m, "active", nothing)
    (a isa AbstractVector && length(a) >= i) ? a[i] === true : true
  end
  # collect with the recorded original position, then restore that order:
  # the WLS objective sums over the vector, so its order is part of the
  # numerical identity of an SE run
  collected = Tuple{Int,Any}[]
  next_pos = Ref(0)
  pos_of(id, i) = begin
    m = get(meta_by_sensor, id, nothing)
    if m !== nothing
      pv = _scf_get(m, "positions", nothing)
      if pv isa AbstractVector && length(pv) >= i
        return _scf_int(pv[i], "measurements.rows.positions")
      end
    end
    next_pos[] += 1
    return 1_000_000 + next_pos[]
  end
  keep!(id, i, m) = (id, i) in skip || push!(collected, (pos_of(id, i), m))
  for row in data.sym_voltage_sensor
    id = row.id
    node = row.measured_object
    ubase = u_rated_by_id[node]
    idx = net.busDict[bus_by_id[node]]
    i = 1
    if _scf_has(row.u_measured)
      keep!(id, i, Measurement(typ = VmMeas, value = value_of(id, i, row.u_measured / ubase), sigma = sigma_of(id, i, (_scf_has(row.u_sigma) ? row.u_sigma : 0.01 * ubase) / ubase), active = active(id, i), busIdx = idx, id = sid(id, i, "Vm_$(node)")))
      i += 1
    end
    if _scf_has(row.u_angle_measured)
      keep!(id, i, Measurement(typ = VaMeas, value = value_of(id, i, rad2deg(row.u_angle_measured)), sigma = sigma_of(id, i, 0.02), active = active(id, i), busIdx = idx, id = sid(id, i, "Va_$(node)")))
    end
  end
  for row in data.sym_power_sensor
    id = row.id
    obj = row.measured_object
    term = _scf_has(row.measured_terminal_type) ? row.measured_terminal_type : "node"
    sigma = (_scf_has(row.power_sigma) ? row.power_sigma : 1.0e6) / 1.0e6
    # PGM's optional per-quantity sigmas take precedence over the shared one
    psig = _scf_has(row.p_sigma) ? row.p_sigma : nothing
    qsig = _scf_has(row.q_sigma) ? row.q_sigma : nothing
    i = 1
    for (reading, indiv, typ) in ((row.p_measured, psig, term == "node" ? PinjMeas : PflowMeas), (row.q_measured, qsig, term == "node" ? QinjMeas : QflowMeas))
      _scf_has(reading) || continue
      value = value_of(id, i, reading / 1.0e6)
      sigma = indiv === nothing ? sigma : indiv / 1.0e6
      if term == "node"
        mid = sid(id, i, "$(typ)_$(obj)")
        # the ZI id rule of formats without a source field holds for 1.0
        # files only; sensors a 1.2 file names in no row carry default ids
        keep!(id, i, Measurement(typ = typ, value = value, sigma = sigma_of(id, i, sigma), active = active(id, i), busIdx = net.busDict[bus_by_id[obj]], id = mid, zero_injection = meta && _zero_injection_written_row(typ, value, mid)))
      else
        keep!(id, i, Measurement(typ = typ, value = value, sigma = sigma_of(id, i, sigma), active = active(id, i), branchIdx = branch_index_by_id[obj], direction = term == "branch_to" ? :to : :from, id = sid(id, i, "$(typ)_$(obj)")))
      end
      i += 1
    end
  end
  idiv = _scf_current_divisor(spar.format_version)
  for row in data.sym_current_sensor
    id = row.id
    obj = row.measured_object
    term = _scf_has(row.measured_terminal_type) ? row.measured_terminal_type : "branch_from"
    dir = term == "branch_to" ? :to : :from
    i = 1
    if _scf_has(row.i_measured)
      keep!(id, i, Measurement(typ = ImagMeas, value = value_of(id, i, row.i_measured / idiv), sigma = sigma_of(id, i, _scf_has(row.i_sigma) ? row.i_sigma / idiv : 1.0), active = active(id, i), branchIdx = branch_index_by_id[obj], direction = dir, id = sid(id, i, "Imag_$(obj)")))
      i += 1
    end
    if _scf_has(row.i_angle_measured)
      keep!(id, i, Measurement(typ = IaMeas, value = value_of(id, i, rad2deg(row.i_angle_measured)), sigma = sigma_of(id, i, rad2deg(_scf_has(row.i_angle_sigma) ? row.i_angle_sigma : deg2rad(0.05))), active = active(id, i), branchIdx = branch_index_by_id[obj], direction = dir, id = sid(id, i, "Ia_$(obj)")))
    end
  end
  # shunt-bay rows: internal units, no PGM sensor behind them
  for r in (meta ? _scf_get(spar.measurements, "shunt_bay_rows", []) : [])
    r isa AbstractDict || throw(ArgumentError("SCF: measurements.shunt_bay_rows entries must be objects, got $(repr(r))."))
    tname = String(_scf_get(r, "type", ""))
    tname in _SCF_SHUNT_BAY_TYPES || throw(ArgumentError("SCF: measurements.shunt_bay_rows type $(repr(tname)) is not one of $(join(_SCF_SHUNT_BAY_TYPES, ", "))."))
    typ = tname == "ShuntQ" ? ShuntQMeas : (tname == "Imag" ? ImagMeas : IaMeas)
    node = _scf_int(_scf_get(r, "node", nothing), "measurements.shunt_bay_rows.node")
    haskey(bus_by_id, node) || throw(ArgumentError("SCF: measurements.shunt_bay_rows references node $(node), which the file does not define."))
    p = _scf_int(_scf_get(r, "position", 0), "measurements.shunt_bay_rows.position")
    if p <= 0
      next_pos[] += 1
      p = 1_000_000 + next_pos[]
    end
    m = Measurement(typ = typ, value = _scf_num(_scf_get(r, "value", nothing), "measurements.shunt_bay_rows.value"), sigma = _scf_num(_scf_get(r, "sigma", nothing), "measurements.shunt_bay_rows.sigma"), active = _scf_get(r, "active", true) === true, busIdx = net.busDict[bus_by_id[node]], id = String(_scf_get(r, "id", "")))
    push!(collected, (p, m))
  end
  sort!(collected; by = e -> e[1])
  for (_, m) in collected
    push!(net.measurements, m)
  end
  return nothing
end

# keys a row of measurements.rows may carry (revision 1.2)
const _SCF_MEASUREMENT_ROW_KEYS = ("type", "id", "node", "terminal", "sensor", "value", "sigma", "valid", "zero_injection")

# Revision 1.2: measurements.rows with the semantics of the measurement
# table, one entry per measured value in the order of the measurement vector.
#
#   type      U, P, Q, Va at a node; P, Q, I, Ia at a terminal
#   node      node id (U, Va, and P/Q as the bus balance)
#   terminal  {"branch": id, "side": "from"|"to"} for a branch end, or
#             {"shunt_bay": node id} for the shunt bay of a node (Q summed
#             over its shunts, I and Ia of the bay)
#   sensor    the PGM sensor that carries the reading for a PGM reader
#   value     in kV (U), MW, MVAr, A or degrees; absent: read from the sensor
#   sigma     in the same unit; absent: the sensor sigma
#   valid     false for an inactive measurement; absent: true
#   zero_injection  true for a zero-injection row (P or Q at a node);
#             absent: false. The id decides nothing in 1.2 (the ZI id rule
#             applies to 1.0 and 1.1 files only).
#
# U is divided by u_rated / 1e3 of the node; a number with more than 17
# significant digits arrives as ScfDecimalText (the parser keeps the text of
# exactly these cells) and is divided from its text at 256 bits and rounded
# once, as the measurement table does. Sensors no row references
# are read after the rows with default ids, as in a file without the block.
function _scf_read_measurement_rows!(net::Net, case::SCFCase, bus_by_id, branch_index_by_id, u_rated_by_id)
  data = case.data
  rows = _scf_get(case.sparlectra.measurements, "rows", [])
  rows isa AbstractVector || throw(ArgumentError("SCF: measurements.rows must be a list."))
  volt = Dict(r.id => r for r in data.sym_voltage_sensor)
  power = Dict(r.id => r for r in data.sym_power_sensor)
  current = Dict(r.id => r for r in data.sym_current_sensor)
  idiv = _scf_current_divisor(case.sparlectra.format_version)
  consumed = Set{Tuple{Int,Int}}()
  ctx(k) = "measurements.rows[$(k)]"
  for (k, r) in enumerate(rows)
    r isa AbstractDict || throw(ArgumentError("SCF: $(ctx(k)) must be an object, got $(repr(r))."))
    for key in keys(r)
      String(key) in _SCF_MEASUREMENT_ROW_KEYS || throw(ArgumentError("SCF: unknown key \"$(key)\" in $(ctx(k)) (known: $(join(_SCF_MEASUREMENT_ROW_KEYS, ", ")))."))
    end
    t = String(_scf_get(r, "type", ""))
    (haskey(r, "node") == haskey(r, "terminal")) && throw(ArgumentError("SCF: $(ctx(k)) needs exactly one of node and terminal."))
    # location: node, branch end or shunt bay
    busIdx = nothing
    branchIdx = nothing
    dir = :none
    node = nothing
    bay = false
    if haskey(r, "node")
      node = _scf_int(r["node"], "$(ctx(k)).node")
      haskey(bus_by_id, node) || throw(ArgumentError("SCF: $(ctx(k)) references node $(node), which the file does not define."))
      busIdx = net.busDict[bus_by_id[node]]
      typ = t == "U" ? VmMeas : t == "Va" ? VaMeas : t == "P" ? PinjMeas : t == "Q" ? QinjMeas : throw(ArgumentError("SCF: $(ctx(k)) type $(repr(t)) at a node; expected U, P, Q or Va."))
    else
      term = r["terminal"]
      term isa AbstractDict || throw(ArgumentError("SCF: $(ctx(k)).terminal must be an object, got $(repr(term))."))
      if haskey(term, "shunt_bay") && length(term) == 1
        bay = true
        node = _scf_int(term["shunt_bay"], "$(ctx(k)).terminal.shunt_bay")
        haskey(bus_by_id, node) || throw(ArgumentError("SCF: $(ctx(k)) references node $(node), which the file does not define."))
        busIdx = net.busDict[bus_by_id[node]]
        typ = t == "Q" ? ShuntQMeas : t == "I" ? ImagMeas : t == "Ia" ? IaMeas : throw(ArgumentError("SCF: $(ctx(k)) type $(repr(t)) at a shunt bay; expected Q, I or Ia."))
      elseif haskey(term, "branch") && haskey(term, "side") && length(term) == 2
        bid = _scf_int(term["branch"], "$(ctx(k)).terminal.branch")
        haskey(branch_index_by_id, bid) || throw(ArgumentError("SCF: $(ctx(k)) references branch $(bid), which the file does not define."))
        side = String(term["side"])
        side in ("from", "to") || throw(ArgumentError("SCF: $(ctx(k)).terminal.side must be from or to, got $(repr(side))."))
        branchIdx = branch_index_by_id[bid]
        dir = Symbol(side)
        typ = t == "P" ? PflowMeas : t == "Q" ? QflowMeas : t == "I" ? ImagMeas : t == "Ia" ? IaMeas : throw(ArgumentError("SCF: $(ctx(k)) type $(repr(t)) at a branch end; expected P, Q, I or Ia."))
      else
        throw(ArgumentError("SCF: $(ctx(k)).terminal must be {\"branch\": id, \"side\": \"from\"|\"to\"} or {\"shunt_bay\": node id}, got $(repr(term))."))
      end
    end
    # the reading of the PGM sensor behind the row, in the unit the
    # estimator holds; its location has to be the row's
    sv = nothing
    ss = nothing
    if haskey(r, "sensor")
      bay && throw(ArgumentError("SCF: $(ctx(k)) is a shunt-bay row; PGM has no sensor for it."))
      sid = _scf_int(r["sensor"], "$(ctx(k)).sensor")
      mismatch = () -> throw(ArgumentError("SCF: $(ctx(k)) names sensor $(sid), which does not measure $(t) at the row's location."))
      if typ == VmMeas || typ == VaMeas
        s = get(volt, sid, nothing)
        (s === nothing || s.measured_object != node) && mismatch()
        ubase = u_rated_by_id[node]
        if typ == VmMeas
          _scf_has(s.u_measured) || mismatch()
          i = 1
          sv = s.u_measured / ubase
          ss = (_scf_has(s.u_sigma) ? s.u_sigma : 0.01 * ubase) / ubase
        else
          _scf_has(s.u_angle_measured) || mismatch()
          i = _scf_has(s.u_measured) ? 2 : 1
          sv = rad2deg(s.u_angle_measured)
          ss = 0.02
        end
      elseif typ in (PinjMeas, QinjMeas, PflowMeas, QflowMeas)
        s = get(power, sid, nothing)
        s === nothing && mismatch()
        sterm = _scf_has(s.measured_terminal_type) ? s.measured_terminal_type : "node"
        want = branchIdx === nothing ? "node" : (dir == :to ? "branch_to" : "branch_from")
        (sterm == want && s.measured_object == (branchIdx === nothing ? node : _scf_int(r["terminal"]["branch"], "$(ctx(k)).terminal.branch"))) || mismatch()
        isp = typ == PinjMeas || typ == PflowMeas
        reading = isp ? s.p_measured : s.q_measured
        _scf_has(reading) || mismatch()
        i = (!isp && _scf_has(s.p_measured)) ? 2 : 1
        indiv = isp ? s.p_sigma : s.q_sigma
        sv = reading / 1.0e6
        ss = (_scf_has(indiv) ? indiv : (_scf_has(s.power_sigma) ? s.power_sigma : 1.0e6)) / 1.0e6
      else
        s = get(current, sid, nothing)
        s === nothing && mismatch()
        sterm = _scf_has(s.measured_terminal_type) ? s.measured_terminal_type : "branch_from"
        (sterm == (dir == :to ? "branch_to" : "branch_from") && s.measured_object == _scf_int(r["terminal"]["branch"], "$(ctx(k)).terminal.branch")) || mismatch()
        if typ == ImagMeas
          _scf_has(s.i_measured) || mismatch()
          i = 1
          sv = s.i_measured / idiv
          ss = _scf_has(s.i_sigma) ? s.i_sigma / idiv : 1.0
        else
          _scf_has(s.i_angle_measured) || mismatch()
          i = _scf_has(s.i_measured) ? 2 : 1
          sv = rad2deg(s.i_angle_measured)
          ss = rad2deg(_scf_has(s.i_angle_sigma) ? s.i_angle_sigma : deg2rad(0.05))
        end
      end
      (sid, i) in consumed && throw(ArgumentError("SCF: $(ctx(k)) takes the $(t) reading of sensor $(sid), which an earlier row already took."))
      push!(consumed, (sid, i))
    end
    # a number of the row in the unit of the table, into the estimator's unit
    vn = typ == VmMeas ? u_rated_by_id[node] / 1.0e3 : 1.0
    own = key -> begin
      x = r[key]
      # a long literal arrives with its text (ScfDecimalText); U divides
      # the text like a table cell, everything else takes its Float64
      x isa ScfDecimalText && return typ == VmMeas ? _measurement_table_pu(x.text, x.value, vn) : x.value
      x isa Number || throw(ArgumentError("SCF: $(ctx(k)).$(key) must be a number, got $(repr(x))."))
      typ == VmMeas || return Float64(x)
      return Float64(x) / vn
    end
    value = haskey(r, "value") ? own("value") : (sv === nothing ? throw(ArgumentError("SCF: $(ctx(k)) has neither a value nor a sensor.")) : sv)
    sigma = haskey(r, "sigma") ? own("sigma") : (ss === nothing ? throw(ArgumentError("SCF: $(ctx(k)) has neither a sigma nor a sensor.")) : ss)
    sigma > 0.0 || throw(ArgumentError("SCF: $(ctx(k)).sigma must be positive, got $(sigma)."))
    valid = _scf_get(r, "valid", true)
    valid isa Bool || throw(ArgumentError("SCF: $(ctx(k)).valid must be true or false, got $(repr(valid))."))
    zi = _scf_get(r, "zero_injection", false)
    zi isa Bool || throw(ArgumentError("SCF: $(ctx(k)).zero_injection must be true or false, got $(repr(zi))."))
    (zi && !(typ == PinjMeas || typ == QinjMeas)) && throw(ArgumentError("SCF: $(ctx(k)).zero_injection is set on a $(t) row that is no bus balance; only P and Q at a node can be a zero injection."))
    id = haskey(r, "id") ? String(r["id"]) : _default_measurement_id(typ; busIdx = busIdx, branchIdx = branchIdx, direction = dir)
    push!(net.measurements, Measurement(typ = typ, value = value, sigma = sigma, active = valid, busIdx = busIdx, branchIdx = branchIdx, direction = dir, id = id, zero_injection = zi))
  end
  # sensors no row took (added by hand or by another tool) come after the
  # rows, read as in a file without the block
  _scf_read_measurements_v10!(net, case, bus_by_id, branch_index_by_id, u_rated_by_id; skip = consumed, meta = false)
  return nothing
end

"""
    build_net(case::SCFCase; config = active_sparlectra_config()) -> Net

The one network constructor of the run path over the typed case:
construct through the existing public
constructors, then apply the run configuration's net parameters exactly
once, here. `config` is the EFFECTIVE run configuration (the run path
passes its resolved `ImportedCase` configuration; interactive use takes
the active one).

System frequency: `meta.f_nom` of the file sets `net.nominal_frequency_hz`;
a file without it takes `model.nominal_frequency_hz` of `config`. When the
user set that key explicitly and it differs from the file, the file wins
with one warning that names both values.
"""
function build_net(case::SCFCase; config::SparlectraConfig = active_sparlectra_config())::Net
  net = _scf_net_from_case(case; f_nom_default = config.model.nominal_frequency_hz)
  key_f = config.model.nominal_frequency_hz
  if net._frequency_from_case && "model.nominal_frequency_hz" in config.user_set_keys && net.nominal_frequency_hz != key_f
    @warn "SCF: the case states meta.f_nom = $(net.nominal_frequency_hz) Hz, the configuration sets model.nominal_frequency_hz = $(key_f) Hz; the case value applies."
  end
  _apply_config_net_parameters!(net, config)
  # Safe AFTER construction on this path only: a case file records the model
  # of every shunt it carries (the shunt_state block), so the field is purely
  # the default for shunts added later. On DTF and CGMES the shunts are built
  # FROM this field, which is why it goes into their constructors instead.
  net.bus_shunt_model = normalize_bus_shunt_model(config.model.bus_shunt_model)
  return net
end

"""
    importSCF(file) -> Net

Read a Sparlectra Case Format file (`.scf.json`, issue #342) and build the
network through the existing constructors. Validation is staged: schema
(unknown keys and unsupported components are hard errors), reference
integrity (every id resolves, ids are unique), then the model checks
`validate!` already performs. Wrapper over [`read_scf_json`](@ref) and
[`build_net`](@ref) with the active configuration.
"""
importSCF(file::AbstractString)::Net = build_net(read_scf_json(file))

"""
    _scf_source_reference(file) -> String

The case file the SCF case was exported from (`sparlectra.meta.
source_reference`), empty when unknown. It is what lets a run tell a set
generated for the SOURCE case apart from a genuinely foreign one.
"""
function _scf_source_reference(file::AbstractString)::String
  isfile(file) || return ""
  lowercase(splitext(String(file))[2]) == ".json" || return ""
  try
    meta = get(get(scf_json_parse(read(String(file), String)), "sparlectra", Dict{String,Any}()), "meta", nothing)
    meta isa AbstractDict || return ""
    v = get(meta, "source_reference", "")
    return v isa AbstractString ? String(v) : ""
  catch
    return ""
  end
end

"""
    _scf_measurement_provenance(file) -> Dict{String,Any}

What a case file records about its own measurement set (`sparlectra.
measurements.provenance`): whether the values carry noise, which generator and
seed produced them, the per-row truth values, and the tap deviations the
generator applied. Empty for a case without measurements, without a provenance
block, or for anything that is not an SCF file, so callers can ask
unconditionally.
"""
function _scf_measurement_provenance(file::AbstractString)::Dict{String,Any}
  isfile(file) || return Dict{String,Any}()
  lowercase(splitext(String(file))[2]) == ".json" || return Dict{String,Any}()
  try
    block = get(get(scf_json_parse(read(String(file), String)), "sparlectra", Dict{String,Any}()), "measurements", nothing)
    block isa AbstractDict || return Dict{String,Any}()
    prov = get(block, "provenance", nothing)
    return prov isa AbstractDict ? Dict{String,Any}(String(k) => v for (k, v) in prov) : Dict{String,Any}()
  catch
    return Dict{String,Any}()
  end
end

"""
    scf_case_studies(file) -> NamedTuple

The study definitions a case file carries: `(short_circuit,)`, an empty
dictionary when the file does not define that study. The block states WHAT
to compute; the result contract of the run itself is unchanged. Contingency
studies live in the `scenarios` block ([`scf_case_scenarios`](@ref)); a file
that still carries the `contingencies` block removed in 0.40.0 fails here
with a message naming the replacement.
"""
function scf_case_studies(file::AbstractString)
  root = scf_json_parse(read(String(file), String))
  spar_rev = _scf_get(root, "sparlectra", Dict{String,Any}())
  isempty(spar_rev) || _scf_require_current_revision(spar_rev)
  _scf_validate_studies(root)
  spar = _scf_get(root, "sparlectra", Dict{String,Any}())
  sc = _scf_get(spar, "short_circuit", Dict{String,Any}())
  return (short_circuit = sc isa AbstractDict ? Dict{String,Any}(sc) : Dict{String,Any}(),)
end

"""
    scf_fault_nodes(file) -> Vector{Int}

The node ids a case file's `data.fault` rows point at, in file order. PGM
expresses "evaluate the fault here" with a `fault` component; a
short-circuit run uses these buses when the study block does not name its
own selection, so a case written in PGM vocabulary is runnable as it stands.
"""
function scf_fault_nodes(file::AbstractString)::Vector{Int}
  root = scf_json_parse(read(String(file), String))
  out = Int[]
  for row in _scf_get(_scf_get(root, "data", Dict{String,Any}()), "fault", [])
    row isa AbstractDict || continue
    _scf_int(_scf_get(row, "status", 1), "fault.status") == 0 && continue
    push!(out, _scf_int(_scf_require(row, "fault_object", "a fault entry"), "fault.fault_object"))
  end
  return out
end

"""
    scf_validate_dataset(root) -> Nothing

Run the file-level validation of a parsed case file: schema, study
definitions, reference integrity. Throws `ArgumentError` with the offending
key or id, and returns nothing when the document is a usable case file.

This is the check without building a network, which is what an upload path
needs: an arbitrary JSON must be rejected BEFORE it is stored, and the user
should see the reader's own message rather than a generic refusal.
"""
function scf_validate_dataset(root::AbstractDict)
  _scf_validate_schema(root)
  _scf_validate_studies(root)
  _scf_validate_references(root)
  return nothing
end

"""
    scf_extra_names(file) -> Dict{Int,String}

The reference name of every component id in a case file, taken from
`sparlectra.extra`. Study blocks address components by id, so a runner that
executes a study from the file resolves the ids through this map: the names
are exactly the ones the network object uses, which is what makes a case
list from the file addressable in the runner.
"""
function scf_extra_names(file::AbstractString)::Dict{Int,String}
  root = scf_json_parse(read(String(file), String))
  extra = _scf_get(_scf_get(root, "sparlectra", Dict{String,Any}()), "extra", Dict{String,Any}())
  out = Dict{Int,String}()
  for (id, entry) in extra
    entry isa AbstractDict && haskey(entry, "name") || continue
    out[parse(Int, String(id))] = String(entry["name"])
  end
  return out
end

"""
    scf_case_units(case::SCFCase) -> Symbol

The declared unit system of the electrical parameters:
`:pu` or `:si`. An ABSENT declaration means
`:si`, which is what every file before the declaration meant, so
existing files read unchanged. Unknown values are a hard error. For
`:pu` the conversion base must be complete: `meta.s_base` (VA,
positive and finite) and a positive finite `u_rated` on every node;
per-unit numbers without their base are not interpretable, so a
missing piece is a hard error naming the field. `build_net` reads the
branch and shunt parameters of a `:pu` file from the `pu` block of
their `sparlectra.extra` entry.
"""
function scf_case_units(case::SCFCase)::Symbol
  # a plain PGM dataset has no sparlectra block at all: that is si
  case.sparlectra === nothing && return :si
  raw = get(case.sparlectra.meta, "units", "si")
  raw isa AbstractString || throw(ArgumentError("sparlectra.meta.units must be a string (\"pu\" or \"si\"), got $(typeof(raw))."))
  u = Symbol(lowercase(strip(String(raw))))
  u in (:pu, :si) || throw(ArgumentError("sparlectra.meta.units must be \"pu\" or \"si\", got \"$(raw)\"."))
  if u === :pu
    haskey(case.sparlectra.meta, "s_base") || throw(ArgumentError("units = pu requires sparlectra.meta.s_base (the global power base in VA)."))
    s_base = _scf_num(case.sparlectra.meta["s_base"], "sparlectra.meta.s_base")
    (isfinite(s_base) && s_base > 0.0) || throw(ArgumentError("units = pu requires a positive finite sparlectra.meta.s_base, got $(s_base)."))
    for row in case.data.node
      (isfinite(row.u_rated) && row.u_rated > 0.0) || throw(ArgumentError("units = pu requires a positive finite u_rated on every node; node $(row.id) carries $(row.u_rated)."))
    end
  end
  return u
end

# the three facts build_net decides the frequency from, without building:
# the stated meta.f_nom, else the configured key; used when a line c1 > 0
# is read in an si file (the condition that sets net._frequency_used)
function _scf_case_frequency(case::SCFCase, config::SparlectraConfig)
  stated = case.sparlectra === nothing ? nothing : _scf_meta_f_nom(case.sparlectra.meta)
  used = scf_case_units(case) === :si && any(row -> row.c1 !== nothing && row.c1 != 0.0, case.data.line)
  return (f = something(stated, config.model.nominal_frequency_hz), used = used, from_case = stated !== nothing)
end

"""
    nominal_frequency_source(case::SCFCase, config::SparlectraConfig) -> String

`nominal_frequency_source(net)` for a
case before it is built, for example the base case of a QSTS series:
the same answer the network `build_net(case; config)` returns would give.
"""
function nominal_frequency_source(case::SCFCase, config::SparlectraConfig)::String
  fr = _scf_case_frequency(case, config)
  return _frequency_source_name(fr.used, fr.from_case)
end

"""
    nominal_frequency_line(case::SCFCase, config::SparlectraConfig) -> String

The run header line of `nominal_frequency_line(net)`
for a case before it is built (the base case of a QSTS series).
"""
function nominal_frequency_line(case::SCFCase, config::SparlectraConfig)::String
  fr = _scf_case_frequency(case, config)
  return _frequency_line_text(fr.f, fr.used, fr.from_case)
end

# Files written by Sparlectra releases before this one carried the whole
# form in sparlectra.config, machine-scope keys included; a case file is a
# published format, so those files keep loading: the machine keys are
# dropped and reported. From this release on the exporter never writes
# them, and a file that carries them is refused as before.
const _SCF_MACHINE_KEYS_TOLERATED_BELOW = v"0.11.0"

# the Sparlectra release that wrote the file, from meta.created_by
# ("Sparlectra 0.10.0"); nothing when the file does not say
function _scf_writer_version(spar::AbstractDict)::Union{Nothing,VersionNumber}
  meta = _scf_get(spar, "meta", Dict{String,Any}())
  meta isa AbstractDict || return nothing
  created = String(get(meta, "created_by", ""))
  m = match(r"^Sparlectra\s+(\d+\.\d+\.\d+)", created)
  m === nothing && return nothing
  return try
    VersionNumber(m.captures[1])
  catch
    nothing
  end
end

"""
    scf_case_config(file) -> Dict{String,Any}

The dotted configuration keys a case file carries in
`sparlectra.config`, ready to be merged as `config_overrides` (the case
file sits below API/CLI overrides and above the YAML file in precedence).
Empty when the file carries no configuration.
"""
function scf_case_config(file::AbstractString; dropped::Union{Nothing,Vector{String}} = nothing)::Dict{String,Any}
  root = scf_json_parse(read(String(file), String))
  spar = _scf_get(root, "sparlectra", Dict{String,Any}())
  # consulted BEFORE the import on the run path, so it must refuse a stale
  # revision itself or the run would first act on old configuration keys
  isempty(spar) || _scf_require_current_revision(spar)
  cfg = _scf_get(spar, "config", Dict{String,Any}())
  out = Dict{String,Any}(String(k) => v for (k, v) in cfg)
  # legacy key names first (0.10.x exports carried the model keys under
  # their old prefixes), then the scope check on the current names
  _apply_flat_config_aliases!(out, string("sparlectra.config in ", basename(String(file))))
  # A key outside the case scope is refused, not ignored: a file that states a
  # setting which quietly does not apply is exactly the silent acceptance this
  # format exists to prevent. Older files that carried the whole form need one
  # re-export.
  bad = sort!(String[k for k in keys(out) if !scf_is_case_config_key(k)])
  # a key this version does not know (removed or misspelled) is dropped
  # with a warning, as in a case configuration; known keys outside the case
  # scope stay refused
  bad = _drop_unknown_case_config_keys!(out, bad, string("sparlectra.config in ", basename(String(file))))
  writer = _scf_writer_version(spar)
  if !isempty(bad) && writer !== nothing && writer < _SCF_MACHINE_KEYS_TOLERATED_BELOW
    # an old file: the machine keys are dropped, the caller reports them
    for k in bad
      delete!(out, k)
    end
    dropped === nothing || append!(dropped, bad)
    return out
  end
  isempty(bad) || throw(ArgumentError("SCF: config key(s) $(join(bad, ", ")) are not case scope (logging, benchmarking, parallelism, Web UI and export settings belong in the configuration file). Re-export the case with format revision $(SCF_FORMAT_VERSION), or remove them."))
  return out
end

# Keys of the `estimation` object of a tap_changer_models entry (0.60.0),
# of one changer in it, of its reactance curve and of one table point. An
# unknown key is refused by name, a misspelt field would otherwise vanish.
const _SCF_SE_TAP_SET_KEYS = ["origin", "ratio_n", "shift_n_deg", "changers"]
const _SCF_SE_TAP_CHANGER_KEYS = ["name", "mrid", "kind", "role", "side", "neutral", "low", "high", "step", "step_assumed", "impedance", "additive", "du", "alpha_deg", "dphi_deg", "step_abs", "outside_du", "x_curve", "rx_correction", "table"]
const _SCF_SE_TAP_X_CURVE_KEYS = ["x_min", "x_max", "a_max_deg", "form"]
const _SCF_SE_TAP_POINT_KEYS = ["step", "ratio", "angle_deg", "r", "x", "g", "b"]

# every key of object `d` is in `known`, else an ArgumentError naming it
function _scf_check_keys(d::AbstractDict, known::Vector{String}, context::AbstractString)
  for k in keys(d)
    String(k) in known || throw(ArgumentError("SCF: $(context) has the unknown key \"$(k)\". Known keys: $(join(known, ", "))."))
  end
  return nothing
end

_scf_str(v, context::AbstractString)::String = v isa AbstractString ? String(v) : throw(ArgumentError("SCF: $(context) must be a string, got $(repr(v))."))
_scf_flag(v, context::AbstractString)::Bool = v isa Bool ? v : throw(ArgumentError("SCF: $(context) must be true or false, got $(repr(v))."))

# The tap changer set of the state estimation from the `estimation` object
# of tap_changer_models entry `bid` (written by _scf_se_tap_changers_dict).
# The SETapChanger and SETapChangerSet constructors check the values; their
# ArgumentError comes back with the entry named.
function _scf_read_se_tap_changers(e, bid::Int)::SETapChangerSet
  ctx = "tap_changer_models $(bid) estimation"
  e isa AbstractDict || throw(ArgumentError("SCF: $(ctx) must be an object."))
  _scf_check_keys(e, _SCF_SE_TAP_SET_KEYS, ctx)
  raw = _scf_require(e, "changers", ctx)
  raw isa AbstractVector || throw(ArgumentError("SCF: $(ctx).changers must be a list."))
  chs = SETapChanger[]
  for (k, c) in enumerate(raw)
    cc = "$(ctx).changers[$(k)]"
    c isa AbstractDict || throw(ArgumentError("SCF: $(cc) must be an object."))
    _scf_check_keys(c, _SCF_SE_TAP_CHANGER_KEYS, cc)
    req(key) = _scf_require(c, key, cc)
    opt(key) = haskey(c, key) ? _scf_num(c[key], "$(cc).$(key)") : nothing
    xc = nothing
    form = :sin2
    if haskey(c, "x_curve")
      x = c["x_curve"]
      x isa AbstractDict || throw(ArgumentError("SCF: $(cc).x_curve must be an object."))
      _scf_check_keys(x, _SCF_SE_TAP_X_CURVE_KEYS, "$(cc).x_curve")
      xc = (_scf_num(_scf_require(x, "x_min", "$(cc).x_curve"), "$(cc).x_curve.x_min"), _scf_num(_scf_require(x, "x_max", "$(cc).x_curve"), "$(cc).x_curve.x_max"), _scf_num(_scf_require(x, "a_max_deg", "$(cc).x_curve"), "$(cc).x_curve.a_max_deg"))
      form = Symbol(_scf_str(_scf_require(x, "form", "$(cc).x_curve"), "$(cc).x_curve.form"))
    end
    table = nothing
    if haskey(c, "table")
      c["table"] isa AbstractVector || throw(ArgumentError("SCF: $(cc).table must be a list."))
      table = SETapPoint[]
      for (j, pt) in enumerate(c["table"])
        pc = "$(cc).table[$(j)]"
        pt isa AbstractDict || throw(ArgumentError("SCF: $(pc) must be an object."))
        _scf_check_keys(pt, _SCF_SE_TAP_POINT_KEYS, pc)
        num(key) = _scf_num(_scf_require(pt, key, pc), "$(pc).$(key)")
        push!(table, SETapPoint(_scf_int(_scf_require(pt, "step", pc), "$(pc).step"), num("ratio"), num("angle_deg"), num("r"), num("x"), num("g"), num("b")))
      end
    end
    # no step key: the changer sits at its middle position, the neutral
    # step; step_assumed true: the written step is that middle position
    neutral = _scf_int(req("neutral"), "$(cc).neutral")
    assumed = !haskey(c, "step") || (haskey(c, "step_assumed") && _scf_flag(c["step_assumed"], "$(cc).step_assumed"))
    kw = (
      name = _scf_str(req("name"), "$(cc).name"),
      mrid = _scf_str(req("mrid"), "$(cc).mrid"),
      kind = Symbol(_scf_str(req("kind"), "$(cc).kind")),
      role = Symbol(_scf_str(req("role"), "$(cc).role")),
      side = Symbol(_scf_str(req("side"), "$(cc).side")),
      additive = haskey(c, "additive") ? _scf_flag(c["additive"], "$(cc).additive") : false,
      du = opt("du"),
      alpha_deg = opt("alpha_deg"),
      dphi_deg = opt("dphi_deg"),
      step_abs = opt("step_abs"),
      neutral = neutral,
      low = _scf_int(req("low"), "$(cc).low"),
      high = _scf_int(req("high"), "$(cc).high"),
      step = haskey(c, "step") ? _scf_num(c["step"], "$(cc).step") : Float64(neutral),
      step_assumed = assumed,
      table = table,
      outside_du = opt("outside_du"),
      x_curve = xc,
      x_curve_form = form,
      impedance = Symbol(_scf_str(req("impedance"), "$(cc).impedance")),
      rx_correction = haskey(c, "rx_correction") ? _scf_flag(c["rx_correction"], "$(cc).rx_correction") : false,
    )
    ch = try
      SETapChanger(; kw...)
    catch err
      # the field checks of the constructor (kind, role, side, band, the
      # data a kind needs): same text, with the entry named
      err isa ArgumentError || rethrow()
      throw(ArgumentError("SCF: $(cc): $(err.msg)"))
    end
    push!(chs, ch)
  end
  ratio_n = _scf_num(_scf_require(e, "ratio_n", ctx), "$(ctx).ratio_n")
  shift_n_deg = _scf_num(_scf_require(e, "shift_n_deg", ctx), "$(ctx).shift_n_deg")
  origin = Symbol(_scf_str(_scf_require(e, "origin", ctx), "$(ctx).origin"))
  # only a set an import stored is a branch record; typed and fields sets
  # are derived on demand, and a stored one would skip the write-back of
  # the steps after a tap estimation (_tap_apply_position!)
  origin in (:cgmes, :powsybl) || throw(ArgumentError("SCF: $(ctx).origin must be \"cgmes\" or \"powsybl\", got $(repr(String(origin)))."))
  return try
    SETapChangerSet(ratio_n = ratio_n, shift_n_deg = shift_n_deg, changers = chs, origin = origin)
  catch err
    # the set checks (ratio_n, origin, an additive changer first)
    err isa ArgumentError || rethrow()
    throw(ArgumentError("SCF: $(ctx): $(err.msg)"))
  end
end

# The typed tap-changer models of one sparlectra.tap_changer_models entry
# onto a winding: `ratio` (a PowerTransformerTaps from its constructor
# inputs) and `phase` (a PhaseTapChangerModel with kind, steps, the formula
# parameters and the table for :tabular); a key that is absent stays
# absent on the winding. A model without `step` sits at its middle
# position, its `neutral_step`; `step_assumed: true` marks a written step
# as that middle position. Returns whether the step of the ratio and of
# the phase model is assumed.
function _scf_apply_tap_models!(w::PowerTransformerWinding, d::AbstractDict, bid::Int)
  if haskey(d, "ratio")
    r = d["ratio"]
    r isa AbstractDict || throw(ArgumentError("SCF: tap_changer_models $(bid): ratio must be an object."))
    w.taps = PowerTransformerTaps(
      Vn_kV = _scf_num(_scf_require(r, "vn_kv", "tap_changer_models $(bid) ratio"), "vn_kv"),
      step = haskey(r, "step") ? _scf_int(r["step"], "step") : _scf_int(_scf_require(r, "neutral_step", "tap_changer_models $(bid) ratio"), "neutral_step"),
      lowStep = _scf_int(_scf_require(r, "low_step", "tap_changer_models $(bid) ratio"), "low_step"),
      highStep = _scf_int(_scf_require(r, "high_step", "tap_changer_models $(bid) ratio"), "high_step"),
      neutralStep = _scf_int(_scf_require(r, "neutral_step", "tap_changer_models $(bid) ratio"), "neutral_step"),
      voltageIncrement_kV = _scf_num(_scf_require(r, "voltage_increment_kv", "tap_changer_models $(bid) ratio"), "voltage_increment_kv"),
      neutralU = haskey(r, "neutral_u_kv") ? _scf_num(r["neutral_u_kv"], "neutral_u_kv") : nothing,
      convention = Symbol(String(_scf_get(r, "convention", "neutral_relative"))),
    )
  end
  if haskey(d, "phase")
    p = d["phase"]
    p isa AbstractDict || throw(ArgumentError("SCF: tap_changer_models $(bid): phase must be an object."))
    kind = Symbol(String(_scf_require(p, "kind", "tap_changer_models $(bid) phase")))
    optnum(key) = haskey(p, key) && p[key] !== nothing ? _scf_num(p[key], key) : nothing
    table = nothing
    if haskey(p, "table") && p["table"] !== nothing
      table = TapTablePoint[TapTablePoint(step = _scf_int(_scf_require(pt, "step", "table point"), "step"), ratio = _scf_num(_scf_get(pt, "ratio", 1.0), "ratio"), angle_deg = _scf_num(_scf_get(pt, "angle_deg", 0.0), "angle_deg"), x_pu = haskey(pt, "x_pu") && pt["x_pu"] !== nothing ? _scf_num(pt["x_pu"], "x_pu") : nothing) for pt in p["table"]]
    end
    w.phase_taps = PhaseTapChangerModel(
      kind = kind,
      step = haskey(p, "step") ? _scf_int(p["step"], "step") : _scf_int(_scf_require(p, "neutral_step", "tap_changer_models $(bid) phase"), "neutral_step"),
      lowStep = kind === :tabular ? nothing : _scf_int(_scf_require(p, "low_step", "tap_changer_models $(bid) phase"), "low_step"),
      highStep = kind === :tabular ? nothing : _scf_int(_scf_require(p, "high_step", "tap_changer_models $(bid) phase"), "high_step"),
      neutralStep = _scf_int(_scf_require(p, "neutral_step", "tap_changer_models $(bid) phase"), "neutral_step"),
      voltage_step_increment = optnum("voltage_step_increment"),
      step_phase_shift_increment = optnum("step_phase_shift_increment"),
      winding_connection_angle_deg = optnum("winding_connection_angle_deg"),
      x_min = optnum("x_min"),
      x_max = optnum("x_max"),
      convention = Symbol(String(_scf_get(p, "convention", "reciprocal_from_side"))),
      table = table,
    )
  end
  # the winding ratio and shift the models move from; the k and theta of
  # the data section already carry the steps, so a file that states them
  # is read without applying the steps a second time. Without the keys
  # the winding keeps k and theta as before.
  haskey(d, "base_ratio") && (w.ratio = _scf_num(d["base_ratio"], "tap_changer_models $(bid) base_ratio"))
  haskey(d, "base_shift_deg") && (w.shift_degree = _scf_num(d["base_shift_deg"], "tap_changer_models $(bid) base_shift_deg"))
  assumed(m) = !haskey(m, "step") || (haskey(m, "step_assumed") && _scf_flag(m["step_assumed"], "tap_changer_models $(bid) step_assumed"))
  return (haskey(d, "ratio") && assumed(d["ratio"]), haskey(d, "phase") && assumed(d["phase"]))
end
