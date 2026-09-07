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
const _SCF_KNOWN_SPARLECTRA_KEYS = ("format_version", "meta", "roles", "components", "extra", "measurements", "start_state", "contingencies", "short_circuit", "config", "transformer_types", "scenarios")
# what `sparlectra.measurements` may carry besides the rows
const _SCF_KNOWN_MEASUREMENT_KEYS = ("rows", "provenance")
const _SCF_KNOWN_COMPONENT_KINDS = ("tap_changer", "transformer3w", "sc_source", "controllers", "shunt_state", "matpower_dcline", "for001_contingency")

_scf_get(d::AbstractDict, key::AbstractString, default = nothing) = get(d, key, default)

function _scf_require(d::AbstractDict, key::AbstractString, context::AbstractString)
  haskey(d, key) || throw(ArgumentError("SCF: $(context) is missing the required key \"$(key)\"."))
  return d[key]
end

_scf_num(v, context::AbstractString)::Float64 = v isa Number ? Float64(v) : throw(ArgumentError("SCF: $(context) must be a number, got $(repr(v))."))
_scf_int(v, context::AbstractString)::Int = v isa Integer ? Int(v) : (v isa AbstractFloat && isinteger(v) ? Int(v) : throw(ArgumentError("SCF: $(context) must be an integer, got $(repr(v)).")))

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
  # their exported cases over a renaming they never saw. Any other revision
  # is refused as before.
  if found == "1.1" && SCF_FORMAT_VERSION == "1.0"
    @info "SCF: file states format revision 1.1, which is the pre-release name of 1.0 (same structure); reading it as $(SCF_FORMAT_VERSION). Re-export to update the file."
  else
    found == SCF_FORMAT_VERSION || throw(ArgumentError("SCF: format revision $(found) found, $(SCF_FORMAT_VERSION) expected. Backward compatibility across revisions is not promised; re-export the case with this version of Sparlectra."))
  end
  return nothing
end

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
const _SCF_CONTINGENCY_MODES = ("explicit", "all_branches", "all_branches_plus")
const _SCF_SHORT_CIRCUIT_CASES = ("max", "min")
const _SCF_SHORT_CIRCUIT_SWEEPS = ("all_buses", "explicit")

function _scf_validate_studies(root::AbstractDict)
  spar = _scf_get(root, "sparlectra", Dict{String,Any}())
  cont = _scf_get(spar, "contingencies", nothing)
  if cont isa AbstractDict
    mode = String(_scf_get(cont, "mode", "explicit"))
    mode in _SCF_CONTINGENCY_MODES || throw(ArgumentError("SCF: contingencies.mode must be one of $(join(_SCF_CONTINGENCY_MODES, ", ")); got $(repr(mode))."))
    for case in _scf_get(cont, "cases", [])
      case isa AbstractDict || throw(ArgumentError("SCF: every contingencies.cases entry must be an object."))
      outages = _scf_require(case, "outages", "a contingency case")
      (outages isa AbstractVector && !isempty(outages)) || throw(ArgumentError("SCF: a contingency case needs a non-empty outages list."))
      for o in outages
        _scf_require(o, "component", "a contingency outage")
      end
    end
  end
  # PGM's own fault vocabulary: Sparlectra evaluates the balanced bolted
  # three-phase fault, so anything else in the file would be silently ignored
  for row in _scf_get(_scf_get(root, "data", Dict{String,Any}()), "fault", [])
    row isa AbstractDict || throw(ArgumentError("SCF: every data.fault entry must be an object."))
    _scf_require(row, "fault_object", "a fault entry")
    ft = String(_scf_get(row, "fault_type", "three_phase"))
    ft == "three_phase" || throw(ArgumentError("SCF: fault.fault_type must be three_phase (Sparlectra evaluates the balanced fault); got $(repr(ft))."))
    for key in ("r_f", "x_f")
      v = _scf_get(row, key, 0.0)
      _scf_num(v, "fault.$(key)") == 0.0 || throw(ArgumentError("SCF: fault.$(key) must be 0 (the short-circuit calculation is bolted); got $(v)."))
    end
  end
  meas_block = _scf_get(spar, "measurements", nothing)
  if meas_block isa AbstractDict
    for key in keys(meas_block)
      String(key) in _SCF_KNOWN_MEASUREMENT_KEYS || throw(ArgumentError("SCF: unknown key \"$(key)\" in sparlectra.measurements (known: $(join(_SCF_KNOWN_MEASUREMENT_KEYS, ", ")))."))
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
behavior use [`build_net`](@ref)).
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

function _scf_net_from_case(case::SCFCase)::Net
  ids = _scf_case_id_kinds(case)
  data = case.data
  spar = case.sparlectra === nothing ? SCFSparlectra() : case.sparlectra
  meta = spar.meta
  extra = spar.extra
  names = _scf_case_names(extra, ids)
  s_base = haskey(meta, "s_base") ? _scf_num(meta["s_base"], "meta.s_base") / 1.0e6 : 100.0
  f_nom = haskey(meta, "f_nom") ? _scf_num(meta["f_nom"], "meta.f_nom") : 50.0
  net = Net(name = String(_scf_get(meta, "case_name", "scf_case")), baseMVA = s_base,
    vmin_pu = haskey(meta, "vmin_pu") ? _scf_num(meta["vmin_pu"], "meta.vmin_pu") : 0.9,
    vmax_pu = haskey(meta, "vmax_pu") ? _scf_num(meta["vmax_pu"], "meta.vmax_pu") : 1.1)

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
    # own bus names are a first-class channel (maintainer 2026-09-04):
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
  # format import does (D10: start_state carries the source system's raw
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
  for entry in branch_rows
    row = entry.row
    id = entry.id
    to = row.to_node
    zb = zbase(to)
    if entry.kind == "line"
      b_s = row.c1 === nothing ? 0.0 : row.c1 * w_nom
      tan1 = row.tan1 === nothing ? 0.0 : row.tan1
      addPIModelACLine!(
        net = net,
        fromBus = bus_by_id[row.from_node],
        toBus = bus_by_id[to],
        r_pu = row.r1 / zb,
        x_pu = row.x1 / zb,
        b_pu = b_s * zb,
        # the direct field wins: tan1 * b_s is 0 whenever the writer had to
        # fall back to g1 (conductance without capacitance)
        g_pu = row.g1 === nothing ? tan1 * b_s * zb : row.g1 * zb,
        status = 1,
        ratedS = row.i_n === nothing ? nothing : row.i_n * sqrt(3.0) * u_rated_by_id[to] / 1.0e6,
        from_status = something(row.from_status, 1),
        to_status = something(row.to_status, 1),
      )
    else
      # k and theta describe the LIVE tap; the neutral position comes from
      # the tap_changer entry when the file carries one, and the nameplate
      # is applied afterwards through the shared path
      addPIModelTrafo!(
        net = net,
        fromBus = bus_by_id[row.from_node],
        toBus = bus_by_id[to],
        r_pu = row.r1 / zb,
        x_pu = row.x1 / zb,
        b_pu = row.b1 === nothing ? 0.0 : row.b1 * zb,
        status = 1,
        ratio = something(row.k, 1.0),
        shift_deg = rad2deg(something(row.theta, 0.0)),
        ratedS = row.sn === nothing ? nothing : row.sn / 1.0e6,
        from_status = something(row.from_status, 1),
        to_status = something(row.to_status, 1),
      )
      row.g1 === nothing || (net.branchVec[end].g_pu = row.g1 * zb)
    end
    branch_index_by_id[id] = length(net.branchVec)
    # nameplate data the PGM dataset cannot express: an UNLIMITED rating (no
    # i_n/sn is "no limit", and the namespaced meta says which of the two it
    # is), plus the neutral ratio/shift of a transformer
    be_all = _scf_get(extra, string(id), Dict{String,Any}())
    meta_e = _scf_get(be_all, "meta", Dict{String,Any}())
    if meta_e isa AbstractDict && !isempty(meta_e)
      br = net.branchVec[end]
      haskey(meta_e, "sn_mva") && (br.sn_MVA = scf_number_or_sentinel(meta_e["sn_mva"], "extra.meta.sn_mva"))
      haskey(meta_e, "neutral_ratio") && (br.ratio = _scf_num(meta_e["neutral_ratio"], "extra.meta.neutral_ratio"))
      haskey(meta_e, "neutral_shift_deg") && (br.angle = _scf_num(meta_e["neutral_shift_deg"], "extra.meta.neutral_shift_deg"))
    end
    # the MATPOWER branch metadata record of the direct import (stage 3a):
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
    if haskey(be_all, "transformer_loss_g_pu") && !isempty(net.trafos)
      net.trafos[end].side1.g = _scf_num(be_all["transformer_loss_g_pu"], "extra.transformer_loss_g_pu")
    end
    # DTF branch metadata (stage 3b): the full importer record, tagged
    # encoding, order preserving; export and report layers read it
    if haskey(be_all, "dtf_branch")
      net.matpower_branch_metadata[length(net.branchVec)] = scf_decode_value(be_all["dtf_branch"])
    end
    # the typed phase-tap model of a DTF or CGMES winding (stage 3b)
    if haskey(be_all, "phase_taps") && entry.kind == "generic_branch" && !isempty(net.trafos)
      net.trafos[end].side1.phase_taps = scf_decode_value(be_all["phase_taps"])
    end
  end
  for row in sort(data.link; by = r -> r.id)
    addLink!(net = net, fromBus = bus_by_id[row.from_node], toBus = bus_by_id[row.to_node], status = something(row.from_status, 1))
  end

  # --- tap-changer nameplates (shared application path) -----------------
  comps = spar.components
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
    psi = 0.0
    psi_given = nothing
    dustep = 0.0
    for ctrl in row.controllers
      idx = ctrl.index
      if idx == 1
        tstep = something(ctrl.step, 0.0)
        tmin = Float64(something(ctrl.pos_min, 0.0))
        tmax = Float64(something(ctrl.pos_max, 0.0))
        tcur = Float64(something(ctrl.pos, 0.0))
      elseif idx == 2
        psi = something(ctrl.alpha_deg, 0.0)
        ctrl.alpha_deg === nothing || (psi_given = psi)
        pmin = Float64(something(ctrl.pos_min, 0.0))
        pmax = Float64(something(ctrl.pos_max, 0.0))
        pcur = Float64(something(ctrl.pos, 0.0))
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
    widened = tstep > 0.0 && !(tmin <= 0.0 <= tmax)
    if widened
      @warn "SCF: the neutral tap position lies outside the regulation band of this transformer; the band is restored as stated in the file" branch = row.branch pos_min = tmin pos_max = tmax
    end
    applyTapNameplate!(br; tap_step = tstep, tap_min_step = widened ? min(tmin, 0.0) : tmin, tap_max_step = widened ? max(tmax, 0.0) : tmax, tap_current_step = tcur, phase_step_deg = pstep, phase_min_step = pmin, phase_max_step = pmax, phase_current_step = pcur, psi_deg = psi, phase_du_step = dustep, context = "sparlectra.components.tap_changer for branch $(row.branch)")
    if widened
      br.tap_min = br.ratio / (1.0 + tmax * tstep)
      br.tap_max = br.ratio / (1.0 + tmin * tstep)
    end
    # the nameplate reads psi = 0 as "not given" and substitutes the 90 deg
    # symmetric-PST convention; a file that states 0 means 0
    psi_given === nothing || (br.tap_est_alpha_deg = psi_given)
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
      # P(U)/Q(U) controllers (stage 3a); the explicit flag keeps every
      # existing case file's behavior untouched
      pu_ctrl = nothing
      qu_ctrl = nothing
      if _scf_get(e, "pq_gen_controller", false) === true
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
        isRegulated = isempty(e) ? reg !== nothing : _scf_get(e, "regulated", false) === true,
        # refreshing the bus types per appliance is O(prosumers) EACH TIME, so
        # a case with n appliances pays n^2 (measured: 6.9 s for case13659,
        # 0.5 s with the deferral). The MATPOWER importer defers the same way
        # and refreshes once when the network is complete.
        defer_bus_type_refresh = true,
        participationFactor = get(participation, id, nothing),
      )
      # ratings are not constructor keywords (the PV/APU flags are derived
      # from the bus type and the voltage setpoint); restore them on the
      # created prosumer so the round trip keeps the nameplate
      ps = net.prosumpsVec[end]
      haskey(e, "rated_s_mva") && (ps.ratedS = scf_number_or_sentinel(e["rated_s_mva"], "extra.rated_s_mva"))
      haskey(e, "rated_u_kv") && (ps.ratedU = scf_number_or_sentinel(e["rated_u_kv"], "extra.rated_u_kv"))
      # the APU flag is derived from the bus type AT CONSTRUCTION time, and the
      # reader sets the final bus types first, so it must come from the file
      _scf_get(e, "apu_node", false) === true && (ps.isAPUNode = true)
      # addProsumer! auto-regulates anything that carries a voltage setpoint,
      # so restoring the setpoint of an UNREGULATED machine would silently
      # promote it to a PV generator. The file knows which it was.
      isempty(e) || (ps.isRegulated = _scf_get(e, "regulated", false) === true)
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
    g_pu = something(row.g1, 0.0) * zb_sh
    b_pu = something(row.b1, 0.0) * zb_sh
    # addShunt! takes the nominal P/Q draw at 1 pu voltage in MW/MVar
    # addShunt! forms y = (P + jQ)/baseMVA, so Q carries the susceptance with
    # the SAME sign the file states; negating it turned a capacitor into a
    # reactor on every round trip
    addShunt!(net = net, busName = bus_by_id[node], pShunt = g_pu * s_base, qShunt = b_pu * s_base)
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
    sh = net.shuntVec[shunt_index_by_id[row.shunt]]
    row.status === nothing || (sh.status = row.status)
    if row.model !== nothing
      model = Symbol(row.model)
      model in (:Y, :VoltageDependentInjection, :PQ) || throw(ArgumentError("SCF: shunt_state.model must be Y, VoltageDependentInjection, or PQ; got $(repr(row.model))."))
      sh.model = model
    end
    row.estimate === nothing || (sh.estimate = row.estimate === true)
  end

  # --- feeder data a PGM `source` declares --------------------------------
  # PGM models a source as a voltage source BEHIND its impedance: `sk` and
  # `rx_ratio` are part of the model, not decoration. They arrive here as
  # declared case data in the same shape addExternalGrid! produces, so the
  # short circuit sees the feeder and `power_flow.external_grid` (source:
  # auto) computes with the file's own numbers instead of config defaults.
  # Our own files carry the richer components.sc_source block; that one wins.
  declared_buses = Set(String(f.bus) for f in net.sc_sources.external_network_injections if f.bus !== nothing)
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

  # --- controllers (FACTS and regulation) --------------------------------
  # the entries are the declarative control.controllers schema, so they go
  # through applyConfiguredControllers! unchanged. Transformer references
  # may use the file's HUMAN branch names (extra block); the built net
  # regenerates component names, so those references are translated to
  # branch indices first (demo-case task: a self-built case says
  # trafo = "Moorau_110_20", not a generated name).
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
  _scf_read_measurements!(net, case, bus_by_id, branch_index_by_id, u_rated_by_id, names)

  # --- CGMES source identity (stage 3c) -----------------------------------
  # the structural-key mRID registry travels in the namespaced meta; with
  # it restored, the CGMES exporter and the mRID-addressed measurement
  # resolution work on a case-built net exactly as on the direct import
  ids_meta = _scf_get(meta, "cgmes_ids", nothing)
  if ids_meta isa AbstractDict
    for (k, v) in ids_meta
      net.cgmes_ids[String(k)] = String(v)
    end
  end

  # --- isolated nodes (roles.isolated_nodes, stage 3a) --------------------
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

  # --- MATPOWER dcline records (stage 3a) ---------------------------------
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

  # --- FOR001 contingency labels (stage 3a/3b) ----------------------------
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

# PGM sensors back into Sparlectra measurement rows; measurements.rows
# supplies the original row ids, weights, and active flags
function _scf_read_measurements!(net::Net, case::SCFCase, bus_by_id, branch_index_by_id, u_rated_by_id, names)
  data = case.data
  spar = case.sparlectra === nothing ? SCFSparlectra() : case.sparlectra
  rows = _scf_get(spar.measurements, "rows", [])
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
  for row in data.sym_voltage_sensor
    id = row.id
    node = row.measured_object
    ubase = u_rated_by_id[node]
    idx = net.busDict[bus_by_id[node]]
    i = 1
    if _scf_has(row.u_measured)
      push!(collected, (pos_of(id, i), Measurement(typ = VmMeas, value = value_of(id, i, row.u_measured / ubase), sigma = sigma_of(id, i, (_scf_has(row.u_sigma) ? row.u_sigma : 0.01 * ubase) / ubase), active = active(id, i), busIdx = idx, id = sid(id, i, "Vm_$(node)"))))
      i += 1
    end
    if _scf_has(row.u_angle_measured)
      push!(collected, (pos_of(id, i), Measurement(typ = VaMeas, value = value_of(id, i, rad2deg(row.u_angle_measured)), sigma = sigma_of(id, i, 0.02), active = active(id, i), busIdx = idx, id = sid(id, i, "Va_$(node)"))))
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
        push!(collected, (pos_of(id, i), Measurement(typ = typ, value = value, sigma = sigma_of(id, i, sigma), active = active(id, i), busIdx = net.busDict[bus_by_id[obj]], id = sid(id, i, "$(typ)_$(obj)"))))
      else
        push!(collected, (pos_of(id, i), Measurement(typ = typ, value = value, sigma = sigma_of(id, i, sigma), active = active(id, i), branchIdx = branch_index_by_id[obj], direction = term == "branch_to" ? :to : :from, id = sid(id, i, "$(typ)_$(obj)"))))
      end
      i += 1
    end
  end
  for row in data.sym_current_sensor
    id = row.id
    obj = row.measured_object
    term = _scf_has(row.measured_terminal_type) ? row.measured_terminal_type : "branch_from"
    dir = term == "branch_to" ? :to : :from
    i = 1
    if _scf_has(row.i_measured)
      push!(collected, (pos_of(id, i), Measurement(typ = ImagMeas, value = value_of(id, i, row.i_measured / 1.0e3), sigma = sigma_of(id, i, (_scf_has(row.i_sigma) ? row.i_sigma : 1.0e3) / 1.0e3), active = active(id, i), branchIdx = branch_index_by_id[obj], direction = dir, id = sid(id, i, "Imag_$(obj)"))))
      i += 1
    end
    if _scf_has(row.i_angle_measured)
      push!(collected, (pos_of(id, i), Measurement(typ = IaMeas, value = value_of(id, i, rad2deg(row.i_angle_measured)), sigma = sigma_of(id, i, rad2deg(_scf_has(row.i_angle_sigma) ? row.i_angle_sigma : deg2rad(0.05))), active = active(id, i), branchIdx = branch_index_by_id[obj], direction = dir, id = sid(id, i, "Ia_$(obj)"))))
    end
  end
  sort!(collected; by = e -> e[1])
  for (_, m) in collected
    push!(net.measurements, m)
  end
  return nothing
end

"""
    build_net(case::SCFCase; config = active_sparlectra_config()) -> Net

The one network constructor of the run path over the typed case (design
decisions D2 and D12): construct through the existing public
constructors, then apply the run configuration's net parameters exactly
once, here. `config` is the EFFECTIVE run configuration (the run path
passes its resolved `ImportedCase` configuration; interactive use takes
the active one).
"""
function build_net(case::SCFCase; config::SparlectraConfig = active_sparlectra_config())::Net
  net = _scf_net_from_case(case)
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

The study definitions a case file carries: `(contingencies, short_circuit)`,
each an empty dictionary when the file does not define that study. The
blocks state WHAT to compute; the result contract of the runs themselves is
unchanged.
"""
function scf_case_studies(file::AbstractString)
  root = scf_json_parse(read(String(file), String))
  spar_rev = _scf_get(root, "sparlectra", Dict{String,Any}())
  isempty(spar_rev) || _scf_require_current_revision(spar_rev)
  _scf_validate_studies(root)
  spar = _scf_get(root, "sparlectra", Dict{String,Any}())
  cont = _scf_get(spar, "contingencies", Dict{String,Any}())
  sc = _scf_get(spar, "short_circuit", Dict{String,Any}())
  return (contingencies = cont isa AbstractDict ? Dict{String,Any}(cont) : Dict{String,Any}(), short_circuit = sc isa AbstractDict ? Dict{String,Any}(sc) : Dict{String,Any}())
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

The declared unit system of the electrical parameters
(task_scf_units_v0100): `:pu` or `:si`. An ABSENT declaration means
`:si`, which is what every file before the declaration meant, so
existing files read unchanged. Unknown values are a hard error. For
`:pu` the conversion base must be complete (D2): `meta.s_base` (VA,
positive and finite) and a positive finite `u_rated` on every node;
per-unit numbers without their base are not interpretable, so a
missing piece is a hard error naming the field.
"""
function scf_case_units(case::SCFCase)::Symbol
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

"""
    scf_case_config(file) -> Dict{String,Any}

The dotted configuration keys a case file carries in
`sparlectra.config`, ready to be merged as `config_overrides` (the case
file sits below API/CLI overrides and above the YAML file in precedence).
Empty when the file carries no configuration.
"""
function scf_case_config(file::AbstractString)::Dict{String,Any}
  root = scf_json_parse(read(String(file), String))
  spar = _scf_get(root, "sparlectra", Dict{String,Any}())
  # consulted BEFORE the import on the run path, so it must refuse a stale
  # revision itself or the run would first act on old configuration keys
  isempty(spar) || _scf_require_current_revision(spar)
  cfg = _scf_get(spar, "config", Dict{String,Any}())
  out = Dict{String,Any}(String(k) => v for (k, v) in cfg)
  # A key outside the case scope is refused, not ignored: a file that states a
  # setting which quietly does not apply is exactly the silent acceptance this
  # format exists to prevent. Older files that carried the whole form need one
  # re-export.
  bad = sort!(String[k for k in keys(out) if !scf_is_case_config_key(k)])
  isempty(bad) || throw(ArgumentError("SCF: config key(s) $(join(bad, ", ")) are not case scope (logging, benchmarking, parallelism, Web UI and export settings belong in the configuration file). Re-export the case with format revision $(SCF_FORMAT_VERSION), or remove them."))
  return out
end
