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

# file: src/adapters/scf/scf_export.jl
# purpose: Sparlectra Case Format (SCF) writer: serialize a Net into the
#          PGM-compatible case file (issue #342). The `data` section is a
#          valid power-grid-model input dataset in SI units; everything
#          Sparlectra adds lives under the namespaced `sparlectra` key.
#          Ids are assigned deterministically (component type, then name),
#          so rebuilding the same net produces byte-identical output.

# Format revision. Bumped whenever the written file changes shape, so a reader
# can refuse a file it does not understand instead of misreading it. 1.1: an
# ineffective voltage setpoint is no longer written, an unlimited reactive
# limit is expressed by omission, and fields equal to their documented default
# are left out (see the default catalog in docs/src/scf.md).
const SCF_FORMAT_VERSION = "1.0"

# Configuration groups a case file may NOT carry: they describe the machine
# and the session, not the network. Logging and result-table shape, benchmark
# runs, parallelism, and the MATPOWER export are decided
# by the installation's configuration, so a delivered case cannot change how
# someone else's system reports, measures, or parallelizes.
const _SCF_NON_CASE_CONFIG_PREFIXES = ("output.", "benchmark.", "runtime.", "webui.", "matpower_export.")

"""
    scf_is_case_config_key(key) -> Bool

Whether a dotted configuration key describes the CASE (its model, import
conventions, solver, estimation, and short-circuit math) rather than the
installation. Only case-scope keys travel inside a case file; everything else
stays in the configuration file. The key must also be GUI-editable, which is
the same allowlist every other configuration surface uses.
"""
function scf_is_case_config_key(key::AbstractString)::Bool
  k = String(key)
  any(p -> startswith(k, p), _SCF_NON_CASE_CONFIG_PREFIXES) && return false
  return k in GUI_EDITABLE_CONFIG_KEYS
end
const SCF_PGM_VERSION = "1.0"

# id assignment order; a component's id is its position in this walk, so the
# file is stable against Dict iteration order and against rebuilds
const _SCF_ID_TYPE_ORDER = ("node", "line", "generic_branch", "link", "source", "sym_load", "sym_gen", "shunt", "voltage_regulator", "sym_voltage_sensor", "sym_power_sensor", "sym_current_sensor", "tap_changer", "transformer3w")

## --- unit helpers ----------------------------------------------------------
## Sparlectra works in per-unit on net.baseMVA with kV bus bases; SCF is SI
## (V, W, var, VA, ohm, siemens, farad, A, rad) like PGM.

_scf_zbase(vn_kV::Float64, baseMVA::Float64) = (vn_kV * 1.0e3)^2 / (baseMVA * 1.0e6)
_scf_u_rated(vn_kV::Float64) = vn_kV * 1.0e3

"""
    _scf_branch_vn(net, branch) -> Float64

Voltage base of a branch's impedance in kV. Sparlectra stamps the tap on the
FROM side (`ui /= tap` in the Y-bus), so `r_pu`/`x_pu` are referenced to the
TO side, which is exactly PGM's `generic_branch` convention.
"""
function _scf_branch_vn(net::Net, branch::Branch)::Float64
  ti = Int(branch.toBus)
  (1 <= ti <= length(net.nodeVec)) && return getNodeVn(net.nodeVec[ti])
  return branch.comp.cVN
end

## --- deterministic ids -----------------------------------------------------

struct ScfIdMap
  node::Dict{Int,Int}          # bus index -> id
  branch::Dict{Int,Int}        # branch index -> id
  link::Dict{Int,Int}          # link index -> id
  prosumer::Dict{Int,Int}      # prosumer position -> id
  shunt::Dict{Int,Int}         # shunt position -> id
  regulator::Dict{Int,Int}     # prosumer position -> voltage_regulator id
  tap_changer::Dict{Int,Int}   # branch index -> tap_changer id
  transformer3w::Dict{Int,Int} # star bus index -> transformer3w group id
  names::Dict{Int,String}      # id -> reference name
end

_scf_prosumer_bus(ps::ProSumer)::Int = something(ps.comp.cFrom_bus, 0)

# a branch is a transformer when it carries a turns ratio; that is the same
# marker the solver and the auto power-flow feature extraction use
_scf_is_transformer(br::Branch)::Bool = br.ratio != 0.0

"""
    _scf_bus_name(net, i) -> String

The REFERENCE name of a bus: the `busDict` key. That is the name
measurements, controllers, contingencies, and reports resolve against, and
the name the reader has to restore; the internal component name ("Bus_1_110",
"Aux_12_380") is kept separately in `extra.component_name` when it differs.
"""
function _scf_bus_name(net::Net, i::Int)::String
  for (name, idx) in net.busDict
    idx == i && return String(name)
  end
  return getCompName(net.nodeVec[i].comp)
end

"""
    scf_id_map(net) -> ScfIdMap

Assign SCF ids deterministically: components are walked in a fixed type
order and, inside a type, in lexicographic name order (ties broken by the
internal index). Rebuilding the same net therefore produces the same ids and
the same file bytes.
"""
function scf_id_map(net::Net)::ScfIdMap
  next = 0
  nid = () -> (next += 1)
  names = Dict{Int,String}()
  nodes = Dict{Int,Int}()
  for i in sort(collect(eachindex(net.nodeVec)); by = i -> (_scf_bus_name(net, i), i))
    id = nid()
    nodes[i] = id
    names[id] = _scf_bus_name(net, i)
  end
  lines = Dict{Int,Int}()
  trafos = Dict{Int,Int}()
  lineidx = [i for i in eachindex(net.branchVec) if !_scf_is_transformer(net.branchVec[i])]
  trafoidx = [i for i in eachindex(net.branchVec) if _scf_is_transformer(net.branchVec[i])]
  for i in sort(lineidx; by = i -> (getCompName(net.branchVec[i].comp), i))
    id = nid()
    lines[i] = id
    names[id] = getCompName(net.branchVec[i].comp)
  end
  for i in sort(trafoidx; by = i -> (getCompName(net.branchVec[i].comp), i))
    id = nid()
    trafos[i] = id
    names[id] = getCompName(net.branchVec[i].comp)
  end
  branch = merge(lines, trafos)
  links = Dict{Int,Int}()
  for i in sort(collect(eachindex(net.linkVec)); by = i -> (net.linkVec[i].cName, i))
    id = nid()
    links[i] = id
    names[id] = net.linkVec[i].cName
  end
  # prosumers split into sources (external network injection at the slack),
  # generators, and loads; the type walk keeps ids stable per class
  srcidx = [i for i in eachindex(net.prosumpsVec) if net.prosumpsVec[i].comp.cTyp == ExternalNetworkInjection]
  genidx = [i for i in eachindex(net.prosumpsVec) if isGenerator(net.prosumpsVec[i]) && !(i in srcidx)]
  loadidx = [i for i in eachindex(net.prosumpsVec) if !isGenerator(net.prosumpsVec[i]) && !(i in srcidx)]
  prosumer = Dict{Int,Int}()
  for group in (srcidx, loadidx, genidx)
    for i in sort(group; by = i -> (getCompName(net.prosumpsVec[i].comp), i))
      id = nid()
      prosumer[i] = id
      names[id] = getCompName(net.prosumpsVec[i].comp)
    end
  end
  shunts = Dict{Int,Int}()
  for i in sort(collect(eachindex(net.shuntVec)); by = i -> (getCompName(net.shuntVec[i].comp), i))
    id = nid()
    shunts[i] = id
    names[id] = getCompName(net.shuntVec[i].comp)
  end
  # voltage regulators: one per PV-controlling generator
  regulators = Dict{Int,Int}()
  for i in sort(genidx; by = i -> (getCompName(net.prosumpsVec[i].comp), i))
    ps = net.prosumpsVec[i]
    _scf_is_pv_generator(net, ps) || continue
    id = nid()
    regulators[i] = id
    names[id] = string("VR_", getCompName(ps.comp))
  end
  # tap changers: one per branch with a declared changer
  taps = Dict{Int,Int}()
  for i in sort(trafoidx; by = i -> (getCompName(net.branchVec[i].comp), i))
    br = net.branchVec[i]
    (br.has_ratio_tap || br.has_phase_tap) || continue
    id = nid()
    taps[i] = id
    names[id] = string("TC_", getCompName(br.comp))
  end
  # three-winding groups: one id per star node (an auxiliary node with
  # exactly three incident branches, i.e. the star point of the equivalent)
  t3w = Dict{Int,Int}()
  for i in sort(collect(eachindex(net.nodeVec)); by = i -> (getCompName(net.nodeVec[i].comp), i))
    _scf_is_aux_node(net, i) || continue
    count(k -> Int(net.branchVec[k].fromBus) == i || Int(net.branchVec[k].toBus) == i, eachindex(net.branchVec)) == 3 || continue
    id = nid()
    t3w[i] = id
    names[id] = string("T3W_", _scf_bus_name(net, i))
  end
  return ScfIdMap(nodes, branch, links, prosumer, shunts, regulators, taps, t3w, names)
end

# a generator counts as PV-controlling when its bus is a PV node
function _scf_is_pv_generator(net::Net, ps::ProSumer)::Bool
  bus = _scf_prosumer_bus(ps)
  (1 <= bus <= length(net.nodeVec)) || return false
  return getNodeType(net.nodeVec[bus]) == PV
end

## --- data section (the PGM subset) -----------------------------------------

function _scf_nodes(net::Net, ids::ScfIdMap)
  rows = Vector{Any}()
  for i in eachindex(net.nodeVec)
    nd = net.nodeVec[i]
    push!(rows, Dict{String,Any}("id" => ids.node[i], "u_rated" => _scf_u_rated(getNodeVn(nd))))
  end
  sort!(rows; by = r -> r["id"])
  return rows
end

"""
    _scf_stable(pu, to_si, to_pu) -> Float64

Return the SI number to WRITE for the per-unit value `pu`: the one that the
reader turns back into exactly `pu`. `to_si` is this writer's conversion,
`to_pu` the reader's.

Unit conversion is not exactly invertible in floating point, and the naive
value cost both guarantees at once: the re-read network differed in the last
bit (so the Y-bus was not bit-identical), and the next export wrote that
neighbour, so a file could drift on EVERY cycle instead of being a fixed
point (seen on case14/case118/case300). Checking the two neighbouring
floats finds an exact preimage in practice; if none exists, the value that
maps to itself keeps at least the byte-identical round trip. The correction
is at most one ulp.
"""
function _scf_stable(pu::Float64, to_si, to_pu)::Float64
  si = to_si(pu)
  to_pu(si) == pu && return si
  for candidate in (nextfloat(si), prevfloat(si))
    to_pu(candidate) == pu && return candidate
  end
  for _ in 1:3
    next = to_si(to_pu(si))
    next == si && return si
    si = next
  end
  # No exact preimage anywhere near: then the value must at least be a FIXED
  # POINT of write-read-write, or the file drifts by one ulp on every cycle
  # (measured on case1354pegase, transformer b1). The per-unit value can move
  # by one ulp here, which is the documented tolerance; the file cannot.
  for candidate in (si, nextfloat(si), prevfloat(si))
    to_si(to_pu(candidate)) == candidate && return candidate
  end
  return si
end

"""
    _scf_stable_step(pos, base, step) -> Int

The band position to WRITE. The reader rebuilds the ratio edge as
`base / (1 + pos * step)`, and a second export derives the position from that
edge again. Rounding can land on the other side of a tie there, so the file
would not be a fixed point (measured on case13659pegase: 6674 branches moved
their band edge by one step on re-export). Pick the neighbour that reproduces
itself; the ratio edge stays within the documented one-step snap either way.
"""
function _scf_stable_step(pos::Int, base::Float64, step::Float64)::Int
  step > 0.0 || return pos
  reproduces(p) = begin
    edge = base / (1.0 + p * step)
    edge == 0.0 && return false
    return Int(round((base / edge - 1.0) / step, RoundNearestTiesAway)) == p
  end
  reproduces(pos) && return pos
  for candidate in (pos - 1, pos + 1)
    reproduces(candidate) && return candidate
  end
  return pos
end

function _scf_lines(net::Net, ids::ScfIdMap, f_nom::Float64)
  rows = Vector{Any}()
  w_nom = 2.0 * pi * f_nom
  for i in eachindex(net.branchVec)
    br = net.branchVec[i]
    _scf_is_transformer(br) && continue
    zb = _scf_zbase(_scf_branch_vn(net, br), net.baseMVA)
    # divide by the base rather than multiply by its reciprocal: one rounding
    # less, and the written value is a fixed point of the read/write cycle
    b_s = br.b_pu / zb
    g_s = br.g_pu / zb
    row = Dict{String,Any}(
      "id" => ids.branch[i],
      "from_node" => ids.node[Int(br.fromBus)],
      "to_node" => ids.node[Int(br.toBus)],
      "from_status" => br.from_status,
      "to_status" => br.to_status,
      # base (physical) impedance, never a FACTS-compensated operating point;
      # each value is the fixed point of the reader's conversion (see _scf_stable)
      "r1" => _scf_stable(br.r_base_pu, v -> v * zb, v -> v / zb),
      "x1" => _scf_stable(br.x_base_pu, v -> v * zb, v -> v / zb),
      "c1" => _scf_stable(br.b_pu, v -> (v / zb) / w_nom, v -> (v * w_nom) * zb),
      "tan1" => b_s == 0.0 ? 0.0 : g_s / b_s,
    )
    # tan1 = g/b cannot carry a conductance without capacitance; the
    # direct field keeps it (the reader prefers g1 when present)
    b_s == 0.0 && g_s != 0.0 && (row["g1"] = _scf_stable(br.g_pu, v -> v / zb, v -> v * zb))
    # an unlimited rating (MATPOWER rateA = 0 arrives as Inf) has no PGM
    # representation: no i_n IS the statement "no limit"
    if br.sn_MVA !== nothing && br.sn_MVA > 0.0 && isfinite(br.sn_MVA)
      row["i_n"] = br.sn_MVA * 1.0e6 / (sqrt(3.0) * _scf_u_rated(_scf_branch_vn(net, br)))
    end
    push!(rows, row)
  end
  sort!(rows; by = r -> r["id"])
  return rows
end

function _scf_generic_branches(net::Net, ids::ScfIdMap)
  rows = Vector{Any}()
  for i in eachindex(net.branchVec)
    br = net.branchVec[i]
    _scf_is_transformer(br) || continue
    zb = _scf_zbase(_scf_branch_vn(net, br), net.baseMVA)
    # LIVE tap position: k and theta describe the branch as it computes,
    # the nameplate detail lives in sparlectra.components.tap_changer
    t = calcBranchRatio(br)
    row = Dict{String,Any}(
      "id" => ids.branch[i],
      "from_node" => ids.node[Int(br.fromBus)],
      "to_node" => ids.node[Int(br.toBus)],
      "from_status" => br.from_status,
      "to_status" => br.to_status,
      "r1" => _scf_stable(br.r_base_pu, v -> v * zb, v -> v / zb),
      "x1" => _scf_stable(br.x_base_pu, v -> v * zb, v -> v / zb),
      "g1" => _scf_stable(br.g_pu, v -> v / zb, v -> v * zb),
      "b1" => _scf_stable(br.b_pu, v -> v / zb, v -> v * zb),
      "k" => abs(t),
      "theta" => angle(t),
    )
    # same as i_n on a line: an unlimited rating is expressed by the ABSENCE
    # of sn here, and the namespaced meta keeps the distinction
    (br.sn_MVA === nothing || !isfinite(br.sn_MVA)) || (row["sn"] = br.sn_MVA * 1.0e6)
    push!(rows, row)
  end
  sort!(rows; by = r -> r["id"])
  return rows
end

function _scf_links(net::Net, ids::ScfIdMap)
  rows = Vector{Any}()
  for i in eachindex(net.linkVec)
    lk = net.linkVec[i]
    push!(rows, Dict{String,Any}("id" => ids.link[i], "from_node" => ids.node[lk.fromBus], "to_node" => ids.node[lk.toBus], "from_status" => lk.status, "to_status" => lk.status))
  end
  sort!(rows; by = r -> r["id"])
  return rows
end

"""
    _scf_cgmes_mrids(net) -> (buses, branches, links)

Resolve the CGMES mRIDs of a delivery-sourced net by their structural keys
(`TN|<bus>`, `ACL|<a>|<b>|<k>`, `PT|<a>|<b>|<k>`, `SW|...`), so `external_id`
carries the source-system identity instead of the internal component id.
Empty dictionaries on nets from other formats.
"""
function _scf_cgmes_mrids(net::Net)
  buses = Dict{Int,String}()
  branches = Dict{Int,String}()
  links = Dict{Int,String}()
  isempty(net.cgmes_ids) && return buses, branches, links
  nby = Dict{Int,String}(idx => name for (name, idx) in net.busDict)
  for i in eachindex(net.nodeVec)
    name = get(nby, i, "")
    isempty(name) && continue
    m = get(net.cgmes_ids, CGMESImporter.cgmesKeyTopologicalNode(name), "")
    isempty(m) || (buses[i] = m)
  end
  # parallel index per bus pair, counted the same way the importer did
  cline = Dict{Tuple{String,String},Int}()
  ctrafo = Dict{Tuple{String,String},Int}()
  for (k, br) in enumerate(net.branchVec)
    a = get(nby, Int(br.fromBus), "")
    b = get(nby, Int(br.toBus), "")
    (isempty(a) || isempty(b)) && continue
    if _scf_is_transformer(br)
      kk = CGMESImporter.cgmesNextParallelIndex!(ctrafo, a, b)
      m = get(net.cgmes_ids, CGMESImporter.cgmesKeyPowerTransformer(a, b, kk), "")
    else
      kk = CGMESImporter.cgmesNextParallelIndex!(cline, a, b)
      m = get(net.cgmes_ids, CGMESImporter.cgmesKeyACLineSegment(a, b, kk), "")
    end
    isempty(m) || (branches[k] = m)
  end
  return buses, branches, links
end

"""
    _scf_sc_sources(net, ids) -> Vector

Short-circuit source data (IEC 60909 feeders and machines) as a
`sparlectra.components.sc_source` list. The data is stored on the net by the
CGMES importer and by `addExternalGrid!`; PGM's `source` component carries
only `sk`/`rx_ratio`, so the full set travels here and stays lossless.
"""
function _scf_sc_sources(net::Net, ids::ScfIdMap)
  rows = Vector{Any}()
  sc = net.sc_sources
  sc === nothing && return rows
  nby = Dict{String,Int}(name => idx for (name, idx) in net.busDict)
  for (kind, entries) in (("external_network_injection", sc.external_network_injections), ("synchronous_machine", sc.synchronous_machines), ("asynchronous_machine", sc.asynchronous_machines), ("equivalent_injection", sc.equivalent_injections))
    for e in entries
      row = Dict{String,Any}("kind" => kind)
      # empty fields are written as null rather than dropped: the consumer
      # reads a source record by field, so a record that comes back short a
      # field would fail the run. Writing null keeps the record SHAPE in the
      # file, and the field vocabulary stays in one place (the record itself).
      for key in propertynames(e)
        v = getproperty(e, key)
        v isa AbstractDict && continue
        row[String(key)] = v isa Symbol ? String(v) : scf_infinity_sentinel(v)
      end
      # resolve the bus name onto the file's node id where possible
      bname = get(row, "bus", nothing)
      if bname isa AbstractString
        bidx = get(nby, String(bname), nothing)
        bidx === nothing || (row["node"] = ids.node[bidx])
      end
      push!(rows, row)
    end
  end
  sort!(rows; by = r -> (String(r["kind"]), string(get(r, "node", 0)), string(get(r, "mrid", ""))))
  return rows
end

# The feeder record of a bus, as the (sk in VA, R/X) pair a PGM `source`
# states. Sparlectra keeps this in the short-circuit source data, which both
# the CGMES import and addExternalGrid! fill.
function _scf_feeder_for_bus(net::Net, busIdx::Int)
  sc = net.sc_sources
  sc === nothing && return nothing
  names = Dict{Int,String}(idx => n for (n, idx) in net.busDict)
  bus = get(names, busIdx, nothing)
  bus === nothing && return nothing
  vn = getNodeVn(net.nodeVec[busIdx])
  for f in sc.external_network_injections
    (f.bus !== nothing && String(f.bus) == bus) || continue
    ik = f.maxInitialSymShCCurrent_A
    (ik === nothing || !isfinite(Float64(ik)) || Float64(ik) <= 0.0) && continue
    rx = f.maxR1ToX1Ratio
    z01 = f.maxZ0ToZ1Ratio
    return (sk_va = sqrt(3.0) * vn * 1.0e3 * Float64(ik), rx = (rx === nothing || !isfinite(Float64(rx))) ? 0.1 : Float64(rx), z01 = z01 === nothing ? nothing : Float64(z01))
  end
  return nothing
end

function _scf_appliances(net::Net, ids::ScfIdMap)
  sources = Vector{Any}()
  loads = Vector{Any}()
  gens = Vector{Any}()
  regs = Vector{Any}()
  for i in eachindex(net.prosumpsVec)
    ps = net.prosumpsVec[i]
    bus = _scf_prosumer_bus(ps)
    id = ids.prosumer[i]
    node = get(ids.node, bus, nothing)
    node === nothing && continue
    p_w = (ps.pVal === nothing ? 0.0 : ps.pVal) * 1.0e6
    q_var = (ps.qVal === nothing ? 0.0 : ps.qVal) * 1.0e6
    # PGM's `source` IS the external network injection, so the COMPONENT type
    # decides (proSumptionType only says Injection or Consumption)
    if ps.comp.cTyp == ExternalNetworkInjection
      row = Dict{String,Any}("id" => id, "node" => node, "status" => 1, "u_ref" => ps.vm_pu === nothing ? 1.0 : ps.vm_pu)
      ps.va_deg === nothing || (row["u_ref_angle"] = deg2rad(ps.va_deg))
      # PGM's source is a voltage source BEHIND its impedance, so a complete
      # source states sk and rx_ratio. The values are the feeder's own
      # short-circuit data, which is where Sparlectra keeps them.
      feeder = _scf_feeder_for_bus(net, bus)
      if feeder !== nothing
        row["sk"] = feeder.sk_va
        row["rx_ratio"] = feeder.rx
        feeder.z01 === nothing || (row["z01_ratio"] = feeder.z01)
      end
      push!(sources, row)
    elseif isGenerator(ps)
      push!(gens, Dict{String,Any}("id" => id, "node" => node, "status" => 1, "type" => 0, "p_specified" => p_w, "q_specified" => q_var))
      rid = get(ids.regulator, i, nothing)
      if rid !== nothing
        row = Dict{String,Any}("id" => rid, "regulated_object" => id, "status" => 1, "u_ref" => ps.vm_pu === nothing ? 1.0 : ps.vm_pu)
        # An unlimited Q band is normal (a generator without Q limits, common
        # in the pegase cases) and arrives as +-Inf. JSON cannot carry that,
        # and the dataset says "no limit" by ABSENCE, which is exactly how the
        # reader takes a missing key. Writing it would abort the export.
        (ps.minQ === nothing || !isfinite(ps.minQ)) || (row["q_min"] = ps.minQ * 1.0e6)
        (ps.maxQ === nothing || !isfinite(ps.maxQ)) || (row["q_max"] = ps.maxQ * 1.0e6)
        push!(regs, row)
      end
    else
      push!(loads, Dict{String,Any}("id" => id, "node" => node, "status" => 1, "type" => 0, "p_specified" => p_w, "q_specified" => q_var))
    end
  end
  for v in (sources, loads, gens, regs)
    sort!(v; by = r -> r["id"])
  end
  return sources, loads, gens, regs
end

"""
    _scf_shunt_state(net, ids) -> Vector

Per-shunt state that PGM's `shunt` (g1/b1 only) has no slot for: whether the
shunt is in service, which model it uses (a voltage-dependent injection is
NOT an admittance), and whether its susceptance is released as a
state-estimation state. Only shunts that deviate from the defaults produce a
row, so an ordinary case file stays free of noise.
"""
function _scf_shunt_state(net::Net, ids::ScfIdMap)
  rows = Vector{Any}()
  for i in eachindex(net.shuntVec)
    sh = net.shuntVec[i]
    row = Dict{String,Any}("shunt" => ids.shunt[i])
    sh.status == 1 || (row["status"] = sh.status)
    sh.model === :Y || (row["model"] = String(Symbol(sh.model)))
    sh.estimate && (row["estimate"] = true)
    length(row) == 1 && continue
    push!(rows, row)
  end
  sort!(rows; by = r -> r["shunt"])
  return rows
end

"""
    _scf_faults(net, ids, buses) -> Vector

The buses a short-circuit study evaluates, in PGM's own `fault` vocabulary
(`fault_object`, `fault_type`, `r_f`/`x_f`). Sparlectra computes the
balanced three-phase bolted fault, so the rows state exactly that; the block
is only written when a study names buses, and it is the PGM-readable twin of
`short_circuit.buses`.
"""
function _scf_faults(net::Net, ids::ScfIdMap, buses)
  rows = Vector{Any}()
  known = Set(values(ids.node))
  for b in buses
    # a study block addresses buses by NODE ID (the file's own vocabulary); a
    # bus name is accepted too and resolved through the network
    node = if b isa Integer
      Int(b) in known || throw(ArgumentError("exportSCF: fault bus id $(b) is not a node of this network."))
      Int(b)
    else
      idx = get(net.busDict, String(b), nothing)
      idx === nothing && throw(ArgumentError("exportSCF: unknown fault bus $(repr(b))."))
      ids.node[idx]
    end
    push!(rows, Dict{String,Any}("id" => 0, "fault_object" => node, "status" => 1, "fault_type" => "three_phase", "r_f" => 0.0, "x_f" => 0.0))
  end
  return rows
end

function _scf_shunts(net::Net, ids::ScfIdMap)
  rows = Vector{Any}()
  for i in eachindex(net.shuntVec)
    sh = net.shuntVec[i]
    zb_sh = _scf_zbase(sh.vn_kV, net.baseMVA)
    # The reader's arithmetic is not `si * zb`: it forms the shunt draw in
    # MW/MVar and `addShunt!` divides by baseMVA again. Multiplying and
    # dividing by the same base is not the identity in floating point, so
    # stabilizing against the short formula optimized the wrong function and
    # the file drifted by one ulp per cycle (case1354pegase shunts).
    to_pu = v -> (v * zb_sh * net.baseMVA) / net.baseMVA
    push!(rows, Dict{String,Any}("id" => ids.shunt[i], "node" => ids.node[sh.busIdx], "status" => 1, "g1" => _scf_stable(real(sh.y_pu_shunt), v -> v / zb_sh, to_pu), "b1" => _scf_stable(imag(sh.y_pu_shunt), v -> v / zb_sh, to_pu)))
  end
  sort!(rows; by = r -> r["id"])
  return rows
end

## --- measurements ----------------------------------------------------------
## PGM sensors carry value and sigma; the Sparlectra row detail that PGM has
## no slot for (id string, weight, active flag, unpaired P/Q rows) is listed
## under sparlectra.measurements so the round trip stays lossless.

# what the reader reconstructs from the SI value of a sensor; used to decide
# whether the exact internal value has to be recorded alongside it
function _scf_roundtrip_value(m, net::Net)::Float64
  if m.typ == VmMeas
    bus = m.busIdx
    bus === nothing && return m.value
    # the reader rebuilds its voltage base from the file's u_rated
    # (u_rated / 1e3 * 1e3), which is a different rounding path than ours
    u = _scf_u_rated(getNodeVn(net.nodeVec[bus]))
    return (m.value * u) / u
  elseif m.typ == VaMeas || m.typ == IaMeas
    return rad2deg(deg2rad(m.value))
  elseif m.typ == ImagMeas
    return (m.value * 1.0e3) / 1.0e3
  else
    return (m.value * 1.0e6) / 1.0e6
  end
end

function _scf_sensors(net::Net, ids::ScfIdMap, next_id::Ref{Int})
  volt = Vector{Any}()
  power = Vector{Any}()
  current = Vector{Any}()
  rows = Vector{Any}()
  meas = net.measurements
  isempty(meas) && return volt, power, current, rows
  # position of each row in the measurement vector: the WLS objective sums in
  # that order, so restoring it is what makes an SE round trip bit-identical
  pos_by_id = Dict{String,Int}(m.id => i for (i, m) in enumerate(meas))
  newid = () -> (next_id[] += 1)
  # Group by measured location so P/Q and Vm/Va become one PGM sensor. A
  # location can carry the quantity MORE than once (redundant transducers are
  # normal in telemetry), so every group holds a list and the n-th reading at
  # a location becomes its own PGM sensor - PGM allows several sensors on one
  # measured object. Keying by location alone silently dropped every reading
  # but the last.
  vm = Dict{Int,Vector{Measurement}}()
  va = Dict{Int,Vector{Measurement}}()
  pinj = Dict{Int,Vector{Measurement}}()
  qinj = Dict{Int,Vector{Measurement}}()
  pflow = Dict{Tuple{Int,Symbol},Vector{Measurement}}()
  qflow = Dict{Tuple{Int,Symbol},Vector{Measurement}}()
  imag_ = Dict{Tuple{Int,Symbol},Vector{Measurement}}()
  iang = Dict{Tuple{Int,Symbol},Vector{Measurement}}()
  other = Measurement[]
  add!(d, key, m) = push!(get!(() -> Measurement[], d, key), m)
  # readings at one location, in the order the measurement vector holds them
  nocc(d, key) = length(get(d, key, ()))
  at(d, key, i) = (v = get(d, key, nothing); v === nothing || i > length(v) ? nothing : v[i])
  for m in meas
    if m.typ == VmMeas && m.busIdx !== nothing
      add!(vm, m.busIdx, m)
    elseif m.typ == VaMeas && m.busIdx !== nothing
      add!(va, m.busIdx, m)
    elseif m.typ == PinjMeas && m.busIdx !== nothing
      add!(pinj, m.busIdx, m)
    elseif m.typ == QinjMeas && m.busIdx !== nothing
      add!(qinj, m.busIdx, m)
    elseif m.typ == PflowMeas && m.branchIdx !== nothing
      add!(pflow, (m.branchIdx, m.direction), m)
    elseif m.typ == QflowMeas && m.branchIdx !== nothing
      add!(qflow, (m.branchIdx, m.direction), m)
    elseif m.typ == ImagMeas && m.branchIdx !== nothing
      add!(imag_, (m.branchIdx, m.direction), m)
    elseif m.typ == IaMeas && m.branchIdx !== nothing
      add!(iang, (m.branchIdx, m.direction), m)
    else
      push!(other, m)
    end
  end
  # One PGM sensor can stand for two Sparlectra rows (P and Q, Vm and Va),
  # and PGM carries ONE sigma per sensor in SI. The authoritative internal
  # sigmas (per-unit for Vm, degrees for Va, MW/MVar for powers, kA for
  # currents) therefore travel here, so a round trip is exact; the weight is
  # derived from sigma and is never written.
  rowrec(sensor_id, ms...) = begin
    entry = Dict{String,Any}("sensor" => sensor_id, "sparlectra_ids" => String[m.id for m in ms], "sigmas" => Float64[m.sigma for m in ms])
    all(m -> m.active, ms) || (entry["active"] = Bool[m.active for m in ms])
    # The PGM sensor stays the authority for the measured value. The SI
    # conversion (MW to W, per-unit to V, degrees to radians) is not always
    # bit-reversible, so a value is recorded here ONLY when reading it back
    # would not reproduce the exact internal number; that keeps the round
    # trip exact without duplicating values that convert cleanly.
    vals = Float64[m.value for m in ms]
    rts = Float64[_scf_roundtrip_value(m, net) for m in ms]
    vals == rts || (entry["values"] = vals)
    entry["positions"] = Int[get(pos_by_id, m.id, 0) for m in ms]
    push!(rows, entry)
    return nothing
  end
  # voltage sensors (Vm always, Va when the bus carries a PMU angle)
  for bus in sort!(collect(keys(vm)))
    node = get(ids.node, bus, nothing)
    node === nothing && continue
    u_base = _scf_u_rated(getNodeVn(net.nodeVec[bus]))
    for k = 1:nocc(vm, bus)
      m = vm[bus][k]
      sid = newid()
      row = Dict{String,Any}("id" => sid, "measured_object" => node, "u_measured" => m.value * u_base, "u_sigma" => m.sigma * u_base)
      a = at(va, bus, k)
      if a === nothing
        rowrec(sid, m)
      else
        row["u_angle_measured"] = deg2rad(a.value)
        rowrec(sid, m, a)
      end
      push!(volt, row)
    end
  end
  for bus in sort!(collect(keys(va)))
    node = get(ids.node, bus, nothing)
    node === nothing && continue
    # angle readings that no magnitude sensor at this bus took along
    for k = (nocc(vm, bus) + 1):nocc(va, bus)
      m = va[bus][k]
      sid = newid()
      push!(volt, Dict{String,Any}("id" => sid, "measured_object" => node, "u_measured" => nothing, "u_angle_measured" => deg2rad(m.value), "u_sigma" => nothing))
      rowrec(sid, m)
    end
  end
  # injection power sensors
  for bus in sort!(unique(vcat(collect(keys(pinj)), collect(keys(qinj)))))
    node = get(ids.node, bus, nothing)
    node === nothing && continue
    for k = 1:max(nocc(pinj, bus), nocc(qinj, bus))
      p = at(pinj, bus, k)
      q = at(qinj, bus, k)
      sid = newid()
      row = Dict{String,Any}(
        "id" => sid,
        "measured_object" => node,
        "measured_terminal_type" => "node",
        "p_measured" => p === nothing ? nothing : p.value * 1.0e6,
        "q_measured" => q === nothing ? nothing : q.value * 1.0e6,
        "power_sigma" => (p === nothing ? q.sigma : p.sigma) * 1.0e6,
      )
      # PGM carries one sigma per sensor; when P and Q differ it also accepts
      # the individual ones, so a PGM reader sees the same weighting we do
      if p !== nothing && q !== nothing && p.sigma != q.sigma
        row["p_sigma"] = p.sigma * 1.0e6
        row["q_sigma"] = q.sigma * 1.0e6
      end
      push!(power, row)
      rowrec(sid, filter(!isnothing, [p, q])...)
    end
  end
  # branch flow power sensors
  for key in sort!(unique(vcat(collect(keys(pflow)), collect(keys(qflow)))); by = k -> (k[1], String(k[2])))
    bidx, dir = key
    bid = get(ids.branch, bidx, nothing)
    bid === nothing && continue
    for k = 1:max(nocc(pflow, key), nocc(qflow, key))
      p = at(pflow, key, k)
      q = at(qflow, key, k)
      sid = newid()
      frow = Dict{String,Any}(
        "id" => sid,
        "measured_object" => bid,
        "measured_terminal_type" => dir == :to ? "branch_to" : "branch_from",
        "p_measured" => p === nothing ? nothing : p.value * 1.0e6,
        "q_measured" => q === nothing ? nothing : q.value * 1.0e6,
        "power_sigma" => (p === nothing ? q.sigma : p.sigma) * 1.0e6,
      )
      if p !== nothing && q !== nothing && p.sigma != q.sigma
        frow["p_sigma"] = p.sigma * 1.0e6
        frow["q_sigma"] = q.sigma * 1.0e6
      end
      push!(power, frow)
      rowrec(sid, filter(!isnothing, [p, q])...)
    end
  end
  # current sensors (magnitude in kA internally, angle in degrees)
  for key in sort!(unique(vcat(collect(keys(imag_)), collect(keys(iang)))); by = k -> (k[1], String(k[2])))
    bidx, dir = key
    bid = get(ids.branch, bidx, nothing)
    bid === nothing && continue
    for k = 1:max(nocc(imag_, key), nocc(iang, key))
      im_ = at(imag_, key, k)
      ia = at(iang, key, k)
      sid = newid()
      row = Dict{String,Any}(
        "id" => sid,
        "measured_object" => bid,
        "measured_terminal_type" => dir == :to ? "branch_to" : "branch_from",
        "angle_measurement_type" => "global_angle",
        "i_measured" => im_ === nothing ? nothing : im_.value * 1.0e3,
        "i_sigma" => im_ === nothing ? nothing : im_.sigma * 1.0e3,
      )
      if ia !== nothing
        row["i_angle_measured"] = deg2rad(ia.value)
        row["i_angle_sigma"] = deg2rad(ia.sigma)
      end
      push!(current, row)
      rowrec(sid, filter(!isnothing, [im_, ia])...)
    end
  end
  for v in (volt, power, current)
    sort!(v; by = r -> r["id"])
  end
  sort!(rows; by = r -> r["sensor"])
  return volt, power, current, rows
end

## --- sparlectra section ----------------------------------------------------

function _scf_roles(net::Net, ids::ScfIdMap)
  slack_nodes = Int[ids.node[i] for i in eachindex(net.nodeVec) if getNodeType(net.nodeVec[i]) == Slack]
  sort!(slack_nodes)
  participation = Vector{Any}()
  for i in eachindex(net.prosumpsVec)
    ps = net.prosumpsVec[i]
    ps.participationFactor === nothing && continue
    ps.participationFactor == 0.0 && continue
    push!(participation, Dict{String,Any}("object" => ids.prosumer[i], "factor" => Float64(ps.participationFactor)))
  end
  sort!(participation; by = r -> r["object"])
  roles = Dict{String,Any}()
  if !isempty(participation)
    # `normalize` is the documented default (true) and therefore not written
    roles["slack"] = Dict{String,Any}("mode" => "distributed", "participation" => participation, "nodes" => slack_nodes)
  elseif length(slack_nodes) == 1
    roles["slack"] = Dict{String,Any}("mode" => "single", "nodes" => slack_nodes)
  elseif isempty(slack_nodes)
    roles["slack"] = Dict{String,Any}("mode" => "source")
  else
    # several ideal slacks (multi-island cases): still "single" per island,
    # listed explicitly so the reader restores every one of them
    roles["slack"] = Dict{String,Any}("mode" => "single", "nodes" => slack_nodes)
  end
  # auxiliary nodes: the explicit AuxBus type plus the star points the
  # three-winding equivalent creates (the CGMES importer names them
  # AUX3WT_<transformer> and leaves their component type at busbar)
  aux = Int[ids.node[i] for i in eachindex(net.nodeVec) if _scf_is_aux_node(net, i)]
  isempty(aux) || (roles["aux_nodes"] = sort!(aux))
  iso = Int[ids.node[b] for b in net.isoNodes if haskey(ids.node, b)]
  isempty(iso) || (roles["isolated_nodes"] = sort!(iso))
  return roles
end

const _SCF_STAR_NODE_PREFIX = "AUX3WT_"

# Documented defaults of the namespaced block (see docs/src/scf.md): a field
# whose value equals its default is not written, and the reader restores it.
const _SCF_DEFAULT_BRANCH_KIND = "BranchC"
const _SCF_DEFAULT_TAP_EST_MODE = "none"

# a node is auxiliary when it carries the AuxBus type or is a three-winding
# star point (the importer's naming convention)
function _scf_is_aux_node(net::Net, i::Int)::Bool
  nd = net.nodeVec[i]
  nd.comp.cTyp == AuxBus && return true
  # the three-winding star point: the importer builds it with isAux (giving
  # the component name an "Aux_" prefix) and names the bus AUX3WT_<trafo>
  startswith(getCompName(nd.comp), "Aux_") && return true
  return startswith(_scf_bus_name(net, i), _SCF_STAR_NODE_PREFIX)
end

"""
    _scf_transformer3w(net, ids) -> Vector

Group the three star legs of a three-winding transformer. Sparlectra builds
the same star equivalent PGM documents for `generic_branch`: three legs and
one auxiliary star node. The electrical parameters therefore live exactly
once, in `data.generic_branch`; this block only says which three branches
belong together, which node is their star point, and which end is which.

Leg direction is Sparlectra's own: `from_node` is the star node, `to_node`
the terminal, and each end names both explicitly, so no reader has to guess
an orientation.
"""
function _scf_transformer3w(net::Net, ids::ScfIdMap)
  rows = Vector{Any}()
  for i in eachindex(net.nodeVec)
    haskey(ids.transformer3w, i) || continue
    legs = [k for k in eachindex(net.branchVec) if Int(net.branchVec[k].fromBus) == i || Int(net.branchVec[k].toBus) == i]
    length(legs) == 3 || continue
    ends = Vector{Any}()
    for k in sort(legs; by = k -> -_scf_branch_vn(net, net.branchVec[k]))
      br = net.branchVec[k]
      terminal = Int(br.fromBus) == i ? Int(br.toBus) : Int(br.fromBus)
      entry = Dict{String,Any}("branch" => ids.branch[k], "terminal_node" => ids.node[terminal], "u_rated" => _scf_u_rated(getNodeVn(net.nodeVec[terminal])))
      (br.sn_MVA === nothing || !isfinite(br.sn_MVA)) || (entry["sn"] = br.sn_MVA * 1.0e6)
      push!(ends, entry)
    end
    # role labels follow the terminal voltage order (highest first); they
    # are reporting labels, not an ordering constraint on the model
    for (j, role) in enumerate(("hv", "mv", "lv"))
      ends[j]["role"] = role
    end
    row = Dict{String,Any}("id" => ids.transformer3w[i], "star_node" => ids.node[i], "ends" => ends, "leg_direction" => "star_to_terminal")
    tc = [get(ids.tap_changer, k, nothing) for k in legs]
    filter!(!isnothing, tc)
    isempty(tc) || (row["tap_changer"] = first(tc))
    push!(rows, row)
  end
  sort!(rows; by = r -> r["star_node"])
  return rows
end

function _scf_tap_changers(net::Net, ids::ScfIdMap)
  rows = Vector{Any}()
  for i in eachindex(net.branchVec)
    id = get(ids.tap_changer, i, nothing)
    id === nothing && continue
    br = net.branchVec[i]
    controllers = Vector{Any}()
    if br.has_ratio_tap && br.tap_step > 0.0
      # position on the fraction grid: r1 = base/current - 1 over tap_step
      cur = br.tap_ratio == 0.0 ? br.ratio : br.tap_ratio
      pos = cur == 0.0 ? 0 : Int(round((br.ratio / cur - 1.0) / br.tap_step, RoundNearestTiesAway))
      push!(controllers, Dict{String,Any}(
        "index" => 1,
        "alpha_deg" => 0.0,
        "step" => br.tap_step,
        "pos" => pos,
        "pos_min" => _scf_stable_step(Int(round((br.ratio / br.tap_max - 1.0) / br.tap_step, RoundNearestTiesAway)), br.ratio, br.tap_step),
        "pos_max" => _scf_stable_step(Int(round((br.ratio / br.tap_min - 1.0) / br.tap_step, RoundNearestTiesAway)), br.ratio, br.tap_step),
      ))
    end
    if br.has_phase_tap
      if br.phase_du_step > 0.0
        # Delta-u phase shifter: the position lives on the additional-voltage
        # grid, the shift angle is a consequence (atan), not the grid
        psi = deg2rad(br.tap_est_alpha_deg)
        tbase = (br.ratio == 0.0 ? 1.0 : br.ratio) * cis(deg2rad(br.angle))
        tlive = (br.tap_ratio == 0.0 ? abs(tbase) : br.tap_ratio) * cis(deg2rad(br.phase_shift_deg))
        r2 = real((tbase / tlive - 1.0) * cis(-psi))
        push!(controllers, Dict{String,Any}(
          "index" => 2,
          "alpha_deg" => br.tap_est_alpha_deg,
          "step" => br.phase_du_step,
          "pos" => Int(round(r2 / br.phase_du_step, RoundNearestTiesAway)),
          "pos_min" => Int(round(br.phase_du_min_step, RoundNearestTiesAway)),
          "pos_max" => Int(round(br.phase_du_max_step, RoundNearestTiesAway)),
        ))
      elseif br.phase_step_deg > 0.0
        push!(controllers, Dict{String,Any}(
          "index" => 2,
          "alpha_deg" => br.tap_est_alpha_deg,
          "step_deg" => br.phase_step_deg,
          # The BAND is stored relative to the neutral already
          # (`phase_min_deg = pos_min * step`, see applyTapNameplate!), the
          # live shift is not (`phase_shift_deg = angle + pos * step`).
          # Subtracting the neutral shift from the band as well moved it by
          # angle/step on EVERY write: the file drifted, and with a case whose
          # neutral shift is not zero the band stopped bracketing the neutral,
          # which made the file unreadable (case13659pegase in its own rad /
          # -1.0 import convention).
          "pos" => Int(round((br.phase_shift_deg - br.angle) / br.phase_step_deg, RoundNearestTiesAway)),
          "pos_min" => Int(round(br.phase_min_deg / br.phase_step_deg, RoundNearestTiesAway)),
          "pos_max" => Int(round(br.phase_max_deg / br.phase_step_deg, RoundNearestTiesAway)),
        ))
      end
    end
    isempty(controllers) && continue
    row = Dict{String,Any}(
      "id" => id,
      "branch" => ids.branch[i],
      "side" => "from",
      "ratio_base" => br.ratio,
      "angle_base_deg" => br.angle,
      "controllers" => controllers,
      "control" => Dict{String,Any}("mode" => "fixed"),
    )
    String(br.tap_est_mode) == _SCF_DEFAULT_TAP_EST_MODE || (row["tap_est_mode"] = String(br.tap_est_mode))
    push!(rows, row)
  end
  sort!(rows; by = r -> r["id"])
  return rows
end

# The MATPOWER converter turns a PQ generator's limits into two constant
# controllers (points at 0.0 and 2.0 pu, the value being the machine's own
# setpoint). That pair is written back as the flag it came from, so every
# existing case file stays byte for byte; anything else is a real
# characteristic and is written as one. The check is EXACT on those two
# voltages: a constant two-point controller a user set deliberately at other
# voltages is a characteristic and stays one.
function _scf_matpower_constant(ch::PiecewiseLinearCharacteristic, setpoint_pu::Float64)::Bool
  ch.interpolation === :linear && length(ch.points) == 2 || return false
  (u1, y1), (u2, y2) = ch.points
  return u1 == 0.0 && u2 == 2.0 && y1 == y2 && isapprox(y1, setpoint_pu; rtol = 1.0e-12, atol = 1.0e-15)
end

function _scf_characteristic_dict(ch::PiecewiseLinearCharacteristic, lo, hi, lo_key::AbstractString, hi_key::AbstractString, s_base::Float64)
  d = Dict{String,Any}("points" => [[u, y] for (u, y) in ch.points])
  # the constructor's default is not written, like every other default
  ch.interpolation === :linear || (d["interpolation"] = String(ch.interpolation))
  # limits in MW / MVAr like max_q_mvar; an unlimited side is stated by absence
  (lo === nothing || !isfinite(lo)) || (d[lo_key] = lo * s_base)
  (hi === nothing || !isfinite(hi)) || (d[hi_key] = hi * s_base)
  return d
end

"Voltage-dependent control of a machine into its `extra` entry (see the note on the MATPOWER constant form above)."
function _scf_write_voltage_control!(d::Dict{String,Any}, ps::ProSumer, s_base::Float64)
  qu = ps.quController
  pu = ps.puController
  (qu === nothing && pu === nothing) && return d
  if qu !== nothing && pu !== nothing &&
     _scf_matpower_constant(qu.characteristic, something(ps.qVal, 0.0) / s_base) &&
     _scf_matpower_constant(pu.characteristic, something(ps.pVal, 0.0) / s_base)
    d["pq_gen_controller"] = true
    return d
  end
  qu === nothing || (d["qu_control"] = _scf_characteristic_dict(qu.characteristic, qu.qmin_pu, qu.qmax_pu, "qmin_mvar", "qmax_mvar", s_base))
  pu === nothing || (d["pu_control"] = _scf_characteristic_dict(pu.characteristic, pu.pmin_pu, pu.pmax_pu, "pmin_mw", "pmax_mw", s_base))
  return d
end

function _scf_extra(net::Net, ids::ScfIdMap)
  extra = Dict{String,Any}()
  put(id, d) = (extra[string(id)] = d)
  # CGMES deliveries: external_id must be the delivery's mRID, not the
  # internal component id, or the file loses its link to the source system
  mrid_bus, mrid_branch, _ = _scf_cgmes_mrids(net)
  for i in eachindex(net.nodeVec)
    nd = net.nodeVec[i]
    refname = _scf_bus_name(net, i)
    d = Dict{String,Any}("name" => refname, "external_id" => get(mrid_bus, i, nd.comp.cID), "bus_index" => i, "node_type" => String(Symbol(getNodeType(nd))))
    refname == getCompName(nd.comp) || (d["component_name"] = getCompName(nd.comp))
    # the SOURCE system's bus number (MATPOWER bus id, DTF node number). It
    # is what the generated component names embed, so losing it renamed every
    # component of a re-read case and changed the file's ids.
    haskey(net.busOrigIdxDict, i) && (d["source_index"] = net.busOrigIdxDict[i])
    haskey(mrid_bus, i) && (d["source_id_kind"] = "cgmes_mrid")
    nd._area === nothing || (d["area"] = nd._area)
    nd._lZone === nothing || (d["zone"] = nd._lZone)
    # a bus inherits the network's limits; only its own deviation is written
    (nd._vmin_pu === nothing || nd._vmin_pu == net.vmin_pu) || (d["vmin_pu"] = nd._vmin_pu)
    (nd._vmax_pu === nothing || nd._vmax_pu == net.vmax_pu) || (d["vmax_pu"] = nd._vmax_pu)
    nd.comp.cTyp == AuxBus && (d["aux"] = true)
    put(ids.node[i], d)
  end
  for i in eachindex(net.branchVec)
    br = net.branchVec[i]
    d = Dict{String,Any}("name" => getCompName(br.comp), "external_id" => get(mrid_branch, i, br.comp.cID), "branch_index" => i)
    # informational, and only when it differs from the documented default
    kindname = String(Symbol(br.comp.cTyp))
    kindname == _SCF_DEFAULT_BRANCH_KIND || (d["kind"] = kindname)
    haskey(mrid_branch, i) && (d["source_id_kind"] = "cgmes_mrid")
    # imported nameplate data Sparlectra stores but does not interpret
    # travels verbatim, so a round trip keeps what the source system said
    meta = Dict{String,Any}()
    br.sn_MVA === nothing || (meta["sn_mva"] = scf_infinity_sentinel(br.sn_MVA))
    br.ratio == 0.0 || (meta["neutral_ratio"] = br.ratio)
    br.angle == 0.0 || (meta["neutral_shift_deg"] = br.angle)
    isempty(meta) || (d["meta"] = meta)
    put(ids.branch[i], d)
  end
  for i in eachindex(net.linkVec)
    lk = net.linkVec[i]
    put(ids.link[i], Dict{String,Any}("name" => lk.cName, "external_id" => lk.cID, "link_index" => i))
  end
  for i in eachindex(net.prosumpsVec)
    ps = net.prosumpsVec[i]
    # the COMPONENT type is what a constructor accepts back ("Generator",
    # "EnergyConsumer", ...); the coarse Injection/Consumption class is
    # derived from it and would not round-trip
    d = Dict{String,Any}("name" => getCompName(ps.comp), "external_id" => ps.comp.cID, "prosumption_type" => String(Symbol(ps.comp.cTyp)), "prosumption_class" => String(Symbol(ps.proSumptionType)), "element_index" => i)
    ps.ratedS === nothing || (d["rated_s_mva"] = scf_infinity_sentinel(ps.ratedS))
    ps.ratedU === nothing || (d["rated_u_kv"] = scf_infinity_sentinel(ps.ratedU))
    ps.maxP === nothing || (d["max_p_mw"] = scf_infinity_sentinel(ps.maxP))
    ps.minP === nothing || (d["min_p_mw"] = scf_infinity_sentinel(ps.minP))
    # The reactive band travels here as well, not only in a voltage_regulator
    # row: a machine written as a PGM `source` has no regulator row, and PGM
    # has no place for its Q limits. Without this a regulating generator at
    # its Q limit came back unlimited and a Q-limit run took a different path
    # (measured on case300).
    # An unlimited band is stated by ABSENCE, in the namespaced block exactly
    # as in the dataset: one rule, no sentinel to interpret. Verified that a
    # missing and an infinite limit produce the same qmin_pu/qmax_pu and the
    # same power flow.
    (ps.maxQ === nothing || !isfinite(ps.maxQ)) || (d["max_q_mvar"] = ps.maxQ)
    (ps.minQ === nothing || !isfinite(ps.minQ)) || (d["min_q_mvar"] = ps.minQ)
    # Same reason for the voltage setpoint: a machine without a
    # voltage_regulator row has nowhere else to keep it, and losing it moves
    # the solved voltage at that bus (measured on case300: 0.05 pu).
    # A setpoint is written only where it ACTS: the machine regulates and has
    # no voltage_regulator row of its own to keep it in. An unregulated
    # machine carries the constructor's 1.0, a value with no effect on the
    # calculation, and writing it made every such machine come back as a
    # regulated one (case1354pegase: 260 regulated machines became 933).
    # Slack sources are unaffected: their setpoint is the `source.u_ref`.
    (ps.vm_pu === nothing || !ps.isRegulated || haskey(ids.regulator, i)) || (d["vm_pu"] = ps.vm_pu)
    ps.isAPUNode && (d["apu_node"] = true)
    # Which appliance carries the reference is a property of the machine, not
    # of its bus: deriving it from "sits on a slack node" gave a LOAD at that
    # bus the reference as well (case57).
    ps.referencePri === nothing || (d["reference_pri"] = true)
    ps.isRegulated && (d["regulated"] = true)
    _scf_write_voltage_control!(d, ps, net.baseMVA)
    put(ids.prosumer[i], d)
  end
  for i in eachindex(net.shuntVec)
    sh = net.shuntVec[i]
    put(ids.shunt[i], Dict{String,Any}("name" => getCompName(sh.comp), "external_id" => sh.comp.cID, "shunt_index" => i))
  end
  return extra
end

## --- assembly --------------------------------------------------------------

"""
    net_to_scf(net; kwargs...) -> Dict{String,Any}

Build the SCF root object of `net` (see [`exportSCF`](@ref) for the keyword
arguments). Pure: the network is not modified.
"""
function net_to_scf(
  net::Net;
  case_name::AbstractString = "",
  f_nom::Float64 = 50.0,
  source_format::AbstractString = "",
  source_reference::AbstractString = "",
  intended_calculations::Vector{String} = String["power_flow"],
  notes::AbstractString = "",
  contingencies::Union{Nothing,AbstractDict} = nothing,
  scenarios::Union{Nothing,AbstractDict} = nothing,
  short_circuit::Union{Nothing,AbstractDict} = nothing,
  include_start_state::Bool = false,
  measurement_provenance::AbstractDict = Dict{String,Any}(),
  strict_pgm::Bool = false,
  created_at::AbstractString = "",
)::Dict{String,Any}
  _scf_assert_base_impedances(net)
  ids = scf_id_map(net)
  next_id = Ref(isempty(ids.names) ? 0 : maximum(keys(ids.names)))
  data = Dict{String,Any}()
  nodes = _scf_nodes(net, ids)
  isempty(nodes) || (data["node"] = nodes)
  lines = _scf_lines(net, ids, f_nom)
  isempty(lines) || (data["line"] = lines)
  gbs = _scf_generic_branches(net, ids)
  isempty(gbs) || (data["generic_branch"] = gbs)
  links = _scf_links(net, ids)
  isempty(links) || (data["link"] = links)
  sources, loads, gens, regs = _scf_appliances(net, ids)
  isempty(sources) || (data["source"] = sources)
  isempty(loads) || (data["sym_load"] = loads)
  isempty(gens) || (data["sym_gen"] = gens)
  isempty(regs) || (data["voltage_regulator"] = regs)
  shunts = _scf_shunts(net, ids)
  isempty(shunts) || (data["shunt"] = shunts)
  volt, power, current, meas_rows = _scf_sensors(net, ids, next_id)
  isempty(volt) || (data["sym_voltage_sensor"] = volt)
  isempty(power) || (data["sym_power_sensor"] = power)
  isempty(current) || (data["sym_current_sensor"] = current)

  meta = Dict{String,Any}(
    "case_name" => isempty(case_name) ? net.name : String(case_name),
    "s_base" => net.baseMVA * 1.0e6,
    "f_nom" => f_nom,
    "created_by" => string("Sparlectra ", version()),
    "intended_calculations" => intended_calculations,
  )
  # Voltage limits are part of the case, not of the reader: pegase cases
  # declare 0.8/1.2 per bus, and losing them moved the solved voltages of half
  # the network (case13659pegase). Written only when they differ from the
  # documented default 0.9/1.1.
  net.vmin_pu == 0.9 || (meta["vmin_pu"] = net.vmin_pu)
  net.vmax_pu == 1.1 || (meta["vmax_pu"] = net.vmax_pu)
  isempty(source_format) || (meta["source_format"] = String(source_format))
  isempty(source_reference) || (meta["source_reference"] = String(source_reference))
  isempty(notes) || (meta["notes"] = String(notes))
  isempty(created_at) || (meta["created_at"] = String(created_at))

  spar = Dict{String,Any}("format_version" => SCF_FORMAT_VERSION, "meta" => meta, "roles" => _scf_roles(net, ids), "extra" => _scf_extra(net, ids))
  components = Dict{String,Any}()
  taps = _scf_tap_changers(net, ids)
  isempty(taps) || (components["tap_changer"] = taps)
  # IEC 60909 source data (CGMES deliveries, addExternalGrid!): PGM's source
  # carries only sk/rx_ratio, the full record travels here
  sc = _scf_sc_sources(net, ids)
  isempty(sc) || (components["sc_source"] = sc)
  # three-winding grouping: the electrical values stay in the three legs,
  # this only says which legs belong together and which end is which
  t3 = _scf_transformer3w(net, ids)
  isempty(t3) || (components["transformer3w"] = t3)
  # FACTS and regulation: PGM has no controller model, so they live here in
  # the declarative control.controllers schema, verbatim
  ctrls = _scf_controllers(net)
  isempty(ctrls) || (components["controllers"] = ctrls)
  shunt_state = _scf_shunt_state(net, ids)
  isempty(shunt_state) || (components["shunt_state"] = shunt_state)
  # PGM's own fault vocabulary for the buses a short-circuit study names
  if short_circuit !== nothing && haskey(short_circuit, "buses")
    faults = _scf_faults(net, ids, short_circuit["buses"])
    for row in faults
      next_id[] += 1
      row["id"] = next_id[]
    end
    isempty(faults) || (data["fault"] = faults)
  end
  isempty(components) || (spar["components"] = components)
  if !isempty(meas_rows)
    block = Dict{String,Any}("rows" => meas_rows)
    # where the values came from, and above all WHETHER THEY CARRY NOISE: a
    # noise-free set has J = 0 by construction, which reads like a perfect
    # estimate and is nothing of the sort. The case file says so instead of
    # leaving that to be guessed from a suspiciously small J.
    isempty(measurement_provenance) || (block["provenance"] = Dict{String,Any}(String(k) => v for (k, v) in measurement_provenance))
    spar["measurements"] = block
  end
  # The writer emits no `sparlectra.config` block any more: a case's
  # settings live in its case configuration file (`<stem>.config.yaml`,
  # write_case_config). The reader still accepts the block as the
  # deprecated precedence level directly below that file (D7).
  # study definitions: what to compute, not what came out. The result
  # contract is untouched; these blocks only describe the study.
  contingencies === nothing || (spar["contingencies"] = Dict{String,Any}(String(k) => v for (k, v) in contingencies))
  # the scenario block (scenario task D3): the scenario-aware path writes
  # scenarios; the deprecated contingencies keyword keeps its historical
  # emission for existing callers
  scenarios === nothing || (spar["scenarios"] = Dict{String,Any}(String(k) => v for (k, v) in scenarios))
  short_circuit === nothing || (spar["short_circuit"] = Dict{String,Any}(String(k) => v for (k, v) in short_circuit))
  if include_start_state
    nodes_state = Dict{String,Any}()
    for i in eachindex(net.nodeVec)
      nd = net.nodeVec[i]
      (nd._vm_pu === nothing || nd._va_deg === nothing) && continue
      nodes_state[string(ids.node[i])] = Dict{String,Any}("vm_pu" => nd._vm_pu, "va_deg" => nd._va_deg)
    end
    isempty(nodes_state) || (spar["start_state"] = Dict{String,Any}("source" => "solved_power_flow", "nodes" => nodes_state))
  end

  root = Dict{String,Any}("version" => SCF_PGM_VERSION, "type" => "input", "is_batch" => false, "attributes" => Dict{String,Any}(), "data" => data, "sparlectra" => spar)
  # Strict PGM: the plain input dataset, for a consumer that reads
  # power-grid-model and nothing else. What the namespaced block carries is
  # then genuinely gone from the file, so the dropped detail is NAMED rather
  # than silently omitted, and reading such a file back is a lossy import.
  if strict_pgm
    dropped = String[]
    haskey(spar, "roles") && push!(dropped, "slack roles and participation")
    haskey(spar, "extra") && push!(dropped, "component names and source ids")
    if haskey(spar, "components")
      comps = spar["components"]
      haskey(comps, "tap_changer") && push!(dropped, "tap-changer nameplates")
      haskey(comps, "transformer3w") && push!(dropped, "three-winding grouping")
      haskey(comps, "sc_source") && push!(dropped, "short-circuit source data")
      haskey(comps, "controllers") && push!(dropped, "controllers")
      haskey(comps, "shunt_state") && push!(dropped, "shunt state")
    end
    haskey(spar, "measurements") && push!(dropped, "measurement rows")
    any(ps -> has_qu_controller(ps) || has_pu_controller(ps), net.prosumpsVec) && push!(dropped, "voltage-dependent Q(U)/P(U) control")
    haskey(spar, "start_state") && push!(dropped, "start state")
    (haskey(spar, "contingencies") || haskey(spar, "short_circuit")) && push!(dropped, "study definitions")
    # g1 on a line row is the SCF extension for a conductance tan1 cannot
    # spell (b == 0); PGM's line has no such attribute, so a strict file
    # genuinely loses that conductance
    stripped_g1 = 0
    for r in get(data, "line", Any[])
      haskey(r, "g1") && (delete!(r, "g1"); stripped_g1 += 1)
    end
    stripped_g1 > 0 && push!(dropped, "$(stripped_g1) line shunt conductance(s) without capacitance (g1)")
    # PGM knows no slack FLAG: its reference is a `source`. A slack generator
    # is therefore rewritten as a source (and its sym_gen/voltage_regulator
    # rows drop out, or the node would inject twice) - otherwise the strict
    # file would have no reference at all and neither PGM nor this reader
    # could solve it.
    slack_nodes = Set(Int[Int(x) for x in get(get(spar["roles"], "slack", Dict{String,Any}()), "nodes", [])])
    if !isempty(slack_nodes)
      gens_rows = get(data, "sym_gen", Any[])
      regs_rows = get(data, "voltage_regulator", Any[])
      srcs = get(data, "source", Any[])
      converted = Any[]
      already = Set(Int(r["node"]) for r in srcs)
      for row in gens_rows
        Int(row["node"]) in slack_nodes || continue
        # a node that already carries a source keeps it: PGM's source IS
        # Sparlectra's external network injection, and adding a second
        # reference there would double the infeed
        Int(row["node"]) in already && continue
        reg = findfirst(r -> Int(r["regulated_object"]) == Int(row["id"]), regs_rows)
        u_ref = reg === nothing ? 1.0 : Float64(get(regs_rows[reg], "u_ref", 1.0))
        push!(converted, Dict{String,Any}("id" => row["id"], "node" => row["node"], "status" => row["status"], "u_ref" => u_ref))
      end
      if !isempty(converted)
        converted_ids = Set(Int(r["id"]) for r in converted)
        filter!(r -> !(Int(r["id"]) in converted_ids), gens_rows)
        filter!(r -> !(Int(r["regulated_object"]) in converted_ids), regs_rows)
        append!(srcs, converted)
        sort!(srcs; by = r -> r["id"])
        data["source"] = srcs
        isempty(gens_rows) ? delete!(data, "sym_gen") : (data["sym_gen"] = gens_rows)
        isempty(regs_rows) ? delete!(data, "voltage_regulator") : (data["voltage_regulator"] = regs_rows)
        push!(dropped, "$(length(converted)) slack generator(s) rewritten as PGM sources")
      end
    end
    isempty(dropped) || @warn "exportSCF(strict_pgm = true): writing the plain power-grid-model dataset; this file does NOT carry $(join(dropped, ", "))" case = String(case_name)
    delete!(root, "sparlectra")
  end
  return root
end

# A case file must describe the EQUIPMENT, not a FACTS operating point: a
# series-FACTS control run stamps its compensated impedance onto r_pu/x_pu
# while r_base_pu/x_base_pu keep the physical value. Writing a compensated
# net would silently misrepresent the network, so it is a hard error.
function _scf_assert_base_impedances(net::Net)
  for (i, br) in enumerate(net.branchVec)
    if !isapprox(br.r_pu, br.r_base_pu; atol = 1e-12, rtol = 1e-9) || !isapprox(br.x_pu, br.x_base_pu; atol = 1e-12, rtol = 1e-9)
      throw(ArgumentError("exportSCF: branch $(i) ($(getCompName(br.comp))) carries a FACTS-compensated impedance (r_pu/x_pu differ from r_base_pu/x_base_pu). A case file must describe the equipment, not an operating point; call restoreBaseImpedances!(net) before exporting."))
    end
  end
  return nothing
end

"""
    exportSCF(net; file, kwargs...) -> String

Write `net` as a Sparlectra Case Format file (`.scf.json`, issue #342) and
return the path. The `data` section is a valid power-grid-model input
dataset in SI units; everything Sparlectra adds (slack roles, tap-changer
cascade, names and source ids, measurement detail) lives under the
namespaced `sparlectra` key. Case-scope configuration travels in the case
configuration file next to the case (`write_case_config`), not inside it.

Ids are deterministic (component type, then name), so exporting the same
network twice produces byte-identical files and Git diffs stay readable.

Keyword arguments:
- `case_name`, `source_format`, `source_reference`, `notes`: `meta` fields;
  `case_name` defaults to the network name.
- `f_nom`: nominal frequency in Hz (default 50.0), needed for the line
  capacitance conversion.
- `intended_calculations`: the contract of the file (default
  `["power_flow"]`).
- `measurement_provenance`: what is known about the measurement set (`noise`,
  `generator`, `seed`); it lands in `sparlectra.measurements.provenance` and
  is what tells a reader that a set is ideal rather than measured.
- `include_start_state`: write the current bus voltages as START values
  (never as results). Off by default, because a stale start state on an
  edited net is worse than none.
- `strict_pgm`: write ONLY the power-grid-model input dataset, without the
  namespaced `sparlectra` block, for a consumer that reads PGM and nothing
  else. Everything that block carries is then absent from the file; the
  export names what it dropped in a warning. Off by default.

Throws an `ArgumentError` when the network carries FACTS-compensated
impedances (call `restoreBaseImpedances!` first).
"""
function exportSCF(net::Net; file::AbstractString, kwargs...)::String
  return write_scf_json(net_to_scfcase(net; kwargs...), String(file))
end

"""
    net_to_scfcase(net; kwargs...) -> SCFCase

The typed case of `net` (design decision D1): what [`exportSCF`](@ref)
writes, as the in-memory [`SCFCase`](@ref). Takes the same keyword
arguments as [`exportSCF`](@ref); pure, the network is not modified. The
document assembly stays the deterministic dict builder (`net_to_scf`),
converted once at the end, so the typed form and the file bytes cannot
diverge.
"""
net_to_scfcase(net::Net; kwargs...)::SCFCase = scfcase_from_root(net_to_scf(net; kwargs...))
