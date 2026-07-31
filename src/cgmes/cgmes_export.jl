# Copyright 2023–2026 Udo Schmitz
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

# file: src/cgmes/cgmes_export.jl
# purpose: CGMES 2.4.15 (RDF/XML) exporter for a Sparlectra Net — EQ, TP,
#          and SSH profiles covering buses, AC lines (incl. optional
#          zero-sequence attributes), two-winding transformers (fixed
#          effective ratio via end ratedU, phase shift as a single-step
#          linear phase tap changer), loads, machines, injections, and
#          shunts, plus the SSH operating point.
#
# Object identity: every exported object resolves its mRID through
# `_cgmesId(net, key)` with a structural key (see cgmes_keys.jl). For nets
# imported from CGMES the importer has recorded the original mRIDs under the
# same keys, so the export reuses them; for all other nets deterministic
# uuid5 ids are minted from the key string. Renaming a component therefore
# never changes its exported mRID — only structural changes do.
#
# Not produced (announced per unit via notices where applicable): static VAr
# compensators, bus links, tap-changer/controller machinery beyond the
# flattened ratio and shift, SV state, ZIP container.
#
# Accuracy note: class/attribute names follow CGMES 2.4.15 (cim16). Validate
# against the ENTSO-E test models before production use; individual profiles
# may require additional mandatory fields.

import UUIDs
import Dates

# Fixed namespace for deterministic mRIDs (uuid5). Same structural key
# => identical UUID on every export (roundtrip-/diff-stable).
const CGMES_UUID_NAMESPACE = UUIDs.UUID("6ba7b810-9dad-11d1-80b4-00c04fd430c8")

# RDF_NS / MD_NS / CIM16_NS come from cgmes_schema.jl (shared with the reader).
const EQ_PROFILES = [
  "http://entsoe.eu/CIM/EquipmentCore/3/1",
  "http://entsoe.eu/CIM/EquipmentShortCircuit/3/1",
  "http://entsoe.eu/CIM/EquipmentOperation/3/1",
]
const TP_PROFILES = ["http://entsoe.eu/CIM/Topology/4/1"]
const SSH_PROFILES = ["http://entsoe.eu/CIM/SteadyStateHypothesis/1/1"]
const SV_PROFILES = ["http://entsoe.eu/CIM/StateVariables/4/1"]

"""
    CGMESLineShortCircuit

Zero-sequence/short-circuit data for one line, for the EquipmentShortCircuit
profile. Physical units (Ohm, Siemens, °C). Only written when supplied — no
zero-sequence values are invented from positive-sequence data.
"""
Base.@kwdef struct CGMESLineShortCircuit
  r0_ohm::Float64
  x0_ohm::Float64
  b0ch_S::Float64 = 0.0
  g0ch_S::Float64 = 0.0
  endTemperature_C::Union{Nothing,Float64} = nothing
end

# ---------------------------------------------------------------------------
# Identity resolution
# ---------------------------------------------------------------------------

# Resolve the mRID for a structural key: reuse the recorded id (CGMES import)
# or mint a deterministic uuid5 from the key and record it on the net, so a
# later export of the same net stays identical.
function _cgmesId(net::Sparlectra.Net, key::AbstractString)::String
  return get!(net.cgmes_ids, String(key)) do
    string(UUIDs.uuid5(CGMES_UUID_NAMESPACE, String(key)))
  end
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function xmlEscape(s::AbstractString)
  s = replace(s, "&" => "&amp;")
  s = replace(s, "<" => "&lt;")
  s = replace(s, ">" => "&gt;")
  s = replace(s, "\"" => "&quot;")
  return s
end

fmtVal(x::Real) = string(Float64(x))

# write `<cim:TAG>value</cim:TAG>` only when the value is known
_optNum(io::IO, tag::AbstractString, v::Union{Nothing,Real}) = v === nothing ? nothing : println(io, "    <cim:", tag, ">", fmtVal(v), "</cim:", tag, ">")
_optBool(io::IO, tag::AbstractString, v::Union{Nothing,Bool}) = v === nothing ? nothing : println(io, "    <cim:", tag, ">", v, "</cim:", tag, ">")

# harvested short-circuit source records by canonical mRID, per class — the
# writer matches them against the exported unit ids (which ARE the canonical
# source mRIDs for captured units)
function _scSourceByMrid(sc::Union{Nothing,CGMESShortCircuitData})
  bymrid = kind -> sc === nothing ? Dict{String,NamedTuple}() : Dict{String,NamedTuple}(cgmesCanonicalMrid(r.mrid) => r for r in getfield(sc, kind))
  return (machines = bymrid(:synchronous_machines), enis = bymrid(:external_network_injections), asms = bymrid(:asynchronous_machines))
end

function writeXmlHeader(io::IO, created::Dates.DateTime)
  println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
  # tool provenance; the stamp is the export date (`created` keyword), so a
  # pinned `created` keeps the file byte-reproducible
  println(io, "<!-- Generated by Sparlectra.jl v$(Sparlectra.version()) on $(Dates.format(created, "yyyy-mm-ddTHH:MM:SSZ")) -->")
  println(io, "<rdf:RDF xmlns:cim=\"$(CIM16_NS)\" xmlns:md=\"$(MD_NS)\" xmlns:rdf=\"$(RDF_NS)\">")
end

function writeFullModel(io::IO, modelId::AbstractString, profiles::Vector{String}, dependentOn::Vector{String}, created::Dates.DateTime)
  # One timestamp drives both header fields — with a caller-supplied `created`
  # the whole file becomes byte-reproducible.
  stamp = Dates.format(created, "yyyy-mm-ddTHH:MM:SSZ")
  println(io, "  <md:FullModel rdf:about=\"urn:uuid:$(modelId)\">")
  println(io, "    <md:Model.scenarioTime>$(stamp)</md:Model.scenarioTime>")
  println(io, "    <md:Model.created>$(stamp)</md:Model.created>")
  println(io, "    <md:Model.description>Sparlectra.jl v$(Sparlectra.version()) export</md:Model.description>")
  println(io, "    <md:Model.version>1</md:Model.version>")
  for p in profiles
    println(io, "    <md:Model.profile>$(p)</md:Model.profile>")
  end
  for d in dependentOn
    println(io, "    <md:Model.DependentOn rdf:resource=\"urn:uuid:$(d)\"/>")
  end
  println(io, "    <md:Model.modelingAuthoritySet>http://sparlectra.local/</md:Model.modelingAuthoritySet>")
  println(io, "  </md:FullModel>")
end

# Shared collection context so the EQ, TP, and SSH profiles use the same
# mRIDs. Trafo/prosumer/shunt records are precomputed NamedTuples because all
# three profile writers need consistent slices of the same data.
struct CGMESContext
  netName::String
  idx2busName::Dict{Int,String}
  baseVoltageIds::Dict{Float64,String}    # vn_kV -> mRID BaseVoltage
  voltageLevelIds::Dict{Float64,String}   # vn_kV -> mRID VoltageLevel
  substationId::String
  regionId::String
  subRegionId::String
  topoNodeIds::Dict{Int,String}           # busIdx -> mRID TopologicalNode
  lineTerminalIds::Vector{NTuple{2,String}} # per line (Terminal1, Terminal2)
  lineIds::Vector{String}                 # mRID ACLineSegment per net.linesAC index
  lineContainerIds::Vector{String}        # mRID cim:Line container per line
  trafoRecs::Vector{NamedTuple}           # one record per exported 2W transformer branch
  trafo3Recs::Vector{NamedTuple}          # one record per reassembled 3W transformer (star group)
  starAuxBuses::Set{Int}                  # star buses consumed by trafo3Recs — no TN/SvVoltage
  prosumerRecs::Vector{NamedTuple}        # one record per exported prosumer
  shuntRecs::Vector{NamedTuple}           # one record per exported shunt
  linkRecs::Vector{NamedTuple}            # one record per exported bus link (Breaker)
  notices::Vector{String}                 # what was NOT exported, one line each
  eqModelId::String
  tpModelId::String
  sshModelId::String
  svModelId::String
end

function _exportBusName(idx2busName::Dict{Int,String}, busIdx::Int)
  return get(idx2busName, busIdx, string("Bus_", busIdx))
end

function buildContext(net::Sparlectra.Net)
  idx2busName = Dict{Int,String}(v => k for (k, v) in net.busDict)

  # Every id resolution is collected so the duplicate guard can name BOTH
  # keys that map to a colliding mRID.
  claimed = Pair{String,String}[]
  claim(key) = begin
    id = _cgmesId(net, key)
    push!(claimed, key => id)
    return id
  end

  baseVoltageIds = Dict{Float64,String}()
  voltageLevelIds = Dict{Float64,String}()
  for node in net.nodeVec
    vn = node.comp.cVN
    if !haskey(baseVoltageIds, vn)
      baseVoltageIds[vn] = claim(cgmesKeyBaseVoltage(vn))
      voltageLevelIds[vn] = claim(cgmesKeyVoltageLevel(vn))
    end
  end

  topoNodeIds = Dict{Int,String}()
  for node in net.nodeVec
    busName = _exportBusName(idx2busName, node.busIdx)
    topoNodeIds[node.busIdx] = claim(cgmesKeyTopologicalNode(busName))
  end

  lineTerminalIds = NTuple{2,String}[]
  lineIds = String[]
  lineContainerIds = String[]
  parallel = Dict{Tuple{String,String},Int}()
  for line in net.linesAC
    busA = _exportBusName(idx2busName, line.comp.cFrom_bus)
    busB = _exportBusName(idx2busName, line.comp.cTo_bus)
    k = cgmesNextParallelIndex!(parallel, busA, busB)
    aclKey = cgmesKeyACLineSegment(busA, busB, k)
    push!(lineIds, claim(aclKey))
    push!(lineContainerIds, claim(cgmesKeyLineContainer(busA, busB, k)))
    push!(lineTerminalIds, (claim(cgmesKeyTerminal(aclKey, 1)), claim(cgmesKeyTerminal(aclKey, 2))))
  end

  notices = String[]
  trafoRecs, trafo3Recs, starAuxBuses = _collectTrafoRecs(net, idx2busName, claim, notices)
  prosumerRecs = _collectProsumerRecs(net, idx2busName, claim, notices)
  shuntRecs = _collectShuntRecs(net, idx2busName, claim)
  linkRecs = _collectLinkRecs(net, idx2busName, claim)

  substationId = claim("SUBSTATION")
  regionId = claim("REGION")
  subRegionId = claim("SUBREGION")
  eqModelId = claim(cgmesKeyModel("EQ"))
  tpModelId = claim(cgmesKeyModel("TP"))
  sshModelId = claim(cgmesKeyModel("SSH"))
  svModelId = claim(cgmesKeyModel("SV"))

  _assertUniqueIds(claimed)

  return CGMESContext(net.name, idx2busName, baseVoltageIds, voltageLevelIds, substationId, regionId, subRegionId, topoNodeIds, lineTerminalIds, lineIds, lineContainerIds, trafoRecs, trafo3Recs, starAuxBuses, prosumerRecs, shuntRecs, linkRecs, notices, eqModelId, tpModelId, sshModelId, svModelId)
end

# Bus links export as closed Breakers (the retained-switch model the importer
# maps back to addLink!). An out-of-service link exports open; the importer
# then keeps the buses separate, which is the same electrical result.
function _collectLinkRecs(net::Sparlectra.Net, idx2busName::Dict{Int,String}, claim)
  recs = NamedTuple[]
  counter = Dict{Tuple{String,String},Int}()
  for lnk in net.linkVec
    busA = _exportBusName(idx2busName, Int(lnk.fromBus))
    busB = _exportBusName(idx2busName, Int(lnk.toBus))
    pair = _cgmesTrafoPairKey(busA, busB)
    k = cgmesNextBusEquipmentIndex!(counter, "LNK", pair)
    key = cgmesKeyBusEquipment("LNK", pair, k)
    push!(
      recs,
      (
        id = claim(key),
        terminalIds = (claim(cgmesKeyTerminal(key, 1)), claim(cgmesKeyTerminal(key, 2))),
        fromIdx = Int(lnk.fromBus),
        toIdx = Int(lnk.toBus),
        name = lnk.cName,
        open = lnk.status != 1,
      ),
    )
  end
  return recs
end

# One record per transformer branch. The branch carries what the solver
# actually uses (pu impedance on the to-side base, effective ratio, shift), so
# every transformer variant — physical windings, PI-model imports, star legs
# of a 3W import — exports uniformly. The effective ratio is written as
# ratedU1 = ratio · vn1 with ratedU2 = vn2: the importer's reconstruction
# (ratedU1/ratedU2)/(vn1/vn2) then reproduces exactly the solved ratio, tap
# positions flattened in.
# live solver values of a transformer branch (calcBranchRatio semantics)
_liveRatio(br)::Float64 = br.ratio == 0.0 ? 1.0 : br.tap_ratio
_liveShift(br)::Float64 = br.ratio == 0.0 ? 0.0 : br.phase_shift_deg

# Detect 3W star groups: an AUX3WT star bus with exactly three transformer
# legs (all from = star bus) and nothing else attached. Returns
# star busIdx → leg indices into `branches`.
function _starGroups(net::Sparlectra.Net, idx2busName::Dict{Int,String}, branches)
  legs_at = Dict{Int,Vector{Int}}()
  for (i, br) in enumerate(branches)
    name = uppercase(_exportBusName(idx2busName, Int(br.fromBus)))
    startswith(name, "AUX3WT") && push!(get!(Vector{Int}, legs_at, Int(br.fromBus)), i)
  end
  groups = Dict{Int,Vector{Int}}()
  for (aux, legs) in legs_at
    length(legs) == 3 || continue
    # nothing but the three legs may touch the star bus
    touches = count(br -> Int(br.fromBus) == aux || Int(br.toBus) == aux, net.branchVec)
    touches == 3 || continue
    any(ps -> Int(ps.comp.cFrom_bus) == aux, net.prosumpsVec) && continue
    any(sh -> Int(sh.busIdx) == aux, net.shuntVec) && continue
    groups[aux] = legs
  end
  return groups
end

function _collectTrafoRecs(net::Sparlectra.Net, idx2busName::Dict{Int,String}, claim, notices::Vector{String})
  recs = NamedTuple[]
  recs3 = NamedTuple[]
  starAux = Set{Int}()
  # PowerTransformer-built branches carry the "2WT" kind marker in their
  # component name (the Branch constructor names lines "ACL" and generic
  # PI models "PI"); the comp cTyp is the generic BranchC for all of them
  branch_idx = [bi for (bi, br) in enumerate(net.branchVec) if occursin("_2WT_", br.comp.cName)]
  branches = [net.branchVec[bi] for bi in branch_idx]
  if length(branches) != length(net.trafos)
    # branch/trafo pairing is positional; a mismatch means an unknown
    # construction path — exporting a guess would corrupt the file
    isempty(net.trafos) || push!(notices, "$(length(net.trafos)) transformer(s) not exported (unrecognized branch structure)")
    return recs, recs3, starAux
  end

  # Star groups whose importer reconstruction is EXACT become one 3W
  # transformer; every other group stays in the 2W star representation.
  # Reconstruction check: the importer places its star bus at the end-1 bus
  # voltage base and derives leg-k ratio (U1/Uk)/(vn_1/vn_k) — that forces
  # vn_star == vn_1 and leg-1 ratio 1; the remaining ratios are free via the
  # end ratedU values.
  in_3w = Set{Int}()
  pt3_counter = Dict{Tuple{String,String},Int}()
  for (aux, legs) in sort(collect(_starGroups(net, idx2busName, branches)); by = first)
    lbrs = [branches[i] for i in legs]
    vn_aux = Sparlectra.getNodeVn(net.nodeVec[aux])
    vns = [Sparlectra.getNodeVn(net.nodeVec[br.toBus]) for br in lbrs]
    ratios = [_liveRatio(br) for br in lbrs]
    # the unit name comes from the star bus name (leg components carry only
    # generic names): unique per unit, and a re-import recreates the same
    # star bus name from it
    name = replace(_exportBusName(idx2busName, aux), r"^(AUX3WT_|Aux3WT_)" => "")
    if !(isapprox(vn_aux, vns[1]; rtol = 1e-9) && isapprox(ratios[1], 1.0; rtol = 1e-9))
      push!(notices, "3W transformer $(name): star reconstruction not exact (leg-1 ratio $(ratios[1]), star base $(vn_aux)/$(vns[1]) kV) — exported as star equivalent")
      continue
    end
    side_buses = [_exportBusName(idx2busName, Int(br.toBus)) for br in lbrs]
    k3 = cgmesNextBusEquipmentIndex!(pt3_counter, "PT3", join(sort(side_buses), "|"))
    pt3Key = cgmesKeyPowerTransformer3W(side_buses[1], side_buses[2], side_buses[3], k3)
    ends = NamedTuple[]
    for (e, br) in enumerate(lbrs)
      r, x, b, g = Sparlectra.fromPU_RXBG(r_pu = br.r_pu, x_pu = br.x_pu, g_pu = br.g_pu, b_pu = br.b_pu, v_kv = vns[e], baseMVA = net.baseMVA)
      shift = _liveShift(br)
      push!(
        ends,
        (
          endId = claim(cgmesKeyTransformerEnd(pt3Key, e)),
          terminalId = claim(cgmesKeyTerminal(pt3Key, e)),
          # the importer applies 3W end taps on the to side, so the shift
          # increment enters negated
          ptcId = shift == 0.0 ? nothing : claim(string(pt3Key, "|PTC", e)),
          angle = -shift,
          busIdx = Int(br.toBus),
          vn = vns[e],
          ratedU = e == 1 ? vn_aux : vns[e] / ratios[e],
          r = r,
          x = x,
          g = g,
          b = b,
          ratedS = br.sn_MVA,
          connected = br.status == 1,
          branchIdx = Int(br.branchIdx),
        ),
      )
    end
    push!(recs3, (pt3Id = claim(pt3Key), ends = ends, name = name))
    push!(starAux, aux)
    union!(in_3w, legs)
  end

  pair_counter = Dict{Tuple{String,String},Int}()
  for (i, br) in enumerate(branches)
    i in in_3w && continue
    tf = net.trafos[i]
    busA = _exportBusName(idx2busName, Int(br.fromBus))
    busB = _exportBusName(idx2busName, Int(br.toBus))
    vn1 = Sparlectra.getNodeVn(net.nodeVec[br.fromBus])
    vn2 = Sparlectra.getNodeVn(net.nodeVec[br.toBus])
    name = tf.comp.cName
    # the solver reads the LIVE tap fields (calcBranchRatio), not the nominal
    # construction values — a tap-controller run must export its final state
    eff_ratio = _liveRatio(br)
    shift = _liveShift(br)
    r, x, b, g = Sparlectra.fromPU_RXBG(r_pu = br.r_pu, x_pu = br.x_pu, g_pu = br.g_pu, b_pu = br.b_pu, v_kv = vn2, baseMVA = net.baseMVA)
    k = cgmesNextParallelIndex!(pair_counter, busA, busB)
    ptKey = cgmesKeyPowerTransformer(busA, busB, k)
    ratedS = br.sn_MVA !== nothing ? br.sn_MVA : tf.side1.ratedS
    # A winding with ratio-tap data re-exports its tap machinery (range,
    # neutral, current step); ratedU1 absorbs the exact residual so the live
    # solved ratio survives regardless of the step correction the importer
    # will re-apply (end-1 changers multiply the ratedU ratio, end-2 divide).
    rtc = nothing
    ratedU1 = eff_ratio * vn1
    tapside = tf.side1.taps !== nothing ? (1, tf.side1.taps) : (tf.side2.taps !== nothing ? (2, tf.side2.taps) : nothing)
    if tapside !== nothing
      side, taps = tapside
      corr = 1.0 + (taps.step - taps.neutralStep) * taps.tapStepPercent / 100.0
      if corr > 0.0
        ratedU1 = side == 1 ? ratedU1 / corr : ratedU1 * corr
        rtc = (rtcId = claim(string(ptKey, "|RTC")), side = side, low = taps.lowStep, high = taps.highStep, neutral = taps.neutralStep, step = taps.step, pct = taps.tapStepPercent, neutralU = taps.neutralU)
      end
    end
    push!(
      recs,
      (
        ptId = claim(ptKey),
        endIds = (claim(cgmesKeyTransformerEnd(ptKey, 1)), claim(cgmesKeyTransformerEnd(ptKey, 2))),
        terminalIds = (claim(cgmesKeyTerminal(ptKey, 1)), claim(cgmesKeyTerminal(ptKey, 2))),
        # phase shift travels as a single-step linear phase tap changer on
        # end 1 (step 1, neutral 0, increment = shift) — the importer's end-1
        # tap application adds it back unnegated
        ptcId = shift == 0.0 ? nothing : claim(string(ptKey, "|PTC")),
        angle = shift,
        rtc = rtc,
        fromIdx = Int(br.fromBus),
        toIdx = Int(br.toBus),
        vn1 = vn1,
        vn2 = vn2,
        ratedU1 = ratedU1,
        ratedU2 = vn2,
        r = r,
        x = x,
        g = g,
        b = b,
        ratedS = ratedS,
        name = name,
        connected = br.status == 1,
        branchIdx = Int(br.branchIdx),
      ),
    )
  end
  return recs, recs3, starAux
end

# Concrete export class per prosumer, from the component-type marker the
# importers maintain (programmatic/MATPOWER nets carry the generic
# Generator/Load pair and classify as machine/load).
function _prosumerExportKind(cTyp)::Union{Nothing,Symbol}
  (cTyp == Sparlectra.SynchronousMachine || cTyp == Sparlectra.Generator) && return :sm
  cTyp == Sparlectra.ExternalNetworkInjection && return :eni
  cTyp == Sparlectra.AsynchronousMachine && return :asm
  (cTyp == Sparlectra.Load || cTyp == Sparlectra.EnergyConsumer) && return :ec
  cTyp == Sparlectra.StaticVarCompensator && return :svc
  return nothing
end

function _collectProsumerRecs(net::Sparlectra.Net, idx2busName::Dict{Int,String}, claim, notices::Vector{String})
  recs = NamedTuple[]
  slackset = Set{Int}(net.slackVec)
  counter = Dict{Tuple{String,String},Int}()
  for ps in net.prosumpsVec
    kind = _prosumerExportKind(ps.comp.cTyp)
    busIdx = Int(ps.comp.cFrom_bus)
    bus = _exportBusName(idx2busName, busIdx)
    if kind === nothing
      push!(notices, "$(ps.comp.cTyp) $(ps.comp.cName) at $(bus) not exported")
      continue
    end
    kindTag = kind === :sm ? "SM" : kind === :eni ? "ENI" : kind === :asm ? "ASM" : kind === :svc ? "SVC" : "EC"
    k = cgmesNextBusEquipmentIndex!(counter, kindTag, bus)
    key = cgmesKeyBusEquipment(kindTag, bus, k)
    vn = ps.comp.cVN
    pv = something(ps.pVal, 0.0)
    qv = something(ps.qVal, 0.0)
    # CGMES load convention in SSH: machine-kind units consume with p > 0,
    # so injections stored injection-positive flip sign
    machine_kind = kind === :sm || kind === :eni || kind === :svc
    rc_target_kV = (machine_kind && ps.isRegulated && ps.vm_pu !== nothing) ? ps.vm_pu * vn : nothing
    push!(
      recs,
      (
        kind = kind,
        id = claim(key),
        terminalId = claim(cgmesKeyTerminal(key, 1)),
        guId = kind === :sm ? claim(string(key, "|GU")) : nothing,
        rcId = rc_target_kV === nothing ? nothing : claim(string(key, "|RC")),
        busIdx = busIdx,
        vn = vn,
        name = ps.comp.cName,
        p_ssh = machine_kind ? -pv : pv,
        q_ssh = machine_kind ? -qv : qv,
        referencePriority = ((kind === :sm || kind === :eni) && busIdx in slackset) ? 1 : 0,
        rc_target_kV = rc_target_kV,
        ratedS = ps.ratedS,
        minQ = ps.minQ,
        maxQ = ps.maxQ,
        minP = ps.minP,
        maxP = ps.maxP,
        normalPF = ps.participationFactor,
      ),
    )
  end
  return recs
end

function _collectShuntRecs(net::Sparlectra.Net, idx2busName::Dict{Int,String}, claim)
  recs = NamedTuple[]
  counter = Dict{Tuple{String,String},Int}()
  for sh in net.shuntVec
    busIdx = Int(sh.busIdx)
    bus = _exportBusName(idx2busName, busIdx)
    k = cgmesNextBusEquipmentIndex!(counter, "SH", bus)
    key = cgmesKeyBusEquipment("SH", bus, k)
    vn = sh.vn_kV
    # invert the MATPOWER-style stamping: y_pu = (Gs + jBs)/baseMVA with
    # Gs/Bs in MW/MVAr at 1 pu → per-section admittance in S at vn
    Gs = real(sh.y_pu_shunt) * net.baseMVA
    Bs = imag(sh.y_pu_shunt) * net.baseMVA
    push!(
      recs,
      (
        id = claim(key),
        terminalId = claim(cgmesKeyTerminal(key, 1)),
        busIdx = busIdx,
        vn = vn,
        name = sh.comp.cName,
        gPerSection = vn > 0.0 ? Gs / vn^2 : 0.0,
        bPerSection = vn > 0.0 ? Bs / vn^2 : 0.0,
        sections = sh.status == 1 ? 1 : 0,
      ),
    )
  end
  return recs
end

# Hard guard: two different structural keys must never share one mRID —
# writing such a file would produce invalid RDF (duplicate rdf:ID). Aborts
# before any file is opened.
function _assertUniqueIds(claimed::Vector{Pair{String,String}})
  seen = Dict{String,String}()
  for (key, id) in claimed
    other = get(seen, id, nothing)
    if other !== nothing && other != key
      error("CGMES export: mRID $(id) is assigned to both \"$(other)\" and \"$(key)\" in net.cgmes_ids — export aborted, no files written")
    end
    seen[id] = key
  end
  return nothing
end

# single-step linear phase tap changer carrying a branch phase shift; the
# step position (1) lives in SSH
function _writeLinearPtc(io::IO, ptcId::AbstractString, name::AbstractString, endId::AbstractString, increment::Float64)
  println(io, "  <cim:PhaseTapChangerLinear rdf:ID=\"_$(ptcId)\">")
  println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(name))_PTC</cim:IdentifiedObject.name>")
  println(io, "    <cim:PhaseTapChanger.TransformerEnd rdf:resource=\"#_$(endId)\"/>")
  println(io, "    <cim:TapChanger.lowStep>0</cim:TapChanger.lowStep>")
  println(io, "    <cim:TapChanger.highStep>1</cim:TapChanger.highStep>")
  println(io, "    <cim:TapChanger.neutralStep>0</cim:TapChanger.neutralStep>")
  println(io, "    <cim:TapChanger.normalStep>1</cim:TapChanger.normalStep>")
  println(io, "    <cim:PhaseTapChangerLinear.stepPhaseShiftIncrement>$(fmtVal(increment))</cim:PhaseTapChangerLinear.stepPhaseShiftIncrement>")
  println(io, "  </cim:PhaseTapChangerLinear>")
end

# Always deliver line parameters in physical units (Ohm/S). PI-model lines
# carry p.u. values and are converted back.
function lineParamsOhm(net::Sparlectra.Net, line)
  if line._isPIModel
    vn = line.comp.cVN
    g_pu = isnothing(line.g) ? 0.0 : line.g
    b_pu = isnothing(line.b) ? 0.0 : line.b
    r, x, b, g = Sparlectra.fromPU_RXBG(r_pu = line.r, x_pu = line.x, g_pu = g_pu, b_pu = b_pu, v_kv = vn, baseMVA = net.baseMVA)
    return r, x, b, g
  else
    r, x, b, g = Sparlectra.getLineRXBG(line)
    return r, x, (isnothing(b) ? 0.0 : b), (isnothing(g) ? 0.0 : g)
  end
end

# ---------------------------------------------------------------------------
# EQ profile
# ---------------------------------------------------------------------------

function writeEQFile(net::Sparlectra.Net, ctx::CGMESContext, path::AbstractString, created::Dates.DateTime; sc_line_data::Dict{Int,CGMESLineShortCircuit} = Dict{Int,CGMESLineShortCircuit}(), sc_source::Union{Nothing,CGMESShortCircuitData} = nothing)
  sc_units = _scSourceByMrid(sc_source)
  open(path, "w") do io
    writeXmlHeader(io, created)
    writeFullModel(io, ctx.eqModelId, EQ_PROFILES, String[], created)

    # Regions / substation / voltage levels (minimal container hierarchy).
    println(io, "  <cim:GeographicalRegion rdf:ID=\"_$(ctx.regionId)\">")
    println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(ctx.netName))_Region</cim:IdentifiedObject.name>")
    println(io, "  </cim:GeographicalRegion>")

    println(io, "  <cim:SubGeographicalRegion rdf:ID=\"_$(ctx.subRegionId)\">")
    println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(ctx.netName))_SubRegion</cim:IdentifiedObject.name>")
    println(io, "    <cim:SubGeographicalRegion.Region rdf:resource=\"#_$(ctx.regionId)\"/>")
    println(io, "  </cim:SubGeographicalRegion>")

    println(io, "  <cim:Substation rdf:ID=\"_$(ctx.substationId)\">")
    println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(ctx.netName))_Substation</cim:IdentifiedObject.name>")
    println(io, "    <cim:Substation.Region rdf:resource=\"#_$(ctx.subRegionId)\"/>")
    println(io, "  </cim:Substation>")

    for (vn, bvId) in sort(collect(ctx.baseVoltageIds); by = first)
      println(io, "  <cim:BaseVoltage rdf:ID=\"_$(bvId)\">")
      println(io, "    <cim:IdentifiedObject.name>BV_$(vn)_kV</cim:IdentifiedObject.name>")
      println(io, "    <cim:BaseVoltage.nominalVoltage>$(fmtVal(vn))</cim:BaseVoltage.nominalVoltage>")
      println(io, "  </cim:BaseVoltage>")

      vlId = ctx.voltageLevelIds[vn]
      println(io, "  <cim:VoltageLevel rdf:ID=\"_$(vlId)\">")
      println(io, "    <cim:IdentifiedObject.name>VL_$(vn)_kV</cim:IdentifiedObject.name>")
      println(io, "    <cim:VoltageLevel.Substation rdf:resource=\"#_$(ctx.substationId)\"/>")
      println(io, "    <cim:VoltageLevel.BaseVoltage rdf:resource=\"#_$(bvId)\"/>")
      println(io, "  </cim:VoltageLevel>")
    end

    # ACLineSegments + terminals.
    for (i, line) in enumerate(net.linesAC)
      lname = line.comp.cName
      vn = line.comp.cVN
      r, x, b, g = lineParamsOhm(net, line)

      println(io, "  <cim:Line rdf:ID=\"_$(ctx.lineContainerIds[i])\">")
      println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(lname))_Line</cim:IdentifiedObject.name>")
      println(io, "    <cim:Line.Region rdf:resource=\"#_$(ctx.subRegionId)\"/>")
      println(io, "  </cim:Line>")

      println(io, "  <cim:ACLineSegment rdf:ID=\"_$(ctx.lineIds[i])\">")
      println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(lname))</cim:IdentifiedObject.name>")
      println(io, "    <cim:Equipment.EquipmentContainer rdf:resource=\"#_$(ctx.lineContainerIds[i])\"/>")
      println(io, "    <cim:ConductingEquipment.BaseVoltage rdf:resource=\"#_$(ctx.baseVoltageIds[vn])\"/>")
      println(io, "    <cim:Conductor.length>$(fmtVal(line.length))</cim:Conductor.length>")
      println(io, "    <cim:ACLineSegment.r>$(fmtVal(r))</cim:ACLineSegment.r>")
      println(io, "    <cim:ACLineSegment.x>$(fmtVal(x))</cim:ACLineSegment.x>")
      println(io, "    <cim:ACLineSegment.bch>$(fmtVal(b))</cim:ACLineSegment.bch>")
      println(io, "    <cim:ACLineSegment.gch>$(fmtVal(g))</cim:ACLineSegment.gch>")

      # Short-circuit/zero-sequence data (EquipmentShortCircuit) — only when supplied.
      if haskey(sc_line_data, i)
        sc = sc_line_data[i]
        println(io, "    <cim:ACLineSegment.r0>$(fmtVal(sc.r0_ohm))</cim:ACLineSegment.r0>")
        println(io, "    <cim:ACLineSegment.x0>$(fmtVal(sc.x0_ohm))</cim:ACLineSegment.x0>")
        println(io, "    <cim:ACLineSegment.b0ch>$(fmtVal(sc.b0ch_S))</cim:ACLineSegment.b0ch>")
        println(io, "    <cim:ACLineSegment.g0ch>$(fmtVal(sc.g0ch_S))</cim:ACLineSegment.g0ch>")
        if !isnothing(sc.endTemperature_C)
          println(io, "    <cim:ACLineSegment.shortCircuitEndTemperature>$(fmtVal(sc.endTemperature_C))</cim:ACLineSegment.shortCircuitEndTemperature>")
        end
      end
      println(io, "  </cim:ACLineSegment>")

      # Terminals (EQ side: equipment association + sequence order).
      for (seq, termId) in enumerate(ctx.lineTerminalIds[i])
        println(io, "  <cim:Terminal rdf:ID=\"_$(termId)\">")
        println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(lname))_T$(seq)</cim:IdentifiedObject.name>")
        println(io, "    <cim:Terminal.ConductingEquipment rdf:resource=\"#_$(ctx.lineIds[i])\"/>")
        println(io, "    <cim:ACDCTerminal.sequenceNumber>$(seq)</cim:ACDCTerminal.sequenceNumber>")
        println(io, "  </cim:Terminal>")
      end
    end

    # Transformers: two ends per unit; the full impedance sits on end 2 (the
    # branch model's to-side base), end 1 carries the effective ratio via its
    # ratedU. That is exactly the layout the importer referral reproduces.
    for rec in ctx.trafoRecs
      println(io, "  <cim:PowerTransformer rdf:ID=\"_$(rec.ptId)\">")
      println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))</cim:IdentifiedObject.name>")
      println(io, "    <cim:Equipment.EquipmentContainer rdf:resource=\"#_$(ctx.substationId)\"/>")
      println(io, "  </cim:PowerTransformer>")
      for e in 1:2
        vn_e = e == 1 ? rec.vn1 : rec.vn2
        println(io, "  <cim:PowerTransformerEnd rdf:ID=\"_$(rec.endIds[e])\">")
        println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))_E$(e)</cim:IdentifiedObject.name>")
        println(io, "    <cim:PowerTransformerEnd.PowerTransformer rdf:resource=\"#_$(rec.ptId)\"/>")
        println(io, "    <cim:TransformerEnd.Terminal rdf:resource=\"#_$(rec.terminalIds[e])\"/>")
        println(io, "    <cim:TransformerEnd.endNumber>$(e)</cim:TransformerEnd.endNumber>")
        println(io, "    <cim:TransformerEnd.BaseVoltage rdf:resource=\"#_$(ctx.baseVoltageIds[vn_e])\"/>")
        println(io, "    <cim:PowerTransformerEnd.ratedU>$(fmtVal(e == 1 ? rec.ratedU1 : rec.ratedU2))</cim:PowerTransformerEnd.ratedU>")
        rec.ratedS === nothing || println(io, "    <cim:PowerTransformerEnd.ratedS>$(fmtVal(rec.ratedS))</cim:PowerTransformerEnd.ratedS>")
        println(io, "    <cim:PowerTransformerEnd.r>$(fmtVal(e == 1 ? 0.0 : rec.r))</cim:PowerTransformerEnd.r>")
        println(io, "    <cim:PowerTransformerEnd.x>$(fmtVal(e == 1 ? 0.0 : rec.x))</cim:PowerTransformerEnd.x>")
        println(io, "    <cim:PowerTransformerEnd.g>$(fmtVal(e == 1 ? 0.0 : rec.g))</cim:PowerTransformerEnd.g>")
        println(io, "    <cim:PowerTransformerEnd.b>$(fmtVal(e == 1 ? 0.0 : rec.b))</cim:PowerTransformerEnd.b>")
        println(io, "  </cim:PowerTransformerEnd>")
        println(io, "  <cim:Terminal rdf:ID=\"_$(rec.terminalIds[e])\">")
        println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))_T$(e)</cim:IdentifiedObject.name>")
        println(io, "    <cim:Terminal.ConductingEquipment rdf:resource=\"#_$(rec.ptId)\"/>")
        println(io, "    <cim:ACDCTerminal.sequenceNumber>$(e)</cim:ACDCTerminal.sequenceNumber>")
        println(io, "  </cim:Terminal>")
      end
      if rec.ptcId !== nothing
        _writeLinearPtc(io, rec.ptcId, rec.name, rec.endIds[1], rec.angle)
      end
      if rec.rtc !== nothing
        rt = rec.rtc
        println(io, "  <cim:RatioTapChanger rdf:ID=\"_$(rt.rtcId)\">")
        println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))_RTC</cim:IdentifiedObject.name>")
        println(io, "    <cim:RatioTapChanger.TransformerEnd rdf:resource=\"#_$(rec.endIds[rt.side])\"/>")
        println(io, "    <cim:TapChanger.lowStep>$(rt.low)</cim:TapChanger.lowStep>")
        println(io, "    <cim:TapChanger.highStep>$(rt.high)</cim:TapChanger.highStep>")
        println(io, "    <cim:TapChanger.neutralStep>$(rt.neutral)</cim:TapChanger.neutralStep>")
        println(io, "    <cim:TapChanger.normalStep>$(rt.neutral)</cim:TapChanger.normalStep>")
        println(io, "    <cim:TapChanger.neutralU>$(fmtVal(rt.neutralU))</cim:TapChanger.neutralU>")
        println(io, "    <cim:RatioTapChanger.stepVoltageIncrement>$(fmtVal(rt.pct))</cim:RatioTapChanger.stepVoltageIncrement>")
        println(io, "  </cim:RatioTapChanger>")
      end
    end

    # Three-winding transformers: one unit, three ends, each with its own
    # impedance and ratedU (see _collectTrafoRecs for the reconstruction
    # contract with the importer).
    for rec in ctx.trafo3Recs
      println(io, "  <cim:PowerTransformer rdf:ID=\"_$(rec.pt3Id)\">")
      println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))</cim:IdentifiedObject.name>")
      println(io, "    <cim:Equipment.EquipmentContainer rdf:resource=\"#_$(ctx.substationId)\"/>")
      println(io, "  </cim:PowerTransformer>")
      for (e, endrec) in enumerate(rec.ends)
        println(io, "  <cim:PowerTransformerEnd rdf:ID=\"_$(endrec.endId)\">")
        println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))_E$(e)</cim:IdentifiedObject.name>")
        println(io, "    <cim:PowerTransformerEnd.PowerTransformer rdf:resource=\"#_$(rec.pt3Id)\"/>")
        println(io, "    <cim:TransformerEnd.Terminal rdf:resource=\"#_$(endrec.terminalId)\"/>")
        println(io, "    <cim:TransformerEnd.endNumber>$(e)</cim:TransformerEnd.endNumber>")
        haskey(ctx.baseVoltageIds, endrec.vn) && println(io, "    <cim:TransformerEnd.BaseVoltage rdf:resource=\"#_$(ctx.baseVoltageIds[endrec.vn])\"/>")
        println(io, "    <cim:PowerTransformerEnd.ratedU>$(fmtVal(endrec.ratedU))</cim:PowerTransformerEnd.ratedU>")
        endrec.ratedS === nothing || println(io, "    <cim:PowerTransformerEnd.ratedS>$(fmtVal(endrec.ratedS))</cim:PowerTransformerEnd.ratedS>")
        println(io, "    <cim:PowerTransformerEnd.r>$(fmtVal(endrec.r))</cim:PowerTransformerEnd.r>")
        println(io, "    <cim:PowerTransformerEnd.x>$(fmtVal(endrec.x))</cim:PowerTransformerEnd.x>")
        println(io, "    <cim:PowerTransformerEnd.g>$(fmtVal(endrec.g))</cim:PowerTransformerEnd.g>")
        println(io, "    <cim:PowerTransformerEnd.b>$(fmtVal(endrec.b))</cim:PowerTransformerEnd.b>")
        println(io, "  </cim:PowerTransformerEnd>")
        println(io, "  <cim:Terminal rdf:ID=\"_$(endrec.terminalId)\">")
        println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))_T$(e)</cim:IdentifiedObject.name>")
        println(io, "    <cim:Terminal.ConductingEquipment rdf:resource=\"#_$(rec.pt3Id)\"/>")
        println(io, "    <cim:ACDCTerminal.sequenceNumber>$(e)</cim:ACDCTerminal.sequenceNumber>")
        println(io, "  </cim:Terminal>")
        if endrec.ptcId !== nothing
          _writeLinearPtc(io, endrec.ptcId, string(rec.name, "_E", e), endrec.endId, endrec.angle)
        end
      end
    end

    # Prosumers: loads, machines (with GeneratingUnit and optional voltage
    # RegulatingControl), external network injections, asynchronous machines.
    for rec in ctx.prosumerRecs
      vnKnown = haskey(ctx.baseVoltageIds, rec.vn)
      if rec.kind === :ec
        println(io, "  <cim:EnergyConsumer rdf:ID=\"_$(rec.id)\">")
        println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))</cim:IdentifiedObject.name>")
        vnKnown && println(io, "    <cim:Equipment.EquipmentContainer rdf:resource=\"#_$(ctx.voltageLevelIds[rec.vn])\"/>")
        println(io, "  </cim:EnergyConsumer>")
      elseif rec.kind === :sm
        println(io, "  <cim:SynchronousMachine rdf:ID=\"_$(rec.id)\">")
        println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))</cim:IdentifiedObject.name>")
        vnKnown && println(io, "    <cim:ConductingEquipment.BaseVoltage rdf:resource=\"#_$(ctx.baseVoltageIds[rec.vn])\"/>")
        println(io, "    <cim:RotatingMachine.GeneratingUnit rdf:resource=\"#_$(rec.guId)\"/>")
        rec.ratedS === nothing || println(io, "    <cim:RotatingMachine.ratedS>$(fmtVal(rec.ratedS))</cim:RotatingMachine.ratedS>")
        rec.minQ === nothing || println(io, "    <cim:SynchronousMachine.minQ>$(fmtVal(rec.minQ))</cim:SynchronousMachine.minQ>")
        rec.maxQ === nothing || println(io, "    <cim:SynchronousMachine.maxQ>$(fmtVal(rec.maxQ))</cim:SynchronousMachine.maxQ>")
        rec.rcId === nothing || println(io, "    <cim:RegulatingCondEq.RegulatingControl rdf:resource=\"#_$(rec.rcId)\"/>")
        # harvested short-circuit source attributes ride along by mRID so a
        # re-imported delivery keeps its Ik'' evaluation inputs
        sm_sc = get(sc_units.machines, rec.id, nothing)
        if sm_sc !== nothing
          # the prosumer model does not retain ratedS/ratedU — without them
          # the machine impedance x''d·U²/S is unusable on re-import
          rec.ratedS === nothing && _optNum(io, "RotatingMachine.ratedS", sm_sc.ratedS_MVA)
          _optNum(io, "SynchronousMachine.satDirectSubtransX", sm_sc.satDirectSubtransX_pu)
          _optNum(io, "SynchronousMachine.satDirectTransX", sm_sc.satDirectTransX_pu)
          _optNum(io, "SynchronousMachine.r0", sm_sc.r0_pu)
          _optNum(io, "SynchronousMachine.x0", sm_sc.x0_pu)
          _optNum(io, "SynchronousMachine.r2", sm_sc.r2_pu)
          _optNum(io, "SynchronousMachine.x2", sm_sc.x2_pu)
          _optBool(io, "SynchronousMachine.earthing", sm_sc.earthing)
          _optNum(io, "RotatingMachine.ratedU", sm_sc.ratedU_kV)
        end
        println(io, "  </cim:SynchronousMachine>")
        println(io, "  <cim:GeneratingUnit rdf:ID=\"_$(rec.guId)\">")
        println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))_GU</cim:IdentifiedObject.name>")
        rec.minP === nothing || println(io, "    <cim:GeneratingUnit.minOperatingP>$(fmtVal(rec.minP))</cim:GeneratingUnit.minOperatingP>")
        rec.maxP === nothing || println(io, "    <cim:GeneratingUnit.maxOperatingP>$(fmtVal(rec.maxP))</cim:GeneratingUnit.maxOperatingP>")
        rec.normalPF === nothing || println(io, "    <cim:GeneratingUnit.normalPF>$(fmtVal(rec.normalPF))</cim:GeneratingUnit.normalPF>")
        println(io, "  </cim:GeneratingUnit>")
      elseif rec.kind === :eni
        println(io, "  <cim:ExternalNetworkInjection rdf:ID=\"_$(rec.id)\">")
        println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))</cim:IdentifiedObject.name>")
        vnKnown && println(io, "    <cim:ConductingEquipment.BaseVoltage rdf:resource=\"#_$(ctx.baseVoltageIds[rec.vn])\"/>")
        rec.minP === nothing || println(io, "    <cim:ExternalNetworkInjection.minP>$(fmtVal(rec.minP))</cim:ExternalNetworkInjection.minP>")
        rec.maxP === nothing || println(io, "    <cim:ExternalNetworkInjection.maxP>$(fmtVal(rec.maxP))</cim:ExternalNetworkInjection.maxP>")
        rec.minQ === nothing || println(io, "    <cim:ExternalNetworkInjection.minQ>$(fmtVal(rec.minQ))</cim:ExternalNetworkInjection.minQ>")
        rec.maxQ === nothing || println(io, "    <cim:ExternalNetworkInjection.maxQ>$(fmtVal(rec.maxQ))</cim:ExternalNetworkInjection.maxQ>")
        rec.rcId === nothing || println(io, "    <cim:RegulatingCondEq.RegulatingControl rdf:resource=\"#_$(rec.rcId)\"/>")
        eni_sc = get(sc_units.enis, rec.id, nothing)
        if eni_sc !== nothing
          _optNum(io, "ExternalNetworkInjection.maxInitialSymShCCurrent", eni_sc.maxInitialSymShCCurrent_A)
          _optNum(io, "ExternalNetworkInjection.minInitialSymShCCurrent", eni_sc.minInitialSymShCCurrent_A)
          _optNum(io, "ExternalNetworkInjection.maxR1ToX1Ratio", eni_sc.maxR1ToX1Ratio)
          _optNum(io, "ExternalNetworkInjection.minR1ToX1Ratio", eni_sc.minR1ToX1Ratio)
          _optNum(io, "ExternalNetworkInjection.maxR0ToX0Ratio", eni_sc.maxR0ToX0Ratio)
          _optNum(io, "ExternalNetworkInjection.maxZ0ToZ1Ratio", eni_sc.maxZ0ToZ1Ratio)
          _optBool(io, "ExternalNetworkInjection.ikSecond", eni_sc.ikSecond)
          _optNum(io, "ExternalNetworkInjection.governorSCD", eni_sc.governorSCD)
        end
        println(io, "  </cim:ExternalNetworkInjection>")
      elseif rec.kind === :asm
        println(io, "  <cim:AsynchronousMachine rdf:ID=\"_$(rec.id)\">")
        println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))</cim:IdentifiedObject.name>")
        vnKnown && println(io, "    <cim:ConductingEquipment.BaseVoltage rdf:resource=\"#_$(ctx.baseVoltageIds[rec.vn])\"/>")
        asm_sc = get(sc_units.asms, rec.id, nothing)
        if asm_sc !== nothing
          _optNum(io, "AsynchronousMachine.iaIrRatio", asm_sc.iaIrRatio)
          _optNum(io, "AsynchronousMachine.rxLockedRotorRatio", asm_sc.rxLockedRotorRatio)
          _optNum(io, "AsynchronousMachine.efficiency", asm_sc.efficiency_percent)
          _optNum(io, "AsynchronousMachine.ratedMechanicalPower", asm_sc.ratedMechanicalPower_MW)
          _optNum(io, "AsynchronousMachine.polePairNumber", asm_sc.polePairNumber)
          _optNum(io, "RotatingMachine.ratedPowerFactor", asm_sc.ratedPowerFactor)
          _optNum(io, "RotatingMachine.ratedS", asm_sc.ratedS_MVA)
          _optNum(io, "RotatingMachine.ratedU", asm_sc.ratedU_kV)
        end
        println(io, "  </cim:AsynchronousMachine>")
      elseif rec.kind === :svc
        println(io, "  <cim:StaticVarCompensator rdf:ID=\"_$(rec.id)\">")
        println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))</cim:IdentifiedObject.name>")
        vnKnown && println(io, "    <cim:ConductingEquipment.BaseVoltage rdf:resource=\"#_$(ctx.baseVoltageIds[rec.vn])\"/>")
        # ratings are reactances (Q = vn²/X) — the exact inverse of the
        # importer's Q-limit derivation
        (rec.maxQ !== nothing && rec.maxQ != 0.0) && _optNum(io, "StaticVarCompensator.capacitiveRating", rec.vn^2 / rec.maxQ)
        (rec.minQ !== nothing && rec.minQ != 0.0) && _optNum(io, "StaticVarCompensator.inductiveRating", rec.vn^2 / rec.minQ)
        rec.rcId === nothing || println(io, "    <cim:RegulatingCondEq.RegulatingControl rdf:resource=\"#_$(rec.rcId)\"/>")
        println(io, "  </cim:StaticVarCompensator>")
      end
      if rec.rc_target_kV !== nothing
        # voltage regulation at the unit's own terminal; targetValue lives in SSH
        println(io, "  <cim:RegulatingControl rdf:ID=\"_$(rec.rcId)\">")
        println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))_RC</cim:IdentifiedObject.name>")
        println(io, "    <cim:RegulatingControl.mode rdf:resource=\"$(CIM16_NS)RegulatingControlModeKind.voltage\"/>")
        println(io, "    <cim:RegulatingControl.Terminal rdf:resource=\"#_$(rec.terminalId)\"/>")
        println(io, "  </cim:RegulatingControl>")
      end
      println(io, "  <cim:Terminal rdf:ID=\"_$(rec.terminalId)\">")
      println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))_T1</cim:IdentifiedObject.name>")
      println(io, "    <cim:Terminal.ConductingEquipment rdf:resource=\"#_$(rec.id)\"/>")
      println(io, "    <cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber>")
      println(io, "  </cim:Terminal>")
    end

    # Bus links: closed Breakers between the two topological nodes.
    for rec in ctx.linkRecs
      println(io, "  <cim:Breaker rdf:ID=\"_$(rec.id)\">")
      println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))</cim:IdentifiedObject.name>")
      println(io, "    <cim:Switch.open>$(rec.open)</cim:Switch.open>")
      println(io, "  </cim:Breaker>")
      for (seq, termId) in enumerate(rec.terminalIds)
        println(io, "  <cim:Terminal rdf:ID=\"_$(termId)\">")
        println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))_T$(seq)</cim:IdentifiedObject.name>")
        println(io, "    <cim:Terminal.ConductingEquipment rdf:resource=\"#_$(rec.id)\"/>")
        println(io, "    <cim:ACDCTerminal.sequenceNumber>$(seq)</cim:ACDCTerminal.sequenceNumber>")
        println(io, "  </cim:Terminal>")
      end
    end

    # Shunts: one-section linear compensators; the switched state (sections)
    # lives in SSH.
    for rec in ctx.shuntRecs
      vnKnown = haskey(ctx.baseVoltageIds, rec.vn)
      println(io, "  <cim:LinearShuntCompensator rdf:ID=\"_$(rec.id)\">")
      println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))</cim:IdentifiedObject.name>")
      vnKnown && println(io, "    <cim:ConductingEquipment.BaseVoltage rdf:resource=\"#_$(ctx.baseVoltageIds[rec.vn])\"/>")
      println(io, "    <cim:LinearShuntCompensator.bPerSection>$(fmtVal(rec.bPerSection))</cim:LinearShuntCompensator.bPerSection>")
      println(io, "    <cim:LinearShuntCompensator.gPerSection>$(fmtVal(rec.gPerSection))</cim:LinearShuntCompensator.gPerSection>")
      println(io, "    <cim:ShuntCompensator.normalSections>1</cim:ShuntCompensator.normalSections>")
      println(io, "    <cim:ShuntCompensator.maximumSections>1</cim:ShuntCompensator.maximumSections>")
      println(io, "  </cim:LinearShuntCompensator>")
      println(io, "  <cim:Terminal rdf:ID=\"_$(rec.terminalId)\">")
      println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(rec.name))_T1</cim:IdentifiedObject.name>")
      println(io, "    <cim:Terminal.ConductingEquipment rdf:resource=\"#_$(rec.id)\"/>")
      println(io, "    <cim:ACDCTerminal.sequenceNumber>1</cim:ACDCTerminal.sequenceNumber>")
      println(io, "  </cim:Terminal>")
    end

    println(io, "</rdf:RDF>")
  end
  return path
end

# ---------------------------------------------------------------------------
# TP profile
# ---------------------------------------------------------------------------

function writeTPFile(net::Sparlectra.Net, ctx::CGMESContext, path::AbstractString, created::Dates.DateTime)
  open(path, "w") do io
    writeXmlHeader(io, created)
    writeFullModel(io, ctx.tpModelId, TP_PROFILES, [ctx.eqModelId], created)

    # TopologicalNodes. Star buses consumed by a 3W reassembly stay internal
    # to the exported transformer and get no node of their own.
    for node in net.nodeVec
      node.busIdx in ctx.starAuxBuses && continue
      busName = _exportBusName(ctx.idx2busName, node.busIdx)
      vn = node.comp.cVN
      tnId = ctx.topoNodeIds[node.busIdx]
      println(io, "  <cim:TopologicalNode rdf:ID=\"_$(tnId)\">")
      println(io, "    <cim:IdentifiedObject.name>$(xmlEscape(busName))</cim:IdentifiedObject.name>")
      println(io, "    <cim:TopologicalNode.BaseVoltage rdf:resource=\"#_$(ctx.baseVoltageIds[vn])\"/>")
      println(io, "    <cim:TopologicalNode.ConnectivityNodeContainer rdf:resource=\"#_$(ctx.voltageLevelIds[vn])\"/>")
      println(io, "  </cim:TopologicalNode>")
    end

    # Terminal -> TopologicalNode association.
    for (i, line) in enumerate(net.linesAC)
      fromIdx = line.comp.cFrom_bus
      toIdx = line.comp.cTo_bus
      busIdxs = (fromIdx, toIdx)
      for (seq, termId) in enumerate(ctx.lineTerminalIds[i])
        busIdx = busIdxs[seq]
        if isnothing(busIdx) || !haskey(ctx.topoNodeIds, busIdx)
          @warn "CGMES-TP: line $(line.comp.cName) terminal $(seq) has no valid bus index — skipped"
          continue
        end
        println(io, "  <cim:Terminal rdf:about=\"#_$(termId)\">")
        println(io, "    <cim:Terminal.TopologicalNode rdf:resource=\"#_$(ctx.topoNodeIds[busIdx])\"/>")
        println(io, "    <cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected>")
        println(io, "  </cim:Terminal>")
      end
    end

    # Transformer terminals: an out-of-service unit exports as disconnected
    # on both ends, which is how the importer reads branch status back.
    for rec in ctx.trafoRecs
      for (e, busIdx) in ((1, rec.fromIdx), (2, rec.toIdx))
        haskey(ctx.topoNodeIds, busIdx) || continue
        println(io, "  <cim:Terminal rdf:about=\"#_$(rec.terminalIds[e])\">")
        println(io, "    <cim:Terminal.TopologicalNode rdf:resource=\"#_$(ctx.topoNodeIds[busIdx])\"/>")
        println(io, "    <cim:ACDCTerminal.connected>$(rec.connected)</cim:ACDCTerminal.connected>")
        println(io, "  </cim:Terminal>")
      end
    end
    for rec in ctx.trafo3Recs
      for endrec in rec.ends
        haskey(ctx.topoNodeIds, endrec.busIdx) || continue
        println(io, "  <cim:Terminal rdf:about=\"#_$(endrec.terminalId)\">")
        println(io, "    <cim:Terminal.TopologicalNode rdf:resource=\"#_$(ctx.topoNodeIds[endrec.busIdx])\"/>")
        println(io, "    <cim:ACDCTerminal.connected>$(endrec.connected)</cim:ACDCTerminal.connected>")
        println(io, "  </cim:Terminal>")
      end
    end

    # Single-terminal equipment (prosumers and shunts).
    for rec in ctx.prosumerRecs
      haskey(ctx.topoNodeIds, rec.busIdx) || continue
      println(io, "  <cim:Terminal rdf:about=\"#_$(rec.terminalId)\">")
      println(io, "    <cim:Terminal.TopologicalNode rdf:resource=\"#_$(ctx.topoNodeIds[rec.busIdx])\"/>")
      println(io, "    <cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected>")
      println(io, "  </cim:Terminal>")
    end
    for rec in ctx.shuntRecs
      haskey(ctx.topoNodeIds, rec.busIdx) || continue
      println(io, "  <cim:Terminal rdf:about=\"#_$(rec.terminalId)\">")
      println(io, "    <cim:Terminal.TopologicalNode rdf:resource=\"#_$(ctx.topoNodeIds[rec.busIdx])\"/>")
      println(io, "    <cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected>")
      println(io, "  </cim:Terminal>")
    end
    for rec in ctx.linkRecs
      for (seq, busIdx) in ((1, rec.fromIdx), (2, rec.toIdx))
        haskey(ctx.topoNodeIds, busIdx) || continue
        println(io, "  <cim:Terminal rdf:about=\"#_$(rec.terminalIds[seq])\">")
        println(io, "    <cim:Terminal.TopologicalNode rdf:resource=\"#_$(ctx.topoNodeIds[busIdx])\"/>")
        println(io, "    <cim:ACDCTerminal.connected>true</cim:ACDCTerminal.connected>")
        println(io, "  </cim:Terminal>")
      end
    end

    println(io, "</rdf:RDF>")
  end
  return path
end

# ---------------------------------------------------------------------------
# SSH profile
# ---------------------------------------------------------------------------

function writeSSHFile(ctx::CGMESContext, path::AbstractString, created::Dates.DateTime)
  open(path, "w") do io
    writeXmlHeader(io, created)
    writeFullModel(io, ctx.sshModelId, SSH_PROFILES, [ctx.eqModelId], created)

    for rec in ctx.prosumerRecs
      if rec.kind === :ec
        println(io, "  <cim:EnergyConsumer rdf:about=\"#_$(rec.id)\">")
        println(io, "    <cim:EnergyConsumer.p>$(fmtVal(rec.p_ssh))</cim:EnergyConsumer.p>")
        println(io, "    <cim:EnergyConsumer.q>$(fmtVal(rec.q_ssh))</cim:EnergyConsumer.q>")
        println(io, "  </cim:EnergyConsumer>")
      elseif rec.kind === :sm
        println(io, "  <cim:SynchronousMachine rdf:about=\"#_$(rec.id)\">")
        println(io, "    <cim:RotatingMachine.p>$(fmtVal(rec.p_ssh))</cim:RotatingMachine.p>")
        println(io, "    <cim:RotatingMachine.q>$(fmtVal(rec.q_ssh))</cim:RotatingMachine.q>")
        rec.referencePriority == 0 || println(io, "    <cim:SynchronousMachine.referencePriority>$(rec.referencePriority)</cim:SynchronousMachine.referencePriority>")
        println(io, "  </cim:SynchronousMachine>")
      elseif rec.kind === :eni
        println(io, "  <cim:ExternalNetworkInjection rdf:about=\"#_$(rec.id)\">")
        println(io, "    <cim:ExternalNetworkInjection.p>$(fmtVal(rec.p_ssh))</cim:ExternalNetworkInjection.p>")
        println(io, "    <cim:ExternalNetworkInjection.q>$(fmtVal(rec.q_ssh))</cim:ExternalNetworkInjection.q>")
        rec.referencePriority == 0 || println(io, "    <cim:ExternalNetworkInjection.referencePriority>$(rec.referencePriority)</cim:ExternalNetworkInjection.referencePriority>")
        println(io, "  </cim:ExternalNetworkInjection>")
      elseif rec.kind === :asm
        println(io, "  <cim:AsynchronousMachine rdf:about=\"#_$(rec.id)\">")
        println(io, "    <cim:RotatingMachine.p>$(fmtVal(rec.p_ssh))</cim:RotatingMachine.p>")
        println(io, "    <cim:RotatingMachine.q>$(fmtVal(rec.q_ssh))</cim:RotatingMachine.q>")
        println(io, "  </cim:AsynchronousMachine>")
      elseif rec.kind === :svc
        println(io, "  <cim:StaticVarCompensator rdf:about=\"#_$(rec.id)\">")
        println(io, "    <cim:StaticVarCompensator.q>$(fmtVal(rec.q_ssh))</cim:StaticVarCompensator.q>")
        println(io, "  </cim:StaticVarCompensator>")
      end
      if rec.rc_target_kV !== nothing
        println(io, "  <cim:RegulatingControl rdf:about=\"#_$(rec.rcId)\">")
        println(io, "    <cim:RegulatingControl.enabled>true</cim:RegulatingControl.enabled>")
        println(io, "    <cim:RegulatingControl.targetValue>$(fmtVal(rec.rc_target_kV))</cim:RegulatingControl.targetValue>")
        println(io, "  </cim:RegulatingControl>")
      end
    end

    for rec in ctx.shuntRecs
      println(io, "  <cim:LinearShuntCompensator rdf:about=\"#_$(rec.id)\">")
      println(io, "    <cim:ShuntCompensator.sections>$(rec.sections)</cim:ShuntCompensator.sections>")
      println(io, "  </cim:LinearShuntCompensator>")
    end

    for rec in ctx.trafoRecs
      if rec.ptcId !== nothing
        println(io, "  <cim:PhaseTapChangerLinear rdf:about=\"#_$(rec.ptcId)\">")
        println(io, "    <cim:TapChanger.step>1</cim:TapChanger.step>")
        println(io, "  </cim:PhaseTapChangerLinear>")
      end
      if rec.rtc !== nothing
        println(io, "  <cim:RatioTapChanger rdf:about=\"#_$(rec.rtc.rtcId)\">")
        println(io, "    <cim:TapChanger.step>$(rec.rtc.step)</cim:TapChanger.step>")
        println(io, "  </cim:RatioTapChanger>")
      end
    end
    for rec in ctx.trafo3Recs
      for endrec in rec.ends
        endrec.ptcId === nothing && continue
        println(io, "  <cim:PhaseTapChangerLinear rdf:about=\"#_$(endrec.ptcId)\">")
        println(io, "    <cim:TapChanger.step>1</cim:TapChanger.step>")
        println(io, "  </cim:PhaseTapChangerLinear>")
      end
    end

    println(io, "</rdf:RDF>")
  end
  return path
end

"""
    cgmesLineShortCircuitData(result::CGMESImportResult) -> Dict{Int,CGMESLineShortCircuit}

Build the `sc_line_data` input of [`writeCGMESFiles`](@ref) from the
zero-sequence line attributes harvested during a CGMES import, keyed by
`net.linesAC` index. Lines are matched through the structural identity keys
recorded on the network, so parallel lines keep their own data. Only lines
with a complete zero-sequence pair (`r0` and `x0`) are included — nothing is
invented for the rest.
"""
function cgmesLineShortCircuitData(result::CGMESImportResult)::Dict{Int,CGMESLineShortCircuit}
  out = Dict{Int,CGMESLineShortCircuit}()
  harvest = Dict{String,NamedTuple}()
  for rec in result.shortcircuit.ac_line_segments
    harvest[cgmesCanonicalMrid(rec.mrid)] = rec
  end
  isempty(harvest) && return out
  net = result.net
  idx2busName = Dict{Int,String}(v => k for (k, v) in net.busDict)
  parallel = Dict{Tuple{String,String},Int}()
  for (i, line) in enumerate(net.linesAC)
    busA = _exportBusName(idx2busName, line.comp.cFrom_bus)
    busB = _exportBusName(idx2busName, line.comp.cTo_bus)
    k = cgmesNextParallelIndex!(parallel, busA, busB)
    mrid = get(net.cgmes_ids, cgmesKeyACLineSegment(busA, busB, k), nothing)
    mrid === nothing && continue
    rec = get(harvest, mrid, nothing)
    rec === nothing && continue
    (rec.r0_ohm === nothing || rec.x0_ohm === nothing) && continue
    out[i] = CGMESLineShortCircuit(
      r0_ohm = rec.r0_ohm,
      x0_ohm = rec.x0_ohm,
      b0ch_S = something(rec.b0ch_S, 0.0),
      g0ch_S = something(rec.g0ch_S, 0.0),
      endTemperature_C = rec.shortCircuitEndTemperature_C,
    )
  end
  return out
end

# ---------------------------------------------------------------------------
# SV profile
# ---------------------------------------------------------------------------

# SvVoltage/SvPowerFlow objects are snapshot artifacts of one export — they
# carry no identity worth preserving across roundtrips, so their ids are
# minted deterministically from the referenced object's id instead of going
# through net.cgmes_ids.
_svObjectId(kind::AbstractString, refId::AbstractString)::String = string(UUIDs.uuid5(CGMES_UUID_NAMESPACE, string(kind, "|", refId)))

function _writeSvPowerFlow(io::IO, termId::AbstractString, p::Float64, q::Float64)
  println(io, "  <cim:SvPowerFlow rdf:ID=\"_$(_svObjectId("SVPF", termId))\">")
  println(io, "    <cim:SvPowerFlow.p>$(fmtVal(p))</cim:SvPowerFlow.p>")
  println(io, "    <cim:SvPowerFlow.q>$(fmtVal(q))</cim:SvPowerFlow.q>")
  println(io, "    <cim:SvPowerFlow.Terminal rdf:resource=\"#_$(termId)\"/>")
  println(io, "  </cim:SvPowerFlow>")
end

"""
Write the SV profile from the network's CURRENT voltage state: after a solve
this is the solution, right after an import it is the start state (for a
CGMES import: the delivery's own SV values). SvPowerFlow rows follow the
same conventions `compareWithSV` reads back — loads and units in load
convention, branch-terminal flows from the branch model at the current
voltages — so a re-import validates cleanly against this profile.
"""
function writeSVFile(net::Sparlectra.Net, ctx::CGMESContext, path::AbstractString, created::Dates.DateTime)
  iso = Set(net.isoNodes)
  V = ComplexF64[something(nd._vm_pu, 1.0) * cis(deg2rad(something(nd._va_deg, 0.0))) for nd in net.nodeVec]
  open(path, "w") do io
    writeXmlHeader(io, created)
    writeFullModel(io, ctx.svModelId, SV_PROFILES, [ctx.tpModelId, ctx.sshModelId], created)

    for node in net.nodeVec
      node.busIdx in ctx.starAuxBuses && continue
      tnId = get(ctx.topoNodeIds, node.busIdx, nothing)
      tnId === nothing && continue
      vm = something(node._vm_pu, 1.0)
      va = something(node._va_deg, 0.0)
      println(io, "  <cim:SvVoltage rdf:ID=\"_$(_svObjectId("SVV", tnId))\">")
      println(io, "    <cim:SvVoltage.v>$(fmtVal(vm * node.comp.cVN))</cim:SvVoltage.v>")
      println(io, "    <cim:SvVoltage.angle>$(fmtVal(va))</cim:SvVoltage.angle>")
      println(io, "    <cim:SvVoltage.TopologicalNode rdf:resource=\"#_$(tnId)\"/>")
      println(io, "  </cim:SvVoltage>")
    end

    # branch-terminal flows from the branch model at the current voltages —
    # the same expressions compareWithSV evaluates on the read side
    branch_flows = function (br)
      ys = inv(br.r_pu + im * br.x_pu)
      ysh2 = (br.g_pu + im * br.b_pu) / 2
      tr = br.ratio == 0.0 ? 1.0 + 0im : br.tap_ratio * cis(deg2rad(br.phase_shift_deg))
      Vf, Vt = V[br.fromBus], V[br.toBus]
      Sfrom = Vf * conj(((ys + ysh2) / abs2(tr)) * Vf - (ys / conj(tr)) * Vt) * net.baseMVA
      Sto = Vt * conj((ys + ysh2) * Vt - (ys / tr) * Vf) * net.baseMVA
      return Sfrom, Sto
    end
    line_branches = [br for br in net.branchVec if occursin("_ACL_", br.comp.cName)]
    if length(line_branches) == length(net.linesAC)
      for (i, br) in enumerate(line_branches)
        (br.status == 0 || br.fromBus in iso || br.toBus in iso) && continue
        Sfrom, Sto = branch_flows(br)
        _writeSvPowerFlow(io, ctx.lineTerminalIds[i][1], real(Sfrom), imag(Sfrom))
        _writeSvPowerFlow(io, ctx.lineTerminalIds[i][2], real(Sto), imag(Sto))
      end
    end
    for rec in ctx.trafoRecs
      br = net.branchVec[rec.branchIdx]
      (br.status == 0 || br.fromBus in iso || br.toBus in iso) && continue
      Sfrom, Sto = branch_flows(br)
      _writeSvPowerFlow(io, rec.terminalIds[1], real(Sfrom), imag(Sfrom))
      _writeSvPowerFlow(io, rec.terminalIds[2], real(Sto), imag(Sto))
    end
    # 3W ends: each exported terminal sits on a leg's side bus — its flow is
    # the leg's to-side flow (the star bus is internal)
    for rec in ctx.trafo3Recs
      for endrec in rec.ends
        br = net.branchVec[endrec.branchIdx]
        (br.status == 0 || br.fromBus in iso || br.toBus in iso) && continue
        _, Sto = branch_flows(br)
        _writeSvPowerFlow(io, endrec.terminalId, real(Sto), imag(Sto))
      end
    end

    # loads and shunts: model values at the current voltages
    for rec in ctx.prosumerRecs
      (rec.kind === :ec || rec.kind === :asm) || continue
      rec.busIdx in iso && continue
      _writeSvPowerFlow(io, rec.terminalId, rec.p_ssh, rec.q_ssh)
    end
    for rec in ctx.shuntRecs
      rec.busIdx in iso && continue
      vm2 = abs2(V[rec.busIdx])
      _writeSvPowerFlow(io, rec.terminalId, rec.gPerSection * rec.sections * rec.vn^2 * vm2, -rec.bPerSection * rec.sections * rec.vn^2 * vm2)
    end

    # units: the bus total is what the model demands at the current state
    # (-(net injection + local load)); the first unit on a bus absorbs the
    # residual against the scheduled values, so per-bus sums are exact even
    # when the solver moved Q or the slack picked up P
    units = [rec for rec in ctx.prosumerRecs if rec.kind === :sm || rec.kind === :eni]
    if !isempty(units)
      n_bus = length(net.nodeVec)
      Yred = Sparlectra.createYBUS(net = net, sparse = true)
      Y = size(Yred, 1) == n_bus ? Yred : Sparlectra._expand_ybus_for_isolated_nodes(Yred, n_bus, net.isoNodes)
      S = V .* conj.(Y * V) .* net.baseMVA
      load_sum = zeros(ComplexF64, n_bus)
      for rec in ctx.prosumerRecs
        (rec.kind === :ec || rec.kind === :asm) || continue
        load_sum[rec.busIdx] += complex(rec.p_ssh, rec.q_ssh)
      end
      scheduled = Dict{Int,ComplexF64}()
      for rec in units
        scheduled[rec.busIdx] = get(scheduled, rec.busIdx, zero(ComplexF64)) + complex(rec.p_ssh, rec.q_ssh)
      end
      first_at_bus = Set{Int}()
      for rec in units
        rec.busIdx in iso && continue
        pf = complex(rec.p_ssh, rec.q_ssh)
        if !(rec.busIdx in first_at_bus)
          push!(first_at_bus, rec.busIdx)
          total = -(S[rec.busIdx] + load_sum[rec.busIdx])
          pf += total - scheduled[rec.busIdx]
        end
        _writeSvPowerFlow(io, rec.terminalId, real(pf), imag(pf))
      end
    end

    println(io, "</rdf:RDF>")
  end
  return path
end

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

"""
    writeCGMESFiles(net; path::AbstractString = pwd(),
                    sc_line_data::Dict{Int,CGMESLineShortCircuit} = Dict(),
                    sc_source::Union{Nothing,CGMESShortCircuitData} = nothing,
                    created::Dates.DateTime = Dates.now(),
                    notices::Union{Nothing,Vector{String}} = nothing,
                    zip::Bool = false)

Exports a Sparlectra `Net` as CGMES-2.4.15 profile files: EQ + TP (buses, AC
lines, transformers, loads, machines, external network injections,
asynchronous machines, shunts), SSH (operating points p/q, slack
`referencePriority`, voltage-regulation targets, shunt sections), and SV
(the network's current voltage state with per-terminal `SvPowerFlow` rows —
after a solve: the solution). `sc_line_data` supplies optional zero-sequence
data per line index (order of `net.linesAC`); `sc_source` writes the
machine/injection/motor short-circuit attributes of a CGMES import harvest
back onto the units by mRID. `zip = true` additionally packs everything into
a re-importable `<name>_CGMES.zip`. `created` sets the
`md:Model.scenarioTime`/`md:Model.created` header stamp — pass a fixed value
to make the output byte-reproducible.

Object mRIDs come from `net.cgmes_ids` (structural keys, see
`cgmes_keys.jl`): ids recorded by the CGMES importer are reused, missing ids
are minted deterministically (uuid5 over the key) and recorded. Two different
keys resolving to the same mRID abort the export before any file is written.

Works for any `Net` — built programmatically, or imported from CGMES,
MATPOWER, or DTF — since the export reads directly from the network model.

Returns a vector of the written file paths.

Transformer phase shifts travel as single-step linear phase tap changers, so
the roundtrip preserves them exactly. What the profiles cannot carry yet is
reported per unit: static VAr compensators and bus links. Pass a `notices`
vector to collect these lines programmatically; without one they are emitted
as warnings. Ratio tap-changer machinery is flattened into the fixed
effective ratio.

Every file carries the tool provenance — a `Generated by Sparlectra.jl
v<version> on <stamp>` header comment and a matching `md:Model.description`
— where the stamp is the `created` timestamp.
"""
function writeCGMESFiles(net::Sparlectra.Net; path::AbstractString = pwd(), sc_line_data::Dict{Int,CGMESLineShortCircuit} = Dict{Int,CGMESLineShortCircuit}(), sc_source::Union{Nothing,CGMESShortCircuitData} = nothing, created::Dates.DateTime = Dates.now(), notices::Union{Nothing,Vector{String}} = nothing, zip::Bool = false)
  ctx = buildContext(net)
  if notices === nothing
    for n in ctx.notices
      @warn "CGMES export: $(n)"
    end
  else
    append!(notices, ctx.notices)
  end

  eqPath = joinpath(path, string(net.name, "_EQ.xml"))
  tpPath = joinpath(path, string(net.name, "_TP.xml"))
  sshPath = joinpath(path, string(net.name, "_SSH.xml"))
  svPath = joinpath(path, string(net.name, "_SV.xml"))

  writeEQFile(net, ctx, eqPath, created; sc_line_data = sc_line_data, sc_source = sc_source)
  writeTPFile(net, ctx, tpPath, created)
  writeSSHFile(ctx, sshPath, created)
  writeSVFile(net, ctx, svPath, created)
  files = [eqPath, tpPath, sshPath, svPath]

  if zip
    # one re-importable delivery file (collectCGMESFiles reads ZIPs directly)
    zipPath = joinpath(path, string(net.name, "_CGMES.zip"))
    open(zipPath, "w") do io
      ZipArchives.ZipWriter(io) do w
        for f in files
          ZipArchives.zip_newfile(w, basename(f))
          write(w, read(f))
        end
      end
    end
    push!(files, zipPath)
  end

  return files
end
