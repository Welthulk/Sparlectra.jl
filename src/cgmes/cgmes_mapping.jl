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

# file: src/cgmes/cgmes_mapping.jl
# purpose: layer 5 of the CGMES importer — map CIM objects onto the
# Sparlectra Net (Stage 1: bus-branch, lines, 2W/3W transformers with fixed
# SSH tap positions, shunts, loads, machines, equivalent injections, slack
# selection D-2). Short-circuit source data is always harvested into
# `CGMESShortCircuitData` (concept §7.7) — read, not evaluated.

"""
Typed short-circuit source data harvested during import (concept §7.7).
Each entry carries the CIM mRID, the object name and — where resolvable —
the Sparlectra bus name. Values stay in CGMES units (Ω, S, A, pu on machine
base); nothing is evaluated in Stage 1.
"""
struct CGMESShortCircuitData
  external_network_injections::Vector{NamedTuple}
  synchronous_machines::Vector{NamedTuple}
  ac_line_segments::Vector{NamedTuple}
  transformer_ends::Vector{NamedTuple}
  equivalent_injections::Vector{NamedTuple}
end

"""
Result of `importCGMES`: the `Net`, the merged `CGMESStore`, the bus-branch
topology, the harvested short-circuit data, the selected slack bus and the
importer's skip/notice messages. `branch_side_of_terminal` (CGMES Terminal
mRID → `(branchIdx, :from|:to)`) and `skipped_equipment` (mRIDs that were not
mapped) provide the provenance `compareWithSV` needs for its `SvPowerFlow`
comparison.
"""
struct CGMESImportResult
  net::Sparlectra.Net
  store::CGMESStore
  topo::CGMESTopology
  shortcircuit::CGMESShortCircuitData
  slack_bus::String
  messages::Vector{String}
  branch_side_of_terminal::Dict{String,Tuple{Int,Symbol}}
  skipped_equipment::Set{String}
end

# mutable importer bookkeeping threaded through the mapping helpers.
# vset_min_pu/vset_max_pu ride along here rather than through every mapping
# signature; see _voltageSetpoint for what the band does.
# split_sides: boundary CN → (source file → per-side bus name), filled by
# _detectNonCancellingBoundarySides! before the element mapping runs.
struct _MapCtx
  messages::Vector{String}
  branch_side::Dict{String,Tuple{Int,Symbol}}
  skipped::Set{String}
  tap_control::Bool
  ignore_connected::Bool
  vset_min_pu::Float64
  vset_max_pu::Float64
  split_sides::Dict{String,Dict{String,String}}
end
_MapCtx(tap_control::Bool = false, ignore_connected::Bool = false; vset_min_pu::Float64 = CGMES_VSET_MIN_PU, vset_max_pu::Float64 = CGMES_VSET_MAX_PU) =
  _MapCtx(String[], Dict{String,Tuple{Int,Symbol}}(), Set{String}(), tap_control, ignore_connected, vset_min_pu, vset_max_pu, Dict{String,Dict{String,String}}())

# SSH `Terminal.connected` through the ctx switch: `ignore_connected = true`
# treats everything as connected (diagnostic mode for snapshots whose SSH
# flags contradict their SV state)
_conn(ctx::_MapCtx, flag::Bool) = ctx.ignore_connected || flag

# SSH `Equipment.inService` — CGMES 3.0 moved the operational status onto the
# equipment itself; 2.4.15 deliveries do not carry it (absent → in service).
# This is NOT redundant with Terminal.connected: ReliCapGrid parks
# out-of-service generators with connected=true, inService=false, no SvVoltage
# and a placeholder regulation target. Ignoring the flag imported 413 MW of
# phantom generation into Svedala alone — the single largest reason the solved
# corridor flows missed the delivery's SV state. Honors ignore_connected for
# consistency with _conn (the override means "treat everything as live").
_inService(ctx::_MapCtx, o::CIMObject)::Bool = ctx.ignore_connected || something(boolval(o, :inService), true)

# --- tap changers (fixed operating point in Stage 1) ------------------------
#
# Step source order: SvTapStep position (the solved state — present in every
# conformity data set and in real exchanges) → SSH `step` (initial operating
# point) → neutralStep. The MicroGrid BaseCase SV was produced with active tap
# control, so SvTapStep differs from SSH for regulated units; a fixed-tap
# import must use the solved position to reproduce SV. Stage 2 controllers
# will start from SSH instead.

"""SvTapStep positions: tap-changer mRID → solved position."""
function _svTapPositions(store::CGMESStore)::Dict{String,Float64}
  out = Dict{String,Float64}()
  for sv in objectsOf(store, :SvTapStep)
    tc = get(sv.refs, :TapChanger, nothing)
    tc === nothing && continue
    p = num(sv, :position)
    p === nothing && continue
    out[tc] = p
  end
  return out
end

_tapStep(tc::CIMObject, svsteps::Dict{String,Float64}) = something(get(svsteps, tc.mrid, nothing), num(tc, :step), num(tc, :neutralStep, 0.0))

function _ratioTapCorrection(rtc::CIMObject, svsteps::Dict{String,Float64})::Float64
  step = _tapStep(rtc, svsteps)
  neutral = num(rtc, :neutralStep, 0.0)
  incr = num(rtc, :stepVoltageIncrement, 0.0)
  return 1.0 + (step - neutral) * incr / 100.0
end

function _phaseTapRatioShift(ptc::CIMObject, svsteps::Dict{String,Float64})::Tuple{Float64,Float64}
  step = Int(round(_tapStep(ptc, svsteps)))
  neutral = Int(round(num(ptc, :neutralStep, 0.0)))
  low = Int(round(num(ptc, :lowStep, Float64(min(step, neutral)))))
  high = Int(round(num(ptc, :highStep, Float64(max(step, neutral)))))
  if ptc.class == :PhaseTapChangerLinear
    shift = (step - neutral) * num(ptc, :stepPhaseShiftIncrement, 0.0)
    return (1.0, shift)
  end
  kind = ptc.class == :PhaseTapChangerSymmetrical ? :symmetrical : :asymmetrical
  vincr = num(ptc, :voltageStepIncrement)
  model = Sparlectra.PhaseTapChangerModel(;
    kind = kind,
    step = step,
    lowStep = low,
    highStep = high,
    neutralStep = neutral,
    # CGMES voltageStepIncrement is in percent; PhaseTapChangerModel expects
    # the per-step fraction (see calcPhaseTapFraction)
    voltage_step_increment = vincr === nothing ? nothing : vincr / 100.0,
    winding_connection_angle_deg = num(ptc, :windingConnectionAngle),
    # the branch's complex ratio t applies the regulating vector directly on
    # the tapped winding (validated against MicroGrid SV, SvTapStep=16)
    convention = :direct_regulating_vector,
  )
  res = Sparlectra.calcPhaseTapAngleRatio(model)
  return (res.effective_ratio, res.effective_shift_deg)
end

# --- Stage-2 tap control (importCGMES(..., tap_control = true)) -------------
#
# With `tap_control = true` the importer starts from the SSH tap positions
# and attaches Sparlectra outer-loop controllers (`addPowerTransformerControl!`)
# for every tap changer whose `controlEnabled` flag AND TapChangerControl
# `enabled` flag are set. The reference SvTapStep positions are then found by
# the control loop (run via `run_sparlectra`), not copied from SV.

"""Branch tap-range fields and controller for one controlled transformer end."""
function _attachTapControl!(net, store::CGMESStore, topo::CGMESTopology, ctx::_MapCtx, svsteps::Dict{String,Float64}, e::CIMObject, branch_idx::Int, on_from_side::Bool, trafo_name::String)
  for tc in _tapChangersOfEnd(store, e.mrid)
    something(boolval(tc, :controlEnabled), false) || continue
    tcc = ref(store, tc, :TapChangerControl)
    tcc === nothing && continue
    something(boolval(tcc, :enabled), false) || continue
    mode = something(enumval(tcc, :mode), "")
    br = net.branchVec[findfirst(b -> b.branchIdx == branch_idx, net.branchVec)]
    low = Int(round(num(tc, :lowStep, 0.0)))
    high = Int(round(num(tc, :highStep, 0.0)))
    neutral = Int(round(num(tc, :neutralStep, 0.0)))
    name = something(str(tc, :name), tc.mrid)

    if tc.class == :RatioTapChanger && endswith(mode, ".voltage")
      # translate step range into branch tap-ratio limits around the branch's
      # nominal ratio (current ratio with the initial correction removed)
      corr0 = _ratioTapCorrection(tc, svsteps)
      rho_nom = on_from_side ? br.tap_ratio / corr0 : br.tap_ratio * corr0
      corr_low = 1.0 + (low - neutral) * num(tc, :stepVoltageIncrement, 0.0) / 100.0
      corr_high = 1.0 + (high - neutral) * num(tc, :stepVoltageIncrement, 0.0) / 100.0
      r1 = on_from_side ? rho_nom * corr_low : rho_nom / corr_high
      r2 = on_from_side ? rho_nom * corr_high : rho_nom / corr_low
      br.has_ratio_tap = true
      br.tap_min = min(r1, r2)
      br.tap_max = max(r1, r2)
      br.tap_step = abs(rho_nom * num(tc, :stepVoltageIncrement, 0.0) / 100.0)
      # regulated bus from the TapChangerControl terminal
      tinfo = _regulatedBus(store, topo, tcc)
      tinfo === nothing && (push!(ctx.messages, "notice: tap controller $(name) — control terminal unresolved, skipped"); continue)
      target = num(tcc, :targetValue, 0.0)
      target <= 0.0 && (push!(ctx.messages, "notice: tap controller $(name) — no positive voltage target, skipped"); continue)
      deadband = num(tcc, :targetDeadband, 0.0)
      Sparlectra.addPowerTransformerControl!(
        net;
        trafo = string(branch_idx),
        mode = :voltage,
        target_bus = tinfo.bus,
        target_vm_pu = target / tinfo.vn_kV,
        deadband_vm_pu = max(deadband / tinfo.vn_kV / 2.0, 1e-3),
        control_ratio = true,
        control_phase = false,
      )
      push!(ctx.messages, "tap control: $(trafo_name)/$(name) → voltage $(round(target / tinfo.vn_kV; digits = 4)) pu at $(tinfo.bus)")
    elseif tc.class != :RatioTapChanger && endswith(mode, ".activePower")
      # phase tap: translate step range into branch phase limits around the
      # branch's base shift (initial phase-tap contribution removed)
      pr0, ps0 = _phaseTapRatioShift(tc, svsteps)
      base_shift = on_from_side ? br.phase_shift_deg - ps0 : br.phase_shift_deg + ps0
      _, ps_low = _phaseTapRatioShift_atstep(tc, low)
      _, ps_high = _phaseTapRatioShift_atstep(tc, high)
      s1 = on_from_side ? base_shift + ps_low : base_shift - ps_low
      s2 = on_from_side ? base_shift + ps_high : base_shift - ps_high
      br.has_phase_tap = true
      br.phase_min_deg = min(s1, s2)
      br.phase_max_deg = max(s1, s2)
      br.phase_step_deg = high > low ? abs(s2 - s1) / (high - low) : 1.0
      # monitored branch side from the control terminal
      tccterm = get(tcc.refs, :Terminal, nothing)
      side = tccterm === nothing ? nothing : get(ctx.branch_side, tccterm, nothing)
      from_name = _busDictName(net, br.fromBus)
      to_name = _busDictName(net, br.toBus)
      (from_name === nothing || to_name === nothing) && (push!(ctx.messages, "notice: tap controller $(name) — bus name lookup failed, skipped"); continue)
      # the controller monitors the branch in stored (from → to) orientation;
      # a control terminal on the to side flips the target sign (the loss
      # asymmetry between the two ends is far below the deadband)
      target_branch = (from_name, to_name)
      target = num(tcc, :targetValue, 0.0)
      (side !== nothing && side[2] == :to) && (target = -target)
      deadband = max(num(tcc, :targetDeadband, 1.0), 0.5)
      Sparlectra.addPowerTransformerControl!(
        net;
        trafo = string(branch_idx),
        mode = :branch_active_power,
        target_branch = target_branch,
        p_target_mw = target,
        deadband_p_mw = deadband / 2.0,
        control_ratio = false,
        control_phase = true,
      )
      push!(ctx.messages, "tap control: $(trafo_name)/$(name) → active power $(target) MW on $(target_branch[1]) → $(target_branch[2])")
    else
      push!(ctx.messages, "notice: tap controller $(name) — unsupported combination $(tc.class)/$(mode), skipped")
    end
  end
  return nothing
end

"""Phase-tap (ratio, shift) of `ptc` evaluated at an explicit step."""
function _phaseTapRatioShift_atstep(ptc::CIMObject, step::Int)::Tuple{Float64,Float64}
  fake = Dict{String,Float64}(ptc.mrid => Float64(step))
  return _phaseTapRatioShift(ptc, fake)
end

# Real deliveries contain connected component(s) without any source: every
# generator/injection in them is switched off, or the feeding branches are
# open (RealGrid: 13 pockets totalling 22 buses). Solving those is impossible
# — they have no reference — so they are de-energized deliberately: loads
# zeroed, branches out of service, buses isolated, one message per component.
# Without this the solver aborts on "island without reference".
"""
De-energize every connected component that contains no source.

A CGMES snapshot regularly carries network parts that are switched off in
reality: open bays, stub connections, feeders parked out of service. They still
appear as buses with loads attached. Left alone, such a component has load but
no generation and no reference — an unsolvable subsystem that would fail the
power flow for a reason that has nothing to do with the interesting network.

The pass therefore zeroes the injections, opens the branches and marks the
buses `Isolated`. This changes the imported model, so every affected component
is reported.

Connectivity is evaluated over in-service branches AND closed links — the
electrical view, since a component fed through a closed switch is energized.
"""
function _deEnergizeSourcelessComponents!(net, ctx::_MapCtx)
  n = length(net.nodeVec)
  n == 0 && return nothing
  parent = collect(1:n)
  find(i) = (while parent[i] != i
    parent[i] = parent[parent[i]]
    i = parent[i]
  end; i)
  function unite(a, b)
    ra, rb = find(a), find(b)
    ra == rb && return nothing
    parent[max(ra, rb)] = min(ra, rb)
    return nothing
  end
  for br in net.branchVec
    br.status == 1 && unite(br.fromBus, br.toBus)
  end
  for l in net.linkVec
    l.status == 1 && unite(Int(l.fromBus), Int(l.toBus))
  end

  # Group buses by component root and record whether any source sits inside.
  # "Source" means generator-like (machines, external network injections) --
  # loads do not energize anything.
  has_source = Dict{Int,Bool}()
  members = Dict{Int,Vector{Int}}()
  for i = 1:n
    r = find(i)
    push!(get!(Vector{Int}, members, r), i)
    haskey(has_source, r) || (has_source[r] = false)
  end
  for ps in net.prosumpsVec
    Sparlectra.isGenerator(ps.proSumptionType) || continue
    has_source[find(ps.comp.cFrom_bus)] = true
  end

  for (r, buses) in members
    has_source[r] && continue
    # A single bus with nothing attached is a harmless leftover (open bay,
    # unconnected boundary node) -- skip it silently instead of drowning the
    # message log. Anything larger, or a lone bus that does carry an injection,
    # is reported.
    length(buses) == 1 && isempty(filter(ps -> ps.comp.cFrom_bus == buses[1], net.prosumpsVec)) && continue
    busset = Set(buses)
    # Zero the injections rather than deleting the prosumers: the components
    # stay visible in reports, they just no longer contribute to the balance.
    nload = 0
    for ps in net.prosumpsVec
      ps.comp.cFrom_bus in busset || continue
      ps.pVal = 0.0
      ps.qVal = 0.0
      nload += 1
    end
    # Any branch touching the component goes out of service. Testing only one
    # end is enough -- a branch with one end inside cannot have the other end in
    # an energized component, or the component would have been connected.
    for br in net.branchVec
      (br.fromBus in busset || br.toBus in busset) && (br.status = 0)
    end
    for i in buses
      i in net.isoNodes && continue
      Sparlectra.setNodeType!(net.nodeVec[i], "Isolated")
      push!(net.isoNodes, i)
    end
    push!(ctx.messages, "notice: de-energized a component of $(length(buses)) bus(es) around $(net.nodeVec[first(buses)].comp.cName) — no source (slack/generator) inside; $(nload) injection(s) zeroed")
  end
  isempty(net.isoNodes) || sort!(unique!(net.isoNodes))
  return nothing
end

# A voltage controller whose target bus is held by a generator (slack or PV)
# cannot influence that voltage — the machine fixes it. Such controllers occur
# in real data sets (MicroGrid BE-TR2_3 regulates the slack busbar); leaving
# them active would burn outer iterations and report a "regulation" that is
# none. Bus types are known only after refreshBusTypesFromProsumers!, so this
# runs as a post-mapping pass.
function _disableIneffectiveTapControllers!(net, ctx::_MapCtx)
  for ctrl in Sparlectra._tap_controllers(net)
    ctrl.enabled || continue
    ctrl.mode in (:voltage, :voltage_and_branch_active_power) || continue
    ctrl.target_bus === nothing && continue
    idx = get(net.busDict, ctrl.target_bus, nothing)
    idx === nothing && continue
    ntype = net.nodeVec[idx]._nodeType
    (ntype == Sparlectra.Slack || ntype == Sparlectra.PV) || continue
    ctrl.enabled = false
    push!(ctx.messages, "notice: tap controller on $(ctrl.target_bus) disabled — target bus is $(ntype) (voltage held by a generator)")
  end
  return nothing
end

"""Bus-dictionary name of node index `idx` (the name `geNetBusIdx` accepts)."""
function _busDictName(net, idx::Int)::Union{Nothing,String}
  for (k, v) in net.busDict
    v == idx && return k
  end
  return nothing
end

"""Regulated bus (name, vn) of a TapChangerControl, or `nothing`."""
function _regulatedBus(store::CGMESStore, topo::CGMESTopology, tcc::CIMObject)
  tmrid = get(tcc.refs, :Terminal, nothing)
  tmrid === nothing && return nothing
  t = get(store.objects, tmrid, nothing)
  t === nothing && return nothing
  tn = tnOfTerminal(t)
  (tn === nothing || !haskey(topo.bus_name, tn)) && return nothing
  return (bus = topo.bus_name[tn], vn_kV = topo.vn_kV[tn])
end

# tap changers attached to a transformer end (RatioTapChanger.TransformerEnd)
function _tapChangersOfEnd(store::CGMESStore, endmrid::String)
  out = CIMObject[]
  for cls in (:RatioTapChanger, :PhaseTapChangerLinear, :PhaseTapChangerSymmetrical, :PhaseTapChangerAsymmetrical, :PhaseTapChangerTabular)
    for tc in objectsOf(store, cls)
      get(tc.refs, :TransformerEnd, "") == endmrid && push!(out, tc)
    end
  end
  return out
end

# fold the fixed tap corrections of transformer end `e` into (ratio, shift).
# `on_from_side`: the end sits on the branch's from bus → ratio multiplies;
# otherwise it divides (referred through the winding).
function _applyEndTaps(store::CGMESStore, svsteps::Dict{String,Float64}, e::CIMObject, ratio::Float64, shift::Float64, on_from_side::Bool, messages::Vector{String})::Tuple{Float64,Float64}
  for tc in _tapChangersOfEnd(store, e.mrid)
    if tc.class == :RatioTapChanger
      corr = _ratioTapCorrection(tc, svsteps)
      ratio = on_from_side ? ratio * corr : ratio / corr
    elseif tc.class == :PhaseTapChangerTabular
      push!(messages, "skip: PhaseTapChangerTabular $(something(str(tc, :name), tc.mrid)) (Stage 4)")
    else
      pratio, pshift = _phaseTapRatioShift(tc, svsteps)
      ratio = on_from_side ? ratio * pratio : ratio / pratio
      shift += on_from_side ? pshift : -pshift
    end
  end
  return (ratio, shift)
end

# --- bus creation (lazy: only TNs referenced by mapped equipment) -----------

function _svVoltageByTN(store::CGMESStore)::Dict{String,Tuple{Float64,Float64}}
  out = Dict{String,Tuple{Float64,Float64}}()
  for sv in objectsOf(store, :SvVoltage)
    tn = get(sv.refs, :TopologicalNode, nothing)
    tn === nothing && continue
    v = num(sv, :v)
    a = num(sv, :angle, 0.0)
    v === nothing && continue
    out[tn] = (v, a)
  end
  return out
end

function _ensureBus!(net::Sparlectra.Net, created::Set{String}, topo::CGMESTopology, svmap::Dict{String,Tuple{Float64,Float64}}, tn::String; isAux::Bool = false, busname::Union{Nothing,String} = nothing)
  # busname overrides the TN's bus for per-side split buses (see
  # _detectNonCancellingBoundarySides!): same TN attributes, different identity.
  bus = something(busname, topo.bus_name[tn])
  bus in created && return bus
  vn = topo.vn_kV[tn]
  vm_pu, va_deg = 1.0, 0.0
  if haskey(svmap, tn) && vn > 0.0
    v, a = svmap[tn]
    vm_pu = v / vn
    va_deg = a
  end
  Sparlectra.addBus!(net = net, busName = bus, vn_kV = vn, vm_pu = vm_pu, va_deg = va_deg, isAux = isAux)
  push!(created, bus)
  return bus
end

"""
Detect assembled boundary nodes whose equivalent injections do **not** cancel
and assign each contributing area (source file) its own side bus.

An AC X-node with both sides present carries two declarations of the *same*
exchange — they cancel, the equivalents are discarded, and the node joins the
areas galvanically. When the pair does NOT cancel, the two injections are two
different devices (in practice: the HVDC converters of a DC border crossing;
the residue is the DC loss). Joining the sides galvanically is then wrong
twice over: it creates an AC corridor where none exists, and it violates KCL
against the delivery's own SV state (ReliCapGrid `BP_SD-EH_DC1`: 94.8 MW flow
into the node on one side, 72.6 MW out on the other — impossible through one
node). Each side therefore keeps its line and its converter injection on its
own bus.

The side of a piece of equipment is the file that defined it — the only
criterion that exists, since both sides' terminals reference the same CN (and
here even the same TN). The first source (sorted) keeps the base bus name,
further sides get `<bus>@<k>`.
"""
function _detectNonCancellingBoundarySides!(ctx::_MapCtx, store::CGMESStore, topo::CGMESTopology)
  by_cn = Dict{String,Vector{Tuple{String,Float64}}}()   # cn → [(source, p)]
  for ei in objectsOf(store, :EquivalentInjection)
    _inService(ctx, ei) || continue
    ts = get(topo.terminals, ei.mrid, nothing)
    (ts === nothing || isempty(ts)) && continue
    cn = get(ts[1].refs, :ConnectivityNode, nothing)
    (cn === nothing || !(cn in store.boundary)) && continue
    push!(get!(Vector{Tuple{String,Float64}}, by_cn, cn), (ei.source, num(ei, :p, 0.0)))
  end
  for (cn, members) in by_cn
    length(members) >= 2 || continue
    _eiPairCancels([m[2] for m in members]) && continue
    sources = sort(unique(first.(members)))
    length(sources) >= 2 || continue
    # find the base bus name via any terminal on this CN that resolves to a TN
    base = nothing
    for t in objectsOf(store, :Terminal)
      get(t.refs, :ConnectivityNode, nothing) == cn || continue
      tn = tnOfTerminal(t)
      (tn === nothing || !haskey(topo.bus_name, tn)) && continue
      base = topo.bus_name[tn]
      break
    end
    base === nothing && continue
    sides = Dict{String,String}()
    for (k, src) in enumerate(sources)
      sides[src] = k == 1 ? base : string(base, "@", k)
    end
    ctx.split_sides[cn] = sides
    push!(
      ctx.messages,
      "notice: boundary node $(base) — the two sides' equivalents do not cancel (device injections, e.g. HVDC converters); node split per side: $(join([string(basename_short(s), " → ", b) for (s, b) in sort(collect(sides))], ", "))",
    )
  end
  return nothing
end

# short human-readable tag for a source file name in messages
basename_short(s::AbstractString) = first(splitext(basename(String(s))))

"""
Bus for one terminal of `eq`, honoring the per-side split of non-cancelling
boundary nodes: when the terminal's CN is split, the equipment lands on its
own area's side bus; everything else falls back to the TN bus.
"""
function _busForTerminal!(net, created, topo::CGMESTopology, svmap, ctx::_MapCtx, eq::CIMObject, seq::Int, tn::String)::String
  if !isempty(ctx.split_sides)
    ts = get(topo.terminals, eq.mrid, nothing)
    if ts !== nothing && length(ts) >= seq
      cn = get(ts[seq].refs, :ConnectivityNode, nothing)
      if cn !== nothing && haskey(ctx.split_sides, cn)
        side = get(ctx.split_sides[cn], eq.source, nothing)
        side === nothing || return _ensureBus!(net, created, topo, svmap, tn; busname = side)
      end
    end
  end
  return _ensureBus!(net, created, topo, svmap, tn)
end

# --- element mapping --------------------------------------------------------

# Retained switches in the bus-branch view: a Switch-family object whose two
# terminals sit on different TopologicalNodes couples those buses (e.g. the
# NL busbar coupler in MicroGrid). Closed → zero-impedance `addLink!`;
# open → ignored. SSH `open` overrides the EQ default `normalOpen`.
const _SWITCH_CLASSES = (:Switch, :Breaker, :Disconnector, :LoadBreakSwitch, :Jumper, :ProtectedSwitch)

function _mapSwitches!(net, store, topo, created, svmap, ctx::_MapCtx)
  for cls in _SWITCH_CLASSES
    for sw in objectsOf(store, cls)
      t1 = busOfEquipment(topo, sw, 1)
      t2 = busOfEquipment(topo, sw, 2)
      (t1 === nothing || t2 === nothing) && continue
      t1.tn == t2.tn && continue
      name = something(str(sw, :name), sw.mrid)
      # an out-of-service switch cannot conduct, whatever its open flag says
      isopen = something(boolval(sw, :open), boolval(sw, :normalOpen), false) || !_inService(ctx, sw)
      if isopen || !(_conn(ctx, t1.connected) && _conn(ctx, t2.connected))
        push!(ctx.messages, "notice: open $(cls) $(name) between $(t1.bus) and $(t2.bus) — buses stay separate")
        continue
      end
      b1 = _ensureBus!(net, created, topo, svmap, t1.tn)
      b2 = _ensureBus!(net, created, topo, svmap, t2.tn)
      Sparlectra.addLink!(net = net, fromBus = b1, toBus = b2, status = 1)
      push!(ctx.messages, "notice: closed $(cls) $(name) — bus link $(b1) = $(b2)")
    end
  end
end

# record the (branchIdx, side) provenance of the just-added branch for the
# two terminals of a line-like equipment (sequence 1 = from, 2 = to)
function _recordBranchTerminals!(ctx::_MapCtx, net, topo::CGMESTopology, eq::CIMObject)
  ts = get(topo.terminals, eq.mrid, nothing)
  ts === nothing && return nothing
  idx = net.branchVec[end].branchIdx
  length(ts) >= 1 && (ctx.branch_side[ts[1].mrid] = (idx, :from))
  length(ts) >= 2 && (ctx.branch_side[ts[2].mrid] = (idx, :to))
  return nothing
end

function _mapLines!(net, store, topo, created, svmap, baseMVA, ctx::_MapCtx)
  for cls in (:ACLineSegment, :SeriesCompensator)
    for line in objectsOf(store, cls)
      t1 = busOfEquipment(topo, line, 1)
      t2 = busOfEquipment(topo, line, 2)
      name = something(str(line, :name), line.mrid)
      if t1 === nothing || t2 === nothing
        push!(ctx.messages, "skip: $(cls) $(name) — terminal without TopologicalNode")
        push!(ctx.skipped, line.mrid)
        continue
      end
      # Terminal.connected handling. Both ends open → the branch is out of
      # service. Exactly ONE end open → the branch carries no longitudinal
      # current, but it still exists: its half charging admittance keeps
      # acting on the closed side. Dropping it outright would shift the
      # reactive balance (RealGrid: 511 such lines against 228 fully open,
      # concentrated on cable-rich sub-transmission levels), so the branch is
      # replaced by that half admittance as a shunt at the closed bus.
      # Equipment.inService == false counts as both ends open: no longitudinal
      # current AND no charging contribution from a de-commissioned line.
      nopen = _inService(ctx, line) ? (_conn(ctx, t1.connected) ? 0 : 1) + (_conn(ctx, t2.connected) ? 0 : 1) : 2
      if nopen == 1
        closed = _conn(ctx, t1.connected) ? t1 : t2
        bus = _ensureBus!(net, created, topo, svmap, closed.tn)
        if cls == :ACLineSegment
          g_half = num(line, :gch, 0.0) / 2.0
          b_half = num(line, :bch, 0.0) / 2.0
          if g_half != 0.0 || b_half != 0.0
            Sparlectra.addShuntMatpower!(net = net, busName = bus, Gs = g_half * closed.vn_kV^2, Bs = b_half * closed.vn_kV^2)
          end
        end
        push!(ctx.messages, "notice: $(cls) $(name) has one open terminal — half charging admittance kept as shunt at $(bus)")
        continue
      end
      status = nopen == 0 ? 1 : 0
      b1 = _busForTerminal!(net, created, topo, svmap, ctx, line, 1, t1.tn)
      b2 = _busForTerminal!(net, created, topo, svmap, ctx, line, 2, t2.tn)
      # branch model convention (calcAdmittance): series/shunt admittance on
      # the TO-side voltage base, complex ratio t at the from side
      r_pu, x_pu, b_pu, g_pu = Sparlectra.toPU_RXBG(
        r = num(line, :r, 0.0),
        x = num(line, :x, 0.0),
        g = cls == :ACLineSegment ? num(line, :gch, 0.0) : 0.0,
        b = cls == :ACLineSegment ? num(line, :bch, 0.0) : 0.0,
        v_kv = t2.vn_kV,
        baseMVA = baseMVA,
      )
      if t1.vn_kV == t2.vn_kV
        Sparlectra.addPIModelACLine!(net = net, fromBus = b1, toBus = b2, r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, g_pu = g_pu, status = status)
        _recordBranchTerminals!(ctx, net, topo, line)
      else
        # boundary line across different nominal voltages (e.g. 380 kV area
        # bus vs 400 kV boundary node): identical physical conductor, so the
        # 2WT formula with ratedU1 == ratedU2 degenerates to ratio = vn2/vn1
        # on a PI branch (impedance on the to-side base).
        ratio = t2.vn_kV / t1.vn_kV
        from = Sparlectra.geNetBusIdx(net = net, busName = b1)
        to = Sparlectra.geNetBusIdx(net = net, busName = b2)
        Sparlectra._addPIModelTrafo_by_idx!(net = net, from = from, to = to, r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, g_pu = g_pu, status = status, ratio = ratio, shift_deg = 0.0)
        _recordBranchTerminals!(ctx, net, topo, line)
        push!(ctx.messages, "notice: $(cls) $(name) spans $(t1.vn_kV)/$(t2.vn_kV) kV — mapped as boundary-line ratio branch ($(round(ratio; digits=4)))")
      end
    end
  end
end

function _endsOfTransformer(store::CGMESStore, pt::CIMObject)::Vector{CIMObject}
  ends = [e for e in objectsOf(store, :PowerTransformerEnd) if get(e.refs, :PowerTransformer, "") == pt.mrid]
  sort!(ends; by = e -> something(num(e, :endNumber), 1.0))
  return ends
end

function _terminalInfo(store::CGMESStore, topo::CGMESTopology, e::CIMObject)
  tmrid = get(e.refs, :Terminal, nothing)
  tmrid === nothing && return nothing
  t = get(store.objects, tmrid, nothing)
  t === nothing && return nothing
  tn = tnOfTerminal(t)
  (tn === nothing || !haskey(topo.bus_name, tn)) && return nothing
  return (bus = topo.bus_name[tn], vn_kV = topo.vn_kV[tn], connected = terminalConnected(t), tn = tn)
end

function _map2WTrafo!(net, store, topo, created, svmap, svsteps, baseMVA, pt, ends, ctx::_MapCtx)
  e1, e2 = ends[1], ends[2]
  i1 = _terminalInfo(store, topo, e1)
  i2 = _terminalInfo(store, topo, e2)
  name = something(str(pt, :name), pt.mrid)
  if i1 === nothing || i2 === nothing
    push!(ctx.messages, "skip: PowerTransformer $(name) — end terminal without TopologicalNode")
    push!(ctx.skipped, pt.mrid)
    return
  end
  U1 = something(num(e1, :ratedU), i1.vn_kV)
  U2 = something(num(e2, :ratedU), i2.vn_kV)
  # refer both end impedances to end 2: the branch model (calcAdmittance)
  # expects series/shunt admittance on the TO-side base, ratio at the from side
  z = Sparlectra.calc2WTEndsReferredRXGB(;
    r1 = num(e1, :r, 0.0), x1 = num(e1, :x, 0.0), g1 = num(e1, :g, 0.0), b1 = num(e1, :b, 0.0),
    r2 = num(e2, :r, 0.0), x2 = num(e2, :x, 0.0), g2 = num(e2, :g, 0.0), b2 = num(e2, :b, 0.0),
    U1 = U1, U2 = U2,
  )
  r_pu, x_pu, b_pu, g_pu = Sparlectra.toPU_RXBG(r = z.r, x = z.x, g = z.g, b = z.b, v_kv = i2.vn_kV, baseMVA = baseMVA)
  ratio = (U1 / U2) / (i1.vn_kV / i2.vn_kV)
  shift = 0.0
  ratio, shift = _applyEndTaps(store, svsteps, e1, ratio, shift, true, ctx.messages)
  ratio, shift = _applyEndTaps(store, svsteps, e2, ratio, shift, false, ctx.messages)
  b1 = _ensureBus!(net, created, topo, svmap, i1.tn)
  b2 = _ensureBus!(net, created, topo, svmap, i2.tn)
  status = (_conn(ctx, i1.connected) && _conn(ctx, i2.connected) && _inService(ctx, pt)) ? 1 : 0
  from = Sparlectra.geNetBusIdx(net = net, busName = b1)
  to = Sparlectra.geNetBusIdx(net = net, busName = b2)
  ratedS = num(e1, :ratedS)
  Sparlectra._addPIModelTrafo_by_idx!(net = net, from = from, to = to, r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, g_pu = g_pu, status = status, ratedU = U1, ratedS = ratedS, ratio = ratio, shift_deg = shift)
  idx = net.branchVec[end].branchIdx
  tm1 = get(e1.refs, :Terminal, nothing)
  tm2 = get(e2.refs, :Terminal, nothing)
  tm1 === nothing || (ctx.branch_side[tm1] = (idx, :from))
  tm2 === nothing || (ctx.branch_side[tm2] = (idx, :to))
  if ctx.tap_control
    _attachTapControl!(net, store, topo, ctx, svsteps, e1, idx, true, name)
    _attachTapControl!(net, store, topo, ctx, svsteps, e2, idx, false, name)
  end
end

function _map3WTrafo!(net, store, topo, created, svmap, svsteps, baseMVA, pt, ends, ctx::_MapCtx)
  name = something(str(pt, :name), pt.mrid)
  infos = [_terminalInfo(store, topo, e) for e in ends]
  if any(i -> i === nothing, infos)
    push!(ctx.messages, "skip: 3W PowerTransformer $(name) — end terminal without TopologicalNode")
    push!(ctx.skipped, pt.mrid)
    return
  end
  e1 = ends[1]
  i1 = infos[1]
  U1 = something(num(e1, :ratedU), i1.vn_kV)
  vn_aux = i1.vn_kV

  # star (AUX) bus at the HV-side voltage base; start voltage from the HV
  # bus SV value — the star point has no SvVoltage of its own
  auxname = "AUX3WT_" * replace(name, r"\s+" => "_")
  if !Sparlectra.hasBusInNet(net = net, busName = auxname)
    vm0, va0 = haskey(svmap, i1.tn) ? (svmap[i1.tn][1] / i1.vn_kV, svmap[i1.tn][2]) : (1.0, 0.0)
    Sparlectra.addBus!(net = net, busName = auxname, vn_kV = vn_aux, vm_pu = vm0, va_deg = va0, isAux = true)
  end
  aux = Sparlectra.geNetBusIdx(net = net, busName = auxname)

  for (k, e) in enumerate(ends)
    ik = infos[k]
    Uk = something(num(e, :ratedU), ik.vn_kV)
    # each leg runs AUX → side bus; end impedances are already on their own
    # end base = the leg's TO side, which is what calcAdmittance expects
    r_pu, x_pu, b_pu, g_pu = Sparlectra.toPU_RXBG(r = num(e, :r, 0.0), x = num(e, :x, 0.0), g = num(e, :g, 0.0), b = num(e, :b, 0.0), v_kv = ik.vn_kV, baseMVA = baseMVA)
    ratio = (U1 / Uk) / (vn_aux / ik.vn_kV)
    shift = 0.0
    # the physical end sits on the branch's "to" side (AUX → side bus)
    ratio, shift = _applyEndTaps(store, svsteps, e, ratio, shift, false, ctx.messages)
    bk = _ensureBus!(net, created, topo, svmap, ik.tn)
    to = Sparlectra.geNetBusIdx(net = net, busName = bk)
    status = (_conn(ctx, ik.connected) && _inService(ctx, pt)) ? 1 : 0
    Sparlectra._addPIModelTrafo_by_idx!(net = net, from = aux, to = to, r_pu = r_pu, x_pu = x_pu, b_pu = b_pu, g_pu = g_pu, status = status, ratedU = Uk, ratedS = num(e, :ratedS), ratio = ratio, shift_deg = shift)
    tm = get(e.refs, :Terminal, nothing)
    tm === nothing || (ctx.branch_side[tm] = (net.branchVec[end].branchIdx, :to))
    if ctx.tap_control
      for tc in _tapChangersOfEnd(store, e.mrid)
        something(boolval(tc, :controlEnabled), false) && push!(ctx.messages, "notice: controlled tap changer on 3W transformer $(name) — controllers on star-equivalent legs are not wired in Stage 2")
      end
    end
  end
end

function _mapTransformers!(net, store, topo, created, svmap, svsteps, baseMVA, ctx::_MapCtx)
  for pt in objectsOf(store, :PowerTransformer)
    ends = _endsOfTransformer(store, pt)
    if length(ends) == 2
      _map2WTrafo!(net, store, topo, created, svmap, svsteps, baseMVA, pt, ends, ctx)
    elseif length(ends) == 3
      _map3WTrafo!(net, store, topo, created, svmap, svsteps, baseMVA, pt, ends, ctx)
    else
      push!(ctx.messages, "skip: PowerTransformer $(something(str(pt, :name), pt.mrid)) with $(length(ends)) ends")
    end
  end
end

function _mapShunts!(net, store, topo, created, svmap, ctx::_MapCtx)
  for sh in objectsOf(store, :LinearShuntCompensator)
    t = busOfEquipment(topo, sh, 1)
    name = something(str(sh, :name), sh.mrid)
    if t === nothing
      push!(ctx.messages, "skip: LinearShuntCompensator $(name) — no TopologicalNode")
      push!(ctx.skipped, sh.mrid)
      continue
    end
    if !_inService(ctx, sh)
      push!(ctx.messages, "skip: LinearShuntCompensator $(name) — out of service (SSH inService=false)")
      push!(ctx.skipped, sh.mrid)
      continue
    end
    if !_conn(ctx, t.connected)
      push!(ctx.messages, "skip: LinearShuntCompensator $(name) — terminal disconnected")
      push!(ctx.skipped, sh.mrid)
      continue
    end
    sections = something(num(sh, :sections), num(sh, :normalSections, 0.0))
    b = num(sh, :bPerSection, 0.0) * sections
    g = num(sh, :gPerSection, 0.0) * sections
    bus = _ensureBus!(net, created, topo, svmap, t.tn)
    # MATPOWER-style shunt: MW/MVAr injected at 1.0 pu of the bus vn
    Sparlectra.addShuntMatpower!(net = net, busName = bus, Gs = g * t.vn_kV^2, Bs = b * t.vn_kV^2)
  end
end

function _mapLoads!(net, store, topo, created, svmap, ctx::_MapCtx)
  for cls in (:EnergyConsumer, :ConformLoad, :NonConformLoad, :StationSupply)
    for load in objectsOf(store, cls)
      t = busOfEquipment(topo, load, 1)
      name = something(str(load, :name), load.mrid)
      if t === nothing
        push!(ctx.messages, "skip: $(cls) $(name) — no TopologicalNode")
        push!(ctx.skipped, load.mrid)
        continue
      end
      if !_inService(ctx, load)
        push!(ctx.messages, "skip: $(cls) $(name) — out of service (SSH inService=false)")
        push!(ctx.skipped, load.mrid)
        continue
      end
      if !_conn(ctx, t.connected)
        push!(ctx.messages, "skip: $(cls) $(name) — terminal disconnected")
        push!(ctx.skipped, load.mrid)
        continue
      end
      p = num(load, :p, 0.0)
      q = num(load, :q, 0.0)
      bus = _ensureBus!(net, created, topo, svmap, t.tn)
      Sparlectra.addProsumer!(net = net, busName = bus, type = "ENERGYCONSUMER", p = p, q = q, defer_bus_type_refresh = true)
    end
  end
end

# regulating-control voltage setpoint for a machine/injection.
# Returns `(vm_pu, voltage_controlled)`:
# - local target → (target/vn, true)
# - remote target (Stage 1 has no remote target-bus control) → (nothing, true);
#   the caller keeps the unit voltage-controlled at its own bus, holding the
#   bus's SV/start voltage instead of a stale fixed SSH q
# - no active voltage control → (nothing, false)
# Default plausibility band for a voltage RegulatingControl target, in p.u. of
# the regulated bus's nominal voltage. Deliberately far wider than any real
# setpoint: across every cached ENTSO-E/ReliCapGrid delivery (1187 targets) the
# actual range is 0.92 .. 1.15 p.u. The band exists to catch placeholder values,
# not to enforce an operating range.
#
# Svedala ships five generator nodes with targetValue = 0.001 kV (≈ 5e-5 p.u.)
# — the same five that carry no SvVoltage, i.e. units the exporting tool did not
# include in its solved state. Taken literally they produce a PV/reference bus
# at essentially zero volts, and the power flow then "converges" onto a
# collapsed all-zero solution that is formally correct and physically nonsense.
#
# These are DEFAULTS, not fixed law: the numbers come from the deliveries we
# have seen, so they are overridable per import (`vset_min_pu`/`vset_max_pu`,
# also as `cgmes_import.*` config keys) instead of silently imposing a policy on
# a delivery that legitimately regulates outside the band.
const CGMES_VSET_MIN_PU = 0.5
const CGMES_VSET_MAX_PU = 1.5

"""
Voltage setpoint of a regulating unit, as (vset_pu, is_voltage_regulating).

`vset_pu === nothing` together with `true` means "regulates voltage, but no
usable local setpoint" — the caller then holds the unit PV at its own bus
voltage (see `_pvVoltage`). That is the fallback for remote controls and for
implausible target values.
"""
function _voltageSetpoint(
  store::CGMESStore,
  topo::CGMESTopology,
  obj::CIMObject,
  own_tn::String,
  messages::Vector{String};
  vset_min_pu::Float64 = CGMES_VSET_MIN_PU,
  vset_max_pu::Float64 = CGMES_VSET_MAX_PU,
)::Tuple{Union{Nothing,Float64},Bool}
  rc = ref(store, obj, :RegulatingControl)
  rc === nothing && return (nothing, false)
  mode = something(enumval(rc, :mode), "")
  endswith(mode, ".voltage") || return (nothing, false)
  something(boolval(rc, :enabled), true) || return (nothing, false)
  target = num(rc, :targetValue)
  (target === nothing || target <= 0.0) && return (nothing, false)
  tmrid = get(rc.refs, :Terminal, nothing)
  tmrid === nothing && return (nothing, false)
  t = get(store.objects, tmrid, nothing)
  t === nothing && return (nothing, false)
  tn = tnOfTerminal(t)
  (tn === nothing || !haskey(topo.vn_kV, tn)) && return (nothing, false)
  if tn != own_tn
    push!(messages, "notice: $(something(str(obj, :name), obj.mrid)) has a remote voltage RegulatingControl (target bus $(get(topo.bus_name, tn, tn))) — held PV at its own bus in Stage 1")
    return (nothing, true)
  end
  vn = topo.vn_kV[tn]
  # A zero/absent nominal voltage makes the p.u. conversion meaningless.
  vn > 0.0 || return (nothing, true)
  vset = target / vn
  if vset < vset_min_pu || vset > vset_max_pu
    push!(
      messages,
      "warning: $(something(str(obj, :name), obj.mrid)) declares an implausible voltage target ($(target) kV = $(round(vset; sigdigits = 3)) pu at $(get(topo.bus_name, tn, tn)), nominal $(vn) kV) — ignored, unit held PV at the bus voltage derived from the nominal data",
    )
    return (nothing, true)
  end
  return (vset, true)
end

# limit hull over both sign-convention readings (see the machine qMin/qMax
# comment in _mapInjections!): [min(a,-b), max(b,-a)]
_limitHullMin(a::Union{Nothing,Float64}, b::Union{Nothing,Float64}) = (a === nothing || b === nothing) ? nothing : min(a, -b)
_limitHullMax(a::Union{Nothing,Float64}, b::Union{Nothing,Float64}) = (a === nothing || b === nothing) ? nothing : max(b, -a)

# PV voltage for a unit: local target if given, else the SV/start voltage of
# its own bus (remote-controlled units, slack fallback)
function _pvVoltage(net, bus::String, vset::Union{Nothing,Float64})::Float64
  vset !== nothing && return vset
  node = net.nodeVec[Sparlectra.geNetBusIdx(net = net, busName = bus)]
  return something(node._vm_pu, 1.0)
end

# slack candidate bookkeeping for decision D-2
#
# The tuples are laid out so that plain ascending `sort!` implements the ranking:
# priority ascending (lower referencePriority = stronger, 0.0 is reserved for the
# SV angle reference and therefore always wins), then -rating so the larger unit
# comes first, then bus name as a deterministic last resort.
mutable struct _SlackScan
  by_priority::Vector{Tuple{Float64,Float64,String}}   # (priority, -rating, bus)
  enis::Vector{Tuple{Float64,String}}                  # (maxP/rating, bus)
  machines::Vector{Tuple{Float64,String}}              # (ratedS, bus)
end
_SlackScan() = _SlackScan([], [], [])

"""
Pick the primary slack bus from the collected candidates.

Fallback chain, strongest evidence first:
 1. `referencePriority` markers — including the SV angle reference, which is
    injected at priority 0.0 because it is the reference the exporting tool
    actually used;
 2. the only / the largest `ExternalNetworkInjection` — a network equivalent is
    the natural reference;
 3. the largest `SynchronousMachine`.

Errors out when a delivery offers none of the three: an angle reference cannot
be invented.
"""
function _chooseSlack(scan::_SlackScan, messages::Vector{String})::String
  if !isempty(scan.by_priority)
    # Ascending sort ranks by the tuple layout documented at _SlackScan.
    sort!(scan.by_priority)
    bus = scan.by_priority[1][3]
    top = scan.by_priority[1][1]
    ties = count(t -> t[1] == top, scan.by_priority)
    # Priority 0.0 is the sentinel for "came from the SV angle reference", not
    # a value read out of the delivery.
    top == 0.0 && push!(messages, "slack: taken from the SV angle reference")
    if ties > 1 && top != 0.0
      # CGMES leaves the winner undefined when several units share the
      # highest referencePriority; we break the tie by rating, then by name,
      # which need not match the exporting tool's choice (RealGrid: two
      # machines at priority 1, its documentation names the other one).
      push!(messages, "notice: $(ties) units share referencePriority $(top) — slack chosen by largest rating; the source tool may have picked a different one")
    end
    push!(messages, "slack: referencePriority marker → $(bus)")
    return bus
  end
  # An ExternalNetworkInjection represents the rest of the grid, so it is the
  # natural reference when no explicit priority is given. The single-ENI case is
  # separated only to produce a clearer message.
  if length(scan.enis) == 1
    push!(messages, "slack: single ExternalNetworkInjection → $(scan.enis[1][2])")
    return scan.enis[1][2]
  end
  if length(scan.enis) > 1
    sort!(scan.enis; by = t -> -t[1])
    push!(messages, "slack: largest ExternalNetworkInjection → $(scan.enis[1][2])")
    return scan.enis[1][2]
  end
  isempty(scan.machines) && error("CGMES import: no slack candidate found (no referencePriority, no ExternalNetworkInjection, no SynchronousMachine)")
  sort!(scan.machines; by = t -> -t[1])
  push!(messages, "slack: largest SynchronousMachine → $(scan.machines[1][2])")
  return scan.machines[1][2]
end

"""
Decide the slack bus set. With `multi_slack` every **electrical** island — in
the sense of `electricalIslandComponents`, i.e. closed links count as
connections — receives at most one reference: the candidate with the highest
`referencePriority` wins, further SV angle references in the same island are
not promoted and keep the bus type their own injection produces (PV for a
regulating machine, PQ otherwise — their voltage setpoint comes from the SV
state anyway). Two references in one island would be physically wrong. On
clean data the guard is a no-op — Svedala's two references sit in two
electrical islands once `Equipment.inService` is honored — but it protects
against deliveries whose island declarations and switch states disagree.
"""
function _selectSlackBuses(net, ctx::_MapCtx, primary::String, angle_ref_buses::Set{String}, pending, multi_slack::Bool)::Set{String}
  chosen = Set{String}([primary])
  # Only angle references that actually carry an injection are candidates: a
  # reference node without a source cannot hold a voltage and would produce a
  # slack bus with nothing behind it.
  candidates = [b for b in angle_ref_buses if b != primary && any(inj -> inj.bus == b, pending)]
  isempty(candidates) && return chosen
  if !multi_slack
    push!(ctx.messages, "notice: the SV declares further angle reference(s) ($(join(sort(candidates), ", "))) — not promoted to slack; pass multi_slack = true to give separate electrical islands their own reference")
    return chosen
  end

  # One reference per electrical island. The electrical view is essential here:
  # counting only branches, retained switches would split one island into
  # several and we would hand out a second slack inside a single island.
  island_of = Sparlectra.electricalIslandOfBus(net)
  taken = Set{Int}()
  pi = get(island_of, get(net.busDict, primary, 0), 0)
  # Island 0 means "bus unknown/isolated" — nothing to reserve in that case.
  pi == 0 || push!(taken, pi)
  # Sorted so the outcome does not depend on Set iteration order.
  for b in sort(candidates)
    isl = get(island_of, get(net.busDict, b, 0), 0)
    if isl == 0 || isl in taken
      push!(ctx.messages, "notice: angle reference $(b) shares its electrical island with $(primary) — not promoted to slack (a second reference in one island would be physically wrong)")
      continue
    end
    push!(taken, isl)
    push!(chosen, b)
    push!(ctx.messages, "slack: additional reference $(b) for electrical island $(isl)")
  end
  return chosen
end

"""
Angle-reference nodes declared by the SV profile: `TopologicalIsland` names
one `AngleRefTopologicalNode` per island. That is the exporting tool's own
reference choice — the most reliable slack information a delivery can carry,
and the only one that scales to multi-island models (ReliCapGrid/Svedala
ships two islands).

!!! warning "Two different island notions"
    A CIM `TopologicalIsland` is the exporting tool's partitioning and need
    **not** agree with the AC island Sparlectra sees after evaluating the
    switch states — which is why `_selectSlackBuses` re-checks the references
    against `electricalIslandOfBus` instead of trusting the CIM island count.
    Instructive history: Svedala *appeared* to disagree (both references in
    one electrical island), but that was an importer artifact — switches with
    SSH `Equipment.inService = false` were treated as closed and merged the
    two islands. With `inService` honored, the electrical view agrees with the
    delivery's two `TopologicalIsland`s exactly. The cross-check stays: it
    costs little and catches deliveries where the disagreement is real.
"""
function _svAngleRefNodes(store::CGMESStore)::Set{String}
  out = Set{String}()
  for isl in objectsOf(store, :TopologicalIsland)
    tn = get(isl.refs, :AngleRefTopologicalNode, nothing)
    tn === nothing || push!(out, tn)
  end
  return out
end

# Do the EquivalentInjections collected on one boundary node describe the
# same exchange seen from both sides? True only when both sides are present
# AND their p values (approximately) cancel — the defining property of an
# assembled AC X-node pair. A non-cancelling pair is two different devices
# (HVDC converters; the residue is the DC loss) and must be kept. The
# tolerance only absorbs rounding in the delivery (largest observed residue of
# a genuine pair: 1.6 % / 0.025 MW; smallest genuine non-pair: 23.2 % / 22 MW).
function _eiPairCancels(ps::Vector{Float64})::Bool
  length(ps) >= 2 || return false
  pmax = maximum(abs, ps)
  return abs(sum(ps)) <= max(0.1 * pmax, 1.0)
end

# Boundary-node identity of an EquivalentInjection's terminal (see the
# ei_p_bnode comment in _mapInjections!): the terminal's ConnectivityNode if
# that CN comes from a boundary file (CGMES 3.0 style), else its
# TopologicalNode if that TN does (2.4.15 TP_BD style), else `nothing` — the
# EI does not sit on a boundary node at all (e.g. an internal equivalent).
function _eiBoundaryNode(store::CGMESStore, topo::CGMESTopology, ei::CIMObject)::Union{Nothing,String}
  ts = get(topo.terminals, ei.mrid, nothing)
  (ts === nothing || isempty(ts)) && return nothing
  t = ts[1]
  cn = get(t.refs, :ConnectivityNode, nothing)
  (cn !== nothing && cn in store.boundary) && return cn
  tn = tnOfTerminal(t)
  (tn !== nothing && tn in store.boundary) && return tn
  return nothing
end

function _mapInjections!(net, store, topo, created, svmap, ctx::_MapCtx; multi_slack::Bool = true)
  scan = _SlackScan()
  angle_refs = _svAngleRefNodes(store)
  angle_ref_buses = Set{String}(topo.bus_name[tn] for tn in angle_refs if haskey(topo.bus_name, tn))
  isempty(angle_ref_buses) || push!(ctx.messages, "notice: SV declares $(length(angle_ref_buses)) angle-reference node(s) — used as slack candidate(s): $(join(sort(collect(angle_ref_buses)), ", "))")
  pending = NamedTuple[]

  for sm in objectsOf(store, :SynchronousMachine)
    t = busOfEquipment(topo, sm, 1)
    name = something(str(sm, :name), sm.mrid)
    if t === nothing
      push!(ctx.messages, "skip: SynchronousMachine $(name) — no TopologicalNode")
      push!(ctx.skipped, sm.mrid)
      continue
    end
    if !_conn(ctx, t.connected)
      push!(ctx.messages, "skip: SynchronousMachine $(name) — terminal disconnected")
      push!(ctx.skipped, sm.mrid)
      continue
    end
    # ReliCapGrid parks out-of-service units exactly like this: connected=true,
    # inService=false, no SvVoltage, placeholder regulation target. Importing
    # them added 413 MW of phantom generation in Svedala alone.
    if !_inService(ctx, sm)
      push!(ctx.messages, "skip: SynchronousMachine $(name) — out of service (SSH inService=false)")
      push!(ctx.skipped, sm.mrid)
      continue
    end
    bus = _ensureBus!(net, created, topo, svmap, t.tn)
    gu = ref(store, sm, :GeneratingUnit)
    ratedS = num(sm, :ratedS, 0.0)
    # machines without usable scalar minQ/maxQ get wide symmetric limits so
    # voltage regulation is not spuriously Q-pinned. This covers both absent
    # values and the degenerate 0/0 pair that data sets write when the real
    # limits live in an InitialReactiveCapabilityCurve (not evaluated in
    # Stage 1) — MicroGrid BE-G1 is such a case.
    qwide = max(ratedS, net.baseMVA)
    qlo = _limitHullMin(num(sm, :minQ), num(sm, :maxQ))
    qhi = _limitHullMax(num(sm, :minQ), num(sm, :maxQ))
    if qlo === nothing || qhi === nothing || (qhi - qlo) < 1e-6
      qlo = -qwide
      qhi = qwide
      push!(ctx.messages, "notice: SynchronousMachine $(name) has no usable scalar minQ/maxQ — using ±$(qwide) MVAr (capability curves are not evaluated in Stage 1)")
    end
    prio = num(sm, :referencePriority, 0.0)
    # an SV angle reference outranks referencePriority: it is what the
    # exporting tool actually used
    bus in angle_ref_buses && push!(scan.by_priority, (0.0, -ratedS, bus))
    prio >= 1.0 && push!(scan.by_priority, (prio, -ratedS, bus))
    push!(scan.machines, (ratedS, bus))
    vset, vctrl = _voltageSetpoint(store, topo, sm, t.tn, ctx.messages; vset_min_pu = ctx.vset_min_pu, vset_max_pu = ctx.vset_max_pu)
    push!(
      pending,
      (
        bus = bus,
        type = "SYNCHRONOUSMACHINE",
        p = -num(sm, :p, 0.0),                       # CGMES machine sign: p < 0 = injection
        q = -num(sm, :q, 0.0),
        pMin = gu === nothing ? nothing : num(gu, :minOperatingP),
        pMax = gu === nothing ? nothing : num(gu, :maxOperatingP),
        # The ENTSO-E test sets are inconsistent about the sign convention of
        # machine minQ/maxQ (MicroGrid NL: load convention like p/q — the SV
        # q lies in [-maxQ, -minQ]; SmallGrid: injection convention). Stage 1
        # uses the hull of both readings so no unit is wrongly Q-pinned;
        # strict limit fidelity is a Stage-2 topic.
        qMin = qlo,
        qMax = qhi,
        vm_pu = vctrl ? _pvVoltage(net, bus, vset) : nothing,
        isRegulated = vctrl,
      ),
    )
  end

  for eni in objectsOf(store, :ExternalNetworkInjection)
    t = busOfEquipment(topo, eni, 1)
    name = something(str(eni, :name), eni.mrid)
    if t === nothing || !_conn(ctx, t.connected)
      push!(ctx.messages, "skip: ExternalNetworkInjection $(name) — missing/disconnected terminal")
      push!(ctx.skipped, eni.mrid)
      continue
    end
    # must run before the slack-candidate bookkeeping below: an out-of-service
    # network equivalent is no reference candidate
    if !_inService(ctx, eni)
      push!(ctx.messages, "skip: ExternalNetworkInjection $(name) — out of service (SSH inService=false)")
      push!(ctx.skipped, eni.mrid)
      continue
    end
    bus = _ensureBus!(net, created, topo, svmap, t.tn)
    prio = num(eni, :referencePriority, 0.0)
    maxP = abs(something(num(eni, :maxP), num(eni, :p, 0.0)))
    bus in angle_ref_buses && push!(scan.by_priority, (0.0, -maxP, bus))
    prio >= 1.0 && push!(scan.by_priority, (prio, -maxP, bus))
    push!(scan.enis, (maxP, bus))
    vset, vctrl = _voltageSetpoint(store, topo, eni, t.tn, ctx.messages; vset_min_pu = ctx.vset_min_pu, vset_max_pu = ctx.vset_max_pu)
    push!(
      pending,
      (
        bus = bus,
        type = "EXTERNALNETWORKINJECTION",
        p = -num(eni, :p, 0.0),                      # load convention → injection
        q = -num(eni, :q, 0.0),
        pMin = _limitHullMin(num(eni, :minP), num(eni, :maxP)),
        pMax = _limitHullMax(num(eni, :minP), num(eni, :maxP)),
        qMin = _limitHullMin(num(eni, :minQ), num(eni, :maxQ)),
        qMax = _limitHullMax(num(eni, :minQ), num(eni, :maxQ)),
        vm_pu = vctrl ? _pvVoltage(net, bus, vset) : nothing,
        isRegulated = vctrl,
      ),
    )
  end

  # EquivalentInjections per boundary node: an assembled model carries one EI
  # per SIDE on each X-node (e.g. BE-Inj-X… and NL-Inj-X…). Two or more EIs
  # on one boundary node therefore mean both sides are present and the
  # equivalents must not be applied on top of the real flows. A dangling
  # line keeps its single EI. (A count of real terminals is NOT a valid
  # criterion — SmallGrid has X-nodes with two parallel dangling circuits
  # and one EI.)
  #
  # Node identity is version-dependent. A CGMES 2.4.15 boundary set ships
  # TopologicalNodes (TP_BD), so both sides' EIs meet on one boundary TN. A
  # CGMES 3.0 boundary carries only ConnectivityNodes — each area's own TP
  # assigns its own TN to the shared CN, and the two sides need not even agree
  # on one TN (ReliCapGrid's Espheim–Portheim border: two TNs, one CN). The
  # shared identity is therefore the terminal's boundary CN when there is one,
  # the boundary TN otherwise.
  #
  # The p values are collected (not just counted) because presence of both
  # sides alone is not sufficient to discard: an assembled AC X-node pair is
  # two views of the SAME exchange and cancels exactly (measured: MicroGrid
  # 0.0 %, ReliCapGrid AC borders 0.0–1.6 % of max|p|), while the pair at a
  # DC crossing consists of the two CONVERTER injections whose difference is
  # the DC loss (ReliCapGrid BP_SD-EH_DC1: +94.57 / −72.60, 23.2 %). Those
  # must survive — they are devices, not a double-counted equivalent.
  ei_p_bnode = Dict{String,Vector{Float64}}()
  # StaticVarCompensator: a controllable reactive source with P = 0 — the
  # MATPOWER-compatible model (generator with Pg = 0). CGMES capacitiveRating /
  # inductiveRating are REACTANCES in Ohm (RealGrid: ±202.5 Ω @ 225 kV ↔
  # 250/-100 MVAr), so the Q limits come from Q = vn²/X. The droop `slope` is
  # deliberately ignored (standard practice), with a notice.
  for svc in objectsOf(store, :StaticVarCompensator)
    t = busOfEquipment(topo, svc, 1)
    name = something(str(svc, :name), svc.mrid)
    if t === nothing || !_conn(ctx, t.connected)
      push!(ctx.messages, "skip: StaticVarCompensator $(name) — missing/disconnected terminal")
      push!(ctx.skipped, svc.mrid)
      continue
    end
    if !_inService(ctx, svc)
      push!(ctx.messages, "skip: StaticVarCompensator $(name) — out of service (SSH inService=false)")
      push!(ctx.skipped, svc.mrid)
      continue
    end
    bus = _ensureBus!(net, created, topo, svmap, t.tn)
    qwide = 10.0 * net.baseMVA
    xcap = num(svc, :capacitiveRating)
    xind = num(svc, :inductiveRating)
    q_from_x(x) = (x === nothing || x == 0.0) ? nothing : t.vn_kV^2 / x
    qcap = q_from_x(xcap)
    qind = q_from_x(xind)
    qlo = qcap === nothing || qind === nothing ? -qwide : min(qcap, qind)
    qhi = qcap === nothing || qind === nothing ? qwide : max(qcap, qind)
    if max(abs(qlo), abs(qhi)) > qwide
      push!(ctx.messages, "notice: StaticVarCompensator $(name) — implausible Q range $(round(qlo; digits=1))..$(round(qhi; digits=1)) MVAr from ratings, clamped to ±$(qwide)")
      qlo = max(qlo, -qwide)
      qhi = min(qhi, qwide)
    end
    slope = num(svc, :slope, 0.0)
    slope != 0.0 && push!(ctx.messages, "notice: StaticVarCompensator $(name) — droop slope $(slope) ignored (constant-voltage model)")
    vset, vctrl = _voltageSetpoint(store, topo, svc, t.tn, ctx.messages; vset_min_pu = ctx.vset_min_pu, vset_max_pu = ctx.vset_max_pu)
    if !vctrl && endswith(something(enumval(svc, :sVCControlMode), ""), ".voltage") && something(boolval(svc, :controlEnabled), true)
      vsp = num(svc, :voltageSetPoint)
      (vsp !== nothing && vsp > 0.0) && ((vset, vctrl) = (vsp / t.vn_kV, true))
    end
    # Consistency guard against snapshots whose control flags contradict their
    # own solved state (RealGrid: RC enabled with target 1.0 pu while the SV
    # voltages sit at 1.008–1.031 and SSH q = 0 — the reference clearly did
    # not let the SVCs regulate; importing them as PV at 1.0 fights the whole
    # 225 kV profile and breaks convergence). If the SV start voltage of the
    # bus deviates from the target by more than 0.005 pu, fall back to PQ
    # with the SSH operating point.
    if vctrl && vset !== nothing
      node = net.nodeVec[Sparlectra.geNetBusIdx(net = net, busName = bus)]
      startvm = node._vm_pu
      if startvm !== nothing && abs(startvm - vset) > 0.005
        push!(ctx.messages, "notice: StaticVarCompensator $(name) — voltage target $(round(vset; digits=4)) pu contradicts the SV state ($(round(startvm; digits=4)) pu); imported as PQ with the SSH operating point")
        vctrl = false
      end
    end
    push!(
      pending,
      (
        bus = bus,
        type = "STATICVARCOMPENSATOR",
        p = 0.0,
        q = -num(svc, :q, 0.0),                      # load convention → injection
        pMin = 0.0,
        pMax = 0.0,
        qMin = qlo,
        qMax = qhi,
        vm_pu = vctrl ? _pvVoltage(net, bus, vset) : nothing,
        isRegulated = vctrl,
      ),
    )
  end

  # HVDC converter stations, Stage-0 model: each converter becomes a fixed
  # injection at its AC connection point with the SSH operating point — the
  # industry-standard load-flow treatment of a point-to-point link (MATPOWER's
  # dcline does the same with two bounded dummy generators). The one property
  # that matters for the power flow is that HVDC has NO angle coupling: the
  # transfer is a control setpoint, not the result of an angle difference. The
  # DC middle (DCLineSegment/DCNode) is deliberately NOT mapped, so two areas
  # joined only through HVDC correctly remain separate electrical islands with
  # their own angle references. The loss shows up as the difference between
  # the two converters' setpoints (ReliCapGrid GA–BH: −100.0 / +109.1 MW).
  for cls in (:VsConverter, :CsConverter)
    for conv in objectsOf(store, cls)
      t = busOfEquipment(topo, conv, 1)
      name = something(str(conv, :name), conv.mrid)
      if t === nothing || !_conn(ctx, t.connected)
        push!(ctx.messages, "skip: $(cls) $(name) — missing/disconnected AC terminal")
        push!(ctx.skipped, conv.mrid)
        continue
      end
      if !_inService(ctx, conv)
        push!(ctx.messages, "skip: $(cls) $(name) — out of service (SSH inService=false)")
        push!(ctx.skipped, conv.mrid)
        continue
      end
      p_ssh = num(conv, :p, 0.0)
      q_ssh = num(conv, :q, 0.0)
      bus = _ensureBus!(net, created, topo, svmap, t.tn)
      push!(
        pending,
        (
          bus = bus,
          type = "EXTERNALNETWORKINJECTION",
          p = -p_ssh,                              # load convention → injection
          q = -q_ssh,
          pMin = -p_ssh,
          pMax = -p_ssh,
          qMin = -q_ssh,
          qMax = -q_ssh,
          vm_pu = nothing,
          isRegulated = false,
        ),
      )
      push!(ctx.messages, "notice: $(cls) $(name) — Stage-0 HVDC: fixed PCC injection at $(bus) (p=$(round(-p_ssh; digits = 1)) MW, q=$(round(-q_ssh; digits = 1)) MVAr from SSH; DC side not mapped, no angle coupling)")
    end
  end

  # Collect on the raw terminals, not via busOfEquipment: an EI whose TN does
  # not resolve (3.0 split-TN border) must still be counted, otherwise its
  # counterpart on the other side sees a count of 1 and survives.
  for ei in objectsOf(store, :EquivalentInjection)
    # an out-of-service EI is no evidence that its side is present
    _inService(ctx, ei) || continue
    bnode = _eiBoundaryNode(store, topo, ei)
    bnode === nothing && continue
    push!(get!(Vector{Float64}, ei_p_bnode, bnode), num(ei, :p, 0.0))
  end

  for ei in objectsOf(store, :EquivalentInjection)
    t = busOfEquipment(topo, ei, 1)
    name = something(str(ei, :name), ei.mrid)
    if t === nothing
      push!(ctx.messages, "skip: EquivalentInjection $(name) — no TopologicalNode")
      push!(ctx.skipped, ei.mrid)
      continue
    end
    if !_conn(ctx, t.connected)
      push!(ctx.messages, "skip: EquivalentInjection $(name) — terminal disconnected")
      push!(ctx.skipped, ei.mrid)
      continue
    end
    if !_inService(ctx, ei)
      push!(ctx.messages, "skip: EquivalentInjection $(name) — out of service (SSH inService=false)")
      push!(ctx.skipped, ei.mrid)
      continue
    end
    bnode = _eiBoundaryNode(store, topo, ei)
    if bnode !== nothing && _eiPairCancels(get(ei_p_bnode, bnode, Float64[]))
      push!(ctx.messages, "skip: EquivalentInjection $(name) — assembled boundary node $(topo.bus_name[t.tn]) has both sides present")
      push!(ctx.skipped, ei.mrid)
      continue
    end
    bus = _busForTerminal!(net, created, topo, svmap, ctx, ei, 1, t.tn)
    # decision D-5: regulationStatus + positive target → PV-like boundary
    vset = nothing
    if something(boolval(ei, :regulationStatus), false)
      target = num(ei, :regulationTarget, 0.0)
      target > 0.0 && (vset = target / t.vn_kV)
    end
    push!(
      pending,
      (
        bus = bus,
        type = "EXTERNALNETWORKINJECTION",
        p = -num(ei, :p, 0.0),
        q = -num(ei, :q, 0.0),
        pMin = nothing,
        pMax = nothing,
        qMin = nothing,
        qMax = nothing,
        vm_pu = vset,
        isRegulated = vset !== nothing,
      ),
    )
  end

  slack_bus = _chooseSlack(scan, ctx.messages)
  # One slack per island: every SV angle-reference bus that carries a unit
  # becomes a reference, not just the primary one. A multi-island delivery
  # (ReliCapGrid/Svedala: two islands) is otherwise unsolvable — each island
  # needs its own angle reference.
  slack_buses = _selectSlackBuses(net, ctx, slack_bus, angle_ref_buses, pending, multi_slack)
  marked = Set{String}()
  for inj in pending
    is_slack = inj.bus in slack_buses && !(inj.bus in marked)
    is_slack && push!(marked, inj.bus)
    vm = inj.vm_pu
    if is_slack && vm === nothing
      # slack needs a voltage target: use the SV/start value of its bus
      node = net.nodeVec[Sparlectra.geNetBusIdx(net = net, busName = inj.bus)]
      vm = something(node._vm_pu, 1.0)
    end
    Sparlectra.addProsumer!(
      net = net,
      busName = inj.bus,
      type = inj.type,
      p = inj.p,
      q = inj.q,
      pMin = inj.pMin,
      pMax = inj.pMax,
      qMin = inj.qMin,
      qMax = inj.qMax,
      vm_pu = vm,
      referencePri = is_slack ? inj.bus : nothing,
      isRegulated = inj.isRegulated || is_slack,
      defer_bus_type_refresh = true,
    )
  end
  return slack_bus
end

# --- short-circuit harvest (§7.7 — read, not evaluated) ---------------------

function _busOrNothing(topo::CGMESTopology, obj::CIMObject)
  ts = get(topo.terminals, obj.mrid, nothing)
  ts === nothing && return nothing
  tn = tnOfTerminal(ts[1])
  tn === nothing && return nothing
  return get(topo.bus_name, tn, nothing)
end

function collectShortCircuitData(store::CGMESStore, topo::CGMESTopology)::CGMESShortCircuitData
  enis = NamedTuple[]
  for o in objectsOf(store, :ExternalNetworkInjection)
    push!(
      enis,
      (
        mrid = o.mrid,
        name = str(o, :name),
        bus = _busOrNothing(topo, o),
        maxInitialSymShCCurrent_A = num(o, :maxInitialSymShCCurrent),
        minInitialSymShCCurrent_A = num(o, :minInitialSymShCCurrent),
        maxR1ToX1Ratio = num(o, :maxR1ToX1Ratio),
        minR1ToX1Ratio = num(o, :minR1ToX1Ratio),
        maxR0ToX0Ratio = num(o, :maxR0ToX0Ratio),
        maxZ0ToZ1Ratio = num(o, :maxZ0ToZ1Ratio),
        ikSecond = boolval(o, :ikSecond),
        governorSCD = num(o, :governorSCD),
      ),
    )
  end
  machines = NamedTuple[]
  for o in objectsOf(store, :SynchronousMachine)
    push!(
      machines,
      (
        mrid = o.mrid,
        name = str(o, :name),
        bus = _busOrNothing(topo, o),
        satDirectSubtransX_pu = num(o, :satDirectSubtransX),
        satDirectTransX_pu = num(o, :satDirectTransX),
        r0_pu = num(o, :r0),
        x0_pu = num(o, :x0),
        r2_pu = num(o, :r2),
        x2_pu = num(o, :x2),
        earthing = boolval(o, :earthing),
        ratedS_MVA = num(o, :ratedS),
        ratedU_kV = num(o, :ratedU),
      ),
    )
  end
  lines = NamedTuple[]
  for o in objectsOf(store, :ACLineSegment)
    if num(o, :r0) !== nothing || num(o, :x0) !== nothing || num(o, :shortCircuitEndTemperature) !== nothing
      push!(
        lines,
        (mrid = o.mrid, name = str(o, :name), r0_ohm = num(o, :r0), x0_ohm = num(o, :x0), b0ch_S = num(o, :b0ch), g0ch_S = num(o, :g0ch), shortCircuitEndTemperature_C = num(o, :shortCircuitEndTemperature)),
      )
    end
  end
  ptes = NamedTuple[]
  for o in objectsOf(store, :PowerTransformerEnd)
    if num(o, :r0) !== nothing || num(o, :x0) !== nothing || boolval(o, :grounded) !== nothing
      push!(
        ptes,
        (mrid = o.mrid, name = str(o, :name), r0_ohm = num(o, :r0), x0_ohm = num(o, :x0), grounded = boolval(o, :grounded), rground_ohm = num(o, :rground), xground_ohm = num(o, :xground)),
      )
    end
  end
  eis = NamedTuple[]
  for o in objectsOf(store, :EquivalentInjection)
    if num(o, :r) !== nothing || num(o, :x) !== nothing || num(o, :r0) !== nothing
      push!(
        eis,
        (mrid = o.mrid, name = str(o, :name), bus = _busOrNothing(topo, o), r_ohm = num(o, :r), x_ohm = num(o, :x), r0_ohm = num(o, :r0), x0_ohm = num(o, :x0), r2_ohm = num(o, :r2), x2_ohm = num(o, :x2)),
      )
    end
  end
  return CGMESShortCircuitData(enis, machines, lines, ptes, eis)
end

# --- public API -------------------------------------------------------------

"""
    importCGMES(; path, baseMVA=100.0, require_boundary=true, name=nothing) -> CGMESImportResult

Import a CGMES delivery (EQ+SSH+TP, plus boundary set) into a Sparlectra
`Net` (Stage 1: bus-branch topology from the TP profile, fixed SSH tap
positions). Short-circuit source data is harvested into the result's
`shortcircuit` container (read, not evaluated). Start voltages are taken
from the SV profile where present.
"""
function importCGMES(;
  path,
  baseMVA::Float64 = 100.0,
  require_boundary::Bool = true,
  tap_control::Bool = false,
  ignore_connected::Bool = false,
  multi_slack::Bool = true,
  vset_min_pu::Float64 = CGMES_VSET_MIN_PU,
  vset_max_pu::Float64 = CGMES_VSET_MAX_PU,
  name::Union{Nothing,String} = nothing,
)::CGMESImportResult
  vset_min_pu <= vset_max_pu || throw(ArgumentError("importCGMES: vset_min_pu ($(vset_min_pu)) must not exceed vset_max_pu ($(vset_max_pu))"))
  store = loadCGMES(path)
  ctx = _MapCtx(tap_control, ignore_connected; vset_min_pu = vset_min_pu, vset_max_pu = vset_max_pu)
  ignore_connected && push!(ctx.messages, "notice: ignore_connected = true — SSH Terminal.connected flags are overridden, everything is treated as in service")
  messages = ctx.messages
  unresolved = unresolvedReferences(store)
  if !isempty(unresolved)
    topo_refs = count(u -> u.key in (:TopologicalNode, :ConnectivityNode), unresolved)
    if topo_refs > 0 && require_boundary
      error("CGMES import: $(topo_refs) unresolved topology references — boundary set missing? (pass the boundary ZIP/folder as additional path, or require_boundary=false)")
    end
    push!(ctx.messages, "notice: $(length(unresolved)) unresolved references (non-fatal)")
  end
  for f in store.files
    f.skipped && push!(ctx.messages, "file skipped: $(f.name) — $(f.skip_reason)")
  end

  topo = buildTopology(store)
  svmap = _svVoltageByTN(store)
  netname = something(name, "cgmes_import")
  net = Sparlectra.Net(name = netname, baseMVA = baseMVA)
  created = Set{String}()

  # must run before any element mapping: line endpoints and equivalent
  # injections at non-cancelling boundary nodes land on per-side buses
  _detectNonCancellingBoundarySides!(ctx, store, topo)

  _mapLines!(net, store, topo, created, svmap, baseMVA, ctx)
  # with tap control the SSH positions are the start point and the outer-loop
  # controllers find the solved positions themselves; without it the SvTapStep
  # positions reproduce the SV state directly (Stage 1)
  svsteps = tap_control ? Dict{String,Float64}() : _svTapPositions(store)
  _mapTransformers!(net, store, topo, created, svmap, svsteps, baseMVA, ctx)
  _mapSwitches!(net, store, topo, created, svmap, ctx)
  _mapShunts!(net, store, topo, created, svmap, ctx)
  _mapLoads!(net, store, topo, created, svmap, ctx)
  slack_bus = _mapInjections!(net, store, topo, created, svmap, ctx; multi_slack = multi_slack)
  Sparlectra.refreshBusTypesFromProsumers!(net)
  # Real deliveries carry TopologicalNodes whose only equipment is a load or a
  # switch that never becomes a branch (RealGrid: 198 such stubs). Mark them so
  # the solver skips them instead of reporting single-bus islands without a
  # reference; the bus stays in the model for provenance.
  _deEnergizeSourcelessComponents!(net, ctx)
  n_iso = length(Sparlectra.markIsolatedBuses!(net = net))
  n_iso > 0 && push!(ctx.messages, "notice: $(n_iso) isolated bus(es) without any in-service branch — marked isolated, excluded from the power flow")
  tap_control && _disableIneffectiveTapControllers!(net, ctx)

  sc = collectShortCircuitData(store, topo)
  isempty(sc.ac_line_segments) || push!(ctx.messages, "short-circuit data: $(length(sc.ac_line_segments)) lines (read, not evaluated)")
  isempty(sc.synchronous_machines) || push!(ctx.messages, "short-circuit data: $(length(sc.synchronous_machines)) machines (read, not evaluated)")
  isempty(sc.external_network_injections) || push!(ctx.messages, "short-circuit data: $(length(sc.external_network_injections)) external injections (read, not evaluated)")

  # Everything else is collected silently and inspected via result.messages.
  # "warning:" is reserved for data defects where the importer had to substitute
  # a value to keep the model usable — those must not stay hidden in a vector
  # nobody reads, so they also go to the log.
  _logCGMESWarnings(ctx.messages)

  return CGMESImportResult(net, store, topo, sc, slack_bus, ctx.messages, ctx.branch_side, ctx.skipped)
end

"""
Emit the collected `warning:` messages through the logging system.

The mapping functions only collect strings; a data defect that made the
importer substitute a derived value has to reach the user's log as well, not
just `result.messages`. Grouped into a single `@warn` so a delivery with many
defective objects does not flood the log.
"""
function _logCGMESWarnings(messages::Vector{String})
  warnings = [m for m in messages if startswith(m, "warning:")]
  isempty(warnings) && return nothing
  @warn "CGMES import: $(length(warnings)) data defect(s) — substituted values, see result.messages" details = join(warnings, "\n")
  return nothing
end

"""
    createNetFromCGMES(; path, baseMVA=100.0, require_boundary=true, name=nothing) -> Net

Thin wrapper around [`importCGMES`](@ref) returning only the `Net`.
"""
createNetFromCGMES(; kwargs...) = importCGMES(; kwargs...).net
