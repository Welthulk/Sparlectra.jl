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

# file: src/upfc_control.jl
# purpose: UPFC as a stationary quadrature composite (issue #325): one call
#          registers the SSSC series converter (line-flow target, injected-
#          voltage limit) and the STATCOM shunt converter (remote voltage
#          target, current-based limit) as one named device pair. This is a
#          registration helper on top of the existing controllers; it adds
#          no new solver mechanism and no new controller struct.

"""
    addUpfcControl!(net; fromBus, toBus, shunt_bus, target_bus, target_vm_pu,
                    p_target_mw, v_inj_max_pu, s_max_mva | i_max_ka, ...)
        -> (name, series, shunt)

Register a UPFC in the STATIONARY QUADRATURE model: an SSSC series converter
on the line branch `fromBus` to `toBus` (steering the branch active power to
`p_target_mw`, limited by the injectable series voltage `v_inj_max_pu`) plus
a STATCOM shunt converter at `shunt_bus` (holding the voltage at `target_bus`
to `target_vm_pu`, limited by the converter rating `s_max_mva` or `i_max_ka`),
paired under one composite name.

# Model and its honest limitation

A real UPFC couples its two converters through the DC-link active-power
balance. In the stationary model the series converter injects its voltage in
quadrature with the line current, so it exchanges (approximately) no active
power with the line, the DC link carries about zero, and the coupling
degenerates: what remains is exactly the SSSC plus STATCOM pair this function
registers. The composite therefore has NO series active-power injection: the
phase-shifter degree of freedom a real UPFC feeds through its DC link is out
of scope, and independent P and Q steering of the line via an injected
voltage of arbitrary phase stays unavailable. For stationary P steering
beyond the quadrature reach, the tap/phase-shift path on the same corridor
remains the answer. Theory in [FACTS Devices](@ref facts_devices).

# Arguments
- `net::Net`: the network.
- `fromBus::String`, `toBus::String`: terminals of the controlled line branch
  in the stored orientation (also the measurement direction of the flow).
- `shunt_bus::String`: bus carrying the shunt converter's machine; must be
  one of `fromBus`/`toBus` (a UPFC sits at one end of its own corridor).
- `target_bus::String`: the regulated (remote, PQ) bus of the shunt side.
- `target_vm_pu::Float64`: voltage target at `target_bus` in p.u.
- `p_target_mw::Float64`: series-side active-power target in MW,
  `fromBus` to `toBus` direction.
- `v_inj_max_pu::Float64`: maximum injectable series voltage in p.u.
  (the SSSC limit, must be positive).
- `s_max_mva`, `i_max_ka`: the shunt converter rating; exactly one of the
  two, positive (the STATCOM limit, see [`addMachineVoltageControl!`](@ref)).
- `deadband_vm_pu = 1e-3`, `deadband_p_mw = 0.5`: per-side convergence bands.
- `prosumer_index`: picks the machine when several generators sit at
  `shunt_bus`.
- `name`: composite name, default `UPFC_<fromBus>_<toBus>`; the two
  controllers are named `<name>_series` and `<name>_shunt`.
- `max_outer_iters::Int = 20`, `enabled::Bool = true`: forwarded to both
  sides.

# Behavior
Registration is all-or-nothing: composite-level validation runs before
anything is registered, and if the shunt-side registration fails after the
series side succeeded, the series controller is removed again and the error
rethrown. The two controllers keep their own result rows (one per actuator,
`at_limit` per converter side) with the device strings
`"UPFC series (VSC pair, stationary quadrature model)"` and
`"UPFC shunt (VSC pair, stationary quadrature model)"`.

Returns `(name = <composite>, series = <SeriesReactanceControl>,
shunt = <MachineVoltageControl>)`.
"""
function addUpfcControl!(
  net::Net;
  fromBus::String,
  toBus::String,
  shunt_bus::String,
  p_target_mw::Float64,
  v_inj_max_pu::Float64,
  model::Symbol = :quadrature,
  target_bus::Union{Nothing,String} = nothing,
  target_vm_pu::Union{Nothing,Float64} = nothing,
  q_target_mvar::Union{Nothing,Float64} = nothing,
  q_shunt_mvar::Float64 = 0.0,
  series_phase::Symbol = :free,
  s_max_mva::Union{Nothing,Float64} = nothing,
  i_max_ka::Union{Nothing,Float64} = nothing,
  deadband_vm_pu::Float64 = 1e-3,
  deadband_p_mw::Float64 = 0.5,
  deadband_q_mvar::Float64 = 0.5,
  prosumer_index::Union{Nothing,Int} = nothing,
  name::Union{Nothing,String} = nothing,
  max_outer_iters::Int = 20,
  enabled::Bool = true,
)
  model in (:quadrature, :full) || error("UpfcControl: model must be :quadrature (default, #325 SSSC+STATCOM composite) or :full (#326 DC-link-coupled model), got $(model).")
  # the full model owns an independent Q DOF; the quadrature composite has
  # none, so q_target_mvar is required for :full and rejected for :quadrature
  if model === :full
    q_target_mvar === nothing && error("UpfcControl: model = :full needs q_target_mvar (the independent reactive line-flow target); the quadrature composite has no independent Q degree of freedom.")
    # the full model's shunt runs on a reactive setpoint (q_shunt_mvar), not a
    # voltage target; closed-loop shunt voltage regulation is deferred (#217)
    (target_bus === nothing && target_vm_pu === nothing) || error("UpfcControl: model = :full uses q_shunt_mvar for the shunt reactive setpoint, not target_bus/target_vm_pu (closed-loop shunt voltage regulation is a follow-up); drop those keywords.")
    cname = something(name, string("UPFC_", fromBus, "_", toBus))
    upfc = addUpfcFullControl!(
      net;
      fromBus = fromBus,
      toBus = toBus,
      shunt_bus = shunt_bus,
      p_target_mw = p_target_mw,
      q_target_mvar = q_target_mvar,
      q_shunt_mvar = q_shunt_mvar,
      v_inj_max_pu = v_inj_max_pu,
      s_max_mva = s_max_mva,
      i_max_ka = i_max_ka,
      series_phase = series_phase,
      deadband_p_mw = deadband_p_mw,
      deadband_q_mvar = deadband_q_mvar,
      prosumer_index = prosumer_index,
      name = cname,
      max_outer_iters = max_outer_iters,
      enabled = enabled,
    )
    upfc.upfc_group = cname
    return (name = cname, upfc = upfc)
  end
  q_target_mvar === nothing || error("UpfcControl: q_target_mvar is only for model = :full; the quadrature composite steers one line quantity, not an independent P/Q pair.")
  # the quadrature composite's shunt is a remote voltage controller (STATCOM),
  # so it needs target_bus/target_vm_pu
  target_bus === nothing && error("UpfcControl: model = :quadrature needs target_bus (the shunt's remote voltage bus).")
  target_vm_pu === nothing && error("UpfcControl: model = :quadrature needs target_vm_pu (the shunt's voltage target).")
  # composite-level validation BEFORE any registration, so a rejected call
  # leaves the net untouched (no half-registered device)
  (s_max_mva === nothing) == (i_max_ka === nothing) && error("UpfcControl: pass exactly one shunt converter rating, s_max_mva or i_max_ka.")
  s_max_mva !== nothing && s_max_mva <= 0.0 && error("UpfcControl: s_max_mva must be positive, got $(s_max_mva)")
  i_max_ka !== nothing && i_max_ka <= 0.0 && error("UpfcControl: i_max_ka must be positive, got $(i_max_ka)")
  v_inj_max_pu > 0.0 || error("UpfcControl: v_inj_max_pu must be positive, got $(v_inj_max_pu)")
  shunt_bus in (fromBus, toBus) || error("UpfcControl: shunt_bus $(shunt_bus) is not an end of the corridor; a UPFC sits at one end of its own branch, admissible buses: $(fromBus), $(toBus).")
  # branch existence pre-check (detailed orientation/line-vs-trafo errors are
  # produced by addSeriesReactanceControl! below, still before anything is
  # registered; the series side registers first)
  if getNetBranch(net = net, fromBus = fromBus, toBus = toBus) === nothing && getNetBranch(net = net, fromBus = toBus, toBus = fromBus) === nothing
    error("UpfcControl: no branch between $(fromBus) and $(toBus).")
  end

  cname = something(name, string("UPFC_", fromBus, "_", toBus))

  # series converter first; its add function validates orientation, branch
  # kind, and the injected-voltage limit and throws before any registration
  series = addSeriesReactanceControl!(
    net;
    fromBus = fromBus,
    toBus = toBus,
    p_target_mw = p_target_mw,
    v_inj_max_pu = v_inj_max_pu,
    deadband_p_mw = deadband_p_mw,
    max_outer_iters = max_outer_iters,
    enabled = enabled,
    name = string(cname, "_series"),
  )
  # shunt converter second; on failure remove the series controller again so
  # the all-or-nothing rule holds
  local shunt::MachineVoltageControl
  try
    addMachineVoltageControl!(
      net;
      bus = shunt_bus,
      target_bus = target_bus,
      target_vm_pu = target_vm_pu,
      s_max_mva = s_max_mva,
      i_max_ka = i_max_ka,
      deadband_vm_pu = deadband_vm_pu,
      prosumer_index = prosumer_index,
      max_outer_iters = max_outer_iters,
      enabled = enabled,
      name = string(cname, "_shunt"),
    )
    shunt = last(net.machineControls)
  catch
    filter!(c -> c !== series, net.machineControls)
    rethrow()
  end

  series.upfc_group = cname
  shunt.upfc_group = cname
  return (name = cname, series = series, shunt = shunt)
end
