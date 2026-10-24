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

# file: src/adapters/scf/scf_controllers.jl
# purpose: controllers (FACTS and regulation) in the Sparlectra Case Format
#          (issue #342). PGM has no controller model at all, so they live in
#          the namespaced block, and they reuse the declarative
#          `control.controllers` schema VERBATIM: same type names, same
#          keyword names. The writer turns the instantiated controllers of a
#          net back into that schema, the reader hands the entries to
#          applyConfiguredControllers!, so there is exactly one construction
#          path and no second vocabulary to keep in sync.

# entries are written sorted by name so the file stays diff-stable
function _scf_controllers(net::Net)
  out = Dict{String,Any}()
  for ctrl in collect_outer_controllers(net)
    entry = _scf_controller_entry(ctrl)
    entry === nothing && continue
    name, payload = entry
    # the rank is written only when set, so files without ranks keep their bytes
    prio = control_priority(ctrl)
    prio == 0 || (payload["priority"] = prio)
    # The controller's own shunt has no shunt row (scf_id_map), so its
    # state (in service, model, estimated) travels here, in the key names
    # of a components.shunt_state row; the reader applies it after the
    # controller has rebuilt the shunt.
    if ctrl isa ShuntVoltageControl
      fields = _scf_shunt_state_fields(_find_controlled_shunt(net, ctrl))
      isempty(fields) || (payload["shunt_state"] = fields)
    end
    key = name
    n = 1
    while haskey(out, key)
      n += 1
      key = string(name, "_", n)
    end
    out[key] = payload
  end
  return out
end

# nothing for controller kinds the declarative schema cannot express (none
# today; the fallback keeps a future controller from silently corrupting a
# case file)
function _scf_controller_entry(ctrl)
  if ctrl isa PowerTransformerControl
    d = Dict{String,Any}("type" => "power_transformer", "trafo" => ctrl.trafo, "mode" => String(ctrl.mode))
    isempty(ctrl.followers) || (d["followers"] = String[f for f in ctrl.followers])
    ctrl.target_bus === nothing || (d["target_bus"] = ctrl.target_bus)
    ctrl.target_branch === nothing || (d["target_branch"] = String[ctrl.target_branch[1], ctrl.target_branch[2]])
    ctrl.target_vm_pu === nothing || (d["target_vm_pu"] = ctrl.target_vm_pu)
    ctrl.p_target_mw === nothing || (d["p_target_mw"] = ctrl.p_target_mw)
    ctrl.q_target_mvar === nothing || (d["q_target_mvar"] = ctrl.q_target_mvar)
    d["control_ratio"] = ctrl.control_ratio
    d["control_phase"] = ctrl.control_phase
    d["is_discrete"] = ctrl.is_discrete
    d["deadband_vm_pu"] = ctrl.deadband_vm_pu
    d["deadband_p_mw"] = ctrl.deadband_p_mw
    d["voltage_error_metric"] = String(ctrl.voltage_error_metric)
    d["max_outer_iters"] = ctrl.max_outer_iters
    d["enabled"] = ctrl.enabled
    # the control range only when set, so files without one keep their bytes
    ctrl.tap_min === nothing || (d["tap_min"] = ctrl.tap_min)
    ctrl.tap_max === nothing || (d["tap_max"] = ctrl.tap_max)
    return (string("tap_", ctrl.trafo), d)
  elseif ctrl isa MachineVoltageControl
    d = Dict{String,Any}("type" => "machine_voltage", "bus" => ctrl.bus, "target_bus" => ctrl.target_bus, "target_vm_pu" => ctrl.target_vm_pu, "deadband_vm_pu" => ctrl.deadband_vm_pu, "max_outer_iters" => ctrl.max_outer_iters, "enabled" => ctrl.enabled, "name" => ctrl.name)
    # :current mode is the STATCOM form (rating instead of a fixed Q box)
    if ctrl.limit_mode === :current && ctrl.s_max_mva !== nothing
      d["s_max_mva"] = ctrl.s_max_mva
    else
      d["qmin_mvar"] = ctrl.qmin_mvar
      d["qmax_mvar"] = ctrl.qmax_mvar
    end
    return (ctrl.name, d)
  elseif ctrl isa ShuntVoltageControl
    d = Dict{String,Any}("type" => "shunt_voltage", "bus" => ctrl.bus, "target_vm_pu" => ctrl.target_vm_pu, "bs_min_mvar" => ctrl.bs_min_mvar, "bs_max_mvar" => ctrl.bs_max_mvar, "deadband_vm_pu" => ctrl.deadband_vm_pu, "max_outer_iters" => ctrl.max_outer_iters, "enabled" => ctrl.enabled, "name" => ctrl.name)
    ctrl.step_mvar === nothing || (d["step_mvar"] = ctrl.step_mvar)
    # the controller's shunt is not written as a shunt row (scf_id_map), so
    # its present susceptance travels as the start value; zero is the
    # default and stays unwritten
    bs = _shunt_start_value(ctrl.bs_mvar, ctrl.bs_min_mvar, ctrl.bs_max_mvar, ctrl.step_mvar)
    bs == 0.0 || (d["bs_start_mvar"] = bs)
    return (ctrl.name, d)
  elseif ctrl isa SeriesReactanceControl
    d = Dict{String,Any}("type" => "series_reactance", "from_bus" => ctrl.fromBus, "to_bus" => ctrl.toBus, "p_target_mw" => ctrl.p_target_mw, "deadband_p_mw" => ctrl.deadband_p_mw, "max_outer_iters" => ctrl.max_outer_iters, "enabled" => ctrl.enabled, "name" => ctrl.name)
    # the two limit forms are exclusive: TCSC window versus SSSC voltage
    if ctrl.limit_mode === :injected_voltage && ctrl.v_inj_max_pu !== nothing
      d["v_inj_max_pu"] = ctrl.v_inj_max_pu
    else
      d["x_min_pu"] = ctrl.x_min_pu
      d["x_max_pu"] = ctrl.x_max_pu
    end
    return (ctrl.name, d)
  elseif ctrl isa HvdcPairControl
    d = Dict{String,Any}("type" => "hvdc_pair", "from_bus" => ctrl.from_bus, "to_bus" => ctrl.to_bus, "p_transfer_mw" => ctrl.p_transfer_mw, "loss_mw" => ctrl.loss_mw, "loss_fraction" => ctrl.loss_fraction, "deadband_vm_pu" => ctrl.deadband_vm_pu, "max_outer_iters" => ctrl.max_outer_iters, "enabled" => ctrl.enabled, "name" => ctrl.name)
    for (key, value) in ("p_rating_mw" => ctrl.p_rating_mw, "from_q_mvar" => ctrl.from_q_mvar, "to_q_mvar" => ctrl.to_q_mvar, "from_vset_pu" => ctrl.from_vset_pu, "to_vset_pu" => ctrl.to_vset_pu)
      value === nothing || (d[key] = value)
    end
    # an unlimited converter Q band is stated by ABSENCE, the same rule the
    # prosumer export documents: JSON has no Inf, and reading a missing key
    # back yields the same unlimited band (found on sp_case188, whose pair
    # runs with default bands)
    for (key, value) in ("from_qmin_mvar" => ctrl.from_qmin_mvar, "from_qmax_mvar" => ctrl.from_qmax_mvar, "to_qmin_mvar" => ctrl.to_qmin_mvar, "to_qmax_mvar" => ctrl.to_qmax_mvar)
      (value === nothing || !isfinite(value)) || (d[key] = value)
    end
    return (ctrl.name, d)
  elseif ctrl isa UpfcFullControl
    d = Dict{String,Any}(
      "type" => "upfc",
      "from_bus" => ctrl.fromBus,
      "to_bus" => ctrl.toBus,
      "shunt_bus" => ctrl.shunt_bus,
      "p_target_mw" => ctrl.p_target_mw,
      "q_target_mvar" => ctrl.q_target_mvar,
      "q_shunt_mvar" => ctrl.q_shunt_mvar,
      "v_inj_max_pu" => ctrl.v_inj_max_pu,
      "s_max_mva" => ctrl.s_max_mva,
      "series_phase" => String(ctrl.series_phase),
      "deadband_p_mw" => ctrl.deadband_p_mw,
      "deadband_q_mvar" => ctrl.deadband_q_mvar,
      "max_outer_iters" => ctrl.max_outer_iters,
      "enabled" => ctrl.enabled,
      "name" => ctrl.name,
      "model" => "full",
    )
    return (ctrl.name, d)
  end
  @warn "SCF export: controller type has no case-format mapping and is not written" type = typeof(ctrl)
  return nothing
end

"""
    _scf_apply_controllers!(net, entries; branch_by_name)

Instantiate the case file's controllers on `net`. The entries are the
declarative `control.controllers` schema, so they go through
`applyConfiguredControllers!` unchanged: one construction path, one
validation, one vocabulary. `branch_by_name` maps the file's branch names
(extra block) to built branch indices, so a `trafo`/`followers` reference
resolves to the branch of that id even where several branches share the
name. Other references (index strings, component ids) are passed through
untouched.
"""
function _scf_apply_controllers!(net::Net, entries::AbstractDict; branch_by_name::AbstractDict{String,Int} = Dict{String,Int}())
  isempty(entries) && return 0
  translate = ref -> begin
    s = String(ref)
    haskey(branch_by_name, s) ? string(branch_by_name[s]) : s
  end
  # the file stores the YAML mapping form (name => entry); ControlConfig
  # holds the entry list, so the key becomes the entry's name
  list = Any[]
  # shunt state of a shunt voltage controller's own shunt, by bus (one
  # such controller per bus); not a controller keyword, so it is taken out
  # before the entry goes to applyConfiguredControllers!
  shunt_states = Pair{String,SCFShuntStateRow}[]
  for key in sort!(String[String(k) for k in keys(entries)])
    body = entries[key]
    body isa AbstractDict || throw(ArgumentError("SCF: sparlectra.components.controllers.$(key) must be an object with a type key."))
    entry = Dict{String,Any}(String(k) => v for (k, v) in body)
    st = pop!(entry, "shunt_state", nothing)
    if st !== nothing
      get(entry, "type", nothing) == "shunt_voltage" || throw(ArgumentError("SCF: sparlectra.components.controllers.$(key).shunt_state is only valid for a shunt_voltage controller."))
      st isa AbstractDict || throw(ArgumentError("SCF: sparlectra.components.controllers.$(key).shunt_state must be an object."))
      for k in keys(st)
        String(k) in ("status", "model", "estimate") || throw(ArgumentError("SCF: sparlectra.components.controllers.$(key).shunt_state has the unknown key $(repr(String(k))); allowed are status, model, estimate."))
      end
      row = Dict{String,Any}(String(k) => v for (k, v) in st)
      row["shunt"] = 0
      push!(shunt_states, string(get(entry, "bus", "")) => _scfcase_shunt_state(row))
    end
    haskey(entry, "name") || (entry["name"] = key)
    haskey(entry, "trafo") && (entry["trafo"] = translate(entry["trafo"]))
    if haskey(entry, "followers") && entry["followers"] isa AbstractVector
      entry["followers"] = [translate(f) for f in entry["followers"]]
    end
    push!(list, entry)
  end
  n = applyConfiguredControllers!(net, ControlConfig(enabled = true, controllers = list))
  for (bus, row) in shunt_states
    ctrl = only(c for c in _shunt_controllers(net) if c.bus == bus)
    _scf_apply_shunt_state!(_find_controlled_shunt(net, ctrl), row)
  end
  return n
end
