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

# file: src/controller_config.jl
# purpose: instantiate outer-loop controllers declared under
#          control.controllers (issue #305). The YAML schema is a NAMED
#          MAPPING (one child mapping per controller, the key is the
#          controller name) because the minimal YAML reader deliberately
#          has no block-sequence support; the free-form child keys follow
#          the power_flow.distributed_slack.weights precedent. Entries are
#          structurally validated at config-load time; network references
#          are validated at apply time by the device add* functions.

# One row per supported type: the discriminator string, the keyword names
# the corresponding add* function accepts (beyond net), the required
# subset, and per-key value converters. Keeping this table next to the
# dispatcher makes "unknown key" and "missing key" errors self-updating.
const _CONTROLLER_TYPE_SPECS = Dict{String,NamedTuple}(
  "power_transformer" => (
    required = ("trafo", "mode"),
    optional = ("target_bus", "target_branch", "target_vm_pu", "p_target_mw", "q_target_mvar", "control_ratio", "control_phase", "is_discrete", "deadband_vm_pu", "deadband_p_mw", "voltage_error_metric", "max_outer_iters", "enabled"),
    symbols = ("mode", "voltage_error_metric"),
    supports_name = false,
  ),
  "machine_voltage" => (
    # s_max_mva / i_max_ka switch the controller into STATCOM mode (#297
    # Draft A); the add function enforces exclusivity with qmin/qmax
    required = ("bus", "target_bus", "target_vm_pu"),
    optional = ("qmin_mvar", "qmax_mvar", "s_max_mva", "i_max_ka", "deadband_vm_pu", "prosumer_index", "max_outer_iters", "enabled", "name"),
    symbols = (),
    supports_name = true,
  ),
  "shunt_voltage" => (
    required = ("bus", "target_vm_pu", "bs_min_mvar", "bs_max_mvar"),
    optional = ("bs_start_mvar", "deadband_vm_pu", "max_outer_iters", "enabled", "name"),
    symbols = (),
    supports_name = true,
  ),
  "series_reactance" => (
    # x_min_pu/x_max_pu (TCSC window) or v_inj_max_pu (SSSC, #297 Draft F);
    # the add function enforces exactly one of the two limit forms
    required = ("from_bus", "to_bus", "p_target_mw"),
    optional = ("x_min_pu", "x_max_pu", "v_inj_max_pu", "deadband_p_mw", "max_outer_iters", "enabled", "name"),
    symbols = (),
    supports_name = true,
  ),
  "hvdc_pair" => (
    # p_transfer_mw is required by the device function in setpoint mode and
    # forbidden in island_feed, so the spec keeps it optional and lets
    # addHvdcPairControl! enforce the mode-dependent rule
    required = ("from_bus", "to_bus"),
    optional = ("mode", "p_transfer_mw", "deadband_p_mw", "loss_mw", "loss_fraction", "p_rating_mw", "from_q_mvar", "to_q_mvar", "from_vset_pu", "to_vset_pu", "from_qmin_mvar", "from_qmax_mvar", "to_qmin_mvar", "to_qmax_mvar", "deadband_vm_pu", "max_outer_iters", "enabled", "name", "from_prosumer", "to_prosumer"),
    symbols = (),
    supports_name = true,
  ),
)

# Normalize the raw control.controllers value into a Vector{Dict{String,Any}}
# with a "name" key per entry. Accepts the canonical named mapping, an
# already-normalized vector of mappings (programmatic construction), and the
# empty placeholders the minimal YAML reader produces ("{}", "", nothing, []).
function _normalize_controller_entries(raw)::Vector{Dict{String,Any}}
  (raw === nothing || (raw isa AbstractString && strip(raw) in ("", "{}", "[]"))) && return Dict{String,Any}[]
  entries = Dict{String,Any}[]
  if raw isa AbstractDict
    for name in sort!(collect(keys(raw)); by = string)
      body = raw[name]
      body isa AbstractDict || throw(ArgumentError("control.controllers.$(name) must be a mapping with a type key."))
      entry = Dict{String,Any}(string(k) => v for (k, v) in body)
      haskey(entry, "name") || (entry["name"] = string(name))
      push!(entries, entry)
    end
    return entries
  end
  if raw isa AbstractVector
    for (i, body) in enumerate(raw)
      body isa AbstractDict || throw(ArgumentError("control.controllers[$(i)] must be a mapping with a type key."))
      push!(entries, Dict{String,Any}(string(k) => v for (k, v) in body))
    end
    return entries
  end
  throw(ArgumentError("control.controllers must be a mapping of named controller entries (block style), got $(typeof(raw))."))
end

# Structural validation at config-load time: known type, no unknown keys,
# required keys present. Reference/limit validation stays in the add*
# functions, which see the net.
function _validate_controller_entries(entries::Vector{Dict{String,Any}})
  for entry in entries
    label = haskey(entry, "name") ? string("control.controllers.", entry["name"]) : "control.controllers entry"
    haskey(entry, "type") || throw(ArgumentError("$(label): missing the type key; supported types: $(join(sort!(collect(keys(_CONTROLLER_TYPE_SPECS))), ", "))."))
    typ = string(entry["type"])
    spec = get(_CONTROLLER_TYPE_SPECS, typ, nothing)
    spec === nothing && throw(ArgumentError("$(label): unknown controller type \"$(typ)\"; supported types: $(join(sort!(collect(keys(_CONTROLLER_TYPE_SPECS))), ", "))."))
    allowed = ("type", "name", spec.required..., spec.optional...)
    for key in keys(entry)
      key in allowed || throw(ArgumentError("$(label) (type=$(typ)): unknown key \"$(key)\"; allowed keys: $(join(sort!(collect(setdiff(Set(allowed), Set(("type",))))), ", "))."))
    end
    for key in spec.required
      haskey(entry, key) || throw(ArgumentError("$(label) (type=$(typ)): required key \"$(key)\" is missing."))
    end
  end
  return nothing
end

_controller_cfg_float(v, label) = v isa Real ? Float64(v) : throw(ArgumentError("$(label) must be a number, got $(typeof(v))"))
_controller_cfg_int(v, label) = v isa Integer ? Int(v) : throw(ArgumentError("$(label) must be an integer, got $(typeof(v))"))
_controller_cfg_bool(v, label) = v isa Bool ? v : throw(ArgumentError("$(label) must be true or false, got $(typeof(v))"))

# Identity used for the skip-if-present rule: applying the same declaration
# twice to one net must be a no-op (run_sparlectra(net = ...) may execute
# repeatedly on the same Net), while a CONFLICTING second controller still
# fails loudly inside the add* function.
function _configured_controller_exists(net::Net, typ::String, entry::Dict{String,Any})::Bool
  if typ == "machine_voltage"
    return any(c -> c isa MachineVoltageControl && c.bus == string(entry["bus"]), net.machineControls)
  elseif typ == "shunt_voltage"
    return any(c -> c isa ShuntVoltageControl && c.bus == string(entry["bus"]), net.machineControls)
  elseif typ == "series_reactance"
    # NOTE: the struct fields are fromBus/toBus (bugfix: from_bus/to_bus
    # threw a FieldError here, so the idempotency check crashed instead of
    # skipping on the second apply of a series_reactance entry)
    return any(c -> c isa SeriesReactanceControl && c.fromBus == string(entry["from_bus"]) && c.toBus == string(entry["to_bus"]), net.machineControls)
  elseif typ == "hvdc_pair"
    return any(c -> c isa HvdcPairControl && c.from_bus == string(entry["from_bus"]) && c.to_bus == string(entry["to_bus"]), net.machineControls)
  elseif typ == "power_transformer"
    return any(c -> c.trafo == string(entry["trafo"]), _tap_controllers(net))
  end
  return false
end

"""
    applyConfiguredControllers!(net::Net, control_cfg::ControlConfig) -> Int

Instantiate the outer-loop controllers declared under `control.controllers`
onto `net` by calling the matching device function
(`addPowerTransformerControl!`, `addMachineVoltageControl!`,
`addShuntVoltageControl!`, `addSeriesReactanceControl!`). Returns the number
of controllers added.

Entries whose controlled element already carries a controller of the same
type are skipped, so repeated runs on the same `Net` stay idempotent; edit
programmatically or rebuild the net to change an already-applied
controller. Structural errors (unknown type or key, missing required key)
and device errors (unknown bus/branch/transformer, invalid limits) throw an
`ArgumentError` naming the entry.

The run pipeline calls this automatically before the outer control loop
when `control.enabled` is true and entries exist; call it directly to apply
a configuration to a programmatically built net.
"""
function applyConfiguredControllers!(net::Net, control_cfg::ControlConfig)::Int
  entries = _normalize_controller_entries(control_cfg.controllers)
  isempty(entries) && return 0
  _validate_controller_entries(entries)
  added = 0
  for entry in entries
    typ = string(entry["type"])
    name = get(entry, "name", nothing)
    label = name === nothing ? "control.controllers entry (type=$(typ))" : "control.controllers.$(name) (type=$(typ))"
    if _configured_controller_exists(net, typ, entry)
      @debug "configured controller already present, skipped" label
      continue
    end
    try
      if typ == "power_transformer"
        target_branch = if haskey(entry, "target_branch")
          tb = entry["target_branch"]
          (tb isa AbstractVector && length(tb) == 2) || throw(ArgumentError("target_branch must be a two-element list [from_bus, to_bus]"))
          (string(tb[1]), string(tb[2]))
        else
          nothing
        end
        addPowerTransformerControl!(
          net;
          trafo = string(entry["trafo"]),
          mode = Symbol(string(entry["mode"])),
          target_bus = haskey(entry, "target_bus") ? string(entry["target_bus"]) : nothing,
          target_branch = target_branch,
          target_vm_pu = haskey(entry, "target_vm_pu") ? _controller_cfg_float(entry["target_vm_pu"], "target_vm_pu") : nothing,
          p_target_mw = haskey(entry, "p_target_mw") ? _controller_cfg_float(entry["p_target_mw"], "p_target_mw") : nothing,
          q_target_mvar = haskey(entry, "q_target_mvar") ? _controller_cfg_float(entry["q_target_mvar"], "q_target_mvar") : nothing,
          control_ratio = haskey(entry, "control_ratio") ? _controller_cfg_bool(entry["control_ratio"], "control_ratio") : true,
          control_phase = haskey(entry, "control_phase") ? _controller_cfg_bool(entry["control_phase"], "control_phase") : false,
          is_discrete = haskey(entry, "is_discrete") ? _controller_cfg_bool(entry["is_discrete"], "is_discrete") : true,
          deadband_vm_pu = haskey(entry, "deadband_vm_pu") ? _controller_cfg_float(entry["deadband_vm_pu"], "deadband_vm_pu") : 1e-3,
          deadband_p_mw = haskey(entry, "deadband_p_mw") ? _controller_cfg_float(entry["deadband_p_mw"], "deadband_p_mw") : 0.5,
          voltage_error_metric = haskey(entry, "voltage_error_metric") ? Symbol(string(entry["voltage_error_metric"])) : :vm,
          max_outer_iters = haskey(entry, "max_outer_iters") ? _controller_cfg_int(entry["max_outer_iters"], "max_outer_iters") : 20,
          enabled = haskey(entry, "enabled") ? _controller_cfg_bool(entry["enabled"], "enabled") : true,
        )
      elseif typ == "machine_voltage"
        addMachineVoltageControl!(
          net;
          bus = string(entry["bus"]),
          target_bus = string(entry["target_bus"]),
          target_vm_pu = _controller_cfg_float(entry["target_vm_pu"], "target_vm_pu"),
          qmin_mvar = haskey(entry, "qmin_mvar") ? _controller_cfg_float(entry["qmin_mvar"], "qmin_mvar") : nothing,
          qmax_mvar = haskey(entry, "qmax_mvar") ? _controller_cfg_float(entry["qmax_mvar"], "qmax_mvar") : nothing,
          s_max_mva = haskey(entry, "s_max_mva") ? _controller_cfg_float(entry["s_max_mva"], "s_max_mva") : nothing,
          i_max_ka = haskey(entry, "i_max_ka") ? _controller_cfg_float(entry["i_max_ka"], "i_max_ka") : nothing,
          deadband_vm_pu = haskey(entry, "deadband_vm_pu") ? _controller_cfg_float(entry["deadband_vm_pu"], "deadband_vm_pu") : 1e-3,
          prosumer_index = haskey(entry, "prosumer_index") ? _controller_cfg_int(entry["prosumer_index"], "prosumer_index") : nothing,
          name = haskey(entry, "name") ? string(entry["name"]) : nothing,
          max_outer_iters = haskey(entry, "max_outer_iters") ? _controller_cfg_int(entry["max_outer_iters"], "max_outer_iters") : 20,
          enabled = haskey(entry, "enabled") ? _controller_cfg_bool(entry["enabled"], "enabled") : true,
        )
      elseif typ == "shunt_voltage"
        addShuntVoltageControl!(
          net;
          bus = string(entry["bus"]),
          target_vm_pu = _controller_cfg_float(entry["target_vm_pu"], "target_vm_pu"),
          bs_min_mvar = _controller_cfg_float(entry["bs_min_mvar"], "bs_min_mvar"),
          bs_max_mvar = _controller_cfg_float(entry["bs_max_mvar"], "bs_max_mvar"),
          bs_start_mvar = haskey(entry, "bs_start_mvar") ? _controller_cfg_float(entry["bs_start_mvar"], "bs_start_mvar") : 0.0,
          deadband_vm_pu = haskey(entry, "deadband_vm_pu") ? _controller_cfg_float(entry["deadband_vm_pu"], "deadband_vm_pu") : 1e-3,
          max_outer_iters = haskey(entry, "max_outer_iters") ? _controller_cfg_int(entry["max_outer_iters"], "max_outer_iters") : 20,
          enabled = haskey(entry, "enabled") ? _controller_cfg_bool(entry["enabled"], "enabled") : true,
          name = haskey(entry, "name") ? string(entry["name"]) : nothing,
        )
      elseif typ == "series_reactance"
        addSeriesReactanceControl!(
          net;
          fromBus = string(entry["from_bus"]),
          toBus = string(entry["to_bus"]),
          p_target_mw = _controller_cfg_float(entry["p_target_mw"], "p_target_mw"),
          x_min_pu = haskey(entry, "x_min_pu") ? _controller_cfg_float(entry["x_min_pu"], "x_min_pu") : nothing,
          x_max_pu = haskey(entry, "x_max_pu") ? _controller_cfg_float(entry["x_max_pu"], "x_max_pu") : nothing,
          v_inj_max_pu = haskey(entry, "v_inj_max_pu") ? _controller_cfg_float(entry["v_inj_max_pu"], "v_inj_max_pu") : nothing,
          deadband_p_mw = haskey(entry, "deadband_p_mw") ? _controller_cfg_float(entry["deadband_p_mw"], "deadband_p_mw") : 0.5,
          max_outer_iters = haskey(entry, "max_outer_iters") ? _controller_cfg_int(entry["max_outer_iters"], "max_outer_iters") : 20,
          enabled = haskey(entry, "enabled") ? _controller_cfg_bool(entry["enabled"], "enabled") : true,
          name = haskey(entry, "name") ? string(entry["name"]) : nothing,
        )
      elseif typ == "hvdc_pair"
        addHvdcPairControl!(
          net;
          from_bus = string(entry["from_bus"]),
          to_bus = string(entry["to_bus"]),
          mode = haskey(entry, "mode") ? Symbol(string(entry["mode"])) : :setpoint,
          p_transfer_mw = haskey(entry, "p_transfer_mw") ? _controller_cfg_float(entry["p_transfer_mw"], "p_transfer_mw") : nothing,
          deadband_p_mw = haskey(entry, "deadband_p_mw") ? _controller_cfg_float(entry["deadband_p_mw"], "deadband_p_mw") : 1e-3,
          loss_mw = haskey(entry, "loss_mw") ? _controller_cfg_float(entry["loss_mw"], "loss_mw") : 0.0,
          loss_fraction = haskey(entry, "loss_fraction") ? _controller_cfg_float(entry["loss_fraction"], "loss_fraction") : 0.0,
          p_rating_mw = haskey(entry, "p_rating_mw") ? _controller_cfg_float(entry["p_rating_mw"], "p_rating_mw") : nothing,
          from_q_mvar = haskey(entry, "from_q_mvar") ? _controller_cfg_float(entry["from_q_mvar"], "from_q_mvar") : nothing,
          to_q_mvar = haskey(entry, "to_q_mvar") ? _controller_cfg_float(entry["to_q_mvar"], "to_q_mvar") : nothing,
          from_vset_pu = haskey(entry, "from_vset_pu") ? _controller_cfg_float(entry["from_vset_pu"], "from_vset_pu") : nothing,
          to_vset_pu = haskey(entry, "to_vset_pu") ? _controller_cfg_float(entry["to_vset_pu"], "to_vset_pu") : nothing,
          from_qmin_mvar = haskey(entry, "from_qmin_mvar") ? _controller_cfg_float(entry["from_qmin_mvar"], "from_qmin_mvar") : nothing,
          from_qmax_mvar = haskey(entry, "from_qmax_mvar") ? _controller_cfg_float(entry["from_qmax_mvar"], "from_qmax_mvar") : nothing,
          to_qmin_mvar = haskey(entry, "to_qmin_mvar") ? _controller_cfg_float(entry["to_qmin_mvar"], "to_qmin_mvar") : nothing,
          to_qmax_mvar = haskey(entry, "to_qmax_mvar") ? _controller_cfg_float(entry["to_qmax_mvar"], "to_qmax_mvar") : nothing,
          deadband_vm_pu = haskey(entry, "deadband_vm_pu") ? _controller_cfg_float(entry["deadband_vm_pu"], "deadband_vm_pu") : 1e-3,
          max_outer_iters = haskey(entry, "max_outer_iters") ? _controller_cfg_int(entry["max_outer_iters"], "max_outer_iters") : 20,
          enabled = haskey(entry, "enabled") ? _controller_cfg_bool(entry["enabled"], "enabled") : true,
          name = haskey(entry, "name") ? string(entry["name"]) : nothing,
          from_prosumer = haskey(entry, "from_prosumer") ? _controller_cfg_int(entry["from_prosumer"], "from_prosumer") : nothing,
          to_prosumer = haskey(entry, "to_prosumer") ? _controller_cfg_int(entry["to_prosumer"], "to_prosumer") : nothing,
        )
      end
    catch err
      err isa InterruptException && rethrow()
      throw(ArgumentError("$(label): $(sprint(showerror, err))"))
    end
    added += 1
  end
  return added
end
