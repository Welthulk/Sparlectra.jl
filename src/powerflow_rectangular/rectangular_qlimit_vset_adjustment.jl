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
#
# This file is included inside module Sparlectra. Do not add a module wrapper here.
#
# Rectangular power-flow Q-limit voltage-setpoint adjustment helpers.

# Supports rectangular Q-limit control mode `:adjust_vset`.
function _build_vset_adjust_controllers(net::Net)
  controllers = Dict{Int,NamedTuple{(:prosumer_idx, :config),Tuple{Int,VoltageAdjustConfig}}}()

  for (ps_idx, ps) in enumerate(net.prosumpsVec)
    isGenerator(ps) || continue
    bus = getPosumerBusIndex(ps)

    has_any = _has_vset_adjust_config(ps)
    if has_any
      cfg = _resolve_vset_adjust_config(ps)
      all_defined = !isnothing(cfg)
      all_defined || error("Bus $(_bus_label(net, bus)): invalid voltage adjustment config at prosumer $ps_idx. vstep_pu, tap_steps_down and tap_steps_up must be provided together.")
      cfg.vstep_pu > 0.0 || error("Bus $(_bus_label(net, bus)): invalid vstep_pu=$(cfg.vstep_pu). Must be > 0.")
      cfg.tap_steps_down >= 0 || error("Bus $(_bus_label(net, bus)): invalid tap_steps_down=$(cfg.tap_steps_down). Must be ≥ 0.")
      cfg.tap_steps_up >= 0 || error("Bus $(_bus_label(net, bus)): invalid tap_steps_up=$(cfg.tap_steps_up). Must be ≥ 0.")
      haskey(controllers, bus) && error("Bus $(_bus_label(net, bus)): multiple prosumers define voltage adjustment data. Only one controller per bus is allowed.")

      controllers[bus] = (prosumer_idx = ps_idx, config = cfg)
    else
      partially_set = !isnothing(ps.vstep_pu) || !isnothing(ps.tap_steps_down) || !isnothing(ps.tap_steps_up)
      partially_set && error("Bus $(_bus_label(net, bus)): incomplete voltage adjustment config at prosumer $ps_idx. vstep_pu, tap_steps_down and tap_steps_up must be provided together.")
    end
  end

  return controllers
end
