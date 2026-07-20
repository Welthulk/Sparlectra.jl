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

# Date: 2026-04-24
# file: examples/tap_control_demo_grid.jl
# purpose: lightweight three-controller demo (OLTC + PST + Schraegregler) via run_sparlectra(net=...), reporting controller/trace rows from latest_control_result(net)

using Sparlectra
using Printf

include(joinpath(@__DIR__, "internal", "example_header.jl"))

function demo_yaml_path(args::AbstractVector{String} = ARGS)
  if length(args) >= 1
    return abspath(args[1])
  end
  from_env = strip(get(ENV, "SPARLECTRA_TAP_DEMO_YAML", ""))
  if !isempty(from_env)
    return abspath(from_env)
  end
  local_path = joinpath(@__DIR__, "tap_control_demo_grid.yaml")
  if isfile(local_path)
    return local_path
  end
  return joinpath(@__DIR__, "tap_control_demo_grid.yaml.example")
end

function load_demo_settings(args::AbstractVector{String} = ARGS)
  path = demo_yaml_path(args)
  return path, Sparlectra.load_yaml_dict(path)
end

function load_demo_sparlectra_config(demo_settings::Dict{String,Any}, demo_path::AbstractString)
  raw_path = get(demo_settings, "sparlectra_config", nothing)
  if raw_path isa AbstractString && !isempty(strip(raw_path))
    cfg_path = isabspath(raw_path) ? raw_path : normpath(joinpath(dirname(demo_path), raw_path))
    return cfg_path, Sparlectra.load_sparlectra_config(cfg_path; reload = true)
  end
  return "default discovery", Sparlectra.load_sparlectra_config(; reload = true)
end

function build_demo_net(demo_settings::Dict{String,Any})
  net = Net(name = "tap_control_demo_grid", baseMVA = 100.0)
  for bus in ("Slack", "Grid_A", "Grid_B", "MV_A", "MV_B", "Load_LV", "Load_MV")
    addBus!(net = net, busName = bus, vn_kV = 110.0)
  end

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.01, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "Load_LV", type = "ENERGYCONSUMER", p = -45.0, q = -12.0)
  addProsumer!(net = net, busName = "Load_MV", type = "ENERGYCONSUMER", p = -28.0, q = -8.0)

  addPIModelTrafo!(net = net, fromBus = "Slack", toBus = "Grid_A", r_pu = 0.006, x_pu = 0.05, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelTrafo!(net = net, fromBus = "Grid_A", toBus = "Grid_B", r_pu = 0.004, x_pu = 0.035, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelTrafo!(net = net, fromBus = "MV_A", toBus = "MV_B", r_pu = 0.004, x_pu = 0.03, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)

  addPIModelACLine!(net = net, fromBus = "Grid_B", toBus = "MV_A", r_pu = 0.015, x_pu = 0.09, b_pu = 0.01, status = 1)
  addPIModelACLine!(net = net, fromBus = "Grid_A", toBus = "MV_A", r_pu = 0.02, x_pu = 0.11, b_pu = 0.01, status = 1)
  addPIModelACLine!(net = net, fromBus = "Grid_B", toBus = "MV_B", r_pu = 0.018, x_pu = 0.10, b_pu = 0.01, status = 1)
  addPIModelACLine!(net = net, fromBus = "MV_B", toBus = "Load_LV", r_pu = 0.015, x_pu = 0.07, b_pu = 0.01, status = 1)
  addPIModelACLine!(net = net, fromBus = "MV_A", toBus = "Load_MV", r_pu = 0.02, x_pu = 0.08, b_pu = 0.01, status = 1)

  t_oltc = getNetBranch(net = net, fromBus = "Slack", toBus = "Grid_A")
  t_pst = getNetBranch(net = net, fromBus = "Grid_A", toBus = "Grid_B")
  t_schraeg = getNetBranch(net = net, fromBus = "MV_A", toBus = "MV_B")
  t_oltc.comp.cName = "T_OLTC"
  t_pst.comp.cName = "T_PST"
  t_schraeg.comp.cName = "T_SCHRAEG"

  oltc = get(demo_settings, "oltc", Dict{String,Any}())
  pst = get(demo_settings, "pst", Dict{String,Any}())
  schraeg = get(demo_settings, "schraeg", Dict{String,Any}())

  t_oltc.has_ratio_tap = true
  t_oltc.tap_min = get(oltc, "tap_min", 0.95)
  t_oltc.tap_max = get(oltc, "tap_max", 1.08)
  t_oltc.tap_step = get(oltc, "tap_step", 0.0125)

  t_pst.has_phase_tap = true
  t_pst.phase_min_deg = get(pst, "phase_min_deg", -12.0)
  t_pst.phase_max_deg = get(pst, "phase_max_deg", 12.0)
  t_pst.phase_step_deg = get(pst, "phase_step_deg", 1.0)

  t_schraeg.has_ratio_tap = true
  t_schraeg.tap_min = get(schraeg, "tap_min", 0.95)
  t_schraeg.tap_max = get(schraeg, "tap_max", 1.08)
  t_schraeg.tap_step = get(schraeg, "tap_step", 0.0125)
  t_schraeg.has_phase_tap = true
  t_schraeg.phase_min_deg = get(schraeg, "phase_min_deg", -10.0)
  t_schraeg.phase_max_deg = get(schraeg, "phase_max_deg", 10.0)
  t_schraeg.phase_step_deg = get(schraeg, "phase_step_deg", 1.0)

  return net, (t_oltc = t_oltc, t_pst = t_pst, t_schraeg = t_schraeg)
end

function add_demo_controllers!(net::Net, trafos, demo_settings::Dict{String,Any})
  # Demo controller setpoints come from the example YAML.
  # Controller objects are still created programmatically.
  # This is not the central control.controllers schema.
  oltc = get(demo_settings, "oltc", Dict{String,Any}())
  pst = get(demo_settings, "pst", Dict{String,Any}())
  schraeg = get(demo_settings, "schraeg", Dict{String,Any}())

  addTapController!(
    net;
    trafo = "T_OLTC",
    mode = :voltage,
    target_bus = String(get(oltc, "target_bus", "Load_LV")),
    target_vm_pu = get(oltc, "target_vm_pu", 1.0),
    control_ratio = true,
    control_phase = false,
    is_discrete = true,
    deadband_vm_pu = get(oltc, "deadband_vm_pu", 5e-3),
    max_outer_iters = get(oltc, "max_outer_iters", 12),
  )

  addTapController!(
    net;
    trafo = "T_PST",
    mode = :branch_active_power,
    target_branch = (String(get(pst, "target_branch_from", "Grid_B")), String(get(pst, "target_branch_to", "MV_B"))),
    p_target_mw = get(pst, "p_target_mw", 20.0),
    control_ratio = false,
    control_phase = true,
    is_discrete = true,
    deadband_p_mw = get(pst, "deadband_p_mw", 1.0),
    max_outer_iters = get(pst, "max_outer_iters", 12),
  )

  addTapController!(
    net;
    trafo = "T_SCHRAEG",
    mode = :voltage_and_branch_active_power,
    target_bus = String(get(schraeg, "target_bus", "Load_MV")),
    target_vm_pu = get(schraeg, "target_vm_pu", 1.0),
    target_branch = (String(get(schraeg, "target_branch_from", "MV_A")), String(get(schraeg, "target_branch_to", "MV_B"))),
    p_target_mw = get(schraeg, "p_target_mw", 10.0),
    control_ratio = true,
    control_phase = true,
    is_discrete = true,
    deadband_vm_pu = get(schraeg, "deadband_vm_pu", 5e-3),
    deadband_p_mw = get(schraeg, "deadband_p_mw", 1.0),
    max_outer_iters = get(schraeg, "max_outer_iters", 20),
  )
end

function print_compact_result(result::ControlRunResult)
  println("Control status: ", result.status)
  println("Outer iterations: ", result.outer_iterations)
  println("PF solves: ", result.powerflow_solves)
  println("Controllers:")
  for row in result.controllers
    @printf("  %-18s mode=%-28s status=%-12s tap=%.4f phase=%.2f°\n", row.controller_name, row.mode, row.status, row.tap_ratio, row.phase_shift_deg)
  end
  println("Trace:")
  println("  rows: ", length(result.trace))
  if !isempty(result.trace)
    println("  first: ", first(result.trace))
    println("  last:  ", last(result.trace))
  end
end

function replace_config_field(config::T, name::Symbol, value) where T
  return T(; (field => (field === name ? value : getfield(config, field)) for field in fieldnames(T))...)
end

function tap_demo_run_config(sparlectra_config::SparlectraConfig, show_classic::Bool)
  show_classic && return sparlectra_config
  compact_output = replace_config_field(sparlectra_config.output, :logfile_results, :off)
  return replace_config_field(sparlectra_config, :output, compact_output)
end

function main(args::AbstractVector{String} = ARGS)
  print_example_banner("examples/tap_control_demo_grid.jl", "lightweight three-controller demo (OLTC + PST + Schraegregler) via run_sparlectra(net=...), reporting controller/trace rows from latest_control_result(net)")
  demo_path, demo_settings = load_demo_settings(args)
  config_label, sparlectra_config = load_demo_sparlectra_config(demo_settings, demo_path)
  show_classic = get(ENV, "SPARLECTRA_TAP_DEMO_CLASSIC", "1") != "0"
  run_config = tap_demo_run_config(sparlectra_config, show_classic)
  net, trafos = build_demo_net(demo_settings)

  initial_result = run_sparlectra(net = net, config = run_config)
  initial_result.numerical_converged || error("Uncontrolled PF failed: $(initial_result.reason_text)")
  vm_lv0 = get_bus_vm_pu(net, "Load_LV")
  vm_mv0 = get_bus_vm_pu(net, "Load_MV")

  add_demo_controllers!(net, trafos, demo_settings)
  # Framework output is selected by run_config.output.
  raw = get(ENV, "SPARLECTRA_TAP_DEMO_RAW", "0") == "1"

  controlled_result = run_sparlectra(net = net, config = run_config)
  controlled_result.numerical_converged || error("Controlled PF failed: $(controlled_result.reason_text)")

  vm_lv = get_bus_vm_pu(net, "Load_LV")
  vm_mv = get_bus_vm_pu(net, "Load_MV")
  p_gb_mb = get_branch_p_from_to_mw(net, "Grid_B", "MV_B")
  p_mva_mvb = get_branch_p_from_to_mw(net, "MV_A", "MV_B")

  result = latest_control_result(net)
  result === nothing && error("Expected control result on net.control_result")

  println("Demo YAML: ", demo_path)
  println("Central config: ", config_label)
  println("Classic output: ", show_classic ? "enabled" : "disabled")
  println("Uncontrolled:")
  @printf("  Load_LV Vm: %.6f pu\n", vm_lv0)
  @printf("  Load_MV Vm: %.6f pu\n", vm_mv0)
  println("Controlled:")
  @printf("  Load_LV Vm: %.6f pu\n", vm_lv)
  @printf("  Load_MV Vm: %.6f pu\n", vm_mv)
  @printf("  P(Grid_B -> MV_B): %.4f MW\n", p_gb_mb)
  @printf("  P(MV_A -> MV_B): %.4f MW\n", p_mva_mvb)

  print_compact_result(result)
  println("Tip: set SPARLECTRA_TAP_DEMO_CLASSIC=0 to disable framework result tables for compact-only output.")
  println("Tip: set SPARLECTRA_TAP_DEMO_RAW=1 for raw ControlRunResult rows.")

  if raw
    println("Raw controller rows:")
    for row in result.controllers
      println("  ", row)
    end
    println("Raw trace rows:")
    for row in result.trace
      println("  ", row)
    end
  end
end

run_example(main, ARGS)

