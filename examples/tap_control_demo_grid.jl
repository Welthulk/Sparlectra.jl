using Sparlectra
using Printf

function demo_yaml_path()
  if length(ARGS) >= 1
    return abspath(ARGS[1])
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

function load_demo_settings()
  path = demo_yaml_path()
  return path, Sparlectra.load_yaml_dict(path)
end

function load_demo_sparlectra_config(demo_settings::Dict{String,Any}, demo_path::AbstractString)
  raw_path = get(demo_settings, "sparlectra_config", nothing)
  if raw_path isa AbstractString && !isempty(strip(raw_path))
    cfg_path = isabspath(raw_path) ? raw_path : normpath(joinpath(dirname(demo_path), raw_path))
    return Sparlectra.load_sparlectra_config(cfg_path; reload = true)
  end
  return Sparlectra.load_sparlectra_config(; reload = true)
end

function build_demo_net(demo_settings::Dict{String,Any})
  net = Net(name = "tap_control_demo_grid", baseMVA = 100.0)
  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "Mid", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = -70.0, q = -20.0)

  tap = get(demo_settings, "tap", Dict{String,Any}())
  tap_min = get(tap, "tap_min", 0.95)
  tap_max = get(tap, "tap_max", 1.05)
  tap_step = get(tap, "tap_step", 0.0125)

  addPIModelTrafo!(net = net, fromBus = "Slack", toBus = "Mid", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "Mid", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)

  trafo = getNetBranch(net = net, fromBus = "Slack", toBus = "Mid")
  trafo.has_ratio_tap = true
  trafo.tap_min = tap_min
  trafo.tap_max = tap_max
  trafo.tap_step = tap_step
  return net, trafo
end

function add_demo_controllers!(net::Net, trafo, demo_settings::Dict{String,Any})
  # This demo reads controller setpoints from examples/tap_control_demo_grid.yaml.
  # This is example-specific input. It is not the central control.controllers schema,
  # which is reserved for future YAML-based controller instantiation.
  controller = get(demo_settings, "controller", Dict{String,Any}())
  target_bus = String(get(controller, "target_bus", "Load"))
  addTapController!(
    net;
    trafo = string(trafo.branchIdx),
    mode = :voltage,
    target_bus = target_bus,
    target_vm_pu = get(controller, "target_vm_pu", 0.98),
    control_ratio = true,
    control_phase = false,
    is_discrete = true,
    deadband_vm_pu = get(controller, "deadband_vm_pu", 5e-3),
    max_outer_iters = get(controller, "max_outer_iters", 8),
  )
end

function print_control_result(result::ControlRunResult)
  println("Control status: ", result.status)
  println("Outer iterations: ", result.outer_iterations)
  println("PF solves: ", result.powerflow_solves)
  println("Controller rows:")
  for row in result.controllers
    println("  ", row)
  end
  println("Trace rows:")
  for row in result.trace
    println("  ", row)
  end
end

function main()
  demo_path, demo_settings = load_demo_settings()
  sparlectra_config = load_demo_sparlectra_config(demo_settings, demo_path)

  net, trafo = build_demo_net(demo_settings)

  _, erg0, _ = run_acpflow(net = net; config = sparlectra_config, show_results = false)
  erg0 == 0 || error("Uncontrolled PF failed with erg=$erg0")
  vm_uncontrolled = get_bus_vm_pu(net, "Load")

  add_demo_controllers!(net, trafo, demo_settings)

  show_classic = get(ENV, "SPARLECTRA_TAP_DEMO_CLASSIC", "0") == "1"
  _, erg, _ = run_acpflow(net = net; config = sparlectra_config, show_results = show_classic)
  erg == 0 || error("Controlled PF failed with erg=$erg")
  vm_controlled = get_bus_vm_pu(net, "Load")

  result = latest_control_result(net)
  result === nothing && error("Expected control result on net.control_result")

  println("Demo YAML: ", demo_path)
  @printf("Uncontrolled Load Vm: %.6f pu\n", vm_uncontrolled)
  @printf("Controlled Load Vm:   %.6f pu\n", vm_controlled)
  print_control_result(result)
end

Base.invokelatest(main)
