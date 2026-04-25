using Sparlectra
using Printf
using Dates

const DEFAULT_CFG = Dict{String,Any}(
  "max_ite" => 40,
  "pf_tol" => 1e-9,
  "pf_method" => "rectangular",
  "t2_target_bus" => "B9",
  "t2_target_vm_pu" => 1.23,
  "pst1_target_from_bus" => "B1",
  "pst1_target_to_bus" => "B2",
  "pst1_target_p_mw" => 7.0,
  "st1_target_bus" => "B5",
  "st1_target_vm_pu" => 1.07,
  "st1_target_from_bus" => "B4",
  "st1_target_to_bus" => "B5",
  "st1_target_p_mw" => -60.0,
  # Transformer tap parameterization (editable from YAML)
  "t2_tap_min" => 0.90,
  "t2_tap_max" => 1.10,
  "t2_tap_step" => 0.00625,
  "pst1_phase_min_deg" => -15.0,
  "pst1_phase_max_deg" => 15.0,
  "pst1_phase_step_deg" => 1.0,
  "st1_tap_min" => 0.95,
  "st1_tap_max" => 1.05,
  "st1_tap_step" => 0.00625,
  "st1_phase_min_deg" => -10.0,
  "st1_phase_max_deg" => 10.0,
  "st1_phase_step_deg" => 0.5,
  "t2_max_outer_iters" => 40,
  "t2_deadband_vm_pu" => 1e-3,
  "t2_voltage_error_metric" => "vm",
  "pst1_max_outer_iters" => 60,
  "st1_max_outer_iters" => 80,
  "st1_deadband_vm_pu" => 1e-3,
  "st1_voltage_error_metric" => "vm",
)

# -----------------------------------------------------------------------------
# YAML config helpers (simple subset)
# -----------------------------------------------------------------------------
function _parse_yaml_scalar(raw::AbstractString)
  s = strip(raw)
  isempty(s) && return nothing

  if (startswith(s, "\"") && endswith(s, "\"")) || (startswith(s, "'") && endswith(s, "'"))
    return s[2:end-1]
  end

  ls = lowercase(s)
  ls == "true" && return true
  ls == "false" && return false
  ls == "null" && return nothing

  iv = tryparse(Int, s)
  !isnothing(iv) && return iv
  fv = tryparse(Float64, s)
  !isnothing(fv) && return fv
  return s
end

function load_yaml_config(path::AbstractString)
  isempty(path) && return Dict{String,Any}()
  isfile(path) || error("YAML config file not found: $path")

  cfg = Dict{String,Any}()
  for line in eachline(path)
    stripped = strip(line)
    isempty(stripped) && continue
    startswith(stripped, "#") && continue
    occursin(":", stripped) || continue

    key, value_raw = split(stripped, ":"; limit = 2)
    key = strip(key)
    value_raw = strip(split(value_raw, "#"; limit = 2)[1])
    cfg[key] = _parse_yaml_scalar(value_raw)
  end
  return cfg
end

function _yaml_path_from_inputs()
  !isempty(ARGS) && return ARGS[1]
  env_path = get(ENV, "SPARLECTRA_TAP_DEMO_YAML", "")
  !isempty(env_path) && return env_path

  local_default = joinpath(@__DIR__, "tap_control_demo_grid.yaml")
  isfile(local_default) && return local_default

  local_example = joinpath(@__DIR__, "tap_control_demo_grid.yaml.example")
  isfile(local_example) && return local_example
  return ""
end

function _cfg_value(cfg::Dict{String,Any}, key::String)
  return get(cfg, key, DEFAULT_CFG[key])
end

function _target_branch(cfg::Dict{String,Any}, prefix::String)
  return (String(_cfg_value(cfg, "$(prefix)_from_bus")), String(_cfg_value(cfg, "$(prefix)_to_bus")))
end

"""
    build_tap_control_demo_grid(cfg)

Builds a three-voltage-level network with a meshed 220 kV level, a meshed 110 kV
level, and a radial 20 kV level. The grid contains one OLTC transformer, one PST,
and one coupled ratio+phase regulating transformer. Tap and phase limits/steps
are read from `cfg` so they can be edited in YAML.
"""
function build_tap_control_demo_grid(cfg::Dict{String,Any})
  net = Net(name = "tap_control_demo_grid", baseMVA = 100.0)

  # 220 kV mesh
  addBus!(net = net, busName = "B1", vn_kV = 220.0)
  addBus!(net = net, busName = "B2", vn_kV = 220.0)
  addBus!(net = net, busName = "B3", vn_kV = 220.0)

  # 110 kV mesh
  addBus!(net = net, busName = "B4", vn_kV = 110.0)
  addBus!(net = net, busName = "B5", vn_kV = 110.0)
  addBus!(net = net, busName = "B6", vn_kV = 110.0)

  # 20 kV radial feeder
  addBus!(net = net, busName = "B7", vn_kV = 20.0)
  addBus!(net = net, busName = "B8", vn_kV = 20.0)
  addBus!(net = net, busName = "B9", vn_kV = 20.0)

  # Injections and loads
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "B1")
  addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 95.0, q = 12.0)
  addProsumer!(net = net, busName = "B5", type = "ENERGYCONSUMER", p = -90.0, q = -28.0)
  addProsumer!(net = net, busName = "B8", type = "ENERGYCONSUMER", p = -24.0, q = -8.0)
  addProsumer!(net = net, busName = "B9", type = "ENERGYCONSUMER", p = -42.0, q = -14.0)

  # 220 kV mesh branches (PST on B1->B2)
  addPIModelTrafo!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.004, x_pu = 0.05, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B3", r_pu = 0.012, x_pu = 0.10, b_pu = 0.003, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.010, x_pu = 0.09, b_pu = 0.003, status = 1)

  # 220/110 kV interconnection
  addPIModelTrafo!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.006, x_pu = 0.07, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)

  # 110 kV mesh with complex tap transformer on B4->B5
  addPIModelTrafo!(net = net, fromBus = "B4", toBus = "B5", r_pu = 0.005, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B4", toBus = "B6", r_pu = 0.014, x_pu = 0.11, b_pu = 0.004, status = 1)
  addPIModelACLine!(net = net, fromBus = "B6", toBus = "B5", r_pu = 0.013, x_pu = 0.10, b_pu = 0.004, status = 1)

  # 110/20 kV interconnection (OLTC on B6->B7)
  addPIModelTrafo!(net = net, fromBus = "B6", toBus = "B7", r_pu = 0.008, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)

  # 20 kV radial feeder
  addPIModelACLine!(net = net, fromBus = "B7", toBus = "B8", r_pu = 0.018, x_pu = 0.10, b_pu = 0.001, status = 1)
  addPIModelACLine!(net = net, fromBus = "B8", toBus = "B9", r_pu = 0.025, x_pu = 0.12, b_pu = 0.001, status = 1)

  pst1 = getNetBranch(net = net, fromBus = "B1", toBus = "B2")
  st1 = getNetBranch(net = net, fromBus = "B4", toBus = "B5")
  t2 = getNetBranch(net = net, fromBus = "B6", toBus = "B7")

  # T2 OLTC (ratio only)
  t2.has_ratio_tap = true
  t2.has_phase_tap = false
  t2.tap_min = Float64(_cfg_value(cfg, "t2_tap_min"))
  t2.tap_max = Float64(_cfg_value(cfg, "t2_tap_max"))
  t2.tap_step = Float64(_cfg_value(cfg, "t2_tap_step"))
  t2.phase_shift_deg = 0.0
  t2.angle = 0.0

  # PST1 (phase only)
  pst1.has_ratio_tap = false
  pst1.has_phase_tap = true
  pst1.tap_ratio = 1.0
  pst1.ratio = 1.0
  pst1.phase_shift_deg = 0.0
  pst1.angle = 0.0
  pst1.phase_min_deg = Float64(_cfg_value(cfg, "pst1_phase_min_deg"))
  pst1.phase_max_deg = Float64(_cfg_value(cfg, "pst1_phase_max_deg"))
  pst1.phase_step_deg = Float64(_cfg_value(cfg, "pst1_phase_step_deg"))

  # ST1 Schrägregler (ratio + phase)
  st1.has_ratio_tap = true
  st1.has_phase_tap = true
  st1.tap_ratio = 1.0
  st1.ratio = 1.0
  st1.tap_min = Float64(_cfg_value(cfg, "st1_tap_min"))
  st1.tap_max = Float64(_cfg_value(cfg, "st1_tap_max"))
  st1.tap_step = Float64(_cfg_value(cfg, "st1_tap_step"))
  st1.phase_shift_deg = 0.0
  st1.angle = 0.0
  st1.phase_min_deg = Float64(_cfg_value(cfg, "st1_phase_min_deg"))
  st1.phase_max_deg = Float64(_cfg_value(cfg, "st1_phase_max_deg"))
  st1.phase_step_deg = Float64(_cfg_value(cfg, "st1_phase_step_deg"))

  return net
end

function _controller_branches(net::Net)
  return (
    pst1 = getNetBranch(net = net, fromBus = "B1", toBus = "B2"),
    st1 = getNetBranch(net = net, fromBus = "B4", toBus = "B5"),
    t2 = getNetBranch(net = net, fromBus = "B6", toBus = "B7"),
  )
end

"""
    configure_tap_controllers!(net, cfg)

Configures one voltage, one phase-power, and one coupled complex tap controller.
All controller setpoints are loaded from YAML or fallback defaults.
"""
function configure_tap_controllers!(net::Net, cfg::Dict{String,Any})
  if !isdefined(Sparlectra, :addTapController!)
    error("Sparlectra.addTapController! is not available. Please update this example to the current API name.")
  end

  empty!(net.tapControllers)
  c = _controller_branches(net)

  addTapController!(net;
    trafo = string(c.t2.branchIdx),
    mode = :voltage,
    target_bus = String(_cfg_value(cfg, "t2_target_bus")),
    target_vm_pu = Float64(_cfg_value(cfg, "t2_target_vm_pu")),
    control_ratio = true,
    control_phase = false,
    is_discrete = true,
    deadband_vm_pu = Float64(_cfg_value(cfg, "t2_deadband_vm_pu")),
    voltage_error_metric = Symbol(_cfg_value(cfg, "t2_voltage_error_metric")),
    max_outer_iters = Int(_cfg_value(cfg, "t2_max_outer_iters")),
  )

  addTapController!(net;
    trafo = string(c.pst1.branchIdx),
    mode = :branch_active_power,
    target_branch = _target_branch(cfg, "pst1_target"),
    p_target_mw = Float64(_cfg_value(cfg, "pst1_target_p_mw")),
    control_ratio = false,
    control_phase = true,
    is_discrete = true,
    deadband_p_mw = 0.5,
    max_outer_iters = Int(_cfg_value(cfg, "pst1_max_outer_iters")),
  )

  addTapController!(net;
    trafo = string(c.st1.branchIdx),
    mode = :voltage_and_branch_active_power,
    target_bus = String(_cfg_value(cfg, "st1_target_bus")),
    target_vm_pu = Float64(_cfg_value(cfg, "st1_target_vm_pu")),
    target_branch = _target_branch(cfg, "st1_target"),
    p_target_mw = Float64(_cfg_value(cfg, "st1_target_p_mw")),
    control_ratio = true,
    control_phase = true,
    is_discrete = true,
    deadband_vm_pu = Float64(_cfg_value(cfg, "st1_deadband_vm_pu")),
    voltage_error_metric = Symbol(_cfg_value(cfg, "st1_voltage_error_metric")),
    deadband_p_mw = 0.5,
    max_outer_iters = Int(_cfg_value(cfg, "st1_max_outer_iters")),
  )

  return nothing
end

"""
    print_tap_control_targets(net, cfg)

Prints configured target values together with achieved branch/bus values and
current tap states for T2, PST1, and ST1.
"""
function print_tap_control_targets(io::IO, net::Net, cfg::Dict{String,Any})
  c = _controller_branches(net)

  t2_target_bus = String(_cfg_value(cfg, "t2_target_bus"))
  t2_target_vm_pu = Float64(_cfg_value(cfg, "t2_target_vm_pu"))
  pst1_target_branch = _target_branch(cfg, "pst1_target")
  pst1_target_p_mw = Float64(_cfg_value(cfg, "pst1_target_p_mw"))
  st1_target_bus = String(_cfg_value(cfg, "st1_target_bus"))
  st1_target_vm_pu = Float64(_cfg_value(cfg, "st1_target_vm_pu"))
  st1_target_branch = _target_branch(cfg, "st1_target")
  st1_target_p_mw = Float64(_cfg_value(cfg, "st1_target_p_mw"))

  println(io, "Tap Control Targets")
  println(io, "-------------------")
  println(io, "Power sign convention: P is positive in the configured target branch direction (from -> to).")
  println(io, "T2 OLTC:")
  @printf(io, "  target bus       : %s\n", t2_target_bus)
  @printf(io, "  target Vm        : %.4f pu\n", t2_target_vm_pu)
  @printf(io, "  achieved Vm      : %.4f pu\n", get_bus_vm_pu(net, t2_target_bus))
  @printf(io, "  tap ratio        : %.5f\n", c.t2.tap_ratio)
  @printf(io, "  tap range        : %.5f .. %.5f\n", c.t2.tap_min, c.t2.tap_max)
  @printf(io, "  tap step         : %.5f\n", c.t2.tap_step)

  println(io, "\nPST1:")
  @printf(io, "  target branch    : %s -> %s\n", pst1_target_branch[1], pst1_target_branch[2])
  @printf(io, "  target P         : %.3f MW\n", pst1_target_p_mw)
  @printf(io, "  achieved P       : %.3f MW\n", get_branch_p_from_to_mw(net, pst1_target_branch[1], pst1_target_branch[2]))
  @printf(io, "  phase shift      : %.5f deg\n", c.pst1.phase_shift_deg)

  println(io, "\nST1 Schraegregler:")
  @printf(io, "  target bus       : %s\n", st1_target_bus)
  @printf(io, "  target Vm        : %.4f pu\n", st1_target_vm_pu)
  @printf(io, "  achieved Vm      : %.4f pu\n", get_bus_vm_pu(net, st1_target_bus))
  @printf(io, "  target branch    : %s -> %s\n", st1_target_branch[1], st1_target_branch[2])
  @printf(io, "  target P         : %.3f MW\n", st1_target_p_mw)
  @printf(io, "  achieved P       : %.3f MW\n", get_branch_p_from_to_mw(net, st1_target_branch[1], st1_target_branch[2]))
  @printf(io, "  tap ratio        : %.5f\n", c.st1.tap_ratio)
  @printf(io, "  phase shift      : %.5f deg\n", c.st1.phase_shift_deg)

  return nothing
end

function _print_section(io::IO, title::String)
  println(io, "\n======================================================================")
  println(io, title)
  println(io, "======================================================================")
end

function _next_versioned_logfile()
  outdir = joinpath(@__DIR__, "_out")
  mkpath(outdir)
  date_tag = Dates.format(Dates.today(), "yyyy-mm-dd")
  prefix = "run_tap_control_demo_grid_$(date_tag)_v"
  existing = filter(name -> startswith(name, prefix) && endswith(name, ".log"), readdir(outdir))
  version = 1
  if !isempty(existing)
    versions = Int[]
    for name in existing
      m = match(r"_v(\d+)\.log$", name)
      isnothing(m) || push!(versions, parse(Int, m.captures[1]))
    end
    !isempty(versions) && (version = maximum(versions) + 1)
  end
  return joinpath(outdir, "$(prefix)$(lpad(version, 2, '0')).log")
end

function _build_classic_report_string(net::Net, ite::Int, tol::Float64, method::Symbol)
  return mktempdir() do tmpdir
    printACPFlowResults(net, 0.0, ite, tol, true, tmpdir; converged = true, solver = method)
    read(joinpath(tmpdir, "result_$(net.name).txt"), String)
  end
end

function _print_dual(io_log::IO, msg::AbstractString)
  print(stdout, msg)
  print(io_log, msg)
end

function main()
  logfile = _next_versioned_logfile()
  open(logfile, "w") do io_log
    _print_dual(io_log, "Tap Control Demo Grid\n")
    _print_dual(io_log, "Logfile: $logfile\n")

    yaml_path = _yaml_path_from_inputs()
    cfg = merge(copy(DEFAULT_CFG), load_yaml_config(yaml_path))
    if !isempty(yaml_path)
      _print_dual(io_log, "Using YAML config: $yaml_path\n")
    else
      _print_dual(io_log, "Using built-in defaults (no YAML file found)\n")
    end

    net = build_tap_control_demo_grid(cfg)

    max_ite = Int(_cfg_value(cfg, "max_ite"))
    pf_tol = Float64(_cfg_value(cfg, "pf_tol"))
    pf_method = Symbol(_cfg_value(cfg, "pf_method"))

    _print_section(stdout, "Uncontrolled power flow")
    _, erg0, _ = run_net_acpflow(net = net, max_ite = max_ite, tol = pf_tol, verbose = 0, method = pf_method, show_results = false)
    erg0 == 0 || error("Uncontrolled power flow failed with erg=$(erg0)")
    print_tap_control_targets(stdout, net, cfg)

    _print_section(io_log, "Uncontrolled power flow")
    print_tap_control_targets(io_log, net, cfg)

    configure_tap_controllers!(net, cfg)

    _print_section(stdout, "Configured tap controller setpoints")
    print_tap_control_targets(stdout, net, cfg)
    _print_section(io_log, "Configured tap controller setpoints")
    print_tap_control_targets(io_log, net, cfg)

    _print_section(stdout, "Controlled power flow with outer-loop tap control")
    ite, erg1, _ = run_net_acpflow(net = net, max_ite = max_ite, tol = pf_tol, verbose = 0, method = pf_method, show_results = false)
    erg1 == 0 || error("Controlled power flow failed with erg=$(erg1)")
    print_tap_control_targets(stdout, net, cfg)
    printTapControllerSummary(stdout, net)

    _print_section(io_log, "Controlled power flow with outer-loop tap control")
    print_tap_control_targets(io_log, net, cfg)
    printTapControllerSummary(io_log, net)

    classic_report = _build_classic_report_string(net, ite, pf_tol, pf_method)
    _print_section(stdout, "Final classic network report")
    print(stdout, classic_report)
    _print_section(io_log, "Final classic network report")
    print(io_log, classic_report)

    report = buildACPFlowReport(net; ct = 0.0, ite = ite, tol = pf_tol, converged = true, solver = pf_method)
    _print_section(stdout, "Tap controller DataFrame-compatible rows")
    for row in report.transformer_controls
      println(stdout, row)
    end
    _print_section(io_log, "Tap controller DataFrame-compatible rows")
    for row in report.transformer_controls
      println(io_log, row)
    end
    return net
  end
end

Base.invokelatest(main)
