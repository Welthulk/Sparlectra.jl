# Example: Voltage-dependent Q(U) and P(U) control in the rectangular solver.

using Sparlectra
using Printf

const DEFAULT_EXAMPLE_CFG = Dict{String,Any}("plot_curve" => true, "qu_points" => [(103.4, 35.0), (110.0, 0.0), (116.6, -25.0)], "pu_points" => [(103.4, 25.0), (110.0, 10.0), (116.6, 0.0)])

function _write_default_yaml_example(path::AbstractString)
  mkpath(dirname(path))
  write(
    path,
    """
# Voltage dependent P(U)/Q(U) control demo.
# Usage:
#   julia --project=. src/examples/example_voltage_dependent_control_rectangular.jl
#   julia --project=. src/examples/example_voltage_dependent_control_rectangular.jl path/to/custom.yaml
#   SPARLECTRA_VCTRL_YAML=path/to/custom.yaml julia --project=. src/examples/example_voltage_dependent_control_rectangular.jl
plot_curve: true
qu_points:
  - [103.4, 35.0]
  - [110.0, 0.0]
  - [116.6, -25.0]
pu_points:
  - [103.4, 25.0]
  - [110.0, 10.0]
  - [116.6, 0.0]
""",
  )
end

_parse_pair_line(line::AbstractString) = begin
  m = match(r"^\s*-\s*\[\s*([^\],]+)\s*,\s*([^\],]+)\s*\]\s*$", line)
  isnothing(m) && return nothing
  return (parse(Float64, m.captures[1]), parse(Float64, m.captures[2]))
end

function load_voltage_ctrl_yaml(path::AbstractString)
  cfg = Dict{String,Any}()
  pairs = Tuple{Float64,Float64}[]
  list_key = ""
  for raw in eachline(path)
    line = strip(split(raw, "#"; limit = 2)[1])
    isempty(line) && continue
    if line == "qu_points:" || line == "pu_points:"
      if !isempty(list_key)
        cfg[list_key] = copy(pairs)
        empty!(pairs)
      end
      list_key = replace(line, ":" => "")
      continue
    end
    if startswith(line, "-")
      p = _parse_pair_line(line)
      isnothing(p) && error("Invalid YAML list entry in $path: $line")
      push!(pairs, p)
      continue
    end
    if !isempty(list_key)
      cfg[list_key] = copy(pairs)
      empty!(pairs)
      list_key = ""
    end
    kv = split(line, ":"; limit = 2)
    length(kv) == 2 || error("Invalid YAML line in $path: $line")
    key = strip(kv[1])
    val = lowercase(strip(kv[2]))
    if key == "plot_curve"
      cfg[key] = (val == "true" || val == "1" || val == "yes" || val == "on")
    end
  end
  if !isempty(list_key)
    cfg[list_key] = copy(pairs)
  end
  return cfg
end

function _resolve_yaml_path()
  length(ARGS) >= 1 && isfile(ARGS[1]) && return ARGS[1]
  env_path = get(ENV, "SPARLECTRA_VCTRL_YAML", "")
  !isempty(env_path) && isfile(env_path) && return env_path

  local_cfg = joinpath(@__DIR__, "example_voltage_dependent_control_rectangular.yaml")
  local_example = joinpath(@__DIR__, "example_voltage_dependent_control_rectangular.yaml.example")
  if !isfile(local_example)
    _write_default_yaml_example(local_example)
  end
  return isfile(local_cfg) ? local_cfg : local_example
end

function _maybe_plot_characteristic(ch::PiecewiseLinearCharacteristic; enable_plot::Bool = false, out_file::String = joinpath(pkgdir(Sparlectra), "src", "examples", "_out", "voltage_control_characteristic.png"))
  enable_plot || return

  if isnothing(Base.find_package("Plots"))
    @info "Skipping curve plot: package Plots.jl is not available in the current/global Julia environment."
    @info "Install globally if desired: julia -e 'using Pkg; Pkg.add(\"Plots\"); Pkg.add(\"GR\")'"
    return
  end

  Base.invokelatest(() -> (@eval import Plots))
  Base.invokelatest(() -> (@eval import GR))

  plots_mod = getfield(Main, :Plots)
  xs = range(ch.points[1][1], ch.points[end][1]; length = 200)
  ys = [first(evaluate_characteristic(ch, x)) for x in xs]
  xk = [p[1] for p in ch.points]
  yk = [p[2] for p in ch.points]

  mkpath(dirname(out_file))
  plt = Base.invokelatest(getfield(plots_mod, :plot), xs, ys; lw = 2, xlabel = "Voltage [pu]", ylabel = "P [pu]", label = "Characteristic", title = "Voltage-dependent characteristic")
  Base.invokelatest(getfield(plots_mod, :scatter!), plt, xk, yk; ms = 4, label = "Support points")
  Base.invokelatest(getfield(plots_mod, :savefig), plt, out_file)
  @info "Characteristic curve plot written to $(abspath(out_file))"
end

function _run_kink_check()
  kink_curve = make_characteristic([(104.5, 20.0), (107.5, 14.0), (110.0, 10.0), (113.0, 14.0), (115.5, 20.0)]; voltage_unit = :kV, value_unit = :MW, vn_kV = 110.0, sbase_MVA = 100.0, interpolation = :polynomial)

  _, slope_left = evaluate_characteristic(kink_curve, 1.0 - 1e-8)
  value_kink, slope_kink = evaluate_characteristic(kink_curve, 1.0)
  _, slope_right = evaluate_characteristic(kink_curve, 1.0 + 1e-8)

  @assert isapprox(value_kink, 0.1; atol = 1e-12)
  @assert isfinite(slope_left)
  @assert isfinite(slope_kink)
  @assert isfinite(slope_right)
  @assert isapprox(slope_kink, slope_left; atol = 5e-6)
  @assert isapprox(slope_kink, slope_right; atol = 5e-6)

  return (curve = kink_curve, slope_left = slope_left, slope_kink = slope_kink, slope_right = slope_right)
end

function run_example_voltage_dependent_control_rectangular(; verbose::Int = 1, plot_curve::Union{Bool,Nothing} = nothing)
  yaml_path = _resolve_yaml_path()
  yaml_cfg = load_voltage_ctrl_yaml(yaml_path)
  cfg = merge(copy(DEFAULT_EXAMPLE_CFG), yaml_cfg)
  used_plot_curve = isnothing(plot_curve) ? Bool(cfg["plot_curve"]) : Bool(plot_curve)

  kink_check = Base.invokelatest(_run_kink_check)
  Base.invokelatest(_maybe_plot_characteristic, kink_check.curve; enable_plot = used_plot_curve)

  net = Net(name = "voltage_dependent_control", baseMVA = 100.0)

  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "Prosumer", vn_kV = 110.0)

  addACLine!(net = net, fromBus = "Slack", toBus = "Prosumer", length = 30.0, r = 0.05, x = 0.45)

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")

  qu_curve = make_characteristic(cfg["qu_points"]; voltage_unit = :kV, value_unit = :MVAr, vn_kV = 110.0, sbase_MVA = net.baseMVA)
  pu_curve = make_characteristic(cfg["pu_points"]; voltage_unit = :kV, value_unit = :MW, vn_kV = 110.0, sbase_MVA = net.baseMVA)

  addProsumer!(net = net, busName = "Prosumer", type = "SYNCHRONOUSMACHINE", p = 10.0, q = 0.0, qu_controller = QUController(qu_curve; qmin_MVAr = -50.0, qmax_MVAr = 50.0, sbase_MVA = net.baseMVA), pu_controller = PUController(pu_curve; pmin_MW = 0.0, pmax_MW = 50.0, sbase_MVA = net.baseMVA))

  addProsumer!(net = net, busName = "Prosumer", type = "ENERGYCONSUMER", p = 45.0, q = 18.0)

  ite, erg, etime = run_net_acpflow(net = net, max_ite = 40, tol = 1e-9, verbose = verbose, opt_sparse = true, method = :rectangular, show_results = false)
  converged = (erg == 0)

  if converged
    V = buildVoltageVector(net)
    Sspec, _, _ = buildControlledSVec(net, V)
    calcNetLosses!(net)

    vm = abs(V[2])
    p_ctrl_pu, _ = evaluate_controller(net.prosumpsVec[2].puController, vm)
    q_ctrl_pu, _ = evaluate_controller(net.prosumpsVec[2].quController, vm)

    println("\nConverged in $ite iterations")
    @printf("Bus %-10s |V| = %.5f pu, angle = %.3f deg\n", "Slack", abs(V[1]), rad2deg(angle(V[1])))
    @printf("Bus %-10s |V| = %.5f pu, angle = %.3f deg\n", "Prosumer", abs(V[2]), rad2deg(angle(V[2])))
    @printf("Controlled prosumer setpoints at Vm=%.5f pu: P=%.3f MW, Q=%.3f MVar\n", vm, p_ctrl_pu * net.baseMVA, q_ctrl_pu * net.baseMVA)
    @printf("Net specified injection at Prosumer bus: P=%.3f MW, Q=%.3f MVar\n", real(Sspec[2]) * net.baseMVA, imag(Sspec[2]) * net.baseMVA)
    @printf("Characteristic slope check (dP/dU): left=%.6f, center=%.6f, right=%.6f\n", kink_check.slope_left, kink_check.slope_kink, kink_check.slope_right)
    println("\nDetailed results:")
    printACPFlowResults(net, etime, ite, 1e-9; converged = true, solver = :rectangular)
  else
    println("Power flow did not converge.")
  end

  return net
end

# IMPORTANT for Julia 1.12 / Revise world-age safety:
# Default: run immediately (also when included from REPL), like a script.
if get(ENV, "SPARLECTRA_SUITE_NO_AUTORUN", "0") != "1"
  main_fn = getfield(@__MODULE__, :main)
  Base.invokelatest(main_fn)
end

Base.invokelatest(() -> run_example_voltage_dependent_control_rectangular())
