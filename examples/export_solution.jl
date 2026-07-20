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

# Date: 2026-05-18
# file: examples/export_solution.jl
# purpose: runs Sparlectra's internal solver, exports a solver-agnostic PFModel/PFSolution, and optionally compares against an external-solver run via runpf_external!

"""
export_solution_for_external_solver.jl

Run Sparlectra's internal solver and export a solver-agnostic PFModel + PFSolution
for external solver experiments. Optionally also run an external-solver path
(via `runpf_external!`) to test the interface and compare against the internal
reference.
"""

using Sparlectra
using Dates
using Printf

include(joinpath(@__DIR__, "internal", "example_header.jl"))

const EXPORT_OUTDIR = joinpath(@__DIR__, "_out", "export_solution")

_angle_delta_deg(calc_deg::Real, ref_deg::Real)::Float64 = mod(Float64(calc_deg) - Float64(ref_deg) + 180.0, 360.0) - 180.0
_case_tag(casefile::AbstractString) = replace(splitext(basename(casefile))[1], r"[^A-Za-z0-9_.-]" => "_")
_csv_quote(s::AbstractString) = "\"" * replace(s, "\"" => "\"\"") * "\""

"Build PF-ordered complex voltage vector from a solved Net using PFModel mapping."
function build_V_pf_from_net(net::Net, model::PFModel)::Vector{ComplexF64}
  V_pf = Vector{ComplexF64}(undef, length(model.busIdx_net))
  @inbounds for (k, busIdx) in enumerate(model.busIdx_net)
    node = net.nodeVec[busIdx]
    V_pf[k] = ComplexF64(node._vm_pu * cosd(node._va_deg), node._vm_pu * sind(node._va_deg))
  end
  return V_pf
end

function _parse_cli_args(args::Vector{String})
  case_override = ""
  run_internal = true
  run_external = true
  show_model = true
  show_solution = true
  show_export_summary = true
  output_dir_override = ""

  for arg in args
    if arg == "--no-internal"
      run_internal = false
    elseif arg == "--no-external"
      run_external = false
    elseif arg == "--hide-model"
      show_model = false
    elseif arg == "--hide-solution"
      show_solution = false
    elseif arg == "--quiet-export"
      show_export_summary = false
    elseif startswith(arg, "--output-dir=")
      output_dir_override = strip(arg[length("--output-dir=")+1:end])
      isempty(output_dir_override) && throw(ArgumentError("`--output-dir=` requires a non-empty path."))
    elseif startswith(arg, "--")
      throw(ArgumentError("Unsupported flag: $(arg)"))
    elseif isempty(case_override)
      case_override = arg
    else
      throw(ArgumentError("Only one positional case override is supported; got extra argument: $(arg)"))
    end
  end

  return (; case_override, run_internal, run_external, show_model, show_solution, show_export_summary, output_dir_override)
end

function _ensure_export_outdir(resolved_case::AbstractString; output_dir_override::AbstractString = "")
  base_dir = isempty(output_dir_override) ? EXPORT_OUTDIR : output_dir_override
  mkpath(base_dir)
  timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
  outdir = joinpath(base_dir, "$(_case_tag(resolved_case))_$(timestamp)")
  mkpath(outdir)
  return outdir
end

function _write_solution_csv(path::AbstractString, net::Net, model::PFModel, sol::PFSolution)
  open(path, "w") do io
    println(io, "pf_index,bus_index,bus_name,vm_pu,va_deg,vr_pu,vi_pu")
    for (k, busIdx) in enumerate(model.busIdx_net)
      node = net.nodeVec[busIdx]
      V = sol.V[k]
      vm = abs(V)
      va = rad2deg(angle(V))
      @printf(io, "%d,%d,%s,%.16g,%.16g,%.16g,%.16g\n", k, busIdx, _csv_quote(node.comp.cName), vm, va, real(V), imag(V))
    end
  end
  return path
end

function _write_comparison_csv(path::AbstractString, net::Net, model_ref::PFModel, sol_ref::PFSolution, sol_ext::PFSolution)
  open(path, "w") do io
    println(io, "pf_index,bus_index,bus_name,vm_internal_pu,vm_external_pu,dvm_pu,va_internal_deg,va_external_deg,dva_deg")
    for (k, busIdx) in enumerate(model_ref.busIdx_net)
      node = net.nodeVec[busIdx]
      Vref = sol_ref.V[k]
      Vext = sol_ext.V[k]
      vm_ref = abs(Vref)
      vm_ext = abs(Vext)
      va_ref = rad2deg(angle(Vref))
      va_ext = rad2deg(angle(Vext))
      dva = _angle_delta_deg(va_ext, va_ref)
      @printf(io, "%d,%d,%s,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g\n", k, busIdx, _csv_quote(node.comp.cName), vm_ref, vm_ext, vm_ext - vm_ref, va_ref, va_ext, dva)
    end
  end
  return path
end

struct DummyExternalSolver <: Sparlectra.AbstractExternalSolver
  casefile::String
  config::Sparlectra.SparlectraConfig
end

function Sparlectra.solvePf(solver::DummyExternalSolver, model::Sparlectra.PFModel; tol = 1e-8, kwargs...)
  local_case = isfile(solver.casefile) ? abspath(solver.casefile) : Sparlectra.FetchMatpowerCase.ensure_casefile(solver.casefile)
  net_tmp = createNetFromMatPowerFile(filename = local_case, log = false, flatstart = solver.config.powerflow.start_mode.flatstart)
  runpf!(net_tmp; config = solver.config)
  model_tmp = buildPfModel(net_tmp; flatstart = false, include_limits = false, verbose = 0)
  V = build_V_pf_from_net(net_tmp, model_tmp)
  r = mismatchInf(model_tmp, V)
  return PFSolution(V = V, converged = (r <= tol), iters = 0, residual_inf = r)
end

function main(args = ARGS)
  print_example_banner("examples/export_solution.jl", "runs Sparlectra's internal solver, exports a solver-agnostic PFModel/PFSolution, and optionally compares against an external-solver run via runpf_external!")
  cli = _parse_cli_args(collect(args))

  config_file = Sparlectra.configuration_path_from_inputs(env_var = "SPARLECTRA_CONFIGURATION_YAML", fallback_paths = [Sparlectra.USER_SPARLECTRA_CONFIG_PATH])

  cfg = isempty(config_file) ? Sparlectra.load_sparlectra_config(; reload = true) : Sparlectra.load_sparlectra_config(config_file; reload = true)

  resolved_case = if !isempty(cli.case_override)
    cli.case_override
  elseif !isempty(cfg.matpower.case)
    cfg.matpower.case
  else
    throw(ArgumentError("No MATPOWER case provided. Pass a casefile as first positional argument or set `matpower.case` in configuration."))
  end

  local_case = Sparlectra.FetchMatpowerCase.ensure_casefile(resolved_case)
  output_dir = _ensure_export_outdir(resolved_case; output_dir_override = cli.output_dir_override)
  summary_file = joinpath(output_dir, "summary.txt")
  matpower_export_file = nothing
  internal_solution_file = nothing
  external_solution_file = nothing
  comparison_file = nothing

  show_results_internal = false
  show_net_verbose = false

  println("Case and configuration:")
  println("  config file            = ", isempty(config_file) ? "default/none" : config_file)
  println("  casefile               = ", resolved_case)
  println("  local case path        = ", local_case)
  println("  powerflow method       = ", cfg.powerflow.method)
  println("  flatstart              = ", cfg.powerflow.start_mode.flatstart)
  println("  tol                    = ", cfg.powerflow.tol)
  println("  run_internal           = ", cli.run_internal)
  println("  run_external           = ", cli.run_external)
  println("  show_model/solution    = ", cli.show_model, " / ", cli.show_solution)
  println("")

  model_ref = nothing
  sol_ref = nothing
  net_ext = nothing
  sol_ext = nothing

  if cli.run_internal
    println("=== Internal solver run (reference) ===")
    solved_result = run_sparlectra(casefile = basename(local_case), path = dirname(local_case), config = cfg)
    solved_net_ref = solved_result.net

    if show_net_verbose
      showNet(solved_net_ref, verbose = true)
    end

    model_ref = buildPfModel(solved_net_ref; flatstart = false, include_limits = true, verbose = (show_net_verbose ? 1 : 0))

    V_ref = build_V_pf_from_net(solved_net_ref, model_ref)
    res_ref = mismatchInf(model_ref, V_ref)

    sol_ref = PFSolution(V = V_ref, converged = isfinite(res_ref) && (res_ref <= cfg.powerflow.tol), iters = -1, residual_inf = res_ref, meta = (solver = :internal, flatstart = cfg.powerflow.start_mode.flatstart))

    if cli.show_export_summary
      println("\nExported INTERNAL reference:")
      println("  n(active)           = ", length(model_ref.busIdx_net))
      println("  slack_pf            = ", model_ref.slack_idx)
      println("  residual ||F(V)||∞  = ", sol_ref.residual_inf)
      println("  converged (by tol)  = ", sol_ref.converged)
    end

    println("")
    internal_solution_file = _write_solution_csv(joinpath(output_dir, "internal_solution.csv"), solved_net_ref, model_ref, sol_ref)
    matpower_export_file = joinpath(output_dir, "$(_case_tag(resolved_case))_export.m")
    writeMatpowerCasefile(solved_net_ref, matpower_export_file)
  end

  if cli.run_external
    println("=== External solver path run (interface test) ===")
    net_ext = createNetFromMatPowerFile(filename = local_case, log = false, flatstart = cfg.powerflow.start_mode.flatstart)

    solver = DummyExternalSolver(local_case, cfg)

    iters_ext, status_ext, sol_ext = runpf_external!(net_ext, solver; tol = cfg.powerflow.tol, flatstart = cfg.powerflow.start_mode.flatstart, include_limits = false, verbose = 0, show_model = cli.show_model, show_solution = cli.show_solution)

    println("External path finished:")
    println("  status              = ", status_ext == 0 ? "converged" : "not converged")
    println("  iters/order         = ", iters_ext)
    println("  residual ||F(V)||∞  = ", sol_ext.residual_inf)
    println("")
    model_ext = buildPfModel(net_ext; flatstart = cfg.powerflow.start_mode.flatstart, include_limits = false, verbose = 0)
    external_solution_file = _write_solution_csv(joinpath(output_dir, "external_solution.csv"), net_ext, model_ext, sol_ext)
  end

  max_dVm = nothing
  max_dVa = nothing
  if cli.run_internal && cli.run_external
    println("=== Compare (internal reference vs external) ===")

    model_cmp = buildPfModel(net_ext; flatstart = cfg.powerflow.start_mode.flatstart, include_limits = false, verbose = 0)

    if length(model_ref.busIdx_net) == length(model_cmp.busIdx_net) && model_ref.busIdx_net == model_cmp.busIdx_net
      max_dVm = maximum(abs.(abs.(sol_ext.V) .- abs.(sol_ref.V)))
      max_dVa = maximum(abs.(_angle_delta_deg.(rad2deg.(angle.(sol_ext.V)), rad2deg.(angle.(sol_ref.V)))))

      println("Max |ΔVm| (p.u.)  = ", max_dVm)
      println("Max |ΔVa| (deg)   = ", max_dVa)
      println("Ref residual      = ", sol_ref.residual_inf)
      println("Ext residual      = ", sol_ext.residual_inf)
      comparison_file = _write_comparison_csv(joinpath(output_dir, "comparison.csv"), net_ext, model_cmp, sol_ref, sol_ext)
    else
      println("Comparison skipped: PF ordering differs between internal/export and external run.")
      println("  internal busIdx_net = ", model_ref.busIdx_net)
      println("  external  busIdx_net = ", model_cmp.busIdx_net)
    end

    println("")
  end

  open(summary_file, "w") do io
    println(io, "Export Solution Summary")
    println(io, "=======================")
    println(io)
    println(io, "Timestamp        : ", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
    println(io, "Config file      : ", isempty(config_file) ? "default/none" : config_file)
    println(io, "Case             : ", resolved_case)
    println(io, "Local case path  : ", local_case)
    println(io, "PF method        : ", cfg.powerflow.method)
    println(io, "Flatstart        : ", cfg.powerflow.start_mode.flatstart)
    println(io, "Tolerance        : ", cfg.powerflow.tol)
    println(io, "Internal enabled : ", cli.run_internal)
    println(io, "External enabled : ", cli.run_external)
    println(io, "MATPOWER export  : ", isnothing(matpower_export_file) ? "n/a (internal disabled)" : matpower_export_file)
    println(io)
    println(io, "Internal reference")
    println(io, "------------------")
    if isnothing(sol_ref)
      println(io, "Converged        : n/a")
      println(io, "Residual inf     : n/a")
    else
      println(io, "Converged        : ", sol_ref.converged)
      println(io, "Residual inf     : ", sol_ref.residual_inf)
    end
    println(io)
    println(io, "External path")
    println(io, "-------------")
    if isnothing(sol_ext)
      println(io, "Converged        : n/a")
      println(io, "Residual inf     : n/a")
    else
      println(io, "Converged        : ", sol_ext.converged)
      println(io, "Residual inf     : ", sol_ext.residual_inf)
    end
    println(io)
    println(io, "Comparison")
    println(io, "----------")
    println(io, "Max |dVm| pu     : ", isnothing(max_dVm) ? "n/a" : string(max_dVm))
    println(io, "Max |dVa| deg    : ", isnothing(max_dVa) ? "n/a" : string(max_dVa))
  end

  println("Export files written to:")
  println("  output directory       : ", output_dir)
  !isnothing(matpower_export_file) && println("  matpower case export   : ", matpower_export_file)
  println("  summary                : ", summary_file)
  !isnothing(internal_solution_file) && println("  internal solution      : ", internal_solution_file)
  !isnothing(external_solution_file) && println("  external solution      : ", external_solution_file)
  !isnothing(comparison_file) && println("  comparison             : ", comparison_file)
  println("Done.")
  return (output_dir = output_dir, summary_file = summary_file, matpower_export_file = matpower_export_file, internal_solution_file = internal_solution_file, external_solution_file = external_solution_file, comparison_file = comparison_file)
end

run_example(main, ARGS)

