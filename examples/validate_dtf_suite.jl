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

# Unified DTF validation suite.
#
# Place this file in examples/ next to dtf_for002_validation_utils.jl.
# It combines:
#   - DTF import audit
#   - native DTF/FOR002 base-case validation
#   - native DTF/FOR002 outage validation
#   - DTF -> MATPOWER -> Sparlectra roundtrip validation
#
# Default input:  <repository>/data/DTF/FOR001*.DAT and FOR002*.DAT
# Default output: examples/_out/dtf_validation_suite/

using Dates
using Printf
using Sparlectra


module NativeBaseValidation
using Sparlectra
using Printf

include(joinpath(@__DIR__, "dtf_for002_validation_utils.jl"))

# Developer notes:
# - Validates the native Testnetz13 DTF base case against FOR002.
# - The intentional execution path is DTFImporter.read_dtf -> DTFImporter.build_net -> runpf!.
# - MATPOWER import/export and the generated FOR001 builder are intentionally not used.
# - FOR002 is treated as a legacy textual reference report.
# - State residuals are diagnostic, not pass/fail thresholds.
# - Console output is intentionally concise; detailed rows are written to
#   CSV/Markdown and returned only in detailed mode.

struct DTFFor002ValidationResult
  output_dir::String
  converged::Bool
  iterations::Int
  final_mismatch::Float64
  metrics::NamedTuple
  written_files::Vector{String}
  max_branch_d_p_MW::Float64
  max_branch_d_q_MVar::Float64
  max_qgen_d_MVar::Float64
  max_state_residual_p_MW::Float64
  max_state_residual_q_MVar::Float64
end

function _assert_finite_dtf_model(case)
  problems = String[]
  for br in case.branches
    all(isfinite, (br.r_ohm, br.x_ohm, br.g_s, br.b_s)) || push!(problems, "branch $(br.index) $(br.from)->$(br.to) has non-finite r/x/g/b")
  end
  for bus in case.buses
    vn = case.nominal_voltages_kv[bus.voltage_level_index]
    vm = bus.start_kv > 0 ? bus.start_kv / vn : 1.0
    isfinite(vn) && vn > 0 || push!(problems, "bus $(bus.name) has invalid nominal voltage $(vn)")
    isfinite(vm) || push!(problems, "bus $(bus.name) has non-finite start vm_pu")
  end
  isempty(problems) || error("DTF parser/model audit failed before solve:\n" * join(problems, "\n"))
  return nothing
end

function _parse_bool(s::AbstractString)
  v = lowercase(strip(s));
  v in ("true", "1", "yes") && return true;
  v in ("false", "0", "no") && return false
  throw(ArgumentError("invalid boolean: $s"))
end

function parse_cli(args)
  # The validation data are stored in <repository>/data/DTF.  Resolve the
  # paths relative to this script so that the program works independently of
  # the current working directory.
  dtf_data_dir = normpath(joinpath(@__DIR__, "..", "data", "DTF"))
  dtf_default = joinpath(dtf_data_dir, "FOR001.DAT")
  for002_default = joinpath(dtf_data_dir, "FOR002.DAT")

  opt = Dict{String,Any}(
    "dtf-file" => dtf_default,
    "for002-file" => for002_default,
    "output-dir" => joinpath(@__DIR__, "_out", "dtf_for002_native_validation"),
    "tol" => 1e-8,
    "max-iter" => 50,
    "method" => "rectangular",
    "write-csv" => true,
    "write-markdown" => true,
    "quiet" => false,
    "strict" => false,
    "details" => false,
    "print-summary" => true,
    "legacy-voltage-level-collapse-230kv" => false,
    "transformer-ratio-mode" => "neutral_one",
  )
  for arg in args
    arg == "--quiet" && (opt["quiet"] = true; continue)
    startswith(arg, "--") || throw(ArgumentError("unsupported argument: $arg"))
    k, v = split(arg[3:end], "="; limit = 2)
    if k in ("tol",)
      opt[k] = parse(Float64, v)
    elseif k in ("max-iter",)
      opt[k] = parse(Int, v)
    elseif k in ("write-csv", "write-markdown", "strict", "quiet", "details", "print-summary", "legacy-voltage-level-collapse-230kv")
      opt[k] = _parse_bool(v)
    else
      opt[k] = v
    end
  end
  return opt
end

function _method_symbol(s::AbstractString)
  s == "default" && return nothing
  s in ("rectangular", "polar") || throw(ArgumentError("--method must be rectangular, polar, or default"))
  return Symbol(s)
end

function _transformer_ratio_mode_symbol(s::AbstractString)
  normalized = replace(lowercase(strip(s)), "-" => "_")
  normalized in ("neutral_one", "winding_over_network") || throw(ArgumentError("--transformer-ratio-mode must be neutral_one or winding_over_network"))
  return Symbol(normalized)
end

function _slack_bus_name(net, case)
  for (idx, n) in enumerate(net.nodeVec)
    if uppercase(Sparlectra.toString(n._nodeType)) == "SLACK"
      return idx <= length(case.buses) ? case.buses[idx].name : n.comp.cName
    end
  end
  return ""
end

function _maxabs(rows, field)
  vals = [abs(Float64(getproperty(r, field))) for r in rows if !(getproperty(r, field) isa Missing) && isfinite(Float64(getproperty(r, field)))]
  return isempty(vals) ? NaN : maximum(vals)
end
_meanabs(rows, field) = (vals = [abs(Float64(getproperty(r, field))) for r in rows if !(getproperty(r, field) isa Missing) && isfinite(Float64(getproperty(r, field)))]; isempty(vals) ? NaN : sum(vals) / length(vals))

function _top(rows, field; n = 10)
  return first(sort(rows, by = r -> -abs(Float64(getproperty(r, field)))), min(n, length(rows)))
end

function _print_summary(result::DTFFor002ValidationResult)
  m = result.metrics
  println("Native DTF/FOR002 validation")
  println("Output directory: ", result.output_dir)
  println("Converged: ", result.converged)
  println("Iterations: ", result.iterations)
  println("Final mismatch: ", result.final_mismatch)
  println("Max |dV|: ", m.max_abs_d_vm_kV, " kV / ", m.max_abs_d_vm_pu, " pu")
  println("Max |dVa|: ", m.max_abs_d_va_deg, " deg")
  println("Max branch |dP|: ", result.max_branch_d_p_MW, " MW")
  println("Max branch |dQ|: ", result.max_branch_d_q_MVar, " MVar")
  println("Max |dQgen| solved-vs-FOR002: ", result.max_qgen_d_MVar, " MVar")
  println("Max state residual |P|: ", result.max_state_residual_p_MW, " MW")
  println("Max state residual |Q|: ", result.max_state_residual_q_MVar, " MVar")
  println("Written files:")
  for path in result.written_files
    println("  - ", basename(path))
  end
  return nothing
end

function run_validation(args = ARGS; return_details::Bool = false)
  opt = parse_cli(args)
  return_details = return_details || opt["details"]
  mkpath(opt["output-dir"])
  case = Sparlectra.DTFImporter.read_dtf(opt["dtf-file"]; strict = opt["strict"])
  _assert_finite_dtf_model(case)
  transformer_ratio_mode = _transformer_ratio_mode_symbol(opt["transformer-ratio-mode"])
  net = Sparlectra.DTFImporter.build_net(case; legacy_voltage_level_collapse_230kv = opt["legacy-voltage-level-collapse-230kv"], transformer_ratio_mode = transformer_ratio_mode)
  ref = parse_for002_ground_load_flow(opt["for002-file"])
  method = _method_symbol(opt["method"])
  iters, status = method === nothing ? Sparlectra.runpf!(net, opt["max-iter"], opt["tol"], 0) : Sparlectra.runpf!(net, opt["max-iter"], opt["tol"], 0; method = method)
  Sparlectra.calcNetLosses!(net)
  final_mismatch = try
    model = Sparlectra.buildPfModel(net; flatstart = false, include_limits = false, verbose = 0)
    V = [ComplexF64(n._vm_pu * cosd(n._va_deg), n._vm_pu * sind(n._va_deg)) for n in net.nodeVec]
    Sparlectra.mismatchInf(model, V)
  catch
    NaN
  end
  converged = status == 0

  bus_by_name = Dict(_norm_name(b.name) => b for b in case.buses)
  import_rows = [(
    baseMVA = case.baseMVA,
    bus_count = length(case.buses),
    branch_count = length(case.branches),
    line_count = count(b -> b.kind != 'T', case.branches),
    transformer_count = count(b -> b.kind == 'T', case.branches),
    transformer_control_count = length(case.transformer_controls),
    outage_count = length(case.outages),
    nominal_voltages_kv = join(case.nominal_voltages_kv, ";"),
    dtf_slack_bus = case.size.slack,
    net_slack_bus = _slack_bus_name(net, case),
    notes = "native DTFImporter path; ground-load-flow only",
    transformer_ratio_mode = transformer_ratio_mode,
  )]

  branch_kcl_p, branch_kcl_q = _branch_kcl_arrays(net)
  bus_rows = NamedTuple[];
  gen_rows = NamedTuple[]
  for db in case.buses
    node = net.nodeVec[db.index]
    key = _norm_name(db.name);
    haskey(ref.buses, key) || continue
    fb = ref.buses[key];
    flags = _bus_type_flags(db.bus_type);
    gl = _bus_generation_load(net, node.busIdx)
    push!(
      bus_rows,
      (
        bus_name = db.name,
        dtf_bus_type = db.bus_type,
        vn_kV = node.comp.cVN,
        for002_vm_kV = fb.v_kV,
        model_vm_kV = node.comp.cVN * node._vm_pu,
        d_vm_kV = node.comp.cVN * node._vm_pu - fb.v_kV,
        for002_vm_pu = fb.v_kV / node.comp.cVN,
        model_vm_pu = node._vm_pu,
        d_vm_pu = node._vm_pu - fb.v_kV / node.comp.cVN,
        for002_va_deg = fb.va_deg,
        model_va_deg = node._va_deg,
        d_va_deg = node._va_deg - fb.va_deg,
        is_slack = flags.is_slack,
        is_pv = flags.is_pv,
        is_pq = flags.is_pq,
      ),
    )
    # Slack/PV generator P/Q are solved quantities, not the specified input
    # values. For Slack Q, infer the solved value from branch KCL plus load.
    # PQ generator Q is expected to stay at its fixed specified value.
    solved_pg = (flags.is_slack || flags.is_pv) ? branch_kcl_p[node.busIdx] + gl.pl : gl.pg_res
    solved_qg = (flags.is_slack || flags.is_pv) ? branch_kcl_q[node.busIdx] + gl.ql : gl.qg_res
    push!(
      gen_rows,
      (
        bus_name = db.name,
        dtf_bus_type = db.bus_type,
        model_bus_type_effective = Sparlectra.toString(node._nodeType),
        for002_pg_MW = fb.p_gen_MW,
        model_pg_specified_MW = gl.pg_spec,
        model_pg_result_MW = solved_pg,
        d_pg_result_vs_for002_MW = solved_pg - fb.p_gen_MW,
        for002_qg_MVar = fb.q_gen_MVar,
        model_qg_specified_MVar = gl.qg_spec,
        model_qg_result_MVar = solved_qg,
        d_qg_result_vs_for002_MVar = solved_qg - fb.q_gen_MVar,
        qmin_MVar = something(db.qmin_mvar, missing),
        qmax_MVar = something(db.qmax_mvar, missing),
        is_regulating = gl.is_regulating || flags.is_pv,
        is_slack = flags.is_slack,
        has_generator = gl.has_generator,
        has_load = gl.has_load,
        p_load_MW = gl.pl,
        q_load_MVar = gl.ql,
        for002_p_load_MW = fb.p_load_MW,
        for002_q_load_MVar = fb.q_load_MVar,
        for002_p_net_MW = fb.p_gen_MW - fb.p_load_MW,
        for002_q_net_MVar = fb.q_gen_MVar - fb.q_load_MVar,
        model_p_net_specified_MW = gl.pg_spec - gl.pl,
        model_q_net_specified_MVar = gl.qg_spec - gl.ql,
        model_p_net_result_MW = solved_pg - gl.pl,
        model_q_net_result_MVar = solved_qg - gl.ql,
        notes = flags.is_slack ? "slack bus: solved P/Q inferred from branch KCL plus load; specified Q is input/schedule" :
                (gl.is_regulating || flags.is_pv ? "PV/regulating bus: result P/Q inferred from branch KCL plus load; specified Q is input/schedule" : (gl.has_generator ? "PQ generator bus: specified Q is expected to remain fixed" : "no model generator")),
      ),
    )
  end

  rd = branch_ref_dict(ref);
  branch_rows = NamedTuple[]
  for (idx, br) in enumerate(net.branchVec)
    meta = net.matpower_branch_metadata[idx]
    nr = uppercase(strip(String(meta.parallel_id)))
    dtfb = case.branches[idx];
    f = dtfb.from;
    t = dtfb.to
    rf = get(rd, (_norm_name(f), _norm_name(t), nr), missing);
    rt = get(rd, (_norm_name(t), _norm_name(f), nr), missing)
    pf = isnothing(br.fBranchFlow) || isnothing(br.fBranchFlow.pFlow) ? missing : br.fBranchFlow.pFlow;
    qf = isnothing(br.fBranchFlow) || isnothing(br.fBranchFlow.qFlow) ? missing : br.fBranchFlow.qFlow
    pt = isnothing(br.tBranchFlow) || isnothing(br.tBranchFlow.pFlow) ? missing : br.tBranchFlow.pFlow;
    qt = isnothing(br.tBranchFlow) || isnothing(br.tBranchFlow.qFlow) ? missing : br.tBranchFlow.qFlow
    push!(
      branch_rows,
      (
        branch_index = idx,
        branch_label = meta.orig_name,
        branch_kind = string(meta.dtf_kind),
        voltage_level_index = meta.voltage_level_index,
        parallel_id = meta.parallel_id,
        from_bus = f,
        to_bus = t,
        u_ref_kV = meta.u_ref_kV,
        r_pu = br.r_pu,
        x_pu = br.x_pu,
        b_pu = br.b_pu,
        ratio = br.ratio,
        transformer_ratio_mode = hasproperty(meta, :transformer_ratio_mode) ? meta.transformer_ratio_mode : missing,
        base_ratio_used = hasproperty(meta, :base_ratio_used) ? meta.base_ratio_used : missing,
        winding_over_network_base_ratio = hasproperty(meta, :winding_over_network_base_ratio) ? meta.winding_over_network_base_ratio : missing,
        nominal_unregulated_kv = hasproperty(meta, :nominal_unregulated_kv) ? something(meta.nominal_unregulated_kv, missing) : missing,
        nominal_regulated_kv = hasproperty(meta, :nominal_regulated_kv) ? something(meta.nominal_regulated_kv, missing) : missing,
        tap_fraction = hasproperty(meta, :tap_fraction) ? meta.tap_fraction : missing,
        skew_angle_deg = hasproperty(meta, :skew_angle_deg) ? meta.skew_angle_deg : missing,
        effective_ratio = hasproperty(meta, :effective_ratio) ? meta.effective_ratio : missing,
        phase_shift_deg = hasproperty(meta, :phase_shift_deg) ? meta.phase_shift_deg : missing,
        for002_p_from_MW = rf === missing ? missing : rf.p_MW,
        model_p_from_MW = pf,
        d_p_from_MW = rf === missing || pf === missing ? missing : pf - rf.p_MW,
        for002_q_from_MVar = rf === missing ? missing : rf.q_MVar,
        model_q_from_MVar = qf,
        d_q_from_MVar = rf === missing || qf === missing ? missing : qf - rf.q_MVar,
        for002_p_to_MW = rt === missing ? missing : rt.p_MW,
        model_p_to_MW = pt,
        d_p_to_MW = rt === missing || pt === missing ? missing : pt - rt.p_MW,
        for002_q_to_MVar = rt === missing ? missing : rt.q_MVar,
        model_q_to_MVar = qt,
        d_q_to_MVar = rt === missing || qt === missing ? missing : qt - rt.q_MVar,
      ),
    )
  end
  residual_rows = for002_state_residual_rows(net, case, ref)
  kcl_rows = branch_kcl_rows(net, case, ref)
  q_diag_rows = q_semantics_diagnostic_rows(gen_rows, kcl_rows, residual_rows)
  p_loss, q_loss = Sparlectra.getTotalLosses(net = net)
  metrics_rows = [(
    converged = converged,
    iterations = iters,
    final_mismatch = final_mismatch,
    max_abs_d_vm_kV = _maxabs(bus_rows, :d_vm_kV),
    max_abs_d_vm_pu = _maxabs(bus_rows, :d_vm_pu),
    max_abs_d_va_deg = _maxabs(bus_rows, :d_va_deg),
    max_abs_d_pg_MW = _maxabs(gen_rows, :d_pg_result_vs_for002_MW),
    max_abs_d_qg_MVar = _maxabs(gen_rows, :d_qg_result_vs_for002_MVar),
    max_abs_branch_d_p_MW = max(_maxabs(branch_rows, :d_p_from_MW), _maxabs(branch_rows, :d_p_to_MW)),
    max_abs_branch_d_q_MVar = max(_maxabs(branch_rows, :d_q_from_MVar), _maxabs(branch_rows, :d_q_to_MVar)),
    max_abs_state_residual_p_MW = _maxabs(residual_rows, :d_p_MW),
    max_abs_state_residual_q_MVar = _maxabs(residual_rows, :d_q_MVar),
    max_abs_branch_kcl_vs_for002_q_MVar = _maxabs(kcl_rows, :d_branch_kcl_vs_for002_q_MVar),
    mean_abs_state_residual_p_MW = _meanabs(residual_rows, :d_p_MW),
    mean_abs_state_residual_q_MVar = _meanabs(residual_rows, :d_q_MVar),
  )]

  written_files = String[]
  if opt["write-csv"]
    for (filename, rows) in [
      ("dtf_import_summary.csv", import_rows),
      ("dtf_bus_comparison.csv", bus_rows),
      ("dtf_generator_comparison.csv", gen_rows),
      ("dtf_bus_kcl_comparison.csv", kcl_rows),
      ("dtf_q_semantics_diagnostics.csv", q_diag_rows),
      ("dtf_branch_comparison.csv", branch_rows),
      ("dtf_state_residual.csv", residual_rows),
      ("dtf_validation_metrics.csv", metrics_rows),
    ]
      path = joinpath(opt["output-dir"], filename)
      write_csv(path, rows)
      push!(written_files, path)
    end
  end
  if opt["write-markdown"]
    summary_path = joinpath(opt["output-dir"], "dtf_for002_validation_summary.md")
    open(summary_path, "w") do io
      println(io, "# Native DTF/FOR002 validation summary\n")
      println(io, "This validation uses `DTFImporter.read_dtf` -> `DTFImporter.build_net` and does not use MATPOWER import, MATPOWER export, or the generated FOR001 builder. Outages are parsed but not executed.\n")
      println(io, "- DTF file: `", opt["dtf-file"], "`")
      println(io, "- FOR002 file: `", opt["for002-file"], "`")
      println(io, "- transformer ratio mode: `", transformer_ratio_mode, "`")
      println(io, "- baseMVA: ", case.baseMVA, "; buses: ", length(case.buses), "; branches: ", length(case.branches), "; lines: ", count(b -> b.kind != 'T', case.branches), "; transformers: ", count(b -> b.kind == 'T', case.branches))
      println(io, "- transformer controls: ", length(case.transformer_controls), "; outages parsed: ", length(case.outages), "; nominal voltages kV: ", join(case.nominal_voltages_kv, ", "))
      println(io, "- DTF slack bus: ", case.size.slack, "; Sparlectra slack bus: ", _slack_bus_name(net, case))
      println(io, "- solver method: ", opt["method"], "; converged: ", converged, "; iterations: ", iters, "; final mismatch: ", final_mismatch)
      println(io, "- total generation MW/MVar: ", sum(r.model_pg_result_MW for r in gen_rows), " / ", sum(r.model_qg_result_MVar for r in gen_rows))
      println(io, "- total load MW/MVar: ", sum(fb.p_load_MW for fb in values(ref.buses)), " / ", sum(fb.q_load_MVar for fb in values(ref.buses)), "; losses MW/MVar: ", p_loss, " / ", q_loss, "\n")
      println(io, "## What are state residuals?\n")
      println(
        io,
        "State residuals force FOR002 printed voltage magnitudes/angles into the native Sparlectra Ybus. The resulting bus injections are compared with the FOR002 printed bus table. This is more sensitive than solved branch-flow comparisons because FOR002 values may be rounded and transformer-adjacent nodes react strongly to small voltage/angle differences. These residuals are diagnostic, not hard pass/fail criteria yet; branch-flow deviations and solved generator/slack comparisons are currently stronger validation signals.\n",
      )
      for (title, rows, field, label) in [
        ("Top 10 bus voltage deviations", bus_rows, :d_vm_kV, "kV"),
        ("Top 10 branch P deviations", branch_rows, :d_p_from_MW, "MW"),
        ("Top 10 branch Q deviations", branch_rows, :d_q_from_MVar, "MVar"),
        ("Top 10 state residual P deviations", residual_rows, :d_p_MW, "MW"),
        ("Top 10 state residual Q deviations", residual_rows, :d_q_MVar, "MVar"),
      ]
        println(io, "## ", title, "\n")
        for r in _top([x for x in rows if !(getproperty(x, field) isa Missing)], field)
          name = hasproperty(r, :bus_name) ? r.bus_name : hasproperty(r, :branch_label) ? r.branch_label : string(r)
          println(io, "- ", name, ": ", getproperty(r, field), " ", label)
        end
        println(io)
      end
      println(io, "## Before/after residual interpretation\n")
      println(io, "Old converted-model FOR002-state residuals included examples such as ALPHA S1 dP around -845 MW. Native DTF residuals for selected buses are:")
      for target in ["ALPHA S1", "BETA1 S1", "DELTA1S1", "BETA2 S1", "DELTA2S1", "WEILERS1"]
        row = findfirst(r -> _norm_name(r.bus_name) == _norm_name(target), residual_rows)
        row === nothing ? println(io, "- ", target, ": not parsed") : println(io, "- ", residual_rows[row].bus_name, ": dP=", residual_rows[row].d_p_MW, " MW, dQ=", residual_rows[row].d_q_MVar, " MVar")
      end
      println(io, "\nThe selected native residuals are materially smaller than the old ALPHA S1 -845 MW example when their absolute dP values are far below that magnitude; no hard pass/fail threshold is encoded yet.")
      max_qg = first(_top(gen_rows, :d_qg_result_vs_for002_MVar; n = 1))
      max_state_q = first(_top(residual_rows, :d_q_MVar; n = 1))
      max_kcl_q = first(_top(kcl_rows, :d_branch_kcl_vs_for002_q_MVar; n = 1))
      slack = gen_rows[findfirst(r -> r.is_slack, gen_rows)]
      println(io, "\n## Q / generator / bus-injection semantics\n")
      println(io, "- max |dQgen|: ", abs(max_qg.d_qg_result_vs_for002_MVar), " MVar at ", max_qg.bus_name)
      println(io, "- max |state residual Q|: ", abs(max_state_q.d_q_MVar), " MVar at ", max_state_q.bus_name)
      println(io, "- max |branch KCL vs FOR002 Q|: ", abs(max_kcl_q.d_branch_kcl_vs_for002_q_MVar), " MVar at ", max_kcl_q.bus_name)
      println(io, "- slack bus Q comparison (", slack.bus_name, "): FOR002 Qgen=", slack.for002_qg_MVar, " MVar; model specified Qgen=", slack.model_qg_specified_MVar, " MVar; model solved/result Qgen=", slack.model_qg_result_MVar, " MVar; dQ=", slack.d_qg_result_vs_for002_MVar, " MVar")
      largest_class = max_qg.is_slack ? "slack" : (max_qg.is_regulating ? "PV/regulating" : (max_qg.has_generator ? "PQ generator" : "non-generator"))
      transformer_adjacent = _norm_name(max_qg.bus_name) in Set(_norm_name.(["BETA1 S1", "BETA2 S1", "DELTA1S1", "DELTA2S1"]))
      println(io, "- largest Q-generator discrepancy class: ", largest_class, transformer_adjacent ? " and transformer-adjacent" : "")
      println(
        io,
        "\nBranch Q deviations remain small relative to the large bus/generator-Q diagnostics. This indicates that the branch-flow model is consistent with FOR002, while FOR002 bus-table and generator-Q reporting semantics still need interpretation. The KCL CSV documents that branch endpoint flows include branch charging/taps as represented by solved branch flows; explicit bus shunts are not added separately.",
      )
    end
    pushfirst!(written_files, summary_path)
  end
  result = DTFFor002ValidationResult(
    opt["output-dir"],
    converged,
    iters,
    final_mismatch,
    metrics_rows[1],
    written_files,
    metrics_rows[1].max_abs_branch_d_p_MW,
    metrics_rows[1].max_abs_branch_d_q_MVar,
    metrics_rows[1].max_abs_d_qg_MVar,
    metrics_rows[1].max_abs_state_residual_p_MW,
    metrics_rows[1].max_abs_state_residual_q_MVar,
  )
  (!opt["quiet"] && opt["print-summary"]) && _print_summary(result)
  # The default return value is lightweight for scripts/CLI use; detailed mode
  # returns diagnostic rows for tests and focused investigations.
  return return_details ?
         (
    case = case,
    net = net,
    metrics = metrics_rows[1],
    residuals = residual_rows,
    buses = bus_rows,
    branches = branch_rows,
    generator_rows = gen_rows,
    kcl_rows = kcl_rows,
    q_diagnostics = q_diag_rows,
    output_dir = opt["output-dir"],
    converged = converged,
    iterations = iters,
    final_mismatch = final_mismatch,
    transformer_ratio_mode = transformer_ratio_mode,
    written_files = written_files,
  ) : result
end


end # module NativeBaseValidation

module NativeOutageValidation
using Sparlectra
using Printf

include(joinpath(@__DIR__, "dtf_for002_validation_utils.jl"))

# Developer notes:
# - Validates native Testnetz13 DTF outage cards against FOR002 outage blocks.
# - The intentional execution path is DTFImporter.read_dtf -> DTFImporter.build_net -> runpf!.
# - MATPOWER import/export and the generated FOR001 builder are intentionally not used.
# - FOR002 is treated as a legacy textual reference report.
# - State residuals are diagnostic, not pass/fail thresholds.
# - Console output is intentionally concise; detailed rows are written to
#   CSV/Markdown and returned only in detailed mode.

struct DTFFor002OutageValidationResult
  output_dir::String
  parsed_outages::Int
  parsed_for002_outage_blocks::Int
  matching_rows::Vector{NamedTuple}
  metrics_rows::Vector{NamedTuple}
  written_files::Vector{String}
end

function _parse_bool(s::AbstractString)
  v = lowercase(strip(s))
  v in ("true", "1", "yes") && return true
  v in ("false", "0", "no") && return false
  throw(ArgumentError("invalid boolean: $s"))
end

function parse_cli(args)
  dtf_default = isfile(joinpath(@__DIR__, "..", "test", "fixtures", "dtf", "FOR001.DAT")) ? joinpath(@__DIR__, "..", "test", "fixtures", "dtf", "FOR001.DAT") : joinpath(@__DIR__, "FOR001.DAT")
  opt = Dict{String,Any}(
    "dtf-file" => normpath(dtf_default),
    "for002-file" => joinpath(@__DIR__, "FOR002.DAT"),
    "output-dir" => joinpath(@__DIR__, "_out", "dtf_for002_native_outages"),
    "tol" => 1e-8,
    "max-iter" => 50,
    "method" => "rectangular",
    "write-csv" => true,
    "write-markdown" => true,
    "quiet" => false,
    "strict" => false,
    "details" => false,
    "print-summary" => true,
  )
  for arg in args
    arg == "--quiet" && (opt["quiet"] = true; continue)
    startswith(arg, "--") || throw(ArgumentError("unsupported argument: $arg"))
    k, v = split(arg[3:end], "="; limit = 2)
    if k == "tol"
      ;
      opt[k] = parse(Float64, v)
    elseif k == "max-iter"
      ;
      opt[k] = parse(Int, v)
    elseif k in ("write-csv", "write-markdown", "strict", "quiet", "details", "print-summary")
      ;
      opt[k] = _parse_bool(v)
    else
      ;
      opt[k] = v
    end
  end
  return opt
end

function _method_symbol(s::AbstractString)
  s == "default" && return nothing
  s in ("rectangular", "polar") || throw(ArgumentError("--method must be rectangular, polar, or default"))
  return Symbol(s)
end

function _maxabs(rows, field)
  vals = Float64[]
  for r in rows
    v = getproperty(r, field)
    (v isa Missing || v === nothing) && continue
    fv = Float64(v);
    isfinite(fv) && push!(vals, abs(fv))
  end
  return isempty(vals) ? NaN : maximum(vals)
end

function _top(rows, field; n = 10)
  usable = [r for r in rows if !(getproperty(r, field) isa Missing)]
  return first(sort(usable, by = r -> -abs(Float64(getproperty(r, field)))), min(n, length(usable)))
end

_outage_label(o) = string(o.kind, o.voltage_level_index, o.parallel_id, " ", o.from, " -> ", o.to)

function _find_branch(case, outage)
  matches = Int[]
  for (i, b) in enumerate(case.branches)
    # Be strict: same branch kind, voltage level, parallel id, and endpoints.
    b.kind == outage.kind || continue
    b.voltage_level_index == outage.voltage_level_index || continue
    uppercase(strip(b.parallel_id)) == uppercase(strip(outage.parallel_id)) || continue
    _norm_name(b.from) == _norm_name(outage.from) || continue
    _norm_name(b.to) == _norm_name(outage.to) || continue
    push!(matches, i)
  end
  return matches
end

function _find_for002_scenario(scenarios, outage)
  candidates = Int[]
  for (i, s) in enumerate(scenarios)
    # FOR002 headings usually omit the DTF scenario/branch kind, so only the
    # available endpoint and parallel-id fields participate in text matching.
    _norm_name(s.from_bus === nothing ? "" : s.from_bus) == _norm_name(outage.from) || continue
    _norm_name(s.to_bus === nothing ? "" : s.to_bus) == _norm_name(outage.to) || continue
    uppercase(strip(s.parallel_id === nothing ? "" : s.parallel_id)) == uppercase(strip(outage.parallel_id)) || continue
    push!(candidates, i)
  end
  return candidates
end

function _apply_single_branch_outage!(net, idx::Int)
  before = [br.status for br in net.branchVec]
  before[idx] == 1 || throw(ArgumentError("matched branch $idx was not initially in service"))
  # Apply exactly one native branch status change and verify no other branch
  # moved; this protects the single-branch DTF outage-card interpretation.
  net.branchVec[idx].status = 0
  after = [br.status for br in net.branchVec]
  changed = findall(i -> before[i] != after[i], eachindex(before))
  length(changed) == 1 && only(changed) == idx || throw(ArgumentError("branch outage changed $(length(changed)) branches, expected exactly one"))
  net.branchVec[idx].status == 0 || throw(ArgumentError("matched branch $idx is still in service after outage"))
  return nothing
end

function _final_mismatch(net)
  try
    model = Sparlectra.buildPfModel(net; flatstart = false, include_limits = false, verbose = 0)
    V = [ComplexF64(n._vm_pu * cosd(n._va_deg), n._vm_pu * sind(n._va_deg)) for n in net.nodeVec]
    return Sparlectra.mismatchInf(model, V)
  catch
    return NaN
  end
end

function _comparison_rows(case, net, ref, outage_index, outage_label, outaged_idx)
  branch_kcl_p, branch_kcl_q = _branch_kcl_arrays(net)
  bus_rows = NamedTuple[];
  gen_rows = NamedTuple[]
  for db in case.buses
    node = net.nodeVec[db.index];
    key = _norm_name(db.name);
    haskey(ref.buses, key) || continue
    fb = ref.buses[key];
    flags = _bus_type_flags(db.bus_type);
    gl = _bus_generation_load(net, node.busIdx)
    push!(
      bus_rows,
      (
        outage_index = outage_index,
        outage_label = outage_label,
        bus_name = db.name,
        dtf_bus_type = db.bus_type,
        vn_kV = node.comp.cVN,
        for002_vm_kV = fb.v_kV,
        model_vm_kV = node.comp.cVN*node._vm_pu,
        d_vm_kV = node.comp.cVN*node._vm_pu-fb.v_kV,
        for002_vm_pu = fb.v_kV/node.comp.cVN,
        model_vm_pu = node._vm_pu,
        d_vm_pu = node._vm_pu-fb.v_kV/node.comp.cVN,
        for002_va_deg = fb.va_deg,
        model_va_deg = node._va_deg,
        d_va_deg = node._va_deg-fb.va_deg,
        is_slack = flags.is_slack,
        is_pv = flags.is_pv,
        is_pq = flags.is_pq,
      ),
    )
    solved_pg = (flags.is_slack || flags.is_pv) ? branch_kcl_p[node.busIdx] + gl.pl : gl.pg_res
    solved_qg = (flags.is_slack || flags.is_pv) ? branch_kcl_q[node.busIdx] + gl.ql : gl.qg_res
    push!(
      gen_rows,
      (
        outage_index = outage_index,
        outage_label = outage_label,
        bus_name = db.name,
        dtf_bus_type = db.bus_type,
        model_bus_type_effective = Sparlectra.toString(node._nodeType),
        for002_pg_MW = fb.p_gen_MW,
        model_pg_specified_MW = gl.pg_spec,
        model_pg_result_MW = solved_pg,
        d_pg_result_vs_for002_MW = solved_pg-fb.p_gen_MW,
        for002_qg_MVar = fb.q_gen_MVar,
        model_qg_specified_MVar = gl.qg_spec,
        model_qg_result_MVar = solved_qg,
        d_qg_result_vs_for002_MVar = solved_qg-fb.q_gen_MVar,
        qmin_MVar = something(db.qmin_mvar, missing),
        qmax_MVar = something(db.qmax_mvar, missing),
        is_regulating = gl.is_regulating || flags.is_pv,
        is_slack = flags.is_slack,
        has_generator = gl.has_generator,
        has_load = gl.has_load,
        p_load_MW = gl.pl,
        q_load_MVar = gl.ql,
        notes = flags.is_slack ? "slack bus: solved P/Q inferred from branch KCL plus load" : (gl.is_regulating || flags.is_pv ? "PV/regulating bus: result P/Q inferred from branch KCL plus load" : ""),
      ),
    )
  end
  rd = branch_ref_dict(ref);
  branch_rows = NamedTuple[]
  for (idx, br) in enumerate(net.branchVec)
    meta = net.matpower_branch_metadata[idx];
    dtfb = case.branches[idx];
    nr = uppercase(strip(String(meta.parallel_id)))
    rf = get(rd, (_norm_name(dtfb.from), _norm_name(dtfb.to), nr), missing);
    rt = get(rd, (_norm_name(dtfb.to), _norm_name(dtfb.from), nr), missing)
    pf = br.fBranchFlow === nothing || br.fBranchFlow.pFlow === nothing ? missing : br.fBranchFlow.pFlow;
    qf = br.fBranchFlow === nothing || br.fBranchFlow.qFlow === nothing ? missing : br.fBranchFlow.qFlow
    pt = br.tBranchFlow === nothing || br.tBranchFlow.pFlow === nothing ? missing : br.tBranchFlow.pFlow;
    qt = br.tBranchFlow === nothing || br.tBranchFlow.qFlow === nothing ? missing : br.tBranchFlow.qFlow
    push!(
      branch_rows,
      (
        outage_index = outage_index,
        outage_label = outage_label,
        branch_index = idx,
        branch_label = meta.orig_name,
        branch_kind = string(meta.dtf_kind),
        voltage_level_index = meta.voltage_level_index,
        parallel_id = meta.parallel_id,
        from_bus = dtfb.from,
        to_bus = dtfb.to,
        is_outaged_branch = idx==outaged_idx,
        branch_status = br.status,
        u_ref_kV = meta.u_ref_kV,
        r_pu = br.r_pu,
        x_pu = br.x_pu,
        b_pu = br.b_pu,
        ratio = br.ratio,
        for002_p_from_MW = rf===missing ? missing : rf.p_MW,
        model_p_from_MW = pf,
        d_p_from_MW = rf===missing || pf===missing ? missing : pf-rf.p_MW,
        for002_q_from_MVar = rf===missing ? missing : rf.q_MVar,
        model_q_from_MVar = qf,
        d_q_from_MVar = rf===missing || qf===missing ? missing : qf-rf.q_MVar,
        for002_p_to_MW = rt===missing ? missing : rt.p_MW,
        model_p_to_MW = pt,
        d_p_to_MW = rt===missing || pt===missing ? missing : pt-rt.p_MW,
        for002_q_to_MVar = rt===missing ? missing : rt.q_MVar,
        model_q_to_MVar = qt,
        d_q_to_MVar = rt===missing || qt===missing ? missing : qt-rt.q_MVar,
      ),
    )
  end
  kcl_base = branch_kcl_rows(net, case, ref)
  kcl_rows = [
    (
      outage_index = outage_index,
      outage_label = outage_label,
      bus_name = r.bus_name,
      dtf_bus_type = r.dtf_bus_type,
      is_slack = r.is_slack,
      is_pv = r.is_pv,
      is_pq = r.is_pq,
      for002_p_net_MW = r.for002_p_net_MW,
      for002_q_net_MVar = r.for002_q_net_MVar,
      model_result_p_net_MW = r.model_result_p_net_MW,
      model_result_q_net_MVar = r.model_result_q_net_MVar,
      branch_kcl_p_net_MW = r.branch_kcl_p_net_MW,
      branch_kcl_q_net_MVar = r.branch_kcl_q_net_MVar,
      d_model_result_vs_for002_p_MW = r.d_model_result_vs_for002_p_MW,
      d_model_result_vs_for002_q_MVar = r.d_model_result_vs_for002_q_MVar,
      d_branch_kcl_vs_for002_p_MW = r.d_branch_kcl_vs_for002_p_MW,
      d_branch_kcl_vs_for002_q_MVar = r.d_branch_kcl_vs_for002_q_MVar,
      notes = r.notes,
    ) for r in kcl_base
  ]
  res_rows = [
    (
      outage_index = outage_index,
      outage_label = outage_label,
      bus_name = r.bus_name,
      forced_for002_vm_pu = r.forced_for002_vm_pu,
      forced_for002_va_deg = r.forced_for002_va_deg,
      ybus_calc_p_MW = r.ybus_calc_p_MW,
      ybus_calc_q_MVar = r.ybus_calc_q_MVar,
      for002_bus_net_p_MW = r.for002_bus_net_p_MW,
      for002_bus_net_q_MVar = r.for002_bus_net_q_MVar,
      d_p_MW = r.d_p_MW,
      d_q_MVar = r.d_q_MVar,
    ) for r in for002_state_residual_rows(net, case, ref)
  ]
  return bus_rows, gen_rows, branch_rows, kcl_rows, res_rows
end

function _write_markdown(path, opt, case, for002_scenarios, matching_rows, metrics_rows, bus_rows, branch_rows, residual_rows)
  open(path, "w") do io
    println(io, "# Native DTF/FOR002 outage validation summary\n")
    println(io, "This validation uses `DTFImporter.read_dtf` -> `DTFImporter.build_net` for each outage and does not use MATPOWER import/export or the generated FOR001 builder.\n")
    println(io, "- DTF file: `", opt["dtf-file"], "`")
    println(io, "- FOR002 file: `", opt["for002-file"], "`")
    println(io, "- parsed DTF outages: ", length(case.outages))
    println(io, "- parsed FOR002 outage blocks: ", length(for002_scenarios), "\n")
    println(io, "## What are state residuals?\n")
    println(io, "State residuals force FOR002 printed voltage magnitudes/angles into the native Sparlectra outage Ybus. The resulting bus injections are compared with the FOR002 printed bus table. This is more sensitive than solved branch-flow comparisons because FOR002 values may be rounded and transformer-adjacent nodes react strongly to small voltage/angle differences. These residuals are diagnostic, not hard pass/fail criteria yet; branch-flow deviations and solved generator/slack comparisons are currently stronger validation signals.\n")
    println(io, "## Matching summary\n")
    for r in matching_rows
      println(io, "- [", r.outage_index, "] ", r.outage_label, ": branch=", r.matched_branch_index, ", FOR002=", r.matched_for002_scenario_index, ", status=", r.match_status)
    end
    for m in metrics_rows
      println(io, "\n## Outage ", m.outage_index, ": ", m.outage_label, "\n")
      println(io, "- converged: ", m.converged, "; iterations: ", m.iterations, "; final mismatch: ", m.final_mismatch)
      println(io, "- max |dV|: ", m.max_abs_d_vm_kV, " kV / ", m.max_abs_d_vm_pu, " pu; max |dVa|: ", m.max_abs_d_va_deg, " deg")
      println(io, "- max branch |dP|/|dQ|: ", m.max_abs_branch_d_p_MW, " MW / ", m.max_abs_branch_d_q_MVar, " MVar")
      for (title, rows, field, unit) in [
        ("Top voltage deviations", bus_rows, :d_vm_kV, "kV"),
        ("Top branch P deviations", branch_rows, :d_p_from_MW, "MW"),
        ("Top branch Q deviations", branch_rows, :d_q_from_MVar, "MVar"),
        ("Top state residual P deviations", residual_rows, :d_p_MW, "MW"),
        ("Top state residual Q deviations", residual_rows, :d_q_MVar, "MVar"),
      ]
        println(io, "\n### ", title, "\n")
        subset = [r for r in rows if r.outage_index == m.outage_index && !(getproperty(r, field) isa Missing)]
        for r in _top(subset, field; n = 5)
          name = hasproperty(r, :bus_name) ? r.bus_name : r.branch_label
          println(io, "- ", name, ": ", getproperty(r, field), " ", unit)
        end
      end
    end
    println(io, "\nState-residual rows force FOR002 outage voltages into the native outage Y-bus; remaining differences are diagnostic and not a pass/fail gate.")
  end
end

function _print_summary(result, opt)
  # Keep normal CLI output short; CSV/Markdown and detailed mode carry row-level
  # diagnostics for developers.
  println("Native DTF/FOR002 outage validation")
  println("DTF file: ", opt["dtf-file"])
  println("FOR002 file: ", opt["for002-file"])
  println("Output directory: ", result.output_dir)
  println("Parsed outages: ", result.parsed_outages)
  println("Matched FOR002 outage blocks: ", count(r -> r.matched_for002_scenario_index !== missing, result.matching_rows))
  println()
  by_index = Dict(r.outage_index => r for r in result.metrics_rows)
  for r in result.matching_rows
    println("[", r.outage_index, "] ", r.outage_label)
    println("    matched branch index: ", r.matched_branch_index)
    if haskey(by_index, r.outage_index)
      m = by_index[r.outage_index]
      println("    converged: ", m.converged)
      println("    iterations: ", m.iterations)
      println("    final mismatch: ", m.final_mismatch)
      println("    max |dV|: ", m.max_abs_d_vm_kV, " kV / ", m.max_abs_d_vm_pu, " pu")
      println("    max |dVa|: ", m.max_abs_d_va_deg, " deg")
      println("    max branch |dP|: ", m.max_abs_branch_d_p_MW, " MW")
      println("    max branch |dQ|: ", m.max_abs_branch_d_q_MVar, " MVar")
    else
      println("    match status: ", r.match_status)
      println("    notes: ", r.notes)
    end
    println()
  end
  println("Written files:")
  for path in result.written_files
    println("  - ", basename(path))
  end
end

function run_validation(args = ARGS; return_details::Bool = false)
  opt = parse_cli(args);
  return_details = return_details || opt["details"]
  mkpath(opt["output-dir"])
  case = Sparlectra.DTFImporter.read_dtf(opt["dtf-file"])
  for002_scenarios = parse_for002_outage_scenarios(opt["for002-file"])
  method = _method_symbol(opt["method"])
  matching_rows = NamedTuple[];
  metrics_rows = NamedTuple[];
  all_bus = NamedTuple[];
  all_gen = NamedTuple[];
  all_branch = NamedTuple[];
  all_kcl = NamedTuple[];
  all_res = NamedTuple[];
  outage_results = NamedTuple[]
  # Execute only DTF-listed outage cards. FOR002 may contain additional legacy
  # report blocks, but they are reference material until FOR001 requests them.
  for outage in case.outages
    label = _outage_label(outage);
    branch_matches = _find_branch(case, outage);
    scen_matches = _find_for002_scenario(for002_scenarios, outage)
    branch_idx = length(branch_matches) == 1 ? only(branch_matches) : missing
    scen_idx = length(scen_matches) == 1 ? only(scen_matches) : missing
    status = branch_idx === missing ? "no_unique_native_branch" : (scen_idx === missing ? "no_unique_for002_block" : "matched")
    notes = "native_branch_candidates=$(join(branch_matches, ';')); for002_candidates=$(join(scen_matches, ';'))"
    push!(
      matching_rows,
      (
        outage_index = outage.index,
        outage_label = label,
        outage_kind = string(outage.kind),
        voltage_level_index = outage.voltage_level_index,
        parallel_id = outage.parallel_id,
        from_bus = outage.from,
        to_bus = outage.to,
        matched_branch_index = branch_idx,
        matched_for002_scenario_index = scen_idx,
        matched_for002_heading = scen_idx===missing ? missing : for002_scenarios[scen_idx].raw_heading,
        match_status = status,
        notes = notes,
      ),
    )
    status == "matched" || continue
    # Build a fresh Net per outage so branch-status mutations and solved state
    # from one contingency cannot leak into the next.
    net = Sparlectra.DTFImporter.build_net(case)
    _apply_single_branch_outage!(net, branch_idx)
    iters, pfstatus = method === nothing ? Sparlectra.runpf!(net, opt["max-iter"], opt["tol"], 0) : Sparlectra.runpf!(net, opt["max-iter"], opt["tol"], 0; method = method)
    Sparlectra.calcNetLosses!(net)
    ref = _scenario_from_outage(for002_scenarios[scen_idx])
    bus_rows, gen_rows, branch_rows, kcl_rows, res_rows = _comparison_rows(case, net, ref, outage.index, label, branch_idx)
    append!(all_bus, bus_rows);
    append!(all_gen, gen_rows);
    append!(all_branch, branch_rows);
    append!(all_kcl, kcl_rows);
    append!(all_res, res_rows)
    metric = (
      outage_index = outage.index,
      outage_label = label,
      converged = pfstatus==0,
      iterations = iters,
      final_mismatch = _final_mismatch(net),
      max_abs_d_vm_kV = _maxabs(bus_rows, :d_vm_kV),
      max_abs_d_vm_pu = _maxabs(bus_rows, :d_vm_pu),
      max_abs_d_va_deg = _maxabs(bus_rows, :d_va_deg),
      max_abs_branch_d_p_MW = max(_maxabs(branch_rows, :d_p_from_MW), _maxabs(branch_rows, :d_p_to_MW)),
      max_abs_branch_d_q_MVar = max(_maxabs(branch_rows, :d_q_from_MVar), _maxabs(branch_rows, :d_q_to_MVar)),
      max_abs_d_pg_MW = _maxabs(gen_rows, :d_pg_result_vs_for002_MW),
      max_abs_d_qg_MVar = _maxabs(gen_rows, :d_qg_result_vs_for002_MVar),
      max_abs_state_residual_p_MW = _maxabs(res_rows, :d_p_MW),
      max_abs_state_residual_q_MVar = _maxabs(res_rows, :d_q_MVar),
      max_abs_branch_kcl_vs_for002_p_MW = _maxabs(kcl_rows, :d_branch_kcl_vs_for002_p_MW),
      max_abs_branch_kcl_vs_for002_q_MVar = _maxabs(kcl_rows, :d_branch_kcl_vs_for002_q_MVar),
    )
    push!(metrics_rows, metric)
    push!(outage_results, (outage = outage, net = net, metrics = metric, bus_rows = bus_rows, generator_rows = gen_rows, branch_rows = branch_rows, kcl_rows = kcl_rows, residual_rows = res_rows))
  end
  written_files = String[]
  if opt["write-csv"]
    for (filename, rows) in [
      ("dtf_outage_matching.csv", matching_rows),
      ("dtf_outage_metrics.csv", metrics_rows),
      ("dtf_outage_bus_comparison.csv", all_bus),
      ("dtf_outage_generator_comparison.csv", all_gen),
      ("dtf_outage_branch_comparison.csv", all_branch),
      ("dtf_outage_bus_kcl_comparison.csv", all_kcl),
      ("dtf_outage_state_residual.csv", all_res),
    ]
      path = joinpath(opt["output-dir"], filename);
      write_csv(path, rows);
      push!(written_files, path)
    end
  end
  if opt["write-markdown"]
    path = joinpath(opt["output-dir"], "dtf_outage_validation_summary.md")
    _write_markdown(path, opt, case, for002_scenarios, matching_rows, metrics_rows, all_bus, all_branch, all_res)
    pushfirst!(written_files, path)
  end
  result = DTFFor002OutageValidationResult(opt["output-dir"], length(case.outages), length(for002_scenarios), matching_rows, metrics_rows, written_files)
  (!opt["quiet"] && opt["print-summary"]) && _print_summary(result, opt)
  return return_details ?
         (
    case = case,
    for002_outage_scenarios = for002_scenarios,
    matching_rows = matching_rows,
    outage_results = outage_results,
    metrics_rows = metrics_rows,
    bus_rows = all_bus,
    generator_rows = all_gen,
    branch_rows = all_branch,
    kcl_rows = all_kcl,
    residual_rows = all_res,
    output_dir = opt["output-dir"],
    written_files = written_files,
  ) : result
end


end # module NativeOutageValidation

module MatpowerRoundtripValidation
using Logging
using Printf
using Sparlectra

include(joinpath(@__DIR__, "dtf_for002_validation_utils.jl"))

Base.@kwdef struct DTFMatpowerExportValidationResult
  output_dir::String
  scenario_results::Vector{NamedTuple}
  metrics_rows::Vector{NamedTuple}
  bus_rows::Vector{NamedTuple}
  branch_rows::Vector{NamedTuple}
  generator_rows::Vector{NamedTuple}
  written_files::Vector{String}
end

function parse_cli(args)
  opt = Dict{String,Any}(
    "dtf-file" => joinpath(@__DIR__, "..", "data", "DTF", "FOR001.DAT"),
    "for002-file" => joinpath(@__DIR__, "..", "data", "DTF", "FOR002.DAT"),
    "output-dir" => joinpath(@__DIR__, "_out", "dtf_matpower_export_testnetz13"),
    "tol" => 1e-8,
    "max-iter" => 50,
    "method" => "default",
    "write-csv" => true,
    "write-markdown" => true,
    "write-matpower" => true,
    "run-outages" => true,
    "quiet" => false,
    "print-summary" => true,
    "details" => false,
    "strict" => false,
  )
  for a in args
    startswith(a, "--") || continue
    k, v = split(a[3:end], "="; limit = 2)
    if haskey(opt, k)
      opt[k] = v in ("true", "false") ? v == "true" : (k in ("tol",) ? parse(Float64, v) : (k in ("max-iter",) ? parse(Int, v) : v))
    end
  end
  return opt
end

_method_symbol(s) = s == "default" ? nothing : Symbol(s)
_outage_label(o) = Sparlectra.DTFImporter.outage_label(o)
_final_mismatch(net) =
  try
    model = Sparlectra.buildPfModel(net; flatstart = false, include_limits = false, verbose = 0)
    V = [ComplexF64(n._vm_pu * cosd(n._va_deg), n._vm_pu * sind(n._va_deg)) for n in net.nodeVec]
    Sparlectra.mismatchInf(model, V)
  catch
    NaN
  end
_statuses(net) = [br.status for br in net.branchVec]
_maxabs(rows, field) = (v = [abs(Float64(getproperty(r, field))) for r in rows if !(getproperty(r, field) isa Missing) && !isnan(Float64(getproperty(r, field)))]; isempty(v) ? 0.0 : maximum(v))
_bname(net, i) = net.nodeVec[i].comp.cName
function _branch_name(net, i)
  meta = get(net.matpower_branch_metadata, i, nothing)
  return meta !== nothing && hasproperty(meta, :orig_name) && meta.orig_name !== nothing ? String(meta.orig_name) : net.branchVec[i].comp.cName
end

function _branch_kind(net, i)
  meta = get(net.matpower_branch_metadata, i, nothing)
  if meta !== nothing && hasproperty(meta, :dtf_kind)
    return string(meta.dtf_kind)
  elseif meta !== nothing && hasproperty(meta, :orig_kind) && meta.orig_kind !== nothing
    k = uppercase(String(meta.orig_kind))
    return k in ("T", "TRAFO", "TRANSFORMER") ? "T" : "L"
  end
  br = net.branchVec[i]
  return (br.ratio != 1.0 || br.angle != 0.0) ? "T" : "L"
end

function _shunt_totals(net)
  gs = bs = 0.0
  for sh in net.shuntVec
    sh.status == 1 || continue
    gs += real(sh.y_pu_shunt) * net.baseMVA
    bs += imag(sh.y_pu_shunt) * net.baseMVA
  end
  return (Gs = gs, Bs = bs)
end

function _metadata_summary(matpower_path)
  mpc = Sparlectra.MatpowerIO.read_case_m(matpower_path; legacy_compat = false)
  contingencies = mpc.for001_contingencies === nothing ? String[] : mpc.for001_contingencies
  return (
    bus_name_count = mpc.bus_name === nothing ? 0 : length(mpc.bus_name),
    branch_name_count = mpc.branch_name === nothing ? 0 : length(mpc.branch_name),
    branch_kind_count = mpc.branch_kind === nothing ? 0 : length(mpc.branch_kind),
    for001_contingency_count = length(contingencies),
    for001_contingencies = contingencies,
    tap_line_values = mpc.branch_kind === nothing ? Float64[] : [mpc.branch[i, 9] for i in axes(mpc.branch, 1) if uppercase(mpc.branch_kind[i]) == "L"],
    tap_transformer_values = mpc.branch_kind === nothing ? Float64[] : [mpc.branch[i, 9] for i in axes(mpc.branch, 1) if uppercase(mpc.branch_kind[i]) == "T"],
    metadata_complete = mpc.bus_name !== nothing && length(mpc.bus_name) == 13 && mpc.branch_name !== nothing && length(mpc.branch_name) == 27 && mpc.branch_kind !== nothing && length(mpc.branch_kind) == 27 && length(contingencies) >= 2,
  )
end

function _run_pf!(net, opt)
  method = _method_symbol(opt["method"])
  r = method === nothing ? Sparlectra.runpf!(net, opt["max-iter"], opt["tol"], 0) : Sparlectra.runpf!(net, opt["max-iter"], opt["tol"], 0; method = method)
  Sparlectra.calcNetLosses!(net)
  return (iterations = r[1], converged = r[2] == 0, final_mismatch = _final_mismatch(net))
end

function _roundtrip(native_net, matpower_path, opt)
  # Existing Sparlectra MATPOWER exporter; this is intentionally Net -> MATPOWER,
  # not a DTF-specific exporter.
  Logging.with_logger(Logging.NullLogger()) do
    Sparlectra.writeMatpowerCasefile(native_net, matpower_path)
  end
  rt = Logging.with_logger(Logging.NullLogger()) do
    # DTF PQ generators are fixed injections. Disable MATPOWER PQ generator
    # controller reinterpretation for this roundtrip so the imported Net
    # preserves native DTF generator semantics instead of voltage-dependent
    # controller behavior.
    Sparlectra.createNetFromMatPowerFile(filename = matpower_path, log = false, flatstart = false, enable_pq_gen_controllers = false, apply_bus_names = true, apply_branch_names = true, apply_branch_kind = true)
  end
  pf = _run_pf!(rt, opt)
  return rt, pf
end

function _bus_rows(scenario, outage_index, outage_label, native, rt)
  rows = NamedTuple[]
  for i in eachindex(native.nodeVec)
    n = native.nodeVec[i];
    r = rt.nodeVec[i]
    push!(
      rows,
      (
        scenario = scenario,
        outage_index = outage_index,
        outage_label = outage_label,
        bus_index_native = i,
        bus_index_roundtrip = i,
        bus_name_native = n.comp.cName,
        bus_name_roundtrip = r.comp.cName,
        baseKV_native = n.comp.cVN,
        baseKV_roundtrip = r.comp.cVN,
        vm_pu_native = n._vm_pu,
        vm_pu_roundtrip = r._vm_pu,
        d_vm_pu = r._vm_pu - n._vm_pu,
        vm_kV_native = n.comp.cVN*n._vm_pu,
        vm_kV_roundtrip = r.comp.cVN*r._vm_pu,
        d_vm_kV = r.comp.cVN*r._vm_pu - n.comp.cVN*n._vm_pu,
        va_deg_native = n._va_deg,
        va_deg_roundtrip = r._va_deg,
        d_va_deg = r._va_deg - n._va_deg,
        bus_type_native = Sparlectra.toString(n._nodeType),
        bus_type_roundtrip = Sparlectra.toString(r._nodeType),
        notes = "",
      ),
    )
  end
  return rows
end

function _branch_rows(scenario, outage_index, outage_label, native, rt, outaged_idx)
  rows = NamedTuple[]
  rt_by_name = Dict(_branch_name(rt, i) => i for i in eachindex(rt.branchVec))
  for i in eachindex(native.branchVec)
    b = native.branchVec[i]
    native_name = _branch_name(native, i)
    j = get(rt_by_name, native_name, i)
    rb = rt.branchVec[j]
    pf = b.fBranchFlow === nothing ? missing : b.fBranchFlow.pFlow;
    qf = b.fBranchFlow === nothing ? missing : b.fBranchFlow.qFlow
    pt = b.tBranchFlow === nothing ? missing : b.tBranchFlow.pFlow;
    qt = b.tBranchFlow === nothing ? missing : b.tBranchFlow.qFlow
    rpf = rb.fBranchFlow === nothing ? missing : rb.fBranchFlow.pFlow;
    rqf = rb.fBranchFlow === nothing ? missing : rb.fBranchFlow.qFlow
    rpt = rb.tBranchFlow === nothing ? missing : rb.tBranchFlow.pFlow;
    rqt = rb.tBranchFlow === nothing ? missing : rb.tBranchFlow.qFlow
    isout = outaged_idx !== nothing && i == outaged_idx
    push!(
      rows,
      (
        scenario = scenario,
        outage_index = outage_index,
        outage_label = outage_label,
        branch_index_native = i,
        branch_index_roundtrip = j,
        branch_name_native = native_name,
        branch_name_roundtrip = _branch_name(rt, j),
        from_bus_native = _bname(native, b.fromBus),
        to_bus_native = _bname(native, b.toBus),
        from_bus_roundtrip = _bname(rt, rb.fromBus),
        to_bus_roundtrip = _bname(rt, rb.toBus),
        r_pu_native = b.r_pu,
        r_pu_roundtrip = rb.r_pu,
        d_r_pu = rb.r_pu-b.r_pu,
        x_pu_native = b.x_pu,
        x_pu_roundtrip = rb.x_pu,
        d_x_pu = rb.x_pu-b.x_pu,
        b_pu_native = b.b_pu,
        b_pu_roundtrip = rb.b_pu,
        d_b_pu = rb.b_pu-b.b_pu,
        g_pu_native = b.g_pu,
        g_pu_roundtrip = rb.g_pu,
        d_g_pu = rb.g_pu-b.g_pu,
        branch_kind_native = _branch_kind(native, i),
        branch_kind_roundtrip = _branch_kind(rt, j),
        exported_tap = _branch_kind(native, i) == "T" ? b.ratio : 0.0,
        ratio_native = b.ratio,
        ratio_roundtrip = rb.ratio,
        d_ratio = rb.ratio-b.ratio,
        status_native = b.status,
        status_roundtrip = rb.status,
        p_from_native_MW = pf,
        p_from_roundtrip_MW = rpf,
        d_p_from_MW = (pf isa Missing || rpf isa Missing || isout) ? missing : rpf-pf,
        q_from_native_MVar = qf,
        q_from_roundtrip_MVar = rqf,
        d_q_from_MVar = (qf isa Missing || rqf isa Missing || isout) ? missing : rqf-qf,
        p_to_native_MW = pt,
        p_to_roundtrip_MW = rpt,
        d_p_to_MW = (pt isa Missing || rpt isa Missing || isout) ? missing : rpt-pt,
        q_to_native_MVar = qt,
        q_to_roundtrip_MVar = rqt,
        d_q_to_MVar = (qt isa Missing || rqt isa Missing || isout) ? missing : rqt-qt,
        is_outaged_branch = isout,
        notes = isout ? "outaged branch; flow deltas excluded" : "",
      ),
    )
  end
  return rows
end

function _gen_rows(scenario, outage_index, outage_label, native, rt)
  rows = NamedTuple[]
  for i in eachindex(native.nodeVec)
    ng = _bus_generation_load(native, i);
    rg = _bus_generation_load(rt, i)
    (ng.has_generator || ng.has_load || rg.has_generator || rg.has_load) || continue
    nn = native.nodeVec[i];
    rn = rt.nodeVec[i]
    push!(
      rows,
      (
        scenario = scenario,
        outage_index = outage_index,
        outage_label = outage_label,
        bus_name_native = nn.comp.cName,
        bus_name_roundtrip = rn.comp.cName,
        bus_index_native = i,
        bus_index_roundtrip = i,
        pg_native_MW = ng.pg_res,
        pg_roundtrip_MW = rg.pg_res,
        d_pg_MW = rg.pg_res-ng.pg_res,
        qg_native_MVar = ng.qg_res,
        qg_roundtrip_MVar = rg.qg_res,
        d_qg_MVar = rg.qg_res-ng.qg_res,
        p_load_native_MW = ng.pl,
        p_load_roundtrip_MW = rg.pl,
        d_p_load_MW = rg.pl-ng.pl,
        q_load_native_MVar = ng.ql,
        q_load_roundtrip_MVar = rg.ql,
        d_q_load_MVar = rg.ql-ng.ql,
        is_slack_native = nn._nodeType == Sparlectra.Slack,
        is_slack_roundtrip = rn._nodeType == Sparlectra.Slack,
        is_regulating_native = ng.is_regulating,
        is_regulating_roundtrip = rg.is_regulating,
        notes = "",
      ),
    )
  end
  return rows
end

function _scenario(case, opt, scenario, outage, outaged_idx)
  native = Sparlectra.DTFImporter.build_net(case)
  if outaged_idx !== nothing
    Sparlectra.DTFImporter.apply_single_branch_outage!(native, outaged_idx)
  end
  npf = _run_pf!(native, opt)
  native_shunts = _shunt_totals(native)
  mp = joinpath(opt["output-dir"], scenario * ".m")
  rt, rpf = opt["write-matpower"] ? _roundtrip(native, mp, opt) : error("write-matpower=false is not supported for roundtrip validation")
  roundtrip_shunts = _shunt_totals(rt)
  metadata = _metadata_summary(mp)
  buses = _bus_rows(scenario, outage === nothing ? missing : outage.index, outage === nothing ? "" : _outage_label(outage), native, rt)
  branches = _branch_rows(scenario, outage === nothing ? missing : outage.index, outage === nothing ? "" : _outage_label(outage), native, rt, outaged_idx)
  gens = _gen_rows(scenario, outage === nothing ? missing : outage.index, outage === nothing ? "" : _outage_label(outage), native, rt)
  metric = (
    scenario = scenario,
    outage_index = outage === nothing ? missing : outage.index,
    outage_label = outage === nothing ? "" : _outage_label(outage),
    native_converged = npf.converged,
    roundtrip_converged = rpf.converged,
    native_iterations = npf.iterations,
    roundtrip_iterations = rpf.iterations,
    native_final_mismatch = npf.final_mismatch,
    roundtrip_final_mismatch = rpf.final_mismatch,
    max_abs_d_vm_pu = _maxabs(buses, :d_vm_pu),
    max_abs_d_vm_kV = _maxabs(buses, :d_vm_kV),
    max_abs_d_va_deg = _maxabs(buses, :d_va_deg),
    max_abs_branch_d_p_MW = max(_maxabs(branches, :d_p_from_MW), _maxabs(branches, :d_p_to_MW)),
    max_abs_branch_d_q_MVar = max(_maxabs(branches, :d_q_from_MVar), _maxabs(branches, :d_q_to_MVar)),
    max_abs_branch_r_pu = _maxabs(branches, :d_r_pu),
    max_abs_branch_x_pu = _maxabs(branches, :d_x_pu),
    max_abs_branch_b_pu = _maxabs(branches, :d_b_pu),
    max_abs_ratio = _maxabs(branches, :d_ratio),
    max_abs_branch_g_pu_native = maximum([abs(br.g_pu) for br in native.branchVec]; init = 0.0),
    count_nonzero_branch_g_pu_native = count(br -> abs(br.g_pu) > 1e-12, native.branchVec),
    native_total_bus_Gs = native_shunts.Gs,
    native_total_bus_Bs = native_shunts.Bs,
    exported_total_bus_Gs = native_shunts.Gs,
    exported_total_bus_Bs = native_shunts.Bs,
    roundtrip_total_bus_Gs = roundtrip_shunts.Gs,
    roundtrip_total_bus_Bs = roundtrip_shunts.Bs,
    bus_count_native = length(native.nodeVec),
    bus_count_roundtrip = length(rt.nodeVec),
    branch_count_native = length(native.branchVec),
    branch_count_roundtrip = length(rt.branchVec),
    generator_count_native = length(native.prosumpsVec),
    generator_count_roundtrip = length(rt.prosumpsVec),
    status_match = _statuses(native) == _statuses(rt),
    notes = metadata.metadata_complete ? "existing exporter: writeMatpowerCasefile; established metadata present" : "existing exporter: writeMatpowerCasefile; metadata incomplete",
  )
  return (native = native, roundtrip = rt, matpower_file = mp, metadata = metadata, metric = metric, bus_rows = buses, branch_rows = branches, generator_rows = gens)
end

function _write_summary(path, opt, results)
  open(path, "w") do io
    println(io, "# DTF MATPOWER export validation\n")
    println(io, "This diagnostic uses `DTFImporter.read_dtf` -> `DTFImporter.build_net` -> `Sparlectra.Net` -> existing `writeMatpowerCasefile` -> `createNetFromMatPowerFile`. It does not implement a DTF-specific MATPOWER exporter.\n")
    println(io, "- native DTF input: `", opt["dtf-file"], "`")
    println(io, "- output directory: `", opt["output-dir"], "`")
    println(io, "- metadata fields: `mpc.bus_name`, `mpc.branch_name`, `mpc.branch_kind`, and `mpc.for001_contingencies` are exported when available; solving does not require them.\n")
    println(io, "- roundtrip import option: `enable_pq_gen_controllers=false`, because native DTF PQ generators are fixed injections and enabling MATPOWER PQ generator controllers changes the solved roundtrip semantics.\n")
    println(io, "- TAP convention: line rows export TAP = 0.0; transformer rows export their explicit ratio. This keeps MATPOWER line rows from being re-imported as nominal-tap transformer rows.\n")
    for r in results
      m = r.metric
      println(io, "## ", m.scenario)
      println(io, "- exported MATPOWER: `", r.matpower_file, "`")
      println(io, "- metadata counts: bus_name=", r.metadata.bus_name_count, ", branch_name=", r.metadata.branch_name_count, ", branch_kind=", r.metadata.branch_kind_count, ", for001_contingencies=", r.metadata.for001_contingency_count)
      println(io, "- exported line TAP values: `", join(r.metadata.tap_line_values, ", "), "`")
      println(io, "- exported transformer TAP values: `", join(r.metadata.tap_transformer_values, ", "), "`")
      println(io, "- metadata complete: ", r.metadata.metadata_complete)
      if !r.metadata.metadata_complete
        println(io, "- metadata limitation: one or more established MATPOWER metadata fields were not preserved with the expected Testnetz13 counts")
      end
      println(io, "- native converged: ", m.native_converged)
      println(io, "- roundtrip converged: ", m.roundtrip_converged)
      println(io, "- max |dV| pu: ", m.max_abs_d_vm_pu)
      println(io, "- max branch |dP| MW: ", m.max_abs_branch_d_p_MW)
      println(io, "- max branch |dQ| MVar: ", m.max_abs_branch_d_q_MVar)
      println(io, "- max native branch g_pu: ", m.max_abs_branch_g_pu_native, " (count > 1e-12: ", m.count_nonzero_branch_g_pu_native, ")")
      if m.count_nonzero_branch_g_pu_native == 0
        println(io, "- branch shunt conductance: all native branch `g_pu` values are zero/negligible")
      else
        affected = [string(row.branch_index_native, ":", row.branch_name_native, "=", row.g_pu_native) for row in r.branch_rows if abs(row.g_pu_native) > 1e-12]
        println(io, "- branch shunt conductance limitation: standard MATPOWER branch data cannot preserve branch `g_pu` exactly; equivalent endpoint bus `Gs` is exported through MATPOWER bus shunts for Ybus preservation. Affected branches: ", join(affected, "; "))
      end
      println(io, "- bus shunt totals native/exported/roundtrip Gs: ", m.native_total_bus_Gs, " / ", m.exported_total_bus_Gs, " / ", m.roundtrip_total_bus_Gs)
      println(io, "- bus shunt totals native/exported/roundtrip Bs: ", m.native_total_bus_Bs, " / ", m.exported_total_bus_Bs, " / ", m.roundtrip_total_bus_Bs)
      println(io, "- status match: ", m.status_match, "\n")
    end
  end
end

function _print_summary(result, opt)
  println("Native DTF -> existing MATPOWER export validation")
  println("DTF file: ", opt["dtf-file"])
  println("Output directory: ", result.output_dir)
  println("Existing exporter: writeMatpowerCasefile")
  for r in result.scenario_results
    m = r.metric
    println(m.scenario == "base" ? "Base case:" : "Outage [$(m.outage_index)] $(m.outage_label):")
    println("    native converged: ", m.native_converged)
    println("    roundtrip converged: ", m.roundtrip_converged)
    println("    max |dV|: ", m.max_abs_d_vm_pu, " pu / ", m.max_abs_d_vm_kV, " kV")
    println("    max branch |dP|: ", m.max_abs_branch_d_p_MW, " MW")
    println("    max branch |dQ|: ", m.max_abs_branch_d_q_MVar, " MVar")
  end
  println("Written files:")
  foreach(p -> println("    ", p), result.written_files)
end

function run_validation(args = ARGS; return_details::Bool = false)
  opt = parse_cli(args);
  return_details |= opt["details"]
  mkpath(opt["output-dir"])
  case = Sparlectra.DTFImporter.read_dtf(opt["dtf-file"])
  results = NamedTuple[]
  push!(results, _scenario(case, opt, "base", nothing, nothing))
  if opt["run-outages"]
    for outage in case.outages
      matches = Sparlectra.DTFImporter.find_outage_branch_indices(case, outage)
      length(matches) == 1 || throw(ArgumentError(Sparlectra.DTFImporter.outage_match_diagnostic(case, outage, matches)))
      push!(results, _scenario(case, opt, "outage_$(outage.index)", outage, only(matches)))
    end
  end
  metrics = [r.metric for r in results];
  buses = reduce(vcat, [r.bus_rows for r in results]; init = NamedTuple[]);
  branches = reduce(vcat, [r.branch_rows for r in results]; init = NamedTuple[]);
  gens = reduce(vcat, [r.generator_rows for r in results]; init = NamedTuple[])
  written = String[]
  if opt["write-csv"]
    for (fn, rows) in [("dtf_matpower_export_metrics.csv", metrics), ("dtf_matpower_export_bus_comparison.csv", buses), ("dtf_matpower_export_branch_comparison.csv", branches), ("dtf_matpower_export_generator_comparison.csv", gens)]
      p = joinpath(opt["output-dir"], fn);
      write_csv(p, rows);
      push!(written, p)
    end
  end
  if opt["write-markdown"]
    p = joinpath(opt["output-dir"], "dtf_matpower_export_summary.md");
    _write_summary(p, opt, results);
    push!(written, p)
  end
  append!(written, [r.matpower_file for r in results if isfile(r.matpower_file)])
  result = DTFMatpowerExportValidationResult(opt["output-dir"], results, metrics, buses, branches, gens, written)
  !opt["quiet"] && opt["print-summary"] && _print_summary(result, opt)
  if opt["strict"]
    all(m -> m.native_converged && m.roundtrip_converged && m.status_match, metrics) || error("roundtrip validation failed strict checks")
  end
  return return_details ? result : (output_dir = result.output_dir, metrics_rows = result.metrics_rows, written_files = result.written_files)
end


end # module MatpowerRoundtripValidation

module ImportAudit
using Printf
using Sparlectra

const CASES = ["A" => "", "B" => "B", "C" => "C", "D" => "D", "E" => "E"]

finite_range(xs) = isempty(xs) ? "n/a" : @sprintf("[%g, %g], finite=%s", minimum(xs), maximum(xs), all(isfinite, xs))

function audit_case(io, case_id, path; strict::Bool = true)
  case = Sparlectra.DTFImporter.read_dtf(path; strict = strict)
  println(io, "## Case ", case_id)
  println(io, "- file: ", path)
  println(io, "- parsed branch count: ", length(case.branches))
  println(io, "- parsed bus count: ", length(case.buses))
  println(io, "- parsed transformer control count: ", length(case.transformer_controls))
  println(io, "- parsed outage/trailing-section count: ", length(case.outages), " outages; ", length(case.trailing_records), " trailing records")
  println(io, "- slack bus: ", case.size.slack)
  println(io, "- PV buses: ", join((b.name for b in case.buses if b.bus_type == 1), ", "))
  println(io, "- PQ buses: ", join((b.name for b in case.buses if b.bus_type in (0, 3)), ", "))
  println(io, "- branch r ohm: ", finite_range([b.r_ohm for b in case.branches]))
  println(io, "- branch x ohm: ", finite_range([b.x_ohm for b in case.branches]))
  println(io, "- branch g S: ", finite_range([b.g_s for b in case.branches]))
  println(io, "- branch b S: ", finite_range([b.b_s for b in case.branches]))
  println(io, "\n### Transformer controls")
  for c in case.transformer_controls
    br_idx = findfirst(b -> b.kind == 'T' && b.from == c.from && b.to == c.to && b.parallel_id == c.parallel_id, case.branches)
    br = br_idx === nothing ? nothing : case.branches[br_idx]
    tap = br === nothing ? nothing : Sparlectra.DTFImporter._dtf_effective_transformer_tap(case, br, c,
      first(b for b in case.buses if b.name == c.from),
      first(b for b in case.buses if b.name == c.to))
    compat_tap = br === nothing ? nothing : Sparlectra.DTFImporter._dtf_effective_transformer_tap(case, br, c,
      first(b for b in case.buses if b.name == c.from),
      first(b for b in case.buses if b.name == c.to);
      transformer_ratio_mode = :winding_over_network)
    println(io, "- source line ", c.index, ": ", c.from, " -> ", c.to, " parallel=", c.parallel_id,
      " raw=", repr(c.raw),
      " regulated_side=", c.regulated_side,
      " unregulated_side=", c.unregulated_side,
      " phase_skew_flag=", c.phase_shifter_flag,
      " nominal_unregulated_kV=", c.nominal_unregulated_kv,
      " nominal_regulated_kV=", c.nominal_regulated_kv,
      " longitudinal_range_percent=", c.longitudinal_range_percent,
      " actual_tap_step=", c.actual_tap_step,
      " max_tap_step=", c.max_tap_step,
      " skew_angle_deg=", c.added_voltage_angle_deg,
      " quadrature_range_percent=", c.quadrature_range_percent,
      " quadrature_max_steps=", c.quadrature_max_steps,
      " quadrature_actual_step=", c.quadrature_actual_step,
      tap === nothing ? "" : string(
        " computed_model=", tap.model,
        " transformer_ratio_mode=", tap.transformer_ratio_mode,
        " base_ratio_used=", tap.base_ratio_used,
        " winding_over_network_base_ratio=", tap.winding_over_network_base_ratio,
        " tap_fraction=", tap.tap_fraction,
        " effective_complex=", tap.effective_complex,
        " effective_ratio=", tap.ratio,
        " effective_shift_deg=", tap.shift_deg,
        " sign_reciprocal_convention=", tap.convention,
        compat_tap === nothing ? "" : string(" compatibility_winding_over_network_ratio=", compat_tap.ratio)))
  end
  println(io, "\n### Node start voltages")
  for b in case.buses
    vn = case.nominal_voltages_kv[b.voltage_level_index]
    vm = b.start_kv > 0 ? b.start_kv / vn : 1.0
    println(io, "- ", b.name, ": start_kV=", b.start_kv, " vn_kV=", vn, " vm_pu=", vm, " finite=", all(isfinite, (b.start_kv, vn, vm)))
  end
  println(io, "\n### Post-bus/trailing records")
  if isempty(case.trailing_records)
    println(io, "- none")
  else
    for r in case.trailing_records
      println(io, "- raw=", repr(r.raw), "; interpreted kind=", r.interpreted_kind)
    end
  end
  println(io)
  return nothing
end


end # module ImportAudit


const _SUITE_MODES = ("audit", "base", "outages", "matpower")

function _parse_bool(value::AbstractString)
  normalized = lowercase(strip(value))
  normalized in ("true", "1", "yes", "on") && return true
  normalized in ("false", "0", "no", "off") && return false
  throw(ArgumentError("invalid boolean value: $value"))
end

function _default_suite_options()
  repo = normpath(joinpath(@__DIR__, ".."))
  return Dict{String,Any}(
    "mode" => "all",
    "case" => "all",
    "data-dir" => joinpath(repo, "data", "DTF"),
    "output-dir" => joinpath(@__DIR__, "_out", "dtf_validation_suite"),
    "tol" => 1e-8,
    "max-iter" => 50,
    "method" => "rectangular",
    "write-csv" => true,
    "write-markdown" => true,
    "write-matpower" => true,
    "run-outages" => true,
    "quiet" => false,
    "strict" => false,
    "parser-strict" => false,
    "continue-on-error" => true,
    "legacy-voltage-level-collapse-230kv" => false,
    "transformer-ratio-mode" => "neutral_one",
  )
end

function _print_help()
  println("""
Unified DTF validation suite

Usage:
  julia --project=. examples/validate_dtf_suite.jl [options]

Main options:
  --mode=all|audit|base|outages|matpower
      Multiple modes may be comma-separated, for example --mode=audit,base.
  --case=all|A|B|C|D|E
      Multiple cases may be comma-separated. Case A maps to FOR001.DAT/FOR002.DAT.
  --data-dir=<path>
  --output-dir=<path>

Solver options:
  --tol=1e-8
  --max-iter=50
  --method=rectangular|polar|default
  --transformer-ratio-mode=neutral_one|winding_over_network
  --legacy-voltage-level-collapse-230kv=true|false

Output and control:
  --write-csv=true|false
  --write-markdown=true|false
  --write-matpower=true|false
  --run-outages=true|false
  --parser-strict=true|false
  --strict=true|false
  --continue-on-error=true|false
  --quiet
  --help

File pairing:
  FOR001.DAT  <-> FOR002.DAT   (case A)
  FOR001B.DAT <-> FOR002B.DAT  (case B)
  Additional matching suffixes are detected automatically.
""")
end

function parse_suite_cli(args = ARGS)
  opt = _default_suite_options()
  boolean_keys = Set([
    "write-csv",
    "write-markdown",
    "write-matpower",
    "run-outages",
    "quiet",
    "strict",
    "parser-strict",
    "continue-on-error",
    "legacy-voltage-level-collapse-230kv",
  ])
  for arg in args
    if arg in ("--help", "-h")
      _print_help()
      return nothing
    elseif arg == "--quiet"
      opt["quiet"] = true
      continue
    elseif arg == "--strict"
      opt["strict"] = true
      continue
    elseif arg == "--no-continue-on-error"
      opt["continue-on-error"] = false
      continue
    end

    startswith(arg, "--") || throw(ArgumentError("unsupported argument: $arg"))
    parts = split(arg[3:end], "="; limit = 2)
    length(parts) == 2 || throw(ArgumentError("expected --key=value, got: $arg"))
    key, value = parts
    haskey(opt, key) || throw(ArgumentError("unknown option --$key"))

    if key == "tol"
      opt[key] = parse(Float64, value)
    elseif key == "max-iter"
      opt[key] = parse(Int, value)
    elseif key in boolean_keys
      opt[key] = _parse_bool(value)
    else
      opt[key] = value
    end
  end

  opt["data-dir"] = normpath(abspath(opt["data-dir"]))
  opt["output-dir"] = normpath(abspath(opt["output-dir"]))
  opt["tol"] > 0 || throw(ArgumentError("--tol must be positive"))
  opt["max-iter"] > 0 || throw(ArgumentError("--max-iter must be positive"))
  return opt
end

function _selected_modes(value::AbstractString)
  raw = lowercase.(strip.(split(value, ",")))
  "all" in raw && return collect(_SUITE_MODES)
  aliases = Dict("native" => "base", "outage" => "outages", "roundtrip" => "matpower")
  modes = String[]
  for token in raw
    isempty(token) && continue
    mode = get(aliases, token, token)
    mode in _SUITE_MODES || throw(ArgumentError("unknown mode '$token'; allowed: all, " * join(_SUITE_MODES, ", ")))
    mode in modes || push!(modes, mode)
  end
  isempty(modes) && throw(ArgumentError("no validation mode selected"))
  return modes
end

function _selected_cases(value::AbstractString)
  tokens = uppercase.(strip.(split(value, ",")))
  "ALL" in tokens && return nothing
  selected = Set(filter(token -> !isempty(token), tokens))
  isempty(selected) && throw(ArgumentError("no case selected"))
  return selected
end

function _discover_for_files(data_dir::AbstractString)
  isdir(data_dir) || throw(ArgumentError("DTF data directory does not exist: $data_dir"))
  for001 = Dict{String,String}()
  for002 = Dict{String,String}()

  for name in readdir(data_dir)
    path = joinpath(data_dir, name)
    isfile(path) || continue

    m1 = match(r"^FOR001(.*)\.DAT$"i, name)
    if m1 !== nothing
      suffix = uppercase(String(something(m1.captures[1], "")))
      haskey(for001, suffix) && throw(ArgumentError("duplicate FOR001 suffix '$suffix': $(for001[suffix]) and $path"))
      for001[suffix] = path
      continue
    end

    m2 = match(r"^FOR002(.*)\.DAT$"i, name)
    if m2 !== nothing
      suffix = uppercase(String(something(m2.captures[1], "")))
      haskey(for002, suffix) && throw(ArgumentError("duplicate FOR002 suffix '$suffix': $(for002[suffix]) and $path"))
      for002[suffix] = path
    end
  end

  suffixes = sort!(collect(union(keys(for001), keys(for002))); by = s -> (isempty(s) ? 0 : 1, s))
  return [(suffix = s, case_id = isempty(s) ? "A" : s, dtf_file = get(for001, s, missing), for002_file = get(for002, s, missing)) for s in suffixes]
end

_safe_case_dir(case_id::AbstractString) = "case_" * replace(case_id, r"[^A-Za-z0-9_.-]+" => "_")

function _csv_cell(value)
  value isa Missing && return ""
  value === nothing && return ""
  text = string(value)
  if occursin(',', text) || occursin('"', text) || occursin('\n', text) || occursin('\r', text)
    return "\"" * replace(text, "\"" => "\"\"") * "\""
  end
  return text
end

function _write_csv(path::AbstractString, rows)
  mkpath(dirname(path))
  open(path, "w") do io
    isempty(rows) && return
    columns = propertynames(first(rows))
    println(io, join(string.(columns), ","))
    for row in rows
      println(io, join((_csv_cell(getproperty(row, column)) for column in columns), ","))
    end
  end
  return path
end

function _finite_max(values)
  usable = Float64[]
  for value in values
    (value isa Missing || value === nothing) && continue
    number = try
      Float64(value)
    catch
      continue
    end
    isfinite(number) && push!(usable, abs(number))
  end
  return isempty(usable) ? missing : maximum(usable)
end

function _md(value)
  value isa Missing && return ""
  value === nothing && return ""
  return replace(string(value), "|" => "\\|", "\n" => " ", "\r" => " ")
end

function _summary_row(;
  case_id,
  suffix,
  mode,
  status,
  dtf_file = missing,
  for002_file = missing,
  output_dir = missing,
  converged = missing,
  scenario_count = missing,
  max_abs_d_vm_kV = missing,
  max_abs_d_va_deg = missing,
  max_abs_d_p_MW = missing,
  max_abs_d_q_MVar = missing,
  message = "",
)
  return (
    case_id = case_id,
    suffix = suffix,
    mode = mode,
    status = status,
    converged = converged,
    scenario_count = scenario_count,
    max_abs_d_vm_kV = max_abs_d_vm_kV,
    max_abs_d_va_deg = max_abs_d_va_deg,
    max_abs_d_p_MW = max_abs_d_p_MW,
    max_abs_d_q_MVar = max_abs_d_q_MVar,
    dtf_file = dtf_file,
    for002_file = for002_file,
    output_dir = output_dir,
    message = message,
  )
end

function _common_validation_args(opt, dtf_file, for002_file, output_dir)
  return String[
    "--dtf-file=$(dtf_file)",
    "--for002-file=$(for002_file)",
    "--output-dir=$(output_dir)",
    "--tol=$(opt["tol"])",
    "--max-iter=$(opt["max-iter"])",
    "--method=$(opt["method"])",
    "--write-csv=$(opt["write-csv"])",
    "--write-markdown=$(opt["write-markdown"])",
    "--quiet=true",
    "--print-summary=false",
  ]
end

function _run_import_audit(case_id, dtf_file, output_dir, opt)
  mkpath(output_dir)
  path = joinpath(output_dir, "dtf_import_audit.md")
  open(path, "w") do io
    println(io, "# DTF import audit\n")
    println(io, "- generated: ", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
    println(io, "- parser strict mode: ", opt["parser-strict"], "\n")
    ImportAudit.audit_case(io, case_id, dtf_file; strict = opt["parser-strict"])
  end
  return path
end

function _run_base_validation(dtf_file, for002_file, output_dir, opt)
  args = _common_validation_args(opt, dtf_file, for002_file, output_dir)
  push!(args, "--strict=$(opt["parser-strict"])")
  push!(args, "--legacy-voltage-level-collapse-230kv=$(opt["legacy-voltage-level-collapse-230kv"])")
  push!(args, "--transformer-ratio-mode=$(opt["transformer-ratio-mode"])")
  return NativeBaseValidation.run_validation(args; return_details = false)
end

function _run_outage_validation(dtf_file, for002_file, output_dir, opt)
  args = _common_validation_args(opt, dtf_file, for002_file, output_dir)
  push!(args, "--strict=false")
  return NativeOutageValidation.run_validation(args; return_details = false)
end

function _run_matpower_validation(dtf_file, for002_file, output_dir, opt)
  opt["write-matpower"] || throw(ArgumentError("MATPOWER roundtrip validation requires --write-matpower=true"))
  args = _common_validation_args(opt, dtf_file, for002_file, output_dir)
  append!(args, [
    "--write-matpower=true",
    "--run-outages=$(opt["run-outages"])",
    "--strict=false",
    "--details=true",
  ])
  return MatpowerRoundtripValidation.run_validation(args; return_details = true)
end

function _write_suite_markdown(path, opt, inventory, rows, modes)
  open(path, "w") do io
    println(io, "# Unified DTF validation suite\n")
    println(io, "- generated: ", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
    println(io, "- data directory: `", opt["data-dir"], "`")
    println(io, "- selected modes: `", join(modes, ", "), "`")
    println(io, "- selected cases: `", opt["case"], "`")
    println(io, "- method: `", opt["method"], "`; tolerance: `", opt["tol"], "`; max iterations: `", opt["max-iter"], "`")
    println(io, "- transformer ratio mode: `", opt["transformer-ratio-mode"], "`")
    println(io, "- parser strict mode: `", opt["parser-strict"], "`\n")

    println(io, "## Input inventory\n")
    println(io, "| Case | Suffix | FOR001 | FOR002 |")
    println(io, "|---|---:|---|---|")
    for item in inventory
      println(io, "| ", _md(item.case_id), " | ", _md(item.suffix), " | ", _md(item.dtf_file), " | ", _md(item.for002_file), " |")
    end

    println(io, "\n## Results\n")
    println(io, "| Case | Mode | Status | Converged | Scenarios | max |dV| kV | max |dVa| deg | max |dP| MW | max |dQ| MVar |")
    println(io, "|---|---|---|---:|---:|---:|---:|---:|---:|")
    for row in rows
      println(
        io,
        "| ", _md(row.case_id),
        " | ", _md(row.mode),
        " | ", _md(row.status),
        " | ", _md(row.converged),
        " | ", _md(row.scenario_count),
        " | ", _md(row.max_abs_d_vm_kV),
        " | ", _md(row.max_abs_d_va_deg),
        " | ", _md(row.max_abs_d_p_MW),
        " | ", _md(row.max_abs_d_q_MVar),
        " |",
      )
    end

    failures = [row for row in rows if !(row.status in ("ok", "ok_no_outages"))]
    if !isempty(failures)
      println(io, "\n## Diagnostics\n")
      for row in failures
        println(io, "- **", row.case_id, " / ", row.mode, " / ", row.status, "**: ", isempty(row.message) ? "no additional message" : row.message)
      end
    end
  end
  return path
end

function _handle_error!(rows, item, mode, mode_output, err, opt)
  message = sprint(showerror, err)
  push!(
    rows,
    _summary_row(
      case_id = item.case_id,
      suffix = item.suffix,
      mode = mode,
      status = "error",
      dtf_file = item.dtf_file,
      for002_file = item.for002_file,
      output_dir = mode_output,
      message = message,
    ),
  )
  opt["quiet"] || println(stderr, "[", item.case_id, "/", mode, "] ERROR: ", message)
  opt["continue-on-error"] || throw(err)
  return nothing
end

function run_suite(args = ARGS)
  opt = parse_suite_cli(args)
  opt === nothing && return nothing

  modes = _selected_modes(opt["mode"])
  selected_cases = _selected_cases(opt["case"])
  inventory = _discover_for_files(opt["data-dir"])
  isempty(inventory) && throw(ArgumentError("no FOR001*.DAT or FOR002*.DAT files found in $(opt["data-dir"])"))

  if selected_cases !== nothing
    inventory = [item for item in inventory if uppercase(item.case_id) in selected_cases]
    isempty(inventory) && throw(ArgumentError("none of the requested cases were found: $(join(sort!(collect(selected_cases)), ", "))"))
  end

  mkpath(opt["output-dir"])
  rows = NamedTuple[]

  for item in inventory
    case_output = joinpath(opt["output-dir"], _safe_case_dir(item.case_id))
    opt["quiet"] || println("\nCase ", item.case_id, ": ", item.dtf_file, " / ", item.for002_file)

    if item.dtf_file isa Missing
      for mode in modes
        push!(
          rows,
          _summary_row(
            case_id = item.case_id,
            suffix = item.suffix,
            mode = mode,
            status = "skipped_missing_for001",
            for002_file = item.for002_file,
            output_dir = joinpath(case_output, mode),
            message = "No matching FOR001 file was found.",
          ),
        )
      end
      continue
    end

    for mode in modes
      mode_output = joinpath(case_output, mode)
      opt["quiet"] || print("  ", rpad(mode, 10), " ... ")

      if mode in ("base", "outages") && item.for002_file isa Missing
        push!(
          rows,
          _summary_row(
            case_id = item.case_id,
            suffix = item.suffix,
            mode = mode,
            status = "skipped_missing_for002",
            dtf_file = item.dtf_file,
            output_dir = mode_output,
            message = "No matching FOR002 reference file was found.",
          ),
        )
        opt["quiet"] || println("SKIPPED (missing FOR002)")
        continue
      end

      try
        if mode == "audit"
          path = _run_import_audit(item.case_id, item.dtf_file, mode_output, opt)
          push!(
            rows,
            _summary_row(
              case_id = item.case_id,
              suffix = item.suffix,
              mode = mode,
              status = "ok",
              dtf_file = item.dtf_file,
              for002_file = item.for002_file,
              output_dir = mode_output,
              scenario_count = 1,
              message = path,
            ),
          )

        elseif mode == "base"
          result = _run_base_validation(item.dtf_file, item.for002_file, mode_output, opt)
          metrics = result.metrics
          status = result.converged ? "ok" : "not_converged"
          push!(
            rows,
            _summary_row(
              case_id = item.case_id,
              suffix = item.suffix,
              mode = mode,
              status = status,
              dtf_file = item.dtf_file,
              for002_file = item.for002_file,
              output_dir = mode_output,
              converged = result.converged,
              scenario_count = 1,
              max_abs_d_vm_kV = metrics.max_abs_d_vm_kV,
              max_abs_d_va_deg = metrics.max_abs_d_va_deg,
              max_abs_d_p_MW = result.max_branch_d_p_MW,
              max_abs_d_q_MVar = result.max_branch_d_q_MVar,
              message = "iterations=$(result.iterations), final_mismatch=$(result.final_mismatch)",
            ),
          )

        elseif mode == "outages"
          result = _run_outage_validation(item.dtf_file, item.for002_file, mode_output, opt)
          metrics = result.metrics_rows
          all_matched = all(row -> row.match_status == "matched", result.matching_rows)
          all_converged = all(row -> row.converged, metrics)
          status = isempty(metrics) ? "ok_no_outages" : (all_matched && all_converged ? "ok" : (!all_matched ? "outage_match_failed" : "not_converged"))
          push!(
            rows,
            _summary_row(
              case_id = item.case_id,
              suffix = item.suffix,
              mode = mode,
              status = status,
              dtf_file = item.dtf_file,
              for002_file = item.for002_file,
              output_dir = mode_output,
              converged = isempty(metrics) ? missing : all_converged,
              scenario_count = length(metrics),
              max_abs_d_vm_kV = _finite_max([m.max_abs_d_vm_kV for m in metrics]),
              max_abs_d_va_deg = _finite_max([m.max_abs_d_va_deg for m in metrics]),
              max_abs_d_p_MW = _finite_max([m.max_abs_branch_d_p_MW for m in metrics]),
              max_abs_d_q_MVar = _finite_max([m.max_abs_branch_d_q_MVar for m in metrics]),
              message = "parsed_outages=$(result.parsed_outages), parsed_for002_blocks=$(result.parsed_for002_outage_blocks)",
            ),
          )

        elseif mode == "matpower"
          for002_arg = item.for002_file isa Missing ? joinpath(opt["data-dir"], "FOR002.DAT") : item.for002_file
          result = _run_matpower_validation(item.dtf_file, for002_arg, mode_output, opt)
          metrics = result.metrics_rows
          all_converged = all(row -> row.native_converged && row.roundtrip_converged, metrics)
          all_status = all(row -> row.status_match, metrics)
          status = all_converged && all_status ? "ok" : (!all_converged ? "not_converged" : "status_mismatch")
          push!(
            rows,
            _summary_row(
              case_id = item.case_id,
              suffix = item.suffix,
              mode = mode,
              status = status,
              dtf_file = item.dtf_file,
              for002_file = item.for002_file,
              output_dir = mode_output,
              converged = all_converged,
              scenario_count = length(metrics),
              max_abs_d_vm_kV = _finite_max([m.max_abs_d_vm_kV for m in metrics]),
              max_abs_d_va_deg = _finite_max([m.max_abs_d_va_deg for m in metrics]),
              max_abs_d_p_MW = _finite_max([m.max_abs_branch_d_p_MW for m in metrics]),
              max_abs_d_q_MVar = _finite_max([m.max_abs_branch_d_q_MVar for m in metrics]),
              message = "roundtrip branch status match=$(all_status)",
            ),
          )
        end
        opt["quiet"] || println(last(rows).status)
      catch err
        _handle_error!(rows, item, mode, mode_output, err, opt)
      end
    end
  end

  written_files = String[]
  if opt["write-csv"]
    push!(written_files, _write_csv(joinpath(opt["output-dir"], "dtf_validation_suite_summary.csv"), rows))
    inventory_rows = [
      (
        case_id = item.case_id,
        suffix = item.suffix,
        dtf_file = item.dtf_file,
        for002_file = item.for002_file,
        has_for001 = !(item.dtf_file isa Missing),
        has_for002 = !(item.for002_file isa Missing),
      ) for item in inventory
    ]
    push!(written_files, _write_csv(joinpath(opt["output-dir"], "dtf_validation_suite_inventory.csv"), inventory_rows))
  end
  if opt["write-markdown"]
    push!(
      written_files,
      _write_suite_markdown(
        joinpath(opt["output-dir"], "dtf_validation_suite_summary.md"),
        opt,
        inventory,
        rows,
        modes,
      ),
    )
  end

  failed_rows = [row for row in rows if !(row.status in ("ok", "ok_no_outages"))]
  if !opt["quiet"]
    println("\nUnified DTF validation suite")
    println("Output directory: ", opt["output-dir"])
    println("Checks: ", length(rows))
    println("Successful: ", length(rows) - length(failed_rows))
    println("Non-successful: ", length(failed_rows))
    for path in written_files
      println("  - ", path)
    end
  end

  result = (
    output_dir = opt["output-dir"],
    inventory = inventory,
    summary_rows = rows,
    written_files = written_files,
    successful = isempty(failed_rows),
  )

  if opt["strict"] && !isempty(failed_rows)
    details = join(("$(row.case_id)/$(row.mode): $(row.status)" for row in failed_rows), "; ")
    error("DTF validation suite strict checks failed: $details")
  end

  return result
end

function _running_as_script()
  return !isempty(PROGRAM_FILE) && abspath(PROGRAM_FILE) == abspath(@__FILE__)
end

if _running_as_script()
  Base.invokelatest(run_suite, ARGS)
  nothing
end
