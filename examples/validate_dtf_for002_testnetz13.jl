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

using Sparlectra
using Printf

include(joinpath(@__DIR__, "dtf_for002_validation_utils.jl"))

# Developer notes:
# - Validates the native Testnetz13 DTF/FOR001 base case against FOR002.
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
  dtf_default = isfile(joinpath(@__DIR__, "..", "test", "fixtures", "dtf", "FOR001.DAT")) ? joinpath(@__DIR__, "..", "test", "fixtures", "dtf", "FOR001.DAT") : joinpath(@__DIR__, "FOR001.DAT")
  opt = Dict{String,Any}(
    "dtf-file" => normpath(dtf_default),
    "for002-file" => joinpath(@__DIR__, "FOR002.DAT"),
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
      println(io, "State residuals force FOR002 printed voltage magnitudes/angles into the native Sparlectra Ybus. The resulting bus injections are compared with the FOR002 printed bus table. This is more sensitive than solved branch-flow comparisons because FOR002 values may be rounded and transformer-adjacent nodes react strongly to small voltage/angle differences. These residuals are diagnostic, not hard pass/fail criteria yet; branch-flow deviations and solved generator/slack comparisons are currently stronger validation signals.\n")
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
  opt["strict"] && !converged && exit(1)
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

if get(ENV, "SPARLECTRA_FOR002_VALIDATION_NO_MAIN", "0") != "1"
  Base.invokelatest(run_validation, ARGS; return_details = false)
  nothing
end
