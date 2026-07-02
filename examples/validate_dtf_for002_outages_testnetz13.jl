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
  opt["strict"] && any(r -> r.match_status != "matched", matching_rows) && exit(1)
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

function _running_as_script()
  return true
  return !isempty(PROGRAM_FILE) && abspath(PROGRAM_FILE) == abspath(@__FILE__)
end

if get(ENV, "SPARLECTRA_FOR002_OUTAGE_VALIDATION_NO_MAIN", "0") != "1" && _running_as_script()
  Base.invokelatest(run_validation, ARGS; return_details = false)
  nothing
end
