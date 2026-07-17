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

# Internal DTF validation module: DTF -> existing MATPOWER export/import
# roundtrip validation. Extracted from validate_dtf_suite.jl; used by the
# suite runner and directly runnable as its own CLI entry point.

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
    "dtf-file" => joinpath(@__DIR__, "..", "..", "data", "DTF", "FOR001.DAT"),
    "for002-file" => joinpath(@__DIR__, "..", "..", "data", "DTF", "FOR002.DAT"),
    "output-dir" => joinpath(@__DIR__, "..", "_out", "dtf_matpower_export_testnetz13"),
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
    "tap-changer-model" => "ideal",
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
function _tap_changer_model_symbol(s::AbstractString)
  normalized = replace(lowercase(strip(s)), "-" => "_")
  normalized in ("ideal", "impedance_correction") || throw(ArgumentError("--tap-changer-model must be ideal or impedance_correction"))
  return Symbol(normalized)
end
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
  native = Sparlectra.DTFImporter.build_net(case; tap_changer_model = _tap_changer_model_symbol(opt["tap-changer-model"]))
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
    println(io, "- tap-changer model (native build): `", opt["tap-changer-model"], "`")
    println(io, "- metadata fields: `mpc.bus_name`, `mpc.branch_name`, `mpc.branch_kind`, and `mpc.for001_contingencies` are exported when available; solving does not require them.\n")
    println(io, "- roundtrip import option: `enable_pq_gen_controllers=false`, because native DTF PQ generators are fixed injections and enabling MATPOWER PQ generator controllers changes the solved roundtrip semantics.\n")
    println(io, "- TAP convention: line rows export TAP = 0.0; transformer rows export their explicit ratio. This keeps MATPOWER line rows from being re-imported as nominal-tap transformer rows.\n")
    for r in results
      m = r.metric
      println(io, "## ", m.scenario)
      println(io, "- exported MATPOWER: `", r.matpower_file, "`")
      println(io, "- metadata counts: bus_name=", r.metadata.bus_name_count, ", branch_name=", r.metadata.branch_name_count, ", branch_kind=", r.metadata.branch_kind_count, ", for001_contingencies=", r.metadata.for001_contingency_count)
      println(io, "- exported line TAP values: `", join(_fmt_num.(r.metadata.tap_line_values), ", "), "`")
      println(io, "- exported transformer TAP values: `", join(_fmt_num.(r.metadata.tap_transformer_values), ", "), "`")
      println(io, "- metadata complete: ", r.metadata.metadata_complete)
      if !r.metadata.metadata_complete
        println(io, "- metadata limitation: one or more established MATPOWER metadata fields were not preserved with the expected Testnetz13 counts")
      end
      println(io, "- native converged: ", m.native_converged)
      println(io, "- roundtrip converged: ", m.roundtrip_converged)
      println(io, "- max |dV| pu: ", _fmt_num(m.max_abs_d_vm_pu))
      println(io, "- max branch |dP| MW: ", _fmt_num(m.max_abs_branch_d_p_MW))
      println(io, "- max branch |dQ| MVar: ", _fmt_num(m.max_abs_branch_d_q_MVar))
      println(io, "- max native branch g_pu: ", _fmt_num(m.max_abs_branch_g_pu_native), " (count > 1e-12: ", m.count_nonzero_branch_g_pu_native, ")")
      if m.count_nonzero_branch_g_pu_native == 0
        println(io, "- branch shunt conductance: all native branch `g_pu` values are zero/negligible")
      else
        affected = [string(row.branch_index_native, ":", row.branch_name_native, "=", _fmt_num(row.g_pu_native)) for row in r.branch_rows if abs(row.g_pu_native) > 1e-12]
        println(io, "- branch shunt conductance limitation: standard MATPOWER branch data cannot preserve branch `g_pu` exactly; equivalent endpoint bus `Gs` is exported through MATPOWER bus shunts for Ybus preservation. Affected branches: ", join(affected, "; "))
      end
      println(io, "- bus shunt totals native/exported/roundtrip Gs: ", _fmt_num(m.native_total_bus_Gs), " / ", _fmt_num(m.exported_total_bus_Gs), " / ", _fmt_num(m.roundtrip_total_bus_Gs))
      println(io, "- bus shunt totals native/exported/roundtrip Bs: ", _fmt_num(m.native_total_bus_Bs), " / ", _fmt_num(m.exported_total_bus_Bs), " / ", _fmt_num(m.roundtrip_total_bus_Bs))
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

_running_as_cli_script() = !isempty(PROGRAM_FILE) && abspath(PROGRAM_FILE) == abspath(@__FILE__)

if _running_as_cli_script()
  Base.invokelatest(run_validation, ARGS)
end

end # module MatpowerRoundtripValidation
