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

using Printf
using Sparlectra

ENV["SPARLECTRA_FOR002_VALIDATION_NO_MAIN"] = "1"
include(joinpath(@__DIR__, "validate_dtf_for002_testnetz13.jl"))

const CASES = [
  (id = "A", suffix = "", description = "base case"),
  (id = "B", suffix = "B", description = "NORD PV bus at 230.5 kV"),
  (id = "C", suffix = "C", description = "DELTA fixed tap step 7"),
  (id = "D", suffix = "D", description = "slack moved to ALPHA"),
  (id = "E", suffix = "E", description = "DELTA tap step 7 plus 60 degree skew-angle regulator"),
]

missing_to_float(x) = x === missing || x === nothing ? NaN : Float64(x)
fmt(x; digits = 3) = isfinite(Float64(x)) ? string(round(Float64(x); digits)) : "n/a"

function firstrow(rows, pred)
  idx = findfirst(pred, rows)
  return idx === nothing ? nothing : rows[idx]
end

function top_rows(rows, field; n = 10)
  present = [r for r in rows if !(getproperty(r, field) isa Missing) && isfinite(Float64(getproperty(r, field)))]
  return first(sort(present; by = r -> -abs(Float64(getproperty(r, field)))), min(n, length(present)))
end

function transformer_loss_rows(details)
  net = details.net
  rows = NamedTuple[]
  for (idx, br) in enumerate(net.branchVec)
    meta = get(net.matpower_branch_metadata, idx, nothing)
    meta === nothing && continue
    hasproperty(meta, :dtf_kind) && uppercase(string(meta.dtf_kind)) == "T" || continue
    pf = br.fBranchFlow === nothing || br.fBranchFlow.pFlow === nothing ? 0.0 : Float64(br.fBranchFlow.pFlow)
    pt = br.tBranchFlow === nothing || br.tBranchFlow.pFlow === nothing ? 0.0 : Float64(br.tBranchFlow.pFlow)
    qf = br.fBranchFlow === nothing || br.fBranchFlow.qFlow === nothing ? 0.0 : Float64(br.fBranchFlow.qFlow)
    qt = br.tBranchFlow === nothing || br.tBranchFlow.qFlow === nothing ? 0.0 : Float64(br.tBranchFlow.qFlow)
    g_pu = hasproperty(meta, :active_no_load_g_pu) ? Float64(meta.active_no_load_g_pu) : 0.0
    gs_each = g_pu * net.baseMVA / 2
    vf2 = net.nodeVec[br.fromBus]._vm_pu^2
    vt2 = net.nodeVec[br.toBus]._vm_pu^2
    no_load = gs_each * (vf2 + vt2)
    push!(rows, (
      branch_index = idx,
      branch_label = meta.orig_name,
      from_bus = meta.orig_name,
      transformer_series_p_loss_MW = pf + pt,
      transformer_no_load_g_p_loss_MW = no_load,
      transformer_total_p_loss_MW = pf + pt + no_load,
      transformer_series_q_loss_MVar = qf + qt,
      transformer_shunt_b_q_loss_MVar = missing,
      transformer_total_q_loss_MVar = missing,
      g_pu = g_pu,
      gs_each_MW_at_1pu = gs_each,
      tap_ratio = hasproperty(meta, :tap_ratio) ? meta.tap_ratio : missing,
      actual_tap_step = hasproperty(meta, :actual_tap_step) ? something(meta.actual_tap_step, missing) : missing,
      phase_shift_deg = hasproperty(meta, :phase_shift_deg) ? meta.phase_shift_deg : missing,
      dtf_control_model = hasproperty(meta, :dtf_control_model) ? meta.dtf_control_model : missing,
      tap_fraction = hasproperty(meta, :tap_fraction) ? meta.tap_fraction : missing,
      skew_angle_deg = hasproperty(meta, :skew_angle_deg) ? meta.skew_angle_deg : missing,
      effective_complex_tap = hasproperty(meta, :effective_complex_tap) ? meta.effective_complex_tap : missing,
      dtf_tap_convention = hasproperty(meta, :dtf_tap_convention) ? meta.dtf_tap_convention : missing,
    ))
  end
  return rows
end

function write_named_csv(path, rows)
  write_csv(path, rows)
  return path
end

function case_note(case_id, details, summary, loss_rows)
  if case_id == "A"
    rows400 = [r for r in details.buses if r.vn_kV >= 350]
    mean400 = isempty(rows400) ? NaN : sum(r.d_vm_kV for r in rows400) / length(rows400)
    return abs(mean400) < 1.0 ? "400-kV mean voltage bias is below 1 kV; old ~2 kV low bias is not visible as a uniform offset." : "400-kV mean voltage bias remains visible at $(fmt(mean400)) kV."
  elseif case_id == "B"
    nord = firstrow(details.generator_rows, r -> _norm_name(r.bus_name) == _norm_name("NORD S1"))
    bus = firstrow(details.buses, r -> _norm_name(r.bus_name) == _norm_name("NORD S1"))
    return nord === nothing || bus === nothing ? "NORD row missing." : "NORD effective type=$(nord.model_bus_type_effective), V=$(fmt(bus.model_vm_kV)) kV, Q=$(fmt(nord.model_qg_result_MVar)) MVar."
  elseif case_id == "C"
    vals = String[]
    for name in ["ALPHA S1", "BETA1 S1", "DELTA1S1", "DELTA2S1"]
      row = firstrow(details.buses, r -> _norm_name(r.bus_name) == _norm_name(name))
      row !== nothing && push!(vals, "$(row.bus_name) dV=$(fmt(row.d_vm_kV)) kV")
    end
    return join(vals, "; ")
  elseif case_id == "D"
    slack = firstrow(details.generator_rows, r -> r.is_slack)
    bstadt = firstrow(details.generator_rows, r -> _norm_name(r.bus_name) == _norm_name("BSTADTS1"))
    return "slack=$(slack === nothing ? "missing" : slack.bus_name); BSTADTS1 P/Q=$(bstadt === nothing ? "missing" : fmt(bstadt.model_pg_result_MW) * "/" * fmt(bstadt.model_qg_result_MVar))."
  elseif case_id == "E"
    beta = [r for r in details.branches if occursin("BETA1", r.branch_label) && occursin("BETA2", r.branch_label)]
    delta = [r for r in details.branches if occursin("DELTA1", r.branch_label) && occursin("DELTA2", r.branch_label)]
    bp = isempty(beta) ? "missing" : join(fmt.(getproperty.(beta, :model_p_from_MW)), ", ")
    dp = isempty(delta) ? "missing" : join(fmt.(getproperty.(delta, :model_p_from_MW)), ", ")
    return "BETA transformer P(from) MW=$(bp); DELTA transformer P(from) MW=$(dp)."
  end
  return ""
end

function main()
  repo = normpath(joinpath(@__DIR__, ".."))
  data_dir = joinpath(repo, "data", "DTF")
  outdir = joinpath(@__DIR__, "_out", "schaefer_dtf_validation")
  mkpath(outdir)
  missing = String[]
  for c in CASES
    for stem in ("FOR001", "FOR002")
      path = joinpath(data_dir, stem * c.suffix * ".DAT")
      isfile(path) || push!(missing, path)
    end
  end
  if !isempty(missing)
    error("Missing DTF/FOR002 files: " * join(missing, ", "))
  end

  summary_rows = NamedTuple[]
  all_loss_rows = NamedTuple[]
  md_path = joinpath(outdir, "schaefer_dtf_validation_summary.md")
  open(md_path, "w") do io
    println(io, "# Schäfer DTF/FOR002 validation summary\n")
    println(io, "Generated with native `DTFImporter.read_dtf` -> `DTFImporter.build_net` -> `runpf!`; this is a diagnostic report, not an electrical-model fix.\n")
    for c in CASES
      case_out = joinpath(outdir, "case_$(c.id)")
      dtf = joinpath(data_dir, "FOR001$(c.suffix).DAT")
      for002 = joinpath(data_dir, "FOR002$(c.suffix).DAT")
      details = run_validation(["--dtf-file=$(dtf)", "--for002-file=$(for002)", "--output-dir=$(case_out)", "--write-csv=true", "--write-markdown=true", "--quiet=true"]; return_details = true)
      Sparlectra.calcNetLosses!(details.net)
      ref = parse_for002_ground_load_flow(for002)
      p_loss, q_loss = Sparlectra.getTotalLosses(net = details.net)
      slack = firstrow(details.generator_rows, r -> r.is_slack)
      loss_rows = transformer_loss_rows(details)
      append!(all_loss_rows, [(case_id = c.id, r...) for r in loss_rows])
      series_p = sum(r.transformer_series_p_loss_MW for r in loss_rows)
      no_load_p = sum(r.transformer_no_load_g_p_loss_MW for r in loss_rows)
      total_p = sum(r.transformer_total_p_loss_MW for r in loss_rows)
      series_q = sum(r.transformer_series_q_loss_MVar for r in loss_rows)
      note = case_note(c.id, details, nothing, loss_rows)
      row = (
        case_id = c.id,
        for001_file = relpath(dtf, repo),
        for002_file = relpath(for002, repo),
        converged = details.converged,
        iterations = details.iterations,
        final_mismatch = details.final_mismatch,
        slack_bus = slack === nothing ? "" : slack.bus_name,
        slack_p_sparlectra_MW = slack === nothing ? NaN : slack.model_pg_result_MW,
        slack_q_sparlectra_MVar = slack === nothing ? NaN : slack.model_qg_result_MVar,
        slack_p_for002_MW = slack === nothing ? NaN : slack.for002_pg_MW,
        slack_q_for002_MVar = slack === nothing ? NaN : slack.for002_qg_MVar,
        slack_dp_MW = slack === nothing ? NaN : slack.d_pg_result_vs_for002_MW,
        slack_dq_MVar = slack === nothing ? NaN : slack.d_qg_result_vs_for002_MVar,
        total_loss_p_sparlectra_MW = p_loss,
        total_loss_q_sparlectra_MVar = q_loss,
        total_loss_p_for002_MW = something(ref.total_p_loss_MW, NaN),
        total_loss_q_for002_MVar = something(ref.total_q_loss_MVar, NaN),
        total_loss_dp_MW = p_loss - something(ref.total_p_loss_MW, NaN),
        total_loss_dq_MVar = q_loss - something(ref.total_q_loss_MVar, NaN),
        max_abs_dV_kV = details.metrics.max_abs_d_vm_kV,
        max_abs_dV_pu = details.metrics.max_abs_d_vm_pu,
        max_abs_dVa_deg = details.metrics.max_abs_d_va_deg,
        max_abs_branch_dP_MW = details.metrics.max_abs_branch_d_p_MW,
        max_abs_branch_dQ_MVar = details.metrics.max_abs_branch_d_q_MVar,
        max_abs_bus_dP_MW = details.metrics.max_abs_state_residual_p_MW,
        max_abs_bus_dQ_MVar = details.metrics.max_abs_state_residual_q_MVar,
        transformer_series_p_loss_MW = series_p,
        transformer_no_load_g_p_loss_MW = no_load_p,
        transformer_total_p_loss_MW = total_p,
        transformer_series_q_loss_MVar = series_q,
        transformer_shunt_b_q_loss_MVar = missing,
        transformer_total_q_loss_MVar = missing,
        notes = note,
      )
      push!(summary_rows, row)
      write_named_csv(joinpath(outdir, "schaefer_dtf_case_$(c.id).csv"), [row])

      println(io, "## Case $(c.id): $(c.description)\n")
      println(io, "- converged: `$(details.converged)`, iterations: `$(details.iterations)`, final mismatch: `$(details.final_mismatch)`")
      println(io, "- slack: `$(row.slack_bus)` Sparlectra P/Q=$(fmt(row.slack_p_sparlectra_MW))/$(fmt(row.slack_q_sparlectra_MVar)); FOR002 P/Q=$(fmt(row.slack_p_for002_MW))/$(fmt(row.slack_q_for002_MVar)); delta=$(fmt(row.slack_dp_MW))/$(fmt(row.slack_dq_MVar))")
      println(io, "- losses Sparlectra P/Q=$(fmt(row.total_loss_p_sparlectra_MW))/$(fmt(row.total_loss_q_sparlectra_MVar)); FOR002 P/Q=$(fmt(row.total_loss_p_for002_MW))/$(fmt(row.total_loss_q_for002_MVar)); delta=$(fmt(row.total_loss_dp_MW))/$(fmt(row.total_loss_dq_MVar))")
      println(io, "- max deviations: |dV|=$(fmt(row.max_abs_dV_kV)) kV ($(fmt(row.max_abs_dV_pu)) pu), |dVa|=$(fmt(row.max_abs_dVa_deg)) deg, branch |dP|=$(fmt(row.max_abs_branch_dP_MW)) MW, branch |dQ|=$(fmt(row.max_abs_branch_dQ_MVar)) MVar")
      println(io, "- transformer active losses: series=$(fmt(series_p)) MW, no-load G/shunt=$(fmt(no_load_p)) MW, total=$(fmt(total_p)) MW; total active network loss=$(fmt(p_loss)) MW")
      println(io, "- note: $(note)\n")
      println(io, "### Top voltage deviations\n")
      for r in top_rows(details.buses, :d_vm_kV)
        println(io, "- $(r.bus_name): model $(fmt(r.model_vm_kV)) kV vs FOR002 $(fmt(r.for002_vm_kV)) kV; dV=$(fmt(r.d_vm_kV)) kV; dVa=$(fmt(r.d_va_deg)) deg")
      end
      println(io, "\n### Top branch P deviations\n")
      for r in top_rows(details.branches, :d_p_from_MW)
        println(io, "- $(r.branch_label): model $(fmt(missing_to_float(r.model_p_from_MW))) MW vs FOR002 $(fmt(missing_to_float(r.for002_p_from_MW))) MW; dP=$(fmt(missing_to_float(r.d_p_from_MW))) MW")
      end
      println(io, "\n### Top branch Q deviations\n")
      for r in top_rows(details.branches, :d_q_from_MVar)
        println(io, "- $(r.branch_label): model $(fmt(missing_to_float(r.model_q_from_MVar))) MVar vs FOR002 $(fmt(missing_to_float(r.for002_q_from_MVar))) MVar; dQ=$(fmt(missing_to_float(r.d_q_from_MVar))) MVar")
      end
      println(io)
    end
  end
  write_named_csv(joinpath(outdir, "schaefer_dtf_validation_summary.csv"), summary_rows)
  write_named_csv(joinpath(outdir, "transformer_loss_summary.csv"), all_loss_rows)
  println("Wrote ", md_path)
  println("Wrote ", joinpath(outdir, "schaefer_dtf_validation_summary.csv"))
  return nothing
end

Base.invokelatest(main)
