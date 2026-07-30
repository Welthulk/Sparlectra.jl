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

# file: examples/experimental/val_realgrid_remeasure.jl
# purpose: RealGrid-Neumessung nach den vier Hebeln, die der Forensik vom
# 2026-07-27 fehlten: P2 Tabular-PST (SV-Start-Mismatch 197 → 5 pu), #192
# distributed slack (Referenz wurde mit Slack-Verteilung gelöst), Q-Limit-
# Fixes (Referenz wurde MIT Limits gelöst), P3 machine_control (die 10 in der
# RealGrid-Doku genannten RegulatingControls haben Off-Nominal-TARGETS; per
# Store-Messung hat keine Maschine ein FERNES RC-Terminal — der Lauf zeigt,
# wie viele Controller real ankommen). Alte Baseline: max Δvm 1.01,
# 241 Busse mit |Δvm| > 0.01, Kollaps-Busse im 63-kV-Level.
#
#   julia --project=. examples/experimental/val_realgrid_remeasure.jl

using Sparlectra
using Printf
using Statistics

const CASEDIR = joinpath(dirname(dirname(@__DIR__)), "data", "CGMES", "cases")
mkpath(CASEDIR)
const COLLAPSE_PROBE = "1_38470111"   # schlimmster Kollaps-Bus der Forensik

function _import_rg(; machine_control::Bool = false)
  zip = Sparlectra.CGMESImporter.fetchCGMESTestSet("realgrid"; outdir = CASEDIR)
  # Referenz laut RealGrid-Doku: feste Taps → tap_control = false
  return importCGMES(path = zip, tap_control = false, machine_control = machine_control)
end

function _voltage_stats(cmp)
  over = count(r -> abs(r.dvm) > 0.01, cmp.rows)
  @printf("  voltages : rms dvm %.4f  max dvm %.4f  |dvm|>0.01: %d von %d   [alt: max 1.01 / 241 Busse]\n",
    cmp.rms_dvm, cmp.max_dvm, over, cmp.n)
  @printf("             rms dva %.2f°  max dva %.2f°\n", cmp.rms_dva, cmp.max_dva)
  worst = sort(collect(cmp.rows); by = r -> -abs(r.dvm))
  println("  worst 5 (vm | sv_vm):")
  for r in first(worst, min(5, length(worst)))
    @printf("      %-24s %.4f | %.4f  (Δ %+.4f)\n", r.bus, r.vm_pu, r.sv_vm_pu, r.dvm)
  end
  probe = findfirst(r -> r.bus == COLLAPSE_PROBE, cmp.rows)
  if probe !== nothing
    r = cmp.rows[probe]
    @printf("  probe    : %s  vm %.4f (sv %.4f) → %s\n", COLLAPSE_PROBE, r.vm_pu, r.sv_vm_pu,
      r.vm_pu < 0.5 ? "KOLLABIERT" : "ok")
  else
    println("  probe    : ", COLLAPSE_PROBE, " nicht in den Vergleichszeilen (isoliert/deenergisiert?)")
  end
end

function _run(label; machine_control = false, qlimits = false, dslack = false)
  println("\n", "="^76, "\n", label, "\n", "="^76)
  result = _import_rg(machine_control = machine_control)
  net = result.net

  pf = PowerFlowConfig(
    max_iter = 60,
    tol = 1e-6,
    islands_enabled = true,
    qlimits = Sparlectra.QLimitConfig(ignore_q_limits = !qlimits),
    distributed_slack = dslack ? Sparlectra.DistributedSlackConfig(enabled = true, p_mode = :pg_weighted, fallback = :ref_only) : Sparlectra.DistributedSlackConfig(),
  )

  ctrls = collect_outer_controllers(net)
  machine_control && println("  machine_control: ", count(c -> c isa MachineVoltageControl, ctrls),
    " Controller angekommen, ", count(m -> occursin("not attachable", m) || occursin("skipped", m), result.messages), " Fallbacks")
  local ite, erg
  if isempty(ctrls)
    ite, erg = runpf!(net; config = pf, verbose = 0)
  else
    cres = run_control!(net; controllers = ctrls, pf_config = pf, control_config = ControlConfig(max_outer_iterations = 20), verbose = 0)
    ite = cres.last_pf_iterations
    erg = cres.status == :pf_failed ? 1 : 0
    println("  control  : status=", cres.status, " outer=", cres.outer_iterations, " pf_solves=", cres.powerflow_solves)
  end
  println("  converged: ", erg == 0, "  (", ite, " it)")
  erg == 0 || return nothing

  cmp = compareWithSV(result)
  _voltage_stats(cmp)

  f = cmp.flows
  f.n > 0 && @printf("  flows    : rms dp %.3f MW  max dp %.1f MW  (n=%d)\n", f.rms_dp, f.max_dp, f.n)

  try
    st = Sparlectra.rectangular_pf_status(net)
    if hasproperty(st, :distributed_slack_active) && st.distributed_slack_active
      @printf("  λ_P      : %+.1f MW über %d Teilnehmer\n",
        st.distributed_slack_lambda_p_mw, st.distributed_slack_participants)
    end
  catch
  end
  return nothing
end

function main()
  # 1) Baseline neu unter aktuellem Code — isoliert den P2-Effekt gegen die alte Forensik
  _run("1) klassisch, Q-Limits AUS — neue Baseline (P2-Effekt)")

  # 2) + Q-Limits (Referenz wurde MIT Limits gelöst; Active-Set jetzt mit Post-Convergence-Round)
  _run("2) + Q-Limits AN"; qlimits = true)

  # 3) + verteilter Slack (Referenz: Zwischenslack + Verteilung über alle Busse)
  _run("3) + distributed slack pg_weighted"; qlimits = true, dslack = true)

  # 4) + Fernregelung (P3; die vset-Plausibilitätsbande greift weiter — das
  #    Skript meldet, wie viele Controller real ankommen)
  _run("4) + machine_control (P3)"; machine_control = true, qlimits = true, dslack = true)
end

Base.invokelatest(main)
