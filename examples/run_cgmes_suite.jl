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

# file: examples/run_cgmes_suite.jl
# purpose: instructive CGMES import walkthrough on the ENTSO-E conformity
# test sets. NOT part of the regular test suite — run it interactively to
# see what the importer does and how to use it:
#
#   julia --project=. examples/run_cgmes_suite.jl
#
# The ENTSO-E test-set package (~22 MB) is fetched on demand into
# <repo>/data/CGMES (gitignored; see cgmes_fetch_testsets.jl). It shows,
# per test configuration:
#   1. summarizeCGMES  — diagnose the delivery before building anything
#   2. importCGMES     — build the Net (+ short-circuit data container)
#   3. runpf!          — solve
#   4. compareWithSV   — validate against the shipped SV reference solution

using Sparlectra
using Printf

# --- ensure the ENTSO-E test sets are cached --------------------------------

const CGMES_CACHE = get(ENV, "SPARLECTRA_CGMES_CACHE", joinpath(dirname(@__DIR__), "data", "CGMES"))

if !isdir(joinpath(CGMES_CACHE, "extracted"))
  println("ENTSO-E test sets not cached yet — fetching (~22 MB, one-time) ...")
  include(joinpath(@__DIR__, "experimental", "cgmes_fetch_testsets.jl"))
  main(String[])
end

const EXTRACTED = joinpath(CGMES_CACHE, "extracted")
mgdir(d) = joinpath(EXTRACTED, "MicroGrid", "BaseCase_BC", "CGMES_v2.4.15_MicroGridTestConfiguration_" * d * "_v2")
sgdir(d) = joinpath(EXTRACTED, "SmallGrid", "BusBranch", d)

# name => (paths, one-line description)
const CASES = [
  "MicroGrid BE" => ([mgdir("BC_BE"), mgdir("BD")], "Belgian half of the ENTSO-E MicroGrid: 3× 2W + 1× 3W transformer, ratio taps, one asymmetrical PST, boundary set with EquivalentInjections."),
  "MicroGrid NL" => ([mgdir("BC_NL"), mgdir("BD")], "Dutch half: three machines with local voltage RegulatingControls and a closed busbar coupler (retained switch → bus link)."),
  "MicroGrid Assembled" => ([mgdir("BC_Assembled"), mgdir("BD")], "Both halves in one model: the X-node EquivalentInjections of BOTH sides are skipped, real cross-border flows replace them."),
  "SmallGrid BusBranch" => ([sgdir("CGMES_v2.4.15_SmallGridTestConfiguration_BaseCase_Complete_v3.0.0"), sgdir("CGMES_v2.4.15_SmallGridTestConfiguration_Boundary_v3.0.0")], "118-bus bus-branch benchmark with 176 lines and a full SV reference state."),
]

# --- the walkthrough --------------------------------------------------------

for (name, (paths, description)) in CASES
  println("\n", "="^78)
  println(name, " — ", description)
  println("="^78)

  # 1. Diagnose first (works on broken/incomplete deliveries too).
  summary = summarizeCGMES(path = paths)
  show(stdout, MIME"text/plain"(), summary)

  # 2. Import. `importCGMES` returns the Net plus the merged store, the
  #    topology, the always-harvested short-circuit data (§7.7) and the
  #    importer's notices; `createNetFromCGMES` is the Net-only shorthand.
  res = importCGMES(path = paths, name = name)
  println("\nimporter notices:")
  for m in res.messages
    println("   ", m)
  end
  sc = res.shortcircuit
  @printf("short-circuit data harvested (read, not evaluated): %d lines, %d machines, %d ENI, %d transformer ends, %d equivalent injections\n",
    length(sc.ac_line_segments), length(sc.synchronous_machines), length(sc.external_network_injections), length(sc.transformer_ends), length(sc.equivalent_injections))
  @printf("net: %d buses, %d branches, %d prosumers, %d shunts — slack: %s\n",
    length(res.net.nodeVec), length(res.net.branchVec), length(res.net.prosumpsVec), length(res.net.shuntVec), res.slack_bus)

  # 3. Solve.
  ite, erg = runpf!(res.net, 30, 1e-8, 0)
  if erg != 0
    println("power flow did NOT converge (erg=", erg, ") — inspect the notices above")
    continue
  end
  println("power flow converged in ", ite, " iterations")

  # 4. Validate against the SV profile shipped with the data set — the
  #    importer's numeric acceptance oracle.
  cmp = compareWithSV(res)
  @printf("SV comparison over %d buses:  max|dVm| = %.2e pu   rms = %.2e   max|dVa| = %.4f°\n", cmp.n, cmp.max_dvm, cmp.rms_dvm, cmp.max_dva)
  worst = sort(cmp.rows; by = r -> -abs(r.dvm))[1:min(3, end)]
  for r in worst
    @printf("   worst: %-18s vm=%.5f  sv=%.5f  dva=%+.4f°\n", r.bus, r.vm_pu, r.sv_vm_pu, r.dva)
  end

  # 5. Flow-level comparison against SvPowerFlow (per terminal; injections,
  #    shunts and loads in these data sets — branch terminals in real
  #    exchanges). Deviations here separate model errors from data-level
  #    SSH↔SV inconsistencies.
  f = cmp.flows
  @printf("SvPowerFlow comparison over %d terminals:  max|dP| = %.3f MW   max|dQ| = %.3f MVAr\n", f.n, f.max_dp, f.max_dq)
  for r in sort(f.rows; by = x -> -max(abs(x.dp), abs(x.dq)))[1:min(3, end)]
    @printf("   worst: %-10s %-18s @%-16s dp=%+.3f MW  dq=%+.3f MVAr\n", r.kind, r.name, r.bus, r.dp, r.dq)
  end
end

# --- Stage 2: let the tap controllers find the operating point --------------
#
# `tap_control = true` starts from the SSH tap positions instead of copying
# the solved SvTapStep positions, and attaches Sparlectra outer-loop
# controllers for every tap changer whose CGMES control flags are set. The
# run then goes through `run_sparlectra` (control framework), not `runpf!`.

println("\n", "="^78)
println("Stage 2 — tap controllers from CGMES TapChangerControl (MicroGrid BE)")
println("="^78)

res = importCGMES(path = [mgdir("BC_BE"), mgdir("BD")], tap_control = true, name = "MicroGrid BE (controlled)")
for m in res.messages
  (startswith(m, "tap control") || occursin("disabled — target bus", m)) && println("   ", m)
end

cfg = SparlectraConfig(powerflow = PowerFlowConfig(max_iter = 30, tol = 1e-8), output = OutputConfig(logfile_results = :off), control = ControlConfig())
runres = run_sparlectra(net = res.net, config = cfg)
cres = latest_control_result(res.net)
@printf("control loop: %s after %d outer iterations (%d power-flow solves)\n", cres.status, cres.outer_iterations, cres.powerflow_solves)
for row in buildTapControllerReportRows(res.net)
  @printf("   %-22s %-22s status=%-10s tap=%.5f  phase=%+.3f°\n", row.controller_name, row.mode, row.status, row.tap_ratio, row.phase_shift_deg)
end
println("""
   Note: the controllers do NOT reproduce the SvTapStep positions exactly, and
   they are not expected to — this data set's CGMES target deadbands are wide
   (PST 35 MW, OLTC 0.5 kV on 10.5 kV), so a controller legitimately stops one
   step short of the reference tool's position. What matters is that each
   target is met inside its own deadband.""")

println("\nDone. Next steps to explore:")
println("  - summarizeCGMES(path = <your delivery>) on any CGMES 2.4.15 folder/ZIP")
println("  - res.shortcircuit — typed source data for future IEC 60909 work (issue #277)")
println("  - task_plan_cgmes_import.md — staging plan (node-breaker, controllers, 3.0 …)")
