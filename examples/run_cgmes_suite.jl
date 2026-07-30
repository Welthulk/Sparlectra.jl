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

# --- full sweep: compute and evaluate EVERY cached/fetchable test set -------
#
# One row per test configuration: import → runpf! (islands enabled, default
# Q-limit handling) → compareWithSV. This is the "does the case solve, and how
# far is it from its own SV state" overview across the whole collection —
# conformity sets, RealGrid, and the ReliCapGrid/Svedala CGMES 3.0 family.
# Import errors and non-convergence become table rows, not aborts.

println("\n", "="^78)
println("Full sweep — every test set: import, solve, SV deviation")
println("="^78)

const CASES_DIR = joinpath(CGMES_CACHE, "cases")
mkpath(CASES_DIR)
_alias(a) = () -> Sparlectra.CGMESImporter.fetchCGMESTestSet(a; outdir = CASES_DIR)
_dirs(parts...) = () -> [joinpath(EXTRACTED, p...) for p in parts]

# label => path provider. Aliases resolve offline out of the cache; the
# non-alias conformity sets are addressed by their extracted directories.
# MicroGrid Type1–Type5 variants exist as well; the sweep covers each distinct
# grid once (the variants alter tap/regulation details of the same grid).
const SWEEP = [
  "MicroGrid BE" => _alias("microgrid_be"),
  "MicroGrid NL" => _alias("microgrid_nl"),
  "MicroGrid Assembled" => _alias("microgrid_assembled"),
  "SmallGrid BusBranch" => _alias("smallgrid"),
  "SmallGrid NodeBreaker" => _alias("smallgrid_nb"),
  "MiniGrid BusBranch" => _dirs(("MiniGrid", "BusBranch", "CGMES_v2.4.15_MiniGridTestConfiguration_BaseCase_v3"), ("MiniGrid", "BusBranch", "CGMES_v2.4.15_MiniGridTestConfiguration_Boundary_v3")),
  "MiniGrid NodeBreaker" => _dirs(("MiniGrid", "NodeBreaker", "CGMES_v2.4.15_MiniGridTestConfiguration_BaseCase_Complete_v3"), ("MiniGrid", "NodeBreaker", "CGMES_v2.4.15_MiniGridTestConfiguration_Boundary_v3")),
  "FullGrid BusBranch" => _alias("fullgrid"),
  "FullGrid NodeBreaker" => _alias("fullgrid_nb"),
  "RealGrid" => _alias("realgrid"),
  "PST PTChLin PTE1" => _dirs(("PST_PTChLin_PTE1_PSEI",)),
  "PST PTChLin PTE2" => _dirs(("PST_PTChLin_PTE2_PSEI",)),
  "PST PTChTab PTE2" => _dirs(("PST_PTChTab_PTE2_PSEI",)),
  "TransformerLineTest" => _dirs(("TransformerLineTest",)),
  "ExplicitLoadFlow" => _dirs(("ENTSOE_ExplicitLoadFlowCalculation", "ENTSOE_CGMES_v2.4_ExplicitLoadFlowCalculation"),),
  "ReliCap Svedala" => _alias("svedala"),
  "ReliCap Svedala+Neighbours" => _alias("svedala_neighbours"),
  "ReliCap Belgovia" => _alias("belgovia"),
  "ReliCap Britheim" => _alias("britheim"),
  "ReliCap Espheim" => _alias("espheim"),
  "ReliCap Galia" => _alias("galia"),
  "ReliCap Nordheim" => _alias("nordheim"),
  "ReliCap Portheim" => _alias("portheim"),
  "ReliCap CGM (all areas)" => _alias("relicapgrid_cgm"),
]

# first ~100 chars of the error with newlines collapsed — enough to identify
# the failure ("no slack bus", "island 1 ... failed: max_iter reached", …)
_err_snippet(err) = replace(sprint(showerror, err), r"\s+" => " ")[1:min(100, end)]
_nan_row(label, status, note; version = "-", buses = 0, iso = 0, it = 0, warnings = 0) =
  (label = label, version = version, buses = buses, iso = iso, status = status, it = it,
    max_dvm = NaN, rms_dvm = NaN, max_dva = NaN, fn = 0, rms_dp = NaN, max_dp = NaN, warnings = warnings, note = note)

sweep_rows = NamedTuple[]
for (label, provider) in SWEEP
  # import and solve fail independently and honestly: an unresolvable delivery
  # is IMPORT ERR, a solver exception (e.g. "no slack bus" on an all-isolated
  # single-converter area) is SOLVE ERR, a clean erg != 0 is NO CONV.
  # (`imp`, not `res`: the walkthrough above owns a global `res`, and a bare
  # assignment here would hit Julia's soft-scope ambiguity warning.)
  imp = try
    importCGMES(path = provider(), name = label)
  catch err
    push!(sweep_rows, _nan_row(label, "IMPORT ERR", _err_snippet(err)))
    @printf("%-28s %s\n", label, "IMPORT ERR")
    continue
  end
  net = imp.net
  warnings = count(m -> startswith(m, "warning:"), imp.messages)
  # Solve from the delivery's own SV/SSH start first; when that fails, retry
  # once from a flat start on a FRESH import (the failed attempt leaves the
  # net in its diverged state). FullGrid is the case this exists for: its SV
  # profile is internally inconsistent (a 14.5° angle jump across a 0.3 Ω
  # line, ≈400 pu start residual), while the network itself — placeholder
  # shunt/tap guarded away — solves cleanly from flat.
  _solve_stats(r, ite, tag) = begin
    cmp = compareWithSV(r)
    (label = label, version = r.store.version, buses = length(r.net.nodeVec), iso = length(r.net.isoNodes), status = tag, it = ite,
      max_dvm = cmp.max_dvm, rms_dvm = cmp.rms_dvm, max_dva = cmp.max_dva,
      fn = cmp.flows.n, rms_dp = cmp.flows.rms_dp, max_dp = cmp.flows.max_dp, warnings = warnings, note = "")
  end
  row = try
    ite, erg = runpf!(net; config = PowerFlowConfig(max_iter = 60, tol = 1e-8), verbose = 0)
    erg == 0 ? _solve_stats(imp, ite, "ok") : nothing
  catch err
    _nan_row(label, "SOLVE ERR", _err_snippet(err); version = imp.store.version, buses = length(net.nodeVec), iso = length(net.isoNodes), warnings = warnings)
  end
  if row === nothing
    row = try
      imp2 = importCGMES(path = provider(), name = label)
      ite2, erg2 = runpf!(imp2.net, 60, 1e-8, 0; opt_flatstart = true)
      erg2 == 0 ? _solve_stats(imp2, ite2, "ok (flat)") :
        _nan_row(label, "NO CONV", "SV and flat start both fail"; version = imp2.store.version, buses = length(imp2.net.nodeVec), iso = length(imp2.net.isoNodes), it = ite2, warnings = warnings)
    catch err
      _nan_row(label, "SOLVE ERR", _err_snippet(err); version = imp.store.version, buses = length(net.nodeVec), iso = length(net.isoNodes), warnings = warnings)
    end
  end
  push!(sweep_rows, row)
  @printf("%-28s %s\n", row.label, row.status)
end

println()
@printf("%-28s %-7s %6s %4s  %-10s %3s  %9s %9s %8s  %6s %9s %8s  %4s  %s\n",
  "case", "version", "buses", "iso", "status", "it", "max|dvm|", "rms dvm", "max|dva|", "n_flow", "rms dp", "max dp", "warn", "note")
println("-"^150)
for r in sweep_rows
  @printf("%-28s %-7s %6d %4d  %-10s %3d  %9s %9s %8s  %6d %9s %8s  %4d  %s\n",
    r.label, r.version, r.buses, r.iso, r.status, r.it,
    isnan(r.max_dvm) ? "-" : @sprintf("%.5f", r.max_dvm), isnan(r.rms_dvm) ? "-" : @sprintf("%.5f", r.rms_dvm),
    isnan(r.max_dva) ? "-" : @sprintf("%.2f°", r.max_dva),
    r.fn, isnan(r.rms_dp) ? "-" : @sprintf("%.3f MW", r.rms_dp), isnan(r.max_dp) ? "-" : @sprintf("%.1f MW", r.max_dp), r.warnings, r.note)
end
println("""

Reading the table:
  - max|dvm|/rms dvm/max|dva| compare the solved state against the delivery's
    own SV profile (excluded: isolated/de-energized buses).
  - "ok (flat)" = the SV-based start failed but a flat start solves; the SV
    columns then measure the SV profile's own inconsistency, not solver error
    (FullGrid: a 14.5° angle jump across a 0.3 Ω line is IN the shipped SV).
  - FullGrid is the import/export COMPLETENESS set (every CGMES class at
    least once, X.99 placeholder values); the importer's plausibility guards
    warn about and skip the placeholder shunt/tap-row, after which the
    network itself solves.
  - warn counts "warning:" importer messages (substituted values).""")

println("\nDone")
