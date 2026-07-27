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
#
# Example: PMU voltage-angle measurements in WLS state estimation
#
# Didactic storyline (four scenarios on the same case):
#   A) SCADA-only measurement set (noisy Vm / injections / flows) — the
#      baseline angle accuracy without any direct angle information.
#   B) SCADA + PMU phasors (Vm + Va, very small sigma) whose time base
#      coincides with the slack reference — angles snap to the truth.
#   C) The same PMUs, but their GPS time base is shifted by +2° against the
#      slack reference. With pmu_ref_offset = :auto the estimator carries a
#      common reference-angle offset α as an ADDITIONAL STATE and recovers
#      both the network angles and α.
#   D) The same shifted PMUs with pmu_ref_offset = :off — the offset stays
#      unmodeled, the objective explodes and the estimate is biased. This is
#      the failure mode the offset state exists to prevent.

# Date: 2026-07-27
# file: examples/state_estimation/state_estimation_pmu_angles.jl
# purpose: shows how PMU voltage-angle measurements enter the WLS estimator via the reference-angle offset state α (pmu_ref_offset), including the failure mode when the offset is not modeled

using Sparlectra
using Printf
using Dates

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))
using Random

const OUTDIR = joinpath(@__DIR__, "..", "_out")
const CASEFILE = "case9.m"

# PMU buses: a small subset is enough — each PMU pins one bus angle, the
# offset state α costs exactly one of those pins.
const PMU_BUS_IDXS = [1, 4, 7]

# The "true" misalignment between the PMU GPS time base and the slack
# reference in scenario C/D.
const PMU_REF_SHIFT_DEG = 2.0

function _angle_rmse_deg(Vref::Vector{ComplexF64}, Vest::Vector{ComplexF64})
  s = 0.0
  for i in eachindex(Vref)
    s += (rad2deg(angle(Vref[i]) - angle(Vest[i])))^2
  end
  return sqrt(s / length(Vref))
end

function _scada_measurements(net, rng)
  # Deliberately modest SCADA basis (no branch flows, coarse power accuracy):
  # observable, but with soft angle information — so the PMU contribution in
  # scenarios B/C stands out.
  std = measurementStdDevs(vm = 0.01, pinj = 3.0, qinj = 3.0)
  return generateMeasurementsFromPF(net; includeVm = true, includePinj = true, includeQinj = true, includePflow = false, includeQflow = false, noise = true, stddev = std, rng = rng)
end

function _add_pmu_measurements!(meas, net, rng; refShiftDeg::Float64)
  # PMU accuracy per IEEE C37.118 (TVE < 1 %): Vm ~0.002 p.u., Va ~0.02°.
  V = buildVoltageVector(net)
  σvm = 0.002
  σva = 0.02
  for i in PMU_BUS_IDXS
    push!(meas, Measurement(typ = Sparlectra.VmMeas, value = abs(V[i]) + randn(rng) * σvm, sigma = σvm, busIdx = i, id = "PMU_Vm_bus_$(i)"))
    push!(meas, Measurement(typ = Sparlectra.VaMeas, value = rad2deg(angle(V[i])) + refShiftDeg + randn(rng) * σva, sigma = σva, busIdx = i, id = "PMU_Va_bus_$(i)"))
  end
  return meas
end

function _report_scenario(io, name, net, meas, Vref; pmuRefOffset::Symbol = :auto)
  gobs = evaluate_global_observability(net, meas; pmuRefOffset = pmuRefOffset)
  se = runse!(net, meas; maxIte = 30, tol = 1e-10, flatstart = true, updateNet = false, pmuRefOffset = pmuRefOffset)
  rmse = _angle_rmse_deg(Vref, se.voltages)
  offset = se.vaRefOffsetDeg
  @printf(io, "%-34s %5d %8d %11s %12.5f %14.4e %12s\n", name, length(meas), gobs.n_states, string(gobs.quality), rmse, se.objectiveJ, offset === nothing ? "-" : @sprintf("%.4f", offset))
  return se
end

function _run_pmu_angle_example(io::IO)
  mkpath(OUTDIR)
  local_case = joinpath(Sparlectra.MPOWER_DIR, CASEFILE)
  case_path = if isfile(local_case)
    local_case
  else
    try
      Sparlectra.FetchMatpowerCase.ensure_casefile(CASEFILE; outdir = Sparlectra.MPOWER_DIR, to_jl = false, overwrite = false)
    catch err
      error("Could not obtain MATPOWER case $(CASEFILE). Please place it in $(Sparlectra.MPOWER_DIR) or allow network download. Original error: $(err)")
    end
  end

  net = createNetFromMatPowerFile(filename = case_path)
  it_pf, erg_pf = runpf!(net, 40, 1e-10, 0; method = :rectangular)
  erg_pf == 0 || error("Power flow did not converge")
  println(io, "PF reference solved in $(it_pf) iterations — used as ground truth.\n")
  Vref = buildVoltageVector(net)

  println(io, "PMU buses: ", PMU_BUS_IDXS, "; PMU sigma: Vm 0.002 p.u., Va 0.02 deg")
  println(io, "Scenario C/D: PMU time base shifted by +$(PMU_REF_SHIFT_DEG) deg against the slack reference.\n")

  @printf(io, "%-34s %5s %8s %11s %12s %14s %12s\n", "Scenario", "m", "n_state", "quality", "Va-RMSE_deg", "objective J", "alpha_deg")
  println(io, "-"^102)

  # A) SCADA baseline.
  measA = _scada_measurements(net, MersenneTwister(42))
  _report_scenario(io, "A SCADA only", net, measA, Vref)

  # B) SCADA + PMUs aligned with the slack reference.
  measB = _add_pmu_measurements!(_scada_measurements(net, MersenneTwister(42)), net, MersenneTwister(7); refShiftDeg = 0.0)
  _report_scenario(io, "B SCADA + PMU (aligned)", net, measB, Vref)

  # C) SCADA + PMUs on a shifted time base, offset modeled (:auto).
  measC = _add_pmu_measurements!(_scada_measurements(net, MersenneTwister(42)), net, MersenneTwister(7); refShiftDeg = PMU_REF_SHIFT_DEG)
  seC = _report_scenario(io, "C shifted PMU, offset state :auto", net, measC, Vref)

  # D) Same measurements, offset NOT modeled (:off).
  seD = _report_scenario(io, "D shifted PMU, no offset (:off)", net, measC, Vref; pmuRefOffset = :off)

  println(io)
  println(io, "Interpretation:")
  println(io, "  A->B: three accurate PMU angle pins sharpen the SCADA-only angle estimate.")
  @printf(io, "  C   : the offset state recovers alpha = %.4f deg (true shift %.1f deg);\n", seC.vaRefOffsetDeg, PMU_REF_SHIFT_DEG)
  println(io, "        network angles remain slack-referenced and unbiased.")
  @printf(io, "  D   : without the offset state the same data yields J = %.3e (vs %.3e in C)\n", seD.objectiveJ, seC.objectiveJ)
  println(io, "        and biased angles - the classic hybrid-SE pitfall.")
  println(io)

  # Bad-data view of scenario D: the unmodeled offset shows up as a block of
  # suspicious PMU angle residuals.
  diagD = validate_measurements(net, measC; pmuRefOffset = :off)
  println(io, "Diagnostics of scenario D (unmodeled reference shift):")
  print_se_diagnostics(diagD; io = io, topN = 6)

  return nothing
end

function main()
  print_example_banner("examples/state_estimation/state_estimation_pmu_angles.jl", "PMU voltage-angle measurements in WLS state estimation: reference-angle offset state alpha, aligned vs. shifted PMU time base, and the unmodeled-offset failure mode")
  mkpath(OUTDIR)
  timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
  logfile = joinpath(OUTDIR, "run_pmu_angles_$(CASEFILE)_$(timestamp).log")

  open(logfile, "w") do io
    redirect_stdout(io) do
      redirect_stderr(io) do
        _run_pmu_angle_example(io)
      end
    end
  end

  println("PMU angle state-estimation example finished. Log written to: $logfile")
  return logfile
end

run_example(main)
