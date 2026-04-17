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
# Example: First classical WLS state estimation workflow
#
# Steps:
# 1) Build/import network
# 2) Solve AC power flow
# 3) Generate synthetic measurements from solved state
# 4) Run WLS state estimation
# 5) Compare estimated voltages with PF reference

# file: src/examples/state_estimation_wls.jl

using Sparlectra
using Printf
using Dates
using Random

const OUTDIR = joinpath(@__DIR__, "_out")
const CASEFILE = "case9.m"

@inline function _measurement_unit(typ::Sparlectra.MeasurementType)
  if typ == Sparlectra.VmMeas
    return "kV"
  elseif typ == Sparlectra.PinjMeas || typ == Sparlectra.PflowMeas
    return "MW"
  elseif typ == Sparlectra.QinjMeas || typ == Sparlectra.QflowMeas
    return "MVAr"
  end
  return "-"
end

function _run_state_estimation_example(io::IO)
  # Ensure output directory exists and fetch the MATPOWER case if needed.
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

  # Build net from the case file and solve the AC power flow.
  net = createNetFromMatPowerFile(filename = case_path)
  it_pf, erg_pf = runpf!(net, 40, 1e-10, 0; method = :rectangular, opt_sparse = true)
  erg_pf == 0 || error("Power flow did not converge")

  # Keep PF voltages as "ground truth" for the synthetic measurement study.
  Vref = buildVoltageVector(net)

  # Configure synthetic measurement standard deviations.
  std = measurementStdDevs(vm = 0.01, pinj = 1.5, qinj = 1.5, pflow = 0.7, qflow = 0.9)

  # Generate all currently supported measurement types with Gaussian noise.
  rng = MersenneTwister(42)
  setMeasurementsFromPF!(net; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = true, stddev = std, rng = rng)
  meas = Measurement[m for m in net.measurements]

  # Reset to a simple flat start before running the state estimator.
  for n in net.nodeVec
    if getNodeType(n) != :Slack
      n._vm_pu = 1.0
      n._va_deg = 0.0
    end
  end

  # Run first WLS state estimation implementation.
  se = runse!(net; maxIte = 12, tol = 1e-6, flatstart = false, jacEps = 1e-6, updateNet = true)

  # Compare estimated voltages against PF reference values.
  Vse = se.voltages
  vm_err = maximum(abs.(abs.(Vref) .- abs.(Vse)))
  va_err = maximum(abs.(rad2deg.(angle.(Vref)) .- rad2deg.(angle.(Vse))))

  # Build measurement predictions for PF reference and SE estimate for tabular output.
  Ybus = createYBUS(net = net)
  Sref_MVA = calc_injections(Ybus, Vref) .* net.baseMVA
  Sse_MVA = calc_injections(Ybus, Vse) .* net.baseMVA

  # Residual diagnostics in physical units and 3σ logic for bad-data localization.
  max_abs_norm_res = 0.0
  out_of_3sigma = 0
  for m in meas
    h_se = Sparlectra._measurement_prediction(m, net, Vse, Sse_MVA)
    residual = h_se - m.value
    abs_norm_res = abs(residual) / m.sigma
    max_abs_norm_res = max(max_abs_norm_res, abs_norm_res)
    out_of_3sigma += abs_norm_res > 3.0
  end

  # Write a compact summary to the selected output stream.
  @printf(io, "Case: %s\n", CASEFILE)
  @printf(io, "PF iterations: %d\n", it_pf)
  @printf(io, "SE converged: %s in %d iterations\n", string(se.converged), se.iterations)
  @printf(io, "Residual norm ||r||2 (in measurement units): %.6e\n", se.residualNorm)
  @printf(io, "WLS objective J = r'Wr: %.6e (dof=%.0f)\n", se.objectiveJ, Float64(se.dof))
  @printf(io, "J in 3σ-band: %s\n", string(se.jWithin3Sigma))
  @printf(io, "Max |r/σ|: %.4f\n", max_abs_norm_res)
  @printf(io, "Measurements outside 3σ: %d / %d\n", out_of_3sigma, length(meas))
  @printf(io, "Max |ΔVm| (p.u.): %.6e\n", vm_err)
  @printf(io, "Max |ΔVa| (deg): %.6e\n", va_err)

  # Bus-wise voltage table in kV.
  println(io, "\nBus voltage comparison (kV)")
  println(io, "-------------------------------------------------------------------------------------------------")
  @printf(io, "%5s %12s %14s %14s %14s %14s\n", "Bus", "Vbase kV", "PF V", "Meas V", "SE V", "SE-PF")
  println(io, "-------------------------------------------------------------------------------------------------")
  for i in eachindex(net.nodeVec)
    vm_meas = NaN
    for m in meas
      if m.typ == Sparlectra.VmMeas && m.busIdx == i
        vm_meas = m.value
        break
      end
    end
    vbase_kV = getNodeVn(net.nodeVec[i])
    v_pf_kV = abs(Vref[i]) * vbase_kV
    v_se_kV = abs(Vse[i]) * vbase_kV
    v_meas_kV = isnan(vm_meas) ? NaN : vm_meas * vbase_kV
    @printf(io, "%5d %12.3f %14.5f %14.5f %14.5f %14.5e\n", i, vbase_kV, v_pf_kV, v_meas_kV, v_se_kV, v_se_kV - v_pf_kV)
  end

  # Injection table (bus measurements) in MW / MVAr.
  println(io, "\nBus injection comparison (MW / MVAr)")
  println(io, "---------------------------------------------------------------------------------------------------------------")
  @printf(io, "%5s %12s %12s %12s %12s %12s %12s %12s %12s\n", "Bus", "PF P", "Meas P", "SE P", "SE-PF P", "PF Q", "Meas Q", "SE Q", "SE-PF Q")
  println(io, "---------------------------------------------------------------------------------------------------------------")
  for i in eachindex(net.nodeVec)
    p_meas = NaN
    q_meas = NaN
    for m in meas
      if m.busIdx == i
        if m.typ == Sparlectra.PinjMeas
          p_meas = m.value
        elseif m.typ == Sparlectra.QinjMeas
          q_meas = m.value
        end
      end
    end
    p_pf = real(Sref_MVA[i])
    q_pf = imag(Sref_MVA[i])
    p_se = real(Sse_MVA[i])
    q_se = imag(Sse_MVA[i])
    @printf(io, "%5d %12.4f %12.4f %12.4f %12.4e %12.4f %12.4f %12.4f %12.4e\n", i, p_pf, p_meas, p_se, p_se - p_pf, q_pf, q_meas, q_se, q_se - q_pf)
  end

  # Measurement-wise table: PF prediction, measured value, SE prediction, residual and 3σ indicator.
  println(io, "\nMeasurement table with 3σ check")
  println(io, "----------------------------------------------------------------------------------------------------------------------------------------------------------------")
  @printf(io, "%4s %-10s %-20s %7s %12s %12s %12s %12s %12s %9s %8s\n", "#", "Type", "Target", "Unit", "PF", "Measured", "SE", "SE-PF", "SE-Meas", "|r|/σ", "3σ")
  println(io, "----------------------------------------------------------------------------------------------------------------------------------------------------------------")
  for (k, m) in enumerate(meas)
    target = if !isnothing(m.busIdx)
      "bus=$(m.busIdx)"
    else
      "br=$(m.branchIdx),$(m.direction)"
    end

    h_pf = Sparlectra._measurement_prediction(m, net, Vref, Sref_MVA)
    h_se = Sparlectra._measurement_prediction(m, net, Vse, Sse_MVA)
    unit = _measurement_unit(m.typ)

    scale = 1.0
    if m.typ == Sparlectra.VmMeas
      bus_idx = something(m.busIdx, 0)
      bus_idx > 0 || error("Vm measurement missing busIdx for unit conversion")
      scale = getNodeVn(net.nodeVec[bus_idx])
    end

    h_pf_u = h_pf * scale
    z_u = m.value * scale
    h_se_u = h_se * scale
    se_pf_u = h_se_u - h_pf_u
    r_u = h_se_u - z_u
    sigma_u = m.sigma * scale
    abs_norm_r = abs(r_u) / sigma_u
    sigma_flag = abs_norm_r <= 3.0 ? "in" : "out"

    @printf(io, "%4d %-10s %-20s %7s %12.4f %12.4f %12.4f %12.4f %12.4f %9.3f %8s\n", k, string(m.typ), target, unit, h_pf_u, z_u, h_se_u, se_pf_u, r_u, abs_norm_r, sigma_flag)
  end

  return se
end

function run_state_estimation_example()
  mkpath(OUTDIR)
  timestamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
  logfile = joinpath(OUTDIR, "run_$(CASEFILE)_$(timestamp).log")

  open(logfile, "w") do io
    # Redirect both stdout and stderr to keep all example output in one file.
    redirect_stdout(io) do
      redirect_stderr(io) do
        _run_state_estimation_example(io)
      end
    end
  end

  println("State estimation example finished. Log written to: $logfile")
  return logfile
end

run_state_estimation_example()
