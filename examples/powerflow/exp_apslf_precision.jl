# Copyright 2023-2026 Udo Schmitz
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
# file: examples/powerflow/exp_apslf_precision.jl
# purpose: APSLF on case6495rte with power_flow.apslf.precision float64 and
#          double_double, without NR polish, against Newton and against the
#          values measured when AnalyticLoadFlow 0.11.0 was released

# case6495rte stops above the tolerance in Float64 from order 32 on: the Padé
# value at s = 1 stagnates with a pole inside |s| = 1, and AnalyticLoadFlow
# reports the diagnosis precision_limit. power_flow.apslf.precision =
# double_double computes the series and the Padé evaluation in double-double
# arithmetic; without NR polish it needs order 48 or 56 on this case.
#
# Run from the Sparlectra checkout (the case is downloaded into the MATPOWER
# case cache on the first run):
#   julia --project=. examples/powerflow/exp_apslf_precision.jl
#   julia --project=. examples/powerflow/exp_apslf_precision.jl /path/to/case6495rte.m
#
# Mode :direct (the default), Q limits off. The "expected" columns are the
# results measured with AnalyticLoadFlow 0.11.0; another version or platform
# may differ in the last digits of the mismatch.

using Sparlectra
using Printf

# (precision, order, nr_polish, expected converged, expected mismatch)
const ROWS = [
  (:float64, 24, false, false, 5.79e-3),
  (:float64, 32, false, false, 1.66e-6),
  (:float64, 40, false, false, 2.21e-4),
  (:float64, 48, false, false, 1.45e-6),
  (:float64, 24, true, true, 7.96e-12),
  (:double_double, 32, false, false, 1.26e-6),
  (:double_double, 48, false, true, 7.10e-10),
  (:double_double, 56, false, true, 3.89e-11),
]

function run_case(casefile::AbstractString, pf::PowerFlowConfig)
  quiet = OutputConfig(logfile_results = :off, console_summary = false, startup_latency_hint = false)
  return run_sparlectra(casefile = casefile, config = SparlectraConfig(powerflow = pf, output = quiet))
end

function diagnosis_line(r)
  st = Sparlectra.rectangular_pf_status(r.net)
  return hasproperty(st, :apslf_diagnosis_line) ? string(st.apslf_diagnosis_line) : "-"
end

function main(args = ARGS)
  casefile = isempty(args) ? Sparlectra.ensure_casefile("case6495rte.m"; to_jl = false) : args[1]
  isfile(casefile) || error("case file not found: $(casefile)")
  noq = Sparlectra.QLimitConfig(ignore_q_limits = true)

  println("AnalyticLoadFlow: ", pkgversion(Sparlectra.AnalyticLoadFlow))
  println("Sparlectra:       ", pkgversion(Sparlectra))
  println("case:             ", casefile)
  println()

  # Newton is the voltage reference
  ref = run_case(casefile, PowerFlowConfig(solver = :rectangular, qlimits = noq))
  ref.final_converged || @warn "Newton reference did not converge" ref.reason
  Vref = Sparlectra.buildVoltageVector(ref.net)
  @printf("Newton reference: converged=%s  mismatch=%.2e  %.2f s\n\n", ref.final_converged, ref.final_mismatch, ref.elapsed_s)

  @printf("%-14s %-6s %-7s %-10s %-10s %-10s %-10s %-10s %-7s\n", "precision", "order", "polish", "converged", "expected", "mismatch", "expected", "max|dV|", "time s")
  println("-"^92)
  for (prec, order, polish, exp_conv, exp_mis) in ROWS
    pf = PowerFlowConfig(solver = :apslf, apslf = Sparlectra.ApslfConfig(order = order, nr_polish = polish, precision = prec), qlimits = noq)
    r = run_case(casefile, pf)
    V = Sparlectra.buildVoltageVector(r.net)
    dv = all(isfinite, V) ? maximum(abs.(V .- Vref)) : NaN
    @printf("%-14s %-6d %-7s %-10s %-10s %-10.2e %-10.2e %-10.2e %-7.2f\n", prec, order, polish ? "yes" : "no", r.final_converged, exp_conv, r.final_mismatch, exp_mis, dv, r.elapsed_s)
    r.final_converged || println("    ", diagnosis_line(r))
  end
  println()
  println("Expected diagnoses: float64 order 24 suggests order 32; float64 orders 32,")
  println("40 and 48 report precision_limit; double_double order 32 suggests order 48.")
  return nothing
end

Base.invokelatest(main)
