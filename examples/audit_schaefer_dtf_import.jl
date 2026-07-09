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

const CASES = ["A" => "", "B" => "B", "C" => "C", "D" => "D", "E" => "E"]

finite_range(xs) = isempty(xs) ? "n/a" : @sprintf("[%g, %g], finite=%s", minimum(xs), maximum(xs), all(isfinite, xs))

function audit_case(io, case_id, path)
  case = Sparlectra.DTFImporter.read_dtf(path; strict = true)
  println(io, "## Case ", case_id)
  println(io, "- file: ", path)
  println(io, "- parsed branch count: ", length(case.branches))
  println(io, "- parsed bus count: ", length(case.buses))
  println(io, "- parsed transformer control count: ", length(case.transformer_controls))
  println(io, "- parsed outage/trailing-section count: ", length(case.outages), " outages; ", length(case.trailing_records), " trailing records")
  println(io, "- slack bus: ", case.size.slack)
  println(io, "- PV buses: ", join((b.name for b in case.buses if b.bus_type == 1), ", "))
  println(io, "- PQ buses: ", join((b.name for b in case.buses if b.bus_type in (0, 3)), ", "))
  println(io, "- branch r ohm: ", finite_range([b.r_ohm for b in case.branches]))
  println(io, "- branch x ohm: ", finite_range([b.x_ohm for b in case.branches]))
  println(io, "- branch g S: ", finite_range([b.g_s for b in case.branches]))
  println(io, "- branch b S: ", finite_range([b.b_s for b in case.branches]))
  println(io, "\n### Transformer controls")
  for c in case.transformer_controls
    println(io, "- source line ", c.index, ": ", c.from, " -> ", c.to, " parallel=", c.parallel_id,
      " nominal_unregulated_kV=", c.nominal_unregulated_kv,
      " nominal_regulated_kV=", c.nominal_regulated_kv,
      " longitudinal_range_percent=", c.longitudinal_range_percent,
      " actual_tap_step=", c.actual_tap_step,
      " max_tap_step=", c.max_tap_step,
      " added_voltage_angle_deg=", c.added_voltage_angle_deg)
  end
  println(io, "\n### Node start voltages")
  for b in case.buses
    vn = case.nominal_voltages_kv[b.voltage_level_index]
    vm = b.start_kv > 0 ? b.start_kv / vn : 1.0
    println(io, "- ", b.name, ": start_kV=", b.start_kv, " vn_kV=", vn, " vm_pu=", vm, " finite=", all(isfinite, (b.start_kv, vn, vm)))
  end
  println(io, "\n### Post-bus/trailing records")
  if isempty(case.trailing_records)
    println(io, "- none")
  else
    for r in case.trailing_records
      println(io, "- raw=", repr(r.raw), "; interpreted kind=", r.interpreted_kind)
    end
  end
  println(io)
  return nothing
end

function main()
  repo = normpath(joinpath(@__DIR__, ".."))
  outdir = joinpath(@__DIR__, "_out", "schaefer_dtf_validation")
  mkpath(outdir)
  out = joinpath(outdir, "dtf_import_audit.md")
  open(out, "w") do io
    println(io, "# Schäfer DTF/FOR001 import audit\n")
    println(io, "Standalone post-bus `A` cards in cases B-E are interpreted as variant branch-echo markers; the following branch-like payload line is preserved as a branch echo and is not added to the electrical model.\n")
    for (case_id, suffix) in CASES
      audit_case(io, case_id, joinpath(repo, "data", "DTF", "FOR001$(suffix).DAT"))
    end
  end
  println("Wrote ", out)
  return nothing
end

Base.invokelatest(main)
