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

# file: examples/experimental/cgmes_export_demo.jl
# purpose: demonstrates the experimental Stage-1 CGMES export on a small
#          programmatically built network.
#
# Experimental example tooling. Not part of the stable Sparlectra API, not
# wired into the Web UI, and may change without notice.

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))
include(joinpath(@__DIR__, "cgmes_export.jl"))

function main(; output_dir::AbstractString = joinpath(@__DIR__, "..", "_out", "cgmes_export"))
  print_example_banner("examples/experimental/cgmes_export_demo.jl", "demonstrates the experimental Stage-1 CGMES export on a small programmatically built network")

  net = Net(name = "cgmes_demo", baseMVA = 100.0)

  addBus!(net = net, busName = "B1", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "B2", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")

  # Lines with physical parameters (Ohm/km) ...
  addACLine!(net = net, fromBus = "B1", toBus = "B2", length = 25.0, r = 0.2, x = 0.39)
  addACLine!(net = net, fromBus = "B2", toBus = "B3", length = 15.0, r = 0.2, x = 0.39)
  # ... and one PI-model line (p.u.) — converted back to Ohm during export.
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B1", r_pu = 0.02, x_pu = 0.10, b_pu = 0.01, status = 1)

  ok, msg = validate!(net = net)
  ok || error("Invalid network: $msg")

  # Optional zero-sequence/short-circuit data per line index (order of
  # net.linesAC). Only supply what is actually known — missing entries are
  # not written.
  sc = Dict{Int,CGMESLineShortCircuit}(
    1 => CGMESLineShortCircuit(r0_ohm = 15.0, x0_ohm = 29.3, b0ch_S = 0.0, endTemperature_C = 160.0),
    2 => CGMESLineShortCircuit(r0_ohm = 9.0, x0_ohm = 17.6),
  )

  mkpath(output_dir)
  files = writeCGMESFiles(net; path = output_dir, sc_line_data = sc)
  return files
end

files = run_example(main)

for f in files
  println("written: ", f)
end
