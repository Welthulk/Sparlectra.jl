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

# file: examples/cgmes/cgmes_export_physical.jl
# purpose: builds a small 110 kV net from line data in ohm/km and nF/km,
#          writes it as a CGMES delivery with writeCGMESFiles and reads it
#          back: the exported ACLineSegment values are the physical ones and
#          the re-import solves to the same voltages.

using Sparlectra

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

# Line data per km as a line catalog gives it.
const LINES = [
  # from, to, length km, r ohm/km, x ohm/km, c nF/km
  ("A", "B", 25.0, 0.12, 0.39, 9.5),
  ("B", "C", 15.0, 0.12, 0.39, 9.5),
  ("C", "A", 30.0, 0.06, 0.30, 12.0),
]

function build_net()
  net = Net(name = "cgmes_physical", baseMVA = 100.0)
  for bus in ("A", "B", "C")
    addBus!(net = net, busName = bus, vn_kV = 110.0)
  end
  addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "A")
  addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 40.0, q = 12.0)
  addProsumer!(net = net, busName = "C", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
  for (from, to, km, r, x, c) in LINES
    # b from the capacitance at the system frequency of the net
    # (net.nominal_frequency_hz, a new Net takes model.nominal_frequency_hz)
    addACLinePerKm!(net = net, fromBus = from, toBus = to, length_km = km, r_ohm_per_km = r, x_ohm_per_km = x, c_nf_per_km = c)
  end
  return net
end

function main(; output_dir::AbstractString = joinpath(@__DIR__, "..", "_out", "cgmes_export_physical"))
  print_example_banner("examples/cgmes/cgmes_export_physical.jl", "CGMES export of a net built from line data in ohm/km and nF/km, with the re-import")
  net = build_net()
  f_hz = net.nominal_frequency_hz
  ite, erg = runpf!(net, 30, 1e-10, 0)
  erg == 0 || error("power flow did not converge")
  println("power flow: converged in ", ite, " iterations")

  mkpath(output_dir)
  files = writeCGMESFiles(net; path = output_dir)

  # the delivery carries the physical line totals
  back = importCGMES(path = output_dir, name = "cgmes_physical_back")
  println()
  println("line      r ohm in / out      x ohm in / out      bch uS in / out")
  names = Sparlectra._bus_name_by_idx(back.net)
  for seg in Sparlectra.CGMESImporter.objectsOf(back.store, :ACLineSegment)
    br = back.net.branchVec[last(back.element_by_mrid[seg.mrid])]
    from, to, km, r, x, c = only(l for l in LINES if Set((l[1], l[2])) == Set((names[br.fromBus], names[br.toBus])))
    r_out = Sparlectra.CGMESImporter.num(seg, :r)
    x_out = Sparlectra.CGMESImporter.num(seg, :x)
    b_out = Sparlectra.CGMESImporter.num(seg, :bch)
    println(rpad("$(from)-$(to)", 10), rpad("$(round(r * km; digits = 4)) / $(round(r_out; digits = 4))", 20),
            rpad("$(round(x * km; digits = 4)) / $(round(x_out; digits = 4))", 20),
            round(2π * f_hz * c * 1e-9 * km * 1e6; digits = 3), " / ", round(b_out * 1e6; digits = 3))
  end

  # the re-import solves to the voltages of the source net
  last(runpf!(back.net, 30, 1e-10, 0)) == 0 || error("re-import did not converge")
  dvm = maximum(abs(back.net.nodeVec[back.net.busDict[bus]]._vm_pu - net.nodeVec[net.busDict[bus]]._vm_pu) for bus in ("A", "B", "C"))
  println()
  println("max |dvm| source against re-import: ", dvm, " pu")
  return files
end

files = run_example(main)

for f in files
  println("written: ", f)
end
