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

# file: tools/build_dtf_demo.jl
# purpose: write the shipped DTF demo deck data/dtf_demo/sp_dtf5.DAT, a
#          self-built five-bus network in the fixed-column FOR001 layout
#          (two voltage levels, a ring of 110 kV lines, one transformer
#          with a longitudinal tap control, two outage records). Every
#          number is invented for this deck; nothing is taken from a deck
#          of another origin. The cards are written column by column, so a
#          rerun leaves the file byte-identical.
#
# usage: julia --project=. tools/build_dtf_demo.jl

using Sparlectra
using Printf

const _DTF_DEMO_PATH = normpath(joinpath(@__DIR__, "..", "data", "dtf_demo", "sp_dtf5.DAT"))

# columns 1 to 3: kind, voltage-level index, parallel identifier; the bus
# names stand in columns 6 to 13 and 16 to 23
_dtf_head(kind::Char, level::Int, from::String, to::String) = @sprintf("%c%d   %-8s  %-8s", kind, level, from, to)

# branch card: series resistance and reactance in ohm, shunt conductance
# and susceptance in siemens, current limit in kA
_dtf_branch(kind::Char, level::Int, from::String, to::String, r::Float64, x::Float64, g::Float64, b::Float64, imax::Float64) = string(_dtf_head(kind, level, from, to), @sprintf(" %9.4f %9.4f %11.4e %11.4e %7.3f", r, x, g, b, imax))

# transformer control card: five-column fields from column 26 on (winding
# voltages, longitudinal range in percent, angle of the added voltage,
# highest step, present step)
_dtf_control(from::String, to::String, u_unregulated::Float64, u_regulated::Float64, range_percent::Float64, angle_deg::Float64, max_step::Int, step::Int) =
  string(@sprintf("     %-8s  %-8s  ", from, to), @sprintf("%5.1f%5.1f%5.1f%5.1f%5d%5d", u_unregulated, u_regulated, range_percent, angle_deg, max_step, step))

# bus card: type (0 load bus, 1 voltage-controlled, 2 slack), level index,
# name in columns 6 to 13, then start voltage in kV, angle, load P and Q,
# generation P and Q, reactive limits
_dtf_bus(type::Int, level::Int, name::String, u_kv::Float64, pd::Float64, qd::Float64, pg::Float64, qg::Float64, qmin::Float64, qmax::Float64) = string(@sprintf("%d%d   %-8s", type, level, name), @sprintf(" %8.3f %7.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f", u_kv, 0.0, pd, qd, pg, qg, qmin, qmax))

function dtf_demo_cards()::Vector{String}
  branches = [
    # a ring of 110 kV overhead lines (0.10 + j0.40 ohm/km, 2.8e-6 S/km)
    _dtf_branch('L', 1, "NORD", "WEST", 2.0, 8.0, 0.0, 5.6e-5, 0.500),
    _dtf_branch('L', 1, "NORD", "OST", 2.5, 10.0, 0.0, 7.0e-5, 0.500),
    _dtf_branch('L', 1, "WEST", "SUED", 3.0, 12.0, 0.0, 8.4e-5, 0.400),
    _dtf_branch('L', 1, "OST", "SUED", 2.0, 8.0, 0.0, 5.6e-5, 0.400),
    _dtf_branch('L', 1, "WEST", "OST", 3.5, 14.0, 0.0, 9.8e-5, 0.400),
    # 40 MVA, 110/20 kV, 12 percent short-circuit voltage, referred to 110 kV
    _dtf_branch('T', 1, "SUED", "STADT", 1.5, 36.3, 1.0e-6, -5.0e-6, 0.210),
  ]
  controls = [_dtf_control("SUED", "STADT", 110.0, 20.0, 16.0, 0.0, 8, 2)]
  buses = [
    _dtf_bus(2, 1, "NORD", 112.2, 0.0, 0.0, 0.0, 0.0, -60.0, 60.0),
    _dtf_bus(0, 1, "WEST", 110.0, 30.0, 10.0, 0.0, 0.0, 0.0, 0.0),
    _dtf_bus(1, 1, "OST", 112.0, 0.0, 0.0, 35.0, 0.0, -20.0, 25.0),
    _dtf_bus(0, 1, "SUED", 110.0, 25.0, 8.0, 0.0, 0.0, 0.0, 0.0),
    _dtf_bus(0, 2, "STADT", 20.0, 18.0, 6.0, 0.0, 0.0, 0.0, 0.0),
  ]
  outages = [_dtf_head('L', 1, "WEST", "OST"), _dtf_head('L', 1, "NORD", "WEST")]
  return vcat(
    ["  50  1.0E-08"],
    ["SPARLECTRA DTF DEMO DECK SP_DTF5", "SELF-BUILT FIVE-BUS NETWORK, INVENTED DATA", "110 KV RING, ONE 110/20 KV TRANSFORMER WITH TAP CONTROL", "TWO OUTAGE RECORDS BETWEEN AUSFALL AND ENDE"],
    ["110.0 20.0"],
    [@sprintf("%5d%5d%5d%5d  %s", length(buses), length(branches), 0, length(controls), "NORD")],
    branches,
    controls,
    buses,
    ["AUSFALL"],
    outages,
    ["ENDE"],
  )
end

function main()
  mkpath(dirname(_DTF_DEMO_PATH))
  open(_DTF_DEMO_PATH, "w") do io
    for card in dtf_demo_cards()
      println(io, card)
    end
  end
  # the deck is only worth shipping when it reads and solves
  case = Sparlectra.DTFImporter.read_dtf(_DTF_DEMO_PATH)
  net = Sparlectra.DTFImporter.build_net(case)
  ite, erg = runpf!(net, 50, 1e-8, 0)
  erg == 0 || error("the demo deck does not solve (erg = $(erg))")
  calcNetLosses!(net)
  @printf("%s: %d buses, %d branches, %d control(s), %d outage record(s); solved in %d iterations, losses %.3f MW\n", basename(_DTF_DEMO_PATH), length(case.buses), length(case.branches), length(case.transformer_controls), length(case.outages), ite, net.totalLosses[end][1])
  for nd in net.nodeVec
    @printf("  %-8s %7.4f pu %8.3f deg\n", nd.comp.cName, nd._vm_pu, nd._va_deg)
  end
  return 0
end

Base.invokelatest(main)
