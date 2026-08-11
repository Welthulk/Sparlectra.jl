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

# Date: 2026-08-09
# file: examples/powerflow/exp_external_grid_comparison.jl
# purpose: 8-bus comparison of grid-connection modeling (issue #299): ideal
#          slack at two alternative buses, non-ideal external-grid source
#          (internal impedance) at the same two buses, and the distributed
#          active-power slack — tabulated per bus so the differences between
#          the representations are directly visible.

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

# Meshed 110 kV ring (B1..B8) with two chords, two PV generators and four
# loads. Scheduled generation (100 MW) deliberately undershoots the load
# (160 MW), so the grid connection has to import a visible amount of power —
# that import is exactly what the comparison tables trace.
#
#   B1 -- B2 -- B3 -- B4 -- B5
#    |     |     |          |
#   B8 -- B7 ----+--- B6 ---+   (chords: B2–B7, B3–B6)
function build_eg8(; connection::String, mode::Symbol, sk_MVA::Float64 = 2000.0, rx::Float64 = 0.1)
  net = Net(name = "eg8_$(connection)_$(mode)", baseMVA = 100.0)
  for b in ("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.010, x_pu = 0.060, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B3", r_pu = 0.015, x_pu = 0.080, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B4", r_pu = 0.020, x_pu = 0.090, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B4", toBus = "B5", r_pu = 0.012, x_pu = 0.070, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B5", toBus = "B6", r_pu = 0.015, x_pu = 0.075, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B6", toBus = "B7", r_pu = 0.018, x_pu = 0.085, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B7", toBus = "B8", r_pu = 0.010, x_pu = 0.055, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B8", toBus = "B1", r_pu = 0.011, x_pu = 0.065, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B2", toBus = "B7", r_pu = 0.020, x_pu = 0.100, b_pu = 0.02, status = 1)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B6", r_pu = 0.022, x_pu = 0.110, b_pu = 0.02, status = 1)

  addProsumer!(net = net, busName = "B3", type = "GENERATOR", p = 60.0, vm_pu = 1.01, qMin = -60.0, qMax = 60.0)
  addProsumer!(net = net, busName = "B6", type = "GENERATOR", p = 40.0, vm_pu = 1.00, qMin = -40.0, qMax = 40.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 45.0, q = 12.0)
  addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 50.0, q = 15.0)
  addProsumer!(net = net, busName = "B7", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
  addProsumer!(net = net, busName = "B8", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)

  # The grid connection: mode :slack is the ideal representation (the
  # connection bus becomes REF), mode :source the non-ideal one (hidden
  # internal bus behind z = Un²/Sk — the terminal bus stays an ordinary bus).
  addExternalGrid!(net = net, busName = connection, vm_pu = 1.02, sk_max_MVA = sk_MVA, rx_max = rx, internal_impedance = (mode === :source))
  return net
end

# Per-bus absorbed active power P_calc(V) − P_spec in MW at the solved state:
# zero everywhere except where a free injection covers imbalance plus losses.
# This is the single most telling row-wise difference between the scenarios.
function absorbed_MW(net)
  Y = Sparlectra.createYBUS(net = net, sparse = true, printYBUS = false)
  V = [n._vm_pu * cis(deg2rad(n._va_deg)) for n in net.nodeVec]
  resid = real.(Sparlectra.calc_injections(Y, V)) .- real.(Sparlectra.buildComplexSVec(net))
  return resid .* net.baseMVA
end

function run_scenario(label::String; connection::String, mode::Symbol, distributed::Bool = false)
  net = build_eg8(connection = connection, mode = mode)
  ite, erg = runpf!(net, 30, 1e-8, 0; distributed_slack_enabled = distributed, distributed_slack_p_mode = :pg_weighted)
  erg == 0 || error("scenario $(label): power flow did not converge (erg = $(erg))")
  names = [name for (name, _) in sort(collect(net.busDict); by = last)]
  return (
    label = label,
    net = net,
    ite = ite,
    bus_names = names,
    vm = [n._vm_pu for n in net.nodeVec],
    va = [n._va_deg for n in net.nodeVec],
    absorbed = absorbed_MW(net),
  )
end

_fmt(x; digits = 4) = x === nothing ? "—" : string(round(x; digits = digits))

# One row per physical bus (B1..B8) plus one synthetic row for the hidden
# internal source bus of the :source scenarios. Values are looked up by bus
# name so scenarios with different bus counts align correctly.
function print_comparison(io::IO, title::String, scenarios, field::Symbol; digits::Int = 4)
  println(io)
  println(io, title)
  colw = 14
  print(io, rpad("bus", 10))
  for s in scenarios
    print(io, lpad(s.label, colw))
  end
  println(io)
  rows = ["B$(i)" for i in 1:8]
  push!(rows, "internal src")
  for row in rows
    print(io, rpad(row, 10))
    for s in scenarios
      idx = if row == "internal src"
        findfirst(endswith("__extgrid_int"), s.bus_names)
      else
        findfirst(==(row), s.bus_names)
      end
      val = idx === nothing ? nothing : getfield(s, field)[idx]
      print(io, lpad(_fmt(val; digits = digits), colw))
    end
    println(io)
  end
end

function main()
  print_example_banner(
    "examples/powerflow/exp_external_grid_comparison.jl",
    "8-bus comparison: ideal slack vs. non-ideal external-grid source at two alternative connection buses, plus the distributed slack",
  )

  scenarios = [
    run_scenario("Slack@B1"; connection = "B1", mode = :slack),
    run_scenario("Slack@B5"; connection = "B5", mode = :slack),
    run_scenario("Source@B1"; connection = "B1", mode = :source),
    run_scenario("Source@B5"; connection = "B5", mode = :source),
    run_scenario("DistSlack@B1"; connection = "B1", mode = :slack, distributed = true),
  ]

  println("Scenarios (all solved with the rectangular Newton-Raphson, tol 1e-8):")
  for s in scenarios
    println("  ", rpad(s.label, 14), " buses=", length(s.bus_names), "  iterations=", s.ite)
  end
  println()
  println("Slack@Bx    : ideal slack (addExternalGrid!, Sk'' = 2000 MVA carried as data only)")
  println("Source@Bx   : non-ideal source (internal_impedance = true, z_pu = baseMVA/Sk'' = 0.05, R/X = 0.1)")
  println("DistSlack@B1: ideal slack at B1 + distributed_slack_enabled (pg_weighted: G3 α=0.6, G6 α=0.4)")

  print_comparison(stdout, "Voltage magnitude Vm [pu]:", scenarios, :vm)
  print_comparison(stdout, "Voltage angle Va [deg] (note: Source scenarios anchor the angle at the internal source bus):", scenarios, :va)
  print_comparison(stdout, "Absorbed active power P_calc − P_spec [MW] (who covers imbalance + losses):", scenarios, :absorbed; digits = 2)

  println()
  println("Reading aid:")
  println("  • Slack@B1 vs Slack@B5: same physics, different boundary condition — the whole")
  println("    voltage profile and the import path change with the slack choice.")
  println("  • Slack@Bx vs Source@Bx: the source's internal impedance lets the connection-bus")
  println("    voltage droop below the 1.02 pu setpoint and shifts all angles; the absorbed")
  println("    power moves to the hidden internal bus (the actual reference).")
  println("  • DistSlack@B1: the imbalance no longer lands on B1 alone — G3/G6 pick up their")
  println("    α-shares (0.6/0.4) while the voltage profile stays close to Slack@B1.")

  for s in scenarios
    total = sum(s.absorbed)
    println()
    println(s.label, ": total absorbed = ", round(total; digits = 2), " MW (scheduled imbalance 60 MW + losses ", round(total - 60.0; digits = 2), " MW)")
  end

  return scenarios
end

run_example(main)
