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

# Date: 2026-07-29
# file: examples/powerflow/exp_distributed_slack_modes.jl
# purpose: compares the classical single-slack solve with the distributed
# active-power slack (issue #192) in pg_weighted and imported mode on a small
# MATPOWER case carrying the APF column

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

# Small 4-bus case with a deliberate 20 MW shortfall (70 MW load vs 50 MW
# scheduled PV generation) and the optional 21st gen column (APF): the two PV
# units carry imported participation factors 0.75 / 0.25 — deliberately
# different from their 30/20 Pg ratio, so pg_weighted and imported mode give
# visibly different shares.
function _write_demo_case(dir::AbstractString)
  path = joinpath(dir, "case4_distributed_slack.m")
  write(path, """
function mpc = case4_distributed_slack
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
1 3 0 0 0 0 1 1.00 0 110 1 1.1 0.9;
2 2 0 0 0 0 1 1.01 0 110 1 1.1 0.9;
3 2 0 0 0 0 1 1.00 0 110 1 1.1 0.9;
4 1 70 20 0 0 1 1.00 0 110 1 1.1 0.9;
];
mpc.gen = [
1 0 0 300 -300 1.00 100 1 250 -250 0 0 0 0 0 0 0 0 0 0 0;
2 30 0 50 -50 1.01 100 1 100 0 0 0 0 0 0 0 0 0 0 0 0.75;
3 20 0 50 -50 1.00 100 1 100 0 0 0 0 0 0 0 0 0 0 0 0.25;
];
mpc.branch = [
1 2 0.01 0.08 0.02 0 0 0 0 0 1 -360 360;
2 3 0.02 0.10 0.02 0 0 0 0 0 1 -360 360;
3 4 0.01 0.06 0.02 0 0 0 0 0 1 -360 360;
1 4 0.03 0.12 0.02 0 0 0 0 0 1 -360 360;
];
""")
  return path
end

# The interesting quantity per run: which bus ends up delivering how much
# active power beyond its schedule.
function _print_bus_p(net)
  Y = Sparlectra.createYBUS(net = net, sparse = true, printYBUS = false)
  V = [n._vm_pu * cis(deg2rad(n._va_deg)) for n in net.nodeVec]
  extra = (real.(Sparlectra.calc_injections(Y, V)) .- real.(Sparlectra.buildComplexSVec(net))) .* net.baseMVA
  for (i, n) in enumerate(net.nodeVec)
    abs(extra[i]) < 0.005 && continue
    println("    bus ", i, " (", Sparlectra.getCompName(n.comp), ")  +", round(extra[i]; digits = 2), " MW beyond schedule")
  end
end

function main()
  print_example_banner("examples/powerflow/exp_distributed_slack_modes.jl", "compares the classical single slack with the distributed active-power slack (pg_weighted and imported APF weights)")
  case_path = _write_demo_case(mktempdir())

  println("=== 1) Classical single slack (default) ===")
  net = createNetFromMatPowerFile(filename = case_path)
  _, erg = runpf_rectangular!(net, 30, 1e-8, 0)
  println("  converged: ", erg == 0)
  println("  the REF bus alone absorbs the imbalance:")
  _print_bus_p(net)

  println()
  println("=== 2) Distributed slack, pg_weighted (shares 30:20 = 0.6/0.4) ===")
  net = createNetFromMatPowerFile(filename = case_path)
  _, erg = runpf_rectangular!(net, 30, 1e-8, 1; distributed_slack_enabled = true, distributed_slack_p_mode = :pg_weighted)
  println("  converged: ", erg == 0)
  _print_bus_p(net)

  println()
  println("=== 3) Distributed slack, imported (MATPOWER APF 0.75/0.25) ===")
  net = createNetFromMatPowerFile(filename = case_path)
  _, erg = runpf_rectangular!(net, 30, 1e-8, 1; distributed_slack_enabled = true, distributed_slack_p_mode = :imported)
  println("  converged: ", erg == 0)
  _print_bus_p(net)
  st = Sparlectra.rectangular_pf_status(net)
  println()
  println("  result metadata: lambda_P = ", round(st.distributed_slack_lambda_p_mw; digits = 3), " MW over ", st.distributed_slack_participants, " participants (mode ", st.distributed_slack_mode, ")")

  return nothing
end

run_example(main)
