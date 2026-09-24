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

# Date: 2026-09-24
# file: examples/powerflow/example_qlimit_reenable_voltage_rule.jl
# purpose: shows the voltage-side PQ->PV release of the active-set Q-limit
#          mode on the Zeng/Chiang 14-bus case: with the release switched off
#          the run ends on a non-physical solution (three machines at Qmax
#          with the voltage above the setpoint, reported by the Q-V check),
#          with the rule it lands on the physical solution that the classic
#          one-at-a-time mode finds

using Sparlectra
using Printf

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

const CASEFILE = "case14.m"

# Zeng, Chiang, Neves, Alberto, IJEPES 147 (2023) 108905, example 1: IEEE 14
# with every load at 125 percent, the reactive bands of their table I, bus 8
# as the angle reference and generator 1 dispatched at 306.06 MW. The case
# is built here from the plain case14.m so the example needs no extra file.
const LOAD_FACTOR = 1.25
const PG_BUS1_MW = 306.06
const REFERENCE_BUS = 8
const Q_BANDS_MVAR = Dict(1 => (0.0, 50.0), 2 => (-40.0, 50.0), 3 => (0.0, 40.0), 6 => (-6.0, 24.0), 8 => (-6.0, 300.0))
const REPORT_BUSES = (1, 2, 3, 6, 8)

# active-set settings shared by every run; only the voltage margin differs
const Q_HYST_PU = 0.01
const COOLDOWN_ITERS = 1
const TOL = 1e-8
const MAX_ITER = 40

function build_zeng_case()::Net
  mpc = Sparlectra.MatpowerIO.read_case(ensure_casefile(CASEFILE); legacy_compat = true)
  bus = copy(mpc.bus)
  gen = copy(mpc.gen)
  # MATPOWER columns: bus PD/QD are 3/4 and the type is 2 (1 PQ, 2 PV, 3 ref);
  # gen PG is 2, QMAX/QMIN are 4/5
  bus[:, 3] .*= LOAD_FACTOR
  bus[:, 4] .*= LOAD_FACTOR
  for r in axes(bus, 1)
    b = Int(bus[r, 1])
    b == 1 && (bus[r, 2] = 2.0)
    b == REFERENCE_BUS && (bus[r, 2] = 3.0)
  end
  for r in axes(gen, 1)
    b = Int(gen[r, 1])
    b == 1 && (gen[r, 2] = PG_BUS1_MW)
    qmin, qmax = Q_BANDS_MVAR[b]
    gen[r, 4] = qmax
    gen[r, 5] = qmin
  end
  zeng = Sparlectra.MatpowerIO.MatpowerCase("case14_zeng_p306", mpc.baseMVA, bus, gen, mpc.branch, mpc.gencost, mpc.bus_name, mpc.branch_name, mpc.branch_kind, mpc.for001_contingencies, mpc.dcline, mpc.sparlectra)
  return Sparlectra.createNetFromMatPowerCase(mpc = zeng)
end

# One solve on a fresh copy of the case. Returns what the report compares:
# the converged flag, the iteration count, the switching and release counts,
# the clamped buses, the Q-V rows and the bus voltages and generator Q.
function solve_variant(label::String, mode::Symbol, reenable_v_hyst_pu::Float64; verbose::Int = 0)
  net = build_zeng_case()
  net.q_hyst_pu = Q_HYST_PU
  net.cooldown_iters = COOLDOWN_ITERS
  net.reenable_v_hyst_pu = reenable_v_hyst_pu
  # the release lines are worth seeing, the rest of the verbose solver
  # output is not: capture stdout and keep the PQ->PV lines only
  capture = tempname()
  ite, erg = open(capture, "w") do io
    redirect_stdout(io) do
      runpf!(net, MAX_ITER, TOL, verbose; qlimit_enforcement_mode = mode)
    end
  end
  release_lines = filter(l -> occursin("PQ->PV", l), readlines(capture))
  rm(capture; force = true)
  st = Sparlectra.rectangular_pf_status(net)
  qv = Sparlectra.qvCharacteristicViolations(net)
  vm = Dict(b => Float64(net.nodeVec[b]._vm_pu) for b in REPORT_BUSES)
  qg = machine_q_mvar(net)
  hits = [(iter = e.iter, bus = e.bus, side = e.side) for e in net.qLimitLog]
  return (label = label, mode = mode, converged = erg == 0, iterations = ite,
    reason = st.reason, active_set_ok = st.active_set_converged,
    switches = st.pv_pq_switching_events, releases = st.qlimit_reenable_events,
    clamped = sort!(collect(keys(net.qLimitEvents))), sides = Dict(net.qLimitEvents),
    hits = hits, qv = qv, vm = vm, qg = qg, release_lines = release_lines)
end

# The reactive output of the machine at each reported bus, from the solved
# voltages: the bus injection (V conj(Y V)) plus the load the bus carries.
# The node field `_qƩGen` is not used here because the two enforcement modes
# finalize a clamped bus differently (the clamp value in one, the injection
# in the other), and the table must compare like with like.
function machine_q_mvar(net::Net)::Dict{Int,Float64}
  Y = Sparlectra.createYBUS(net = net, sparse = true, printYBUS = false)
  S = Sparlectra.calc_injections(Y, buildVoltageVector(net))
  return Dict(b => imag(S[b]) * net.baseMVA + Float64(something(net.nodeVec[b]._qƩLoad, 0.0)) for b in REPORT_BUSES)
end

function report_run(r)
  println(r.label)
  status = r.converged ? "converged" : string("not accepted (", r.reason, ")")
  println("  mode: ", r.mode, ", status: ", status, ", iterations: ", r.iterations)
  println("  PV->PQ switching events: ", r.switches, ", PQ->PV release events: ", r.releases, " (one event per iteration with a release)")
  println("  Q-limit hits (iteration, bus, side): ", join((string("(", h.iter, ", ", h.bus, ", ", h.side, ")") for h in r.hits), " "))
  println("  clamped at the end: ", isempty(r.clamped) ? "none" : join((string(b, " (", r.sides[b], ")") for b in r.clamped), ", "))
  significant = filter(row -> row.significant, r.qv)
  if isempty(significant)
    println("  Q-V check: no non-physical generator state")
  else
    println("  Q-V check, machines at a limit with the voltage on the wrong side of the setpoint:")
    for row in significant
      @printf("    bus %2d at Q%s: Vm = %.4f against Vset = %.4f, Q = %.2f MVAr\n", row.bus, row.side, row.vm_pu, row.vset_pu, row.q_MVAr)
    end
  end
  for line in r.release_lines
    println("    ", line)
  end
end

function main()
  print_example_banner("examples/powerflow/example_qlimit_reenable_voltage_rule.jl",
    "voltage-side PQ->PV release of the active-set Q-limit mode on the Zeng/Chiang 14-bus case")

  println("Case: IEEE 14 (Zeng, Chiang, Neves, Alberto 2023, example 1): loads at $(LOAD_FACTOR) times,")
  println("reactive bands of table I, bus $(REFERENCE_BUS) as reference, Pg(bus 1) = $(PG_BUS1_MW) MW.")
  println("Settings: enforcement_mode active_set, q_hyst_pu = $(Q_HYST_PU), cooldown_iters = $(COOLDOWN_ITERS), tol = $(TOL).")
  println()

  # A: the margin is so large that the voltage rule never fires (the state
  # before 0.17.1); B: the default margin; C: the classic mode as control
  a = solve_variant("Run A: active_set, reenable_v_hyst_pu = 1.0 (release switched off)", :active_set, 1.0)
  report_run(a)
  println()
  b = solve_variant("Run B: active_set, reenable_v_hyst_pu = 1e-4 (default, voltage rule active)", :active_set, 1e-4; verbose = 1)
  report_run(b)
  println()
  c = solve_variant("Run C: classic_one_at_a_time (control)", :classic_one_at_a_time, 1e-4)
  report_run(c)
  println()

  println("Bus voltages and generator reactive power at the end of each run:")
  println("  bus │     Vm A │     Vm B │     Vm C │    Qg A [MVAr] │    Qg B [MVAr] │    Qg C [MVAr]")
  for bus in REPORT_BUSES
    @printf("  %3d │ %8.4f │ %8.4f │ %8.4f │ %14.2f │ %14.2f │ %14.2f\n", bus, a.vm[bus], b.vm[bus], c.vm[bus], a.qg[bus], b.qg[bus], c.qg[bus])
  end
  max_dvm = maximum(abs(b.vm[bus] - c.vm[bus]) for bus in REPORT_BUSES)
  @printf("max|dVm| between B and C over the reported buses: %.2e pu (for information: two switching\n", max_dvm)
  println("strategies need not end on the same point, what matters is whether the point is physical)")
  println()
  physical(r) = isempty(filter(row -> row.significant, r.qv))
  println("Physical solution (Q-V check clean): A ", physical(a), ", B ", physical(b), ", C ", physical(c))
  if b.releases >= 1 && physical(b)
    println("Run B released clamped machines and ends on a physical solution; run A kept them clamped and does not.")
  else
    println("Run B did not end on a physical solution; see the release lines and the clamped set above.")
  end
  return nothing
end

run_example(main)
