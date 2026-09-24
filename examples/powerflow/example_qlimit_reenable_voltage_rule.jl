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
#          mode on the Zeng/Chiang 14-bus case against the two classic modes:
#          with the release switched off the run ends on a non-physical
#          solution (three machines at Qmax with the voltage above the
#          setpoint), with the rule it lands on the physical solution that
#          classic one-at-a-time finds; classic simultaneous clamps every
#          violation at once and ends on the non-physical point as well.
#          One table per run with the switching path, the end state and the
#          physical verdict per bus, then the comparison and the solver times.

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

# One solve on a fresh copy of the case. The solver runs with verbose = 1
# into a capture file; its `PQ->PV Bus n ... released` lines are the record
# of the releases, its warnings are shown under the run. Returns what the
# table shows: status, iteration count, the hit log, the released buses, the
# clamped set, the Q-V rows, and per bus the setpoint, voltage, machine Q
# and band.
function solve_variant(label::String, mode::Symbol, reenable_v_hyst_pu::Float64)
  net = build_zeng_case()
  net.q_hyst_pu = Q_HYST_PU
  net.cooldown_iters = COOLDOWN_ITERS
  net.reenable_v_hyst_pu = reenable_v_hyst_pu
  capture = tempname()
  seconds = 0.0
  ite, erg = open(capture, "w") do io
    redirect_stdio(stdout = io, stderr = io) do
      # the whole solve including the switching logic; the file capture adds
      # nothing measurable, the compile time is paid by the warm-up in main
      result = nothing
      seconds = @elapsed result = runpf!(net, MAX_ITER, TOL, 1; qlimit_enforcement_mode = mode)
      result
    end
  end
  captured = readlines(capture)
  rm(capture; force = true)
  release_lines = filter(l -> occursin("PQ->PV", l), captured)
  # the solver's own warnings (a run it does not accept) belong under the
  # run they concern, not in front of its heading
  notes = [strip(replace(l, r"^.*Warning: " => "")) for l in captured if occursin("Warning:", l)]
  released = Int[]
  for line in release_lines
    m = match(r"PQ->PV Bus (\d+)", line)
    m === nothing || push!(released, parse(Int, m.captures[1]))
  end
  st = Sparlectra.rectangular_pf_status(net)
  qv = Sparlectra.qvCharacteristicViolations(net)
  qmin_pu, qmax_pu = Sparlectra.getQLimits_pu(net)
  vm = Dict(b => Float64(net.nodeVec[b]._vm_pu) for b in REPORT_BUSES)
  vset = Dict(b => generator_setpoint(net, b) for b in REPORT_BUSES)
  band = Dict(b => (qmin_pu[b] * net.baseMVA, qmax_pu[b] * net.baseMVA) for b in REPORT_BUSES)
  qg = machine_q_mvar(net)
  hits = [(iter = e.iter, bus = e.bus, side = e.side) for e in net.qLimitLog]
  return (label = label, mode = mode, converged = erg == 0, iterations = ite, reason = st.reason, seconds = seconds,
    released = unique(released), clamped = Dict(net.qLimitEvents), hits = hits, qv = qv, notes = notes,
    vm = vm, vset = vset, band = band, qg = qg)
end

# the voltage setpoint of the machine at `bus` (NaN when no machine regulates there)
function generator_setpoint(net::Net, bus::Int)::Float64
  for ps in net.prosumpsVec
    isGenerator(ps) && getPosumerBusIndex(ps) == bus && ps.vm_pu !== nothing && return Float64(ps.vm_pu)
  end
  return NaN
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

# end state and history of one bus: what the hit log and the release lines
# say, in words a reader can check against the table
function bus_story(r, bus::Int)
  bus == REFERENCE_BUS && return ("slack", "")
  side = get(r.clamped, bus, nothing)
  was_released = bus in r.released
  state = side === nothing ? (was_released ? "PV, released" : "PV") : string("PQ at Q", side)
  events = String[]
  for (k, h) in enumerate(filter(h -> h.bus == bus, r.hits))
    push!(events, string(k == 1 ? "hit Q" : "hit again Q", h.side, " it ", h.iter))
    k == 1 && was_released && push!(events, "released back to PV")
  end
  return (state, join(events, ", "))
end

# the Q-V verdict of one bus: a machine at a reactive limit with the voltage
# on the wrong side of its setpoint is a non-physical state (it could not
# hold that voltage at that limit); everything else is physical
function physical_verdict(r, bus::Int)::String
  for row in r.qv
    row.bus == bus && row.significant && return string("NO, Vm ", row.side == :max ? "above" : "below", " Vset at Q", row.side)
  end
  return "yes"
end

# the switching path of one bus in words: none, PV->PQ (clamped and kept),
# PV->PQ->PV (released back), PV->PQ->PV->PQ (released and clamped again)
function switching_path(r, bus::Int)::String
  hits = filter(h -> h.bus == bus, r.hits)
  isempty(hits) && return "none"
  parts = ["PV"]
  for (k, _) in enumerate(hits)
    push!(parts, "PQ")
    k == 1 && bus in r.released && push!(parts, "PV")
  end
  return join(parts, "->")
end

function report_run(r)
  println(r.label)
  status = r.converged ? "converged" : string("not accepted (", r.reason, ")")
  @printf("  mode: %s, status: %s, iterations: %d, solver time: %.4f s\n", r.mode, status, r.iterations, r.seconds)
  for note in r.notes
    println("  solver note: ", note)
  end
  println("  bus │   Vset │     Vm │  Qg [MVAr] │  band [MVAr] │ switching      │ end state     │ physical                     │ history")
  for bus in REPORT_BUSES
    state, story = bus_story(r, bus)
    lo, hi = r.band[bus]
    @printf("  %3d │ %6.4f │ %6.4f │ %10.2f │ [%4.0f, %4.0f] │ %-14s │ %-13s │ %-28s │ %s\n", bus, r.vset[bus], r.vm[bus], r.qg[bus], lo, hi, switching_path(r, bus), state, physical_verdict(r, bus), story)
  end
end

function main()
  print_example_banner("examples/powerflow/example_qlimit_reenable_voltage_rule.jl",
    "voltage-side PQ->PV release of the active-set Q-limit mode on the Zeng/Chiang 14-bus case")

  println("Case: IEEE 14 (Zeng, Chiang, Neves, Alberto 2023, example 1): loads at $(LOAD_FACTOR) times,")
  println("reactive bands of table I, bus $(REFERENCE_BUS) as reference, Pg(bus 1) = $(PG_BUS1_MW) MW.")
  println("Settings: enforcement_mode active_set, q_hyst_pu = $(Q_HYST_PU), cooldown_iters = $(COOLDOWN_ITERS), tol = $(TOL).")
  println()

  # one throwaway solve per enforcement mode pays the compile time (the
  # classic outer loop is its own code path), so the times below are solver time
  solve_variant("warm-up", :active_set, 1e-4)
  solve_variant("warm-up", :classic_one_at_a_time, 1e-4)
  solve_variant("warm-up", :classic_simultaneous, 1e-4)
  # A: the margin is so large that the voltage rule never fires (the state
  # before 0.17.1); B: the default margin; C and D: the classic modes as control
  a = solve_variant("Run A: active_set, reenable_v_hyst_pu = 1.0 (release switched off)", :active_set, 1.0)
  report_run(a)
  println()
  b = solve_variant("Run B: active_set, reenable_v_hyst_pu = 1e-4 (default, voltage rule active)", :active_set, 1e-4)
  report_run(b)
  println()
  c = solve_variant("Run C: classic_one_at_a_time (control, the largest violation per outer pass)", :classic_one_at_a_time, 1e-4)
  report_run(c)
  println()
  d = solve_variant("Run D: classic_simultaneous (control, all violations at once)", :classic_simultaneous, 1e-4)
  report_run(d)
  println()

  println("Bus voltages and machine reactive power at the end of each run:")
  println("  bus │   Vset │   Vm A │   Vm B │   Vm C │   Vm D │ Qg A [MVAr] │ Qg B [MVAr] │ Qg C [MVAr] │ Qg D [MVAr] │ phys A │ phys B │ phys C │ phys D")
  short(r, bus) = physical_verdict(r, bus) == "yes" ? "yes" : "NO"
  for bus in REPORT_BUSES
    @printf("  %3d │ %6.4f │ %6.4f │ %6.4f │ %6.4f │ %6.4f │ %11.2f │ %11.2f │ %11.2f │ %11.2f │ %-6s │ %-6s │ %-6s │ %-6s\n", bus, b.vset[bus], a.vm[bus], b.vm[bus], c.vm[bus], d.vm[bus], a.qg[bus], b.qg[bus], c.qg[bus], d.qg[bus], short(a, bus), short(b, bus), short(c, bus), short(d, bus))
  end
  println()
  println("Runs (solver time is the whole solve including the switching logic, after a warm-up):")
  println("  run │ mode                  │ status                                          │ iterations │ solver time [s]")
  for (name, r) in (("A", a), ("B", b), ("C", c), ("D", d))
    status = r.converged ? "converged" : string("not accepted (", r.reason, ")")
    @printf("  %3s │ %-21s │ %-47s │ %10d │ %15.4f\n", name, r.mode, status, r.iterations, r.seconds)
  end
  return nothing
end

run_example(main)
