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

# Date: 2026-09-29
# file: examples/others/exp_tap_sensitivity.jl
# purpose: the sensitivity of every bus voltage to the tap of one
#          transformer, read from the rectangular Newton-Raphson Jacobian
#          at the solution: dV/dtau = -J^-1 * dF/dtau (and the same for the
#          phase shift). One linear solve gives the answer for all buses.
#          From it the reach of the transformer follows: the buses a tap
#          move of a few steps changes by more than a threshold. The
#          prediction is checked against a power flow solved again with
#          the tap moved (the method of exp_tap_influence_zone.jl). Three
#          shipped SCF cases: a transformer on a radial stub (sp_case14),
#          a 380/110 kV coupler into a meshed 110 kV zone (sp_case188) and
#          a 220/110 kV coupler (sp_case60). The console output is compact
#          (checks, reach as rings, voltage levels, walls); the full result
#          (every bus, F_u, sensitivities) is written as a Markdown report
#          to examples/_out/tap_sensitivity/tap_sensitivity.md.
#
# approach
#
#   1. The Jacobian J = dF/dV is the sensitivity of the power equations to
#      the voltages. Its entries between the two buses of a transformer
#      contain the tap: with the from-side complex tap t = tau*exp(j*phi)
#          Y_ff = (y + y0_from) / |t|^2,  Y_ft = -y / conj(t),
#          Y_tf = -y / t,                 Y_tt = y + y0_to.
#      J has no column for the tap: tau is no variable of the Newton
#      system, Sparlectra moves it in an outer loop around the power flow.
#
#   2. How a Jacobian block changes over the Newton iterations (see
#      exp_jacobian_block.jl) is not a sensitivity. It is the second
#      derivative of F times the Newton step, dJ ~ H*dV, and it settles
#      because the step shrinks quadratically.
#
#   3. The sensitivity comes from the implicit function theorem. At the
#      solution F(x, u) = 0 with x the rectangular voltages and u the
#      control (tap_ratio tau, or phase_shift_deg phi); moving u and
#      keeping the equations at zero gives
#          J * dx/du + F_u = 0   =>   dx = -J^-1 * (F_u * du).
#      F_u = dF/du is the partial derivative of the mismatch with the
#      voltages held fixed: `tapMismatchDerivative`. F has the layout
#      [dP_1, dQ|dV_1, dP_2, ...] of the non-slack buses and is
#      F = S_calc - S_spec with S_calc = V .* conj(Y*V), so F_u = dS_calc/du.
#      Only the admittance entries of the transformer carry u, so F_u has
#      at most four entries: dP and dQ at the two transformer buses (a PV
#      bus has |V| - Vset as second equation, which does not depend on u;
#      the slack bus has no rows). With t = tau*exp(j*phi), ys the series
#      and y0f the from-side shunt admittance of the branch:
#          dY_ff/dtau = -2 (ys + y0f) / tau^3
#          dY_ft/dtau =  ys e^(+j phi) / tau^2    dY_tf/dtau = ys e^(-j phi) / tau^2
#          dY_ft/dphi = -j ys e^(+j phi) / tau    dY_tf/dphi = +j ys e^(-j phi) / tau
#      (Y_tt does not depend on u, Y_ff not on phi; phi in radians, the code
#      multiplies by pi/180 for degrees), and
#          dS_f = V_f conj(dY_ff V_f + dY_ft V_t),   dS_t = V_t conj(dY_tf V_f),
#          F_u = [Re dS, Im dS] at the two buses.
#      The formula is the default. Central differences of the mismatch with
#      the admittance matrix rebuilt at u +- h are its check; if the two
#      disagree, the differences are used and the example says so.
#
#   4. The model is what the solver factored. `runpf!` contracts closed
#      bus links (impedance-less couplers, sp_case60 and sp_case188 have
#      one each) into a working net; the Ybus of the uncontracted net
#      would miss the coupler. The example takes that working net from
#      Sparlectra._merged_pf_net (an internal function) and builds the
#      canonical model with buildPfModel: active buses only, Ybus, bus
#      types, setpoints, injections and the solved voltages, with
#      model.busIdx_net mapping a model index back to the bus of the net.
#      The mismatch at the solution is printed; if it is not near zero the
#      model does not reproduce the solved state (Q-limit switching, an
#      unsupported control) and the sensitivities are not to be used.
#
#   5. In rectangular coordinates the unknowns are Vr and Vi, so the result
#      is a complex dV per bus. Magnitude and angle follow from
#          d|V| = Re(conj(V) dV) / |V|,  d(delta) = Im(conj(V) dV) / |V|^2.
#      A PV bus has |V|^2 - Vset^2 as its second equation, which does not
#      depend on tau: its d|V|/dtau is exactly zero, and the bus answers
#      with reactive power. A nonzero value at a PV bus would mean the row
#      layout is mixed up.
#
#   6. Steps. The SCF import and the scenario engine use the grid
#          tap_ratio = ratio / (1 + n * tap_step),
#      `ratio` the neutral position, n the step number. A positive step
#      lowers tap_ratio and raises the voltage of the low side. The tap
#      move of n steps from the present position is
#          dtau = tap_ratio_at(n0 + n) - tap_ratio,
#      and the predicted voltage change is dV/dtau * dtau.
#
#   7. Check: the power flow is solved again on a copy with tap_ratio set
#      to the new value, and the actual change of |V| is put next to the
#      prediction, bus by bus. The zone taken from the prediction (one
#      solve of J) and the zone taken from the second power flow are
#      compared as sets. They differ only by the second-order term, which
#      grows with the square of the tap move.
#
#   8. Reach. J^-1 is dense, so every bus of the island gets a nonzero
#      dV/dtau and the reach is infinite on paper. The question "I move
#      the tap by x steps, which buses change enough to matter" needs a
#      cut: the predicted change against an absolute threshold in pu (a
#      controller deadband, the accuracy of a meter; 0.001 pu here).
#      Below the cut a bus is reached in theory only. The set belongs to
#      one operating point: a PV bus that hits its Q limit and switches to
#      PQ removes a wall and the set grows.
#
#   9. Radial vs. meshed. On a stub the tap sets the voltage of the stub.
#      Inside a mesh, parallel paths hold the two sides together; the tap
#      drives a reactive current around the loop, which lifts |V| on one
#      side and lowers it on the other by amounts the loop impedances
#      decide. d|V|/dtau is then smaller per step, spread over more buses
#      and of both signs. A coupler from the 380 kV (or 220 kV) level into
#      a meshed 110 kV zone is the case in point: the zone hangs on the
#      upper level at several places, so the tap also shifts flow between
#      those infeed points. Voltage-holding buses (slack, PV) are the
#      walls that cut the reach.
#
#  10. Scope. The power flow is `runpf!`, without the outer controllers
#      (tap controllers, STATCOM, SSSC, HVDC setpoint, PST). The
#      sensitivity belongs to that operating point, which is also the one
#      such a controller starts from.

using Sparlectra
using LinearAlgebra
using SparseArrays
using Printf
using Dates

include(joinpath(@__DIR__, "example_header.jl"))

const OUT_DIR = joinpath(dirname(@__DIR__), "_out", "tap_sensitivity")

# Three shipped SCF cases (data/scf). All of them run by default; select one
# with SPARLECTRA_TAP_SENS_SCENARIO=<name> or main(; scenarios = [:sp_case188]).
#
#   sp_case14   the transformer Moorau_110_20 feeds a radial 20 kV stub.
#               The reach is obvious: the stub as a whole.
#   sp_case188  380 kV double ring, four meshed 110 kV zones. The coupler
#               Waldbrok_380_110 (Bornbrok_380 -> Waldbrok_110, 200 MVA)
#               feeds a 110 kV bus whose neighbours are Oberloh_110 and
#               Lindbach_110; Oberloh_110 hangs on Holtau_380 through the
#               PST, so the coupler lies in a loop over both voltage levels.
#   sp_case60   220 kV double ring feeding the meshed 110 kV zone A. The
#               coupler Westdorf_220_110 (140 MVA) feeds Westdorf_110 with
#               the neighbours Falkmoor_110 and Bornloh_110; an SVC holds
#               the Bornloh_110 pocket and should act as a wall.
const SCENARIO_ORDER = [:sp_case14, :sp_case188, :sp_case60]
const SCENARIOS = Dict{Symbol,NamedTuple}(
    :sp_case14 => (trafoBus=("Moorau_110", "Moorau_20"), targetBus="Moorau_20"),
    :sp_case188 => (trafoBus=("Bornbrok_380", "Waldbrok_110"), targetBus="Waldbrok_110"),
    :sp_case60 => (trafoBus=("Westdorf_220", "Westdorf_110"), targetBus="Westdorf_110"),
)

# ---------------------------------------------------------------------------
# tap grid (the one of the SCF import): tap_ratio = ratio / (1 + n * tap_step)
# ---------------------------------------------------------------------------

tap_step_number(br) = (br.ratio / br.tap_ratio - 1.0) / br.tap_step
tap_ratio_at(br, n::Real) = br.ratio / (1.0 + n * br.tap_step)
tap_step_band(br) = ((br.ratio / br.tap_max - 1.0) / br.tap_step, (br.ratio / br.tap_min - 1.0) / br.tap_step)

# The transformer in service between two buses. `reps` maps a bus of the
# uncontracted net onto its representative when the net is the working net.
function find_transformer(net::Net, busA::String, busB::String; reps=nothing)
    a = geNetBusIdx(net=net, busName=busA)
    b = geNetBusIdx(net=net, busName=busB)
    if reps !== nothing
        a = reps[a]
        b = reps[b]
    end
    pair = minmax(a, b)
    for br in net.branchVec
        br.status == 1 && br.has_ratio_tap && minmax(Int(br.fromBus), Int(br.toBus)) == pair && return br
    end
    error("no transformer with a ratio tap in service between $(busA) and $(busB)")
end

# ---------------------------------------------------------------------------
# the model at the solution
# ---------------------------------------------------------------------------

# The state the solver worked on, at the solution: closed bus links
# contracted as in runpf!, then the canonical model of the active buses
# (buildPfModel refreshes the bus types from the prosumers itself).
function solved_state(net::Net)
    wnet, reps, merged = Sparlectra._merged_pf_net(net)
    model = Sparlectra.buildPfModel(wnet; flatstart=false, include_limits=false)
    return (; wnet, reps, merged, model)
end

jacobian_at(m, V) = Sparlectra.build_rectangular_jacobian_pq_pv_sparse(m.Ybus, V, m.busType, m.Vset, m.slack_idx)
mismatch_at(m, V; Ybus=m.Ybus) = Sparlectra.mismatch_rectangular(Ybus, V, m.Sspec, m.busType, m.Vset, m.slack_idx)

# Names of the model buses, in the order of the model.
function bus_labels(net::Net, m)
    names = Dict{Int,String}(idx => bus for (bus, idx) in net.busDict)
    return String[names[i] for i in m.busIdx_net]
end

# The admittance matrix of the working net with one tap field of the
# transformer moved to `value`; the branch is restored afterwards.
function ybus_with(wnet::Net, br, field::Symbol, value::Float64)
    saved = getfield(br, field)
    try
        setfield!(br, field, value)
        return createYBUS(net=wnet, sparse=true)
    finally
        setfield!(br, field, saved)
    end
end

"""
    tapMismatchDerivative(st, br, field; method = :analytic, h = 1e-6) -> Vector{Float64}

`F_u = dF/du`, the derivative of the rectangular mismatch vector `F` with
respect to a tap control `u` of the transformer `br`, the voltages held at
the solution of the state `st` (see `solved_state`). `field` is `:tap_ratio`
(per unit of the ratio) or `:phase_shift_deg` (per degree). The result has
the layout of `F`: `[dP_1, dQ|dV_1, dP_2, ...]` over the non-slack buses of
the model, nonzero only at the two transformer buses.

`method = :analytic` stamps the derivative of the Ybus entries of the
transformer (see the approach comment, point 3); `method = :differences`
takes central differences of the mismatch with the admittance matrix
rebuilt at `u ± h`.
"""
function tapMismatchDerivative(st, br, field::Symbol; method::Symbol=:analytic, h::Float64=1e-6)
    field in (:tap_ratio, :phase_shift_deg) || error("tapMismatchDerivative: field must be :tap_ratio or :phase_shift_deg, got $(field)")
    method === :analytic && return _tap_mismatch_derivative_analytic(st, br, field)
    method === :differences && return _tap_mismatch_derivative_differences(st, br, field, h)
    error("tapMismatchDerivative: method must be :analytic or :differences, got $(method)")
end

function _tap_mismatch_derivative_differences(st, br, field::Symbol, h::Float64)
    m = st.model
    base = getfield(br, field)
    Fup = mismatch_at(m, m.V0; Ybus=ybus_with(st.wnet, br, field, base + h))
    Fdn = mismatch_at(m, m.V0; Ybus=ybus_with(st.wnet, br, field, base - h))
    return (Fup .- Fdn) ./ (2h)
end

function _tap_mismatch_derivative_analytic(st, br, field::Symbol)
    m = st.model
    n = length(m.V0)
    Sparlectra._branch_terminal_state(br) === :closed || error("the transformer is open at a terminal; its tap derivative is not defined")
    kf = findfirst(==(Int(br.fromBus)), m.busIdx_net)
    kt = findfirst(==(Int(br.toBus)), m.busIdx_net)
    (kf === nothing || kt === nothing) && error("a transformer terminal is not an active bus of the model")

    ys = Sparlectra.calcBranchYser(br)
    y0f = Sparlectra._branch_y0_from(br)
    t = Sparlectra.calcBranchRatio(br)      # tau * exp(j*phi)
    tau = abs(t)
    ephi = t / tau
    if field === :tap_ratio
        dYff = -2.0 * (ys + y0f) / tau^3
        dYft = ys * ephi / tau^2
        dYtf = ys * conj(ephi) / tau^2
    else
        c = π / 180                            # per degree
        dYff = 0.0 + 0.0im
        dYft = -im * c * ys * ephi / tau
        dYtf = im * c * ys * conj(ephi) / tau
    end
    Vf = m.V0[kf]
    Vt = m.V0[kt]
    dSf = Vf * conj(dYff * Vf + dYft * Vt)
    dSt = Vt * conj(dYtf * Vf)               # dY_tt/du = 0

    F = zeros(Float64, 2 * (n - 1))
    for (k, dS) in ((kf, dSf), (kt, dSt))
        k == m.slack_idx && continue           # the slack bus has no rows
        pos = k < m.slack_idx ? k : k - 1      # position among the non-slack buses
        F[2*pos-1] += real(dS)             # dP
        m.busType[k] == :PQ && (F[2*pos] += imag(dS))   # dQ; the PV row |V| - Vset does not depend on the tap
    end
    return F
end

"""
    tap_sensitivity(st, br, field; h = 1e-6, rtol = 1e-6)

dV/dtau (`field = :tap_ratio`) or dV/dphi (`field = :phase_shift_deg`, per
degree) for every bus of the model, from one solve with the Jacobian at the
solution:

    J * dx = -F_u

`F_u` comes from `tapMismatchDerivative`: the analytic formula, checked
against central differences; if they differ by more than `rtol` of the
largest entry, the differences are used (`source = :differences`).
Returned per model bus: `dVm` (magnitude, pu per unit of the field) and
`dVa_deg` (angle, degrees per unit of the field), plus `Fu`, `deviation`
(analytic against differences) and `source`. A PV bus keeps its magnitude,
so its `dVm` is zero.
"""
function tap_sensitivity(st, br, field::Symbol; h::Float64=1e-6, rtol::Float64=1e-6)
    m = st.model
    Fu_analytic = tapMismatchDerivative(st, br, field; method=:analytic)
    Fu_differences = tapMismatchDerivative(st, br, field; method=:differences, h=h)
    deviation = maximum(abs, Fu_analytic .- Fu_differences)
    agree = deviation <= rtol * max(maximum(abs, Fu_differences), 1e-12)
    Fu = agree ? Fu_analytic : Fu_differences
    J = jacobian_at(m, m.V0)
    dx = -(J \ Fu)
    n = length(m.V0)
    non_slack = [k for k in 1:n if k != m.slack_idx]
    dV = zeros(ComplexF64, n)
    for (k, bus) in enumerate(non_slack)
        dV[bus] = dx[k] + im * dx[(n-1)+k]
    end
    dVm = [real(conj(m.V0[b]) * dV[b]) / abs(m.V0[b]) for b in 1:n]
    dVa_deg = [rad2deg(imag(conj(m.V0[b]) * dV[b]) / abs2(m.V0[b])) for b in 1:n]
    return (dVm=dVm, dVa_deg=dVa_deg, Fu=Fu, deviation=deviation, source=agree ? :analytic : :differences)
end

# The nonzero entries of F_u with row, equation and bus.
function print_Fu(labels::Vector{String}, m, Fu::Vector{Float64}, unit::String)
    n = length(labels)
    non_slack = [k for k in 1:n if k != m.slack_idx]
    println("F_u = dF/d", unit, ", nonzero entries (row of F, equation, bus number, bus):")
    for r in eachindex(Fu)
        abs(Fu[r]) > 1e-9 || continue
        k = non_slack[(r+1)÷2]
        kind = isodd(r) ? "dP" : (m.busType[k] == :PV ? "dV" : "dQ")
        @printf("  row %3d  %s  %4d  %-14s %+.5f\n", r, kind, m.busIdx_net[k], labels[k], Fu[r])
    end
    return nothing
end

# The power flow solved again on a copy of the net with the tap ratio of the
# transformer set to `tap_ratio`.
function resolve_with(net::Net, busA::String, busB::String, tap_ratio::Float64)
    net2 = deepcopy(net)
    br2 = find_transformer(net2, busA, busB)
    br2.tap_ratio = tap_ratio
    ite, erg = runpf!(net2, 40, 1e-10, 0; method=:rectangular)
    erg == 0 || error("power flow with the moved tap did not converge (erg = $(erg))")
    return (state=solved_state(net2), iterations=ite)
end

"""
    affected_buses(labels, m, predicted, actual; threshold_pu)

The buses sorted by the size of the predicted change of |V| (`predicted`,
from the sensitivity) with the change of the second power flow (`actual`)
next to it, and the two zones: the buses whose predicted change, and whose
actual change, reaches `threshold_pu`. The threshold is the cut that makes
the reach finite: J^-1 is dense, every bus of the island gets a nonzero
value, but below the threshold no controller or meter would notice it. PV
buses hold |V| and never show up, whatever the tap does.
"""
function affected_buses(labels::Vector{String}, m, predicted::Vector{Float64}, actual::Vector{Float64}; threshold_pu::Float64)
    n = length(labels)
    order = sort(collect(1:n); by=k -> -abs(predicted[k]))
    peak = maximum(abs, predicted)
    rows = [(bus=k, nr=m.busIdx_net[k], name=labels[k], type=k == m.slack_idx ? "Slack" : string(m.busType[k]), predicted=predicted[k], actual=actual[k], share=peak > 0 ? abs(predicted[k]) / peak : 0.0) for k in order]
    by_prediction = String[r.name for r in rows if abs(r.predicted) >= threshold_pu]
    by_resolve = String[r.name for r in rows if abs(r.actual) >= threshold_pu]
    return (rows=rows, by_prediction=by_prediction, by_resolve=by_resolve)
end

# The model as the Jacobian sees it: type and position k of every bus.
function print_network(labels::Vector{String}, m, marked)
    n = length(labels)
    width = maximum(length, labels)
    pos(k) = k < m.slack_idx ? k : k - 1
    println("network (number = bus index of the net, k = position of the bus in the Jacobian, * = a transformer bus):")
    for k in 1:n
        label = k == m.slack_idx ? "Slack      " : @sprintf("%-5s k=%-2d ", string(m.busType[k]), pos(k))
        println("  ", labels[k] in marked ? "*" : " ", " ", @sprintf("%4d", m.busIdx_net[k]), "  ", rpad(labels[k], width), "  ", label)
    end
    println("  rows of J: 2k-1 (dP), 2k (dQ or dV); columns: k (Vr), ", n - 1, "+k (Vi)")
    return nothing
end

names_line(names::Vector{String}; limit::Int=12) = join(first(names, limit), ", ") * (length(names) > limit ? ", ..." : "")

"""
    hop_distances(st, br) -> Vector{Int}

Distance of every model bus from the terminals of the transformer `br`, in
branches in service (the terminals are 0, a neighbour of a terminal 1, ...;
-1 for a bus that is not connected to it).
"""
function hop_distances(st, br)
    m = st.model
    n = length(m.V0)
    pos = Dict{Int,Int}(b => k for (k, b) in enumerate(m.busIdx_net))
    adj = [Int[] for _ in 1:n]
    for b in st.wnet.branchVec
        b.status == 1 || continue
        (haskey(pos, Int(b.fromBus)) && haskey(pos, Int(b.toBus))) || continue
        a = pos[Int(b.fromBus)]
        c = pos[Int(b.toBus)]
        a == c && continue
        push!(adj[a], c)
        push!(adj[c], a)
    end
    dist = fill(-1, n)
    queue = Int[pos[Int(br.fromBus)], pos[Int(br.toBus)]]
    for k in queue
        dist[k] = 0
    end
    head = 1
    while head <= length(queue)
        k = queue[head]
        head += 1
        for o in adj[k]
            dist[o] == -1 || continue
            dist[o] = dist[k] + 1
            push!(queue, o)
        end
    end
    return dist
end

# Nominal voltage in kV of every model bus.
bus_levels_kv(st) = Float64[Sparlectra.getNodeVn(st.wnet.nodeVec[i]) for i in st.model.busIdx_net]

# ---------------------------------------------------------------------------
# the reach as data: rings, voltage levels, walls
# ---------------------------------------------------------------------------

# The reach as rings around the transformer: per distance the number of
# buses, how many are in the reach, and the bus (number and name) with the
# largest change.
# With `fold`, all rings beyond the first one without a bus in the reach
# are one row. Returns the rows and the last ring with a bus in the reach.
function ring_rows(labels::Vector{String}, numbers::Vector{Int}, predicted::Vector{Float64}, in_zone::Vector{Bool}, hops::Vector{Int}; fold::Bool=true)
    maxhop = maximum(hops)
    last_in = any(in_zone) ? maximum(hops[k] for k in eachindex(hops) if in_zone[k]) : 0
    show_to = fold ? min(maxhop, last_in + 1) : maxhop
    rows = NamedTuple[]
    for h in 0:show_to
        ks = [k for k in eachindex(hops) if hops[k] == h]
        isempty(ks) && continue
        top = ks[argmax([abs(predicted[k]) for k in ks])]
        push!(rows, (hops=string(h), buses=length(ks), in_reach=count(k -> in_zone[k], ks), top_value=predicted[top], top_name=labels[top], top_nr=numbers[top]))
    end
    rest = [k for k in eachindex(hops) if hops[k] > show_to]
    if !isempty(rest)
        top = rest[argmax([abs(predicted[k]) for k in rest])]
        push!(rows, (hops=string(">", show_to), buses=length(rest), in_reach=count(k -> in_zone[k], rest), top_value=predicted[top], top_name=labels[top], top_nr=numbers[top]))
    end
    return rows, last_in
end

level_rows(levels::Vector{Float64}, in_zone::Vector{Bool}) = [(kv=kv, in_reach=count(k -> in_zone[k], findall(==(kv), levels)), buses=count(==(kv), levels)) for kv in sort(unique(levels); rev=true)]

# Slack and PV buses connected to the transformer, nearest first: (hops, name, bus number).
wall_rows(labels::Vector{String}, m, hops::Vector{Int}) = sort([(hops[k], labels[k], m.busIdx_net[k]) for k in eachindex(hops) if m.busType[k] != :PQ && hops[k] >= 0])

function print_rings(rows, lonely::Int)
    println("  by distance from the transformer (hops over the branches in service):")
    println("    hops  buses  in reach  largest d|V| pu     nr  bus")
    for r in rows
        @printf("    %4s  %5d  %8d  %+.5f  %5d  %s\n", r.hops, r.buses, r.in_reach, r.top_value, r.top_nr, r.top_name)
    end
    lonely > 0 && println("    not connected to the transformer: ", lonely, " buses")
    return nothing
end

function print_levels(rows)
    println("  by voltage level:")
    for r in rows
        @printf("    %6.0f kV  %3d of %-3d buses in reach\n", r.kv, r.in_reach, r.buses)
    end
    return nothing
end

function print_walls(walls)
    if isempty(walls)
        println("  walls: no slack or PV bus connected to the transformer")
    else
        println("  walls, slack and PV buses that hold |V| (hops): ", names_line(String["$(w[2]) ($(w[1]))" for w in walls]; limit=8))
    end
    return nothing
end

# ---------------------------------------------------------------------------
# the Markdown report
# ---------------------------------------------------------------------------

# One scenario as a Markdown section; `d` is the NamedTuple run_scenario
# assembles. The pipe character is avoided in the table headers (dVm and dVa
# stand for the change of the voltage magnitude and angle).
function scenario_markdown(d)
    io = IOBuffer()
    println(io, "## ", d.name, ": ", d.trafo)
    println(io)
    println(io, "| item | value |")
    println(io, "|---|---|")
    @printf(io, "| tap move | %+.0f step(s), tap_ratio %.5f to %.5f |\n", d.moved, d.tap_old, d.tap_new)
    @printf(io, "| buses in the model | %d%s |\n", d.n, d.nmerged > 0 ? @sprintf(" (%d merged by closed bus links)", d.nmerged) : "")
    @printf(io, "| mismatch at the solution | %.2e |\n", d.res0)
    @printf(io, "| F_u, analytic against central differences | %.2e |\n", d.Fu_deviation)
    @printf(io, "| %s, dVm | predicted %+.5f pu, second power flow %+.5f pu |\n", d.target, d.pred_t, d.act_t)
    @printf(io, "| largest deviation over all buses | %.2e pu at %s |\n", d.worst_dev, d.worst_name)
    if abs(d.per_step) > 1e-12
        @printf(io, "| one step at %s | %+.5f pu; %+.4f pu needs about %+.1f steps |\n", d.target, d.per_step, d.dVm_wanted, d.dVm_wanted / d.per_step)
    end
    @printf(io, "| reach at %.4f pu | %d of %d buses (%.0f %%); second power flow %s |\n", d.threshold_pu, d.zone_pred, d.n, 100 * d.zone_pred / d.n, d.same_set ? "finds the same set" : "differs")
    println(io)

    println(io, "### Reach by distance from the transformer")
    println(io)
    println(io, "| hops | buses | in reach | largest dVm pu | nr | bus |")
    println(io, "|---:|---:|---:|---:|---:|---|")
    for r in d.rings
        @printf(io, "| %s | %d | %d | %+.5f | %d | %s |\n", r.hops, r.buses, r.in_reach, r.top_value, r.top_nr, r.top_name)
    end
    println(io)

    println(io, "### Reach by voltage level")
    println(io)
    println(io, "| kV | buses | in reach |")
    println(io, "|---:|---:|---:|")
    for r in d.level_rows
        @printf(io, "| %.0f | %d | %d |\n", r.kv, r.buses, r.in_reach)
    end
    println(io)

    println(io, "### Walls (slack and PV buses that hold the voltage magnitude)")
    println(io)
    if isempty(d.walls)
        println(io, "No slack or PV bus is connected to the transformer.")
    else
        println(io, join(String["`$(w[2])` (nr $(w[3]), $(w[1]) hops)" for w in d.walls], ", "))
    end
    println(io)

    println(io, "### F_u = dF/dtap_ratio, nonzero entries")
    println(io)
    println(io, "| row of F | equation | nr | bus | value |")
    println(io, "|---:|---|---:|---|---:|")
    non_slack = [k for k in eachindex(d.labels) if k != d.m.slack_idx]
    for r in eachindex(d.Fu)
        abs(d.Fu[r]) > 1e-9 || continue
        k = non_slack[(r+1)÷2]
        kind = isodd(r) ? "dP" : (d.m.busType[k] == :PV ? "dV" : "dQ")
        @printf(io, "| %d | %s | %d | %s | %+.5f |\n", r, kind, d.m.busIdx_net[k], d.labels[k], d.Fu[r])
    end
    println(io)

    in_reach = [r for r in d.reach.rows if abs(r.predicted) >= d.threshold_pu]
    println(io, "### Buses in the reach (", length(in_reach), ")")
    println(io)
    println(io, "| nr | bus | kV | type | hops | predicted dVm | second PF dVm | share of max |")
    println(io, "|---:|---|---:|---|---:|---:|---:|---:|")
    for r in in_reach
        @printf(io, "| %d | %s | %.0f | %s | %d | %+.5f | %+.5f | %.1f %% |\n", r.nr, r.name, d.levels[r.bus], r.type, d.hops[r.bus], r.predicted, r.actual, 100 * r.share)
    end
    println(io)

    println(io, "<details><summary>All ", length(d.reach.rows), " buses with the sensitivities</summary>")
    println(io)
    println(io, "| nr | bus | kV | type | hops | Vm pu | dVm/dtau | dVa/dtau deg | dVm/dphi | dVa/dphi deg | predicted dVm | second PF dVm |")
    println(io, "|---:|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in d.reach.rows
        k = r.bus
        @printf(io, "| %d | %s | %.0f | %s | %d | %.5f | %.4f | %.4f | %.5f | %.5f | %+.5f | %+.5f |\n", r.nr, r.name, d.levels[k], r.type, d.hops[k], abs(d.m.V0[k]), d.s_tau.dVm[k], d.s_tau.dVa_deg[k], d.s_phi.dVm[k], d.s_phi.dVa_deg[k], r.predicted, r.actual)
    end
    println(io)
    println(io, "</details>")
    println(io)
    return String(take!(io))
end

# The report: summary table, a short legend, then one section per scenario.
function write_report(results, sections::Vector{String}, tap_steps::Int, threshold_pu::Float64)
    mkpath(OUT_DIR)
    path = joinpath(OUT_DIR, "tap_sensitivity.md")
    io = IOBuffer()
    println(io, "# Tap sensitivity: reach of a transformer")
    println(io)
    @printf(io, "Written by `examples/others/exp_tap_sensitivity.jl` on %s. Tap move: %d step(s) per scenario, threshold %.4f pu.\n\n", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM"), tap_steps, threshold_pu)
    println(io, "## Summary")
    println(io)
    println(io, "| scenario | transformer | buses | reach J | reach 2nd PF | share | hops | wall | max dVm pu | deviation pu |")
    println(io, "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in results
        @printf(io, "| %s | %s | %d | %d | %d | %.0f %% | %d | %s | %.5f | %.2e |\n", r.scenario, r.trafo, r.buses, r.zone, r.zone_resolve, 100 * r.zone / r.buses, r.last_hop, r.wall_hop < 0 ? "-" : string(r.wall_hop), r.max_dvm, r.deviation)
    end
    println(io)
    println(io, "- **reach J**: buses whose voltage magnitude changes by at least the threshold, predicted from one solve with the Jacobian (dx = -J \\ (F_u * du)).")
    println(io, "- **reach 2nd PF**: the same, taken from a second power flow with the tap moved.")
    println(io, "- **hops**: last ring around the transformer (in branches) that still has a bus in the reach.")
    println(io, "- **wall**: distance of the nearest slack or PV bus, which holds the voltage magnitude and cuts the reach.")
    println(io, "- **deviation**: largest difference between prediction and second power flow, the second-order term.")
    println(io, "- dVm is the change of the voltage magnitude in pu, dVa of the angle in degrees; tau is tap_ratio, phi is phase_shift_deg.")
    println(io, "- nr is the bus index of the net (`net.busDict`); a bus merged into a neighbour by a closed bus link appears under the number of its representative.")
    println(io)
    for sec in sections
        print(io, sec)
    end
    write(path, String(take!(io)))
    return path
end

"""
    run_scenario(name; tap_steps = 2, threshold_pu = 0.001, dVm_wanted = 0.01, verbose = false)

One scenario of `SCENARIOS`: solves the SCF case, builds the model at the
solution, computes dV/dtau and dV/dphi of the transformer, checks the
prediction of a tap move of `tap_steps` steps against a second power flow,
and prints the reach (buses at or above `threshold_pu`) as rings around the
transformer, per voltage level, and with the walls. `verbose` adds the
details to the console (network, F_u, sensitivity table, one line per bus).
Returns a NamedTuple with the numbers of the summary and the Markdown
section (`md`), or `nothing` if the model does not reproduce the solution.
"""
function run_scenario(name::Symbol; tap_steps::Int=2, threshold_pu::Float64=0.001, dVm_wanted::Float64=0.01, verbose::Bool=false)
    sc = SCENARIOS[name]
    casefile = joinpath(pkgdir(Sparlectra), "data", "scf", string(name, ".scf.json"))
    println("-"^80)
    net = Sparlectra.importSCF(casefile)
    ite, erg = runpf!(net, 40, 1e-10, 0; method=:rectangular)
    erg == 0 || error("power flow did not converge (erg = $(erg))")

    st = solved_state(net)
    m = st.model
    n = length(m.V0)
    labels = bus_labels(st.wnet, m)
    nmerged = count(i -> st.reps[i] != i, eachindex(st.reps))
    res0 = maximum(abs, mismatch_at(m, m.V0))
    if res0 > 1e-6
        @printf("%s: the model does not reproduce the solved state (max mismatch %.2e); no sensitivities for this case (Q-limit switching or a control the model does not carry)\n\n", name, res0)
        return nothing
    end

    br = find_transformer(st.wnet, sc.trafoBus...; reps=st.reps)
    br.tap_step > 0.0 || error("transformer $(sc.trafoBus) has no ratio tap step")
    n0 = tap_step_number(br)
    nmin, nmax = tap_step_band(br)
    n1 = clamp(n0 + tap_steps, nmin, nmax)
    moved = n1 - n0
    isapprox(moved, tap_steps; atol=1e-6) || @warn "tap move clamped to the step band" requested = tap_steps used = moved
    tap_new = tap_ratio_at(br, n1)
    dtau = tap_new - br.tap_ratio
    dtau1 = tap_ratio_at(br, n0 + 1) - br.tap_ratio

    s_tau = tap_sensitivity(st, br, :tap_ratio)
    s_phi = tap_sensitivity(st, br, :phase_shift_deg)

    # check: the tap is moved, the power flow is solved again
    re = resolve_with(net, sc.trafoBus[1], sc.trafoBus[2], tap_new)
    m2 = re.state.model
    m2.busIdx_net == m.busIdx_net || error("the active buses differ between the two power flows")
    actual = abs.(m2.V0) .- abs.(m.V0)
    predicted = s_tau.dVm .* dtau
    worst = argmax(abs.(actual .- predicted))
    width = max(4, maximum(length, labels))
    t = findfirst(==(st.reps[geNetBusIdx(net=net, busName=sc.targetBus)]), m.busIdx_net)
    t === nothing && error("target bus $(sc.targetBus) is not an active bus of the model")
    reach = affected_buses(labels, m, predicted, actual; threshold_pu=threshold_pu)
    in_zone = Bool[abs(predicted[k]) >= threshold_pu for k in 1:n]
    hops = hop_distances(st, br)
    levels = bus_levels_kv(st)
    per_step = s_tau.dVm[t] * dtau1
    Fu_deviation = max(s_tau.deviation, s_phi.deviation)
    same_set = isempty(symdiff(Set(reach.by_prediction), Set(reach.by_resolve)))
    rings, last_hop = ring_rows(labels, m.busIdx_net, predicted, in_zone, hops; fold=true)
    levelrows = level_rows(levels, in_zone)
    walls = wall_rows(labels, m, hops)
    wall_hop = isempty(walls) ? -1 : walls[1][1]

    # ---- compact result on the console ----------------------------------------
    @printf("%s: %s -T- %s, %+.0f step(s) (tap_ratio %.5f -> %.5f), %d buses in the model", name, sc.trafoBus[1], sc.trafoBus[2], moved, br.tap_ratio, tap_new, n)
    println(nmerged > 0 ? @sprintf(" (%d merged by closed bus links)", nmerged) : "")
    @printf("  checks: mismatch at the solution %.2e, F_u analytic against differences %.2e\n", res0, Fu_deviation)
    @printf("          %s d|V| predicted %+.5f, second power flow %+.5f pu, largest deviation over all buses %.2e pu\n", sc.targetBus, predicted[t], actual[t], abs(actual[worst] - predicted[worst]))
    (s_tau.source === :analytic && s_phi.source === :analytic) || println("          the analytic F_u disagrees with the differences; the differences are used")
    m2.busType == m.busType || println("          a bus changed its type in the second power flow, the linear prediction does not cover that")
    if abs(per_step) > 1e-12
        @printf("  one step moves %s by %+.5f pu; %+.4f pu there needs about %+.1f steps\n", sc.targetBus, per_step, dVm_wanted, dVm_wanted / per_step)
    else
        println("  ", sc.targetBus, " does not react to the tap (a PV bus, or not reached through this transformer)")
    end
    println()
    @printf("  reach at %.4f pu: %d of %d buses (%.0f %%)\n", threshold_pu, count(in_zone), n, 100 * count(in_zone) / n)
    println("  ", same_set ? "the second power flow finds the same set" : "the second power flow differs in: " * names_line(sort!(collect(symdiff(Set(reach.by_prediction), Set(reach.by_resolve))))))
    print_rings(rings, count(==(-1), hops))
    print_levels(levelrows)
    print_walls(walls)
    println()

    # ---- details on the console -----------------------------------------------
    if verbose
        small = n <= 30   # bus lists are printed in full for small cases only
        println("  details (solved in ", ite, " iteration(s))")
        if small
            print_network(labels, m, collect(sc.trafoBus))
            println()
        end
        @printf("one step = %+.5f in tap_ratio (a positive step lowers tap_ratio and raises the low side); step %+.0f of %+.0f..%+.0f, tap_step = %.5f, phase_shift = %.3f deg\n", dtau1, n0, nmin, nmax, br.tap_step, br.phase_shift_deg)
        print_Fu(labels, m, s_tau.Fu, "tap_ratio")
        println()
        println("sensitivity of the bus voltages, from one solve with J at the solution")
        println("    nr  bus          type    |V| pu    d|V|/dtau  ddelta/dtau   d|V|/dphi  ddelta/dphi")
        println("                                           pu/1     deg/1        pu/deg      deg/deg")
        shown = small ? collect(1:n) : sort(sort(collect(1:n); by=k -> -abs(s_tau.dVm[k]))[1:min(15, n)])
        small || println("  (the 15 buses with the largest d|V|/dtau of ", n, ")")
        for k in shown
            typ = k == m.slack_idx ? "Slack" : string(m.busType[k])
            @printf("  %4d  %s %-5s  %8.5f  %10.4f  %10.4f   %10.5f  %10.5f\n", m.busIdx_net[k], rpad(labels[k], width), typ, abs(m.V0[k]), s_tau.dVm[k], s_tau.dVa_deg[k], s_phi.dVm[k], s_phi.dVa_deg[k])
        end
        println()
        println("buses by the predicted change of |V| (* by prediction, + by second power flow)")
        println("        nr  bus          type   predicted     actual   share of max")
        zone_size = max(length(reach.by_prediction), length(reach.by_resolve))
        listed = small ? length(reach.rows) : min(length(reach.rows), max(10, zone_size + 3))
        for r in reach.rows[1:listed]
            mark = string(abs(r.predicted) >= threshold_pu ? "*" : " ", abs(r.actual) >= threshold_pu ? "+" : " ")
            bar = repeat("#", round(Int, 20 * r.share))
            @printf("  %s %4d  %s %-5s  %+.5f  %+.5f   %5.1f %%  %s\n", mark, r.nr, rpad(r.name, width), r.type, r.predicted, r.actual, 100 * r.share, bar)
        end
        listed < length(reach.rows) && println("  ... ", length(reach.rows) - listed, " more buses below")
        println()
    end

    # ---- Markdown section (always the full detail) ------------------------------
    trafo = string(sc.trafoBus[1], " -T- ", sc.trafoBus[2])
    md = scenario_markdown((; name, trafo, moved, tap_old=br.tap_ratio, tap_new, n, nmerged, res0, Fu_deviation, target=sc.targetBus, pred_t=predicted[t], act_t=actual[t], worst_dev=abs(actual[worst] - predicted[worst]), worst_name=labels[worst], per_step, dVm_wanted, threshold_pu, zone_pred=length(reach.by_prediction), same_set, rings=ring_rows(labels, m.busIdx_net, predicted, in_zone, hops; fold=false)[1], level_rows=levelrows, walls, Fu=s_tau.Fu, labels, m, levels, hops, s_tau, s_phi, reach))

    return (scenario=name, trafo=trafo, buses=n, moved=moved, zone=length(reach.by_prediction), zone_resolve=length(reach.by_resolve), last_hop=last_hop, wall_hop=wall_hop, max_dvm=maximum(abs, predicted), deviation=abs(actual[worst] - predicted[worst]), md=md)
end

function print_summary(results, threshold_pu::Float64)
    println("="^80)
    @printf("summary: the same tap move in every scenario, threshold %.4f pu\n", threshold_pu)
    println("  scenario     transformer                    buses  reach J|PF   share  hops  wall  max d|V| pu  deviation pu")
    for r in results
        wall = r.wall_hop < 0 ? "-" : string(r.wall_hop)
        @printf("  %-11s  %-28s %6d  %4d|%-4d  %5.0f%%  %4d  %4s  %11.5f  %12.2e\n", string(r.scenario), r.trafo, r.buses, r.zone, r.zone_resolve, 100 * r.zone / r.buses, r.last_hop, wall, r.max_dvm, r.deviation)
    end
    println()
    println("Reading the numbers:")
    println("  reach   buses whose |V| changes by at least the threshold (J = predicted from one solve, PF = second power flow)")
    println("  hops    last ring around the transformer that still has a bus in the reach")
    println("  wall    distance of the nearest slack or PV bus, which holds |V| and cuts the reach")
    println("  deviation   largest difference between prediction and second power flow, the second-order term")
    return nothing
end

# All scenarios by default; SPARLECTRA_TAP_SENS_SCENARIO=<name> selects one.
function selected_scenarios()
    choice = get(ENV, "SPARLECTRA_TAP_SENS_SCENARIO", "all")
    choice == "all" && return copy(SCENARIO_ORDER)
    haskey(SCENARIOS, Symbol(choice)) || error("unknown scenario $(choice); one of $(join(SCENARIO_ORDER, ", ")) or all")
    return [Symbol(choice)]
end

"""
    main(; scenarios, tap_steps, threshold_pu, dVm_wanted, verbose)

Runs every scenario of `scenarios` (default: all three, or the one named in
SPARLECTRA_TAP_SENS_SCENARIO), prints a compact result per scenario and a
summary table, and writes the full result as a Markdown report to
`examples/_out/tap_sensitivity/tap_sensitivity.md`. `verbose` (or
SPARLECTRA_TAP_SENS_VERBOSE=1) adds the details to the console. A scenario
that fails prints its error and the others still run. Returns `nothing`, so
a run from the REPL does not echo a result tuple; call `run_scenario`
directly for the numbers.
"""
function main(; scenarios::Vector{Symbol}=selected_scenarios(), tap_steps::Int=2, threshold_pu::Float64=0.001, dVm_wanted::Float64=0.01, verbose::Bool=get(ENV, "SPARLECTRA_TAP_SENS_VERBOSE", "0") == "1")
    print_example_banner("examples/others/exp_tap_sensitivity.jl", "the sensitivity of the bus voltages to a transformer tap, from the Jacobian at the solution")
    results = NamedTuple[]
    sections = String[]
    for name in scenarios
        try
            r = run_scenario(name; tap_steps=tap_steps, threshold_pu=threshold_pu, dVm_wanted=dVm_wanted, verbose=verbose)
            if r === nothing
                push!(sections, string("## ", name, "\n\nThe model does not reproduce the solved state; no sensitivities for this case.\n\n"))
            else
                push!(results, r)
                push!(sections, r.md)
            end
        catch err
            println("scenario ", name, " failed:")
            showerror(stdout, err, catch_backtrace())
            println()
            println()
            push!(sections, string("## ", name, "\n\nFailed: `", sprint(showerror, err), "`\n\n"))
        end
    end
    print_summary(results, threshold_pu)
    path = write_report(results, sections, tap_steps, threshold_pu)
    println()
    println("report: ", relpath(path, dirname(dirname(@__DIR__))))
    return nothing
end

run_example(main)