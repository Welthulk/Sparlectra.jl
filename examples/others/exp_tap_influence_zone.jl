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

# Date: 2026-09-23
# file: examples/others/exp_tap_influence_zone.jl
# purpose: moves one transformer tap by a given number of steps, solves the power flow before and after, and lists the buses whose voltage magnitude or angle changes beyond a threshold (influence zone of the tap)

using Sparlectra
using Printf

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

# Settings
const CASE_NAME = "sp_case14"        # SCF case under data/scf
const TRAFO_BUSES = ("Moorau_110", "Moorau_20")  # terminal buses; nothing = first transformer with a ratio tap
const DELTA_STEPS = 2                # tap move in steps, positive or negative
const VM_TOL_PU = 0.01               # voltage magnitude threshold
const VA_TOL_DEG = 0.2               # voltage angle threshold
const MAX_ITE = 30
const TOL = 1e-8

"""
    tapStep(br) -> Float64

Current tap step of a ratio tap changer on the cascade grid
`tap_ratio = ratio / (1 + n * tap_step)`, with `ratio` the neutral position.
"""
tapStep(br) = (br.ratio / br.tap_ratio - 1.0) / br.tap_step

"""
    tapStepRange(br) -> (nmin, nmax)

Step band of a ratio tap changer, derived from the ratio bounds.
"""
tapStepRange(br) = ((br.ratio / br.tap_max - 1.0) / br.tap_step, (br.ratio / br.tap_min - 1.0) / br.tap_step)

"""
    setTapStep!(br, n)

Sets the live tap ratio of `br` to step `n`; the neutral ratio stays untouched.
"""
function setTapStep!(br, n::Real)
    br.tap_ratio = br.ratio / (1.0 + n * br.tap_step)
    return br
end

# The SCF import regenerates branch component names, so the transformer is
# selected by its terminal bus names (either orientation).
function findTrafo(net, buses)
    pair = buses === nothing ? nothing : Set(geNetBusIdx(net=net, busName=b) for b in buses)
    for br in net.branchVec
        br.has_ratio_tap || continue
        (pair === nothing || Set((br.fromBus, br.toBus)) == pair) && return br
    end
    error("no transformer with a ratio tap changer found (buses = $(buses))")
end

busName(net, idx) = first(k for (k, v) in net.busDict if v == idx)

# hop distance from the transformer terminals over closed branches
function hopDistance(net, br0)
    adj = Dict{Int,Vector{Int}}()
    for br in net.branchVec
        br.status == 1 || continue
        push!(get!(adj, br.fromBus, Int[]), br.toBus)
        push!(get!(adj, br.toBus, Int[]), br.fromBus)
    end
    dist = Dict{Int,Int}(br0.fromBus => 0, br0.toBus => 0)
    queue = [br0.fromBus, br0.toBus]
    while !isempty(queue)
        b = popfirst!(queue)
        for nb in get(adj, b, Int[])
            haskey(dist, nb) && continue
            dist[nb] = dist[b] + 1
            push!(queue, nb)
        end
    end
    return dist
end

wrapDeg(x) = mod(x + 180.0, 360.0) - 180.0

function main()
    print_example_banner("examples/others/exp_tap_influence_zone.jl", "moves one transformer tap by a given number of steps, solves the power flow before and after, and lists the buses whose voltage magnitude or angle changes beyond a threshold (influence zone of the tap)")

    path = joinpath(pkgdir(Sparlectra), "data", "scf", "$(CASE_NAME).scf.json")
    net1 = importSCF(path)
    net2 = deepcopy(net1)

    # state 1: tap as imported
    br1 = findTrafo(net1, TRAFO_BUSES)
    n1 = tapStep(br1)
    nmin, nmax = tapStepRange(br1)

    # state 2: tap moved by DELTA_STEPS, clamped to the step band
    n2 = clamp(n1 + DELTA_STEPS, nmin, nmax)
    n2 == n1 + DELTA_STEPS || @warn "tap move clamped to the step band" requested = n1 + DELTA_STEPS used = n2
    br2 = net2.branchVec[findfirst(b -> b.branchIdx == br1.branchIdx, net2.branchVec)]
    setTapStep!(br2, n2)

    ite1, erg1 = runpf!(net1, MAX_ITE, TOL, 0)
    ite2, erg2 = runpf!(net2, MAX_ITE, TOL, 0)
    (erg1 == 0 && erg2 == 0) || error("power flow did not converge (state 1: $(erg1), state 2: $(erg2))")

    dist = hopDistance(net1, br1)
    rows = NamedTuple[]
    for (k, (a, b)) in enumerate(zip(net1.nodeVec, net2.nodeVec))
        dvm = b._vm_pu - a._vm_pu
        dva = wrapDeg(b._va_deg - a._va_deg)
        push!(rows, (bus=a.comp.cName, hops=get(dist, a.busIdx, -1), vm1=a._vm_pu, dvm=dvm, dva=dva,
            in_zone=abs(dvm) > VM_TOL_PU || abs(dva) > VA_TOL_DEG))
    end
    sort!(rows; by=r -> -abs(r.dvm))

    return (trafo=string(busName(net1, br1.fromBus), " -> ", busName(net1, br1.toBus)), n1=n1, n2=n2, nmin=nmin, nmax=nmax,
        ratio1=br1.tap_ratio, ratio2=br2.tap_ratio, ite1=ite1, ite2=ite2, rows=rows)
end

result = run_example(main)

@printf("Transformer %s: step %.0f -> %.0f (band %.0f..%.0f), ratio %.5f -> %.5f\n",
    result.trafo, result.n1, result.n2, result.nmin, result.nmax, result.ratio1, result.ratio2)
@printf("Power flow: %d / %d iterations; thresholds |dVm| > %.3f pu or |dVa| > %.2f deg\n\n",
    result.ite1, result.ite2, VM_TOL_PU, VA_TOL_DEG)
@printf("%-20s %5s %9s %10s %10s  %s\n", "bus", "hops", "Vm1 [pu]", "dVm [pu]", "dVa [deg]", "zone")
for r in result.rows
    @printf("%-20s %5d %9.4f %+10.5f %+10.4f  %s\n", r.bus, r.hops, r.vm1, r.dvm, r.dva, r.in_zone ? "x" : "")
end
zone = [r.bus for r in result.rows if r.in_zone]
println("\nInfluence zone (", length(zone), " of ", length(result.rows), " buses): ", isempty(zone) ? "none" : join(zone, ", "))