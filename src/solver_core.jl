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

# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 2024-06-01
# Purpose: core functions for power flow solver
# file: src/solver_core.jl

# ------------------------------------------------------------
# Helper: bus currents and complex power injections
# S = V .* conj(I), I = Y * V
# ------------------------------------------------------------

using LinearAlgebra


"""
    calc_currents(Y, V) -> I

Compute bus current injections I = Y * V.
"""
function calc_currents(Y::AbstractMatrix{ComplexF64}, V::AbstractVector{ComplexF64})
    return Y * V
end

"""
    calc_injections(Y, V) -> S

Compute complex bus power injections S = V .* conj(Y*V).
"""
function calc_injections(Y::AbstractMatrix{ComplexF64}, V::AbstractVector{ComplexF64})
    I = calc_currents(Y, V)
    return V .* conj.(I)
end


"""
    solve_linear(A, b; allow_pinv=true)

Solve A*x=b. Falls SingularException und allow_pinv=true: pinv(A)*b.
"""
function solve_linear(A, b; allow_pinv::Bool = true)
    try
        return A \ b
    catch e
        if allow_pinv && (e isa LinearAlgebra.SingularException)
            return pinv(Matrix(A)) * b
        end
        rethrow(e)
    end
end

"""
    non_slack_indices(n, slack_idx) -> Vector{Int}
"""
non_slack_indices(n::Int, slack_idx::Int) = [i for i in 1:n if i != slack_idx]

"""
    build_pos_map(non_slack, n) -> pos::Vector{Int}
pos[i]=k, wenn i==non_slack[k], sonst 0.
"""
function build_pos_map(non_slack::Vector{Int}, n::Int)
    pos = zeros(Int, n)
    @inbounds for (k, b) in enumerate(non_slack)
        pos[b] = k
    end
    return pos
end

"""
    slack_elimination_indices(n, slack_idx)

Return (non_slack, row_idx, col_idx) for [P;Q] and [Vr;Vi] layout of size 2n.
"""
function slack_elimination_indices(n::Int, slack_idx::Int)
    non_slack = non_slack_indices(n, slack_idx)
    idx = vcat(non_slack, n .+ non_slack)
    return non_slack, idx, idx
end


"""
    extract_bus_types_and_vset(net) -> (bus_types::Vector{Symbol}, Vset::Vector{Float64}, slack_idx::Int)

bus_types uses :Slack/:PV/:PQ, Vset is vm setpoint (PV) or current vm default.
"""
function extract_bus_types_and_vset(net::Net)
    n = length(net.nodeVec)
    bus_types = Vector{Symbol}(undef, n)
    Vset      = Vector{Float64}(undef, n)

    slack_idx = 0
    @inbounds for (k, node) in enumerate(net.nodeVec)
        nt = getNodeType(node)
        if nt == Slack
            bus_types[k] = :Slack
            slack_idx = k
        elseif nt == PV
            bus_types[k] = :PV
        elseif nt == PQ
            bus_types[k] = :PQ
        else
            error("extract_bus_types_and_vset: unsupported node type at bus $k")
        end
        Vset[k] = isnothing(node._vm_pu) ? 1.0 : node._vm_pu
    end
    slack_idx == 0 && error("extract_bus_types_and_vset: no slack bus found")
    return bus_types, Vset, slack_idx
end

"""
    build_qload_pu(net) -> Qload_pu::Vector{Float64}

Aggregated reactive load per bus in p.u. (sum of non-generator prosmps qVal / baseMVA).
"""
function build_qload_pu(net::Net)
    nb = length(net.nodeVec)
    Qload = zeros(Float64, nb)
    @inbounds for ps in net.prosumpsVec
        isGenerator(ps) && continue
        bus = getPosumerBusIndex(ps)
        (1 <= bus <= nb) || continue
        q = isnothing(ps.qVal) ? 0.0 : ps.qVal
        Qload[bus] += q / net.baseMVA
    end
    return Qload
end

"""
    build_voltage_vector(busVec) -> V::Vector{ComplexF64}
"""
build_voltage_vector(busVec::Vector{BusData}) =
    ComplexF64[bus.vm_pu * exp(im * bus.va_rad) for bus in busVec]


"""
    compute_sbus_and_totals(Y, V, baseMVA) -> (Sbus_pu, p_MW, q_MVar)
"""
function compute_sbus_and_totals(Y, V, baseMVA::Float64)
    Sbus_pu = calc_injections(Y, V)
    p = sum(real.(Sbus_pu)) * baseMVA
    q = sum(imag.(Sbus_pu)) * baseMVA
    return Sbus_pu, p, q
end
