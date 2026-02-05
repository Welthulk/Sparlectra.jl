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

# file: src/limits.jl
# Helper function: ensures that the vector is filled to at least 'bus'
# and fills with a default value (via append!) up to the required length.
function _ensure_bus_index!(v::Vector{T}, bus::Int, default::T) where {T}
  if length(v) < bus
    append!(v, fill(default, bus - length(v)))
  end
  return v
end

"""
    printQLimitLog(net::Net; sort_by=:iter, io::IO=stdout)

Pretty-prints the structured Q-limit log (`net.qLimitLog`) as a small table.

# Keyword arguments
- `sort_by`:  `:iter` (default) or `:bus` — controls sorting.
- `io`: optional output stream (default = `stdout`).

Each line shows:
where `Side` is `:min` or `:max`.
"""
function printQLimitLog(net::Net; sort_by::Symbol = :iter, io::IO = stdout)
  logdata = net.qLimitLog
  if isempty(logdata)
    println(io, "No Q-limit events recorded.")
    return
  end

  # sort by chosen field
  if sort_by == :bus
    sort!(logdata, by = x -> x.bus)
  else
    sort!(logdata, by = x -> x.iter)
  end

  println(io, "─"^34)
  println(io, " Iteration │ Bus │ Side")
  println(io, "─"^34)

  for ev in logdata
    @printf(io, " %9d │ %3d │ %-4s\n", ev.iter, ev.bus, String(ev.side))
  end

  println(io, "─"^34)
  println(io, "Total events: ", length(logdata))
end

Base.@kwdef struct QLimitEvent
  iter::Int
  bus::Int
  side::Symbol   # :min | :max
end

function logQLimitHit!(net::Net, iter::Int, bus::Int, side::Symbol)
  push!(net.qLimitLog, QLimitEvent(iter = iter, bus = bus, side = side))
  net.qLimitEvents[bus] = side
end

"Returns the last iteration number where `bus` hit a Q-limit, or `nothing`."
function lastQLimitIter(net::Net, bus::Int)
  for i = length(net.qLimitLog):-1:1
    if net.qLimitLog[i].bus == bus
      return net.qLimitLog[i].iter
    end
  end
  return nothing
end

function resetQLimitLog!(net::Net)
  empty!(net.qLimitEvents)
  empty!(net.qLimitLog)
  return nothing
end

"""
    getQLimits_pu(net::Net) -> (qmin_pu, qmax_pu)

Return per-bus Q limits in p.u. (build once if empty).
"""
function getQLimits_pu(net::Net)
  if isempty(net.qmin_pu) || isempty(net.qmax_pu)
    buildQLimits!(net)
  end
  _sanitize_q_limits!(net.qmin_pu, net.qmax_pu, length(net.nodeVec))
  return net.qmin_pu, net.qmax_pu
end
"""
    pv_hit_q_limit(net, pv_names)

Gibt `true` zurück, wenn einer der PV-Busse aus `pv_names`
in `net.qLimitEvents` vorkommt.

`pv_names` ist eine Liste von Busnamen (Strings).
"""
function pv_hit_q_limit(net, pv_names)
  pv_idx = map(name -> geNetBusIdx(net = net, busName = name), pv_names)
  return any(haskey(net.qLimitEvents, i) for i in pv_idx)
end



# ------------------------------
# Q-limit helpers (local arrays)
# ------------------------------
@inline function has_q_limits(qmin_pu::AbstractVector, qmax_pu::AbstractVector, i::Int)::Bool
  has_lo = (i <= length(qmin_pu)) && isfinite(qmin_pu[i])
  has_hi = (i <= length(qmax_pu)) && isfinite(qmax_pu[i])
  return has_lo || has_hi
end

@inline function q_limit_band(qmin_pu::AbstractVector, qmax_pu::AbstractVector, i::Int, q_hyst_pu::Float64)::Tuple{Float64,Float64}
  lo = ((i <= length(qmin_pu)) && isfinite(qmin_pu[i])) ? (qmin_pu[i] + q_hyst_pu) : -Inf
  hi = ((i <= length(qmax_pu)) && isfinite(qmax_pu[i])) ? (qmax_pu[i] - q_hyst_pu) :  Inf
  return lo, hi
end

# ------------------------------
# Q-limit sanitizing (drop-in)
# ------------------------------
function _sanitize_q_limits!(qmin_pu::Vector{Float64}, qmax_pu::Vector{Float64}, nb::Int)
  _ensure_bus_index!(qmin_pu, nb, -Inf)
  _ensure_bus_index!(qmax_pu, nb,  Inf)

  @inbounds for i in 1:nb
    # Treat NaN as "no limit"
    if isnan(qmin_pu[i])
      qmin_pu[i] = -Inf
    end
    if isnan(qmax_pu[i])
      qmax_pu[i] =  Inf
    end
  end
  return qmin_pu, qmax_pu
end