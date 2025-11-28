# src/limits.jl
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
