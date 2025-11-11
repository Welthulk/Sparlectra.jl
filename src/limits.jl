# src/limits.jl
# Helper function: ensures that the vector is filled to at least 'bus'
# and fills with a default value (via append!) up to the required length.
function _ensure_bus_index!(v::Vector{T}, bus::Int, default::T) where T
    if length(v) < bus
        append!(v, fill(default, bus - length(v)))
    end
    return v
end

function buildQLimits!(net::Net; reset::Bool=false)
    # Number of buses, if known via nodeVec (otherwise we grow on-demand)
    nbus = length(net.nodeVec)

    # Ensure base size (without field reassignment, only mutations)
    if reset
        empty!(net.qmin_pu); empty!(net.qmax_pu)
        resetQLimitLog!(net)
    end
    if length(net.qmin_pu) < nbus
        append!(net.qmin_pu, fill(Inf,  nbus - length(net.qmin_pu)))
    end
    if length(net.qmax_pu) < nbus
        append!(net.qmax_pu, fill(-Inf, nbus - length(net.qmax_pu)))
    end

    # Aggregation over prosumers/generators
    for ps in net.prosumpsVec
        isGenerator(ps) || continue
        bus = getPosumerBusIndex(ps)

        # Fill on-demand up to bus index (if Bus > nbus)
        _ensure_bus_index!(net.qmin_pu, bus,  Inf)
        _ensure_bus_index!(net.qmax_pu, bus, -Inf)

        qmin_pu = isnothing(ps.minQ) ? -Inf : ps.minQ / net.baseMVA
        qmax_pu = isnothing(ps.maxQ) ?  Inf : ps.maxQ / net.baseMVA

        # Only write if "not yet set" (Sentinel) OR the new value is stricter:
        # - For qmin: smaller = stricter (minimize), Sentinel +Inf => set directly
        cur_qmin = net.qmin_pu[bus]
        net.qmin_pu[bus] = isfinite(cur_qmin) ? min(cur_qmin, qmin_pu) : qmin_pu

        # - For qmax: larger = wider (maximize), Sentinel -Inf => set directly
        cur_qmax = net.qmax_pu[bus]
        net.qmax_pu[bus] = isfinite(cur_qmax) ? max(cur_qmax, qmax_pu) : qmax_pu
    end

    return nothing
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
function printQLimitLog(net::Net; sort_by::Symbol=:iter, io::IO=stdout)
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
    push!(net.qLimitLog, QLimitEvent(iter=iter, bus=bus, side=side))
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