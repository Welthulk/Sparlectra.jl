# src/limits.jl

"""
    _prosumer_bus_index(ps::ProSumer) -> Int

Internal helper: derive the prosumer's bus index from its component.
Prefers `cFrom_bus`, then falls back to `cTo_bus`. Throws a clear error
if neither is available.
"""
function _prosumer_bus_index(ps::ProSumer)::Int
    c = ps.comp
    if hasproperty(c, :cFrom_bus) && getfield(c, :cFrom_bus) !== nothing
        return Int(getfield(c, :cFrom_bus))
    elseif hasproperty(c, :cTo_bus) && getfield(c, :cTo_bus) !== nothing
        return Int(getfield(c, :cTo_bus))
    end
    error("ProSumer: cannot determine bus index (component has neither :cFrom_bus nor :cTo_bus).")
end

# ---------- Storage helpers (immutable Net: never rebind fields) ----------

"""
    _ensure_limit_storage!(net::Net, nbus::Int)

Make sure `net.qmin_pu`/`net.qmax_pu` exist with length `nbus`.
No field rebinds; only mutate/resize/fill.
"""
function _ensure_limit_storage!(net::Net, nbus::Int)
    if length(net.qmin_pu) != nbus
        empty!(net.qmin_pu); empty!(net.qmax_pu)
        append!(net.qmin_pu, fill(-Inf, nbus))
        append!(net.qmax_pu, fill( Inf, nbus))
    else
        fill!(net.qmin_pu, -Inf)
        fill!(net.qmax_pu,  Inf)
    end
    return nothing
end

"""
    _ensure_log_storage!(net::Net)

Ensure logging containers exist (mutating only).
"""
function _ensure_log_storage!(net::Net)
    # Dict exists by type; just ensure it's empty or usable
    if net.qLimitEvents === nothing
        # very defensive; with typed field this should not happen
        empty!(net.qLimitEvents)
    end
    # qLimitLog is a Vector{NamedTuple{(:iter,:bus,:side),...}}; it exists by construction
    return nothing
end

# ---------- Public API ----------

"""
    buildQLimits!(net::Net)

Aggregate generator Q-limits per bus and store them in `net.qmin_pu`/`net.qmax_pu`
(in p.u., base on `net.baseMVA`). Only generators contribute.

- If multiple generators are connected to a bus, take min of all `minQ` and max of all `maxQ`.
- Missing limits remain at ±Inf.
"""
function buildQLimits!(net::Net)
    nbus = length(net.nodeVec)
    _ensure_limit_storage!(net, nbus)

    for ps in net.prosumpsVec
        isGenerator(ps) || continue
        bus = _prosumer_bus_index(ps)

        qmin_pu = isnothing(ps.minQ) ? -Inf : ps.minQ / net.baseMVA
        qmax_pu = isnothing(ps.maxQ) ?  Inf : ps.maxQ / net.baseMVA

        # aggregate
        net.qmin_pu[bus] = min(net.qmin_pu[bus], qmin_pu)
        net.qmax_pu[bus] = max(net.qmax_pu[bus], qmax_pu)
    end

    _ensure_log_storage!(net)
    return nothing
end

"""
    resetQLimitLog!(net::Net)

Clears PV→PQ switch events and the structured log.
"""
function resetQLimitLog!(net::Net)
    empty!(net.qLimitEvents)
    empty!(net.qLimitLog)
    return nothing
end

"""
    logQLimitHit!(net::Net, iter::Int, bus::Int, side::Symbol)

Append a structured Q-limit event: `side` is `:min` or `:max`.
Also mirrors the last state per bus into `net.qLimitEvents[bus]`.
"""
function logQLimitHit!(net::Net, iter::Int, bus::Int, side::Symbol)
    push!(net.qLimitLog, (iter=iter, bus=bus, side=side))
    net.qLimitEvents[bus] = side
    return nothing
end

"""
    get_Q_limits_pu(net::Net) -> (qmin_pu, qmax_pu)

Return per-bus Q limits in p.u. (build once if empty).
"""
function get_Q_limits_pu(net::Net)
    if isempty(net.qmin_pu) || isempty(net.qmax_pu)
        buildQLimits!(net)
    end
    return net.qmin_pu, net.qmax_pu
end

# Legacy compatibility: some callers pass a second parameter (e.g., bus count)
function get_Q_limits_pu(net::Net, ::Int)
    @warn "get_Q_limits_pu: additional argument ignored (compatibility fallback)."
    return get_Q_limits_pu(net)
end

"""
    lastQLimitIter(net::Net, bus::Int) -> Union{Nothing,Int}

Return the iteration index of the most recent Q-limit event for the given bus, if any.
"""
function lastQLimitIter(net::Net, bus::Int)::Union{Nothing,Int}
    for i = length(net.qLimitLog):-1:1
        ev = net.qLimitLog[i]
        if ev.bus == bus
            return ev.iter
        end
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
