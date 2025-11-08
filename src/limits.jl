# src/limits.jl

"""
    _prosumer_bus_index(ps::ProSumer) -> Union{Int,Nothing}

Internal helper: derive the bus index of a `ProSumer` from its component.
Prefers ImpPGMComp's `cFrom_bus`. Tries a few other common names as fallback.
Returns `nothing` if none found (callers decide how to handle it).
Must NOT throw.
"""
function _prosumer_bus_index(ps::ProSumer)::Union{Int,Nothing}
    c = ps.comp
    # Prefer ImpPGMComp field names
    for sym in (:cFrom_bus, :cTo_bus, :fromBus, :from, :busIdx, :bus, :node, :idx)
        if hasproperty(c, sym)
            val = getfield(c, sym)
            return isnothing(val) ? nothing : Int(val)
        end
    end
    return nothing
end

"""
    buildQLimits!(net::Net) -> Nothing

Aggregates reactive power (Q) limits from generator prosumers to bus level and
stores them into `net.qmin_pu` / `net.qmax_pu` (per unit, referred to `net.baseMVA`).

Rules:
- Only generators (`isGenerator(ps)`) contribute.
- Multiple generators at the same bus: take the min of all minQ (p.u.) and
  the max of all maxQ (p.u.).
- Missing limits remain at ±Inf.

`Net` is immutable → mutate arrays/dicts in-place (resize!/fill!/empty!).
"""
function buildQLimits!(net::Net)
    nbus = length(net.nodeVec)

    # Prepare arrays in-place
    resize!(net.qmin_pu, nbus); fill!(net.qmin_pu, -Inf)
    resize!(net.qmax_pu, nbus); fill!(net.qmax_pu,  Inf)

    for ps in net.prosumpsVec
        isGenerator(ps) || continue

        bus = _prosumer_bus_index(ps)
        if bus === nothing
            cname = hasproperty(ps.comp, :cName) ? String(getfield(ps.comp, :cName)) :
                    (hasproperty(ps.comp, :name)  ? String(getfield(ps.comp, :name))  :
                     string(typeof(ps.comp)))
            @warn "buildQLimits!: cannot determine bus index for prosumer component $(cname); skipping its Q-limits"
            continue
        end
        @assert 1 <= bus <= nbus "Prosumer bus index $bus out of range 1:$nbus"

        # Convert MVar -> p.u.
        qmin_pu = isnothing(ps.minQ) ? -Inf : ps.minQ / net.baseMVA
        qmax_pu = isnothing(ps.maxQ) ?  Inf : ps.maxQ / net.baseMVA

        # Aggregate per bus
        net.qmin_pu[bus] = min(net.qmin_pu[bus], qmin_pu)
        net.qmax_pu[bus] = max(net.qmax_pu[bus], qmax_pu)
    end

    # Reset PV→PQ event log (mutate dict, don't replace)
    empty!(net.qLimitEvents)
    return nothing
end

"""
    resetQLimitLog!(net::Net) -> Nothing

Clears the event log for Q-limit driven PV→PQ switches.
"""
function resetQLimitLog!(net::Net)
    empty!(net.qLimitEvents)
    return nothing
end

"""
    get_Q_limits_pu(net::Net) -> (qmin_pu::Vector{Float64}, qmax_pu::Vector{Float64})

Returns the per-bus Q-limits in p.u. Builds them lazily on first use.
"""
function get_Q_limits_pu(net::Net)
    if isempty(net.qmin_pu) || isempty(net.qmax_pu)
        buildQLimits!(net)
    end
    return net.qmin_pu, net.qmax_pu
end

"""
    get_Q_limits_pu(net::Net, _) -> (qmin_pu, qmax_pu)

Compatibility shim: older callers may pass an additional argument
(e.g., a bus count). The extra argument is ignored.
"""
function get_Q_limits_pu(net::Net, ::Any)
    @warn "get_Q_limits_pu: additional argument ignored (compatibility fallback)."
    return get_Q_limits_pu(net)
end
