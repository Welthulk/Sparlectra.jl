# Q-limit utilities (wird im Modul Sparlectra inkludiert)

"""
    buildQLimits!(net::Net)

Aggregiert je Bus die Q-Grenzen (p.u.) aus Prosumer-/Generator-Objekten und
schreibt sie in `net.qmin_pu` und `net.qmax_pu`.

Konvention:
- Fehlende Grenzen → ±Inf
- Mehrere Prosumer am Bus → engste Grenzen (max der Minima, min der Maxima)
"""
function buildQLimits!(net::Net)
    nb = length(net.nodeVec)
    qmin = fill(-Inf, nb)
    qmax = fill(+Inf, nb)
    Sbase = net.baseMVA
    @inbounds for ps in net.prosumpsVec
        b = ps.busIdx
        if !isnothing(ps.minQ); qmin[b] = max(qmin[b], ps.minQ / Sbase) end
        if !isnothing(ps.maxQ); qmax[b] = min(qmax[b], ps.maxQ / Sbase) end
    end
    net.qmin_pu = qmin
    net.qmax_pu = qmax
    return net
end

"""
    resetQLimitLog!(net::Net)

Löscht das Q-Limit-Ereignis-Logging im Netz.
"""
function resetQLimitLog!(net::Net)
    empty!(net.qLimitEvents)
    return net
end
