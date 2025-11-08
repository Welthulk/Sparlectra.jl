"""
    get_Q_limits_pu(net::Net, nb::Int) -> (qmin_pu::Vector{Float64}, qmax_pu::Vector{Float64})

Aggregiert je Bus die Q-Grenzen (p.u.) aus Prosumer/Generator-Objekten.
Mehrere Prosumer am Bus werden durch „engste“ Grenzen kombiniert.
Fehlende Grenzen -> ±Inf.
"""
function get_Q_limits_pu(net::Net, nb::Int)
    qmin = fill(-Inf, nb)
    qmax = fill(+Inf, nb)
    Sbase = net.baseMVA
    @inbounds for ps in net.prosumpsVec
        b = ps.busIdx
        if !isnothing(ps.minQ)
            qmin[b] = max(qmin[b], ps.minQ / Sbase)
        end
        if !isnothing(ps.maxQ)
            qmax[b] = min(qmax[b], ps.maxQ / Sbase)
        end
    end
    return qmin, qmax
end

