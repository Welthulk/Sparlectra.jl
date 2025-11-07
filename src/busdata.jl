# busdata.jl — Datentyp für NR / Lastfluss
# Dieser Typ war zuvor lokal in jacobian.jl; ausgelagert zur gemeinsamen Nutzung.

mutable struct BusData
    idx::Int         # Index des Busses (nach Sortierung)
    vm_pu::Float64   # Spannung in p.u.
    va_rad::Float64  # Winkel in rad
    pƩ::Float64      # Summe Wirkleistung (p.u.)
    qƩ::Float64      # Summe Blindleistung (p.u.)
    _pRes::Float64   # berechnete Wirkleistung (p.u.)
    _qRes::Float64   # berechnete Blindleistung (p.u.)
    type::NodeType

    function BusData(idx::Int, vm_pu::Float64, va_rad::Float64,
                     sumP::Float64, sumQ::Float64, type::Sparlectra.NodeType)
        new(idx, vm_pu, va_rad, sumP, sumQ, 0.0, 0.0, type)
    end
end

function Base.show(io::IO, bus::BusData)
    va_deg = round(rad2deg(bus.va_rad), digits = 3)
    print(io, "BusData($(bus.idx), $(bus.vm_pu), $(va_deg)°), $(bus.pƩ), $(bus.qƩ), $(bus._pRes), $(bus._qRes), $(bus.type))")
end
