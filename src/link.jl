# Author: Udo Schmitz (https://github.com/Welthulk)
# include-file link.jl

"""
    BusLink

Topological, impedance-less connection between two buses. Bus links are not part
of the electrical branch model (YBUS). They are intended for post-power-flow KCL
allocation, e.g. busbar couplers / sectionalizers.
"""
mutable struct BusLink
  linkIdx::Int
  fromBus::Int
  toBus::Int
  status::Int
  pFlow_MW::Union{Nothing,Float64}
  qFlow_MVar::Union{Nothing,Float64}
  iFrom_kA::Union{Nothing,Float64}
  iTo_kA::Union{Nothing,Float64}

  function BusLink(; linkIdx::Int, fromBus::Int, toBus::Int, status::Int = 1)
    @assert fromBus != toBus "fromBus and toBus must differ"
    @assert status in (0, 1) "status must be 0 or 1"
    new(linkIdx, fromBus, toBus, status, nothing, nothing, nothing, nothing)
  end

  function Base.show(io::IO, l::BusLink)
    print(io, "BusLink(")
    print(io, "idx=", l.linkIdx, ", ")
    print(io, "from=", l.fromBus, ", ")
    print(io, "to=", l.toBus, ", ")
    print(io, "status=", l.status)
    if !isnothing(l.pFlow_MW)
      print(io, ", P=", l.pFlow_MW, " MW")
    end
    if !isnothing(l.qFlow_MVar)
      print(io, ", Q=", l.qFlow_MVar, " MVar")
    end
    if !isnothing(l.iFrom_kA)
      print(io, ", Ifrom=", l.iFrom_kA, " kA")
    end
    if !isnothing(l.iTo_kA)
      print(io, ", Ito=", l.iTo_kA, " kA")
    end
    print(io, ")")
  end
end

function setLinkStatus!(link::BusLink, status::Int)
  @assert status in (0, 1) "status must be 0 or 1"
  link.status = status
end

function setLinkFlow!(link::BusLink, pFlow_MW::Float64, qFlow_MVar::Float64)
  link.pFlow_MW = pFlow_MW
  link.qFlow_MVar = qFlow_MVar
end

function setLinkCurrent!(link::BusLink, iFrom_kA::Float64, iTo_kA::Float64)
  link.iFrom_kA = iFrom_kA
  link.iTo_kA = iTo_kA
end
