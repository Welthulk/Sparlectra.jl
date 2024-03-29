# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file branch.jl

# helper
struct BranchFlow
  vm_pu::Union{Nothing,Float64} # voltage magnitude
  va_deg::Union{Nothing,Float64} # voltage angle
  pFlow::Union{Nothing,Float64} # active power flow
  qFlow::Union{Nothing,Float64} # reactive power flow
  function BranchFlow(vm_pu::Union{Nothing,Float64} = nothing, va_deg::Union{Nothing,Float64} = nothing, pFlow::Union{Nothing,Float64} = nothing, qFlow::Union{Nothing,Float64} = nothing)
    new(vm_pu, va_deg, pFlow, qFlow)
  end
  function Base.show(io::IO, b::BranchFlow)
    print(io, "BranchFlow( ")
    print(io, "vm: ", b.vm_pu, ", ")
    print(io, "va: ", b.va_deg, ", ")
    print(io, "pFlow: ", b.pFlow, ", ")
    print(io, "qFlow: ", b.qFlow, ", ")
    println(io, ")")
  end
end
"""
Purpose: Branch to connect two nodes and save pq-flow-data
"""
mutable struct Branch
  comp::AbstractComponent
  fromBus::Integer        # Bus number of the from bus
  toBus::Integer          # Bus number of the to bus    
  _from::Integer          # original bus number of the from bus
  _to::Integer            # original bus number of the to bus
  fromNodeID::String      # Node ID starting node (Seite 1)
  toNodeID::String        # Node ID ending node (Seite 2)
  fBranchFlow::Union{Nothing,BranchFlow} # flow from fromNodeID to toNodeID
  tBranchFlow::Union{Nothing,BranchFlow} # flow from toNodeID to fromNodeID
  r_pu::Float64           # resistance
  x_pu::Float64           # reactance
  b_pu::Float64           # total line charging susceptance, kapazitiver Anteil
  g_pu::Float64           # total line charging conductance, ohmscher Anteil
  ratio::Float64          # transformer off nominal turns ratio
  angle::Float64          # transformer off nominal phase shift angle
  status::Integer         # 1 = in service, 0 = out of service
  isParallel::Bool        # is a parallel branch? (true/false)  

  function Branch(
    branchC::Component,
    fromBus::Integer,
    toBus::Integer,
    fromNodeID::String,
    toNodeID::String,
    r_pu::Float64,
    x_pu::Float64,
    b_pu::Float64,
    g_pu::Float64,
    ratio::Float64,
    angle::Float64,
    status::Integer,
    fBracnhFlow::Union{Nothing,BranchFlow} = nothing,
    tBracnhFlow::Union{Nothing,BranchFlow} = nothing,
    isParallel::Bool = false,
  )
    new(branchC, fromBus, toBus, fromBus, toBus, fromNodeID, toNodeID, fBracnhFlow, tBracnhFlow, r_pu, x_pu, b_pu, g_pu, ratio, angle, status, isParallel)
  end

  function Branch(
    branchC::Component,
    fromBus::Integer,
    toBus::Integer,
    fromOrigBus::Integer,
    toOrigBus::Integer,
    fromNodeID::String,
    toNodeID::String,
    r_pu::Float64,
    x_pu::Float64,
    b_pu::Float64,
    g_pu::Float64,
    ratio::Float64,
    angle::Float64,
    status::Integer,
    fBracnhFlow::Union{Nothing,BranchFlow} = nothing,
    tBracnhFlow::Union{Nothing,BranchFlow} = nothing,
    isParallel::Bool = false,
  )
    new(branchC, fromBus, toBus, fromOrigBus, toOrigBus, fromNodeID, toNodeID, fBracnhFlow, tBracnhFlow, r_pu, x_pu, b_pu, g_pu, ratio, angle, status, isParallel)
  end

  function Base.show(io::IO, b::Branch)
    print(io, "Branch( ")
    print(io, b.comp, ", ")
    print(io, "FromBus: ", b.fromBus, " ($(b._from))", ", ")
    print(io, "ToBus: ", b.toBus, " ($(b._to))", ", ")
    print(io, "FromNodeID: ", b.fromNodeID, ", ")
    print(io, "ToNodeID: ", b.toNodeID, ", ")

    if (!isnothing(b.fBranchFlow))
      print(io, "BranchFlow (from): ", b.fBranchFlow, ", ")
    end

    if (!isnothing(b.tBranchFlow))
      print(io, "BranchFlow (to): ", b.tBranchFlow, ", ")
    end

    print(io, "r_pu: ", b.r_pu, ", ")
    print(io, "x_pu: ", b.x_pu, ", ")
    print(io, "b_pu: ", b.b_pu, ", ")
    print(io, "g_pu: ", b.g_pu, ", ")
    print(io, "ratio: ", b.ratio, ", ")
    print(io, "angle: ", b.angle, ", ")
    print(io, "status: ", b.status, ", ")
    print(io, "parallel: ", b.isParallel, ")")
  end
end

# helper
function setBranchFlow!(tfBranchFlow::BranchFlow, fBranchFlow::BranchFlow, branch::Branch)
  branch.tBranchFlow = tfBranchFlow
  branch.fBranchFlow = fBranchFlow
end

# helper
function setBranchStatus!(service::Bool, branch::Branch)
  if (service)
    branch.status = 1
  else
    branch.status = 0
  end
end
