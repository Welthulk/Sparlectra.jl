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

mutable struct Branch
  comp::AbstractComponent
  fromBus::Integer        
  toBus::Integer          
  r_pu::Float64                          # resistance
  x_pu::Float64                          # reactance
  b_pu::Float64                          # total line charging susceptance
  g_pu::Float64                          # total line charging conductance
  ratio::Float64                         # transformer off nominal turns ratio
  angle::Float64                         # transformer off nominal phase shift angle
  status::Integer                        # 1 = in service, 0 = out of service
  sn_MVA::Union{Nothing,Float64}         # nominal power of the branch = rateA
  fBranchFlow::Union{Nothing,BranchFlow} # flow from fromNodeID to toNodeID
  tBranchFlow::Union{Nothing,BranchFlow} # flow from toNodeID to fromNodeID

  function Branch(; from::Int, to::Int, baseMVA::Float64, branch::AbstractBranch, id::Int, status::Integer=1,  ratio::Union{Nothing,Float64}=nothing,  side::Union{Nothing,Int}=nothing, vn_kV::Union{Nothing,Float64}=nothing )    
    if isa(branch, ACLineSegment)
      @assert !isnothing(vn_kV) "vn_kV must be set for an ACLineSegment"
      c = getBranchComp(vn_kV, from, to, id, "ACLine")
      r, x, b, g = getRXBG(branch)
      baseZ = (vn_kV)^2 / baseMVA
      r_pu = r / baseZ
      x_pu = x / baseZ
      b_pu = b * baseZ
      g_pu = g * baseZ
      if isnothing(ratio)
        ratio = 0.0
      end
      new(c, from, to, r_pu, x_pu, b_pu, g_pu, ratio, 0.0, status, nothing, nothing, nothing)
    elseif isa(branch, PowerTransformer)      
      if (isnothing(side) && branch.isBiWinder)        
        side = getSideNumber2WT(branch)
      else
        error("side must be set for a PowerTransformer")      
      end
      
      w = (side in [1, 2, 3]) ? (side == 1 ? branch.side1 : (side == 2 ? branch.side2 : branch.side3)) : error("wrong value for 'side'")
      if isnothing(vn_kV)
        vn_kV = w.Vn 
      end   
      c = getBranchComp(vn_kV, from, to, id, "Transformer")
      sn_MVA = w.ratedS
      r, x, b, g = getRXBG(w)
      baseZ = (vn_kV)^2 / baseMVA
      r_pu = r / baseZ
      x_pu = x / baseZ
      b_pu = b * baseZ
      g_pu = g * baseZ
      if isnothing(ratio)
        ratio = 0.0
      end

      new(c, from, to, r_pu, x_pu, b_pu, g_pu, ratio, w.shift_degree, status, sn_MVA, nothing, nothing)
    end
  end

  function Base.show(io::IO, b::Branch)
    print(io, "Branch( ")
    print(io, b.comp, ", ")
    print(io, "fromBus: ", b.fromBus, " ($(b._from))", ", ")
    print(io, "toBus: ", b.toBus, " ($(b._to))", ", ")

    print(io, "r_pu: ", b.r_pu, ", ")
    print(io, "x_pu: ", b.x_pu, ", ")
    print(io, "b_pu: ", b.b_pu, ", ")
    print(io, "g_pu: ", b.g_pu, ", ")
    print(io, "ratio: ", b.ratio, ", ")
    print(io, "angle: ", b.angle, ", ")
    print(io, "status: ", b.status, ", ")
    print(io, "parallel: ", b.isParallel, ", ")

    if !isnothing(b.sn_MVA)
      print(io, "sn_MVA: ", b.sn_MVA, ", ")
    end

    if (!isnothing(b.fBranchFlow))
      print(io, "BranchFlow (from): ", b.fBranchFlow, ", ")
    end

    if (!isnothing(b.tBranchFlow))
      print(io, "BranchFlow (to): ", b.tBranchFlow, ", ")
    end

    println(io, ")")
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

function getBranchComp(Vn_kV::Float64, from::Int, to::Int, idx::Int, kind::String)
  cTyp = toComponentTyp("Branch")
  name = "B_$(kind)_$(string(convert(Int,trunc(Vn_kV))))_$(Int(from))_$(Int(to))"
  cID = "#" * name * "#" * string(idx)
  return ImpPGMComp(cID, name, cTyp, Vn_kV, from, to)
end
