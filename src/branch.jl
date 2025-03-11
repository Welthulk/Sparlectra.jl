# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file branch.jl

# helper
"""
    BranchFlow

A structure representing the flow in a branch of a power system.

# Fields
- `vm_pu::Union{Nothing,Float64}`: The voltage magnitude in per unit.
- `va_deg::Union{Nothing,Float64}`: The voltage angle in degrees.
- `pFlow::Union{Nothing,Float64}`: The active power flow.
- `qFlow::Union{Nothing,Float64}`: The reactive power flow.

# Constructors
- `BranchFlow(vm_pu::Union{Nothing,Float64} = nothing, va_deg::Union{Nothing,Float64} = nothing, pFlow::Union{Nothing,Float64} = nothing, qFlow::Union{Nothing,Float64} = nothing)`: Creates a new `BranchFlow` instance.

# Methods
- `Base.show(io::IO, b::BranchFlow)`: Prints the `BranchFlow` instance.

# Example
```julia
BranchFlow(vm_pu = 1.0, va_deg = 0.0, pFlow = 100.0, qFlow = 50.0)
```
"""
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
    BranchModel

A structure representing a branch model in a power system.

# Fields
- `r_pu::Float64`: The per unit resistance of the branch.
- `x_pu::Float64`: The per unit reactance of the branch.
- `b_pu::Float64`: The per unit total line charging susceptance of the branch.
- `g_pu::Float64`: The per unit total line charging conductance of the branch.
- `ratio::Float64`: The transformer off nominal turns ratio.
- `angle::Float64`: The transformer off nominal phase shift angle.
- `sn_MVA::Union{Nothing,Float64}`: The nominal power of the branch = rateA.

# Constructors
- `BranchModel(; r_pu::Float64, x_pu::Float64, b_pu::Float64, g_pu::Float64, ratio::Float64, angle::Float64, sn_MVA::Union{Nothing,Float64} = nothing)`: Creates a new `BranchModel` instance.

# Example
```julia
BranchModel(r_pu = 0.01, x_pu = 0.1, b_pu = 0.02, g_pu = 0.02, ratio = 1.0, angle = 0.0, sn_MVA = 100.0)
```
"""
struct BranchModel
  r_pu::Float64
  x_pu::Float64
  b_pu::Float64
  g_pu::Float64
  ratio::Float64
  angle::Float64
  sn_MVA::Union{Nothing,Float64}
end

"""
    Branch

A mutable structure representing a branch in a power system.

# Fields
- `comp::AbstractComponent`: The component of the branch.
- `fromBus::Integer`: The index of the bus where the branch starts.
- `toBus::Integer`: The index of the bus where the branch ends.
- `r_pu::Float64`: The per unit resistance of the branch.
- `x_pu::Float64`: The per unit reactance of the branch.
- `b_pu::Float64`: The per unit total line charging susceptance of the branch.
- `g_pu::Float64`: The per unit total line charging conductance of the branch.
- `ratio::Float64`: The transformer off nominal turns ratio.
- `angle::Float64`: The transformer off nominal phase shift angle.
- `status::Integer`: The status of the branch. 1 = in service, 0 = out of service.
- `sn_MVA::Union{Nothing,Float64}`: The nominal power of the branch = rateA.
- `fBranchFlow::Union{Nothing,BranchFlow}`: The flow from fromNodeID to toNodeID.
- `tBranchFlow::Union{Nothing,BranchFlow}`: The flow from toNodeID to fromNodeID.
- `pLosses::Union{Nothing,Float64}`: The active power losses.
- `qLosses::Union{Nothing,Float64}`: The reactive power losses.

# Constructors
- `Branch(; from::Int, to::Int, baseMVA::Float64, branch::AbstractBranch, id::Int, status::Integer = 1, ratio::Union{Nothing,Float64} = nothing, side::Union{Nothing,Int} = nothing, vn_kV::Union{Nothing,Float64} = nothing,
                    fromOid::Union{Nothing,Int} = nothing, toOid::Union{Nothing,Int} = nothing)`: Creates a new `Branch` instance.

# Methods
- `Base.show(io::IO, b::Branch)`: Prints the `Branch` instance.
"""
mutable struct Branch <: AbstractBranch
  comp::AbstractComponent
  branchIdx::Int
  fromBus::Integer
  toBus::Integer
  r_pu::Float64                          # resistance
  x_pu::Float64                          # reactance
  b_pu::Float64                          # total line charging susceptance
  g_pu::Float64                          # total line charging conductance
  ratio::Float64                         # nominal turns ratio
  angle::Float64                         # nominal phase shift angle in degrees
  status::Integer                        # 1 = in service, 0 = out of service
  sn_MVA::Union{Nothing,Float64}         # nominal power of the branch = rateA
  fBranchFlow::Union{Nothing,BranchFlow} # flow from fromNodeID to toNodeID
  tBranchFlow::Union{Nothing,BranchFlow} # flow from toNodeID to fromNodeID
  pLosses::Union{Nothing,Float64}        # active power losses
  qLosses::Union{Nothing,Float64}        # reactive power losses

  function Branch(;
    branchIdx::Int,
    from::Int,
    to::Int,
    baseMVA::Float64,
    branch::AbstractBranch,
    id::Int,
    status::Integer = 1,
    ratio::Union{Nothing,Float64} = nothing,
    side::Union{Nothing,Int} = nothing,
    vn_kV::Union{Nothing,Float64} = nothing,
    fromOid::Union{Nothing,Int} = nothing,
    toOid::Union{Nothing,Int} = nothing,
  )
    if isa(branch, ACLineSegment) # Line
      @assert !isnothing(vn_kV) "vn_kV must be set for an ACLineSegment"
      if isnothing(ratio)
        ratio = 0.0
      end
      if !isnothing(fromOid) && !isnothing(toOid)
        c = getBranchComp(vn_kV, fromOid, toOid, id, "ACL")
      else
        c = getBranchComp(vn_kV, from, to, id, "ACL")
      end
      c = getBranchComp(vn_kV, from, to, id, "ACL")

      r_pu, x_pu, b_pu, g_pu = getLineRXBG_pu(branch, vn_kV, baseMVA)
      new(c, branchIdx, from, to, r_pu, x_pu, b_pu, g_pu, 0.0, 0.0, status, branch.ratedS, nothing, nothing, nothing, nothing)
    elseif isa(branch, PowerTransformer) # Transformer     
      if (isnothing(side) && branch.isBiWinder)
        side = getSideNumber2WT(branch)
      elseif (isnothing(side) && !branch.isBiWinder)
        error("side must be set for a PowerTransformer")
      end

      c = if !isnothing(fromOid) && !isnothing(toOid)
        getBranchComp(vn_kV, fromOid, toOid, id, "2WT")
      else
        getBranchComp(vn_kV, from, to, id, "2WT")
      end

      w = (side in [1, 2, 3]) ? (side == 1 ? branch.side1 : (side == 2 ? branch.side2 : branch.side3)) : error("wrong value for 'side'")
      vn_kV = isnothing(vn_kV) ? w.Vn : vn_kV
      sn_MVA = getWindingRatedS(w)
      r_pu, x_pu, b_pu, g_pu = getTrafoRXBG_pu(w, vn_kV, baseMVA)

      ratio = isnothing(ratio) ? 1.0 : ratio
      @assert ratio != 0.0 "ratio must not be 0.0 for transformers"
      angle = isnothing(w.shift_degree) ? 0.0 : w.shift_degree

      new(c, branchIdx, from, to, r_pu, x_pu, b_pu, g_pu, ratio, angle, status, sn_MVA, nothing, nothing, nothing, nothing)
    elseif isa(branch, BranchModel) # PI-Model
      @assert !isnothing(vn_kV) "vn_kV must be set for PI-Model"

      if !isnothing(fromOid) && !isnothing(toOid)
        c = getBranchComp(vn_kV, fromOid, toOid, id, "PI")
      else
        c = getBranchComp(vn_kV, from, to, id, "PI")
      end

      new(c, branchIdx, from, to, branch.r_pu, branch.x_pu, branch.b_pu, branch.g_pu, branch.ratio, branch.angle, status, branch.sn_MVA, nothing, nothing, nothing, nothing)
    else
      error("Branch type not supported")
    end
  end

  function Base.show(io::IO, b::Branch)
    print(io, "Branch( ")
    print(io, b.comp, ", ")
    print(io, "branchIdx: ", b.branchIdx, ", ")
    print(io, "fromBus: ", b.fromBus, ", ")
    print(io, "toBus: ", b.toBus, ", ")

    print(io, "r_pu: ", b.r_pu, ", ")
    print(io, "x_pu: ", b.x_pu, ", ")
    print(io, "b_pu: ", b.b_pu, ", ")
    print(io, "g_pu: ", b.g_pu, ", ")
    print(io, "ratio: ", b.ratio, ", ")
    print(io, "angle: ", b.angle, ", ")
    print(io, "status: ", b.status, ", ")
    if !isnothing(b.sn_MVA)
      print(io, "sn_MVA: ", b.sn_MVA, ", ")
    end
    if (!isnothing(b.fBranchFlow))
      print(io, "BranchFlow (from): ", b.fBranchFlow, ", ")
    end
    if (!isnothing(b.tBranchFlow))
      print(io, "BranchFlow (to): ", b.tBranchFlow, ", ")
    end
    if (!isnothing(b.pLosses))
      print(io, "pLosses: ", b.pLosses, ", ")
    end
    if (!isnothing(b.qLosses))
      print(io, "qLosses: ", b.qLosses, ", ")
    end

    println(io, ")")
  end
end

function calcAdmittance(branch::Branch, u_rated::Float64, s_rated::Float64)::Tuple{ComplexF64,ComplexF64,ComplexF64,ComplexF64}
  # Series Admittance ys
  ys = calcBranchYser(branch)
  # Shunt Admittance ysh
  ysh = calcBranchYshunt(branch)
  # calc complex ratio
  t = calcBranchRatio(branch)
  # Calculate Y_from_from, Y_from_to, Y_to_from, Y_to_to
  Y_11 = (ys + 0.5 * ysh) / abs2(t)
  Y_12 = -1.0 * ys / conj(t)
  Y_21 = -1.0 * ys / t
  Y_22 = ys + 0.5 * ysh
  return (Y_11, Y_12, Y_21, Y_22)
end

# helper
function setBranchFlow!(branch::Branch, tfBranchFlow::BranchFlow, fBranchFlow::BranchFlow)
  branch.tBranchFlow = tfBranchFlow
  branch.fBranchFlow = fBranchFlow
end

# helper
function setBranchStatus!(branch::Branch, service::Bool)
  if (service)
    branch.status = 1
  else
    branch.status = 0
  end
end

function getBranchFlow(branch::Branch, from::Node, to::Node)
  if (branch.fromBus == from.busIdx && branch.toBus == to.busIdx)
    return branch.fBranchFlow
  elseif (branch.fromBus == to.busIdx && branch.toBus == from.busIdx)
    return branch.tBranchFlow
  else
    error("Nodes do not match the branch")
  end
end

function getBranchIdx(branch::Branch)
  return branch.branchIdx
end

function getBranchLosses(branch::Branch)
  return branch.pLosses, branch.qLosses
end

function setBranchLosses!(branch::Branch, pLosses::Float64, qLosses::Float64)
  branch.pLosses = pLosses
  branch.qLosses = qLosses
end

function calcBranchYser(branch::Branch)::ComplexF64
  return inv((branch.r_pu + branch.x_pu * im))
end

function calcBranchYshunt(branch::Branch)::ComplexF64
  return (branch.g_pu + branch.b_pu * im)
end

function calcBranchRatio(branch::Branch)::ComplexF64
  t = 1.0 + 0.0 * im
  ratio = (isnothing(branch.ratio) || branch.ratio == 0.0) ? 1.0 : branch.ratio
  shift = isnothing(branch.angle) ? 0.0 : branch.angle
  t = (shift != 0.0) ? calcComplexRatio(tapRatio = ratio, angleInDegrees = shift) : 1.0 + 0.0im
  return t
end

function getBranchComp(Vn_kV::Float64, from::Int, to::Int, idx::Int, kind::String)
  cTyp = toComponentTyp("Branch")
  name = "B_$(kind)_$(string(convert(Int,trunc(Vn_kV))))_$(Int(from))_$(Int(to))"
  cID = "#" * name * "#" * string(idx)
  return ImpPGMComp(cID, name, cTyp, Vn_kV, from, to)
end
