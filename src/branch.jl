# Copyright 2023–2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# file: src/branch.jl
# purpose: branch model types (Branch, BranchModel, BranchFlow) and
#          admittance, ratio, flow, and loss helpers for network branches

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
  tap_ratio::Float64
  phase_shift_deg::Float64
  has_ratio_tap::Bool
  has_phase_tap::Bool
  tap_min::Float64
  tap_max::Float64
  tap_step::Float64
  phase_min_deg::Float64
  phase_max_deg::Float64
  phase_step_deg::Float64
  # per-terminal service state (r0.10.0): 1 = closed, 0 = open. The
  # aggregate `status` stays the user-facing switch and is kept consistent
  # by the setters (status = 1 iff both terminals closed); consumers must
  # read the state through _branch_terminal_state, never the raw fields.
  from_status::Int
  to_status::Int
  # open-end voltage of a one-sided open branch (result, not a bus
  # voltage): filled by calcNetLosses! from the pi-model voltage divider,
  # nothing while closed or fully open
  open_end_vm_pu::Union{Nothing,Float64}
  open_end_va_deg::Union{Nothing,Float64}

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
    angle::Union{Nothing,Float64} = nothing,
    values_are_pu::Bool=false,
    from_status::Union{Nothing,Integer} = nothing,
    to_status::Union{Nothing,Integer} = nothing,
  )
    # terminal flags default to the aggregate; the stored aggregate is then
    # recomputed so that status = 1 iff both terminals are closed
    fs = Int(something(from_status, status))
    ts = Int(something(to_status, status))
    fs in (0, 1) && ts in (0, 1) || error("Branch: from_status/to_status must be 0 or 1")
    status = (fs == 1 && ts == 1) ? 1 : 0
    if isa(branch, ACLineSegment) # Line
      @assert !isnothing(vn_kV) "vn_kV must be set for an ACLineSegment"
      if isnothing(ratio)
        # to distinguish line from transformer
        ratio  = 0.0
      end
      if !isnothing(fromOid) && !isnothing(toOid)
        c = getBranchComp(vn_kV, fromOid, toOid, id, "ACL")
      else
        c = getBranchComp(vn_kV, from, to, id, "ACL")
      end
      if values_are_pu
        r_pu, x_pu, b_pu, g_pu = getLineRXBG(branch)
      else 
        r_pu, x_pu, b_pu, g_pu = getLineRXBG_pu(branch, vn_kV, baseMVA)        
      end    
      new(c, branchIdx, from, to, r_pu, x_pu, b_pu, g_pu, 0.0, 0.0, status, branch.ratedS, nothing, nothing, nothing, nothing, 1.0, 0.0, false, false, 0.9, 1.1, 0.00625, -30.0, 30.0, 1.25, fs, ts, nothing, nothing)
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
      
      if values_are_pu
        r_pu, x_pu, b_pu, g_pu = getTrafoRXBG(w)
      else 
        r_pu, x_pu, b_pu, g_pu = getTrafoRXBG_pu(w, vn_kV, baseMVA)        
      end
      

      ratio = isnothing(ratio) ? 1.0 : ratio
      @assert ratio != 0.0 "ratio must not be 0.0 for transformers"
      #angle = isnothing(w.shift_degree) ? 0.0 : w.shift_degree
      angle = isnothing(angle) ? (isnothing(w.shift_degree) ? 0.0 : w.shift_degree) : angle
      tap_min = 0.9
      tap_max = 1.1
      tap_step = 0.00625
      if !isnothing(w.taps)
        tap_min, tap_max, tap_step = calcRatioTapRange(w.taps)
      end

      new(c, branchIdx, from, to, r_pu, x_pu, b_pu, g_pu, ratio, angle, status, sn_MVA, nothing, nothing, nothing, nothing, ratio, angle, true, true, tap_min, tap_max, tap_step, -30.0, 30.0, 1.25, fs, ts, nothing, nothing)
    elseif isa(branch, BranchModel) # PI-Model
      @assert !isnothing(vn_kV) "vn_kV must be set for PI-Model"

      if !isnothing(fromOid) && !isnothing(toOid)
        c = getBranchComp(vn_kV, fromOid, toOid, id, "PI")
      else
        c = getBranchComp(vn_kV, from, to, id, "PI")
      end

      is_tap = branch.ratio != 0.0
      initial_ratio = is_tap ? branch.ratio : 1.0
      initial_angle = is_tap ? branch.angle : 0.0
      new(c, branchIdx, from, to, branch.r_pu, branch.x_pu, branch.b_pu, branch.g_pu, branch.ratio, branch.angle, status, branch.sn_MVA, nothing, nothing, nothing, nothing, initial_ratio, initial_angle, is_tap, is_tap, 0.9, 1.1, 0.00625, -30.0, 30.0, 1.25, fs, ts, nothing, nothing)
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
    if _branch_terminal_state(b) != :closed
      print(io, "terminal_state: ", _branch_terminal_state(b), ", ")
    end
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
    print(io, "tap_ratio: ", b.tap_ratio, ", ")
    print(io, "phase_shift_deg: ", b.phase_shift_deg, ", ")

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

"""
    _open_terminal_yin(branch) -> ComplexF64

Exact pi-model input admittance of a one-sided open branch seen from its
closed bus, as the Schur complement of the two-port from `calcAdmittance`:
with the TO end open `Y_in = Y11 - Y12*Y21/Y22`, with the FROM end open
`Y_in = Y22 - Y21*Y12/Y11`. The two-port already carries the complex ratio,
so lines and transformers (off-nominal ratio, phase shift) are covered
uniformly without a case distinction. For `|Y_s| >> |Y_0|` this approaches
the FULL line charging `g + jb` (not half of it).
"""
function _open_terminal_yin(branch::Branch)::ComplexF64
  Y11, Y12, Y21, Y22 = calcAdmittance(branch, 1.0, 1.0)
  st = _branch_terminal_state(branch)
  if st == :open_to
    abs(Y22) < 1e-12 && return 0.0 + 0.0im
    return Y11 - Y12 * Y21 / Y22
  elseif st == :open_from
    abs(Y11) < 1e-12 && return 0.0 + 0.0im
    return Y22 - Y21 * Y12 / Y11
  end
  error("_open_terminal_yin: branch $(branch.branchIdx) is $(st), not one-sided open")
end

"""
    _open_end_voltage(branch, u_closed::ComplexF64) -> ComplexF64

Voltage at the open terminal of a one-sided open branch from the pi-model
voltage divider (zero current at the open end): with TO open
`U_open = -Y21/Y22 * U_from`, with FROM open `U_open = -Y12/Y11 * U_to`.
Reproduces the Ferranti rise (`|U_open| > |U_closed|` for b > 0) without
adding a node to the solved system.
"""
function _open_end_voltage(branch::Branch, u_closed::ComplexF64)::ComplexF64
  Y11, Y12, Y21, Y22 = calcAdmittance(branch, 1.0, 1.0)
  st = _branch_terminal_state(branch)
  if st == :open_to
    abs(Y22) < 1e-12 && return u_closed
    return -Y21 / Y22 * u_closed
  elseif st == :open_from
    abs(Y11) < 1e-12 && return u_closed
    return -Y12 / Y11 * u_closed
  end
  error("_open_end_voltage: branch $(branch.branchIdx) is $(st), not one-sided open")
end

# helper
function setBranchFlow!(branch::Branch, tfBranchFlow::BranchFlow, fBranchFlow::BranchFlow)
  branch.tBranchFlow = tfBranchFlow
  branch.fBranchFlow = fBranchFlow
end

"""
    _branch_terminal_state(br::Branch) -> Symbol

Single source of truth for the service state of a branch: `:closed` (both
terminals in service), `:open_from` / `:open_to` (exactly one terminal
open, the branch reduces to its pi-model Schur complement at the closed
bus), or `:open` (out of service: aggregate `status == 0` with untouched
terminal flags, or both flags open). Every consumer reads the state
through this helper, never the raw fields.
"""
function _branch_terminal_state(br::Branch)::Symbol
  # aggregate forced open without touching the terminal flags (legacy
  # direct writes of br.status): both flags still 1 means fully open
  (br.status == 0 && br.from_status == 1 && br.to_status == 1) && return :open
  br.from_status == 0 && br.to_status == 0 && return :open
  br.from_status == 0 && return :open_from
  br.to_status == 0 && return :open_to
  return :closed
end

"""
    setBranchStatus!(branch, service::Bool)

User-facing aggregate switch: sets the aggregate `status` and BOTH
terminal flags consistently (in service = all closed, out of service =
all open).
"""
function setBranchStatus!(branch::Branch, service::Bool)
  v = service ? 1 : 0
  branch.status = v
  branch.from_status = v
  branch.to_status = v
  return branch
end

"""
    setBranchTerminalStatus!(branch; from = nothing, to = nothing)

Open or close individual branch terminals (`true` = closed, `false` =
open; `nothing` leaves a terminal unchanged) and recompute the aggregate
`status` (1 iff both terminals closed). A branch open at exactly one
terminal stays in the model as its exact pi reduction at the closed bus,
see the "One-sided open branches" section of the branch-model docs.
"""
function setBranchTerminalStatus!(branch::Branch; from::Union{Nothing,Bool} = nothing, to::Union{Nothing,Bool} = nothing)
  from === nothing || (branch.from_status = from ? 1 : 0)
  to === nothing || (branch.to_status = to ? 1 : 0)
  branch.status = (branch.from_status == 1 && branch.to_status == 1) ? 1 : 0
  return branch
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
  ratio = (branch.ratio == 0.0) ? 1.0 : branch.tap_ratio
  shift = (branch.ratio == 0.0) ? 0.0 : branch.phase_shift_deg

  # use ratio even if shift is 0
  if ratio != 1.0 || shift != 0.0
    return calcComplexRatio(tapRatio = ratio, angleInDegrees = shift)
  else
    return 1.0 + 0.0im
  end
end

function getBranchComp(Vn_kV::Float64, from::Int, to::Int, idx::Int, kind::String)
  cTyp = toComponentTyp("Branch")
  name = "B_$(kind)_$(string(convert(Int,trunc(Vn_kV))))_$(Int(from))_$(Int(to))"
  cID = "#" * name * "#" * string(idx)
  return ImpPGMComp(cID, name, cTyp, Vn_kV, from, to)
end
