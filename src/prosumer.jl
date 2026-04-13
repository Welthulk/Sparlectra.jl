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
# file: src/prosumer.jl

abstract type AbstractVoltageDependentController end

"""
    PiecewiseLinearCharacteristic(points)

Piecewise linear characteristic `y = f(u)` represented by ordered `(u, y)` points.
Outside the point range, values are clamped to the edge points and the derivative is `0.0`.
"""
struct PiecewiseLinearCharacteristic
  points::Vector{Tuple{Float64,Float64}}
  interpolation::Symbol
  spline_second_derivatives::Vector{Float64}
  polynomial_nodes::Vector{Float64}
  polynomial_coefficients::Vector{Float64}

  function PiecewiseLinearCharacteristic(points::Vector{Tuple{Float64,Float64}}; interpolation::Symbol = :linear)
    length(points) >= 2 || error("PiecewiseLinearCharacteristic requires at least 2 points.")
    u_prev = points[1][1]
    for i in 2:length(points)
      u = points[i][1]
      u > u_prev || error("PiecewiseLinearCharacteristic points must have strictly increasing voltage values.")
      u_prev = u
    end

    if interpolation === :linear
      return new(points, interpolation, Float64[], Float64[], Float64[])
    elseif interpolation === :spline
      if length(points) == 2
        return new(points, :linear, Float64[], Float64[], Float64[])
      end
      return new(points, interpolation, _natural_spline_second_derivatives(points), Float64[], Float64[])
    elseif interpolation === :polynomial
      if length(points) == 2
        return new(points, :linear, Float64[], Float64[], Float64[])
      end
      nodes, coeffs = _newton_polynomial_coefficients(points)
      return new(points, interpolation, Float64[], nodes, coeffs)
    else
      error("PiecewiseLinearCharacteristic: unsupported interpolation $(interpolation). Use :linear, :spline, or :polynomial.")
    end
  end
end

function _natural_spline_second_derivatives(points::Vector{Tuple{Float64,Float64}})
  n = length(points)
  h = Vector{Float64}(undef, n - 1)
  for i in 1:(n-1)
    h[i] = points[i+1][1] - points[i][1]
  end

  lower = Vector{Float64}(undef, n - 3)
  diag = Vector{Float64}(undef, n - 2)
  upper = Vector{Float64}(undef, n - 3)
  rhs = Vector{Float64}(undef, n - 2)

  for i in 1:(n-2)
    him1 = h[i]
    hi = h[i+1]
    diag[i] = 2.0 * (him1 + hi)
    rhs[i] = 6.0 * (((points[i+2][2] - points[i+1][2]) / hi) - ((points[i+1][2] - points[i][2]) / him1))
    if i <= n - 3
      upper[i] = hi
      lower[i] = hi
    end
  end

  m = zeros(Float64, n)
  m[2:n-1] = Tridiagonal(lower, diag, upper) \ rhs
  return m
end

function _newton_polynomial_coefficients(points::Vector{Tuple{Float64,Float64}})
  n = length(points)
  nodes = [p[1] for p in points]
  coeffs = [p[2] for p in points]

  for j in 2:n
    for i in n:-1:j
      coeffs[i] = (coeffs[i] - coeffs[i-1]) / (nodes[i] - nodes[i-j+1])
    end
  end

  return nodes, coeffs
end

function _evaluate_newton_polynomial(nodes::Vector{Float64}, coeffs::Vector{Float64}, x::Float64)
  n = length(coeffs)
  value = coeffs[n]
  slope = 0.0

  for i in (n-1):-1:1
    slope = slope * (x - nodes[i]) + value
    value = value * (x - nodes[i]) + coeffs[i]
  end

  return value, slope
end

struct QUController <: AbstractVoltageDependentController
  characteristic::PiecewiseLinearCharacteristic
  qmin_pu::Union{Nothing,Float64}
  qmax_pu::Union{Nothing,Float64}
end

struct PUController <: AbstractVoltageDependentController
  characteristic::PiecewiseLinearCharacteristic
  pmin_pu::Union{Nothing,Float64}
  pmax_pu::Union{Nothing,Float64}
end

@inline _to_pu_power(value::Float64, sbase_MVA::Float64) = value / sbase_MVA
@inline _to_pu_voltage(value::Float64, vn_kV::Float64) = value / vn_kV

function _convert_power_limit_to_pu(name::String, value::Union{Nothing,Float64}, unit::Symbol, sbase_MVA::Union{Nothing,Float64})
  isnothing(value) && return nothing
  if unit === :pu
    return value
  elseif (unit === :MW) || (unit === :MVAr)
    isnothing(sbase_MVA) && error("$(name): sbase_MVA is required when using $(unit) limits.")
    return _to_pu_power(value, sbase_MVA)
  else
    error("$(name): unsupported unit $(unit). Use :pu, :MW, or :MVAr.")
  end
end

"""
    QUController(characteristic; qmin_pu=nothing, qmax_pu=nothing,
                 qmin_MVAr=nothing, qmax_MVAr=nothing, sbase_MVA=nothing)

Convenience constructor for `QUController` that accepts limits either in p.u. or in MVAr.
"""
function QUController(characteristic::PiecewiseLinearCharacteristic; qmin_pu::Union{Nothing,Float64} = nothing, qmax_pu::Union{Nothing,Float64} = nothing, qmin_MVAr::Union{Nothing,Float64} = nothing, qmax_MVAr::Union{Nothing,Float64} = nothing, sbase_MVA::Union{Nothing,Float64} = nothing)
  (!isnothing(qmin_pu) && !isnothing(qmin_MVAr)) && error("QUController: specify either qmin_pu or qmin_MVAr, not both.")
  (!isnothing(qmax_pu) && !isnothing(qmax_MVAr)) && error("QUController: specify either qmax_pu or qmax_MVAr, not both.")
  qmin = isnothing(qmin_MVAr) ? qmin_pu : _convert_power_limit_to_pu("QUController", qmin_MVAr, :MVAr, sbase_MVA)
  qmax = isnothing(qmax_MVAr) ? qmax_pu : _convert_power_limit_to_pu("QUController", qmax_MVAr, :MVAr, sbase_MVA)
  return QUController(characteristic, qmin, qmax)
end

"""
    PUController(characteristic; pmin_pu=nothing, pmax_pu=nothing,
                 pmin_MW=nothing, pmax_MW=nothing, sbase_MVA=nothing)

Convenience constructor for `PUController` that accepts limits either in p.u. or in MW.
"""
function PUController(characteristic::PiecewiseLinearCharacteristic; pmin_pu::Union{Nothing,Float64} = nothing, pmax_pu::Union{Nothing,Float64} = nothing, pmin_MW::Union{Nothing,Float64} = nothing, pmax_MW::Union{Nothing,Float64} = nothing, sbase_MVA::Union{Nothing,Float64} = nothing)
  (!isnothing(pmin_pu) && !isnothing(pmin_MW)) && error("PUController: specify either pmin_pu or pmin_MW, not both.")
  (!isnothing(pmax_pu) && !isnothing(pmax_MW)) && error("PUController: specify either pmax_pu or pmax_MW, not both.")
  pmin = isnothing(pmin_MW) ? pmin_pu : _convert_power_limit_to_pu("PUController", pmin_MW, :MW, sbase_MVA)
  pmax = isnothing(pmax_MW) ? pmax_pu : _convert_power_limit_to_pu("PUController", pmax_MW, :MW, sbase_MVA)
  return PUController(characteristic, pmin, pmax)
end

@inline function _apply_limits(value_pu::Float64, slope_pu::Float64, min_pu::Union{Nothing,Float64}, max_pu::Union{Nothing,Float64})
  if !isnothing(min_pu) && value_pu < min_pu
    return min_pu, 0.0
  end
  if !isnothing(max_pu) && value_pu > max_pu
    return max_pu, 0.0
  end
  return value_pu, slope_pu
end

"""
    evaluate_characteristic(ch, u_pu) -> (value, slope)

Evaluate piecewise linear characteristic at `u_pu`.

For interior breakpoints (turning points), the implementation uses the segment
found first during the left-to-right scan, i.e. the slope of the segment on
the left side of the breakpoint. Outside the point range, the value is clamped
and the returned slope is `0.0`.
"""
function evaluate_characteristic(ch::PiecewiseLinearCharacteristic, u_pu::Float64)
  pts = ch.points
  u1 = pts[1][1]
  uN = pts[end][1]

  if u_pu <= u1
    return pts[1][2], 0.0
  elseif u_pu >= uN
    return pts[end][2], 0.0
  end

  if ch.interpolation === :polynomial
    return _evaluate_newton_polynomial(ch.polynomial_nodes, ch.polynomial_coefficients, u_pu)
  end

  for i in 1:(length(pts)-1)
    uk, yk = pts[i]
    uk1, yk1 = pts[i+1]
    if uk <= u_pu <= uk1
      if ch.interpolation === :spline
        h = uk1 - uk
        m = ch.spline_second_derivatives
        a = (uk1 - u_pu) / h
        b = (u_pu - uk) / h
        value = a * yk + b * yk1 + (((a^3 - a) * m[i] + (b^3 - b) * m[i+1]) * h^2 / 6.0)
        slope = (yk1 - yk) / h + (((-3.0 * a^2 + 1.0) * m[i] + (3.0 * b^2 - 1.0) * m[i+1]) * h / 6.0)
        return value, slope
      else
        slope = (yk1 - yk) / (uk1 - uk)
        value = yk + slope * (u_pu - uk)
        return value, slope
      end
    end
  end

  return pts[end][2], 0.0
end

function evaluate_controller(ctrl::QUController, u_pu::Float64)
  value, slope = evaluate_characteristic(ctrl.characteristic, u_pu)
  return _apply_limits(value, slope, ctrl.qmin_pu, ctrl.qmax_pu)
end

function evaluate_controller(ctrl::PUController, u_pu::Float64)
  value, slope = evaluate_characteristic(ctrl.characteristic, u_pu)
  return _apply_limits(value, slope, ctrl.pmin_pu, ctrl.pmax_pu)
end

"""
    make_characteristic(points; voltage_unit=:pu, value_unit=:pu, vn_kV=nothing, sbase_MVA=nothing, interpolation=:linear)

Create a piecewise linear characteristic from points in p.u. or physical units.

- `voltage_unit = :pu` expects voltage in per-unit.
- `voltage_unit = :kV` expects voltage in kV and requires `vn_kV`.
- `value_unit = :pu` expects power output in per-unit.
- `value_unit = :MW` or `:MVAr` expects physical power values and requires `sbase_MVA`.
- `interpolation = :linear` uses piecewise linear interpolation.
- `interpolation = :spline` uses natural cubic spline interpolation through all points.
  If only two points are provided, linear interpolation is used automatically.
- `interpolation = :polynomial` uses one global polynomial through all points.
  If only two points are provided, linear interpolation is used automatically.
"""
function make_characteristic(points::Vector{Tuple{Float64,Float64}}; voltage_unit::Symbol = :pu, value_unit::Symbol = :pu, vn_kV::Union{Nothing,Float64} = nothing, sbase_MVA::Union{Nothing,Float64} = nothing, interpolation::Symbol = :linear)
  converted = Tuple{Float64,Float64}[]
  for (u, y) in points
    u_pu = if voltage_unit === :pu
      u
    elseif voltage_unit === :kV
      isnothing(vn_kV) && error("make_characteristic: vn_kV is required when voltage_unit=:kV.")
      _to_pu_voltage(u, vn_kV)
    else
      error("make_characteristic: unsupported voltage_unit $(voltage_unit). Use :pu or :kV.")
    end

    y_pu = if value_unit === :pu
      y
    elseif (value_unit === :MW) || (value_unit === :MVAr)
      isnothing(sbase_MVA) && error("make_characteristic: sbase_MVA is required when value_unit=$(value_unit).")
      _to_pu_power(y, sbase_MVA)
    else
      error("make_characteristic: unsupported value_unit $(value_unit). Use :pu, :MW, or :MVAr.")
    end
    push!(converted, (u_pu, y_pu))
  end
  return PiecewiseLinearCharacteristic(converted; interpolation = interpolation)
end

# Data type to describe producers and consumers
"""
    ProSumer

A mutable structure representing a prosumer in a power system. A prosumer is an entity that either produces or consumes power.

# Fields
- `comp::AbstractComponent`: The component of the prosumer.
- `busIdx::Int`: The index of the bus where the prosumer is connected.
- `pGen::Float64`: The active power generation of the prosumer.
- `qGen::Float64`: The reactive power generation of the prosumer.
- `pLoad::Float64`: The active power consumption of the prosumer.
- `qLoad::Float64`: The reactive power consumption of the prosumer.
- `status::Int`: The status of the prosumer. 1 = in service, 0 = out of service.

# Constructors
- `ProSumer(comp::AbstractComponent, busIdx::Int, pGen::Float64, qGen::Float64, pLoad::Float64, qLoad::Float64, status::Int)`: Creates a new `ProSumer` instance.

# Methods
- `Base.show(io::IO, prosumer::ProSumer)`: Prints the `ProSumer` instance.

# Example
```julia
prosumer = ProSumer(comp, 1, 100.0, 50.0, 80.0, 40.0, 1)
```
"""
mutable struct ProSumer
  comp::AbstractComponent
  ratedS::Union{Nothing,Float64}
  ratedU::Union{Nothing,Float64}
  qPercent::Union{Nothing,Float64}
  pVal::Union{Nothing,Float64}
  qVal::Union{Nothing,Float64}
  maxP::Union{Nothing,Float64}
  minP::Union{Nothing,Float64}
  maxQ::Union{Nothing,Float64}
  minQ::Union{Nothing,Float64}
  ratedPowerFactor::Union{Nothing,Float64}
  referencePri::Union{Nothing,Integer}
  vm_pu::Union{Nothing,Float64}
  va_deg::Union{Nothing,Float64}
  vstep_pu::Union{Nothing,Float64}
  tap_steps_down::Union{Nothing,Int}
  tap_steps_up::Union{Nothing,Int}
  isRegulated::Bool
  proSumptionType::ProSumptionType
  isAPUNode::Bool
  quController::Union{Nothing,QUController}
  puController::Union{Nothing,PUController}
  qGenRepl::Union{Nothing,Float64}
  pRes::Union{Nothing,Float64}
  qRes::Union{Nothing,Float64}

  function ProSumer(;
    vn_kv::Float64,
    busIdx::Int,
    oID::Int,
    type::ProSumptionType,
    ratedS::Union{Nothing,Float64} = nothing,
    ratedU::Union{Nothing,Float64} = nothing,
    qPercent::Union{Nothing,Float64} = nothing,
    p::Union{Nothing,Float64} = nothing,
    q::Union{Nothing,Float64} = nothing,
    maxP::Union{Nothing,Float64} = nothing,
    minP::Union{Nothing,Float64} = nothing,
    maxQ::Union{Nothing,Float64} = nothing,
    minQ::Union{Nothing,Float64} = nothing,
    ratedPowerFactor::Union{Nothing,Float64} = nothing,
    referencePri::Union{Nothing,Integer} = nothing,
    vm_pu::Union{Nothing,Float64} = nothing,
    va_deg::Union{Nothing,Float64} = nothing,
    vstep_pu::Union{Nothing,Float64} = nothing,
    tap_steps_down::Union{Nothing,Int} = nothing,
    tap_steps_up::Union{Nothing,Int} = nothing,
    isRegulated::Bool = false,
    isAPUNode::Bool = false,
    quController::Union{Nothing,QUController} = nothing,
    puController::Union{Nothing,PUController} = nothing,
  )
    comp = getProSumPGMComp(vn_kv, busIdx, isGenerator(type), oID)

    if isnothing(vm_pu)
      vm_pu = 1.0
    end

    if isnothing(va_deg)
      va_deg = 0.0
    end

    new(comp, ratedS, ratedU, qPercent, p, q, maxP, minP, maxQ, minQ, ratedPowerFactor, referencePri, vm_pu, va_deg, vstep_pu, tap_steps_down, tap_steps_up, isRegulated, type, isAPUNode, quController, puController, nothing, nothing, nothing)
  end

  function Base.show(io::IO, prosumption::ProSumer)
    print(io, "ProSumption( ")
    print(io, prosumption.comp, ", ")

    if (!isnothing(prosumption.ratedS))
      print(io, "ratedS: ", prosumption.ratedS, " MVA, ")
    end

    if (!isnothing(prosumption.ratedU))
      print(io, "ratedU: ", prosumption.ratedU, " kV, ")
    end

    if (!isnothing(prosumption.qPercent))
      print(io, "qPercent: ", prosumption.qPercent, " %, ")
    end

    if (!isnothing(prosumption.pVal))
      print(io, "pVal: ", prosumption.pVal, " MW, ")
    end

    if (!isnothing(prosumption.qVal))
      print(io, "qVal: ", prosumption.qVal, " MVar, ")
    end

    if (!isnothing(prosumption.maxP))
      print(io, "maxP: ", prosumption.maxP, " MW, ")
    end

    if (!isnothing(prosumption.minP))
      print(io, "minP: ", prosumption.minP, " MW, ")
    end

    if (!isnothing(prosumption.maxQ))
      print(io, "maxQ: ", prosumption.maxQ, " MVar, ")
    end

    if (!isnothing(prosumption.minQ))
      print(io, "minQ: ", prosumption.minQ, " MVar, ")
    end

    if (!isnothing(prosumption.ratedPowerFactor))
      print(io, "ratedPowerFactor: ", prosumption.ratedPowerFactor, ", ")
    end

    if (!isnothing(prosumption.referencePri))
      print(io, "referencePri: ", prosumption.referencePri, ", ")
    end

    if (!isnothing(prosumption.vm_pu))
      print(io, "vm_pu: ", prosumption.vm_pu, ", ")
    end

    if (!isnothing(prosumption.va_deg))
      print(io, "va_deg: ", prosumption.va_deg, ", ")
    end
    if (!isnothing(prosumption.vstep_pu))
      print(io, "vstep_pu: ", prosumption.vstep_pu, ", ")
    end
    if (!isnothing(prosumption.tap_steps_down))
      print(io, "tap_steps_down: ", prosumption.tap_steps_down, ", ")
    end
    if (!isnothing(prosumption.tap_steps_up))
      print(io, "tap_steps_up: ", prosumption.tap_steps_up, ", ")
    end
    if prosumption.isRegulated
      print(io, "isRegulated: true, ")
    end
    if !isnothing(prosumption.quController)
      print(io, "quController: yes, ")
    end
    if !isnothing(prosumption.puController)
      print(io, "puController: yes, ")
    end

    if (!isnothing(prosumption.qGenRepl))
      print(io, "qGenRepl: ", prosumption.qGenRepl, ", ")
    end

    if (!isnothing(prosumption.pRes))
      print(io, "pRes: ", prosumption.pRes, " MW, ")
    end
    if (!isnothing(prosumption.qRes))
      print(io, "qRes: ", prosumption.qRes, " MVar, ")
    end

    print(io, "proSumptionType: ", prosumption.proSumptionType, ")")
  end
end

function setPQResult!(ps::ProSumer, p::Float64, q::Float64)
  ps.pRes = p
  ps.qRes = q
end

function getPosumerBusIndex(ps::ProSumer)::Int
  c = ps.comp
  if hasproperty(c, :cFrom_bus) && getfield(c, :cFrom_bus) !== nothing
    return Int(getfield(c, :cFrom_bus))
  elseif hasproperty(c, :cTo_bus) && getfield(c, :cTo_bus) !== nothing
    return Int(getfield(c, :cTo_bus))
  end
  error("ProSumer: cannot determine bus index (component has neither :cFrom_bus nor :cTo_bus).")
end

function isAPUNode(o::ProSumer)
  return o.isAPUNode
end

function isGenerator(c::Sparlectra.ProSumptionType)::Bool
  if c == Injection
    return true
  else
    return false
  end
end # isGenerator

function isGenerator(o::ProSumer)::Bool
  return isGenerator(o.proSumptionType)
end # isGenerator

@inline has_qu_controller(ps::ProSumer)::Bool = !isnothing(ps.quController)
@inline has_pu_controller(ps::ProSumer)::Bool = !isnothing(ps.puController)

function getProSumPGMComp(Vn::Float64, from::Int, isGen::Bool, id::Int)
  cTyp = isGen ? toComponentTyp("GENERATOR") : toComponentTyp("LOAD")
  cName = isGen ? "Gen_$(string(convert(Int,trunc(Vn))))" : "Ld_$(string(round(Vn,digits=1)))"
  cID = "#$cName\\_$from\\_#$(string(id))"
  return ImpPGMComp(cID, cName, cTyp, Vn, from, from)
end

function setQGenReplacement!(o::ProSumer, q::Float64)
  o.qGenRepl = q
end

function getQGenReplacement(o::ProSumer)::Union{Nothing,Float64}
  return o.qGenRepl
end

function updatePQ!(o::ProSumer, p::Union{Nothing,Float64}, q::Union{Nothing,Float64})
  if !isnothing(p)
    o.pVal = p
  end
  if !isnothing(q)
    o.qVal = q
  end
end

function isSlack(o::ProSumer)
  if !isnothing(o.referencePri) && o.referencePri > 0
    return true
  else
    return false
  end
end

function isRegulating(o::ProSumer)::Bool
  has_controller = !isnothing(o.vstep_pu) || !isnothing(o.tap_steps_down) || !isnothing(o.tap_steps_up)
  return o.isRegulated || has_controller
end

function toProSumptionType(o::ComponentTyp)::ProSumptionType
  if o == Generator || o == ExternalNetworkInjection || o == SynchronousMachine
    return Injection
  elseif o == AsynchronousMachine || o == Load || o == EnergyConsumer
    return Consumption
  else
    return UnknownP
  end
end

function toProSumptionType(o::String)::ProSumptionType
  val = uppercase(o)
  if val == "GENERATOR" || val == "EXTERNALNETWORKINJECTION" || val == "SYNCHRONOUSMACHINE"
    return Injection
  elseif val == "ASYNCHRONOUSMACHINE" || val == "LOAD" || val == "ENERGYCONSUMER"
    return Consumption
  else
    return UnknownP
  end
end

function toString(o::ProSumptionType)::String
  if o == Injection
    return "Injection"
  elseif o == Consumption
    return "Consumption"
  else
    return "UnknownP"
  end
end
