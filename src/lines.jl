# Copyright 2023-2026 Udo Schmitz
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
# file: src/lines.jl
# purpose: ACLineSegment branch type: physical line totals in ohm and siemens
#          or per-unit PI-model data, plus accessor and show methods

"""
    ACLineSegment

An AC line between two buses of one voltage level.

A physical line (built by [`addACLine!`](@ref) or [`addACLinePerKm!`](@ref))
stores the totals of the whole line: `r` and `x` in ohm, `b` and `g` in
siemens. A PI-model line (built by [`addPIModelACLine!`](@ref)) stores per-unit
values on the bus voltage and the network base. [`getLineRXBG_pu`](@ref)
returns per unit for both kinds.

# Fields
- `comp::AbstractComponent`: the component record.
- `length_km::Union{Nothing,Float64}`: line length in km, information only
  (used by `show` and as CGMES `Conductor.length` in the export); it never
  scales `r`, `x`, `b` or `g`.
- `r::Float64`, `x::Float64`: series resistance and reactance, ohm (physical)
  or per unit (PI model).
- `b::Float64`, `g::Float64`: total shunt susceptance and conductance,
  siemens (physical) or per unit (PI model).
- `ratedS::Union{Nothing,Float64}`: rated power in MVA.
- `_isPIModel::Bool`: `true` when the values are per unit.

# Constructor
    ACLineSegment(; vn_kv, from, to, r, x, b = 0.0, g = 0.0, length_km = nothing, ratedS = nothing, isPIModel)

`isPIModel` has no default: the caller states the unit of the four values.
"""
mutable struct ACLineSegment <: AbstractBranch
  comp::AbstractComponent
  length_km::Union{Nothing,Float64}
  r::Float64
  x::Float64
  b::Float64
  g::Float64
  ratedS::Union{Nothing,Float64}
  _isPIModel::Bool
  function ACLineSegment(; vn_kv::Float64, from::Int, to::Int, r::Float64, x::Float64, b::Float64 = 0.0, g::Float64 = 0.0, length_km::Union{Nothing,Float64} = nothing, ratedS::Union{Nothing,Float64} = nothing, isPIModel::Bool)
    c = getLineImpPGMComp(vn_kv, from, to)
    new(c, length_km, r, x, b, g, ratedS, isPIModel)
  end

  function Base.show(io::IO, acseg::ACLineSegment)
    unit_z, unit_y = acseg._isPIModel ? ("pu", "pu") : ("Ohm", "S")
    print(io, "ACLineSegment( ")
    print(io, acseg.comp)
    print(io, " ")
    isnothing(acseg.length_km) || print(io, "length_km: $(acseg.length_km), ")
    print(io, "r: $(acseg.r) $(unit_z), ")
    print(io, "x: $(acseg.x) $(unit_z), ")
    print(io, "b: $(acseg.b) $(unit_y), ")
    print(io, "g: $(acseg.g) $(unit_y), ")
    if !isnothing(acseg.ratedS)
      print(io, "ratedS: $(acseg.ratedS), ")
    end
    print(io, ")")
  end
end

"""
    getLineRXBG(o::ACLineSegment) -> NTuple{4,Float64}

The stored values `(r, x, b, g)` of the line: totals in ohm and siemens for a
physical line, per unit for a PI-model line (see [`isLinePIModel`](@ref)).
Use [`getLineRXBG_pu`](@ref) for per-unit values of either kind.
"""
function getLineRXBG(o::ACLineSegment)::NTuple{4,Float64}
  return (o.r, o.x, o.b, o.g)
end

"""
    getLineRXBG_pu(o, vn_kV, baseMVA) -> (r, x, b, g)

The line parameters in per unit on `vn_kV` and `baseMVA`. A physical line is
converted with [`toPU_RXBG`](@ref); a PI-model line already holds per unit
and is returned unchanged.
"""
function getLineRXBG_pu(o::ACLineSegment, vn_kV::Float64, baseMVA::Float64)::NTuple{4,Float64}
  o._isPIModel && return getLineRXBG(o)
  return toPU_RXBG(r = o.r, x = o.x, g = o.g, b = o.b, v_kv = vn_kV, baseMVA = baseMVA)
end

"""
    get_line_parameters(line::ACLineSegment)::Dict{Symbol,Any}

Returns a dictionary of the parameters of an AC line segment. If a parameter is `nothing`, it is replaced with `0.0`.

# Arguments
- `line::ACLineSegment`: The AC line segment.

# Returns
- `parameters::Dict{Symbol,Any}`: A dictionary where the keys are the parameter names and the values are the parameter values.

# Example
```julia
get_line_parameters(acLineSegment)
```
"""
function get_line_parameters(line::ACLineSegment)
  parameters = Dict{Symbol,Any}()

  for field in fieldnames(ACLineSegment)
    value = getproperty(line, field)
    parameters[field] = value === nothing ? 0.0 : value
  end

  return parameters
end

"""
    isLinePIModel(line) -> Bool

Whether the line was built directly from per-unit PI-model parameters.
"""
function isLinePIModel(line::ACLineSegment)
  return line._isPIModel
end

function getLineBusID(Vn::Float64, from::Int, to::Int)
  name = "ACL_$(string(round(Vn,digits=1)))"
  id = "#$name\\_$from\\_$to"
  return name, id
end

function getLineImpPGMComp(Vn::Float64, from::Int, to::Int)
  cName, cID = getLineBusID(Vn, from, to)
  return ImpPGMComp(cID, cName, toComponentTyp("ACLINESEGMENT"), Vn, from, to)
end
