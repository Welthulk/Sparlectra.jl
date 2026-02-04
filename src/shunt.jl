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
# file: src/shunt.jl
"""
    Shunt

A mutable structure representing a shunt in a power system.

# Fields
- `comp::AbstractComponent`: The component of the shunt.
- `busIdx::Int`: The index of the bus.
- `p_shunt::Float64`: The active power of the shunt.
- `q_shunt::Float64`: The reactive power of the shunt.
- `G_shunt::Float64`: The conductance of the shunt.
- `B_shunt::Float64`: The susceptance of the shunt.
- `y_pu_shunt::ComplexF64`: The shunt admittance in per unit.
- `status::Int`: The status of the shunt. 1 = in service, 0 = out of service.

# Constructors
- `Shunt(; fromBus::Int, id::Int, base_MVA::Float64, vn_kV_shunt::Float64, p_shunt::Union{Nothing,Float64} = nothing, q_shunt::Union{Nothing,Float64} = nothing, g_shunt::Union{Nothing,Float64} = nothing, b_shunt::Union{Nothing,Float64} = nothing, ratio::Float64 = 1.0, status::Int = 1)`: Creates a new `Shunt` instance.

# Methods
- `Base.show(io::IO, shunt::Shunt)`: Prints the `Shunt` instance.
"""
mutable struct Shunt
  comp::AbstractComponent
  vn_kV::Float64           # rated voltage in kV
  baseMVA::Float64         # base MVA
  busIdx::Int              # short cut to node    
  p_shunt::Float64         # aktive power
  q_shunt::Float64         # reactive power
  G_shunt::Float64         # conductance
  B_shunt::Float64         # susceptance      
  y_pu_shunt::ComplexF64   # shunt admittance in p.u.
  status::Int              # 1 = in service, 0 = out of service

  function Shunt(; fromBus::Int, id::Int, base_MVA::Float64, vn_kV_shunt::Float64, p_shunt::Union{Nothing,Float64} = nothing, q_shunt::Union{Nothing,Float64} = nothing, g_shunt::Union{Nothing,Float64} = nothing, b_shunt::Union{Nothing,Float64} = nothing, ratio::Float64 = 1.0, status::Int = 1)
    comp = getShuntPGMComp(vn_kV_shunt, fromBus, id)
    busIdx = fromBus
    if (isnothing(p_shunt) && isnothing(q_shunt)) && (isnothing(g_shunt) && isnothing(b_shunt))
      error("Either p_shunt and q_shunt or g_shunt and b_shunt must be given")
    end

    if isnothing(p_shunt) && isnothing(q_shunt)
      p_shunt = g_shunt * vn_kV_shunt^2 * ratio
      q_shunt = b_shunt * vn_kV_shunt^2 * ratio
      y_pu_shunt = Complex(g_shunt, b_shunt) * base_MVA
      new(comp, vn_kV_shunt, base_MVA, busIdx, p_shunt, q_shunt, g_shunt, b_shunt, y_pu_shunt, status)
    else
      g_shunt = p_shunt / vn_kV_shunt^2
      b_shunt = q_shunt / vn_kV_shunt^2
      y_pu_shunt = Complex(p_shunt, q_shunt) / base_MVA
      new(comp, vn_kV_shunt, base_MVA, busIdx, p_shunt, q_shunt, g_shunt, b_shunt, y_pu_shunt, status)
    end
  end

  function Base.show(io::IO, shunt::Shunt)
    print(io, "Shunt( ")
    print(io, shunt.comp, ", ")
    print(io, "bus: ", shunt.busIdx, ", ")
    print(io, "p_shunt: ", shunt.p_shunt, " MW, ")
    print(io, "q_shunt: ", shunt.q_shunt, " MVar, ")
    print(io, "G_shunt: ", shunt.G_shunt, " S, ")
    print(io, "B_shunt: ", shunt.B_shunt, " S, ")
    print(io, "y_pu_shunt: ", shunt.y_pu_shunt, " p.u., ")
    print(io, "status: ", shunt.status, ")")
  end
end

function getGBShunt(o::Shunt)::Tuple{Float64,Float64}
  return (o.G_shunt, o.B_shunt)
end

function getPQShunt(o::Shunt)::Tuple{Float64,Float64}
  return (o.p_shunt, o.q_shunt)
end

function updatePQShunt!(o::Shunt, p::Float64, q::Float64)
  o.p_shunt = p
  o.q_shunt = q
  o.G_shunt = p / o.vn_kV^2
  o.B_shunt = q / o.vn_kV^2
  o.y_pu_shunt = Complex(o.p_shunt, o.q_shunt) / o.baseMVA
end

function getShuntPGMComp(Vn::Float64, from::Int, id::Int)
  cTyp = toComponentTyp("SHUNT")
  cName = "Sh_$(from)_$(string(round(Vn,digits=1)))"
  cID = "#$cName\\_#$(string(id))"
  return ImpPGMComp(cID, cName, cTyp, Vn, from, from)
end
