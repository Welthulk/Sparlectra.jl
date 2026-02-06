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
- `vn_kV::Float64`: The nominal voltage of the shunt in kV.
- `baseMVA::Float64`: The base MVA for per unit calculations.
- `busIdx::Int`: The index of the bus.
- `p_shunt::Float64`: The active power of the shunt (MW) - results after solve.
- `q_shunt::Float64`: The reactive power of the shunt (MVar) - results after solve.
- `G_shunt::Float64`: The conductance of the shunt (S) - optional input/debug.
- `B_shunt::Float64`: The susceptance of the shunt (S) - optional input/debug.
- `y_pu_shunt::ComplexF64`: The shunt admittance in per unit (solver-relevant for Y-shunts).
- `model::Symbol`: The shunt model - `:Y` for admittance or `:PQ` for constant power.
- `status::Int`: The status of the shunt. 1 = in service, 0 = out of service.

# Constructors
- `Shunt(; fromBus::Int, id::Int, base_MVA::Float64, vn_kV_shunt::Float64, p_shunt::Union{Nothing,Float64} = nothing, q_shunt::Union{Nothing,Float64} = nothing, g_shunt::Union{Nothing,Float64} = nothing, b_shunt::Union{Nothing,Float64} = nothing, y_pu::Union{Nothing,ComplexF64} = nothing, values_are_pu::Bool = false, model::Symbol = :Y, ratio::Float64 = 1.0, status::Int = 1)`: Creates a new `Shunt` instance. Provide either `y_pu` OR (`g_shunt`, `b_shunt`) OR (`p_shunt`, `q_shunt`).

# Methods
- `Base.show(io::IO, shunt::Shunt)`: Prints the `Shunt` instance.
"""

mutable struct Shunt
  comp::AbstractComponent
  vn_kV::Float64
  baseMVA::Float64
  busIdx::Int

  # Report / results (MW/MVar) after solve
  p_shunt::Float64
  q_shunt::Float64

  # Optional (input helper / debug)
  G_shunt::Float64
  B_shunt::Float64

  # Solver-relevant for Y-shunts
  y_pu_shunt::ComplexF64

  # Semantics
  model::Symbol      # :Y  or :PQ   (start with :Y if you want)
  status::Int

  function Shunt(; fromBus::Int, id::Int, base_MVA::Float64, vn_kV_shunt::Float64,
                 p_shunt::Union{Nothing,Float64}=nothing,
                 q_shunt::Union{Nothing,Float64}=nothing,
                 g_shunt::Union{Nothing,Float64}=nothing,
                 b_shunt::Union{Nothing,Float64}=nothing,
                 y_pu::Union{Nothing,ComplexF64}=nothing,
                 values_are_pu::Bool=false,
                 model::Symbol=:Y,
                 ratio::Float64=1.0,
                 status::Int=1)

    comp   = getShuntPGMComp(vn_kV_shunt, fromBus, id)
    busIdx = fromBus

    # --- build y_pu_shunt depending on what was provided ---
    if y_pu !== nothing
      y_pu_shunt = ComplexF64(y_pu)
      G = real(y_pu_shunt); B = imag(y_pu_shunt)

    elseif (g_shunt !== nothing || b_shunt !== nothing)
      @assert (g_shunt !== nothing && b_shunt !== nothing) "Provide both g_shunt and b_shunt"
      if values_are_pu
        # g/b already in pu-admittance on system base
        y_pu_shunt = ComplexF64(g_shunt, b_shunt)
      else
        # g/b in Siemens -> convert to pu: y_pu = y_phys / y_base, y_base = baseMVA / V^2
        # => y_pu = (g + j b) * V^2 / baseMVA
        y_pu_shunt = ComplexF64(g_shunt, b_shunt) * (vn_kV_shunt^2) / base_MVA * ratio
      end
      G = real(y_pu_shunt); B = imag(y_pu_shunt)

    elseif (p_shunt !== nothing || q_shunt !== nothing)
      @assert (p_shunt !== nothing && q_shunt !== nothing) "Provide both p_shunt and q_shunt"
      if model == :Y
        # interpret p/q as MATPOWER-style (MW/MVar at V=1pu) => y_pu = (P + jQ)/baseMVA
        y_pu_shunt = ComplexF64(p_shunt, q_shunt) / base_MVA
        G = real(y_pu_shunt); B = imag(y_pu_shunt)
      else
        # model == :PQ (constant power) -> NOT a Y-shunt; keep y_pu_shunt = 0 here
        y_pu_shunt = 0.0 + 0.0im
        G = 0.0; B = 0.0
      end
    else
      error("Provide either y_pu OR (g_shunt,b_shunt) OR (p_shunt,q_shunt).")
    end

    # Results init (will be overwritten by post-solve update)
    p0 = 0.0
    q0 = 0.0

    return new(comp, vn_kV_shunt, base_MVA, busIdx, p0, q0, G, B, y_pu_shunt, model, status)
 end

 function Base.show(io::IO, shunt::Shunt)
    print(io, "Shunt( ")
    print(io, shunt.comp, ", ")
    print(io, "bus: ", shunt.busIdx, ", ")
    print(io, "model: ", shunt.model, ", ")
    print(io, "status: ", shunt.status, ", ")
    print(io, "p_shunt: ", shunt.p_shunt, " MW, ")
    print(io, "q_shunt: ", shunt.q_shunt, " MVar ")
    print(io, "G_shunt: ", shunt.G_shunt, " S, ")
    print(io, "B_shunt: ", shunt.B_shunt, " S ")
    print(IO, "y_pu_shunt: ", shunt.y_pu_shunt)
    print(io, ")")
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
