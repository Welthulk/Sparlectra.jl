# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file shunt.jl

mutable struct Shunt
  comp::AbstractComponent
  busIdx::Int              # short cut to node    
  p_shunt::Float64         # aktive power
  q_shunt::Float64         # reactive power
  G_shunt::Float64         # conductance
  B_shunt::Float64         # susceptance      
  y_pu_shunt::ComplexF64   # shunt admittance in p.u.
  status::Int              # 1 = in service, 0 = out of service

  function Shunt(; fromBus::Int, id::Int, base_MVA::Float64, Vn_kV_shunt::Float64, p_shunt::Union{Nothing,Float64} = nothing, q_shunt::Union{Nothing,Float64} = nothing, g_shunt::Union{Nothing,Float64} = nothing, b_shunt::Union{Nothing,Float64} = nothing, ratio::Float64 = 1.0, status::Int = 1)
    comp = getShuntPGMComp(Vn_kV_shunt, fromBus, id)
    busIdx = fromBus
    if (isnothing(p_shunt) && isnothing(q_shunt)) && (isnothing(g_shunt) && isnothing(b_shunt))
      error("Either p_shunt and q_shunt or g_shunt and b_shunt must be given")
    end

    if isnothing(p_shunt) && isnothing(q_shunt)
      p_shunt = g_shunt * Vn_kV_shunt^2 * ratio
      q_shunt = b_shunt * Vn_kV_shunt^2 * ratio
      y_pu_shunt = Complex(p_shunt, q_shunt) / base_MVA
      new(comp, busIdx, p_shunt, q_shunt, g_shunt, b_shunt, y_pu_shunt, status)
    else
      g_shunt = p_shunt / Vn_kV_shunt^2
      b_shunt = q_shunt / Vn_kV_shunt^2
      y_pu_shunt = Complex(g_shunt, q_shunt) / base_MVA
      new(comp, busIdx, p_shunt, q_shunt, g_shunt, b_shunt, y_pu_shunt, status)
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

function getShuntPGMComp(Vn::Float64, from::Int, id::Int)
  cTyp = toComponentTyp("SHUNT")
  cName = "Sh_$(from)_$(string(round(Vn,digits=1)))"
  cID = "#$cName\\_#$(string(id))"
  return ImpPGMComp(cID, cName, cTyp, Vn, from, from)
end
