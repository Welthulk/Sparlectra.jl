# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file shunt.jl

mutable struct Shunt
  comp::AbstractComponent  
  busIdx::Int              # short cut to node    
  p_shunt::Float64         # aktive power
  q_shunt::Float64         # reactive power      
  y_pu_shunt::ComplexF64   # shunt admittance in p.u.
  status::Int              # 1 = in service, 0 = out of service

  function Shunt(comp::AbstractComponent, busIdx::Int, p_shunt::Float64, q_shunt::Float64, y_pu_shunt::ComplexF64, status::Int = 1)
    new(comp, busIdx, p_shunt, q_shunt, y_pu_shunt, status)
  end

  function Shunt(;comp::ImpPGMComp, base_MVA::Float64, Vn_kV_shunt::Float64, g_shunt::Float64, b_shunt::Float64, ratio::Float64=1.0, status::Int = 1)
    busIdx = comp.cFrom_bus
    p_shunt = g_shunt * Vn_kV_shunt^2*ratio
    q_shunt = b_shunt * Vn_kV_shunt^2*ratio
    Sh_ref = Complex(p_shunt, q_shunt)
    y_pu_shunt = Sh_ref / base_MVA

    new(comp, busIdx, p_shunt, q_shunt, y_pu_shunt, status)
  end

  function Base.show(io::IO, shunt::Shunt)
    print(io, "Shunt( ")
    print(io, shunt.comp, ", ")    
    print(io, "bus: ", shunt.busIdx, ", ")
    print(io, "p_shunt: ", shunt.p_shunt, " MW, ")
    print(io, "q_shunt: ", shunt.q_shunt, " MVar, ")
    print(io, "y_pu_shunt: ", shunt.y_pu_shunt, " p.u., ")
    print(io, "status: ", shunt.status, ")")
  end
end

function getShuntPGMComp(Vn::Float64, from::Int, id::Int)
  cTyp = toComponentTyp("SHUNT")
  cName = "Sh_$(string(round(Vn,digits=1)))"
  cID = "#$cName\\_$from\\_#$(string(id))"
  return ImpPGMComp(cID, cName, cTyp, Vn, from, from)  
end

