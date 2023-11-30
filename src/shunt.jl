# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file shunt.jl

struct Shunt
  comp::Component
  nodeID::String           # ID for mapping to node
  busIdx::Int              # short cut to node    
  p_shunt::Float64         # aktive power
  q_shunt::Float64         # reactive power      
  y_pu_shunt::ComplexF64   # shunt admittance in p.u.
  status::Int              # 1 = in service, 0 = out of service

  function Shunt(comp::Component, nodeID::String, busIdx::Int, p_shunt::Float64, q_shunt::Float64, y_pu_shunt::ComplexF64, status::Int = 1)
    new(comp, nodeID, busIdx, p_shunt, q_shunt, y_pu_shunt, status)
  end

  function Base.show(io::IO, shunt::Shunt)
    print(io, "Shunt( ")
    print(io, shunt.comp, ", ")
    print(io, "NodeID: ", shunt.nodeID, ", ")
    print(io, "bus: ", shunt.busIdx, ", ")
    print(io, "p_shunt: ", shunt.p_shunt, " MW, ")
    print(io, "q_shunt: ", shunt.q_shunt, " MVar, ")
    print(io, "y_pu_shunt: ", shunt.y_pu_shunt, " p.u., ")
    print(io, "status: ", shunt.status, ")")
  end
end
