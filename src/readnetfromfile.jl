# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 01.8.2023
# include-file readnetfromfile.jl
# Bus
struct Bus
  busIdx::Int64
  name::String
  id::String
  wt3id::String
  vn_kv::Float64
  type::Int64
  pgmID::Union{Nothing,Int64} # original id uses for pgm import

  function Bus(busIdx::Int64, name::String, id::String, wt3id::String, vn_kv::Float64, type::Int64)
    new(busIdx, name, id, wt3id, vn_kv, type, nothing)
  end

  function Bus(busIdx::Int64, name::String, id::String, wt3id::String, vn_kv::Float64, type::Int64, pgmID::Int64)
    new(busIdx, name, id, wt3id, vn_kv, type, pgmID)
  end

  function Base.show(io::IO, bus::Bus)
    if bus.pgmID === nothing
      print(io, "Bus(", bus.busIdx, ", ", bus.name, ", ", bus.id, ", ", bus.wt3id, ", ", bus.vn_kv, ", ", bus.type, ")")
    else
      print(io, "Bus(", bus.busIdx, ", ", bus.name, ", ", bus.id, ", ", bus.wt3id, ", ", bus.vn_kv, ", ", bus.type, ", ", bus.pgmID, ")")
    end    
  end
end

function readTapMod(tap::Dict{String,Any})
  side = tap["tap_side"]
  tap_min = Int(tap["tap_min"])
  tap_max = Int(tap["tap_max"])
  tap_neutral = Int(tap["tap_neutral"])
  tap_step_percent = float(tap["tap_step_percent"])
  shift_degree = float(tap["shift_degree"])
  neutralU_ratio = float(tap["neutralU_ratio"])

  return side, tap_min, tap_max, tap_neutral, tap_step_percent, shift_degree, neutralU_ratio
end

function checkBusNumber(bus::Int64, busVec::Vector{Bus})::Bool
  for b in busVec
    if b.busIdx == bus
      return true
    end
  end
  @warn "bus $bus not found!"
  return false
end

include("jsonimport.jl")
include("pgmimport.jl ")


