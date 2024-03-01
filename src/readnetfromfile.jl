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
  vm_pu::Union{Nothing,Float64}
  va_deg::Union{Nothing,Float64}

  function Bus(busIdx::Int64, vn_kv::Float64, type::Int64, vm_pu::Float64=1.0, va_deg::Float64=0.0)    
   name = "Bus_$(string(busIdx))_$(string(convert(Int,trunc(vn_kv))))"
   id = "#$name"
    
   new(busIdx, name, id, "", vn_kv, type, vm_pu, va_deg)
  end
  function Bus(busIdx::Int64, name::String, id::String, wt3id::String, vn_kv::Float64, type::Int64)
    new(busIdx, name, id, wt3id, vn_kv, type, nothing, nothing)
  end

  function Bus(busIdx::Int64, name::String, id::String, wt3id::String, vn_kv::Float64, type::Int64, vm_pu::Union{Nothing,Float64}, va_deg::Union{Nothing,Float64})
    new(busIdx, name, id, wt3id, vn_kv, type, vm_pu, va_deg)
  end

  function Base.show(io::IO, bus::Bus)
    print(io, "Bus(", bus.busIdx, ", ", bus.name, ", ", bus.id, ", ", bus.wt3id, ", ", bus.vn_kv, ", ", bus.type, ")")
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


