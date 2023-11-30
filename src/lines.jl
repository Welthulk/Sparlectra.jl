# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file lines.jl

mutable struct ACLineSegment
  comp::Component
  length::Float64
  r::Float64
  x::Float64
  b::Union{Nothing,Float64}
  g::Union{Nothing,Float64}
  maxIka::Union{Nothing,Float64}
  c_nf_per_km::Union{Nothing,Float64}

  function ACLineSegment(c::Component, length::Float64, r::Float64, x::Float64, b::Float64)
    new(c, length, r, x, b, nothing, nothing, nothing)
  end

  function ACLineSegment(id::String, name::String, Vn::Float64, length::Float64, r::Float64, x::Float64, b::Union{Nothing,Float64} = nothing, g::Union{Nothing,Float64} = nothing, maxIka::Union{Nothing,Float64} = nothing, c_nf_per_km::Union{Nothing,Float64} = nothing)
    comp = Component(id, name, ResDataTypes.LineC, Vn)
    new(comp, length, r, x, b, g, maxIka, c_nf_per_km)
  end

  function Base.show(io::IO, acseg::ACLineSegment)
    print(io, "ACLineSegment( ")
    print(io, acseg.comp)
    print(io, " ")
    print(io, "length: $(acseg.length), ")
    print(io, "r: $(acseg.r), ")
    print(io, "x: $(acseg.x), ")
    if !isnothing(acseg.b)
      print(io, "b: $(acseg.b), ")
    end
    if !isnothing(acseg.g)
      print(io, "g: $(acseg.g), ")
    end
    if !isnothing(acseg.c_nf_per_km)
      print(io, "c_nf_per_km: $(acseg.c_nf_per_km), ")
    end
    if !isnothing(acseg.maxIka)
      print(io, "maxIka: $(acseg.maxIka), ")
    end
    print(io, ")")
  end
end

