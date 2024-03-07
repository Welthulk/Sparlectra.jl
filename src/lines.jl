# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file lines.jl

mutable struct ACLineSegment <: AbstractBranch
  comp::AbstractComponent
  length::Float64
  r::Float64
  x::Float64
  b::Union{Nothing,Float64}
  g::Union{Nothing,Float64}
  c_nf_per_km::Union{Nothing,Float64}
  tanδ::Union{Nothing,Float64}

  function ACLineSegment(c::AbstractComponent, length::Float64, r::Float64, x::Float64, c_nf_per_km::Union{Nothing,Float64} = nothing, tanδ::Union{Nothing,Float64} = nothing)
    g = 0.0
    b = 0.0
    if !isnothing(c_nf_per_km) && !isnothing(tanδ)
      y1_shunt_ = 2.0 * pi * 50.0 * c_nf_per_km * 1e-9 * (tanδ + im * 1.0)
      g = real(y1_shunt_)
      b = imag(y1_shunt_)
    end
    new(c, length, r, x, b, g, c_nf_per_km, tanδ)
  end

  function ACLineSegment(id::String, name::String, Vn::Float64, length::Float64, r::Float64, x::Float64, b::Union{Nothing,Float64} = nothing, g::Union{Nothing,Float64} = nothing, c_nf_per_km::Union{Nothing,Float64} = nothing)
    comp = Component(id, name, ResDataTypes.LineC, Vn)
    new(comp, length, r, x, b, g, c_nf_per_km, nothing)
  end

  function Base.show(io::IO, acseg::ACLineSegment)
    print(io, "ACLineSegment( ")
    print(io, acseg.comp)
    print(io, " ")
    print(io, "length: $(acseg.length), ")
    print(io, "r: $(acseg.r), ")
    print(io, "x: $(acseg.x), ")
    if !isnothing(acseg.g)
      print(io, "g: $(acseg.g), ")
    end
    if !isnothing(acseg.b)
      print(io, "b: $(acseg.b), ")
    end
    if !isnothing(acseg.c_nf_per_km)
      print(io, "c_nf_per_km: $(acseg.c_nf_per_km), ")
    end
    if !isnothing(acseg.tanδ)
      print(io, "tanδ: $(acseg.tanδ), ")
    end
    print(io, ")")
  end
end

function getRXBG(o::ACLineSegment)::Tuple{Float64,Float64,Union{Nothing,Float64},Union{Nothing,Float64}}
  return (o.r, o.x, o.b, o.g)
end

function get_line_parameters(line::ACLineSegment)
  parameters = Dict{Symbol,Any}()

  for field in fieldnames(ACLineSegment)
    value = getproperty(line, field)
    parameters[field] = value === nothing ? 0.0 : value
  end

  return parameters
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
