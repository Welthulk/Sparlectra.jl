# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file lines.jl

mutable struct ACLineSegment
  comp::AbstractComponent
  length::Float64
  r::Float64
  x::Float64
  b::Union{Nothing,Float64}
  g::Union{Nothing,Float64}
  c_nf_per_km::Union{Nothing,Float64}
  tanδ::Union{Nothing,Float64}

  function ACLineSegment(c::AbstractComponent, length::Float64, r::Float64, x::Float64, b::Union{Nothing,Float64} = nothing, g::Union{Nothing,Float64} = nothing, c_nf_per_km::Union{Nothing,Float64} = nothing, tanδ::Union{Nothing,Float64} = nothing)
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
    if !isnothing(acseg.b)
      print(io, "b: $(acseg.b), ")
    end
    if !isnothing(acseg.g)
      print(io, "g: $(acseg.g), ")
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

function getRXBG(o::ACLineSegment)::Tuple{Float64, Float64, Union{Nothing,Float64}, Union{Nothing,Float64}}
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
