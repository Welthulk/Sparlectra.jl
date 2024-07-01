# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# include-file lines.jl

"""
    ACLineSegment

A mutable structure representing an AC line segment in a power system.

# Fields
- `comp::AbstractComponent`: The component of the AC line segment.
- `length::Float64`: The length of the AC line segment.
- `r::Float64`: The resistance of the AC line segment.
- `x::Float64`: The reactance of the AC line segment.
- `b::Union{Nothing,Float64}`: The susceptance of the AC line segment. It can be `Nothing` or a `Float64` value.
- `g::Union{Nothing,Float64}`: The conductance of the AC line segment. It can be `Nothing` or a `Float64` value.
- `c_nf_per_km::Union{Nothing,Float64}`: The capacitance per kilometer of the AC line segment. It can be `Nothing` or a `Float64` value.
- `tanδ::Union{Nothing,Float64}`: The tangent of the loss angle of the AC line segment. It can be `Nothing` or a `Float64` value.
- `ratedS::Union{Nothing,Float64}`: The rated power of the AC line segment. It can be `Nothing` or a `Float64` value.
- `paramsBasedOnLength::Bool`: A boolean indicating whether the parameters are based on the length of the AC line segment, if false, the parameters (r,x,g,b) have to multiply by length.
- `_isPIModel::Bool`: A boolean indicating whether the AC line segment is a PI model, all parameters are p.u. .

# Constructors
- `ACLineSegment(; vn_kv::Float64, from::Int, to::Int, length::Float64, r::Float64, x::Float64, b::Union{Nothing,Float64} = nothing, c_nf_per_km::Union{Nothing,Float64} = nothing, 
                           tanδ::Union{Nothing,Float64} = nothing, ratedS::Union{Nothing,Float64} = nothing, paramsBasedOnLength=true)`: Creates a new `ACLineSegment` instance.

# Methods
- `Base.show(io::IO, acseg::ACLineSegment)`: Prints the `ACLineSegment` instance.
"""
mutable struct ACLineSegment <: AbstractBranch
  comp::AbstractComponent
  length::Float64
  r::Float64
  x::Float64
  b::Union{Nothing,Float64}
  g::Union{Nothing,Float64}
  c_nf_per_km::Union{Nothing,Float64}
  tanδ::Union{Nothing,Float64}
  ratedS::Union{Nothing,Float64}
  paramsBasedOnLength::Bool
  _isPIModel::Bool
  isLongLine::Bool
  #! format: off
  function ACLineSegment(; vn_kv::Float64,  from::Int,  to::Int,  length::Float64, r::Float64, x::Float64, b::Union{Nothing,Float64} = nothing,
                           c_nf_per_km::Union{Nothing,Float64} = nothing, tanδ::Union{Nothing,Float64} = nothing, ratedS::Union{Nothing,Float64} = nothing,
                           paramsBasedOnLength::Bool = true, isPIModel::Bool = false, isLongLine::Bool = false)
  #! format: on
    c = getLineImpPGMComp(vn_kv, from, to)
    g = 0.0

    if isPIModel
      new(c, length, r, x, b, g, nothing, nothing, ratedS, paramsBasedOnLength, true)
    else
      if !isnothing(b)
        b = b
      elseif !isnothing(c_nf_per_km) && !isnothing(tanδ)
        y1_shunt_ = 2.0 * pi * 50.0 * c_nf_per_km * 1e-9 * (tanδ + im * 1.0)
        g = real(y1_shunt_)
        b = imag(y1_shunt_)
      else
        b = 0.0
      end

      new(c, length, r, x, b, g, c_nf_per_km, tanδ, ratedS, paramsBasedOnLength, isPIModel, isLongLine)
    end
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
    if !isnothing(acseg.ratedS)
      print(io, "ratedS: $(acseg.ratedS), ")
    end
    print(io, ")")
  end
end

function getABCDParms(branch::ACLineSegment, u_rated::Float64, s_rated::Float64)
  r,x,g,b = getRXBG(branch)  
  r_pu, x_pu, g_pu, b_pu = toPU_RXGB(r = r, x = x, g = g, b = b, v_kv = u_rated, baseMVA = s_rated)
  if branch.isLongLine
    Z = Complex(r_pu, x_pu)
    Y = Complex(g_pu, b_pu)
    
    γ = sqrt(Z * Y)          # Propagation constant    
    Z_c = sqrt(Z / Y)        # Characteristic impedance
    
    A = cosh(γ)
    B = Z_c * sinh(γ)
    C = sinh(γ) / Z_c
    D = cosh(γ)
  else    
    Z = Complex(r_pu, x_pu)
    Y = Complex(g_pu, b_pu)
        
    A = 1 + Z * Y / 2
    B = Z
    C = Y * (1 + Z * Y / 4)
    D = 1 + Z * Y / 2
  end
  return A, B, C, D    
end


"""
    getRXBG(o::ACLineSegment)::Tuple{Float64,Float64,Union{Nothing,Float64},Union{Nothing,Float64}}

Returns the resistance, reactance, susceptance, and conductance of an AC line segment. 
If the parameters are based on length, they are multiplied by the length of the line segment.

# Arguments
- `o::ACLineSegment`: The AC line segment.

# Returns
- `r::Float64`: The resistance of the AC line segment.
- `x::Float64`: The reactance of the AC line segment.
- `b::Union{Nothing,Float64}`: The susceptance of the AC line segment. It can be `Nothing` or a `Float64` value.
- `g::Union{Nothing,Float64}`: The conductance of the AC line segment. It can be `Nothing` or a `Float64` value.

# Example
```julia
getRXBG(acLineSegment)
```
"""
function getRXBG(o::ACLineSegment)::Tuple{Float64,Float64,Union{Nothing,Float64},Union{Nothing,Float64}}
  if o.paramsBasedOnLength || o._isPIModel
    return (o.r, o.x, o.b, o.g)
  else
    return (o.r * o.length, o.x * o.length, o.b * o.length, o.g * o.length)
  end
end

"""
    get_line_parameters(line::ACLineSegment)::Dict{Symbol,Any}

Returns a dictionary of the parameters of an AC line segment. If a parameter is `nothing`, it is replaced with `0.0`.

# Arguments
- `line::ACLineSegment`: The AC line segment.

# Returns
- `parameters::Dict{Symbol,Any}`: A dictionary where the keys are the parameter names and the values are the parameter values.

# Example
```julia
get_line_parameters(acLineSegment)
```
"""
function get_line_parameters(line::ACLineSegment)
  parameters = Dict{Symbol,Any}()

  for field in fieldnames(ACLineSegment)
    value = getproperty(line, field)
    parameters[field] = value === nothing ? 0.0 : value
  end

  return parameters
end

function isPIModel(line::ACLineSegment)
  return line._isPIModel
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


#=
using LinearAlgebra

function calculate_PI_ABCD_parameters(r, x, g, b, d)
    # Serienschaltung und Paralleladmittanzen
    Z = Complex(r, x) * d
    Y = Complex(g, b) * d
    
    # Berechnung der ABCD-Parameter
    A = 1 + Z * Y / 2
    B = Z
    C = Y * (1 + Z * Y / 4)
    D = 1 + Z * Y / 2
    
    return A, B, C, D
end

# Beispielparameter für das PI-Ersatzschaltbild
r = 0.1       # Ω/km
x = 0.314     # Ω/km
g = 0.0       # S/km
b = 3.14e-3   # S/km
d = 100.0     # km

A, B, C, D = calculate_PI_ABCD_parameters(r, x, g, b, d)

println("A: ", A)
println("B: ", B)
println("C: ", C)
println("D: ", D)





using LinearAlgebra

function calculate_ABCD_parameters(r, l, g, c, d, f)
    ω = 2π * f
    Z = Complex(r, ω * l)
    Y = Complex(g, ω * c)
    
    γ = sqrt(Z * Y)
    Z_c = sqrt(Z / Y)
    
    A = cosh(γ * d)
    B = Z_c * sinh(γ * d)
    C = sinh(γ * d) / Z_c
    D = cosh(γ * d)
    
    return A, B, C, D
end

# Beispielparameter
r = 0.1
l = 0.001
g = 0.0
c = 0.01e-6
d = 100.0
f = 50.0

A, B, C, D = calculate_ABCD_parameters(r, l, g, c, d, f)

println("A: ", A)
println("B: ", B)
println("C: ", C)
println("D: ", D)

using LinearAlgebra

function calculate_ABCD_parameters(r, x, g, b, d)
    Z = Complex(r, x)
    Y = Complex(g, b)
    
    γ = sqrt(Z * Y)
    Z_c = sqrt(Z / Y)
    
    A = cosh(γ * d)
    B = Z_c * sinh(γ * d)
    C = sinh(γ * d) / Z_c
    D = cosh(γ * d)
    
    return A, B, C, D
end

# Beispielparameter
r = 0.1       # Ω/km
x = 0.314     # Ω/km
g = 0.0       # S/km
b = 3.14e-3   # S/km
d = 100.0     # km

A, B, C, D = calculate_ABCD_parameters(r, x, g, b, d)

println("A: ", A)
println("B: ", B)
println("C: ", C)
println("D: ", D)

=#