# Copyright 2023-2026 Udo Schmitz
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

# file: src/dynamic_data.jl
# purpose: optional dynamic parameters of a generating unit (#506 part 1):
#          DynamicData, the unit type enum and the validation error. A
#          ProSumer refers to one through its `dynamics` field. No
#          calculation reads them; the CGMES island slack choice reads the
#          unit type only.

"""
    DynamicUnitType

How a unit holds the grid frequency:

- `SynchronousUnit`: a synchronous machine, inertia from its rotating mass.
- `GridFollowingUnit`: a converter that follows the grid voltage and needs
  another unit to form it. It cannot hold an island alone.
- `GridFormingUnit`: a converter that forms voltage and frequency itself,
  with a synthetic inertia (`virtual_inertia_s` of [`DynamicData`](@ref)).

SCF writes the values as `"synchronous"`, `"grid_following"`,
`"grid_forming"` ([`dynamic_unit_type_name`](@ref)).
"""
@enum DynamicUnitType SynchronousUnit GridFollowingUnit GridFormingUnit

"""
    dynamic_unit_type_name(t::DynamicUnitType) -> String

The SCF spelling of a unit type: `"synchronous"`, `"grid_following"` or
`"grid_forming"`.
"""
function dynamic_unit_type_name(t::DynamicUnitType)::String
  t == SynchronousUnit && return "synchronous"
  t == GridFollowingUnit && return "grid_following"
  return "grid_forming"
end

"""
    parse_dynamic_unit_type(s::AbstractString) -> DynamicUnitType

The unit type of its SCF spelling ([`dynamic_unit_type_name`](@ref)).
Throws `ArgumentError` naming the allowed values for anything else.
"""
function parse_dynamic_unit_type(s::AbstractString)::DynamicUnitType
  s == "synchronous" && return SynchronousUnit
  s == "grid_following" && return GridFollowingUnit
  s == "grid_forming" && return GridFormingUnit
  throw(ArgumentError("unit_type must be \"synchronous\", \"grid_following\" or \"grid_forming\", got $(repr(s))"))
end

"""
    DynamicDataError(field, value, rule)

Thrown by the [`DynamicData`](@ref) constructor when a value breaks the rule
of its field. `field` names the field, `value` is the rejected value and
`rule` the condition it fails (for example `"must be >= 0"`).
"""
struct DynamicDataError <: Exception
  field::Symbol
  value::Float64
  rule::String
end

Base.showerror(io::IO, e::DynamicDataError) = print(io, "DynamicDataError: ", e.field, " = ", e.value, " ", e.rule)

# one Float64 field: finite, and -0.0 stored as 0.0 (the SCF writer prints
# every zero as 0.0, so a negative zero would not come back bit for bit)
function _dyn_value(field::Symbol, x::Float64)::Float64
  isfinite(x) || throw(DynamicDataError(field, x, "must be finite"))
  return x + 0.0
end
_dyn_value(::Symbol, ::Nothing)::Nothing = nothing

function _dyn_nonnegative(field::Symbol, x::Union{Nothing,Float64})::Union{Nothing,Float64}
  v = _dyn_value(field, x)
  (v === nothing || v >= 0.0) || throw(DynamicDataError(field, v, "must be >= 0"))
  return v
end

function _dyn_positive(field::Symbol, x::Union{Nothing,Float64})::Union{Nothing,Float64}
  v = _dyn_value(field, x)
  (v === nothing || v > 0.0) || throw(DynamicDataError(field, v, "must be > 0 when set"))
  return v
end

"""
    DynamicData(; inertia_s = nothing, damping_pu = 0.0, rating_mva = nothing,
                droop_pu = nothing, deadband_hz = 0.0, primary_reserve_mw = 0.0,
                governor_time_s = nothing, agc_participation = 0.0,
                secondary_reserve_mw = 0.0, unit_type = nothing,
                virtual_inertia_s = nothing, black_start = false,
                min_stable_mw = nothing)

Dynamic parameters of one generating unit, referenced from
`ProSumer.dynamics` (`nothing` when the unit has none). No power flow, state
estimation, short circuit or contingency calculation reads them. The only
reader is the island slack choice of the CGMES import, which takes
`unit_type` into account.

# Fields
- `inertia_s::Union{Nothing,Float64}`: inertia constant H in s, on the rating
  below. `nothing`: the unit has no inertia.
- `damping_pu::Float64`: damping coefficient D in pu/pu. Default 0.0.
- `rating_mva::Union{Nothing,Float64}`: rating in MVA that H and D refer to,
  when it differs from the unit's `ratedS`. `nothing`: `ratedS` applies.
- `droop_pu::Union{Nothing,Float64}`: primary control droop (statics) in pu,
  for example 0.05. `nothing`: no primary control.
- `deadband_hz::Float64`: primary control deadband in Hz. Default 0.0.
- `primary_reserve_mw::Float64`: primary reserve limit in MW, symmetric.
  Default 0.0.
- `governor_time_s::Union{Nothing,Float64}`: governor/turbine time constant
  in s of a first-order lag. `nothing`: not given.
- `agc_participation::Float64`: share of the secondary control signal, 0 to
  1. Default 0.0.
- `secondary_reserve_mw::Float64`: secondary reserve limit in MW,
  symmetric. Default 0.0.
- `unit_type::DynamicUnitType`: [`DynamicUnitType`](@ref). Default
  `SynchronousUnit` when `inertia_s` is set, else `GridFollowingUnit`.
- `virtual_inertia_s::Union{Nothing,Float64}`: synthetic inertia in s of a
  grid-forming unit. `nothing`: none.
- `black_start::Bool`: the unit can energise a dead network. Default false.
- `min_stable_mw::Union{Nothing,Float64}`: minimum stable generation in MW.
  `nothing`: not given.

# Validation
Every number must be finite; inertia, virtual inertia, time constant,
deadband and both reserves must not be negative; droop and rating must be
greater than zero when set; `agc_participation` lies in [0, 1]. A violation
throws [`DynamicDataError`](@ref) naming the field. A negative zero is
stored as 0.0.

# Example
```julia
ps.dynamics = DynamicData(inertia_s = 4.0, droop_pu = 0.05, primary_reserve_mw = 20.0)
```
"""
struct DynamicData
  inertia_s::Union{Nothing,Float64}
  damping_pu::Float64
  rating_mva::Union{Nothing,Float64}
  droop_pu::Union{Nothing,Float64}
  deadband_hz::Float64
  primary_reserve_mw::Float64
  governor_time_s::Union{Nothing,Float64}
  agc_participation::Float64
  secondary_reserve_mw::Float64
  unit_type::DynamicUnitType
  virtual_inertia_s::Union{Nothing,Float64}
  black_start::Bool
  min_stable_mw::Union{Nothing,Float64}

  function DynamicData(inertia_s::Union{Nothing,Float64}, damping_pu::Float64, rating_mva::Union{Nothing,Float64}, droop_pu::Union{Nothing,Float64}, deadband_hz::Float64, primary_reserve_mw::Float64, governor_time_s::Union{Nothing,Float64}, agc_participation::Float64, secondary_reserve_mw::Float64, unit_type::DynamicUnitType, virtual_inertia_s::Union{Nothing,Float64}, black_start::Bool, min_stable_mw::Union{Nothing,Float64})
    agc = _dyn_value(:agc_participation, agc_participation)
    0.0 <= agc <= 1.0 || throw(DynamicDataError(:agc_participation, agc, "must lie in [0, 1]"))
    return new(
      _dyn_nonnegative(:inertia_s, inertia_s),
      _dyn_value(:damping_pu, damping_pu),
      _dyn_positive(:rating_mva, rating_mva),
      _dyn_positive(:droop_pu, droop_pu),
      _dyn_nonnegative(:deadband_hz, deadband_hz),
      _dyn_nonnegative(:primary_reserve_mw, primary_reserve_mw),
      _dyn_nonnegative(:governor_time_s, governor_time_s),
      agc,
      _dyn_nonnegative(:secondary_reserve_mw, secondary_reserve_mw),
      unit_type,
      _dyn_nonnegative(:virtual_inertia_s, virtual_inertia_s),
      black_start,
      _dyn_value(:min_stable_mw, min_stable_mw),
    )
  end
end

# keyword form with the documented defaults; integers and other reals are
# taken as Float64, the unit type follows the inertia when not given
_dyn_float(::Nothing) = nothing
_dyn_float(x::Real) = Float64(x)

function DynamicData(;
  inertia_s::Union{Nothing,Real} = nothing,
  damping_pu::Real = 0.0,
  rating_mva::Union{Nothing,Real} = nothing,
  droop_pu::Union{Nothing,Real} = nothing,
  deadband_hz::Real = 0.0,
  primary_reserve_mw::Real = 0.0,
  governor_time_s::Union{Nothing,Real} = nothing,
  agc_participation::Real = 0.0,
  secondary_reserve_mw::Real = 0.0,
  unit_type::Union{Nothing,DynamicUnitType} = nothing,
  virtual_inertia_s::Union{Nothing,Real} = nothing,
  black_start::Bool = false,
  min_stable_mw::Union{Nothing,Real} = nothing,
)
  ut = unit_type === nothing ? (inertia_s === nothing ? GridFollowingUnit : SynchronousUnit) : unit_type
  return DynamicData(_dyn_float(inertia_s), Float64(damping_pu), _dyn_float(rating_mva), _dyn_float(droop_pu), Float64(deadband_hz), Float64(primary_reserve_mw), _dyn_float(governor_time_s), Float64(agc_participation), Float64(secondary_reserve_mw), ut, _dyn_float(virtual_inertia_s), black_start, _dyn_float(min_stable_mw))
end

function Base.show(io::IO, d::DynamicData)
  print(io, "DynamicData(", dynamic_unit_type_name(d.unit_type))
  for f in fieldnames(DynamicData)
    f === :unit_type && continue
    v = getfield(d, f)
    v === nothing && continue
    print(io, ", ", f, " = ", v)
  end
  print(io, ")")
end
