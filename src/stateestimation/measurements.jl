# Copyright 2023–2026 Udo Schmitz
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

# file: src/stateestimation/measurements.jl
# purpose: state-estimation measurement model: MeasurementType and Measurement,
#          add helpers for Vm/Va/PMU/injection/flow/zero-injection measurements,
#          and synthetic measurement generation from power-flow results

using Random

@inline function _net_measurements(net::Net)::Vector
  return net.measurements
end

"""
    @enum MeasurementType

Supported measurement types for the first WLS state-estimation implementation.

`ImagMeas` is a current-magnitude measurement in ampere, in two variants:
the branch variant (`branchIdx` + `direction`) and the bus-referenced
shunt-bay variant (`busIdx` only, SE phase 2) measuring the current of the
shunt bay at that bus. Current magnitudes never carry observability: the
observability checks exclude them, and the estimator gates them by
iteration count and by the `value >= 3 sigma` rule (see `runse!`).

`ShuntQMeas` is the reactive power of a shunt bay in MVar (`busIdx`
referenced, sign convention identical to the solved `q_shunt` results). It
is the direct measurement that releases a shunt susceptance for estimation
(SE phase 2, case A) and the carrier of the derived pseudo-measurements
from the shunt back-calculation (case B, `deriveShuntPseudoMeasurements!`).

`IaMeas` is a PMU current-phasor ANGLE in degrees at a branch end (or the
shunt bay, same variants as `ImagMeas`), referenced to the PMU time base:
the common reference-angle offset alpha applies exactly as for `VaMeas`.
A current angle is meaningless near zero current, so an `IaMeas` row is
active only while the predicted current magnitude passes
`I >= 3 * sigma_I_ref` (the paired `ImagMeas` sigma at the same end, else
the config floor `state_estimation.ia_current_floor_A`); it shares the
`ImagMeas` iteration gate and, like `ImagMeas`, never carries
observability in this stage (a linear PMU-only estimation that would use
current phasors for observability is future work).
"""
@enum MeasurementType begin
  VmMeas
  PinjMeas
  QinjMeas
  PflowMeas
  QflowMeas
  VaMeas
  ImagMeas
  ShuntQMeas
  IaMeas
end

"""
    Measurement

Generic state-estimation measurement model.

Fields:
- `typ`: Measurement type.
- `value`: Measured value (Vm in p.u., powers in MW/MVar, Va in degrees).
- `sigma`: Standard deviation in measurement units.
- `weight`: Weight used in WLS (`1/sigma^2`).
- `active`: If `false`, measurement is ignored by estimator.
- `busIdx`: Bus index for bus measurements.
- `branchIdx`: Branch index for branch flow measurements.
- `direction`: Branch direction `:from` or `:to`, otherwise `:none`.
- `linkIdx`: Bus-link index for link flow measurements (SE phase 3). A
  link-referenced flow measurement never enters the WLS estimator: the link
  is not part of the Ybus, so it constrains only the post-SE flow
  allocation (`calcLinkFlowsSE!`).
- `id`: Optional measurement identifier.
"""
struct Measurement
  typ::MeasurementType
  value::Float64
  sigma::Float64
  weight::Float64
  active::Bool
  busIdx::Union{Nothing,Int}
  branchIdx::Union{Nothing,Int}
  direction::Symbol
  id::String
  linkIdx::Union{Nothing,Int}
end

function Measurement(; typ::MeasurementType, value::Real, sigma::Real, active::Bool = true, busIdx::Union{Nothing,Int} = nothing, branchIdx::Union{Nothing,Int} = nothing, direction::Symbol = :none, id::AbstractString = "", linkIdx::Union{Nothing,Int} = nothing)
  σ = Float64(sigma)
  σ <= 0.0 && error("Measurement sigma must be > 0.0")
  w = inv(σ * σ)
  return Measurement(typ, Float64(value), σ, w, active, busIdx, branchIdx, direction, String(id), linkIdx)
end

@inline function _default_measurement_id(typ::MeasurementType; busIdx::Union{Nothing,Int} = nothing, branchIdx::Union{Nothing,Int} = nothing, direction::Symbol = :none, linkIdx::Union{Nothing,Int} = nothing)
  if typ == VmMeas
    return "Vm_bus_$(something(busIdx, 0))"
  elseif typ == VaMeas
    return "Va_bus_$(something(busIdx, 0))"
  elseif typ == PinjMeas
    return "Pinj_bus_$(something(busIdx, 0))"
  elseif typ == QinjMeas
    return "Qinj_bus_$(something(busIdx, 0))"
  elseif typ == PflowMeas
    if isnothing(branchIdx) && !isnothing(linkIdx)
      return "Pflow_link_$(linkIdx)"
    end
    return "Pflow_branch_$(something(branchIdx, 0))_$(direction)"
  elseif typ == QflowMeas
    if isnothing(branchIdx) && !isnothing(linkIdx)
      return "Qflow_link_$(linkIdx)"
    end
    return "Qflow_branch_$(something(branchIdx, 0))_$(direction)"
  elseif typ == ImagMeas
    if isnothing(branchIdx) && !isnothing(busIdx)
      return "Imag_shunt_bus_$(busIdx)"
    end
    return "Imag_branch_$(something(branchIdx, 0))_$(direction)"
  elseif typ == ShuntQMeas
    return "ShuntQ_bus_$(something(busIdx, 0))"
  elseif typ == IaMeas
    if isnothing(branchIdx) && !isnothing(busIdx)
      return "Ia_shunt_bus_$(busIdx)"
    end
    return "Ia_branch_$(something(branchIdx, 0))_$(direction)"
  end
  error("Unsupported measurement type")
end

@inline function _resolve_branch_idx(net::Net; branchNr::Union{Nothing,Int} = nothing, fromBus::Union{Nothing,String} = nothing, toBus::Union{Nothing,String} = nothing)
  if !isnothing(branchNr)
    branchNr < 1 && error("branchNr must be > 0")
    branchNr > length(net.branchVec) && error("branchNr $(branchNr) not found in network")
    return branchNr
  end

  isnothing(fromBus) && error("Either branchNr or fromBus must be provided")
  isnothing(toBus) && error("Either branchNr or toBus must be provided")

  brVec = getNetBranchNumberVec(net = net, fromBus = fromBus, toBus = toBus)
  isempty(brVec) && error("No branch found between $(fromBus) and $(toBus)")
  length(brVec) > 1 && error("Multiple branches found between $(fromBus) and $(toBus); please specify branchNr")
  return brVec[1]
end

"""
    addMeasurement!(measurements; typ, value, sigma, active=true, busIdx=nothing, branchIdx=nothing, direction=:none, id="")

Append a state-estimation measurement to `measurements` and return it.
"""
function addMeasurement!(measurements::Vector; typ::MeasurementType, value::Real, sigma::Real, active::Bool = true, busIdx::Union{Nothing,Int} = nothing, branchIdx::Union{Nothing,Int} = nothing, direction::Symbol = :none, id::AbstractString = "", linkIdx::Union{Nothing,Int} = nothing)
  if !isnothing(linkIdx) && (typ != PflowMeas && typ != QflowMeas)
    error("linkIdx is only valid for flow measurements (PflowMeas/QflowMeas)")
  end
  if typ == VmMeas || typ == VaMeas || typ == PinjMeas || typ == QinjMeas
    isnothing(busIdx) && error("Bus measurement requires busIdx")
  elseif (typ == PflowMeas || typ == QflowMeas) && !isnothing(linkIdx)
    # link flow measurement (SE phase 3): allocation input, never a WLS row.
    # Positive from link.fromBus to link.toBus; only :from is accepted (one
    # measured quantity per link and kind is enough at zero impedance).
    isnothing(branchIdx) || error("Link flow measurement must not carry a branchIdx")
    direction == :to && error("Link flow measurement direction :to is not supported; measure positive from link.fromBus to link.toBus (direction :from)")
    direction = :from
  elseif typ == PflowMeas || typ == QflowMeas
    isnothing(branchIdx) && error("Flow measurement requires branchIdx")
    (direction == :from || direction == :to) || error("Flow measurement direction must be :from or :to")
  elseif typ == ImagMeas
    # two variants: branch current (branchIdx + direction) or shunt-bay
    # current (busIdx only, SE phase 2)
    if isnothing(branchIdx) && !isnothing(busIdx)
      direction == :none || error("Shunt-bay current measurement (busIdx) takes no direction")
    else
      isnothing(branchIdx) && error("Current-magnitude measurement requires branchIdx (branch variant) or busIdx (shunt-bay variant)")
      (direction == :from || direction == :to) || error("Current-magnitude measurement direction must be :from or :to")
    end
  elseif typ == ShuntQMeas
    isnothing(busIdx) && error("Shunt reactive-power measurement requires busIdx")
  end

  meas_id = isempty(id) ? _default_measurement_id(typ; busIdx = busIdx, branchIdx = branchIdx, direction = direction, linkIdx = linkIdx) : String(id)
  meas = Measurement(typ = typ, value = value, sigma = sigma, active = active, busIdx = busIdx, branchIdx = branchIdx, direction = direction, id = meas_id, linkIdx = linkIdx)
  push!(measurements, meas)
  return meas
end

function addMeasurement!(net::Net; kwargs...)
  return addMeasurement!(_net_measurements(net); kwargs...)
end

"""
    addVmMeasurement!(measurements; net, busName, value, sigma, active=true, id="")

Append a bus voltage-magnitude measurement identified by `busName`.
"""
function addVmMeasurement!(measurements::Vector; net::Net, busName::String, value::Real, sigma::Real, active::Bool = true, id::AbstractString = "")
  busIdx = geNetBusIdx(net = net, busName = busName)
  return addMeasurement!(measurements; typ = VmMeas, value = value, sigma = sigma, active = active, busIdx = busIdx, id = id)
end

function addVmMeasurement!(net::Net; busName::String, value::Real, sigma::Real, active::Bool = true, id::AbstractString = "")
  return addVmMeasurement!(_net_measurements(net); net = net, busName = busName, value = value, sigma = sigma, active = active, id = id)
end

"""
    addVaMeasurement!(measurements; net, busName, value, sigma, active=true, id="")

Append a bus voltage-angle measurement (PMU synchrophasor) identified by
`busName`. `value` and `sigma` are in degrees, referenced to the common PMU
time base (see the `pmu_ref_offset` handling in `runse!`).
"""
function addVaMeasurement!(measurements::Vector; net::Net, busName::String, value::Real, sigma::Real, active::Bool = true, id::AbstractString = "")
  busIdx = geNetBusIdx(net = net, busName = busName)
  return addMeasurement!(measurements; typ = VaMeas, value = value, sigma = sigma, active = active, busIdx = busIdx, id = id)
end

function addVaMeasurement!(net::Net; busName::String, value::Real, sigma::Real, active::Bool = true, id::AbstractString = "")
  return addVaMeasurement!(_net_measurements(net); net = net, busName = busName, value = value, sigma = sigma, active = active, id = id)
end

"""
    addPmuPhasorMeasurement!(measurements; net, busName, vm_pu, va_deg, sigmaVm=0.002, sigmaVa=0.02, active=true, idPrefix="PMU")

Append a complete PMU voltage-phasor measurement — magnitude and angle — for
one bus and return the pair `(vmMeas, vaMeas)`.

The magnitude enters as an ordinary `VmMeas` (p.u.); PMU accuracy is
expressed solely through the tight default `sigmaVm`. The angle enters as a
`VaMeas` (degrees, referenced to the common PMU time base; see the
`pmu_ref_offset` handling in `runse!`).
"""
function addPmuPhasorMeasurement!(measurements::Vector; net::Net, busName::String, vm_pu::Real, va_deg::Real, sigmaVm::Real = 0.002, sigmaVa::Real = 0.02, active::Bool = true, idPrefix::AbstractString = "PMU")
  busIdx = geNetBusIdx(net = net, busName = busName)
  vm = addMeasurement!(measurements; typ = VmMeas, value = vm_pu, sigma = sigmaVm, active = active, busIdx = busIdx, id = "$(idPrefix)_Vm_bus_$(busIdx)")
  va = addMeasurement!(measurements; typ = VaMeas, value = va_deg, sigma = sigmaVa, active = active, busIdx = busIdx, id = "$(idPrefix)_Va_bus_$(busIdx)")
  return (vm, va)
end

function addPmuPhasorMeasurement!(net::Net; busName::String, vm_pu::Real, va_deg::Real, sigmaVm::Real = 0.002, sigmaVa::Real = 0.02, active::Bool = true, idPrefix::AbstractString = "PMU")
  return addPmuPhasorMeasurement!(_net_measurements(net); net = net, busName = busName, vm_pu = vm_pu, va_deg = va_deg, sigmaVm = sigmaVm, sigmaVa = sigmaVa, active = active, idPrefix = idPrefix)
end

"""
    addPinjMeasurement!(measurements; net, busName, value, sigma, active=true, id="")

Append an active-power injection measurement identified by `busName`.
"""
function addPinjMeasurement!(measurements::Vector; net::Net, busName::String, value::Real, sigma::Real, active::Bool = true, id::AbstractString = "")
  busIdx = geNetBusIdx(net = net, busName = busName)
  return addMeasurement!(measurements; typ = PinjMeas, value = value, sigma = sigma, active = active, busIdx = busIdx, id = id)
end

function addPinjMeasurement!(net::Net; busName::String, value::Real, sigma::Real, active::Bool = true, id::AbstractString = "")
  return addPinjMeasurement!(_net_measurements(net); net = net, busName = busName, value = value, sigma = sigma, active = active, id = id)
end

"""
    addQinjMeasurement!(measurements; net, busName, value, sigma, active=true, id="")

Append a reactive-power injection measurement identified by `busName`.
"""
function addQinjMeasurement!(measurements::Vector; net::Net, busName::String, value::Real, sigma::Real, active::Bool = true, id::AbstractString = "")
  busIdx = geNetBusIdx(net = net, busName = busName)
  return addMeasurement!(measurements; typ = QinjMeas, value = value, sigma = sigma, active = active, busIdx = busIdx, id = id)
end

function addQinjMeasurement!(net::Net; busName::String, value::Real, sigma::Real, active::Bool = true, id::AbstractString = "")
  return addQinjMeasurement!(_net_measurements(net); net = net, busName = busName, value = value, sigma = sigma, active = active, id = id)
end

@inline function _resolve_link_idx(net::Net, linkNr::Int; branchNr = nothing, fromBus = nothing, toBus = nothing)
  (isnothing(branchNr) && isnothing(fromBus) && isnothing(toBus)) || error("linkNr (link flow measurement) is mutually exclusive with the branch keywords branchNr/fromBus/toBus")
  (1 <= linkNr <= length(net.linkVec)) || error("linkNr $(linkNr) not found in network")
  return linkNr
end

"""
    addPflowMeasurement!(measurements; net, value, sigma, direction=:from, branchNr=nothing, fromBus=nothing, toBus=nothing, linkNr=nothing, active=true, id="")

Append an active-power flow measurement identified by `branchNr` or a unique
`fromBus`/`toBus` branch pair, or, with `linkNr`, on a bus link (SE phase 3):
positive from `link.fromBus` to `link.toBus`, `direction = :from` only. A
link flow measurement never enters the WLS estimator; it is consumed by the
post-SE flow allocation (`calcLinkFlowsSE!`).
"""
function addPflowMeasurement!(measurements::Vector; net::Net, value::Real, sigma::Real, direction::Symbol = :from, branchNr::Union{Nothing,Int} = nothing, fromBus::Union{Nothing,String} = nothing, toBus::Union{Nothing,String} = nothing, linkNr::Union{Nothing,Int} = nothing, active::Bool = true, id::AbstractString = "")
  if !isnothing(linkNr)
    lidx = _resolve_link_idx(net, linkNr; branchNr = branchNr, fromBus = fromBus, toBus = toBus)
    return addMeasurement!(measurements; typ = PflowMeas, value = value, sigma = sigma, active = active, direction = direction, id = id, linkIdx = lidx)
  end
  bridx = _resolve_branch_idx(net; branchNr = branchNr, fromBus = fromBus, toBus = toBus)
  return addMeasurement!(measurements; typ = PflowMeas, value = value, sigma = sigma, active = active, branchIdx = bridx, direction = direction, id = id)
end

function addPflowMeasurement!(net::Net; value::Real, sigma::Real, direction::Symbol = :from, branchNr::Union{Nothing,Int} = nothing, fromBus::Union{Nothing,String} = nothing, toBus::Union{Nothing,String} = nothing, linkNr::Union{Nothing,Int} = nothing, active::Bool = true, id::AbstractString = "")
  return addPflowMeasurement!(_net_measurements(net); net = net, value = value, sigma = sigma, direction = direction, branchNr = branchNr, fromBus = fromBus, toBus = toBus, linkNr = linkNr, active = active, id = id)
end

"""
    addQflowMeasurement!(measurements; net, value, sigma, direction=:from, branchNr=nothing, fromBus=nothing, toBus=nothing, linkNr=nothing, active=true, id="")

Append a reactive-power flow measurement identified by `branchNr` or a unique
`fromBus`/`toBus` branch pair, or, with `linkNr`, on a bus link (SE phase 3,
see `addPflowMeasurement!` for the link conventions).
"""
function addQflowMeasurement!(measurements::Vector; net::Net, value::Real, sigma::Real, direction::Symbol = :from, branchNr::Union{Nothing,Int} = nothing, fromBus::Union{Nothing,String} = nothing, toBus::Union{Nothing,String} = nothing, linkNr::Union{Nothing,Int} = nothing, active::Bool = true, id::AbstractString = "")
  if !isnothing(linkNr)
    lidx = _resolve_link_idx(net, linkNr; branchNr = branchNr, fromBus = fromBus, toBus = toBus)
    return addMeasurement!(measurements; typ = QflowMeas, value = value, sigma = sigma, active = active, direction = direction, id = id, linkIdx = lidx)
  end
  bridx = _resolve_branch_idx(net; branchNr = branchNr, fromBus = fromBus, toBus = toBus)
  return addMeasurement!(measurements; typ = QflowMeas, value = value, sigma = sigma, active = active, branchIdx = bridx, direction = direction, id = id)
end

function addQflowMeasurement!(net::Net; value::Real, sigma::Real, direction::Symbol = :from, branchNr::Union{Nothing,Int} = nothing, fromBus::Union{Nothing,String} = nothing, toBus::Union{Nothing,String} = nothing, linkNr::Union{Nothing,Int} = nothing, active::Bool = true, id::AbstractString = "")
  return addQflowMeasurement!(_net_measurements(net); net = net, value = value, sigma = sigma, direction = direction, branchNr = branchNr, fromBus = fromBus, toBus = toBus, linkNr = linkNr, active = active, id = id)
end

"""
    addImagMeasurement!(measurements; net, value, sigma, direction=:from, branchNr=nothing, fromBus=nothing, toBus=nothing, busName=nothing, active=true, id="")

Append a current-magnitude measurement in ampere. Two variants:

- Branch current: identified by `branchNr` or a unique `fromBus`/`toBus`
  branch pair; `direction` selects the measured branch end (`:from`/`:to`).
- Shunt-bay current (SE phase 2): identified by `busName` alone (mutually
  exclusive with the branch keywords); measures the current drawn by the
  shunt at that bus. Errors when the bus carries no shunt.

Current magnitudes are auxiliary measurements: they must only
supplement power measurements, never replace them. The observability checks
exclude them, and the estimator activates them only from the configured
iteration on (`state_estimation.imag_activation_iteration`) and only when
`value >= 3 sigma` holds.
"""
function addImagMeasurement!(measurements::Vector; net::Net, value::Real, sigma::Real, direction::Symbol = :from, branchNr::Union{Nothing,Int} = nothing, fromBus::Union{Nothing,String} = nothing, toBus::Union{Nothing,String} = nothing, busName::Union{Nothing,String} = nothing, active::Bool = true, id::AbstractString = "")
  if !isnothing(busName)
    (isnothing(branchNr) && isnothing(fromBus) && isnothing(toBus)) || error("addImagMeasurement!: busName (shunt-bay variant) is mutually exclusive with the branch keywords branchNr/fromBus/toBus")
    busIdx = geNetBusIdx(net = net, busName = busName)
    haskey(net.shuntDict, busIdx) || error("addImagMeasurement!: no shunt at bus $(busName)")
    return addMeasurement!(measurements; typ = ImagMeas, value = value, sigma = sigma, active = active, busIdx = busIdx, direction = :none, id = id)
  end
  bridx = _resolve_branch_idx(net; branchNr = branchNr, fromBus = fromBus, toBus = toBus)
  return addMeasurement!(measurements; typ = ImagMeas, value = value, sigma = sigma, active = active, branchIdx = bridx, direction = direction, id = id)
end

"""
    addIaMeasurement!(measurements; net, value, sigma, direction=:from, branchNr=nothing, fromBus=nothing, toBus=nothing, busName=nothing, active=true, id="")

Append a PMU current-phasor ANGLE measurement in degrees, PMU-referenced
like `VaMeas` (the common reference-angle offset alpha applies). Same two
location variants as `addImagMeasurement!` (branch end or shunt bay). The
row is active in the solve only while the predicted current magnitude
passes the floor gate (see `IaMeas`); like `ImagMeas` it never carries
observability.
"""
function addIaMeasurement!(measurements::Vector; net::Net, value::Real, sigma::Real, direction::Symbol = :from, branchNr::Union{Nothing,Int} = nothing, fromBus::Union{Nothing,String} = nothing, toBus::Union{Nothing,String} = nothing, busName::Union{Nothing,String} = nothing, active::Bool = true, id::AbstractString = "")
  if !isnothing(busName)
    (isnothing(branchNr) && isnothing(fromBus) && isnothing(toBus)) || error("addIaMeasurement!: busName (shunt-bay variant) is mutually exclusive with the branch keywords branchNr/fromBus/toBus")
    busIdx = geNetBusIdx(net = net, busName = busName)
    haskey(net.shuntDict, busIdx) || error("addIaMeasurement!: no shunt at bus $(busName)")
    return addMeasurement!(measurements; typ = IaMeas, value = value, sigma = sigma, active = active, busIdx = busIdx, direction = :none, id = id)
  end
  bridx = _resolve_branch_idx(net; branchNr = branchNr, fromBus = fromBus, toBus = toBus)
  return addMeasurement!(measurements; typ = IaMeas, value = value, sigma = sigma, active = active, branchIdx = bridx, direction = direction, id = id)
end

function addIaMeasurement!(net::Net; value::Real, sigma::Real, direction::Symbol = :from, branchNr::Union{Nothing,Int} = nothing, fromBus::Union{Nothing,String} = nothing, toBus::Union{Nothing,String} = nothing, busName::Union{Nothing,String} = nothing, active::Bool = true, id::AbstractString = "")
  return addIaMeasurement!(_net_measurements(net); net = net, value = value, sigma = sigma, direction = direction, branchNr = branchNr, fromBus = fromBus, toBus = toBus, busName = busName, active = active, id = id)
end

"""
    addCurrentPhasorMeasurement!(measurements; net, i_A, ia_deg, sigmaI=10.0, sigmaIa=0.05, direction=:from, branchNr=nothing, fromBus=nothing, toBus=nothing, busName=nothing, active=true, idPrefix="PMU")

Append a full PMU current phasor as the `ImagMeas` + `IaMeas` pair (mirror
of `addPmuPhasorMeasurement!` for voltages). Returns `(imag, ia)`.
"""
function addCurrentPhasorMeasurement!(measurements::Vector; net::Net, i_A::Real, ia_deg::Real, sigmaI::Real = 10.0, sigmaIa::Real = 0.05, direction::Symbol = :from, branchNr::Union{Nothing,Int} = nothing, fromBus::Union{Nothing,String} = nothing, toBus::Union{Nothing,String} = nothing, busName::Union{Nothing,String} = nothing, active::Bool = true, idPrefix::AbstractString = "PMU")
  mi = addImagMeasurement!(measurements; net = net, value = i_A, sigma = sigmaI, direction = direction, branchNr = branchNr, fromBus = fromBus, toBus = toBus, busName = busName, active = active, id = "")
  ma = addIaMeasurement!(measurements; net = net, value = ia_deg, sigma = sigmaIa, direction = direction, branchNr = branchNr, fromBus = fromBus, toBus = toBus, busName = busName, active = active, id = "")
  mi2 = Measurement(typ = mi.typ, value = mi.value, sigma = mi.sigma, active = mi.active, busIdx = mi.busIdx, branchIdx = mi.branchIdx, direction = mi.direction, id = string(idPrefix, "_", mi.id))
  ma2 = Measurement(typ = ma.typ, value = ma.value, sigma = ma.sigma, active = ma.active, busIdx = ma.busIdx, branchIdx = ma.branchIdx, direction = ma.direction, id = string(idPrefix, "_", ma.id))
  measurements[end-1] = mi2
  measurements[end] = ma2
  return (mi2, ma2)
end

function addCurrentPhasorMeasurement!(net::Net; kwargs...)
  return addCurrentPhasorMeasurement!(_net_measurements(net); net = net, kwargs...)
end

function addImagMeasurement!(net::Net; value::Real, sigma::Real, direction::Symbol = :from, branchNr::Union{Nothing,Int} = nothing, fromBus::Union{Nothing,String} = nothing, toBus::Union{Nothing,String} = nothing, busName::Union{Nothing,String} = nothing, active::Bool = true, id::AbstractString = "")
  return addImagMeasurement!(_net_measurements(net); net = net, value = value, sigma = sigma, direction = direction, branchNr = branchNr, fromBus = fromBus, toBus = toBus, busName = busName, active = active, id = id)
end

"""
    addShuntQMeasurement!(measurements; net, busName, value, sigma, active=true, id="")

Append a shunt-bay reactive-power measurement in MVar for the shunt at
`busName` (sign convention identical to the solved `q_shunt` results).
Errors when the bus carries no shunt. This is the direct measurement that
releases a shunt susceptance for estimation (SE phase 2, case A; see
`setShuntEstimation!`).
"""
function addShuntQMeasurement!(measurements::Vector; net::Net, busName::String, value::Real, sigma::Real, active::Bool = true, id::AbstractString = "")
  busIdx = geNetBusIdx(net = net, busName = busName)
  haskey(net.shuntDict, busIdx) || error("addShuntQMeasurement!: no shunt at bus $(busName)")
  return addMeasurement!(measurements; typ = ShuntQMeas, value = value, sigma = sigma, active = active, busIdx = busIdx, id = id)
end

function addShuntQMeasurement!(net::Net; busName::String, value::Real, sigma::Real, active::Bool = true, id::AbstractString = "")
  return addShuntQMeasurement!(_net_measurements(net); net = net, busName = busName, value = value, sigma = sigma, active = active, id = id)
end

"""
    measurementStdDevs(; vm=0.005, pinj=1.0, qinj=1.0, pflow=1.0, qflow=1.0, va=0.02, imag=10.0, shuntq=1.0)

Create default standard-deviation map for synthetic measurement generation.
`va` is the PMU voltage-angle standard deviation in degrees; typical PMU
accuracy is 0.01 to 0.05° (IEEE C37.118 TVE < 1 %). `imag` is the
current-magnitude standard deviation in ampere; the 10 A default matches a
class 1 instrument on a measuring range of a few hundred ampere (typical
HV line current transformer). `shuntq` is the shunt-bay reactive-power
standard deviation in MVar.
"""
function measurementStdDevs(; vm::Float64 = 0.005, pinj::Float64 = 1.0, qinj::Float64 = 1.0, pflow::Float64 = 1.0, qflow::Float64 = 1.0, va::Float64 = 0.02, imag::Float64 = 10.0, shuntq::Float64 = 1.0, ia::Float64 = 0.05)
  return Dict(VmMeas => vm, PinjMeas => pinj, QinjMeas => qinj, PflowMeas => pflow, QflowMeas => qflow, VaMeas => va, ImagMeas => imag, ShuntQMeas => shuntq, IaMeas => ia)
end

"""
    measurementSigmaFloors(; kwargs...) -> Dict{MeasurementType,Float64}

Per-type sigma floors for `generateMeasurementsFromPF(relativeSigma = true)`:
a relative sigma is `max(fraction * |value|, floor)`, so a near-zero reading
never gets a near-zero sigma (weight `1/sigma^2` would explode). Units are
the measurement units: pu (`vm`), MW/MVar (`pinj`, `qinj`, `pflow`,
`qflow`, `shuntq`), ampere (`imag`), degrees (`va`, also the absolute sigma
in relative mode).

The power floors model the RANGE term of a transducer class (error =
percent of reading plus a range contribution) and double as a numerical
guard: a floor that is too small turns every zero-injection bus into a
near-constraint and visibly slows the flat-start Gauss-Newton (MiniGrid:
32 iterations at 0.01 MW versus 11 at the 0.05 MW default).
"""
function measurementSigmaFloors(; vm::Float64 = 1e-4, pinj::Float64 = 0.05, qinj::Float64 = 0.05, pflow::Float64 = 0.05, qflow::Float64 = 0.05, va::Float64 = 0.01, imag::Float64 = 0.1, shuntq::Float64 = 0.05, ia::Float64 = 0.01)
  return Dict(VmMeas => vm, PinjMeas => pinj, QinjMeas => qinj, PflowMeas => pflow, QflowMeas => qflow, VaMeas => va, ImagMeas => imag, ShuntQMeas => shuntq, IaMeas => ia)
end

@inline _bus_power_value(x::Union{Nothing,Float64}) = isnothing(x) ? 0.0 : x

"""
    findPassiveBuses(net; atol=1e-9, includeSlack=false) -> Vector{Int}

Return bus indices that have no generation, no load, and no shunt contribution
within the given tolerance `atol`.

This is useful for state-estimation workflows where passive / transit buses are
often modeled through zero-injection pseudo-measurements.
"""
function findPassiveBuses(net::Net; atol::Float64 = 1e-9, includeSlack::Bool = false)
  passive = Int[]
  for i in eachindex(net.nodeVec)
    node = net.nodeVec[i]
    if !includeSlack && getNodeType(node) == Slack
      continue
    end

    pinj = _bus_power_value(node._pƩGen) - _bus_power_value(node._pƩLoad) - _bus_power_value(node._pShunt)
    qinj = _bus_power_value(node._qƩGen) - _bus_power_value(node._qƩLoad) - _bus_power_value(node._qShunt)
    if abs(pinj) <= atol && abs(qinj) <= atol
      push!(passive, i)
    end
  end
  return passive
end

"""
    ZERO_INJECTION_SIGMA

Default sigma of a zero-injection pseudo-measurement.

A passive bus injects nothing, so the row is a CONSTRAINT rather than a
reading, and the sigma only says how hard it is enforced. The former 1e-6
means a weight of 1e12, and on a large network that is not "tight", it is
crushing: on a 25000-bus case with 27268 such rows the estimation did not
converge in 30 iterations, reported an objective of 3.8e15 per degree of
freedom and losses 2069 percent above the reference. The same set with 1e-3
converges in 7 iterations to J/dof 4.97 and losses within 0.02 percent.

1e-3 still binds the bus far tighter than any real transducer (the tightest
real rows in that set carry 0.01), while leaving the normal equations
solvable.
"""
const ZERO_INJECTION_SIGMA = 1.0e-3

"""
    addZeroInjectionMeasurements!(measurements; net, sigma=ZERO_INJECTION_SIGMA, busNames=nothing, busIdxs=nothing, active=true, idPrefix="ZI") -> Vector{Measurement}

Append active- and reactive-power zero-injection pseudo-measurements for the
selected buses and return the newly added measurements.

Selection rules:
- If `busIdxs` is provided, those indices are used.
- Else if `busNames` is provided, names are resolved to indices.
- Else passive buses are detected automatically via `findPassiveBuses(net)`.

These pseudo-measurements are the current way to encode equality constraints
`P_inj = 0` and `Q_inj = 0` in the WLS estimator.
"""
function addZeroInjectionMeasurements!(measurements::Vector; net::Net, sigma::Real = ZERO_INJECTION_SIGMA, busNames::Union{Nothing,Vector{String}} = nothing, busIdxs::Union{Nothing,Vector{Int}} = nothing, active::Bool = true, idPrefix::AbstractString = "ZI")
  selected = if !isnothing(busIdxs)
    copy(busIdxs)
  elseif !isnothing(busNames)
    [geNetBusIdx(net = net, busName = busName) for busName in busNames]
  else
    findPassiveBuses(net)
  end

  isempty(selected) && return Measurement[]

  added = Measurement[]
  for busIdx in selected
    1 <= busIdx <= length(net.nodeVec) || error("Bus index $(busIdx) out of bounds")
    push!(added, addMeasurement!(measurements; typ = PinjMeas, value = 0.0, sigma = sigma, active = active, busIdx = busIdx, id = "$(idPrefix)_PINJ_bus_$(busIdx)"))
    push!(added, addMeasurement!(measurements; typ = QinjMeas, value = 0.0, sigma = sigma, active = active, busIdx = busIdx, id = "$(idPrefix)_QINJ_bus_$(busIdx)"))
  end
  return added
end

function addZeroInjectionMeasurements!(net::Net; sigma::Real = ZERO_INJECTION_SIGMA, busNames::Union{Nothing,Vector{String}} = nothing, busIdxs::Union{Nothing,Vector{Int}} = nothing, active::Bool = true, idPrefix::AbstractString = "ZI")
  return addZeroInjectionMeasurements!(_net_measurements(net); net = net, sigma = sigma, busNames = busNames, busIdxs = busIdxs, active = active, idPrefix = idPrefix)
end

# Complex shunt-bay power in MVA at busIdx: S = |V|^2 * conj(G + jB) * base,
# the exact Ybus stamping term (same expression as updateShuntPowers!). The
# susceptance defaults to the model value; `shuntB` (busIdx -> B in pu)
# overrides it for shunts whose B is an estimator state (SE phase 2).
@inline function _shunt_bay_S_MVA(net::Net, busIdx::Int, V::Vector{ComplexF64}, shuntB::Union{Nothing,Dict{Int,Float64}})
  haskey(net.shuntDict, busIdx) || error("Shunt measurement at bus $(busIdx): no shunt registered at this bus")
  sh = net.shuntVec[net.shuntDict[busIdx]]
  b = shuntB !== nothing && haskey(shuntB, busIdx) ? shuntB[busIdx] : imag(sh.y_pu_shunt)
  y = complex(real(sh.y_pu_shunt), b)
  return abs2(V[busIdx]) * conj(y) * net.baseMVA
end

## complex three-phase end current in ampere at the referenced branch end
## (:from/:to) or shunt bay (busIdx variant), shared by the ImagMeas
## magnitude and the IaMeas angle: I = 1000/(sqrt(3) Vn_kV) * conj(S_MVA / V_pu)
## (S = V conj(I) => I = conj(S/V); |I| matches the established magnitude
## formula bitwise since |conj(S/V)| = |S|/|V|). Returns 0 + 0im on a dead
## bus voltage.
function _end_current_A(net::Net, meas::Measurement, V::Vector{ComplexF64}, shuntB::Union{Nothing,Dict{Int,Float64}}, tapOverlay = nothing)::ComplexF64
  bridx = something(meas.branchIdx, 0)
  if bridx < 1
    # shunt-bay variant (SE phase 2): the current drawn by the shunt at the
    # referenced bus, from the same analytic bay power as ShuntQMeas
    i = something(meas.busIdx, 0)
    i < 1 && error("current measurement requires branchIdx or busIdx")
    abs(V[i]) <= eps(Float64) && return 0.0 + 0.0im
    s_MVA = _shunt_bay_S_MVA(net, i, V, shuntB)
    return 1000.0 / (sqrt(3.0) * getNodeVn(net.nodeVec[i])) * conj(s_MVA / V[i])
  end
  return _end_current_A(net, _prediction_branch(net, bridx, tapOverlay), meas.direction, V)
end

## branch-end variant of the shared current helper
function _end_current_A(net::Net, br, direction::Symbol, V::Vector{ComplexF64})::ComplexF64
  if direction == :from
    s_MVA = branchFlow_pu(br, br.fromBus, br.toBus, 1, V) * net.baseMVA
    endBus = Int(br.fromBus)
  elseif direction == :to
    s_MVA = branchFlow_pu(br, br.toBus, br.fromBus, 2, V) * net.baseMVA
    endBus = Int(br.toBus)
  else
    error("current measurement direction must be :from or :to")
  end
  abs(V[endBus]) <= eps(Float64) && return 0.0 + 0.0im
  return 1000.0 / (sqrt(3.0) * getNodeVn(net.nodeVec[endBus])) * conj(s_MVA / V[endBus])
end

## resolve the branch a row evaluates against: the tap overlay's scratch
## branch (released trafo at the cascade position) or the model branch
@inline function _prediction_branch(net::Net, bridx::Int, tapOverlay)
  tapOverlay !== nothing && haskey(tapOverlay, bridx) && return tapOverlay[bridx].branch
  return net.branchVec[bridx]
end

@inline function _measurement_prediction(meas::Measurement, net::Net, V::Vector{ComplexF64}, Sbus_MVA::Vector{ComplexF64}, vaOffsetRad::Float64 = 0.0, shuntB::Union{Nothing,Dict{Int,Float64}} = nothing, tapOverlay = nothing)
  if meas.typ == VmMeas
    i = something(meas.busIdx, 0)
    i < 1 && error("Vm measurement missing busIdx")
    return abs(V[i])
  elseif meas.typ == VaMeas
    i = something(meas.busIdx, 0)
    i < 1 && error("Va measurement missing busIdx")
    # PMU angle model: z = θ_i(slack-referenced) + α, α = slack angle in the
    # PMU time base. α is 0 unless the estimator carries the offset state.
    return rad2deg(angle(V[i]) + vaOffsetRad)
  elseif meas.typ == PinjMeas
    i = something(meas.busIdx, 0)
    i < 1 && error("Pinj measurement missing busIdx")
    return real(Sbus_MVA[i])
  elseif meas.typ == QinjMeas
    i = something(meas.busIdx, 0)
    i < 1 && error("Qinj measurement missing busIdx")
    return imag(Sbus_MVA[i])
  elseif meas.typ == PflowMeas || meas.typ == QflowMeas
    bridx = something(meas.branchIdx, 0)
    bridx < 1 && error("Flow measurement missing branchIdx")
    br = _prediction_branch(net, bridx, tapOverlay)
    if meas.direction == :from
      s = branchFlow_pu(br, br.fromBus, br.toBus, 1, V) * net.baseMVA
    elseif meas.direction == :to
      s = branchFlow_pu(br, br.toBus, br.fromBus, 2, V) * net.baseMVA
    else
      error("Flow measurement direction must be :from or :to")
    end
    return meas.typ == PflowMeas ? real(s) : imag(s)
  elseif meas.typ == ShuntQMeas
    i = something(meas.busIdx, 0)
    i < 1 && error("ShuntQ measurement missing busIdx")
    return imag(_shunt_bay_S_MVA(net, i, V, shuntB))
  elseif meas.typ == ImagMeas
    return abs(_end_current_A(net, meas, V, shuntB, tapOverlay))
  elseif meas.typ == IaMeas
    # PMU current-phasor angle: network-referenced phase of the same complex
    # end current the magnitude uses, plus the PMU reference offset alpha
    # (identical convention to VaMeas). Near zero current the angle is
    # undefined; the activity gate keeps such rows out of the solve.
    return rad2deg(angle(_end_current_A(net, meas, V, shuntB, tapOverlay)) + vaOffsetRad)
  else
    error("Unsupported measurement type")
  end
end

"""
    generateMeasurementsFromPF(net; kwargs...) -> Vector{Measurement}

Generate synthetic measurements from the current solved network state.

Keyword options:
- `includeVm`, `includePinj`, `includeQinj`, `includePflow`, `includeQflow`
- `includeImag`: add branch current-magnitude measurements in ampere for both
  ends (default `false`; auxiliary only, see `addImagMeasurement!`)
- `includeShuntQ`: add shunt-bay reactive-power measurements in MVar for every
  in-service shunt (default `false`; sign convention of the `q_shunt` results)
- `includeVa`: add PMU voltage-angle measurements (degrees, default `false`)
- `vaBusIdxs`: restrict Va measurements to these bus indices (default: all buses)
- `vaRefOffsetDeg`: common angle offset added to all generated Va values,
  emulating a PMU time base that differs from the slack reference
- `noise`: add Gaussian noise if `true`
- `stddev`: dictionary from `MeasurementType => sigma`
- `relativeSigma`: interpret the `stddev` values as FRACTIONS of the true
  measured value (`0.01` = 1 percent of reading) instead of absolute units.
  This keeps one accuracy setting meaningful across voltage levels: 1
  percent of a 400 MW flow and of a 4 MW flow both get a class-appropriate
  sigma, where an absolute 2 MW would be tight at 400 kV and absurd at
  30 kV. Per row `sigma = max(fraction * |value|, sigmaFloor[type])`; the
  floor guards near-zero readings (a zero flow must not get weight
  infinity). `VaMeas` stays absolute in degrees even in relative mode (an
  angle passes through zero, a fraction of it is meaningless).
- `sigmaFloor`: the per-type floors for `relativeSigma` (default
  `measurementSigmaFloors()`)
- `rng`: random number generator
"""
function generateMeasurementsFromPF(
  net::Net;
  includeVm::Bool = true,
  includePinj::Bool = true,
  includeQinj::Bool = true,
  includePflow::Bool = true,
  includeQflow::Bool = true,
  includeImag::Bool = false,
  includeIa::Bool = false,
  includeShuntQ::Bool = false,
  includeVa::Bool = false,
  vaBusIdxs::Union{Nothing,Vector{Int}} = nothing,
  vaRefOffsetDeg::Float64 = 0.0,
  noise::Bool = false,
  stddev::Dict{MeasurementType,Float64} = measurementStdDevs(),
  relativeSigma::Bool = false,
  sigmaFloor::Dict{MeasurementType,Float64} = measurementSigmaFloors(),
  rng::AbstractRNG = Random.default_rng(),
)
  # per-row sigma: absolute from stddev, or (relative mode) a fraction of
  # the reading with the per-type floor; the angle types Va and Ia stay
  # absolute in both modes (an angle passes through zero)
  rowσ(typ::MeasurementType, value::Float64) = (relativeSigma && !(typ in (VaMeas, IaMeas))) ? max(stddev[typ] * abs(value), sigmaFloor[typ]) : stddev[typ]
  includeIa && !includeImag && @warn "generateMeasurementsFromPF: current ANGLES without current magnitudes (includeIa without includeImag); the IaMeas activity gate then falls back to the configured current floor"
  V = buildVoltageVector(net)
  # createYBUS drops isolated nodes (its rows are iso-shifted), so the
  # injection computation must run in that compressed space; V stays in full
  # bus numbering for the voltage and branch-flow rows. Isolated buses get
  # no measurements at all: they carry no defined state.
  iso_sorted = sort(net.isoNodes)
  is_iso(b::Int) = insorted(b, iso_sorted)
  iso_shift(b::Int) = b - (searchsortedfirst(iso_sorted, b) - 1)
  Vc = isempty(iso_sorted) ? V : V[[i for i in eachindex(V) if !is_iso(i)]]
  Sbus_pu = calc_injections(createYBUS(net = net), Vc)
  Sbus_MVA = Sbus_pu .* net.baseMVA

  m = Measurement[]
  nbus = length(net.nodeVec)

  vaSelected = isnothing(vaBusIdxs) ? nothing : Set(vaBusIdxs)

  for i = 1:nbus
    is_iso(i) && continue
    if includeVm
      σ = rowσ(VmMeas, abs(V[i]))
      z = abs(V[i]) + (noise ? randn(rng) * σ : 0.0)
      push!(m, Measurement(typ = VmMeas, value = z, sigma = σ, busIdx = i, id = "Vm_bus_$(i)"))
    end
    if includeVa && (isnothing(vaSelected) || i in vaSelected)
      σ = stddev[VaMeas]
      z = rad2deg(angle(V[i])) + vaRefOffsetDeg + (noise ? randn(rng) * σ : 0.0)
      push!(m, Measurement(typ = VaMeas, value = z, sigma = σ, busIdx = i, id = "Va_bus_$(i)"))
    end
    if includePinj
      σ = rowσ(PinjMeas, real(Sbus_MVA[iso_shift(i)]))
      z = real(Sbus_MVA[iso_shift(i)]) + (noise ? randn(rng) * σ : 0.0)
      push!(m, Measurement(typ = PinjMeas, value = z, sigma = σ, busIdx = i, id = "Pinj_bus_$(i)"))
    end
    if includeQinj
      σ = rowσ(QinjMeas, imag(Sbus_MVA[iso_shift(i)]))
      z = imag(Sbus_MVA[iso_shift(i)]) + (noise ? randn(rng) * σ : 0.0)
      push!(m, Measurement(typ = QinjMeas, value = z, sigma = σ, busIdx = i, id = "Qinj_bus_$(i)"))
    end
  end

  for br in net.branchVec
    # no flow measurements on branches that are fully open or open at one
    # terminal (r0.9.10): the estimator's measurement model rejects the
    # partial case, see the runse! scope guard
    _branch_terminal_state(br) == :closed || continue
    # a nominally closed branch hanging on an isolated bus carries no flow
    (is_iso(Int(br.fromBus)) || is_iso(Int(br.toBus))) && continue
    if includePflow
      sfrom = branchFlow_pu(br, br.fromBus, br.toBus, 1, V) * net.baseMVA
      sto = branchFlow_pu(br, br.toBus, br.fromBus, 2, V) * net.baseMVA
      σf = rowσ(PflowMeas, real(sfrom))
      σt = rowσ(PflowMeas, real(sto))
      push!(m, Measurement(typ = PflowMeas, value = real(sfrom) + (noise ? randn(rng) * σf : 0.0), sigma = σf, branchIdx = br.branchIdx, direction = :from, id = "Pflow_branch_$(br.branchIdx)_from"))
      push!(m, Measurement(typ = PflowMeas, value = real(sto) + (noise ? randn(rng) * σt : 0.0), sigma = σt, branchIdx = br.branchIdx, direction = :to, id = "Pflow_branch_$(br.branchIdx)_to"))
    end
    if includeQflow
      sfrom = branchFlow_pu(br, br.fromBus, br.toBus, 1, V) * net.baseMVA
      sto = branchFlow_pu(br, br.toBus, br.fromBus, 2, V) * net.baseMVA
      σf = rowσ(QflowMeas, imag(sfrom))
      σt = rowσ(QflowMeas, imag(sto))
      push!(m, Measurement(typ = QflowMeas, value = imag(sfrom) + (noise ? randn(rng) * σf : 0.0), sigma = σf, branchIdx = br.branchIdx, direction = :from, id = "Qflow_branch_$(br.branchIdx)_from"))
      push!(m, Measurement(typ = QflowMeas, value = imag(sto) + (noise ? randn(rng) * σt : 0.0), sigma = σt, branchIdx = br.branchIdx, direction = :to, id = "Qflow_branch_$(br.branchIdx)_to"))
    end
    if includeImag || includeIa
      # shared complex end current (same helper as the predictions), per end
      for dir in (:from, :to)
        Ic = _end_current_A(net, br, dir, V)
        if includeImag
          i_A = abs(Ic)
          σ = rowσ(ImagMeas, i_A)
          push!(m, Measurement(typ = ImagMeas, value = i_A + (noise ? randn(rng) * σ : 0.0), sigma = σ, branchIdx = br.branchIdx, direction = dir, id = "Imag_branch_$(br.branchIdx)_$(dir)"))
        end
        if includeIa
          σa = stddev[IaMeas]   # absolute degrees, like VaMeas
          ang = rad2deg(angle(Ic)) + vaRefOffsetDeg
          push!(m, Measurement(typ = IaMeas, value = ang + (noise ? randn(rng) * σa : 0.0), sigma = σa, branchIdx = br.branchIdx, direction = dir, id = "Ia_branch_$(br.branchIdx)_$(dir)"))
        end
      end
    end
  end

  if includeShuntQ
    for sh in net.shuntVec
      sh.status == 1 || continue
      is_iso(sh.busIdx) && continue
      i = sh.busIdx
      # same analytic bay power as the ShuntQMeas prediction (and as
      # updateShuntPowers!), so the roundtrip is exact
      q_MVar = imag(abs2(V[i]) * conj(sh.y_pu_shunt) * net.baseMVA)
      σ = rowσ(ShuntQMeas, q_MVar)
      push!(m, Measurement(typ = ShuntQMeas, value = q_MVar + (noise ? randn(rng) * σ : 0.0), sigma = σ, busIdx = i, id = "ShuntQ_bus_$(i)"))
    end
  end

  return m
end

"""
    deriveShuntPseudoMeasurements!(net; sigmaFloor=nothing) -> Vector{Measurement}

Shunt back-calculation (SE phase 2, case B): a measurement-preprocessing step
that runs BEFORE `runse!` and derives a `ShuntQMeas` pseudo-measurement (id
prefix `SHDERIV`) for every in-service shunt whose bay has

- an active shunt-bay current measurement (bus-referenced `ImagMeas`),
- NO active direct `ShuntQMeas` (a direct measurement takes precedence,
  case A), and
- an active `VmMeas` at the bus. The voltage enters `Q = -B V^2` squared, so
  a measured voltage is a hard prerequisite; without one the
  shunt is skipped with a warning, never computed from nominal voltage.

Derivation per shunt: `U_kV = vn_kV * vm_meas`, `S_MVA = sqrt(3) U_kV I_A /
1000`, magnitude `sqrt(max(S^2 - P^2, 0))` with the sign taken from the model
susceptance (reactor versus capacitor, `q_shunt` convention).

Deliberate 0.10.0 restriction: no bay active-power measurement type exists,
so the bay `P` is always taken as 0 and the propagated sigma is doubled to
cover the neglected active part. For reactors and capacitor banks `P` is
small against `Q` and the error is secondary; for filter circuits with a
substantial active component the derived value is biased, use a direct
`ShuntQMeas` (case A) there instead.

The sigma comes from first-order error propagation over `sigma_I` and
`sigma_Vm`; `sigmaFloor` optionally sets a minimum sigma so the derived
pseudo-measurement cannot outweigh real measurements. Returns the added
measurements.
"""
function deriveShuntPseudoMeasurements!(net::Net; sigmaFloor::Union{Nothing,Float64} = nothing)
  meas = _net_measurements(net)
  added = Measurement[]
  for sh in net.shuntVec
    sh.status == 1 || continue
    i = sh.busIdx
    # direct Q measurement present -> case A territory, nothing to derive
    any(m -> m.active && m.typ == ShuntQMeas && m.busIdx == i, meas) && continue
    imeas = findfirst(m -> m.active && m.typ == ImagMeas && m.branchIdx === nothing && m.busIdx == i, meas)
    imeas === nothing && continue
    vmeas = findfirst(m -> m.active && m.typ == VmMeas && m.busIdx == i, meas)
    if vmeas === nothing
      @warn "deriveShuntPseudoMeasurements!: shunt at bus $(i) skipped (no active Vm measurement at the bus; the voltage enters Q quadratically and nominal voltage is not an acceptable substitute)"
      continue
    end
    b_model = imag(sh.y_pu_shunt)
    if b_model == 0.0
      @warn "deriveShuntPseudoMeasurements!: shunt at bus $(i) skipped (model susceptance is zero, no sign for the derived Q)"
      continue
    end
    mi = meas[imeas]
    mv = meas[vmeas]
    vn_kV = getNodeVn(net.nodeVec[i])
    U_kV = vn_kV * mv.value
    S_MVA = sqrt(3.0) * U_kV * mi.value / 1000.0
    # no bay P measurement type exists yet: P = 0, sigma widened by 2 below
    Q_mag = S_MVA
    # sign from the model susceptance, q_shunt convention Q = -|V|^2 B base
    q_derived = -sign(b_model) * Q_mag
    # first-order propagation: dS/dI = sqrt(3) U / 1000, dS/dVm = sqrt(3) vn I / 1000
    σ = sqrt((sqrt(3.0) * U_kV / 1000.0 * mi.sigma)^2 + (sqrt(3.0) * vn_kV * mi.value / 1000.0 * mv.sigma)^2)
    σ = 2.0 * σ   # widened: bay P unmeasured (losses neglected in Q = sqrt(S^2 - P^2))
    if sigmaFloor !== nothing
      σ = max(σ, sigmaFloor)
    end
    push!(added, addMeasurement!(meas; typ = ShuntQMeas, value = q_derived, sigma = σ, busIdx = i, id = "SHDERIV_ShuntQ_bus_$(i)"))
  end
  return added
end

"""
    addMeasurementNoise!(net::Net; stddev, relativeSigma = true, sigmaFloor, rng, sigmas_from_rows = true) -> Vector{Measurement}

Perturb the measurement values `net` already holds with Gaussian noise, in
place, and return the perturbed rows.

This is the counterpart to generating a set: it needs no power flow and no
truth state, because it works on the VALUES that are there. That matters for
a set delivered with a case file, which is typically noise-free: an
estimation on ideal values returns `J = 0` by construction, and the only way
to a realistic run used to be regenerating the whole set from a fresh solve.

Each row is drawn as `value + sigma * randn(rng)`. `sigmas_from_rows = true`
(the default) uses each row's OWN sigma, which is what the set already
declares about its accuracy; with `false` the sigma comes from `stddev`
(relative to the reading when `relativeSigma`, with the per-type floor). A
row with a virtual sigma (a zero-injection constraint, `sigma <= 1e-6`) is
never touched: it states a physical law, not a reading.
"""
function addMeasurementNoise!(
  net::Net;
  stddev::Dict{MeasurementType,Float64} = measurementStdDevs(),
  relativeSigma::Bool = true,
  sigmaFloor::Dict{MeasurementType,Float64} = measurementSigmaFloors(),
  rng::AbstractRNG = Random.default_rng(),
  sigmas_from_rows::Bool = true,
)::Vector{Measurement}
  isempty(net.measurements) && throw(ArgumentError("addMeasurementNoise!: the network carries no measurements to perturb."))
  for (k, m) in enumerate(net.measurements)
    # virtual rows encode a constraint (zero injection), not a reading
    m.sigma <= 1e-6 && continue
    σ = if sigmas_from_rows
      m.sigma
    elseif relativeSigma && !(m.typ in (VaMeas, IaMeas))
      max(stddev[m.typ] * abs(m.value), sigmaFloor[m.typ])
    else
      stddev[m.typ]
    end
    net.measurements[k] = Measurement(typ = m.typ, value = m.value + σ * randn(rng), sigma = σ, active = m.active,
                                      busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id, linkIdx = m.linkIdx)
  end
  return Measurement[m for m in net.measurements]
end

"""
    setMeasurementsFromPF!(net; kwargs...)

Replace the network's measurements with a set generated from its solved
power flow (see `generateMeasurementsFromPF`).
"""
function setMeasurementsFromPF!(net::Net; kwargs...)
  measurements = generateMeasurementsFromPF(net; kwargs...)
  empty!(net.measurements)
  append!(net.measurements, measurements)
  return Measurement[m for m in net.measurements]
end

# ---------------------------------------------------------------------------
# Measurement CSV v1 (SE phase 5)
# ---------------------------------------------------------------------------

const _MEASUREMENT_CSV_VERSION = "sparlectra-measurements v1"
const _MEASUREMENT_CSV_HEADER = "type,bus,from_bus,to_bus,branch_nr,link_nr,direction,value,sigma,active,id"
## first-cut v1 header without the branch_nr column: still accepted on read
## (branch resolution falls back to the bus pair, rejecting parallels)
const _MEASUREMENT_CSV_HEADER_NO_BRANCH_NR = "type,bus,from_bus,to_bus,link_nr,direction,value,sigma,active,id"

@inline _csv_field(s) = strip(s)

"""
    writeMeasurementsCSV(net; file) -> NamedTuple

Write `net.measurements` as a measurement CSV v1 file:

    # sparlectra-measurements v1
    type,bus,from_bus,to_bus,link_nr,direction,value,sigma,active,id

One row per measurement; exactly one location group per row (`bus` for bus
measurements including the shunt-bay `ImagMeas`, `from_bus`/`to_bus` plus
`branch_nr` for branch measurements, `link_nr` for link flow measurements).
Buses are written by NAME (matched whitespace-insensitively on read) or,
with `busReference = :mrid`, by their preserved CGMES mRID (ENTSO-E UUID;
name fallback where none is recorded; MATPOWER/DTF nets have no mRIDs and
always use bus numbers/names). The reader resolves names, component
ids/names, and mRIDs. Values use round-trip-exact Float64 formatting
(decimal point, UTF-8), so `readMeasurementsCSV!` restores a bitwise-equal
measurement vector. Returns `(count = number of rows written,)`.

`branch_nr` carries the branch index and disambiguates parallel branches
between the same bus pair (routine in CGMES nets); the reader verifies it
against the named endpoints. Files written without the column (first-cut
v1 header) are still read, with pair-based resolution that rejects
parallels (same rule as `addPflowMeasurement!`). `headerComments` lines are
written as `#` comments after the version line (free-form metadata, e.g.
the transformer tap positions a generated set is based on); the reader
skips them.
"""
function writeMeasurementsCSV(net::Net; file::AbstractString, headerComments::Vector{String} = String[], busReference::Symbol = :name)
  busReference in (:name, :mrid) || error("writeMeasurementsCSV: busReference must be :name or :mrid")
  name_by_idx = _bus_name_by_idx(net)
  # :mrid (CGMES only): reference buses by their preserved ENTSO-E UUID;
  # falls back to the name where no mRID is recorded. MATPOWER/DTF nets
  # have no mRIDs and always write bus numbers/names.
  function busname(i)
    nm = get(name_by_idx, i, string(i))
    busReference == :name && return nm
    m = _bus_mrid(net, i; name_by_idx = name_by_idx)
    return isempty(m) ? nm : m
  end
  n = 0
  open(file, "w") do io
    println(io, "# ", _MEASUREMENT_CSV_VERSION)
    # free-form metadata (e.g. the transformer tap positions the set was
    # generated from); the reader skips comment lines
    for c in headerComments
      println(io, "# ", c)
    end
    println(io, _MEASUREMENT_CSV_HEADER)
    for m in net.measurements
      bus = ""
      fromB = ""
      toB = ""
      branchNr = ""
      linkNr = ""
      if m.linkIdx !== nothing
        linkNr = string(m.linkIdx)
      elseif m.branchIdx !== nothing
        br = net.branchVec[m.branchIdx]
        fromB = busname(Int(br.fromBus))
        toB = busname(Int(br.toBus))
        branchNr = string(m.branchIdx)
      elseif m.busIdx !== nothing
        bus = busname(m.busIdx)
      end
      dir = m.direction == :none ? "" : String(m.direction)
      println(io, string(m.typ), ",", bus, ",", fromB, ",", toB, ",", branchNr, ",", linkNr, ",", dir, ",", repr(m.value), ",", repr(m.sigma), ",", m.active ? "true" : "false", ",", m.id)
      n += 1
    end
  end
  return (count = n,)
end

## CGMES mRIDs (ENTSO-E rdf:ID UUIDs, canonicalized without the leading
## underscore) survive the import in net.cgmes_ids under structural keys;
## empty string when the net is not CGMES-sourced or the object was not
## recorded. MATPOWER/DTF nets have no mRIDs; bus numbers stay the
## reference there.
function _bus_mrid(net::Net, busIdx::Int; name_by_idx::AbstractDict = _bus_name_by_idx(net))::String
  nm = get(name_by_idx, busIdx, "")
  isempty(nm) && return ""
  return get(net.cgmes_ids, CGMESImporter.cgmesKeyTopologicalNode(nm), "")
end

## branchIdx -> PowerTransformer mRID for every transformer branch, using
## the exporter's first-seen parallel counting so the structural keys match
## the ones the importer recorded
function _transformer_mrids(net::Net)::Dict{Int,String}
  out = Dict{Int,String}()
  isempty(net.cgmes_ids) && return out
  nby = _bus_name_by_idx(net)
  counter = Dict{Tuple{String,String},Int}()
  for (k, br) in enumerate(net.branchVec)
    br.ratio != 0.0 || continue
    a = get(nby, Int(br.fromBus), string(Int(br.fromBus)))
    b = get(nby, Int(br.toBus), string(Int(br.toBus)))
    kk = CGMESImporter.cgmesNextParallelIndex!(counter, a, b)
    m = get(net.cgmes_ids, CGMESImporter.cgmesKeyPowerTransformer(a, b, kk), "")
    isempty(m) || (out[k] = m)
  end
  return out
end

## CSV bus-name resolution: the reader strips every field, but bus names
## (CGMES) can legally carry leading/trailing whitespace in busDict. Exact
## key first, then a whitespace-normalized unique match; ambiguity errors
## instead of guessing.
function _measurement_csv_bus_key(net::Net, name::AbstractString)::String
  haskey(net.busDict, name) && return String(name)
  hits = [String(k) for k in keys(net.busDict) if strip(k) == name]
  length(hits) == 1 && return hits[1]
  length(hits) > 1 && error("bus name '$(name)' is ambiguous after whitespace normalization ($(length(hits)) matches)")
  # component-id fallback: CGMES-style files may reference buses by the
  # component name or id instead of the busDict key (the two diverge on
  # CGMES imports); a unique match resolves, ambiguity errors
  name_by_idx = _bus_name_by_idx(net)
  cid_hits = Int[]
  for (i, nd) in enumerate(net.nodeVec)
    (getCompName(nd.comp) == name || nd.comp.cID == name) && push!(cid_hits, i)
  end
  length(cid_hits) == 1 && haskey(name_by_idx, cid_hits[1]) && return name_by_idx[cid_hits[1]]
  length(cid_hits) > 1 && error("bus reference '$(name)' matches $(length(cid_hits)) component ids/names")
  # CGMES mRID fallback: TN|<busname> => mRID lives in net.cgmes_ids; a file
  # written with busReference = :mrid references buses by exactly that UUID
  for (k, v) in net.cgmes_ids
    if v == name && startswith(k, "TN|")
      busname = String(k[4:end])
      haskey(net.busDict, busname) && return busname
    end
  end
  error("Bus $(name) not found in the network")
end

## enum name -> MeasurementType, strict (unknown names reject the file)
function _measurement_type_from_name(s::AbstractString)
  for t in instances(MeasurementType)
    string(t) == s && return t
  end
  return nothing
end

"""
    readMeasurementsCSV!(net; file, replace=true) -> NamedTuple

Read a measurement CSV v1 file (see `writeMeasurementsCSV`) into
`net.measurements`. The import is ATOMIC: the whole file is parsed and
validated first, and any error (unknown version comment, unknown type name,
malformed row, unresolvable bus/branch/link) aborts with a line-precise
message (`file:line: reason`) leaving `net.measurements` untouched; there is
never a partial import. `replace = true` (default) replaces the stored
measurements, `replace = false` appends.

Returns `(counts = Dict{MeasurementType,Int}, total, skipped)`; `skipped`
counts empty and comment lines after the version line.
"""
function readMeasurementsCSV!(net::Net; file::AbstractString, replace::Bool = true)
  lines = readlines(file)
  isempty(lines) && error("$(file): empty file, expected the version comment '# $(_MEASUREMENT_CSV_VERSION)'")
  vline = strip(lines[1])
  expected = "# " * _MEASUREMENT_CSV_VERSION
  vline == expected || error("$(file):1: unknown measurement file version ('$(vline)'; expected '$(expected)')")

  parsed = Measurement[]
  counts = Dict{MeasurementType,Int}()
  skipped = 0
  headerSeen = false
  hasBranchNr = true
  for (ln, raw) in enumerate(lines)
    ln == 1 && continue
    line = strip(raw)
    if isempty(line) || startswith(line, "#")
      skipped += 1
      continue
    end
    if !headerSeen
      if line == _MEASUREMENT_CSV_HEADER
        hasBranchNr = true
      elseif line == _MEASUREMENT_CSV_HEADER_NO_BRANCH_NR
        hasBranchNr = false
      else
        error("$(file):$(ln): unexpected header '$(line)' (expected '$(_MEASUREMENT_CSV_HEADER)')")
      end
      headerSeen = true
      continue
    end
    nfields = hasBranchNr ? 11 : 10
    fields = split(line, ","; limit = nfields)
    length(fields) == nfields || error("$(file):$(ln): expected $(nfields) comma-separated fields, got $(length(fields))")
    branchStr = ""
    if hasBranchNr
      tstr, bus, fromB, toB, branchStr, linkStr, dirStr, valStr, sigStr, actStr, id = map(_csv_field, fields)
    else
      tstr, bus, fromB, toB, linkStr, dirStr, valStr, sigStr, actStr, id = map(_csv_field, fields)
    end

    typ = _measurement_type_from_name(tstr)
    typ === nothing && error("$(file):$(ln): unknown measurement type '$(tstr)'")
    value = tryparse(Float64, valStr)
    value === nothing && error("$(file):$(ln): value '$(valStr)' is not a number")
    sigma = tryparse(Float64, sigStr)
    sigma === nothing && error("$(file):$(ln): sigma '$(sigStr)' is not a number")
    active = actStr in ("true", "1") ? true : (actStr in ("false", "0") ? false : nothing)
    active === nothing && error("$(file):$(ln): active '$(actStr)' must be true/false")

    hasBus = !isempty(bus)
    hasBranch = !isempty(fromB) || !isempty(toB) || !isempty(branchStr)
    hasLink = !isempty(linkStr)
    (hasBus + hasBranch + hasLink) == 1 || error("$(file):$(ln): exactly one location group required (bus | from_bus+to_bus | link_nr)")

    m = try
      if hasLink
        linkNr = tryparse(Int, linkStr)
        linkNr === nothing && error("link_nr '$(linkStr)' is not an integer")
        (1 <= linkNr <= length(net.linkVec)) || error("link_nr $(linkNr) not found in network")
        Measurement(typ = typ, value = value, sigma = sigma, active = active, direction = :from, id = id, linkIdx = linkNr)
      elseif hasBranch
        (isempty(fromB) || isempty(toB)) && error("branch measurement needs both from_bus and to_bus")
        fromKey = _measurement_csv_bus_key(net, fromB)
        toKey = _measurement_csv_bus_key(net, toB)
        bridx = if isempty(branchStr)
          # no branch_nr (first-cut v1 file): pair resolution, parallels reject
          _resolve_branch_idx(net; fromBus = fromKey, toBus = toKey)
        else
          # branch_nr disambiguates parallel branches; the named endpoints
          # must match that branch (either orientation) so a renumbered net
          # cannot silently misplace the measurement
          bn = tryparse(Int, branchStr)
          bn === nothing && error("branch_nr '$(branchStr)' is not an integer")
          (1 <= bn <= length(net.branchVec)) || error("branch_nr $(bn) not found in network")
          br = net.branchVec[bn]
          fi = geNetBusIdx(net = net, busName = fromKey)
          ti = geNetBusIdx(net = net, busName = toKey)
          ((Int(br.fromBus) == fi && Int(br.toBus) == ti) || (Int(br.fromBus) == ti && Int(br.toBus) == fi)) || error("branch_nr $(bn) connects other buses than $(fromB) and $(toB)")
          bn
        end
        dir = dirStr == "from" ? :from : (dirStr == "to" ? :to : error("direction '$(dirStr)' must be from or to for a branch measurement"))
        Measurement(typ = typ, value = value, sigma = sigma, active = active, branchIdx = bridx, direction = dir, id = id)
      else
        isempty(dirStr) || error("direction must be empty for a bus measurement")
        busIdx = geNetBusIdx(net = net, busName = _measurement_csv_bus_key(net, bus))
        Measurement(typ = typ, value = value, sigma = sigma, active = active, busIdx = busIdx, id = id)
      end
    catch e
      msg = e isa ErrorException ? e.msg : sprint(showerror, e)
      error("$(file):$(ln): $(msg)")
    end
    push!(parsed, m)
    counts[typ] = get(counts, typ, 0) + 1
  end
  headerSeen || error("$(file): missing header line '$(_MEASUREMENT_CSV_HEADER)'")

  # A zero-injection row is a constraint, and older sets wrote it with sigma
  # 1e-6, i.e. weight 1e12. On a large network that does not bind the bus, it
  # crushes the normal equations: a 25000-bus set with 27268 such rows did not
  # converge in 30 iterations and reported an objective of 3.8e15 per degree of
  # freedom. The same set at ZERO_INJECTION_SIGMA converges in 7 iterations to
  # J/dof 4.97 (and to 1.00 once the documented tap deviations are released).
  # Reading therefore lifts a constraint row to the current floor and says so;
  # ordinary telemetry is untouched, its tightest rows are orders of magnitude
  # above this.
  lifted = 0
  for (i, m) in enumerate(parsed)
    m.sigma < ZERO_INJECTION_SIGMA || continue
    parsed[i] = Measurement(typ = m.typ, value = m.value, sigma = ZERO_INJECTION_SIGMA, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id, active = m.active)
    lifted += 1
  end
  lifted > 0 && @info "measurement set: $(lifted) constraint row(s) had a sigma below $(ZERO_INJECTION_SIGMA) and were read at that value (a tighter weight makes the normal equations unsolvable at scale)"

  # atomic take-over only after the whole file validated
  if replace
    empty!(net.measurements)
  end
  append!(net.measurements, parsed)

  # A generated set holds each physical quantity once. Two readings of the same
  # quantity are legal in real telemetry (redundant transducers), but a set that
  # doubles EVERY quantity is the signature of a corrupted file, and it inflates
  # J without any single residual looking wrong - which then reads as a topology
  # error. Name it here instead of letting the estimate mislead.
  # checked on the RESULT, not on the parsed rows: with `replace = false` the
  # doubling arises from the combination, which is exactly how a set gets
  # corrupted (generated rows appended onto rows the case already carried)
  duplicates = _duplicate_measured_quantities(net.measurements)
  if !isempty(duplicates)
    @warn "$(basename(file)): $(sum(values(duplicates))) measurement(s) repeat a quantity already measured; a set that repeats every quantity is usually a corrupted file and inflates J." duplicates
  end
  return (counts = counts, total = length(parsed), skipped = skipped, duplicates = duplicates)
end

"""
    _duplicate_measured_quantities(ms) -> Dict{String,Int}

Count, per measurement type, how many entries repeat a quantity (same type, same
bus/branch/link, same direction) that an earlier entry already covers. The first
occurrence of each quantity is not counted, so an empty result means the set
measures every quantity at most once.

Only ACTIVE rows count: a deactivated row contributes nothing to the objective,
so it never competes with the reading it sits next to.
"""
function _duplicate_measured_quantities(ms::AbstractVector)
  seen = Set{Tuple{MeasurementType,Union{Nothing,Int},Union{Nothing,Int},Union{Nothing,Int},Symbol}}()
  dups = Dict{String,Int}()
  for m in ms
    m.active || continue
    key = (m.typ, m.busIdx, m.branchIdx, m.linkIdx, m.direction)
    if key in seen
      k = string(m.typ)
      dups[k] = get(dups, k, 0) + 1
    else
      push!(seen, key)
    end
  end
  return dups
end
