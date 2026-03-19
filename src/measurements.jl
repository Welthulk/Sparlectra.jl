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

# file: src/measurements.jl

using Random

@inline function _net_measurements(net::Net)::Vector
  return net.measurements
end

"""
    @enum MeasurementType

Supported measurement types for the first WLS state-estimation implementation.
"""
@enum MeasurementType begin
  VmMeas
  PinjMeas
  QinjMeas
  PflowMeas
  QflowMeas
end

"""
    Measurement

Generic state-estimation measurement model.

Fields:
- `typ`: Measurement type.
- `value`: Measured value (Vm in p.u., powers in MW/MVar).
- `sigma`: Standard deviation in measurement units.
- `weight`: Weight used in WLS (`1/sigma^2`).
- `active`: If `false`, measurement is ignored by estimator.
- `busIdx`: Bus index for bus measurements.
- `branchIdx`: Branch index for branch flow measurements.
- `direction`: Branch direction `:from` or `:to`, otherwise `:none`.
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
end

function Measurement(; typ::MeasurementType, value::Real, sigma::Real, active::Bool = true, busIdx::Union{Nothing,Int} = nothing, branchIdx::Union{Nothing,Int} = nothing, direction::Symbol = :none, id::AbstractString = "")
  σ = Float64(sigma)
  σ <= 0.0 && error("Measurement sigma must be > 0.0")
  w = inv(σ * σ)
  return Measurement(typ, Float64(value), σ, w, active, busIdx, branchIdx, direction, String(id))
end

@inline function _default_measurement_id(typ::MeasurementType; busIdx::Union{Nothing,Int} = nothing, branchIdx::Union{Nothing,Int} = nothing, direction::Symbol = :none)
  if typ == VmMeas
    return "Vm_bus_$(something(busIdx, 0))"
  elseif typ == PinjMeas
    return "Pinj_bus_$(something(busIdx, 0))"
  elseif typ == QinjMeas
    return "Qinj_bus_$(something(busIdx, 0))"
  elseif typ == PflowMeas
    return "Pflow_branch_$(something(branchIdx, 0))_$(direction)"
  elseif typ == QflowMeas
    return "Qflow_branch_$(something(branchIdx, 0))_$(direction)"
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
function addMeasurement!(measurements::Vector; typ::MeasurementType, value::Real, sigma::Real, active::Bool = true, busIdx::Union{Nothing,Int} = nothing, branchIdx::Union{Nothing,Int} = nothing, direction::Symbol = :none, id::AbstractString = "")
  if typ == VmMeas || typ == PinjMeas || typ == QinjMeas
    isnothing(busIdx) && error("Bus measurement requires busIdx")
  elseif typ == PflowMeas || typ == QflowMeas
    isnothing(branchIdx) && error("Flow measurement requires branchIdx")
    (direction == :from || direction == :to) || error("Flow measurement direction must be :from or :to")
  end

  meas_id = isempty(id) ? _default_measurement_id(typ; busIdx = busIdx, branchIdx = branchIdx, direction = direction) : String(id)
  meas = Measurement(typ = typ, value = value, sigma = sigma, active = active, busIdx = busIdx, branchIdx = branchIdx, direction = direction, id = meas_id)
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

"""
    addPflowMeasurement!(measurements; net, value, sigma, direction=:from, branchNr=nothing, fromBus=nothing, toBus=nothing, active=true, id="")

Append an active-power flow measurement identified by `branchNr` or a unique
`fromBus`/`toBus` branch pair.
"""
function addPflowMeasurement!(measurements::Vector; net::Net, value::Real, sigma::Real, direction::Symbol = :from, branchNr::Union{Nothing,Int} = nothing, fromBus::Union{Nothing,String} = nothing, toBus::Union{Nothing,String} = nothing, active::Bool = true, id::AbstractString = "")
  bridx = _resolve_branch_idx(net; branchNr = branchNr, fromBus = fromBus, toBus = toBus)
  return addMeasurement!(measurements; typ = PflowMeas, value = value, sigma = sigma, active = active, branchIdx = bridx, direction = direction, id = id)
end

function addPflowMeasurement!(net::Net; value::Real, sigma::Real, direction::Symbol = :from, branchNr::Union{Nothing,Int} = nothing, fromBus::Union{Nothing,String} = nothing, toBus::Union{Nothing,String} = nothing, active::Bool = true, id::AbstractString = "")
  return addPflowMeasurement!(_net_measurements(net); net = net, value = value, sigma = sigma, direction = direction, branchNr = branchNr, fromBus = fromBus, toBus = toBus, active = active, id = id)
end

"""
    addQflowMeasurement!(measurements; net, value, sigma, direction=:from, branchNr=nothing, fromBus=nothing, toBus=nothing, active=true, id="")

Append a reactive-power flow measurement identified by `branchNr` or a unique
`fromBus`/`toBus` branch pair.
"""
function addQflowMeasurement!(measurements::Vector; net::Net, value::Real, sigma::Real, direction::Symbol = :from, branchNr::Union{Nothing,Int} = nothing, fromBus::Union{Nothing,String} = nothing, toBus::Union{Nothing,String} = nothing, active::Bool = true, id::AbstractString = "")
  bridx = _resolve_branch_idx(net; branchNr = branchNr, fromBus = fromBus, toBus = toBus)
  return addMeasurement!(measurements; typ = QflowMeas, value = value, sigma = sigma, active = active, branchIdx = bridx, direction = direction, id = id)
end

function addQflowMeasurement!(net::Net; value::Real, sigma::Real, direction::Symbol = :from, branchNr::Union{Nothing,Int} = nothing, fromBus::Union{Nothing,String} = nothing, toBus::Union{Nothing,String} = nothing, active::Bool = true, id::AbstractString = "")
  return addQflowMeasurement!(_net_measurements(net); net = net, value = value, sigma = sigma, direction = direction, branchNr = branchNr, fromBus = fromBus, toBus = toBus, active = active, id = id)
end

"""
    measurementStdDevs(; vm=0.005, pinj=1.0, qinj=1.0, pflow=1.0, qflow=1.0)

Create default standard-deviation map for synthetic measurement generation.
"""
function measurementStdDevs(; vm::Float64 = 0.005, pinj::Float64 = 1.0, qinj::Float64 = 1.0, pflow::Float64 = 1.0, qflow::Float64 = 1.0)
  return Dict(VmMeas => vm, PinjMeas => pinj, QinjMeas => qinj, PflowMeas => pflow, QflowMeas => qflow)
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
    addZeroInjectionMeasurements!(measurements; net, sigma=1e-6, busNames=nothing, busIdxs=nothing, active=true, idPrefix="ZI") -> Vector{Measurement}

Append active- and reactive-power zero-injection pseudo-measurements for the
selected buses and return the newly added measurements.

Selection rules:
- If `busIdxs` is provided, those indices are used.
- Else if `busNames` is provided, names are resolved to indices.
- Else passive buses are detected automatically via `findPassiveBuses(net)`.

These pseudo-measurements are the current way to encode equality constraints
`P_inj = 0` and `Q_inj = 0` in the WLS estimator.
"""
function addZeroInjectionMeasurements!(measurements::Vector; net::Net, sigma::Real = 1e-6, busNames::Union{Nothing,Vector{String}} = nothing, busIdxs::Union{Nothing,Vector{Int}} = nothing, active::Bool = true, idPrefix::AbstractString = "ZI")
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

function addZeroInjectionMeasurements!(net::Net; sigma::Real = 1e-6, busNames::Union{Nothing,Vector{String}} = nothing, busIdxs::Union{Nothing,Vector{Int}} = nothing, active::Bool = true, idPrefix::AbstractString = "ZI")
  return addZeroInjectionMeasurements!(_net_measurements(net); net = net, sigma = sigma, busNames = busNames, busIdxs = busIdxs, active = active, idPrefix = idPrefix)
end

@inline function _measurement_prediction(meas::Measurement, net::Net, V::Vector{ComplexF64}, Sbus_MVA::Vector{ComplexF64})
  if meas.typ == VmMeas
    i = something(meas.busIdx, 0)
    i < 1 && error("Vm measurement missing busIdx")
    return abs(V[i])
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
    br = net.branchVec[bridx]
    if meas.direction == :from
      s = branchFlow_pu(br, br.fromBus, br.toBus, 1, V) * net.baseMVA
    elseif meas.direction == :to
      s = branchFlow_pu(br, br.toBus, br.fromBus, 2, V) * net.baseMVA
    else
      error("Flow measurement direction must be :from or :to")
    end
    return meas.typ == PflowMeas ? real(s) : imag(s)
  else
    error("Unsupported measurement type")
  end
end

"""
    generateMeasurementsFromPF(net; kwargs...) -> Vector{Measurement}

Generate synthetic measurements from the current solved network state.

Keyword options:
- `includeVm`, `includePinj`, `includeQinj`, `includePflow`, `includeQflow`
- `noise`: add Gaussian noise if `true`
- `stddev`: dictionary from `MeasurementType => sigma`
- `rng`: random number generator
"""
function generateMeasurementsFromPF(net::Net; includeVm::Bool = true, includePinj::Bool = true, includeQinj::Bool = true, includePflow::Bool = true, includeQflow::Bool = true, noise::Bool = false, stddev::Dict{MeasurementType,Float64} = measurementStdDevs(), rng::AbstractRNG = Random.default_rng())
  V = buildVoltageVector(net)
  Sbus_pu = calc_injections(createYBUS(net = net), V)
  Sbus_MVA = Sbus_pu .* net.baseMVA

  m = Measurement[]
  nbus = length(net.nodeVec)

  for i = 1:nbus
    if includeVm
      σ = stddev[VmMeas]
      z = abs(V[i]) + (noise ? randn(rng) * σ : 0.0)
      push!(m, Measurement(typ = VmMeas, value = z, sigma = σ, busIdx = i, id = "Vm_bus_$(i)"))
    end
    if includePinj
      σ = stddev[PinjMeas]
      z = real(Sbus_MVA[i]) + (noise ? randn(rng) * σ : 0.0)
      push!(m, Measurement(typ = PinjMeas, value = z, sigma = σ, busIdx = i, id = "Pinj_bus_$(i)"))
    end
    if includeQinj
      σ = stddev[QinjMeas]
      z = imag(Sbus_MVA[i]) + (noise ? randn(rng) * σ : 0.0)
      push!(m, Measurement(typ = QinjMeas, value = z, sigma = σ, busIdx = i, id = "Qinj_bus_$(i)"))
    end
  end

  for br in net.branchVec
    br.status == 0 && continue
    if includePflow
      σ = stddev[PflowMeas]
      sfrom = branchFlow_pu(br, br.fromBus, br.toBus, 1, V) * net.baseMVA
      sto = branchFlow_pu(br, br.toBus, br.fromBus, 2, V) * net.baseMVA
      push!(m, Measurement(typ = PflowMeas, value = real(sfrom) + (noise ? randn(rng) * σ : 0.0), sigma = σ, branchIdx = br.branchIdx, direction = :from, id = "Pflow_branch_$(br.branchIdx)_from"))
      push!(m, Measurement(typ = PflowMeas, value = real(sto) + (noise ? randn(rng) * σ : 0.0), sigma = σ, branchIdx = br.branchIdx, direction = :to, id = "Pflow_branch_$(br.branchIdx)_to"))
    end
    if includeQflow
      σ = stddev[QflowMeas]
      sfrom = branchFlow_pu(br, br.fromBus, br.toBus, 1, V) * net.baseMVA
      sto = branchFlow_pu(br, br.toBus, br.fromBus, 2, V) * net.baseMVA
      push!(m, Measurement(typ = QflowMeas, value = imag(sfrom) + (noise ? randn(rng) * σ : 0.0), sigma = σ, branchIdx = br.branchIdx, direction = :from, id = "Qflow_branch_$(br.branchIdx)_from"))
      push!(m, Measurement(typ = QflowMeas, value = imag(sto) + (noise ? randn(rng) * σ : 0.0), sigma = σ, branchIdx = br.branchIdx, direction = :to, id = "Qflow_branch_$(br.branchIdx)_to"))
    end
  end

  return m
end

function setMeasurementsFromPF!(net::Net; kwargs...)
  measurements = generateMeasurementsFromPF(net; kwargs...)
  empty!(net.measurements)
  append!(net.measurements, measurements)
  return Measurement[m for m in net.measurements]
end
