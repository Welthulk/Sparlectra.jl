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

"""
    measurementStdDevs(; vm=0.005, pinj=1.0, qinj=1.0, pflow=1.0, qflow=1.0)

Create default standard-deviation map for synthetic measurement generation.
"""
function measurementStdDevs(; vm::Float64 = 0.005, pinj::Float64 = 1.0, qinj::Float64 = 1.0, pflow::Float64 = 1.0, qflow::Float64 = 1.0)
  return Dict(VmMeas => vm, PinjMeas => pinj, QinjMeas => qinj, PflowMeas => pflow, QflowMeas => qflow)
end

@inline function _branch_flow_pu(branch::Branch, from::Int, to::Int, tapSide::Int, V::Vector{ComplexF64})
  @assert tapSide == 1 || tapSide == 2
  if branch.status == 0
    return 0.0 + 0.0im
  end

  ui = V[from]
  uj = V[to]

  ratio = (branch.ratio != 0.0) ? branch.ratio : 1.0
  angle = (branch.ratio != 0.0) ? branch.angle : 0.0
  tap = calcComplexRatio(tapRatio = ratio, angleInDegrees = angle)

  if tapSide == 1
    ui /= tap
  else
    uj /= tap
  end

  Yik = inv(branch.r_pu + im * branch.x_pu)
  Y0ik = 0.5 * (branch.g_pu + im * branch.b_pu)
  return abs(ui)^2 * conj(Y0ik + Yik) - ui * conj(uj) * conj(Yik)
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
      s = _branch_flow_pu(br, br.fromBus, br.toBus, 1, V) * net.baseMVA
    elseif meas.direction == :to
      s = _branch_flow_pu(br, br.toBus, br.fromBus, 2, V) * net.baseMVA
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
      sfrom = _branch_flow_pu(br, br.fromBus, br.toBus, 1, V) * net.baseMVA
      sto = _branch_flow_pu(br, br.toBus, br.fromBus, 2, V) * net.baseMVA
      push!(m, Measurement(typ = PflowMeas, value = real(sfrom) + (noise ? randn(rng) * σ : 0.0), sigma = σ, branchIdx = br.branchIdx, direction = :from, id = "Pflow_branch_$(br.branchIdx)_from"))
      push!(m, Measurement(typ = PflowMeas, value = real(sto) + (noise ? randn(rng) * σ : 0.0), sigma = σ, branchIdx = br.branchIdx, direction = :to, id = "Pflow_branch_$(br.branchIdx)_to"))
    end
    if includeQflow
      σ = stddev[QflowMeas]
      sfrom = _branch_flow_pu(br, br.fromBus, br.toBus, 1, V) * net.baseMVA
      sto = _branch_flow_pu(br, br.toBus, br.fromBus, 2, V) * net.baseMVA
      push!(m, Measurement(typ = QflowMeas, value = imag(sfrom) + (noise ? randn(rng) * σ : 0.0), sigma = σ, branchIdx = br.branchIdx, direction = :from, id = "Qflow_branch_$(br.branchIdx)_from"))
      push!(m, Measurement(typ = QflowMeas, value = imag(sto) + (noise ? randn(rng) * σ : 0.0), sigma = σ, branchIdx = br.branchIdx, direction = :to, id = "Qflow_branch_$(br.branchIdx)_to"))
    end
  end

  return m
end
