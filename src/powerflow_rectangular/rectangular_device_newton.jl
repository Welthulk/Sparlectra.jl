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

# file: src/powerflow_rectangular/rectangular_device_newton.jl
# purpose: the seam for a Newton loop that runs on the GPU (mismatch,
#          Jacobian build, factorization, solve and voltage update in device
#          memory, ext/SparlectraCUDAExt.jl). Study switch only, no
#          configuration key: without _with_device_newton every run takes
#          the CPU loop of runpf_rectangular! unchanged.

# Study switch: true inside _with_device_newton. The CPU loop of
# runpf_rectangular! hands the solve to the device loop only when it is
# set and the request check below accepts the solve.
const _DEVICE_NEWTON = Base.ScopedValues.ScopedValue(false)

# run `f` with the device Newton loop allowed for every GPU solve inside it
_with_device_newton(@nospecialize(f)) = Base.ScopedValues.with(f, _DEVICE_NEWTON => true)

"""
    DeviceNewtonInput

Everything the device Newton loop needs from `runpf_rectangular!`. The
arrays `S`, `bus_types`, `Vset`, `dPinj_dVm`, `dQinj_dVm` and `workspace`
are the solve's own: the device loop changes them in place exactly where
the CPU loop changes them (Q-limit switches, the injection rebuild of
voltage-dependent controllers), so the code after the loop and a CPU loop
that continues the solve read the same state. `V` is the start iterate,
`V0[slack_idx]` the fixed reference voltage. `qlimit_args` carries the
arguments of `_handle_rectangular_qlimit_iteration!` after `history`, in
call order. `linear_ctx` is the GPU context of the solve; its counters and
timings are updated by the loop.
"""
struct DeviceNewtonInput
  net::Net
  Ybus::SparseMatrixCSC{ComplexF64,Int}
  V::Vector{ComplexF64}
  V0::Vector{ComplexF64}
  S::Vector{ComplexF64}
  bus_types::Vector{Symbol}
  Vset::Vector{Float64}
  slack_idx::Int
  maxiter::Int
  tol::Float64
  damp::Float64
  autodamp_min::Float64
  newton_update::Symbol
  has_vdep_control::Bool
  isolated_rows::Vector{Int}
  dPinj_dVm::Vector{Float64}
  dQinj_dVm::Vector{Float64}
  qlimits_enabled::Bool
  qlimit_args::NamedTuple
  workspace::RectangularIterationWorkspace
  performance_profile::Any
  linear_ctx::AbstractGpuNewtonContext
  check_cancel::Any
  set_phase::Any
  verbose::Int
end

"""
    DeviceNewtonTimings

Times of one device Newton solve in seconds. Device phases (`mismatch_s`,
`build_s`, `refactor_s`, `solve_s`, `update_s`) come from CUDA event pairs
read after the loop; `map_s`, `upload_s`, `analysis_s` and
`factor_first_s` are host clocks around synchronized work before the first
iteration; `d2h_s` and `h2d_s` are the host copies the loop issues itself
after the pre-flight (`d2h_bytes`, `h2d_bytes` and `d2h_count`, `h2d_count`
count them; copies inside libcudss are not counted), `host_qlimit_s` the
CPU Q-limit handler, `host_injections_s` the rebuild of voltage-dependent
injections, `core_s` the whole device loop.
`first_call` is true for the first device solve of the process (kernel
compilation included).
"""
mutable struct DeviceNewtonTimings
  first_call::Bool
  map_s::Float64
  upload_s::Float64
  analysis_s::Float64
  factor_first_s::Float64
  mismatch_s::Float64
  build_s::Float64
  refactor_s::Float64
  solve_s::Float64
  update_s::Float64
  d2h_s::Float64
  h2d_s::Float64
  host_qlimit_s::Float64
  host_injections_s::Float64
  core_s::Float64
  launches_per_iteration::Int
  d2h_bytes::Int
  h2d_bytes::Int
  d2h_count::Int
  h2d_count::Int
end

DeviceNewtonTimings() = DeviceNewtonTimings(false, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0, 0)

"""
    DeviceNewtonOutcome

What the device loop hands back to `runpf_rectangular!`: the iterate `V`,
the loop locals the code after the Newton loop reads, and `next_it`, the
first iteration the CPU loop still has to run (`maxiter + 1` when the
device loop finished the solve; smaller after a GPU error or a failed GPU
step, when the CPU continues from the iterate read back from the device,
or from the start iterate if the CUDA context is lost). The
history vectors and `step_diagnostics` belong to the device loop and are
appended to the solve's own vectors by the caller.
"""
struct DeviceNewtonOutcome
  V::Vector{ComplexF64}
  next_it::Int
  iters::Int
  converged::Bool
  rejection_reason::Symbol
  history::Vector{Float64}
  pre_switch_history::Vector{Float64}
  type_switched::Vector{Bool}
  step_diagnostics::Vector{NamedTuple}
  best_finite_mismatch::Float64
  best_finite_iteration::Int
  best_finite_voltage::Union{Nothing,Vector{ComplexF64}}
  last_finite_iteration::Int
  last_finite_voltage::Union{Nothing,Vector{ComplexF64}}
  first_nonfinite_iteration::Int
  first_nonfinite_voltage::Union{Nothing,Vector{ComplexF64}}
  qlimit_active_set_changes::Int
  qlimit_reenable_events::Int
  qlimit_switch_iterations::Int
  newton_steps::Int
end

# Whether the device loop may take this solve, and if not, why. The GPU
# path covers the default Newton: fixed step, no start projection, no
# distributed slack, no Jacobian reuse, Q limits by PV->PQ switching.
# Everything else stays on the CPU loop (Jacobian and mismatch on the CPU,
# linear solve as the compute backend decided).
function _device_newton_request(linear_ctx; power_mode::Bool, autodamp::Bool, merit_enabled::Bool, trust_region_enabled::Bool, start_projection::Bool, distributed_slack::Bool, jacobian_reuse::Bool, qlimit_mode::Symbol, isolated_rows::Vector{Int})::NamedTuple{(:use, :reason),Tuple{Bool,String}}
  linear_ctx isa AbstractGpuNewtonContext || return (use = false, reason = "the linear solve does not run on the GPU")
  linear_ctx.cpu === nothing || return (use = false, reason = "the GPU context already fell back to the CPU")
  backend = _GPU_BACKEND[]
  (backend !== nothing && hasproperty(backend, :newton)) || return (use = false, reason = "the CUDA extension has no device Newton loop")
  power_mode && return (use = false, reason = "power mode keeps its own backend")
  (autodamp || merit_enabled) && return (use = false, reason = "damping (autodamp or merit line search) runs on the CPU")
  trust_region_enabled && return (use = false, reason = "the trust region runs on the CPU")
  start_projection && return (use = false, reason = "the start projection runs on the CPU")
  distributed_slack && return (use = false, reason = "distributed slack is not in the device map")
  jacobian_reuse && return (use = false, reason = "Jacobian reuse (chord steps) is not in the device loop")
  qlimit_mode === :switch_to_pq || return (use = false, reason = "qlimit_mode $(qlimit_mode) is not in the device loop")
  isempty(isolated_rows) || return (use = false, reason = "isolated buses are not in the device map")
  return (use = true, reason = "")
end

"""
    _run_device_newton(input::DeviceNewtonInput) -> Union{Nothing,DeviceNewtonOutcome}

Run the Newton loop of one solve on the GPU through the CUDA extension.
`nothing` means the device loop did not start (a GPU error before the
first iteration with `gpu_fallback_to_cpu` on; the reason is in the
context's `fallback_reason`) and no solver state was changed, so the CPU
loop runs the whole solve. A GPU error with fallback off is rethrown.
"""
function _run_device_newton(input::DeviceNewtonInput)
  return _GPU_BACKEND[].newton(input)
end

# One Newton step on the CPU from `V` with the current state, for a step
# the GPU did not finish: a GPU error after the Q-limit handling ran, a
# GPU error before the step was checked, or a non-finite step. `linear_ctx`
# has its `cpu` context set already, so the solve goes to UMFPACK; a singular Jacobian
# there is the CPU loop's :singular_newton_step. Returns the new iterate
# (slack reset as in the CPU loop) or `nothing` for a singular step.
function _device_newton_cpu_step(input::DeviceNewtonInput, V::Vector{ComplexF64}, step_diagnostics, active_set_changed::Bool)
  try
    Vnew = complex_newton_step_rectangular(
      input.Ybus,
      V,
      input.S;
      slack_idx = input.slack_idx,
      damp = input.damp,
      autodamp = false,
      autodamp_min = input.autodamp_min,
      newton_update = input.newton_update,
      bus_types = input.bus_types,
      Vset = input.Vset,
      dPinj_dVm = input.dPinj_dVm,
      dQinj_dVm = input.dQinj_dVm,
      performance_profile = input.performance_profile,
      step_diagnostics = step_diagnostics,
      active_set_changed = active_set_changed,
      linear_ctx = input.linear_ctx,
    )
    Vnew[input.slack_idx] = input.V0[input.slack_idx]
    return Vnew
  catch step_error
    _is_rectangular_linear_step_failure(step_error) || rethrow(step_error)
    input.verbose > 0 && @warn "Rectangular Newton step failed because the linear Jacobian solve was singular; returning non-convergence." exception = (typeof(step_error), sprint(showerror, step_error))
    return nothing
  end
end
