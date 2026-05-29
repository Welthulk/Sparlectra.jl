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
#
# This file is included inside module Sparlectra. Do not add a module wrapper here.
#
# Standalone rectangular array-level Newton-Raphson solver.

# Date: 29.5.2026
# file: src/powerflow_rectangular/rectangular_standalone_solver.jl

function run_complex_nr_rectangular(
  Ybus,
  V0,
  S;
  slack_idx::Int = 1,
  maxiter::Int = 20,
  tol::Float64 = 1e-8,
  verbose::Bool = false,
  damp::Float64 = 1.0,
  autodamp::Bool = false,
  autodamp_min::Float64 = 0.05,
  wrong_branch_detection::Symbol = :warn,
  wrong_branch_rescue::Bool = false,
  wrong_branch_min_vm_pu::Float64 = 0.70,
  wrong_branch_max_vm_pu::Float64 = 1.30,
  wrong_branch_max_angle_spread_deg::Float64 = 180.0,
  wrong_branch_max_branch_angle_deg::Float64 = 90.0,
  wrong_branch_min_low_vm_count::Int = 1,
  wrong_branch_rescue_max_attempts::Int = 2,
  bus_types::Vector{Symbol},
  Vset::Vector{Float64},
  dPinj_dVm::Vector{Float64} = zeros(Float64, length(V0)),
  dQinj_dVm::Vector{Float64} = zeros(Float64, length(V0)),
  performance_profile = nothing,
)
  V = copy(V0)
  history = Float64[]

  for iter = 1:maxiter
    F = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
    max_mis = maximum(abs.(F))
    push!(history, max_mis)

    if max_mis <= tol
      return V, true, iter, history
    end

    V = complex_newton_step_rectangular(Ybus, V, S; slack_idx = slack_idx, damp = damp, autodamp = autodamp, autodamp_min = autodamp_min, bus_types = bus_types, Vset = Vset, dPinj_dVm = dPinj_dVm, dQinj_dVm = dQinj_dVm, performance_profile = performance_profile)
  end

  return V, false, maxiter, history
end
