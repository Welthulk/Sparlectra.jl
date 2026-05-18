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

# file: src/utilities.jl
# ------------------------------
# Utility: set a full row of a CSC matrix to zero
# ------------------------------

function zero_row!(J::AbstractMatrix, r::Int)
  @views J[r, :] .= zero(eltype(J))
  return J
end

function print_jacobian(J::AbstractMatrix; label::AbstractString = "Jacobian", digits::Int = 3, states_per_bus::Int = 2)
  nrows, ncols = size(J)
  @assert nrows == ncols "Jacobian must be square."
  @assert nrows % states_per_bus == 0 "Matrix size not compatible with states_per_bus."

  nbuses = nrows ÷ states_per_bus

  println(label, ": size = ", nrows, " × ", ncols, " (", nbuses, " buses, ", states_per_bus, " states per bus)")

  println("\nJacobian (rounded, bus-block view):")

  # Column header
  print("i\\j\t")
  for b = 1:nbuses
    for s = 1:states_per_bus
      print("B", b, "_", s, '\t')
    end
  end
  println()

  # Rows
  for bi = 1:nbuses
    for si = 1:states_per_bus
      row = (bi - 1) * states_per_bus + si
      print("B", bi, "_", si, '\t')

      for bj = 1:nbuses
        for sj = 1:states_per_bus
          col = (bj - 1) * states_per_bus + sj
          val = round(J[row, col]; digits = digits)
          print(val, '\t')
        end
      end
      println()
    end
  end
end

function _perf_profile_enabled(profile)::Bool
  profile isa AbstractDict || return false
  return Bool(get(profile, :enabled, false))
end

function _perf_profile_wants_allocations(profile)::Bool
  _perf_profile_enabled(profile) || return false
  return Bool(get(profile, :show_allocations, false))
end

function _perf_profile_add!(profile, phase::Symbol, elapsed_s::Real, bytes::Integer = 0)
  _perf_profile_enabled(profile) || return nothing
  timings = get!(profile, :timings) do
    Dict{Symbol,NamedTuple{(:calls,:elapsed_s,:bytes),Tuple{Int,Float64,Int}}}()
  end
  prev = get(timings, phase, (calls = 0, elapsed_s = 0.0, bytes = 0))
  timings[phase] = (calls = prev.calls + 1, elapsed_s = prev.elapsed_s + Float64(elapsed_s), bytes = prev.bytes + Int(bytes))
  return nothing
end

function _perf_profile_push_iteration!(profile, row)
  _perf_profile_enabled(profile) || return nothing
  push!(get!(profile, :iterations, NamedTuple[]), row)
  return nothing
end

function _perf_profile_time!(f::F, profile, phase::Symbol) where {F}
  if !_perf_profile_enabled(profile)
    return f()
  elseif _perf_profile_wants_allocations(profile)
    timed = @timed f()
    _perf_profile_add!(profile, phase, timed.time, timed.bytes)
    return timed.value
  else
    t0 = time_ns()
    value = f()
    _perf_profile_add!(profile, phase, (time_ns() - t0) / 1e9, 0)
    return value
  end
end
