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

# file: src/examples/h_matrix_observability_demo.jl

"""
Demo for small measurement Jacobians H (m×n).

This example reuses public observability helpers from `state_estimation.jl`:
- `evaluate_observability_matrix`
- `evaluate_local_observability_matrix`
- `numerical_observable` / `structural_observable`
- `numerical_row_redundant` / `structural_row_redundant`

Run from project root:
  julia --project=. src/examples/h_matrix_observability_demo.jl
"""

using Sparlectra
using Printf

abstract type Meas end
struct FlowMeas <: Meas
  a::Int
  b::Int
end
struct InjMeas <: Meas
  k::Int
end

"""Toy spanning-tree-like reachability with optional forbidden measurement index."""
function can_build_tree(nbus::Int, slack::Int, meas::Vector{Meas}; forbidden::Int = 0)
  reached = falses(nbus)
  reached[slack] = true
  inj_available = falses(nbus)

  changed = true
  while changed
    changed = false

    for (idx, m) in enumerate(meas)
      idx == forbidden && continue
      if m isa InjMeas
        k = (m::InjMeas).k
        if reached[k] && !inj_available[k]
          inj_available[k] = true
          changed = true
        end
      end
    end

    for (idx, m) in enumerate(meas)
      idx == forbidden && continue
      if m isa FlowMeas
        a = (m::FlowMeas).a
        b = (m::FlowMeas).b
        if reached[a] ⊻ reached[b]
          reached[a] = true
          reached[b] = true
          changed = true
        end
      end
    end

    if !changed
      for k in eachindex(reached)
        if reached[k] && inj_available[k]
          for (idx, m) in enumerate(meas)
            idx == forbidden && continue
            if m isa FlowMeas
              a = (m::FlowMeas).a
              b = (m::FlowMeas).b
              if a == k && !reached[b]
                reached[b] = true
                inj_available[k] = false
                changed = true
                break
              elseif b == k && !reached[a]
                reached[a] = true
                inj_available[k] = false
                changed = true
                break
              end
            end
          end
        end
        changed && break
      end
    end
  end

  return all(reached)
end

"""Toy graph redundancy by disabling one measurement."""
function graph_row_redundant(nbus::Int, slack::Int, meas::Vector{Meas}, i::Int)
  @assert 1 <= i <= length(meas)
  return can_build_tree(nbus, slack, meas; forbidden = i)
end

function print_row_flags(title::String, flags::Vector{Bool})
  println(title)
  for (i, f) in enumerate(flags)
    @printf("  row %2d : %s\n", i, f ? "REDUNDANT" : "CRITICAL")
  end
  return nothing
end

function analyze_matrix(name::String, H::AbstractMatrix{<:Real}; tol = nothing)
  m, n = size(H)
  println("============================================================")
  println("Matrix: $name  (m=$m, n=$n)")
  println("H =")
  display(H)

  obs = evaluate_observability_matrix(H; tol = tol)
  println("\nGlobal observability:")
  @printf("  Structural observable : %s\n", obs.structural_observable ? "YES" : "NO")
  @printf("  Numerical observable  : %s\n", obs.numerical_observable ? "YES" : "NO")

  s_red = [structural_row_redundant(H, i) for i in axes(H, 1)]
  n_red = [numerical_row_redundant(H, i; tol = tol) for i in axes(H, 1)]
  println("\nRow redundancy classification:")
  print_row_flags("Structural (matching) redundancy:", s_red)
  print_row_flags("Numerical (rank) redundancy:", n_red)
  return nothing
end

function analyze_local(name::String, H::AbstractMatrix{<:Real}, cols::Vector{Int}; tol = nothing)
  obs = evaluate_local_observability_matrix(H, cols; tol = tol)
  Hloc = H[obs.rows, cols]

  println("------------------------------------------------------------")
  println("Local subproblem: $name")
  @printf("  Selected columns: %s\n", string(cols))
  @printf("  Rows touching cols: %s\n", string(obs.rows))
  println("  H_local =")
  display(Hloc)

  println("  Local observability:")
  @printf("    Structural: %s\n", obs.structural_observable ? "YES" : "NO")
  @printf("    Numerical : %s\n", obs.numerical_observable ? "YES" : "NO")

  s_red = [structural_row_redundant(Hloc, i) for i in axes(Hloc, 1)]
  n_red = [numerical_row_redundant(Hloc, i; tol = tol) for i in axes(Hloc, 1)]
  print_row_flags("  Structural local row redundancy:", s_red)
  print_row_flags("  Numerical  local row redundancy:", n_red)
  return nothing
end

function main()
  tol = 1e-10

  # H_A: clearly observable with intentional duplicate information.
  # Expected global result:
  # - Observable (structural + numerical), because rows 1..3 already span all 3 states.
  # - Row 3 is CRITICAL (it is the only row touching state x3).
  # - Rows 1,2,4,5 are REDUNDANT (alternative rows still keep full column rank/coverage).
  H_A = [
    1.0 0.0 0.0
    0.0 1.0 0.0
    0.0 0.0 1.0
    1.0 1.0 0.0
    1.0 1.0 0.0
  ]

  # H_B: minimal square identity case.
  # Expected global result:
  # - Observable.
  # - Every row is CRITICAL, because removing any row immediately loses one state equation.
  H_B = [
    1.0 0.0 0.0
    0.0 1.0 0.0
    0.0 0.0 1.0
  ]

  # H_C: same structural pattern as H_A for rows 1..5, but rows 4 and 5 are almost identical.
  # Expected global result with tol=1e-10:
  # - Still observable.
  # - Similar critical/redundant pattern as H_A in this tolerance regime.
  # Why this is useful: by changing tol, numerical redundancy can flip while structural
  # redundancy stays unchanged (sparsity-only check vs value-dependent check).
  ϵ = 1e-12
  H_C = [
    1.0 0.0 0.0
    0.0 1.0 0.0
    0.0 0.0 1.0
    1.0 1.0 0.0
    1.0 1.0+ϵ 0.0
  ]

  # H_D: sparse 6x5 case with mixed redundancy.
  # Expected global result:
  # - Observable.
  # - Row 3 is CRITICAL (only direct support for state x3).
  # - Row 4 tends to be CRITICAL (needed to anchor x4/x5 coupling path).
  # - Additional rows (especially row 6) provide alternative support and can be redundant.
  # Structural and numerical classification can differ for some rows (e.g. row 5), which is
  # exactly the didactic point of running both checks.
  H_D = [
    1.0 0.0 0.0 0.0 1.0
    0.0 1.0 0.0 0.0 0.0
    0.0 0.0 1.0 0.0 0.0
    0.0 0.0 0.0 1.0 1.0
    1.0 1.0 0.0 0.0 0.0
    1.0 1.0 0.0 0.0 1.0
  ]

  # --- Toy graph example E --------------------------------------------------
  # nbus/slack/measE describe a tiny measurement set in graph form:
  # - buses: 1..4
  # - slack bus: 1 (start node for reachability expansion)
  # - measurements:
  #   1) Flow 1-2, 2) Flow 2-3, 3) Flow 3-4, 4) Injection at bus 2, 5) Flow 1-3
  #
  # Intuition / expected toy-graph redundancy:
  # - Meas 3 (Flow 3-4) is CRITICAL: without it bus 4 cannot be reached at all.
  # - Meas 1,2,5 are mostly REDUNDANT because the triangle-like alternative path
  #   (1-3 plus 2-3) can still connect buses 2 and 3 from slack 1.
  # - Meas 4 (Inj at 2) acts as an optional "joker" in this toy rule and is typically redundant.
  nbus = 4
  slack = 1
  measE = Meas[FlowMeas(1, 2), FlowMeas(2, 3), FlowMeas(3, 4), InjMeas(2), FlowMeas(1, 3)]

  # H_E is an incidence-like Jacobian-style matrix tied to measE row-by-row.
  # Column convention: x = [x1, x2, x3, x4] (bus-related state placeholders).
  # Row mapping:
  #   row 1 <-> FlowMeas(1,2): +x1 -x2
  #   row 2 <-> FlowMeas(2,3): +x2 -x3
  #   row 3 <-> FlowMeas(3,4): +x3 -x4   (typically critical for bus/state 4)
  #   row 4 <-> InjMeas(2):    +x2       (local support at bus/state 2)
  #   row 5 <-> FlowMeas(1,3): +x1 -x3   (alternative path information)
  #
  # This is intentionally illustrative (not a full physical SE Jacobian), but it allows
  # direct comparison between graph intuition (measE) and matrix-based checks (H_E).
  H_E = [
    1.0 -1.0 0.0 0.0
    0.0 1.0 -1.0 0.0
    0.0 0.0 1.0 -1.0
    0.0 1.0 0.0 0.0
    1.0 0.0 -1.0 0.0
  ]

  analyze_matrix("A (redundant duplicate rows)", H_A; tol = tol)
  analyze_local("A local around states 1..2", H_A, [1, 2]; tol = tol)

  analyze_matrix("B (minimal, all critical)", H_B; tol = tol)
  analyze_local("B local around states 2..3", H_B, [2, 3]; tol = tol)

  analyze_matrix("C (structural OK, numerically fragile)", H_C; tol = tol)
  analyze_local("C local around states 1..2", H_C, [1, 2]; tol = tol)

  analyze_matrix("D (sparse, one extra row)", H_D; tol = tol)
  analyze_local("D local around states 1..2", H_D, [1, 2]; tol = tol)

  println("============================================================")
  println("Graph / spanning-tree-style observability (toy) for E")
  println("Measurements (indexed):")
  for (i, m) in enumerate(measE)
    if m isa FlowMeas
      mm = m::FlowMeas
      @printf("  %2d: FlowMeas(%d,%d)\n", i, mm.a, mm.b)
    else
      mm = m::InjMeas
      @printf("  %2d: InjMeas(%d)\n", i, mm.k)
    end
  end

  global_ok = can_build_tree(nbus, slack, measE)
  @printf("\nGlobal (toy) reachable from slack=%d: %s\n", slack, global_ok ? "YES" : "NO")

  flags = [graph_row_redundant(nbus, slack, measE, i) for i in eachindex(measE)]
  print_row_flags("Toy graph redundancy:", flags)

  println("\nFor comparison, run structural/numerical on the illustrative H_E:")
  analyze_matrix("E (incidence-like H for the toy graph)", H_E; tol = tol)
  analyze_local("E local around buses 2..4", H_E, [2, 3, 4]; tol = tol)

  return nothing
end

main()
