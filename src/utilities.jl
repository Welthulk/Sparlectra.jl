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
