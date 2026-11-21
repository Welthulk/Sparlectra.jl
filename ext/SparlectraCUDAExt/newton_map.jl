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
#
# file: ext/SparlectraCUDAExt/newton_map.jl
# purpose: the index map from Ybus nonzeros to positions in the CSR value
#          vector of the Jacobian, built on the CPU once per Ybus pattern
#          and slack, uploaded once. The pattern is the CPU structural
#          pattern with every non-slack bus PQ; a PV bus keeps its full Q
#          row with explicit zeros off the diagonal, so a PV/PQ switch
#          changes values only and cuDSS analyses once per solve.

# Layout (n buses, slack s, non-slack rank pos[i], nvar = 2(n-1)):
#   columns  Vr(j) = pos[j], Vi(j) = (n-1) + pos[j]
#   rows     P(i) = 2 pos[i] - 1, Q(i) = 2 pos[i]   (as mismatch_rectangular)
#   C_i      sorted({j in row i of Ybus, j != s} union {i}), m_i = |C_i|
#   row P(i) holds Vr(c) for c in C_i, then Vi(c) for c in C_i; row Q(i)
#   the same columns. base_i is the 1-based start of row P(i); row Q(i)
#   starts at base_i + 2 m_i.
#   Ybus entry k = (i, j), i, j != s, j != i: jpos[k] = base_i + rank0(j in C_i)
#   is (P(i), Vr(j)); + m_i is (P(i), Vi(j)); + 2 m_i (Q(i), Vr(j)); + 3 m_i
#   (Q(i), Vi(j)). dpos[i] the same for the diagonal block. rank0 is the
#   0-based rank, so base_i + rank0 is 1-based.
struct HostNewtonMap
  n::Int
  slack::Int
  nvar::Int
  nnzY::Int
  nnzJ::Int
  yrowptr::Vector{Int32}
  ycol::Vector{Int32}
  yrow::Vector{Int32}
  yval::Vector{ComplexF64}
  jpos::Vector{Int32}
  mcnt::Vector{Int32}
  dpos::Vector{Int32}
  ydiag::Vector{Int32}
  pos::Vector{Int32}
  nonslack::Vector{Int32}
  jrowptr::Vector{Int32}
  jcol::Vector{Int32}
  max_row_len::Int
end

function HostNewtonMap(Ybus::SparseMatrixCSC{ComplexF64}, slack::Int)
  n = size(Ybus, 1)
  size(Ybus, 2) == n || throw(ArgumentError("Ybus is not square"))
  # CSR of Ybus: the CSC of its transpose. transpose, never the adjoint
  # (Ybus' would conjugate the values).
  Yt = copy(transpose(Ybus))
  yrowptr = Yt.colptr
  ycol = Yt.rowval
  yval = Yt.nzval
  nnzY = length(yval)
  yrow = Vector{Int32}(undef, nnzY)
  for i = 1:n, k = yrowptr[i]:(yrowptr[i+1]-1)
    yrow[k] = i
  end
  nonslack = [i for i = 1:n if i != slack]
  pos = zeros(Int32, n)
  for (r, i) in enumerate(nonslack)
    pos[i] = r
  end
  nvar = 2 * (n - 1)
  mcnt = zeros(Int32, n)
  dpos = zeros(Int32, n)
  ydiag = zeros(Int32, n)
  jpos = zeros(Int32, nnzY)
  jrowptr = Vector{Int32}(undef, nvar + 1)
  jcol = Int32[]
  max_row_len = 0
  base = 1
  Ci = Int[]
  for i in nonslack
    lo = yrowptr[i]
    hi = yrowptr[i+1] - 1
    max_row_len = max(max_row_len, hi - lo + 1)
    empty!(Ci)
    for k = lo:hi
      j = ycol[k]
      j == slack && continue
      push!(Ci, j)
      j == i && (ydiag[i] = k)
    end
    ydiag[i] == 0 && push!(Ci, i)
    sort!(Ci)
    mi = length(Ci)
    mcnt[i] = mi
    p = pos[i]
    jrowptr[2p-1] = base
    jrowptr[2p] = base + 2mi
    for _ = 1:2
      for c in Ci
        push!(jcol, pos[c])
      end
      for c in Ci
        push!(jcol, (n - 1) + pos[c])
      end
    end
    for k = lo:hi
      j = ycol[k]
      (j == slack || j == i) && continue
      jpos[k] = base + (searchsortedfirst(Ci, j) - 1)
    end
    dpos[i] = base + (searchsortedfirst(Ci, i) - 1)
    base += 4mi
  end
  jrowptr[nvar+1] = base
  nnzJ = base - 1
  # the device code indexes scenario s at (s-1)*nnzJ + p in Int64, but the
  # map itself is Int32
  nnzJ < typemax(Int32) || throw(ArgumentError("Jacobian too large for Int32 indices"))
  return HostNewtonMap(n, slack, nvar, nnzY, nnzJ, Int32.(yrowptr), Int32.(ycol), yrow, yval, jpos, mcnt, dpos, ydiag, pos, Int32.(nonslack), jrowptr, jcol, max_row_len)
end

# the map in device memory
struct DeviceNewtonMap
  host::HostNewtonMap
  yrowptr::CuVector{Int32}
  ycol::CuVector{Int32}
  yrow::CuVector{Int32}
  yval::CuVector{ComplexF64}
  jpos::CuVector{Int32}
  mcnt::CuVector{Int32}
  dpos::CuVector{Int32}
  ydiag::CuVector{Int32}
  pos::CuVector{Int32}
  nonslack::CuVector{Int32}
  jrowptr::CuVector{Int32}
  jcol::CuVector{Int32}
end

DeviceNewtonMap(h::HostNewtonMap) = DeviceNewtonMap(h, CuVector(h.yrowptr), CuVector(h.ycol), CuVector(h.yrow), CuVector(h.yval), CuVector(h.jpos), CuVector(h.mcnt), CuVector(h.dpos), CuVector(h.ydiag), CuVector(h.pos), CuVector(h.nonslack), CuVector(h.jrowptr), CuVector(h.jcol))

# bus type codes of the device mask
const _BT_SLACK = UInt8(0)
const _BT_PQ = UInt8(1)
const _BT_PV = UInt8(2)

function _bus_type_code(t::Symbol)
  t === :PQ && return _BT_PQ
  t === :PV && return _BT_PV
  t === :Slack && return _BT_SLACK
  throw(ArgumentError("unsupported bus type $(t)"))
end

_bus_type_codes!(dst::AbstractVector{UInt8}, types::AbstractVector{Symbol}) = (dst .= _bus_type_code.(types); dst)
