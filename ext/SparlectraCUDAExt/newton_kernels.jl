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
# file: ext/SparlectraCUDAExt/newton_kernels.jl
# purpose: the kernels of the device Newton loop: current and mismatch,
#          Jacobian values straight into the CSR value vector cuDSS
#          factorizes, the voltage update. Batch kernels: blockIdx().y is
#          the scenario, a single run is N = 1. Each kernel repeats the CPU
#          arithmetic operation by operation and in the CPU's order (see the
#          comments with the CPU source), with the non-contracting helpers
#          of device_arith.jl, so mismatch and Jacobian equal the CPU values
#          bit for bit. The update kernel is the exception: it uses Base
#          complex division and cis (sincos of libdevice) and is not
#          bitwise equal to the CPU.

const _THREADS = 256

# Shared-memory max over the block, NaN-propagating like Julia's max (CUDA.jl
# lowers max(::Float64, ::Float64) to llvm.maximum). Every thread of the
# block must call it. Float32 for phase 1 of the two-phase batch.
@inline function _block_max(val::R) where {R<:Union{Float32,Float64}}
  sh = CuStaticSharedArray(R, _THREADS)
  t = threadIdx().x
  @inbounds sh[t] = val
  sync_threads()
  stride = _THREADS ÷ 2
  while stride > 0
    if t <= stride
      @inbounds sh[t] = max(sh[t], sh[t+stride])
    end
    sync_threads()
    stride ÷= 2
  end
  return @inbounds sh[1]
end

# Mismatch rows of bus i from its current I_i, as mismatch_rectangular:
#   S_calc = V .* conj.(I)           Base *(z, w): re = zr*wr - zi*wi,
#                                                  im = zr*wi + zi*wr
#   PQ: dP = re(S_calc) - re(S), dQ = im(S_calc) - im(S)
#   PV: dP, dV = abs(V_i) - Vset_i   (abs(V) is hypot)
# writes F, rhs = T(-F) and returns max(abs(dP), abs(second row)). The
# kernels are generic in the real type of V (R): Float64 as above, Float32
# in phase 1 of the two-phase batch (zero(R) and one(R) give the Float64
# literals of the CPU code for R = Float64).
@inline function _mismatch_rows!(F, rhs, i, p, s, Ir::R, Ii::R, V, S, Vset, btype) where {R<:Union{Float32,Float64}}
  v = @inbounds V[i, s]
  a = real(v)
  b = imag(v)
  sc_re = _sb(_m(a, Ir), _m(b, -Ii))
  sc_im = _ad(_m(a, -Ii), _m(b, Ir))
  sp = @inbounds S[i, s]
  dP = _sb(sc_re, real(sp))
  second = if @inbounds(btype[i, s]) == _BT_PV
    _sb(_hypot_cr(a, b), @inbounds(Vset[i, s]))
  else
    _sb(sc_im, imag(sp))
  end
  T = eltype(rhs)
  @inbounds F[2p-1, s] = dP
  @inbounds F[2p, s] = second
  @inbounds rhs[2p-1, s] = T(-dP)
  @inbounds rhs[2p, s] = T(-second)
  return max(abs(dP), abs(second))
end

# I = Ybus * V and the mismatch, one thread per bus and scenario.
# The CPU computes I in SparseArrays._spmatmul! (stdlib v1.13, linalg.jl
# lines 168 to 190): C[i] = muladd(Y_ik, V_k, C[i]) over the entries of the
# row in ascending column order, from C = 0. The complex muladd
# (base/complex.jl lines 293 to 297) is
#   re = yr*xr - (yi*xi - cr),  im = yr*xi + (yi*xr + ci)
# and is compiled without fma on the CPU (SLP-vectorized vmulpd/vsubpd/
# vaddpd). benchmark/gpu/check_jacobian_gpu.jl checks on every run that
# Sparlectra's CPU SpMV gives the bits of this unfused association and
# fails if it gives the fused one.
function _k_current_mismatch!(I, F, rhs, partial, yrowptr, ycol, yval, V, S, Vset, btype, pos, active)
  s = blockIdx().y
  @inbounds active[s] == 0x00 && return nothing
  i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
  n = size(V, 1)
  R = real(eltype(V))
  val = zero(R)
  if i <= n
    acc_re = zero(R)
    acc_im = zero(R)
    @inbounds for k = yrowptr[i]:(yrowptr[i+1]-Int32(1))
      y = yval[k]
      g = real(y)
      h = imag(y)
      vj = V[ycol[k], s]
      c = real(vj)
      d = imag(vj)
      acc_re = _sb(_m(g, c), _sb(_m(h, d), acc_re))
      acc_im = _ad(_m(g, d), _ad(_m(h, c), acc_im))
    end
    @inbounds I[i, s] = Complex{R}(acc_re, acc_im)
    p = @inbounds pos[i]
    if p != 0
      val = _mismatch_rows!(F, rhs, i, p, s, acc_re, acc_im, V, S, Vset, btype)
    end
  end
  m = _block_max(val)
  threadIdx().x == 1 && (@inbounds partial[blockIdx().x, s] = m)
  return nothing
end

# The mismatch from the current I of the same iterate, after a Q-limit
# switch changed S and the bus types (V and I are unchanged). Only the
# scenarios flagged in `recompute`.
function _k_mismatch_from_current!(F, rhs, partial, I, V, S, Vset, btype, pos, recompute)
  s = blockIdx().y
  @inbounds recompute[s] == 0x00 && return nothing
  i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
  n = size(V, 1)
  val = zero(real(eltype(V)))
  if i <= n
    p = @inbounds pos[i]
    if p != 0
      cur = @inbounds I[i, s]
      val = _mismatch_rows!(F, rhs, i, p, s, real(cur), imag(cur), V, S, Vset, btype)
    end
  end
  m = _block_max(val)
  threadIdx().x == 1 && (@inbounds partial[blockIdx().x, s] = m)
  return nothing
end

# norm[s] = max over the block maxima of scenario s, one block per scenario;
# written to stat[s]
function _k_reduce_max!(stat, partial, nblocks, active)
  s = blockIdx().x
  @inbounds active[s] == 0x00 && return nothing
  t = threadIdx().x
  val = zero(eltype(stat))
  b = t
  while b <= nblocks
    val = max(val, @inbounds partial[b, s])
    b += _THREADS
  end
  m = _block_max(val)
  t == 1 && (@inbounds stat[s] = m)
  return nothing
end

# Off-diagonal Jacobian entries, one thread per Ybus entry k = (i, j),
# j != i, i and j non-slack, and scenario. The CPU
# (_assemble_rectangular_jacobian!) computes
#   S_dVr = V_i * conj(Y_ij)
#   S_dVi = im * (-V_i * conj(Y_ij))     (-V_i) * conj(Y), then Complex{Bool}
#                                        times it: false*x = copysign(0, x),
#                                        true*x = x (base/bool.jl)
# and writes re/im of both to the P row, and to the Q row for PQ. A PV row
# holds explicit zeros here (fixed pattern).
function _k_jacobian_offdiag!(nz, yrow, jpos, mcnt, yval, V, btype, active, nnzJ)
  s = blockIdx().y
  @inbounds active[s] == 0x00 && return nothing
  k = (blockIdx().x - 1) * blockDim().x + threadIdx().x
  k > length(jpos) && return nothing
  p = @inbounds jpos[k]
  p == 0 && return nothing
  i = @inbounds yrow[k]
  mi = Int64(@inbounds mcnt[i])
  v = @inbounds V[i, s]
  a = real(v)
  b = imag(v)
  y = @inbounds yval[k]
  g = real(y)
  h = imag(y)
  p_re = _sb(_m(a, g), _m(b, -h))
  p_im = _ad(_m(a, -h), _m(b, g))
  q_re = _sb(_m(-a, g), _m(-b, -h))
  q_im = _ad(_m(-a, -h), _m(-b, g))
  dvi_re = _sb(copysign(zero(q_re), q_re), q_im)
  dvi_im = _ad(copysign(zero(q_im), q_im), q_re)
  T = eltype(nz)
  o = (Int64(s) - 1) * Int64(nnzJ) + Int64(p)
  @inbounds nz[o] = T(p_re)
  @inbounds nz[o+mi] = T(dvi_re)
  if @inbounds(btype[i, s]) == _BT_PV
    @inbounds nz[o+2mi] = zero(T)
    @inbounds nz[o+3mi] = zero(T)
  else
    @inbounds nz[o+2mi] = T(p_im)
    @inbounds nz[o+3mi] = T(dvi_im)
  end
  return nothing
end

# The 2 x 2 diagonal block of a non-slack bus, one thread per non-slack bus
# and scenario. With a stored Y_ii the CPU adds to the off-diagonal formula
#   S_dVr += conj(I_i),   S_dVi += im * conj(I_i)
# (without one it adds neither). Then the chain terms of voltage-dependent
# injections, emitted for every bus in structural mode and summed with the
# S part by `sparse` (two summands, so the order does not matter):
#   vm_safe = vm > vm_eps ? vm : vm_eps,  dvm = (a, b) / vm_safe
#   P row: + (-dPinj) * dvm;  PQ Q row: + (-dQinj) * dvm
# A PV Q row is d|V|/dVr, d|V|/dVi = (a, b) / vm_div with
# vm_div = vm > 0 ? vm : vm_eps (another threshold than vm_safe, as on the
# CPU).
function _k_jacobian_diag!(nz, I, V, yval, ydiag, dpos, mcnt, dPinj, dQinj, btype, nonslack, active, vm_eps, nnzJ)
  s = blockIdx().y
  @inbounds active[s] == 0x00 && return nothing
  t = (blockIdx().x - 1) * blockDim().x + threadIdx().x
  t > length(nonslack) && return nothing
  i = @inbounds nonslack[t]
  mi = Int64(@inbounds mcnt[i])
  v = @inbounds V[i, s]
  a = real(v)
  b = imag(v)
  vm = _hypot_cr(a, b)
  vm_safe = vm > vm_eps ? vm : vm_eps
  dvm_dvr = a / vm_safe
  dvm_dvi = b / vm_safe
  dP = @inbounds dPinj[i, s]
  c_pr = _m(-dP, dvm_dvr)
  c_pi = _m(-dP, dvm_dvi)
  pv = @inbounds(btype[i, s]) == _BT_PV
  k = @inbounds ydiag[i]
  T = eltype(nz)
  o = (Int64(s) - 1) * Int64(nnzJ) + Int64(@inbounds dpos[i])
  if k != 0
    cur = @inbounds I[i, s]
    ir = real(cur)
    ii = imag(cur)
    y = @inbounds yval[k]
    g = real(y)
    h = imag(y)
    p_re = _sb(_m(a, g), _m(b, -h))
    p_im = _ad(_m(a, -h), _m(b, g))
    q_re = _sb(_m(-a, g), _m(-b, -h))
    q_im = _ad(_m(-a, -h), _m(-b, g))
    dvi_re = _sb(copysign(zero(q_re), q_re), q_im)
    dvi_im = _ad(copysign(zero(q_im), q_im), q_re)
    dvr_re = _ad(p_re, ir)
    dvr_im = _ad(p_im, -ii)
    t_re = _sb(copysign(zero(ir), ir), -ii)
    t_im = _ad(copysign(zero(ii), -ii), ir)
    dvi_re = _ad(dvi_re, t_re)
    dvi_im = _ad(dvi_im, t_im)
    @inbounds nz[o] = T(_ad(dvr_re, c_pr))
    @inbounds nz[o+mi] = T(_ad(dvi_re, c_pi))
    if !pv
      dQ = @inbounds dQinj[i, s]
      @inbounds nz[o+2mi] = T(_ad(dvr_im, _m(-dQ, dvm_dvr)))
      @inbounds nz[o+3mi] = T(_ad(dvi_im, _m(-dQ, dvm_dvi)))
    end
  else
    @inbounds nz[o] = T(c_pr)
    @inbounds nz[o+mi] = T(c_pi)
    if !pv
      dQ = @inbounds dQinj[i, s]
      @inbounds nz[o+2mi] = T(_m(-dQ, dvm_dvr))
      @inbounds nz[o+3mi] = T(_m(-dQ, dvm_dvi))
    end
  end
  if pv
    vm_div = vm > zero(vm) ? vm : vm_eps
    @inbounds nz[o+2mi] = T(a / vm_div)
    @inbounds nz[o+3mi] = T(b / vm_div)
  end
  return nothing
end

# abs(V) as the CPU computes it, for the proof
function _k_vm!(vm, V)
  s = blockIdx().y
  i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
  i > size(V, 1) && return nothing
  v = @inbounds V[i, s]
  @inbounds vm[i, s] = _hypot_cr(real(v), imag(v))
  return nothing
end

# The voltage update of _apply_rectangular_delta!, one thread per non-slack
# bus and scenario. A non-finite solution entry gives a non-finite V entry
# and so a non-finite mismatch norm in the next iteration, which is how the
# loop detects a failed step. The iterate before the step is kept in Vlast
# by the loop.
#   rectangular: V += alpha * dx
#   polar:       q = dV / V_b, V_b * (1 + real(q)) * cis(imag(q)), the linear
#                step for V_b == 0
# Products and sums use the non-contracting helpers as everywhere else; the
# complex division dV / V_b (Base code) and cis (sincos of libdevice) do
# not give the CPU's bits, so this kernel is exempt from the static FMA
# check and the proof measures its difference to the CPU update.
function _k_update!(V, sol, alpha, nonslack, active, polar::Bool)
  s = blockIdx().y
  @inbounds active[s] == 0x00 && return nothing
  t = (blockIdx().x - 1) * blockDim().x + threadIdx().x
  nn = length(nonslack)
  t > nn && return nothing
  i = @inbounds nonslack[t]
  R = real(eltype(V))
  x = R(@inbounds sol[t, s])
  y = R(@inbounds sol[nn+t, s])
  al = @inbounds alpha[s]
  vb = @inbounds V[i, s]
  dv = Complex{R}(_m(al, x), _m(al, y))
  if polar && vb != zero(vb)
    q = dv / vb
    w = cis(imag(q))
    r = _ad(one(R), real(q))
    # Complex * Real, then Complex * Complex as in Base
    zr = _m(r, real(vb))
    zi = _m(r, imag(vb))
    @inbounds V[i, s] = Complex{R}(_sb(_m(zr, real(w)), _m(zi, imag(w))), _ad(_m(zr, imag(w)), _m(zi, real(w))))
  else
    @inbounds V[i, s] = Complex{R}(_ad(real(vb), real(dv)), _ad(imag(vb), imag(dv)))
  end
  return nothing
end

# What the CPU Q-limit decision reads of the iterate, single scenario:
# out[1] = the mismatch norm stat[1], out[1 + k] for request k = req[k]:
#   req[k] = b > 0: imag(S_calc_b) = imag(V_b * conj(I_b)), the operations
#                   of the mismatch rows (bitwise the CPU calc_injections)
#   req[k] = -b:    abs(V_b) as the CPU computes it
# so one download per iteration carries the norm and the Q-limit values.
function _k_qlimit_sample!(out, req, nreq, I, V, stat)
  k = (blockIdx().x - 1) * blockDim().x + threadIdx().x
  k == 1 && (@inbounds out[1] = stat[1])
  k > nreq && return nothing
  r = Int(@inbounds req[k])
  if r > 0
    v = @inbounds V[r, 1]
    cur = @inbounds I[r, 1]
    @inbounds out[1+k] = _ad(_m(real(v), -imag(cur)), _m(imag(v), real(cur)))
  else
    v = @inbounds V[-r, 1]
    @inbounds out[1+k] = _hypot_cr(real(v), imag(v))
  end
  return nothing
end

# A Q-limit switch on the device, single scenario: buf holds four values
# per switched bus (bus, type code, real S, imag S), written into btype and
# S.
function _k_switch_scatter!(S, btype, buf, m)
  k = (blockIdx().x - 1) * blockDim().x + threadIdx().x
  k > m && return nothing
  o = 4 * (k - 1)
  bus = unsafe_trunc(Int, @inbounds buf[o+1])
  @inbounds btype[bus, 1] = unsafe_trunc(UInt8, @inbounds buf[o+2])
  @inbounds S[bus, 1] = ComplexF64(@inbounds(buf[o+3]), @inbounds(buf[o+4]))
  return nothing
end

_blocks(n::Integer) = cld(max(Int(n), 1), _THREADS)
