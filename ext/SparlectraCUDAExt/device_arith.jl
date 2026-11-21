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
# file: ext/SparlectraCUDAExt/device_arith.jl
# purpose: arithmetic for the device Newton kernels that gives the CPU's
#          bits. GPUCompiler hands the NVPTX back end fma_contraction = true
#          (GPUCompiler 2.22.2, src/ptx.jl lines 337 to 345), so a plain
#          a*b + c in a kernel becomes fma.rn.f64, which the CPU code does
#          not do. The .rn intrinsics of CUDA.jl (CUDACore 6.5.1,
#          src/device/intrinsics/math.jl lines 573 to 615) are never fused.

# Every Float64 (or Float32) +, - and * of the mismatch and Jacobian kernels
# goes through these three. Division and sqrt stay plain: they lower to
# div.rn and sqrt.rn (IEEE, correctly rounded) unless the afn flag is set,
# which nothing here does. Never use muladd in device code: CUDACore maps it
# to fma on the device.
@inline _m(x::T, y::T) where {T<:Union{Float32,Float64}} = CUDA.mul_rn(x, y)
@inline _ad(x::T, y::T) where {T<:Union{Float32,Float64}} = CUDA.add_rn(x, y)
@inline _sb(x::T, y::T) where {T<:Union{Float32,Float64}} = CUDA.sub_rn(x, y)

# Constants of Base._hypot, evaluated on the host once (sqrt is correctly
# rounded on both sides, so the values are the CPU's).
const _HYPOT_WIDE = sqrt(eps(Float64) / 2)
const _HYPOT_SCALE = eps(Float64) * sqrt(floatmin(Float64))
const _HYPOT_BIG = sqrt(floatmax(Float64) / 2)
const _HYPOT_SMALL = sqrt(floatmin(Float64))
const _HYPOT_INV_SCALE = inv(_HYPOT_SCALE)

# abs(complex(x, y)) as the CPU computes it: a copy of Base._hypot(x, y)
# (base/math.jl, Julia 1.13.3), the branch for hardware fma, which an x86-64
# CPU with FMA takes (the proof script prints
# Core.Intrinsics.have_fma(Float64) of the CPU it runs on). CUDA.jl
# overrides Base.hypot with __nv_hypot, which is not correctly rounded, so
# abs and hypot are never called in a kernel. The muladd of Base is written
# as fma_rn: on the CPU a real muladd compiles to vfmadd. The fma branch is
# correctly rounded, so the result does not depend on that choice. Every
# other operation is the non-contracting form.
@inline function _hypot_cr(x::Float64, y::Float64)
  ax = abs(x)
  ay = abs(y)
  (isinf(ax) || isinf(ay)) && return Inf
  if ay > ax
    ax, ay = ay, ax
  end
  # widely varying operands, also ay == 0
  ay <= _m(ax, _HYPOT_WIDE) && return ax
  scale = _HYPOT_SCALE
  if ax > _HYPOT_BIG
    ax = _m(ax, scale)
    ay = _m(ay, scale)
    scale = inv(scale)
  elseif ay < _HYPOT_SMALL
    # Base divides by scale, a power of two (2^-563): the same as the exact
    # product with 2^563. LLVM rewrites the division (also div_rn) into a
    # plain mul.f64, which the static check of the proof script does not
    # accept, so the product is written out with mul_rn.
    ax = _m(ax, _HYPOT_INV_SCALE)
    ay = _m(ay, _HYPOT_INV_SCALE)
  else
    scale = 1.0
  end
  h = sqrt(CUDA.fma_rn(ax, ax, _m(ay, ay)))
  hsquared = _m(h, h)
  axsquared = _m(ax, ax)
  corr = _sb(_ad(CUDA.fma_rn(-ay, ay, _sb(hsquared, axsquared)), CUDA.fma_rn(h, h, -hsquared)), CUDA.fma_rn(ax, ax, -axsquared))
  h = _sb(h, corr / _m(2.0, h))
  return _m(h, scale)
end

# Float32, for phase 1 of the two-phase batch, which has no CPU result to
# match: the plain formula without the scaling branches (bus voltages in pu
# are far from overflow and underflow), never fused.
@inline _hypot_cr(x::Float32, y::Float32) = sqrt(_ad(_m(x, x), _m(y, y)))
