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

# file: src/stateestimation/tap_changer_data.jl
# purpose: tap changer data of the tap estimation (0.60.0): the effective
#          changer set of a transformer branch (stored by an import, or
#          derived from the typed models or the branch fields), the map from
#          the step states s to the branch ratio and shift, its tangent
#          dt/ds, the impedance factors per step with their slopes, and the
#          start position of a set on the live branch. The per-changer
#          factor and slope live in equicircuit.jl (_tap_changer_factor).
#          A derived set is built per call and never cached on the branch:
#          island estimations run on parallel tasks and must not write
#          shared state.

# A phase tap changer counts as DECLARED nameplate data (not a constructor
# default) when a nameplate row applied it (phase_tap_declared, set by
# applyTapNameplate! for MATPOWER mpc.sparlectra.tap_changers and SCF
# tap_changer rows) or when a caller set its direction psi in
# tap_est_alpha_deg by hand. The row decides, not psi: a degree-form changer
# of an SCF file with alpha_deg 0 is declared (0.60.0; before, psi 0 kept it
# out). Only declared phase changers join the mass release and the
# generator's deviation targeting. Moved here from the SE service (0.60.0),
# which imports it.
_declared_phase_tap(br::Branch)::Bool = br.has_phase_tap && (br.phase_step_deg > 0.0 || br.phase_du_step > 0.0) && (br.phase_tap_declared || br.tap_est_alpha_deg != 0.0)

# angle difference in degrees folded into [-180, 180]
_tap_wrap_deg(d::Float64)::Float64 = d - 360.0 * round(d / 360.0)

# The tap changers of branch `br` for the tap estimation, `nothing` for a
# line. In this order: the set an import stored on the branch; the set
# derived from the typed models of a resolved winding (origin :typed); the
# set derived from the branch fields (origin :fields). An empty set means the
# branch has no tap changer and is never released.
function _effective_tap_changers(br::Branch)::Union{Nothing,SETapChangerSet}
  br.ratio == 0.0 && return nothing
  br.se_tap_changers === nothing || return br.se_tap_changers
  (br.taps_derived && br.tap_winding !== nothing) && return _typed_tap_changers(br, br.tap_winding)
  return _field_tap_changers(br)
end

# Set of a winding with typed models, mirroring resolve_winding_taps: the
# neutral point (w.ratio, w.shift_degree), then the ratio model, which
# multiplies (side :from), then the phase model. A ratio model without a
# step increment moves nothing and is left out (the old estimator also
# ignored it, _tap_r1_model).
function _typed_tap_changers(br::Branch, w::PowerTransformerWinding)::SETapChangerSet
  chs = SETapChanger[]
  t = w.taps
  if t !== nothing && t.tapStepPercent != 0.0
    lo, hi = minmax(t.lowStep, t.highStep)
    push!(chs, SETapChanger(name = "ratio", kind = :ratio, role = :ratio, side = :from, du = t.tapStepPercent / 100.0, step_assumed = br.ratio_step_assumed, neutral = t.neutralStep, low = lo, high = hi, step = t.step))
  end
  m = w.phase_taps
  m === nothing || push!(chs, _typed_phase_changer(br, m))
  return SETapChangerSet(ratio_n = something(w.ratio, 1.0), shift_n_deg = something(w.shift_degree, 0.0), changers = chs, origin = :typed)
end

# One typed phase model as a changer. The reciprocal convention divides the
# branch ratio (side :to), the direct one multiplies; a table is taken
# verbatim for both conventions (calcPhaseTapAngleRatio), so it is :from.
# An asymmetrical model along the real axis (psi 0 or 180, the DTF in-phase
# changer) is a ratio changer with the signed increment u cos(psi). The
# reactance curve keeps the resolver's form: tan2 for asymmetrical, sin2 for
# symmetrical, a_max at highStep; `rx_correction` follows the branch's
# model.tap_changer_model, which the resolver applies to the asymmetrical
# kind only.
function _typed_phase_changer(br::Branch, m::PhaseTapChangerModel)::SETapChanger
  if m.kind === :tabular
    allx = all(p -> p.x_pu !== nothing, m.table)
    pts = SETapPoint[SETapPoint(step = p.step, ratio = p.ratio, angle_deg = p.angle_deg, x = allx ? p.x_pu : 1.0) for p in m.table]
    return SETapChanger(name = "phase", kind = :tabular, role = :phase, side = :from, step_assumed = br.phase_step_assumed, neutral = m.neutralStep, low = m.lowStep, high = m.highStep, step = m.step, table = pts, impedance = allx ? :table : :none)
  end
  side = m.convention === :reciprocal_from_side ? :to : :from
  u = something(m.voltage_step_increment, 0.0)
  xc = nothing
  if m.x_min !== nothing && m.x_max !== nothing
    amax = calcPhaseTapAngleRatio(m; step = m.highStep).effective_shift_deg
    # without an angle range the resolver's curve divides by zero
    amax != 0.0 && (xc = (m.x_min, m.x_max, amax))
  end
  imp = xc === nothing ? :none : :x_curve
  if m.kind === :symmetrical
    return SETapChanger(name = "phase", kind = :symmetrical, role = :phase, side = side, du = u, step_assumed = br.phase_step_assumed, neutral = m.neutralStep, low = m.lowStep, high = m.highStep, step = m.step, x_curve = xc, x_curve_form = :sin2, impedance = imp)
  end
  rx = br.tap_correction === :impedance_correction
  psi = something(m.winding_connection_angle_deg, 0.0)
  if abs(sind(psi)) <= 1e-9
    return SETapChanger(name = "ratio", kind = :ratio, role = :ratio, side = side, du = u * cosd(psi), step_assumed = br.phase_step_assumed, neutral = m.neutralStep, low = m.lowStep, high = m.highStep, step = m.step, x_curve = xc, x_curve_form = :tan2, impedance = imp, rx_correction = rx)
  end
  return SETapChanger(name = "phase", kind = :asymmetrical, role = :phase, side = side, du = u, alpha_deg = psi, step_assumed = br.phase_step_assumed, neutral = m.neutralStep, low = m.lowStep, high = m.highStep, step = m.step, x_curve = xc, x_curve_form = :tan2, impedance = imp, rx_correction = rx)
end

# Set from the branch fields, neutral point (ratio, angle) as
# applyTapNameplate! keeps it. Only declared changers join (0.60.0): the
# ratio changer when a nameplate row declared it (ratio_tap_declared), the
# phase changer when declared (_declared_phase_tap). The constructor gives
# every transformer has_ratio_tap and a default grid, which is no changer
# model: a transformer that only carries a ratio (plain MATPOWER without
# mpc.sparlectra.tap_changers, a DTF transformer without changer data, a
# hand-built branch without applyTapNameplate!) gets an empty set and is
# never released. The Delta-u form divides (side :to), the degree form adds
# degrees (side :from, kind :linear). A
# branch TapStepImpedance gives the changers of the roles it has values for
# the impedance source :tsi.
function _field_tap_changers(br::Branch)::SETapChangerSet
  chs = SETapChanger[]
  tsi = br.tap_step_impedance
  imp_ratio = (tsi !== nothing && !isempty(tsi.ratio_values)) ? :tsi : :none
  imp_phase = (tsi !== nothing && !isempty(tsi.phase_values)) ? :tsi : :none
  declared = _declared_phase_tap(br)
  du_changer = declared && br.phase_du_step > 0.0
  psi = br.tap_est_alpha_deg
  # with a Delta-u changer tap_ratio carries both factors; the ratio
  # changer's own position comes from the split
  ratio_live, r2 = du_changer ? _tap_du_positions(br, psi) : (br.tap_ratio, 0.0)
  (br.has_ratio_tap && br.ratio_tap_declared) && push!(chs, _field_ratio_changer(br, ratio_live, imp_ratio))
  if du_changer
    du = br.phase_du_step
    lo, hi = minmax(Int(round(br.phase_du_min_step)), Int(round(br.phase_du_max_step)))
    push!(chs, SETapChanger(name = "phase", kind = :asymmetrical, role = :phase, side = :to, du = du, alpha_deg = psi, step_assumed = br.phase_step_assumed, neutral = 0, low = lo, high = hi, step = r2 / du, impedance = imp_phase))
  elseif declared
    dphi = br.phase_step_deg
    # the degree band is stored relative to neutral (applyTapNameplate!)
    lo, hi = minmax(Int(round(br.phase_min_deg / dphi)), Int(round(br.phase_max_deg / dphi)))
    push!(chs, SETapChanger(name = "phase", kind = :linear, role = :phase, side = :from, dphi_deg = dphi, step_assumed = br.phase_step_assumed, neutral = 0, low = lo, high = hi, step = (br.phase_shift_deg - br.angle) / dphi, impedance = imp_phase))
  end
  return SETapChangerSet(ratio_n = br.ratio, shift_n_deg = br.angle, changers = chs, origin = :fields)
end

# Positions (ratio_live, r2) of a declared Delta-u changer from the live
# branch through _du_tap_split. A ratio changer and a Delta-u changer along
# the real axis cannot be split when t has a phase part (_du_tap_split
# throws there); the set then starts both at neutral and the start position
# search (_tap_s_from_branch) reports the mismatch. The test below is the
# one _du_tap_split makes before it throws.
function _tap_du_positions(br::Branch, psi_deg::Float64)::Tuple{Float64,Float64}
  if br.has_ratio_tap && br.tap_step > 0.0 && abs(sin(deg2rad(psi_deg))) < 1e-9
    base = br.ratio
    tlive = (br.tap_ratio == 0.0 ? base : br.tap_ratio) * cis(deg2rad(br.phase_shift_deg))
    v = base * cis(deg2rad(br.angle)) / tlive
    abs(imag(v)) <= 1e-12 * abs(v) || return (base, 0.0)
  end
  return _du_tap_split(br, psi_deg)
end

# The ratio changer of the branch fields on its grid, positions counted from
# `ratio` (neutral 0):
# - :fraction (nameplate default, constructor): factor 1 + s du dividing the
#   ratio (side :to), du = tap_step; the band from tap_min/tap_max as the old
#   estimator read it (0.9/1.1 where they are not set);
# - :linear without points: additive, R = ratio + s tap_step;
# - :linear with points (CGMES RatioTapChangerTable): a table of the
#   positions the tap controller reaches (_tap_table_positions), ratio as
#   factor of `ratio`.
function _field_ratio_changer(br::Branch, ratio_live::Float64, imp::Symbol)::SETapChanger
  base = br.ratio
  if br.tap_ratio_grid === :linear
    step_abs = _tap_step1(br)
    if br.tap_ratio_points !== nothing
      pos = _tap_table_positions(br)
      if length(pos) > 1
        pts = SETapPoint[SETapPoint(step = n, ratio = r / base) for (n, r) in pos]
        return SETapChanger(name = "ratio", kind = :tabular, role = :ratio, side = :from, step_assumed = br.ratio_step_assumed, neutral = 0, low = pos[1][1], high = pos[end][1], step = _tap_table_position_of(pts, ratio_live / base), table = pts, impedance = imp)
      end
    end
    lo, hi = minmax(ceil(Int, (br.tap_min - base) / step_abs - 1e-9), floor(Int, (br.tap_max - base) / step_abs + 1e-9))
    return SETapChanger(name = "ratio", kind = :ratio, role = :ratio, side = :from, additive = true, step_abs = step_abs, du = step_abs / base, step_assumed = br.ratio_step_assumed, neutral = 0, low = lo, high = hi, step = (ratio_live - base) / step_abs, impedance = imp)
  end
  du = _tap_step1(br)
  tmin = br.tap_min > 0.0 ? br.tap_min : 0.9
  tmax = br.tap_max > 0.0 ? br.tap_max : 1.1
  # the highest step gives the smallest ratio
  lo, hi = minmax(ceil(Int, (base / tmax - 1.0) / du - 1e-9), floor(Int, (base / tmin - 1.0) / du + 1e-9))
  return SETapChanger(name = "ratio", kind = :ratio, role = :ratio, side = :to, du = du, step_assumed = br.ratio_step_assumed, neutral = 0, low = lo, high = hi, step = (base / ratio_live - 1.0) / du, impedance = imp)
end

# continuous step of the factor `f` on a table with ascending ratios, linear
# between the points, the end intervals continued beyond them
function _tap_table_position_of(pts::Vector{SETapPoint}, f::Float64)::Float64
  k = findfirst(p -> p.ratio >= f, pts)
  k = k === nothing ? length(pts) : max(k, 2)
  a, b = pts[k - 1], pts[k]
  return a.step + (b.step - a.step) * (f - a.ratio) / (b.ratio - a.ratio)
end

# Branch ratio and shift (degrees) of `set` at the step states `s` (one per
# changer, relative to neutral): the import's polar accumulation in import
# order. An additive changer adds its steps to the ratio; every other one
# multiplies (side :from) or divides (side :to) by its factor magnitude and
# adds or subtracts its angle.
function _tap_set_ratio_shift(set::SETapChangerSet, s::AbstractVector{<:Real})::Tuple{Float64,Float64}
  length(s) == length(set.changers) || throw(ArgumentError("_tap_set_ratio_shift: $(length(s)) step states for $(length(set.changers)) changers."))
  R = set.ratio_n
  Phi = set.shift_n_deg
  for (c, ch) in enumerate(set.changers)
    if ch.additive
      R += s[c] * ch.step_abs
    else
      m, a = _tap_changer_factor(ch, s[c])
      if ch.side === :from
        R *= m
        Phi += a
      else
        R /= m
        Phi -= a
      end
    end
  end
  return (R, Phi)
end

# (R, Phi_deg, dR/ds_c, dPhi_deg/ds_c) of `set` at `s`. A multiplicative
# changer with sign e (+1 :from, -1 :to) gives dR = e R dm/m and
# dPhi = e da. The additive changer is the first of its set and the later
# changers multiply onto R_a = R_n + s step_abs, so dR = step_abs R / R_a.
function _tap_set_partials(set::SETapChangerSet, s::AbstractVector{<:Real}, c::Int)::NTuple{4,Float64}
  R, Phi = _tap_set_ratio_shift(set, s)
  ch = set.changers[c]
  ch.additive && return (R, Phi, ch.step_abs * R / (set.ratio_n + s[c] * ch.step_abs), 0.0)
  m, _ = _tap_changer_factor(ch, s[c])
  dm, da = _tap_changer_slope(ch, s[c])
  e = ch.side === :from ? 1.0 : -1.0
  return (R, Phi, e * R * dm / m, e * da)
end

# dt/ds_c of the branch tap t = R e^(j Phi) at `s`:
# e^(j Phi) (dR + j R dPhi), dPhi in radians
function _tap_set_tangent(set::SETapChangerSet, s::AbstractVector{<:Real}, c::Int)::ComplexF64
  R, Phi, dR, dPhi = _tap_set_partials(set, s, c)
  return cis(deg2rad(Phi)) * (dR + im * R * deg2rad(dPhi))
end

# whether changer `ch` moves the impedance with its step: per-step data
# (`impedance` source) or the analytic correction; which of them applies is
# decided per transformer by _tap_reactance_source
_tap_changer_has_impedance(ch::SETapChanger)::Bool = ch.impedance !== :none || ch.rx_correction

# whether `set` carries per-step impedance data: a table with r and x, the
# reactance curve from xMin/xMax, or TapStepImpedance factors
_tap_set_has_table(set::SETapChangerSet)::Bool = any(ch -> ch.impedance !== :none, set.changers)

# Source of the reactance of a released transformer under the impedance
# rule `mode` (state_estimation.tap_estimation.impedance_correction), one
# rule per transformer (0.60.0):
# per-step data always follow the step, whatever `mode` says; without them
# `auto` takes the analytic correction (rx_correction,
# |1 + f e^(j psi)|^2) where a changer has it and `off` keeps the start
# impedance. Returns `:table`, `:analytic` or `:nominal` (the impedance
# does not move). Where the analytic factor applies next to per-step data
# says _tap_changer_impedance.
function _tap_reactance_source(set::SETapChangerSet, mode::Symbol)::Symbol
  _tap_set_has_table(set) && return :table
  (mode === :auto && any(ch -> ch.rx_correction, set.changers)) && return :analytic
  return :nominal
end

# The set as a run under `mode` uses it for the impedance. Under `off` the
# analytic correction is not computed, so a typed set (origin :typed)
# loses its rx_correction flag; its x curve or x table stays and follows
# the step. Only typed sets need this: a stored set with per-step data
# drops the flag in _tap_changer_impedance under both values, and one
# without them is :nominal under `off`. A typed set is never written back
# as a stored set (_tap_apply_position!), so the flag is not lost on the
# branch.
function _tap_set_under(set::SETapChangerSet, mode::Symbol)::SETapChangerSet
  (mode === :off && set.origin === :typed && any(ch -> ch.rx_correction, set.changers)) || return set
  chs = SETapChanger[ch.rx_correction ? SETapChanger(ch.name, ch.mrid, ch.kind, ch.role, ch.side, ch.additive, ch.du, ch.alpha_deg, ch.dphi_deg, ch.step_abs, ch.neutral, ch.low, ch.high, ch.step, ch.table, ch.outside_du, ch.x_curve, ch.x_curve_form, ch.impedance, false, ch.step_assumed) : ch for ch in set.changers]
  return SETapChangerSet(set.ratio_n, set.shift_n_deg, chs, set.origin)
end

# Impedance factors (r, x, g, b) of `set` at `s`: the product of the
# factors of every changer with per-step data. The estimator scales the
# start impedance by F(s) / F(s0). `tsi` is the branch TapStepImpedance for
# changers with source :tsi, whose pairs were recorded with the other
# changers at their start position: the actuator value of a :tsi changer is
# taken with the others at `s0`.
function _tap_impedance_factors(set::SETapChangerSet, s::AbstractVector{<:Real}; s0::AbstractVector{<:Real} = s, tsi::Union{Nothing,TapStepImpedance} = nothing)::NTuple{4,Float64}
  F = (1.0, 1.0, 1.0, 1.0)
  for c in eachindex(set.changers)
    _tap_changer_has_impedance(set.changers[c]) || continue
    F = F .* first(_tap_changer_impedance(set, c, s, s0, tsi))
  end
  return F
end

# dF/ds_c of _tap_impedance_factors at `s`; the factors of the other
# changers do not depend on s_c (a :tsi changer reads its actuator with the
# others at s0)
function _tap_impedance_slope(set::SETapChangerSet, s::AbstractVector{<:Real}, c::Int; s0::AbstractVector{<:Real} = s, tsi::Union{Nothing,TapStepImpedance} = nothing)::NTuple{4,Float64}
  _tap_changer_has_impedance(set.changers[c]) || return (0.0, 0.0, 0.0, 0.0)
  dF = last(_tap_changer_impedance(set, c, s, s0, tsi))
  for k in eachindex(set.changers)
    (k == c || !_tap_changer_has_impedance(set.changers[k])) && continue
    dF = dF .* first(_tap_changer_impedance(set, k, s, s0, tsi))
  end
  return dF
end

# (factors, d factors / ds_c) of changer c. The source gives r, x, g, b.
# The analytic correction (rx_correction) scales r always and x unless a
# reactance curve gives x, as the power flow resolver does for a typed
# model (resolve_winding_taps: x_alpha replaces the corrected x). On a set
# from an import (CGMES, PowSyBl, a stored set) with per-step data the
# flag is dropped: the table, the reactance curve or the TapStepImpedance
# factors give the impedance alone.
# A typed set keeps the resolver's form, so the estimation under `auto`
# equals the power flow; under `off` _tap_set_under has removed the flag.
function _tap_changer_impedance(set::SETapChangerSet, c::Int, s::AbstractVector{<:Real}, s0::AbstractVector{<:Real}, tsi::Union{Nothing,TapStepImpedance})::Tuple{NTuple{4,Float64},NTuple{4,Float64}}
  ch = set.changers[c]
  sc = Float64(s[c])
  F = (1.0, 1.0, 1.0, 1.0)
  dF = (0.0, 0.0, 0.0, 0.0)
  src = ch.impedance
  if src === :table
    F, dF = _tap_table_impedance(ch, sc)
  elseif src === :x_curve
    fx, dfx = _tap_x_curve(ch, sc)
    F = (1.0, fx, 1.0, 1.0)
    dF = (0.0, dfx, 0.0, 0.0)
  elseif src === :tsi
    F, dF = _tap_tsi_impedance(set, c, s, s0, tsi)
  end
  if ch.rx_correction && (set.origin === :typed || !_tap_set_has_table(set))
    q, dq = _tap_rx_factor(ch, sc)
    fr, dfr = F[1] * q, dF[1] * q + F[1] * dq
    fx, dfx = src === :x_curve ? (F[2], dF[2]) : (F[2] * q, dF[2] * q + F[2] * dq)
    F = (fr, fx, F[3], F[4])
    dF = (dfr, dfx, dF[3], dF[4])
  end
  return (F, dF)
end

# |1 + s du e^(j alpha)|^2 and its slope 2 Re(conj(v) v'), the factor of
# calcTapImpedanceCorrectionFactor (alpha 0 for an in-phase changer)
function _tap_rx_factor(ch::SETapChanger, s::Float64)::Tuple{Float64,Float64}
  ch.du === nothing && return (1.0, 0.0)
  alpha = ch.kind === :asymmetrical ? ch.alpha_deg : 0.0
  dv = ch.du * cis(deg2rad(alpha))
  v = 1.0 + s * dv
  return (abs2(v), 2.0 * real(conj(v) * dv))
end

# Reactance curve (x_min, x_max, a_max) over the changer's own angle a(s):
# x = x_min + (x_max - x_min) g^2 with g = sin(a/2)/sin(a_max/2) (:sin2,
# CGMES import) or tan(a)/tan(a_max) (:tan2, typed asymmetrical); returns x
# and dx/ds. At a whole step x / x_min is the import's factor.
function _tap_x_curve(ch::SETapChanger, s::Float64)::Tuple{Float64,Float64}
  xmin, xmax, amax = ch.x_curve
  a, da = ch.additive ? (0.0, 0.0) : (last(_tap_changer_factor(ch, s)), last(_tap_changer_slope(ch, s)))
  if ch.x_curve_form === :sin2
    g = sind(a / 2.0) / sind(amax / 2.0)
    dg = (pi / 360.0) * cosd(a / 2.0) / sind(amax / 2.0)
  else
    g = tand(a) / tand(amax)
    dg = (pi / 180.0) / (cosd(a)^2 * tand(amax))
  end
  return (xmin + (xmax - xmin) * g^2, (xmax - xmin) * 2.0 * g * dg * da)
end

# r, x, g, b of a table changer at s, linear between the points with the
# slope of the current interval. Outside a ratio table (outside_du set) the
# nearest end point, slope 0, like the import's _tapTablePointAtStep; a
# phase table continues its end interval, each factor floored at half its
# smallest table value.
function _tap_table_impedance(ch::SETapChanger, s::Float64)::Tuple{NTuple{4,Float64},NTuple{4,Float64}}
  pts = ch.table
  n = length(pts)
  zero4 = (0.0, 0.0, 0.0, 0.0)
  vals(p) = (p.r, p.x, p.g, p.b)
  n == 1 && return (vals(pts[1]), zero4)
  p = ch.neutral + s
  if ch.outside_du !== nothing
    p >= pts[end].step && return (vals(pts[end]), zero4)
    p < pts[1].step && return (vals(pts[1]), zero4)
  end
  k = _tap_table_interval(pts, p)
  a, b = vals(pts[k]), vals(pts[k + 1])
  h = Float64(pts[k + 1].step - pts[k].step)
  F = ntuple(i -> a[i] + (b[i] - a[i]) * (p - pts[k].step) / h, 4)
  dF = ntuple(i -> (b[i] - a[i]) / h, 4)
  floors = ntuple(i -> 0.5 * minimum(q -> vals(q)[i], pts), 4)
  return (ntuple(i -> max(F[i], floors[i]), 4), ntuple(i -> F[i] < floors[i] ? 0.0 : dF[i], 4))
end

# TapStepImpedance factors of changer c: its actuator value (branch ratio
# for :ratio, shift in degrees for :phase) with the other changers at s0,
# piecewise linear over the recorded pairs in step order, the end interval
# continued with each factor floored at half its smallest value; the slope
# is dF/dv dv/ds_c.
function _tap_tsi_impedance(set::SETapChangerSet, c::Int, s::AbstractVector{<:Real}, s0::AbstractVector{<:Real}, tsi::Union{Nothing,TapStepImpedance})::Tuple{NTuple{4,Float64},NTuple{4,Float64}}
  ones4 = (1.0, 1.0, 1.0, 1.0)
  zero4 = (0.0, 0.0, 0.0, 0.0)
  tsi === nothing && return (ones4, zero4)
  ch = set.changers[c]
  vals, facs = ch.role === :ratio ? (tsi.ratio_values, tsi.ratio_factors) : (tsi.phase_values, tsi.phase_factors)
  isempty(vals) && return (ones4, zero4)
  length(vals) == 1 && return (facs[1], zero4)
  sx = Float64[Float64(x) for x in s0]
  sx[c] = Float64(s[c])
  R, Phi, dR, dPhi = _tap_set_partials(set, sx, c)
  v, dv = ch.role === :ratio ? (R, dR) : (Phi, dPhi)
  n = length(vals)
  k = findfirst(j -> min(vals[j], vals[j + 1]) <= v <= max(vals[j], vals[j + 1]), 1:n - 1)
  k === nothing && (k = abs(v - vals[1]) <= abs(v - vals[n]) ? 1 : n - 1)
  h = vals[k + 1] - vals[k]
  h == 0.0 && return (facs[k], zero4)
  w = (v - vals[k]) / h
  F = ntuple(i -> facs[k][i] + (facs[k + 1][i] - facs[k][i]) * w, 4)
  dF = ntuple(i -> (facs[k + 1][i] - facs[k][i]) / h * dv, 4)
  floors = ntuple(i -> 0.5 * minimum(f -> f[i], facs), 4)
  return (ntuple(i -> max(F[i], floors[i]), 4), ntuple(i -> F[i] < floors[i] ? 0.0 : dF[i], 4))
end

# band of a changer in steps relative to neutral (lower end first): the
# changer's low..high, for a table without the ratio formula outside (a
# phase table) its table range
function _tap_band(ch::SETapChanger)::Tuple{Float64,Float64}
  (ch.kind === :tabular && ch.outside_du === nothing) && return (Float64(ch.table[1].step - ch.neutral), Float64(ch.table[end].step - ch.neutral))
  return (Float64(ch.low - ch.neutral), Float64(ch.high - ch.neutral))
end

# Size of one step of a changer in the pu of the branch tap (sigma; the
# estimator scales its tap columns and the step test by it): |du| for the
# formula kinds (an additive changer without du: step_abs / ratio_n), the
# angle per step in radians for :linear, for a table the mean over its
# intervals of |d ratio| + |d angle| (radians) per step; at least 1e-6.
function _tap_scale(ch::SETapChanger; ratio_n::Float64 = 1.0)::Float64
  k = ch.kind
  if k in (:ratio, :asymmetrical, :symmetrical)
    sc = ch.du === nothing ? abs(ch.step_abs) / ratio_n : abs(ch.du)
  elseif k === :linear
    sc = deg2rad(abs(ch.dphi_deg))
  else
    pts = ch.table
    if length(pts) > 1
      sc = sum((abs(pts[j + 1].ratio - pts[j].ratio) + deg2rad(abs(pts[j + 1].angle_deg - pts[j].angle_deg))) / (pts[j + 1].step - pts[j].step) for j in 1:length(pts) - 1) / (length(pts) - 1)
    else
      sc = ch.outside_du === nothing ? 0.0 : abs(ch.outside_du)
    end
  end
  return max(sc, 1e-6)
end

# release mode of a transformer from the roles of its effective set:
# :both, :ratio, :pst, or nothing (line, or no tap changer)
function _tap_release_mode(br::Branch)::Union{Nothing,Symbol}
  set = _effective_tap_changers(br)
  (set === nothing || isempty(set.changers)) && return nothing
  hr = any(ch -> ch.role === :ratio, set.changers)
  hp = any(ch -> ch.role === :phase, set.changers)
  return hr && hp ? :both : (hr ? :ratio : :pst)
end

# True when `net` has at least one closed transformer that is no machine
# transformer and carries a DECLARED tap changer (the Web UI checks
# "estimate taps" by default for such a case), that is a transformer the
# SE service releases (_tap_release_mode): a non-empty set an import
# stored, a typed winding with a changer, or a declared changer of the
# branch fields (_field_tap_changers).
function _case_has_tap_changers(net::Net)::Bool
  for (k, br) in enumerate(net.branchVec)
    br.ratio == 0.0 && continue
    _branch_terminal_state(br) == :closed || continue
    _tap_release_mode(br) === nothing && continue
    _is_machine_transformer(net, k) && continue
    return true
  end
  return false
end

# index of the changer a documented deviation moves: the first ratio
# changer, else the first phase changer, 0 for an empty set. Shared by the
# measurement generator and the auto release so both act on the same one.
function _tap_deviation_changer(set::SETapChangerSet)::Int
  k = findfirst(ch -> ch.role === :ratio, set.changers)
  k === nothing && (k = findfirst(ch -> ch.role === :phase, set.changers))
  return k === nothing ? 0 : k
end

# Which changers of `set` raise the rank of the real 2 x k matrix of
# (dt/ds_c)/t = (dR/R, dPhi) at `s`, walked in set order (true = separable
# from the ones before). The ratio plane has two real directions, so a
# third changer never raises it; two changers along the same direction (two
# ratio changers, ratio plus asymmetrical at alpha 0 or 180, two pure angle
# changers) do not either.
function _tap_separable_mask(set::SETapChangerSet, s::AbstractVector{<:Real})::Vector{Bool}
  out = falses(length(set.changers))
  basis = nothing
  for c in eachindex(set.changers)
    R, _, dR, dPhi = _tap_set_partials(set, s, c)
    col = (dR / R, deg2rad(dPhi))
    ncol = hypot(col[1], col[2])
    if basis === nothing
      ncol > 1e-12 || continue
      basis = col
      out[c] = true
    else
      cross = basis[1] * col[2] - basis[2] * col[1]
      if abs(cross) > 1e-9 * hypot(basis[1], basis[2]) * ncol
        out[c] = true
        break
      end
    end
  end
  return collect(out)
end

# Start position of `set` on the live branch: the stored steps when they
# reproduce tap_ratio and phase_shift_deg to 1e-10, else a damped
# Gauss-Newton search on the steps of every changer (at most 50 iterations,
# each step at most one band width), else `nothing` (the estimator freezes
# the branch with :model_mismatch). The search runs only when every changer
# is separable; inseparable steps (two ratio changers on one branch) cannot
# be read from one complex number. Needed because tap controllers, a
# write-back or user edits move the branch after the import.
function _tap_s_from_branch(set::SETapChangerSet, br::Branch)::Union{Nothing,Vector{Float64}}
  s = Float64[ch.step - ch.neutral for ch in set.changers]
  function resid(x)
    R, Phi = _tap_set_ratio_shift(set, x)
    return [R / br.tap_ratio - 1.0, deg2rad(_tap_wrap_deg(Phi - br.phase_shift_deg))]
  end
  r = resid(s)
  maximum(abs, r) <= 1e-10 && return s
  k = length(s)
  (k == 0 || !all(_tap_separable_mask(set, s))) && return nothing
  span = Float64[max(ch.high - ch.low, 1) for ch in set.changers]
  for _ in 1:50
    J = zeros(2, k)
    for c in 1:k
      _, _, dR, dPhi = _tap_set_partials(set, s, c)
      J[1, c] = dR / br.tap_ratio
      J[2, c] = deg2rad(dPhi)
    end
    G = J' * J
    gmax = maximum(G[i, i] for i in 1:k)
    gmax > 0.0 || return nothing
    dx = -((G + 1e-9 * gmax * I) \ (J' * r))
    s .+= clamp.(dx, -span, span)
    r = resid(s)
    maximum(abs, r) <= 1e-10 && return s
  end
  return nothing
end
