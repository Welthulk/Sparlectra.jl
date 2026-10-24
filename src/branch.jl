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

# Author: Udo Schmitz (https://github.com/Welthulk)
# Date: 10.05.2023
# file: src/branch.jl
# purpose: branch model types (Branch, BranchModel, BranchFlow) and
#          admittance, ratio, flow, and loss helpers for network branches

# helper
"""
    BranchFlow

A structure representing the flow in a branch of a power system.

# Fields
- `vm_pu::Union{Nothing,Float64}`: The voltage magnitude in per unit.
- `va_deg::Union{Nothing,Float64}`: The voltage angle in degrees.
- `pFlow::Union{Nothing,Float64}`: The active power flow.
- `qFlow::Union{Nothing,Float64}`: The reactive power flow.

# Constructors
- `BranchFlow(vm_pu::Union{Nothing,Float64} = nothing, va_deg::Union{Nothing,Float64} = nothing, pFlow::Union{Nothing,Float64} = nothing, qFlow::Union{Nothing,Float64} = nothing)`: Creates a new `BranchFlow` instance.

# Methods
- `Base.show(io::IO, b::BranchFlow)`: Prints the `BranchFlow` instance.

# Example
```julia
BranchFlow(vm_pu = 1.0, va_deg = 0.0, pFlow = 100.0, qFlow = 50.0)
```
"""
struct BranchFlow
  vm_pu::Union{Nothing,Float64} # voltage magnitude
  va_deg::Union{Nothing,Float64} # voltage angle
  pFlow::Union{Nothing,Float64} # active power flow
  qFlow::Union{Nothing,Float64} # reactive power flow

  function BranchFlow(vm_pu::Union{Nothing,Float64} = nothing, va_deg::Union{Nothing,Float64} = nothing, pFlow::Union{Nothing,Float64} = nothing, qFlow::Union{Nothing,Float64} = nothing)
    new(vm_pu, va_deg, pFlow, qFlow)
  end

  function Base.show(io::IO, b::BranchFlow)
    print(io, "BranchFlow( ")
    print(io, "vm: ", b.vm_pu, ", ")
    print(io, "va: ", b.va_deg, ", ")
    print(io, "pFlow: ", b.pFlow, ", ")
    print(io, "qFlow: ", b.qFlow, ", ")
    println(io, ")")
  end
end
"""
    BranchModel

A structure representing a branch model in a power system.

# Fields
- `r_pu::Float64`: The per unit resistance of the branch.
- `x_pu::Float64`: The per unit reactance of the branch.
- `b_pu::Float64`: The per unit total line charging susceptance of the branch.
- `g_pu::Float64`: The per unit total line charging conductance of the branch.
- `ratio::Float64`: The transformer off nominal turns ratio.
- `angle::Float64`: The transformer off nominal phase shift angle.
- `sn_MVA::Union{Nothing,Float64}`: The nominal power of the branch = rateA.
- `g_from_pu`, `b_from_pu`, `g_to_pu`, `b_to_pu::Float64`: the charging
  admittance per terminal; the totals are always their sums, the
  symmetric half split is the default.

# Constructors
- `BranchModel(; r_pu, x_pu, b_pu, g_pu, ratio, angle, sn_MVA = nothing, g_from_pu = nothing, b_from_pu = nothing, g_to_pu = nothing, b_to_pu = nothing)`:
  the four terminal values are given together or not at all; without them
  the totals are split in halves.

# Example
```julia
BranchModel(r_pu = 0.01, x_pu = 0.1, b_pu = 0.02, g_pu = 0.02, ratio = 1.0, angle = 0.0, sn_MVA = 100.0)
```
"""
struct BranchModel
  r_pu::Float64
  x_pu::Float64
  b_pu::Float64
  g_pu::Float64
  ratio::Float64
  angle::Float64
  sn_MVA::Union{Nothing,Float64}
  # per-terminal charging admittance (0.20.0); the totals above are their sums
  g_from_pu::Float64
  b_from_pu::Float64
  g_to_pu::Float64
  b_to_pu::Float64
end

function BranchModel(; r_pu::Float64, x_pu::Float64, b_pu::Float64, g_pu::Float64, ratio::Float64, angle::Float64, sn_MVA::Union{Nothing,Float64} = nothing, g_from_pu = nothing, b_from_pu = nothing, g_to_pu = nothing, b_to_pu = nothing)
  split = _branch_shunt_split_kw("BranchModel", g_from_pu, b_from_pu, g_to_pu, b_to_pu)
  gf, bf, gt, bt = _split_branch_shunt(g_pu, b_pu, split)
  return BranchModel(r_pu, x_pu, bf + bt, gf + gt, ratio, angle, sn_MVA, gf, bf, gt, bt)
end

# The four terminal values of a branch shunt from the totals and an optional
# explicit split; `nothing` is the symmetric half at each end (the MATPOWER
# model and every builder that passes totals only).
function _split_branch_shunt(g_pu, b_pu, split::Union{Nothing,NTuple{4,Float64}})::NTuple{4,Float64}
  split === nothing || return split
  g = Float64(something(g_pu, 0.0))
  b = Float64(something(b_pu, 0.0))
  return (0.5 * g, 0.5 * b, 0.5 * g, 0.5 * b)
end

# Keyword validation shared by every branch builder: all four terminal
# values or none. A partial set is a caller error naming the branch, so a
# forgotten to-side value cannot silently become a symmetric half.
function _branch_shunt_split_kw(name::AbstractString, g_from_pu, b_from_pu, g_to_pu, b_to_pu)::Union{Nothing,NTuple{4,Float64}}
  given = (g_from_pu !== nothing, b_from_pu !== nothing, g_to_pu !== nothing, b_to_pu !== nothing)
  all(given) && return (Float64(g_from_pu), Float64(b_from_pu), Float64(g_to_pu), Float64(b_to_pu))
  any(given) && throw(ArgumentError("$(name): g_from_pu, b_from_pu, g_to_pu and b_to_pu must be given together (all four) or not at all"))
  return nothing
end

"""
    TapStepImpedance

Impedance of a transformer branch per tap step, for tap changer tables
whose points carry r, x, g, b deviations (CGMES `TapChangerTablePoint`).
`r_pu`, `x_pu` and the four terminal admittances are the branch values
without the factors of the tables listed here. `ratio_values` holds the
branch `tap_ratio` at each step of the ratio table and `ratio_factors` the
factors `(r, x, g, b)` of that step; `phase_values` and `phase_factors` do
the same for `phase_shift_deg`. An empty pair means no table on that
actuator. A tap controller that moves the branch calls
`_apply_tap_step_impedance!`, which picks the step nearest to the live
`tap_ratio` and `phase_shift_deg`.
"""
struct TapStepImpedance
  r_pu::Float64
  x_pu::Float64
  g_from_pu::Float64
  b_from_pu::Float64
  g_to_pu::Float64
  b_to_pu::Float64
  ratio_values::Vector{Float64}
  ratio_factors::Vector{NTuple{4,Float64}}
  phase_values::Vector{Float64}
  phase_factors::Vector{NTuple{4,Float64}}
end

"""
    SETapPoint

One point of a tap changer table as the state estimation reads it: the
absolute `step`, the changer's own factor at that step (`ratio`,
`angle_deg`), and the impedance at the step as factors `r`, `x`, `g`, `b`
(1 + percent/100 of a CGMES or IIDM point; the absolute `x` in pu for a
typed table with `x_pu`; 1.0 without data). Only ratios of the impedance
numbers enter the estimator.
"""
struct SETapPoint
  step::Int
  ratio::Float64
  angle_deg::Float64
  r::Float64
  x::Float64
  g::Float64
  b::Float64
end

SETapPoint(; step::Integer, ratio::Real, angle_deg::Real = 0.0, r::Real = 1.0, x::Real = 1.0, g::Real = 1.0, b::Real = 1.0) = SETapPoint(Int(step), Float64(ratio), Float64(angle_deg), Float64(r), Float64(x), Float64(g), Float64(b))

"""
    SETapChanger

Data of one tap changer of a transformer branch for the tap estimation
(state estimation only; power flow, controllers and export do not read it).
`kind` is the changer type: `:ratio` (in-phase, factor `1 + s du`),
`:asymmetrical` (factor `1 + s du e^(j alpha)`, a quadrature changer is
alpha 90), `:symmetrical` (angle `2 atan(s du / 2)`, magnitude 1),
`:tabular` (table points, linear in between) or `:linear` (angle
`s dphi`, magnitude 1). `role` is `:ratio` or `:phase` and decides which
release moves it. `side` is `:from` when the factor multiplies the branch
ratio and `:to` when it divides it; `additive` marks the ratio changer of a
linear nameplate grid, `R = R_n + s step_abs`, which must be the first
changer of its set. Positions are absolute steps (`neutral`, `low`, `high`,
the import position `step`); the estimator's state is `s = step - neutral`.
`outside_du` continues a ratio table with `1 + s outside_du` beyond its
points. `x_curve = (x_min, x_max, a_max_deg)` is a reactance curve over the
changer's own angle in the form `x_curve_form` (`:sin2` of the CGMES
import, `:tan2` of a typed asymmetrical model). `impedance` names the
per-step impedance data (`:none`, `:table`, `:x_curve`, `:tsi` for the
branch `TapStepImpedance`), `rx_correction` the impedance correction of
`model.tap_changer_model = impedance_correction`. `step_assumed` marks a
changer whose import carried no step: `step` is then the middle position
(the neutral step, the band middle without one) and the report row of the
tap estimation says so. Optional numbers are `nothing`. The keyword
constructor checks the fields and throws an `ArgumentError` naming the
field.
"""
struct SETapChanger
  name::String
  mrid::String
  kind::Symbol
  role::Symbol
  side::Symbol
  additive::Bool                              # :ratio only: R = R_n + s * step_abs (linear nameplate grid)
  du::Union{Nothing,Float64}                  # increment per step as a fraction of the winding voltage, signed
  alpha_deg::Union{Nothing,Float64}           # direction of the additional voltage (:asymmetrical), from the import, never a state
  dphi_deg::Union{Nothing,Float64}            # angle per step (:linear)
  step_abs::Union{Nothing,Float64}            # ratio increment per step (additive :ratio)
  neutral::Int
  low::Int
  high::Int
  step::Float64                               # import position, absolute (a CGMES SvTapStep may be fractional)
  table::Union{Nothing,Vector{SETapPoint}}    # :tabular: points, strictly ascending in step
  outside_du::Union{Nothing,Float64}          # :tabular from a ratio table: 1 + s * outside_du beyond the points
  x_curve::Union{Nothing,NTuple{3,Float64}}   # (x_min, x_max, a_max_deg) of the reactance curve
  x_curve_form::Symbol                        # :sin2 (CGMES import) or :tan2 (typed asymmetrical model)
  impedance::Symbol                           # :none, :table, :x_curve, :tsi
  rx_correction::Bool                         # typed model under model.tap_changer_model = :impedance_correction
  step_assumed::Bool                          # the import carried no step, `step` is the middle position
end

const SE_TAP_CHANGER_KINDS = [:ratio, :asymmetrical, :symmetrical, :tabular, :linear]
const SE_TAP_IMPEDANCE_SOURCES = [:none, :table, :x_curve, :tsi]

_se_tap_opt(v) = v === nothing ? nothing : Float64(v)

function SETapChanger(;
  name::AbstractString = "",
  mrid::AbstractString = "",
  kind::Symbol,
  role::Symbol,
  side::Symbol,
  additive::Bool = false,
  du::Union{Nothing,Real} = nothing,
  alpha_deg::Union{Nothing,Real} = nothing,
  dphi_deg::Union{Nothing,Real} = nothing,
  step_abs::Union{Nothing,Real} = nothing,
  neutral::Integer,
  low::Integer,
  high::Integer,
  step::Real,
  table::Union{Nothing,Vector{SETapPoint}} = nothing,
  outside_du::Union{Nothing,Real} = nothing,
  x_curve::Union{Nothing,NTuple{3,Real}} = nothing,
  x_curve_form::Symbol = :sin2,
  impedance::Symbol = :none,
  rx_correction::Bool = false,
  step_assumed::Bool = false,
)
  what = isempty(name) ? "SETapChanger" : "SETapChanger $(name)"
  kind in SE_TAP_CHANGER_KINDS || throw(ArgumentError("$(what): kind must be one of $(join(SE_TAP_CHANGER_KINDS, ", ")), got $(repr(kind))."))
  role in (:ratio, :phase) || throw(ArgumentError("$(what): role must be :ratio or :phase, got $(repr(role))."))
  side in (:from, :to) || throw(ArgumentError("$(what): side must be :from or :to, got $(repr(side))."))
  kind === :ratio && role !== :ratio && throw(ArgumentError("$(what): role must be :ratio for kind :ratio."))
  kind in (:asymmetrical, :symmetrical, :linear) && role !== :phase && throw(ArgumentError("$(what): role must be :phase for kind $(repr(kind))."))
  low <= high || throw(ArgumentError("$(what): low ($(low)) must not exceed high ($(high))."))
  isfinite(step) || throw(ArgumentError("$(what): step must be finite, got $(step)."))
  additive && kind !== :ratio && throw(ArgumentError("$(what): additive is allowed only for kind :ratio."))
  additive && step_abs === nothing && throw(ArgumentError("$(what): step_abs is required for an additive changer."))
  kind === :ratio && !additive && du === nothing && throw(ArgumentError("$(what): du is required for kind :ratio."))
  kind === :asymmetrical && du === nothing && throw(ArgumentError("$(what): du is required for kind :asymmetrical."))
  kind === :asymmetrical && alpha_deg === nothing && throw(ArgumentError("$(what): alpha_deg is required for kind :asymmetrical."))
  kind === :symmetrical && du === nothing && throw(ArgumentError("$(what): du is required for kind :symmetrical."))
  kind === :linear && dphi_deg === nothing && throw(ArgumentError("$(what): dphi_deg is required for kind :linear."))
  if kind === :tabular
    (table === nothing || isempty(table)) && throw(ArgumentError("$(what): table is required and must not be empty for kind :tabular."))
    all(k -> table[k].step < table[k + 1].step, 1:length(table) - 1) || throw(ArgumentError("$(what): table steps must be strictly ascending."))
    all(p -> p.ratio > 0.0, table) || throw(ArgumentError("$(what): table ratios must be positive."))
  end
  outside_du !== nothing && kind !== :tabular && throw(ArgumentError("$(what): outside_du is allowed only for kind :tabular."))
  x_curve_form in (:sin2, :tan2) || throw(ArgumentError("$(what): x_curve_form must be :sin2 or :tan2, got $(repr(x_curve_form))."))
  impedance in SE_TAP_IMPEDANCE_SOURCES || throw(ArgumentError("$(what): impedance must be one of $(join(SE_TAP_IMPEDANCE_SOURCES, ", ")), got $(repr(impedance))."))
  impedance === :table && table === nothing && throw(ArgumentError("$(what): impedance :table needs a table."))
  impedance === :x_curve && x_curve === nothing && throw(ArgumentError("$(what): impedance :x_curve needs x_curve."))
  xc = x_curve === nothing ? nothing : (Float64(x_curve[1]), Float64(x_curve[2]), Float64(x_curve[3]))
  return SETapChanger(String(name), String(mrid), kind, role, side, additive, _se_tap_opt(du), _se_tap_opt(alpha_deg), _se_tap_opt(dphi_deg), _se_tap_opt(step_abs), Int(neutral), Int(low), Int(high), Float64(step), table, _se_tap_opt(outside_du), xc, x_curve_form, impedance, rx_correction, step_assumed)
end

"""
    SETapChangerSet

The tap changers of one transformer branch for the tap estimation: the
branch ratio `ratio_n` and shift `shift_n_deg` without the changers'
factors (the branch at the neutral steps of the formula kinds), the
changers in the composition order of the import, and the `origin`
(`:cgmes`, `:powsybl` for a set an import stored on the branch; `:typed`,
`:fields` for a set derived on demand from the typed models or the branch
fields). An empty `changers` vector stored by an import declares "no tap
changer": the transformer is never released. The keyword constructor
checks the fields (`ArgumentError` naming the field), among them that an
additive changer comes first.
"""
struct SETapChangerSet
  ratio_n::Float64
  shift_n_deg::Float64
  changers::Vector{SETapChanger}
  origin::Symbol
end

function SETapChangerSet(; ratio_n::Real, shift_n_deg::Real = 0.0, changers::Vector{SETapChanger} = SETapChanger[], origin::Symbol)
  (isfinite(ratio_n) && ratio_n > 0.0) || throw(ArgumentError("SETapChangerSet: ratio_n must be positive and finite, got $(ratio_n)."))
  isfinite(shift_n_deg) || throw(ArgumentError("SETapChangerSet: shift_n_deg must be finite, got $(shift_n_deg)."))
  origin in (:cgmes, :powsybl, :typed, :fields) || throw(ArgumentError("SETapChangerSet: origin must be :cgmes, :powsybl, :typed or :fields, got $(repr(origin))."))
  for (k, ch) in enumerate(changers)
    ch.additive && k != 1 && throw(ArgumentError("SETapChangerSet: changers[$(k)] ($(ch.name)) is additive and must be the first changer of its set."))
  end
  return SETapChangerSet(Float64(ratio_n), Float64(shift_n_deg), changers, origin)
end

"""
    Branch

A mutable structure representing a branch in a power system.

# Fields
- `comp::AbstractComponent`: The component of the branch.
- `fromBus::Integer`: The index of the bus where the branch starts.
- `toBus::Integer`: The index of the bus where the branch ends.
- `r_pu::Float64`: The per unit resistance of the branch.
- `x_pu::Float64`: The per unit reactance of the branch.
- `b_pu::Float64`: The per unit total line charging susceptance of the branch (sum of the two terminal values).
- `g_pu::Float64`: The per unit total line charging conductance of the branch (sum of the two terminal values).
- `g_from_pu`, `b_from_pu`, `g_to_pu`, `b_to_pu::Float64`: the charging admittance per terminal; the from arm sits behind the ideal transformer. Change them through `set_branch_shunt!`, never one total alone.
- `ratio::Float64`: The transformer off nominal turns ratio.
- `angle::Float64`: The transformer off nominal phase shift angle.
- `status::Integer`: The status of the branch. 1 = in service, 0 = out of service.
- `sn_MVA::Union{Nothing,Float64}`: The nominal power of the branch = rateA.
- `fBranchFlow::Union{Nothing,BranchFlow}`: The flow from fromNodeID to toNodeID.
- `tBranchFlow::Union{Nothing,BranchFlow}`: The flow from toNodeID to fromNodeID.
- `pLosses::Union{Nothing,Float64}`: The active power losses.
- `qLosses::Union{Nothing,Float64}`: The reactive power losses.
- `tap_step_impedance::Union{Nothing,TapStepImpedance}`: impedance per tap step of a tap changer table, applied when a tap controller moves the branch.
- `se_tap_changers::Union{Nothing,SETapChangerSet}`: the tap changers of a transformer as an import stored them for the tap estimation (state estimation only); `nothing` means the set is derived from the typed models or the branch fields.
- `phase_tap_declared::Bool`: a nameplate row declared the phase tap changer of the branch fields (`applyTapNameplate!`); the constructor default phase grid is not declared.
- `ratio_tap_declared::Bool`: a nameplate row declared the ratio tap changer of the branch fields (`applyTapNameplate!` with `tap_step > 0`); the constructor default ratio grid is not declared.
- `ratio_step_assumed::Bool`, `phase_step_assumed::Bool`: the import carried no current step for the ratio changer (phase changer) of the branch fields or the typed models, so it sits at its middle position (the neutral step, the band middle without one); the tap estimation names it in its report row.

# Constructors
- `Branch(; from::Int, to::Int, baseMVA::Float64, branch::AbstractBranch, id::Int, status::Integer = 1, ratio::Union{Nothing,Float64} = nothing, side::Union{Nothing,Int} = nothing, vn_kV::Union{Nothing,Float64} = nothing,
                    fromOid::Union{Nothing,Int} = nothing, toOid::Union{Nothing,Int} = nothing)`: Creates a new `Branch` instance.

# Methods
- `Base.show(io::IO, b::Branch)`: Prints the `Branch` instance.
"""
mutable struct Branch <: AbstractBranch
  comp::AbstractComponent
  branchIdx::Int
  fromBus::Integer
  toBus::Integer
  r_pu::Float64                          # resistance (live; a series-FACTS control run stamps its operating point here)
  x_pu::Float64                          # reactance (live; a series-FACTS control run stamps its operating point here)
  # physical equipment impedance (#329): equals r_pu/x_pu until a series-FACTS
  # controller (TCSC/SSSC, full UPFC) writes a compensated operating point onto
  # the live fields during a power-flow control run. Short circuit and the
  # CGMES/MATPOWER exports read these base values so they see the equipment, not
  # the operating point; the power flow keeps reading the live r_pu/x_pu.
  r_base_pu::Float64
  x_base_pu::Float64
  b_pu::Float64                          # total line charging susceptance (b_from_pu + b_to_pu)
  g_pu::Float64                          # total line charging conductance (g_from_pu + g_to_pu)
  # per-terminal charging admittance (0.20.0): the from arm sits behind the
  # ideal transformer and is stamped through |t|^2, the to arm on the to
  # bus. The symmetric half split is the default and the MATPOWER case;
  # PowSyBl and CGMES place a transformer's magnetizing admittance on one
  # end. set_branch_shunt! is the one writer that keeps the totals in step.
  g_from_pu::Float64
  b_from_pu::Float64
  g_to_pu::Float64
  b_to_pu::Float64
  ratio::Float64                         # nominal turns ratio
  angle::Float64                         # nominal phase shift angle in degrees
  status::Integer                        # 1 = in service, 0 = out of service
  sn_MVA::Union{Nothing,Float64}         # nominal power of the branch = rateA
  fBranchFlow::Union{Nothing,BranchFlow} # flow from fromNodeID to toNodeID
  tBranchFlow::Union{Nothing,BranchFlow} # flow from toNodeID to fromNodeID
  pLosses::Union{Nothing,Float64}        # active power losses
  qLosses::Union{Nothing,Float64}        # reactive power losses
  tap_ratio::Float64
  phase_shift_deg::Float64
  has_ratio_tap::Bool
  has_phase_tap::Bool
  tap_min::Float64
  tap_max::Float64
  tap_step::Float64
  phase_min_deg::Float64
  phase_max_deg::Float64
  phase_step_deg::Float64
  # per-terminal service state (r0.9.10): 1 = closed, 0 = open. The
  # aggregate `status` stays the user-facing switch and is kept consistent
  # by the setters (status = 1 iff both terminals closed); consumers must
  # read the state through _branch_terminal_state, never the raw fields.
  from_status::Int
  to_status::Int
  # open-end voltage of a one-sided open branch (result, not a bus
  # voltage): filled by calcNetLosses! from the pi-model voltage divider,
  # nothing while closed or fully open
  open_end_vm_pu::Union{Nothing,Float64}
  open_end_va_deg::Union{Nothing,Float64}
  # tap estimation release (0.10.0): :none | :ratio | :pst | :both, which
  # changers of the branch's changer set become step states. The direction
  # psi of a nameplate Delta-u changer is nameplate data (never estimated).
  # Set through setTapEstimation!, honored by the SE overlay.
  tap_est_mode::Symbol
  tap_est_alpha_deg::Float64
  # additional-voltage phase stepper (Delta-u PST, 0.10.0): the mechanical
  # grid is the ADDITIONAL-VOLTAGE amplitude per step (n * phase_du_step
  # along psi = tap_est_alpha_deg, the factor 1 + n du e^(j psi) that divides
  # the branch ratio); the shift angle follows as atan and is NOT the grid.
  # 0 = degree-grid mode (the classic phase_step_deg fields apply). Set by
  # the mpc.sparlectra.tap_changers import; the tap estimation reads it as
  # an asymmetrical changer and fixes it on whole Delta-u steps.
  phase_du_step::Float64
  phase_du_min_step::Float64
  phase_du_max_step::Float64
  # typed tap model (0.20.0): `taps_derived` is true when ratio, angle,
  # impedance correction and the tap grid were derived from the models of
  # `tap_winding` by resolve_branch_taps!, so a controller moves the model
  # step and the resolver rewrites the fields; false keeps the legacy
  # degree/ratio grid. `tap_correction` is the model.tap_changer_model
  # option the resolver applies to r_pu/x_pu from the equipment base.
  taps_derived::Bool
  tap_winding::Union{Nothing,PowerTransformerWinding}
  tap_correction::Symbol
  # discrete ratios of a tabular ratio tap changer on the legacy ratio grid
  # (CGMES RatioTapChangerTable), sorted ascending: a discrete ratio
  # controller moves from point to point instead of by `tap_step`, which is
  # then only the mean step. `nothing` keeps the uniform grid.
  tap_ratio_points::Union{Nothing,Vector{Float64}}
  # per-step impedance of tap changer tables with r/x/g/b deviations,
  # nothing for every branch without such a table (see TapStepImpedance)
  tap_step_impedance::Union{Nothing,TapStepImpedance}
  # how a ratio tap position maps onto the ratio, counted in steps from
  # `ratio`: :fraction is the grid ratio / (1 + n * tap_step) of nameplate
  # data (applyTapNameplate!) and the constructor default; :linear is
  # ratio + n * tap_step (or the table points), the grid of a CGMES
  # RatioTapChanger on the tapped from side, which the tap controller also
  # steps on. The tap estimation derives its ratio changer from this grid.
  tap_ratio_grid::Symbol
  # tap changers for the tap estimation (0.60.0), stored by an import that
  # knows them (CGMES, PowSyBl, SCF `estimation` key); an empty set means
  # "no tap changer". nothing: derived on demand from the typed models or
  # the fields above (_effective_tap_changers). Power flow, controllers and
  # export do not read it.
  se_tap_changers::Union{Nothing,SETapChangerSet}
  # true when a nameplate row (applyTapNameplate!: MATPOWER
  # mpc.sparlectra.tap_changers, SCF tap_changer controller 2) declared the
  # phase changer of the fields above; the constructor's default phase grid
  # is not declared and never released for the tap estimation (0.60.0). A
  # degree-form changer is declared by its row, whatever its psi.
  phase_tap_declared::Bool
  # the same for the ratio changer of the fields: true when a nameplate row
  # with tap_step > 0 applied it. The constructor sets has_ratio_tap and a
  # default grid for every transformer, so has_ratio_tap alone cannot tell
  # a declared changer from that grid. Only a declared ratio changer joins
  # the changer set of the fields (_field_tap_changers), so a transformer
  # without one is never released for the tap estimation (0.60.0).
  ratio_tap_declared::Bool
  # true when the import carried no current step for the ratio changer
  # (phase changer) of the fields above or of the typed models: it sits at
  # its middle position (applyTapNameplate! with a NaN step, SCF
  # tap_changer_models without step, a DTF control without Stufe). Read by
  # the tap estimation for its report row (0.60.0); a step written back off
  # that position clears it.
  ratio_step_assumed::Bool
  phase_step_assumed::Bool

  function Branch(;
    branchIdx::Int,
    from::Int,
    to::Int,
    baseMVA::Float64,
    branch::AbstractBranch,
    id::Int,
    status::Integer = 1,
    ratio::Union{Nothing,Float64} = nothing,
    side::Union{Nothing,Int} = nothing,
    vn_kV::Union{Nothing,Float64} = nothing,
    fromOid::Union{Nothing,Int} = nothing,
    toOid::Union{Nothing,Int} = nothing,
    angle::Union{Nothing,Float64} = nothing,
    values_are_pu::Bool=false,
    from_status::Union{Nothing,Integer} = nothing,
    to_status::Union{Nothing,Integer} = nothing,
    shunt_split::Union{Nothing,NTuple{4,Float64}} = nothing,
  )
    # terminal flags default to the aggregate; the stored aggregate is then
    # recomputed so that status = 1 iff both terminals are closed
    fs = Int(something(from_status, status))
    ts = Int(something(to_status, status))
    fs in (0, 1) && ts in (0, 1) || error("Branch: from_status/to_status must be 0 or 1")
    status = (fs == 1 && ts == 1) ? 1 : 0
    if isa(branch, ACLineSegment) # Line
      @assert !isnothing(vn_kV) "vn_kV must be set for an ACLineSegment"
      if isnothing(ratio)
        # to distinguish line from transformer
        ratio  = 0.0
      end
      if !isnothing(fromOid) && !isnothing(toOid)
        c = getBranchComp(vn_kV, fromOid, toOid, id, "ACL")
      else
        c = getBranchComp(vn_kV, from, to, id, "ACL")
      end
      # the line knows the unit of its values (ohm/siemens or per unit), so
      # values_are_pu is not consulted here; it only matters for transformers
      r_pu, x_pu, b_pu, g_pu = getLineRXBG_pu(branch, vn_kV, baseMVA)
      # (g_from, b_from, g_to, b_to): explicit split or the symmetric half;
      # the stored totals are the sums either way
      gf, bf, gt, bt = _split_branch_shunt(g_pu, b_pu, shunt_split)
      new(c, branchIdx, from, to, r_pu, x_pu, r_pu, x_pu, bf + bt, gf + gt, gf, bf, gt, bt, 0.0, 0.0, status, branch.ratedS, nothing, nothing, nothing, nothing, 1.0, 0.0, false, false, 0.9, 1.1, 0.00625, -30.0, 30.0, 1.25, fs, ts, nothing, nothing, :none, 0.0, 0.0, 0.0, 0.0, false, nothing, :ideal, nothing, nothing, :fraction, nothing, false, false, false, false)
    elseif isa(branch, PowerTransformer) # Transformer     
      if (isnothing(side) && branch.isBiWinder)
        side = getSideNumber2WT(branch)
      elseif (isnothing(side) && !branch.isBiWinder)
        error("side must be set for a PowerTransformer")
      end

      c = if !isnothing(fromOid) && !isnothing(toOid)
        getBranchComp(vn_kV, fromOid, toOid, id, "2WT")
      else
        getBranchComp(vn_kV, from, to, id, "2WT")
      end

      w = (side in [1, 2, 3]) ? (side == 1 ? branch.side1 : (side == 2 ? branch.side2 : branch.side3)) : error("wrong value for 'side'")
      vn_kV = isnothing(vn_kV) ? w.Vn : vn_kV
      sn_MVA = getWindingRatedS(w)
      
      if values_are_pu
        r_pu, x_pu, b_pu, g_pu = getTrafoRXBG(w)
      else 
        r_pu, x_pu, b_pu, g_pu = getTrafoRXBG_pu(w, vn_kV, baseMVA)        
      end
      

      ratio = isnothing(ratio) ? 1.0 : ratio
      @assert ratio != 0.0 "ratio must not be 0.0 for transformers"
      #angle = isnothing(w.shift_degree) ? 0.0 : w.shift_degree
      angle = isnothing(angle) ? (isnothing(w.shift_degree) ? 0.0 : w.shift_degree) : angle
      tap_min = 0.9
      tap_max = 1.1
      tap_step = 0.00625
      if !isnothing(w.taps)
        tap_min, tap_max, tap_step = calcRatioTapRange(w.taps)
      end

      gf, bf, gt, bt = _split_branch_shunt(g_pu, b_pu, shunt_split)
      new(c, branchIdx, from, to, r_pu, x_pu, r_pu, x_pu, bf + bt, gf + gt, gf, bf, gt, bt, ratio, angle, status, sn_MVA, nothing, nothing, nothing, nothing, ratio, angle, true, true, tap_min, tap_max, tap_step, -30.0, 30.0, 1.25, fs, ts, nothing, nothing, :none, 0.0, 0.0, 0.0, 0.0, false, nothing, :ideal, nothing, nothing, :fraction, nothing, false, false, false, false)
    elseif isa(branch, BranchModel) # PI-Model
      @assert !isnothing(vn_kV) "vn_kV must be set for PI-Model"

      if !isnothing(fromOid) && !isnothing(toOid)
        c = getBranchComp(vn_kV, fromOid, toOid, id, "PI")
      else
        c = getBranchComp(vn_kV, from, to, id, "PI")
      end

      is_tap = branch.ratio != 0.0
      initial_ratio = is_tap ? branch.ratio : 1.0
      initial_angle = is_tap ? branch.angle : 0.0
      # the model carries its own split; an explicit keyword overrides it
      gf, bf, gt, bt = shunt_split === nothing ? (branch.g_from_pu, branch.b_from_pu, branch.g_to_pu, branch.b_to_pu) : shunt_split
      new(c, branchIdx, from, to, branch.r_pu, branch.x_pu, branch.r_pu, branch.x_pu, bf + bt, gf + gt, gf, bf, gt, bt, branch.ratio, branch.angle, status, branch.sn_MVA, nothing, nothing, nothing, nothing, initial_ratio, initial_angle, is_tap, is_tap, 0.9, 1.1, 0.00625, -30.0, 30.0, 1.25, fs, ts, nothing, nothing, :none, 0.0, 0.0, 0.0, 0.0, false, nothing, :ideal, nothing, nothing, :fraction, nothing, false, false, false, false)
    else
      error("Branch type not supported")
    end
  end

  function Base.show(io::IO, b::Branch)
    print(io, "Branch( ")
    print(io, b.comp, ", ")
    print(io, "branchIdx: ", b.branchIdx, ", ")
    print(io, "fromBus: ", b.fromBus, ", ")
    print(io, "toBus: ", b.toBus, ", ")

    print(io, "r_pu: ", b.r_pu, ", ")
    print(io, "x_pu: ", b.x_pu, ", ")
    print(io, "b_pu: ", b.b_pu, ", ")
    print(io, "g_pu: ", b.g_pu, ", ")
    if !has_symmetric_shunt(b)
      print(io, "shunt from/to: (", b.g_from_pu, " + ", b.b_from_pu, "im) / (", b.g_to_pu, " + ", b.b_to_pu, "im), ")
    end
    b.taps_derived && print(io, "taps: derived from model, ")
    print(io, "ratio: ", b.ratio, ", ")
    print(io, "angle: ", b.angle, ", ")
    print(io, "status: ", b.status, ", ")
    if _branch_terminal_state(b) != :closed
      print(io, "terminal_state: ", _branch_terminal_state(b), ", ")
    end
    if !isnothing(b.sn_MVA)
      print(io, "sn_MVA: ", b.sn_MVA, ", ")
    end
    if (!isnothing(b.fBranchFlow))
      print(io, "BranchFlow (from): ", b.fBranchFlow, ", ")
    end
    if (!isnothing(b.tBranchFlow))
      print(io, "BranchFlow (to): ", b.tBranchFlow, ", ")
    end
    if (!isnothing(b.pLosses))
      print(io, "pLosses: ", b.pLosses, ", ")
    end
    if (!isnothing(b.qLosses))
    print(io, "qLosses: ", b.qLosses, ", ")
    end
    print(io, "tap_ratio: ", b.tap_ratio, ", ")
    print(io, "phase_shift_deg: ", b.phase_shift_deg, ", ")

    println(io, ")")
  end
end

"""
    calcAdmittance(branch, u_rated, s_rated) -> (yaa, yab, yba, ybb)

The four admittance-matrix entries of the branch on the given base,
including ratio and phase shift.
"""
function calcAdmittance(branch::Branch, u_rated::Float64, s_rated::Float64)::Tuple{ComplexF64,ComplexF64,ComplexF64,ComplexF64}
  # Series Admittance ys
  ys = calcBranchYser(branch)
  # the two shunt arms (0.20.0): the from arm sits BEHIND the ideal
  # transformer and is seen through |t|^2 like the series admittance, the
  # to arm sits on the to bus. With the symmetric split both arms equal
  # 0.5 * (g + jb) and the four entries are the classic MATPOWER stamp:
  #   Y_ff = (ys + y0_from) / |t|^2     Y_ft = -ys / conj(t)
  #   Y_tf = -ys / t                    Y_tt =  ys + y0_to
  y0_from = _branch_y0_from(branch)
  y0_to = _branch_y0_to(branch)
  # calc complex ratio
  t = calcBranchRatio(branch)
  Y_11 = (ys + y0_from) / abs2(t)
  Y_12 = -1.0 * ys / conj(t)
  Y_21 = -1.0 * ys / t
  Y_22 = ys + y0_to
  return (Y_11, Y_12, Y_21, Y_22)
end

"""
    _open_terminal_yin(branch) -> ComplexF64

Exact pi-model input admittance of a one-sided open branch seen from its
closed bus, as the Schur complement of the two-port from `calcAdmittance`:
with the TO end open `Y_in = Y11 - Y12*Y21/Y22`, with the FROM end open
`Y_in = Y22 - Y21*Y12/Y11`. The two-port already carries the complex ratio,
so lines and transformers (off-nominal ratio, phase shift) are covered
uniformly without a case distinction. For `|Y_s| >> |Y_0|` this approaches
the FULL line charging `g + jb` (not half of it).
"""
function _open_terminal_yin(branch::Branch)::ComplexF64
  Y11, Y12, Y21, Y22 = calcAdmittance(branch, 1.0, 1.0)
  st = _branch_terminal_state(branch)
  if st == :open_to
    abs(Y22) < 1e-12 && return 0.0 + 0.0im
    return Y11 - Y12 * Y21 / Y22
  elseif st == :open_from
    abs(Y11) < 1e-12 && return 0.0 + 0.0im
    return Y22 - Y21 * Y12 / Y11
  end
  error("_open_terminal_yin: branch $(branch.branchIdx) is $(st), not one-sided open")
end

"""
    _open_end_voltage(branch, u_closed::ComplexF64) -> ComplexF64

Voltage at the open terminal of a one-sided open branch from the pi-model
voltage divider (zero current at the open end): with TO open
`U_open = -Y21/Y22 * U_from`, with FROM open `U_open = -Y12/Y11 * U_to`.
Reproduces the Ferranti rise (`|U_open| > |U_closed|` for b > 0) without
adding a node to the solved system.
"""
function _open_end_voltage(branch::Branch, u_closed::ComplexF64)::ComplexF64
  Y11, Y12, Y21, Y22 = calcAdmittance(branch, 1.0, 1.0)
  st = _branch_terminal_state(branch)
  if st == :open_to
    abs(Y22) < 1e-12 && return u_closed
    return -Y21 / Y22 * u_closed
  elseif st == :open_from
    abs(Y11) < 1e-12 && return u_closed
    return -Y12 / Y11 * u_closed
  end
  error("_open_end_voltage: branch $(branch.branchIdx) is $(st), not one-sided open")
end

# helper
"""
    setBranchFlow!(branch, tfBranchFlow, fBranchFlow)

Store the solved flows of both branch ends.
"""
function setBranchFlow!(branch::Branch, tfBranchFlow::BranchFlow, fBranchFlow::BranchFlow)
  branch.tBranchFlow = tfBranchFlow
  branch.fBranchFlow = fBranchFlow
end

"""
    _branch_terminal_state(br::Branch) -> Symbol

Single source of truth for the service state of a branch: `:closed` (both
terminals in service), `:open_from` / `:open_to` (exactly one terminal
open, the branch reduces to its pi-model Schur complement at the closed
bus), or `:open` (out of service: aggregate `status == 0` with untouched
terminal flags, or both flags open). Every consumer reads the state
through this helper, never the raw fields.
"""
function _branch_terminal_state(br::Branch)::Symbol
  # aggregate forced open without touching the terminal flags (legacy
  # direct writes of br.status): both flags still 1 means fully open
  (br.status == 0 && br.from_status == 1 && br.to_status == 1) && return :open
  br.from_status == 0 && br.to_status == 0 && return :open
  br.from_status == 0 && return :open_from
  br.to_status == 0 && return :open_to
  return :closed
end

"""
    setBranchStatus!(branch, service::Bool)

User-facing aggregate switch: sets the aggregate `status` and BOTH
terminal flags consistently (in service = all closed, out of service =
all open).
"""
function setBranchStatus!(branch::Branch, service::Bool)
  v = service ? 1 : 0
  branch.status = v
  branch.from_status = v
  branch.to_status = v
  return branch
end

"""
    setBranchTerminalStatus!(branch; from = nothing, to = nothing)

Open or close individual branch terminals (`true` = closed, `false` =
open; `nothing` leaves a terminal unchanged) and recompute the aggregate
`status` (1 iff both terminals closed). A branch open at exactly one
terminal stays in the model as its exact pi reduction at the closed bus,
see the "One-sided open branches" section of the branch-model docs.
"""
function setBranchTerminalStatus!(branch::Branch; from::Union{Nothing,Bool} = nothing, to::Union{Nothing,Bool} = nothing)
  from === nothing || (branch.from_status = from ? 1 : 0)
  to === nothing || (branch.to_status = to ? 1 : 0)
  branch.status = (branch.from_status == 1 && branch.to_status == 1) ? 1 : 0
  return branch
end

"""
    getBranchFlow(branch, from, to) -> BranchFlow

The stored flow record for the given orientation of the branch.
"""
function getBranchFlow(branch::Branch, from::Node, to::Node)
  if (branch.fromBus == from.busIdx && branch.toBus == to.busIdx)
    return branch.fBranchFlow
  elseif (branch.fromBus == to.busIdx && branch.toBus == from.busIdx)
    return branch.tBranchFlow
  else
    error("Nodes do not match the branch")
  end
end

"""
    getBranchIdx(branch) -> Int

The position of the branch in the network's branch vector.
"""
function getBranchIdx(branch::Branch)
  return branch.branchIdx
end

"""
    getBranchLosses(branch) -> (p, q)

The stored losses of the branch in MW/MVar.
"""
function getBranchLosses(branch::Branch)
  return branch.pLosses, branch.qLosses
end

"""
    setBranchLosses!(branch, pLosses, qLosses)

Store the solved losses of the branch in MW/MVar.
"""
function setBranchLosses!(branch::Branch, pLosses::Float64, qLosses::Float64)
  branch.pLosses = pLosses
  branch.qLosses = qLosses
end

"""
    calcBranchYser(branch) -> ComplexF64

The series admittance of the branch in per unit.
"""
function calcBranchYser(branch::Branch)::ComplexF64
  return inv((branch.r_pu + branch.x_pu * im))
end

"""
    calcBranchYserBase(branch)

Series admittance from the PHYSICAL equipment impedance (`r_base_pu`/`x_base_pu`,
issue #329). Short circuit and the CGMES/MATPOWER exports use this instead of
`calcBranchYser` so a series-FACTS operating point stamped onto the live
`r_pu`/`x_pu` does not leak into fault or interchange data. The power flow keeps
using `calcBranchYser` (the live, compensated value).
"""
function calcBranchYserBase(branch::Branch)::ComplexF64
  return inv((branch.r_base_pu + branch.x_base_pu * im))
end

# The branchIdx of every branch whose negative series resistance or reactance
# is taken over unchanged by the exports and the short circuit:
#  - a branch built from a CGMES EquivalentBranch (a network equivalent, no
#    equipment; net.equivalent_branch_metadata),
#  - a leg of a three-winding transformer in its star equivalent
#    (_three_winding_groups: an auxiliary star bus whose branches are all
#    transformers). The star split of the pair impedances
#    z1 = (z12 + z13 - z23)/2 and so on gives a negative leg r or x for real
#    transformers; only the pair sums are physical.
# A two-winding transformer or a line keeps the check.
function _negative_r_allowed_branches(net)::Set{Int}
  allowed = Set{Int}(keys(net.equivalent_branch_metadata))
  for g in _three_winding_groups(net), k in g.legs
    push!(allowed, Int(net.branchVec[k].branchIdx))
  end
  return allowed
end

"""
    assertPhysicalBranchImpedances(net, context)

Defensive check (#329): refuse to proceed when any branch carries a NEGATIVE
BASE series resistance. A physical line/transformer has `r >= 0`. Since short
circuit and the CGMES/MATPOWER exports read the base (equipment) impedance
`r_base_pu`/`x_base_pu`, a series-FACTS operating point stamped onto the live
`r_pu`/`x_pu` (e.g. a full UPFC's `Re(z_add) < 0`, #326) no longer reaches
them; this assertion is therefore expected NOT to fire in the normal FACTS plus
export/SC workflow. It stays as a last-resort guard against a corrupted base
model. `context` names the caller (e.g. "CGMES export", "MATPOWER export").
Two kinds of branch are exempt (`_negative_r_allowed_branches`): a branch
built from a CGMES `EquivalentBranch` (`net.equivalent_branch_metadata`), a
network equivalent where a negative r or x is normal in Ward and Kron
reductions, and a leg of a three-winding star equivalent, where a single leg
can be negative and only the sums of two legs are physical.
"""
function assertPhysicalBranchImpedances(net, context::AbstractString)
  allowed = _negative_r_allowed_branches(net)
  for br in net.branchVec
    if br.r_base_pu < 0.0 && !(Int(br.branchIdx) in allowed)
      error("$(context): branch $(getCompName(br.comp)) has a negative base series resistance (r_base = $(br.r_base_pu) pu). The physical equipment impedance must be non-negative; check the branch model, not a FACTS operating point.")
    end
  end
  return nothing
end

"""
    restoreBaseImpedances!(net)

Reset every branch's live series impedance (`r_pu`/`x_pu`) to its physical base
(`r_base_pu`/`x_base_pu`), discarding any series-FACTS operating point that a
power-flow control run stamped in place (#329). Non-FACTS branches are
unchanged (their base equals the live value). Use this to return a net to its
equipment model, for example before a second, independent study on the same
net. Returns `net`.
"""
function restoreBaseImpedances!(net)
  for br in net.branchVec
    br.r_pu = br.r_base_pu
    br.x_pu = br.x_base_pu
  end
  return net
end

"""
    calcBranchYshunt(branch) -> ComplexF64

The total shunt admittance of the branch in per unit.
"""
function calcBranchYshunt(branch::Branch)::ComplexF64
  return (branch.g_pu + branch.b_pu * im)
end

# The shunt arm at one terminal of a branch, as a complex admittance in pu.
_branch_y0_from(branch::Branch)::ComplexF64 = branch.g_from_pu + im * branch.b_from_pu
_branch_y0_to(branch::Branch)::ComplexF64 = branch.g_to_pu + im * branch.b_to_pu

"""
    has_symmetric_shunt(branch::Branch; atol = 1e-12) -> Bool

Whether the two terminal shunt arms of `branch` are equal, the symmetric
pi model that MATPOWER writes as one `BR_B`.
"""
has_symmetric_shunt(branch::Branch; atol::Float64 = 1e-12)::Bool = isapprox(branch.g_from_pu, branch.g_to_pu; atol = atol) && isapprox(branch.b_from_pu, branch.b_to_pu; atol = atol)

"""
    set_branch_shunt!(branch::Branch; g_from_pu, b_from_pu, g_to_pu, b_to_pu) -> Branch

Set the charging admittance of both terminals of `branch` in pu and keep
the totals `g_pu`, `b_pu` equal to their sums. This is the ONE writer of
these fields after construction; no code writes a total alone.
"""
function set_branch_shunt!(branch::Branch; g_from_pu::Real, b_from_pu::Real, g_to_pu::Real, b_to_pu::Real)::Branch
  branch.g_from_pu = Float64(g_from_pu)
  branch.b_from_pu = Float64(b_from_pu)
  branch.g_to_pu = Float64(g_to_pu)
  branch.b_to_pu = Float64(b_to_pu)
  branch.g_pu = branch.g_from_pu + branch.g_to_pu
  branch.b_pu = branch.b_from_pu + branch.b_to_pu
  return branch
end

# Factors of the step whose actuator value lies nearest to `v`; (1, 1, 1, 1)
# without a table on that actuator.
function _tap_step_factors(values::Vector{Float64}, factors::Vector{NTuple{4,Float64}}, v::Float64)::NTuple{4,Float64}
  isempty(values) && return (1.0, 1.0, 1.0, 1.0)
  best = 1
  for k in eachindex(values)
    abs(values[k] - v) < abs(values[best] - v) && (best = k)
  end
  return factors[best]
end

# Set the impedance of a branch with a TapStepImpedance to the step of its
# live tap_ratio and phase_shift_deg. The equipment base r_base_pu/x_base_pu
# follows as well: the table gives the impedance of the transformer at that
# step, not an operating point. No-op without a table.
function _apply_tap_step_impedance!(br::Branch)::Branch
  t = br.tap_step_impedance
  t === nothing && return br
  ra = _tap_step_factors(t.ratio_values, t.ratio_factors, br.tap_ratio)
  ph = _tap_step_factors(t.phase_values, t.phase_factors, br.phase_shift_deg)
  fr, fx, fg, fb = ra[1] * ph[1], ra[2] * ph[2], ra[3] * ph[3], ra[4] * ph[4]
  br.r_pu = t.r_pu * fr
  br.x_pu = t.x_pu * fx
  br.r_base_pu = br.r_pu
  br.x_base_pu = br.x_pu
  set_branch_shunt!(br; g_from_pu = t.g_from_pu * fg, b_from_pu = t.b_from_pu * fb, g_to_pu = t.g_to_pu * fg, b_to_pu = t.b_to_pu * fb)
  return br
end

"""
    set_branch_shunt_total!(branch::Branch; g_pu, b_pu) -> Branch

Set the total charging admittance of `branch` with the symmetric split,
half at each terminal (the MATPOWER form).
"""
set_branch_shunt_total!(branch::Branch; g_pu::Real, b_pu::Real)::Branch = set_branch_shunt!(branch; g_from_pu = 0.5 * g_pu, b_from_pu = 0.5 * b_pu, g_to_pu = 0.5 * g_pu, b_to_pu = 0.5 * b_pu)

"""
    calcBranchRatio(branch) -> ComplexF64

The complex winding ratio of the branch (magnitude and phase shift).
"""
function calcBranchRatio(branch::Branch)::ComplexF64
  ratio = (branch.ratio == 0.0) ? 1.0 : branch.tap_ratio
  shift = (branch.ratio == 0.0) ? 0.0 : branch.phase_shift_deg

  # use ratio even if shift is 0
  if ratio != 1.0 || shift != 0.0
    return calcComplexRatio(tapRatio = ratio, angleInDegrees = shift)
  else
    return 1.0 + 0.0im
  end
end

function getBranchComp(Vn_kV::Float64, from::Int, to::Int, idx::Int, kind::String)
  cTyp = toComponentTyp("Branch")
  name = "B_$(kind)_$(string(convert(Int,trunc(Vn_kV))))_$(Int(from))_$(Int(to))"
  cID = "#" * name * "#" * string(idx)
  return ImpPGMComp(cID, name, cTyp, Vn_kV, from, to)
end

# The element a branch was built from: :line (ACLineSegment), :transformer
# (PowerTransformer) or :pi (BranchModel). getBranchComp writes the kind
# into the generated name (B_ACL_, B_2WT_, B_PI_), and the component type is
# the generic Branch for all three. A branch that carries a name of its own
# (an SCF case restores extra.name) is told by the line marker instead: only
# the line constructor stores ratio 0.0, a transformer asserts ratio != 0.0.
function _branch_build_kind(br::Branch)::Symbol
  name = br.comp.cName
  occursin("_ACL_", name) && return :line
  occursin("_2WT_", name) && return :transformer
  occursin("_PI_", name) && return :pi
  return br.ratio == 0.0 ? :line : :transformer
end

"""
    applyTapNameplate!(branch; tap_step, tap_min_step, tap_max_step,
                       tap_current_step, phase_step_deg, phase_min_step,
                       phase_max_step, phase_current_step, psi_deg = 0.0,
                       phase_du_step = 0.0, context = "tap changer")

Apply nameplate tap-changer data to a transformer branch: the branch's own
`ratio`/`angle` stay the NEUTRAL position, and the current steps move the
LIVE `tap_ratio`/`phase_shift_deg` off it on the mechanical grids the tap
estimation reads as its changers.

- Ratio changer: `tap = neutral / (1 + n * tap_step)` (fraction grid,
  `tap_ratio_grid = :fraction`); with `tap_ratio_grid = :linear` the steps
  lie on the grid of a CGMES changer, `tap = neutral + n * tap_step`
  (`_ratio_grid_shift`), and `ratio` is the position the steps count from.
  A `tap_step` of 0 explicitly declares "no ratio tap changer" and keeps the
  transformer out of the estimator's mass release.
- Phase changer, additional-voltage form (`phase_du_step > 0`): the
  mechanical grid is the additional-voltage amplitude `r2 = n *
  phase_du_step` along the nameplate direction `psi_deg`; the shift angle
  follows from it and is NOT the grid. The factor divides the LIVE
  ratio of the ratio changer (the neutral ratio without one):
  `t = tap_ratio * cis(angle) / (1 + r2 * cis(psi))`. With both changers
  `psi_deg` must not be 180 degrees (0 reads as 90 here), else the two
  positions cannot be told apart from `t` (`ArgumentError`). The SCF import
  builds a stated `alpha_deg` 0 as 0 (a missing key reads as 90) and
  refuses such a pair at 0 as well, unless the Delta-u position is 0.
- Phase changer, degree form (`phase_step_deg > 0`): additive degrees on the
  neutral shift; the band is stored relative to neutral.

`phase_step_deg` and `phase_du_step` are exclusive (a changer has one
mechanical grid). A `NaN` current step means the file carried no step: the
changer starts at its middle position, the neutral position 0 the steps
count from, and `ratio_step_assumed` (`phase_step_assumed`) is set. A phase
changer applied here is declared nameplate data
(`phase_tap_declared`), which releases it for the tap estimation whatever
its psi; a ratio changer applied here (`tap_step > 0`) is declared in the
same way (`ratio_tap_declared`). The tap estimation releases declared
changers only, never the constructor's default grid. `context` prefixes the error messages, so a MATPOWER
`mpc.sparlectra.tap_changers` row and an SCF `tap_changer` entry report in
their own vocabulary. This is the single tap-nameplate application path;
every importer routes through it.
"""
function applyTapNameplate!(
  br::Branch;
  # zero means "this changer is not present", which is what a transformer
  # with only a ratio tap (or only a phase tap) declares; importers pass the
  # full set, a hand-built net passes what it has
  tap_step::Float64 = 0.0,
  tap_min_step::Float64 = 0.0,
  tap_max_step::Float64 = 0.0,
  tap_current_step::Float64 = 0.0,
  phase_step_deg::Float64 = 0.0,
  phase_min_step::Float64 = 0.0,
  phase_max_step::Float64 = 0.0,
  phase_current_step::Float64 = 0.0,
  psi_deg::Float64 = 0.0,
  phase_du_step::Float64 = 0.0,
  tap_ratio_grid::Symbol = :fraction,
  context::AbstractString = "tap changer",
)
  # MATPOWER convention: psi 0 means "not given" and reads as the 90 deg
  # symmetric PST. The SCF import, where a stated 0 means 0, calls
  # _apply_tap_nameplate! with the direction as stated.
  psi = psi_deg != 0.0 ? psi_deg : 90.0
  # with a ratio changer the Delta-u factor sits on its live ratio; along
  # the real axis the two positions cannot be told apart from t
  (tap_step > 0.0 && phase_du_step > 0.0 && abs(sin(deg2rad(psi))) < 1e-9) && throw(ArgumentError("$(context): a ratio changer and a Delta-u changer with psi_deg = $(psi) act along the same real direction; their positions cannot be told apart."))
  return _apply_tap_nameplate!(br, psi; tap_step, tap_min_step, tap_max_step, tap_current_step, phase_step_deg, phase_min_step, phase_max_step, phase_current_step, phase_du_step, tap_ratio_grid, context)
end

# The body of applyTapNameplate! with the direction psi_deg used as given,
# so 0 stays 0. The collinear check of a ratio and Delta-u pair is the
# caller's: applyTapNameplate! refuses the pair along the real axis, the SCF
# import keeps a stated 0 at Delta-u position 0 (factor 1, exact split).
function _apply_tap_nameplate!(
  br::Branch,
  psi_deg::Float64;
  tap_step::Float64,
  tap_min_step::Float64,
  tap_max_step::Float64,
  tap_current_step::Float64,
  phase_step_deg::Float64,
  phase_min_step::Float64,
  phase_max_step::Float64,
  phase_current_step::Float64,
  phase_du_step::Float64,
  tap_ratio_grid::Symbol,
  context::AbstractString,
)
  br.ratio != 0.0 || throw(ArgumentError("$(context): not a transformer branch (its neutral ratio is 0)."))
  # a NaN current step: the file carried no step, so the changer starts at
  # its middle position, the neutral position 0 the steps count from
  ratio_assumed = isnan(tap_current_step)
  phase_assumed = isnan(phase_current_step)
  ratio_assumed && (tap_current_step = 0.0)
  phase_assumed && (phase_current_step = 0.0)
  br.ratio_step_assumed = tap_step > 0.0 && ratio_assumed
  br.phase_step_assumed = (phase_step_deg > 0.0 || phase_du_step > 0.0) && phase_assumed
  tap_ratio_grid in (:fraction, :linear) || throw(ArgumentError("$(context): tap_ratio_grid must be :fraction or :linear, got $(repr(tap_ratio_grid))."))
  base = br.ratio
  if tap_step > 0.0
    (tap_min_step <= 0.0 <= tap_max_step) || throw(ArgumentError("$(context): the neutral position 0 must lie inside [tap_min_step, tap_max_step]."))
    (tap_min_step <= tap_current_step <= tap_max_step) || throw(ArgumentError("$(context): tap_current_step $(tap_current_step) outside [$(tap_min_step), $(tap_max_step)]."))
    br.has_ratio_tap = true
    br.ratio_tap_declared = true
    br.tap_step = tap_step
    br.tap_ratio_grid = tap_ratio_grid
    if tap_ratio_grid === :linear
      # the grid of a CGMES changer, the one the tap controller and the tap
      # estimation step on (whole tap_steps, or the table points if any)
      br.tap_min = _ratio_grid_shift(br, base, tap_min_step; window = (-Inf, Inf))
      br.tap_max = _ratio_grid_shift(br, base, tap_max_step; window = (-Inf, Inf))
      br.tap_ratio = _ratio_grid_shift(br, base, tap_current_step; window = (-Inf, Inf))
    else
      # ratio bounds from the step band on the fraction grid: the highest step
      # gives the smallest multiplier
      br.tap_min = base / (1.0 + tap_max_step * tap_step)
      br.tap_max = base / (1.0 + tap_min_step * tap_step)
      br.tap_ratio = base / (1.0 + tap_current_step * tap_step)
    end
  else
    br.has_ratio_tap = false
    br.ratio_tap_declared = false
  end
  (phase_step_deg > 0.0 && phase_du_step > 0.0) && throw(ArgumentError("$(context): phase_step_deg and phase_du_step are exclusive (a changer has ONE mechanical grid)."))
  if phase_du_step > 0.0
    (phase_min_step <= 0.0 <= phase_max_step) || throw(ArgumentError("$(context): the neutral position 0 must lie inside [phase_min_step, phase_max_step]."))
    (phase_min_step <= phase_current_step <= phase_max_step) || throw(ArgumentError("$(context): phase_current_step $(phase_current_step) outside [$(phase_min_step), $(phase_max_step)]."))
    br.has_phase_tap = true
    br.phase_tap_declared = true
    br.phase_du_step = phase_du_step
    br.phase_du_min_step = phase_min_step
    br.phase_du_max_step = phase_max_step
    br.phase_step_deg = 0.0
    # psi carried in tap_est_alpha_deg: inert while tap_est_mode is :none,
    # and the marker that says this phase changer is DECLARED nameplate
    # data (plain constructor defaults never join the mass release)
    br.tap_est_alpha_deg = psi_deg
    # with a ratio changer the Delta-u factor sits on its live ratio; built
    # on the neutral ratio, the ratio position was lost
    tbase = (tap_step > 0.0 ? br.tap_ratio : base) * cis(deg2rad(br.angle))
    tlive = tbase / (1.0 + phase_current_step * phase_du_step * cis(deg2rad(br.tap_est_alpha_deg)))
    br.tap_ratio = abs(tlive)
    br.phase_shift_deg = rad2deg(angle(tlive))
    # informative degree band (atan of the band ends); the tap estimation
    # works on the Delta-u steps, not on these
    br.phase_min_deg = br.angle - rad2deg(angle(1.0 + phase_max_step * phase_du_step * cis(deg2rad(br.tap_est_alpha_deg))))
    br.phase_max_deg = br.angle - rad2deg(angle(1.0 + phase_min_step * phase_du_step * cis(deg2rad(br.tap_est_alpha_deg))))
  elseif phase_step_deg > 0.0
    (phase_min_step <= 0.0 <= phase_max_step) || throw(ArgumentError("$(context): the neutral position 0 must lie inside [phase_min_step, phase_max_step]."))
    (phase_min_step <= phase_current_step <= phase_max_step) || throw(ArgumentError("$(context): phase_current_step $(phase_current_step) outside [$(phase_min_step), $(phase_max_step)]."))
    br.has_phase_tap = true
    br.phase_tap_declared = true
    br.phase_step_deg = phase_step_deg
    # the band is stored relative to neutral; the tap estimation reads it
    # as the step band of a linear changer (_field_tap_changers)
    br.phase_min_deg = phase_min_step * phase_step_deg
    br.phase_max_deg = phase_max_step * phase_step_deg
    br.phase_shift_deg = br.angle + phase_current_step * phase_step_deg
    br.tap_est_alpha_deg = psi_deg
  else
    br.has_phase_tap = false
    br.phase_tap_declared = false
  end
  return br
end

# Middle position of a tap changer band without a stated neutral step (a
# changer without a transferred step sits at its middle position, the
# neutral step, else the band middle). With an even number of steps the
# lower of the two middle steps.
_tap_band_middle(low::Real, high::Real)::Int = fld(Int(round(low)) + Int(round(high)), 2)

# Import lines for tap changers whose file carried no step: one line per
# changer, or one line listing all of them when there are more than
# _TAP_STEP_ASSUMED_LIST_MAX. `items` are (name, middle step) pairs; the
# importer emits the lines on its own channel (@info, the CGMES messages,
# the PowSyBl report).
const _TAP_STEP_ASSUMED_LIST_MAX = 10
function _tap_step_assumed_lines(items::AbstractVector)::Vector{String}
  isempty(items) && return String[]
  fmt(x) = isinteger(x) ? string(Int(x)) : string(x)
  length(items) <= _TAP_STEP_ASSUMED_LIST_MAX && return String["tap changer $(n): step not transferred, middle position $(fmt(s)) assumed" for (n, s) in items]
  return String[string(length(items), " tap changers: step not transferred, middle position assumed: ", join((string(n, " (", fmt(s), ")") for (n, s) in items), ", "))]
end

# The changers of branch `br` that still sit at the assumed middle position,
# as a text for the export warnings ("ratio", "phase", "ratio and phase"),
# or nothing. Both places an import marks are read: the branch flags of the
# nameplate and the typed models, and the changers of a stored set.
function _tap_step_assumed_text(br::Branch)::Union{Nothing,String}
  set = br.se_tap_changers
  ratio = br.ratio_step_assumed || (set !== nothing && any(ch -> ch.step_assumed && ch.role === :ratio, set.changers))
  phase = br.phase_step_assumed || (set !== nothing && any(ch -> ch.step_assumed && ch.role === :phase, set.changers))
  (ratio || phase) || return nothing
  return ratio && phase ? "ratio and phase" : (ratio ? "ratio" : "phase")
end

# The export warning for a transformer whose step was not transferred: the
# writers that have no mark for it (MATPOWER, CGMES) write the assumed
# middle position as the step (or in the ratio and shift), and the file
# then states it as a real one.
_tap_step_assumed_export_text(name::AbstractString, which::AbstractString)::String = "transformer $(name): $(which) step not transferred, the export writes the assumed middle position"

# Take the mark of an assumed step off the `which` changer (:ratio or
# :phase) of branch `br`: the branch flag and the changers of that role in
# a stored set. A tap controller or a QSTS carry that moves the changer
# calls it, so the step it reaches is a step, not the middle position, and
# the writers do not mark it (the tap estimation does the same in
# _tap_apply_position!). The stored set keeps its steps.
function _clear_tap_step_assumed!(br::Branch, which::Symbol)
  which === :ratio ? (br.ratio_step_assumed = false) : (br.phase_step_assumed = false)
  set = br.se_tap_changers
  (set !== nothing && any(ch -> ch.step_assumed && ch.role === which, set.changers)) || return nothing
  chs = SETapChanger[(ch.step_assumed && ch.role === which) ? SETapChanger(ch.name, ch.mrid, ch.kind, ch.role, ch.side, ch.additive, ch.du, ch.alpha_deg, ch.dphi_deg, ch.step_abs, ch.neutral, ch.low, ch.high, ch.step, ch.table, ch.outside_du, ch.x_curve, ch.x_curve_form, ch.impedance, ch.rx_correction, false) : ch for ch in set.changers]
  br.se_tap_changers = SETapChangerSet(set.ratio_n, set.shift_n_deg, chs, set.origin)
  return nothing
end

# Live positions of a Delta-u phase changer (phase_du_step > 0) from the
# branch tap, the inverse of applyTapNameplate!:
#   tap_ratio * cis(phase_shift_deg) = ratio_live * cis(angle) / (1 + r2 * cis(psi))
# Returns (ratio_live, r2): the live ratio of the ratio changer and the
# additional-voltage amplitude. Without a ratio changer ratio_live is the
# neutral ratio and r2 the projection on cis(psi), the formula the writers
# used before, so those files stay byte-identical. With both changers
# v = ratio * cis(angle) / t = (ratio / ratio_live) * (1 + r2 * cis(psi))
# splits along the real axis and cis(psi). `psi_deg` comes from the caller
# because the MATPOWER writer maps 0 to 90 itself for a single changer.
# With both changers psi 0 or 180 is collinear. The importers never build
# such a pair with a nonzero Delta-u position (applyTapNameplate! refuses
# psi 180, the SCF import a stated 0 off position 0), and
# setTapEstimation! keeps the stored psi, so a psi of 0 here comes from an
# SCF pair stated at 0 and position 0 or from a caller that set the field.
# Then a real v means r2 = 0 exactly, and the split is the ratio alone; a v
# with a phase part cannot be split without the direction and throws.
function _du_tap_split(br::Branch, psi_deg::Float64)::Tuple{Float64,Float64}
  psi = deg2rad(psi_deg)
  base = br.ratio == 0.0 ? 1.0 : br.ratio
  tbase = base * cis(deg2rad(br.angle))
  tlive = (br.tap_ratio == 0.0 ? base : br.tap_ratio) * cis(deg2rad(br.phase_shift_deg))
  v = tbase / tlive
  (br.has_ratio_tap && br.tap_step > 0.0) || return (base, real((v - 1.0) * cis(-psi)))
  sa = sin(psi)
  if abs(sa) < 1e-9
    abs(imag(v)) <= 1e-12 * abs(v) && return (base / real(v), 0.0)
    throw(ArgumentError("branch $(br.fromBus)-$(br.toBus): a ratio changer and a Delta-u changer at psi = $(psi_deg) deg with a Delta-u position off 0: the two positions cannot be told apart. Set the nameplate psi of the Delta-u changer (tap_est_alpha_deg) or re-import the case before writing it."))
  end
  u = imag(v) / sa                 # u = (ratio / ratio_live) * r2
  k = real(v) - u * cos(psi)       # k = ratio / ratio_live
  return (base / k, u / k)
end
