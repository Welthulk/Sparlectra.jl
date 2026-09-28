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

# file: src/tap_changer_kinds.jl
# purpose: the one mapping table from the tap-changer names of the source
#          formats (CGMES classes, PowSyBl tables, DTF regulator types,
#          MATPOWER nameplate rows) to Sparlectra's tap-changer kinds; the
#          importers and exporters call tap_changer_kind instead of
#          spelling the names themselves.

"""
    TAP_CHANGER_KIND_TABLE

One row per source name: `source`, `name`, the Sparlectra `kind`
(`:symmetrical`, `:asymmetrical`, `:tabular` for a `PhaseTapChangerModel`,
`:ratio` for a `PowerTransformerTaps`, `:grid` for the branch-level grid
without a model), the `model` family (`:phase`, `:ratio`, `:grid`), the
default winding angle `psi_deg` of an asymmetrical changer (`nothing`
when the nameplate states it), and a note. CGMES `PhaseTapChangerLinear`
is a `:tabular` model with generated points (there is no fourth kind);
the DTF regulator types all map to `:asymmetrical` with their nameplate
angle, the longitudinal regulator included (angle 0), because that is the
regulating vector the DTF importer reproduces.
"""
const TAP_CHANGER_KIND_TABLE = [
  (source = :cgmes, name = "PhaseTapChangerSymmetrical", kind = :symmetrical, model = :phase, psi_deg = nothing, note = ""),
  (source = :cgmes, name = "PhaseTapChangerAsymmetrical", kind = :asymmetrical, model = :phase, psi_deg = nothing, note = "psi from windingConnectionAngle"),
  (source = :cgmes, name = "PhaseTapChangerLinear", kind = :tabular, model = :phase, psi_deg = nothing, note = "points generated from stepPhaseShiftIncrement over lowStep..highStep"),
  (source = :cgmes, name = "PhaseTapChangerTabular", kind = :tabular, model = :phase, psi_deg = nothing, note = ""),
  (source = :cgmes, name = "RatioTapChanger", kind = :ratio, model = :ratio, psi_deg = nothing, note = "with or without a table"),
  (source = :powsybl, name = "phase_tap_changer_steps", kind = :tabular, model = :phase, psi_deg = nothing, note = "follow-up, the importer uses the _at_current_tap columns"),
  (source = :powsybl, name = "ratio_tap_changer_steps", kind = :ratio, model = :ratio, psi_deg = nothing, note = "follow-up"),
  (source = :dtf, name = "Querregler", kind = :asymmetrical, model = :phase, psi_deg = 90.0, note = "quadrature booster"),
  (source = :dtf, name = "Schraegregler", kind = :asymmetrical, model = :phase, psi_deg = nothing, note = "psi from the nameplate: 0, 30, 60, 90"),
  (source = :dtf, name = "Laengsregler", kind = :asymmetrical, model = :phase, psi_deg = 0.0, note = "longitudinal regulator, regulating vector 1 + f"),
  (source = :matpower, name = "tap_changers", kind = :grid, model = :grid, psi_deg = nothing, note = "branch grid, no model"),
]

"""
    tap_changer_kind(source::Symbol, name::AbstractString; angle = nothing) -> NamedTuple

The Sparlectra kind of a source's tap changer from
[`TAP_CHANGER_KIND_TABLE`](@ref): `(kind, model, psi_deg)`. `angle`
supplies the winding angle of an asymmetrical changer when the nameplate
states it (CGMES `windingConnectionAngle`, the DTF angle); without it the
table's default applies. Throws an `ArgumentError` for a name the table
does not carry, so an importer cannot silently invent a kind.
"""
function tap_changer_kind(source::Symbol, name::AbstractString; angle = nothing)
  for row in TAP_CHANGER_KIND_TABLE
    (row.source === source && row.name == String(name)) || continue
    psi = row.kind === :asymmetrical ? (angle === nothing ? row.psi_deg : Float64(angle)) : nothing
    return (kind = row.kind, model = row.model, psi_deg = psi)
  end
  throw(ArgumentError("tap_changer_kind: no mapping for source $(source) and name $(repr(String(name)))"))
end
