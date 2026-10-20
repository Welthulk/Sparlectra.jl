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

# file: src/qsts/qsts_chart.jl
# purpose: the time series chart of one quantity of one element of a QSTS
#          result as a self-contained SVG (no JavaScript, a tooltip per point
#          through <title>, the band dashed); one renderer for the service
#          artefacts, the Web UI view and the workshop.

"""
    QstsChart

The chart of one quantity of one element over the steps of a
[`QstsResult`](@ref), built by [`qstsChart`](@ref): `element_id`,
`quantity`, `labels` (axis label per step), `values` (`NaN` where the step
has no value), `band` (`(low, high)` with `NaN` for an open side, `nothing`
without a band) and `svg`, the chart as SVG text. `show(io,
MIME"image/svg+xml"(), c)` writes the SVG, so notebooks and Documenter
display it directly.
"""
struct QstsChart
  element_id::String
  quantity::String
  labels::Vector{String}
  values::Vector{Float64}
  band::Union{Nothing,Tuple{Float64,Float64}}
  svg::String
end

Base.show(io::IO, ::MIME"image/svg+xml", c::QstsChart) = print(io, c.svg)
Base.show(io::IO, c::QstsChart) = print(io, "QstsChart(", repr(c.element_id), ", ", repr(c.quantity), ", ", count(!isnan, c.values), " of ", length(c.values), " steps)")

# the band :auto gives a quantity: the voltage band and the loading limit of
# the run's qsts section, nothing for the other quantities
function _qsts_chart_band(r::QstsResult, quantity::AbstractString, band)::Union{Nothing,Tuple{Float64,Float64}}
  band === nothing && return nothing
  band isa Tuple && return (Float64(band[1]), Float64(band[2]))
  band === :auto || throw(ArgumentError(string("qstsChart: band must be :auto, nothing or a (low, high) tuple; got ", repr(band))))
  quantity == "vm_pu" && return (r.qsts.vm_min_pu, r.qsts.vm_max_pu)
  quantity == "loading_pct" && return (NaN, r.qsts.loading_max_pct)
  return nothing
end

# XML text escape for names and labels in the SVG
_qsts_svg_escape(s::AbstractString)::String = replace(String(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")

_qsts_svg_num(x::Real)::String = @sprintf("%.6g", x)

# The SVG text of one series: the renderer behind qstsChart, also called by
# the Web UI chart route, which holds the series of its view artifact and no
# QstsResult. `b` is the band as qstsChart resolved it.
function _qsts_chart_svg(element_id::AbstractString, quantity::AbstractString, labels::Vector{String}, values::Vector{Float64}, b::Union{Nothing,Tuple{Float64,Float64}},
                         width::Integer, height::Integer)::String
  n = length(values)
  W = Float64(width)
  H = Float64(height)
  left, right, top, bottom = 64.0, 16.0, 28.0, 36.0
  finite = filter(!isnan, values)
  bounds = copy(finite)
  b === nothing || append!(bounds, filter(!isnan, [b[1], b[2]]))
  vlo, vhi = isempty(bounds) ? (0.0, 1.0) : extrema(bounds)
  spread = vhi == vlo ? 1.0 : vhi - vlo
  # 5 percent margin above and below the values and the band
  lo = (vhi == vlo ? vlo - 0.5 : vlo) - 0.05 * spread
  hi = (vhi == vlo ? vhi + 0.5 : vhi) + 0.05 * spread
  xof(k) = n <= 1 ? left + (W - left - right) / 2 : left + (k - 1) * (W - left - right) / (n - 1)
  yof(v) = top + (hi - v) * (H - top - bottom) / (hi - lo)
  title = isempty(element_id) ? string("step ", quantity) : string(element_id, " ", quantity)
  io = IOBuffer()
  println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"", Int(width), "\" height=\"", Int(height), "\" viewBox=\"0 0 ", Int(width), " ", Int(height), "\" font-family=\"sans-serif\" font-size=\"11\">")
  println(io, "<title>", _qsts_svg_escape(title), "</title>")
  println(io, "<rect x=\"0\" y=\"0\" width=\"", Int(width), "\" height=\"", Int(height), "\" fill=\"white\"/>")
  println(io, "<text x=\"", _qsts_svg_num(left), "\" y=\"16\" font-size=\"13\">", _qsts_svg_escape(title), "</text>")
  # axes and their labels: value range left, first and last step below
  println(io, "<line x1=\"", _qsts_svg_num(left), "\" y1=\"", _qsts_svg_num(top), "\" x2=\"", _qsts_svg_num(left), "\" y2=\"", _qsts_svg_num(H - bottom), "\" stroke=\"#444\"/>")
  println(io, "<line x1=\"", _qsts_svg_num(left), "\" y1=\"", _qsts_svg_num(H - bottom), "\" x2=\"", _qsts_svg_num(W - right), "\" y2=\"", _qsts_svg_num(H - bottom), "\" stroke=\"#444\"/>")
  println(io, "<text x=\"", _qsts_svg_num(left - 4), "\" y=\"", _qsts_svg_num(top + 4), "\" text-anchor=\"end\">", _qsts_svg_num(hi), "</text>")
  println(io, "<text x=\"", _qsts_svg_num(left - 4), "\" y=\"", _qsts_svg_num(H - bottom), "\" text-anchor=\"end\">", _qsts_svg_num(lo), "</text>")
  if n >= 1
    println(io, "<text x=\"", _qsts_svg_num(xof(1)), "\" y=\"", _qsts_svg_num(H - bottom + 16), "\" text-anchor=\"start\">", _qsts_svg_escape(labels[1]), "</text>")
    println(io, "<text x=\"", _qsts_svg_num(xof(n)), "\" y=\"", _qsts_svg_num(H - bottom + 16), "\" text-anchor=\"end\">", _qsts_svg_escape(labels[n]), "</text>")
  end
  # the band, dashed, one line per closed side
  if b !== nothing
    for (side, v) in (("low", b[1]), ("high", b[2]))
      isnan(v) && continue
      println(io, "<line class=\"band\" x1=\"", _qsts_svg_num(left), "\" y1=\"", _qsts_svg_num(yof(v)), "\" x2=\"", _qsts_svg_num(W - right), "\" y2=\"", _qsts_svg_num(yof(v)),
              "\" stroke=\"#c33\" stroke-dasharray=\"6 4\"><title>", side, " ", _qsts_svg_num(v), "</title></line>")
    end
  end
  # the line, broken at steps without a value, then one point per value
  segment = String[]
  flush_segment() = (length(segment) >= 2 && println(io, "<polyline fill=\"none\" stroke=\"#1f5fa8\" stroke-width=\"1.5\" points=\"", join(segment, " "), "\"/>"); empty!(segment))
  for k in 1:n
    if isnan(values[k])
      flush_segment()
    else
      push!(segment, string(_qsts_svg_num(xof(k)), ",", _qsts_svg_num(yof(values[k]))))
    end
  end
  flush_segment()
  for k in 1:n
    isnan(values[k]) && continue
    println(io, "<circle cx=\"", _qsts_svg_num(xof(k)), "\" cy=\"", _qsts_svg_num(yof(values[k])), "\" r=\"2.5\" fill=\"#1f5fa8\"><title>step ", k,
            isempty(labels[k]) ? "" : string(" (", _qsts_svg_escape(labels[k]), ")"), ": ", _qsts_svg_num(values[k]), "</title></circle>")
  end
  print(io, "</svg>\n")
  return String(take!(io))
end

"""
    qstsChart(r::QstsResult, element_id, quantity; band = :auto, width = 720, height = 260) -> QstsChart

The time series of one quantity of one element (an empty `element_id`
selects a step quantity, see [`qstsValues`](@ref)) as an SVG chart: one
point per step with a value, joined by a line that breaks at a step without
a value, each point with a tooltip (`<title>`: step, label, value), the
axis labelled with the first and the last step label and the value range.
`band = :auto` draws the voltage band `qsts.vm_min_pu`/`qsts.vm_max_pu` for
`vm_pu` and the limit `qsts.loading_max_pct` for `loading_pct` as dashed
lines; a `(low, high)` tuple draws that band (`NaN` leaves a side open),
`nothing` none. The SVG holds no script and no time stamp, so equal results
give equal bytes. Throws `ArgumentError` for an unknown element or quantity
and for a band of another form.
"""
function qstsChart(r::QstsResult, element_id::AbstractString, quantity::AbstractString; band = :auto, width::Integer = 720, height::Integer = 260)::QstsChart
  values = qstsValues(r, element_id, quantity)
  b = _qsts_chart_band(r, quantity, band)
  labels = String[s.label for s in r.steps]
  return QstsChart(String(element_id), String(quantity), labels, values, b, _qsts_chart_svg(element_id, quantity, labels, values, b, width, height))
end
