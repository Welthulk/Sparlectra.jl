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

# file: src/adapters/powsybl/powsybl_report.jl
# purpose: the import report of the PowSyBl builder: counts per element
#          type (read, built, skipped), one entry per skipped element with
#          its reason, the slack choice per synchronous component, and the
#          notices (remote regulation held local, contracted switches,
#          limit types skipped). The console block follows the CGMES
#          importer's style: count lines first, then the decisions, then
#          one line per skipped element.

"""
    PowsyblImportReport

What the builder did with a bundle. `counts` maps an element kind
(`"bus"`, `"line"`, `"generator"`, ...) to `(read, built, skipped)`;
`skipped` lists every element that was not built with its id, kind and
reason; `slack` records the slack decision per synchronous component;
`messages` holds the notices. `bus_substation` and `bus_view` keep the
bus metadata the `Net` has no field for (substation id and bus-view id
per bus-breaker bus), `temporary_limits` the temporary current limits per
branch element id as `(duration_s, value_A)` pairs.
"""
Base.@kwdef struct PowsyblImportReport
  counts::Dict{String,NamedTuple{(:read, :built, :skipped),NTuple{3,Int}}} = Dict{String,NamedTuple{(:read, :built, :skipped),NTuple{3,Int}}}()
  skipped::Vector{NamedTuple{(:id, :kind, :reason),Tuple{String,String,String}}} = NamedTuple{(:id, :kind, :reason),Tuple{String,String,String}}[]
  slack::Vector{NamedTuple{(:component, :generator, :bus, :reason),Tuple{Int,String,String,String}}} = NamedTuple{(:component, :generator, :bus, :reason),Tuple{Int,String,String,String}}[]
  messages::Vector{String} = String[]
  bus_substation::Dict{String,String} = Dict{String,String}()
  bus_view::Dict{String,String} = Dict{String,String}()
  temporary_limits::Dict{String,Vector{Tuple{Int,Float64}}} = Dict{String,Vector{Tuple{Int,Float64}}}()
end

# The order the console block lists the element kinds in.
const _POWSYBL_REPORT_KINDS = String["bus", "switch", "line", "2wt", "3wt", "generator", "load", "shunt", "svc", "dangling_line", "tie_line", "hvdc"]

function _powsybl_count!(report::PowsyblImportReport, kind::AbstractString; read::Int = 0, built::Int = 0, skipped::Int = 0)
  c = get(report.counts, String(kind), (read = 0, built = 0, skipped = 0))
  report.counts[String(kind)] = (read = c.read + read, built = c.built + built, skipped = c.skipped + skipped)
  return nothing
end

function _powsybl_skip!(report::PowsyblImportReport, kind::AbstractString, id::AbstractString, reason::AbstractString)
  push!(report.skipped, (id = String(id), kind = String(kind), reason = String(reason)))
  _powsybl_count!(report, kind; skipped = 1)
  return nothing
end

_powsybl_notice!(report::PowsyblImportReport, text::AbstractString) = (push!(report.messages, "notice: " * String(text)); nothing)

"""
    format_powsybl_report(report::PowsyblImportReport) -> String

The console block of an import: one count line per element kind (read,
built, skipped), the slack decision per synchronous component, the
notices, then one line per skipped element with id, kind and reason.
"""
function format_powsybl_report(report::PowsyblImportReport)::String
  io = IOBuffer()
  println(io, "PowSyBl import")
  kinds = vcat([k for k in _POWSYBL_REPORT_KINDS if haskey(report.counts, k)], sort([k for k in keys(report.counts) if !(k in _POWSYBL_REPORT_KINDS)]))
  width = maximum(length, kinds; init = 4)
  for kind in kinds
    c = report.counts[kind]
    println(io, "  ", rpad(kind, width), "  read ", lpad(c.read, 5), "  built ", lpad(c.built, 5), "  skipped ", lpad(c.skipped, 5))
  end
  for s in report.slack
    println(io, "  slack component ", s.component, ": ", isempty(s.generator) ? "none" : s.generator * " at bus " * s.bus, " (", s.reason, ")")
  end
  for m in report.messages
    println(io, "  ", m)
  end
  for s in report.skipped
    println(io, "  skip: ", s.kind, " ", s.id, ": ", s.reason)
  end
  return String(take!(io))
end
