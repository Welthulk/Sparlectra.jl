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

# file: tools/gen_smallgrid_measurements.jl
# purpose: build the synthetic measurement-table fixture
#          test/fixtures/measurements/smallgrid_synthetic.table.csv from the
#          ENTSO-E SmallGrid (NodeBreaker BaseCase, test/fixtures/cgmes_smallgrid):
#          two snapshots generated from the SV state with
#          generateMeasurementTableFromSV, plus a few hand-made rows that
#          break one rule each, so the import report has something to count.
#          Deterministic: equal seeds give an identical file.
#
# usage: julia --project=. tools/gen_smallgrid_measurements.jl

using Sparlectra

const _ROOT = normpath(joinpath(@__DIR__, ".."))
const _CASE = joinpath(_ROOT, "test", "fixtures", "cgmes_smallgrid")
const _OUT = joinpath(_ROOT, "test", "fixtures", "measurements", "smallgrid_synthetic.table.csv")
const _T1 = "2026-01-01T12:00Z"
const _T2 = "2026-01-01T13:00Z"

function main()
  files = sort!(filter(f -> endswith(f, ".zip"), readdir(_CASE; join = true)))
  result = redirect_stdout(() -> importCGMES(path = files), devnull)
  # snapshot 1: branch ends and nodes with 0.5 percent gross errors of 20
  # sigma; snapshot 2: branch ends, nodes and injections without gross errors
  s1 = generateMeasurementTableFromSV(result; timestamp = _T1, bad_data_fraction = 0.005, bad_data_sigmas = 20.0, seed = 1)
  s2 = generateMeasurementTableFromSV(result; timestamp = _T2, include_injections = true, seed = 2)
  Row = eltype(s1)
  rows = Row[]
  # one healthy branch-end row of snapshot 1 is switched off: it must come
  # back as an inactive measurement, not as a problem
  index = cgmesMeasurementKeys(result)
  flipped = false
  for r in s1
    if !flipped && r.source == "synthetic" && !isempty(r.terminal_mrid) && index.terminals[r.terminal_mrid].kind === :branch
      push!(rows, merge(r, (valid = false,)))
      flipped = true
    else
      push!(rows, r)
    end
  end
  # rows that break one rule each (snapshot 1)
  breaker = first(sort!([k for (k, t) in index.terminals if t.kind === :not_modelled && startswith(t.detail, "Breaker is not")]))
  node = first(sort!([k for (k, b) in index.nodes if b > 0]))
  push!(rows, Row((_T1, "", "00000000-0000-0000-0000-000000000000", "", "P", 1.0, 0.1, "handmade;unknown_key", true)))
  push!(rows, Row((_T1, "", breaker, "", "P", 1.0, 0.1, "handmade;not_modelled", true)))
  push!(rows, Row((_T1, "", "", node, "I", 1.0, 0.1, "handmade;unknown_type", true)))
  append!(rows, s2)
  comments = [
    "synthetic measurements of the ENTSO-E SmallGrid (NodeBreaker BaseCase), built by tools/gen_smallgrid_measurements.jl",
    "snapshot $(_T1): branch ends and nodes, seed 1, noise sigma P/Q max(1 % of reading, 0.1), U 0.2 % of reading, 0.5 % gross errors of 20 sigma (source synthetic;bad_data)",
    "snapshot $(_T1): one healthy branch row with valid=false, three handmade rows that break one rule each (source handmade;<reason>)",
    "snapshot $(_T2): branch ends, nodes and injections (SvPowerFlow), seed 2, same noise, no gross errors",
  ]
  writeMeasurementTable(_OUT, rows; comments = comments)
  println("wrote ", length(rows), " rows to ", _OUT)
end

Base.invokelatest(main)
