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

# file: tools/build_powsybl_sc_demo.jl
# purpose: builds data/powsybl/ieee14_sc.xiidm, the IEEE 14-bus network
#          with short-circuit data on its generators, from the shipped
#          ieee14.xiidm (which stays untouched: it is the OpenLoadFlow
#          reference). The IEEE case defines neither rated powers nor
#          reactances; the values below are ASSUMED for demonstration
#          (x''d = 0.2 pu on the machine base, ratings sized to the
#          scheduled output and the reactive range) and are written into
#          the IIDM file as the attribute ratedS and the extension
#          generatorShortCircuit, the form a PowSyBl export carries them.
#          Run: julia --project=. tools/build_powsybl_sc_demo.jl

using Sparlectra

const ROOT = dirname(@__DIR__)
const SOURCE = joinpath(ROOT, "data", "powsybl", "ieee14.xiidm")
const TARGET = joinpath(ROOT, "data", "powsybl", "ieee14_sc.xiidm")

# generator id => (rated power in MVA, nominal voltage of its level in kV)
const RATINGS = ["B1-G" => (300.0, 135.0), "B2-G" => (80.0, 135.0), "B3-G" => (60.0, 135.0), "B6-G" => (40.0, 12.0), "B8-G" => (40.0, 20.0)]
const XDPP_PU = 0.2

function main()
  text = read(SOURCE, String)
  extensions = IOBuffer()
  for (id, (rated, vn)) in RATINGS
    occursin("<iidm:generator id=\"$(id)\"", text) || error("generator $(id) not found in $(SOURCE)")
    text = replace(text, "<iidm:generator id=\"$(id)\"" => "<iidm:generator id=\"$(id)\" ratedS=\"$(rated)\"")
    x_ohm = round(XDPP_PU * vn^2 / rated; digits = 6)
    println(extensions, "    <iidm:extension id=\"$(id)\">")
    println(extensions, "        <gsc:generatorShortCircuit xmlns:gsc=\"http://www.powsybl.org/schema/iidm/ext/generator_short_circuit/1_0\" directSubtransX=\"$(x_ohm)\"/>")
    println(extensions, "    </iidm:extension>")
  end
  text = replace(text, "</iidm:network>" => String(take!(extensions)) * "</iidm:network>")
  write(TARGET, text)
  # the file must read: the reader names what it cannot map
  tables = Sparlectra.read_iidm_tables(TARGET)
  println("wrote ", TARGET, ": ", length(tables.generators.id), " generators, ", count(isfinite, tables.generators.direct_subtrans_x), " with short-circuit data")
  return nothing
end

Base.invokelatest(main)
