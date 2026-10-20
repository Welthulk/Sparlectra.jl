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
# file: examples/others/exp_dynamic_data.jl
# purpose: DynamicData on the generators of sp_case14: the optional dynamic
#          parameters (inertia, droop, reserves, unit type), the validation
#          error and the SCF round trip (format 1.3). No calculation reads
#          the fields. A CGMES delivery with a DY profile fills them with
#          cgmes_import.read_dynamics = true.
#          julia --project=. examples/others/exp_dynamic_data.jl

using Sparlectra

include(joinpath(@__DIR__, "..", "others", "example_header.jl"))

function main()
  print_example_banner("examples/others/exp_dynamic_data.jl", "DynamicData: dynamic parameters of a unit and their SCF round trip")
  net = importSCF(joinpath(dirname(dirname(pathof(Sparlectra))), "data", "scf", "sp_case14.scf.json"))
  gens = [ps for ps in net.prosumpsVec if Sparlectra.isGenerator(ps)]

  # 1. a synchronous machine with primary control, and a grid-following
  #    converter; the unit type follows the inertia when it is not given
  gens[1].dynamics = DynamicData(inertia_s = 5.0, damping_pu = 1.0, droop_pu = 0.05, primary_reserve_mw = 20.0, governor_time_s = 0.5)
  gens[2].dynamics = DynamicData(agc_participation = 0.4, secondary_reserve_mw = 10.0)
  for ps in gens[1:2]
    println(Sparlectra.getCompName(ps.comp), ": ", ps.dynamics)
  end

  # 2. a value against the rule of its field is a named error
  try
    DynamicData(droop_pu = 0.0)
  catch err
    # the one failure this call is built to show
    err isa DynamicDataError || rethrow()
    println("refused: ", sprint(showerror, err))
  end

  # 3. SCF carries the records bit for bit; a unit without them writes none
  file = exportSCF(net; file = joinpath(mktempdir(), "sp_case14_dyn.scf.json"))
  back = importSCF(file)
  same = all(k -> back.prosumpsVec[k].dynamics == net.prosumpsVec[k].dynamics, eachindex(net.prosumpsVec))
  println("SCF round trip: ", same ? "every record equal" : "DIFFERENT", ", units with DynamicData: ", count(ps -> ps.dynamics !== nothing, back.prosumpsVec))
  return nothing
end

Base.invokelatest(main)
