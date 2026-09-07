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

# Date: 2026-09-03
# file: examples/others/exp_demo_case_tour.jl
# purpose: load the shipped demo case sp_case14 (data/scf, self-built, no
#          downloads) and run all five run kinds on it: power flow, N-1
#          contingencies, the file's own scenarios, state estimation on the
#          measurements the file carries, and an IEC 60909 short circuit

using Sparlectra
using Printf

include(joinpath(@__DIR__, "example_header.jl"))

function demo_case_tour(; case_name::String = "sp_case14")
  print_example_banner("examples/others/exp_demo_case_tour.jl", "five run kinds on the shipped demo case $(case_name)")
  path = joinpath(pkgdir(Sparlectra), "data", "scf", "$(case_name).scf.json")

  # 1) power flow: the shipped file carries a solved start_state, so the
  # solve settles immediately
  res = Sparlectra.import_case(path, Sparlectra.SparlectraConfig(Dict()))
  net = res.net
  ite, erg = runpf!(net; verbose = 0)
  calcNetLosses!(net)
  p_losses, _ = Sparlectra.getTotalLosses(net = net)
  @printf("power flow      : converged=%s iterations=%d losses=%.3f MW\n", erg == 0, ite, p_losses)

  # 2) N-1 contingencies over every in-service branch
  n1net = Sparlectra.import_case(path, Sparlectra.SparlectraConfig(Dict())).net
  n1 = runContingencies!(n1net, generateN1Branches(n1net))
  islanded = count(r -> !r.converged, n1)
  @printf("N-1             : %d cases, %d converged, %d islanded/failed\n", length(n1), count(r -> r.converged, n1), islanded)

  # 3) the file's own explicit scenarios (the demo cases ship a scenarios
  # block; the same list drives the Web UI's "case file scenarios" source)
  icase = read_scf_json(path)
  iset = Sparlectra.scf_case_scenarios(icase)
  scen_net = Sparlectra.import_case(path, Sparlectra.SparlectraConfig(Dict())).net
  scen = runScenarios!(scen_net, Sparlectra.ScenarioSet(scenarios = iset.scenarios, mode = :explicit); index = Sparlectra.ScenarioIndex(icase))
  for r in scen
    load_txt = (r.max_branch_loading_pct !== nothing && isfinite(r.max_branch_loading_pct)) ? @sprintf("max loading %.1f%%", r.max_branch_loading_pct) : "islanding case"
    @printf("scenario        : %-24s converged=%-5s %s\n", r.name, r.converged, load_txt)
  end

  # 4) state estimation on the measurements the file itself carries
  se_net = Sparlectra.import_case(path, Sparlectra.SparlectraConfig(Dict())).net
  se = runse!(se_net, Vector{Sparlectra.Measurement}(se_net.measurements), Sparlectra.StateEstimationConfig())
  @printf("state estimation: converged=%s J=%.2f dof=%d (J/dof %.2f)\n", se.converged, se.objectiveJ, se.dof, se.objectiveJ / se.dof)

  # 5) short circuit at the case's study buses, from the file's own feeder
  # record (the external grid doubles as the IEC 60909 source)
  sc = runShortCircuit!(net, net.sc_sources; buses = ["Ostheim_110", "Moorau_20"], case = :max)
  for row in sc.rows
    @printf("short circuit   : bus %-4s Ik'' = %.3f kA (Sk'' = %.1f MVA)\n", row.bus, row.ik_kA, row.sk_MVA)
  end
  return nothing
end

Base.invokelatest(demo_case_tour)
