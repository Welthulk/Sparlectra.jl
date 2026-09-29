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

# file: docs/lit/warmup/workshop_scenarios.jl
# purpose: compile warm-up of the scenarios workshop
#          (docs/lit/workshop_scenarios.jl): runs every path the chapters
#          use once on the small shipped sp_case9, so the chapters do not
#          stall on first-call compilation. Included by the first code cell
#          of the notebook after it has set up the study case. Not part of
#          the library.

"""
    warmup()

Run every compile-heavy path of the scenarios workshop once on the 9-bus
case and print how long it took. The output of the paths is discarded.
"""
function warmup()
  println("warm-up: compiles the paths of this workshop once; how long it takes depends on the machine (a Colab session is several times slower than a desktop)")
  ## every path of the chapters once on the 9-bus case, output discarded:
  ## the MATPOWER import, the SCF export and read, the N-1 expansion,
  ## hand-written scenarios (including an islanding double outage and a
  ## scenario that does not converge, the two rows of the study case that
  ## otherwise compile in the screening chapter), runScenarios! with and
  ## without screening and writing the scenarios block back
  t_paths = @elapsed redirect_stdout(devnull) do
    wdir = mktempdir()
    wm = joinpath(dirname(dirname(pathof(Sparlectra))), "data", "mpower", "sp_case9.m")
    wnet = createNetFromMatPowerFile(filename=wm, flatstart=false, enable_pq_gen_controllers=true, bus_shunt_model=:admittance, matpower_shift_sign=1.0, matpower_shift_unit=:deg, matpower_ratio=:normal, tap_changer_model=:ideal)
    wscf = joinpath(wdir, "sp_case9_warmup.scf.json")
    exportSCF(wnet; file=wscf, case_name="sp_case9_warmup", source_reference="sp_case9.m")
    wcase = Sparlectra.read_scf_json(wscf)
    windex = ScenarioIndex(wcase)
    expand_scenarios(ScenarioSet(mode=:n1_all), wnet, windex)
    expand_scenarios(ScenarioSet(mode=:n1_branches), wnet, windex)
    wset = ScenarioSet(scenarios=[
        ## a double outage that splits off an island without a reference
        ## (the reference takeover the study case hits)
        Scenario(name="double outage", weight=2.0, ops=[PatchOp(op=:status, target=:branch, id=wcase.data.line[1].id, value=0.0), PatchOp(op=:status, target=:branch, id=wcase.data.line[2].id, value=0.0)]),
        ## a load scaling the solver cannot carry (the non-converged row)
        Scenario(name="overload", ops=[PatchOp(op=:scale, target=:load, id=first(id for (id, k) in windex.kind_by_id if k === :load), factor=20.0)]),
        Scenario(name="scale", ops=[PatchOp(op=:scale, target=:load, id=first(id for (id, k) in windex.kind_by_id if k === :load), factor=1.2)]),
        Scenario(name="setpoint", ops=[PatchOp(op=:set, target=:generator, id=first(id for (id, k) in windex.kind_by_id if k === :generator), field=:p, value=25.0)]),
    ])
    validate_scenarios(wset, windex)
    ## the result rows printed the way the chapters print them
    for r in runScenarios!(wnet, wset; index=windex)
      println(r.name, ": converged = ", r.converged, ", vmin = ", round(r.min_vm_pu; digits=4), " pu, iterations = ", r.iterations)
    end
    runScenarios!(wnet, ScenarioSet(mode=:n1_all); index=windex)
    runScenarios!(wnet, ScenarioSet(mode=:n1_all); index=windex, screening_mode=:flag)
    wcase.sparlectra.scenarios = scenario_set_dict(wset)
    Sparlectra.write_scf_json(wcase, wscf)
    scf_case_scenarios(Sparlectra.read_scf_json(wscf))
  end
  println("warm-up: ", round(t_paths; digits=2), " s (import, SCF, N-1 expansion, scenarios, screening); everything warm")
  return nothing
end
