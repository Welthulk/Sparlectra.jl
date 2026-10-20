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

# file: docs/lit/warmup/workshop_qsts.jl
# purpose: compile warm-up of the QSTS workshop (docs/lit/workshop_qsts.jl):
#          runs every path the chapters use once on the first two steps of
#          the shipped series sp_case14_day, so the chapters do not stall on
#          first-call compilation. Included by the first code cell of the
#          notebook. Not part of the library.

"""
    warmup()

Run every compile-heavy path of the QSTS workshop once on two steps of the
example series and print how long it took. The output of the paths is
discarded.
"""
function warmup()
  println("warm-up: compiles the paths of this workshop once; how long it takes depends on the machine (a Colab session is several times slower than a desktop)")
  ## every path of the chapters once on two steps, output discarded: the
  ## series read, mode off, mode A with and without controllers, mode B,
  ## the accessors, the chart and the result files
  t_paths = @elapsed redirect_stdout(devnull) do
    wdata = joinpath(dirname(dirname(pathof(Sparlectra))), "data", "qsts")
    wfull = read_scf_series(joinpath(wdata, "sp_case14_day.qsts.zip"))
    wseries = SCFSeries(name = "warm-up", base = wfull.base, steps = wfull.steps[1:2], ops = wfull.ops[1:2], step_minutes = 15)
    wcontrol = joinpath(wdata, "sp_case14_day.control.yaml")
    try
      runQsts(wseries; config = SparlectraConfig())
    catch err
      err isa ArgumentError || rethrow()
    end
    a = runQsts(wseries; config = SparlectraConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "snapshots"))))
    qstsValues(a, "8", "vm_pu")
    qstsValues(a, "26", "p_from_mw")
    qstsValues(a, "", "total_p_loss_mw")
    ac = runQsts(wseries; config = load_sparlectra_config(wcontrol; reload = true, overrides = Dict("qsts" => Dict("mode" => "snapshots", "snapshot_controls" => true))))
    qstsValues(ac, "30", "tap_ratio")
    b = runQsts(wseries; config = load_sparlectra_config(wcontrol; reload = true, overrides = Dict("qsts" => Dict("mode" => "sequential"))))
    qstsValues(b, "30#OLTC", "position")
    qstsSwitching(b)
    qstsMetrics(b)
    qstsChart(b, "14", "vm_pu"; band = (0.99, 1.01))
    writeQstsResults(mktempdir(), b; format = "technical")
  end
  println("warm-up: ", round(t_paths; digits=2), " s (series, mode A, mode B, accessors, chart, files); everything warm")
  return nothing
end
