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

# file: docs/lit/warmup/workshop_tour_cgmes.jl
# purpose: compile warm-up of the foreign-formats tour
#          (docs/lit/workshop_tour_cgmes.jl): runs every path the chapters
#          use once (CGMES summary, analysis, import, SV comparison, export,
#          the IIDM reader and the DTF deck), so the chapters do not stall on
#          first-call compilation. Included by the first code cell of the
#          notebook, which defines the data paths used here (microgrid_be,
#          microgrid_bd, minigrid_nb, minigrid_bd, powsybl_dir,
#          powsybl_reference, dtf_deck).
#          Not part of the library.

"""
    warmup()

Run every compile-heavy path of the foreign-formats tour once and print how
long the first solve and the rest took. The output of the paths is
discarded.
"""
function warmup()
  println("warm-up: compiles the paths of this workshop once; how long it takes depends on the machine (a Colab session is several times slower than a desktop)")
  wnet = Net(name = "warmup", baseMVA = 100.0)
  addBus!(net = wnet, busName = "A", vn_kV = 110.0)
  addBus!(net = wnet, busName = "B", vn_kV = 110.0)
  addProsumer!(net = wnet, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = wnet, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
  addPIModelACLine!(net = wnet, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  t_first = @elapsed runpf!(wnet, 10, 1e-8, 0; islands_enabled = true)
  t_second = @elapsed runpf!(wnet, 10, 1e-8, 0; islands_enabled = true)
  calcNetLosses!(wnet)
  println("warm-up solve: first ", round(t_first; digits = 2), " s (compiles), second ", round(t_second * 1000; digits = 2), " ms")

  ## every further path of the chapters once, output discarded
  t_paths = @elapsed redirect_stdout(devnull) do
    print(summarizeCGMES(path = [microgrid_be, microgrid_bd]))
    print(analyzeCGMES(path = microgrid_be))
    wres = importCGMES(path = [microgrid_be, microgrid_bd], name = "warmup_microgrid")
    wite, werg = runpf!(wres.net, 30, 1e-8, 0; islands_enabled = true)
    calcNetLosses!(wres.net)
    printACPFlowResults(wres.net, 0.0, wite, 1e-8)
    compareWithSV(wres)
    wdir = mktempdir()
    writeCGMESFiles(wres.net; path = wdir)
    importCGMES(path = wdir, name = "warmup_roundtrip")
    wfiles = [f for f in readdir(minigrid_nb; join = true) if endswith(f, ".xml") && !occursin("_TP", basename(f)) && !occursin("_SV", basename(f))]
    wnotp = importCGMES(path = vcat(wfiles, minigrid_bd), name = "warmup_no_tp")
    ## the flat-start solve of the two node-breaker twins and their comparison
    wnotp.net.flatstart = true
    runpf!(wnotp.net, 40, 1e-8, 0; islands_enabled = true)
    wva = sort([n._vm_pu for n in wnotp.net.nodeVec])
    println("max |Vm| difference: ", maximum(abs.(wva .- wva)))
    wtables = Sparlectra.read_iidm_tables(joinpath(powsybl_dir, "micro_grid_be.xiidm"))
    wpnet, wreport = build_net_from_powsybl(wtables, PowsyblAdapterOptions(remote_regulation = :remote))
    println(format_powsybl_report(wreport))
    run_sparlectra(net = wpnet, config = SparlectraConfig(
      powerflow = PowerFlowConfig(max_iter = 60, tol = 1e-9, islands_enabled = true, qlimits = Sparlectra.QLimitConfig(hysteresis_pu = 1e-6), distributed_slack = Sparlectra.DistributedSlackConfig(enabled = true, p_mode = :imported)),
      output = OutputConfig(logfile_results = :off, startup_latency_hint = false),
      control = ControlConfig(),
    ))
    ## the OpenLoadFlow reference file of the comparison
    wref, whdr = Sparlectra.DelimitedFiles.readdlm(powsybl_reference, ',', String; quotes = true, header = true)
    wcols = Dict(Symbol(strip(h)) => k for (k, h) in enumerate(vec(whdr)))
    parse(Float64, wref[1, wcols[:v_pu]])
    wcase = Sparlectra.DTFImporter.read_dtf(dtf_deck)
    wdnet = Sparlectra.DTFImporter.build_net(wcase)
    runpf!(wdnet, 50, 1e-8, 0)
    if !isempty(wcase.outages)
      wm = Sparlectra.DTFImporter.find_outage_branch_indices(wcase, first(wcase.outages))
      println("outage resolves to ", wm, ": ", Sparlectra.DTFImporter.outage_match_diagnostic(wcase, first(wcase.outages), wm))
      if length(wm) == 1
        wonet = Sparlectra.DTFImporter.build_net(wcase)
        Sparlectra.DTFImporter.apply_single_branch_outage!(wonet, wm[1])
        markIsolatedBuses!(net = wonet, log = false)
        runpf!(wonet, 50, 1e-8, 0; islands_enabled = true)
        println(maximum(abs(wonet.nodeVec[i]._vm_pu - wdnet.nodeVec[i]._vm_pu) for i in eachindex(wdnet.nodeVec)))
      end
    end
  end
  println("further paths  : ", round(t_paths; digits = 2), " s (CGMES summary, analysis, import, export, IIDM, DTF); everything warm")
  return nothing
end
