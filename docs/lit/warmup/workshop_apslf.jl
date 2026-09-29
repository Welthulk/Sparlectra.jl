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

# file: docs/lit/warmup/workshop_apslf.jl
# purpose: compile warm-up of the APSLF workshop (docs/lit/workshop_apslf.jl):
#          runs every path the parts use once on throwaway nets, so the
#          parts do not stall on first-call compilation. Included by the
#          first code cell of the notebook, which defines the helpers used
#          here (build_ring7, bus_vm, bus_va, scaled_ring7, clamped, quiet,
#          cfg_nr, cfg_apslf) and loads Printf. Not part of the library.

"""
    warmup()

Run every compile-heavy path of the APSLF workshop once and print how long
it took. The output of the paths is discarded.
"""
function warmup()
  println("warm-up: compiles the paths of this workshop once; how long it takes depends on the machine (a Colab session is several times slower than a desktop)")
  t_paths = @elapsed redirect_stdout(devnull) do
    w_nr = run_sparlectra(net=build_ring7("warmup NR"), config=cfg_nr)
    w_ap = run_sparlectra(net=scaled_ring7(1.5), config=cfg_apslf)
    println(maximum(abs.(bus_vm(w_nr.net) .- bus_vm(w_ap.net))), bus_va(w_ap.net))
    println(Sparlectra.rectangular_pf_status(w_ap.net).apslf_convergence_line)
    w_cfg = SparlectraConfig(powerflow=PowerFlowConfig(solver=:apslf, apslf=Sparlectra.ApslfConfig(order=8, convergence_radius=false)), output=quiet)
    run_sparlectra(net=build_ring7("warmup order"), config=w_cfg)
    w_hyb = SparlectraConfig(powerflow=PowerFlowConfig(solver=:rectangular, apslf_start=Sparlectra.ApslfStartConfig(enabled=true)), output=quiet)
    run_sparlectra(net=build_ring7("warmup hybrid"), config=w_hyb)
    w_case = joinpath(dirname(dirname(pathof(Sparlectra))), "data", "mpower", "sp_case9.m")
    w9 = run_sparlectra(casefile=basename(w_case), path=dirname(w_case), config=cfg_apslf)
    println(clamped(w9.net))
    w14 = joinpath(dirname(dirname(pathof(Sparlectra))), "data", "scf", "sp_case14.scf.json")
    run_sparlectra(net=importSCF(w14), config=w_hyb)
    @printf("%.2e\n", 0.0)
    ## the order table of Part 1 (low and high order, its table format)
    for w_order in (4, 40)
      w_cfg_o = SparlectraConfig(powerflow=PowerFlowConfig(solver=:apslf, apslf=Sparlectra.ApslfConfig(order=w_order)), output=quiet)
      w_r = run_sparlectra(net=build_ring7("warmup order $(w_order)"), config=w_cfg_o)
      @printf("  %3d    %.3e          %s\n", w_order, w_r.final_mismatch, w_r.final_converged)
    end
    ## the loadability sweep of Part 2: both solvers past the fold (neither
    ## converges), the NamedTuple rows and the table format
    w_rows = NamedTuple[]
    for w_lambda in (1.0, 4.0)
      w_rap = run_sparlectra(net=scaled_ring7(w_lambda), config=cfg_apslf)
      w_rnr = run_sparlectra(net=scaled_ring7(w_lambda), config=cfg_nr)
      w_st = Sparlectra.rectangular_pf_status(w_rap.net)
      w_vmin = w_rap.final_converged ? minimum(bus_vm(w_rap.net)) : NaN
      push!(w_rows, (lambda=w_lambda, dmin=w_st.apslf_convergence_radius, level=w_st.apslf_convergence_level, apslf=w_rap.final_converged, nr=w_rnr.final_converged, vmin=w_vmin))
    end
    for w_row in w_rows
      @printf("  %.1f   %7.3f   %-5s   %-5s   %-5s  %.3f\n", w_row.lambda, w_row.dmin, w_row.level, w_row.apslf, w_row.nr, w_row.vmin)
    end
  end
  println("warm: both solvers, series options, loadability sweep, hybrid start, case import ", round(t_paths; digits = 2), " s (first calls compile); everything warm")
  return nothing
end
