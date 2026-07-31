# Copyright 2023–2026 Udo Schmitz
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

# file: examples/others/exp_ac_rescue_dc_fallback.jl
# purpose: non-convergence handling — the power_flow.rescue strategy ladder
#          recovers a failed AC solve, and power_flow.dc.fallback leaves a
#          usable DC state when AC has no solution at all.

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

"""
    main()

Two runs on a two-bus feeder. First a solvable case whose start state is
poisoned so the plain Newton solve fails: with `rescue = true` the solver
retries from the original start with the strategy ladder (alternate start,
autodamp, DC-seeded projection) and reports which strategy converged.
Second a genuinely infeasible case (load far beyond the transfer limit):
the ladder is exhausted and `dc.fallback = true` runs the standalone DC
power flow, so the net still carries angles and branch P flows — while the
AC status honestly stays non-converged.
"""
function main()
  print_example_banner("examples/others/exp_ac_rescue_dc_fallback.jl", "non-convergence handling: rescue strategy ladder plus standalone-DC fallback")

  build = function (; p_mw::Float64)
    rnet = Net(name = "rescue_demo", baseMVA = 100.0)
    addBus!(net = rnet, busName = "Grid", vn_kV = 110.0)
    addBus!(net = rnet, busName = "Load", vn_kV = 110.0)
    addPIModelACLine!(net = rnet, fromBus = "Grid", toBus = "Load", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
    addProsumer!(net = rnet, busName = "Grid", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Grid")
    addProsumer!(net = rnet, busName = "Load", type = "ENERGYCONSUMER", p = p_mw, q = 50.0)
    return rnet
  end

  # Solvable, but the stored start state is garbage and the budget is tight.
  rescued = build(p_mw = 300.0)
  Sparlectra.setVmVa!(node = rescued.nodeVec[2], vm_pu = 0.01, va_deg = 179.0)
  cfg = Sparlectra.PowerFlowConfig(max_iter = 5, tol = 1e-8, rescue = true, start_mode = Sparlectra.StartModeConfig(flatstart = false))
  profile = Dict{Symbol,Any}()
  ite, erg = runpf!(rescued, cfg; performance_profile = profile)
  rescue_result = (converged = erg == 0, iterations = ite, strategy = get(profile, :ac_rescue_strategy, :none))

  # No AC solution exists; the DC fallback still delivers a flow picture.
  infeasible = build(p_mw = 5000.0)
  cfg_fb = Sparlectra.PowerFlowConfig(max_iter = 10, tol = 1e-8, rescue = true, dc = Sparlectra.DcPowerFlowConfig(fallback = true))
  profile_fb = Dict{Symbol,Any}()
  _, erg_fb = runpf!(infeasible, cfg_fb; performance_profile = profile_fb)
  fallback_result = (ac_converged = erg_fb == 0, dc_fallback = get(profile_fb, :dc_fallback_applied, false), va_load_deg = infeasible.nodeVec[2]._va_deg)

  return (rescue = rescue_result, fallback = fallback_result)
end

result = run_example(main)
println()
println("rescue:   converged=", result.rescue.converged, " via strategy '", result.rescue.strategy, "' after ", result.rescue.iterations, " iteration(s)")
println("fallback: AC converged=", result.fallback.ac_converged, ", DC fallback applied=", result.fallback.dc_fallback, ", Load angle ", round(result.fallback.va_load_deg; digits = 2), "°")
