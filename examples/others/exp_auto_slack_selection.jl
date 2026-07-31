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

# file: examples/others/exp_auto_slack_selection.jl
# purpose: automatic slack selection (power_flow.auto_slack / ensureSlack!)
#          on a network whose data registers no voltage reference.

using Sparlectra

include(joinpath(@__DIR__, "..", "internal", "example_header.jl"))

"""
    main()

A small grid with two generators and a load, none of them marked as slack —
the situation of an edited or partial data set. The first solve attempt
aborts with the no-slack error. With `auto_slack = true` the solver promotes
the strongest candidate itself (`ensureSlack!`): external network injections
outrank generators, within each group the largest unit wins (`ratedS`, then
`maxP`, then dispatch). The promotion is logged and the run converges.
"""
function main()
  print_example_banner("examples/others/exp_auto_slack_selection.jl", "automatic slack selection for cases whose data registers no voltage reference")

  build = function ()
    aunet = Net(name = "auto_slack_demo", baseMVA = 100.0)
    for bus in ("GenA", "GenB", "Load")
      addBus!(net = aunet, busName = bus, vn_kV = 110.0)
    end
    addPIModelACLine!(net = aunet, fromBus = "GenA", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
    addPIModelACLine!(net = aunet, fromBus = "GenB", toBus = "Load", r_pu = 0.02, x_pu = 0.12, b_pu = 0.01, status = 1)
    # Neither generator carries referencePri — no slack is registered.
    addProsumer!(net = aunet, busName = "GenA", type = "GENERATOR", p = 30.0, q = 5.0, pMax = 80.0)
    addProsumer!(net = aunet, busName = "GenB", type = "GENERATOR", p = 40.0, q = 5.0, pMax = 200.0)
    addProsumer!(net = aunet, busName = "Load", type = "ENERGYCONSUMER", p = 60.0, q = 15.0)
    return aunet
  end

  # Without the option the run aborts — keep the message, it names the fix.
  failed_message = try
    runpf!(build(), 30, 1e-8, 0)
    ""
  catch err
    sprint(showerror, err)
  end

  # With the option the strongest generator (GenB, 200 MW) becomes the slack.
  net = build()
  ite, erg = runpf!(net, 30, 1e-8, 0; auto_slack = true)
  slack_bus = getCompName(net.nodeVec[only(net.slackVec)].comp)
  return (failed_message = failed_message, iterations = ite, converged = erg == 0, slack_bus = slack_bus)
end

result = run_example(main)
println("without auto_slack: ", result.failed_message)
println()
println("with auto_slack:    converged=", result.converged, " after ", result.iterations, " iteration(s), slack bus: ", result.slack_bus)
