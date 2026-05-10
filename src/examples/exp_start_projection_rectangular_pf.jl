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

using Sparlectra

"""
    main()

Demonstrate APSLF-style start projection before the rectangular Newton power
flow. The solver sanitizes the raw seed, tries a DC-angle start, and scans
convex blends between both starts before entering Newton iterations.
"""
function main()
  net = Net(name = "start_projection_demo", baseMVA = 100.0)

  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)
  addACLine!(net = net, fromBus = "Slack", toBus = "Load", length = 1.0, r = 0.01, x = 0.10)
  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = 80.0, q = 25.0)

  iterations, status, elapsed = run_net_acpflow(
    net = net,
    max_ite = 30,
    tol = 1e-8,
    verbose = 1,
    method = :rectangular,
    opt_flatstart = true,
    start_projection = true,
    start_projection_try_dc_start = true,
    start_projection_try_blend_scan = true,
    start_projection_blend_lambdas = [0.25, 0.5, 0.75],
    start_projection_dc_angle_limit_deg = 60.0,
    show_results = true,
  )

  println("start-projected rectangular PF: status=$(status), iterations=$(iterations), elapsed=$(elapsed) s")
  return status == 0
end

Base.invokelatest(getfield(@__MODULE__, :main))
