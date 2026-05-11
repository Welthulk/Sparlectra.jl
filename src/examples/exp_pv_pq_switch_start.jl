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

function _build_demo_net()
  net = Net(name = "pv_pq_switch_start_demo", baseMVA = 100.0)
  addBus!(net = net, busName = "SL", vn_kV = 110.0)
  addBus!(net = net, busName = "PV", vn_kV = 110.0)
  addBus!(net = net, busName = "PQ", vn_kV = 110.0)
  addACLine!(net = net, fromBus = "SL", toBus = "PV", length = 8.0, r = 0.08, x = 0.20, c_nf_per_km = 8.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "PV", toBus = "PQ", length = 8.0, r = 0.10, x = 0.25, c_nf_per_km = 8.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "SL", toBus = "PQ", length = 12.0, r = 0.12, x = 0.30, c_nf_per_km = 8.0, tanδ = 0.0)
  addProsumer!(net = net, busName = "SL", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "SL")
  addProsumer!(net = net, busName = "PV", type = "SYNCHRONOUSMACHINE", p = 70.0, q = 33.2, vm_pu = 1.025, qMin = -15.0, qMax = 15.0)
  addProsumer!(net = net, busName = "PQ", type = "ENERGYCONSUMER", p = 100.0, q = 30.0)
  return net
end

function main()
  net_iter = _build_demo_net()
  _, status_iter = runpf!(net_iter, 30, 1e-8, 1; method = :rectangular, qlimit_start_iter = 4)
  println("iteration-start status: ", status_iter, ", Q-limit events: ", length(net_iter.qLimitLog))

  net_auto = _build_demo_net()
  _, status_auto = runpf!(net_auto, 30, 1e-8, 1; method = :rectangular, qlimit_start_mode = :auto_q_delta, qlimit_auto_q_delta_pu = 0.1)
  println("auto-Q-delta status: ", status_auto, ", Q-limit events: ", length(net_auto.qLimitLog))
end

Base.invokelatest(main)
