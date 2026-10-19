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

# file: test/test_legacy_voltage_dependent_control.jl
# purpose: legacy profile: the testsets of test/test_voltage_dependent_control.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_voltage_dependent_control.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_voltage_dependent_control_tests()
  # cut from testset: Voltage dependent P(U)/Q(U) control
    # A Q(U) unit makes the solver rebuild S every iteration; the rebuild
    # used to drop the Q-limit clamp, so the regulating unit ran at its
    # scheduled 3 MVAr instead of its 5 MVAr limit. Q(U) unit on the clamped
    # bus itself and on another bus, and the classic outer loop, which
    # clamps by fixing the unit's qVal. The clamped bus reports its physical
    # generator Q (injection plus load), not the band value at 1.04 pu.
    # Last case: the bus is on lock_pv_to_pq_buses with a schedule of
    # -20 MVAr inside the band; the regulating unit keeps its schedule (it
    # ran at -10 when the lock counted the Q(U) unit with its qVal of 2 MVAr
    # instead of its -8 MVAr at 1.04 pu).
    @testset "Q-limit clamp survives the Q(U) rebuild" begin (function ()
      for (qu_bus, mode, lock, q_reg, q_reg_expected) in (("PV", :active_set, false, 3.0, 5.0), ("Q", :active_set, false, 3.0, 5.0), ("PV", :classic_simultaneous, false, 3.0, 5.0), ("PV", :active_set, true, -20.0, -20.0))
        net = Net(name = "qu_clamp_$(qu_bus)_$(mode)_$(lock)", baseMVA = 100.0)
        addBus!(net = net, busName = "Slack", vn_kV = 110.0)
        addBus!(net = net, busName = "PV", vn_kV = 110.0)
        addBus!(net = net, busName = "Q", vn_kV = 110.0)
        addPIModelACLine!(net = net, fromBus = "Slack", toBus = "PV", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = net, fromBus = "Slack", toBus = "Q", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
        addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "Slack")
        addProsumer!(net = net, busName = "PV", type = "SYNCHRONOUSMACHINE", p = 20.0, q = q_reg, vm_pu = 1.04, isRegulated = true, qMin = -30.0, qMax = 5.0)
        addProsumer!(net = net, busName = "PV", type = "ENERGYCONSUMER", p = 30.0, q = 20.0)
        qu = QUController(make_characteristic([(0.95, 0.10), (1.05, -0.10)]); qmin_MVAr = -50.0, qmax_MVAr = 50.0, sbase_MVA = 100.0)
        addProsumer!(net = net, busName = qu_bus, type = "GENERATOR", p = 5.0, q = 2.0, qMin = -50.0, qMax = 50.0, qu_controller = qu)

        b = geNetBusIdx(net = net, busName = "PV")
        ite, erg = runpf!(net, 30, 1e-10, 0; method = :rectangular, qlimit_enforcement_mode = mode, lock_pv_to_pq_buses = lock ? [b] : Int[])
        @test erg == 0
        @test haskey(net.qLimitEvents, b)
        V = buildVoltageVector(net)
        qgen_phys = imag(Sparlectra.calc_injections(Sparlectra.createYBUS(net = net), V)[b]) * net.baseMVA + 20.0
        q_qu = qu_bus == "PV" ? first(evaluate_controller(qu, abs(V[b]))) * net.baseMVA : 0.0
        # the regulating unit sits at its limit, or at its schedule when locked
        @test isapprox(qgen_phys - q_qu, q_reg_expected; atol = 1e-6)
        @test isapprox(net.nodeVec[b]._qƩGen, qgen_phys; atol = 1e-6)
      end
    end)() end

  return true
end
