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
#
# file: test/test_island_diagnostics.jl
# purpose: regression tests for AC island diagnostics reporting -- the failure
#          message must show the failing island's own iteration count and
#          stage, and islands that were never attempted must not inherit the
#          failed island's solve statistics in ac_island_solver_summary.csv.
#          Repro source: a 6209-bus CGMES delivery (1 main island + 158
#          single-bus islands) where the run-level message claimed
#          "iterations=0 / stage=before_nr" while ac_island_1_solver.log
#          recorded 80 Newton iterations, and every single-bus island row
#          duplicated island 1's mismatch/iteration/switching values.

using Test
using Sparlectra

# Per-island status as the rectangular solver stores it in
# performance_profile[:ac_island_solver_statuses] for a failed island. Values
# mirror the observed CGMES run (80 NR iterations, active-set instability).
function _island_diag_failed_status()
  return (
    island_id = 1,
    status = :not_converged,
    final_converged = false,
    numerical_converged = false,
    iterations = 80,
    initial_mismatch = 2220.614510584437,
    final_mismatch = 240324.3267233851,
    best_mismatch = 731.175415721738,
    reason = :type_identification_divergence,
    stage = :newton_iteration,
    pv_pq_switching_events = 435,
    qlimit_active_set_changes = 13,
    mismatch_history = [2220.614510584437, 731.175415721738, 240324.3267233851],
    exception_type = "",
    exception_message = "",
    stacktrace_top = "",
  )
end

# Two-bus island with one branch plus one isolated bus. The isolated bus forms
# a second diagnostics island that the solver never attempts -- the miniature
# version of the CGMES delivery's 158 single-bus islands.
function _island_diag_testnet()
  net = Net(name = "island_diag", baseMVA = 100.0)
  for b in ("A1", "A2", "B1")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.02, status = 1)
  addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "A1")
  addProsumer!(net = net, busName = "A2", type = "ENERGYCONSUMER", p = 50.0, q = 10.0)
  return net
end

# Split a summary CSV into a header-name => column-index map plus data rows
# (delimiter from the header line, see run_csv).
function _island_diag_read_summary(path::AbstractString)
  csv = run_csv(path)
  col = Dict(name => i for (i, name) in enumerate(csv.header))
  return col, csv.rows
end

function run_island_diagnostics_tests()

  @testset "no-slack failure names the isolation cause" begin (function ()
    # A branch-less delivery (e.g. one side of a DC border crossing) leaves
    # every bus isolated; the old error was a bare "no slack bus found".
    # The message must now name the actual cause: nothing is energized.
    onebus = Net(name = "island_diag_onebus", baseMVA = 100.0)
    addBus!(net = onebus, busName = "B1", vn_kV = 110.0)
    addProsumer!(net = onebus, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
    # CGMES import marks isolation before solving; mirror that state here.
    markIsolatedBuses!(net = onebus)
    err = @test_throws ErrorException runpf!(onebus, 30, 1e-8, 0)
    @test occursin("all 1 bus(es) of this network are isolated", sprint(showerror, err.value))

    # A net without any registered voltage reference keeps the plain
    # "no slack" wording, now with a hint on how to register one.
    noslack = Net(name = "island_diag_noslack", baseMVA = 100.0)
    addBus!(net = noslack, busName = "C1", vn_kV = 110.0)
    addBus!(net = noslack, busName = "C2", vn_kV = 110.0)
    addPIModelACLine!(net = noslack, fromBus = "C1", toBus = "C2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
    addProsumer!(net = noslack, busName = "C2", type = "ENERGYCONSUMER", p = 10.0, q = 2.0)
    err2 = @test_throws ErrorException runpf!(noslack, 30, 1e-8, 0)
    @test occursin("no slack bus registered", sprint(showerror, err2.value))
  end)() end

  @testset "auto slack promotes the strongest candidate" begin (function ()
    # Build a solvable two-bus net whose generator is not marked as slack.
    # NOTE: private inner name — an anonymous fixture assigning to a name
    # that also exists in the enclosing testset would rebind that local.
    mkautonet = function ()
      anet = Net(name = "island_diag_autoslack", baseMVA = 100.0)
      addBus!(net = anet, busName = "G1", vn_kV = 110.0)
      addBus!(net = anet, busName = "L1", vn_kV = 110.0)
      addPIModelACLine!(net = anet, fromBus = "G1", toBus = "L1", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addProsumer!(net = anet, busName = "G1", type = "GENERATOR", p = 20.0, q = 5.0, pMax = 150.0)
      addProsumer!(net = anet, busName = "L1", type = "ENERGYCONSUMER", p = 10.0, q = 2.0)
      return anet
    end

    # Keyword path: promotion happens, run converges, slack is registered.
    auto1 = mkautonet()
    erg = last(runpf!(auto1, 30, 1e-8, 0; auto_slack = true))
    @test erg == 0
    @test auto1.slackVec == [1]
    @test getNodeType(auto1.nodeVec[1]) == Sparlectra.Slack

    # Config path: power_flow.auto_slack reaches the solver.
    auto2 = mkautonet()
    cfg = Sparlectra.PowerFlowConfig(auto_slack = true)
    erg2 = last(runpf!(auto2, cfg))
    @test erg2 == 0
    @test auto2.slackVec == [1]

    # An external network injection outranks a larger generator.
    # addProsumer! collapses the type string to Injection/Consumption, so the
    # precise CGMES class marker is set directly, as the importer does.
    auto3 = mkautonet()
    addBus!(net = auto3, busName = "X1", vn_kV = 110.0)
    addPIModelACLine!(net = auto3, fromBus = "L1", toBus = "X1", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
    addProsumer!(net = auto3, busName = "X1", type = "EXTERNALNETWORKINJECTION", p = 0.0, q = 0.0, pMax = 50.0)
    auto3.prosumpsVec[end].comp.cTyp = Sparlectra.ExternalNetworkInjection
    picked = ensureSlack!(auto3; log = false)
    @test picked == 3

    # A registered slack is left alone: ensureSlack! is a no-op.
    auto4 = mkautonet()
    addProsumer!(net = auto4, busName = "G1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "G1")
    @test ensureSlack!(auto4; log = false) === nothing
  end)() end

  @testset "an island without unknowns: a lone generator and a lone load" begin (function ()
    # A PV generator left alone on its AC island (its only line out) is
    # the island's reference and nothing else: no unknown, Newton has
    # nothing to solve. The loop used to reduce over the empty mismatch
    # vector and the whole run failed with an ArgumentError (found on
    # case_ACTIVSg2000, a generator behind its outaged unit transformer).
    # Power mode builds the light final diagnostics, which took the worst
    # row of the same empty vector (argmax over nothing); both paths run.
    for pm in (false, true)
      net = Net(name = "pv_alone", baseMVA = 100.0)
      for b in ("S", "L", "G")
        addBus!(net = net, busName = b, vn_kV = 110.0)
      end
      addPIModelACLine!(net = net, fromBus = "S", toBus = "L", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = net, fromBus = "L", toBus = "G", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addProsumer!(net = net, busName = "S", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "S")
      addProsumer!(net = net, busName = "L", type = "ENERGYCONSUMER", p = 30.0, q = 10.0)
      addProsumer!(net = net, busName = "G", type = "SYNCHRONOUSMACHINE", p = 20.0, vm_pu = 1.02, va_deg = 0.0)
      Sparlectra.setBranchStatus!(net.branchVec[2], false)
      prof = Dict{Symbol,Any}()
      _, erg = redirect_stdout(devnull) do
        runpf!(net, 30, 1e-8, 0; islands_enabled = true, power_mode = pm, performance_profile = prof)
      end
      @test erg == 0
      lone = [st for st in values(prof[:ac_island_solver_statuses]) if get(st, :no_unknowns, false) === true]
      @test length(lone) == 1
      @test lone[1].status === :converged
      @test lone[1].iterations == 0
      @test net.nodeVec[3]._vm_pu == 1.02
    end

    # The N-1 engine on the shipped sp_case60: an outage that cuts off a
    # generator bus and one that cuts off a pure load bus, both without
    # error and with their own note (the generator out of service, the
    # load not supplied; a bus without infeed is never promoted)
    net60 = Sparlectra.importSCF(joinpath(dirname(@__DIR__), "data", "scf", "sp_case60.scf.json"))
    cases = [c for c in generateN1Branches(net60) if c.name in ("Ostsee_Westhof", "Altfeld_Neudorf")]
    @test length(cases) == 2
    res = Dict(r.name => r for r in redirect_stdout(() -> runContingencies!(net60, cases), devnull))
    gen_out = res["Ostsee_Westhof"]
    load_out = res["Altfeld_Neudorf"]
    @test gen_out.converged && load_out.converged
    @test occursin("generating unit(s)", something(gen_out.error, "")) && gen_out.shed_load_mw == 0.0
    @test occursin("load disconnected", something(load_out.error, "")) && load_out.shed_load_mw > 0.0
  end)() end

end
