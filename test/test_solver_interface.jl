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

# file: test/test_solver_interface.jl

using Test
using Sparlectra



"""
    test_external_solver_interface() -> Bool

Minimal unit test for the external solver interface:
- PFModel construction
- canonical mismatch evaluation
- solution write-back into Net

Returns `true` if all checks pass.
"""
function test_external_solver_interface()::Bool

  @testset "External solver interface (PFModel/PFSolution)" begin

    # Existing 3-bus test network from testgrid.jl
    net = createTest3BusNet()

    # Build canonical PF model
    model = buildPfModel(
      net;
      opt_sparse = true,
      flatstart = true,
      include_limits = false,
      verbose = 0,
    )

    # --- Structural checks ---------------------------------------------------
    @test size(model.Ybus, 1) == length(model.busIdx_net)
    @test size(model.Ybus, 2) == length(model.busIdx_net)

    @test length(model.busIdx_net) == 3
    @test model.busIdx_net == [1, 2, 3]
    @test model.busType == [:PQ, :PV, :Slack]
    @test model.slack_idx == 3

    # --- Initial state -------------------------------------------------------
    @test length(model.V0) == 3
    @test all(isfinite.(abs.(model.V0)))

    # --- Canonical mismatch --------------------------------------------------
    r0 = mismatchInf(model, model.V0)
    @test isfinite(r0)
    @test r0 ≥ 0.0

    # --- Apply solution back into Net ---------------------------------------
    sol = PFSolution(
      V = copy(model.V0),
      converged = true,
      iters = 0,
      residual_inf = r0,
    )

    applyPfSolution!(
      net,
      model,
      sol;
      write_pq_results = false,
      verbose = 0,
    )

    for (k, busIdx) in enumerate(model.busIdx_net)
      node = net.nodeVec[busIdx]
      @test isapprox(node._vm_pu, abs(model.V0[k]); atol = 1e-12)
      @test isapprox(node._va_deg, rad2deg(angle(model.V0[k])); atol = 1e-10)
    end

  end

  return true
end

function run_solver_interface_tests()
  # Covers solver-interface integration and option behavior:
  # external model API, Q-limit reporting/autocorrection, PV->PQ locking, and final-limit reporting.
  @testset "Solver interface" begin
    @test test_external_solver_interface() == true
    # Checks tabular Q-limit output truncation plus sign-validation/autocorrect behavior.
    @testset "Q-limit reporting and validation options" begin
      net = createTest3BusNet()
      net.nodeVec[1]._nodeType = Sparlectra.PV
      net.nodeVec[2]._nodeType = Sparlectra.PV
      net.nodeVec[3]._nodeType = Sparlectra.PV
      empty!(net.qmin_pu)
      append!(net.qmin_pu, [-0.1, 0.05, -0.2])
      empty!(net.qmax_pu)
      append!(net.qmax_pu, [0.2, 0.1, 0.3])

      io = IOBuffer()
      printPVQLimitsTable(net; io = io, max_rows = 2)
      printed = String(take!(io))
      @test occursin("more PV rows omitted", printed)

      res = validate_q_limit_signs!(net.qmin_pu, net.qmax_pu; io = io, autocorrect = true, warn = false)
      @test res.flagged >= 1
      @test res.corrected >= 1
      @test net.qmin_pu[2] <= 0.0
      @test net.qmax_pu[2] >= 0.0
    end

    # Verifies that configured lock_pv_to_pq_buses keeps selected PV buses from switching to PQ.
    @testset "PV->PQ lock option" begin
      net_unlocked = createTest3BusNet()
      setQLimits!(net = net_unlocked, qmin_MVar = -1.0, qmax_MVar = 1.0, busName = "STATION1")
      _, erg_unlocked = runpf!(net_unlocked, 20, 1e-6, 0; method = :rectangular)
      @test erg_unlocked == 0
      @test getNodeType(net_unlocked.nodeVec[2]) == Sparlectra.PQ

      net_locked = createTest3BusNet()
      setQLimits!(net = net_locked, qmin_MVar = -1.0, qmax_MVar = 1.0, busName = "STATION1")
      _, erg_locked = runpf!(net_locked, 20, 1e-6, 0; method = :rectangular, lock_pv_to_pq_buses = [2])
      @test erg_locked == 0
      @test getNodeType(net_locked.nodeVec[2]) == Sparlectra.PV
    end

    # Ensures final-limit validation remains robust when q-generation data is partially missing.
    @testset "Final limit validation tolerates missing qgen" begin
      net = createTest3BusNet()
      net.nodeVec[3]._qƩGen = nothing
      io = IOBuffer()
      res = printFinalLimitValidation(net; q_headroom = 0.20, io = io)
      out = String(take!(io))
      @test haskey(res, :q_violations)
      @test haskey(res, :v_violations)
      @test occursin("Final limit validation:", out)
    end
  end
end
