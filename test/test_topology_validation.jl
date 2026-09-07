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

# file: test/test_topology_validation.jl
# purpose: topology validation (0.10.0): stage-1 pre-check matrix
#          incl. false-positive guards, the advisory contract, the stage-2
#          fingerprint classification, the stage-3 hypothesis test with an
#          untouched input net, and the singular-normal-equations
#          regression from parallel released taps.

using Sparlectra
using Test
using Logging
using LinearAlgebra

# meshed two-level net plus a linked busbar section: every stage-1 check
# has something to bite on (branches, a link, complete injection coverage)
function _topo_net(; link_status::Int = 1)
  net = Net(name = "topo_val", baseMVA = 100.0)
  for (b, vn) in (("H1", 110.0), ("H2", 110.0), ("B3", 110.0), ("L1", 20.0), ("L2", 20.0))
    addBus!(net = net, busName = b, vn_kV = vn)
  end
  addProsumer!(net = net, busName = "H1", type = "EXTERNALNETWORKINJECTION", referencePri = "H1", vm_pu = 1.02, va_deg = 0.0)
  addProsumer!(net = net, busName = "L1", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
  addProsumer!(net = net, busName = "L2", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
  addPIModelACLine!(net = net, fromBus = "H1", toBus = "H2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "L1", toBus = "L2", r_pu = 0.02, x_pu = 0.1, b_pu = 0.0, status = 1)
  addPIModelTrafo!(net = net, fromBus = "H1", toBus = "L1", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelTrafo!(net = net, fromBus = "H2", toBus = "L2", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addLink!(net = net, fromBus = "H2", toBus = "B3", status = link_status)
  for br in net.branchVec
    br.ratio == 0.0 && continue
    br.has_ratio_tap = true
    br.tap_step = 0.00625
  end
  return net
end

function _topo_measurements(net)
  ite, erg = runpf!(net, 40, 1e-12, 0; method = :rectangular)
  @test erg == 0
  return generateMeasurementsFromPF(net; includeImag = true, noise = false)
end

_meas_with(m::Measurement; value = m.value, active = m.active) = Measurement(typ = m.typ, value = value, sigma = m.sigma, active = active, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id, linkIdx = m.linkIdx)

function test_topology_stage1()::Bool
  @testset "Topology stage 1 precheck matrix" begin
    # clean reference: a consistent set fires NO finding, and the summary
    # says explicitly what was checked (the always-logged line)
    net = _topo_net()
    meas = _topo_measurements(net)
    pre = validate_topology(net, meas)
    @test isempty(pre.findings)
    @test occursin("no findings", pre.summary)
    @test pre.n_checked_branches == 4

    # 1) open element with measured flow: the status contradiction
    net = _topo_net()
    meas = _topo_measurements(net)
    setBranchStatus!(net.branchVec[2], false)
    pre = validate_topology(net, meas)
    @test any(f -> f.kind == :open_element_with_flow && occursin("branch 2", f.location), pre.findings)
    # the estimation still runs (advisory) and carries the findings
    res = runse!(net, meas; maxIte = 40, tol = 1e-8)
    @test res.topologyFindings !== nothing
    @test any(f -> f.kind == :open_element_with_flow, res.topologyFindings)

    # 2) closed branch reading dead while its neighbourhood is loaded
    net = _topo_net()
    meas = _topo_measurements(net)
    dead = [m.branchIdx == 2 && m.typ in (Sparlectra.PflowMeas, Sparlectra.QflowMeas, Sparlectra.ImagMeas) ? _meas_with(m; value = 0.001) : m for m in meas]
    pre = validate_topology(net, dead)
    @test any(f -> f.kind == :closed_element_without_flow && occursin("branch 2", f.location), pre.findings)
    @test all(f -> f.severity in (:warning,), [f for f in pre.findings if f.kind == :closed_element_without_flow])

    # false-positive guard: when the WHOLE neighbourhood reads dead the
    # branch is legitimately unloaded, no finding
    alldead = [m.typ in (Sparlectra.PflowMeas, Sparlectra.QflowMeas, Sparlectra.ImagMeas) ? _meas_with(m; value = 0.0001) : m for m in meas]
    pre = validate_topology(net, alldead)
    @test !any(f -> f.kind == :closed_element_without_flow, pre.findings)

    # 3) closed link with disagreeing voltage measurements on both sides
    net = _topo_net()
    meas = _topo_measurements(net)
    b3 = 3   # bus B3
    mism = [m.typ == Sparlectra.VmMeas && m.busIdx == b3 ? _meas_with(m; value = m.value + 0.05) : m for m in meas]
    pre = validate_topology(net, mism)
    @test any(f -> f.kind == :closed_link_voltage_mismatch && occursin("link 1", f.location), pre.findings)
    # an OPEN link may disagree freely (real separation): no finding
    net = _topo_net()
    measO = _topo_measurements(net)
    net.linkVec[1].status = 0
    mismO = [m.typ == Sparlectra.VmMeas && m.busIdx == b3 ? _meas_with(m; value = m.value + 0.05) : m for m in measO]
    pre = validate_topology(net, mismO)
    @test !any(f -> f.kind == :closed_link_voltage_mismatch, pre.findings)

    # 4) node balance at a completely measured node
    net = _topo_net()
    meas = _topo_measurements(net)
    l2 = 5   # bus L2
    kcl = [m.typ == Sparlectra.PinjMeas && m.busIdx == l2 ? _meas_with(m; value = m.value + 20.0) : m for m in meas]
    pre = validate_topology(net, kcl)
    @test any(f -> f.kind == :kcl_violation && occursin("L2", f.location), pre.findings)
    # partially measured node: drop one flow row at L2, no verdict
    part = [m.typ == Sparlectra.PflowMeas && m.branchIdx == 2 && m.direction == :to ? _meas_with(m; active = false) : m for m in kcl]
    pre = validate_topology(net, part)
    @test !any(f -> f.kind == :kcl_violation, pre.findings)
  end

  return true
end

function test_topology_stage2_fingerprint()::Bool
  @testset "Topology stage 2 fingerprint classification" begin
    # truth: line L1-L2 is OPEN; the model believes it closed. The
    # measurements come from the true (open) state, so the estimator can
    # neither absorb nor eliminate its way out: elimination exhausts,
    # the band stays :high, and the suspects cluster at the line's stations.
    # the open element is the TRANSFORMER H2-L2: with it gone the whole
    # low-voltage load flows through H1-L1 and the measured state is far
    # from anything the closed model can explain (an open branch whose
    # terminals sit at similar voltages anyway would be undetectable)
    tnet = _topo_net()
    setBranchStatus!(tnet.branchVec[4], false)
    meas = _topo_measurements(tnet)
    net = _topo_net()   # model: branch 4 closed
    diag = runse_diagnostics(net, meas; max_eliminations = 2, maxIte = 100, tol = 1e-8)
    @test diag.stop_reason == :max_eliminations
    @test diag.topology_findings !== nothing
    f = first(something(diag.topology_findings, NamedTuple[]))
    @test f.kind == :topology_error_suspected_at_station
    @test occursin("H2", f.location) || occursin("L2", f.location)
    @test f.severity == :strong

    # control: a curable gross error NEVER produces a topology finding,
    # the ordinary elimination handles it (full-fingerprint requirement)
    net = _topo_net()
    meas = _topo_measurements(net)
    gi = findfirst(m -> m.typ == Sparlectra.PflowMeas, meas)
    meas[gi] = _meas_with(meas[gi]; value = meas[gi].value + 10.0 * meas[gi].sigma)
    diag = runse_diagnostics(net, meas; max_eliminations = 3, maxIte = 40, tol = 1e-8)
    @test diag.stop_reason != :max_eliminations
    @test diag.topology_findings === nothing
  end

  return true
end

function test_topology_stage3_hypotheses()::Bool
  @testset "Topology stage 3 hypothesis test" begin
    tnet = _topo_net()
    setBranchStatus!(tnet.branchVec[4], false)
    meas = _topo_measurements(tnet)
    net = _topo_net()   # model: branch 4 closed, truth open

    # input snapshot for the bitwise-untouched assertion
    stat0 = [(br.status, br.from_status, br.to_status, br.tap_ratio) for br in net.branchVec]
    link0 = [l.status for l in net.linkVec]
    vm0 = [nd._vm_pu for nd in net.nodeVec]

    rep = test_topology_hypotheses(net, meas; max_candidates = 8, maxIte = 100, tol = 1e-8)
    @test rep.base_verdict == :high
    @test rep.n_candidates >= 1
    top = first(rep.recommendations)
    @test top.verdict == :hypothesis_supported
    @test top.kind == :branch && top.idx == 4
    @test top.hypothesis == "actually open"
    @test top.j_after < top.j_before
    @test top.elapsed_s >= 0.0
    # the small net is honestly ambiguous (opening the quiet H1-H2 line
    # also lands in the band); the report says so instead of hiding it,
    # and the J ranking still puts the true hypothesis first
    @test rep.ambiguous == true
    # ranked recommendations only, nothing switched: the input is untouched
    @test [(br.status, br.from_status, br.to_status, br.tap_ratio) for br in net.branchVec] == stat0
    @test [l.status for l in net.linkVec] == link0
    @test [nd._vm_pu for nd in net.nodeVec] == vm0

    # island case: the toggled candidate sits in island b of a two-island
    # net; the run works island-wise and the input again stays untouched
    function _two_island()
      n = Net(name = "topo_islands", baseMVA = 100.0)
      for sfx in ("a", "b")
        for (b, vn) in (("H1", 110.0), ("H2", 110.0), ("L1", 20.0), ("L2", 20.0))
          addBus!(net = n, busName = b * sfx, vn_kV = vn)
        end
        addProsumer!(net = n, busName = "H1" * sfx, type = "EXTERNALNETWORKINJECTION", referencePri = "H1" * sfx, vm_pu = 1.02, va_deg = 0.0)
        addProsumer!(net = n, busName = "L1" * sfx, type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
        addProsumer!(net = n, busName = "L2" * sfx, type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
        addPIModelACLine!(net = n, fromBus = "H1" * sfx, toBus = "H2" * sfx, r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = n, fromBus = "L1" * sfx, toBus = "L2" * sfx, r_pu = 0.02, x_pu = 0.1, b_pu = 0.0, status = 1)
        addPIModelTrafo!(net = n, fromBus = "H1" * sfx, toBus = "L1" * sfx, r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
        addPIModelTrafo!(net = n, fromBus = "H2" * sfx, toBus = "L2" * sfx, r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
      end
      return n
    end
    tn2 = _two_island()
    setBranchStatus!(tn2.branchVec[8], false)   # trafo H2b-L2b, island b
    runpf!(tn2, 40, 1e-12, 0; islands_enabled = true)
    meas2 = generateMeasurementsFromPF(tn2; includeImag = true, noise = false)
    net2 = _two_island()
    stat2 = [(br.status, br.from_status, br.to_status) for br in net2.branchVec]
    rep2 = test_topology_hypotheses(net2, meas2; candidates = [(kind = :branch, idx = 8)], maxIte = 100, tol = 1e-8)
    @test rep2.n_candidates == 1
    @test first(rep2.recommendations).verdict == :hypothesis_supported
    @test [(br.status, br.from_status, br.to_status) for br in net2.branchVec] == stat2
  end

  return true
end

function test_topology_singular_normal_equations()::Bool
  @testset "Singular WLS normal equations regression" begin
    # solve_linear contract: with allow_pinv and a raised svd_max_n a
    # singular dense system degrades to a least-squares solution instead of
    # throwing (the service crash on collectively dependent tap columns);
    # the small default cap still protects the power-flow callers
    A = zeros(100, 100)
    for i = 1:99
      A[i, i] = 1.0
    end
    A[:, 100] = A[:, 99]   # exact column duplicate: singular
    b = ones(100)
    x = Sparlectra.solve_linear(A, b; allow_pinv = true, svd_max_n = 20_000)
    @test all(isfinite, x)
    @test_throws Exception Sparlectra.solve_linear(A, b; allow_pinv = true)

    # SE level: two PARALLEL released taps between the same buses pass the
    # per-column observability test individually but are collectively
    # dependent; the run must converge on a least-squares split instead of
    # crashing (the sum of both steps carries the information)
    function _par_net()
      n = Net(name = "topo_partap", baseMVA = 100.0)
      for (b, vn) in (("H1", 110.0), ("H2", 110.0), ("L1", 20.0), ("L2", 20.0))
        addBus!(net = n, busName = b, vn_kV = vn)
      end
      addProsumer!(net = n, busName = "H1", type = "EXTERNALNETWORKINJECTION", referencePri = "H1", vm_pu = 1.02, va_deg = 0.0)
      addProsumer!(net = n, busName = "L1", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
      addProsumer!(net = n, busName = "L2", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
      addPIModelACLine!(net = n, fromBus = "H1", toBus = "H2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = n, fromBus = "L1", toBus = "L2", r_pu = 0.02, x_pu = 0.1, b_pu = 0.0, status = 1)
      addPIModelTrafo!(net = n, fromBus = "H1", toBus = "L1", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
      addPIModelTrafo!(net = n, fromBus = "H1", toBus = "L1", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
      for br in n.branchVec
        br.ratio == 0.0 && continue
        br.has_ratio_tap = true
        br.tap_step = 0.00625
      end
      return n
    end
    tnet = _par_net()
    br = tnet.branchVec[3]
    ct = Sparlectra._cascade_tap(br, 2 * 0.00625, 0.0, 0.0)
    br.tap_ratio = ct.tap_ratio
    br.phase_shift_deg = ct.phase_shift_deg
    ite, erg = runpf!(tnet, 40, 1e-12, 0; method = :rectangular)
    @test erg == 0
    meas = generateMeasurementsFromPF(tnet; includeImag = true, noise = false)
    net = _par_net()
    setTapEstimation!(net; trafo = 3, mode = :ratio)
    setTapEstimation!(net; trafo = 4, mode = :ratio)
    res = runse!(net, meas; maxIte = 40, tol = 1e-8)
    @test res.converged
    @test res.tapEstimates !== nothing && length(res.tapEstimates) == 2
    # per-branch flow rows follow the LIVE tap (branchFlow_pu fix), so the
    # parallel split is IDENTIFIED, not a least-squares choice: trafo 3
    # carries the full deviation, its partner stays at neutral, and the
    # difference-direction prior stays out of measured pairs
    s = sum(t.electrical_step_1 for t in res.tapEstimates)
    @test abs(s - 2.0) < 0.05
    fixedByBranch = Dict(t.branch => t.fixed_step_1 for t in res.tapEstimates)
    @test fixedByBranch[3] == 2
    @test fixedByBranch[4] == 0
    # regression (maintainer run 1908605e): robust weighting on parallel
    # released taps used to excite the unobservable difference direction
    # into divergence; with the targeted prior plus the weight freeze the
    # robust run converges too
    net = _par_net()
    setTapEstimation!(net; trafo = 3, mode = :ratio)
    setTapEstimation!(net; trafo = 4, mode = :ratio)
    resR = runse!(net, meas; maxIte = 50, tol = 1e-8, robust = true)
    @test resR.converged
    @test resR.iterations < 50
  end

  return true
end

function run_topology_validation_tests()
  @testset "Topology validation" begin
    tests = [
      ("Stage 1 precheck matrix", test_topology_stage1),
      ("Stage 2 fingerprint classification", test_topology_stage2_fingerprint),
      ("Stage 3 hypothesis test", test_topology_stage3_hypotheses),
      ("Singular normal equations", test_topology_singular_normal_equations),
    ]
    for (name, testfn) in tests
      @testset "$name" begin
        @test _se_run_quiet(testfn) == true
      end
    end
  end
end
