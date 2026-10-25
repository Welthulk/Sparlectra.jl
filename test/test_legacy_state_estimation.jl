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

# file: test/test_legacy_state_estimation.jl
# purpose: legacy profile: the testsets of test/test_state_estimation.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_state_estimation.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function test_state_estimation_link_allocation_w2()::Bool
  # W2 allocation. Measurement-free runs equal the KCL
  # allocation bitwise (shared core); a link measurement steers the split of
  # a zero-impedance ring while the component KCL still holds.
  @testset "State estimation W2 link allocation" begin (function ()
    # ring of three link-joined buses, fed and loaded through real branches
    net = Net(name = "w2_ring", baseMVA = 100.0)
    for b in ("S", "R1", "R2", "R3", "L2", "L3")
      addBus!(net = net, busName = b, vn_kV = 110.0)
    end
    addPIModelACLine!(net = net, fromBus = "S", toBus = "R1", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = net, fromBus = "R2", toBus = "L2", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = net, fromBus = "R3", toBus = "L3", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
    addProsumer!(net = net, busName = "S", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "S")
    addProsumer!(net = net, busName = "L2", type = "ENERGYCONSUMER", p = 12.0, q = 4.0)
    addProsumer!(net = net, busName = "L3", type = "ENERGYCONSUMER", p = 18.0, q = 5.0)
    addLink!(net = net, fromBus = "R1", toBus = "R2", status = 1)
    addLink!(net = net, fromBus = "R2", toBus = "R3", status = 1)
    addLink!(net = net, fromBus = "R3", toBus = "R1", status = 1)
    ok, msg = validate!(net = net)
    @test ok
    ite, erg = runpf!(net, 40, 1e-10, 0; method = :rectangular)
    @test erg == 0
    @test ite > 0
    calcNetLosses!(net)

    # measurement-free: bitwise regression against the KCL allocation
    calcLinkFlowsKCL!(net)
    fK = [(l.pFlow_MW, l.qFlow_MVar) for l in net.linkVec]
    rep0 = calcLinkFlowsSE!(net)
    fS = [(l.pFlow_MW, l.qFlow_MVar) for l in net.linkVec]
    @test fS == fK
    @test all(r.source == :kcl for r in rep0)
    @test all(isnan(r.p_meas_residual) for r in rep0)

    # one P measurement on link 1 pins the ring circulation: the measured
    # link matches z within sigma and the nodal balances still hold
    z = fK[1][1] + 5.0
    σ = 0.01
    addPflowMeasurement!(net; value = z, sigma = σ, linkNr = 1)
    rep1 = calcLinkFlowsSE!(net)
    @test rep1[1].source == :measured_ls
    @test abs(net.linkVec[1].pFlow_MW - z) <= 3.0 * σ
    @test abs(rep1[1].p_meas_residual) <= 3.0 * σ
    @test rep1[2].source == :kcl && rep1[3].source == :kcl
    # component KCL invariance: the ring's null space is the pure circulation
    # (all three links oriented around the ring), so the measurement may only
    # ADD one common circulation to the minimum-norm split. Equal deltas on
    # all three links mean every nodal balance is untouched.
    fM = [l.pFlow_MW for l in net.linkVec]
    d = [fM[j] - fK[j][1] for j in 1:3]
    @test isapprox(d[1], d[2]; atol = 1e-8)
    @test isapprox(d[2], d[3]; atol = 1e-8)
    @test abs(d[1]) > 1.0   # the measurement really moved the circulation

    # a single link between two buses: the nodal balance fixes its flow, so
    # a gross measurement on it (100 MW off, sigma 2 MW) must not move it,
    # the error shows in the residual. With the two balance rows at weight 1
    # against 1/sigma^2 = 0.25 the measurement would pull the flow by
    # 100 * 0.25 / 2.25 = 11.1 MW.
    net1 = Net(name = "w2_single", baseMVA = 100.0)
    for b in ("S", "A", "B", "L")
      addBus!(net = net1, busName = b, vn_kV = 110.0)
    end
    addPIModelACLine!(net = net1, fromBus = "S", toBus = "A", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = net1, fromBus = "B", toBus = "L", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
    addProsumer!(net = net1, busName = "S", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "S")
    addProsumer!(net = net1, busName = "L", type = "ENERGYCONSUMER", p = 50.0, q = 10.0)
    addLink!(net = net1, fromBus = "A", toBus = "B", status = 1)
    ok1, _ = validate!(net = net1)
    @test ok1
    _, erg1 = runpf!(net1, 40, 1e-10, 0; method = :rectangular)
    @test erg1 == 0
    calcNetLosses!(net1)
    calcLinkFlowsKCL!(net1)
    fK1 = net1.linkVec[1].pFlow_MW
    @test abs(fK1 - 50.0) < 1.0   # the load flow through the link
    addPflowMeasurement!(net1; value = fK1 - 100.0, sigma = 2.0, linkNr = 1)
    rep = calcLinkFlowsSE!(net1)
    @test rep[1].source == :measured_ls
    @test abs(net1.linkVec[1].pFlow_MW - fK1) < 1e-6
    @test abs(rep[1].p_meas_residual - 100.0) < 1e-6
  end)() end

  return true
end

function test_state_estimation_islands()::Bool
  # Island-wise SE (0.10.0): a two-island net (no connection between the
  # groups) is partitioned and every measured island is estimated with its
  # own reference; the merged SEResult carries per-island rows. Also covers
  # the snapshot-vs-state chain distinction on a measurement/model
  # discrepancy, which the noise-free phase-5 chain test could not see
  # (there takeover == model and the two modes trivially coincide).
  @testset "State estimation island-wise" begin (function ()
    net = _se_two_island_net()
    ite, erg = runpf!(net, 40, 1e-10, 0; islands_enabled = true)
    @test erg == 0
    vm_true = [n._vm_pu for n in net.nodeVec]
    meas = generateMeasurementsFromPF(net; noise = false)

    res = with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = false) do
      runse!(net, meas)
    end
    @test res.converged
    @test res.islands !== nothing
    @test length(res.islands) == 2
    @test all(i.estimated for i in res.islands)
    @test all(i.converged for i in res.islands)
    # per-island accounting present (dof per island, own band verdict)
    @test all(i.dof > 0 for i in res.islands)
    @test all(hasproperty(i, :z_wh) for i in res.islands)
    @test res.dof == sum(i.dof for i in res.islands)
    # merged voltages are in ORIGINAL bus numbering and match the PF truth
    @test maximum(abs.(abs.(res.voltages) .- vm_true)) < 1e-6

    # a shunt measurement in the second island (0.30.2): the island subnet
    # looked shunts up in the dictionary of the whole net, keyed by the
    # net-wide bus numbers, and failed with "no shunt registered at this bus"
    netS = _se_two_island_net()
    addShunt!(net = netS, busName = "B2", pShunt = 0.0, qShunt = 5.0)
    runpf!(netS, 40, 1e-10, 0; islands_enabled = true)
    vm_shunt = [n._vm_pu for n in netS.nodeVec]
    measS = generateMeasurementsFromPF(netS; noise = false, includeShuntQ = true)
    @test any(m -> m.typ == Sparlectra.ShuntQMeas, measS)
    resS = with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = false) do
      runse!(netS, measS)
    end
    @test resS.converged
    @test maximum(abs.(abs.(resS.voltages) .- vm_shunt)) < 1e-6

    # the printed diagnostics carry the per-island band test (no island can
    # hide inside the summed chi-square)
    rep = with_state_estimation_config(max_iter = 30, tol = 1e-8) do
      validate_measurements(net, meas)
    end
    @test rep.islands !== nothing && count(i -> i.measured, rep.islands) == 2
    io = IOBuffer()
    print_se_diagnostics(io, rep; topN = 3)
    txt = String(take!(io))
    @test occursin("Islands", txt)
    @test occursin("island 1", txt) && occursin("island 2", txt)

    # unmeasured island: estimated where measured, skipped elsewhere
    bIdxs = Set([net.busDict["B1"], net.busDict["B2"]])
    measA = Measurement[m for m in meas if (m.busIdx !== nothing && !(m.busIdx in bIdxs)) || (m.branchIdx !== nothing && !(Int(net.branchVec[m.branchIdx].fromBus) in bIdxs))]
    resA = with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = false) do
      runse!(net, measA)
    end
    @test resA.converged
    @test count(i -> i.estimated, resA.islands) == 1
    skipped = only([i for i in resA.islands if !i.estimated])
    @test skipped.reason == :no_measurements

    # snapshot-vs-state distinction (regression for the inert phase-5
    # takeover): create a measurement/model discrepancy AFTER the
    # measurement set was taken, then chain both modes
    net2 = _se_two_island_net()
    runpf!(net2, 40, 1e-10, 0; islands_enabled = true)
    meas2 = generateMeasurementsFromPF(net2; noise = false)
    # model drifts after the measurements were taken: +8 MW load at A3
    addProsumer!(net = net2, busName = "A3", type = "ENERGYCONSUMER", p = 8.0, q = 0.0)
    with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = true) do
      runse!(net2, meas2)
    end
    nPros = length(net2.prosumpsVec)
    rState = runpf_from_se!(net2, 40, 1e-10, 0; mode = :se_state, islands_enabled = true)
    vmState = [n._vm_pu for n in net2.nodeVec]
    rSnap = runpf_from_se!(net2, 40, 1e-10, 0; mode = :se_snapshot, islands_enabled = true)
    vmSnap = [n._vm_pu for n in net2.nodeVec]
    @test rState.converged && rSnap.converged
    # state mode balances the drift through the slack; snapshot reproduces
    # the estimated operating point (zero mismatch at the start)
    @test abs(rState.slack_pickup_mw) > 4.0
    # the island-wise PF aggregates the per-island iteration counts, so the
    # zero-mismatch snapshot start allows one iteration per island
    @test rSnap.iterations <= 2
    @test rSnap.iterations < rState.iterations
    @test abs(rSnap.slack_pickup_mw) < 1e-6
    @test maximum(abs.(vmState .- vmSnap)) > 1e-5
    # the temporary delta prosumers are gone again
    @test length(net2.prosumpsVec) == nPros
  end)() end

  return true
end

function test_state_estimation_tap_pf_equivalence()::Bool
  # PF-equivalence gate of the step states: at the start position s0 of the
  # changer set the unstamped-Ybus + overlay predictions reproduce the
  # STAMPED power flow to 1e-10 for every changer case (kinds a to f, the
  # nameplate, typed and stored data paths), and setTapEstimation! returns
  # that s0. This is the correctness anchor of the whole tap overlay.
  @testset "State estimation tap PF equivalence gate" begin (function ()
    for c in _tap_changer_cases()
      net = c.net()
      c.build(net, c.n, true)
      ite, erg = runpf!(net, 40, 1e-12, 0; method = :rectangular)
      @test erg == 0
      V = buildVoltageVector(net)
      meas = generateMeasurementsFromPF(net; includeImag = true, noise = false)
      z = [m.value for m in meas]
      rel = setTapEstimation!(net; trafo = 3, mode = c.mode)
      @test rel.s0 !== nothing && isapprox(rel.s0[1], c.n; atol = 1e-9)
      nbus = length(net.nodeVec)
      tapMap, out, _ = Sparlectra._build_tap_state_map(net, (nbus - 1) + nbus + 1, :auto)
      @test tapMap !== nothing && isempty(out)
      Ybus = createYBUS(net = net)
      Sparlectra._tap_unstamp!(Ybus, net, tapMap)
      x = Sparlectra._initial_state_vector(net, Sparlectra._find_slack_idx(net); flatstart = false, tapSInit = Sparlectra._tap_s_init(tapMap))
      h = Sparlectra._predict_measurements(meas, net, V, Ybus; tapOverlay = Sparlectra._tap_overlay(x, tapMap, net))
      # 1e-10: the start position is read back from the live branch, one
      # floating-point roundtrip (about 5e-12 on the MVA scale)
      ok = maximum(abs.(h .- z)) <= 1e-10
      @test ok
      ok || println("      PF gate failed: ", c.name, " max ", maximum(abs.(h .- z)))
    end

    # regression (case57 corner finding at the parallel-trafo buses 4/18):
    # a LIVE tap off neutral (tap_ratio moved by a generator deviation or a
    # tap-controller step) must keep generated measurements self-consistent.
    # branchFlow_pu used to read the NEUTRAL ratio/angle while the Ybus
    # stamps the live tap, so the measured KCL missed at the trafo corners.
    knet = _tap_gate_net()
    kbr = knet.branchVec[3]
    kbr.tap_ratio = kbr.ratio / (1.0 + 2.0 * kbr.tap_step)   # 2 mechanical steps off neutral
    _, kerg = runpf!(knet, 40, 1e-10, 0)
    @test kerg == 0
    kmeas = generateMeasurementsFromPF(knet; noise = false)
    for bus in eachindex(knet.nodeVec)
      kpinj = only(m.value for m in kmeas if m.typ == Sparlectra.PinjMeas && m.busIdx == bus)
      kqinj = only(m.value for m in kmeas if m.typ == Sparlectra.QinjMeas && m.busIdx == bus)
      kpfl = 0.0
      kqfl = 0.0
      for m in kmeas
        m.typ in (Sparlectra.PflowMeas, Sparlectra.QflowMeas) || continue
        kbrm = knet.branchVec[m.branchIdx]
        ((m.direction == :from && Int(kbrm.fromBus) == bus) || (m.direction == :to && Int(kbrm.toBus) == bus)) || continue
        m.typ == Sparlectra.PflowMeas ? (kpfl += m.value) : (kqfl += m.value)
      end
      @test isapprox(kpinj, kpfl; atol = 1e-8)
      @test isapprox(kqinj, kqfl; atol = 1e-8)
    end
  end)() end

  return true
end

function test_state_estimation_tap_roundtrip()::Bool
  # Tap-estimation roundtrips: the TRUE network sits n steps off the
  # model's neutral position; measurements come from the true power flow
  # (noise-free). The estimator must recover the electrical step, the
  # mandatory fixation run must land on the exact step, and after the
  # fixation the regulator is no state any more: J after must be
  # numerically zero and the dof gains the freed state(s) back.
  @testset "State estimation tap roundtrip and fixation" begin (function ()
    step = 0.00625

    function truthmeas(n)
      tnet = _tap_gate_net()
      br = tnet.branchVec[3]
      # the fraction grid of the constructor: ratio / (1 + n step)
      br.tap_ratio = br.ratio / (1.0 + n * step)
      return _tap_truth_meas(tnet)
    end
    run(net, meas) = with_state_estimation_config(() -> runse!(net, meas); max_iter = 40, tol = 1e-10)

    # 1) :ratio, truth 2 steps up on the fraction grid
    meas = truthmeas(2)
    net = _tap_gate_net()
    setTapEstimation!(net; trafo = 3, mode = :ratio)
    res = run(net, meas)
    @test res.converged
    @test res.iterations < 40
    @test res.tapEstimates !== nothing && length(res.tapEstimates) == 1
    te = res.tapEstimates[1]
    @test te.branch == 3 && te.mode == :ratio && te.kind === :ratio && te.side === :to && te.fixed
    @test te.mrid == ""                       # not a CGMES case: no mRID
    @test abs(te.s_est - 2.0) < 0.05 && isapprox(te.du_est, te.s_est * step; rtol = 1e-14)
    @test te.step_fixed == 2
    # the model started at neutral, so the row says where the tap came from
    @test te.step_model == 0 && te.deviation_steps == 2
    @test !te.out_of_range
    tf = res.tapFixation
    @test tf !== nothing && tf.fixed
    @test tf.dof_after == tf.dof_before + 1 && tf.k_removed == 1
    @test tf.j_after < 1e-12
    @test !tf.j_rise_exceeds_quantile && !tf.sequential && tf.j_simultaneous == tf.j_after
    @test res.objectiveJ == tf.j_after && res.dof == tf.dof_after
    # Delta-J notice with k and the 95 % quantile of chi2(1) (3.841)
    @test occursin("k = 1 step state(s) removed, 95 % quantile of chi2(1) = 3.841", tf.delta_j_notice)

    # 2) :pst, truth 2 steps of a declared degree-form changer (kind f, no
    # direction to give); until 0.60.0 this ran on the undeclared degree
    # grid of the constructor with alpha 90 set by hand
    degree(n) = (tn = _tap_gate_net(); applyTapNameplate!(tn.branchVec[3]; phase_step_deg = 1.25, phase_min_step = -24.0, phase_max_step = 24.0, phase_current_step = Float64(n), psi_deg = 90.0); tn)
    meas = _tap_truth_meas(degree(2))
    net = degree(0)
    rel = setTapEstimation!(net; trafo = 3, mode = :pst)
    @test rel.alpha_deg === nothing && rel.s0 == [0.0]
    res = run(net, meas)
    @test res.converged && res.iterations < 40
    te = res.tapEstimates[1]
    @test te.mode == :pst && te.kind === :linear && te.du_est === nothing
    @test abs(te.s_est - 2.0) < 0.05
    @test te.step_fixed == 2
    tf = res.tapFixation
    @test tf.dof_after == tf.dof_before + 1
    @test tf.j_after < 1e-12

    # 3) :both, declared pair: truth 2 fraction steps and -1 Delta-u step
    # along psi 30 (both changers divide the branch ratio)
    meas = _tap_truth_meas(_tap_gate_pair_net(2, -1; psi = 30.0))
    net = _tap_gate_pair_net(0, 0; psi = 30.0)
    setTapEstimation!(net; trafo = 3, mode = :both, alpha_deg = 30.0)
    res = run(net, meas)
    @test res.converged && res.iterations < 40
    @test [t.kind for t in res.tapEstimates] == [:ratio, :asymmetrical] && all(t -> t.mode == :both, res.tapEstimates)
    @test abs(res.tapEstimates[1].s_est - 2.0) < 1e-6 && abs(res.tapEstimates[2].s_est - (-1.0)) < 1e-6
    @test [t.step_fixed for t in res.tapEstimates] == [2, -1]
    tf = res.tapFixation
    @test tf.dof_after == tf.dof_before + 2 && tf.k_removed == 2   # both regulators leave the state
    @test tf.j_after < 1e-10
    # enabled = false and mode = :ratio keep the stored direction, so a
    # release after them estimates the same as the first release (until
    # 0.60.0 both set alpha to 0, and this release was refused as collinear)
    te = res.tapEstimates
    net = _tap_gate_pair_net(0, 0; psi = 30.0)
    setTapEstimation!(net; trafo = 3, mode = :both, alpha_deg = 30.0)
    setTapEstimation!(net; trafo = 3, enabled = false)
    setTapEstimation!(net; trafo = 3, mode = :ratio)
    @test net.branchVec[3].tap_est_alpha_deg == 30.0
    setTapEstimation!(net; trafo = 3, mode = :both)
    res2 = run(net, meas)
    te2 = res2.tapEstimates
    @test res2.converged && [t.step_fixed for t in te2] == [2, -1]
    @test isapprox(te2[1].s_est, te[1].s_est; atol = 1e-9) && isapprox(te2[2].s_est, te[2].s_est; atol = 1e-9)
    # alpha comes from the import: a different one, one on a changer
    # without direction, and a phase release without phase data are refused
    @test_throws ArgumentError setTapEstimation!(net; trafo = 3, mode = :both, alpha_deg = 60.0)
    @test_throws ArgumentError setTapEstimation!(degree(0); trafo = 3, mode = :pst, alpha_deg = 90.0)
    @test_throws ArgumentError setTapEstimation!(_tap_gate_net(); trafo = 3, mode = :pst)

    # 4) off-grid truth (2.4 fraction steps): the electrical estimate sits
    # between steps, the fixation rounds to 2, and the residual J after is
    # honestly nonzero (the model cannot represent a half-step position)
    # The J rise above the 95 % quantile of chi2(1) gives a warning and the
    # flag; with one step state there is no sequential fixation
    meas = truthmeas(2.4)
    net = _tap_gate_net()
    setTapEstimation!(net; trafo = 3, mode = :ratio)
    res = @test_logs (:warn, r"above the 95 % quantile 3\.841; one step state") match_mode = :any run(net, meas)
    @test res.converged
    @test res.tapFixation.j_rise_exceeds_quantile && !res.tapFixation.sequential
    te = res.tapEstimates[1]
    @test abs(te.s_est - 2.4) < 0.05
    @test te.step_fixed == 2
    tf = res.tapFixation
    @test tf.j_before < 1e-12                 # continuous fit absorbs it fully
    @test tf.j_after > 1e-6                   # fixed step cannot
    @test res.objectiveJ == tf.j_after
    # noise-free sigmas overstate the errors, so the off-grid remainder
    # stays INSIDE the band here: no verdict flip, no off-grid note
    @test tf.offgrid_residual == false

    # 5) :offgrid_tap_residual: noisy set, truth half a
    # step off the grid, tightened sigmas. The continuous fit passes the
    # band, the fixation must round half a step away and fails it :high;
    # the note marks that flip as off-grid truth, NOT bad data.
    tightσ = measurementStdDevs(vm = 0.0025, pinj = 0.5, qinj = 0.5, pflow = 0.5, qflow = 0.5, imag = 5.0)
    tnet = _tap_gate_net()
    br = tnet.branchVec[3]
    br.tap_ratio = br.ratio / (1.0 + 2.5 * step)
    ite, erg = runpf!(tnet, 40, 1e-12, 0; method = :rectangular)
    @test erg == 0
    meas = generateMeasurementsFromPF(tnet; includeImag = true, noise = true, stddev = tightσ, rng = Random.MersenneTwister(7))
    net = _tap_gate_net()
    setTapEstimation!(net; trafo = 3, mode = :ratio)
    res = @test_logs (:warn, r"above the 95 % quantile 3\.841; one step state") match_mode = :any run(net, meas)
    @test res.converged
    tf = res.tapFixation
    @test tf.fixed && tf.j_rise_exceeds_quantile
    @test Sparlectra._band_test_verdict(tf.j_before, tf.dof_before).reason in (:ok, :low)
    @test Sparlectra._band_test_verdict(tf.j_after, tf.dof_after).reason == :high
    @test tf.offgrid_residual == true

    # 6) the fixation uses the steps of the changer, so a truth on them
    # comes back as its own step with J after zero: a CGMES changer
    # (tap_ratio_grid :linear, ratio + n * tap_step) at both range ends, a
    # CGMES table with uneven points, nameplate data with a neutral ratio
    # off 1 at its range end, and a typed ratio model. On the fraction grid
    # these came back as other steps (13 and 19 for the CGMES ends, the
    # nameplate end clamped to a band computed as if the ratio were 1).
    gridnet = function (kind::Symbol)
      if kind === :typed
        return _tap_gate_net(taps = PowerTransformerTaps(Vn_kV = 110.0, step = 0, lowStep = -16, highStep = 16, neutralStep = 0, voltageIncrement_kV = 1.1))
      end
      n = _tap_gate_net()
      br = n.branchVec[3]
      if kind === :nameplate
        br.ratio = 0.95
        br.tap_ratio = 0.95
        applyTapNameplate!(br; tap_step = 0.0125, tap_min_step = -8.0, tap_max_step = 8.0)
      else
        # the fields the CGMES import sets for a 1 % changer at its neutral
        br.tap_ratio_grid = :linear
        br.tap_step = 0.01
        br.tap_min, br.tap_max = 0.84, 1.16
        if kind === :table
          br.tap_ratio_points = [0.9, 0.93, 0.97, 1.0, 1.02, 1.05, 1.1]
          br.tap_min, br.tap_max = 0.9, 1.1
        end
      end
      return n
    end
    for (kind, n, truth) in ((:linear, 16, 1.16), (:linear, -16, 0.84), (:table, 2, 1.05), (:table, -3, 0.9), (:nameplate, -8, 0.95 / 0.9), (:typed, 16, 1.16))
      tnet = gridnet(kind)
      tnet.branchVec[3].tap_ratio = truth
      meas = _tap_truth_meas(tnet)
      net = gridnet(kind)
      setTapEstimation!(net; trafo = 3, mode = :ratio)
      res = run(net, meas)
      @test res.converged
      te = res.tapEstimates[1]
      @test abs(te.s_est - n) < 0.05
      @test te.step_fixed == n
      @test !te.out_of_range
      @test res.tapFixation.j_after < 1e-10
    end
  end)() end

  return true
end

function test_state_estimation_tap_guards()::Bool
  # Release guards: a bridge transformer without
  # a voltage pin on its cut-off side freezes structurally; machine (GSU)
  # transformers are detected for the mass-release skip; partial freezes
  # re-pack the state columns; a fully frozen release is bitwise identical
  # to no release at all.
  @testset "State estimation tap release guards" begin (function ()
    step = 0.00625

    function _radial_net(; p = 8.0, q = 2.0)
      net = _tap_gate_net()
      addBus!(net = net, busName = "L3", vn_kV = 10.0)
      addProsumer!(net = net, busName = "L3", type = "ENERGYCONSUMER", p = p, q = q)
      addPIModelTrafo!(net = net, fromBus = "L2", toBus = "L3", r_pu = 0.004, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
      net.branchVec[5].has_ratio_tap = true
      net.branchVec[5].ratio_tap_declared = true
      net.branchVec[5].tap_step = step
      return net
    end

    # 1) machine-transformer detection: generator bus hanging on one trafo
    gsu = _tap_gate_net()
    addBus!(net = gsu, busName = "G1", vn_kV = 10.0)
    addProsumer!(net = gsu, busName = "G1", type = "SYNCHRONOUSMACHINE", p = -20.0, q = -5.0)
    addPIModelTrafo!(net = gsu, fromBus = "H2", toBus = "G1", r_pu = 0.003, x_pu = 0.1, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
    @test Sparlectra._is_machine_transformer(gsu, 5)
    @test !Sparlectra._is_machine_transformer(gsu, 3)   # loop transformer
    @test !Sparlectra._is_machine_transformer(gsu, 4)

    # 2) radial trafo WITH a voltage pin on the far side: estimable, the
    # guard must not overblock (truth 2 steps, full measurement set)
    tnet = _radial_net()
    br = tnet.branchVec[5]
    br.tap_ratio = br.ratio / (1.0 + 2 * step)
    meas = _tap_truth_meas(tnet)
    net = _radial_net()
    setTapEstimation!(net; trafo = 5, mode = :ratio)
    res = with_state_estimation_config(max_iter = 40, tol = 1e-10) do
      runse!(net, meas)
    end
    @test res.converged
    te = res.tapEstimates[1]
    @test te.reason == :none
    @test te.step_fixed == 2
    @test res.tapFixation.j_after < 1e-12

    # 3) same release WITHOUT the far-side voltage measurement. At the flat
    # start no current flows, the start test keeps the tap out
    # (radial_no_voltage_pin); at the converged state the losses of the
    # 8 MW load tie the far-side voltage level to the measured flows, the
    # recheck at the estimate (0.60.0) finds the column determinable and
    # estimates it again. With exact rows the estimate is the truth step 2,
    # but sigma_s from the row sigmas is several steps, above the bound of
    # 0.5 for a re-estimate: the changer keeps its transferred step 0 with
    # :weak_at_estimate. Until the bound it was fixed at 2; until the recheck
    # ran on a run without a fixation the tap stayed frozen here.
    idxL3 = 5   # bus L3 is the fifth bus of _radial_net
    measNoVm = [m for m in meas if !(m.typ == Sparlectra.VmMeas && m.busIdx == idxL3)]
    net = _radial_net()
    setTapEstimation!(net; trafo = 5, mode = :ratio)
    resE = with_state_estimation_config(max_iter = 40, tol = 1e-10) do
      runse!(net, measNoVm)
    end
    rowE = only(resE.tapEstimates)
    @test resE.converged && rowE.status === :not_estimated && rowE.reason === :weak_at_estimate && rowE.step_fixed == 0
    @test (est = match(r"estimated again at step (\S+) with sigma_s (\S+) steps", rowE.notice); est !== nothing && abs(parse(Float64, est[1]) - 2.0) < 1e-3 && parse(Float64, est[2]) >= 0.5)
    # without load nothing flows at the estimate either, the tap is not
    # determinable there: noisy rows (seed 11) give the estimate a small
    # current, the recheck finds a rank, the joint estimate does not
    # converge and the state before it is put back.
    # The frozen run is BITWISE identical to a run without any release
    # (stamp restored, also after that failed joint estimate)
    tnet = _radial_net(p = 0.0, q = 0.0)
    br = tnet.branchVec[5]
    br.tap_ratio = br.ratio / (1.0 + 2 * step)
    measNoVm = [m for m in _tap_truth_meas(tnet; noise = true) if !(m.typ == Sparlectra.VmMeas && m.busIdx == idxL3)]
    net = _radial_net(p = 0.0, q = 0.0)
    setTapEstimation!(net; trafo = 5, mode = :ratio)
    resF = with_state_estimation_config(max_iter = 40, tol = 1e-10) do
      runse!(net, measNoVm)
    end
    @test resF.converged
    @test resF.tapEstimates !== nothing && length(resF.tapEstimates) == 1
    @test resF.tapEstimates[1].reason == :radial_no_voltage_pin
    @test occursin("bridge transformer", resF.tapEstimates[1].notice)
    # the voltage hint of a changer the rank test kept out: measured minus
    # estimated |V| at the terminal with a Vm
    # row, here only L2 (bus 4) is left
    idxL2 = 4
    hint = match(r"voltage hint: measured minus estimated \|V\| at L2 (\S+) pu", resF.tapEstimates[1].notice)
    mL2 = only(m for m in measNoVm if m.typ == Sparlectra.VmMeas && m.busIdx == idxL2)
    @test hint !== nothing && abs(parse(Float64, hint[1]) - (mL2.value - abs(resF.voltages[idxL2]))) <= 5e-3 * abs(mL2.value - abs(resF.voltages[idxL2])) + 1e-15
    @test resF.tapEstimates[1].fixed == false
    @test resF.tapFixation === nothing   # no tap state entered the solve
    ref = with_state_estimation_config(max_iter = 40, tol = 1e-10) do
      runse!(_radial_net(p = 0.0, q = 0.0), measNoVm)
    end
    @test resF.objectiveJ == ref.objectiveJ
    @test resF.dof == ref.dof
    @test resF.voltages == ref.voltages

    # 4) freezes per regulator re-pack the state columns (unit level): the
    # phase regulator of a :both release loses its column, the following
    # :ratio release moves up; freezing the rest restores the stamps bitwise
    net = _tap_gate_pair_net(0, 0; psi = 30.0)
    net.branchVec[4].tap_step = step
    net.branchVec[4].ratio_tap_declared = true
    setTapEstimation!(net; trafo = 3, mode = :both)
    setTapEstimation!(net; trafo = 4, mode = :ratio)
    col0 = (length(net.nodeVec) - 1) + length(net.nodeVec) + 1
    map0, _, _ = Sparlectra._build_tap_state_map(net, col0, :auto)
    @test map0 !== nothing && map0.ncols == 3
    Ybus = createYBUS(net = net)
    Yref = copy(Matrix(Ybus))
    Sparlectra._tap_unstamp!(Ybus, net, map0)
    map0.trafos[1].frozen[2] = :not_observable
    m1, moved = Sparlectra._tap_repack(map0, net, Ybus)
    @test m1.ncols == 2 && isempty(moved)
    @test m1.trafos[1].cols == [col0, 0] && m1.trafos[2].cols == [col0 + 1]
    m1.trafos[1].frozen[1] = :not_observable
    m1.trafos[2].frozen[1] = :not_observable
    m2, moved2 = Sparlectra._tap_repack(m1, net, Ybus)
    @test m2 === nothing && length(moved2) == 2
    @test Matrix(Ybus) == Yref
  end)() end

  return true
end

function test_state_estimation_tap_determinability()::Bool
  # Determinability guard of the step states, on the radial branch 5 of
  # _tap_radial_net (bus R hangs on L2):
  # with angle measurements on both sides a ratio plus phase regulator pair
  # comes back exactly (A2b); without them the phase regulator of a tree
  # edge cannot be determined and is frozen with a reason and a notice
  # instead of diverging (A3). Truths with the importer's own functions on
  # an in-memory CGMES store (cgmes_mapping.jl) or through applyTapNameplate!.
  @testset "State estimation tap determinability guard" begin (function ()
    S = Sparlectra
    C = S.CGMESImporter
    rtc = _tap_cim("_g_rtc", :RatioTapChanger, [:neutralStep => "0", :lowStep => "-10", :highStep => "10", :stepVoltageIncrement => "1"])
    q90 = _tap_cim("_g_q90", :PhaseTapChangerAsymmetrical, [:neutralStep => "0", :lowStep => "-10", :highStep => "10", :voltageStepIncrement => "1", :windingConnectionAngle => "90"])
    tab27 = _tap_cim("_g_tab", :PhaseTapChangerTabular, [:neutralStep => "14", :lowStep => "1", :highStep => "27"]; refs = Dict(:PhaseTapChangerTable => "_g_tt"))
    # a synthetic 27-point phase table: constant ratio, angle bending a little
    tangle(k) = 1.1 * (k - 14) + 0.01 * (k - 14)^2
    tpts = [_tap_cim("_g_tt_$k", :PhaseTapChangerTablePoint, [:step => string(k), :ratio => "1.0", :angle => string(tangle(k))]; refs = Dict(:PhaseTapChangerTable => "_g_tt")) for k in 1:27]
    store = _tap_cim_store(vcat([rtc, q90, tab27], tpts))
    rch(st) = S.SETapChanger(name = "RTC", kind = :ratio, role = :ratio, side = :from, du = 0.01, neutral = 0, low = -10, high = 10, step = st)
    qch(st) = S.SETapChanger(name = "PTC", kind = :asymmetrical, role = :phase, side = :from, du = 0.01, alpha_deg = 90.0, neutral = 0, low = -10, high = 10, step = st)
    tch(st) = S.SETapChanger(name = "PTT", kind = :tabular, role = :phase, side = :from, neutral = 14, low = 1, high = 27, step = st, table = [S.SETapPoint(step = k, ratio = 1.0, angle_deg = tangle(k)) for k in 1:27])
    # the import multiplies the from-side factors of both changers onto t_n
    # (_applyEndTaps): the RTC formula and the phase changer of the store
    function pair_net(phase::Symbol, nr::Int, np::Int)
      net = _tap_radial_net()
      br = net.branchVec[5]
      obj, ch = phase === :q90 ? (q90, qch) : (tab27, tch)
      neutral = phase === :q90 ? 0 : 14
      br.se_tap_changers = S.SETapChangerSet(ratio_n = 1.0, changers = [rch(nr), ch(neutral + np)], origin = :cgmes)
      m, a = C._phaseTapRatioShift_atstep(store, obj, neutral + np)
      br.tap_ratio = C._ratioTapCorrectionAt(rtc, nothing, Float64(nr))[1] * m
      br.phase_shift_deg = a
      return net
    end
    rowof(res, kind) = only(t for t in res.tapEstimates if t.kind === kind)

    # A2b) T4-like (ratio 1 % plus quadrature 1 %) and T5-like (ratio plus a
    # 27-point phase table), radial, Va on both sides: both steps come back
    for (phase, kind, np) in ((:q90, :asymmetrical, -4), (:tab, :tabular, 6))
      meas = _tap_truth_meas(pair_net(phase, 3, np); includeVa = true)
      net = pair_net(phase, 0, 0)
      setTapEstimation!(net; trafo = 5, mode = :both)
      res = _tap_se_run(net, meas)
      @test res.converged
      @test all(t -> t.reason === :none && t.note === :none, res.tapEstimates)
      @test abs(rowof(res, :ratio).s_est - 3.0) < 1e-6 && rowof(res, :ratio).step_fixed == 3
      @test abs(rowof(res, kind).s_est - np) < 1e-6 && rowof(res, kind).deviation_steps == np
      @test res.tapFixation.j_after < 1e-9
    end
    # the T4-like pair with one row whose weight lies 1e3 above every other
    # row (Vm at the slack bus H1): both regulators are still estimated and
    # come back, neither is frozen. The earlier column rule estimated this
    # case as well (tap columns a factor 1.73 below the strongest Vm column,
    # it froze at 1e6), so the case pins the behaviour; the regression
    # against the old rule is i2 below
    meas =_tap_truth_meas(pair_net(:q90, 3, -4); includeVa = true)
    isheavy(m) = m.typ == S.VmMeas && m.busIdx == 1
    wmax = maximum(1 / m.sigma^2 for m in meas if !isheavy(m))
    meas = [isheavy(m) ? Measurement(typ = S.VmMeas, value = m.value, sigma = 1 / sqrt(1e3 * wmax), busIdx = 1, id = "Vm_H1_heavy") : m for m in meas]
    net = pair_net(:q90, 0, 0)
    setTapEstimation!(net; trafo = 5, mode = :both)
    res = _tap_se_run(net, meas)
    @test res.converged && count(isheavy, meas) == 1
    @test all(t -> t.reason === :none, res.tapEstimates)
    @test abs(rowof(res, :ratio).s_est - 3.0) < 1e-6 && rowof(res, :ratio).step_fixed == 3
    @test abs(rowof(res, :asymmetrical).s_est + 4.0) < 1e-6 && rowof(res, :asymmetrical).step_fixed == -4

    # A3 a) T4-like without Va: the phase regulator is frozen (no angle
    # anchor), the ratio regulator estimated with note :conditional (the
    # frozen quadrature changer moves the magnitude too), no divergence
    meas = _tap_truth_meas(pair_net(:q90, 3, -4))
    net = pair_net(:q90, 0, 0)
    setTapEstimation!(net; trafo = 5, mode = :both)
    res = _tap_se_run(net, meas)
    @test res.converged
    @test rowof(res, :asymmetrical).reason === :no_angle_anchor && occursin("no angle measurement", rowof(res, :asymmetrical).notice)
    @test rowof(res, :ratio).reason === :none && rowof(res, :ratio).note === :conditional && !isempty(rowof(res, :ratio).notice)
    @test rowof(res, :ratio).step_fixed == 3
    # h) T5-like without Va: frozen as well, but no :conditional, the table
    # ratio is constant
    meas = _tap_truth_meas(pair_net(:tab, 3, 0))
    net = pair_net(:tab, 0, 0)
    setTapEstimation!(net; trafo = 5, mode = :both)
    res = _tap_se_run(net, meas)
    @test res.converged && rowof(res, :tabular).reason === :no_angle_anchor
    @test rowof(res, :ratio).note === :none && rowof(res, :ratio).step_fixed == 3

    # An elimination whose control solve converges only in the fallback
    # without tap estimation does not count (c758bc91, the fallback of
    # 0.60.0 since). On the shipped sp_case5 the bridge transformer br5 gets
    # a Delta-u changer of 5 % at psi 60 (truth step -6, start 0), no Va,
    # and the Vm row of its far bus a gross error of -0.05 pu; the run with
    # the released step does not converge (the fallback itself is tested in
    # test_state_estimation_tap_fallback). A second Vm meter at that bus
    # reads 0.07 pu above the truth: the base run converges with the tap
    # released, the top suspect is that meter, and the control solve without
    # it converges only in the fallback. The row goes back to active and the
    # loop stops with not_converged.
    case5 = joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.scf.json")
    m5 = redirect_stdout(() -> S.import_case(case5, S.resolve_config(S.DEFAULT_SPARLECTRA_CONFIG_PATH, case5).config; run_kind = :state_estimation), devnull)
    empty!(m5.net.measurements)
    du5(n) = (tn = deepcopy(m5.net); applyTapNameplate!(tn.branchVec[5]; phase_du_step = 0.05, phase_min_step = -10.0, phase_max_step = 10.0, phase_current_step = Float64(n), psi_deg = 60.0); tn)
    tn = du5(-6)
    _, erg = redirect_stdout(() -> runpf!(tn; config = S._copy_powerflow_with(m5.config.powerflow; tol = 1e-10, max_iter = 40)), devnull)
    @test erg == 0
    lv = Int(tn.branchVec[5].toBus)
    meas = [m.typ == S.VmMeas && m.busIdx == lv ? Measurement(typ = S.VmMeas, value = m.value - 0.05, sigma = m.sigma, busIdx = lv, id = m.id * "_gross") : m for m in generateMeasurementsFromPF(tn; includeImag = true, noise = true, rng = Random.MersenneTwister(1))]
    lvm = only(m for m in meas if m.typ == S.VmMeas && m.busIdx == lv)
    measc = [meas; Measurement(typ = S.VmMeas, value = lvm.value + 0.05 + 0.07, sigma = lvm.sigma, busIdx = lv, id = "Vm_counter")]
    diag5(n, ms, f) = S.with_sparlectra_config(m5.config) do
      with_state_estimation_config(() -> f(n, ms); max_iter = 40, tol = 1e-6, update_net = false)
    end
    net = du5(0)
    setTapEstimation!(net; trafo = 5, mode = :pst)
    dc5 = @test_logs (:warn, r"falling back") match_mode = :any diag5(net, measc, runse_diagnostics)
    @test dc5.diagnostics.converged && isempty(S._se_report_tap_fallback(dc5.diagnostics))
    @test dc5.reverted !== nothing && dc5.reverted.id == "Vm_counter" && isempty(dc5.eliminations) && dc5.stop_reason === :not_converged
    measc[end] = S._set_measurement_active(measc[end], false)
    ctl = @test_logs (:warn, r"falling back") match_mode = :any diag5(net, measc, validate_measurements)
    @test ctl.converged && S._se_report_tap_fallback(ctl) == [1]

    # b) a symmetrical changer alone on the radial edge without Va: frozen at
    # its symmetric start, and the run is bitwise the run without a release
    sym_net() = (n = _tap_radial_net(); n.branchVec[5].se_tap_changers = S.SETapChangerSet(ratio_n = 1.0, changers = [S.SETapChanger(name = "PS", kind = :symmetrical, role = :phase, side = :from, du = 0.02, neutral = 0, low = -8, high = 8, step = 0)], origin = :cgmes); n)
    meas = _tap_truth_meas(sym_net())
    net = sym_net()
    setTapEstimation!(net; trafo = 5, mode = :pst)
    res = _tap_se_run(net, meas)
    ref = _tap_se_run(sym_net(), meas)
    @test only(res.tapEstimates).reason === :symmetric_start && res.tapFixation === nothing
    @test res.voltages == ref.voltages && res.objectiveJ == ref.objectiveJ
    # the same without the Vm row of the far bus R (bus 5): the bridge rule
    # names only changers that move the magnitude, so the pure angle
    # changer keeps symmetric_start and the notice of a pure angle
    # regulator (the bridge cause and its text "the tap would absorb the
    # downstream voltage" were not true for it)
    measNoVm = [m for m in meas if !(m.typ == S.VmMeas && m.busIdx == 5)]
    net = sym_net()
    setTapEstimation!(net; trafo = 5, mode = :pst)
    t = only(_tap_se_run(net, measNoVm).tapEstimates)
    @test t.status === :not_estimated && t.reason === :symmetric_start && occursin("pure angle regulator", t.notice) && !occursin("downstream voltage", t.notice)

    # c) a Delta-u changer at psi 90 on the radial edge: at its start 0 the
    # magnitude does not move (symmetric start, frozen); from step 2 the
    # start test releases it, the estimate is the truth 3, and the recheck
    # at the estimate (0.60.0) withdraws it: |t| at the mirror
    # step -3 is the same and no angle row tells the two apart
    # (:mirror_step_fits), the transferred step 2 stays. Before the recheck it
    # was fixed at 3 with note :sign_from_start
    du90(n) = (tn = _tap_radial_net(); applyTapNameplate!(tn.branchVec[5]; phase_du_step = 0.01, phase_min_step = -10.0, phase_max_step = 10.0, phase_current_step = Float64(n), psi_deg = 90.0); tn)
    meas = _tap_truth_meas(du90(3))
    net = du90(0)
    setTapEstimation!(net; trafo = 5, mode = :pst)
    @test only(_tap_se_run(net, meas).tapEstimates).reason === :symmetric_start
    # the same start without the rows of branch 5 and the injections at its
    # buses: the noise test freezes it (zero column), and the notice carries
    # no note of the structural pass (":sign_from_start" was left on it
    # before, its text contradicted the freeze)
    blind = [m for m in meas if !(m.branchIdx == 5) && !(m.typ in (S.PinjMeas, S.QinjMeas) && m.busIdx in (4, 5))]
    net = du90(0)
    setTapEstimation!(net; trafo = 5, mode = :pst)
    t = only(_tap_se_run(net, blind).tapEstimates)
    @test t.reason === :not_observable && t.note === :none && occursin("do not see", t.notice) && !occursin("sign of its step", t.notice)
    net = du90(2)
    setTapEstimation!(net; trafo = 5, mode = :pst)
    res = _tap_se_run(net, meas)
    row = only(res.tapEstimates)
    @test res.converged && row.status === :not_estimated && row.reason === :mirror_step_fits && row.step_fixed == 2 && abs(row.s_est - 3.0) < 1e-6

    # d) the same changer at psi 60 starts with a moving magnitude: released
    # and exact without any angle measurement
    du60(n) = (tn = _tap_radial_net(); applyTapNameplate!(tn.branchVec[5]; phase_du_step = 0.01, phase_min_step = -10.0, phase_max_step = 10.0, phase_current_step = Float64(n), psi_deg = 60.0); tn)
    meas = _tap_truth_meas(du60(4))
    net = du60(0)
    setTapEstimation!(net; trafo = 5, mode = :pst)
    res = _tap_se_run(net, meas)
    @test res.converged && only(res.tapEstimates).reason === :none && abs(only(res.tapEstimates).s_est - 4.0) < 1e-6

    # e) two parallel Delta-u changers into the radial bus, no Va: at psi 90
    # the lower branch index (the tree edge) is frozen at its symmetric
    # start, the other is estimated only relative to that frozen one (note
    # :relative_to_frozen): not fixed, its continuous value stays in the
    # fit (before, it was fixed with :relative_only); at
    # psi 60 both are released, the higher index with :relative_only
    function par(psi::Float64)
      n = _tap_radial_net()
      addPIModelTrafo!(net = n, fromBus = "L2", toBus = "R", r_pu = 0.004, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
      for k in (5, 6)
        applyTapNameplate!(n.branchVec[k]; phase_du_step = 0.01, phase_min_step = -10.0, phase_max_step = 10.0, phase_current_step = 0.0, psi_deg = psi)
      end
      return n
    end
    for (psi, frozen5) in ((90.0, :symmetric_start), (60.0, :none))
      meas = _tap_truth_meas(par(psi))
      net = par(psi)
      setTapEstimation!(net; trafo = 5, mode = :pst)
      setTapEstimation!(net; trafo = 6, mode = :pst)
      res = _tap_se_run(net, meas)
      r5 = only(t for t in res.tapEstimates if t.branch == 5)
      r6 = only(t for t in res.tapEstimates if t.branch == 6)
      @test res.converged && r5.reason === frozen5 && r6.reason === :none
      if psi == 90.0
        @test r6.note === :relative_to_frozen && !r6.fixed && abs(r6.s_fixed) < 1e-6 && occursin("relative to the frozen phase regulator", r6.notice)
        @test res.tapFixation.k_removed == 0 && res.dof == res.tapFixation.dof_before
      else
        @test r6.note === :relative_only && r6.fixed && r6.step_fixed == 0
      end
    end

    # f) an Ia row at the slack-side terminal of the PST and a Va row on its
    # far side share the PMU time base. The rank test counts a current row
    # only when its predicted current at the start passes 3 sigma of its
    # magnitude (the Ia gate of the WLS loop): from a flat start the PST at
    # s = 0 between equal voltages carries no current, the Ia row does not
    # count and at the start the column lies in the span of the angle states (not
    # estimated, symmetric_start). From a warm start (the voltages of a power
    # flow of the start model) the Ia row counts and anchors the edge:
    # estimated, exact. Before the rank test the structural rule counted the
    # Ia row in both cases.
    meas = _tap_truth_meas(du90(3); includeVa = true, includeIa = true)
    rbus = 5
    meas = [m for m in meas if !(m.typ == S.VaMeas && m.busIdx != rbus) && !(m.typ == S.IaMeas && !(m.branchIdx == 5 && m.direction == :from))]
    @test count(m -> m.typ == S.VaMeas, meas) == 1 && count(m -> m.typ == S.IaMeas, meas) == 1
    net = du90(0)
    setTapEstimation!(net; trafo = 5, mode = :pst)
    res = _tap_se_run(net, meas)
    # (0.60.0) the recheck at the estimate: kept out at the flat start, the
    # run has no tap state, at the converged state the Ia row counts and the
    # changer is estimated again, on the truth 3, but with sigma_s of
    # several steps at the default sigmas (1 % per step): not taken over
    # (:weak_at_estimate), the transferred step 0 stays. The same path with
    # tight rows is taken over (recheck test of the se group). Until the
    # bound it was fixed at 3; until the recheck ran on a run without a
    # fixation it stayed symmetric_start
    rowR = only(res.tapEstimates)
    @test res.converged && rowR.status === :not_estimated && rowR.reason === :weak_at_estimate && rowR.step_fixed == 0.0
    @test (est = match(r"estimated again at step (\S+) with sigma_s (\S+) steps", rowR.notice); est !== nothing && abs(parse(Float64, est[1]) - 3.0) < 1e-3 && parse(Float64, est[2]) >= 0.5)
    net = du90(0)
    runpf!(net, 40, 1e-12, 0; method = :rectangular)
    setTapEstimation!(net; trafo = 5, mode = :pst)
    res = _tap_se_run(net, meas; flatstart = false)
    @test res.converged && only(res.tapEstimates).reason === :none && abs(only(res.tapEstimates).s_est - 3.0) < 1e-6

    # g) a regulator no row sees (no flow on its branch, no injection at its
    # buses): its column is zero within its own FD noise, it is not
    # estimated (:not_observable, notice "do not see"); i) the same with one
    # injection row back at a huge sigma (1e9): the row sees the regulator,
    # it is estimated (the former weighted column rule froze it),
    # and the hint "weak column" reads the unweighted sensitivity, so the
    # huge sigma gives no hint any more. The sz3 pattern with the weights
    # of a study case is the regression test below (it replaces the case
    # with a Vm row at sigma 1e-7 next to that row at sigma 1e3, a weight
    # ratio of 1e20)
    gmeas = _tap_truth_meas((tn = _tap_gate_net(); tn.branchVec[3].tap_ratio = 1.0 / (1.0 + 2 * 0.00625); tn))
    blind = [m for m in gmeas if !(m.branchIdx == 3) && !(m.typ in (S.PinjMeas, S.QinjMeas) && m.busIdx in (1, 3))]
    qinj3(sigma) = Measurement(typ = S.QinjMeas, value = only(m.value for m in gmeas if m.typ == S.QinjMeas && m.busIdx == 3), sigma = sigma, busIdx = 3, id = "Qinj_L1_wide")
    for (meas, frozen) in ((blind, :not_observable), (vcat(blind, [qinj3(1e9)]), :none))
      net = _tap_gate_net()
      setTapEstimation!(net; trafo = 3, mode = :ratio)
      res = _tap_se_run(net, meas)
      t = only(res.tapEstimates)
      @test res.converged && t.reason === frozen && !occursin("weak column", t.notice)
      if frozen === :none
        @test abs(t.s_est - 2.0) < 1e-6 && t.step_fixed == 2
      else
        @test occursin("do not see", t.notice)
      end
    end
    # The constellation of a study case: a passive bus Z on H2 over a
    # link-like line (x 1e-5 pu), its P and Q injection rows zero-injection
    # pseudo measurements with weight 5e5 against weight 1 of the ordinary
    # P and Q rows (sigma 1 MW, Mvar), far from the ratio changer on branch
    # 3 (truth 2 steps, start 0). The Q row of Z puts about 1e7 Mvar per pu
    # into the Vm column of Z, as a link aggregate Q row of the study case
    # put 2.97e6 into the Vm column of its bus. The former column rule
    # froze the changer: weighted column over its step 3331 against jac_eps
    # times the largest weighted Vm column 7071 (study case 832 to 873
    # against 2109); this test fails with that rule (reason not_observable,
    # step 0) and passes here. The rank test does not read the weights: the
    # changer is estimated, exact, and gets no weak-column hint.
    function study_link_net()
      n = _tap_gate_net()
      addBus!(net = n, busName = "Z", vn_kV = 110.0)
      addPIModelACLine!(net = n, fromBus = "H2", toBus = "Z", r_pu = 0.0, x_pu = 1e-5, b_pu = 0.0, status = 1)
      return n
    end
    tn = study_link_net()
    tn.branchVec[3].tap_ratio = 1.0 / (1.0 + 2 * 0.00625)
    runpf!(tn, 40, 1e-12, 0; method = :rectangular)
    zbus = only(i for (i, nd) in enumerate(tn.nodeVec) if getCompName(nd.comp) == "Z")
    meas = generateMeasurementsFromPF(tn; stddev = measurementStdDevs(vm = 0.004, pinj = 1.0, qinj = 1.0, pflow = 1.0, qflow = 1.0, imag = 5.0), relativeSigma = false)
    meas = [m.typ in (S.PinjMeas, S.QinjMeas) && m.busIdx == zbus ? Measurement(typ = m.typ, value = m.value, sigma = 1 / sqrt(5e5), busIdx = zbus, id = m.id * "_zi") : m for m in meas]
    @test count(m -> m.busIdx == zbus && 1 / m.sigma^2 > 4.9e5, meas) == 2
    net = study_link_net()
    setTapEstimation!(net; trafo = 3, mode = :ratio)
    res = _tap_se_run(net, meas)
    t = only(res.tapEstimates)
    @test res.converged && t.reason === :none && t.step_fixed == 2 && abs(t.s_est - 2.0) < 1e-6 && isempty(t.notice)
    # N1 on the bound itself: one tap column of the gate net on a row shaped
    # like a zero injection (value 0, terms of size 1e3 seen through a Va
    # column), entries a tenth and a hundred times the FD noise bound
    # n = 10 eps T sigma_c / jac_eps with T = |z| + sum (Vm |H_Vm| + |H_Va|)
    # = 1e3; a bound from |z| alone would be a thousand times lower and leave
    # the tenth released (arithmetic, not run against that code). Neither
    # test reads a weight. The hint compares the unweighted column over its
    # step (100 n / sigma_c = 2.2e-4) with jac_eps times the largest Vm
    # column (bus 1): 1e-3 with an entry 1e3 gives the hint, 1e-6 with an
    # entry 1 does not
    function colmap()
      n = _tap_gate_net()
      setTapEstimation!(n; trafo = 3, mode = :ratio)
      nb = length(n.nodeVec)
      return nb, first(S._build_tap_state_map(n, 2 * nb, :auto))
    end
    for (mult, vmcol, frozen, hint) in ((0.1, 1e3, :not_observable, false), (100.0, 1e3, :none, true), (100.0, 1.0, :none, false))
      nb, map = colmap()
      ts = only(map.trafos)
      noise = 10 * eps(Float64) * 1e3 * ts.sigma[1] / 1e-6
      H = zeros(2, 2 * nb)
      H[1, nb] = vmcol
      H[2, 1] = 1e3
      H[2, ts.cols[1]] = mult * noise
      S._tap_column_freezes!(map, H, [1.0, 0.0], ones(2 * nb), nb, 1e-6)
      hints = S._tap_weak_hints!(map, H, nb, 1e-6)
      @test ts.frozen[1] === frozen && (length(hints) == 1) == hint && (ts.weak[1] > 0.0) == hint
    end
    # the noise test is per row: a row whose tap entry is half its noise
    # bound does not freeze a column that another row sees at 100 times its
    # bound, both rows with T = 1e3; a pooled test weighted 1e6 on the first
    # row and 1 on the second would freeze it (sum w H^2 = 2.6e5 n^2 against
    # sum w n^2 = 1e6 n^2)
    let (nb, map) = colmap()
      ts = only(map.trafos)
      noise = 10 * eps(Float64) * 1e3 * ts.sigma[1] / 1e-6
      H = zeros(3, 2 * nb)
      H[1, nb] = 1.0
      H[2, 1] = 1e3
      H[2, ts.cols[1]] = 0.5 * noise
      H[3, 1] = 1e3
      H[3, ts.cols[1]] = 100.0 * noise
      S._tap_column_freezes!(map, H, [1.0, 0.0, 0.0], ones(2 * nb), nb, 1e-6)
      @test ts.frozen[1] === :none
    end

    # j) two changers along one direction on one branch: two RTCs (both ends
    # of a CGMES delivery), and a ratio changer next to an asymmetrical one
    # at alpha 0; the second of each pair is frozen (:not_separable)
    for (chs, mode) in (([rch(0), S.SETapChanger(name = "RTC2", kind = :ratio, role = :ratio, side = :to, du = 0.01, neutral = 0, low = -10, high = 10, step = 0)], :ratio), ([rch(0), S.SETapChanger(name = "PTC0", kind = :asymmetrical, role = :phase, side = :from, du = 0.01, alpha_deg = 0.0, neutral = 0, low = -10, high = 10, step = 0)], :both))
      sep() = (n = _tap_gate_net(); n.branchVec[3].se_tap_changers = S.SETapChangerSet(ratio_n = 1.0, changers = chs, origin = :cgmes); n)
      meas = _tap_truth_meas(sep())
      net = sep()
      setTapEstimation!(net; trafo = 3, mode = mode)
      res = _tap_se_run(net, meas)
      @test res.converged && [t.reason for t in res.tapEstimates] == [:none, :not_separable]
    end

    # k) a branch off its changer curves (its live shift of 2 deg lies on no
    # changer of its set: a ratio row without a phase row, the angle base 0),
    # released through the service mass release: the run completes, the row
    # says :model_mismatch with a notice (until 0.60.0 the r0 split threw
    # and the run was lost without SE metadata)
    mktempdir() do d
      txt = read(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"), String)
      # the phase row of tap changer 46 (branch 30) goes, the live angle of
      # branch 30 becomes 2 deg (theta in radians)
      ind(n) = " "^n
      ctrl2 = join([ind(12) * "},", ind(12) * "{", ind(14) * "\"alpha_deg\": 0.0,", ind(14) * "\"index\": 2,", ind(14) * "\"pos\": 0,", ind(14) * "\"pos_max\": 24,", ind(14) * "\"pos_min\": -24,", ind(14) * "\"step_deg\": 1.25", ind(12) * "}", ind(10) * "],", ind(10) * "\"id\": 46,"], "\n")
      theta0 = join([ind(8) * "\"theta\": 0.0,", ind(8) * "\"to_node\": 8,"], "\n")
      @test occursin(ctrl2, txt) && occursin(theta0, txt)
      txt = replace(txt, ctrl2 => join([ind(12) * "}", ind(10) * "],", ind(10) * "\"id\": 46,"], "\n"))
      txt = replace(txt, theta0 => join([ind(8) * "\"theta\": $(deg2rad(2.0)),", ind(8) * "\"to_node\": 8,"], "\n"))
      case = joinpath(d, "sp_case14_offcurve.scf.json")
      write(case, txt)
      set = joinpath(d, "offcurve.measurements.csv")
      SparlectraApp._se_generate_measurement_set(case, set, SparlectraApp.MeasurementGeneratorOptions(; noise = false, gross_k = 0.0, tap_steps = 0.0, include_i = false, sigma_u_pct = 0.5, sigma_i_pct = 1.0, sigma_p_pct = 1.0, sigma_q_pct = 1.0, sigma_ia_deg = 0.0, seed = 3))
      res = redirect_stdout(devnull) do
        SparlectraApp._run_state_estimation_service(case, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, joinpath(d, "run"), "offcurve", set; max_iter = 50)
      end
      dd = SparlectraApp.to_dict(res)
      @test dd["status"] == "succeeded"
      rows = [r for r in dd["metadata"]["se_tap_estimates"] if r["reason"] == "model_mismatch"]
      @test length(rows) == 1 && occursin("do not lie on the curves", rows[1]["notice"])
      @test occursin("notice: tap regulator", read(joinpath(d, "run", "run.log"), String))
    end

    # l) a phase table whose ratio has its minimum at the start point: a
    # vertex, the one-sided slopes of |t| change sign there. Until the rank
    # test the structural rule froze it (:symmetric_start). The step column
    # follows the tangent of the table to the right of the start (ratio
    # slope 0.002 per step), which the far-side angles do not explain, so
    # the rank test keeps it: estimated from its start, it stays on the side
    # of the start (note :sign_from_start), the truth at the start comes
    # back, and the vertex notice says the step is not reliable
    vtab() = (n = _tap_radial_net(); n.branchVec[5].se_tap_changers = S.SETapChangerSet(ratio_n = 1.0, changers = [S.SETapChanger(name = "PTV", kind = :tabular, role = :phase, side = :from, neutral = 0, low = -4, high = 4, step = 0, table = [S.SETapPoint(step = k, ratio = 1.0 + 0.002 * k^2, angle_deg = 3.0 * k) for k in -4:4])], origin = :cgmes); n)
    meas = _tap_truth_meas(vtab())
    net = vtab()
    setTapEstimation!(net; trafo = 5, mode = :pst)
    res = _tap_se_run(net, meas)
    row = only(res.tapEstimates)
    @test res.converged && row.reason === :none && row.note === :sign_from_start && row.step_fixed == 0
    @test occursin("its step is not reliable", row.notice)

    # a Delta-u changer at psi 90 started at step 2 with the
    # truth at 0, the vertex of |t|: released (its start is no vertex) and
    # estimated; once with exact measurements, once with noise. Since the
    # recheck at the estimate (0.60.0) none of these runs keeps
    # the step: at the vertex 0 the step column has no rank
    # (:not_observable_at_estimate), at any other fixed step the mirror step
    # fits as well (:mirror_step_fits); the transferred step 2 stays. Until
    # the recheck they were fixed (exact at 0, seed 6 at 4, seed 5 at 0) with
    # the vertex notice.
    # Exact: near the vertex h moves with s^2, Gauss-Newton halves s per
    # iteration (measured 28 iterations to the stop test on s du, against 3
    # away from the vertex), the estimate goes to 0
    meas = _tap_truth_meas(du90(0))
    net = du90(2)
    setTapEstimation!(net; trafo = 5, mode = :pst)
    res = _tap_se_run(net, meas)
    row = only(res.tapEstimates)
    @test res.converged && !row.fixed && row.reason === :not_observable_at_estimate && row.step_fixed == 2 && abs(row.s_est) < 0.01
    @test res.tapFixation.iterations_estimation <= 40
    # noise, tight sigmas (vm 0.0005 pu, P and Q 0.1 MW, I 1 A), seed 6:
    # converged, the estimate rounds to 4 (measured s 3.95), far from the
    # truth 0; J at the fixed step 4 and at the mirror step -4 are equal
    # (27.21), so the recheck withdraws it
    tight = measurementStdDevs(vm = 0.0005, pinj = 0.1, qinj = 0.1, pflow = 0.1, qflow = 0.1, imag = 1.0)
    tn = du90(0)
    runpf!(tn, 40, 1e-12, 0; method = :rectangular)
    meas = generateMeasurementsFromPF(tn; includeImag = true, noise = true, stddev = tight, rng = Random.MersenneTwister(6))
    net = du90(2)
    setTapEstimation!(net; trafo = 5, mode = :pst)
    res = _tap_se_run(net, meas)
    row = only(res.tapEstimates)
    @test res.converged && !row.fixed && row.reason === :mirror_step_fits && row.step_fixed == 2 && round(row.s_est) == 4 && occursin("at the mirror step -4", row.notice)
    # noise, the generator's default sigmas, seed 5: the data put the best
    # |t| below the turning point, the step state moves half the way to it,
    # then onto it once the remaining magnitude effect is below the FD
    # tolerance (_tap_keep_side_pins!): converged in 6 iterations, s 0,
    # the vertex, withdrawn by the recheck. Before the side rule: no
    # convergence in 40 iterations, the state jumped between -3 and 2
    meas = _tap_truth_meas(du90(0); noise = true, rng = Random.MersenneTwister(5))
    net = du90(2)
    setTapEstimation!(net; trafo = 5, mode = :pst)
    res = _tap_se_run(net, meas)
    row = only(res.tapEstimates)
    @test res.converged && !row.fixed && row.reason === :not_observable_at_estimate && row.step_fixed == 2 && row.s_est == 0.0
    @test res.tapFixation.iterations_estimation <= 10
    # seeds 1 to 12: every run converges and stays on the side of the start
    # (s >= 0) with J at most that of the run without the release. Before
    # the side rule seeds 2, 5, 7 and 9 did not converge and seeds 3 and 8
    # ended at s -2.69 and -1.85. Three seeds end above the band end 10
    # (s 10.55, 12.23, 12.52); each such run logs exactly one
    # out-of-band warning naming the mirror step, the others none. Every run
    # is withdrawn by the recheck at the estimate (seeds 2, 5, 7 and 9 at the
    # vertex 0, the others by the mirror step)
    for seed in 1:12
      meas = _tap_truth_meas(du90(0); noise = true, rng = Random.MersenneTwister(seed))
      net = du90(2)
      setTapEstimation!(net; trafo = 5, mode = :pst)
      oorLog = Test.TestLogger(min_level = Logging.Warn)
      res = Logging.with_logger(() -> _tap_se_run(net, meas), oorLog)
      row = only(res.tapEstimates)
      ref = with_state_estimation_config(() -> runse!(du90(2), meas); max_iter = 40, tol = 1e-10)
      @test res.converged && !row.fixed && row.reason in (:mirror_step_fits, :not_observable_at_estimate) && row.step_fixed == 2 && row.s_est >= 0.0 && res.tapFixation.j_before <= ref.objectiveJ + 1e-6
      @test [r.message for r in oorLog.logs] == (row.out_of_range ? [string("runse!: ", S._tap_out_of_range_line(row))] : String[])
      @test !row.out_of_range || occursin("lies outside the band -10 to 10, |t| has the same value at the mirror step", row.notice)
    end

    # h) mesh edges only relative to a frozen phase regulator, on the
    # shipped sp_case60: br70 to br74 each carry a ratio and a degree phase
    # changer (1.25 deg per step, br72 1.0 deg), br70, br71 and br72 run in
    # parallel. Without Va the phase regulators of the tree edges br70, br73
    # and br74 freeze (:no_angle_anchor); br71 and br72 close loops over
    # br70, so only their angle against the frozen br70 is measured. Truth:
    # every changer 2 steps up, exact measurements. br70 stays at 0, so
    # br72 gives 2 - 2 * 1.25 / 1.0 = -0.5 and br71 gives 0. Before, br72
    # was fixed at -1 (fixed = true, note :relative_only) and J rose from
    # 1.9e-23 to 19.3 at dof 357; now both stay states, the
    # fit keeps J at zero and every ratio step is fixed at the truth
    case60 = joinpath(dirname(@__DIR__), "data", "scf", "sp_case60.scf.json")
    m60 = redirect_stdout(() -> S.import_case(case60, S.resolve_config(S.DEFAULT_SPARLECTRA_CONFIG_PATH, case60).config; run_kind = :state_estimation), devnull)
    empty!(m60.net.measurements)
    tn = deepcopy(m60.net)
    for k in 70:74
      br = tn.branchVec[k]
      set = S._effective_tap_changers(br)
      s0 = S._tap_s_from_branch(set, br)
      S._tap_apply_position!(br, set, round.(s0) .+ 2.0; s_from = s0, impedance = m60.config.state_estimation.tap_estimation.impedance_correction)
    end
    _, erg = redirect_stdout(() -> runpf!(tn; config = S._copy_powerflow_with(m60.config.powerflow; tol = 1e-10, max_iter = 40)), devnull)
    @test erg == 0
    meas = generateMeasurementsFromPF(tn; stddev = measurementStdDevs(vm = 0.002, pinj = 0.5, qinj = 0.5, pflow = 0.5, qflow = 0.5), relativeSigma = false)
    net = deepcopy(m60.net)
    for k in 70:74
      setTapEstimation!(net; trafo = k, mode = :both)
    end
    res = S.with_sparlectra_config(m60.config) do
      with_state_estimation_config(() -> runse!(net, meas); max_iter = 50, tol = 1e-6, update_net = false)
    end
    phase = Dict(t.branch => t for t in res.tapEstimates if t.role === :phase)
    @test res.converged && [phase[k].reason for k in 70:74] == [:no_angle_anchor, :none, :none, :no_angle_anchor, :no_angle_anchor]
    name70 = getCompName(net.branchVec[70].comp)
    for k in (71, 72)
      @test phase[k].note === :relative_to_frozen && !phase[k].fixed && occursin(name70, phase[k].notice)
    end
    @test abs(phase[72].s_est + 0.5) < 1e-6 && abs(phase[72].s_fixed + 0.5) < 1e-6 && abs(phase[71].s_fixed) < 1e-6
    @test all(t -> t.fixed && t.step_fixed == t.neutral_step + 2, (t for t in res.tapEstimates if t.role === :ratio))
    tf = res.tapFixation
    @test tf.k_removed == 5 && tf.dof_after == tf.dof_before + 5 && tf.j_after < 1e-9 && !tf.j_rise_exceeds_quantile && !tf.sequential
    # the diagnostics Jacobian carries the same extended state as the run,
    # br71 and br72 stay states there too: the trace of the residual
    # sensitivity matrix is m - n (wii = w_i Omega_ii), so sum(wii) is the
    # dof of the run. Without the two columns the sum would be dof + 2.
    d60 = S.with_sparlectra_config(m60.config) do
      with_state_estimation_config(() -> validate_measurements(deepcopy(net), meas); max_iter = 50, tol = 1e-6, update_net = false)
    end
    @test d60.result.dof == res.dof && length(d60.residuals) == length(res.residuals)
    @test abs(sum(d60.residual_sensitivities) - res.dof) < 1e-3
    @test maximum(abs.(d60.residuals .- res.residuals)) < 1e-6
    # rn and wii of one row against the run itself: a gross error e of 10
    # sigma on the active flow of br71 or br72 moves the run residual by
    # wii * e (wii = dr/dz), so wii = r/e and rn = r/(sigma sqrt(wii)) come
    # from the run alone. The diagnostics match to wii 7e-6 and rn 4e-5;
    # without the tap columns wii was 0.729 against 0.670 and rn 7.84
    # against 8.18 on br71.
    for id in ("Pflow_branch_71_from", "Pflow_branch_72_from")
      i = only(k for k in eachindex(meas) if meas[k].id == id)
      mg = meas[i]
      e = 10.0 * mg.sigma
      measg = copy(meas)
      measg[i] = Measurement(typ = mg.typ, value = mg.value + e, sigma = mg.sigma, branchIdx = mg.branchIdx, direction = mg.direction, id = mg.id)
      resg, dg = S.with_sparlectra_config(m60.config) do
        with_state_estimation_config(max_iter = 50, tol = 1e-6, update_net = false) do
          (runse!(deepcopy(net), measg), validate_measurements(deepcopy(net), measg))
        end
      end
      @test resg.converged && all(t -> t.step_fixed == t.neutral_step + 2, (t for t in resg.tapEstimates if t.role === :ratio))
      rk = only(x for x in dg.measurement_ranking if x.measurement_index == i)
      rr = resg.residuals[rk.local_index]
      @test abs(rk.residual - rr) < 1e-9
      @test abs(rk.wii - rr / e) < 1e-3
      @test abs(rk.normalized_residual - rr / (mg.sigma * sqrt(rr / e))) < 1e-2
    end

    # i) the frozen one is frozen outside the angle graph: branch 4, parallel
    # to branch 3, carries a stored phase table, its live shift lies on no
    # point (:model_mismatch), and the model shift 3.5 deg differs from the
    # truth 3.0 deg. Branch 3 is a typed Delta-u winding at psi 90 (truth
    # step 2, start 0), so only its angle against branch 4 is measured:
    # s_est 1.127, not fixed, the notice names branch 4. Before, branch 4
    # tied its buses into one class, branch 3 was an edge inside one
    # class and was rounded and fixed (code reading, not run). The
    # diagnostics stamp keeps the continuous position: J of the diagnostics
    # residuals equals the run J (0.0190); with the winding steps rounded
    # it was 58.2
    function mm(step3::Int, shift4::Float64)
      n = _tap_gate_net(phase_taps = PhaseTapChangerModel(kind = :asymmetrical, step = step3, lowStep = -10, highStep = 10, neutralStep = 0, voltage_step_increment = 0.01, winding_connection_angle_deg = 90.0))
      b4 = n.branchVec[4]
      b4.se_tap_changers = S.SETapChangerSet(ratio_n = 1.0, changers = [S.SETapChanger(name = "PT4", kind = :tabular, role = :phase, side = :from, neutral = 0, low = -4, high = 4, step = 0, table = [S.SETapPoint(step = k, ratio = 1.0 + 0.002 * k^2, angle_deg = 3.0 * k) for k in -4:4])], origin = :cgmes)
      b4.phase_shift_deg = shift4
      b4.angle = shift4
      return n
    end
    meas = _tap_truth_meas(mm(2, 3.0))
    net = mm(0, 3.5)
    setTapEstimation!(net; trafo = 3, mode = :pst)
    setTapEstimation!(net; trafo = 4, mode = :pst)
    res, d = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_net = false) do
      (runse!(deepcopy(net), meas), validate_measurements(deepcopy(net), meas))
    end
    r3 = only(t for t in res.tapEstimates if t.branch == 3)
    r4 = only(t for t in res.tapEstimates if t.branch == 4)
    @test res.converged && r4.reason === :model_mismatch && r3.reason === :none
    @test r3.note === :relative_to_frozen && !r3.fixed && occursin(getCompName(net.branchVec[4].comp), r3.notice)
    @test abs(r3.s_est - 1.1272) < 1e-3 && abs(r3.s_fixed - r3.s_est) < 1e-6 && res.tapFixation.k_removed == 0
    @test length(d.residuals) == length(meas) && isapprox(sum(abs2, d.residuals ./ [m.sigma for m in meas]), res.objectiveJ; rtol = 1e-6)
    # branch 3 stays a state in the diagnostics Jacobian as in the run:
    # sum(wii) = m - n = the dof of the run (dof + 1 without the column),
    # and the normalized residuals are the run's
    # residuals over sqrt(Omega_ii) = sigma_i sqrt(wii)
    @test d.result.dof == res.dof && abs(sum(d.residual_sensitivities) - res.dof) < 1e-3
    @test maximum(abs.(d.residuals .- res.residuals)) < 1e-6
    sig = [m.sigma for m in meas]
    loc = [i for i in eachindex(meas) if d.residual_sensitivities[i] > 1e-6]
    @test !isempty(loc) && maximum(abs.(d.normalized_residuals[loc] .- res.residuals[loc] ./ (sig[loc] .* sqrt.(d.residual_sensitivities[loc])))) < 1e-6
    # wii of one row against the run itself, as on sp_case60: a gross error
    # e of 10 sigma on the active flow of branch 3 moves the run residual by
    # wii * e, so wii = (r_e - r_0) / e comes from the run alone
    i3 = only(k for k in eachindex(meas) if meas[k].id == "Pflow_branch_3_from")
    m3 = meas[i3]
    e3 = 10.0 * m3.sigma
    measg = copy(meas)
    measg[i3] = Measurement(typ = m3.typ, value = m3.value + e3, sigma = m3.sigma, branchIdx = m3.branchIdx, direction = m3.direction, id = m3.id)
    resg, dg = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_net = false) do
      (runse!(deepcopy(net), measg), validate_measurements(deepcopy(net), measg))
    end
    rk = only(x for x in dg.measurement_ranking if x.measurement_index == i3)
    wrun = (resg.residuals[rk.local_index] - res.residuals[rk.local_index]) / e3
    @test resg.converged && abs(rk.wii - wrun) < 1e-3
  end)() end
  return true
end

# #381: the estimator reads its settings from the registry; a deviating run
# installs them for its duration and the registry is the same afterwards
function test_state_estimation_config_argument()::Bool
  @testset "SE settings from the registry, scoped overrides (#381)" begin (function ()
    scf_dir = joinpath(dirname(@__DIR__), "data", "scf")
    load5() = begin
      net = importSCF(joinpath(scf_dir, "sp_case5.scf.json"))
      readMeasurementTable!(net; file = joinpath(scf_dir, "sp_case5.measurements.csv"), zero_injection = false)
      net
    end
    before = state_estimation_config()
    # two runs in one session with different settings: one iteration cannot
    # converge, the full budget does; nothing is passed to runse! itself
    r1 = with_state_estimation_config(() -> runse!(load5()); max_iter = 1, tol = 1e-10)
    r2 = with_state_estimation_config(() -> runse!(load5()); max_iter = 40, tol = 1e-10)
    @test !r1.converged
    @test r2.converged
    @test r1.iterations < r2.iterations
    # the registry is the same object afterwards, also after an error inside
    @test state_estimation_config() === before
    @test_throws ErrorException with_state_estimation_config(() -> error("boom"); max_iter = 1)
    @test state_estimation_config() === before
    @test_throws ArgumentError with_state_estimation_config(() -> nothing; no_such_setting = 1)
    # a case sidecar with state_estimation.* reaches the run through the
    # resolved configuration installed for the run, the way the service
    # does it; no set_sparlectra_config! involved
    dir = mktempdir()
    for f in ("sp_case5.scf.json", "sp_case5.measurements.csv", "sp_case5.measurements.baddata.csv")
      cp(joinpath(scf_dir, f), joinpath(dir, f))
    end
    # (case-scope keys only, issue #377: the elimination budget and threshold)
    write(joinpath(dir, "sp_case5.config.yaml"), "config_version: 1\nscope: case\ncase: sp_case5.scf.json\nstate_estimation:\n  max_eliminations: 0\n  k_eliminate: 100.0\n")
    resolved = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, joinpath(dir, "sp_case5.scf.json")).config
    @test resolved.state_estimation.max_eliminations == 0
    @test resolved.state_estimation.k_eliminate == 100.0
    net = importSCF(joinpath(dir, "sp_case5.scf.json"))
    readMeasurementTable!(net; file = joinpath(dir, "sp_case5.measurements.baddata.csv"), zero_injection = false)
    dside = with_sparlectra_config(() -> runse_diagnostics(net), resolved)
    @test isempty(dside.eliminations)
    @test state_estimation_config() === before
    # without a wrapper the run behaves as before (registry defaults)
    r0 = runse!(load5())
    rg = with_sparlectra_config(() -> runse!(load5()), Sparlectra.ACTIVE_SPARLECTRA_CONFIG[])
    @test r0.converged == rg.converged && r0.iterations == rg.iterations
    @test isapprox(r0.objectiveJ, rg.objectiveJ; rtol = 1e-12)
    # the diagnostics and the precheck read the same registry
    d = with_state_estimation_config(() -> runse_diagnostics(load5()); k_eliminate = 100.0, max_eliminations = 1)
    @test d.stop_reason !== nothing
    pre = with_state_estimation_config(() -> validate_topology(load5()); topology_dead_flow_k = 3.0)
    @test pre.n_checked_branches > 0
  end)() end
  return true
end

function run_legacy_state_estimation_tests()
  # the same wrapper as the list in run_state_estimation_tests: each function
  # runs with the expected warnings of the group and must return true
  tests = [
    ("W2 link allocation", test_state_estimation_link_allocation_w2),
    ("Island-wise SE", test_state_estimation_islands),
    ("Tap PF equivalence gate", test_state_estimation_tap_pf_equivalence),
    ("Tap roundtrip and fixation", test_state_estimation_tap_roundtrip),
    ("Tap release guards", test_state_estimation_tap_guards),
    ("Tap determinability guard", test_state_estimation_tap_determinability),
    ("Run configuration as an argument (#381)", test_state_estimation_config_argument),
  ]
  for (name, testfn) in tests
    @testset "$name" begin (function ()
      @test _se_run_quiet(testfn) == true
    end)() end
  end
  return true
end
