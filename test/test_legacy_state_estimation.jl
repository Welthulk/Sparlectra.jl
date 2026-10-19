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
  # PF-equivalence gate for the tap cascade: at (r1_0, r2_0) the unstamped-Ybus + cascade-overlay predictions
  # reproduce the STAMPED power flow to 1e-12 for all three release
  # variants. This is the correctness anchor of the whole tap overlay.
  @testset "State estimation tap PF equivalence gate" begin (function ()
    for (mode, alpha, r1, r2) in ((:ratio, 0.0, 0.0125, 0.0), (:pst, 90.0, 0.0, 0.031), (:both, 30.0, 0.02, -0.015))
      net = _tap_gate_net()
      br = net.branchVec[3]
      ct = Sparlectra._cascade_tap(br, r1, r2, alpha)
      br.tap_ratio = ct.tap_ratio
      br.phase_shift_deg = ct.phase_shift_deg
      ite, erg = runpf!(net, 40, 1e-12, 0; method = :rectangular)
      @test erg == 0
      V = buildVoltageVector(net)
      meas = generateMeasurementsFromPF(net; includeImag = true, noise = false)
      z = [m.value for m in meas]
      rel = setTapEstimation!(net; trafo = 3, mode = mode, alpha_deg = mode == :ratio ? nothing : alpha)
      @test abs(rel.r1_0 - r1) < 1e-12 && abs(rel.r2_0 - r2) < 1e-12
      nbus = length(net.nodeVec)
      tapMap = Sparlectra._build_tap_state_map(net, (nbus - 1) + nbus + 1)
      @test tapMap !== nothing
      Ybus = createYBUS(net = net)
      Sparlectra._tap_unstamp!(Ybus, net, tapMap)
      x = Sparlectra._initial_state_vector(net, Sparlectra._find_slack_idx(net); flatstart = false, tapRInit = Sparlectra._tap_r_init(tapMap))
      h = Sparlectra._predict_measurements(meas, net, V, Ybus; tapOverlay = Sparlectra._tap_overlay(x, tapMap, net))
      # 1e-10: flows follow the LIVE tap since the branchFlow_pu fix, so the
      # r0 reconstruction (inverting tap_ratio back to a fraction) adds one
      # floating-point roundtrip (~5e-12 on MVA scale) that the old
      # neutral-ratio path never saw; still far below any physical tolerance
      @test maximum(abs.(h .- z)) <= 1e-10
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
  # Tap-estimation roundtrips: the TRUE network sits n mechanical steps off the
  # model's neutral position; measurements come from the true power flow
  # (noise-free). The estimator must recover the electrical step, the
  # mandatory fixation run must land on the exact mechanical step, and
  # after the fixation the tap is no state any more: J after must be
  # numerically zero and the dof gains the freed state(s) back.
  @testset "State estimation tap roundtrip and fixation" begin (function ()
    step = 0.00625
    pstep = 1.25

    function truthmeas(r1, r2, alpha)
      tnet = _tap_gate_net()
      br = tnet.branchVec[3]
      ct = Sparlectra._cascade_tap(br, r1, r2, alpha)
      br.tap_ratio = ct.tap_ratio
      br.phase_shift_deg = ct.phase_shift_deg
      ite, erg = runpf!(tnet, 40, 1e-12, 0; method = :rectangular)
      @test erg == 0
      return generateMeasurementsFromPF(tnet; includeImag = true, noise = false)
    end

    # 1) :ratio, truth 2 steps up on the fraction grid
    meas = truthmeas(2 * step, 0.0, 0.0)
    net = _tap_gate_net()
    setTapEstimation!(net; trafo = 3, mode = :ratio)
    res = with_state_estimation_config(max_iter = 40, tol = 1e-10) do
      runse!(net, meas)
    end
    @test res.converged
    @test res.iterations < 40
    @test res.tapEstimates !== nothing && length(res.tapEstimates) == 1
    te = res.tapEstimates[1]
    @test te.branch == 3 && te.mode == :ratio && te.fixed
    @test te.mrid == ""                       # not a CGMES case: no mRID
    @test abs(te.electrical_step_1 - 2.0) < 0.05
    @test te.fixed_step_1 == 2 && te.fixed_step_2 == 0
    # the model started at neutral, so the row says where the tap came from
    @test abs(te.model_step_1) < 1e-9 && abs(te.model_step_2) < 1e-9
    @test !te.out_of_range
    tf = res.tapFixation
    @test tf !== nothing && tf.fixed
    @test tf.dof_after == tf.dof_before + 1
    @test tf.j_after < 1e-12
    @test res.objectiveJ == tf.j_after && res.dof == tf.dof_after

    # 2) :pst, truth 2 steps on the shift grid (direction 90 degrees)
    meas = truthmeas(0.0, _tap_r2_for_shift(2 * pstep, 90.0), 90.0)
    net = _tap_gate_net()
    setTapEstimation!(net; trafo = 3, mode = :pst, alpha_deg = 90.0)
    res = with_state_estimation_config(max_iter = 40, tol = 1e-10) do
      runse!(net, meas)
    end
    @test res.converged && res.iterations < 40
    te = res.tapEstimates[1]
    @test te.mode == :pst && te.fixed_step_1 == 0
    @test abs(te.electrical_step_2 - 2.0) < 0.05
    @test te.fixed_step_2 == 2
    tf = res.tapFixation
    @test tf.dof_after == tf.dof_before + 1
    @test tf.j_after < 1e-12

    # 3) :both, truth 2 fraction steps and -1 shift step (direction 30)
    meas = truthmeas(2 * step, _tap_r2_for_shift(-pstep, 30.0), 30.0)
    net = _tap_gate_net()
    setTapEstimation!(net; trafo = 3, mode = :both, alpha_deg = 30.0)
    res = with_state_estimation_config(max_iter = 40, tol = 1e-10) do
      runse!(net, meas)
    end
    @test res.converged && res.iterations < 40
    te = res.tapEstimates[1]
    @test te.mode == :both
    @test abs(te.electrical_step_1 - 2.0) < 0.1
    @test abs(te.electrical_step_2 - (-1.0)) < 0.1
    @test te.fixed_step_1 == 2 && te.fixed_step_2 == -1
    tf = res.tapFixation
    @test tf.dof_after == tf.dof_before + 2   # both regulators leave the state
    @test tf.j_after < 1e-10

    # 4) off-grid truth (2.4 fraction steps): the electrical estimate sits
    # between steps, the fixation rounds to 2, and the residual J after is
    # honestly nonzero (the model cannot represent a half-step position)
    meas = truthmeas(2.4 * step, 0.0, 0.0)
    net = _tap_gate_net()
    setTapEstimation!(net; trafo = 3, mode = :ratio)
    res = with_state_estimation_config(max_iter = 40, tol = 1e-10) do
      runse!(net, meas)
    end
    @test res.converged
    te = res.tapEstimates[1]
    @test abs(te.electrical_step_1 - 2.4) < 0.05
    @test te.fixed_step_1 == 2
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
    ct = Sparlectra._cascade_tap(br, 2.5 * step, 0.0, 0.0)
    br.tap_ratio = ct.tap_ratio
    br.phase_shift_deg = ct.phase_shift_deg
    ite, erg = runpf!(tnet, 40, 1e-12, 0; method = :rectangular)
    @test erg == 0
    meas = generateMeasurementsFromPF(tnet; includeImag = true, noise = true, stddev = tightσ, rng = Random.MersenneTwister(7))
    net = _tap_gate_net()
    setTapEstimation!(net; trafo = 3, mode = :ratio)
    res = with_state_estimation_config(max_iter = 40, tol = 1e-10) do
      runse!(net, meas)
    end
    @test res.converged
    tf = res.tapFixation
    @test tf.fixed
    @test Sparlectra._band_test_verdict(tf.j_before, tf.dof_before).reason in (:ok, :low)
    @test Sparlectra._band_test_verdict(tf.j_after, tf.dof_after).reason == :high
    @test tf.offgrid_residual == true

    # 6) the fixation uses the ratio grid of the branch, so a truth on that
    # grid comes back as its own step with J after zero: a CGMES changer
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
      ite, erg = runpf!(tnet, 40, 1e-12, 0; method = :rectangular)
      @test erg == 0
      meas = generateMeasurementsFromPF(tnet; includeImag = true, noise = false)
      net = gridnet(kind)
      setTapEstimation!(net; trafo = 3, mode = :ratio)
      res = with_state_estimation_config(max_iter = 40, tol = 1e-10) do
        runse!(net, meas)
      end
      @test res.converged
      te = res.tapEstimates[1]
      @test abs(te.electrical_step_1 - n) < 0.05
      @test te.fixed_step_1 == n
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

    function _radial_net()
      net = _tap_gate_net()
      addBus!(net = net, busName = "L3", vn_kV = 10.0)
      addProsumer!(net = net, busName = "L3", type = "ENERGYCONSUMER", p = 8.0, q = 2.0)
      addPIModelTrafo!(net = net, fromBus = "L2", toBus = "L3", r_pu = 0.004, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
      net.branchVec[5].has_ratio_tap = true
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
    ct = Sparlectra._cascade_tap(br, 2 * step, 0.0, 0.0)
    br.tap_ratio = ct.tap_ratio
    br.phase_shift_deg = ct.phase_shift_deg
    ite, erg = runpf!(tnet, 40, 1e-12, 0; method = :rectangular)
    @test erg == 0
    meas = generateMeasurementsFromPF(tnet; includeImag = true, noise = false)
    net = _radial_net()
    setTapEstimation!(net; trafo = 5, mode = :ratio)
    res = with_state_estimation_config(max_iter = 40, tol = 1e-10) do
      runse!(net, meas)
    end
    @test res.converged
    te = res.tapEstimates[1]
    @test te.frozen_reason == :none
    @test te.fixed_step_1 == 2
    @test res.tapFixation.j_after < 1e-12

    # 3) same release WITHOUT the far-side voltage measurement: the
    # structural bridge guard freezes the tap, and the frozen run is
    # BITWISE identical to a run without any release (stamp restored)
    idxL3 = 5   # bus L3 is the fifth bus of _radial_net
    measNoVm = [m for m in meas if !(m.typ == Sparlectra.VmMeas && m.busIdx == idxL3)]
    net = _radial_net()
    setTapEstimation!(net; trafo = 5, mode = :ratio)
    resF = with_state_estimation_config(max_iter = 40, tol = 1e-10) do
      runse!(net, measNoVm)
    end
    @test resF.converged
    @test resF.tapEstimates !== nothing && length(resF.tapEstimates) == 1
    @test resF.tapEstimates[1].frozen_reason == :radial_no_voltage_pin
    @test resF.tapEstimates[1].fixed == false
    @test resF.tapFixation === nothing   # no tap state entered the solve
    ref = with_state_estimation_config(max_iter = 40, tol = 1e-10) do
      runse!(_radial_net(), measNoVm)
    end
    @test resF.objectiveJ == ref.objectiveJ
    @test resF.dof == ref.dof
    @test resF.voltages == ref.voltages

    # 4) partial freeze re-packs the state columns (unit level): a :both
    # release loses its r2 column, a following :ratio release moves up
    net = _tap_gate_net()
    setTapEstimation!(net; trafo = 3, mode = :both, alpha_deg = 30.0)
    setTapEstimation!(net; trafo = 4, mode = :ratio)
    col0 = (length(net.nodeVec) - 1) + length(net.nodeVec) + 1
    map0 = Sparlectra._build_tap_state_map(net, col0)
    @test map0 !== nothing && map0.ncols == 3
    Ybus = createYBUS(net = net)
    Yref = copy(Matrix(Ybus))
    Sparlectra._tap_unstamp!(Ybus, net, map0)
    m1 = Sparlectra._tap_apply_freezes(map0, [false, false], [true, false], net, Ybus)
    @test m1.ncols == 2
    @test m1.r1cols == [col0, col0 + 1] && m1.r2cols == [0, 0]
    # freezing the remaining columns restores the stamps bitwise
    m2 = Sparlectra._tap_apply_freezes(m1, [true, true], [true, true], net, Ybus)
    @test m2 === nothing
    @test Matrix(Ybus) == Yref
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
    ("Run configuration as an argument (#381)", test_state_estimation_config_argument),
  ]
  for (name, testfn) in tests
    @testset "$name" begin (function ()
      @test _se_run_quiet(testfn) == true
    end)() end
  end
  return true
end
