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

# file: test/test_observability.jl
# purpose: tests the observability side of state estimation: global/local
#          observability metrics and the labeled measurement Jacobian, the
#          matrix-only helpers (rank, matching, null-space dark states,
#          criticality), the FD-aware rank tolerance on a linked net, and
#          the Takahashi selected-inverse diagnostics (Omega_ii/wii/rn,
#          FD column coloring, diag(Omega) criticality, #394)

using Test
using Logging
using Sparlectra
using Random

# Runs alone with test/test_runner_helpers.jl (the link-net builder
# `create_se_link_net`, the warning allow-list `SE_EXPECTED_WARNINGS`, its
# runner `_se_run_quiet`, `large_case_path`) and test/testgrid.jl
# (`createTest3BusNet`) included first, as test/runtests.jl does.

"""
    observability_fixture() -> NamedTuple

Builds the networks the observability test sets share ONCE: the solved
3-bus net (`net3`, with the power-flow outcome in `pf3` so the tests can
still assert it), the solved shipped `sp_case14` (`net14`), the solved link
net with the bus link closed (`net_link`), the same net with the link open
(`net_link_open`, judged unsolved as before) and the solved link net with a
shunt on B2 (`net_link_shunt`). Only functions that take the measurement
vector explicitly are called on the shared nets, or the test sets
`net.measurements` themselves right before a convenience-form call, so the
sets stay order independent.
"""
function observability_fixture()
  net3 = createTest3BusNet()
  ite3, erg3 = runpf!(net3, 40, 1e-10, 0; method = :rectangular)
  net14 = importSCF(joinpath(pkgdir(Sparlectra), "data", "scf", "sp_case14.scf.json"))
  runpf!(net14, 40, 1e-8, 0; method = :rectangular)
  net_link = create_se_link_net()
  runpf!(net_link, 40, 1e-10, 0; method = :rectangular)
  net_link_open = create_se_link_net(link_status = 0)
  net_link_shunt = create_se_link_net(shunt_buses = ["B2"])
  runpf!(net_link_shunt, 40, 1e-10, 0; method = :rectangular)
  return (net3 = net3, pf3 = (ite = ite3, erg = erg3), net14 = net14, net_link = net_link, net_link_open = net_link_open, net_link_shunt = net_link_shunt)
end

function test_observability_metrics(fx)::Bool
  # Validates global/local observability metrics for full, reduced, and sparse
  # measurement sets, including expected not-observable behavior for sparse data.
  @testset "State estimation observability metrics" begin
    net = deepcopy(fx.net3)
    # the fixture solved the PF once; its outcome is asserted here, where the
    # measurements are generated from it
    @test fx.pf3.erg == 0
    @test fx.pf3.ite > 0

    std = measurementStdDevs(vm = 1e-5, pinj = 1e-4, qinj = 1e-4, pflow = 1e-4, qflow = 1e-4)
    meas = generateMeasurementsFromPF(net; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false, stddev = std, rng = MersenneTwister(42))

    gobs = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
      evaluate_global_observability(net, meas)
    end
    @test gobs.n_states == 2 * length(net.nodeVec) - 1
    @test gobs.n_measurements == count(m -> m.active, meas)

    # measurement_jacobian: the labeled H behind the observability checks
    # (same evaluation, plus described rows and state-column labels)
    empty!(net.measurements)
    append!(net.measurements, meas)
    mj = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
      measurement_jacobian(net)
    end
    @test size(mj.H) == (gobs.n_measurements, gobs.n_states)
    @test length(mj.rows) == gobs.n_measurements
    @test length(mj.cols) == gobs.n_states
    @test count(c -> startswith(c, "Va("), mj.cols) == length(net.nodeVec) - 1
    @test count(c -> startswith(c, "Vm("), mj.cols) == length(net.nodeVec)
    # a flow row is labeled with its oriented branch, a bus row with its bus
    @test any(r -> r.type == Sparlectra.PflowMeas && occursin("->", r.location), mj.rows)
    @test any(r -> r.type == Sparlectra.VmMeas && !occursin("->", r.location), mj.rows)
    # every row touches at least one state; the rank matches the check above
    @test all(i -> any(j -> abs(mj.H[i, j]) > 1e-12, axes(mj.H, 2)), axes(mj.H, 1))
    @test Sparlectra.numeric_rank(mj.H) == gobs.numerical_rank
    empty!(net.measurements)
    @test gobs.redundancy == gobs.n_measurements - gobs.n_states
    @test isapprox(gobs.redundancy_ratio, gobs.n_measurements / gobs.n_states)
    @test gobs.dof == gobs.redundancy
    @test gobs.structural_observable == true
    @test gobs.numerical_observable == true
    @test gobs.quality in (:good, :critical)

    nbus = length(net.nodeVec)
    local_cols = [1, nbus + 1]
    lobs = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
      evaluate_local_observability(net, meas, local_cols)
    end
    @test lobs.n_states == length(local_cols)
    @test lobs.n_measurements >= lobs.n_states
    @test !isempty(lobs.rows)
    @test lobs.redundancy == lobs.n_measurements - lobs.n_states
    @test lobs.structural_observable == true
    @test lobs.numerical_observable == true

    meas_reduced = copy(meas)
    for i in eachindex(meas_reduced)
      m = meas_reduced[i]
      if m.typ == Sparlectra.PflowMeas || m.typ == Sparlectra.QflowMeas
        meas_reduced[i] = Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = false, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id)
      end
    end

    gobs_reduced = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
      evaluate_global_observability(net, meas_reduced)
    end
    @test gobs_reduced.n_measurements < gobs.n_measurements
    @test gobs_reduced.structural_matching <= gobs.structural_matching
    @test gobs_reduced.numerical_rank <= gobs.numerical_rank

    meas_sparse = copy(meas)
    first_kept = false
    for i in eachindex(meas_sparse)
      m = meas_sparse[i]
      keep = !first_kept
      first_kept = true
      meas_sparse[i] = Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = keep, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id)
    end

    gobs_sparse = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
      evaluate_global_observability(net, meas_sparse)
    end
    @test gobs_sparse.quality == :not_observable
    @test gobs_sparse.numerical_observable == false

    @test_throws ErrorException with_state_estimation_config(() -> evaluate_local_observability(net, meas_sparse, local_cols); flatstart = true, jac_eps = 1e-6)
  end

  # The rank test runs on a column-normalized Jacobian. Without it the
  # relative tolerance is measured against the largest column scale, and on
  # a 25000-bus network a set with 203916 rows for 49999 states came out as
  # not observable (deficit 191 raw). Rows are deliberately not scaled by
  # sigma: rank is invariant under positive row scaling, and dividing by it
  # lifts zero-injection rows (sigma 1e-6) six orders of magnitude above the
  # rest, which made the same case worse (deficit 12308). Column
  # normalization alone leaves it at 0.
  # A reviewer's objection, and it was justified: the interesting failure is
  # not a set without any flow information (there the matching fails
  # obviously), it is one where the Hopcroft-Karp matching covers every
  # state column while information is still missing, so two angles are only
  # jointly determined. An earlier version let the structural result
  # overrule the numerical one and waved exactly those through as
  # observable. Thinning sp_case14 produces them: five of sixty draws had
  # full matching with a rank one or two short. The property below is what
  # must hold - a numerical deficit is never reported as observable, no
  # matter what the structure says.
  @testset "observability: structure never overrules a numerical deficit" begin
    net14 = deepcopy(fx.net14)
    full = generateMeasurementsFromPF(net14; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false)
    rng = MersenneTwister(11)
    seen_gap = false
    for _ = 1:40
      keep = Measurement[m for m in full if rand(rng) < 0.55]
      isempty(keep) && continue
      obs = try
        with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
          evaluate_global_observability(net14, keep)
        end
      catch
        continue
      end
      if obs.numerical_rank < obs.n_states
        @test obs.quality == :not_observable
        obs.structural_matching == obs.n_states && (seen_gap = true)
      end
    end
    # the draws must actually contain the interesting case, otherwise the
    # test above proves nothing
    @test seen_gap
  end

  @testset "observability: rank verdict is scale independent" begin
    netH = deepcopy(fx.net3)
    msH = generateMeasurementsFromPF(netH; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false)
    empty!(netH.measurements)
    append!(netH.measurements, msH)
    obsH = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
      evaluate_global_observability(netH, msH)
    end
    @test obsH.numerical_observable
    @test obsH.numerical_rank == obsH.n_states
    # scaling every measurement value by a constant must not change the
    # verdict: the rank question is about direction, not magnitude
    scaled = [Measurement(typ = m.typ, value = m.value, sigma = m.sigma * 1000, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id) for m in msH]
    obsS = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
      evaluate_global_observability(netH, scaled)
    end
    @test obsS.numerical_rank == obsH.numerical_rank
    # the convenience form forwards rankTolFactor (it silently ignored the
    # keyword before, so the configured value never reached the rank test)
    obsK = with_state_estimation_config(flatstart = true, jac_eps = 1e-6, rank_tol_factor = 0.1) do
      evaluate_global_observability(netH)
    end
    @test obsK.numerical_observable
  end

  return true
end

function test_observability_matrix_helpers()::Bool
  # Unit-checks matrix-only observability helpers (rank/matching/redundancy),
  # including local-column selection and expected error paths.
  @testset "Unobservable state columns (null-space dark states)" begin
    # the column-restricted local test is necessary but NOT sufficient: one
    # flow between two buses makes the 1x1 submatrix for column 1 full
    # rank, yet only the difference x1 - x2 is determined. The global
    # null-space answer names both columns as dark.
    H = [1.0 -1.0]
    sub = evaluate_local_observability_matrix(H, [1])
    @test sub.numerical_observable                         # the misleading positive verdict
    glob = evaluate_observability_matrix(H)
    @test !glob.numerical_observable
    @test glob.unobservable_state_columns == [1, 2]        # the rigorous answer

    # 7-bus workshop network: a spanning-tree flow set plus one Vm anchor
    # is observable (empty dark list); dropping the B6-B7 pair leaves
    # exactly B7's states dark: angle column 6 and magnitude column 13
    net7 = Net(name = "dark_states", baseMVA = 100.0)
    addBus!(net = net7, busName = "B1", vn_kV = 110.0, vm_pu = 1.02, va_deg = 0.0)
    for i = 2:7
      addBus!(net = net7, busName = "B$(i)", vn_kV = 110.0, vm_pu = 1.0, va_deg = 0.0)
    end
    ring = [("B1", "B2"), ("B2", "B3"), ("B3", "B4"), ("B4", "B5"), ("B5", "B6"), ("B6", "B7"), ("B7", "B1"), ("B2", "B5"), ("B3", "B6")]
    for (f, t) in ring
      addPIModelACLine!(net = net7, fromBus = f, toBus = t, r_pu = 0.010, x_pu = 0.080, b_pu = 0.0, status = 1)
    end
    addProsumer!(net = net7, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.02, va_deg = 0.0)
    addProsumer!(net = net7, busName = "B3", type = "GENERATOR", p = 60.0, q = 10.0)
    for (b, p, q) in [("B2", 35.0, 10.0), ("B4", 45.0, 15.0), ("B5", 25.0, 8.0), ("B6", 30.0, 10.0), ("B7", 20.0, 6.0)]
      addProsumer!(net = net7, busName = b, type = "LOAD", p = p, q = q)
    end
    ok7, msg7 = validate!(net = net7)
    ok7 || error("test net invalid: $msg7")
    _, erg7 = runpf!(net7, 40, 1e-10, 0)
    @test erg7 == 0
    calcNetLosses!(net7)

    tree = [("B1", "B2"), ("B2", "B3"), ("B3", "B4"), ("B4", "B5"), ("B5", "B6"), ("B6", "B7")]
    fill_tree! = function (n, pairs)
      empty!(n.measurements)
      for (f, t) in pairs
        addPflowMeasurement!(n; fromBus = f, toBus = t, value = get_branch_p_from_to_mw(n, f, t), sigma = 0.8, direction = :from)
        addQflowMeasurement!(n; fromBus = f, toBus = t, value = get_branch_q_from_to_mvar(n, f, t), sigma = 0.8, direction = :from)
      end
      addVmMeasurement!(n; busName = "B1", value = n.nodeVec[n.busDict["B1"]]._vm_pu, sigma = 0.002)
    end

    fill_tree!(net7, tree)
    intact = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
      evaluate_global_observability(net7)
    end
    @test intact.numerical_observable
    @test intact.unobservable_state_columns == Int[]

    fill_tree!(net7, tree[1:5])
    broken = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
      evaluate_global_observability(net7)
    end
    @test !broken.numerical_observable
    @test broken.unobservable_state_columns == [6, 13]
  end

  @testset "State estimation matrix observability helpers" begin
    H = [
      1.0 0.0 0.0
      0.0 1.0 0.0
      0.0 0.0 1.0
      1.0 1.0 0.0
      1.0 1.0 0.0
    ]

    obs = evaluate_observability_matrix(H)
    @test obs.numerical_observable == true
    @test obs.structural_observable == true
    @test obs.n_states == 3
    @test obs.n_measurements == 5
    @test obs.redundancy == 2

    @test numerical_observable(H) == true
    @test structural_observable(H) == true
    @test numerical_row_redundant(H, 3) == false
    @test structural_row_redundant(H, 3) == false

    lobs = evaluate_local_observability_matrix(H, [1, 2])
    @test lobs.n_states == 2
    @test lobs.n_measurements == 4
    @test !isempty(lobs.rows)

    @test_throws ErrorException evaluate_local_observability_matrix(H, Int[])
    @test_throws ErrorException evaluate_local_observability_matrix(H, [10])
  end

  return true
end

function test_observability_fd_rank_tolerance(fx)::Bool
  # Island-wise SE (0.10.0): two-stage observability.
  # An OPEN link is a real separation; with island-wise estimation each
  # measured island is judged on its own subnet and the merged verdict is
  # :good with the :island_partition note (the pre-island contract returned
  # :not_observable/:structural_islands here). The FD-aware tolerance is
  # still exercised by the local check below, whose column selection spans
  # the island-internal angle-offset null space.
  @testset "State estimation FD-aware rank tolerance" begin
    _lstd() = measurementStdDevs(vm = 1e-4, pinj = 1e-3, qinj = 1e-3, pflow = 1e-3, qflow = 1e-3)
    netC = deepcopy(fx.net_link)
    meas = generateMeasurementsFromPF(netC; noise = false, stddev = _lstd())

    # closed link (contracted): observable, no structural note
    obsC = evaluate_global_observability(netC, meas)
    @test obsC.quality == :good
    @test isempty(obsC.notes)

    # open link: a real separation; both islands are measured, each is
    # observable on its own subnet with its own reference, the merged
    # verdict reports the partition
    netO = deepcopy(fx.net_link_open)
    obsO = evaluate_global_observability(netO, [m for m in meas if m.linkIdx === nothing])
    @test obsO.quality == :good
    @test :island_partition in obsO.notes
    @test obsO.n_measured_islands == 2
    @test obsO.numerical_rank == obsO.n_states
    @test all(i.result.quality == :good for i in obsO.islands if i.measured)

    # an explicit tol still wins over the FD-aware default
    obsT = evaluate_global_observability(netC, meas; tol = 1e-14)
    @test obsT.numerical_rank == obsT.n_states

    # local-check rest case: a column selection spanning the ISLAND-INTERNAL
    # theta pair of the second island (B1a, B2; slack sits in island 1)
    # carries the common angle-offset null space in its own submatrix. The
    # FD-aware tolerance must fail it instead of silently passing on FD
    # noise. (Buses: S=1, B1=2, B1a=3, B2=4; theta columns 1..3 for the
    # non-slack buses, so island 2 is columns [2, 3].)
    obsL = evaluate_local_observability(netO, [m for m in meas if m.linkIdx === nothing], [2, 3])
    @test !obsL.numerical_observable
    @test obsL.quality == :not_observable
  end

  return true
end

function test_observability_takahashi_diagnostics(fx)::Bool
  # Takahashi task (0.10.0): the SE diagnostics Omega_ii/wii/rn via the
  # shared selected inverse agree with the dense pinv path; K requests and
  # guard failures fall back to dense; FD zeros are exact.
  # helper: H, r, w at the solved SE state (optionally with shunt B states,
  # the same construction validate_measurements uses)
  function _diag_inputs(net, meas; robust::Bool = false)
    res = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_net = false, robust = robust) do
      runse!(net, meas)
    end
    @assert res.converged
    prep = Sparlectra._se_prepare(net, meas)
    snet = prep.snet
    am = prep.meas
    nb = length(snet.nodeVec)
    slackIdx = Sparlectra._find_slack_idx(snet)
    releasedRows = res.shuntEstimates === nothing ? NamedTuple[] : NamedTuple[rw for rw in res.shuntEstimates if !rw.frozen]
    bInit = Float64[rw.B_est for rw in releasedRows]
    x = Sparlectra._initial_state_vector(snet, slackIdx; flatstart = false, shuntBInit = bInit)
    # res.voltages already lives in the contracted (snet) index space
    p = 0
    for i = 1:nb
      if i != slackIdx
        p += 1
        x[p] = angle(res.voltages[i])
      end
      x[nb-1+i] = abs(res.voltages[i])
    end
    Ybus = createYBUS(net = snet)
    shuntMap = nothing
    if !isempty(releasedRows)
      sidxs = Int[]
      sbuses = Int[]
      for rw in releasedRows
        sb = prep.busmap[prep.reps[rw.busIdx]]
        push!(sidxs, snet.shuntDict[sb])
        push!(sbuses, sb)
        sh = snet.shuntVec[snet.shuntDict[sb]]
        Ybus[sb, sb] -= sh.y_pu_shunt
      end
      shuntMap = Sparlectra.ShuntStateMap(sidxs, sbuses, (nb - 1) + nb + 1)
    end
    H, h = Sparlectra._measurement_jacobian_fd(am, snet, x, slackIdx, nb, Ybus; shuntMap = shuntMap)
    z = Sparlectra._measurement_vector(am)
    w = Sparlectra._weight_vector(am)
    return H, z - h, w
  end

  function _assert_paths_agree(H, r, w)
    d = Sparlectra._residual_diagnostics(H, r, w; minStates = typemax(Int))
    t = Sparlectra._residual_diagnostics(H, r, w; minStates = 1)
    @test d.omega_path == :dense
    @test t.omega_path == :takahashi
    @test t.omega === nothing
    @test maximum(abs.(d.wii .- t.wii) ./ max.(abs.(d.wii), 1e-12)) < 1e-10
    @test maximum(abs.(d.rn .- t.rn) ./ max.(abs.(d.rn), 1e-8)) < 1e-8
    @test maximum(abs.(d.state_variances .- t.state_variances) ./ abs.(d.state_variances)) < 1e-10
    # suspicion order identical
    @test sortperm(abs.(d.rn); rev = true) == sortperm(abs.(t.rn); rev = true)
  end

  @testset "State estimation Takahashi diagnostics" begin
    stdn = measurementStdDevs(vm = 0.002, pinj = 0.5, qinj = 0.5, pflow = 0.5, qflow = 0.5, shuntq = 0.5)

    # 1) small meshed case
    net1 = deepcopy(fx.net3)
    meas1 = generateMeasurementsFromPF(net1; noise = true, stddev = stdn, rng = MersenneTwister(1))
    H1, r1, w1 = _diag_inputs(net1, meas1)
    _assert_paths_agree(H1, r1, w1)

    # FD sparsity: exact zeros exactly where the structure demands them.
    # A Vm measurement depends only on its own Vm column; the flow on the
    # branch S-A does not touch bus B at all.
    am1 = Sparlectra._active_measurements(meas1)
    vmRow = findfirst(m -> m.typ == Sparlectra.VmMeas && m.busIdx == 2, am1)
    nb = 3
    colVm(i) = (nb - 1) + i
    @test H1[vmRow, colVm(2)] != 0.0
    @test H1[vmRow, colVm(1)] == 0.0 && H1[vmRow, colVm(3)] == 0.0
    @test H1[vmRow, 1] == 0.0 && H1[vmRow, 2] == 0.0   # no theta dependence
    @test count(iszero, H1) > 0

    # 2) released shunt plus link contraction
    net2 = deepcopy(fx.net_link_shunt)
    meas2 = generateMeasurementsFromPF(net2; includeShuntQ = true, noise = true, stddev = stdn, rng = MersenneTwister(2))
    setShuntEstimation!(net2; busName = "B2")
    H2, r2, w2 = _diag_inputs(net2, meas2)
    _assert_paths_agree(H2, r2, w2)

    # 3) robust-mode run (statistics stay on original sigmas either way)
    net3 = deepcopy(fx.net3)
    meas3 = generateMeasurementsFromPF(net3; noise = true, stddev = stdn, rng = MersenneTwister(3))
    b3 = findfirst(m -> m.typ == Sparlectra.PflowMeas, meas3)
    bm = meas3[b3]
    meas3[b3] = Measurement(typ = bm.typ, value = bm.value + 10.0 * bm.sigma, sigma = bm.sigma, active = bm.active, busIdx = bm.busIdx, branchIdx = bm.branchIdx, direction = bm.direction, id = bm.id)
    H3, r3, w3 = _diag_inputs(net3, meas3; robust = true)
    _assert_paths_agree(H3, r3, w3)

    # 4) guard failure falls back to dense with one warning, never a wrong
    # number: a zero column makes G singular (factorization guard)
    Hbad = [1.0 0.0; 1.0 0.0]
    rbad = [0.1, -0.1]
    wbad = [1.0, 1.0]
    fb = Test.@test_logs (:warn, r"falling back to the dense path") match_mode = :any Sparlectra._residual_diagnostics(Hbad, rbad, wbad; minStates = 1)
    @test fb.omega_path == :dense
    # unsymmetric pivot ordering guard on the shared module directly
    fake = (p = [1, 2], q = [2, 1], L = SparseArrays.sparse(1.0 * LinearAlgebra.I, 2, 2), U = SparseArrays.sparse(1.0 * LinearAlgebra.I, 2, 2), Rs = [1.0, 1.0])
    S, ok, info = takahashi_selected_inverse(fake)
    @test !ok && S === nothing && occursin("unsymmetric", info)
    dg, okd, infod = takahashi_diag(fake)
    @test !okd && isempty(dg) && occursin("unsymmetric", infod)

    # 5) K request pins the dense path, results unchanged
    repK = with_state_estimation_config(max_iter = 40, tol = 1e-8, report_residual_correlation = true) do
      validate_measurements(net3, meas3)
    end
    @test repK.omega_path == :dense
    @test all(!isnan(row.max_abs_correlation) for row in repK.measurement_ranking)
    # state variances are reported on both paths
    @test repK.state_variances !== nothing && all(v -> v >= 0.0 || isnan(v), repK.state_variances)

    # 6) the Jacobian arrives sparse, and the LDLt second
    # attempt hands the shared pattern pass a valid symmetric-pivot
    # factorization (this is the route that rescues symmetric G whose
    # UMFPACK auto strategy pivots unsymmetrically; found on
    # case1354pegase in the step-0 baseline)
    @test H3 isa SparseArrays.SparseMatrixCSC{Float64}
    Gsym = SparseArrays.sparse(LinearAlgebra.Symmetric(H3' * (LinearAlgebra.Diagonal(w3) * H3)))
    tak_in = Sparlectra._takahashi_ldlt_input(Gsym)
    @test tak_in !== nothing
    @test tak_in.p == tak_in.q
    @test all(isone, LinearAlgebra.diag(tak_in.L))
    @test maximum(abs.(Matrix(tak_in.L * tak_in.U) .- Matrix(Gsym[tak_in.p, tak_in.p]))) < 1e-12 * max(1.0, maximum(abs, Gsym))
    Sl, okl, infol = takahashi_selected_inverse(tak_in)
    @test (okl, infol) == (true, "ok")
    GI_ref = inv(Matrix(Gsym))
    @test all(abs(Sl[j, j] - GI_ref[j, j]) < 1e-8 for j in axes(Gsym, 1))

    # 6b) FD column coloring: the colored assembly
    # must be BIT-IDENTICAL to the per-column assembly at the same state
    # (same perturbation, same per-entry formula, same write order); any
    # deviation means the coloring changed the numerics and the feature
    # must stop. Probed on the shipped cases and, cache-gated, on
    # case1354pegase.
    fd_probe = function (net)
      prep = Sparlectra._se_prepare(net, Sparlectra.Measurement[mm for mm in net.measurements])
      snet = prep.snet
      meas = prep.meas
      nbus = length(snet.nodeVec)
      slackIdx = Sparlectra._find_slack_idx(snet)
      Ybus = createYBUS(net = snet)
      wva = Sparlectra._va_offset_active(meas, :auto)
      x = Sparlectra._initial_state_vector(snet, slackIdx; flatstart = false, withVaOffset = wva)
      H1, _ = Sparlectra._measurement_jacobian_fd(meas, snet, x, slackIdx, nbus, Ybus; withVaOffset = wva)
      # the PRODUCTION coloring is the structural one (point patterns are
      # not stable, see the report); the colored assembly must use it
      # without a single detector fallback
      col = Sparlectra._fd_structural_coloring(meas, length(x), slackIdx, nbus, Ybus, wva, nothing, nothing, snet)
      @test col !== nothing
      H2, _, used = Sparlectra._measurement_jacobian_fd(meas, snet, x, slackIdx, nbus, Ybus; withVaOffset = wva, coloring = col)
      @test used === true
      # superset inclusion: every numerically occupied
      # position must lie inside the structural pattern, probed at the
      # warm point and at a deterministically shifted one; a position
      # outside means the coupling table is incomplete
      for xp in (x, x .+ 0.01 .* sin.(1.0:length(x)))
        Hp, _ = Sparlectra._measurement_jacobian_fd(meas, snet, xp, slackIdx, nbus, Ybus; withVaOffset = wva)
        rvp = SparseArrays.rowvals(Hp)
        inside = true
        for k in axes(Hp, 2)
          for ptr in SparseArrays.nzrange(Hp, k)
            insorted(rvp[ptr], col.rows_per_col[k]) || (inside = false)
          end
        end
        @test inside
      end
      return H1, H2, col
    end
    for demo in ("sp_case60", "sp_case188")
      netd = Sparlectra.importSCF(joinpath(dirname(@__DIR__), "data", "scf", demo * ".scf.json"))
      H1, H2, col = fd_probe(netd)
      @test (demo, H2.colptr == H1.colptr, H2.rowval == H1.rowval, H2.nzval == H1.nzval) == (demo, true, true, true)
      ncolors = length(col.groups)
      println("      fd coloring ", demo, ": ", ncolors, " colors for ", size(H1, 2), " states")
      @test ncolors < size(H1, 2) / 5
    end
    # size probe: the shipped data/mpower/sp_case1354.m (synthetic 1354-bus
    # grid, not shipped generator) replaces the cache-gated case1354pegase
    big1354 = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "sp_case1354.m"))
    @test isfile(big1354)
    begin
      netp = Sparlectra.createNetFromMatPowerFile(filename = big1354, flatstart = false, bus_shunt_model = :admittance, matpower_shift_sign = 1.0, matpower_shift_unit = :deg, matpower_ratio = :normal, tap_changer_model = :ideal)
      runpf!(netp, 60, 1e-8, 0)
      msp = generateMeasurementsFromPF(netp; includeVm = true, includePinj = true, includeQinj = true, includePflow = false, includeQflow = false, noise = true, rng = Xoshiro(42))
      empty!(netp.measurements)
      append!(netp.measurements, msp)
      H1, H2, col = fd_probe(netp)
      @test (H2.colptr == H1.colptr, H2.rowval == H1.rowval, H2.nzval == H1.nzval) == (true, true, true)
      println("      fd coloring sp_case1354: RAN (", length(col.groups), " colors for ", size(H1, 2), " states)")
      @test length(col.groups) < size(H1, 2) / 5
    end

    # 7) criticality (issue #394): the default reads diag(Omega) from one
    # selected-inverse pass and knows no budget; a 600 x 600 identity has
    # every row critical (each state seen by exactly one row). The former
    # per-row rank tests stay as criticality_method = :rank, budgeted and
    # skipped above the budget with the warning, as before.
    Hwide = SparseArrays.sparse(1.0 * LinearAlgebra.I, 600, 600)
    big = Sparlectra.evaluate_observability_matrix(Hwide)
    @test big.criticality_method === :omega
    @test !big.criticality_skipped
    @test length(big.numerical_critical_measurement_indices) == 600
    @test all(<=(1e-8), big.criticality_wii)
    big_rank = Test.@test_logs (:warn, r"single-row criticality skipped") match_mode = :any Sparlectra.evaluate_observability_matrix(Hwide; criticality_method = :rank)
    @test big_rank.criticality_skipped
    @test isempty(big_rank.numerical_critical_measurement_indices)
    small = Sparlectra.evaluate_observability_matrix(SparseArrays.sparse(1.0 * LinearAlgebra.I, 3, 3); criticality_method = :rank)
    @test !small.criticality_skipped
    # a state column no row touches is a not-observable verdict, not an
    # error (the estimator's tap and shunt release guards rely on it; a
    # thinned set left a released tap without rows, Web UI run 93b08476)
    Hgap = SparseArrays.sparse([1, 2], [1, 2], [1.0, 1.0], 2, 3)
    gap = Sparlectra.evaluate_local_observability_matrix(Hgap, [3])
    @test !gap.numerical_observable
    @test isempty(gap.rows)
    # both methods agree on a real set at the boundary: sp_case5 with one
    # flow row removed, the rows the rank tests call critical are exactly
    # the rows with Omega_ii at zero, and the redundant rows carry wii > 0
    net5 = importSCF(joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.scf.json"))
    readMeasurementsCSV!(net5; file = joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.measurements.csv"))
    kflow = findfirst(m -> m.typ == Sparlectra.PflowMeas, net5.measurements)
    deleteat!(net5.measurements, kflow)
    obs_omega = with_state_estimation_config(() -> evaluate_global_observability(net5); criticality_method = :omega)
    obs_rank = with_state_estimation_config(() -> evaluate_global_observability(net5); criticality_method = :rank)
    @test obs_omega.criticality_method === :omega && obs_rank.criticality_method === :rank
    @test sort(obs_omega.numerical_critical_measurement_indices) == sort(obs_rank.numerical_critical_measurement_indices)
    @test sort(obs_omega.structural_critical_measurement_indices) == sort(obs_rank.structural_critical_measurement_indices)
    @test obs_omega.quality == obs_rank.quality
    @test length(obs_omega.criticality_wii) == length(obs_omega.active_measurement_indices)
    redundant = [k for k in eachindex(obs_omega.active_measurement_indices) if !(obs_omega.active_measurement_indices[k] in obs_omega.numerical_critical_measurement_indices)]
    @test all(k -> obs_omega.criticality_wii[k] > 1e-8, redundant)
  end

  return true
end

function run_observability_tests()
  # Aggregates the observability test sets; the shared nets are built once
  # and handed to every set that needs them.
  @testset "Observability" begin
    fx = observability_fixture()
    tests =
      [("Observability metrics", () -> test_observability_metrics(fx)), ("Matrix observability helpers", test_observability_matrix_helpers), ("FD-aware rank tolerance", () -> test_observability_fd_rank_tolerance(fx)), ("Takahashi diagnostics", () -> test_observability_takahashi_diagnostics(fx))]

    for (name, testfn) in tests
      @testset "$name" begin
        @test _se_run_quiet(testfn) == true
      end
    end
  end
end
