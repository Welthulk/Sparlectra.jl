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
#          FD column coloring, diag(Omega) criticality, #394), the
#          observable islands of the P-theta model with the boundary
#          iteration, the deficits of the P-theta and Q-V models (SVD and
#          gain-matrix iteration, boundary rows, network islands), the
#          estimation of the observable islands one by one
#          (on_unobservable), the measurement placement and the weak
#          class against an SVD reference

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
shunt on B2 (`net_link_shunt`), and for the observable-island sets the
solved shipped IEEE 14-bus topology with its slack on bus 8 (`net14z`,
case14_zeng_p306_activeSet_A), the solved sp_case60 (`net60`) and the solved
two-island two_islands_prio (`net_two`). Only functions that take the measurement
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
  scf = joinpath(pkgdir(Sparlectra), "data", "scf")
  net14z = importSCF(joinpath(scf, "case14_zeng_p306_activeSet_A.scf.json"))
  runpf!(net14z, 40, 1e-8, 0; method = :rectangular)
  net60 = importSCF(joinpath(scf, "sp_case60.scf.json"))
  runpf!(net60, 40, 1e-8, 0; method = :rectangular)
  net_two = redirect_stdout(() -> importSCF(joinpath(scf, "two_islands_prio.scf.json")), devnull)
  redirect_stdout(() -> runpf!(net_two, 40, 1e-8, 0; method = :rectangular, islands_enabled = true), devnull)
  return (net3 = net3, pf3 = (ite = ite3, erg = erg3), net14 = net14, net_link = net_link, net_link_open = net_link_open, net_link_shunt = net_link_shunt, net14z = net14z, net60 = net60, net_two = net_two)
end

function test_observability_metrics(fx)::Bool
  # Validates global/local observability metrics for full, reduced, and sparse
  # measurement sets, including expected not-observable behavior for sparse data.
  @testset "State estimation observability metrics" begin (function ()
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
        meas_reduced[i] = Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = false, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id, origin = m.origin)
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
      meas_sparse[i] = Measurement(typ = m.typ, value = m.value, sigma = m.sigma, active = keep, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id, origin = m.origin)
    end

    gobs_sparse = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
      evaluate_global_observability(net, meas_sparse)
    end
    @test gobs_sparse.quality == :not_observable
    @test gobs_sparse.numerical_observable == false

    @test_throws ErrorException with_state_estimation_config(() -> evaluate_local_observability(net, meas_sparse, local_cols); flatstart = true, jac_eps = 1e-6)
  end)() end

  # The rank test runs on a column-normalized Jacobian. Without it the
  # relative tolerance is measured against the largest column scale, and on
  # a 25000-bus network a set with 203916 rows for 49999 states came out as
  # not observable (deficit 191 raw). Rows are deliberately not scaled by
  # sigma: rank is invariant under positive row scaling, and dividing by it
  # lifts zero-injection rows (sigma 1e-6) six orders of magnitude above the
  # rest, which made the same case worse (deficit 12308). Column
  # normalization alone leaves it at 0.
  # An objection from review, and it was justified: the interesting failure is
  # not a set without any flow information (there the matching fails
  # obviously), it is one where the Hopcroft-Karp matching covers every
  # state column while information is still missing, so two angles are only
  # jointly determined. An earlier version let the structural result
  # overrule the numerical one and waved exactly those through as
  # observable. Thinning sp_case14 produces them: five of sixty draws had
  # full matching with a rank one or two short. The property below is what
  # must hold - a numerical deficit is never reported as observable, no
  # matter what the structure says.
  @testset "observability: structure never overrules a numerical deficit" begin (function ()
    net14 = deepcopy(fx.net14)
    full = generateMeasurementsFromPF(net14; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false)
    rng = MersenneTwister(11)
    seen_gap = false
    for _ = 1:40
      keep = Measurement[m for m in full if rand(rng) < 0.55]
      isempty(keep) && continue
      obs = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
        evaluate_global_observability(net14, keep)
      end
      if obs.numerical_rank < obs.n_states
        @test obs.quality == :not_observable
        obs.structural_matching == obs.n_states && (seen_gap = true)
      end
    end
    # the draws must actually contain the interesting case, otherwise the
    # test above proves nothing
    @test seen_gap
  end)() end

  @testset "observability: rank verdict is scale independent" begin (function ()
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
    scaled = [Measurement(typ = m.typ, value = m.value, sigma = m.sigma * 1000, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id, origin = m.origin) for m in msH]
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
    # a bus coupler entered as a branch (r = 0, x = 1e-5 pu, admittance 1e5
    # against about 600 for the feeding lines, the case141 pattern): its rows
    # dominated both bus columns, the column normalization alone shrank the
    # feeding line's entries below the cut and a fully measured radial
    # net read rank 6 of 7, "not observable"
    netC = Net(name = "coupler", baseMVA = 10.0)
    for i = 1:4
      addBus!(net = netC, busName = "B$(i)", vn_kV = 12.47, vm_pu = 1.0, va_deg = 0.0)
    end
    addPIModelACLine!(net = netC, fromBus = "B1", toBus = "B2", r_pu = 0.0037, x_pu = 0.0016, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = netC, fromBus = "B2", toBus = "B3", r_pu = 0.0037, x_pu = 0.0016, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = netC, fromBus = "B3", toBus = "B4", r_pu = 0.0, x_pu = 1e-5, b_pu = 0.0, status = 1)
    addProsumer!(net = netC, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.0, va_deg = 0.0)
    for (b, p, q) in [("B2", 0.3, 0.1), ("B3", 0.5, 0.2), ("B4", 0.15, 0.05)]
      addProsumer!(net = netC, busName = b, type = "LOAD", p = p, q = q)
    end
    okC, msgC = validate!(net = netC)
    okC || error("test net invalid: $msgC")
    _, ergC = runpf!(netC, 40, 1e-10, 0)
    @test ergC == 0
    msC = generateMeasurementsFromPF(netC; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false)
    obsC = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
      evaluate_global_observability(netC, msC)
    end
    @test obsC.numerical_observable
    @test obsC.numerical_rank == obsC.n_states
  end)() end

  return true
end

function test_observability_matrix_helpers()::Bool
  # Unit-checks matrix-only observability helpers (rank/matching/redundancy),
  # including local-column selection and expected error paths.
  @testset "Unobservable state columns (null-space dark states)" begin (function ()
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
  end)() end

  @testset "State estimation matrix observability helpers" begin (function ()
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
  end)() end

  return true
end

function test_observability_fd_rank_tolerance(fx)::Bool
  # Island-wise SE (0.10.0): two-stage observability.
  # An OPEN link is a real separation; with island-wise estimation each
  # measured island is judged on its own subnet and the merged verdict is
  # :good with the :island_partition note. The FD-aware tolerance is
  # still exercised by the local check below, whose column selection spans
  # the island-internal angle-offset null space.
  @testset "State estimation FD-aware rank tolerance" begin (function ()
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
  end)() end

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
    # the blockwise route (above the dense state cap) computes the same
    # Omega diagonal and the same K maxima from the sparse factorization
    b = Sparlectra._residual_diagnostics_blockwise(H, r, w)
    @test b.omega_path == :blockwise
    @test maximum(abs.(d.wii .- b.wii) ./ max.(abs.(d.wii), 1e-12)) < 1e-10
    @test maximum(abs.(d.state_variances .- b.state_variances) ./ abs.(d.state_variances)) < 1e-10
    kd = Sparlectra._residual_correlation_max(d.omega)   # the dense path always carries the full Omega
    @test maximum(abs.(kd .- Sparlectra._residual_correlation_max_blockwise(H, w))) < 1e-8
  end

  @testset "State estimation Takahashi diagnostics" begin (function ()
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
    meas3[b3] = Measurement(typ = bm.typ, value = bm.value + 10.0 * bm.sigma, sigma = bm.sigma, active = bm.active, busIdx = bm.busIdx, branchIdx = bm.branchIdx, direction = bm.direction, id = bm.id, origin = bm.origin)
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
      # the row-restricted assembly the WLS loop runs after its first
      # iteration (only the rows a color's columns touch are re-evaluated)
      # must give the same H bit for bit
      H3, _, used3 = Sparlectra._measurement_jacobian_fd(meas, snet, x, slackIdx, nbus, Ybus; withVaOffset = wva, coloring = col, check_pattern = false)
      @test used3 === true
      @test (H3.colptr == H2.colptr, H3.rowval == H2.rowval, H3.nzval == H2.nzval) == (true, true, true)
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
      # a P and a Q row on the coupler of the case (0.60.0): link rows in
      # the colored assembly, bit-identical and inside the pattern as well
      couplers = [k for (k, l) in enumerate(netd.linkVec) if l.status == 1]
      @test !isempty(couplers)
      for k in couplers
        addPflowMeasurement!(netd; value = 0.0, sigma = 1.0, linkNr = k)
        addQflowMeasurement!(netd; value = 0.0, sigma = 1.0, linkNr = k)
      end
      H1, H2, col = fd_probe(netd)
      @test (demo, H2.colptr == H1.colptr, H2.rowval == H1.rowval, H2.nzval == H1.nzval) == (demo, true, true, true)
      ncolors = length(col.groups)
      println("      fd coloring ", demo, ": ", ncolors, " colors for ", size(H1, 2), " states")
      @test ncolors < size(H1, 2) / 5
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
    # thinned set left a released tap without rows, a Web UI run)
    Hgap = SparseArrays.sparse([1, 2], [1, 2], [1.0, 1.0], 2, 3)
    gap = Sparlectra.evaluate_local_observability_matrix(Hgap, [3])
    @test !gap.numerical_observable
    @test isempty(gap.rows)
    # both methods agree on a real set at the boundary: sp_case5 with one
    # flow row removed, the rows the rank tests call critical are exactly
    # the rows with Omega_ii at zero, and the redundant rows carry wii > 0
    net5 = importSCF(joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.scf.json"))
    readMeasurementTable!(net5; file = joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.measurements.csv"), zero_injection = false)
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
  end)() end

  return true
end

## Measurement sets of the observable-island tests (concept 1.1.0, block A):
## the full noise-free set of a solved net and filters on it.
_obs_full_set(net) = generateMeasurementsFromPF(net; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false)
_obs_is_injection(m) = m.typ in (Sparlectra.PinjMeas, Sparlectra.QinjMeas)
_obs_cfg(f) = with_state_estimation_config(f; flatstart = true, jac_eps = 1e-6)
_obs_mode(f, mode::Symbol; zero_injection::Symbol = :on) = with_state_estimation_config(f; flatstart = true, jac_eps = 1e-6, zero_injection = zero_injection, observability = Sparlectra.ObservabilityConfig(on_unobservable = mode))
## the schedule and U-order tests without the passive-bus candidates
## (zero_injection = off reports them in `passive_not_inserted`)
_obs_zioff(f) = with_state_estimation_config(f; flatstart = true, jac_eps = 1e-6, zero_injection = :off)

## no observable island may contradict the region: the angle of every bus
## outside the observable island of the reference bus floats against the
## reference, so it has a share in the P-theta null space (a region bus).
## The reverse does not hold: a set of injections without flows spreads the
## null space over buses the iteration keeps in the reference island. Per
## network island the observable islands are numbered in the order of
## their first bus and labeled "k.j" ("k" for one), and the island numbers
## run 1..n_islands.
function _obs_islands_consistent(net, isl)::Bool
  slack = findfirst(nd -> getNodeType(nd) == Sparlectra.Slack, net.nodeVec)
  ref = isl.island_of_bus[slack]
  floating = [b for b in eachindex(isl.island_of_bus) if isl.island_of_bus[b] != 0 && isl.island_of_bus[b] != ref && isl.network_island_of_bus[b] == isl.network_island_of_bus[slack]]
  ok = issubset(floating, isl.region_buses) && sort(unique(filter(!=(0), isl.island_of_bus))) == collect(1:isl.n_islands)
  for row in isl.network_islands
    row.measured || continue
    ids = [id for id in eachindex(isl.island_buses) if isl.network_island_of_bus[first(isl.island_buses[id])] == row.network_island]
    k = row.network_island
    labels = length(ids) == 1 ? [string(k)] : [string(k, ".", j) for j in eachindex(ids)]
    ok &= issorted([first(isl.island_buses[id]) for id in ids]) && isl.island_labels[ids] == labels && row.n_islands == length(ids)
  end
  return ok
end

## rank and dark columns of a set recomputed without the observability
## code: the labeled Jacobian, rows then columns normalized to unit norm
## (a zero row or column stays zero), a dense SVD, the FD-aware rank
## tolerance and the 1.0.0 dark rule (null-vector entries above 1e-8)
function _obs_svd_reference(net, meas)
  netc = deepcopy(net)
  empty!(netc.measurements)
  append!(netc.measurements, meas)
  H = Matrix{Float64}(_obs_cfg(() -> measurement_jacobian(netc)).H)
  for i in axes(H, 1)
    r = LinearAlgebra.norm(view(H, i, :))
    r > 0.0 && (H[i, :] ./= r)
  end
  for j in axes(H, 2)
    c = LinearAlgebra.norm(view(H, :, j))
    c > 0.0 && (H[:, j] ./= c)
  end
  F = LinearAlgebra.svd(H; full = true)
  cfg = _obs_cfg(state_estimation_config)
  tol = cfg.rank_tol_factor * cfg.jac_eps * F.S[1]
  r = count(>(tol), F.S)
  # dark: the row of the orthonormal null basis above tol / sigma_max, at
  # least 1e-8 (concept D1; 1.0.0 took every entry above 1e-8)
  tau = max(1.0e-8, tol / F.S[1])
  N = F.V[:, (r+1):end]
  dark = Int[j for j in axes(N, 1) if LinearAlgebra.norm(view(N, j, :)) > tau]
  return r, dark
end

## the dark buses of the report as state columns (angle of bus b in column
## b, b - 1 behind the slack; magnitude in column nb - 1 + b): the verdict
## field `unobservable_state_columns` and the report read the same null
## bases of the decoupled models with the same tolerance (concept D1, G4')
function _obs_report_columns(net, isl)::Vector{Int}
  nb = length(net.nodeVec)
  slack = findfirst(nd -> getNodeType(nd) == Sparlectra.Slack, net.nodeVec)
  return sort(vcat([b < slack ? b : b - 1 for b in isl.region_buses], [nb - 1 + b for b in isl.dark_magnitude_buses]))
end

function test_observability_islands(fx)::Bool
  # Test 1 of the concept. IEEE 14-bus topology, slack on bus 8, bus 8
  # hangs on bus 7 only (branch 14, 7-8). Removed: both flows of 7-8 and the
  # injections at 7 and 8, i.e. every row that contains the angle of bus 8
  # together with another angle. By hand: every remaining row depends on
  # angle differences inside the other 13 buses only, so those 13 angles
  # can shift together against the reference without changing a
  # measurement. One null direction of the P-theta model, two observable
  # islands {8} and the rest (labels 1.1 and 1.2 of the one network
  # island), branch 7-8 unobservable, the 13 non-slack angles in the region,
  # no magnitude dark (every bus keeps its Vm row, the Q-V model has full
  # rank), nothing to remove (the injections at the ends of 7-8 are already
  # gone). The verdict field `unobservable_state_columns` names the same 13
  # angle columns: 1.0.0 counted every SVD entry above 1e-8 and listed nine
  # magnitude columns as well, whose null-space entries are
  # forward-difference noise (at most 3.7e-8, against 0.21 to 0.35 for the
  # 13 angles, concept D1 and D2).
  @testset "observable islands: case14 split at branch 7-8, slack alone on bus 8" begin (function ()
    net = fx.net14z
    full = _obs_full_set(net)
    meas = filter(m -> m.branchIdx != 14 && !(_obs_is_injection(m) && m.busIdx in (7, 8)), full)
    obs = _obs_cfg(() -> evaluate_global_observability(net, meas))
    isl = obs.observable_islands
    @test obs.quality == :not_observable
    @test isl.n_islands == 2
    @test isl.island_of_bus == [1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1]
    @test isl.network_island_of_bus == fill(1, 14)
    @test isl.island_buses[2] == [8] && isl.island_labels == ["1.1", "1.2"]
    @test isl.unobservable_branches == [14]
    @test isl.unobservable_branch_names == ["7 - 8"]
    @test (isl.rank_deficit, isl.angle_rank_deficit, isl.magnitude_rank_deficit) == (1, 1, 0)
    @test (obs.rank_deficit, obs.angle_rank_deficit, obs.magnitude_rank_deficit) == (1, 1, 0)
    @test isl.region_buses == [1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14]
    @test isempty(isl.dark_magnitude_buses)
    @test isempty(isl.removed_boundary_injections) && isl.passes == 1
    @test isempty(isl.unmeasured_network_islands) && length(isl.network_islands) == 1
    @test _obs_islands_consistent(net, isl)
    @test obs.unobservable_state_columns == _obs_report_columns(net, isl) == collect(1:13)
    # not observable, so the weak class is not determined
    @test obs.weak_observability.basis === :not_determined && !obs.weakly_observable
    # the standalone call gives the same report
    @test _obs_cfg(() -> identify_observable_islands(net, meas)).island_of_bus == isl.island_of_bus
    io = IOBuffer()
    print_observable_islands(io, obs)
    @test occursin("7 - 8", String(take!(io)))
    # counter-check: a zero injection at bus 7 is the balance of the three
    # branches at 7; with the flows on 4-7 and 7-9 measured it fixes the
    # flow on 7-8, and the set is one island
    meas_zi = copy(meas)
    addZeroInjectionMeasurements!(meas_zi; net = net, busIdxs = [7])
    obs_zi = _obs_cfg(() -> evaluate_global_observability(net, meas_zi))
    @test obs_zi.numerical_observable
    @test isempty(obs_zi.observable_islands)
    @test _obs_cfg(() -> identify_observable_islands(net, meas_zi)).n_islands == 1
  end)() end

  # Test 10 of the concept (I6): runse! with on_unobservable =
  # observable_islands on the set of test 1. Observable island 1.2 is bus 8
  # alone with the slack and its Vm row; island 1.1 takes the bus of its
  # best voltage-controlled unit as reference (bus 1, the largest unit of
  # case14), so its angles are relative. The flow on 7-8 follows from no
  # measurement and stays NaN, not zero. Each island equals runse! on its
  # subnet, cut out with the helpers of the network-island path, with the
  # rows that lie inside it; with noise-free rows both reproduce the power
  # flow, because a reference keeps its stored angle. Two variants:
  # without the Vm rows of island 1.1 its magnitudes can shift together
  # (the Q-V null direction), all 13 are dark and leave (concept M1: they
  # are listed as removed, the island has no row, bus 8 is island 1.1 of
  # the cut); with the Q injection at 7 kept and the Vm row, the Q flows
  # at 7 and the Q injections at 4 and 9 (which would tie it through 4-7
  # and 7-9) removed, the magnitude of 7 is seen through that injection
  # only (not dark), but the injection holds the flow on 7-8 and is not
  # used, so island 1.1 fails the rank check (Q-V deficit 1). Under :abort
  # the not-observable set takes the old path (no island rows), an
  # observable set gives the same result bit for bit in both modes.
  @testset "observable islands: estimated one by one" begin (function ()
    net = fx.net14z
    full = _obs_full_set(net)
    meas = filter(m -> m.branchIdx != 14 && !(_obs_is_injection(m) && m.busIdx in (7, 8)), full)
    closed = [k for (k, b) in enumerate(net.branchVec) if Sparlectra._branch_terminal_state(b) == :closed]
    est = deepcopy(net)
    res = _obs_mode(() -> runse!(est, meas), :observable_islands)
    oi = res.observableIslands
    @test res.converged
    @test [(r.label, r.network_island, r.estimated, r.reason) for r in res.islands] == [("1.1", 1, true, :estimated), ("1.2", 1, true, :estimated)]
    @test res.islands[1].buses == [1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14] && res.islands[2].buses == [8]
    @test (res.islands[1].reference_bus, res.islands[1].relative_angles) == (1, true)
    @test (res.islands[2].reference_bus, res.islands[2].relative_angles) == (8, false)
    @test oi.label_of_bus == [fill("1.1", 7); "1.2"; fill("1.1", 6)]
    @test oi.relative_angle == [fill(true, 7); false; fill(true, 6)] && all(oi.estimated)
    @test oi.unobservable_branches == [14] && isempty(oi.dropped_measurements)
    @test isnan(est.branchVec[14].fBranchFlow.pFlow) && isnan(est.branchVec[14].tBranchFlow.qFlow)
    @test all(k -> isfinite(est.branchVec[k].fBranchFlow.pFlow), setdiff(closed, [14]))
    @test all(isnan, est.totalLosses[end])
    @test maximum(abs.(res.voltages .- buildVoltageVector(net))) < 1e-9
    for row in res.islands
      bset = Set(row.buses)
      brs = [k for k in closed if Int(net.branchVec[k].fromBus) in bset && Int(net.branchVec[k].toBus) in bset && !(k in oi.unobservable_branches)]
      slack = [b for b in row.buses if getNodeType(net.nodeVec[b]) == Sparlectra.Slack]
      sub = Sparlectra._se_island_subnet(net, (buses = row.buses, branches = brs, n_bus = length(row.buses), n_ref = length(slack), chosen_ref_bus = isempty(slack) ? row.reference_bus : 0))
      inside = [q for (q, m) in enumerate(meas) if m.branchIdx !== nothing ? m.branchIdx in brs : m.busIdx in bset]
      smeas, _ = Sparlectra._se_remap_island_measurements(sub, meas, inside)
      sres = _obs_cfg(() -> runse!(sub.inet, smeas))
      @test sres.converged && maximum(abs.(sres.voltages .- res.voltages[row.buses])) < 1e-12
      @test (sres.objectiveJ, sres.dof) == (row.objectiveJ, row.dof)
    end
    io = IOBuffer()
    print_observable_islands(io, res)
    @test length(split(strip(String(take!(io))), '\n')) == 3
    # the weak class at the estimate, per estimated island (concept M2):
    # both islands carry every row of their buses, nothing is weak
    @test res.weakObservability.basis === :estimated_state && res.weakObservability.n_buses == 0
    @test [r.label for r in res.weakObservability.network_islands] == ["1.1", "1.2"]

    # no magnitude anchor in island 1.1
    noVm = filter(m -> !(m.typ == Sparlectra.VmMeas && m.busIdx != 8), meas)
    est2 = deepcopy(net)
    res2 = _obs_mode(() -> runse!(est2, noVm), :observable_islands)
    @test [(r.label, r.estimated, r.reason, r.buses) for r in res2.islands] == [("1.1", true, :estimated, [8])]
    @test res2.observableIslands.removed_buses == [(bus = b, reason = :dark_magnitude) for b in res.islands[1].buses]
    @test res2.observableIslands.label_of_bus == [fill("", 7); "1.1"; fill("", 6)]
    @test all(isnan, res2.voltages[res.islands[1].buses]) && all(b -> isnan(est2.nodeVec[b]._vm_pu), res.islands[1].buses)
    @test all(k -> isnan(est2.branchVec[k].fBranchFlow.pFlow), closed)

    # the magnitude of 7 seen through its injection only
    at7 = [k for k in closed if 7 in (Int(net.branchVec[k].fromBus), Int(net.branchVec[k].toBus))]
    q7 = filter(m -> m.branchIdx != 14 && !(_obs_is_injection(m) && m.busIdx == 8) && !(m.typ == Sparlectra.PinjMeas && m.busIdx == 7) && !(m.typ == Sparlectra.VmMeas && m.busIdx == 7) && !(m.typ == Sparlectra.QinjMeas && m.busIdx in (4, 9)) && !(m.typ == Sparlectra.QflowMeas && m.branchIdx in at7), full)
    isl7 = _obs_cfg(() -> identify_observable_islands(net, q7))
    @test isl7.island_labels == ["1.1", "1.2"] && isempty(isl7.dark_magnitude_buses) && isl7.magnitude_rank_deficit == 0
    res3 = _obs_mode(() -> runse!(deepcopy(net), q7), :observable_islands)
    @test [(r.label, r.estimated, r.reason, r.rank_deficit) for r in res3.islands] == [("1.1", false, :not_observable, 1), ("1.2", true, :estimated, 0)]
    @test "Qinj_bus_7" in res3.observableIslands.dropped_measurements

    # abort: the old path, and an observable set is the same in both modes
    ab = _obs_mode(() -> runse!(deepcopy(net), meas), :abort)
    @test ab.islands === nothing && ab.observableIslands === nothing && all(isfinite, ab.voltages)
    a = _obs_mode(() -> runse!(deepcopy(net), full), :abort)
    o = _obs_mode(() -> runse!(deepcopy(net), full), :observable_islands)
    @test isequal(a.voltages, o.voltages) && a.objectiveJ === o.objectiveJ && (a.iterations, a.dof) == (o.iterations, o.dof)
    @test o.islands === nothing && o.observableIslands === nothing
  end)() end

  # Test 2: buses 9 to 14 carry no measurement (no bus row and no flow
  # measured at their end of a branch). Flows measured at the far ends
  # (4-9, 7-9, 6-11, 6-12, 6-13) still fix 9, 11, 12 and 13; buses 10 and
  # 14 float. By hand in the decoupled models: the angles of 10 and 14 have
  # no P row (P-theta deficit 2), their magnitudes no Q or Vm row (Q-V
  # deficit 2), full deficit 4, here equal to the AC deficit and not the
  # island count minus one.
  @testset "observable islands: a larger gap" begin (function ()
    net = fx.net14z
    endbus(m) = Int(m.direction == :to ? net.branchVec[m.branchIdx].toBus : net.branchVec[m.branchIdx].fromBus)
    meas = filter(m -> !(m.busIdx !== nothing && m.busIdx in 9:14) && !(m.branchIdx !== nothing && endbus(m) in 9:14), _obs_full_set(net))
    obs = _obs_cfg(() -> evaluate_global_observability(net, meas))
    isl = obs.observable_islands
    @test isl.island_buses == [[1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13], [10], [14]]
    @test isl.island_labels == ["1.1", "1.2", "1.3"]
    @test isl.unobservable_branches == [16, 17, 18, 20]
    @test (isl.rank_deficit, isl.angle_rank_deficit, isl.magnitude_rank_deficit) == (4, 2, 2)
    @test isl.rank_deficit == obs.n_states - obs.numerical_rank
    @test isl.angle_rank_deficit == isl.n_islands - 1
    @test isl.region_buses == isl.dark_magnitude_buses == [10, 14]
    @test _obs_islands_consistent(net, isl)
    # the dark state columns of the verdict name the same buses: angle
    # column of bus b is b - 1 (slack 8 has none for b > 8), magnitude 13 + b
    @test obs.unobservable_state_columns == _obs_report_columns(net, isl) == [9, 13, 23, 27]
  end)() end

  # Test 3: a boundary injection. Removed: the flows on 4-7, 4-9, 5-6, 6-11,
  # 7-9, 9-14 and the injections at 4, 5, 6, 7, 11, 14. Left are the rigid
  # groups G1 = {1..5}, G2 = {6, 12, 13, 14} and G3 = {9, 10, 11} (measured
  # flows inside them) and {7, 8} with the slack. By hand on the P-theta
  # model: G1 and G2 shift freely by a and b, the P injection at bus 9 sees
  # the branches 4-9, 7-9 and 9-14 and fixes the shift of G3 to
  # (b49 a + b914 b) / (b49 + b79 + b914). Deficit 2, four groups with
  # different shifts in the first pass, the six branches between them
  # unobservable. P9 sits at their end buses and is dropped (Monticelli and
  # Wu); the second pass finds no new unobservable branch, so the observable
  # islands are the four groups. The region is every non-slack bus outside
  # {7, 8}. The AC rank (deficit 1) is weaker than the P-theta model here.
  @testset "observable islands: a boundary injection is dropped" begin (function ()
    net = fx.net14z
    meas = filter(m -> !(m.branchIdx !== nothing && m.branchIdx in (8, 9, 10, 11, 15, 17)) && !(_obs_is_injection(m) && m.busIdx in (4, 5, 6, 7, 11, 14)), _obs_full_set(net))
    obs = _obs_cfg(() -> evaluate_global_observability(net, meas))
    fin = obs.observable_islands
    @test obs.quality == :not_observable && !obs.numerical_observable
    @test (obs.angle_rank_deficit, obs.magnitude_rank_deficit, obs.n_states - obs.numerical_rank) == (2, 0, 1)
    @test fin.island_buses == [[1, 2, 3, 4, 5], [6, 12, 13, 14], [7, 8], [9, 10, 11]]
    @test fin.unobservable_branches == [8, 9, 10, 11, 15, 17]
    @test fin.removed_boundary_injections == ["Pinj_bus_9"] && fin.passes == 2
    @test fin.region_buses == [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14]
    @test (fin.angle_rank_deficit, fin.rank_deficit) == (2, 2)
    @test _obs_islands_consistent(net, fin)
    # the first pass alone finds the same groups and drops nothing
    one = _obs_cfg(() -> identify_observable_islands(net, meas; max_passes = 1))
    @test one.island_buses == fin.island_buses && isempty(one.removed_boundary_injections) && one.region_buses == fin.region_buses
    io = IOBuffer()
    print_observable_islands(io, fin)
    out = String(take!(io))
    @test occursin("observable island 1.4", out) && occursin("removed boundary rows: 1", out)
    @test_throws ErrorException identify_observable_islands(net, meas; max_passes = 0)
  end)() end

  # Test 4: bus 12 keeps its P injection and nothing else: no Q row and no
  # Vm row touches it (flows on 6-12 and 12-13 and the injections at 6 and
  # 13 removed as well). In the P-theta model the P row fixes the angle of
  # bus 12, in the Q-V model its magnitude has no row: one observable
  # island, deficit 0 + 1, the magnitude dark, the angle not (the AC null
  # vector mixes both, 790900bc listed the angle as dark as well).
  @testset "observable islands: a bus with P only is no island" begin (function ()
    net = fx.net14z
    meas = filter(m -> !(m.branchIdx !== nothing && m.branchIdx in (12, 19)) && !(m.busIdx == 12 && m.typ != Sparlectra.PinjMeas) && !(_obs_is_injection(m) && m.busIdx in (6, 13)), _obs_full_set(net))
    obs = _obs_cfg(() -> evaluate_global_observability(net, meas))
    isl = obs.observable_islands
    @test !obs.numerical_observable
    @test isl.n_islands == 1 && isl.island_labels == ["1"]
    @test isempty(isl.unobservable_branches)
    @test (isl.rank_deficit, isl.angle_rank_deficit, isl.magnitude_rank_deficit) == (1, 0, 1)
    @test isempty(isl.region_buses) && isl.dark_magnitude_buses == [12]
    @test _obs_islands_consistent(net, isl)
    @test obs.unobservable_state_columns == [25]
  end)() end

  # Test 4b (concept G4'): bus 12 is seen through one Q flow only, on 6-12
  # read at bus 6 (its bus rows, both P flows, the Q flow on 12-13 and the
  # injections at 6 and 13 removed). By hand: no P row contains the angle of
  # bus 12, the P-theta deficit is 1 (bus 12 its own observable island,
  # 6-12 and 12-13 unobservable); the Q flow b (V6 - V12) with V6 measured
  # fixes V12, the Q-V deficit is 0. Full deficit 1. The AC Jacobian has
  # one null direction as well (two states, one row); 790900bc added the
  # rank of its magnitude part to the P-theta deficit and reported 2.
  @testset "observable islands: a bus seen through one Q flow" begin (function ()
    net = fx.net14z
    meas = filter(m -> m.busIdx != 12 && !(_obs_is_injection(m) && m.busIdx in (6, 13)) && m.branchIdx != 19 && !(m.branchIdx == 12 && !(m.typ == Sparlectra.QflowMeas && m.direction == :from)), _obs_full_set(net))
    obs = _obs_cfg(() -> evaluate_global_observability(net, meas))
    isl = obs.observable_islands
    @test (obs.rank_deficit, obs.angle_rank_deficit, obs.magnitude_rank_deficit, obs.n_states - obs.numerical_rank) == (1, 1, 0, 1)
    @test isl.island_buses == [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14], [12]]
    @test isl.unobservable_branches == [12, 19]
    @test isl.region_buses == [12] && isempty(isl.dark_magnitude_buses)
    @test obs.unobservable_state_columns == [11]
  end)() end

  # Boundary rows on two small nets.
  # Net 1: S - A (lossy), A - X (lossless), A - Y (lossy); Vm at every bus,
  # both flows of S - A and the Q injection at A. The P-theta model has the
  # P flow of S - A as its only row: the angles of X and Y are free, two
  # zero columns, so the observable islands are {S, A}, {X}, {Y} in the
  # first pass, with A - X and A - Y unobservable. The Q injection is no
  # P-theta row and is not dropped.
  # Net 2: S - B1 (lossy), B1 = B1a (closed link), B1 - X (lossless),
  # B1a - Y (lossy); Vm at every bus and the flows of S - B1. A Q flow
  # through the link or the Q injections at B1 and B1a (aggregated into one
  # row of the fused bus) carry no P-theta row: the same three islands,
  # nothing dropped. A P flow through the link enters the P-theta model by
  # its row model: the minimum-norm split of the fused bus, half the
  # balance of each side, so the flows on S - B1, B1 - X and B1a - Y. With
  # S - B1 measured it is one equation for the two free angles of X and Y:
  # deficit 1, X and Y move in a fixed ratio, so the first pass finds three
  # groups with B1 - X and B1a - Y unobservable. The link row contains both
  # flows and is dropped; the second pass gives the same three islands.
  @testset "observable islands: boundary rows and link rows" begin (function ()
    function small_net(name, lines, loads; link = nothing, open = Int[])
      net = Net(name = name, baseMVA = 100.0)
      for b in unique(vcat(["S"], [l[1] for l in lines], [l[2] for l in lines], link === nothing ? String[] : collect(link)))
        addBus!(net = net, busName = b, vn_kV = 110.0)
      end
      for (k, (f, t, r)) in enumerate(lines)
        addPIModelACLine!(net = net, fromBus = f, toBus = t, r_pu = r, x_pu = 0.2, b_pu = 0.0, status = k in open ? 0 : 1)
      end
      addProsumer!(net = net, busName = "S", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "S")
      for (b, p, q) in loads
        addProsumer!(net = net, busName = b, type = "ENERGYCONSUMER", p = p, q = q)
      end
      link === nothing || addLink!(net = net, fromBus = link[1], toBus = link[2], status = 1)
      ok, msg = validate!(net = net)
      ok || error(msg)
      runpf!(net, 40, 1e-10, 0; method = :rectangular)
      return net
    end
    net1 = small_net("noisecol", [("S", "A", 0.02), ("A", "X", 0.0), ("A", "Y", 0.02)], [("X", 10.0, 3.0), ("Y", 20.0, 6.0)])
    meas1 = filter(m -> m.typ == Sparlectra.VmMeas || m.branchIdx == 1 || (m.typ == Sparlectra.QinjMeas && m.busIdx == 2), _obs_full_set(net1))
    one = _obs_cfg(() -> identify_observable_islands(net1, meas1; max_passes = 1))
    @test one.island_buses == [[1, 2], [3], [4]]
    for dmn in (Sparlectra._SE_DENSE_LINALG_MAX_N, 0)
      fin = _obs_cfg(() -> identify_observable_islands(net1, meas1; dense_max_n = dmn))
      @test fin.island_buses == [[1, 2], [3], [4]]
      @test fin.unobservable_branch_names == ["A - X", "A - Y"]
      @test isempty(fin.removed_boundary_injections) && fin.passes == 1
      @test (fin.angle_rank_deficit, fin.magnitude_rank_deficit) == (2, 0)
    end

    net2 = small_net("linkbnd", [("S", "B1", 0.02), ("B1", "X", 0.0), ("B1a", "Y", 0.02)], [("B1a", 5.0, 2.0), ("X", 10.0, 3.0), ("Y", 20.0, 6.0)]; link = ("B1", "B1a"))
    base2 = filter(m -> m.typ == Sparlectra.VmMeas || (m.branchIdx == 1 && m.typ in (Sparlectra.PflowMeas, Sparlectra.QflowMeas)), _obs_full_set(net2))
    with_link = copy(base2)
    Sparlectra.addMeasurement!(with_link; typ = Sparlectra.QflowMeas, value = 0.02, sigma = 0.01, direction = :from, id = "link_Q", linkIdx = 1)
    with_inj = copy(base2)
    Sparlectra.addMeasurement!(with_inj; typ = Sparlectra.QinjMeas, value = -0.02, sigma = 0.01, busIdx = 3, id = "Qinj_B1a")
    Sparlectra.addMeasurement!(with_inj; typ = Sparlectra.QinjMeas, value = 0.0, sigma = 0.01, busIdx = 2, id = "Qinj_B1")
    il = _obs_cfg(() -> identify_observable_islands(net2, with_link))
    ii = _obs_cfg(() -> identify_observable_islands(net2, with_inj))
    i0 = _obs_cfg(() -> identify_observable_islands(net2, base2))
    @test il.island_buses == ii.island_buses == i0.island_buses == [[1, 2, 3], [4], [5]]
    @test il.unobservable_branch_names == ii.unobservable_branch_names == ["B1 - X", "B1a - Y"]
    @test isempty(il.removed_boundary_injections) && isempty(ii.removed_boundary_injections)
    with_p = copy(base2)
    Sparlectra.addMeasurement!(with_p; typ = Sparlectra.PflowMeas, value = 0.1, sigma = 0.01, direction = :from, id = "link_P", linkIdx = 1)
    ip = _obs_cfg(() -> identify_observable_islands(net2, with_p))
    @test ip.island_buses == [[1, 2, 3], [4], [5]] && ip.unobservable_branch_names == ["B1 - X", "B1a - Y"]
    @test ip.removed_boundary_injections == ["link_P"] && ip.passes == 2 && ip.region_buses == [4, 5]
    @test ip.angle_rank_deficit == 1 && i0.angle_rank_deficit == 2
    # the same pattern with the open line S - Y in front of the two branches
    # to X and Y: the working net keeps branches 1, 3 and 4, so its
    # numbering differs from the caller's. The placement names the
    # unobservable branches 3 and 4 like identify_observable_islands (1.1.0
    # mapped them twice and failed with a BoundsError), and two P
    # injections close the two angle gaps.
    net3 = small_net("linkopen", [("S", "B1", 0.02), ("S", "Y", 0.02), ("B1", "X", 0.02), ("B1a", "Y", 0.02)], [("B1a", 5.0, 2.0), ("X", 10.0, 3.0), ("Y", 20.0, 6.0)]; link = ("B1", "B1a"), open = [2])
    base3 = filter(m -> m.typ == Sparlectra.VmMeas || (m.branchIdx == 1 && m.typ in (Sparlectra.PflowMeas, Sparlectra.QflowMeas)), _obs_full_set(net3))
    @test _obs_cfg(() -> identify_observable_islands(net3, base3)).unobservable_branches == [3, 4]
    pl3 = _obs_cfg(() -> suggest_measurement_placement(net3, base3; kind = :injection))
    @test (pl3.unobservable_branches_before, pl3.unobservable_branches_after, length(pl3.steps), pl3.reached) == ([3, 4], Int[], 2, true)
  end)() end

  # Test 5: an observable set computes nothing new. The old fields are
  # compared against the same call with the islands switched off, which
  # shows that the island step has no side effect. Rank and dark columns
  # are also recomputed outside the observability code: the labeled
  # Jacobian of `measurement_jacobian`, rows and columns normalized here,
  # a dense SVD, the rank tolerance rank_tol_factor * jac_eps * sigma_max
  # and the dark rule of 1.0.0 (null-vector entries above 1e-8). The whole
  # result against a run of the 1.0.0 sources was compared outside the
  # suite (block A report).
  @testset "observable islands: an observable set stays as it was" begin (function ()
    for net in (fx.net14z, fx.net60)
      full = _obs_full_set(net)
      on = _obs_cfg(() -> evaluate_global_observability(net, full))
      off = _obs_cfg(() -> evaluate_global_observability(net, full; observable_islands = false))
      @test on.numerical_observable
      @test isempty(on.observable_islands) && isempty(off.observable_islands)
      @test all(k -> k === :observable_islands || isequal(getfield(on, k), getfield(off, k)), keys(off))
      @test (on.numerical_rank, on.unobservable_state_columns) == _obs_svd_reference(net, full) == (on.n_states, Int[])
      @test (on.rank_deficit, on.angle_rank_deficit, on.magnitude_rank_deficit) == (0, 0, 0)
    end
    # on a not-observable set the islands change no old field either
    net = fx.net14z
    meas = filter(m -> m.branchIdx != 14 && !(_obs_is_injection(m) && m.busIdx in (7, 8)), _obs_full_set(net))
    on = _obs_cfg(() -> evaluate_global_observability(net, meas))
    off = _obs_cfg(() -> evaluate_global_observability(net, meas; observable_islands = false))
    @test isempty(off.observable_islands) && !isempty(on.observable_islands)
    @test all(k -> k === :observable_islands || isequal(getfield(on, k), getfield(off, k)), keys(off))
    @test (on.numerical_rank, on.unobservable_state_columns) == _obs_svd_reference(net, meas)
  end)() end

  # Test 6: the sparse path (dense_max_n = 0) against the SVD path. The
  # gain-matrix iteration gives an orthonormal null basis (a matrix with
  # known dependent columns), it counts near-singular directions the
  # column test of the sparse QR keeps (Kahan's matrix), and both paths give
  # the same rank, deficits (AC, P-theta, Q-V), observable islands with
  # their labels, region, unobservable branches and removed rows, the same
  # dark buses and the same weak report, on the case14 and sp_case60 gaps
  # and on fixed removal patterns of sp_case118 and sp_case300. Without a
  # removed row the P-theta deficit is the island count minus one.
  @testset "observable islands: SVD and sparse QR agree" begin (function ()
    rng = Xoshiro(521)
    A = SparseArrays.sprandn(rng, 40, 30, 0.15)
    A[:, 5] = A[:, 3] + A[:, 7]          # columns 3, 5, 7 dependent
    A[:, 12] .= 0.0                        # an unmeasured state
    A[:, 20] = 2.0 .* A[:, 1]              # columns 1, 20 dependent
    SparseArrays.dropzeros!(A)
    d, N = Sparlectra._gain_null_space(A, 1e-10, 1)
    @test d == 3
    @test LinearAlgebra.norm(N' * N - LinearAlgebra.I) <= 1e-12
    @test LinearAlgebra.norm(A * N) <= 1e-12
    @test Sparlectra._null_basis_dark_columns(N, 1e-8) == [1, 3, 5, 7, 12, 20]
    @test evaluate_observability_matrix(A; dense_max_n = 0).unobservable_state_columns == evaluate_observability_matrix(A).unobservable_state_columns == [1, 3, 5, 7, 12, 20]
    # Kahan's matrix K (n = 40, theta = 1.2): upper triangular with the
    # diagonal sin(theta)^(i-1), at least 0.064, so the column test of the
    # sparse QR keeps every column at tol 1e-6, while the smallest singular
    # value is 6.9e-7. With an unmeasured state next to it the set has a
    # deficit, and both paths count rank 39 of 41
    s, c = sin(1.2), cos(1.2)
    K = [i <= j ? s^(i - 1) * (i == j ? 1.0 : -c) : 0.0 for i = 1:40, j = 1:40]
    Ak = SparseArrays.sparse(hcat(vcat(K, zeros(5, 40)), zeros(45)))
    @test LinearAlgebra.rank(LinearAlgebra.qr(Ak; tol = 1e-6)) == 40
    dk, Nk = Sparlectra._gain_null_space(Ak, 1e-6, 1)
    @test 41 - dk == 39 == count(>(1e-6), LinearAlgebra.svdvals(Matrix(Ak)))
    @test all(j -> LinearAlgebra.norm(Ak * Nk[:, j]) <= 1e-6, 1:dk)
    @test evaluate_observability_matrix(Ak; tol = 1e-6, dense_max_n = 0).numerical_rank == evaluate_observability_matrix(Ak; tol = 1e-6).numerical_rank == 39

    z = fx.net14z
    full14 = _obs_full_set(z)
    endbus(m) = Int(m.direction == :to ? z.branchVec[m.branchIdx].toBus : z.branchVec[m.branchIdx].fromBus)
    sets = [
      ("case14 7-8", z, filter(m -> m.branchIdx != 14 && !(_obs_is_injection(m) && m.busIdx in (7, 8)), full14)),
      ("case14 gap", z, filter(m -> !(m.busIdx !== nothing && m.busIdx in 9:14) && !(m.branchIdx !== nothing && endbus(m) in 9:14), full14)),
      ("case14 boundary", z, filter(m -> !(m.branchIdx !== nothing && m.branchIdx in (8, 9, 10, 11, 15, 17)) && !(_obs_is_injection(m) && m.busIdx in (4, 5, 6, 7, 11, 14)), full14)),
      # buses 50 to 54 and every branch at them unmeasured: five islands,
      # one of them the floating group {52, 58, 59, 60}
      ("sp_case60 gap", fx.net60, filter(m -> !(m.busIdx !== nothing && m.busIdx in 50:54) && !(m.branchIdx !== nothing && (Int(fx.net60.branchVec[m.branchIdx].fromBus) in 50:54 || Int(fx.net60.branchVec[m.branchIdx].toBus) in 50:54)), _obs_full_set(fx.net60))),
    ]
    # fixed removal patterns on the shipped sp_case118 and sp_case300, with
    # the SVD as the reference: thinned sets (each row kept with a drawn
    # share of 30 to 60 percent, fixed seeds) and ten-bus regions whose bus
    # rows and branch rows are all removed. Left out: sp_case118 seed 5,
    # where the SVD has one singular value (1.38e-5) just below the
    # tolerance (1.62e-5) and the sparse QR verdict calls the set observable;
    # an observable verdict of the sparse QR is not checked
    mp_net(name) = begin
      net = Sparlectra.createNetFromMatPowerFile(filename = joinpath(pkgdir(Sparlectra), "data", "mpower", name), flatstart = false, bus_shunt_model = :admittance, matpower_shift_sign = 1.0, matpower_shift_unit = :deg, matpower_ratio = :normal, tap_changer_model = :ideal)
      redirect_stdout(() -> runpf!(net, 40, 1e-8, 0), devnull)
      net
    end
    touches(net, m, r) = (m.busIdx !== nothing && m.busIdx in r) || (m.branchIdx !== nothing && (Int(net.branchVec[m.branchIdx].fromBus) in r || Int(net.branchVec[m.branchIdx].toBus) in r))
    for (name, seeds, regions) in (("sp_case118.m", (1, 2, 3, 7), (1:10, 59:68, 98:107)), ("sp_case300.m", (1, 2, 3, 7, 8), (1:10, 40:49, 280:289)))
      net = mp_net(name)
      full = _obs_full_set(net)
      for seed in seeds
        rs = Xoshiro(seed)
        frac = 0.3 + 0.3 * rand(rs)
        push!(sets, ("$(name) seed $(seed)", net, [m for m in full if rand(rs) < frac]))
      end
      for r in regions
        push!(sets, ("$(name) region $(r)", net, filter(m -> !touches(net, m, r), full)))
      end
    end
    # the shipped MATPOWER cases warn by design on these sets: the
    # diag(Omega) criticality of a singular gain matrix falls back to the
    # dense path
    expected = (r"Takahashi selected inverse unavailable",)
    for (label, net, meas) in sets
      a = run_with_expected_warnings(() -> _obs_cfg(() -> evaluate_global_observability(net, meas)), expected)
      b = run_with_expected_warnings(() -> _obs_cfg(() -> evaluate_global_observability(net, meas; dense_max_n = 0)), expected)
      ia = a.observable_islands
      ib = b.observable_islands
      @test (label, a.numerical_rank, a.angle_rank_deficit, a.magnitude_rank_deficit) == (label, b.numerical_rank, b.angle_rank_deficit, b.magnitude_rank_deficit)
      @test (label, ia.n_islands >= 1, a.quality) == (label, true, :not_observable)
      @test (label, ia.island_of_bus, ia.island_labels, ia.unobservable_branches, ia.removed_boundary_injections, ia.region_buses) == (label, ib.island_of_bus, ib.island_labels, ib.unobservable_branches, ib.removed_boundary_injections, ib.region_buses)
      @test (label, ia.rank_deficit, ia.angle_rank_deficit) == (label, ib.rank_deficit, ib.angle_rank_deficit)
      isempty(ia.removed_boundary_injections) && @test (label, ia.angle_rank_deficit) == (label, ia.n_islands - 1)
      @test (label, a.weak_observability) == (label, b.weak_observability)
      @test (label, ia.dark_magnitude_buses) == (label, ib.dark_magnitude_buses)
      @test (label, _obs_islands_consistent(net, ia), _obs_islands_consistent(net, ib)) == (label, true, true)
      # the verdict field and the report read the same null bases with the
      # same tolerance (without couplers state columns and buses map one to
      # one; sp_case60 has couplers), and both paths give the same dark
      # columns (concept D1: one rule for the SVD and the iteration)
      isempty(net.linkVec) && @test (label, a.unobservable_state_columns, b.unobservable_state_columns) == (label, _obs_report_columns(net, ia), _obs_report_columns(net, ib))
      @test (label, b.unobservable_state_columns) == (label, a.unobservable_state_columns)
    end
  end)() end

  # Test 11 (concept N5): network islands and observable islands. In
  # two_islands_prio (network islands {1, 2, 3} and the ring 4-5-6-7,
  # branch 3-4 out of service) the branches 4-5 and 6-7 are opened as well:
  # three network islands {1, 2, 3}, {4, 7} and {5, 6}. Bus 2 floats in the
  # first (its rows, the flows on 1-2 and 2-3 and the injections at 1 and 3
  # removed), the second is measured in full, the third carries no
  # measurement. By hand: in network island 1 neither the angle nor the
  # magnitude of bus 2 appears in a row, P-theta deficit 1 and Q-V deficit
  # 1, observable islands 1.1 = {1, 3} and 1.2 = {2}; network island 2 is
  # one observable island "2" without a deficit; network island 3 is
  # reported as not measured and adds nothing. The full deficit is 2, from
  # network island 1 alone, on the dense and on the sparse path. The
  # measurements come from the power flow before the branches were opened;
  # observability reads only their structure.
  @testset "observable islands: network islands and observable islands" begin (function ()
    net = deepcopy(fx.net_two)
    full = _obs_full_set(net)
    setBranchStatus!(net.branchVec[5], false)
    setBranchStatus!(net.branchVec[7], false)
    meas = filter(m -> !(m.busIdx in (2, 5, 6)) && !(m.branchIdx !== nothing && m.branchIdx in (1, 2, 5, 6, 7)) && !(_obs_is_injection(m) && m.busIdx in (1, 3)), full)
    reports = [_obs_cfg(() -> evaluate_global_observability(net, meas; dense_max_n = dmn)) for dmn in (Sparlectra._SE_DENSE_LINALG_MAX_N, 0)]
    for obs in reports
      isl = obs.observable_islands
      @test :island_partition in obs.notes
      @test (obs.n_measured_islands, obs.n_unmeasured_islands) == (2, 1)
      @test (obs.rank_deficit, obs.angle_rank_deficit, obs.magnitude_rank_deficit) == (2, 1, 1)
      @test isl.island_buses == [[1, 3], [2], [4, 7]] && isl.island_labels == ["1.1", "1.2", "2"]
      @test isl.network_island_of_bus == [1, 1, 1, 2, 3, 3, 2]
      @test isl.island_of_bus == [1, 2, 1, 3, 0, 0, 3]
      @test isl.unmeasured_network_islands == [3]
      @test [(r.network_island, r.measured, r.n_islands, r.angle_rank_deficit, r.magnitude_rank_deficit) for r in isl.network_islands] == [(1, true, 2, 1, 1), (2, true, 1, 0, 0), (3, false, 0, 0, 0)]
      @test isl.unobservable_branches == [1, 2] && isl.region_buses == [2] && isl.dark_magnitude_buses == [2]
      @test _obs_islands_consistent(net, isl)
    end
    @test reports[1].observable_islands.network_islands == reports[2].observable_islands.network_islands
    io = IOBuffer()
    print_observable_islands(io, reports[1])
    out = String(take!(io))
    @test occursin("network island 3 (2 buses): not measured", out) && occursin("observable island 1.1", out)
    # an observable set: no report in the verdict, one observable island
    # per measured network island
    @test isempty(_obs_cfg(() -> evaluate_global_observability(fx.net_two, full)).observable_islands)
    two = _obs_cfg(() -> identify_observable_islands(fx.net_two, full))
    @test two.island_of_bus == [1, 1, 1, 2, 2, 2, 2] && two.island_labels == ["1", "2"]
    # estimated one by one: bus 2 has a dark magnitude and leaves (concept
    # M1), so network island 1 keeps one unit, numbered 1.1 because it lost
    # a bus; network island 2 keeps its promoted reference (bus 4), so no
    # angle is relative; network island 3 is not measured.
    res = _obs_mode(() -> runse!(deepcopy(net), meas), :observable_islands)
    @test [(r.label, r.network_island, r.estimated, r.reason) for r in res.islands] == [("1.1", 1, true, :estimated), ("2", 2, true, :estimated), ("3", 3, false, :no_measurements)]
    @test res.observableIslands.removed_buses == [(bus = 2, reason = :dark_magnitude)]
    @test res.observableIslands.label_of_bus == ["1.1", "", "1.1", "2", "3", "3", "2"]
    @test res.observableIslands.estimated == [true, false, true, true, false, false, true]
    @test !any(res.observableIslands.relative_angle) && isequal(isnan.(res.voltages), .!res.observableIslands.estimated)
    # The SE service on a case file with these network islands. abort:
    # rejected with the islands in the metadata (1.1.0 crashed on the
    # criticality fields the merged verdict lacked); observable_islands:
    # estimated per observable island as runse! above (2 units). With P rows
    # only no magnitude has a row (Q-V deficit 7, every bus dark), no
    # observable island can be estimated, and the service rejects the set
    # like abort instead of failing with a backtrace. On the observable set
    # of two_islands_prio with weak_check_after_estimation = false the run
    # computes no weak class after the estimation.
    service(m, tag, mode; weak_check::Bool = true) = mktempdir() do d
      case_net = deepcopy(m === full ? fx.net_two : net)
      case_net.measurements = m
      case = exportSCF(case_net; file = joinpath(d, "case.scf.json"), intended_calculations = String["power_flow", "state_estimation"])
      cfg = joinpath(d, "se.yaml")
      write(cfg, "config_version: 1\nstate_estimation:\n  observability:\n    on_unobservable: $(mode)\n    weak_check_after_estimation: $(weak_check)\n")
      r = _se_run_quiet(() -> redirect_stdout(() -> SparlectraApp._run_state_estimation_service(case, cfg, joinpath(d, "run"), tag, joinpath(d, "none.csv")), devnull))
      SparlectraApp.to_dict(r)
    end
    ab = service(meas, "isl_abort", "abort")
    @test (ab["status"], ab["reason"], ab["metadata"]["se_islands_total"], ab["metadata"]["se_criticality_skipped"]) == ("failed", "se_not_observable", 3, false)
    oi = service(meas, "isl_oi", "observable_islands")
    @test (oi["status"], oi["metadata"]["se_observable_islands_estimated"]) == ("succeeded", 2)
    ponly = service(filter(m -> m.typ in (Sparlectra.PinjMeas, Sparlectra.PflowMeas), meas), "isl_ponly", "observable_islands")
    @test (ponly["status"], ponly["reason"]) == ("failed", "se_not_observable")
    @test occursin("none of its observable islands can be estimated", ponly["message"])
    nowk = service(full, "isl_noweak", "abort"; weak_check = false)
    @test (nowk["status"], nowk["metadata"]["se_weak_basis"]) == ("succeeded", "not_determined")
  end)() end

  @testset "observable islands: a link cluster in one of two network islands" begin (function ()
    # A1-A2-A3 and B1-B2 with B2a on B2 through a closed link only (no
    # branch, the last bus). The cluster injections of B2/B2a become LINKAGG
    # rows with caller index 0. 1.1.0 failed twice here: the measurement
    # generator read the Y-bus row of B2a eagerly (BoundsError), and the
    # merge of the network islands indexed the kept rows with that 0.
    net = Net(name = "obs_islands_link", baseMVA = 100.0)
    for b in ("A1", "A2", "A3", "B1", "B2", "B2a")
      addBus!(net = net, busName = b, vn_kV = 110.0)
    end
    addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = net, fromBus = "A2", toBus = "A3", r_pu = 0.012, x_pu = 0.09, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = net, fromBus = "A1", toBus = "A3", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
    addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", referencePri = "A1", vm_pu = 1.02, va_deg = 0.0)
    addProsumer!(net = net, busName = "A2", type = "ENERGYCONSUMER", p = 30.0, q = 10.0)
    addProsumer!(net = net, busName = "A3", type = "ENERGYCONSUMER", p = 20.0, q = 6.0)
    addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.0, va_deg = 0.0)
    addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 12.0, q = 4.0)
    addProsumer!(net = net, busName = "B2a", type = "ENERGYCONSUMER", p = 3.0, q = 1.0)
    addLink!(net = net, fromBus = "B2", toBus = "B2a", status = 1)
    @test first(validate!(net = net))
    _, erg = runpf!(net, 40, 1e-10, 0; method = :rectangular, islands_enabled = true)
    @test erg == 0
    meas = generateMeasurementsFromPF(net; noise = false)
    # the cluster rows read generation minus load of their own bus
    pB2a = only(m for m in meas if m.typ == Sparlectra.PinjMeas && m.busIdx == 6)
    @test isapprox(pB2a.value, -3.0; atol = 1e-9)
    obs = _obs_cfg(() -> evaluate_global_observability(net, meas))
    @test obs.n_measured_islands == 2 && obs.quality !== :not_observable
    @test length(obs.active_measurement_indices) == length(obs.criticality_wii)
    @test count(==(0), obs.active_measurement_indices) == 2   # the LINKAGG Pinj and Qinj rows
    @test all(i -> 0 <= i <= length(meas), obs.active_measurement_indices)
  end)() end

  @testset "observable islands: a later run starts from a net with NaN voltages" begin (function ()
    # sp_case14 with the magnitude of the leaf bus 14 dark (its Vm, Qinj and
    # Q flows and the Qinj of its neighbour 7 removed): Q-V deficit 1. The
    # SE service with released taps runs the staged start: a pre-solve per
    # observable island writes NaN for bus 14 into the net, the released run
    # starts from the stored voltages. 1.1.0 took the NaN into the Jacobian
    # and failed with "invalid argument #4 to LAPACK call"; a NaN stored
    # voltage now starts flat.
    net = importSCF(joinpath(pkgdir(Sparlectra), "data", "scf", "sp_case14.scf.json"))
    runpf!(net, 40, 1e-8, 0; method = :rectangular)
    b = 14
    nb = Set{Int}(br.fromBus == b ? br.toBus : br.fromBus for br in net.branchVec if b in (br.fromBus, br.toBus))
    @test nb == Set([7])
    atB(m) = m.busIdx == b || (m.branchIdx !== nothing && b in (net.branchVec[m.branchIdx].fromBus, net.branchVec[m.branchIdx].toBus))
    qv = (Sparlectra.VmMeas, Sparlectra.QinjMeas, Sparlectra.QflowMeas)
    meas = filter(m -> !(atB(m) && m.typ in qv) && !(m.typ == Sparlectra.QinjMeas && m.busIdx in nb), _obs_full_set(net))
    obs = _obs_cfg(() -> evaluate_global_observability(net, meas))
    @test (obs.quality, obs.angle_rank_deficit, obs.magnitude_rank_deficit) == (:not_observable, 0, 1)
    # the library: a second runse! from the stored state of the first
    netL = deepcopy(net)
    oc = Sparlectra.ObservabilityConfig(on_unobservable = :observable_islands)
    r1 = _se_run_quiet(() -> with_state_estimation_config(() -> runse!(netL, meas); flatstart = true, update_net = true, observability = oc))
    @test r1.converged && isnan(netL.nodeVec[b]._vm_pu)
    r2 = _se_run_quiet(() -> with_state_estimation_config(() -> runse!(netL, meas); flatstart = false, update_net = true, observability = oc))
    @test r2.converged && isnan(netL.nodeVec[b]._vm_pu) && !isnan(netL.nodeVec[1]._vm_pu)
    # the service with the staged tap start
    mktempdir() do d
      netS = deepcopy(net)
      netS.measurements = meas
      case = exportSCF(netS; file = joinpath(d, "case.scf.json"), intended_calculations = String["power_flow", "state_estimation"])
      cfg = joinpath(d, "se.yaml")
      write(cfg, "config_version: 1\nstate_estimation:\n  observability:\n    on_unobservable: observable_islands\n")
      r = SparlectraApp.to_dict(_se_run_quiet(() -> redirect_stdout(() -> SparlectraApp._run_state_estimation_service(case, cfg, joinpath(d, "run"), "nan_start", joinpath(d, "none.csv")), devnull)))
      @test (r["status"], r["metadata"]["se_tap_staged_start"], r["metadata"]["se_observable_islands_estimated"]) == ("succeeded", true, 1)
    end
    # a non-finite Jacobian entry is named, not a LAPACK argument error
    t = @test_throws ArgumentError Sparlectra._sigma_max(Sparlectra.sparse([1.0 NaN; 0.0 2.0]))
    @test occursin("non-finite entries (first at row 1, column 2)", t.value.msg)
    # also with an explicit tolerance, where _sigma_max is not called and
    # the NaN reached CHOLMOD and the SVD before
    tp = @test_throws ArgumentError Sparlectra.numeric_rank_pivots(Sparlectra.sparse([1.0 NaN; 0.0 2.0]); tol = 1e-8)
    @test occursin("non-finite entries (first at row 1, column 2)", tp.value.msg)
    # a LINKAGG aggregate (index 0 in active_measurement_indices) has no row
    # to name: the nearly critical list leaves it out instead of "#0"
    obs0 = (criticality_skipped = false, criticality_method = :omega, numerical_critical_measurement_indices = Int[], structural_critical_measurement_indices = Int[], structural_criticality_skipped = false, criticality_wii = [0.1, 0.2], active_measurement_indices = [0, 1])
    netM = deepcopy(net)
    netM.measurements = [Sparlectra.Measurement(typ = Sparlectra.VmMeas, value = 1.0, sigma = 0.01, busIdx = 1, id = "u1")]
    near0 = only(filter(l -> occursin("nearly critical", l), SparlectraApp._se_criticality_lines(netM, obs0)))
    @test occursin("u1 (0.2)", near0) && !occursin("#0", near0) && occursin(": 1: ", near0)
  end)() end

  return true
end

## The weak gaps, the state file and the restoration with scheduled values
## (concept 1.1.0, M1, M2, M4 and block B) on the shipped IEEE 14-bus
## topology with its slack on bus 8 (`net14z`) and on two_islands_prio.
function test_observability_restore(fx)::Bool
  # M1: a bus with a dark magnitude cuts an observable island in two. Bus
  # 7 is the hub of 4-7, 7-8 and 7-9. Removed: the Vm row at 7, the Q
  # injections at 4, 7, 8 and 9 and the Q flows on the three branches at 7
  # (8, 14, 15), i.e. every Q-V row that holds the magnitude of 7. By
  # hand: every P row stays, so the P-theta model has full rank and the set
  # is one observable island; the magnitude of 7 appears in no row (Q-V
  # deficit 1, bus 7 dark). Without bus 7 the branches 8, 14 and 15 are
  # unobservable and the island falls into {1..6, 9..14} and {8}: two
  # observable islands of their own, 1.1 and 1.2 (no third level), bus 7
  # removed with reason :dark_magnitude. The P injections at 4, 8 and 9 hold
  # a flow of those branches and are not used; island 1.2 keeps the Vm row
  # of bus 8 and is estimated with the slack, island 1.1 through its flows
  # with a relative reference.
  @testset "observable islands: a dark magnitude cuts an island in two" begin (function ()
    net = fx.net14z
    full = _obs_full_set(net)
    cut = filter(m -> !(m.typ == Sparlectra.VmMeas && m.busIdx == 7) && !(m.typ == Sparlectra.QinjMeas && m.busIdx in (4, 7, 8, 9)) && !(m.typ == Sparlectra.QflowMeas && m.branchIdx in (8, 14, 15)), full)
    isl = _obs_cfg(() -> identify_observable_islands(net, cut))
    @test (isl.n_islands, isl.dark_magnitude_buses, isl.angle_rank_deficit, isl.magnitude_rank_deficit) == (1, [7], 0, 1)
    res = _obs_mode(() -> runse!(deepcopy(net), cut), :observable_islands)
    oi = res.observableIslands
    @test res.converged
    @test [(r.label, r.estimated, r.reason, r.buses) for r in res.islands] == [("1.1", true, :estimated, [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14]), ("1.2", true, :estimated, [8])]
    @test oi.removed_buses == [(bus = 7, reason = :dark_magnitude)]
    @test oi.unobservable_branches == [8, 14, 15]
    @test oi.label_of_bus == [fill("1.1", 6); ""; "1.2"; fill("1.1", 6)]
    @test res.islands[1].relative_angles && !res.islands[2].relative_angles
    @test issubset(["Pinj_bus_4", "Pinj_bus_8", "Pinj_bus_9"], oi.dropped_measurements)
    @test isnan(res.voltages[7]) && all(b -> isfinite(res.voltages[b]), setdiff(1:14, [7]))
    io = IOBuffer()
    print_observable_islands(io, res)
    @test occursin("buses not estimated (dark magnitude): 1 bus (7)", String(take!(io)))
  end)() end

  # M2: a weak direction at the estimate is a gap. The flows on 7-8 and the
  # injections at 7 are removed and the P injection at 8 carries 50 MW
  # instead of 1 MW. Bus 8 hangs on 7-8 alone, so by hand the active flow
  # on 7-8 is the injection at 8 (losses aside) and follows from that row
  # only: its standard deviation is about 50 MW, above the 20 MW bound, the
  # other flows are measured with 1 MW. The set is observable, so the first
  # run estimates the whole net; at its estimate branch 14 is the weak one
  # and its ends 7 and 8 are the affected buses (concept G3). They leave
  # with reason :weak, the rest {1..6, 9..14} is estimated again as one
  # observable island (1.1, it lost buses; relative angles, the slack 8 is
  # gone), and at that estimate nothing is weak. Under abort runse! only
  # reports (weak_check_after_estimation, default on): branch 14 above the
  # bound, its ends 7 and 8 affected, every voltage estimated.
  @testset "weak directions at the estimate are a gap" begin (function ()
    net = fx.net14z
    full = _obs_full_set(net)
    base = filter(m -> m.branchIdx != 14 && !(_obs_is_injection(m) && m.busIdx == 7), full)
    weakSet = [m.typ == Sparlectra.PinjMeas && m.busIdx == 8 ? Measurement(typ = m.typ, value = m.value, sigma = 50.0, busIdx = 8, id = m.id) : m for m in base]
    @test _obs_cfg(() -> evaluate_global_observability(net, weakSet)).quality === :weakly_observable
    res = _obs_mode(() -> runse!(deepcopy(net), weakSet), :observable_islands)
    oi = res.observableIslands
    @test res.converged && oi !== nothing
    @test oi.removed_buses == [(bus = 7, reason = :weak), (bus = 8, reason = :weak)]
    @test [(r.label, r.estimated, r.buses, r.relative_angles) for r in res.islands] == [("1.1", true, [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14], true)]
    @test res.weakObservability.basis === :estimated_state && res.weakObservability.n_buses == 0
    # abort: one estimate of the whole net, the weak class reported only
    ab = _obs_mode(() -> runse!(deepcopy(net), weakSet), :abort)
    @test ab.observableIslands === nothing && all(isfinite, ab.voltages)
    wa = ab.weakObservability
    @test (wa.basis, wa.weakly_observable, wa.branches, wa.buses) == (:estimated_state, true, [14], [7, 8])
    @test 40.0 < only(wa.branch_flow_sigma_mw) < 55.0
    # A cut estimate that does not converge keeps the converged one before
    # the cut. A gross error of 1000 MW on the flow of branch 9 (the start
    # ladder off, max_iter 30): the whole set absorbs it (17 iterations,
    # still weak at 7 and 8), the set without 7 and 8 does not converge in
    # 30. The result is the first estimate with its weak report, the net
    # holds its voltages (those of the same run under abort), and the flow
    # of branch 14 is left empty.
    gross = [m.id == "Pflow_branch_9_from" ? Measurement(typ = m.typ, value = m.value + 1000.0, sigma = m.sigma, branchIdx = m.branchIdx, direction = m.direction, id = m.id) : m for m in weakSet]
    gcfg(f, mode) = with_state_estimation_config(f; flatstart = true, jac_eps = 1e-6, start_ladder = false, max_iter = 30, observability = Sparlectra.ObservabilityConfig(on_unobservable = mode))
    gnet = deepcopy(net)
    gr = gcfg(() -> runse!(gnet, gross), :observable_islands)
    ga = gcfg(() -> runse!(deepcopy(net), gross), :abort)
    @test (ga.converged, ga.iterations, ga.weakObservability.buses) == (true, 17, [7, 8])
    @test (gr.converged, gr.iterations, gr.observableIslands === nothing, gr.weakObservability.buses) == (true, 17, true, [7, 8])
    @test gr.voltages == ga.voltages && [nd._vm_pu for nd in gnet.nodeVec] == abs.(ga.voltages)
    @test isnan(gnet.branchVec[14].fBranchFlow.pFlow)
  end)() end

  # M4: the state file of an estimate per observable island. The case14
  # split of test 1 estimates both islands with relative angles: read back,
  # every bus is estimated, branch 7-8 between the two islands is not, and
  # the state is partial, so runpf_from_se! refuses it with the message of
  # the in-memory estimate. Without the Vm rows of island 1.1 its 13 buses
  # are not estimated and no closed branch is (each has an end among them).
  # A full estimate reads back as before, without a warning.
  @testset "state file of an estimate per observable island" begin (function ()
    net = fx.net14z
    full = _obs_full_set(net)
    meas = filter(m -> m.branchIdx != 14 && !(_obs_is_injection(m) && m.busIdx in (7, 8)), full)
    pferr(n) = try
      runpf_from_se!(n, 40, 1e-8, 0)
      ""
    catch e
      sprint(showerror, e)
    end
    dir = mktempdir()
    est = deepcopy(net)
    _obs_mode(() -> runse!(est, meas), :observable_islands)
    f = joinpath(dir, "split.csv")
    @test writeSEStateCSV(est; file = f) == (count = 14, not_estimated = 0)
    @test startswith(readlines(f)[1], "# sparlectra-se-state v2")
    back = deepcopy(net)
    rb = @test_logs (:warn, r"0 bus\(es\) and 1 branch\(es\) not estimated") readSEStateCSV!(back; file = f)
    @test (rb.not_estimated_buses, rb.not_estimated_branches, rb.partial) == (Int[], [14], true)
    @test isnan(back.branchVec[14].fBranchFlow.pFlow)
    msg = pferr(est)
    @test occursin("observable islands only", msg) && pferr(back) == msg
    est2 = deepcopy(net)
    _obs_mode(() -> runse!(est2, filter(m -> !(m.typ == Sparlectra.VmMeas && m.busIdx != 8), meas)), :observable_islands)
    f2 = joinpath(dir, "dark.csv")
    @test writeSEStateCSV(est2; file = f2).not_estimated == 13
    back2 = deepcopy(net)
    rb2 = @test_logs (:warn, r"13 bus\(es\) and 20 branch\(es\) not estimated") readSEStateCSV!(back2; file = f2)
    @test rb2.not_estimated_buses == setdiff(1:14, [8]) && rb2.not_estimated_branches == collect(1:20)
    @test all(b -> isnan(back2.nodeVec[b]._vm_pu), rb2.not_estimated_buses) && pferr(back2) == msg
    estF = deepcopy(net)
    _obs_cfg(() -> runse!(estF, full))
    f3 = joinpath(dir, "full.csv")
    writeSEStateCSV(estF; file = f3)
    back3 = deepcopy(net)
    rb3 = @test_logs min_level = Logging.Warn readSEStateCSV!(back3; file = f3)
    @test !rb3.partial && isempty(rb3.not_estimated_buses) && isempty(rb3.not_estimated_branches)
  end)() end

  # Test 8 of the concept: restoration with scheduled values. The case14
  # split of test 1 has a P-theta deficit of 1 (bus 8 floats against the
  # rest) and full Q-V rank. P candidates are the buses without a P
  # injection row, 7 and 8; bus 7 has no load or generator (a passive bus,
  # no schedule), bus 8 the condenser with P = 0. Under zero_injection = off
  # the passive bus is no candidate: by hand the injection at 8 is the flow
  # on 7-8 and ties the angle of 8 to bus 7, so it closes the deficit; its
  # sigma is the floor, the default P sigma of 1 MW (30 % of 0 is 0), and
  # the zero injection at 7 is reported as one that would close it (its
  # row holds the flow on 7-8 as well). At the stored (power-flow) state
  # nothing is weak, so the restoration succeeds (R1), and runse! estimates
  # the whole net with the row. Under zero_injection = on (concept U3) the
  # passive bus comes first: P = 0 at bus 7 with the zero-injection sigma
  # closes the deficit, flagged as a zero-injection row, and the schedule
  # of bus 8 is not needed.
  @testset "restoration: the case14 split closed from the schedule" begin (function ()
    net = fx.net14z
    full = _obs_full_set(net)
    meas = filter(m -> m.branchIdx != 14 && !(_obs_is_injection(m) && m.busIdx in (7, 8)), full)
    rr = _obs_zioff(() -> restore_observability(net, meas))
    @test [(r.id, r.bus, r.quantity, r.value, r.sigma, r.source) for r in rr.restored] == [("restore_Pinj_bus_8", 8, :P, 0.0, 1.0, :schedule)]
    @test isempty(rr.skipped) && (rr.rank_deficit_before, rr.rank_deficit_after) == (1, 0)
    @test rr.observable && rr.weak_judged && rr.success && rr.verdict === :observable
    @test [(r.bus, r.quantity) for r in rr.passive_not_inserted] == [(7, :P)]
    @test only(rr.rows).origin === :pseudo_restore && Sparlectra._elimination_protected(only(rr.rows))
    res = _obs_mode(() -> runse!(deepcopy(net), meas), :scheduled; zero_injection = :off)
    @test res.converged && res.observableIslands === nothing
    @test [r.id for r in res.restoredRows] == ["restore_Pinj_bus_8"]
    @test _obs_mode(() -> runse!(deepcopy(net), meas), :observable_islands).restoredRows === nothing
    rz = _obs_cfg(() -> restore_observability(net, meas))
    zs = Sparlectra.state_estimation_config().zero_injection_sigma_mva
    @test [(r.id, r.bus, r.quantity, r.value, r.sigma, r.source) for r in rz.restored] == [("restore_Pinj_bus_7", 7, :P, 0.0, zs, :passive)]
    @test only(rz.rows).zero_injection && isempty(rz.passive_not_inserted) && rz.verdict === :observable
    # the replacement row is protected, a gross error elsewhere goes: the
    # from-end P flow of branch 1 shifted by 50 MW (50 sigma)
    bad = [m.typ == Sparlectra.PflowMeas && m.branchIdx == 1 && m.direction == :from ? Measurement(typ = m.typ, value = m.value + 50.0, sigma = m.sigma, branchIdx = 1, direction = :from, id = m.id) : m for m in meas]
    badId = only(m.id for m in bad if m.typ == Sparlectra.PflowMeas && m.branchIdx == 1 && m.direction == :from)
    diag = _obs_cfg(() -> runse_diagnostics(net, vcat(bad, rr.rows)))
    ids = [t.id for t in diag.eliminations]
    @test badId in ids && !("restore_Pinj_bus_8" in ids)
    # a redundant candidate: without the P injection at bus 2 its schedule
    # is a candidate before bus 8, but every flow at bus 2 is measured, so
    # its row is the sum of measured flow rows and does not lower the deficit
    no2 = filter(m -> !(m.typ == Sparlectra.PinjMeas && m.busIdx == 2), meas)
    r2 = _obs_zioff(() -> restore_observability(net, no2))
    @test [(r.bus, r.quantity) for r in r2.restored] == [(8, :P)]
    @test r2.skipped == [(bus = 2, network_island = 1, quantity = :P, source = :schedule)]
    # the schedule cannot close the gap (no P set point at bus 8): nothing
    # is taken, and :scheduled runs the rest like :observable_islands
    nosched = deepcopy(net)
    for ps in nosched.prosumpsVec
      Sparlectra.getPosumerBusIndex(ps) == 8 && (ps.pVal = nothing)
    end
    r3 = _obs_zioff(() -> restore_observability(nosched, meas))
    @test isempty(r3.restored) && (r3.rank_deficit_before, r3.rank_deficit_after) == (1, 1) && !r3.observable && !r3.success && r3.verdict === :not_observable
    s3 = _obs_mode(() -> runse!(deepcopy(nosched), meas), :scheduled; zero_injection = :off)
    o3 = _obs_mode(() -> runse!(deepcopy(nosched), meas), :observable_islands; zero_injection = :off)
    @test isempty(s3.restoredRows) && [r.label for r in s3.islands] == ["1.1", "1.2"]
    @test isequal(s3.voltages, o3.voltages)
  end)() end

  # Magnitude gaps closed by U in the order of concept R4, one case per
  # rule. Rules 1 and 2: the split without the Vm row at 8 leaves the
  # magnitude of 8 in no Q-V row (its injections and the flows on 7-8 are
  # gone), Q-V deficit 1; bus 8 carries the condenser with set point 1.09
  # pu (rule 2), with an SV value stored for bus 8 rule 1 wins. Rule 3: the
  # set of the M1 test (magnitude of 7 dark) on a copy whose loads at 4, 8
  # and 9 carry no Q set point, so no Q candidate: bus 7 has no unit and no
  # stored state, so U is the median of the 13 measured Vm rows of its
  # level, sigma their standard deviation (at least 0.005 pu). Rule 4: no
  # Vm row at all and no voltage-controlled unit (every prosumer
  # unregulated on a copy): the Q rows see magnitude differences only, one
  # common shift is free (Q-V deficit 1), and U at the first bus closes it
  # at 1.0 pu with restore_vm_sigma_nominal_pu (0.05). All under
  # zero_injection = off: the passive bus 7 would otherwise come first.
  @testset "restoration: U by the order of the stored state, set point, level median, nominal" begin (function ()
    net = fx.net14z
    full = _obs_full_set(net)
    split = filter(m -> m.branchIdx != 14 && !(_obs_is_injection(m) && m.busIdx in (7, 8)), full)
    no8 = filter(m -> !(m.typ == Sparlectra.VmMeas && m.busIdx == 8), split)
    rr = _obs_zioff(() -> restore_observability(net, no8))
    @test [(r.bus, r.quantity, r.value, r.sigma, r.source) for r in rr.restored] == [(8, :P, 0.0, 1.0, :schedule), (8, :U, 1.09, 0.005, :unit_setpoint)]
    @test rr.network_islands == [(network_island = 1, angle_rank_deficit = 1, magnitude_rank_deficit = 1, angle_rank_deficit_after = 0, magnitude_rank_deficit_after = 0, n_taken = 2)]
    sv = deepcopy(net)
    sv._sv_vm_pu[8] = 1.085
    r1 = _obs_zioff(() -> restore_observability(sv, no8))
    @test [(r.quantity, r.value, r.source) for r in r1.restored] == [(:P, 0.0, :schedule), (:U, 1.085, :stored_state)]
    cut = filter(m -> !(m.typ == Sparlectra.VmMeas && m.busIdx == 7) && !(m.typ == Sparlectra.QinjMeas && m.busIdx in (4, 7, 8, 9)) && !(m.typ == Sparlectra.QflowMeas && m.branchIdx in (8, 14, 15)), full)
    noq = deepcopy(net)
    for ps in noq.prosumpsVec
      Sparlectra.getPosumerBusIndex(ps) in (4, 8, 9) && (ps.qVal = nothing)
    end
    vms = sort([m.value for m in cut if m.typ == Sparlectra.VmMeas])
    @test length(vms) == 13
    mu = sum(vms) / 13
    r3 = _obs_zioff(() -> restore_observability(noq, cut))
    @test [(r.bus, r.quantity, r.source) for r in r3.restored] == [(7, :U, :level_median)]
    @test only(r3.restored).value == vms[7] && only(r3.restored).sigma ≈ max(sqrt(sum(abs2, vms .- mu) / 13), 0.005)
    @test r3.rank_deficit_after == 0
    novm = filter(m -> m.typ != Sparlectra.VmMeas, full)
    flat = deepcopy(net)
    for ps in flat.prosumpsVec
      ps.isRegulated = false
    end
    r4 = _obs_zioff(() -> restore_observability(flat, novm))
    @test [(r.bus, r.quantity, r.value, r.sigma, r.source) for r in r4.restored] == [(1, :U, 1.0, 0.05, :nominal)]
    @test (r4.network_islands[1].magnitude_rank_deficit, r4.rank_deficit_after) == (1, 0)
  end)() end

  # N4, N5: candidates stay inside their network island. The set of the
  # network-island test (network island 1 with a P-theta and a Q-V deficit
  # at bus 2, network island 2 observable, network island 3 not measured):
  # every taken row lies in network island 1, network island 3 gets none
  # although its buses carry schedules, and network island 2 needs none.
  @testset "restoration: candidates never cross a network island" begin (function ()
    net = deepcopy(fx.net_two)
    full = _obs_full_set(net)
    setBranchStatus!(net.branchVec[5], false)
    setBranchStatus!(net.branchVec[7], false)
    meas = filter(m -> !(m.busIdx in (2, 5, 6)) && !(m.branchIdx !== nothing && m.branchIdx in (1, 2, 5, 6, 7)) && !(_obs_is_injection(m) && m.busIdx in (1, 3)), full)
    isl = _obs_cfg(() -> identify_observable_islands(net, meas))
    rr = _obs_cfg(() -> restore_observability(net, meas))
    @test [r.network_island for r in rr.network_islands] == [1, 2]
    @test (rr.network_islands[1].angle_rank_deficit, rr.network_islands[1].magnitude_rank_deficit, rr.network_islands[2].n_taken) == (1, 1, 0)
    @test !isempty(rr.restored) && all(r -> r.network_island == 1 && isl.network_island_of_bus[r.bus] == 1, rr.restored)
    @test all(r -> isl.network_island_of_bus[r.bus] == 1, rr.skipped)
  end)() end

  return true
end

## The weak class (concept K3, G3, U4) against references computed here
## from the labeled Jacobian of `measurement_jacobian` at the stored state.
## The covariance reference: G = H' W H with W = 1/sigma^2, the covariance
## inv(G) dense, or with constraint rows (`constraint = true`: the
## zero-injection rows, weight 0) the upper-left block of the dense inverse
## of the KKT matrix [G C'; C 0]; per in-service branch sqrt(g Cov g') in
## MW, g from `measurement_jacobian` on one P flow row at the from end, per
## bus the root of the diagonal entry of its Vm column in pu.
function _obs_cov_reference(net, meas; constraint::Bool = false)
  stored(f) = with_state_estimation_config(f; flatstart = false, jac_eps = 1e-6)
  netc = deepcopy(net)
  empty!(netc.measurements)
  append!(netc.measurements, meas)
  jm = stored(() -> measurement_jacobian(netc))
  H = Matrix{Float64}(jm.H)
  w = Float64[1.0 / meas[r.index].sigma^2 for r in jm.rows]
  con = constraint ? [i for (i, r) in enumerate(jm.rows) if meas[r.index].zero_injection] : Int[]
  w[con] .= 0.0
  G = H' * (w .* H)
  n = size(G, 1)
  C = H[con, :]
  Cov = isempty(con) ? inv(G) : inv([G C'; C zeros(length(con), length(con))])[1:n, 1:n]
  brs = [k for (k, br) in enumerate(net.branchVec) if br.status == 1]
  netf = deepcopy(net)
  empty!(netf.measurements)
  for k in brs
    push!(netf.measurements, Measurement(typ = Sparlectra.PflowMeas, value = 0.0, sigma = 1.0, branchIdx = k, direction = :from, id = "g$(k)"))
  end
  Gf = Matrix{Float64}(stored(() -> measurement_jacobian(netf)).H)
  nb = length(net.nodeVec)
  flow = [sqrt(max(LinearAlgebra.dot(Gf[k, :], Cov * Gf[k, :]), 0.0)) for k in eachindex(brs)]
  vm = [sqrt(Cov[nb - 1 + b, nb - 1 + b]) for b = 1:nb]
  return (brs = brs, flow = flow, vm = vm, cond_g = LinearAlgebra.cond(G))
end

## The direction reference of the singular fallback (`_weak_directions`):
## W^(1/2) H, a dense SVD, per direction the largest branch-flow standard
## deviation |g_k . v_i| / sigma_i and the largest magnitude one |v_ji| /
## sigma_i. Returns the affected buses (the ends of the branches whose flow
## standard deviation over the weak directions exceeds `wf`, and the buses
## whose magnitude standard deviation exceeds `wv`), those branches with
## their standard deviation, the largest branch-flow standard deviation
## over the weak directions (per branch, as on the covariance path), and the buses of the rule of
## 790900bc (angle and magnitude rows of the weak directions scaled by
## 1/sigma_i above rank_tol_factor * jac_eps).
function _obs_weak_reference(net, meas; wf = 20.0, wv = 0.01)
  stored(f) = with_state_estimation_config(f; flatstart = false, jac_eps = 1e-6)
  netc = deepcopy(net)
  empty!(netc.measurements)
  append!(netc.measurements, meas)
  jm = stored(() -> measurement_jacobian(netc))
  sig = Float64[meas[r.index].sigma for r in jm.rows]
  A = Matrix{Float64}(jm.H) ./ sig
  brs = [k for (k, br) in enumerate(net.branchVec) if br.status == 1]
  netf = deepcopy(net)
  empty!(netf.measurements)
  for k in brs
    push!(netf.measurements, Measurement(typ = Sparlectra.PflowMeas, value = 0.0, sigma = 1.0, branchIdx = k, direction = :from, id = "g$(k)"))
  end
  G = Matrix{Float64}(stored(() -> measurement_jacobian(netf)).H)
  F = LinearAlgebra.svd(A)
  nb = length(net.nodeVec)
  slack = findfirst(nd -> getNodeType(nd) == Sparlectra.Slack, net.nodeVec)
  GV = G * F.V
  s = [maximum(abs, GV[:, i]) / F.S[i] for i in eachindex(F.S)]
  t = [maximum(abs, F.V[nb:(2nb-1), i]) / F.S[i] for i in eachindex(F.S)]
  weak = [i for i in eachindex(F.S) if s[i] > wf || t[i] > wv]
  R = F.V[:, weak] ./ F.S[weak]'
  std = [LinearAlgebra.norm(GV[k, weak] ./ F.S[weak]) for k in eachindex(brs)]
  hit = [k for k in eachindex(brs) if std[k] > wf]
  vmstd = [LinearAlgebra.norm(R[nb - 1 + b, :]) for b = 1:nb]
  buses = sort(unique(vcat([b for b = 1:nb if vmstd[b] > wv], [Int(net.branchVec[brs[k]].fromBus) for k in hit], [Int(net.branchVec[brs[k]].toBus) for k in hit])))
  cfg = _obs_cfg(state_estimation_config)
  tau = cfg.rank_tol_factor * cfg.jac_eps
  oldBuses = [b for b = 1:nb if LinearAlgebra.norm(vcat(b == slack ? Float64[] : R[b < slack ? b : b - 1, :], R[nb - 1 + b, :])) > tau]
  return (buses = buses, branches = brs[hit], std = std[hit], smax = maximum(std; init = 0.0), old_buses = oldBuses)
end

## `_weak_directions` (the fallback on a singular gain matrix) called the
## way `_weak_network_island` calls it, at the stored state
function _obs_weak_directions(net, meas, dense_max_n)
  return _obs_cfg() do
    cfg = Sparlectra.state_estimation_config()
    prep, act, _ = Sparlectra._observability_rows(net, meas, "test")
    snet = prep.snet
    sl = Sparlectra._find_slack_idx(snet)
    Y = Sparlectra.createYBUS(net = snet)
    x = Sparlectra._initial_state_vector(snet, sl; flatstart = false)
    H = SparseArrays.sparse(first(Sparlectra._measurement_jacobian_fd(act, snet, x, sl, length(snet.nodeVec), Y; eps = cfg.jac_eps)))
    Sparlectra._weak_directions(H, act, snet, x, sl, Y, cfg.jac_eps, false, 20.0, 0.01, dense_max_n)
  end
end

## bounds below every standard deviation: every in-service branch is
## listed with its flow standard deviation (the values against a reference)
const _OBS_WEAK_ALL = Sparlectra.SparlectraConfig(Dict{String,Any}("state_estimation" => Dict{String,Any}("observability" => Dict{String,Any}("weak_flow_sigma_mw" => 1.0e-9, "weak_vm_sigma_pu" => 1.0e-12))))

function test_observability_weak(fx)::Bool
  # Weakly observable (concept K3, G2, G3, U4, U5) on the shipped
  # sp_case118 with Vm, P and Q at every bus and no flows, at the stored
  # (power-flow) state of the net. The class comes from the covariance of
  # the estimate: per branch the flow standard deviation sqrt(g G^-1 g'),
  # per bus the magnitude one, from the Takahashi selected inverse. With the
  # default sigmas (1 MW) the largest flow standard deviation is a few MW
  # and the set is observable; with P and Q at 30 MW (pseudo measurements)
  # it rises above the default bound of 20 MW. Every flow and magnitude
  # standard deviation equals the dense inv(G) reference (bounds below all
  # of them list every branch), on the dense and the sparse path, with no
  # covariance entry outside the pattern of the selected inverse; the
  # affected buses are the ends of the branches above the bound and the
  # buses above the magnitude bound.
  @testset "weak class: sp_case118 covariance against a dense reference" begin (function ()
    net = Sparlectra.createNetFromMatPowerFile(filename = joinpath(pkgdir(Sparlectra), "data", "mpower", "sp_case118.m"), flatstart = false, bus_shunt_model = :admittance, matpower_shift_sign = 1.0, matpower_shift_unit = :deg, matpower_ratio = :normal, tap_changer_model = :ideal)
    redirect_stdout(() -> runpf!(net, 40, 1e-8, 0), devnull)
    inj = generateMeasurementsFromPF(net; includeVm = true, includePinj = true, includeQinj = true, includePflow = false, includeQflow = false, noise = false)
    good = _obs_cfg(() -> evaluate_global_observability(net, inj))
    refG = _obs_cov_reference(net, inj)
    @test good.quality in (:good, :critical) && !good.weakly_observable
    @test good.weak_observability.basis === :stored_state
    @test good.weak_observability.n_buses == 0 && good.weak_observability.max_flow_sigma_mw < 20.0
    @test isapprox(good.weak_observability.max_flow_sigma_mw, maximum(refG.flow); rtol = 1e-8)
    @test isapprox(good.weak_observability.max_vm_sigma_pu, maximum(refG.vm); rtol = 1e-8)
    coarse = [Measurement(typ = m.typ, value = m.value, sigma = m.typ in (Sparlectra.PinjMeas, Sparlectra.QinjMeas) ? 30.0 : m.sigma, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id, origin = m.origin) for m in inj]
    ref = _obs_cov_reference(net, coarse)
    hit = [k for k in eachindex(ref.brs) if ref.flow[k] > 20.0]
    buses = sort(unique(vcat([b for b in eachindex(ref.vm) if ref.vm[b] > 0.01], [Int(net.branchVec[ref.brs[k]].fromBus) for k in hit], [Int(net.branchVec[ref.brs[k]].toBus) for k in hit])))
    @test !isempty(hit)
    for dmn in (Sparlectra._SE_DENSE_LINALG_MAX_N, 0)
      o = _obs_cfg(() -> evaluate_global_observability(net, coarse; dense_max_n = dmn))
      w = o.weak_observability
      @test (dmn, o.quality, o.weakly_observable, o.numerical_observable, o.angle_rank_deficit, o.magnitude_rank_deficit) == (dmn, :weakly_observable, true, true, 0, 0)
      @test (dmn, w.basis, w.n_branches, w.branches, w.buses, w.n_buses, w.pattern_misses) == (dmn, :stored_state, length(hit), ref.brs[hit], buses, length(buses), 0)
      @test (dmn, only(w.network_islands).method) == (dmn, :covariance)
      @test isapprox(w.branch_flow_sigma_mw, ref.flow[hit]; rtol = 1e-8)
      @test isapprox(w.max_flow_sigma_mw, maximum(ref.flow); rtol = 1e-8)
      # the old fields of the verdict do not change with the class
      @test (o.numerical_rank, o.unobservable_state_columns, isempty(o.observable_islands)) == (o.n_states, Int[], true)
      # the standalone call at the same state gives the same report
      ws = _obs_cfg(() -> evaluate_weak_observability(net, coarse; dense_max_n = dmn))
      @test ws.weakly_observable && (ws.n_branches, ws.buses, ws.branches) == (w.n_branches, w.buses, w.branches)
    end
    # every branch and bus against the reference
    all118 = with_sparlectra_config(() -> _obs_cfg(() -> evaluate_weak_observability(net, coarse)), _OBS_WEAK_ALL)
    @test all118.branches == ref.brs && isapprox(all118.branch_flow_sigma_mw, ref.flow; rtol = 1e-8)
    @test all118.buses == collect(eachindex(ref.vm)) && isapprox(all118.max_vm_sigma_pu, maximum(ref.vm); rtol = 1e-8)
    # The singular fallback (`_weak_directions`, a gain matrix without
    # covariance) against the SVD reference on the same set, on both paths;
    # its affected buses are a proper subset of the buses the D1 tolerance
    # flagged in 790900bc, and its branches a subset of the covariance ones
    # (the directions add only the weak part of the covariance)
    dref = _obs_weak_reference(net, coarse)
    @test issubset(dref.buses, dref.old_buses) && length(dref.buses) < length(dref.old_buses)
    @test issubset(dref.branches, ref.brs[hit])
    for dmn in (Sparlectra._SE_DENSE_LINALG_MAX_N, 0)
      d = _obs_weak_directions(net, coarse, dmn)
      @test (dmn, d.method, d.buses, d.branches) == (dmn, :directions, dref.buses, dref.branches)
      @test isapprox(d.branch_flow_sigma_mw, dref.std; rtol = 1e-4) && isapprox(d.max_flow_sigma_mw, dref.smax; rtol = 1e-4)
    end
    # a net that stores the flat profile: the class is not determined
    flat = deepcopy(net)
    for nd in flat.nodeVec
      nd._vm_pu = 1.0
      nd._va_deg = 0.0
    end
    wf = _obs_cfg(() -> evaluate_weak_observability(flat, coarse))
    @test (wf.basis, wf.weakly_observable, wf.n_buses) == (:not_determined, false, 0)
    @test _obs_cfg(() -> evaluate_global_observability(flat, coarse)).weak_observability.basis === :not_determined
    # the bounds come from state_estimation.observability: a file value
    # parses, a non-positive one is refused with the key named, and a bound
    # below the largest flow standard deviation of the default set makes
    # that set weakly observable as well
    cfg = Sparlectra.SparlectraConfig(Dict{String,Any}("state_estimation" => Dict{String,Any}("observability" => Dict{String,Any}("weak_flow_sigma_mw" => 2.0, "weak_vm_sigma_pu" => 0.5, "weak_check_after_estimation" => false))))
    @test (cfg.state_estimation.observability.weak_flow_sigma_mw, cfg.state_estimation.observability.weak_vm_sigma_pu, cfg.state_estimation.observability.weak_check_after_estimation) == (2.0, 0.5, false)
    @test Sparlectra.ObservabilityConfig().weak_check_after_estimation
    @test_throws "state_estimation.observability.weak_flow_sigma_mw" Sparlectra.SparlectraConfig(Dict{String,Any}("state_estimation" => Dict{String,Any}("observability" => Dict{String,Any}("weak_flow_sigma_mw" => 0.0))))
    @test_throws "state_estimation.observability.weak_vm_sigma_pu" Sparlectra.SparlectraConfig(Dict{String,Any}("state_estimation" => Dict{String,Any}("observability" => Dict{String,Any}("weak_vm_sigma_pu" => -1.0))))
    tight = with_sparlectra_config(cfg) do
      _obs_cfg(() -> evaluate_global_observability(net, inj))
    end
    @test tight.quality === :weakly_observable && tight.weak_observability.n_buses > 0
    @test tight.weak_observability.max_flow_sigma_mw == good.weak_observability.max_flow_sigma_mw
    # after the estimation (weak_check_after_estimation, default on):
    # runse! takes the class at the estimate it wrote into the net, the
    # same as the standalone call on that net; with the key off the field
    # is nothing
    tight2 = Sparlectra.SparlectraConfig(Dict{String,Any}("state_estimation" => Dict{String,Any}("observability" => Dict{String,Any}("weak_flow_sigma_mw" => 2.0, "weak_vm_sigma_pu" => 0.5))))
    est = deepcopy(net)
    res = with_sparlectra_config(() -> _obs_cfg(() -> runse!(est, inj)), tight2)
    @test res.converged && res.weakObservability !== nothing
    @test res.weakObservability.basis === :estimated_state && res.weakObservability.weakly_observable
    again = with_sparlectra_config(() -> _obs_cfg(() -> evaluate_weak_observability(est, inj)), tight2)
    @test (again.basis, again.n_branches, again.branches) == (:stored_state, res.weakObservability.n_branches, res.weakObservability.branches)
    plain = _obs_cfg(() -> runse!(deepcopy(net), inj))
    @test plain.weakObservability !== nothing && !plain.weakObservability.weakly_observable && plain.weakObservability.basis === :estimated_state
    @test with_sparlectra_config(() -> _obs_cfg(() -> runse!(deepcopy(net), inj)), cfg).weakObservability === nothing
  end)() end

  # The IEEE 14-bus topology (net14z, slack on bus 8) with Vm, P and Q at
  # every bus and no flows; bus 7 is passive, its injection rows are the
  # zero-injection rows (sigma 1e-3 MW). zero_injection = on weights them,
  # the covariance is inv(G); zero_injection = constraint solves them as
  # equality constraints, weight 0 in G, and the covariance of the
  # constrained estimate is the upper-left block of the inverse of the KKT
  # matrix [G C'; C 0] (derivation at `_weak_covariance_factor`). Every flow
  # and magnitude standard deviation equals the dense reference of its
  # mode; the two modes agree to the size of the zero-injection sigma, and
  # without the zero-injection rows the flows at bus 7 are looser.
  @testset "weak class: case14 covariance, weighted and constrained zero injections" begin (function ()
    net = fx.net14z
    base = filter(m -> !(_obs_is_injection(m) && m.busIdx == 7), generateMeasurementsFromPF(net; includeVm = true, includePinj = true, includeQinj = true, includePflow = false, includeQflow = false, noise = false))
    meas = copy(base)
    addZeroInjectionMeasurements!(meas; net = net, busIdxs = [7])
    vals = Dict{Symbol,Any}()
    for mode in (:on, :constraint)
      w = with_sparlectra_config(() -> with_state_estimation_config(() -> evaluate_weak_observability(net, meas); flatstart = true, jac_eps = 1e-6, zero_injection = mode), _OBS_WEAK_ALL)
      ref = _obs_cov_reference(net, meas; constraint = mode === :constraint)
      @test (mode, w.branches, w.pattern_misses) == (mode, ref.brs, 0)
      @test isapprox(w.branch_flow_sigma_mw, ref.flow; rtol = 1e-7)
      @test isapprox(w.max_vm_sigma_pu, maximum(ref.vm); rtol = 1e-7)
      vals[mode] = w.branch_flow_sigma_mw
    end
    @test isapprox(vals[:on], vals[:constraint]; rtol = 1e-4)
    loose = with_sparlectra_config(() -> _obs_cfg(() -> evaluate_weak_observability(net, base)), _OBS_WEAK_ALL)
    at7 = [k for (k, br) in enumerate(net.branchVec) if br.status == 1 && 7 in (Int(br.fromBus), Int(br.toBus))]
    @test all(k -> loose.branch_flow_sigma_mw[k] > vals[:constraint][k], at7)
  end)() end
  return true
end

## Measurement placement (concept test 9, A4, N4) on the shipped IEEE
## 14-bus topology with its slack on bus 8 (`net14z`) and on
## two_islands_prio. Every result is checked by evaluate_global_observability
## on the set plus the suggested rows, not by the placement's own counts.
function test_observability_placement(fx)::Bool
  # Vm rows only: the Q-V model has full rank (every magnitude has its own
  # row), the P-theta model has no row, so all 13 non-slack angles are free
  # (deficit 13 + 0). One P injection lowers the P-theta deficit by at most
  # one, so :injection needs exactly 13 steps. The first step is bus 8: it
  # hangs on 7-8 alone, its injection is the flow on 7-8, and it is the
  # only P row that makes a branch observable on its own (the plain pass
  # orders by that, ties by bus). With flows every candidate adds one P row
  # as well, so 13 is a lower bound for :flow; the count itself is
  # computed. A PMU angle is theta_i plus the PMU offset (pmu_ref_offset =
  # auto, the default): a set without a PMU has no offset column, so the
  # first PMU fixes only the offset and lowers no deficit, every further one
  # lowers it by one, 14 steps. Under pmu_ref_offset = off there is no
  # offset column and a PMU at the slack row adds nothing: 13 steps.
  @testset "placement: case14 with Vm rows only" begin (function ()
    net = fx.net14z
    vm = filter(m -> m.typ == Sparlectra.VmMeas, _obs_full_set(net))
    o0 = _obs_cfg(() -> evaluate_global_observability(net, vm))
    @test (o0.angle_rank_deficit, o0.magnitude_rank_deficit) == (13, 0)
    observable_with(rows) = _obs_cfg(() -> evaluate_global_observability(net, vcat(vm, rows))).numerical_observable
    inj = _obs_cfg(() -> suggest_measurement_placement(net, vm; kind = :injection))
    @test inj.reached && inj.stop_reason === :reached
    @test length(inj.steps) == o0.angle_rank_deficit + o0.magnitude_rank_deficit
    @test all(s -> s.quantity === :P && s.angle_rank_deficit_after == s.angle_rank_deficit_before - 1, inj.steps)
    @test first(inj.steps).bus == 8 && first(inj.steps).unobservable_branches_before == 20
    @test (inj.angle_rank_deficit_after, inj.magnitude_rank_deficit_after) == (0, 0)
    @test observable_with(inj.rows) && !observable_with(inj.rows[1:end-1])
    @test all(r -> startswith(r.id, "placement_Pinj_bus_"), inj.rows)
    fl = _obs_cfg(() -> suggest_measurement_placement(net, vm; kind = :flow))
    @test fl.reached && length(fl.steps) >= 13 && observable_with(fl.rows)
    @test all(s -> s.quantity === :PQ && s.branch > 0 && s.bus == 0, fl.steps) && length(fl.rows) == 2 * length(fl.steps)
    pm = _obs_cfg(() -> suggest_measurement_placement(net, vm; kind = :pmu))
    @test pm.reached && observable_with(pm.rows)
    @test length(pm.steps) == 14
    @test (pm.steps[1].angle_rank_deficit_before, pm.steps[1].angle_rank_deficit_after) == (13, 13)
    @test all(s -> s.angle_rank_deficit_after == s.angle_rank_deficit_before - 1, pm.steps[2:end])
    off = with_state_estimation_config(() -> suggest_measurement_placement(net, vm; kind = :pmu); flatstart = true, jac_eps = 1e-6, pmu_ref_offset = :off)
    @test off.reached && length(off.steps) == 13
    # max_steps stops early and says so; the configured default reaches the
    # call
    short = _obs_cfg(() -> suggest_measurement_placement(net, vm; kind = :injection, max_steps = 3))
    @test (length(short.steps), short.reached, short.stop_reason, short.angle_rank_deficit_after) == (3, false, :max_steps, 10)
    capped = with_state_estimation_config(() -> suggest_measurement_placement(net, vm; kind = :injection); flatstart = true, jac_eps = 1e-6, observability = Sparlectra.ObservabilityConfig(placement_max_steps = 5))
    @test (length(capped.steps), capped.stop_reason) == (5, :max_steps)
    io = IOBuffer()
    print_measurement_placement(io, short)
    @test countlines(seekstart(io)) == 4
  end)() end

  # target :no_critical on an observable set: Vm everywhere and P at every
  # bus but 14. Observable (the 13 P rows fix the 13 angles), but every row
  # is critical (no row has a second). The placement continues until the
  # criticality check of the verdict finds no critical row; the count is
  # computed, the result is checked by evaluate_global_observability.
  @testset "placement: no critical measurement left" begin (function ()
    net = fx.net14z
    full = _obs_full_set(net)
    set = filter(m -> m.typ == Sparlectra.VmMeas || (m.typ == Sparlectra.PinjMeas && m.busIdx != 14), full)
    o0 = _obs_cfg(() -> evaluate_global_observability(net, set))
    @test o0.numerical_observable && o0.quality === :critical
    r = _obs_cfg(() -> suggest_measurement_placement(net, set; kind = :injection, target = :no_critical))
    @test r.reached && !isempty(r.steps) && all(s -> s.phase === :no_critical && s.critical_after < s.critical_before, r.steps)
    o1 = _obs_cfg(() -> evaluate_global_observability(net, vcat(set, r.rows)))
    @test isempty(o1.numerical_critical_measurement_indices) && isempty(o1.structural_critical_measurement_indices)
  end)() end

  # N4, N5: the set of the network-island test (network island 1 with a
  # P-theta and a Q-V deficit at bus 2, network island 2 observable,
  # network island 3 not measured). No step leaves its network island,
  # network islands 2 and 3 get none. :injection and :flow close the gap;
  # :pmu closes the angle and stops on the magnitude of bus 2, which no
  # angle row reaches.
  @testset "placement: no candidate crosses a network island" begin (function ()
    net = deepcopy(fx.net_two)
    full = _obs_full_set(net)
    setBranchStatus!(net.branchVec[5], false)
    setBranchStatus!(net.branchVec[7], false)
    meas = filter(m -> !(m.busIdx in (2, 5, 6)) && !(m.branchIdx !== nothing && m.branchIdx in (1, 2, 5, 6, 7)) && !(_obs_is_injection(m) && m.busIdx in (1, 3)), full)
    isl = _obs_cfg(() -> identify_observable_islands(net, meas))
    island_of(s) = s.branch == 0 ? isl.network_island_of_bus[s.bus] : isl.network_island_of_bus[Int(net.branchVec[s.branch].fromBus)]
    for kind in (:injection, :flow, :pmu)
      r = _obs_cfg(() -> suggest_measurement_placement(net, meas; kind = kind))
      @test (kind, [x.network_island for x in r.network_islands]) == (kind, [1, 2])
      @test (kind, all(s -> s.network_island == 1 && island_of(s) == 1, r.steps)) == (kind, true)
      @test (kind, r.network_islands[2].n_steps) == (kind, 0)
      o = _obs_cfg(() -> evaluate_global_observability(net, vcat(meas, r.rows)))
      @test (kind, r.reached, o.numerical_observable) == (kind, kind !== :pmu, kind !== :pmu)
      kind === :pmu && @test (r.stop_reason, o.angle_rank_deficit, o.magnitude_rank_deficit) == (:no_candidate, 0, 1)
    end
    # a magnitude gap alone (the dark bus 7 of the M1 set, no PMU in the
    # set): a PMU enters only the P-theta model, whose deficit is 0, so
    # :pmu takes no step (1.1.0 took one PMU that only fixed the offset)
    cut = filter(m -> !(m.typ == Sparlectra.VmMeas && m.busIdx == 7) && !(m.typ == Sparlectra.QinjMeas && m.busIdx in (4, 7, 8, 9)) && !(m.typ == Sparlectra.QflowMeas && m.branchIdx in (8, 14, 15)), _obs_full_set(fx.net14z))
    pm = _obs_cfg(() -> suggest_measurement_placement(fx.net14z, cut; kind = :pmu))
    @test (length(pm.steps), pm.reached, pm.stop_reason, pm.angle_rank_deficit_before, pm.magnitude_rank_deficit_after) == (0, false, :no_candidate, 0, 1)
  end)() end

  return true
end

function run_observability_tests()

  @testset "Observability" begin (function ()
    fx = observability_fixture()
    tests =
      [("Observability metrics", () -> test_observability_metrics(fx)), ("Matrix observability helpers", test_observability_matrix_helpers), ("FD-aware rank tolerance", () -> test_observability_fd_rank_tolerance(fx)), ("Takahashi diagnostics", () -> test_observability_takahashi_diagnostics(fx)), ("Observable islands", () -> test_observability_islands(fx)), ("Restoration and weak gaps", () -> test_observability_restore(fx)), ("Measurement placement", () -> test_observability_placement(fx)), ("Weak class", () -> test_observability_weak(fx))]

    for (name, testfn) in tests
      @testset "$name" begin (function ()
        @test _se_run_quiet(testfn) == true
      end)() end
    end
  end)() end
end

function run_observability_rank_tests()
  # the rank read from the LDLt pivots against the decomposition
  @testset "observability: the rank from the LDLt pivots agrees with the decomposition (#399)" begin (function ()
    # the pivot rule reads the rank from the factorization the criticality
    # pass needs anyway; it is selectable and must agree with the SVD/QR
    # verdict on every case here, a disagreement is a finding
    both(H; tol = nothing) = begin
      a = evaluate_observability_matrix(H; tol = tol, rank_method = :decomposition)
      b = evaluate_observability_matrix(H; tol = tol, rank_method = :pivots)
      @test a.rank_method === :decomposition
      @test b.numerical_rank == a.numerical_rank
      @test b.numerical_observable == a.numerical_observable
      b
    end
    # full rank with one critical row
    Hc = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 1.0 0.0; 1.0 1.0 0.0]
    @test both(Hc).rank_method === :pivots
    @test both(Hc).numerical_rank == 3
    # an unmeasured state (a dark column, the unmeasured island in matrix
    # form): the gain matrix is singular, the pivot path states the deficit
    # through the decomposition and says so
    Hd = [1.0 0.0 0.0; 0.0 1.0 0.0; 1.0 1.0 0.0]
    bd = both(Hd)
    @test bd.numerical_rank == 2 && !bd.numerical_observable
    @test bd.rank_method === :decomposition
    # a weight spread of 1e6 in the rows (sigma 1e-4 next to 1.0 scales a
    # row by 1e4) changes neither verdict
    Hw = copy(Hc)
    Hw[1, :] .*= 1.0e4
    Hw[4, :] .*= 1.0e4
    @test both(Hw).rank_method === :pivots
    # FD noise: a column that carries nothing but a 1e-9 entry is dark at
    # the rank tolerance of the estimator (FD noise floor), for both rules
    Hn = [1.0 0.0 1.0e-9; 0.0 1.0 0.0; 1.0 1.0 0.0]
    bn = both(Hn; tol = 1.0e-6)
    @test bn.numerical_rank == 2 && !bn.numerical_observable
    # the fixture networks through the estimator's own path
    fx = observability_fixture()
    for (label, net_fx) in (("3-bus", fx.net3), ("sp_case14", fx.net14), ("link", fx.net_link), ("link with shunt", fx.net_link_shunt))
      # the fixture nets carry no measurements; a full set from the solved
      # state, as the other observability sets attach it
      net_a = deepcopy(net_fx)
      empty!(net_a.measurements)
      append!(net_a.measurements, generateMeasurementsFromPF(net_a; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false))
      net_b = deepcopy(net_a)
      obs_a = with_state_estimation_config(() -> evaluate_global_observability(net_a); rank_method = :decomposition)
      obs_b = with_state_estimation_config(() -> evaluate_global_observability(net_b); rank_method = :pivots)
      @test (label, obs_b.numerical_rank, obs_b.numerical_observable) == (label, obs_a.numerical_rank, obs_a.numerical_observable)
    end
    # the large shipped case: sparse QR against the CHOLMOD pivots
    big = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "sp_case1354.m"))
    netq = Sparlectra.createNetFromMatPowerFile(filename = big, flatstart = false, bus_shunt_model = :admittance, matpower_shift_sign = 1.0, matpower_shift_unit = :deg, matpower_ratio = :normal, tap_changer_model = :ideal)
    runpf!(netq, 60, 1e-8, 0)
    msq = generateMeasurementsFromPF(netq; includeVm = true, includePinj = true, includeQinj = true, includePflow = false, includeQflow = false, noise = true, rng = Xoshiro(399))
    empty!(netq.measurements)
    append!(netq.measurements, msq)
    obs_a = with_state_estimation_config(() -> evaluate_global_observability(netq); rank_method = :decomposition)
    obs_b = with_state_estimation_config(() -> evaluate_global_observability(netq); rank_method = :pivots)
    @test obs_b.numerical_rank == obs_a.numerical_rank
    @test obs_b.numerical_observable == obs_a.numerical_observable
    @test obs_b.rank_method === :pivots
    # the full set (P, Q and U at every bus): observable in both decoupled
    # models, not weakly observable at the default bounds at the stored
    # state (concept U4: the largest branch-flow standard deviation of the
    # covariance is 9.8 MW, below 20 MW)
    @test (obs_a.angle_rank_deficit, obs_a.magnitude_rank_deficit) == (0, 0) && obs_a.numerical_observable
    @test !obs_a.weakly_observable && obs_a.weak_observability.basis === :stored_state && obs_a.weak_observability.n_buses == 0
    @test 9.0 < obs_a.weak_observability.max_flow_sigma_mw < 11.0
    # The patterns of the evidence table (concept K2, E1, G4'): every bus
    # row of a few buses removed, the set has injections and Vm rows and no
    # flows. The P-theta deficit is the number of missing P injections minus
    # one (the injections of the other buses fix every angle difference but
    # one per missing injection, and the balance of the net closes one). The
    # Q-V deficit is 0: the Vm rows fix every other magnitude, and the Q
    # injections of the neighbours fix the missing ones. The AC rank is
    # weaker than the P-theta model: deficit 0 for the first two, the Q
    # rows carry angle information through the coupling (E1: flow errors up
    # to 1000 MW without 100:105, unreported in 1.0.0). The iteration
    # splits the net into many observable islands, most of them single
    # buses; the region is every non-slack bus.
    msfull = Measurement[m for m in netq.measurements]
    for (r, d, dac) in ((700:702, 2, 0), (100:105, 5, 0), (50:60, 10, 5))
      o = with_state_estimation_config(() -> evaluate_global_observability(netq, filter(m -> !(m.busIdx !== nothing && m.busIdx in r), msfull)); rank_method = :decomposition)
      @test (r, o.quality, o.angle_rank_deficit, o.magnitude_rank_deficit, o.rank_deficit, o.n_states - o.numerical_rank) == (r, :not_observable, d, 0, d, dac)
      @test (r, length(o.observable_islands.region_buses), o.observable_islands.n_islands > d) == (r, 1353, true)
    end
    # Restoration with scheduled values on two of the patterns (concept
    # test 8, U2, U3): without 700:702 the deficit 2 closes with the P
    # schedules of the first two candidates by bus, 700 (a generator of 180
    # MW, sigma 30 % = 54 MW) and 701 (a load of 6 MW, 1.8 MW); at the
    # stored state seven branches stay above the flow bound (largest 29.9
    # MW, 124 MW before), so the verdict is weakly observable after
    # restoring. Without 100:105 the buses 100 and 104 carry no load,
    # generator or shunt (passive): their zero injections come first, then
    # the P rows at 101, 102 and 103 close the five, and the set is
    # observable (13.6 MW). Under zero_injection = off the passive rows are
    # only reported, the schedules at 101, 102, 103 and 105 close four of
    # the five, one stays.
    for (r, buses, srcs, verdict) in ((700:702, [700, 701], [:schedule, :schedule], :weakly_observable), (100:105, [100, 104, 101, 102, 103], [:passive, :passive, :schedule, :schedule, :schedule], :observable))
      rr = restore_observability(netq, filter(m -> !(m.busIdx !== nothing && m.busIdx in r), msfull))
      @test (r, [x.bus for x in rr.restored], [x.source for x in rr.restored], all(x -> x.quantity === :P, rr.restored), rr.rank_deficit_after, rr.observable) == (r, buses, srcs, true, 0, true)
      @test (r, rr.verdict, rr.success, isempty(rr.weak.branches)) == (r, verdict, verdict === :observable, verdict === :observable)
      r == 700:702 && @test reduce(vcat, [[x.value, x.sigma] for x in rr.restored]) ≈ [180.0, 54.0, -6.0, 1.8]
    end
    ro = with_state_estimation_config(() -> restore_observability(netq, filter(m -> !(m.busIdx !== nothing && m.busIdx in 100:105), msfull)); zero_injection = :off)
    @test ([x.bus for x in ro.restored], [(x.bus, x.quantity) for x in ro.passive_not_inserted], ro.rank_deficit_after, ro.verdict) == ([101, 102, 103, 105], [(100, :P), (104, :P)], 1, :not_observable)
    # The kept weak estimate (concept G2, M2): without 700:702 runse! under
    # scheduled takes the two schedules, the estimate is weakly observable,
    # and without its affected buses no observable island passes its rank
    # check (the boundary injections leave with them), so the estimate is
    # kept. It carries the weak report (the service maps it to
    # se_result_usable = false), and the flows of the weak branches and the
    # total losses are NaN in the net, never zero. Without 100:105 the same
    # set is not weak at its estimate (13.6 MW at the stored state) and
    # does not reach this path, neither under scheduled nor under
    # observable_islands.
    netk = deepcopy(netq)
    rk = @test_logs (:info, r"no observable island passes its rank check, the estimate is kept") match_mode = :any with_state_estimation_config(() -> runse!(netk, filter(m -> !(m.busIdx !== nothing && m.busIdx in 700:702), msfull)); update_net = true, observability = Sparlectra.ObservabilityConfig(on_unobservable = :scheduled))
    wk = rk.weakObservability
    @test (rk.converged, rk.observableIslands === nothing, length(rk.restoredRows), wk.basis, wk.n_branches, wk.n_buses) == (true, true, 2, :estimated_state, 7, 10)
    @test all(j -> isnan(netk.branchVec[j].fBranchFlow.pFlow) && isnan(netk.branchVec[j].tBranchFlow.qFlow), wk.branches)
    @test all(isnan, netk.totalLosses[end])
    @test all(j -> isfinite(netk.branchVec[j].fBranchFlow.pFlow), [j for j in eachindex(netk.branchVec) if netk.branchVec[j].status == 1 && !(j in wk.branches)])
    # With every flow on a branch outside the removed region measured as
    # well (evidence E2) both models agree: deficits 0, 0 and 1, the full
    # deficit equals the AC deficit. Without 50:60 the one null direction
    # lies on buses 30 and 52, which hang on the region only, and the P
    # injection at bus 30 ties the two: the first pass finds the branches
    # 30-50, 30-52 and 50-52 unobservable, the iteration drops the P
    # injection at bus 30, and the observable islands are the rest (1352
    # buses), {30} and {52}; the region is {30, 52}.
    allq = generateMeasurementsFromPF(netq; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = true, rng = Xoshiro(1))
    inreg(m, r) = (m.busIdx !== nothing && m.busIdx in r) || (m.branchIdx !== nothing && (Int(netq.branchVec[m.branchIdx].fromBus) in r || Int(netq.branchVec[m.branchIdx].toBus) in r))
    for (r, d) in ((700:702, 0), (100:105, 0), (50:60, 1))
      o = with_state_estimation_config(() -> evaluate_global_observability(netq, filter(m -> !inreg(m, r), allq)); rank_method = :decomposition)
      @test (r, o.angle_rank_deficit, o.magnitude_rank_deficit, o.n_states - o.numerical_rank) == (r, d, 0, d)
      d == 0 && @test (r, o.numerical_observable, isempty(o.observable_islands)) == (r, true, true)
      if d > 0
        isl = o.observable_islands
        @test (length.(isl.island_buses), isl.island_buses[2:3], isl.region_buses) == ([1352, 1, 1], [[30], [52]], [30, 52])
        @test (isl.removed_boundary_injections, isl.rank_deficit) == (["Pinj_bus_30"], 1)
      end
    end
    # eleven buses without any measurement, the pattern 50:60
    dark = Set(50:60)
    filter!(m -> !(hasproperty(m, :busIdx) && m.busIdx in dark), netq.measurements)
    # Above the dense threshold (n = 2707) the AC rank and the weak
    # directions come from the gain-matrix iteration; the former warning
    # "dark-state identification skipped" is gone, and no other warning may
    # take its place. The AC deficit is five: the SVD of the same normalized
    # Jacobian (computed outside the suite, it is a dense 2707-column
    # decomposition) has the singular values 2.4e-15, 9.2e-7, 2.3e-6,
    # 3.4e-6 and 4.7e-6 below the tolerance 1.7e-5 and 4.4e-5 as the next
    # one. The column test of the sparse QR alone kept three of them (rank
    # 2705); the iteration on the gain matrix counts all five. The decoupled
    # models (1353 and 1354 columns) stay below the threshold: P-theta
    # deficit ten, Q-V deficit zero, no dark magnitude.
    def_a = @test_logs min_level = Logging.Warn with_state_estimation_config(() -> evaluate_global_observability(netq); rank_method = :decomposition)
    def_b = with_state_estimation_config(() -> evaluate_global_observability(netq); rank_method = :pivots)
    @test !def_a.numerical_observable
    @test def_a.n_states - def_a.numerical_rank == 5
    @test (def_a.angle_rank_deficit, def_a.magnitude_rank_deficit, def_a.rank_deficit) == (10, 0, 10)
    @test def_b.numerical_rank == def_a.numerical_rank
    @test def_b.rank_method === :decomposition
    @test !isempty(def_a.unobservable_state_columns)
    @test def_b.unobservable_state_columns == def_a.unobservable_state_columns
    isl = def_a.observable_islands
    @test (isl.angle_rank_deficit, isl.magnitude_rank_deficit, isl.rank_deficit) == (10, 0, 10)
    @test isempty(isl.dark_magnitude_buses) && length(isl.region_buses) == 1353
    @test _obs_islands_consistent(netq, isl)
    @test def_a.unobservable_state_columns == _obs_report_columns(netq, isl)
    @test (def_b.observable_islands.island_of_bus, def_b.observable_islands.region_buses, def_b.observable_islands.rank_deficit) == (isl.island_of_bus, isl.region_buses, isl.rank_deficit)
    # the SVD as the reference (the dense path with the threshold raised):
    # the same rank, deficits, islands and region
    ref = with_state_estimation_config(() -> evaluate_global_observability(netq; dense_max_n = 100_000); rank_method = :decomposition)
    iref = ref.observable_islands
    @test ref.numerical_rank == def_a.numerical_rank
    @test (iref.island_of_bus, iref.region_buses, iref.rank_deficit, iref.removed_boundary_injections) == (isl.island_of_bus, isl.region_buses, isl.rank_deficit, isl.removed_boundary_injections)
    @test (iref.dark_magnitude_buses, iref.angle_rank_deficit, ref.weak_observability) == (isl.dark_magnitude_buses, isl.angle_rank_deficit, def_a.weak_observability)
  end)() end

  return true
end
