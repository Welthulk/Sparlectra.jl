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
    @testset "Rectangular PF statuses use weak net keys" begin
      @test isempty(Sparlectra._RECTANGULAR_PF_STATUS.entries) || all(entry[2] isa WeakRef for entry in Sparlectra._RECTANGULAR_PF_STATUS.entries)

      net = createTest3BusNet()
      status = (status = :test_status,)
      @test Sparlectra._set_rectangular_pf_status!(net, status) === status
      @test Sparlectra.rectangular_pf_status(net) === status
    end
    @testset "Wrong-branch helper classification" begin
      @test isapprox(Sparlectra._circular_angle_spread_deg([179.0, -179.0]), 2.0; atol = 1e-8)
      @test isapprox(Sparlectra._circular_angle_spread_deg([-170.0, 170.0]), 20.0; atol = 1e-8)
      Vok = ComplexF64[1.0 + 0im, 1.01 + 0.01im, 1.03 - 0.02im]
      bus_types = [:PQ, :PV, :Slack]
      vset = [1.0, 1.01, 1.03]
      net = createTest3BusNet()
      ok = Sparlectra._check_wrong_branch_solution(net, Vok, bus_types, vset, 3; min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 180.0, max_branch_angle_deg = 90.0, min_low_vm_count = 1)
      @test ok.status == :ok

      Vlow = ComplexF64[0.40 + 0im, 1.0 + 0im, 1.03 + 0im]
      low = Sparlectra._check_wrong_branch_solution(net, Vlow, bus_types, vset, 3; min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 180.0, max_branch_angle_deg = 90.0, min_low_vm_count = 1)
      @test low.status == :warn
      @test low.reason == :low_voltage_magnitude

      Vbad = ComplexF64[NaN + 0im, 1.0 + 0im, 1.03 + 0im]
      bad = Sparlectra._check_wrong_branch_solution(net, Vbad, bus_types, vset, 3; min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 180.0, max_branch_angle_deg = 90.0, min_low_vm_count = 1)
      @test bad.status == :fail
      @test bad.reason == :nonfinite_voltage

      branch_net = Net(name = "branch_wrap_test", baseMVA = 100.0)
      addBus!(net = branch_net, busName = "B1", vn_kV = 110.0)
      addBus!(net = branch_net, busName = "B2", vn_kV = 110.0)
      addBus!(net = branch_net, busName = "B3", vn_kV = 110.0)
      addACLine!(net = branch_net, fromBus = "B1", toBus = "B2", length = 1.0, r = 0.0, x = 0.1, c_nf_per_km = 0.0, tanδ = 0.0)
      addACLine!(net = branch_net, fromBus = "B2", toBus = "B3", length = 1.0, r = 0.0, x = 0.1, c_nf_per_km = 0.0, tanδ = 0.0)
      setNetBranchStatus!(net = branch_net, branchNr = 2, status = 0)
      branch_net.branchVec[1].phase_shift_deg = 10.0
      branch_net.branchVec[1].angle = 10.0

      Vwrap = ComplexF64[exp(im * deg2rad(179.0)), exp(im * deg2rad(-179.0)), exp(im * deg2rad(0.0))]
      wrap = Sparlectra._check_wrong_branch_solution(Vwrap, [:PQ, :PQ, :Slack], [1.0, 1.0, 1.0], 3; net = branch_net, min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 360.0, max_branch_angle_deg = 20.0, min_low_vm_count = 1)
      @test wrap.status == :ok
      @test wrap.branch_angle_violation_count == 0
      @test !isnothing(wrap.worst_branch)
      @test wrap.worst_branch.angle_check_basis == :effective_bus_angle_minus_phase_shift
      @test isapprox(wrap.worst_branch.angle_diff_raw_deg, 2.0; atol = 1e-8)
      @test isapprox(wrap.worst_branch.angle_diff_effective_deg, 12.0; atol = 1e-8)

      strict = Sparlectra._check_wrong_branch_solution(Vwrap, [:PQ, :PQ, :Slack], [1.0, 1.0, 1.0], 3; net = branch_net, min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 360.0, max_branch_angle_deg = 11.0, min_low_vm_count = 1)
      @test strict.status == :warn
      @test strict.reason == :branch_angle_exceeded
      @test strict.branch_angle_violation_count == 1

      wrap_spread = Sparlectra._check_wrong_branch_solution(Vwrap, [:PQ, :PQ, :Slack], [1.0, 1.0, 1.0], 3; net = nothing, min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 10.0, max_branch_angle_deg = Inf, min_low_vm_count = 1)
      @test wrap_spread.status == :warn
      @test wrap_spread.reason == :angle_spread_exceeded
      @test isapprox(wrap_spread.angle_spread_deg, 181.0; atol = 1e-8)

      Vwide = ComplexF64[exp(im * deg2rad(120.0)), exp(im * deg2rad(-120.0)), exp(im * deg2rad(0.0))]
      wide_spread = Sparlectra._check_wrong_branch_solution(Vwide, [:PQ, :PQ, :Slack], [1.0, 1.0, 1.0], 3; net = nothing, min_vm_pu = 0.70, max_vm_pu = 1.30, max_angle_spread_deg = 30.0, max_branch_angle_deg = Inf, min_low_vm_count = 1)
      @test wide_spread.status == :warn
      @test wide_spread.reason == :angle_spread_exceeded
      @test isapprox(wide_spread.angle_spread_deg, 240.0; atol = 1e-8)
    end
    @testset "Flat-start voltage setpoints" begin
      net = createTest3BusNet()
      net.nodeVec[1]._vm_pu = 0.94
      net.nodeVec[1]._va_deg = -7.0
      net.nodeVec[2]._vm_pu = 1.04
      net.nodeVec[2]._va_deg = -3.0
      net.nodeVec[3]._vm_pu = 1.06
      net.nodeVec[3]._va_deg = 2.0

      Vflat, slack_idx = Sparlectra.initialVrect(net; flatstart = true)
      @test slack_idx == 3
      @test abs(Vflat[1]) == 1.0
      @test angle(Vflat[1]) == 0.0
      @test abs(Vflat[2]) == 1.04
      @test angle(Vflat[2]) == 0.0
      @test abs(Vflat[3]) == 1.06
      @test isapprox(rad2deg(angle(Vflat[3])), 2.0; atol = 1e-12)

      busVec, polar_slack_idx = Sparlectra.getBusData(net.nodeVec, net.baseMVA, true; net = net)
      @test polar_slack_idx == 3
      @test busVec[1].vm_pu == 1.0
      @test busVec[1].va_rad == 0.0
      @test busVec[2].vm_pu == 1.04
      @test busVec[2].va_rad == 0.0
      @test busVec[3].vm_pu == 1.06
      @test isapprox(rad2deg(busVec[3].va_rad), 2.0; atol = 1e-12)

      model = Sparlectra.buildPfModel(net; flatstart = true, include_limits = false, start_projection = false)
      @test model.busIdx_net == [1, 2, 3]
      @test abs(model.V0[1]) == 1.0
      @test angle(model.V0[1]) == 0.0
      @test abs(model.V0[2]) == 1.04
      @test angle(model.V0[2]) == 0.0
      @test abs(model.V0[3]) == 1.06
      @test isapprox(rad2deg(angle(model.V0[3])), 2.0; atol = 1e-12)

      Vseeded, _ = Sparlectra.initialVrect(net; flatstart = false)
      @test isapprox(abs(Vseeded[1]), 0.94; atol = 1e-12)
      @test isapprox(rad2deg(angle(Vseeded[1])), -7.0; atol = 1e-12)
      @test isapprox(abs(Vseeded[2]), 1.04; atol = 1e-12)
      @test isapprox(rad2deg(angle(Vseeded[2])), -3.0; atol = 1e-12)
      @test isapprox(abs(Vseeded[3]), 1.06; atol = 1e-12)
      @test isapprox(rad2deg(angle(Vseeded[3])), 2.0; atol = 1e-12)

      iso_net = Net(name = "iso_flatstart", baseMVA = 100.0)
      for busName in ("B1", "B2", "B3", "B4")
        addBus!(net = iso_net, busName = busName, vn_kV = 110.0)
      end
      addACLine!(net = iso_net, fromBus = "B1", toBus = "B3", length = 10.0, r = 0.0, x = 0.4, c_nf_per_km = 9.55, tanδ = 0.0)
      addACLine!(net = iso_net, fromBus = "B3", toBus = "B4", length = 10.0, r = 0.0, x = 0.4, c_nf_per_km = 9.55, tanδ = 0.0)
      addACLine!(net = iso_net, fromBus = "B4", toBus = "B1", length = 10.0, r = 0.0, x = 0.4, c_nf_per_km = 9.55, tanδ = 0.0)
      addProsumer!(net = iso_net, busName = "B4", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.06, va_deg = 5.0, referencePri = "B4")
      addProsumer!(net = iso_net, busName = "B3", type = "SYNCHRONOUSMACHINE", p = 20.0, q = 0.0, vm_pu = 1.04)
      addProsumer!(net = iso_net, busName = "B1", type = "ENERGYCONSUMER", p = 30.0, q = 10.0)
      setNodeVoltage!(net = iso_net, busName = "B4", vm_pu = 1.06, va_deg = 5.0)
      setPVBusVset!(net = iso_net, busName = "B3", vm_pu = 1.04)
      markIsolatedBuses!(net = iso_net, log = false)

      iso_model = Sparlectra.buildPfModel(iso_net; flatstart = true, include_limits = false, start_projection = false)
      @test iso_net.isoNodes == [2]
      @test iso_model.busIdx_net == [1, 3, 4]
      @test iso_model.busType == [:PQ, :PV, :Slack]
      @test abs.(iso_model.V0) == [1.0, 1.04, 1.06]
      @test all(isapprox.(rad2deg.(angle.(iso_model.V0)), [0.0, 0.0, 5.0]; atol = 1e-12))
      @test iso_model.Vset == [1.0, 1.04, 1.06]
    end
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
      @test occursin("switch comparisons use per-unit internally", printed)

      switch_io = IOBuffer()
      active_set_q_limits!(
        net,
        10,
        length(net.nodeVec);
        get_qreq_pu = bus -> bus == 2 ? 0.2 : 0.0,
        is_pv = bus -> bus == 2,
        make_pq! = (bus, qclamp, side) -> nothing,
        make_pv! = bus -> nothing,
        qmin_pu = fill(-Inf, length(net.nodeVec)),
        qmax_pu = [Inf, 0.1, Inf],
        pv_orig_mask = trues(length(net.nodeVec)),
        allow_reenable = false,
        q_hyst_pu = 0.0,
        cooldown_iters = 0,
        verbose = 1,
        io = switch_io,
      )
      switch_log = String(take!(switch_io))
      @test occursin("Q=0.200000 pu", switch_log)
      @test occursin("Qmax=0.100000 pu", switch_log)
      @test occursin("MVAr", switch_log)

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
      @test erg_locked == 1
      @test getNodeType(net_locked.nodeVec[2]) == Sparlectra.PV
    end

    @testset "Rectangular Q-limit trace diagnostics" begin
      net = createTest3BusNet()
      setQLimits!(net = net, qmin_MVar = -1.0, qmax_MVar = 1.0, busName = "STATION1")
      trace_path, trace_io = mktemp()
      close(trace_io)
      open(trace_path, "w") do io
        redirect_stdout(io) do
          runpf!(net, 3, 1e-8, 0; method = :rectangular, qlimit_start_iter = 99, qlimit_trace_buses = [2])
        end
      end
      trace_text = read(trace_path, String)
      rm(trace_path; force = true)
      @test occursin("Q-limit trace enabled", trace_text)
      @test occursin("BUS_I=2", trace_text)
      @test occursin("qlimit_start_iter has not been reached", trace_text)
      pq_io = IOBuffer()
      Sparlectra._print_rectangular_qlimit_trace(
        pq_io,
        net,
        4,
        2,
        :PQ,
        ComplexF64[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im],
        [1.0, 1.0, 1.0],
        NaN,
        [-Inf, -0.01, -Inf],
        [Inf, 0.01, Inf];
        q_hyst_pu = 0.0,
        qlimit_start_iter = 2,
        qlimit_start_mode = :iteration_or_auto,
        qlimit_iter_ready = true,
        qlimit_auto_ready = true,
        qlimit_ready = true,
        qlimit_check_active = true,
        converged_before_switching = false,
        qlimit_mode = :switch_to_pq,
        lock_mask = falses(3),
      )
      pq_trace_text = String(take!(pq_io))
      @test occursin("decision=keep PQ", pq_trace_text)
      @test occursin("reason=bus is no longer an active PV/REF bus", pq_trace_text)
      @test !occursin("decision=keep PV reason=bus is not an active PV/REF bus", pq_trace_text)

      ignore_io = IOBuffer()
      Sparlectra._print_rectangular_qlimit_trace(
        ignore_io, net, 1, 2, :PV, ComplexF64[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im], [1.0, 1.0, 1.0], 0.02, [-Inf, -0.01, -Inf], [Inf, 0.01, Inf];
        q_hyst_pu = 0.0, qlimit_start_iter = 1, qlimit_start_mode = :iteration, qlimit_iter_ready = true, qlimit_auto_ready = false, qlimit_ready = true, qlimit_check_active = true, converged_before_switching = false, qlimit_mode = :switch_to_pq, lock_mask = Bool[false, true, false], qlimit_lock_reason = :ignore_q_limits,
      )
      @test occursin("ignore_q_limits=true disables Q-limit switching", String(take!(ignore_io)))

      final_io = IOBuffer()
      res = Sparlectra._print_rectangular_qlimit_summary(
        final_io,
        net,
        ComplexF64[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im],
        ComplexF64[0.0 + 0.0im, 0.0 + 0.02im, 0.0 + 0.0im],
        [:PQ, :PV, :Slack],
        [-Inf, -0.01, -0.05],
        [Inf, 0.01, 0.05],
        zeros(Float64, 3);
        q_hyst_pu = 0.0,
        tolerance_pu = 1e-6,
        max_rows = 30,
      )
      final_text = String(take!(final_io))
      @test res.violating == 1
      @test res.pv_violations == 1
      @test res.ref_violations == 0
      @test occursin("Final PV Q-limit active-set check: FAIL", final_text)
      @test occursin("active PV buses checked", final_text)
      @test occursin("active REF buses checked", final_text)
      @test occursin("PV violations: 1", final_text)
      @test occursin("REF violations: 0", final_text)
      @test occursin("Final PV/REF Q-limit check: FAIL", final_text)
      @test occursin("Qcalc pu", final_text)
      @test occursin("switched", final_text)
      @test occursin("return band", final_text)

      ref_io = IOBuffer()
      ref_res = Sparlectra._print_rectangular_qlimit_summary(
        ref_io,
        net,
        ComplexF64[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im],
        ComplexF64[0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.08im],
        [:PQ, :PV, :Slack],
        [-Inf, -0.01, -0.05],
        [Inf, 0.01, 0.05],
        zeros(Float64, 3);
        q_hyst_pu = 0.0,
        tolerance_pu = 1e-6,
        max_rows = 30,
      )
      ref_text = String(take!(ref_io))
      @test ref_res.pv_violations == 0
      @test ref_res.ref_violations == 1
      @test ref_res.status == :warn
      @test occursin("Final PV Q-limit active-set check: OK", ref_text)
      @test occursin("Final REF Q-limit diagnostic: WARN", ref_text)
      @test occursin("PV violations: 0", ref_text)
      @test occursin("REF violations: 1", ref_text)
      @test occursin("Final PV/REF Q-limit check: WARN", ref_text)

  limited_io = IOBuffer()
  limited_res = Sparlectra._print_rectangular_qlimit_summary(
    limited_io,
    net,
    ComplexF64[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im],
    ComplexF64[0.0 + 0.0im, 0.0 + 0.02im, 0.0 + 0.08im],
    [:PQ, :PV, :Slack],
    [-Inf, -0.01, -0.05],
    [Inf, 0.01, 0.05],
    zeros(Float64, 3);
    q_hyst_pu = 0.0,
    tolerance_pu = 1e-6,
    max_rows = 30,
    max_console_rows = 1,
  )
  limited_text = String(take!(limited_io))
  @test limited_res.violating == 2
  @test occursin("1 more violation rows omitted", limited_text)
  @test occursin("max_rows/max_console_rows", limited_text)

  summary_io = IOBuffer()
  Sparlectra._print_rectangular_convergence_summary(summary_io, (
    numerical_converged = true,
    q_limit_active_set_ok = false,
    final_converged = false,
    reason_text = "remaining PV Q-limit violations",
  ))
  summary_text = String(take!(summary_io))
  @test occursin("numerical_solution=OK", summary_text)
  @test occursin("q_limit_active_set=FAIL", summary_text)
  @test occursin("final_converged=false", summary_text)
  @test occursin("reason=remaining PV Q-limit violations", summary_text)
  @test Sparlectra._rectangular_solver_status_symbol(true, false, false, :remaining_pv_q_limit_violations) == :converged_limits_failed
  @test Sparlectra._rectangular_solver_status_symbol(false, false, false, :singular_newton_step) == :singular_jacobian
  @test Sparlectra._rectangular_solver_status_symbol(false, false, false, :nr_mismatch_not_converged) == :not_converged
end

    @testset "Q-limit guard and active-set status" begin
      net = createTest3BusNet()
      bus = geNetBusIdx(net = net, busName = "STATION1")
      qmin_pu = [-Inf, 0.0, -Inf]
      qmax_pu = [Inf, 0.0, Inf]
      bus_types = [:PQ, :PV, :Slack]
      S = ComplexF64[0.0 + 0.0im, 0.1 + 0.0im, 0.0 + 0.0im]
      guarded = Sparlectra._apply_qlimit_guard_to_rectangular_active_set!(
        net,
        bus_types,
        S,
        zeros(Float64, 3),
        qmin_pu,
        qmax_pu;
        min_q_range_pu = 1e-4,
        zero_range_mode = :lock_pq,
        narrow_range_mode = :prefer_pq,
        log = false,
        verbose = 0,
      )
      @test guarded == [bus]
      @test bus_types[bus] == :PQ
      @test length(net.qLimitLog) == 1

      bus_types[bus] = :PV
      changed, reenabled = Sparlectra.active_set_q_limits!(
        net,
        2,
        3;
        get_qreq_pu = _ -> 0.2,
        is_pv = b -> (bus_types[b] == :PV),
        make_pq! = (b, _qclamp, _side) -> (bus_types[b] = :PQ),
        make_pv! = b -> (bus_types[b] = :PV),
        qmin_pu = qmin_pu,
        qmax_pu = qmax_pu,
        pv_orig_mask = trues(3),
        allow_reenable = false,
        q_hyst_pu = 0.0,
        cooldown_iters = 0,
        qlimit_guard_max_switches = 1,
        qlimit_guard_freeze_after_repeated_switching = true,
      )
      @test !changed
      @test !reenabled
      @test bus_types[bus] == :PV
      @test length(net.qLimitLog) == 1

      status_io = IOBuffer()
      status = (
        numerical_converged = true,
        q_limit_active_set_ok = false,
        final_mismatch = 2.1e-9,
        pv_pq_switching_events = 1842,
        qlimit_active_set_changes = 91,
        qlimit_reenable_events = 14,
        oscillating_buses = 37,
        guarded_narrow_q_pv_buses = 512,
        status = :qlimit_chatter,
      )
      Sparlectra._print_qlimit_active_set_summary(status_io, status)
      status_text = String(take!(status_io))
      @test occursin("NR convergence             : yes", status_text)
      @test occursin("Active-set convergence     : no", status_text)
      @test occursin("Final status               : qlimit_chatter", status_text)

      function zero_range_pv_net()
        guarded_net = createTest3BusNet()
        guarded_bus = geNetBusIdx(net = guarded_net, busName = "STATION1")
        guarded_net.nodeVec[guarded_bus]._nodeType = Sparlectra.PV
        empty!(guarded_net.qmin_pu)
        append!(guarded_net.qmin_pu, [-Inf, 0.0, -Inf])
        empty!(guarded_net.qmax_pu)
        append!(guarded_net.qmax_pu, [Inf, 0.0, Inf])
        return guarded_net
      end

      # Regression: lower-level solver callers must opt in before the narrow-Q guard
      # locks PV buses to PQ during rectangular pre-processing.
      default_net = zero_range_pv_net()
      redirect_stdout(devnull) do
        runpf!(default_net; config = PowerFlowConfig(max_iter = 0))
      end
      @test isempty(default_net.qLimitLog)

      opt_in_net = zero_range_pv_net()
      redirect_stdout(devnull) do
        runpf!(opt_in_net; config = PowerFlowConfig(max_iter = 0, qlimits = QLimitConfig(guard = true)))
      end
      @test length(opt_in_net.qLimitLog) == 1
    end

    @testset "Rectangular performance profile exposes solver control path" begin
      net = createTest3BusNet()
      profile = Dict{Symbol,Any}(:enabled => true, :show_allocations => false, :show_iteration_table => true)
      _, erg = runpf!(net, 20, 1e-8, 0; method = :rectangular, performance_profile = profile)
      @test erg == 0
      timings = profile[:timings]
      for phase in (
        :solver_bus_type_scan,
        :solver_qlimit_extraction,
        :solver_active_set_origin_mask,
        :solver_active_set_setup,
        :iteration_qreq_vector,
        :iteration_control_bookkeeping,
        :iteration_state_update,
        :solver_final_injection_vectors,
        :solver_result_bus_writeback,
        :solver_status_bookkeeping,
      )
        @test haskey(timings, phase)
        @test timings[phase].calls >= 1
      end
    end

    @testset "Wrong-branch detection modes are distinct" begin
      base_cfg = PowerFlowConfig(max_iter = 40, tol = 1e-8, start_mode = StartModeConfig(flatstart = true), wrong_branch_min_vm_pu = 1.20)

      net_warn = createTest3BusNet()
      _, erg_warn = runpf!(net_warn; config = PowerFlowConfig(; NamedTuple{fieldnames(PowerFlowConfig)}(getfield.(Ref(base_cfg), fieldnames(PowerFlowConfig)))..., wrong_branch_detection = :warn))
      st_warn = Sparlectra.rectangular_pf_status(net_warn)
      @test erg_warn == 0
      @test st_warn.wrong_branch_status == :warn

      net_fail = createTest3BusNet()
      _, erg_fail = runpf!(net_fail; config = PowerFlowConfig(; NamedTuple{fieldnames(PowerFlowConfig)}(getfield.(Ref(base_cfg), fieldnames(PowerFlowConfig)))..., wrong_branch_detection = :fail))
      st_fail = Sparlectra.rectangular_pf_status(net_fail)
      @test erg_fail != 0
      @test st_fail.wrong_branch_status == :fail

      net_rescue = createTest3BusNet()
      _, erg_rescue = runpf!(net_rescue; config = PowerFlowConfig(; NamedTuple{fieldnames(PowerFlowConfig)}(getfield.(Ref(base_cfg), fieldnames(PowerFlowConfig)))..., wrong_branch_detection = :rescue))
      st_rescue = Sparlectra.rectangular_pf_status(net_rescue)
      @test erg_rescue != 0
      @test st_rescue.wrong_branch_status == :wrong_branch_rescue_not_implemented
      @test st_rescue.wrong_branch_reason == :rescue_requested_but_not_available
      @test st_rescue.wrong_branch_rescue_attempted === false
    end

    @testset "Typed power-flow config entry points" begin
      @test isdefined(Main, :run_sparlectra)
      @test isdefined(Main, :run_acpflow)
      @test isdefined(Main, :SparlectraRunResult)
      @test isdefined(Sparlectra, :ensure_casefile)
      @test isdefined(Main, :ensure_casefile)
      @test Sparlectra.ensure_casefile === Sparlectra.FetchMatpowerCase.ensure_casefile

      pf_config = PowerFlowConfig(max_iter = 40, tol = 1e-8, sparse = true, start_mode = StartModeConfig(flatstart = true))

      net_direct = createTest3BusNet()
      _, erg_direct = runpf!(net_direct; config = pf_config)
      @test erg_direct == 0

      net_project = createTest3BusNet()
      _, erg_project = runpf!(net_project; config = SparlectraConfig(powerflow = pf_config))
      @test erg_project == 0

      net_rectangular = createTest3BusNet()
      _, erg_rectangular = runpf_rectangular!(net_rectangular, 20, 1e-8, 0; method = :rectangular)
      @test erg_rectangular == 0

      cfg_output_off = SparlectraConfig(powerflow = pf_config, output = OutputConfig(logfile_results = :off))
      net_runner = createTest3BusNet()
      result = run_sparlectra(net = net_runner, config = cfg_output_off)
      @test result isa SparlectraRunResult
      @test result.net === net_runner
      @test result.numerical_converged
      @test result.solution_available
      @test result.control_status === :none
      @test result.final_converged
      @test result.outcome in (:converged, :converged_with_limit_warnings)
      @test result.limit_validation_status in (:ok, :skip)
      @test result.reason isa Symbol
      @test result.reason_text isa String
      @test result.iterations >= 0
      @test result.final_mismatch isa Float64
      @test result.diagnostics isa NamedTuple

      fallback_result = SparlectraRunResult(
        result.net,
        :converged,
        true,
        true,
        :skip,
        true,
        :none,
        "none",
        0,
        0.0,
        nothing,
        NaN,
        :classic,
        :none,
        nothing,
        (converged = false, erg = 1, numerical_converged = false, solution_available = false, final_converged = false, outcome = :diagnostic_override, reason = :diagnostic_override),
      )
      fallback_row = Sparlectra._sparlectra_status_row(fallback_result)
      @test fallback_row.converged === true
      @test fallback_row.erg == 0
      @test fallback_row.numerical_converged === true
      @test fallback_row.solution_available === true
      @test fallback_row.final_converged === true
      @test fallback_row.outcome === :converged
      @test fallback_row.reason === :none
      @test !hasproperty(fallback_row, :q_limit_active_set_ok)

      net_alias = createTest3BusNet()
      alias_result = run_acpflow(net = net_alias, config = cfg_output_off)
      @test alias_result isa SparlectraRunResult
      @test alias_result.net === net_alias
      @test alias_result.outcome == result.outcome
      @test alias_result.method == result.method

      @test_throws TypeError run_sparlectra(net = createTest3BusNet(), config = pf_config)
      @test_throws TypeError run_acpflow(net = createTest3BusNet(), config = pf_config)

      net_with_explicit_cfg = createTest3BusNet()
      explicit_output = mktemp() do path, io
        redirect_stdout(io) do
          @test run_sparlectra(net = net_with_explicit_cfg, config = cfg_output_off).numerical_converged
        end
        flush(io)
        return read(path, String)
      end
      @test !occursin("AC Power Flow Results", explicit_output)

      cfg_output_full = SparlectraConfig(powerflow = pf_config, output = OutputConfig(logfile_results = :full))
      net_with_output = createTest3BusNet()
      full_output = mktemp() do path, io
        redirect_stdout(io) do
          @test run_sparlectra(net = net_with_output, config = cfg_output_full).numerical_converged
        end
        flush(io)
        return read(path, String)
      end
      @test occursin("AC Power Flow Results", full_output)
      @test Sparlectra._sparlectra_result_mode(net_with_output, OutputConfig(logfile_results = :compact, result_table_large_case_threshold_buses = 1, result_table_large_case_mode = :summary)) === :summary
    end

    @testset "Framework status composition" begin
      rejected_outcomes = (:wrong_branch_detected, :angle_spread_exceeded, :branch_angle_exceeded, :wrong_branch_rescue_not_implemented)
      for outcome in rejected_outcomes
        rect_status = (status = outcome, numerical_converged = true, reason = :none)
        @test !Sparlectra._rectangular_solution_available(rect_status)
      end
      rejected_reasons = (:wrong_branch_detected, :angle_spread_exceeded, :branch_angle_exceeded, :wrong_branch_rescue_not_implemented, :rescue_requested_but_not_available)
      for reason in rejected_reasons
        rect_status = (status = :converged, numerical_converged = true, reason = reason)
        @test !Sparlectra._rectangular_solution_available(rect_status)
      end
      @test !Sparlectra._rectangular_solution_available(nothing)
      @test !Sparlectra._rectangular_solution_available((status = :not_converged, numerical_converged = false, reason = :nr_mismatch_not_converged))
      @test Sparlectra._rectangular_solution_available((status = :converged_limits_failed, numerical_converged = true, reason = :remaining_q_limit_violations))

      base = (outcome = :converged, numerical_converged = true, solution_available = true, limit_validation_status = :ok, final_converged = true, reason = :none, reason_text = "none", final_mismatch = 0.0)
      for accepted_status in (:none, :converged, :disabled, :no_active_controllers, :no_controllers)
        accepted = Sparlectra._compose_framework_status(base, accepted_status)
        @test accepted.final_converged
        @test accepted.outcome === :converged
      end

      blocked = Sparlectra._compose_framework_status(base, :blocked)
      @test blocked.numerical_converged
      @test blocked.solution_available
      @test !blocked.final_converged
      @test blocked.outcome === :control_blocked
      @test blocked.reason === :control_blocked
      @test occursin("blocked", blocked.reason_text)

      exhausted = Sparlectra._compose_framework_status(base, :max_outer_iterations)
      @test exhausted.outcome === :control_max_outer_iterations
      @test !exhausted.final_converged

      unknown = Sparlectra._compose_framework_status(base, :custom_terminal_status)
      @test unknown.outcome === :control_not_converged
      @test !unknown.final_converged

      rejected_base = merge(base, (outcome = :wrong_branch_detected, solution_available = false, final_converged = false, reason = :wrong_branch_detected, reason_text = "wrong branch detected"))
      rejected = Sparlectra._compose_framework_status(rejected_base, :pf_failed)
      @test rejected.numerical_converged
      @test !rejected.solution_available
      @test !rejected.final_converged
      @test rejected.outcome === :wrong_branch_detected
      @test rejected.reason === :wrong_branch_detected
      @test occursin("wrong branch", lowercase(rejected.reason_text))

      failed_base = merge(base, (outcome = :not_converged, numerical_converged = false, solution_available = false, final_converged = false, reason = :nr_mismatch_not_converged, reason_text = "NR mismatch did not converge", final_mismatch = Inf))
      failed = Sparlectra._compose_framework_status(failed_base, :pf_failed)
      @test !failed.numerical_converged
      @test !failed.solution_available
      @test !failed.final_converged
      @test failed.outcome === :pf_failed
      @test failed.reason === :pf_failed
    end

    @testset "Configured MATPOWER batch parsing and deterministic local execution" begin
      cfg_cases = SparlectraConfig(Dict("matpower" => Dict("cases" => [" case_a.m ", "case_b.m"])))
      @test cfg_cases.matpower.cases == ["case_a.m", "case_b.m"]
      @test configured_matpower_cases(cfg_cases) == ["case_a.m", "case_b.m"]

      cfg_single = SparlectraConfig(Dict("matpower" => Dict("case" => " case_single.m ")))
      @test configured_matpower_cases(cfg_single) == ["case_single.m"]

      cfg_precedence = SparlectraConfig(Dict("matpower_import" => Dict("case" => "case_single.m", "cases" => ["case_a.m", "case_b.m"])))
      @test configured_matpower_cases(cfg_precedence) == ["case_a.m", "case_b.m"]
      @test_throws ArgumentError SparlectraConfig(Dict("matpower" => Dict("cases" => "case_a.m")))
      @test_throws ArgumentError SparlectraConfig(Dict("matpower" => Dict("cases" => ["case_a.m", " "])))
      @test_throws ArgumentError run_sparlectra_cases(config = SparlectraConfig(matpower = MatpowerImportConfig(case = "")))
      @test_throws ArgumentError run_sparlectra_cases(config = cfg_single, performance_profile = Dict{Symbol,Any}())

      mktempdir() do tmpdir
        for (case_name, load_mw) in [("case_a.m", 40), ("case_b.m", 60)]
          function_name = splitext(case_name)[1]
          write(joinpath(tmpdir, case_name), """
function mpc = $(function_name)
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
1 3 0 0 0 0 1 1.0 0 110 1 1.1 0.9;
2 1 $(load_mw) 15 0 0 1 1.0 0 110 1 1.1 0.9;
];
mpc.gen = [
1 100 0 300 -300 1.02 100 1 300 0;
];
mpc.branch = [
1 2 0.01 0.05 0.0 999 999 999 0 0 1 -360 360;
];
""")
        end
        active_before = active_sparlectra_config()
        results = run_sparlectra_cases(config = cfg_precedence, path = tmpdir)
        @test length(results) == 2
        @test [result.net.name for result in results] == ["case_a", "case_b"]
        @test all(result -> result isa SparlectraRunResult, results)
        @test all(result -> result.numerical_converged, results)
        @test active_sparlectra_config() === active_before
      end
    end

    @testset "Framework runner input validation and removed legacy keywords" begin
      net = createTest3BusNet()
      @test_throws ArgumentError run_sparlectra(net = net, casefile = "case3.m")
      @test_throws ArgumentError run_sparlectra()
      @test_throws ArgumentError run_sparlectra(net = net, path = "unused")
      @test_throws ArgumentError run_acpflow(net = net, casefile = "case3.m")
      @test_throws ArgumentError run_acpflow()
      @test_throws ArgumentError run_acpflow(net = net, path = "unused")
      @test_throws MethodError run_acpflow(net = net, max_ite = 10)
      @test_throws MethodError run_acpflow(net = net, tol = 1e-6)
      @test_throws MethodError run_acpflow(net = net, verbose = 0)
      @test_throws MethodError run_acpflow(casefile = "case14.m", matpower_ratio = :normal)
    end

    @testset "Rectangular damping defaults and validation" begin
      @test Sparlectra.active_sparlectra_config().powerflow.autodamp_min == 0.05

      @test_throws ErrorException Sparlectra._validate_rectangular_damping(1.0, 0.0)
      @test_throws ErrorException Sparlectra._validate_rectangular_damping(0.2, 0.3)
      @test_throws ErrorException Sparlectra._validate_rectangular_damping(0.0, 0.01)
      @test_throws ErrorException Sparlectra._validate_rectangular_damping(1.01, 0.05)
    end

    @testset "Current-iteration rejection log reports candidate guard details" begin
      mktempdir() do tmpdir
        run_guarded_current_iteration_start = getfield(Sparlectra, :_run_guarded_current_iteration_start)
        Ybus = ComplexF64[1 -1; -1 1]
        original = ComplexF64[1.0 + 0.0im, 1.0 + 0.0im]
        S = ComplexF64[0.0 + 0.0im, 10.0 + 0.0im]
        bus_types = [:Slack, :PQ]
        profile = Dict{Symbol, Any}(:output_dir => tmpdir)

        restored, summary = run_guarded_current_iteration_start(
          Ybus,
          original,
          S,
          bus_types,
          [1.0, 1.0],
          1;
          enabled = true,
          max_iter = 1,
          tol = 1.0e-6,
          damping = 1.0,
          accept_only_if_improved = true,
          min_improvement_factor = 0.98,
          vm_min_pu = 0.5,
          vm_max_pu = 1.5,
          max_angle_step_deg = 30.0,
          only_for_large_cases = false,
          large_case_min_buses = 10,
          performance_profile = profile,
        )

        @test restored == original
        @test summary.current_iteration_attempted === true
        @test summary.current_iteration_accepted === false
        @test summary.current_iteration_reason === :voltage_magnitude_guard
        log_text = read(joinpath(tmpdir, "current_iteration_start.log"), String)
        for needle in (
          "current_iteration_attempted: true",
          "current_iteration_accepted: false",
          "current_iteration_reason: voltage_magnitude_guard",
          "candidate_voltage_magnitude_min:",
          "candidate_voltage_magnitude_max: 11.0",
          "candidate_voltage_low_count:",
          "candidate_voltage_high_count: 1",
          "candidate_voltage_worst_high_bus: 2",
          "candidate_voltage_worst_high_value: 11.0",
          "candidate_max_angle_step_deg:",
          "rejected_at_iteration: 1",
          "rejection_stage: candidate_guard",
          "original_start_values_restored: true",
          "restored_voltage_magnitude_min: 1.0",
          "restored_voltage_magnitude_max: 1.0",
        )
          @test occursin(needle, log_text)
        end
      end
    end

    @testset "AC island detection and independent solving" begin
      function two_island_net(; second_ref::Symbol = :slack)
        net = Net(name = "two_island", baseMVA = 100.0)
        for name in ("A1", "A2", "B1", "B2")
          addBus!(net = net, busName = name, vn_kV = 110.0)
        end
        addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.10, b_pu = 0.0, status = 1)
        addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "A1")
        addProsumer!(net = net, busName = "A2", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
        if second_ref === :slack
          addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
        elseif second_ref === :pv
          addProsumer!(net = net, busName = "B1", type = "SYNCHRONOUSMACHINE", p = 12.0, q = 0.0, vm_pu = 1.0, isRegulated = true)
        end
        addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 8.0, q = 2.0)
        refreshBusTypesFromProsumers!(net)
        return net
      end

      net = two_island_net()
      report = Sparlectra.detect_ac_islands(net)
      @test length(report.rows) == 2
      @test all(row -> row.n_branch == 1, report.rows)
      @test_throws ErrorException runpf!(deepcopy(net); config = PowerFlowConfig(max_iter = 40))

      mktempdir() do tmpdir
        solved = deepcopy(net)
        _, erg = runpf!(solved; config = PowerFlowConfig(max_iter = 40, islands_enabled = true), performance_profile = Dict{Symbol,Any}(:output_dir => tmpdir))
        @test erg == 0
        @test isfile(joinpath(tmpdir, "ac_islands.csv"))
        @test !isfile(joinpath(pwd(), "ac_islands.csv"))
        @test count(!isempty, split(read(joinpath(tmpdir, "ac_islands.csv"), String), '\n')) == 3
        @test all(node -> isfinite(node._vm_pu) && isfinite(node._va_deg), solved.nodeVec)
      end

      pv_ref_net = two_island_net(second_ref = :pv)
      pv_report = Sparlectra.detect_ac_islands(pv_ref_net)
      @test pv_report.rows[2].status == "promote_pv_ref"
      _, pv_erg = runpf!(pv_ref_net; config = PowerFlowConfig(max_iter = 40, islands_enabled = true))
      @test pv_erg == 0

      no_ref_net = two_island_net(second_ref = :none)
      no_ref_report = Sparlectra.detect_ac_islands(no_ref_net)
      @test no_ref_report.rows[2].status == "missing_ref"
      @test_throws ErrorException runpf!(no_ref_net; config = PowerFlowConfig(max_iter = 40, islands_enabled = true))
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
      @test occursin("Final Q-limit validation", out)
      @test occursin("Qgen [MVAr]", out) || occursin("Q-limits: no violations", out)
      empty!(net.qmin_pu)
      append!(net.qmin_pu, [-Inf, -0.01, -Inf])
      empty!(net.qmax_pu)
      append!(net.qmax_pu, [Inf, 0.01, Inf])
      net.nodeVec[2]._qƩGen = 5.0
      io = IOBuffer()
      printFinalLimitValidation(net; q_headroom = 0.0, io = io)
      out = String(take!(io))
      @test occursin("Bus │ Qgen [MVAr] │ Qmin [MVAr] │ Qmax [MVAr] │ Violation [MVAr]", out)
      @test occursin("Displayed values are MVAr. Switching comparisons use internal p.u. values.", out)
      io = IOBuffer()
      printFinalLimitValidation(net; q_headroom = 0.20, io = io, converged = false)
      @test occursin("Last-iteration Q-limit diagnostic (NR did not converge; values are not a valid final solution)", String(take!(io)))
    end
  end
end
