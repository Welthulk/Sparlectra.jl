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

      # Regression: direct run_net_acpflow callers must opt in before the narrow-Q guard
      # locks PV buses to PQ during rectangular pre-processing.
      default_net = zero_range_pv_net()
      redirect_stdout(devnull) do
        run_net_acpflow(net = default_net, max_ite = 0, verbose = 0, show_results = false)
      end
      @test isempty(default_net.qLimitLog)

      opt_in_net = zero_range_pv_net()
      redirect_stdout(devnull) do
        run_net_acpflow(net = opt_in_net, max_ite = 0, verbose = 0, show_results = false, qlimit_guard = true)
      end
      @test length(opt_in_net.qLimitLog) == 1
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
