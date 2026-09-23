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

# file: test/test_state_estimation.jl
# purpose: tests WLS state estimation (runse!): measurement add helpers,
#          zero-injection handling, bad-data diagnostics, PMU/current
#          measurements, shunt and link handling, robust R, island-wise SE,
#          measurement CSV, the SE-to-PF chain and tap estimation; the
#          observability side lives in test/test_observability.jl

using Test
using Logging
using Sparlectra
using Random

function test_state_estimation_wls_first_version()::Bool
  # Verifies the baseline WLS state-estimation workflow:
  # PF reference generation, synthetic measurements, convergence, and voltage accuracy.
  @testset "State estimation WLS first version" begin (function ()
    net = createTest3BusNet()

    ite, erg = runpf!(net, 40, 1e-10, 0; method = :rectangular)
    @test erg == 0
    @test ite > 0

    Vref = buildVoltageVector(net)

    std = measurementStdDevs(vm = 1e-5, pinj = 1e-4, qinj = 1e-4, pflow = 1e-4, qflow = 1e-4)
    meas = generateMeasurementsFromPF(net; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false, stddev = std, rng = MersenneTwister(42))

    @test !isempty(meas)
    @test all(m -> m.active, meas)

    result = with_state_estimation_config(max_iter = 12, tol = 1e-8, flatstart = true, jac_eps = 1e-6, update_net = true) do
      runse!(net, meas)
    end

    @test result.converged == true
    @test result.iterations > 0
    @test length(result.voltages) == length(Vref)
    @test result.objectiveJ >= 0.0
    @test result.dof == length(meas) - (2 * length(net.nodeVec) - 1)

    vm_ref = abs.(Vref)
    vm_est = abs.(result.voltages)
    @test maximum(abs.(vm_ref .- vm_est)) < 2e-4

    for i in eachindex(Vref)
      nt = getNodeType(net.nodeVec[i])
      if nt != :Slack
        @test abs(rad2deg(angle(Vref[i])) - rad2deg(angle(result.voltages[i]))) < 1e-2
      end
    end
  end)() end

  return true
end

function test_state_estimation_measurement_add_helpers()::Bool
  # Verifies measurement insertion helper APIs for Vm/Pinj/Qinj/Pflow/Qflow,
  # field mapping, storage in net.measurements, and invalid-branch error handling.
  @testset "State estimation measurement add helpers" begin (function ()
    net = createTest3BusNet()
    @test isempty(net.measurements)

    vm = addVmMeasurement!(net; busName = "ASTADT", value = 1.01, sigma = 0.002)
    pinj = addPinjMeasurement!(net; busName = "STATION1", value = -25.0, sigma = 1.0, id = "PINJ_STATION1")
    qinj = addQinjMeasurement!(net; busName = "STATION1", value = -8.0, sigma = 1.2)
    pflow = addPflowMeasurement!(net; fromBus = "ASTADT", toBus = "STATION1", value = 24.0, sigma = 0.8, direction = :from)
    qflow = addQflowMeasurement!(net; branchNr = 1, value = 6.0, sigma = 0.9, direction = :to, active = false)

    @test length(net.measurements) == 5
    @test vm.typ == Sparlectra.VmMeas
    @test vm.busIdx == Sparlectra.geNetBusIdx(net = net, busName = "ASTADT")
    @test pinj.typ == Sparlectra.PinjMeas
    @test pinj.id == "PINJ_STATION1"
    @test qinj.typ == Sparlectra.QinjMeas
    @test pflow.typ == Sparlectra.PflowMeas
    @test pflow.branchIdx == 1
    @test pflow.direction == :from
    @test qflow.typ == Sparlectra.QflowMeas
    @test qflow.direction == :to
    @test qflow.active == false
    @test net.measurements[end] == qflow

    @test_throws ErrorException addPflowMeasurement!(Measurement[]; net = net, fromBus = "ASTADT", toBus = "UNKNOWN", value = 1.0, sigma = 0.1)
  end)() end

  return true
end

function create_state_estimation_passive_transit_net()::Net
  net = Net(name = "se_passive_transit", baseMVA = 100.0)

  addBus!(net = net, busName = "Slack", vn_kV = 110.0)
  addBus!(net = net, busName = "Transit", vn_kV = 110.0)
  addBus!(net = net, busName = "Load", vn_kV = 110.0)

  addACLine!(net = net, fromBus = "Slack", toBus = "Transit", length = 10.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "Transit", toBus = "Load", length = 8.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)

  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
  # Moderate loading keeps the flat start inside the Gauss-Newton basin of
  # attraction for the undamped WLS iteration. At 40 MW the flat start lies
  # outside the basin and the iteration wanders chaotically before it locks
  # in; the count then depends on BLAS/LAPACK rounding details (45 on Julia
  # 1.12.6, 69 on 1.12.7 locally, >100 on the CI runner). At 10 MW it
  # converges in 5 iterations on every platform tested. The passive-bus and
  # observability properties this fixture exists for do not depend on the
  # loading level.
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = 10.0, q = 4.0)

  return net
end

function test_state_estimation_passive_bus_zero_injection_helpers()::Bool
  # Checks passive-bus detection and automatic zero-injection measurement addition,
  # and confirms improved observability/SE convergence after augmentation.
  @testset "State estimation passive bus zero-injection helpers" begin (function ()
    net = create_state_estimation_passive_transit_net()

    ite, erg = runpf!(net, 40, 1e-10, 0; method = :rectangular)
    @test erg == 0
    @test ite > 0

    passive = findPassiveBuses(net)
    @test passive == [2]

    std = measurementStdDevs(vm = 1e-5, pinj = 1e-5, qinj = 1e-5, pflow = 1e-5, qflow = 1e-5)
    all_meas = setMeasurementsFromPF!(net; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false, stddev = std, rng = MersenneTwister(7))

    base_meas = Measurement[]
    for m in all_meas
      keep = (m.typ == Sparlectra.VmMeas && m.busIdx == 1) ||
             ((m.typ == Sparlectra.PflowMeas || m.typ == Sparlectra.QflowMeas) && m.branchIdx == 1 && m.direction == :from) ||
             ((m.typ == Sparlectra.PflowMeas || m.typ == Sparlectra.QflowMeas) && m.branchIdx == 2 && m.direction == :to)
      keep && push!(base_meas, m)
    end

    @test length(base_meas) == 5
    empty!(net.measurements)
    append!(net.measurements, base_meas)

    base_obs = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
      evaluate_global_observability(net)
    end
    @test base_obs.quality == :critical
    @test base_obs.redundancy == 0
    @test !isempty(base_obs.numerical_critical_measurement_indices)

    added = addZeroInjectionMeasurements!(net; sigma = 1e-6)
    @test length(added) == 2
    @test all(m -> m.busIdx == 2, added)
    @test sort([m.typ for m in added]) == sort([Sparlectra.PinjMeas, Sparlectra.QinjMeas])

    zi_obs = with_state_estimation_config(flatstart = true, jac_eps = 1e-6) do
      evaluate_global_observability(net)
    end
    @test zi_obs.quality == :good
    @test zi_obs.redundancy == 2
    @test isempty(zi_obs.numerical_critical_measurement_indices)
    @test isempty(zi_obs.structural_critical_measurement_indices)

    # 5 iterations measured; budget 20 leaves margin without accepting a
    # chaotic wander (see the fixture comment on the loading level). The
    # solution stays exact (J = 0, noise-free measurement set).
    result = with_state_estimation_config(max_iter = 20, tol = 1e-8, flatstart = true, jac_eps = 1e-6, update_net = false) do
      runse!(deepcopy(net))
    end
    @test result.converged == true
    @test result.objectiveJ < 1e-6
    @test result.residualNorm < 1e-6
  end)() end

  return true
end

function test_state_estimation_bad_data_diagnostics()::Bool
  # Injects a controlled bad measurement and verifies diagnostics:
  # residual ranking, suspect identification, formatted reporting, and rerun improvement.
  @testset "State estimation bad-data diagnostics" begin (function ()
    net = createTest3BusNet()

    ite, erg = runpf!(net, 40, 1e-10, 0; method = :rectangular)
    @test erg == 0
    @test ite > 0

    std = measurementStdDevs(vm = 1e-5, pinj = 1e-4, qinj = 1e-4, pflow = 1e-4, qflow = 1e-4)
    meas = generateMeasurementsFromPF(net; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false, stddev = std, rng = MersenneTwister(42))
    @test !isempty(meas)

    bad_idx = findfirst(m -> m.typ == Sparlectra.VmMeas, meas)
    @test !isnothing(bad_idx)
    idx = something(bad_idx, 0)
    bad_meas = meas[idx]
    meas[idx] = Measurement(
      typ = bad_meas.typ,
      value = bad_meas.value + 0.15,
      sigma = bad_meas.sigma,
      active = bad_meas.active,
      busIdx = bad_meas.busIdx,
      branchIdx = bad_meas.branchIdx,
      direction = bad_meas.direction,
      id = bad_meas.id,
    )

    report = with_state_estimation_config(max_iter = 12, tol = 1e-8, flatstart = true, jac_eps = 1e-6, k_eliminate = 3.0) do
      validate_measurements(net, meas)
    end
    @test report.converged == true
    @test report.objective.dof == report.result.dof
    @test report.global_consistency == false
    @test !isempty(report.measurement_ranking)
    @test report.largest_normalized_residual.measurement_index == idx
    @test any(x -> x.measurement_index == idx, report.suspicious_measurements)

    summary = summarize_se_diagnostics(report)
    @test summary.global_consistency == false
    @test summary.suspicious_count > 0
    @test occursin("Objective J is outside", summary.reason)

    io = IOBuffer()
    print_se_diagnostics(io, report; topN = 5)
    txt = String(take!(io))
    @test occursin("State-estimation diagnostics", txt)
    @test occursin("Global consistency: false", txt)
    @test occursin("Flag", txt)
    @test occursin("BAD", txt)

    io_md = IOBuffer()
    print_se_diagnostics(io_md, report; topN = 5, format = :markdown)
    txt_md = String(take!(io_md))
    @test occursin("## State-estimation diagnostics", txt_md)
    @test occursin("| Idx | ID | Type |", txt_md)
    @test occursin("|", txt_md)

    diag = with_state_estimation_config(max_eliminations = 1, max_iter = 12, tol = 1e-8, flatstart = true, jac_eps = 1e-6, k_eliminate = 3.0) do
      runse_diagnostics(net, meas)
    end
    @test !isnothing(diag.rerun)
    @test diag.rerun.deactivated_measurement_index == idx
    @test diag.rerun.diagnostics.objective.value < diag.diagnostics.objective.value

    io2 = IOBuffer()
    print_se_diagnostics(io2, diag; topN = 5)
    txt2 = String(take!(io2))
    @test occursin("Deactivate-and-rerun", txt2)

    io2_md = IOBuffer()
    print_se_diagnostics(io2_md, diag; topN = 5, format = :markdown)
    txt2_md = String(take!(io2_md))
    @test occursin("### Deactivate-and-rerun", txt2_md)
  end)() end

  return true
end

function test_state_estimation_pmu_va_measurements()::Bool
  # Verifies PMU voltage-angle (VaMeas) support: the common reference-angle
  # offset α is carried as an additional state, estimated network angles stay
  # slack-referenced, and pmuRefOffset=:off treats PMU angles as
  # slack-referenced without the extra state.
  @testset "State estimation PMU angle measurements" begin (function ()
    net = createTest3BusNet()

    ite, erg = runpf!(net, 40, 1e-10, 0; method = :rectangular)
    @test erg == 0
    @test ite > 0
    Vref = buildVoltageVector(net)
    nbus = length(net.nodeVec)

    std = measurementStdDevs(vm = 1e-5, pinj = 1e-4, qinj = 1e-4, pflow = 1e-4, qflow = 1e-4, va = 0.02)
    @test haskey(std, Sparlectra.VaMeas)

    # PMU time base shifted by +5° against the slack reference.
    meas = generateMeasurementsFromPF(net; includeVa = true, vaRefOffsetDeg = 5.0, noise = false, stddev = std)
    vaRows = [m for m in meas if m.typ == Sparlectra.VaMeas]
    @test length(vaRows) == nbus
    @test all(m -> startswith(m.id, "Va_bus_"), vaRows)

    result = with_state_estimation_config(max_iter = 20, tol = 1e-8, flatstart = true, update_net = false) do
      runse!(net, meas)
    end
    @test result.converged
    @test result.vaRefOffsetDeg !== nothing
    @test isapprox(result.vaRefOffsetDeg, 5.0; atol = 1e-6)
    # dof accounts for the extra offset state: n = 2*nbus - 1 + 1
    @test result.dof == length(meas) - 2 * nbus

    # Estimated network angles stay slack-referenced despite the PMU offset.
    for i in eachindex(Vref)
      @test abs(abs(Vref[i]) - abs(result.voltages[i])) < 1e-8
      @test abs(rad2deg(angle(Vref[i]) - angle(result.voltages[i]))) < 1e-6
    end

    # Observability includes the offset column.
    gobs = evaluate_global_observability(net, meas)
    @test gobs.n_states == 2 * nbus
    @test gobs.quality == :good

    # :off mode: no offset state, slack-referenced PMU angles estimate cleanly.
    meas0 = generateMeasurementsFromPF(net; includeVa = true, vaRefOffsetDeg = 0.0, noise = false, stddev = std)
    res0 = with_state_estimation_config(max_iter = 20, tol = 1e-8, update_net = false, pmu_ref_offset = :off) do
      runse!(net, meas0)
    end
    @test res0.converged
    @test res0.vaRefOffsetDeg === nothing
    @test res0.dof == length(meas0) - (2 * nbus - 1)
    for i in eachindex(Vref)
      @test abs(rad2deg(angle(Vref[i]) - angle(res0.voltages[i]))) < 1e-6
    end

    # :off mode with a foreign PMU time base leaves the offset unmodeled and
    # inflates the objective by orders of magnitude.
    resBiased = with_state_estimation_config(max_iter = 20, tol = 1e-8, update_net = false, pmu_ref_offset = :off) do
      runse!(net, meas)
    end
    @test resBiased.objectiveJ > 1e6 * max(result.objectiveJ, eps())

    # Without any Va measurement, no offset state appears even in :auto mode.
    measNoVa = generateMeasurementsFromPF(net; includeVa = false, noise = false, stddev = std)
    resNoVa = with_state_estimation_config(max_iter = 20, tol = 1e-8, update_net = false) do
      runse!(net, measNoVa)
    end
    @test resNoVa.vaRefOffsetDeg === nothing

    # Va measurements restricted to selected buses.
    measSub = generateMeasurementsFromPF(net; includeVa = true, vaBusIdxs = [2], noise = false, stddev = std)
    @test count(m -> m.typ == Sparlectra.VaMeas, measSub) == 1

    # Manual helper resolves bus names and applies the default id.
    va = addVaMeasurement!(net; busName = "ASTADT", value = -2.5, sigma = 0.02)
    @test va.typ == Sparlectra.VaMeas
    @test startswith(va.id, "Va_bus_")
    @test va.busIdx == geNetBusIdx(net = net, busName = "ASTADT")
    @test_throws Exception addMeasurement!(net; typ = Sparlectra.VaMeas, value = 0.0, sigma = 0.02)

    # Combined PMU phasor helper: one call appends the Vm + Va pair with
    # PMU-class default sigmas.
    vmP, vaP = addPmuPhasorMeasurement!(net; busName = "ASTADT", vm_pu = 1.01, va_deg = -2.4)
    @test vmP.typ == Sparlectra.VmMeas
    @test vaP.typ == Sparlectra.VaMeas
    @test vmP.sigma == 0.002
    @test vaP.sigma == 0.02
    @test vmP.busIdx == vaP.busIdx == geNetBusIdx(net = net, busName = "ASTADT")
    @test startswith(vmP.id, "PMU_Vm_bus_")
    @test startswith(vaP.id, "PMU_Va_bus_")

    # Bad-data diagnostics with the offset state present: the estimator
    # converges, and the Wilson-Hilferty band test (SE phase 4) now flags the
    # noise-free synthetic set as :low (J implausibly small, sigmas
    # overestimated) instead of silently passing like the old symmetric test.
    diag = validate_measurements(net, meas)
    @test diag.converged
    @test !diag.global_consistency
    @test diag.objective.reason == :low
  end)() end

  return true
end

function test_state_estimation_imag_measurements()::Bool
  # Verifies the branch current-magnitude measurement (ImagMeas, ampere):
  # generator/prediction roundtrip against an independent current computation,
  # the 3-sigma value gate, the busIdx-only rejection (phase 2 door), and
  # that observability is invariant under ImagMeas (currents never carry it).
  @testset "State estimation current-magnitude measurements (ImagMeas)" begin (function ()
    net = createTest3BusNet()
    ite, erg = runpf!(net, 40, 1e-10, 0; method = :rectangular)
    @test erg == 0
    @test ite > 0
    calcNetLosses!(net)

    # 1) roundtrip: generated value == prediction, and both match the current
    # computed independently from the calcNetLosses! branch flows
    std = measurementStdDevs(vm = 1e-5, pinj = 1e-4, qinj = 1e-4, pflow = 1e-4, qflow = 1e-4, imag = 1e-3)
    meas = generateMeasurementsFromPF(net; includeImag = true, noise = false, stddev = std)
    imags = [m for m in meas if m.typ == Sparlectra.ImagMeas]
    @test length(imags) == 2 * length(net.branchVec)
    V = buildVoltageVector(net)
    Sbus_MVA = Sparlectra.calc_injections(createYBUS(net = net), V) .* net.baseMVA
    for m in imags
      pred = Sparlectra._measurement_prediction(m, net, V, Sbus_MVA)
      @test isapprox(pred, m.value; rtol = 1e-8)
      # independent reference from the stored branch-end flows
      br = net.branchVec[m.branchIdx]
      bf = m.direction == :from ? br.fBranchFlow : br.tBranchFlow
      endBus = m.direction == :from ? Int(br.fromBus) : Int(br.toBus)
      i_ref = 1000.0 * sqrt(bf.pFlow^2 + bf.qFlow^2) / (sqrt(3.0) * getNodeVn(net.nodeVec[endBus]) * abs(V[endBus]))
      @test isapprox(pred, i_ref; rtol = 1e-6)
    end

    # SE with currents converges to the same state as without
    resWith = with_state_estimation_config(max_iter = 20, tol = 1e-8, update_net = false) do
      runse!(net, meas)
    end
    @test resWith.converged

    # 2) value gate: a current below 3 sigma never enters the active set
    weak = copy(meas)
    addImagMeasurement!(weak; net = net, value = 1.0, sigma = 10.0, fromBus = "ASTADT", toBus = "STATION1", id = "Imag_weak")
    gated = Test.@test_logs (:info, r"Imag_weak excluded") match_mode = :any with_state_estimation_config(() -> runse!(net, weak); max_iter = 20, tol = 1e-8, update_net = false)
    @test gated.converged
    report = with_state_estimation_config(max_iter = 20, tol = 1e-8) do
      validate_measurements(net, weak)
    end
    @test length(report.measurement_ranking) == length(meas)   # weak row filtered out
    @test all(row.id != "Imag_weak" for row in report.measurement_ranking)

    # 3) busIdx-only ImagMeas is the shunt-bay variant (SE phase 2): valid
    # without a direction, rejected with one (the phase-1 not-yet-supported
    # error was replaced by the implementation, by design)
    mshunt = addMeasurement!(copy(meas); typ = Sparlectra.ImagMeas, value = 100.0, sigma = 1.0, busIdx = 1)
    @test mshunt.direction == :none
    @test startswith(mshunt.id, "Imag_shunt_bus_")
    err = try
      addMeasurement!(copy(meas); typ = Sparlectra.ImagMeas, value = 100.0, sigma = 1.0, busIdx = 1, direction = :from)
      nothing
    catch e
      e
    end
    @test err !== nothing
    @test occursin("takes no direction", sprint(showerror, err))

    # 4) observability invariance: identical verdict with and without currents
    noImag = [m for m in meas if m.typ != Sparlectra.ImagMeas]
    obsWith = evaluate_global_observability(net, meas)
    obsWithout = evaluate_global_observability(net, noImag)
    @test obsWith.numerical_rank == obsWithout.numerical_rank
    @test obsWith.quality == obsWithout.quality
    @test obsWith.n_measurements == obsWithout.n_measurements
  end)() end

  return true
end

function test_state_estimation_sequential_elimination()::Bool
  # Verifies bad-data localization stage 1: wii exposure, the localization
  # benefit of current measurements, sequential elimination with trace and
  # stop reasons, and the ZIB protection.
  @testset "State estimation sequential elimination and wii" begin (function ()
    std = measurementStdDevs(vm = 0.001, pinj = 0.05, qinj = 0.05, pflow = 0.05, qflow = 0.05, imag = 2.0)

    _fresh_meas(net; withImag::Bool) = generateMeasurementsFromPF(net; includeImag = withImag, noise = true, stddev = std, rng = MersenneTwister(7))

    _inject!(meas, idx, nsigma) = begin
      bm = meas[idx]
      meas[idx] = Measurement(typ = bm.typ, value = bm.value + nsigma * bm.sigma, sigma = bm.sigma, active = bm.active, busIdx = bm.busIdx, branchIdx = bm.branchIdx, direction = bm.direction, id = bm.id)
      bm.id
    end

    # 1) localization benefit: one 10-sigma gross error on a Pflow measurement,
    # identified first both without and with current measurements; the wii of
    # the faulted measurement must not decrease when currents are added
    net = createTest3BusNet()
    ite, erg = runpf!(net, 40, 1e-10, 0; method = :rectangular)
    @test erg == 0
    @test ite > 0

    wiiByCase = Float64[]
    for withImag in (false, true)
      meas = _fresh_meas(net; withImag = withImag)
      badIdx = findfirst(m -> m.typ == Sparlectra.PflowMeas, meas)
      badId = _inject!(meas, badIdx, 10.0)
      diag = with_state_estimation_config(max_eliminations = 3, max_iter = 20, tol = 1e-8) do
        runse_diagnostics(net, meas)
      end
      @test !isempty(diag.eliminations)
      @test diag.eliminations[1].id == badId
      row = only(r for r in diag.diagnostics.measurement_ranking if r.measurement_index == badIdx)
      @test row.wii >= 0.0
      @test row.localizable == (row.wii > 0.3)
      push!(wiiByCase, row.wii)
    end
    @test wiiByCase[2] >= wiiByCase[1] - 1e-6

    # 2) two independent gross errors, both identified, stop :consistent
    meas = _fresh_meas(net; withImag = false)
    idxA = findfirst(m -> m.typ == Sparlectra.PflowMeas, meas)
    idxB = findfirst(m -> m.typ == Sparlectra.VmMeas && m.busIdx == 1, meas)
    idA = _inject!(meas, idxA, 10.0)
    idB = _inject!(meas, idxB, 10.0)
    diag = with_state_estimation_config(max_eliminations = 3, max_iter = 20, tol = 1e-8) do
      runse_diagnostics(net, meas)
    end
    @test diag.stop_reason == :consistent
    @test length(diag.eliminations) == 2
    @test Set(t.id for t in diag.eliminations) == Set([idA, idB])
    @test diag.final_diagnostics.objective.within_3sigma
    # trace rows carry the objective drop
    @test all(t.objective_after < t.objective_before for t in diag.eliminations)

    # state_estimation.max_eliminations = 1 does exactly one pass
    diag1 = with_state_estimation_config(max_eliminations = 1, max_iter = 20, tol = 1e-8) do
      runse_diagnostics(net, meas)
    end
    @test length(diag1.eliminations) == 1
    @test diag1.rerun !== nothing
    @test diag1.rerun.deactivated_measurement_index == diag1.eliminations[1].measurement_index

    # printed report carries the wii column and the elimination trace
    io = IOBuffer()
    print_se_diagnostics(io, diag; topN = 5)
    txt = String(take!(io))
    @test occursin("wii", txt)
    @test occursin("Sequential elimination", txt)
    @test occursin("consistent", txt)

    # optional K-matrix report: MaxK column and warning flag fields
    reportK = with_state_estimation_config(max_iter = 20, tol = 1e-8, report_residual_correlation = true) do
      validate_measurements(net, meas)
    end
    @test reportK.correlation_enabled
    @test all(!isnan(row.max_abs_correlation) for row in reportK.measurement_ranking)
    ioK = IOBuffer()
    print_se_diagnostics(ioK, reportK; topN = 5)
    @test occursin("MaxK", String(take!(ioK)))

    # 3) ZIB protection: gross error on a regular measurement with ZIB
    # pseudo-measurements present; no ZI id is ever eliminated
    tnet = create_state_estimation_passive_transit_net()
    ite2, erg2 = runpf!(tnet, 40, 1e-10, 0; method = :rectangular)
    @test erg2 == 0
    @test ite2 > 0
    tmeas = generateMeasurementsFromPF(tnet; includePinj = false, includeQinj = false, noise = true, stddev = std, rng = MersenneTwister(11))
    addZeroInjectionMeasurements!(tmeas; net = tnet)
    @test any(m -> startswith(m.id, "ZI"), tmeas)
    tbad = findfirst(m -> m.typ == Sparlectra.PflowMeas, tmeas)
    _inject!(tmeas, tbad, 10.0)
    tdiag = with_state_estimation_config(max_eliminations = 3, max_iter = 30, tol = 1e-8) do
      runse_diagnostics(tnet, tmeas)
    end
    @test all(!startswith(t.id, "ZI") for t in tdiag.eliminations)
  end)() end

  return true
end

function create_se_shunt_net(; qShunt::Float64 = 30.0, shunt_model = nothing)::Net
  # meshed 3-bus net with one shunt at LoadB (reactor-style for qShunt > 0)
  net = Net(name = "se_shunt", baseMVA = 100.0)
  for b in ("Slack", "LoadA", "LoadB")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addACLine!(net = net, fromBus = "Slack", toBus = "LoadA", length = 12.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "LoadA", toBus = "LoadB", length = 9.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addACLine!(net = net, fromBus = "Slack", toBus = "LoadB", length = 11.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
  addProsumer!(net = net, busName = "Slack", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack")
  addProsumer!(net = net, busName = "LoadA", type = "ENERGYCONSUMER", p = 35.0, q = 12.0)
  addProsumer!(net = net, busName = "LoadB", type = "ENERGYCONSUMER", p = 28.0, q = 9.0)
  if shunt_model === nothing
    addShunt!(net = net, busName = "LoadB", pShunt = 0.0, qShunt = qShunt)
  else
    addShunt!(net = net, busName = "LoadB", pShunt = 0.0, qShunt = qShunt, bus_shunt_model = shunt_model)
  end
  ok, msg = validate!(net = net)
  ok || error("se_shunt net invalid: $msg")
  return net
end

function test_state_estimation_shunt_estimation()::Bool
  # SE phase 2 case A: released shunt susceptance as an estimator state
  # (recovery, freeze guards, write-back gate, voltage sensitivity) plus the
  # bus-referenced ImagMeas prediction and the injection-mode rejection.
  _se_shunt_std() = measurementStdDevs(vm = 1e-4, pinj = 1e-3, qinj = 1e-3, pflow = 1e-3, qflow = 1e-3, shuntq = 1e-3, imag = 1e-2)

  @testset "State estimation shunt estimation (case A)" begin (function ()
    # --- 1) roundtrip: model B off by 20 percent, direct ShuntQ + Vm recover it
    net = create_se_shunt_net()
    ite, erg = runpf!(net, 40, 1e-10, 0; method = :rectangular)
    @test erg == 0
    @test ite > 0
    meas = generateMeasurementsFromPF(net; includeShuntQ = true, noise = false, stddev = _se_shunt_std())
    sh = net.shuntVec[1]
    btrue = imag(sh.y_pu_shunt)
    sh.y_pu_shunt = complex(real(sh.y_pu_shunt), 1.2 * btrue)   # perturbed model
    y_before = sh.y_pu_shunt                                    # snapshot before the run
    setShuntEstimation!(net; busName = "LoadB")
    res = with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = false) do
      runse!(net, meas)
    end
    @test res.converged
    @test res.shuntEstimates !== nothing
    row = only(res.shuntEstimates)
    @test row.busName == "LoadB"
    @test !row.frozen
    @test isapprox(row.B_model, 1.2 * btrue; rtol = 1e-12)
    @test abs(row.B_est - btrue) / abs(btrue) < 1e-3
    @test isapprox(row.delta, row.B_est - row.B_model; atol = 1e-12)
    # estimation never silently overwrites the model: the full complex model
    # admittance is bitwise unchanged after the updateShunts = false run
    @test sh.y_pu_shunt == y_before
    @test imag(sh.y_pu_shunt) == 1.2 * btrue
    # write-back only behind updateShunts
    with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = false, update_shunts = true) do
      runse!(net, meas)
    end
    @test abs(imag(sh.y_pu_shunt) - btrue) / abs(btrue) < 1e-3
    @test isapprox(sh.B_shunt, imag(sh.y_pu_shunt); rtol = 1e-12)
    dev_with_vm = abs(row.B_est - btrue)

    # --- 2) freeze guard: released without any direct shunt measurement
    net2 = create_se_shunt_net()
    runpf!(net2, 40, 1e-10, 0; method = :rectangular)
    meas2 = generateMeasurementsFromPF(net2; includeShuntQ = false, noise = false, stddev = _se_shunt_std())
    setShuntEstimation!(net2; busName = "LoadB")
    res2 = Test.@test_logs (:warn, r"no active direct shunt measurement") match_mode = :any with_state_estimation_config(() -> runse!(net2, meas2); max_iter = 30, tol = 1e-10, update_net = false)
    @test res2.converged
    frow = only(res2.shuntEstimates)
    @test frow.frozen
    @test frow.B_est == frow.B_model
    # identical estimate to a run without the release
    net2b = create_se_shunt_net()
    runpf!(net2b, 40, 1e-10, 0; method = :rectangular)
    res2b = with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = false) do
      runse!(net2b, meas2)
    end
    @test maximum(abs.(res2.voltages .- res2b.voltages)) <= 1e-12

    # --- 3) voltage sensitivity: no Vm measurement at the shunt bus degrades B
    #        (the concept's Vm recommendation; V enters Q = -B V^2 squared).
    #        The local Vm measurement is the only voltage anchor here, so its
    #        absence must visibly degrade B; averaged over seeds because a
    #        single noise draw is not meaningful.
    meanDev = Dict{Bool,Float64}()
    for dropVm in (false, true)
      acc = 0.0
      nseeds = 10
      for seed = 1:nseeds
        netv = create_se_shunt_net()
        runpf!(netv, 40, 1e-10, 0; method = :rectangular)
        Vv = buildVoltageVector(netv)
        stdn = measurementStdDevs(pinj = 0.5, qinj = 0.5, pflow = 0.5, qflow = 0.5, shuntq = 0.5)
        rngv = MersenneTwister(1000 + seed)
        mv = generateMeasurementsFromPF(netv; includeVm = false, includeShuntQ = true, noise = true, stddev = stdn, rng = rngv)
        shuntBus = netv.shuntVec[1].busIdx
        if !dropVm
          # the recommended companion measurement: Vm at the shunt bus
          σvm = 0.002
          addVmMeasurement!(mv; net = netv, busName = "LoadB", value = abs(Vv[shuntBus]) + randn(rngv) * σvm, sigma = σvm)
        end
        shv = netv.shuntVec[1]
        btruev = imag(shv.y_pu_shunt)
        shv.y_pu_shunt = complex(real(shv.y_pu_shunt), 1.2 * btruev)
        setShuntEstimation!(netv; busName = "LoadB")
        resv = with_state_estimation_config(max_iter = 60, tol = 1e-8, update_net = false) do
          runse!(netv, mv)
        end
        @test resv.converged
        acc += abs(only(resv.shuntEstimates).B_est - btruev)
      end
      meanDev[dropVm] = acc / nseeds
    end
    @test meanDev[true] > meanDev[false]

    # --- 4) bus-referenced ImagMeas: prediction matches the PF shunt current
    net4 = create_se_shunt_net()
    runpf!(net4, 40, 1e-10, 0; method = :rectangular)
    updateShuntPowers!(net = net4)
    sh4 = net4.shuntVec[1]
    V4 = buildVoltageVector(net4)
    i4 = sh4.busIdx
    # independent reference from the post-solve shunt powers
    s_ref = sqrt(sh4.p_shunt^2 + sh4.q_shunt^2)
    i_ref = 1000.0 * s_ref / (sqrt(3.0) * getNodeVn(net4.nodeVec[i4]) * abs(V4[i4]))
    m4 = addImagMeasurement!(Measurement[]; net = net4, value = i_ref, sigma = 1.0, busName = "LoadB")
    Sbus4 = Sparlectra.calc_injections(createYBUS(net = net4), V4) .* net4.baseMVA
    pred = Sparlectra._measurement_prediction(m4, net4, V4, Sbus4)
    @test isapprox(pred, i_ref; rtol = 1e-8)
    # shunt-bay variant rejects a direction and a shuntless bus
    @test_throws ErrorException addImagMeasurement!(Measurement[]; net = net4, value = 1.0, sigma = 1.0, busName = "LoadA")
    @test_throws ErrorException addImagMeasurement!(Measurement[]; net = net4, value = 1.0, sigma = 1.0, busName = "LoadB", branchNr = 1)

    # --- 5) injection-mode shunts are rejected for estimation
    net5 = create_se_shunt_net(shunt_model = "voltage_dependent_injection")
    runpf!(net5, 40, 1e-10, 0; method = :rectangular)
    meas5 = generateMeasurementsFromPF(net5; includeShuntQ = true, noise = false, stddev = _se_shunt_std())
    setShuntEstimation!(net5; busName = "LoadB")
    err5 = try
      with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = false) do
        runse!(net5, meas5)
      end
      nothing
    catch e
      e
    end
    @test err5 !== nothing
    @test occursin("voltage-dependent injection", sprint(showerror, err5))
  end)() end

  return true
end

function test_state_estimation_shunt_back_calculation()::Bool
  # SE phase 2 case B: deriveShuntPseudoMeasurements! (bay current plus Vm to
  # a SHDERIV ShuntQ pseudo-measurement), sign for reactor and capacitor,
  # the hard Vm prerequisite, and the SHDERIV elimination protection.
  @testset "State estimation shunt back-calculation (case B)" begin (function ()
    for (qShunt, kind) in ((30.0, "reactor"), (-25.0, "capacitor"))
      net = create_se_shunt_net(qShunt = qShunt)
      runpf!(net, 40, 1e-10, 0; method = :rectangular)
      updateShuntPowers!(net = net)
      sh = net.shuntVec[1]
      V = buildVoltageVector(net)
      i = sh.busIdx
      # synthetic bay telemetry from the solved state
      i_A = 1000.0 * sqrt(sh.p_shunt^2 + sh.q_shunt^2) / (sqrt(3.0) * getNodeVn(net.nodeVec[i]) * abs(V[i]))
      σ_i = 2.0
      σ_vm = 0.002
      addImagMeasurement!(net; value = i_A, sigma = σ_i, busName = "LoadB")
      addVmMeasurement!(net; busName = "LoadB", value = abs(V[i]), sigma = σ_vm)
      added = deriveShuntPseudoMeasurements!(net)
      @test length(added) == 1
      d = added[1]
      @test startswith(d.id, "SHDERIV")
      @test d.typ == Sparlectra.ShuntQMeas
      # sign matches the q_shunt convention for both device kinds
      @test sign(d.value) == sign(sh.q_shunt)
      # value matches the solved q_shunt within the propagated sigma
      @test abs(d.value - sh.q_shunt) <= d.sigma
      @test d.sigma > 0.0
      kind == "reactor" && @test d.value < 0.0
      kind == "capacitor" && @test d.value > 0.0
    end

    # a direct ShuntQ measurement takes precedence: nothing derived
    netd = create_se_shunt_net()
    runpf!(netd, 40, 1e-10, 0; method = :rectangular)
    updateShuntPowers!(net = netd)
    addShuntQMeasurement!(netd; busName = "LoadB", value = netd.shuntVec[1].q_shunt, sigma = 0.5)
    addImagMeasurement!(netd; value = 100.0, sigma = 2.0, busName = "LoadB")
    addVmMeasurement!(netd; busName = "LoadB", value = 1.0, sigma = 0.002)
    @test isempty(deriveShuntPseudoMeasurements!(netd))

    # missing Vm at the bus: hard prerequisite, skip with warning
    netm = create_se_shunt_net()
    runpf!(netm, 40, 1e-10, 0; method = :rectangular)
    addImagMeasurement!(netm; value = 100.0, sigma = 2.0, busName = "LoadB")
    addedm = Test.@test_logs (:warn, r"no active Vm measurement") match_mode = :any deriveShuntPseudoMeasurements!(netm)
    @test isempty(addedm)

    # SHDERIV pseudo-measurements are protected from sequential elimination
    nete = create_se_shunt_net()
    runpf!(nete, 40, 1e-10, 0; method = :rectangular)
    updateShuntPowers!(net = nete)
    stdn = measurementStdDevs(vm = 0.002, pinj = 0.5, qinj = 0.5, pflow = 0.5, qflow = 0.5)
    mease = generateMeasurementsFromPF(nete; noise = true, stddev = stdn, rng = MersenneTwister(17))
    append!(nete.measurements, mease)
    Ve = buildVoltageVector(nete)
    she = nete.shuntVec[1]
    ie = she.busIdx
    iA = 1000.0 * sqrt(she.p_shunt^2 + she.q_shunt^2) / (sqrt(3.0) * getNodeVn(nete.nodeVec[ie]) * abs(Ve[ie]))
    addImagMeasurement!(nete; value = iA, sigma = 2.0, busName = "LoadB")
    derived = deriveShuntPseudoMeasurements!(nete)
    @test length(derived) == 1
    # corrupt the derived value grossly; the elimination must not remove it
    didx = findfirst(m -> startswith(m.id, "SHDERIV"), nete.measurements)
    dm = nete.measurements[didx]
    nete.measurements[didx] = Measurement(typ = dm.typ, value = dm.value + 20.0 * dm.sigma, sigma = dm.sigma, active = dm.active, busIdx = dm.busIdx, branchIdx = dm.branchIdx, direction = dm.direction, id = dm.id)
    diag = with_state_estimation_config(max_eliminations = 3, max_iter = 30, tol = 1e-8) do
      runse_diagnostics(nete)
    end
    @test all(!startswith(t.id, "SHDERIV") for t in diag.eliminations)
  end)() end

  return true
end

function test_state_estimation_link_contraction()::Bool
  # Link contraction: the estimator runs on the contracted net. Closed link
  # clusters fuse onto their representative (member voltages identical after
  # write-back), open links stay real separations (observability reports the
  # unobservable second island), cluster injections aggregate only as a
  # whole, and link measurements never enter the WLS.
  _lstd() = measurementStdDevs(vm = 1e-4, pinj = 1e-3, qinj = 1e-3, pflow = 1e-3, qflow = 1e-3, shuntq = 1e-3)

  @testset "State estimation link contraction" begin (function ()
    # --- 1) closed link: SE converges on the fused net, members share V
    net = create_se_link_net()
    ite, erg = runpf!(net, 40, 1e-10, 0; method = :rectangular)
    @test erg == 0
    @test ite > 0
    calcNetLosses!(net)
    vref = [n._vm_pu for n in net.nodeVec]
    meas = generateMeasurementsFromPF(net; noise = false, stddev = _lstd())
    res = with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = true) do
      runse!(net, meas)
    end
    @test res.converged
    b1 = geNetBusIdx(net = net, busName = "B1")
    b1a = geNetBusIdx(net = net, busName = "B1a")
    @test net.nodeVec[b1]._vm_pu == net.nodeVec[b1a]._vm_pu
    @test net.nodeVec[b1]._va_deg == net.nodeVec[b1a]._va_deg
    @test abs(net.nodeVec[b1]._vm_pu - vref[b1]) < 1e-6

    # the aggregated cluster injection is reported, not the swallowed members
    report = with_state_estimation_config(max_iter = 30, tol = 1e-8) do
      validate_measurements(net, meas)
    end
    aggRows = [r for r in report.measurement_ranking if startswith(r.id, "LINKAGG")]
    @test length(aggRows) == 2   # Pinj and Qinj aggregate of the B1/B1a cluster
    @test all(r.measurement_index == 0 for r in aggRows)

    # --- 2) aggregation values: sum and root-sum-square sigma
    prep = Sparlectra._se_prepare(net, meas)
    @test prep.has_merges
    aggP = only(m for m in prep.meas if startswith(m.id, "LINKAGG_Pinj"))
    memberP = [m for m in meas if m.typ == Sparlectra.PinjMeas && m.busIdx in (b1, b1a)]
    @test isapprox(aggP.value, sum(m.value for m in memberP); atol = 1e-12)
    @test isapprox(aggP.sigma, sqrt(sum(m.sigma^2 for m in memberP)); atol = 1e-15)
    @test any(n -> n.reason == :aggregated_cluster_injection, prep.notes)

    # partial cluster: only one member measured -> excluded with warning,
    # WLS identical to a run without that measurement
    measPartial = [m for m in meas if !(m.typ == Sparlectra.PinjMeas && m.busIdx == b1a)]
    resPartial = Test.@test_logs (:warn, r"not every cluster member is measured") match_mode = :any with_state_estimation_config(() -> runse!(net, measPartial); max_iter = 30, tol = 1e-10, update_net = false)
    measNone = [m for m in meas if !(m.typ == Sparlectra.PinjMeas && m.busIdx in (b1, b1a))]
    resNone = with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = false) do
      runse!(net, measNone)
    end
    @test resPartial.converged && resNone.converged
    @test maximum(abs.(resPartial.voltages .- resNone.voltages)) <= 1e-12
    @test resPartial.objectiveJ == resNone.objectiveJ

    # --- 3) link measurement never enters the WLS: J and state identical
    # (fresh reference run: the updateNet write-back above changed the flat
    # start's slack seed by the estimation residual, so `res` is not bitwise
    # comparable any more)
    resRef = with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = false) do
      runse!(net, meas)
    end
    measLink = copy(meas)
    addPflowMeasurement!(measLink; net = net, value = 99.0, sigma = 0.01, linkNr = 1)
    resLink = with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = false) do
      runse!(net, measLink)
    end
    @test resLink.objectiveJ == resRef.objectiveJ
    @test maximum(abs.(resLink.voltages .- resRef.voltages)) <= 1e-12

    # --- 4) open link: a real separation, never silently merged. The open
    # link leaves the representative map at identity (no contraction) and the
    # island detection reports the two synchronous islands.
    netO = create_se_link_net(link_status = 0)
    repsO = Sparlectra._active_link_representative_map(netO)
    @test all(repsO[i] == i for i in eachindex(repsO))
    prepO = Sparlectra._se_prepare(netO, [m for m in meas if m.linkIdx === nothing])
    @test !prepO.has_merges
    islands = Sparlectra.detect_ac_islands(netO)
    @test length(islands.rows) == 2

    # --- 5) released shunt on a cluster member resolves to the representative
    netS = create_se_link_net(shunt_buses = ["B1a"])
    runpf!(netS, 40, 1e-10, 0; method = :rectangular)
    measS = generateMeasurementsFromPF(netS; includeShuntQ = true, noise = false, stddev = _lstd())
    setShuntEstimation!(netS; busName = "B1a")
    resS = with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = false) do
      runse!(netS, measS)
    end
    @test resS.converged
    rowS = only(resS.shuntEstimates)
    @test rowS.busIdx == min(geNetBusIdx(net = netS, busName = "B1"), geNetBusIdx(net = netS, busName = "B1a"))
    @test !rowS.frozen

    # two released shunts in one cluster: one B state per fused bus
    netS2 = create_se_link_net(shunt_buses = ["B1", "B1a"])
    runpf!(netS2, 40, 1e-10, 0; method = :rectangular)
    measS2 = generateMeasurementsFromPF(netS2; includeShuntQ = true, noise = false, stddev = _lstd())
    setShuntEstimation!(netS2; busName = "B1")
    setShuntEstimation!(netS2; busName = "B1a")
    errS2 = try
      with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = false) do
        runse!(netS2, measS2)
      end
      nothing
    catch e
      e
    end
    @test errS2 !== nothing
    @test occursin("same (fused) bus", sprint(showerror, errS2))
  end)() end

  return true
end

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
    addACLine!(net = net, fromBus = "S", toBus = "R1", length = 10.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
    addACLine!(net = net, fromBus = "R2", toBus = "L2", length = 8.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
    addACLine!(net = net, fromBus = "R3", toBus = "L3", length = 8.0, r = 0.02, x = 0.2, c_nf_per_km = 0.0, tanδ = 0.0)
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

function test_state_estimation_se_view()::Bool
  # The static SE view classifies controllers (frozen),
  # shunt releases, link clusters, and excluded measurements. The controller
  # row uses a series-reactance controller (the view is type-agnostic via
  # collect_outer_controllers).
  @testset "State estimation se_view" begin (function ()
    net = create_se_link_net(shunt_buses = ["B2"])
    runpf!(net, 40, 1e-10, 0; method = :rectangular)
    meas = generateMeasurementsFromPF(net; includeShuntQ = true, noise = false)
    append!(net.measurements, meas)
    # post-PF additions: se_view itself never solves, so the injection-mode
    # shunt (which the merged PF path rejects) can be registered here
    addLink!(net = net, fromBus = "B1", toBus = "B2", status = 0)
    addShunt!(net = net, busName = "B1", pShunt = 0.0, qShunt = 5.0, bus_shunt_model = "voltage_dependent_injection")
    addSeriesReactanceControl!(net; fromBus = "S", toBus = "B1", p_target_mw = 5.0, x_min_pu = 0.05, x_max_pu = 0.4)
    addPflowMeasurement!(net; value = 1.0, sigma = 0.5, linkNr = 1)
    setShuntEstimation!(net; busName = "B2")
    setShuntEstimation!(net; busName = "B1")

    view = se_view(net)
    @test length(view.controllers) == 1
    @test view.controllers[1].enabled
    @test length(view.shunts) == 2
    statuses = Dict(s.busName => s.status for s in view.shunts)
    @test statuses["B2"] == :released
    @test statuses["B1"] == :injection_mode_rejected
    @test length(view.links.closed_clusters) == 1
    @test length(view.links.open_links) == 1
    cl = view.links.closed_clusters[1]
    @test sort(cl.members) == sort([geNetBusIdx(net = net, busName = "B1"), geNetBusIdx(net = net, busName = "B1a")])
    @test any(n -> n.reason == :link_measurement, view.excluded_measurements)

    io = IOBuffer()
    print_se_view(io, view)
    txt = String(take!(io))
    @test occursin("SE view", txt)
    @test occursin("closed cluster", txt)
    io2 = IOBuffer()
    print_se_view(io2, view; format = :markdown)
    @test occursin("## SE view", String(take!(io2)))
  end)() end

  return true
end

function test_state_estimation_robust_r()::Bool
  # Two-stage robust R modification. The gate is the
  # solve/diagnosis separation: robust changes the estimate, never the
  # original-sigma statistics.
  @testset "State estimation robust R modification" begin (function ()
    # --- stage classification (unit level): t = 2, 4, 8 -> stages 0/1/2
    σ = 0.5
    s0 = Sparlectra._robust_stage(2.0 * σ, σ)
    s1 = Sparlectra._robust_stage(4.0 * σ, σ)
    s2 = Sparlectra._robust_stage(8.0 * σ, σ)
    @test s0 == (0, σ)
    @test s1[1] == 1
    @test isapprox(s1[2], σ * (2.0 * 4.0 / 3.0 - 1.0); rtol = 1e-12)
    @test s2[1] == 2
    @test isapprox(s2[2], 8.0 * σ / 3.0; rtol = 1e-12)
    # virtual measurements are exempt regardless of t
    @test Sparlectra._robust_stage(100.0, 1e-7) == (0, 1e-7)

    # --- integration: one 10 sigma gross error on a meshed case
    std = measurementStdDevs(vm = 0.002, pinj = 0.5, qinj = 0.5, pflow = 0.5, qflow = 0.5)
    net = createTest3BusNet()
    ite, erg = runpf!(net, 40, 1e-10, 0; method = :rectangular)
    @test erg == 0
    @test ite > 0
    Vref = buildVoltageVector(net)

    _meas() = generateMeasurementsFromPF(net; noise = true, stddev = std, rng = MersenneTwister(42))
    meas = _meas()
    badIdx = findfirst(m -> m.typ == Sparlectra.PflowMeas, meas)
    bm = meas[badIdx]
    meas[badIdx] = Measurement(typ = bm.typ, value = bm.value + 10.0 * bm.sigma, sigma = bm.sigma, active = bm.active, busIdx = bm.busIdx, branchIdx = bm.branchIdx, direction = bm.direction, id = bm.id)

    resPlain = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_net = false, robust = false) do
      runse!(net, meas)
    end
    resRobust = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_net = false, robust = true) do
      runse!(net, meas)
    end
    # noise-only reference: the same measurement set without the gross error
    resClean0 = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_net = false, robust = false) do
      runse!(net, _meas())
    end
    @test resPlain.converged && resRobust.converged && resClean0.converged
    devPlain = maximum(abs.(resPlain.voltages .- Vref))
    devRobust = maximum(abs.(resRobust.voltages .- Vref))
    devClean = maximum(abs.(resClean0.voltages .- Vref))
    # robust suppresses the gross error: closer to the true state, and back
    # to the noise floor of the error-free run (within a factor of 2)
    @test devRobust < devPlain
    @test devRobust < 2.0 * devClean
    # the faulted measurement reaches stage 2 in the final iteration
    @test resPlain.robustRows === nothing
    @test resRobust.robustRows !== nothing
    rr = only(r for r in resRobust.robustRows if r.id == bm.id)
    @test rr.stage == 2
    @test rr.t > 6.0
    @test rr.measurement_index == badIdx
    @test rr.sigma_factor > 1.0

    # --- recovery: a consistent set has large flat-start residuals but every
    # measurement returns to stage 0 by the final iteration
    resClean = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_net = false, robust = true) do
      runse!(net, _meas())
    end
    @test resClean.converged
    @test isempty(resClean.robustRows)

    # --- statistics separation (the gate): with robust active, the
    # normalized residual of the faulted measurement still exceeds the
    # threshold (original sigmas) and the ranking names it first
    report = with_state_estimation_config(max_iter = 40, tol = 1e-8, robust = true) do
      validate_measurements(net, meas)
    end
    @test report.robust_rows !== nothing && !isempty(report.robust_rows)
    top = first(report.measurement_ranking)
    @test top.measurement_index == badIdx
    @test top.abs_normalized_residual >= 3.0
    @test top.suspicious

    # --- generalized knees (0.10.0): defaults reproduce the fixed 3/6
    # formulas bitwise, t = k1 is still stage 0, just above k1 the modified
    # sigma is continuous (sigma_mod -> sigma), and the constant stage
    # divides by k1
    @test Sparlectra._robust_stage(4.0 * σ, σ, 3.0, 6.0) === Sparlectra._robust_stage(4.0 * σ, σ)
    k1 = 2.5
    k2 = 5.0
    @test Sparlectra._robust_stage(k1 * σ, σ, k1, k2) == (0, σ)
    scont = Sparlectra._robust_stage((k1 + 1e-9) * σ, σ, k1, k2)
    @test scont[1] == 1
    @test isapprox(scont[2], σ; atol = 1e-6)
    shigh = Sparlectra._robust_stage((k2 + 1.0) * σ, σ, k1, k2)
    @test shigh[1] == 2
    @test isapprox(shigh[2], (k2 + 1.0) * σ / k1; rtol = 1e-12)

    # --- robust_mode selection (0.10.0): staged with default knees is
    # BITWISE the legacy Bool robust = true; an explicit robust_mode wins
    # over the Bool
    resStaged = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_net = false, robust_mode = :staged, robust_k1 = 3.0, robust_k2 = 6.0) do
      runse!(net, meas)
    end
    @test resStaged.voltages == resRobust.voltages
    @test resStaged.objectiveJ == resRobust.objectiveJ
    @test Sparlectra._se_effective_robust_mode(Sparlectra.StateEstimationConfig(robust = true)) === :staged
    @test Sparlectra._se_effective_robust_mode(Sparlectra.StateEstimationConfig(robust = true, robust_mode = :replacement)) === :replacement

    # --- replacement mode (0.10.0): judged on the CONVERGED state, it
    # suppresses exactly the corrupted row (stage 3, fixed suppression
    # sigma), pulls the estimate back toward the truth, and keeps all
    # statistics on the original sigmas. An unreachable threshold leaves
    # the solve bitwise plain.
    resReplHi = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_net = false, robust_mode = :replacement, k_suppress = 1.0e9) do
      runse!(net, meas)
    end
    @test resReplHi.voltages == resPlain.voltages
    resRepl = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_net = false, robust_mode = :replacement, k_suppress = 6.0, suppression_sigma = 2000.0) do
      runse!(net, meas)
    end
    @test resRepl.converged
    @test resRepl.robustRows !== nothing && length(resRepl.robustRows) == 1
    rrep = only(resRepl.robustRows)
    @test rrep.id == bm.id
    @test rrep.stage == 3
    @test rrep.sigma_factor ≈ 2000.0 / bm.sigma
    devRepl = maximum(abs.(resRepl.voltages .- Vref))
    @test devRepl < devPlain
    @test devRepl < 2.0 * devClean
    # J is reported on the ORIGINAL sigmas: removing the bad row's pull
    # makes its residual grow, so J exceeds the plain solve's J
    @test resRepl.objectiveJ > resPlain.objectiveJ
    # virtual rows stay exempt: a ZI pseudo-measurement far off the truth
    # would otherwise be suppressed; the classifier alone proves the guard
    @test Sparlectra._robust_stage(100.0, 1e-7, 3.0, 6.0) == (0, 1e-7)

    # --- k_eliminate exposure: a lowered threshold marks more rows
    # suspicious, a raised one clears the corrupted row from candidacy
    repLoose = with_state_estimation_config(max_iter = 40, tol = 1e-8, k_eliminate = 50.0) do
      validate_measurements(net, meas)
    end
    @test !any(r -> r.suspicious, repLoose.measurement_ranking)
    repTight = with_state_estimation_config(max_iter = 40, tol = 1e-8, k_eliminate = 1.0) do
      validate_measurements(net, meas)
    end
    @test count(r -> r.suspicious, repTight.measurement_ranking) >= count(r -> r.suspicious, report.measurement_ranking)
  end)() end

  return true
end

function test_state_estimation_wilson_hilferty()::Bool
  # The Wilson-Hilferty band test decides for all nu,
  # with the two-sided reasons and the redundancy notes.
  @testset "State estimation Wilson-Hilferty band test" begin (function ()
    # nu = 4: asymmetric acceptance interval. J = 15 exceeds the OLD
    # symmetric upper bound (nu + 3 sqrt(2 nu) = 12.49) but passes WH
    # (upper bound = nu (mu + 3 sigma)^3 = 18.02); the lower bound is above 0.
    v15 = Sparlectra._band_test_verdict(15.0, 4)
    @test v15.passed
    @test abs(v15.z_legacy) > 3.0   # the old test would have failed here
    vlow = Sparlectra._band_test_verdict(0.01, 4)
    @test !vlow.passed
    @test vlow.reason == :low
    vhigh = Sparlectra._band_test_verdict(30.0, 4)
    @test !vhigh.passed
    @test vhigh.reason == :high
    # nu = 0: no redundancy, test skipped and reported as not passed
    v0 = Sparlectra._band_test_verdict(0.0, 0)
    @test !v0.passed
    @test v0.reason == :no_redundancy
    # nu = 10 carries the small-redundancy note; nu = 40 does not
    @test Sparlectra._band_test_verdict(10.0, 10).small_redundancy
    @test !Sparlectra._band_test_verdict(40.0, 40).small_redundancy

    # integration: noise-free synthetic set fails :low, gross error :high
    net = createTest3BusNet()
    runpf!(net, 40, 1e-10, 0; method = :rectangular)
    std = measurementStdDevs(vm = 0.002, pinj = 0.5, qinj = 0.5, pflow = 0.5, qflow = 0.5)
    measClean = generateMeasurementsFromPF(net; noise = false, stddev = std)
    repLow = with_state_estimation_config(max_iter = 30, tol = 1e-8) do
      validate_measurements(net, measClean)
    end
    @test repLow.converged
    @test !repLow.global_consistency
    @test repLow.objective.reason == :low
    @test repLow.objective.small_redundancy   # nu = 16 < 30
    sLow = summarize_se_diagnostics(repLow)
    @test occursin("implausibly small", sLow.reason)
    @test occursin("overestimated", sLow.reason)

    measBad = generateMeasurementsFromPF(net; noise = true, stddev = std, rng = MersenneTwister(7))
    bi = findfirst(m -> m.typ == Sparlectra.VmMeas, measBad)
    bmm = measBad[bi]
    measBad[bi] = Measurement(typ = bmm.typ, value = bmm.value + 10.0 * bmm.sigma, sigma = bmm.sigma, active = bmm.active, busIdx = bmm.busIdx, branchIdx = bmm.branchIdx, direction = bmm.direction, id = bmm.id)
    repHigh = with_state_estimation_config(max_iter = 30, tol = 1e-8) do
      validate_measurements(net, measBad)
    end
    @test !repHigh.global_consistency
    @test repHigh.objective.reason == :high

    # printed report carries the Wilson-Hilferty line
    io = IOBuffer()
    print_se_diagnostics(io, repLow; topN = 3)
    txt = String(take!(io))
    @test occursin("Wilson-Hilferty", txt)
    @test occursin("reason=low", txt)
  end)() end

  return true
end

function test_state_estimation_diagnostics_unification()::Bool
  # Diagnostics share the estimator's extended state
  # definition (B states of released shunts, nu included).
  @testset "State estimation diagnostics unification" begin (function ()
    stdn = measurementStdDevs(vm = 0.002, pinj = 0.5, qinj = 0.5, pflow = 0.5, qflow = 0.5, shuntq = 0.5)
    net = create_se_shunt_net()
    runpf!(net, 40, 1e-10, 0; method = :rectangular)
    meas = generateMeasurementsFromPF(net; includeShuntQ = true, noise = true, stddev = stdn, rng = MersenneTwister(5))
    # gross error elsewhere (a Pflow), released shunt with direct measurement
    badIdx = findfirst(m -> m.typ == Sparlectra.PflowMeas, meas)
    bm = meas[badIdx]
    meas[badIdx] = Measurement(typ = bm.typ, value = bm.value + 10.0 * bm.sigma, sigma = bm.sigma, active = bm.active, busIdx = bm.busIdx, branchIdx = bm.branchIdx, direction = bm.direction, id = bm.id)
    setShuntEstimation!(net; busName = "LoadB")

    report = with_state_estimation_config(max_iter = 40, tol = 1e-8) do
      validate_measurements(net, meas)
    end
    @test report.converged
    # nu counts the B state: m - (2 nbus - 1) - 1
    nbus = length(net.nodeVec)
    m = length(report.measurement_ranking)
    @test report.objective.dof == m - (2 * nbus - 1) - 1
    @test !only(r for r in report.result.shuntEstimates).frozen
    # ranking consistent with the estimator: the gross error tops the list
    @test first(report.measurement_ranking).measurement_index == badIdx

    # frozen case: no direct measurement, diagnostics behave exactly like the
    # estimator (frozen with warning, no B state in nu)
    measNoQ = [mm for mm in meas if mm.typ != Sparlectra.ShuntQMeas]
    repFrozen = Test.@test_logs (:warn, r"no active direct shunt measurement") match_mode = :any with_state_estimation_config(() -> validate_measurements(net, measNoQ); max_iter = 40, tol = 1e-8)
    @test only(r for r in repFrozen.result.shuntEstimates).frozen
    mF = length(repFrozen.measurement_ranking)
    @test repFrozen.objective.dof == mF - (2 * nbus - 1)
  end)() end

  return true
end

function test_state_estimation_measurement_csv()::Bool
  # Measurement CSV v1 core. Lossless roundtrip over all
  # types and location groups, atomic import, line-precise errors, version
  # guard.
  @testset "State estimation measurement CSV v1" begin (function ()
    net = create_se_link_net(shunt_buses = ["B2"])
    runpf!(net, 40, 1e-10, 0; method = :rectangular)
    m = generateMeasurementsFromPF(net; includeImag = true, includeShuntQ = true, includeVa = true, noise = true, rng = MersenneTwister(3))
    append!(net.measurements, m)
    addImagMeasurement!(net; value = 120.0, sigma = 2.0, busName = "B2")   # shunt-bay current
    addPflowMeasurement!(net; value = 12.5, sigma = 0.4, linkNr = 1)        # link flow
    addMeasurement!(net; typ = Sparlectra.VmMeas, value = 1.0, sigma = 0.01, busIdx = 2, active = false, id = "inactive_probe")
    orig = copy(net.measurements)

    dir = mktempdir()
    f = joinpath(dir, "meas.csv")
    w = writeMeasurementsCSV(net; file = f)
    @test w.count == length(orig)
    r = readMeasurementsCSV!(net; file = f, replace = true)
    @test r.total == length(orig)
    # bitwise-equal measurement vectors (all fields, incl. linkIdx and active)
    @test length(net.measurements) == length(orig)
    @test all(a == b for (a, b) in zip(orig, net.measurements))
    @test sum(values(r.counts)) == r.total
    @test haskey(r.counts, Sparlectra.ShuntQMeas) && haskey(r.counts, Sparlectra.ImagMeas)

    # replace = false appends, and appending a set onto itself is exactly the
    # corruption the reader warns about: every quantity is then measured twice
    n0 = length(net.measurements)
    @test_logs (:warn, r"repeat a quantity already measured") readMeasurementsCSV!(net; file = f, replace = false)
    @test length(net.measurements) == 2 * n0
    empty!(net.measurements)
    append!(net.measurements, orig)

    # version guard: v0 and a missing comment reject the file
    bad = joinpath(dir, "bad.csv")
    open(io -> println(io, "# sparlectra-measurements v0"), bad, "w")
    errV = try
      readMeasurementsCSV!(net; file = bad)
      nothing
    catch e
      sprint(showerror, e)
    end
    @test errV !== nothing && occursin("unknown measurement file version", errV)

    # line-precise errors, atomic import: net.measurements stays untouched
    lines = readlines(f)
    lines[5] = replace(lines[5], r"^[A-Za-z]+," => "BogusMeas,")
    broken = joinpath(dir, "broken.csv")
    open(io -> foreach(l -> println(io, l), lines), broken, "w")
    before = copy(net.measurements)
    errL = try
      readMeasurementsCSV!(net; file = broken)
      nothing
    catch e
      sprint(showerror, e)
    end
    @test errL !== nothing
    @test occursin("broken.csv:5", errL)
    @test occursin("unknown measurement type", errL)
    @test net.measurements == before   # nothing imported

    # malformed number, line-precise (field 8 = value in the branch_nr header)
    lines2 = readlines(f)
    parts = split(lines2[6], ",")
    parts[8] = "not_a_number"
    lines2[6] = join(parts, ",")
    broken2 = joinpath(dir, "broken2.csv")
    open(io -> foreach(l -> println(io, l), lines2), broken2, "w")
    errN = try
      readMeasurementsCSV!(net; file = broken2)
      nothing
    catch e
      sprint(showerror, e)
    end
    @test errN !== nothing && occursin("broken2.csv:6", errN) && occursin("not a number", errN)

    # parallel branches (CGMES-typical): branch_nr keeps each measurement on
    # its own branch where pair resolution must reject
    pnet = Net(name = "csv_parallel", baseMVA = 100.0)
    addBus!(net = pnet, busName = "P1", vn_kV = 110.0)
    addBus!(net = pnet, busName = "P2", vn_kV = 110.0)
    addProsumer!(net = pnet, busName = "P1", type = "EXTERNALNETWORKINJECTION", referencePri = "P1", vm_pu = 1.0, va_deg = 0.0)
    addProsumer!(net = pnet, busName = "P2", type = "ENERGYCONSUMER", p = 20.0, q = 5.0)
    addPIModelACLine!(net = pnet, fromBus = "P1", toBus = "P2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = pnet, fromBus = "P1", toBus = "P2", r_pu = 0.02, x_pu = 0.16, b_pu = 0.0, status = 1)
    addPflowMeasurement!(pnet; branchNr = 1, value = 12.0, sigma = 0.5, direction = :from)
    addPflowMeasurement!(pnet; branchNr = 2, value = 8.0, sigma = 0.5, direction = :from)
    porig = copy(pnet.measurements)
    fpar = joinpath(dir, "parallel.csv")
    writeMeasurementsCSV(pnet; file = fpar)
    readMeasurementsCSV!(pnet; file = fpar, replace = true)
    @test [m.branchIdx for m in pnet.measurements] == [1, 2]
    @test all(a == b for (a, b) in zip(porig, pnet.measurements))
    # a branch_nr pointing at a branch with other endpoints is rejected
    plines = readlines(fpar)
    pparts = split(plines[3], ",")
    pparts[5] = "3"
    plines[3] = join(pparts, ",")
    badnr = joinpath(dir, "badnr.csv")
    open(io -> foreach(l -> println(io, l), plines), badnr, "w")
    errB = try
      readMeasurementsCSV!(pnet; file = badnr)
      nothing
    catch e
      sprint(showerror, e)
    end
    @test errB !== nothing && occursin("badnr.csv:3", errB) && occursin("branch_nr", errB)

    # first-cut v1 header (no branch_nr column) is still read via pair
    # resolution: rebuild the roundtrip file without column 5
    legacy = joinpath(dir, "legacy.csv")
    open(legacy, "w") do io
      for (i, l) in enumerate(readlines(f))
        if i == 1
          println(io, l)
        elseif i == 2
          println(io, "type,bus,from_bus,to_bus,link_nr,direction,value,sigma,active,id")
        else
          c = split(l, ","; limit = 11)
          println(io, join(vcat(c[1:4], c[6:end]), ","))
        end
      end
    end
    readMeasurementsCSV!(net; file = legacy, replace = true)
    @test length(net.measurements) == length(orig)
    @test all(a == b for (a, b) in zip(orig, net.measurements))

    # whitespace-decorated bus names (CGMES): the reader strips CSV fields,
    # so resolution must match busDict keys whitespace-insensitively
    @test Sparlectra._measurement_csv_bus_key(net, "B2") == "B2"
    net.busDict["WSPAD "] = net.busDict["B2"]
    @test Sparlectra._measurement_csv_bus_key(net, "WSPAD") == "WSPAD "
    delete!(net.busDict, "WSPAD ")

    # relative sigmas (percent of reading): per row sigma = max(f*|value|,
    # floor), floors bind on tiny readings, Va stays absolute
    rnet = create_se_link_net()
    runpf!(rnet, 40, 1e-10, 0; method = :rectangular)
    fr = measurementStdDevs(vm = 0.005, pinj = 0.02, qinj = 0.01, pflow = 0.02, qflow = 0.01, va = 0.02, imag = 0.02)
    fl = measurementSigmaFloors()
    rmeas = generateMeasurementsFromPF(rnet; includeImag = true, includeVa = true, noise = false, stddev = fr, relativeSigma = true, sigmaFloor = fl)
    for mm in rmeas
      if mm.typ == Sparlectra.VaMeas
        @test mm.sigma == 0.02   # absolute even in relative mode
      else
        @test mm.sigma == max(fr[mm.typ] * abs(mm.value), fl[mm.typ])
      end
    end
    # a dominating floor wins on every row
    bigfl = measurementSigmaFloors(vm = 10.0, pinj = 1e3, qinj = 1e3, pflow = 1e3, qflow = 1e3, imag = 1e6)
    rmeas2 = generateMeasurementsFromPF(rnet; includeImag = true, noise = false, stddev = fr, relativeSigma = true, sigmaFloor = bigfl)
    @test all(mm.sigma == bigfl[mm.typ] for mm in rmeas2)
    # absolute mode is bitwise unchanged by the new keywords' defaults
    a1 = generateMeasurementsFromPF(rnet; noise = false, stddev = fr)
    a2 = generateMeasurementsFromPF(rnet; noise = false, stddev = fr, relativeSigma = false)
    @test all(x == y for (x, y) in zip(a1, a2))

    # header comments round-trip: written after the version line, skipped by
    # the reader, sniffer still recognizes the file
    empty!(rnet.measurements)
    append!(rnet.measurements, a1)
    fcm = joinpath(dir, "commented.csv")
    writeMeasurementsCSV(rnet; file = fcm, headerComments = ["transformer taps at generation:", "transformer branch 1 A-B tap_ratio=1.0"])
    lines_cm = readlines(fcm)
    @test startswith(lines_cm[2], "# transformer taps")
    rcm = readMeasurementsCSV!(rnet; file = fcm, replace = true)
    @test rcm.total == length(a1)
    @test all(x == y for (x, y) in zip(a1, rnet.measurements))

    # component-id fallback: a bus referenced by its comp cID (CGMES-style
    # divergence between busDict key and component id) still resolves
    cid = rnet.nodeVec[geNetBusIdx(net = rnet, busName = "B2")].comp.cID
    @test !isempty(String(cid))
    @test Sparlectra._measurement_csv_bus_key(rnet, String(cid)) == "B2"

    # mRID reference mode (CGMES ENTSO-E UUIDs preserved in net.cgmes_ids):
    # the writer references buses by UUID (name fallback where none is
    # recorded), the reader resolves them back; roundtrip bitwise
    rnet.cgmes_ids["TN|B2"] = "11111111-2222-3333-4444-555555555555"
    @test Sparlectra._bus_mrid(rnet, geNetBusIdx(net = rnet, busName = "B2")) == "11111111-2222-3333-4444-555555555555"
    empty!(rnet.measurements)
    append!(rnet.measurements, a1)
    fmr = joinpath(dir, "mrid.csv")
    writeMeasurementsCSV(rnet; file = fmr, busReference = :mrid)
    @test occursin("11111111-2222-3333-4444-555555555555", read(fmr, String))
    readMeasurementsCSV!(rnet; file = fmr, replace = true)
    @test length(rnet.measurements) == length(a1)
    @test all(x == y for (x, y) in zip(a1, rnet.measurements))
    empty!(rnet.cgmes_ids)
  end)() end

  return true
end

function test_state_estimation_se_chain()::Bool
  # PF started from the SE result. :se_state starts from
  # the estimated voltages (model injections authoritative), :se_snapshot
  # additionally takes the nodal balances (0/1 iterations, slack pickup at
  # tolerance, persistent model untouched). Both error without a preceding SE.
  @testset "State estimation SE-to-PF chain" begin (function ()
    net = create_se_shunt_net()
    iteFlat, ergFlat = runpf!(net, 40, 1e-10, 0; method = :rectangular, opt_flatstart = true)
    @test ergFlat == 0
    std = measurementStdDevs(vm = 1e-4, pinj = 1e-3, qinj = 1e-3, pflow = 1e-3, qflow = 1e-3, shuntq = 1e-3)
    meas = generateMeasurementsFromPF(net; includeShuntQ = true, noise = false, stddev = std)
    res = with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = true) do
      runse!(net, meas)
    end
    @test res.converged

    # :se_state: strictly fewer iterations than the flat start on this case
    r1 = runpf_from_se!(net, 40, 1e-10, 0; mode = :se_state, method = :rectangular)
    @test r1.converged
    @test r1.iterations < iteFlat
    @test r1.mode == :se_state

    # :se_snapshot: immediate convergence, tiny slack pickup, persistent
    # model loads bitwise untouched
    with_state_estimation_config(() -> runse!(net, meas); max_iter = 30, tol = 1e-10, update_net = true)   # re-register after the PF overwrote V
    loadsBefore = [(n._pƩLoad, n._qƩLoad) for n in net.nodeVec]
    r2 = runpf_from_se!(net, 40, 1e-10, 0; mode = :se_snapshot, method = :rectangular)
    @test r2.converged
    @test r2.iterations <= 1
    @test abs(r2.slack_pickup_mw) < 1e-6
    @test abs(r2.slack_pickup_mvar) < 1e-6
    @test loadsBefore == [(n._pƩLoad, n._qƩLoad) for n in net.nodeVec]

    # SE state CSV persistence: a fresh net restores and chains
    f = joinpath(mktempdir(), "se_state.csv")
    with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = true) do
      runse!(net, meas)
    end
    writeSEStateCSV(net; file = f)
    net2 = create_se_shunt_net()
    readSEStateCSV!(net2; file = f)
    r3 = runpf_from_se!(net2, 40, 1e-10, 0; mode = :se_snapshot, method = :rectangular)
    @test r3.converged && r3.iterations <= 1

    # both modes error clearly without a preceding SE result
    net3 = create_se_shunt_net()
    for mode in (:se_state, :se_snapshot)
      err = try
        runpf_from_se!(net3, 40, 1e-10, 0; mode = mode)
        nothing
      catch e
        sprint(showerror, e)
      end
      @test err !== nothing && occursin("no preceding state-estimation", err)
    end

    # measurement/model DISCREPANCY: the case the original noise-free chain
    # test could not distinguish (takeover == model there), which left the
    # phase-5 snapshot takeover inert without any test noticing. With the
    # model drifting after the measurements were taken, the two modes must
    # converge differently: :se_state balances the drift through the slack,
    # :se_snapshot reproduces the estimated point with zero start mismatch.
    net4 = create_se_shunt_net()
    runpf!(net4, 40, 1e-10, 0; method = :rectangular)
    meas4 = generateMeasurementsFromPF(net4; noise = false, stddev = std)
    addProsumer!(net = net4, busName = "LoadA", type = "ENERGYCONSUMER", p = 7.0, q = 0.0)   # post-measurement drift
    with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = true) do
      runse!(net4, meas4)
    end
    nPros4 = length(net4.prosumpsVec)
    rSt = runpf_from_se!(net4, 40, 1e-10, 0; mode = :se_state, method = :rectangular)
    vmSt = [n._vm_pu for n in net4.nodeVec]
    rSn = runpf_from_se!(net4, 40, 1e-10, 0; mode = :se_snapshot, method = :rectangular)
    vmSn = [n._vm_pu for n in net4.nodeVec]
    @test rSt.converged && rSn.converged
    @test abs(rSt.slack_pickup_mw) > 4.0          # the drift lands in the slack
    @test rSn.iterations <= 1                     # zero mismatch at the start
    @test abs(rSn.slack_pickup_mw) < 1e-6
    @test maximum(abs.(vmSt .- vmSn)) > 1e-5      # the modes are measurably different
    @test length(net4.prosumpsVec) == nPros4      # delta prosumers cleaned up
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
    function _two_island_net()
      net = Net(name = "se_islands", baseMVA = 100.0)
      for b in ("A1", "A2", "A3", "B1", "B2")
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
      return net
    end

    net = _two_island_net()
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
    net2 = _two_island_net()
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

function test_state_estimation_ia_measurements()::Bool
  # PMU current-phasor angles (IaMeas, 0.10.0). Alpha recovery from
  # a full phasor set, the paired and unpaired activity gates, observability
  # invariance, the exact angle-wrap helper, and the phasor pair helper.
  @testset "State estimation current-angle measurements (IaMeas)" begin (function ()
    net = Net(name = "ia_test", baseMVA = 100.0)
    for b in ("A", "B", "C")
      addBus!(net = net, busName = b, vn_kV = 110.0)
    end
    addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.02, va_deg = 0.0)
    addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 30.0, q = 10.0)
    addProsumer!(net = net, busName = "C", type = "ENERGYCONSUMER", p = 20.0, q = 6.0)
    addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = net, fromBus = "B", toBus = "C", r_pu = 0.012, x_pu = 0.09, b_pu = 0.0, status = 1)
    addPIModelACLine!(net = net, fromBus = "A", toBus = "C", r_pu = 0.011, x_pu = 0.085, b_pu = 0.0, status = 1)
    ite, erg = runpf!(net, 40, 1e-10, 0; method = :rectangular)
    @test erg == 0
    Vref = buildVoltageVector(net)

    # 1) roundtrip and alpha recovery: current phasors share the PMU time
    # base with the voltage angles; a 3 degree shift comes back as alpha
    meas = generateMeasurementsFromPF(net; includeImag = true, includeIa = true, includeVa = true, vaRefOffsetDeg = 3.0, noise = false)
    nIa = count(m -> m.typ == Sparlectra.IaMeas, meas)
    @test nIa == 2 * length(net.branchVec)
    res = with_state_estimation_config(max_iter = 25, tol = 1e-9, update_net = false) do
      runse!(net, meas)
    end
    @test res.converged
    @test res.vaRefOffsetDeg !== nothing && isapprox(res.vaRefOffsetDeg, 3.0; atol = 1e-6)
    @test res.objectiveJ < 1e-12
    for i in eachindex(Vref)
      @test abs(rad2deg(angle(Vref[i]) - angle(res.voltages[i]))) < 1e-6
    end
    # generated value equals the prediction at the truth (shared helper)
    Sb = Sparlectra.calc_injections(createYBUS(net = net), Vref) .* net.baseMVA
    ia1 = first(m for m in meas if m.typ == Sparlectra.IaMeas)
    @test abs(Sparlectra._measurement_prediction(ia1, net, Vref, Sb, deg2rad(3.0)) - ia1.value) < 1e-9

    # 2) IaMeas ALONE also activates the offset state (with the warning)
    measIaOnly = Test.@test_logs (:warn, r"includeIa without includeImag") match_mode = :any generateMeasurementsFromPF(net; includeImag = false, includeIa = true, includeVa = false, vaRefOffsetDeg = 2.0, noise = false)
    resIa = with_state_estimation_config(max_iter = 25, tol = 1e-9, update_net = false) do
      runse!(net, [m for m in meas if m.typ != Sparlectra.IaMeas && m.typ != Sparlectra.VaMeas] ∪ [m for m in measIaOnly if m.typ == Sparlectra.IaMeas])
    end
    @test resIa.converged
    @test resIa.vaRefOffsetDeg !== nothing && isapprox(resIa.vaRefOffsetDeg, 2.0; atol = 1e-5)

    # 3) paired gate: a weak magnitude drops its angle row with it, dof
    # shrinks by both rows
    weak = copy(meas)
    addImagMeasurement!(weak; net = net, value = 1.0, sigma = 10.0, branchNr = 2, direction = :from, id = "Imag_weak")
    addIaMeasurement!(weak; net = net, value = 10.0, sigma = 0.05, branchNr = 2, direction = :from, id = "Ia_weak")
    resW = Test.@test_logs (:info, r"Ia_weak excluded") match_mode = :any with_state_estimation_config(() -> runse!(net, weak); max_iter = 25, tol = 1e-9, update_net = false)
    @test resW.converged
    @test resW.dof == res.dof   # both extra rows gated away again

    # 4) observability invariance: identical verdict with and without Ia
    noIa = [m for m in meas if m.typ != Sparlectra.IaMeas]
    obsWith = evaluate_global_observability(net, meas)
    obsWithout = evaluate_global_observability(net, noIa)
    @test obsWith.numerical_rank == obsWithout.numerical_rank
    @test obsWith.n_measurements == obsWithout.n_measurements
    @test obsWith.quality == obsWithout.quality

    # 5) exact angle wrap: out-of-range residuals fold once, in-range stay
    # bitwise untouched (also for tiny values), non-angle rows never touched
    rw = [359.0, -359.0, 1e-12, 250.0]
    mw = Measurement[
      Measurement(typ = Sparlectra.IaMeas, value = 0.0, sigma = 1.0, branchIdx = 1, direction = :from),
      Measurement(typ = Sparlectra.VaMeas, value = 0.0, sigma = 1.0, busIdx = 1),
      Measurement(typ = Sparlectra.VaMeas, value = 0.0, sigma = 1.0, busIdx = 2),
      Measurement(typ = Sparlectra.PflowMeas, value = 0.0, sigma = 1.0, branchIdx = 1, direction = :from),
    ]
    Sparlectra._wrap_angle_residuals!(rw, mw)
    @test rw == [-1.0, 1.0, 1e-12, 250.0]

    # 6) phasor pair helper: PMU-prefixed ids, both rows appended
    n0 = length(net.measurements)
    mi, ma = addCurrentPhasorMeasurement!(net; i_A = 120.0, ia_deg = -15.0, branchNr = 1, direction = :from)
    @test length(net.measurements) == n0 + 2
    @test mi.typ == Sparlectra.ImagMeas && ma.typ == Sparlectra.IaMeas
    @test startswith(mi.id, "PMU_Imag_") && startswith(ma.id, "PMU_Ia_")
    empty!(net.measurements)

    # 7) gating notes surface in the diagnostics report and printout
    rep = with_state_estimation_config(max_iter = 25, tol = 1e-8) do
      validate_measurements(net, weak)
    end
    @test any(gn -> gn.reason == :ia_below_current_floor && gn.id == "Ia_weak", rep.gating_notes)
    io = IOBuffer()
    print_se_diagnostics(io, rep; topN = 3)
    @test occursin("Gated measurements", String(take!(io)))
  end)() end

  return true
end

function _tap_gate_net()::Net
  # two voltage levels, meshed through two parallel transformer paths (the
  # closed loop that makes a tap observable), lines on both levels
  net = Net(name = "tap_gate", baseMVA = 100.0)
  for (b, vn) in (("H1", 110.0), ("H2", 110.0), ("L1", 20.0), ("L2", 20.0))
    addBus!(net = net, busName = b, vn_kV = vn)
  end
  addProsumer!(net = net, busName = "H1", type = "EXTERNALNETWORKINJECTION", referencePri = "H1", vm_pu = 1.02, va_deg = 0.0)
  addProsumer!(net = net, busName = "L1", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
  addProsumer!(net = net, busName = "L2", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
  addPIModelACLine!(net = net, fromBus = "H1", toBus = "H2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "L1", toBus = "L2", r_pu = 0.02, x_pu = 0.1, b_pu = 0.0, status = 1)
  addPIModelTrafo!(net = net, fromBus = "H1", toBus = "L1", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  addPIModelTrafo!(net = net, fromBus = "H2", toBus = "L2", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  net.branchVec[3].has_ratio_tap = true
  net.branchVec[3].tap_step = 0.00625
  return net
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

## r2 producing EXACTLY the shift n * phase_step_deg under direction alpha:
## with phi = -shift (the regulating-vector angle) the closed form
## r2 = sin(phi)/sin(alpha - phi) satisfies angle(1 + r2 cis alpha) = phi
function _tap_r2_for_shift(shift_deg::Float64, alpha_deg::Float64)::Float64
  ϕ = -deg2rad(shift_deg)
  return sin(ϕ) / sin(deg2rad(alpha_deg) - ϕ)
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

function test_state_estimation_tap_islands()::Bool
  # Island aggregation of tap estimates. Two
  # disconnected copies of the meshed tap topology, one released trafo per
  # island with different true deviations; the merged result carries one
  # row per island with the island id, the fixation sums add, and the
  # islands stay numerically isolated from each other.
  @testset "State estimation tap estimation across islands" begin (function ()
    step = 0.00625

    function _two_island_tap_net()
      net = Net(name = "tap_islands", baseMVA = 100.0)
      for sfx in ("a", "b")
        for (b, vn) in (("H1", 110.0), ("H2", 110.0), ("L1", 20.0), ("L2", 20.0))
          addBus!(net = net, busName = b * sfx, vn_kV = vn)
        end
        addProsumer!(net = net, busName = "H1" * sfx, type = "EXTERNALNETWORKINJECTION", referencePri = "H1" * sfx, vm_pu = 1.02, va_deg = 0.0)
        addProsumer!(net = net, busName = "L1" * sfx, type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
        addProsumer!(net = net, busName = "L2" * sfx, type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
        addPIModelACLine!(net = net, fromBus = "H1" * sfx, toBus = "H2" * sfx, r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = net, fromBus = "L1" * sfx, toBus = "L2" * sfx, r_pu = 0.02, x_pu = 0.1, b_pu = 0.0, status = 1)
        addPIModelTrafo!(net = net, fromBus = "H1" * sfx, toBus = "L1" * sfx, r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
        addPIModelTrafo!(net = net, fromBus = "H2" * sfx, toBus = "L2" * sfx, r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
      end
      for br in net.branchVec
        br.ratio == 0.0 && continue
        br.has_ratio_tap = true
        br.tap_step = step
      end
      return net
    end

    # branch layout per island block: 1/2 lines, 3/4 trafos (a), 5/6 lines,
    # 7/8 trafos (b). Truth: +2 steps on trafo 3 (island a), -1 on trafo 7.
    tnet = _two_island_tap_net()
    for (k, n) in ((3, 2.0), (7, -1.0))
      br = tnet.branchVec[k]
      ct = Sparlectra._cascade_tap(br, n * step, 0.0, 0.0)
      br.tap_ratio = ct.tap_ratio
      br.phase_shift_deg = ct.phase_shift_deg
    end
    ite, erg = runpf!(tnet, 40, 1e-12, 0; islands_enabled = true)
    @test erg == 0
    meas = generateMeasurementsFromPF(tnet; includeImag = true, noise = false)

    net = _two_island_tap_net()
    setTapEstimation!(net; trafo = 3, mode = :ratio)
    setTapEstimation!(net; trafo = 7, mode = :ratio)
    res = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_net = false) do
      runse!(net, meas)
    end
    @test res.converged
    @test res.islands !== nothing && length(res.islands) == 2
    @test res.tapEstimates !== nothing && length(res.tapEstimates) == 2
    rowA = only(t for t in res.tapEstimates if t.branch == 3)
    rowB = only(t for t in res.tapEstimates if t.branch == 7)
    @test rowA.island != rowB.island
    @test rowA.fixed_step_1 == 2
    @test rowB.fixed_step_1 == -1
    @test rowA.fixed && rowB.fixed
    tf = res.tapFixation
    @test tf.fixed
    @test tf.dof_after == tf.dof_before + 2   # one freed state per island
    @test tf.j_after < 1e-10

    # numerical isolation: releasing ONLY island a's trafo must reproduce
    # island a's voltages bitwise (island b's model discrepancy stays in b)
    netA = _two_island_tap_net()
    setTapEstimation!(netA; trafo = 3, mode = :ratio)
    resA = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_net = false) do
      runse!(netA, meas)
    end
    @test resA.converged
    @test length(resA.tapEstimates) == 1
    @test resA.tapEstimates[1].fixed_step_1 == 2
    @test res.voltages[1:4] == resA.voltages[1:4]
  end)() end

  return true
end

function test_state_estimation_tap_completion()::Bool
  # Tap-estimation completion tests: tap write-back protection, PMU synergy on the
  # released tap, bad data next to a released trafo, the machine-trafo
  # back-calculation, and the se_view release listing.
  @testset "State estimation tap completion" begin (function ()
    step = 0.00625

    function _truth_meas(r1; noise = false, rng = Random.MersenneTwister(1), grossQ = 0.0, includeIa = false)
      tnet = _tap_gate_net()
      br = tnet.branchVec[3]
      ct = Sparlectra._cascade_tap(br, r1, 0.0, 0.0)
      br.tap_ratio = ct.tap_ratio
      br.phase_shift_deg = ct.phase_shift_deg
      ite, erg = runpf!(tnet, 40, 1e-12, 0; method = :rectangular)
      @test erg == 0
      meas = generateMeasurementsFromPF(tnet; includeImag = true, includeIa = includeIa, noise = noise, rng = rng)
      if grossQ != 0.0
        gi = findfirst(m -> m.typ == Sparlectra.QflowMeas && m.branchIdx == 3, meas)
        m0 = meas[gi]
        meas[gi] = Measurement(typ = m0.typ, value = m0.value + grossQ * m0.sigma, sigma = m0.sigma, active = m0.active, busIdx = m0.busIdx, branchIdx = m0.branchIdx, direction = m0.direction, id = m0.id, linkIdx = m0.linkIdx)
      end
      return meas
    end

    # 1) write-back protection: without updateTaps the model stays bitwise
    # untouched; with it the FIXED cascade position lands in the branch and
    # a follow-up run without any release explains the same measurements
    meas = _truth_meas(2 * step)
    net = _tap_gate_net()
    setTapEstimation!(net; trafo = 3, mode = :ratio)
    ratio0 = net.branchVec[3].tap_ratio
    res = with_state_estimation_config(max_iter = 40, tol = 1e-10) do
      runse!(net, meas)
    end
    @test res.converged && res.tapEstimates[1].fixed_step_1 == 2
    @test net.branchVec[3].tap_ratio === ratio0            # bitwise protection
    res = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_taps = true) do
      runse!(net, meas)
    end
    @test res.converged
    expect = net.branchVec[3].ratio / (1.0 + 2 * step)
    @test isapprox(net.branchVec[3].tap_ratio, expect; atol = 1e-14)
    net.branchVec[3].tap_est_mode = :none                  # no release: the written-back model
    res2 = with_state_estimation_config(() -> runse!(net, meas); max_iter = 40, tol = 1e-10)     # explains the set by itself
    @test res2.converged
    @test res2.objectiveJ < 1e-10

    # 2) PMU synergy: on a TWO-REGULATOR release the phase
    # regulator r2 lives in the angle information; an accurate current
    # PHASOR pair on the released trafo shrinks its estimate spread over
    # noisy realizations (a magnitude-only set leaves the angle soft)
    r2true = _tap_r2_for_shift(-1.25, 90.0)   # one shift step down
    function _truth_both()
      tn = _tap_gate_net()
      br = tn.branchVec[3]
      ct = Sparlectra._cascade_tap(br, 2 * step, r2true, 90.0)
      br.tap_ratio = ct.tap_ratio
      br.phase_shift_deg = ct.phase_shift_deg
      ite2, erg2 = runpf!(tn, 40, 1e-12, 0; method = :rectangular)
      @test erg2 == 0
      return tn
    end
    tnB = _truth_both()
    VB = buildVoltageVector(tnB)
    IendB = Sparlectra._end_current_A(tnB, tnB.branchVec[3], :from, VB)
    function _estimates(withPMU::Bool)
      out = Float64[]
      for seed = 1:6
        m = generateMeasurementsFromPF(_truth_both(); includeImag = true, noise = true, rng = Random.MersenneTwister(seed))
        n = _tap_gate_net()
        withPMU && addCurrentPhasorMeasurement!(n; i_A = abs(IendB), ia_deg = rad2deg(angle(IendB)), branchNr = 3, direction = :from, sigmaI = 0.5, sigmaIa = 0.02)
        setTapEstimation!(n; trafo = 3, mode = :both, alpha_deg = 90.0)
        append!(n.measurements, m)
        r = with_state_estimation_config(max_iter = 40, tol = 1e-10) do
          runse!(n)
        end
        @test r.converged
        push!(out, r.tapEstimates[1].electrical_step_2)
      end
      return out
    end
    base = _estimates(false)
    pmu = _estimates(true)
    spread(v) = maximum(abs.(v .- (-1.0)))
    @test spread(pmu) < spread(base)

    # 3) bad data NEXT TO the released trafo: the diagnostics (running on
    # the tap-unified net) name exactly the corrupted row, and the tap
    # estimate stays within half a step of the truth
    meas = _truth_meas(2 * step; noise = true, rng = Random.MersenneTwister(3), grossQ = 12.0)
    net = _tap_gate_net()
    setTapEstimation!(net; trafo = 3, mode = :ratio)
    append!(net.measurements, meas)
    diag = with_state_estimation_config(max_eliminations = 2, max_iter = 40, tol = 1e-10) do
      runse_diagnostics(net)
    end
    @test !isempty(diag.eliminations)
    @test occursin("Qflow", diag.eliminations[1].id)
    res = with_state_estimation_config(max_iter = 40, tol = 1e-10, robust = true) do
      runse!(net)
    end
    @test res.converged
    @test abs(res.tapEstimates[1].electrical_step_1 - 2.0) < 0.5

    # 4) machine-trafo back-calculation: mutually exclusive with a release,
    # needs a registered SE, and recovers the GSU position from the AVR
    # setpoint plus the dispatch against the estimated network voltage
    gnet = _tap_gate_net()
    addBus!(net = gnet, busName = "G1", vn_kV = 10.0)
    addProsumer!(net = gnet, busName = "G1", type = "SYNCHRONOUSMACHINE", p = 20.0, q = 5.0, vm_pu = 1.02)
    addPIModelTrafo!(net = gnet, fromBus = "H2", toBus = "G1", r_pu = 0.001, x_pu = 0.09, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
    gnet.branchVec[5].has_ratio_tap = true
    gnet.branchVec[5].tap_step = step
    # TRUE machine tap: 2 steps; the machine holds 1.02 pu at G1
    brg = gnet.branchVec[5]
    ctg = Sparlectra._cascade_tap(brg, 2 * step, 0.0, 0.0)
    brg.tap_ratio = ctg.tap_ratio
    brg.phase_shift_deg = ctg.phase_shift_deg
    ite, erg = runpf!(gnet, 40, 1e-12, 0; method = :rectangular)
    @test erg == 0
    gmeas = generateMeasurementsFromPF(gnet; includeImag = true, noise = false)
    append!(gnet.measurements, gmeas)
    # the back-calculation never reads the model tap: it scans the tap
    # fraction against the estimated NETWORK voltage plus AVR setpoint,
    # dispatch P, and the MEASURED machine Q (under AVR control Q is not
    # the schedule; the SCADA reading of the truth stands in for it here)
    vg = gnet.nodeVec[5]._vm_pu   # the AVR held terminal voltage of the truth
    Vtrue = buildVoltageVector(gnet)
    yg = Sparlectra.calcAdmittance(brg, brg.comp.cVN, gnet.baseMVA)
    q_scada = imag(Vtrue[5] * conj(yg[3] * Vtrue[2] + yg[4] * Vtrue[5]) * gnet.baseMVA)
    # mutual exclusion and the no-SE precondition give clear errors
    setTapEstimation!(gnet; trafo = 5, mode = :ratio)
    @test_throws ErrorException calcMachineTrafoTapFromSE(gnet; trafo = 5)
    setTapEstimation!(gnet; trafo = 5, enabled = false)
    @test_throws ErrorException calcMachineTrafoTapFromSE(gnet; trafo = 5)
    rse = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_net = true) do
      runse!(gnet)
    end
    @test rse.converged
    bt = calcMachineTrafoTapFromSE(gnet; trafo = 5, v_machine_pu = vg, q_mvar = q_scada)
    @test bt.fixed_step == 2
    @test abs(bt.electrical_step - 2.0) < 0.25
    @test bt.machine_bus == 5
    @test bt.q_residual_mvar < 0.1

    # 5) se_view lists the release with mode and alpha
    vnet = _tap_gate_net()
    setTapEstimation!(vnet; trafo = 3, mode = :both, alpha_deg = 30.0)
    view = se_view(vnet)
    @test length(view.taps) == 1
    @test view.taps[1].branch == 3 && view.taps[1].mode == :both && view.taps[1].alpha_deg == 30.0
    io = IOBuffer()
    print_se_view(io, view)
    out = String(take!(io))
    @test occursin("Tap estimation releases: 1", out)
    @test occursin("mode both, alpha 30.0 deg", out)
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
      readMeasurementsCSV!(net; file = joinpath(scf_dir, "sp_case5.measurements.csv"))
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
    readMeasurementsCSV!(net; file = joinpath(dir, "sp_case5.measurements.baddata.csv"))
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

"""
Targeted critical thinning of the measurement generator. With
`critical_count = k` the generated set ends with at least k critical rows,
stays observable, and the number of removed rows is small: every step
removes the PARTNER of the least redundant row (largest normalized
residual covariance), so one critical row costs about one removal, where
the former greedy sweep removed the least redundant row itself and needed
many removals per critical row. Protected zero-injection rows are never
removed, and the footer comments name the target, the critical rows and
the removed rows. The invalid option combinations throw at construction
of the options, before any case is imported.
"""
function test_measurement_generator_critical_thinning()::Bool
  @testset "Measurement generator critical thinning" begin (function ()
    case = joinpath(dirname(@__DIR__), "data", "scf", "sp_case60.scf.json")
    k = 3
    base = (noise = false, gross_k = 0.0, tap_steps = 0.0, include_i = false, sigma_u_pct = 0.5, sigma_i_pct = 1.0, sigma_p_pct = 1.0, sigma_q_pct = 1.0, sigma_ia_deg = 0.0, seed = 42)
    mktempdir() do d
      out = joinpath(d, "crit.csv")
      g = SparlectraApp._se_generate_measurement_set(case, out, SparlectraApp.MeasurementGeneratorOptions(; base..., critical_count = k))
      @test occursin("critical row(s) after removing", g.critical_note)
      lines = readlines(out)
      target = findfirst(l -> startswith(l, "# critical_target: $(k) reached: "), lines)
      @test target !== nothing
      reached = target === nothing ? 0 : parse(Int, strip(lines[target][length("# critical_target: $(k) reached: ")+1:end]))
      @test reached >= k
      removed_line = findfirst(l -> startswith(l, "# critical_removed: "), lines)
      @test removed_line !== nothing
      removed = removed_line === nothing ? String[] : split(strip(lines[removed_line][length("# critical_removed: ")+1:end]))
      # small: on a fully measured set a critical row costs the local
      # redundancy of its bus (about three rows: the other flow end, the
      # injection, the voltage), where the former greedy sweep needed about
      # eighteen removals per critical row on sp_case188
      @test 1 <= length(removed) <= 4 * k
      @test !any(id -> startswith(id, "ZI"), removed)
      # the written set is what the criticality check reports it to be
      net = Sparlectra.import_case(case, Sparlectra.load_sparlectra_config()).net
      Sparlectra.readMeasurementsCSV!(net; file = out)
      obs = evaluate_global_observability(net)
      @test obs.quality != :not_observable
      @test length(obs.numerical_critical_measurement_indices) >= k
      # the same seed and target regenerate the identical file
      out2 = joinpath(d, "crit2.csv")
      SparlectraApp._se_generate_measurement_set(case, out2, SparlectraApp.MeasurementGeneratorOptions(; base..., critical_count = k))
      @test read(out, String) == read(out2, String)
    end
    # option validation happens at construction, with the user-readable texts
    @test_throws ArgumentError SparlectraApp.MeasurementGeneratorOptions(; base..., critical_count = -1)
    @test_throws ArgumentError SparlectraApp.MeasurementGeneratorOptions(; base..., flow_ends = :none)
    @test_throws ArgumentError SparlectraApp.MeasurementGeneratorOptions(; base..., truth_source = :from_run)
    @test_throws ArgumentError SparlectraApp.MeasurementGeneratorOptions(; base..., truth_source = :from_run, run_id = "r1", run_root = "/tmp", tap_steps = 1.0)
    @test_throws ArgumentError SparlectraApp.MeasurementGeneratorOptions(; base..., passive_sigma = 0.0)
    okopts = SparlectraApp.MeasurementGeneratorOptions(; base...)
    @test okopts.critical_count == 0 && okopts.flow_ends === :both && okopts.passive_as_zi
  end)() end
  return true
end

function run_state_estimation_tests()
  # Aggregates all state-estimation unit tests to keep coverage explicit and ordered.

"""
State estimation on a DTF network case: the format has to travel with the
request. A bare `.DAT` is ambiguous (FOR001 network vs FOR002 reference), so
the detector refuses to guess; the SE service used to detect the format
itself and therefore rejected every DTF case with "Ambiguous .DAT input"
before it ever looked at the measurements (found 2026-09-05 while tracing
the sysimage workload). It now takes `case_format` like the power-flow
service does, and DTF is an accepted SE format.
"""
function test_state_estimation_dtf_service()
  dtf = joinpath(dirname(@__DIR__), "data", "DTF", "FOR001.DAT")
  isfile(dtf) || return true
  cfg = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH
  net = Sparlectra.DTFImporter.build_net(Sparlectra.DTFImporter.read_dtf(dtf))
  runpf!(net, 40, 1e-8, 0; method = :rectangular)
  append!(net.measurements, generateMeasurementsFromPF(net; noise = false))
  mktempdir() do d
    mf = joinpath(d, "dtf.measurements.csv")
    writeMeasurementsCSV(net; file = mf)

    ok = redirect_stdout(devnull) do
      SparlectraApp._run_state_estimation_service(dtf, cfg, joinpath(d, "run_ok"), "dtf_se", mf; case_format = :dtf_for001)
    end
    d_ok = SparlectraApp.to_dict(ok)
    @test d_ok["status"] == "succeeded"
    @test d_ok["metadata"]["run_mode"] == "se"

    # without the format the run must fail with the way out NAMED, not with
    # a bare "unknown format" the caller cannot act on
    amb = redirect_stdout(devnull) do
      SparlectraApp._run_state_estimation_service(dtf, cfg, joinpath(d, "run_amb"), "dtf_se_amb", mf)
    end
    d_amb = SparlectraApp.to_dict(amb)
    @test d_amb["status"] == "failed"
    @test occursin("dtf_for001", d_amb["message"])
  end
  # the import entry point accepts the format and rejects an unsupported one
  @test SparlectraApp._se_import_case(dtf, Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true); requested_format = :dtf_for001).net isa Sparlectra.Net
  return true
end


"""
Tap-estimation fallback: released taps are extra states, and a measurement
set that estimates the voltages cleanly can still be far too thin to pin
them. The run then does not settle at all and the user gets nothing,
although the SAME set works without the taps (seen 2026-09-06 on
case300 with 98 released taps and a CGMES delivery). A non-convergence WITH
released taps therefore freezes the taps back to their model position and
repeats the estimation once, and says so in the log: a silent retry would
hide that the reported tap positions are model values, not estimates.
"""
function test_state_estimation_tap_fallback()::Bool
  @testset "State estimation tap fallback" begin
    case = joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m")
    if !isfile(case)
      println("      SE tap fallback SKIPPED (fixture case missing)")
      return true
    end
    cfg = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH
    mktempdir() do d
      set = joinpath(d, "gen.csv")
      SparlectraApp._se_generate_measurement_set(case, set, SparlectraApp.MeasurementGeneratorOptions(; noise = true, gross_k = 0.0, tap_steps = 0.0,
        include_i = false, sigma_u_pct = 0.5, sigma_i_pct = 1.0, sigma_p_pct = 1.0, sigma_q_pct = 1.0,
        sigma_ia_deg = 0.0, seed = 42))
      # one iteration is not enough for anything, so the first solve cannot
      # converge; the fallback must still deliver a result and say why
      res = redirect_stdout(devnull) do
        SparlectraApp._run_state_estimation_service(case, cfg, joinpath(d, "run"), "tapfb", set;
          tap_estimation = true, max_iter = 1)
      end
      dd = SparlectraApp.to_dict(res)
      # with a cap of one the repeat cannot converge either: what is asserted
      # here is that the fallback RAN and is visible, not that it rescues
      log = read(joinpath(d, "run", "run.log"), String)
      @test occursin("repeating WITHOUT tap estimation", log)
      @test occursin("released transformer tap", log)

      # and the honest case: a normal cap converges, and then no fallback
      # line appears at all
      res_ok = redirect_stdout(devnull) do
        SparlectraApp._run_state_estimation_service(case, cfg, joinpath(d, "run_ok"), "tapok", set;
          tap_estimation = true, max_iter = 50)
      end
      d_ok = SparlectraApp.to_dict(res_ok)
      @test d_ok["status"] == "succeeded"
      log_ok = read(joinpath(d, "run_ok", "run.log"), String)
      @test !occursin("repeating WITHOUT tap estimation", log_ok)
      @test get(d_ok["metadata"], "se_tap_estimation_fallback", false) == false
    end
    # The step limit. A released regulator state is
    # bounded in how far ONE iteration may move it, at a quarter of the
    # changer's declared mechanical travel. Without it the Gauss-Newton step
    # drives r1 toward -1, where the cascade
    # t = t_base/((1+r1)(1+r2 e^{j alpha})) is singular: measured 2026-09-06
    # on a CGMES delivery, one weakly determined transformer ran to -329
    # electrical steps on a band of about -14 to +18 and took the estimation
    # down. Bounding the STATE into that band was measured to be worse than
    # no bound at all (five transformers that had converged stuck to the
    # lower bound and diverged), which is why the limit is on the step.
    begin
      tnet = SparlectraApp._se_import_case_net(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"),
        Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true))
      # ratio != 0.0 is the check setTapEstimation! itself uses; tap_step is
      # NOT a transformer marker, it carries a constructor default on every
      # branch
      trafo = findfirst(b -> b.ratio != 0.0 && b.tap_step > 0.0 && b.tap_min > 0.0 && b.tap_max > 0.0, tnet.branchVec)
      if trafo !== nothing
        setTapEstimation!(tnet; trafo = trafo, mode = :ratio, enabled = true)
        tmap = Sparlectra._build_tap_state_map(tnet, 1)
        @test tmap !== nothing
        limits = Sparlectra._tap_step_limits(tmap, tnet)
        @test !isempty(limits)
        br = tnet.branchVec[trafo]
        full = abs((1.0 / br.tap_min - 1.0) - (1.0 / br.tap_max - 1.0))
        # a quarter of the full travel, and strictly positive so a released
        # regulator can still move
        for (_, lim) in limits
          @test lim > 0.0
          @test lim ≈ 0.25 * full atol = 1e-12
        end
        setTapEstimation!(tnet; trafo = trafo, enabled = false)
      end
    end

    # the fallback must be visible on all THREE
    # surfaces, because a run whose tap positions are MODEL values looks
    # exactly like a successful tap estimation otherwise, and its J measures
    # those model positions. sp_case60 at a cap of 4 is the shipped fixture
    # that produces a SUCCEEDING run with the fallback used: the tap solve
    # needs more iterations than the plain one, so the first attempt fails
    # and the repeat converges (measured 2026-09-06; the same holds for
    # sp_case188, while sp_case14 already converges with the taps).
    mktempdir() do d
      shipped = joinpath(dirname(@__DIR__), "data", "scf", "sp_case60.scf.json")
      set = joinpath(d, "gen60.csv")
      SparlectraApp._se_generate_measurement_set(shipped, set, SparlectraApp.MeasurementGeneratorOptions(; noise = true, gross_k = 0.0, tap_steps = 0.0,
        include_i = false, sigma_u_pct = 0.5, sigma_i_pct = 1.0, sigma_p_pct = 1.0, sigma_q_pct = 1.0,
        sigma_ia_deg = 0.0, seed = 42))
      res = redirect_stdout(devnull) do
        SparlectraApp._run_state_estimation_service(shipped, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH,
          joinpath(d, "fb"), "tapfb60", set; tap_estimation = true, max_iter = 4)
      end
      dd = SparlectraApp.to_dict(res)
      @test dd["status"] == "succeeded"
      @test get(dd["metadata"], "se_tap_estimation_fallback", false) == true

      # surface 3: se_diagnostics.md
      diagmd = read(joinpath(d, "fb", "se_diagnostics.md"), String)
      @test occursin(SparlectraApp._SE_TAP_FALLBACK_NOTE, diagmd)

      # surface 1: the result summary of the Web UI
      summary = SparlectraApp._webui_se_summary(dd)
      @test summary !== nothing
      @test occursin(SparlectraApp._SE_TAP_FALLBACK_NOTE, summary)

      # surface 2: the tap table is REPLACED by the reason, not shown with
      # positions that look estimated
      table = SparlectraApp._webui_se_tap_section(dd)
      @test occursin("did NOT converge", table)
      @test !occursin("<th>Electrical step</th>", table)

      # and the run history must not call a state estimation "rectangular"
      # (seen 2026-09-06): the method comes from the run kind
      @test SparlectraApp._powerflow_run_index_solver(res) == "wls"
    end

    # the helpers the fallback is built from
    net = SparlectraApp._se_import_case_net(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"),
      Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true))
    @test SparlectraApp._any_tap_released(net) == false
    @test SparlectraApp._freeze_all_tap_estimation!(net) == 0
  end
  return true
end

  @testset "State estimation" begin (function ()
    tests = [
      ("WLS", test_state_estimation_wls_first_version),
      ("Measurement add helpers", test_state_estimation_measurement_add_helpers),
      ("Passive bus zero injection helpers", test_state_estimation_passive_bus_zero_injection_helpers),
      ("Bad-data diagnostics", test_state_estimation_bad_data_diagnostics),
      ("PMU angle measurements", test_state_estimation_pmu_va_measurements),
      ("Current-magnitude measurements (ImagMeas)", test_state_estimation_imag_measurements),
      ("Sequential elimination and wii", test_state_estimation_sequential_elimination),
      ("Shunt estimation (case A)", test_state_estimation_shunt_estimation),
      ("Shunt back-calculation (case B)", test_state_estimation_shunt_back_calculation),
      ("Link contraction", test_state_estimation_link_contraction),
      ("W2 link allocation", test_state_estimation_link_allocation_w2),
      ("se_view", test_state_estimation_se_view),
      ("Robust R modification", test_state_estimation_robust_r),
      ("Wilson-Hilferty band test", test_state_estimation_wilson_hilferty),
      ("Diagnostics unification", test_state_estimation_diagnostics_unification),
      ("Measurement CSV v1", test_state_estimation_measurement_csv),
      ("SE-to-PF chain", test_state_estimation_se_chain),
      ("Island-wise SE", test_state_estimation_islands),
      ("Current-angle measurements (IaMeas)", test_state_estimation_ia_measurements),
      ("Tap PF equivalence gate", test_state_estimation_tap_pf_equivalence),
      ("Tap roundtrip and fixation", test_state_estimation_tap_roundtrip),
      ("Tap release guards", test_state_estimation_tap_guards),
      ("Tap estimation across islands", test_state_estimation_tap_islands),
      ("Tap completion (write-back, PMU, machine trafo)", test_state_estimation_tap_completion),
      ("DTF case through the SE service", test_state_estimation_dtf_service),
      ("Tap fallback on non-convergence", test_state_estimation_tap_fallback),
      ("Measurement generator critical thinning", test_measurement_generator_critical_thinning),
      ("Run configuration as an argument (#381)", test_state_estimation_config_argument),
    ]

    for (name, testfn) in tests
      @testset "$name" begin (function ()
        @test _se_run_quiet(testfn) == true
      end)() end
    end
  end)() end
end
