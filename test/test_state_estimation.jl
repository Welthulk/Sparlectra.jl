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

using Test
using Sparlectra
using Random

function test_state_estimation_wls_first_version()::Bool
  @testset "State estimation WLS first version" begin
    net = createTest3BusNet()

    ite, erg = runpf!(net, 40, 1e-10, 0; method = :polar_full, opt_sparse = true)
    @test erg == 0
    @test ite > 0

    Vref = buildVoltageVector(net)

    std = measurementStdDevs(vm = 1e-5, pinj = 1e-4, qinj = 1e-4, pflow = 1e-4, qflow = 1e-4)
    meas = generateMeasurementsFromPF(net; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false, stddev = std, rng = MersenneTwister(42))

    @test !isempty(meas)
    @test all(m -> m.active, meas)

    result = runse!(net, meas; maxIte = 12, tol = 1e-8, flatstart = true, jacEps = 1e-6, updateNet = true)

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
  end

  return true
end

function test_state_estimation_observability_metrics()::Bool
  @testset "State estimation observability metrics" begin
    net = createTest3BusNet()

    ite, erg = runpf!(net, 40, 1e-10, 0; method = :polar_full, opt_sparse = true)
    @test erg == 0
    @test ite > 0

    std = measurementStdDevs(vm = 1e-5, pinj = 1e-4, qinj = 1e-4, pflow = 1e-4, qflow = 1e-4)
    meas = generateMeasurementsFromPF(net; includeVm = true, includePinj = true, includeQinj = true, includePflow = true, includeQflow = true, noise = false, stddev = std, rng = MersenneTwister(42))

    gobs = evaluate_global_observability(net, meas; flatstart = true, jacEps = 1e-6)
    @test gobs.n_states == 2 * length(net.nodeVec) - 1
    @test gobs.n_measurements == count(m -> m.active, meas)
    @test gobs.redundancy == gobs.n_measurements - gobs.n_states
    @test isapprox(gobs.redundancy_ratio, gobs.n_measurements / gobs.n_states)
    @test gobs.dof == gobs.redundancy
    @test gobs.structural_observable == true
    @test gobs.numerical_observable == true
    @test gobs.quality in (:good, :critical)

    nbus = length(net.nodeVec)
    local_cols = [1, nbus + 1]
    lobs = evaluate_local_observability(net, meas, local_cols; flatstart = true, jacEps = 1e-6)
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
        meas_reduced[i] = Measurement(
          typ = m.typ,
          value = m.value,
          sigma = m.sigma,
          active = false,
          busIdx = m.busIdx,
          branchIdx = m.branchIdx,
          direction = m.direction,
          id = m.id,
        )
      end
    end

    gobs_reduced = evaluate_global_observability(net, meas_reduced; flatstart = true, jacEps = 1e-6)
    @test gobs_reduced.n_measurements < gobs.n_measurements
    @test gobs_reduced.structural_matching <= gobs.structural_matching
    @test gobs_reduced.numerical_rank <= gobs.numerical_rank

    meas_sparse = copy(meas)
    first_kept = false
    for i in eachindex(meas_sparse)
      m = meas_sparse[i]
      keep = !first_kept
      first_kept = true
      meas_sparse[i] = Measurement(
        typ = m.typ,
        value = m.value,
        sigma = m.sigma,
        active = keep,
        busIdx = m.busIdx,
        branchIdx = m.branchIdx,
        direction = m.direction,
        id = m.id,
      )
    end

    gobs_sparse = evaluate_global_observability(net, meas_sparse; flatstart = true, jacEps = 1e-6)
    @test gobs_sparse.quality == :not_observable
    @test gobs_sparse.numerical_observable == false

    @test_throws ErrorException evaluate_local_observability(net, meas_sparse, local_cols; flatstart = true, jacEps = 1e-6)
  end

  return true
end

function test_state_estimation_matrix_observability_helpers()::Bool
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

function test_state_estimation_measurement_add_helpers()::Bool
  @testset "State estimation measurement add helpers" begin
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
  end

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
  addProsumer!(net = net, busName = "Load", type = "ENERGYCONSUMER", p = 40.0, q = 15.0)

  return net
end

function test_state_estimation_passive_bus_zero_injection_helpers()::Bool
  @testset "State estimation passive bus zero-injection helpers" begin
    net = create_state_estimation_passive_transit_net()

    ite, erg = runpf!(net, 40, 1e-10, 0; method = :polar_full, opt_sparse = true)
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

    base_obs = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
    @test base_obs.quality == :critical
    @test base_obs.redundancy == 0
    @test !isempty(base_obs.numerical_critical_measurement_indices)

    added = addZeroInjectionMeasurements!(net; sigma = 1e-6)
    @test length(added) == 2
    @test all(m -> m.busIdx == 2, added)
    @test sort([m.typ for m in added]) == sort([Sparlectra.PinjMeas, Sparlectra.QinjMeas])

    zi_obs = evaluate_global_observability(net; flatstart = true, jacEps = 1e-6)
    @test zi_obs.quality == :good
    @test zi_obs.redundancy == 2
    @test isempty(zi_obs.numerical_critical_measurement_indices)
    @test isempty(zi_obs.structural_critical_measurement_indices)

    result = runse!(deepcopy(net); maxIte = 12, tol = 1e-8, flatstart = true, jacEps = 1e-6, updateNet = false)
    @test result.converged == true
  end

  return true
end

function test_state_estimation_bad_data_diagnostics()::Bool
  @testset "State estimation bad-data diagnostics" begin
    net = createTest3BusNet()

    ite, erg = runpf!(net, 40, 1e-10, 0; method = :polar_full, opt_sparse = true)
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

    report = validate_measurements(net, meas; maxIte = 12, tol = 1e-8, flatstart = true, jacEps = 1e-6, normalizedThreshold = 3.0)
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

    diag = runse_diagnostics(net, meas; deactivate_and_rerun = true, maxIte = 12, tol = 1e-8, flatstart = true, jacEps = 1e-6, normalizedThreshold = 3.0)
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
  end

  return true
end

function run_state_estimation_tests()
  @testset "State estimation" begin
    tests = [
      ("WLS", test_state_estimation_wls_first_version),
      ("Observability metrics", test_state_estimation_observability_metrics),
      ("Matrix observability helpers", test_state_estimation_matrix_observability_helpers),
      ("Measurement add helpers", test_state_estimation_measurement_add_helpers),
      ("Passive bus zero injection helpers", test_state_estimation_passive_bus_zero_injection_helpers),
      ("Bad-data diagnostics", test_state_estimation_bad_data_diagnostics),
    ]

    for (name, testfn) in tests
      @testset "$name" begin
        @test testfn() == true
      end
    end
  end
end
