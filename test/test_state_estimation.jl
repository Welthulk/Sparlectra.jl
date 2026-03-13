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
