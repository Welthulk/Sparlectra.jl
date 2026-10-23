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
import ZipArchives

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

  addPIModelACLine!(net = net, fromBus = "Slack", toBus = "Transit", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "Transit", toBus = "Load", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)

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

"""
`state_estimation.zero_injection = constraint` against `on` (weighted ZI
rows at 1 kW): on the shipped sp_case14 and sp_case60 with noisy sets whose
passive buses carry only ZI rows, the same state within 1e-6 pu (measured
4.2e-7 and below), the same iterations
and dof, the constraint rows hold to 1e-6 MW and are marked as constraints
in the ranking, never suspicious. A wrong constraint (ZI at a load bus) plus
a gross error: the bad-data elimination never removes a constraint. A
two-island net with a passive bus per island: each island holds its
constraint and recovers the power flow.
"""
function test_state_estimation_zero_injection_constraint()::Bool
  @testset "Zero injection as equality constraint" begin (function ()
    for name in ("sp_case14", "sp_case60")
      path = joinpath(dirname(@__DIR__), "data", "scf", name * ".scf.json")
      net = SparlectraApp._se_import_case_net(path, Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, path).config)
      _, erg = runpf!(net, 60, 1e-8, 0; method = :rectangular)
      @test erg == 0
      meas = generateMeasurementsFromPF(net; noise = true, relativeSigma = true, rng = MersenneTwister(7))
      passive = Set(findPassiveBuses(net))
      @test !isempty(passive)
      filter!(m -> !(m.typ in (Sparlectra.PinjMeas, Sparlectra.QinjMeas) && m.busIdx in passive), meas)
      @test Sparlectra._add_missing_zero_injections!(meas, net) == 2 * length(passive)
      runs = Dict{Symbol,Any}()
      for mode in (:on, :constraint)
        n2 = deepcopy(net)
        r = with_state_estimation_config(zero_injection = mode) do
          runse!(n2, meas)
        end
        @test r.converged
        runs[mode] = (r = r, vm = Float64[nd._vm_pu for nd in n2.nodeVec])
      end
      @test runs[:on].r.iterations == runs[:constraint].r.iterations
      @test runs[:on].r.dof == runs[:constraint].r.dof
      @test maximum(abs.(runs[:on].vm .- runs[:constraint].vm)) < 1e-6
      d = with_state_estimation_config(zero_injection = :constraint) do
        validate_measurements(net, meas)
      end
      con = [row for row in d.measurement_ranking if row.constraint]
      # every ZI row is a constraint; on sp_case60 the passive buses 7 and 21
      # sit in one closed-link cluster and enter as its aggregate, which
      # carries the zero-injection flag of its members
      aggs = name == "sp_case60" ? ["LINKAGG_Pinj_rep_7", "LINKAGG_Qinj_rep_7"] : String[]
      for a in aggs
        @test only(row for row in d.measurement_ranking if row.id == a).constraint
      end
      @test length(con) == count(row -> startswith(String(row.id), "ZI"), d.measurement_ranking) + length(aggs)
      @test all(row -> abs(row.residual) < 1e-6 && row.normalized_residual == 0.0 && !row.suspicious, con)
      # with on the same rows are measurements
      don = with_state_estimation_config(zero_injection = :on) do
        validate_measurements(net, meas)
      end
      @test !any(row -> row.constraint, don.measurement_ranking)
      # under on the aggregate keeps its id and is protected by the flag
      for a in aggs
        @test Sparlectra._elimination_protected(only(m for m in Sparlectra._se_prepare(net, meas).meas if m.id == a))
      end
    end

    # J bookkeeping: constraint rows are not part of J. After one iteration
    # from the flat start the constraints are not yet met, so their weighted
    # residuals are large; J must equal the sum over the measurement rows.
    path14 = joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")
    net14 = SparlectraApp._se_import_case_net(path14, Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, path14).config)
    _, erg = runpf!(net14, 60, 1e-8, 0; method = :rectangular)
    @test erg == 0
    m14 = generateMeasurementsFromPF(net14; noise = true, relativeSigma = true, rng = MersenneTwister(7))
    p14 = Set(findPassiveBuses(net14))
    filter!(m -> !(m.typ in (Sparlectra.PinjMeas, Sparlectra.QinjMeas) && m.busIdx in p14), m14)
    Sparlectra._add_missing_zero_injections!(m14, net14)
    r1 = with_state_estimation_config(zero_injection = :constraint, max_iter = 1, update_net = false) do
      runse!(deepcopy(net14), m14)
    end
    @test length(r1.residuals) == length(m14)   # no contraction on sp_case14: rows in input order
    jM = sum(m.zero_injection ? 0.0 : r1.residuals[i]^2 / m.sigma^2 for (i, m) in enumerate(m14))
    jC = sum(m.zero_injection ? r1.residuals[i]^2 / m.sigma^2 : 0.0 for (i, m) in enumerate(m14))
    @test jC > 1e-3 * jM   # the check below would see constraint rows in J
    @test isapprox(r1.objectiveJ, jM; rtol = 1e-12)

    # a real injection measurement whose id happens to start with ZI (a
    # measurement table uses the meter mRID as id) stays a weighted row
    mreal = copy(m14)
    k = findfirst(m -> m.typ == Sparlectra.PinjMeas && !m.zero_injection && m.value < -5.0, mreal)
    g = mreal[k]
    mreal[k] = Measurement(typ = g.typ, value = g.value, sigma = g.sigma, busIdx = g.busIdx, id = "ZI_real")
    dreal = with_state_estimation_config(zero_injection = :constraint) do
      validate_measurements(net14, mreal)
    end
    @test !only(row for row in dreal.measurement_ranking if row.id == "ZI_real").constraint
    @test !Sparlectra._elimination_protected(mreal[k])

    # dependent constraints (one balance stated twice): the singular KKT
    # step falls back to least squares and the run says why
    mdup = copy(m14)
    zi1 = first(m for m in m14 if m.zero_injection)
    push!(mdup, Measurement(typ = zi1.typ, value = 0.0, sigma = zi1.sigma, busIdx = zi1.busIdx, id = zi1.id * "_dup", zero_injection = true))
    @test_logs (:warn, r"KKT matrix of the zero-injection constraints was singular, dependent constraints") match_mode = :any with_state_estimation_config(zero_injection = :constraint, update_net = false) do
      runse!(deepcopy(net14), mdup)
    end

    # a wrong constraint (zero injection at a load bus) and one gross error:
    # the elimination works on the telemetry and never on a constraint
    path = joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")
    net = SparlectraApp._se_import_case_net(path, Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, path).config)
    _, erg = runpf!(net, 60, 1e-8, 0; method = :rectangular)
    @test erg == 0
    meas = generateMeasurementsFromPF(net; noise = true, relativeSigma = true, rng = MersenneTwister(11))
    loadBus = first(m.busIdx for m in meas if m.typ == Sparlectra.PinjMeas && m.value < -5.0)
    filter!(m -> !(m.typ in (Sparlectra.PinjMeas, Sparlectra.QinjMeas) && m.busIdx == loadBus), meas)
    addZeroInjectionMeasurements!(meas; net = net, busIdxs = [loadBus])
    k = findfirst(m -> m.typ == Sparlectra.PflowMeas && abs(m.value) > 5.0, meas)
    g = meas[k]
    meas[k] = Measurement(typ = g.typ, value = g.value + 30.0 * g.sigma, sigma = g.sigma, busIdx = g.busIdx, branchIdx = g.branchIdx, direction = g.direction, id = g.id, active = g.active)
    dd = with_state_estimation_config(zero_injection = :constraint) do
      runse_diagnostics(net, meas)
    end
    @test !isempty(dd.eliminations)
    @test !any(t -> startswith(String(t.id), "ZI"), dd.eliminations)
    @test all(row -> !row.suspicious, (row for row in dd.final_diagnostics.measurement_ranking if row.constraint))

    # two islands, one passive transit bus each, constraint in both
    function _two_island_transit_net()
      tnet = Net(name = "se_zi_islands", baseMVA = 100.0)
      for sfx in ("a", "b")
        for b in ("Slack", "Transit", "Load")
          addBus!(net = tnet, busName = b * sfx, vn_kV = 110.0)
        end
        addPIModelACLine!(net = tnet, fromBus = "Slack" * sfx, toBus = "Transit" * sfx, r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
        addPIModelACLine!(net = tnet, fromBus = "Transit" * sfx, toBus = "Load" * sfx, r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
        addProsumer!(net = tnet, busName = "Slack" * sfx, type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "Slack" * sfx)
        addProsumer!(net = tnet, busName = "Load" * sfx, type = "ENERGYCONSUMER", p = sfx == "a" ? 10.0 : 6.0, q = sfx == "a" ? 4.0 : 2.0)
      end
      return tnet
    end
    inet = _two_island_transit_net()
    _, erg = runpf!(inet, 40, 1e-10, 0; method = :rectangular, islands_enabled = true)
    @test erg == 0
    @test findPassiveBuses(inet) == [2, 5]
    truth = Float64[nd._vm_pu for nd in inet.nodeVec]
    std = measurementStdDevs(vm = 1e-5, pinj = 1e-5, qinj = 1e-5, pflow = 1e-5, qflow = 1e-5)
    all_meas = generateMeasurementsFromPF(inet; noise = false, stddev = std)
    imeas = Measurement[m for m in all_meas if (m.typ == Sparlectra.VmMeas && m.busIdx in (1, 4)) ||
                                              (m.typ in (Sparlectra.PflowMeas, Sparlectra.QflowMeas) && ((m.branchIdx in (1, 3) && m.direction == :from) || (m.branchIdx in (2, 4) && m.direction == :to)))]
    @test length(imeas) == 10
    @test Sparlectra._add_missing_zero_injections!(imeas, inet) == 4
    snet = _two_island_transit_net()
    ri = with_state_estimation_config(zero_injection = :constraint, max_iter = 20, tol = 1e-8, jac_eps = 1e-6) do
      runse!(snet, imeas)
    end
    @test ri.converged
    @test ri.islands !== nothing && length(ri.islands) == 2
    @test ri.objectiveJ < 1e-6
    @test maximum(abs.(Float64[nd._vm_pu for nd in snet.nodeVec] .- truth)) < 1e-6
    di = with_state_estimation_config(zero_injection = :constraint, max_iter = 20, tol = 1e-8, jac_eps = 1e-6) do
      validate_measurements(_two_island_transit_net(), imeas)
    end
    icon = [row for row in di.measurement_ranking if row.constraint]
    @test length(icon) == 4
    @test sort(unique(row.island for row in icon)) == [1, 2]
    @test all(row -> abs(row.residual) < 1e-6, icon)
    @test occursin("4 row(s) solved as equality constraints", sprint(io -> print_se_diagnostics(io, di; format = :markdown)))
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
    # nothing was eliminated, so a row over the threshold is "over k", not BAD
    @test occursin("over k", txt)
    @test !occursin("BAD", txt)

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
    # BAD marks exactly the eliminated rows of the ranking; the rows the
    # gross error smeared into stay "over k"
    bad_rows = filter(l -> startswith(l, "| ") && endswith(l, "| BAD |"), split(txt2_md, '\n'))
    @test length(bad_rows) == length(diag.eliminations) == 1
    @test startswith(only(bad_rows), "| $(idx) |")
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
    # The K report is on by default since 0.17.4, so a set beyond the m×m
    # limit must not fail the diagnostics any more: the columns are skipped
    # with a warning, omega stays nothing, the ranking keeps its rows. One
    # nonzero per row makes G diagonal, so the Takahashi path runs cheaply.
    (function ()
      mbig = Sparlectra._SE_FULL_OMEGA_MAX_M + 1
      nbig = 50
      Hbig = SparseArrays.sparse(1:mbig, mod1.(1:mbig, nbig), 1.0, mbig, nbig)
      rdbig = @test_logs (:warn, r"K report .* skipped") Sparlectra._residual_diagnostics(Hbig, zeros(Float64, mbig), ones(Float64, mbig); need_full_omega = true, minStates = 1)
      @test rdbig.omega === nothing
      @test length(rdbig.wii) == mbig
      # beyond the dense state cap with few measurements (case1354pegase:
      # 2707 states failed the whole run with "no affordable route"): no
      # error and no warning, the Takahashi route serves w_ii and the K
      # columns come blockwise; a Takahashi refusal there lands on the
      # blockwise route instead of the error
      nwide = Sparlectra._SE_DENSE_LINALG_MAX_N + 1
      Hwide = SparseArrays.sparse(1:nwide, 1:nwide, 1.0, nwide, nwide)
      rdwide = @test_logs Sparlectra._residual_diagnostics(Hwide, zeros(Float64, nwide), ones(Float64, nwide); need_full_omega = true, minStates = 1)
      @test rdwide.omega === nothing
      @test rdwide.omega_path === :takahashi
      @test length(rdwide.wii) == nwide
      rdblk = Sparlectra._residual_diagnostics(Hwide, zeros(Float64, nwide), ones(Float64, nwide); minStates = typemax(Int))
      @test rdblk.omega_path === :blockwise
      @test rdblk.wii == rdwide.wii
      @test Sparlectra._residual_correlation_max_blockwise(Hwide, ones(Float64, nwide)) == zeros(Float64, nwide)
    end)()
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
  addPIModelACLine!(net = net, fromBus = "Slack", toBus = "LoadA", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "LoadA", toBus = "LoadB", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "Slack", toBus = "LoadB", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
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
    # update_shunts is on by default since 0.17.4; this run pins the switched-off case
    res = with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = false, update_shunts = false) do
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
    # with the write-back switched off the estimation never touches the model:
    # the full complex model admittance is bitwise unchanged after the run
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
    # a :low band is not bad data: the sequential elimination does not start
    # (until 0.32.1 it ran on every band failure and reported the budget or
    # the missing suspects instead of the reason)
    diagLow = with_state_estimation_config(max_iter = 30, tol = 1e-8, max_eliminations = 3) do
      runse_diagnostics(net, measClean)
    end
    @test diagLow.stop_reason == :band_low
    @test isempty(diagLow.eliminations)

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

# two islands without a connection: A1-A2-A3 (meshed) and B1-B2, each with
# its own reference
function _se_two_island_net()
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

# two disconnected copies of the meshed tap topology (branches per island
# block: 1/2 lines, 3/4 trafos (a), 5/6 lines, 7/8 trafos (b)), every
# transformer with the fraction grid of `step`
function _se_two_island_tap_net(step::Float64)::Net
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

# largest deviation in MW/MVAr of the branch flows (both ends) and losses of
# `net` from those of the solved truth net `tnet` (same topology); Inf when
# a value is missing, for example on a net that never got flows. The truth
# flows come from calcNetLosses! on the solved truth voltages (runpf! leaves
# the branch flows empty).
function _se_flow_deviation(net::Net, tnet::Net)::Float64
  _val(x) = x === nothing ? Inf : x
  calcNetLosses!(tnet)
  dev = 0.0
  for (br, tb) in zip(net.branchVec, tnet.branchVec)
    for (a, b) in ((br.fBranchFlow, tb.fBranchFlow), (br.tBranchFlow, tb.tBranchFlow))
      (a === nothing || b === nothing) && return Inf
      dev = max(dev, abs(_val(a.pFlow) - _val(b.pFlow)), abs(_val(a.qFlow) - _val(b.qFlow)))
    end
    dev = max(dev, abs(_val(br.pLosses) - _val(tb.pLosses)), abs(_val(br.qLosses) - _val(tb.qLosses)))
  end
  return dev
end

function test_state_estimation_islands_parallel()::Bool
  # The island-wise SE on worker threads against the serial island loop. Run
  # in the se profile and, repeated with four threads, in the threads profile
  # (tools/run_gates.sh threads).
  @testset "State estimation island-wise, parallel equals serial" begin (function ()
    net = _se_two_island_net()
    ite, erg = runpf!(net, 40, 1e-10, 0; islands_enabled = true)
    @test erg == 0
    meas = generateMeasurementsFromPF(net; noise = false)

    # 0.30.2: islands estimated on their own tasks equal the serial run.
    # The fixture's islands have 3 and 2 buses, so the size threshold is
    # lowered to 1; the serial reference switches runtime.parallel off.
    function _se_islands_with(parallel::Bool)
      active = Sparlectra.active_sparlectra_config()
      pf = Sparlectra._copy_powerflow_with(active.powerflow; islands_parallel_min_buses = 1)
      rt = Sparlectra.RuntimeConfig(parallel = Sparlectra.ParallelRuntimeConfig(enabled = parallel))
      cfg = Sparlectra._sparlectra_config_with(active; powerflow = pf, runtime = rt)
      return Sparlectra.with_sparlectra_config(cfg) do
        with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = false) do
          runse!(net, meas)
        end
      end
    end
    se_serial = _se_islands_with(false)
    for _ = 1:5
      se_par = _se_islands_with(true)
      @test se_par.voltages == se_serial.voltages
      @test se_par.residuals == se_serial.residuals
      @test se_par.objectiveJ == se_serial.objectiveJ
      @test [i.iterations for i in se_par.islands] == [i.iterations for i in se_serial.islands]
    end
    parallel_comparison_report("SE, two islands (5 repeats)")

    # tap regulators on worker threads (0.60.0): island a releases a set
    # derived from the branch fields (fraction grid, truth +2 steps), island
    # b a stored set (CGMES-like ratio changer of 1 %, truth -3 steps);
    # update_taps writes the fixed steps back through the island mirror. Rows
    # and the written-back taps equal the serial run, repeated.
    stored(step) = Sparlectra.SETapChangerSet(ratio_n = 1.0, changers = [Sparlectra.SETapChanger(name = "RTC", mrid = "_rtc_b", kind = :ratio, role = :ratio, side = :from, du = 0.01, neutral = 0, low = -8, high = 8, step = step)], origin = :cgmes)
    function _tap_island_net(truth::Bool)
      n = _se_two_island_tap_net(0.00625)
      a, b = n.branchVec[3], n.branchVec[7]
      truth && (a.tap_ratio = a.ratio / (1.0 + 2 * 0.00625))
      # the CGMES ratio formula 1 + (step - neutral) * stepVoltageIncrement / 100
      b.se_tap_changers = stored(truth ? -3 : 0)
      truth && (b.tap_ratio = 1.0 + (-3) * 0.01)
      return n
    end
    tmeas = _tap_truth_meas(_tap_island_net(true); islands = true)
    function _tap_islands_with(parallel::Bool)
      n = _tap_island_net(false)
      setTapEstimation!(n; trafo = 3, mode = :ratio)
      setTapEstimation!(n; trafo = 7, mode = :ratio)
      active = Sparlectra.active_sparlectra_config()
      pf = Sparlectra._copy_powerflow_with(active.powerflow; islands_parallel_min_buses = 1)
      rt = Sparlectra.RuntimeConfig(parallel = Sparlectra.ParallelRuntimeConfig(enabled = parallel))
      cfg = Sparlectra._sparlectra_config_with(active; powerflow = pf, runtime = rt)
      r = Sparlectra.with_sparlectra_config(cfg) do
        with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = false, update_taps = true) do
          runse!(n, tmeas)
        end
      end
      return r, (n.branchVec[3].tap_ratio, n.branchVec[7].tap_ratio, n.branchVec[7].se_tap_changers.changers[1].step)
    end
    tap_serial, wb_serial = _tap_islands_with(false)
    @test [t.step_fixed for t in tap_serial.tapEstimates] == [2.0, -3.0]
    @test wb_serial[1] == 1.0 / (1.0 + 2 * 0.00625) && isapprox(wb_serial[2], 0.97; rtol = 1e-14) && wb_serial[3] == -3.0
    for _ = 1:5
      tap_par, wb_par = _tap_islands_with(true)
      @test tap_par.tapEstimates == tap_serial.tapEstimates
      @test tap_par.voltages == tap_serial.voltages
      @test wb_par == wb_serial
    end
    parallel_comparison_report("SE with tap regulators, two islands (5 repeats)")
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

function _tap_gate_net(; taps::Union{Nothing,PowerTransformerTaps} = nothing, phase_taps::Union{Nothing,PhaseTapChangerModel} = nothing, model::Symbol = :ideal)::Net
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
  # `taps`, `phase_taps`: typed models on the first transformer, which then
  # takes its changers from the models (`model` = model.tap_changer_model)
  addPIModelTrafo!(net = net, fromBus = "H1", toBus = "L1", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1, taps = taps, phase_taps = phase_taps, tap_changer_model = model)
  addPIModelTrafo!(net = net, fromBus = "H2", toBus = "L2", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  net.branchVec[3].has_ratio_tap = true
  taps === nothing && (net.branchVec[3].tap_step = 0.00625)
  return net
end

# the gate net with a declared nameplate changer pair on branch 3: the
# fraction ratio changer of 0.625 % (-14..17, the constructor band) at step
# `nr` and a Delta-u changer of 1 % along `psi` (-10..10) at step `np`
function _tap_gate_pair_net(nr::Real, np::Real; psi::Float64 = 90.0)::Net
  net = _tap_gate_net()
  applyTapNameplate!(net.branchVec[3]; tap_step = 0.00625, tap_min_step = -14.0, tap_max_step = 17.0, tap_current_step = Float64(nr), phase_du_step = 0.01, phase_min_step = -10.0, phase_max_step = 10.0, phase_current_step = Float64(np), psi_deg = psi)
  return net
end

# a small in-memory CGMES store of tap changer objects; the tests evaluate
# the truth with the importer's own functions on it (cgmes_mapping.jl)
function _tap_cim_store(objs)
  C = Sparlectra.CGMESImporter
  byclass = Dict{Symbol,Vector{String}}()
  for o in objs
    push!(get!(Vector{String}, byclass, o.class), o.mrid)
  end
  return C.CGMESStore(Dict(o.mrid => o for o in objs), byclass, C.CGMESFileInfo[], "2.4.15", Set{String}(), C.CGMESMultiRefs())
end
_tap_cim(m, cls, attrs; refs = Dict{Symbol,String}()) = Sparlectra.CGMESImporter.CIMObject(m, cls, Dict{Symbol,String}(attrs), refs, "u1")

# the meshed gate net plus a radial branch: bus R (20 kV, load) hangs on
# bus L2 through branch 5, a transformer whose changers the guard tests set
function _tap_radial_net()::Net
  net = _tap_gate_net()
  net.branchVec[3].has_ratio_tap = false
  addBus!(net = net, busName = "R", vn_kV = 20.0)
  addProsumer!(net = net, busName = "R", type = "ENERGYCONSUMER", p = 8.0, q = 2.0)
  addPIModelTrafo!(net = net, fromBus = "L2", toBus = "R", r_pu = 0.004, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  net.branchVec[4].has_ratio_tap = false
  net.branchVec[5].has_ratio_tap = false
  return net
end

# exact measurements of `net` after a power flow (both flow ends, currents,
# optionally Va and Ia; `islands` solves a net of several islands)
function _tap_truth_meas(net::Net; includeVa::Bool = false, includeIa::Bool = false, noise::Bool = false, rng = Random.MersenneTwister(11), islands::Bool = false)
  ite, erg = runpf!(net, 40, 1e-12, 0; method = :rectangular, islands_enabled = islands)
  @test erg == 0
  return generateMeasurementsFromPF(net; includeImag = true, includeVa = includeVa, includeIa = includeIa, noise = noise, rng = rng)
end

# one tap estimation with the impedance rule `imp` and the run settings
function _tap_se_run(net::Net, meas; imp::Symbol = :auto, max_iter::Int = 40, kw...)
  return with_state_estimation_config(max_iter = max_iter, tol = 1e-10, update_net = false, tap_estimation = Sparlectra.TapEstimationConfig(impedance_correction = imp); kw...) do
    runse!(net, meas)
  end
end

# Transformer branch with nameplate data (applyTapNameplate!) on the meshed
# tap topology; ratio and angle are the neutral point the steps count from
function _tcd_nameplate_branch(; ratio::Float64 = 1.032, angle::Float64 = 2.0, kw...)::Branch
  br = _tap_gate_net().branchVec[3]
  br.ratio = ratio
  br.tap_ratio = ratio
  br.angle = angle
  br.phase_shift_deg = angle
  applyTapNameplate!(br; kw...)
  return br
end

# Transformer branch whose winding carries typed tap models, resolved at
# build (taps_derived)
function _tcd_typed_branch(; taps = nothing, phase_taps = nothing, model::Symbol = :ideal)::Branch
  net = Net(name = "tcd_typed", baseMVA = 100.0)
  addBus!(net = net, busName = "H", vn_kV = 110.0)
  addBus!(net = net, busName = "L", vn_kV = 20.0)
  addProsumer!(net = net, busName = "H", type = "EXTERNALNETWORKINJECTION", referencePri = "H", vm_pu = 1.0, va_deg = 0.0)
  addPIModelTrafo!(net = net, fromBus = "H", toBus = "L", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.02, shift_deg = 1.5, status = 1, taps = taps, phase_taps = phase_taps, tap_changer_model = model)
  return net.branchVec[1]
end

# central difference of f at s (scalar or complex valued, or a tuple)
_tcd_fd(f, s::Float64; h::Float64 = 1e-4) = (f(s + h) .- f(s - h)) ./ (2h)
_tcd_close(a, b; rtol::Float64 = 1e-7, atol::Float64 = 1e-11) = all(abs.(a .- b) .<= rtol .* abs.(b) .+ atol)

function test_tap_changer_data()::Bool
  # Tap changer data of the tap estimation: the set evaluation
  # s -> (ratio, shift) of every changer kind and side against the
  # importers' own functions at every whole step of a band (the set code is
  # never its own truth), the analytic slopes and tangents against central
  # differences at an off-step s, the impedance slopes against central
  # differences, the derived sets of nameplate rows, typed windings and
  # constructor defaults, and the start position search.
  @testset "Tap changer data: s -> t(s)" begin (function ()
    S = Sparlectra
    rs(set, s) = S._tap_set_ratio_shift(set, s)
    # ratio relative 1e-12, angle absolute 1e-12 deg against a truth (R, Phi)
    agrees(got, R, Phi) = isapprox(got[1], R; rtol = 1e-12) && abs(got[2] - Phi) <= 1e-12

    # --- validation: an ArgumentError names the field
    @test_throws ArgumentError S.SETapChanger(kind = :quadrature, role = :phase, side = :from, neutral = 0, low = -1, high = 1, step = 0)
    @test_throws ArgumentError S.SETapChanger(kind = :asymmetrical, role = :phase, side = :from, du = 0.01, neutral = 0, low = -1, high = 1, step = 0)
    @test_throws ArgumentError S.SETapChanger(kind = :ratio, role = :ratio, side = :to, du = 0.01, neutral = 0, low = 2, high = 1, step = 0)
    thrown = @test_throws ArgumentError S.SETapChanger(kind = :linear, role = :phase, side = :from, neutral = 0, low = -1, high = 1, step = 0)
    @test occursin("dphi_deg", thrown.value.msg)
    add = S.SETapChanger(kind = :ratio, role = :ratio, side = :from, additive = true, step_abs = 0.01, du = 0.01, neutral = 0, low = -2, high = 2, step = 0)
    rat = S.SETapChanger(kind = :ratio, role = :ratio, side = :to, du = 0.01, neutral = 0, low = -2, high = 2, step = 0)
    @test_throws ArgumentError S.SETapChangerSet(ratio_n = 1.0, changers = [rat, add], origin = :cgmes)

    # --- (a) nameplate fraction grid: base / (1 + n tap_step), side :to
    let set = S._effective_tap_changers(_tcd_nameplate_branch(tap_step = 0.0125, tap_min_step = -8.0, tap_max_step = 8.0, tap_current_step = 3.0))
      @test set.origin === :fields && length(set.changers) == 1
      ch = set.changers[1]
      @test (ch.kind, ch.side, ch.low, ch.high) == (:ratio, :to, -8, 8) && isapprox(ch.step, 3.0; atol = 1e-12)
      @test all(agrees(rs(set, [n]), _tcd_nameplate_branch(tap_step = 0.0125, tap_min_step = -8.0, tap_max_step = 8.0, tap_current_step = Float64(n)).tap_ratio, 2.0) for n in -8:8)
    end
    # --- (a) nameplate linear grid: base + n tap_step, additive
    let set = S._effective_tap_changers(_tcd_nameplate_branch(tap_step = 0.01, tap_min_step = -5.0, tap_max_step = 5.0, tap_current_step = -2.0, tap_ratio_grid = :linear))
      ch = set.changers[1]
      @test ch.additive && (ch.low, ch.high) == (-5, 5) && isapprox(ch.step, -2.0; atol = 1e-12)
      @test all(agrees(rs(set, [n]), _tcd_nameplate_branch(tap_step = 0.01, tap_min_step = -5.0, tap_max_step = 5.0, tap_current_step = Float64(n), tap_ratio_grid = :linear).tap_ratio, 2.0) for n in -5:5)
    end
    # --- (b) nameplate pair: fraction ratio changer and Delta-u changer at psi 60, both divide
    pair(nr, np; grid = :fraction) = _tcd_nameplate_branch(tap_step = 0.0125, tap_min_step = -4.0, tap_max_step = 4.0, tap_current_step = Float64(nr), phase_du_step = 0.01, phase_min_step = -10.0, phase_max_step = 10.0, phase_current_step = Float64(np), psi_deg = 60.0, tap_ratio_grid = grid)
    for grid in (:fraction, :linear)
      set = S._effective_tap_changers(pair(3, -7; grid = grid))
      @test [c.kind for c in set.changers] == [:ratio, :asymmetrical]
      @test set.changers[2].side === :to && set.changers[2].alpha_deg == 60.0 && (set.changers[2].low, set.changers[2].high) == (-10, 10)
      @test isapprox(set.changers[1].step, 3.0; atol = 1e-9) && isapprox(set.changers[2].step, -7.0; atol = 1e-9)
      @test S._tap_release_mode(pair(3, -7; grid = grid)) === :both
      ok = true
      for nr in -4:4, np in -10:10
        tb = pair(nr, np; grid = grid)
        ok &= agrees(rs(set, [nr, np]), tb.tap_ratio, tb.phase_shift_deg)
      end
      @test ok
    end
    # --- (b) nameplate Delta-u changer alone (no ratio changer, has_ratio_tap false)
    du_only(np) = _tcd_nameplate_branch(phase_du_step = 0.01, phase_min_step = -10.0, phase_max_step = 10.0, phase_current_step = Float64(np), psi_deg = 60.0)
    let br = du_only(4), set = S._effective_tap_changers(br)
      @test !br.has_ratio_tap && [c.kind for c in set.changers] == [:asymmetrical] && set.changers[1].side === :to
      @test (set.changers[1].low, set.changers[1].high) == (-10, 10) && isapprox(set.changers[1].step, 4.0; atol = 1e-9)
      @test all(agrees(rs(set, [np]), du_only(np).tap_ratio, du_only(np).phase_shift_deg) for np in -10:10)
    end
    # --- (f) nameplate degree form: angle + n step, magnitude untouched
    let set = S._effective_tap_changers(_tcd_nameplate_branch(phase_step_deg = 1.5, phase_min_step = -6.0, phase_max_step = 6.0, phase_current_step = 4.0, psi_deg = 30.0))
      @test [c.kind for c in set.changers] == [:linear] && set.changers[1].side === :from && (set.changers[1].low, set.changers[1].high) == (-6, 6)
      @test all(agrees(rs(set, [n]), 1.032, _tcd_nameplate_branch(phase_step_deg = 1.5, phase_min_step = -6.0, phase_max_step = 6.0, phase_current_step = Float64(n), psi_deg = 30.0).phase_shift_deg) for n in -6:6)
      @test rs(set, [2.5])[1] === 1.032
    end

    # --- constructor defaults: the fraction grid 0.9/1.1 of 0.00625, no undeclared phase changer
    let net = _tap_gate_net()
      set = S._effective_tap_changers(net.branchVec[4])
      @test length(set.changers) == 1 && set.changers[1].side === :to && (set.changers[1].low, set.changers[1].high) == (-14, 17)
      @test S._tap_release_mode(net.branchVec[4]) === :ratio
      @test S._effective_tap_changers(net.branchVec[1]) === nothing && S._tap_release_mode(net.branchVec[1]) === nothing
      net.branchVec[4].has_ratio_tap = false
      @test isempty(S._effective_tap_changers(net.branchVec[4]).changers) && S._tap_release_mode(net.branchVec[4]) === nothing
    end

    # --- typed windings against resolve_winding_taps at every step
    ratio_taps() = PowerTransformerTaps(Vn_kV = 110.0, step = 12, lowStep = 1, highStep = 19, neutralStep = 10, voltageIncrement_kV = 1.65)
    asym(conv) = PhaseTapChangerModel(kind = :asymmetrical, step = 3, lowStep = -10, highStep = 10, neutralStep = 0, voltage_step_increment = 0.012, winding_connection_angle_deg = 60.0, x_min = 0.05, x_max = 0.08, convention = conv)
    symm(conv) = PhaseTapChangerModel(kind = :symmetrical, step = -2, lowStep = -8, highStep = 8, neutralStep = 0, voltage_step_increment = 0.02, x_min = 0.05, x_max = 0.09, convention = conv)
    table = [TapTablePoint(step = k, ratio = 1.0 + 0.002 * k^2, angle_deg = 3.0 * k - 0.1 * k^2, x_pu = 0.06 + 0.001 * k^2) for k in -4:5]
    tabular() = PhaseTapChangerModel(kind = :tabular, step = 2, neutralStep = 0, table = table)
    inphase() = PhaseTapChangerModel(kind = :asymmetrical, step = 4, lowStep = -9, highStep = 9, neutralStep = 0, voltage_step_increment = 0.01, winding_connection_angle_deg = 180.0)
    typed_cases = [
      ("ratio + asymmetrical reciprocal", ratio_taps(), asym(:reciprocal_from_side), :ideal, [:ratio, :asymmetrical], [:from, :to]),
      ("asymmetrical direct, rx correction", nothing, asym(:direct_regulating_vector), :impedance_correction, [:asymmetrical], [:from]),
      ("symmetrical reciprocal", nothing, symm(:reciprocal_from_side), :ideal, [:symmetrical], [:to]),
      ("symmetrical direct", ratio_taps(), symm(:direct_regulating_vector), :ideal, [:ratio, :symmetrical], [:from, :from]),
      ("tabular", nothing, tabular(), :ideal, [:tabular], [:from]),
      ("DTF in-phase psi 180, rx correction", nothing, inphase(), :impedance_correction, [:ratio], [:to]),
    ]
    for (label, t, m, model, kinds, sides) in typed_cases
      br = _tcd_typed_branch(taps = t, phase_taps = m, model = model)
      set = S._effective_tap_changers(br)
      w = br.tap_winding
      @test set.origin === :typed && [c.kind for c in set.changers] == kinds && [c.side for c in set.changers] == sides
      # the start position is the model step and reproduces the branch
      s0 = S._tap_s_from_branch(set, br)
      @test s0 == [c.step - c.neutral for c in set.changers]
      res0 = S.resolve_winding_taps(w; tap_changer_model = model)
      F0 = S._tap_impedance_factors(set, s0)
      ok = true
      rsteps = t === nothing ? [0] : collect(t.lowStep:t.highStep)
      psteps = m.kind === :tabular ? [p.step for p in m.table] : collect(m.lowStep:m.highStep)
      for rstep in rsteps, pstep in psteps
        # the truth moves the model steps the way a tap controller does and
        # resolves the winding (the resolver reads a table's reactance at the
        # model step, not at its phase_step keyword)
        t === nothing || (w.taps.step = rstep)
        w.phase_taps.step = pstep
        res = S.resolve_winding_taps(w; tap_changer_model = model)
        s = Float64[]
        t === nothing || push!(s, rstep - t.neutralStep)
        push!(s, pstep - m.neutralStep)
        ok &= agrees(rs(set, s), res.ratio, res.shift_deg)
        # impedance: r scales with rx_factor, x with the model's X(alpha) where it has one
        F = S._tap_impedance_factors(set, s)
        xr = rr -> rr.x_alpha === nothing ? rr.rx_factor : rr.x_alpha
        ok &= isapprox(F[1] / F0[1], res.rx_factor / res0.rx_factor; rtol = 1e-12) && isapprox(F[2] / F0[2], xr(res) / xr(res0); rtol = 1e-12)
      end
      @test ok
      ok || println("      typed case failed: ", label)
    end

    # --- stored sets (CGMES-like) against the importer's own functions
    # (cgmes_mapping.jl): the truths come from tap changer objects in a small
    # in-memory store, evaluated by _phaseTapRatioShift_atstep (asymmetrical,
    # symmetrical, linear), _ratioTapTable and _ratioTapCorrectionAt (ratio
    # table). The importer applies a changer on the from side as it is and on
    # end 2 as 1 / t with the angle negated (_applyEndTaps).
    C = S.CGMESImporter
    cim(m, cls, attrs; refs = Dict{Symbol,String}()) = C.CIMObject(m, cls, Dict{Symbol,String}(attrs), refs, "u1")
    ptc_asym = cim("_ptc_a", :PhaseTapChangerAsymmetrical, [:neutralStep => "13", :lowStep => "1", :highStep => "25", :voltageStepIncrement => "1.25", :windingConnectionAngle => "90"])
    ptc_sym = cim("_ptc_s", :PhaseTapChangerSymmetrical, [:neutralStep => "0", :lowStep => "-8", :highStep => "8", :voltageStepIncrement => "2"])
    ptc_lin = cim("_ptc_l", :PhaseTapChangerLinear, [:neutralStep => "0", :lowStep => "-10", :highStep => "10", :stepPhaseShiftIncrement => "1.2"])
    rtc_ratio(k) = 1.0 + 0.011 * (k - 3) + 0.0007 * (k - 3)^2
    rtc_obj = cim("_rtc", :RatioTapChanger, [:neutralStep => "3", :lowStep => "-2", :highStep => "8", :stepVoltageIncrement => "1"]; refs = Dict(:RatioTapChangerTable => "_rtct"))
    rtc_rows = [cim("_rtct_$k", :RatioTapChangerTablePoint, [:step => string(k), :ratio => string(rtc_ratio(k))]; refs = Dict(:RatioTapChangerTable => "_rtct")) for k in 1:5]
    # a one-point ratio table for the slope check further down
    rtc1_obj = cim("_rtc1", :RatioTapChanger, [:neutralStep => "0", :lowStep => "-5", :highStep => "5", :stepVoltageIncrement => "1"]; refs = Dict(:RatioTapChangerTable => "_rtct1"))
    rtc1_row = cim("_rtct1_0", :RatioTapChangerTablePoint, [:step => "0", :ratio => "1.004"]; refs = Dict(:RatioTapChangerTable => "_rtct1"))
    objs = vcat([ptc_asym, ptc_sym, ptc_lin, rtc_obj, rtc1_obj], rtc_rows, [rtc1_row])
    byclass = Dict{Symbol,Vector{String}}()
    for o in objs
      push!(get!(Vector{String}, byclass, o.class), o.mrid)
    end
    store = C.CGMESStore(Dict(o.mrid => o for o in objs), byclass, C.CGMESFileInfo[], "2.4.15", Set{String}(), C.CGMESMultiRefs())
    imp(ptc, k) = C._phaseTapRatioShift_atstep(store, ptc, k)
    # PhaseTapChangerAsymmetrical: from side multiplies, end 2 divides and negates the angle
    asym_ch(side) = S.SETapChanger(name = "PTC", kind = :asymmetrical, role = :phase, side = side, du = 0.0125, alpha_deg = 90.0, neutral = 13, low = 1, high = 25, step = 16)
    tn = 1.0526315789
    from_set = S.SETapChangerSet(ratio_n = tn, changers = [asym_ch(:from)], origin = :cgmes)
    to_set = S.SETapChangerSet(ratio_n = tn, changers = [asym_ch(:to)], origin = :cgmes)
    @test all(agrees(rs(from_set, [k - 13]), tn * imp(ptc_asym, k)[1], imp(ptc_asym, k)[2]) for k in 1:25)
    @test all(agrees(rs(to_set, [k - 13]), tn / imp(ptc_asym, k)[1], -imp(ptc_asym, k)[2]) for k in 1:25)
    # RatioTapChangerTable: the point ratio inside the table, outside the
    # linear formula of stepVoltageIncrement
    rtc_pts = [S.SETapPoint(step = k, ratio = rtc_ratio(k), r = 1.0 + 0.01 * k, x = 1.0 + 0.02 * k) for k in 1:5]
    rtc = S.SETapChanger(name = "RTC", kind = :tabular, role = :ratio, side = :from, neutral = 3, low = -2, high = 8, step = 4, table = rtc_pts, outside_du = 0.01, impedance = :table)
    rtc_set = S.SETapChangerSet(ratio_n = 1.0, changers = [rtc], origin = :cgmes)
    rtc_table = C._ratioTapTable(store, rtc_obj)
    rtc_truth(k) = C._ratioTapCorrectionAt(rtc_obj, rtc_table, Float64(k))[1]
    @test length(rtc_table) == 5
    @test all(agrees(rs(rtc_set, [k - 3]), rtc_truth(k), 0.0) for k in -2:8)
    @test S._tap_band(rtc) == (-5.0, 5.0)
    # between the end point and the next whole step the curve is continuous
    @test isapprox(rs(rtc_set, [2.0 + 1e-9])[1], rtc_truth(5); atol = 1e-10) && isapprox(rs(rtc_set, [3.0 - 1e-9])[1], rtc_truth(6); atol = 1e-10)
    # a one-point ratio table: on its point the slope is the one of the
    # interval to the right, linear to the formula value 1 + 0.01 at s = 1,
    # not 0
    one = S.SETapChanger(kind = :tabular, role = :ratio, side = :from, neutral = 0, low = -5, high = 5, step = 0, table = [S.SETapPoint(step = 0, ratio = 1.004)], outside_du = 0.01)
    rtc1_table = C._ratioTapTable(store, rtc1_obj)
    @test all(isapprox(S._tap_changer_factor(one, Float64(k))[1], C._ratioTapCorrectionAt(rtc1_obj, rtc1_table, Float64(k))[1]; rtol = 1e-12) for k in -5:5)
    @test _tcd_close(S._tap_changer_slope(one, 0.0), (1.01 - 1.004, 0.0)) && _tcd_close(S._tap_changer_slope(one, 0.0), S._tap_changer_slope(one, 0.5))
    # PhaseTapChangerLinear (the set's nominal angle -5 deg adds to the
    # changer angle) and Symmetrical, both on the from side
    lin_set = S.SETapChangerSet(ratio_n = 1.0, shift_n_deg = -5.0, changers = [S.SETapChanger(kind = :linear, role = :phase, side = :from, dphi_deg = 1.2, neutral = 0, low = -10, high = 10, step = 0)], origin = :cgmes)
    @test all(agrees(rs(lin_set, [k]), imp(ptc_lin, k)[1], -5.0 + imp(ptc_lin, k)[2]) for k in -10:10)
    sym_set = S.SETapChangerSet(ratio_n = 1.0, changers = [S.SETapChanger(kind = :symmetrical, role = :phase, side = :from, du = 0.02, neutral = 0, low = -8, high = 8, step = 0)], origin = :cgmes)
    @test all(agrees(rs(sym_set, [k]), imp(ptc_sym, k)[1], imp(ptc_sym, k)[2]) for k in -8:8)

    # --- slopes of every kind against central differences, off-step
    xtan = (0.05, 0.07, 9.0)
    phase_pts = [S.SETapPoint(step = k, ratio = 1.0 + 0.003 * k^2, angle_deg = 2.5 * k, r = 1.0 + 0.004 * k^2, x = 1.0 + 0.01 * k^2, g = 1.0 - 0.002 * k, b = 1.0 + 0.003 * k) for k in -3:3]
    kinds = [
      S.SETapChanger(kind = :ratio, role = :ratio, side = :to, du = -0.0121, neutral = 0, low = -9, high = 9, step = 0),
      S.SETapChanger(kind = :asymmetrical, role = :phase, side = :from, du = 0.01, alpha_deg = 60.0, neutral = 0, low = -9, high = 9, step = 0, x_curve = xtan, x_curve_form = :tan2, impedance = :x_curve, rx_correction = true),
      S.SETapChanger(kind = :asymmetrical, role = :phase, side = :to, du = 0.015, alpha_deg = 90.0, neutral = 0, low = -9, high = 9, step = 0),
      S.SETapChanger(kind = :symmetrical, role = :phase, side = :from, du = 0.02, neutral = 0, low = -9, high = 9, step = 0, x_curve = (0.1, 0.13, 20.0), impedance = :x_curve),
      S.SETapChanger(kind = :tabular, role = :phase, side = :from, neutral = 0, low = -3, high = 3, step = 0, table = phase_pts, impedance = :table),
      S.SETapChanger(kind = :linear, role = :phase, side = :to, dphi_deg = 2.0, neutral = 0, low = -9, high = 9, step = 0, x_curve = (0.1, 0.13, 18.0), impedance = :x_curve),
      S.SETapChanger(kind = :ratio, role = :ratio, side = :from, du = 0.01, neutral = 0, low = -9, high = 9, step = 0, rx_correction = true),
    ]
    for ch in kinds, s in (2.37, -1.61, 4.2)
      @test _tcd_close(S._tap_changer_slope(ch, s), _tcd_fd(x -> S._tap_changer_factor(ch, x), s))
      single = S.SETapChangerSet(ratio_n = 1.03, shift_n_deg = 1.0, changers = [ch], origin = :cgmes)
      t1 = x -> (v = rs(single, [x]); v[1] * cis(deg2rad(v[2])))
      @test _tcd_close(S._tap_set_tangent(single, [s], 1), _tcd_fd(t1, s))
      @test _tcd_close(S._tap_impedance_slope(single, [s], 1), _tcd_fd(x -> S._tap_impedance_factors(single, [x]), s))
    end
    # composed sets: the tangent and impedance slope of each changer with
    # the others held (additive followed by a Delta-u changer, fraction
    # pair, curve plus table)
    sets = [
      S._effective_tap_changers(pair(2, -3; grid = :linear)),
      S._effective_tap_changers(pair(2, -3)),
      S.SETapChangerSet(ratio_n = 0.98, shift_n_deg = -2.0, changers = [kinds[2], kinds[5]], origin = :cgmes),
    ]
    for set in sets, c in 1:2
      sbase = [1.37, -2.41]
      at = x -> (v = copy(sbase); v[c] = x; v)
      t2 = x -> (v = rs(set, at(x)); v[1] * cis(deg2rad(v[2])))
      @test _tcd_close(S._tap_set_tangent(set, sbase, c), _tcd_fd(t2, sbase[c]))
      @test _tcd_close(S._tap_impedance_slope(set, sbase, c), _tcd_fd(x -> S._tap_impedance_factors(set, at(x)), sbase[c]))
    end

    # --- TapStepImpedance (:tsi): the recorded factor at each whole step,
    # the slope against central differences with the other changer held at s0
    let br = _tcd_nameplate_branch(tap_step = 0.0125, tap_min_step = -4.0, tap_max_step = 4.0, tap_current_step = 1.0, phase_step_deg = 1.5, phase_min_step = -6.0, phase_max_step = 6.0, phase_current_step = 0.0, psi_deg = 30.0)
      steps = -4:4
      vals = [1.032 / (1.0 + n * 0.0125) for n in steps]
      facs = [(1.0 + 0.01 * n, 1.0 + 0.02 * n + 0.001 * n^2, 1.0, 1.0 - 0.005 * n) for n in steps]
      br.tap_step_impedance = S.TapStepImpedance(br.r_pu, br.x_pu, 0.0, 0.0, 0.0, 0.0, vals, facs, Float64[], NTuple{4,Float64}[])
      set = S._effective_tap_changers(br)
      tsi = br.tap_step_impedance
      @test [c.impedance for c in set.changers] == [:tsi, :none]
      s0 = [1.0, 0.0]
      @test all(_tcd_close(S._tap_impedance_factors(set, [Float64(n), 0.0]; s0 = s0, tsi = tsi), facs[n + 5]; rtol = 1e-12, atol = 1e-12) for n in steps)
      # the actuator is read with the phase changer at s0: moving it does not move the factors
      @test S._tap_impedance_factors(set, [2.3, 3.0]; s0 = s0, tsi = tsi) == S._tap_impedance_factors(set, [2.3, 0.0]; s0 = s0, tsi = tsi)
      @test _tcd_close(S._tap_impedance_slope(set, [2.3, 0.0], 1; s0 = s0, tsi = tsi), _tcd_fd(x -> S._tap_impedance_factors(set, [x, 0.0]; s0 = s0, tsi = tsi), 2.3))
      @test S._tap_impedance_slope(set, [2.3, 0.0], 2; s0 = s0, tsi = tsi) == (0.0, 0.0, 0.0, 0.0)
    end

    # --- start position search on a moved branch
    let br = _tap_gate_net().branchVec[3]
      br.ratio = tn
      br.se_tap_changers = from_set
      # the stored step reproduces the branch: returned as it is
      br.tap_ratio, br.phase_shift_deg = tn * imp(ptc_asym, 16)[1], imp(ptc_asym, 16)[2]
      @test S._tap_s_from_branch(S._effective_tap_changers(br), br) == [3.0]
      # moved to step 20 after the import (a controller, a write-back)
      br.tap_ratio, br.phase_shift_deg = tn * imp(ptc_asym, 20)[1], imp(ptc_asym, 20)[2]
      s = S._tap_s_from_branch(from_set, br)
      @test s !== nothing && isapprox(s[1], 7.0; atol = 1e-8)
      # off the curve: a pure quadrature changer cannot give this branch
      br.tap_ratio = 1.2
      @test S._tap_s_from_branch(from_set, br) === nothing
      # two ratio changers of one branch cannot be told apart from t
      two = S.SETapChangerSet(ratio_n = 1.0, changers = [rat, S.SETapChanger(kind = :ratio, role = :ratio, side = :from, du = 0.01, neutral = 0, low = -2, high = 2, step = 0)], origin = :cgmes)
      @test S._tap_separable_mask(two, [0.0, 0.0]) == [true, false]
      br.tap_ratio, br.phase_shift_deg = 1.0 / 1.01, 0.0
      @test S._tap_s_from_branch(two, br) === nothing
    end

    # --- degree-form phase changers of an SCF file with alpha_deg 0 are
    # declared by their row (0.60.0): the changer list decides the release,
    # not psi (until 0.60.0 psi 0 kept them out). Every transformer of the
    # Zeng/Chiang case carries a ratio row and a degree row at alpha 0.
    let net = importSCF(joinpath(dirname(@__DIR__), "data", "scf", "case14_zeng_p306_activeSet_A.scf.json"))
      trafos = [br for br in net.branchVec if br.ratio != 0.0]
      @test !isempty(trafos) && all(br -> br.tap_est_alpha_deg == 0.0 && br.phase_step_deg == 1.25, trafos)
      @test all(br -> S._tap_release_mode(br) === :both && [c.kind for c in S._effective_tap_changers(br).changers] == [:ratio, :linear], trafos)
      # the constructor's own phase grid is not declared
      @test S._tap_release_mode(_tap_gate_net().branchVec[3]) === :ratio
    end

    # --- small helpers
    @test S._tap_band(kinds[5]) == (-3.0, 3.0) && S._tap_band(asym_ch(:from)) == (-12.0, 12.0)
    @test S._tap_scale(kinds[1]) == 0.0121 && S._tap_scale(kinds[6]) == deg2rad(2.0)
    @test S._tap_deviation_changer(sets[3]) == 1 && S._tap_deviation_changer(S.SETapChangerSet(ratio_n = 1.0, changers = [kinds[3], kinds[1]], origin = :cgmes)) == 2
    @test S._tap_deviation_changer(S.SETapChangerSet(ratio_n = 1.0, origin = :cgmes)) == 0
  end)() end
  return true
end

# The changer cases of the tap estimation tests, each on
# branch 3 of the meshed gate net: `build(net, n, zstep)` puts the branch at
# step n (relative to neutral) through the importer's own path (nameplate
# rows through applyTapNameplate!, typed models through
# resolve_branch_taps!, CGMES-like stored sets with the ratio, angle and
# impedance of the CGMES functions on an in-memory store), with the
# impedance of step n (`zstep`) or of the start step 0. `imp` says whether
# the case has per-step impedance data, so auto and off differ.
function _tap_changer_cases()
  C = Sparlectra.CGMESImporter
  S = Sparlectra
  ptc(m, cls, attrs) = _tap_cim(m, cls, attrs)
  asym60 = ptc("_pa60", :PhaseTapChangerAsymmetrical, [:neutralStep => "0", :lowStep => "-10", :highStep => "10", :voltageStepIncrement => "1.0", :windingConnectionAngle => "60", :xMin => "10", :xMax => "13"])
  asym90 = ptc("_pa90", :PhaseTapChangerAsymmetrical, [:neutralStep => "0", :lowStep => "-10", :highStep => "10", :voltageStepIncrement => "1.0", :windingConnectionAngle => "90", :xMin => "10", :xMax => "14"])
  sym = ptc("_ps", :PhaseTapChangerSymmetrical, [:neutralStep => "0", :lowStep => "-8", :highStep => "8", :voltageStepIncrement => "2", :xMin => "10", :xMax => "12"])
  lin = ptc("_pl", :PhaseTapChangerLinear, [:neutralStep => "0", :lowStep => "-10", :highStep => "10", :stepPhaseShiftIncrement => "1.2", :xMin => "10", :xMax => "12.5"])
  tab = _tap_cim("_pt", :PhaseTapChangerTabular, [:neutralStep => "0", :lowStep => "-4", :highStep => "5"]; refs = Dict(:PhaseTapChangerTable => "_ptt"))
  # a non-uniform table: ratio and angle bend with the step, r and x move
  tabpts = [_tap_cim("_ptt_$k", :PhaseTapChangerTablePoint, [:step => string(k), :ratio => string(1.0 + 0.002 * k^2), :angle => string(3.0 * k - 0.1 * k^2), :r => string(2.0 * k), :x => string(3.0 * k + 0.5 * k^2)]; refs = Dict(:PhaseTapChangerTable => "_ptt")) for k in -4:5]
  rtc = ptc("_rtc", :RatioTapChanger, [:neutralStep => "10", :lowStep => "1", :highStep => "21", :stepVoltageIncrement => "1.25"])
  store = _tap_cim_store(vcat([asym60, asym90, sym, lin, tab, rtc], tabpts))
  curve(o) = (c = C._phaseTapXCurve(store, o); (c.xmin, c.xmax, c.amax))
  function stored!(br, ch, obj, n, zstep, neutral)
    br.se_tap_changers = S.SETapChangerSet(ratio_n = br.ratio, changers = [ch(neutral + n)], origin = :cgmes)
    m, a = obj.class == :RatioTapChanger ? (C._ratioTapCorrectionAt(obj, nothing, Float64(neutral + n))[1], 0.0) : C._phaseTapRatioShift_atstep(store, obj, neutral + n)
    br.tap_ratio = br.ratio * m
    br.phase_shift_deg = br.angle + a
    zstep || return br
    if obj.class == :PhaseTapChangerTabular
      f, f0 = C._tapTablePointFactors(store, obj, Float64(n)), C._tapTablePointFactors(store, obj, 0.0)
      br.r_pu *= f[1] / f0[1]
      br.x_pu *= f[2] / f0[2]
    elseif obj.class != :RatioTapChanger
      br.x_pu *= C._phaseTapXCurveFactor(store, obj, Float64(n)) / C._phaseTapXCurveFactor(store, obj, 0.0)
    end
    return br
  end
  asymch(o, alpha) = st -> S.SETapChanger(name = o.mrid, kind = :asymmetrical, role = :phase, side = :from, du = 0.01, alpha_deg = alpha, neutral = 0, low = -10, high = 10, step = st, x_curve = curve(o), impedance = :x_curve)
  tabch = st -> S.SETapChanger(name = "_pt", kind = :tabular, role = :phase, side = :from, neutral = 0, low = -4, high = 5, step = st, table = [S.SETapPoint(step = k, ratio = 1.0 + 0.002 * k^2, angle_deg = 3.0 * k - 0.1 * k^2, r = 1.0 + 2.0 * k / 100, x = 1.0 + (3.0 * k + 0.5 * k^2) / 100) for k in -4:5], impedance = :table)
  typed_asym() = PhaseTapChangerModel(kind = :asymmetrical, step = 0, lowStep = -10, highStep = 10, neutralStep = 0, voltage_step_increment = 0.01, winding_connection_angle_deg = 60.0, convention = :reciprocal_from_side)
  cases = [
    (name = "a nameplate fraction grid", mode = :ratio, n = 3, imp = false, net = () -> _tap_gate_net(), build = (net, n, z) -> applyTapNameplate!(net.branchVec[3]; tap_step = 0.0125, tap_min_step = -8.0, tap_max_step = 8.0, tap_current_step = Float64(n))),
    (name = "a CGMES ratio changer", mode = :ratio, n = -4, imp = false, net = () -> _tap_gate_net(), build = (net, n, z) -> stored!(net.branchVec[3], st -> S.SETapChanger(name = "_rtc", kind = :ratio, role = :ratio, side = :from, du = 0.0125, neutral = 10, low = 1, high = 21, step = st), rtc, n, z, 10)),
    (name = "b nameplate Delta-u psi 60", mode = :pst, n = 4, imp = false, net = () -> _tap_gate_net(), build = (net, n, z) -> applyTapNameplate!(net.branchVec[3]; phase_du_step = 0.01, phase_min_step = -10.0, phase_max_step = 10.0, phase_current_step = Float64(n), psi_deg = 60.0)),
    (name = "b CGMES asymmetrical 60 with x curve", mode = :pst, n = 4, imp = true, net = () -> _tap_gate_net(), build = (net, n, z) -> stored!(net.branchVec[3], asymch(asym60, 60.0), asym60, n, z, 0)),
    (name = "c nameplate Delta-u psi 90", mode = :pst, n = -3, imp = false, net = () -> _tap_gate_net(), build = (net, n, z) -> applyTapNameplate!(net.branchVec[3]; phase_du_step = 0.01, phase_min_step = -10.0, phase_max_step = 10.0, phase_current_step = Float64(n), psi_deg = 90.0)),
    (name = "c CGMES quadrature with x curve", mode = :pst, n = -3, imp = true, net = () -> _tap_gate_net(), build = (net, n, z) -> stored!(net.branchVec[3], asymch(asym90, 90.0), asym90, n, z, 0)),
    (name = "d CGMES symmetrical with x curve", mode = :pst, n = 5, imp = true, net = () -> _tap_gate_net(), build = (net, n, z) -> stored!(net.branchVec[3], st -> S.SETapChanger(name = "_ps", kind = :symmetrical, role = :phase, side = :from, du = 0.02, neutral = 0, low = -8, high = 8, step = st, x_curve = curve(sym), impedance = :x_curve), sym, n, z, 0)),
    (name = "e CGMES phase table with r and x", mode = :pst, n = 3, imp = true, net = () -> _tap_gate_net(), build = (net, n, z) -> stored!(net.branchVec[3], tabch, tab, n, z, 0)),
    (name = "f nameplate degree form", mode = :pst, n = -5, imp = false, net = () -> _tap_gate_net(), build = (net, n, z) -> applyTapNameplate!(net.branchVec[3]; phase_step_deg = 1.5, phase_min_step = -8.0, phase_max_step = 8.0, phase_current_step = Float64(n), psi_deg = 30.0)),
    (name = "f CGMES linear with x curve", mode = :pst, n = 6, imp = true, net = () -> _tap_gate_net(), build = (net, n, z) -> stored!(net.branchVec[3], st -> S.SETapChanger(name = "_pl", kind = :linear, role = :phase, side = :from, dphi_deg = 1.2, neutral = 0, low = -10, high = 10, step = st, x_curve = curve(lin), impedance = :x_curve), lin, n, z, 0)),
    # typed asymmetrical model under model.tap_changer_model = impedance_correction:
    # the power flow scales r and x with |1 + f e^(j psi)|^2 at every step
    (name = "b typed asymmetrical, impedance correction", mode = :pst, n = 7, imp = true, net = () -> _tap_gate_net(phase_taps = typed_asym(), model = :impedance_correction), build = function (net, n, z)
      br = net.branchVec[3]
      r0, x0 = br.r_pu, br.x_pu
      br.tap_winding.phase_taps.step = n
      Sparlectra.resolve_branch_taps!(br, br.tap_winding; tap_changer_model = :impedance_correction)
      z || (br.r_pu = r0; br.x_pu = x0)
      return br
    end),
  ]
  return cases
end

function test_tap_estimation_changer_types()::Bool
  # Tap estimation per changer type: every kind a to f, from a model at the
  # start step 0 to the true step with exact measurements on both sides:
  # the estimate comes back to 1e-6 steps, the fixation lands on the true
  # step, J below 1e-9 before and after. Run with impedance_correction auto
  # and off (there is no on, auto takes its role); the truth of the off run
  # keeps the impedance of the start step. A case with
  # per-step impedance data runs off once more against the truth with the
  # step impedance: the step still comes back, J is no longer zero. The
  # nameplate fraction case is a regression guard (the fraction cascade of
  # 0.10.0 recovered it as well).
  @testset "Tap estimation per changer type" begin (function ()
    for c in _tap_changer_cases()
      for imp in (:auto, :off)
        truth = c.net()
        c.build(truth, c.n, imp === :auto)
        meas = _tap_truth_meas(truth)
        model = c.net()
        c.build(model, 0, true)
        rel = setTapEstimation!(model; trafo = 3, mode = c.mode)
        @test rel.s0 == [0.0]
        res = _tap_se_run(model, meas; imp = imp)
        row = only(res.tapEstimates)
        ok = res.converged && abs(row.s_est - c.n) <= 1e-6 && row.step_fixed == row.neutral_step + c.n && row.deviation_steps == c.n && res.tapFixation.j_before < 1e-9 && res.tapFixation.j_after < 1e-9 && row.frozen_reason === :none
        @test ok
        ok || println("      tap changer case failed: ", c.name, " (", imp, "): s ", row.s_est, ", J ", res.tapFixation.j_before, " / ", res.tapFixation.j_after)
      end
      c.imp || continue
      truth = c.net()
      c.build(truth, c.n, true)
      meas = _tap_truth_meas(truth)
      model = c.net()
      c.build(model, 0, true)
      setTapEstimation!(model; trafo = 3, mode = c.mode)
      # the impedance left out of the model shows as a J rise of the
      # fixation; a rise above the 95 % quantile of chi2(1) warns and sets
      # the flag (seen on two of these cases), and the run warns about
      # nothing else
      logger = Test.TestLogger(min_level = Logging.Warn)
      res = Logging.with_logger(() -> _tap_se_run(model, meas; imp = :off), logger)
      warned = [r.message for r in logger.logs if r.level >= Logging.Warn]
      row = only(res.tapEstimates)
      ok = res.converged && row.step_fixed == row.neutral_step + c.n && res.tapFixation.j_after > 1e-9
      ok &= all(m -> occursin("above the 95 % quantile 3.841; one step state", m), warned) && length(warned) == (res.tapFixation.j_rise_exceeds_quantile ? 1 : 0)
      @test ok
      ok || println("      tap changer case failed: ", c.name, " (off against the step impedance): step ", row.step_fixed, ", J ", res.tapFixation.j_after)
    end
    # update_taps on a typed model: with auto the branch gets r and x
    # of the fixed step, with off it keeps r and x of the run (the resolver
    # of the typed winding rewrites them, the write-back puts them back)
    ct = only(cc for cc in _tap_changer_cases() if cc.name == "b typed asymmetrical, impedance correction")
    for imp in (:auto, :off)
      truth = ct.net()
      ct.build(truth, ct.n, imp === :auto)
      meas = _tap_truth_meas(truth)
      model = ct.net()
      ct.build(model, 0, true)
      r0, x0 = model.branchVec[3].r_pu, model.branchVec[3].x_pu
      setTapEstimation!(model; trafo = 3, mode = ct.mode)
      res = _tap_se_run(model, meas; imp = imp, update_taps = true)
      br = model.branchVec[3]
      rx = imp === :auto ? (truth.branchVec[3].r_pu, truth.branchVec[3].x_pu) : (r0, x0)
      ok = res.converged && res.tapFixation.fixed && br.tap_winding.phase_taps.step == ct.n && abs(br.r_pu - rx[1]) <= 1e-12 && abs(br.x_pu - rx[2]) <= 1e-12 && (imp === :off || abs(rx[1] - r0) > 1e-6)
      @test ok
      ok || println("      typed write-back (", imp, "): step ", br.tap_winding.phase_taps.step, ", r ", br.r_pu, " x ", br.x_pu, " expected ", rx)
    end
  end)() end
  return true
end

function test_tap_estimation_jacobian()::Bool
  # The tap column of the FD Jacobian (the chain rule: dh/dt by forward
  # difference along the analytic tangent dt/ds) against a central
  # difference of the exact prediction over s, with
  # impedance_correction auto (the column includes the impedance part), at
  # an off-step s, for an asymmetrical changer with a reactance curve (b)
  # and a phase table with r and x (e).
  @testset "Tap estimation Jacobian column" begin (function ()
    S = Sparlectra
    cases = _tap_changer_cases()
    for (name, soff) in (("b CGMES asymmetrical 60 with x curve", 1.37), ("e CGMES phase table with r and x", 0.61))
      c = only(cc for cc in cases if cc.name == name)
      truth = c.net()
      c.build(truth, c.n, true)
      meas = _tap_truth_meas(truth)
      net = c.net()
      c.build(net, 0, true)
      runpf!(net, 40, 1e-12, 0; method = :rectangular)
      setTapEstimation!(net; trafo = 3, mode = c.mode)
      nbus = length(net.nodeVec)
      slack = S._find_slack_idx(net)
      Ybus = createYBUS(net = net)
      map, _, _ = S._build_tap_state_map(net, (nbus - 1) + nbus + 1, :auto)
      @test map !== nothing && map.trafos[1].impedance
      S._tap_unstamp!(Ybus, net, map)
      x = S._initial_state_vector(net, slack; flatstart = false, tapSInit = [soff])
      col = map.trafos[1].cols[1]
      H, _, _ = S._measurement_jacobian_fd(meas, net, x, slack, nbus, Ybus; eps = 1e-6, tapMap = map)
      hat(sv) = (xv = copy(x); xv[col] = sv; S._predict_measurements(meas, net, S._state_to_voltage(xv, slack, nbus), Ybus; tapOverlay = S._tap_overlay(xv, map, net)))
      fd = (hat(soff + 1e-4) .- hat(soff - 1e-4)) ./ 2e-4
      # the column is a forward difference of step jac_eps / du = 1e-4
      # steps: measured 9.9e-5 relative for b (as large as the exact forward
      # difference of that step, the curvature of h, not the tangent) and
      # 2.2e-6 for e
      err = norm(Vector(H[:, col]) .- fd) / norm(fd)
      @test err < 1e-3
      err < 1e-3 || println("      Jacobian column ", name, ": relative error ", err)
    end
  end)() end
  return true
end

function test_tap_estimation_t2_line()::Bool
  # An asymmetrical changer of 1 % at psi 60 on a radial transformer with
  # t_n 1.0320366, true step 10, Vm on both sides and no angle measurement.
  # The factor of the tapped winding is the straight line 1 + s du e^(j psi):
  # at step 10 the branch angle is 4.715004 deg; the r1/r2 cascade of 0.10.0
  # put the same |t| on a circle at 5.841 deg.
  @testset "Tap estimation T2 line" begin (function ()
    S = Sparlectra
    C = S.CGMESImporter
    t2 = _tap_cim("_t2", :PhaseTapChangerAsymmetrical, [:neutralStep => "0", :lowStep => "-16", :highStep => "16", :voltageStepIncrement => "1", :windingConnectionAngle => "60"])
    store = _tap_cim_store([t2])
    tn = 1.0320366
    function t2net(n::Int)
      net = Net(name = "t2", baseMVA = 100.0)
      addBus!(net = net, busName = "H", vn_kV = 110.0)
      addBus!(net = net, busName = "L", vn_kV = 110.0)
      addProsumer!(net = net, busName = "H", type = "EXTERNALNETWORKINJECTION", referencePri = "H", vm_pu = 1.02, va_deg = 0.0)
      addProsumer!(net = net, busName = "L", type = "ENERGYCONSUMER", p = 40.0, q = 12.0)
      addPIModelTrafo!(net = net, fromBus = "H", toBus = "L", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = tn, shift_deg = 0.0, status = 1)
      br = net.branchVec[1]
      br.se_tap_changers = S.SETapChangerSet(ratio_n = tn, changers = [S.SETapChanger(name = "T2", kind = :asymmetrical, role = :phase, side = :from, du = 0.01, alpha_deg = 60.0, neutral = 0, low = -16, high = 16, step = n)], origin = :cgmes)
      # the importer: t_n times the factor of the PhaseTapChangerAsymmetrical
      m, a = C._phaseTapRatioShift_atstep(store, t2, n)
      br.tap_ratio = tn * m
      br.phase_shift_deg = a
      return net
    end
    meas = _tap_truth_meas(t2net(10))
    net = t2net(0)
    setTapEstimation!(net; trafo = 1, mode = :pst)
    res = _tap_se_run(net, meas)
    row = only(res.tapEstimates)
    @test res.converged && row.frozen_reason === :none
    @test row.step_fixed == 10 && abs(row.s_est - 10.0) < 1e-6
    @test abs(row.angle_est_deg - 4.715004) < 1e-6
    @test isapprox(row.angle_est_deg, C._phaseTapRatioShift_atstep(store, t2, 10)[2]; atol = 1e-6)
    # a tree edge without an angle anchor, but the magnitude of the factor
    # moves at psi 60: released, and at step 10 no vertex notice
    @test row.notice == ""
  end)() end
  return true
end

function test_tap_estimation_microgrid_be()::Bool
  # MicroGrid BE: the phase shifter BE-TR2_1
  # (PhaseTapChangerAsymmetrical, 1.25 %, psi 90, steps 1 to 25, neutral 13,
  # SvTapStep 16) on the changer set the import stores. Truth: the same
  # delivery with its SvTapStep at 2, re-imported; exact measurements of
  # that state on the normal import (step 16), BE-TR2_1 released alone: the
  # estimate comes back to step 2 on the absolute band 1 to 25. The
  # relative reading of 0.10.0 gave the band 4 to 28 (the SV step minus
  # neutral moved it by 3 steps) and clamped this truth at 4 (spec probe,
  # not rerun). Every stored set gives the ratio and shift of its branch at
  # the import steps, and every transformer branch without a changer
  # carries an empty set. Needs the MicroGrid package in the CGMES cache:
  # RAN or SKIPPED is printed.
  @testset "Tap estimation MicroGrid BE" begin (function ()
    S = Sparlectra
    C = S.CGMESImporter
    zip = joinpath(C.cgmesTestSetCacheDir(), "cases", "cgmes_microgrid_be.zip")
    if !isfile(zip)
      println("      tap estimation MicroGrid BE: SKIPPED (", zip, " not in the CGMES cache)")
      return
    end
    println("      tap estimation MicroGrid BE: RAN")
    store = C.loadCGMES(zip)
    ptc = only(tc for tc in C.objectsOf(store, :PhaseTapChangerAsymmetrical) if C.str(tc, :name) == "BE-TR2_1")
    model = redirect_stdout(() -> importCGMES(path = deepcopy(store), name = "be"), devnull)
    k = only(i for (i, b) in enumerate(model.net.branchVec) if b.se_tap_changers !== nothing && any(c -> c.mrid == ptc.mrid, b.se_tap_changers.changers))
    ch = only(model.net.branchVec[k].se_tap_changers.changers)
    @test (ch.kind, ch.side, ch.du, ch.alpha_deg, ch.neutral, ch.low, ch.high, ch.step) == (:asymmetrical, :from, 0.0125, 90.0, 13, 1, 25, 16.0)
    @test S._tap_band(ch) == (-12.0, 12.0)
    trafos = [b for b in model.net.branchVec if b.ratio != 0.0]
    @test all(b -> b.se_tap_changers !== nothing && b.se_tap_changers.origin === :cgmes, trafos)
    @test count(b -> isempty(b.se_tap_changers.changers), trafos) == 7
    ok = true
    for b in trafos
      set = b.se_tap_changers
      R, P = S._tap_set_ratio_shift(set, Float64[c.step - c.neutral for c in set.changers])
      ok &= isapprox(R, b.tap_ratio; rtol = 1e-12) && abs(P - b.phase_shift_deg) <= 1e-12
    end
    @test ok
    truth_store = deepcopy(store)
    sv = only(o for o in C.objectsOf(truth_store, :SvTapStep) if get(o.refs, :TapChanger, "") == ptc.mrid)
    @test C.num(sv, :position) == 16.0
    sv.attrs[:position] = "2"
    truth = redirect_stdout(() -> importCGMES(path = truth_store, name = "be_step2"), devnull)
    meas = _tap_truth_meas(truth.net)
    net = deepcopy(model.net)
    @test setTapEstimation!(net; trafo = k, mode = :pst).s0 == [3.0]
    res = _tap_se_run(net, meas)
    row = only(res.tapEstimates)
    @test res.converged && row.frozen_reason === :none
    @test abs(row.s_est + 11.0) < 1e-6 && row.step_fixed == 2.0 && row.deviation_steps == -14.0
    @test res.tapFixation.fixed && res.tapFixation.j_after < 1e-9

    # strongly coupled steps: br9 and br12 both leave bus 3 (br12 into the
    # star point of a three-winding transformer). Truth: the four regulators
    # (br8 phase, br9, br10, br12 ratio) 2 steps above their SV steps, noisy
    # measurements of seed 240, no Va. The continuous estimate s9 -0.744,
    # s12 2.530 (J 50.1) rounded at once gives (-1, 3) and J 3761, a rise far
    # above the 95 % quantile 9.49 of chi2(4); the truth (-1, 2) has J 58.17.
    # The rise gives a warning and the flag, the steps are fixed one at a
    # time (the one closest to its physical step first, the rest estimated
    # again) and the smaller J is kept: every step at the truth, J equal to
    # a plain run with the true taps (before the sequential fixation the
    # rounded result stayed, with an info line only)
    cfg = S.resolve_config(S.DEFAULT_SPARLECTRA_CONFIG_PATH, zip).config
    mg = redirect_stdout(() -> S.import_case(zip, cfg; run_kind = :state_estimation), devnull)
    empty!(mg.net.measurements)
    regs = [8, 9, 10, 12]
    tn = deepcopy(mg.net)
    for k in regs
      br = tn.branchVec[k]
      set = S._effective_tap_changers(br)
      s0 = S._tap_s_from_branch(set, br)
      S._tap_apply_position!(br, set, round.(s0) .+ 2.0; s_from = s0, impedance = mg.config.state_estimation.tap_estimation.impedance_correction)
    end
    _, erg = redirect_stdout(() -> runpf!(tn; config = S._copy_powerflow_with(mg.config.powerflow; tol = 1e-10, max_iter = 40)), devnull)
    @test erg == 0
    meas = generateMeasurementsFromPF(tn; noise = true, stddev = measurementStdDevs(vm = 0.002, pinj = 0.5, qinj = 0.5, pflow = 0.5, qflow = 0.5), relativeSigma = false, rng = Random.MersenneTwister(240))
    net = deepcopy(mg.net)
    for k in regs
      setTapEstimation!(net; trafo = k, mode = S._tap_release_mode(net.branchVec[k]))
    end
    seRun(n) = S.with_sparlectra_config(mg.config) do
      with_state_estimation_config(() -> runse!(n, meas); max_iter = 50, tol = 1e-6, update_net = false)
    end
    res = @test_logs (:warn, r"above the 95 % quantile 9\.488; sequential fixation") match_mode = :any seRun(net)
    tf = res.tapFixation
    ref = seRun(deepcopy(tn))
    @test res.converged && tf.j_rise_exceeds_quantile && tf.sequential && tf.j_simultaneous > 3000.0
    @test all(t -> t.fixed && t.deviation_steps == 2.0, res.tapEstimates)
    @test abs(tf.j_after - ref.objectiveJ) < 1e-6 && tf.dof_after == ref.dof
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

    _two_island_tap_net() = _se_two_island_tap_net(step)

    # branch layout per island block: 1/2 lines, 3/4 trafos (a), 5/6 lines,
    # 7/8 trafos (b). Truth: +2 steps on trafo 3 (island a), -1 on trafo 7.
    tnet = _two_island_tap_net()
    for (k, n) in ((3, 2.0), (7, -1.0))
      br = tnet.branchVec[k]
      # the fraction grid of the constructor: ratio / (1 + n step)
      br.tap_ratio = br.ratio / (1.0 + n * step)
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
    @test rowA.step_fixed == 2
    @test rowB.step_fixed == -1
    @test rowA.fixed && rowB.fixed
    # every row field travels through the island merge (the field list
    # there used to drop the model steps)
    @test rowA.step_model == 0 && rowA.deviation_steps == 2 && rowA.notice == ""
    @test res.tapFixation.iterations_estimation >= res.tapFixation.iterations_fixation >= 1
    tf = res.tapFixation
    @test tf.fixed
    @test tf.dof_after == tf.dof_before + 2   # one freed state per island
    @test tf.j_after < 1e-10

    # since 0.60.0, with update_net and update_taps the island path mirrors
    # the fixed steps onto the parent net and recomputes the parent flows
    # and losses with them. The parent first carries the flows of a power
    # flow on the model taps; after the run its flows equal the truth.
    netW = _two_island_tap_net()
    ite, erg = runpf!(netW, 40, 1e-12, 0; islands_enabled = true)
    @test erg == 0
    calcNetLosses!(netW)
    setTapEstimation!(netW; trafo = 3, mode = :ratio)
    setTapEstimation!(netW; trafo = 7, mode = :ratio)
    resW = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_net = true, update_taps = true) do
      runse!(netW, meas)
    end
    @test resW.converged && resW.islands !== nothing && length(resW.islands) == 2
    @test _se_flow_deviation(netW, tnet) < 1e-6

    # numerical isolation: releasing ONLY island a's trafo must reproduce
    # island a's voltages bitwise (island b's model discrepancy stays in b)
    netA = _two_island_tap_net()
    setTapEstimation!(netA; trafo = 3, mode = :ratio)
    resA = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_net = false) do
      runse!(netA, meas)
    end
    @test resA.converged
    @test length(resA.tapEstimates) == 1
    @test resA.tapEstimates[1].step_fixed == 2
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

    function _truth_net(r1)
      tnet = _tap_gate_net()
      br = tnet.branchVec[3]
      # the fraction grid of the constructor: ratio / (1 + n step)
      br.tap_ratio = br.ratio / (1.0 + r1)
      ite, erg = runpf!(tnet, 40, 1e-12, 0; method = :rectangular)
      @test erg == 0
      return tnet
    end
    function _truth_meas(r1; noise = false, rng = Random.MersenneTwister(1), grossQ = 0.0, includeIa = false)
      tnet = _truth_net(r1)
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
    @test res.converged && res.tapEstimates[1].step_fixed == 2
    @test net.branchVec[3].tap_ratio === ratio0            # bitwise protection
    res = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_taps = true) do
      runse!(net, meas)
    end
    @test res.converged
    expect = net.branchVec[3].ratio / (1.0 + 2 * step)
    @test isapprox(net.branchVec[3].tap_ratio, expect; atol = 1e-14)
    # since 0.60.0 the flows and losses on the net use the written-back tap,
    # so they equal the truth power flow (before, calcNetLosses! ran with
    # the old tap and the voltages of the fixed-tap solve)
    @test _se_flow_deviation(net, _truth_net(2 * step)) < 1e-6
    net.branchVec[3].tap_est_mode = :none                  # no release: the written-back model
    res2 = with_state_estimation_config(() -> runse!(net, meas); max_iter = 40, tol = 1e-10)     # explains the set by itself
    @test res2.converged
    @test res2.objectiveJ < 1e-10

    # 2) PMU synergy: on a TWO-REGULATOR release (declared nameplate pair:
    # fraction ratio changer and Delta-u changer at psi 90) the phase
    # regulator lives in the angle information; an accurate current PHASOR
    # pair on the released trafo shrinks its estimate spread over noisy
    # realizations (a magnitude-only set leaves the angle soft). Until
    # 0.60.0 the truth sat on the degree grid of an undeclared branch,
    # which a :both release now refuses.
    function _truth_both()
      tn = _tap_gate_pair_net(2, -1)
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
        n = _tap_gate_pair_net(0, 0)
        withPMU && addCurrentPhasorMeasurement!(n; i_A = abs(IendB), ia_deg = rad2deg(angle(IendB)), branchNr = 3, direction = :from, sigmaI = 0.5, sigmaIa = 0.02)
        setTapEstimation!(n; trafo = 3, mode = :both, alpha_deg = 90.0)
        append!(n.measurements, m)
        r = with_state_estimation_config(max_iter = 40, tol = 1e-10) do
          runse!(n)
        end
        @test r.converged
        push!(out, only(t for t in r.tapEstimates if t.role === :phase).s_est)
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
    @test abs(res.tapEstimates[1].s_est - 2.0) < 0.5

    # 4) machine-trafo back-calculation: mutually exclusive with a release,
    # needs a registered SE, and recovers the GSU position from the AVR
    # setpoint plus the dispatch against the estimated network voltage
    gnet = _tap_gate_net()
    addBus!(net = gnet, busName = "G1", vn_kV = 10.0)
    addProsumer!(net = gnet, busName = "G1", type = "SYNCHRONOUSMACHINE", p = 20.0, q = 5.0, vm_pu = 1.02)
    addPIModelTrafo!(net = gnet, fromBus = "H2", toBus = "G1", r_pu = 0.001, x_pu = 0.09, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
    gnet.branchVec[5].has_ratio_tap = true
    gnet.branchVec[5].tap_step = step
    # TRUE machine tap: 2 steps on the fraction grid; the machine holds
    # 1.02 pu at G1
    brg = gnet.branchVec[5]
    brg.tap_ratio = brg.ratio / (1.0 + 2 * step)
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
    @test abs(bt.electrical_step - 2.0) < 0.25 && bt.s_est == bt.electrical_step
    @test bt.machine_bus == 5
    @test bt.q_residual_mvar < 0.1

    # 5) se_view lists the release with its mode and one entry per released
    # regulator (declared Delta-u changer at psi 30 next to the ratio changer)
    vnet = _tap_gate_pair_net(0, 0; psi = 30.0)
    setTapEstimation!(vnet; trafo = 3, mode = :both, alpha_deg = 30.0)
    view = se_view(vnet)
    @test length(view.taps) == 1
    @test view.taps[1].branch == 3 && view.taps[1].mode == :both
    @test [g.kind for g in view.taps[1].regulators] == [:ratio, :asymmetrical] && view.taps[1].regulators[2].alpha_deg == 30.0
    io = IOBuffer()
    print_se_view(io, view)
    out = String(take!(io))
    @test occursin("Tap estimation releases: 1", out)
    @test occursin("mode both; ratio (a) du 0.00625 side to; asymmetrical (b) alpha 30.0 du 0.01 side to", out)
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
of the options, before any case is imported. The relative sigmas of
`generateMeasurementsFromPF` are checked here too: per row
sigma = max(f*|value|, floor), Va stays absolute, a dominating floor wins.
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
      Sparlectra.readMeasurementTable!(net; file = out, zero_injection = false)
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
    # passive_as_zi not stated: state_estimation.zero_injection of the case
    # decides at generation (default on, the former literal true)
    @test okopts.critical_count == 0 && okopts.flow_ends === :both && okopts.passive_as_zi === nothing
    # the key is case scope: set in the case configuration file (what the
    # Web UI checkbox writes), it switches the generator's ZI rows off
    mktempdir() do d
      copy_case = joinpath(d, basename(case))
      cp(case, copy_case)
      for want in (true, false)
        want || Sparlectra.write_case_config(copy_case, Dict{String,Any}("state_estimation.zero_injection" => "off"))
        out = joinpath(d, "zi_$(want).csv")
        SparlectraApp._se_generate_measurement_set(copy_case, out, okopts)
        @test any(l -> occursin("ZI_PINJ_bus_", l), readlines(out)) == want
      end
    end
    # a noise-free set fits the model it was generated from. sp_case188 has
    # a coupler between Falkwik_380 (160 MW generation) and the passive
    # section Falkwik_380b: the injection rows of the section with the
    # generator carried the Y-bus value (103.8 MW, the link flow left out)
    # while the passive section became a zero-injection row, so the
    # cluster sum missed 56.2 MW and J was 3859 (dof 1054) instead of 0.
    # With the distributed slack on, Falkwik_380 delivers its schedule plus
    # its share (156.1 MW); the device sum took the schedule alone, missed
    # the cluster's Y-bus total, fell back to the Y-bus values with a
    # warning and J was 3757 (dof 1054)
    for ds in (false, true)
      mktempdir() do d
        c188 = joinpath(dirname(@__DIR__), "data", "scf", "sp_case188.scf.json")
        if ds
          cp(c188, joinpath(d, basename(c188)))
          c188 = joinpath(d, basename(c188))
          Sparlectra.write_case_config(c188, Dict{String,Any}("power_flow.tol" => 1.0e-8, "power_flow.distributed_slack.enabled" => true))
        end
        set = joinpath(d, "nf188.csv")
        gen_log = Test.TestLogger(min_level = Logging.Warn)
        Logging.with_logger(gen_log) do
          SparlectraApp._se_generate_measurement_set(c188, set, okopts)
        end
        @test !any(r -> occursin("do not add up", string(r.message)), gen_log.logs)
        res = redirect_stdout(devnull) do
          SparlectraApp._run_state_estimation_service(c188, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, joinpath(d, "run"), "nf188", set; flatstart = true)
        end
        dres = SparlectraApp.to_dict(res)
        @test dres["status"] == "succeeded"
        jm = match(r"converged in \d+ iteration\(s\); J = ([0-9.eE+-]+), dof", read(joinpath(d, "run", "run.log"), String))
        @test jm !== nothing && parse(Float64, jm.captures[1]) < 1e-6
      end
    end
    # relative sigmas (percent of reading): per row sigma = max(f*|value|,
    # floor), floors bind on tiny readings, Va stays absolute (moved from
    # the removed measurement CSV v1 testset)
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
  end)() end
  return true
end

# The net of the tap changer set test of the CGMES import (A4b in
# test/test_cgmes_importer.jl): a chain A..H of seven transformers with a
# typed asymmetrical phase changer (step 16 of 1..25, neutral 13, 1.25 %,
# psi 90), a symmetrical one with an xMin/xMax curve, a phase table
# (-3..3 at `tab_step`), a RatioTapChanger, one transformer without changer,
# a fixed shift (the export's single-step PhaseTapChangerLinear) and a
# transformer tapped on both ends; the ENI at A makes A-B a machine
# transformer. Written by writeCGMESFiles and zipped into `dir`.
function _generator_a4b_delivery(dir::AbstractString, name::AbstractString; tab_step::Int)::String
  net = Net(name = name, baseMVA = 100.0)
  for b in ("A", "B", "C", "D", "E", "F", "G", "H")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
  addProsumer!(net = net, busName = "H", type = "ENERGYCONSUMER", p = 10.0, q = 2.0)
  asym() = PhaseTapChangerModel(kind = :asymmetrical, step = 16, lowStep = 1, highStep = 25, neutralStep = 13, voltage_step_increment = 0.0125, winding_connection_angle_deg = 90.0, convention = :direct_regulating_vector)
  sym = PhaseTapChangerModel(kind = :symmetrical, step = 3, lowStep = -8, highStep = 8, neutralStep = 0, voltage_step_increment = 0.02, x_min = 0.05, x_max = 0.06, convention = :direct_regulating_vector)
  tab = PhaseTapChangerModel(kind = :tabular, step = tab_step, neutralStep = 0, convention = :direct_regulating_vector, table = [TapTablePoint(step = k, ratio = 1.0 + 0.001 * k^2, angle_deg = 2.0 * k) for k in -3:3])
  rtc = PowerTransformerTaps(Vn_kV = 110.0, step = 2, lowStep = -5, highStep = 5, neutralStep = 0, voltageIncrement_kV = 1.1)
  tr(f, t; shift = 0.0, phase_taps = nothing, taps = nothing) = addPIModelTrafo!(net = net, fromBus = f, toBus = t, r_pu = 0.002, x_pu = 0.05, b_pu = 0.0, status = 1, ratio = 1.0, shift_deg = shift, phase_taps = phase_taps, taps = taps)
  tr("A", "B"; phase_taps = asym())
  tr("B", "C"; phase_taps = sym)
  tr("C", "D"; phase_taps = tab)
  tr("D", "E"; taps = rtc)
  tr("E", "F")
  tr("F", "G"; shift = -3.0)
  tr("G", "H"; phase_taps = asym())
  net.trafos[end].side2.taps = PowerTransformerTaps(Vn_kV = 110.0, step = -1, lowStep = -4, highStep = 4, neutralStep = 0, voltageIncrement_kV = 1.65)
  return _generator_zip_delivery(net, dir, name)
end

# writeCGMESFiles of `net` into dir/name, zipped to dir/name.zip
function _generator_zip_delivery(net::Net, dir::AbstractString, name::AbstractString)::String
  files = joinpath(dir, name)
  mkpath(files)
  writeCGMESFiles(net; path = files)
  z = joinpath(dir, name * ".zip")
  ZipArchives.ZipWriter(z) do w
    for file in sort(readdir(files; join = true))
      ZipArchives.zip_newfile(w, basename(file))
      write(w, read(file))
    end
  end
  return z
end

"""
The measurement generator moves the tap changer the tap estimation
estimates: on a CGMES delivery the deviation goes onto the
stored changer set of each drawn transformer (the first ratio changer, else
the first phase changer), its ratio, angle and, under the default
`impedance_correction = auto`, the reactance of the x curve follow the
changer data, and the default run with every tap released recovers it: the
phase regulators fix two steps off their model step and J after the fixation
is zero (noise-free set). Under `off` the truth keeps the start impedance,
the set differs from the `auto` one and the `off` run recovers it as well.
The taps table of the set lists the transformers with a changer only, in
absolute steps. A table changer has no position off its points: a
deviation there is refused. A delivery whose transformers have no changer
(the last pool) gets no deviation and a note.
"""
function test_measurement_generator_tap_changer_sets()::Bool
  @testset "Measurement generator tap deviation on the changer set" begin (function ()
    opts(; kw...) = SparlectraApp.MeasurementGeneratorOptions(; noise = false, gross_k = 0.0, include_i = true, sigma_u_pct = 0.5, sigma_i_pct = 1.0, sigma_p_pct = 1.0, sigma_q_pct = 1.0, sigma_ia_deg = 0.05, kw...)
    cfg = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH
    tap_rows(path) = [split(l[3:end], ",") for l in readlines(path) if occursin(r"^# \d+,", l)]
    data_rows(path) = [l for l in readlines(path) if !startswith(l, "#")]
    mktempdir() do d
      # the phase table sits at step 2 of -3..3: two steps up is no point
      z2 = _generator_a4b_delivery(d, "a4b_tab2"; tab_step = 2)
      @test_throws ArgumentError SparlectraApp._se_generate_measurement_set(z2, joinpath(d, "refused.csv"), opts(; tap_steps = 2.0, tap_count = 7))
      # table at step 1: every drawn changer has room for two steps
      z = _generator_a4b_delivery(d, "a4b"; tab_step = 1)
      runs = Dict{Symbol,String}()
      for imp in (:auto, :off)
        case = z
        if imp === :off
          case = joinpath(d, "a4b_off.zip")
          cp(z, case)
          Sparlectra.write_case_config(case, Dict{String,Any}("state_estimation.tap_estimation.impedance_correction" => "off"))
        end
        set = joinpath(d, "a4b_$(imp).csv")
        runs[imp] = set
        g = SparlectraApp._se_generate_measurement_set(case, set, opts(; tap_steps = 2.0, tap_count = 7))
        # A-B is a machine transformer (ENI at A), E-F and F-G have no
        # changer: four transformers are drawn, the deviation of the branch
        # tapped on both ends goes onto its ratio changer
        @test occursin("phase deviation 2.0 step(s) (0.02 pu each, symmetrical) on PST branch 2", g.tap_note)
        @test occursin("phase deviation 2.0 table step(s) on PST branch 3", g.tap_note)
        @test occursin("tap deviation 2.0 steps (~1.96% ratio) on transformer branch 4", g.tap_note)
        @test occursin("on transformer branch 7", g.tap_note)
        @test occursin("limited to 4 eligible transformer(s), max 7 requested", g.tap_note)
        # taps table: no row for the changerless E-F and F-G; absolute
        # steps after the deviation (A-B stays at 16, the RTC of G-H moves
        # from -1 to 1), step column du, neutral ratio of the set (0.985 on
        # G-H, the export's end-2 correction)
        rows = tap_rows(set)
        @test [(r[1], r[6], r[7], r[8], r[9], r[11]) for r in rows] == [("1", "1.0", "0.0125", "16.0", "16", "0.0"), ("2", "1.0", "0.02", "5.0", "5", "2.0"), ("3", "1.0", "0.03790658503988655", "3.0", "3", "2.0"), ("4", "1.0", "0.01", "4.0", "4", "2.0"), ("7", "0.985", "0.015", "1.0", "1", "2.0")]
        res = redirect_stdout(devnull) do
          SparlectraApp._run_state_estimation_service(case, cfg, joinpath(d, "run_$(imp)"), "a5b_$(imp)", set; tap_estimation = true, max_iter = 50)
        end
        dd = SparlectraApp.to_dict(res)
        @test dd["status"] == "succeeded"
        md = dd["metadata"]
        @test md["se_tap_impedance_correction"] == String(imp)
        @test md["se_tap_j_after"] < 1e-9
        regs = Dict((t["branch"], t["regulator"]) => t for t in md["se_tap_estimates"] if get(t, "source", "") == "estimated")
        # the phase regulators of B-C and C-D: role phase, fixed two steps
        # off the model step; the ratio regulators of D-E and G-H as well,
        # the phase changer of G-H stays at its model step
        for (b, role, dev) in ((2, "phase", 2), (3, "phase", 2), (4, "ratio", 2), (7, "phase", 0))
          r = regs[(b, 1)]
          @test (r["role"], r["deviation_steps"], r["fixed"], r["frozen_reason"]) == (role, dev, true, "none")
        end
        @test (regs[(7, 2)]["role"], regs[(7, 2)]["deviation_steps"]) == ("ratio", 2)
      end
      # the x curve of B-C moves the reactance of the auto truth only
      @test data_rows(runs[:auto]) != data_rows(runs[:off])

      # a delivery whose transformers have no changer: only the last pool,
      # each drawn transformer is skipped with a note, the set has no
      # deviation and its taps block says why
      net = Net(name = "plain", baseMVA = 100.0)
      for b in ("A", "B", "C")
        addBus!(net = net, busName = b, vn_kV = 110.0)
      end
      addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", referencePri = "A", vm_pu = 1.0, va_deg = 0.0)
      addProsumer!(net = net, busName = "C", type = "ENERGYCONSUMER", p = 10.0, q = 2.0)
      addPIModelTrafo!(net = net, fromBus = "A", toBus = "B", r_pu = 0.002, x_pu = 0.05, b_pu = 0.0, status = 1, ratio = 1.0, shift_deg = 0.0)
      addPIModelTrafo!(net = net, fromBus = "B", toBus = "C", r_pu = 0.002, x_pu = 0.05, b_pu = 0.0, status = 1, ratio = 1.0, shift_deg = -3.0)
      zp = _generator_zip_delivery(net, d, "plain")
      setp = joinpath(d, "plain.csv")
      gp = SparlectraApp._se_generate_measurement_set(zp, setp, opts(; tap_steps = 2.0, tap_count = 2))
      @test occursin("transformer branch 1 (A-B) has no tap changer, no deviation", gp.tap_note)
      @test occursin("transformer branch 2 (B-C) has no tap changer, no deviation", gp.tap_note)
      @test isempty(tap_rows(setp))
      @test any(==("# case has no tap changers (no tap steps to estimate)"), readlines(setp))
    end
  end)() end
  return true
end

function run_state_estimation_tests()
  # Aggregates all state-estimation unit tests to keep coverage explicit and ordered.

"""
State estimation on a DTF network case, on the shipped demo deck. The SE
service takes `case_format` like the power-flow service does, and DTF is
an accepted SE format. Since the format detection judges a deck by its
content, the run needs no named format either: the same request without
`case_format` gives the same result. (Before that a `.DAT` file without a
named format was refused, and the SE service, which once detected the
format itself, rejected every DTF case before it looked at the
measurements.)
"""
function test_state_estimation_dtf_service()
  dtf = joinpath(dirname(@__DIR__), "data", "dtf_demo", "sp_dtf5.DAT")
  cfg = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH
  net = Sparlectra.DTFImporter.build_net(Sparlectra.DTFImporter.read_dtf(dtf))
  runpf!(net, 40, 1e-8, 0; method = :rectangular)
  append!(net.measurements, generateMeasurementsFromPF(net; noise = false))
  mktempdir() do d
    mf = joinpath(d, "dtf.measurements.csv")
    writeMeasurementTable(mf, Sparlectra._measurement_table_from_measurements(net, net.measurements, nameMeasurementKeys(net)).rows)

    ok = redirect_stdout(devnull) do
      SparlectraApp._run_state_estimation_service(dtf, cfg, joinpath(d, "run_ok"), "dtf_se", mf; case_format = :dtf_for001)
    end
    d_ok = SparlectraApp.to_dict(ok)
    @test d_ok["status"] == "succeeded"
    @test d_ok["metadata"]["run_mode"] == "se"

    # without the format the deck is recognised by its content: the same
    # run, the same objective
    auto = redirect_stdout(devnull) do
      SparlectraApp._run_state_estimation_service(dtf, cfg, joinpath(d, "run_auto"), "dtf_se_auto", mf)
    end
    d_auto = SparlectraApp.to_dict(auto)
    @test d_auto["status"] == "succeeded"
    @test d_auto["metadata"]["se_objective"] == d_ok["metadata"]["se_objective"]

    # a .DAT file that is no deck fails with a message that names the file
    # and the reason
    report = joinpath(d, "report.DAT")
    write(report, "plain unsupported data\n")
    refused = redirect_stdout(devnull) do
      SparlectraApp._run_state_estimation_service(report, cfg, joinpath(d, "run_refused"), "dtf_se_refused", mf)
    end
    d_refused = SparlectraApp.to_dict(refused)
    @test d_refused["status"] == "failed"
    @test occursin("report.DAT is not a DTF network deck", d_refused["message"])
  end
  # the import entry point accepts the format and rejects an unsupported one
  @test SparlectraApp._se_import_case(dtf, Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true); requested_format = :dtf_for001).net isa Sparlectra.Net
  return true
end



# Both MiniGrid regression sets below run through one service call: the case
# zip is rebuilt from the ENTSO-E conformity cache (bus-branch base case v3
# with its boundary, packed as `cgmes_minigrid.zip`, the name the sets are
# bound to; never downloaded here), the set goes through the SE service with
# tap estimation, and the expected warnings of these sets (zero-injection
# rows at the two isolated buses of the base case) are captured. A warning
# matching one of `expected` is captured too, and at least one such warning
# per pattern must come. Everything else is re-emitted for the runner's
# warning check. `f` receives the result dict, run.log and
# se_diagnostics.md; without the cache the run is reported as SKIPPED and
# `f` is not called.
function _minigrid_fallback_run(f::Function, set::AbstractString, name::AbstractString; staged::Bool = false, expected::Vector{String} = String[])
  root = joinpath(Sparlectra.CGMESImporter.cgmesTestSetCacheDir(), "extracted", "MiniGrid", "BusBranch")
  base = joinpath(root, "CGMES_v2.4.15_MiniGridTestConfiguration_BaseCase_v3")
  boundary = joinpath(root, "CGMES_v2.4.15_MiniGridTestConfiguration_Boundary_v3")
  if !(isdir(base) && isdir(boundary))
    println("      SE tap fallback diagnostics (MiniGrid, $(basename(set))): SKIPPED (conformity package not in the local cache)")
    return
  end
  println("      SE tap fallback diagnostics (MiniGrid, $(basename(set))): RAN")
  mktempdir() do d
    z = joinpath(d, "cgmes_minigrid.zip")
    ZipArchives.ZipWriter(z) do w
      for file in vcat(sort(readdir(base; join = true)), sort(readdir(boundary; join = true)))
        ZipArchives.zip_newfile(w, basename(file))
        write(w, read(file))
      end
    end
    out = joinpath(d, "run")
    logger = Test.TestLogger(min_level = Logging.Warn)
    res = Logging.with_logger(logger) do
      redirect_stdout(devnull) do
        SparlectraApp._run_state_estimation_service(z, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, out, name, set; tap_estimation = true, max_iter = 50, tap_staged_start = staged)
      end
    end
    for p in expected
      @test any(record -> occursin(p, record.message), logger.logs)
    end
    for record in logger.logs
      occursin("excluded (its bus is isolated)", record.message) && continue
      any(p -> occursin(p, record.message), expected) && continue
      @warn record.message
    end
    f(SparlectraApp.to_dict(res), read(joinpath(out, "run.log"), String), read(joinpath(out, "se_diagnostics.md"), String))
  end
  return
end

# The reported block of se_diagnostics.md describes the frozen-tap run the
# status line reports: same run id, case, set, iterations, J and dof as
# run.log and the metadata. `case` is the case name the block shows.
function _assert_tap_fallback_report(dd, log::AbstractString, diag::AbstractString, name::AbstractString, set::AbstractString; case::AbstractString = "cgmes_minigrid.zip")
  @test dd["status"] == "succeeded"
  @test get(dd["metadata"], "se_tap_estimation_fallback", false) == true
  @test occursin("repeating WITHOUT tap estimation", log)
  jm = match(r"converged in (\d+) iteration\(s\); J = ([0-9.eE+-]+), dof = (\d+)", log)
  @test jm !== nothing
  jm === nothing && return
  @test occursin("Tap estimation fallback", diag)
  @test occursin("## Reported result", diag)
  @test occursin("**Converged:** true in $(jm.captures[1]) iteration(s)", diag)
  @test occursin("**Objective J:** $(jm.captures[2]) (dof=$(jm.captures[3])), every tap frozen at its model position", diag)
  @test occursin("- **Run:** $(name)", diag)
  @test occursin("- **Case:** $(case)", diag)
  @test occursin("- **Measurement set:** $(basename(set))", diag)
  @test occursin("- **J/dof:** $(round(parse(Float64, jm.captures[2]) / parse(Int, jm.captures[3]); digits = 2))", diag)
  @test occursin("## Diagnostics pass (every row active, before any elimination)", diag)
  # the diagnostics pass printed below the block describes the same frozen-tap run
  @test !occursin("**Converged:** false", diag)
  @test isapprox(dd["metadata"]["se_objective"], parse(Float64, jm.captures[2]); atol = 1e-6)   # run.log rounds to six digits
  @test dd["metadata"]["se_dof"] == parse(Int, jm.captures[3])
  # the verdict file counts what the run did: every row a diagnostics pass
  # removed (the pass on the released taps and the repeat on the frozen
  # ones) is eliminated, inactive only where the table says valid=false,
  # and its suspicious and eliminated measurements are the counts of the
  # metadata (until the fix the repeat's own trace replaced the first
  # pass, which then showed as inactive, and se_suspicious counted the
  # pass before the elimination)
  vrows = [split(l, ",") for l in readlines(joinpath(dd["output_dir"], "measurements_verdict.csv"))[2:end]]
  tab = Dict(r.line => r for r in Sparlectra.readMeasurementTable(joinpath(dd["output_dir"], "measurements.csv")))
  @test all(v -> occursin("inactive", v[end]) == !tab[parse(Int, v[1])].valid, vrows)
  @test length(unique(v[4] for v in vrows if v[9] == "true")) == dd["metadata"]["se_eliminations"]
  @test length(unique(v[4] for v in vrows if occursin("suspicious", v[end]))) == dd["metadata"]["se_suspicious"]
  return
end

"""
The diagnostics artifact describes the run the status line reports. On the
CGMES MiniGrid with the fixture set `cgmes_minigrid_tap_fallback` (noisy,
with zero-injection constraints; its truth is an adopted SE state of an
earlier run, not a fresh solve, which is why the KCL precheck flags its
passive buses: regression data, not a plausible measurement campaign) the
released-tap attempt does not converge and the frozen-tap repeat does. The
artifact used to describe the failed attempt (the diverged J at the dof
with tap states, the frozen run's J on the page); now the diagnostics
repeat on the frozen taps and the artifact opens with the reported
numbers. This set confirms
all three of its eliminations and stops at `max_eliminations`. On the set
`cgmes_minigrid_tap_rejected` the staged estimate was rejected
(`worse_fit`) before the sequential fixation; it is now kept
with the fixation flags. The rejection itself is checked on an off-grid
sp_case14 in `test_state_estimation_tap_fallback`.
Needs the conformity package in the local cache: RAN or SKIPPED is printed.
"""
function test_state_estimation_tap_fallback_diagnostics_consistent()::Bool
  @testset "Tap fallback: diagnostics artifact matches the reported run (MiniGrid)" begin
    # a set built from a run outside the band carries the generator's warning
    @test occursin("outside the band", SparlectraApp._se_truth_run_warning(Dict{String,Any}("se_band_reason" => "high")))
    @test occursin("frozen at their model positions", SparlectraApp._se_truth_run_warning(Dict{String,Any}("se_band_reason" => "ok", "se_tap_estimation_fallback" => true)))
    @test SparlectraApp._se_truth_run_warning(Dict{String,Any}("se_band_reason" => "ok", "numerical_status" => "converged")) == ""
    set = joinpath(@__DIR__, "fixtures", "measurements", "cgmes_minigrid_tap_fallback.measurements.csv")
    # the COLD mass release (no staged start) does not converge here and
    # falls back; the fixture pins that path's report
    _minigrid_fallback_run(set, "minigrid_fb") do dd, log, diag
      _assert_tap_fallback_report(dd, log, diag, "minigrid_fb", set)
      @test occursin("- **Eliminated rows:** 3 ", diag)
      @test occursin("(stop: max_eliminations)", diag)
      @test !occursin("Reverted:", diag)
      @test !occursin("Not eliminated", diag)
    end
    # the staged start on this set keeps its estimate: the import stores
    # the changers of the delivery (0.60.0), so only the two 380 kV legs of
    # the three-winding transformers with a RatioTapChanger are released
    # (1 % steps, 1 to 33, neutral 17), no longer the eight transformer
    # branches of the constructor grid; their fixed steps give J 10752.1 at
    # dof 75 against 46692.2 at dof 78 with the model positions (with the
    # eight releases the estimate was rejected, J 48298.6 at dof 75)
    # The fixation of each solve here (the estimate and the three
    # eliminations) raises J by about 9 against the 95 % quantile 5.99 of
    # chi2(2) on this set (J/dof 140 to 620), so each warns
    _minigrid_fallback_run(set, "minigrid_fb_staged"; staged = true, expected = ["above the 95 % quantile 5.991"]) do dd, log, _
      @test dd["status"] == "succeeded"
      @test get(dd["metadata"], "se_tap_staged_start", false) == true
      @test get(dd["metadata"], "se_tap_estimation_fallback", false) == false
      @test occursin("staged start, one solve with the 2 released tap(s) frozen", log)
    end
    # the set `cgmes_minigrid_tap_rejected` (generated on this case by the SE
    # measurement generator: fresh solve, gaussian noise of seed 165, sigma
    # U 0.5 %, P and Q 1 %, no tap deviation and no gross error, so the
    # model steps are the truth): the two released legs leave bus 1 into
    # the star points of two three-winding transformers and are strongly
    # coupled. Their estimates sit at -0.498 and -0.503 steps; rounded at
    # once to 0 and -1 they gave J 68770 at dof 82 against 73.44 with the
    # model positions, and the staged estimate was rejected (worse_fit)
    # before the sequential fixation. Now the rise above the 95 % quantile of
    # chi2(2) warns and the steps are fixed one at a time: -0.503 first to
    # -1, the other leg follows to -1. The pair (16, 16) one step below the
    # truth gives J 73.43, the data barely see the common step of the two
    # legs; the estimate is kept with the flags in the metadata and the
    # notice in run.log. Seeds 1 to 150 of the same generator options give
    # no rejection either; the worse_fit fallback is checked on an off-grid
    # sp_case14 in test_state_estimation_tap_fallback
    rejected = joinpath(@__DIR__, "fixtures", "measurements", "cgmes_minigrid_tap_rejected.measurements.csv")
    _minigrid_fallback_run(rejected, "minigrid_rej_staged"; staged = true, expected = ["above the 95 % quantile 5.991; sequential fixation"]) do dd, log, _
      @test dd["status"] == "succeeded"
      @test get(dd["metadata"], "se_tap_staged_start", false) == true
      @test get(dd["metadata"], "se_tap_estimation_fallback", false) == false
      @test get(dd["metadata"], "se_tap_j_rise_exceeds_quantile", false) == true && get(dd["metadata"], "se_tap_sequential_fixation", false) == true
      @test occursin("kept the sequential result", log)
      legs = [t for t in dd["metadata"]["se_tap_estimates"] if get(t, "source", "") == "estimated"]
      @test length(legs) == 2 && all(t -> t["fixed"] == true && abs(t["deviation_steps"]) <= 1, legs) && legs[1]["step_fixed"] == legs[2]["step_fixed"]
      # J of the kept pair against the solve at the model steps of the staged start
      jm = match(r"staged start, one solve with the 2 released tap\(s\) frozen \(\d+ iteration\(s\), J = ([0-9.eE+-]+)", log)
      @test jm !== nothing && dd["metadata"]["se_tap_j_after"] <= parse(Float64, jm.captures[1]) * (1.0 + 1e-6)
    end
  end
  return true
end

"""
An elimination counts only when its control solve converges. The fixture
set `cgmes_minigrid_elimination_reverted` (same case, same truth source as
the set above, a later generation) ends its third elimination in a
non-converged control solve: that removal is reverted (the row is active
again), the artifact names it as "Reverted" in the pass and as "Not
eliminated" in the reported block, the sequence stops with
`not_converged`, and only the two confirmed eliminations leave the
reported run. Needs the conformity package in the local cache: RAN or
SKIPPED is printed.
"""
function test_state_estimation_elimination_reverted()::Bool
  @testset "Tap fallback: an unconfirmed elimination is reverted (MiniGrid)" begin
    set = joinpath(@__DIR__, "fixtures", "measurements", "cgmes_minigrid_elimination_reverted.measurements.csv")
    _minigrid_fallback_run(set, "minigrid_rev") do dd, log, diag
      _assert_tap_fallback_report(dd, log, diag, "minigrid_rev", set)
      @test occursin("- **Eliminated rows:** 2 ", diag)
      @test occursin("- **Not eliminated:** `", diag)
      @test occursin("(stop: not_converged)", diag)
      @test occursin("Reverted:", diag)
    end
    # the staged start (the service default) converges on the same set, and
    # the fixed taps fit better than the model positions (measured J 19177
    # at dof 75 against 39952 at dof 78), so the estimate is kept. With the
    # changers the import stores (0.60.0) only the two 380 kV legs with a
    # RatioTapChanger are released; with the eight transformer branches of
    # the constructor grid released the kept estimate had J 34251. Until
    # a5121638 this set was rejected (J 46180 against 39952): the 380 kV
    # legs fixed on the constructor grid of 0.625 % instead of their
    # delivered 1 % steps
    _minigrid_fallback_run(set, "minigrid_rev_staged"; staged = true) do dd, log, _
      @test dd["status"] == "succeeded"
      @test get(dd["metadata"], "se_tap_staged_start", false) == true
      @test get(dd["metadata"], "se_tap_estimation_fallback", false) == false
      @test occursin("tap estimation: staged start", log)
      @test !occursin("repeating WITHOUT tap estimation", log)
    end
  end
  return true
end


"""
Tap-estimation fallback: released taps are extra states, and a measurement
set that estimates the voltages cleanly can still be far too thin to pin
them. The run then does not settle at all and the user gets nothing,
although the SAME set works without the taps (seen on
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
      # the diagnostics artifact states the reported numbers first, the
      # same J and dof run.log carries (the diagnostics pass's own J can
      # differ when it eliminated rows)
      diag_ok = read(joinpath(d, "run_ok", "se_diagnostics.md"), String)
      @test occursin("## Reported result", diag_ok)
      jm = match(r"converged in \d+ iteration\(s\); J = ([0-9.eE+-]+), dof = (\d+)", log_ok)
      @test jm !== nothing
      @test occursin("**Objective J:** $(jm.captures[1]) (dof=$(jm.captures[2]))", diag_ok)
      @test occursin("**Converged:** true in", diag_ok)
    end
    # The step limit. A released regulator state is
    # bounded in how far ONE iteration may move it, at a quarter of the
    # changer's declared mechanical travel. Without it the Gauss-Newton step
    # drives r1 toward -1, where the cascade
    # t = t_base/((1+r1)(1+r2 e^{j alpha})) is singular: measured
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
        tmap, _, _ = Sparlectra._build_tap_state_map(tnet, 1, :auto)
        @test tmap !== nothing
        limits = Sparlectra._tap_step_limits(tmap)
        @test !isempty(limits)
        br = tnet.branchVec[trafo]
        # the step band of the nameplate row in steps of the fraction grid:
        # tap_min = ratio / (1 + pos_max step), tap_max = ratio / (1 + pos_min step)
        full = round((br.ratio / br.tap_min - 1.0) / br.tap_step) - round((br.ratio / br.tap_max - 1.0) / br.tap_step)
        # a quarter of the full travel in steps (until 0.60.0 in units of
        # the fraction r1), and strictly positive so a released regulator
        # can still move
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
    # those model positions. Under a flat start the staged start (taps
    # frozen for one solve, then released from that state) rescues the
    # shipped cases at a cap of 4, where the cold tap solve used to fail
    # (sp_case60 and sp_case188 needed the fallback before; case1354pegase
    # with 240 released taps needed 63 iterations cold and 4 staged). The
    # fallback is still reached from a configured WARM start, which skips
    # the staging: sp_case188 at a cap of 2 converges plain but not with
    # the taps, so the repeat succeeds (measured). The cap was 3 until
    # 0.50.0, while the generator's injection rows at the Falkwik coupler
    # did not fit the model; with a fitting set the taps converge at 3
    # (measured: no fallback at 3, fallback not_converged at 2).
    mktempdir() do d
      shipped = joinpath(dirname(@__DIR__), "data", "scf", "sp_case188.scf.json")
      set = joinpath(d, "gen188.csv")
      SparlectraApp._se_generate_measurement_set(shipped, set, SparlectraApp.MeasurementGeneratorOptions(; noise = true, gross_k = 0.0, tap_steps = 0.0,
        include_i = false, sigma_u_pct = 0.5, sigma_i_pct = 1.0, sigma_p_pct = 1.0, sigma_q_pct = 1.0,
        sigma_ia_deg = 0.0, seed = 42))
      # the staged start itself: flat start, cap 4, the mass release of the
      # 13 ratio changers converges WITHOUT the fallback. Since 0.60.0 the
      # degree rows of the SCF tap changers (1.25 deg, -24..24, alpha 0) are
      # declared changers and join the release (13 transformers, 26 states),
      # pinned on the shipped case; the ratio-only release the shipped case
      # had until 0.60.0 runs on a copy without the degree rows (measured: no
      # fallback at cap 4). With the degree rows the staged estimate was
      # rejected (worse_fit) while the mesh-edge degree regulators were
      # rounded although only relative to frozen ones and the steps were
      # rounded at once (J 1043.0 at dof 1052
      # against 1001.6 at dof 1057 with the model positions). Now five
      # regulators stay states, the rise above the 95 % quantile of chi2(13)
      # warns, the steps fixed one at a time give J 996.3, and the estimate
      # is kept
      ratio_only = joinpath(d, "sp_case188.scf.json")
      write(ratio_only, replace(read(shipped, String), r",\s*\{\s*\"alpha_deg\": 0\.0,\s*\"index\": 2,[^}]*\}" => ""))
      for (case, tag, rise) in ((ratio_only, "st", false), (shipped, "stdeg", true))
        logger = Test.TestLogger(min_level = Logging.Warn)
        res_st = Logging.with_logger(logger) do
          redirect_stdout(devnull) do
            SparlectraApp._run_state_estimation_service(case, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH,
              joinpath(d, tag), "tapst188", set; tap_estimation = true, max_iter = 4, flatstart = true)
          end
        end
        rises = [r.message for r in logger.logs if r.level >= Logging.Warn && occursin("above the 95 % quantile", r.message)]
        @test isempty(rises) == !rise && all(m -> occursin("kept the sequential result", m), rises)
        # every other warning goes to the capture of the test set
        for r in logger.logs
          (r.level >= Logging.Warn && !(r.message in rises)) && @warn r.message
        end
        dst = SparlectraApp.to_dict(res_st)
        @test dst["status"] == "succeeded"
        @test get(dst["metadata"], "se_tap_staged_start", false) == true
        @test get(dst["metadata"], "se_tap_estimation_fallback", false) == false
        @test get(dst["metadata"], "se_tap_j_rise_exceeds_quantile", false) == rise && get(dst["metadata"], "se_tap_sequential_fixation", false) == rise
      end
      res = redirect_stdout(devnull) do
        SparlectraApp._run_state_estimation_service(shipped, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH,
          joinpath(d, "fb"), "tapfb188", set; tap_estimation = true, max_iter = 2, flatstart = false)
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
      # (as seen): the method comes from the run kind
      @test SparlectraApp._powerflow_run_index_solver(res) == "wls"
    end

    # The acceptance rule of the staged start (worse_fit): the fixed steps
    # are kept only if they fit at least as well as the model positions.
    # Trigger: sp_case14 with the ratio changer of branch 31 at the model
    # position 4.4, between two physical steps, and exact measurements of
    # that state. The estimate comes back to 4.4, the fixation can only reach
    # step 4: J 0.233 against 4.4e-16 at the model positions (dof 83 both),
    # so the estimate is rejected and the run falls back with its own reason
    # and note on every surface. The rise stays below the 95 % quantile 5.99
    # of chi2(2), so no warning. A noisy set of the same state was kept
    # (seed 42: J 77.94 fixed against 78.09 at the model positions), and the
    # MiniGrid set cgmes_minigrid_tap_rejected no longer triggers the rule
    # since the sequential fixation
    mktempdir() do d
      offgrid = joinpath(d, "sp_case14_offgrid.scf.json")
      txt = read(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"), String)
      write(offgrid, replace(txt, r"(\"branch\": 31,[^\]]*?\"pos\": )4," => s"\g<1>4.4,"))
      @test read(offgrid, String) != txt
      set = joinpath(d, "gen14.csv")
      SparlectraApp._se_generate_measurement_set(offgrid, set, SparlectraApp.MeasurementGeneratorOptions(; noise = false, gross_k = 0.0, tap_steps = 0.0,
        include_i = false, sigma_u_pct = 0.5, sigma_i_pct = 1.0, sigma_p_pct = 1.0, sigma_q_pct = 1.0,
        sigma_ia_deg = 0.0, seed = 42))
      res = redirect_stdout(devnull) do
        SparlectraApp._run_state_estimation_service(offgrid, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH,
          joinpath(d, "rej"), "offgrid_rej", set; tap_estimation = true, max_iter = 50)
      end
      dd = SparlectraApp.to_dict(res)
      log = read(joinpath(d, "rej", "run.log"), String)
      diag = read(joinpath(d, "rej", "se_diagnostics.md"), String)
      _assert_tap_fallback_report(dd, log, diag, "offgrid_rej", set; case = basename(offgrid))
      @test get(dd["metadata"], "se_tap_staged_start", false) == true
      @test get(dd["metadata"], "se_tap_estimation_fallback_reason", "") == "worse_fit"
      @test occursin("tap estimation rejected: the 2 estimated tap(s)", log)
      @test occursin(SparlectraApp._SE_TAP_REJECTED_NOTE, diag)
      @test occursin(SparlectraApp._SE_TAP_REJECTED_NOTE, SparlectraApp._webui_se_summary(dd))
      @test occursin(SparlectraApp._webui_escape(SparlectraApp._SE_TAP_REJECTED_NOTE), SparlectraApp._webui_se_tap_section(dd))
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
      ("Zero injection as equality constraint", test_state_estimation_zero_injection_constraint),
      ("Bad-data diagnostics", test_state_estimation_bad_data_diagnostics),
      ("PMU angle measurements", test_state_estimation_pmu_va_measurements),
      ("Current-magnitude measurements (ImagMeas)", test_state_estimation_imag_measurements),
      ("Sequential elimination and wii", test_state_estimation_sequential_elimination),
      ("Shunt estimation (case A)", test_state_estimation_shunt_estimation),
      ("Shunt back-calculation (case B)", test_state_estimation_shunt_back_calculation),
      ("Link contraction", test_state_estimation_link_contraction),
      ("se_view", test_state_estimation_se_view),
      ("Robust R modification", test_state_estimation_robust_r),
      ("Wilson-Hilferty band test", test_state_estimation_wilson_hilferty),
      ("Diagnostics unification", test_state_estimation_diagnostics_unification),
      ("SE-to-PF chain", test_state_estimation_se_chain),
      ("Island-wise SE, parallel against serial", test_state_estimation_islands_parallel),
      ("Current-angle measurements (IaMeas)", test_state_estimation_ia_measurements),
      ("Tap changer data: s -> t(s)", test_tap_changer_data),
      ("Tap estimation per changer type", test_tap_estimation_changer_types),
      ("Tap estimation Jacobian column", test_tap_estimation_jacobian),
      ("Tap estimation T2 line", test_tap_estimation_t2_line),
      ("Tap estimation MicroGrid BE", test_tap_estimation_microgrid_be),
      ("Tap estimation across islands", test_state_estimation_tap_islands),
      ("Tap completion (write-back, PMU, machine trafo)", test_state_estimation_tap_completion),
      ("DTF case through the SE service", test_state_estimation_dtf_service),
      ("Tap fallback on non-convergence", test_state_estimation_tap_fallback),
      ("Tap fallback: diagnostics artifact matches the reported run (MiniGrid)", test_state_estimation_tap_fallback_diagnostics_consistent),
      ("Tap fallback: an unconfirmed elimination is reverted (MiniGrid)", test_state_estimation_elimination_reverted),
      ("Measurement generator critical thinning", test_measurement_generator_critical_thinning),
      ("Measurement generator tap deviation on the changer set", test_measurement_generator_tap_changer_sets),
    ]

    for (name, testfn) in tests
      @testset "$name" begin (function ()
        @test _se_run_quiet(testfn) == true
      end)() end
    end
  end)() end
end

# The island-wise SE on worker threads against the serial loop, alone: the
# threads profile (tools/run_gates.sh threads) repeats it with four threads.
function run_state_estimation_parallel_tests()
  @testset "State estimation, parallel against serial" begin (function ()
    @test _se_run_quiet(test_state_estimation_islands_parallel) == true
  end)() end
end
