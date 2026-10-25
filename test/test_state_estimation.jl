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
    # resolution: the same rows rounded to the step in table units (U in kV
    # through the bus Vn, P in MW), each within half a step of the unrounded
    # value, the sigma unchanged, other types untouched
    res = Dict(Sparlectra.VmMeas => 0.1, Sparlectra.PflowMeas => 0.5)
    mr = generateMeasurementsFromPF(net; noise = false, stddev = std, resolution = res)
    @test length(mr) == length(meas)
    for (a, b) in zip(meas, mr)
      @test a.id == b.id && a.sigma == b.sigma
      step = get(res, a.typ, 0.0)
      scale = a.typ == Sparlectra.VmMeas ? net.nodeVec[a.busIdx].comp.cVN : 1.0
      k = b.value * scale / step
      @test step == 0.0 ? a.value == b.value : (abs(k - round(k)) < 1e-9 && abs(b.value - a.value) * scale <= step / 2 + 1e-9)
    end

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

## The start ladder (state_estimation.start_ladder, 0.60.0): a WLS from the
## flat start that does not converge is solved once more from the voltages
## stored in the net. sp_case14 with its shipped set (measured): 4
## iterations from the flat start, 3 from the stored state of the case, so a
## cap of 3 fails flat and converges on the second try. The second try
## restores the stored voltages the failed try overwrote (update_net) and
## is the run with flatstart = false bit for bit; a net without stored
## voltages gets no second try. Also the abort point of the rank test of the
## tap columns, the loop a Web UI abort of case6495rte could not stop.
function test_state_estimation_start_ladder()::Bool
  @testset "State estimation start ladder" begin (function ()
    data = joinpath(dirname(@__DIR__), "data", "scf")
    net0 = Sparlectra.import_case(joinpath(data, "sp_case14.scf.json"), Sparlectra.SparlectraConfig()).net
    readMeasurementTable!(net0; file = joinpath(data, "sp_case14.measurements.csv"))
    @test Sparlectra._se_has_stored_state(net0)
    solve(net; kw...) = with_state_estimation_config(() -> runse!(net); merge((max_iter = 3,), kw)...)
    lad = solve(deepcopy(net0))
    @test lad.converged && lad.start === :stored_after_flat
    warm = solve(deepcopy(net0); flatstart = false)
    @test warm.converged && warm.start === :stored
    @test lad.voltages == warm.voltages && lad.objectiveJ == warm.objectiveJ && lad.iterations == warm.iterations
    off = solve(deepcopy(net0); start_ladder = false)
    @test !off.converged && off.start === :flat
    flat = solve(deepcopy(net0); max_iter = 50)
    @test flat.converged && flat.start === :flat
    # no stored voltages: the second try would start where the first did
    bare = deepcopy(net0)
    for nd in bare.nodeVec
      nd._vm_pu = 1.0
      nd._va_deg = 0.0
    end
    @test !Sparlectra._se_has_stored_state(bare)
    nb = solve(bare)
    @test !nb.converged && nb.start === :flat

    # a set abort token stops the rank test of the tap columns
    H = SparseArrays.sparse([1.0 0.0 1.0; 0.0 1.0 1.0; 1.0 1.0 0.0; 0.5 0.0 1.0])
    @test Sparlectra._tap_rank_kept(H, [1, 2], [3]) == [true]
    tok = Threads.Atomic{Bool}(true)
    err = fetch(@async begin
      Sparlectra.sparlectra_arm_abort_token!(tok)
      try
        Sparlectra._tap_rank_kept(H, [1, 2], [3])
      catch e
        e
      finally
        Sparlectra.sparlectra_disarm_abort_token!()
      end
    end)
    @test err isa Sparlectra.PowerFlowAborted
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
      # SEResult.conditioning: the gain matrix under on, the scaled KKT
      # matrix under constraint, and the KKT matrix is not worse conditioned
      # than the gain matrix with rows at 1e-3 (1-norm estimates, measured
      # 2.6e4 against 6.9e9 on sp_case14, 2.4e7 against 2.3e11 on sp_case60)
      con_c = runs[:constraint].r.conditioning
      on_c = runs[:on].r.conditioning
      @test (on_c.zero_injection, on_c.matrix, con_c.zero_injection, con_c.matrix) == (:on, :gain, :constraint, :kkt)
      @test isfinite(con_c.condition_number) && con_c.condition_number <= on_c.condition_number
      # the multipliers: one entry per passive bus, the fused coupler bus of
      # sp_case60 once as its LINKAGG pair; none under on
      @test runs[:on].r.zeroInjectionMultipliers === nothing
      zm = runs[:constraint].r.zeroInjectionMultipliers
      @test length(zm) == length(passive) - (isempty(aggs) ? 0 : 1)
      @test count(z -> z.kind === :link_aggregate, zm) == (isempty(aggs) ? 0 : 1)
      @test all(z -> isfinite(z.lambda_p) && isfinite(z.lambda_q) && !isempty(z.bus), zm)
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

    # on with the sigma tending to 0 (1e-7 MW per row, the rows added under
    # state_estimation.zero_injection_sigma_mva) gives the state of
    # constraint to 1e-6 pu (measured 1.2e-7, after 43 iterations against
    # 6), and constraint holds the balance exactly: the injection at every
    # passive bus from the Ybus and the estimated voltages is below 1e-9 pu
    # (measured 4.2e-13). On sp_case60 the gain matrix at 1e-7 is singular
    # to working precision (estimate 1.9e31) and on does not converge.
    m7 = Measurement[m for m in m14 if !m.zero_injection]
    with_state_estimation_config(zero_injection_sigma_mva = 1e-7) do
      @test Sparlectra._add_missing_zero_injections!(m7, net14) == 2 * length(p14)
    end
    @test all(m -> m.sigma == 1e-7, (m for m in m7 if m.zero_injection))
    v7 = Dict{Symbol,Vector{ComplexF64}}()
    for mode in (:on, :constraint)
      r7 = with_state_estimation_config(zero_injection = mode, zero_injection_sigma_mva = 1e-7, max_iter = 100) do
        runse!(deepcopy(net14), m7)
      end
      @test r7.converged
      v7[mode] = r7.voltages
    end
    @test maximum(abs.(v7[:on] .- v7[:constraint])) < 1e-6
    s7 = v7[:constraint] .* conj.(createYBUS(net = net14) * v7[:constraint])
    @test maximum(abs.(s7[collect(p14)])) < 1e-9

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
    # the multiplier is the pull of the measurements against a balance:
    # without the gross error and with the passive buses as constraints
    # too, the wrong constraint at the load bus carries the largest one
    # (measured |lambda_q| 1.45 per Mvar against at most 0.77 elsewhere)
    mw = copy(meas)
    mw[k] = g
    pw = Set(findPassiveBuses(net))
    filter!(m -> !(m.typ in (Sparlectra.PinjMeas, Sparlectra.QinjMeas) && m.busIdx in pw), mw)
    Sparlectra._add_missing_zero_injections!(mw, net)
    rw = with_state_estimation_config(zero_injection = :constraint) do
      runse!(deepcopy(net), mw)
    end
    @test rw.converged
    zw = rw.zeroInjectionMultipliers
    @test length(zw) == length(pw) + 1
    @test zw[argmax([max(abs(z.lambda_p), abs(z.lambda_q)) for z in zw])].bus == getCompName(net.nodeVec[loadBus].comp)

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

    # Angle reference: the estimate keeps the slack at 0;
    # the warm start, the angle predictions, the result and the write-back
    # are in the model reference. The same net with the slack at -5 deg: the
    # power flow solution is the one above turned by -5 deg (no shifter).
    # Noise-free sets in the model reference (offset 0) and with a PMU time
    # base 5 deg ahead: every run ends on the power flow angles with the slack
    # at its stored -5 deg, alpha is the time base offset (0 or 5, not the
    # slack angle), J stays at the noise level under :off, and the warm start
    # from the converged state needs one iteration. Before the correction
    # alpha was -5, the net came back turned by +5 deg, J under :off was
    # 1.9e5 and the warm start took 6 iterations.
    tn = createTest3BusNet()
    slk = findfirst(nd -> Sparlectra.getNodeType(nd) == Sparlectra.Slack, tn.nodeVec)
    setVmVa!(node = tn.nodeVec[slk], vm_pu = something(tn.nodeVec[slk]._vm_pu, 1.0), va_deg = -5.0)
    _, ergT = runpf!(tn, 40, 1e-10, 0; method = :rectangular)
    vaPF = Float64[nd._va_deg for nd in tn.nodeVec]
    @test ergT == 0 && abs(vaPF[slk] + 5.0) < 1e-12 && maximum(abs.(vaPF .- (rad2deg.(angle.(Vref)) .- 5.0))) < 1e-8
    for (offDeg, mode, alpha) in ((0.0, :auto, 0.0), (5.0, :auto, 5.0), (0.0, :off, nothing)), fs in (true, false)
      measT = generateMeasurementsFromPF(tn; includeVa = true, vaRefOffsetDeg = offDeg, noise = false, stddev = std)
      n2 = deepcopy(tn)
      rT = with_state_estimation_config(max_iter = 20, tol = 1e-8, flatstart = fs, update_net = true, pmu_ref_offset = mode) do
        runse!(n2, measT)
      end
      vaSE = Float64[nd._va_deg for nd in n2.nodeVec]
      got = (rT.converged, rT.objectiveJ < 1e-8, vaSE[slk] == vaPF[slk], maximum(abs.(vaSE .- vaPF)) < 1e-6, maximum(abs.(rad2deg.(angle.(rT.voltages)) .- vaPF)) < 1e-6,
             alpha === nothing ? rT.vaRefOffsetDeg === nothing : isapprox(rT.vaRefOffsetDeg, alpha; atol = 1e-6), fs || offDeg != 0.0 || rT.iterations == 1)
      @test (offDeg, mode, fs, got) == (offDeg, mode, fs, (true, true, true, true, true, true, true))
    end
    # the same on a CGMES case with an SV angle reference: MicroGrid BE,
    # power flow, noise-free set with Va, warm start with update_net
    zipBE = joinpath(Sparlectra.CGMESImporter.cgmesTestSetCacheDir(), "cases", "cgmes_microgrid_be.zip")
    if isfile(zipBE)
      println("      SE angle reference MicroGrid BE: RAN")
      cfgBE = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, zipBE).config
      be = redirect_stdout(() -> Sparlectra.import_case(zipBE, cfgBE; run_kind = :state_estimation), devnull)
      empty!(be.net.measurements)
      _, ergBE = redirect_stdout(() -> runpf!(be.net; config = Sparlectra._copy_powerflow_with(be.config.powerflow; tol = 1e-10, max_iter = 40)), devnull)
      sBE = findfirst(nd -> Sparlectra.getNodeType(nd) == Sparlectra.Slack, be.net.nodeVec)
      vaBE = Float64[something(nd._va_deg, 0.0) for nd in be.net.nodeVec]
      @test ergBE == 0 && abs(vaBE[sBE]) > 0.1
      measBE = generateMeasurementsFromPF(be.net; includeVa = true, noise = false)
      nBE = deepcopy(be.net)
      rBE = Sparlectra.with_sparlectra_config(be.config) do
        with_state_estimation_config(() -> runse!(nBE, measBE); max_iter = 20, tol = 1e-6, flatstart = false, update_net = true)
      end
      vaSE = Float64[something(nd._va_deg, 0.0) for nd in nBE.nodeVec]
      @test (rBE.converged, rBE.iterations, nBE.nodeVec[sBE]._va_deg == vaBE[sBE], maximum(abs.(vaSE .- vaBE)) < 1e-6, isapprox(rBE.vaRefOffsetDeg, 0.0; atol = 1e-6)) == (true, 1, true, true, true)
    else
      println("      SE angle reference MicroGrid BE: SKIPPED (", zipBE, " not in the CGMES cache)")
    end

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

    # 4) range check of the voltage measurements on the
    # shipped sp_case14: a voltage measurement outside 0.5 to 1.5 pu at an
    # energized bus is replaced by a pseudo measurement (imported state,
    # else 1.0 pu, sigma times bad_data.substitute_sigma_factor), never
    # eliminated; at a bus the switching state does not energize, 0 kV is
    # right and the bus with its branches leaves the estimation
    (function ()
      std14 = measurementStdDevs(vm = 0.005, pinj = 0.5, qinj = 0.5, pflow = 0.5, qflow = 0.5)
      net14 = load_fixture_net("sp_case14")
      ite14, erg14 = runpf!(net14, 40, 1e-10, 0; method = :rectangular)
      @test erg14 == 0
      @test ite14 > 0
      meas14 = generateMeasurementsFromPF(net14; noise = true, stddev = std14, rng = MersenneTwister(7))
      vm_at(ms, b) = findfirst(m -> m.typ == Sparlectra.VmMeas && m.busIdx == b, ms)
      with_value(m, v) = Measurement(typ = m.typ, value = v, sigma = m.sigma, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id, linkIdx = m.linkIdx)
      # 0 kV at Ilmrode_110 (bus 9) with an imported state of 1.02 pu, a
      # non-finite value at Moorau_20 (bus 10), 1.6 pu at Hassel_20 (bus 11):
      # all three buses are energized
      k9, k10, k11 = vm_at(meas14, 9), vm_at(meas14, 10), vm_at(meas14, 11)
      bad = copy(meas14)
      bad[k9] = with_value(meas14[k9], 0.0)
      bad[k10] = with_value(meas14[k10], NaN)
      bad[k11] = with_value(meas14[k11], 1.6)
      net14._sv_vm_pu[9] = 1.02
      quiet_diag(n, ms) = with_state_estimation_config(() -> runse_diagnostics(n, ms); update_net = false)
      diag = @test_logs (:warn, r"Vm_bus_9 at bus Ilmrode_110 is 0\.0 pu.*replaced by 1\.02 pu \(imported state\).*\(range_check_substituted\)") (:warn, r"Vm_bus_10 .* NaN pu.*replaced by 1\.0 pu \(nominal voltage\)") (:warn, r"Vm_bus_11 .* 1\.6 pu.*range_check_substituted") match_mode = :any quiet_diag(net14, bad)
      @test (diag.range_check.substituted, diag.range_check.deenergized) == (3, 0)
      @test [(r.measurement_index, r.substitute, r.source) for r in diag.range_check.rows] == [(k9, 1.02, :imported_state), (k10, 1.0, :nominal), (k11, 1.0, :nominal)]
      # sigma times the default factor 3
      @test all(r.sigma_substitute == 3.0 * r.sigma for r in diag.range_check.rows)
      @test diag.final_diagnostics.converged
      @test diag.diagnostics.range_check.substituted == 3
      # every bus keeps its voltage measurement: in the ranking, flagged, not
      # eliminated, and the caller's vector is unchanged
      flagged = Dict(r.measurement_index => r.range_check for r in diag.final_diagnostics.measurement_ranking if r.measurement_index >= 1)
      @test [get(flagged, k, :missing) for k in (k9, k10, k11)] == fill(:range_check_substituted, 3)
      @test !any(e -> e.measurement_index in (k9, k10, k11), diag.eliminations)
      @test bad[k9].value == 0.0 && isnan(bad[k10].value)
      io = IOBuffer()
      print_se_diagnostics(io, diag; topN = 30)
      @test occursin("Range check of the voltage measurements: 3 substituted, 0 deenergized", String(take!(io)))
      # runse! runs the same check before its WLS
      res = @test_logs (:warn, r"runse!: range check: voltage measurement Vm_bus_9") match_mode = :any with_state_estimation_config(() -> runse!(net14, bad); update_net = false)
      @test res.converged
      # 0 kV at Torfloh_20 (bus 13) behind its open branch: no source, not
      # energized; the bus row, its injection and the flow of its branch
      # leave the estimation (deenergized), the rest converges
      dead = load_fixture_net("sp_case14")
      kbr = findfirst(br -> Int(br.fromBus) == 13 || Int(br.toBus) == 13, dead.branchVec)
      Sparlectra.setBranchStatus!(dead.branchVec[kbr], false)
      # without the isolation marking bus 13 is an island without a source
      @test !(13 in Sparlectra._se_energized_buses(dead))
      @test 12 in Sparlectra._se_energized_buses(dead)
      Sparlectra.markIsolatedBuses!(net = dead)
      @test runpf!(dead, 40, 1e-10, 0; method = :rectangular)[2] == 0
      dmeas = generateMeasurementsFromPF(dead; noise = true, stddev = std14, rng = MersenneTwister(7))
      push!(dmeas, Measurement(typ = Sparlectra.VmMeas, value = 0.0, sigma = 0.005, busIdx = 13, id = "Vm_bus_13"))
      push!(dmeas, Measurement(typ = Sparlectra.PinjMeas, value = 0.0, sigma = 0.5, busIdx = 13, id = "Pinj_bus_13"))
      push!(dmeas, Measurement(typ = Sparlectra.PflowMeas, value = 0.0, sigma = 0.5, branchIdx = kbr, direction = :from, id = "Pflow_open_branch"))
      n0 = length(dmeas)
      ddiag = @test_logs (:info, r"Vm_bus_13 at bus Torfloh_20 is 0\.0 pu.*not energized.*3 measurement\(s\) inactive \(deenergized\)") match_mode = :any quiet_diag(dead, dmeas)
      @test (ddiag.range_check.substituted, ddiag.range_check.deenergized) == (0, 1)
      @test only(ddiag.range_check.rows).deactivated == [n0 - 2, n0 - 1, n0]
      @test ddiag.final_diagnostics.converged
      @test !any(r -> r.measurement_index > n0 - 3, ddiag.final_diagnostics.measurement_ranking)
    end)()

    # 5) the same check in the SE service: run.log line, metadata counter and
    # the verdict reason, on the shipped sp_case5 with its U row at
    # Lindhof_110 set to 0 kV
    (function ()
      dir = test_scratch_path()
      mkpath(dir)
      src = joinpath(dirname(@__DIR__), "data", "scf")
      cp(joinpath(src, "sp_case5.scf.json"), joinpath(dir, "sp_case5.scf.json"))
      lines = readlines(joinpath(src, "sp_case5.measurements.csv"))
      lk = findfirst(l -> startswith(l, ",,,Lindhof_110,U,"), lines)
      lines[lk] = ",,,Lindhof_110,U,0.0,0.55,converted;v1,true"
      # the truth comment of a generated set, so the run writes se_deltas.csv
      # (truth: the original reading 105.7238 kV on 110 kV)
      push!(lines, "# truth_value,U@Lindhof_110,$(105.72380696737058 / 110.0)")
      write(joinpath(dir, "zero_kv.csv"), join(lines, "\n") * "\n")
      out = joinpath(dir, "run")
      service() = redirect_stdout(() -> SparlectraApp._run_state_estimation_service(joinpath(dir, "sp_case5.scf.json"), Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, out, "range_check", joinpath(dir, "zero_kv.csv")), devnull)
      r = @test_logs (:warn, r"range_check_substituted") match_mode = :any service()
      d = SparlectraApp.to_dict(r)
      @test d["status"] == "succeeded"
      @test (d["metadata"]["se_range_check_substituted"], d["metadata"]["se_range_check_deenergized"]) == (1, 0)
      @test occursin(r"range check: voltage measurement \S+ at bus Lindhof_110 is 0\.0 pu.*replaced by 1\.0 pu \(nominal voltage\).*\(range_check_substituted\)", read(joinpath(out, "run.log"), String))
      v = run_csv(joinpath(out, "measurements_verdict.csv"))
      hit = [row for row in v.rows if occursin("range_check_substituted", row[end])]
      @test [(row[2], row[3]) for row in hit] == [("Lindhof_110", "U")]
      # se_deltas.csv: the row keeps its raw reading as measured, but its
      # residual belongs to the substitute 1.0 pu solved with sigma times 3,
      # so estimated = 1.0 - residual (about 0.96 pu, not about 0) and
      # t_est = |1.0 - estimated| / (3 sigma)
      dl = run_csv(joinpath(out, "se_deltas.csv"))
      drow = only(row for row in dl.rows if row[2] == "U@Lindhof_110")
      @test dl.number(drow[7]) == 0.0
      est = dl.number(drow[9])
      @test 0.9 < est < 1.1
      @test isapprox(dl.number(drow[14]), abs(1.0 - est) / (3.0 * dl.number(drow[13])); rtol = 1e-5)
    end)()
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
  # whole, and link measurements are WLS rows (a single link, and a coupler
  # ring with circulation states).
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

    # --- 3) link measurement as a WLS row (0.60.0), the study case: a
    # single link, its P measurement about 100 MW off the balance flow at
    # sigma 2 MW. A capacitor at B1a puts a shunt term into the row. The
    # truth is the KCL at B1 from the power flow: the line S-B1 arrives at
    # B1, the load there is 10 MW + 3 Mvar, the link takes the rest.
    netL = create_se_link_net(shunt_buses = ["B1a"])
    runpf!(netL, 40, 1e-10, 0; method = :rectangular)
    calcNetLosses!(netL)
    b1L = geNetBusIdx(net = netL, busName = "B1")
    brSL = only(br for br in netL.branchVec if Int(br.toBus) == b1L)
    fP = -brSL.tBranchFlow.pFlow - 10.0
    fQ = -brSL.tBranchFlow.qFlow - 3.0
    measL = generateMeasurementsFromPF(netL; noise = false, stddev = measurementStdDevs(vm = 2e-3, pinj = 1.0, qinj = 1.0, pflow = 1.0, qflow = 1.0))
    resL0 = with_state_estimation_config(max_iter = 30, tol = 1e-8, update_net = false) do
      runse!(netL, measL)
    end
    # a consistent link row is one more row of the WLS and fits
    measC = copy(measL)
    addPflowMeasurement!(measC; net = netL, value = fP, sigma = 2.0, linkNr = 1, id = "P_link")
    resC = with_state_estimation_config(max_iter = 30, tol = 1e-8, update_net = false) do
      runse!(netL, measC)
    end
    @test resC.converged
    @test resC.dof == resL0.dof + 1
    @test abs(resC.residuals[end]) < 1e-6
    # 100 MW off: the row has a normalized residual and is localizable
    measB = copy(measL)
    addPflowMeasurement!(measB; net = netL, value = fP - 100.0, sigma = 2.0, linkNr = 1, id = "P_link")
    repB = with_state_estimation_config(max_iter = 30, tol = 1e-8) do
      validate_measurements(netL, measB)
    end
    topB = first(repB.measurement_ranking)
    @test topB.id == "P_link"
    @test topB.measurement_index == length(measB)
    @test topB.suspicious && topB.localizable
    @test abs(topB.residual) > 50.0
    # the bad-data loop eliminates it first
    diagB = with_state_estimation_config(max_iter = 30, tol = 1e-8) do
      runse_diagnostics(netL, measB)
    end
    @test !isempty(diagB.eliminations)
    @test diagB.eliminations[1].id == "P_link"
    # after the elimination the flow stays on the balance value, P and Q
    # (the Q flow carries the capacitor at B1a)
    measF = copy(measB)
    measF[end] = Sparlectra._set_measurement_active(measF[end], false)
    append!(netL.measurements, measF)
    with_state_estimation_config(max_iter = 30, tol = 1e-8, update_net = true) do
      runse!(netL, measF)
    end
    rowsF = calcLinkFlowsSE!(netL)
    @test abs(rowsF[1].pFlow_MW - fP) < 1e-6
    @test abs(rowsF[1].qFlow_MVar - fQ) < 1e-6

    # --- 3b) coupler ring: R1 -> R2 -> R3 -> R1, fed and loaded through
    # lines. The balances fix the flows up to one circulation; P rows on all
    # three links with 5 MW added around the ring give one circulation
    # state, and the flows after the estimate carry it. A 30 MW error on
    # link 1 (sigma 1 MW) is the largest normalized residual, the loop
    # eliminates it and the other two rows keep the circulation.
    # the ring net; `island2` adds a second AC island T-U with its own
    # reference (case 3c), `pendant` a fourth link R3 -> R4 outside the ring
    # with a load at L4 behind R4 (case 3d)
    function _ring_net(; pL2 = 12.0, island2 = false, pendant = false)
      n = Net(name = "se_link_ring", baseMVA = 100.0)
      for b in (island2 ? ("S", "R1", "R2", "R3", "L2", "L3", "T", "U") : ("S", "R1", "R2", "R3", "L2", "L3"))
        addBus!(net = n, busName = b, vn_kV = 110.0)
      end
      if pendant
        addBus!(net = n, busName = "R4", vn_kV = 110.0)
        addBus!(net = n, busName = "L4", vn_kV = 110.0)
        addPIModelACLine!(net = n, fromBus = "R4", toBus = "L4", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
        addProsumer!(net = n, busName = "L4", type = "ENERGYCONSUMER", p = 9.0, q = 3.0)
      end
      addPIModelACLine!(net = n, fromBus = "S", toBus = "R1", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = n, fromBus = "R2", toBus = "L2", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = n, fromBus = "R3", toBus = "L3", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
      addProsumer!(net = n, busName = "S", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "S")
      addProsumer!(net = n, busName = "L2", type = "ENERGYCONSUMER", p = pL2, q = 4.0)
      addProsumer!(net = n, busName = "L3", type = "ENERGYCONSUMER", p = 18.0, q = 5.0)
      addLink!(net = n, fromBus = "R1", toBus = "R2", status = 1)
      addLink!(net = n, fromBus = "R2", toBus = "R3", status = 1)
      addLink!(net = n, fromBus = "R3", toBus = "R1", status = 1)
      pendant && addLink!(net = n, fromBus = "R3", toBus = "R4", status = 1)
      if island2
        addPIModelACLine!(net = n, fromBus = "T", toBus = "U", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
        addProsumer!(net = n, busName = "T", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "T")
        addProsumer!(net = n, busName = "U", type = "ENERGYCONSUMER", p = 15.0, q = 5.0)
      end
      island2 ? runpf!(n, 40, 1e-10, 0; method = :rectangular, islands_enabled = true) : runpf!(n, 40, 1e-10, 0; method = :rectangular)
      return n
    end
    netR = _ring_net()
    calcNetLosses!(netR)
    calcLinkFlowsKCL!(netR)
    fR = [l.pFlow_MW for l in netR.linkVec] .+ 5.0
    measR = generateMeasurementsFromPF(netR; noise = false, stddev = measurementStdDevs(vm = 2e-3, pinj = 1.0, qinj = 1.0, pflow = 1.0, qflow = 1.0))
    for j in 1:3
      addPflowMeasurement!(measR; net = netR, value = fR[j], sigma = 1.0, linkNr = j, id = "P_link_$(j)")
    end
    append!(netR.measurements, measR)
    resR = with_state_estimation_config(max_iter = 30, tol = 1e-8, update_net = true) do
      runse!(netR, measR)
    end
    @test resR.converged
    @test resR.linkCirculations !== nothing && length(resR.linkCirculations) == 1
    rowsR = calcLinkFlowsSE!(netR)
    @test all(abs(rowsR[j].pFlow_MW - fR[j]) < 1e-6 for j in 1:3)
    @test all(r.source == :estimated for r in rowsR)
    i1 = findfirst(m -> m.id == "P_link_1", measR)
    measRB = copy(measR)
    measRB[i1] = Measurement(typ = measR[i1].typ, value = measR[i1].value + 30.0, sigma = 1.0, linkIdx = 1, direction = :from, id = "P_link_1")
    repR = with_state_estimation_config(max_iter = 30, tol = 1e-8) do
      validate_measurements(netR, measRB)
    end
    @test first(repR.measurement_ranking).id == "P_link_1"
    diagR = with_state_estimation_config(max_iter = 30, tol = 1e-8) do
      runse_diagnostics(netR, measRB)
    end
    @test diagR.eliminations[1].id == "P_link_1"
    measRF = copy(measRB)
    measRF[i1] = Sparlectra._set_measurement_active(measRF[i1], false)
    empty!(netR.measurements)
    append!(netR.measurements, measRF)
    with_state_estimation_config(max_iter = 30, tol = 1e-8, update_net = true) do
      runse!(netR, measRF)
    end
    rowsRF = calcLinkFlowsSE!(netR)
    @test all(abs(rowsRF[j].pFlow_MW - fR[j]) < 1e-6 for j in 1:3)
    # a power flow after the estimate replaces its state: the circulation
    # of the estimate no longer applies, and without link measurements the
    # flows are exactly the KCL split of the new state
    empty!(netR.measurements)
    runpf!(netR, 40, 1e-10, 0; method = :rectangular, opt_flatstart = true)
    calcNetLosses!(netR)
    calcLinkFlowsKCL!(netR)
    fKR = [(l.pFlow_MW, l.qFlow_MVar) for l in netR.linkVec]
    rowsPF = calcLinkFlowsSE!(netR)
    @test [(r.pFlow_MW, r.qFlow_MVar) for r in rowsPF] == fKR
    @test all(r.source == :kcl for r in rowsPF)

    # --- 3c) the same ring in a two-island net, measured at another
    # operating point (L2 load 32 instead of 12 MW) than the last power flow
    # of the net. The island path writes the estimate back and recomputes
    # the branch flows and shunt powers, so the link flows after the
    # estimate are the KCL split at the measured point plus the estimated
    # circulation, not the split of the old power flow.
    netI = _ring_net(island2 = true)
    netIB = _ring_net(pL2 = 32.0, island2 = true)
    calcNetLosses!(netIB)
    calcLinkFlowsKCL!(netIB)
    fI = [l.pFlow_MW for l in netIB.linkVec] .+ 5.0
    measI = generateMeasurementsFromPF(netIB; noise = false, stddev = measurementStdDevs(vm = 2e-3, pinj = 1.0, qinj = 1.0, pflow = 1.0, qflow = 1.0))
    for j in 1:3
      addPflowMeasurement!(measI; net = netIB, value = fI[j], sigma = 1.0, linkNr = j, id = "P_link_$(j)")
    end
    append!(netI.measurements, measI)
    resI = with_state_estimation_config(max_iter = 30, tol = 1e-8, update_net = true) do
      runse!(netI, measI)
    end
    @test resI.converged
    @test length(resI.islands) == 2
    rowsI = calcLinkFlowsSE!(netI)
    @test all(abs(rowsI[j].pFlow_MW - fI[j]) < 1e-6 for j in 1:3)

    # --- 3d) a P meter only on the pendant link R3 -> R4 of the ring net:
    # the link lies in no loop, so the row gets no circulation state (its
    # row of nullspace(A) is rounding, 4.8e-18, and a rank cut relative to
    # the largest singular value counted it as a loop, with a zero H column
    # and one dof too few). The row adds exactly one dof, and the flow after
    # the estimate is the balance flow to L4.
    netD = _ring_net(pendant = true)
    calcNetLosses!(netD)
    calcLinkFlowsKCL!(netD)
    fD4 = netD.linkVec[4].pFlow_MW
    @test abs(fD4 - 9.0) < 1.0
    measD0 = generateMeasurementsFromPF(netD; noise = false, stddev = measurementStdDevs(vm = 2e-3, pinj = 1.0, qinj = 1.0, pflow = 1.0, qflow = 1.0))
    measD = copy(measD0)
    addPflowMeasurement!(measD; net = netD, value = fD4, sigma = 1.0, linkNr = 4, id = "P_link_4")
    resD0 = with_state_estimation_config(max_iter = 30, tol = 1e-8, update_net = false) do
      runse!(netD, measD0)
    end
    append!(netD.measurements, measD)
    resD = with_state_estimation_config(max_iter = 30, tol = 1e-8, update_net = true) do
      runse!(netD, measD)
    end
    @test resD0.converged && resD.converged
    @test resD.linkCirculations === nothing
    @test resD.dof == resD0.dof + 1
    rowsD = calcLinkFlowsSE!(netD)
    @test abs(rowsD[4].pFlow_MW - fD4) < 1e-6

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

    # update_shunts on a cluster with a link Q row: the estimated B reaches
    # the shunt draw calcLinkFlowsSE! reads. The model B is halved after the
    # measurements were taken, the estimate finds the true B again, and the
    # link Q flow is the one of the measured state (as in 3a).
    netU = create_se_link_net(shunt_buses = ["B1a"])
    runpf!(netU, 40, 1e-10, 0; method = :rectangular)
    calcNetLosses!(netU)
    b1U = geNetBusIdx(net = netU, busName = "B1")
    brSU = only(br for br in netU.branchVec if Int(br.toBus) == b1U)
    fQU = -brSU.tBranchFlow.qFlow - 3.0
    measU = generateMeasurementsFromPF(netU; includeShuntQ = true, noise = false, stddev = _lstd())
    addQflowMeasurement!(measU; net = netU, value = fQU, sigma = 2.0, linkNr = 1, id = "Q_link")
    append!(netU.measurements, measU)
    shU = only(netU.shuntVec)
    bTrue = imag(shU.y_pu_shunt)
    shU.y_pu_shunt = complex(real(shU.y_pu_shunt), 0.5 * bTrue)
    shU.B_shunt = 0.5 * bTrue
    setShuntEstimation!(netU; busName = "B1a")
    resU = with_state_estimation_config(max_iter = 30, tol = 1e-10, update_net = true, update_shunts = true) do
      runse!(netU, measU)
    end
    @test resU.converged
    @test abs(only(resU.shuntEstimates).B_est - bTrue) < 1e-6
    rowsU = calcLinkFlowsSE!(netU)
    @test abs(rowsU[1].qFlow_MVar - fQU) < 1e-6
    @test abs(rowsU[1].q_meas_residual) < 1e-6

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

    # --- 6) link-fed passive bus: M has couplers only, no branch, so the
    # branch-only marking of a CGMES import types it Isolated. It is a member of its fused node: its
    # zero-injection rows belong to the cluster, otherwise the estimator
    # drops the ZI rows of every member as a partial sum. Truth: M is
    # passive, the link flows into M add up to 0. :chain is A = M = C with
    # lines S-A and C-L (the estimator split it at M into two islands, J
    # 2.6e4 on exact data); :parallel is A = M by two couplers with lines
    # S-A, A-C, C-L (the link flows into M did not add up to 0 on noisy
    # data). dof: states 2 * buses - 1 of the contracted net, rows Vm at
    # all five buses, Pinj/Qinj at S and L, P/Q flows at both ends of each
    # line, the LINKAGG pair, the ZI pair of C where C is outside the
    # cluster. :chain 5 + 4 + 8 + 2 - 5 = 14, :parallel 5 + 4 + 12 + 2 + 2 - 7 = 18.
    for (topo, dofExpected) in ((:chain, 14), (:parallel, 18))
      netF = Net(name = "se_link_fed", baseMVA = 100.0)
      for b in ("S", "A", "M", "C", "L")
        addBus!(net = netF, busName = b, vn_kV = 110.0)
      end
      addPIModelACLine!(net = netF, fromBus = "S", toBus = "A", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
      topo == :parallel && addPIModelACLine!(net = netF, fromBus = "A", toBus = "C", r_pu = 0.01, x_pu = 0.1, b_pu = 0.0, status = 1)
      addPIModelACLine!(net = netF, fromBus = "C", toBus = "L", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
      addProsumer!(net = netF, busName = "S", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "S")
      addProsumer!(net = netF, busName = "L", type = "ENERGYCONSUMER", p = 30.0, q = 10.0)
      addLink!(net = netF, fromBus = "A", toBus = "M", status = 1)
      addLink!(net = netF, fromBus = "M", toBus = topo == :chain ? "C" : "A", status = 1)
      @test last(runpf!(netF, 40, 1e-10, 0; method = :rectangular)) == 0
      calcNetLosses!(netF)
      # the marking of a CGMES import counts branches only
      markIsolatedBuses!(net = netF)
      bM = geNetBusIdx(net = netF, busName = "M")
      @test netF.isoNodes == [bM]
      measF = generateMeasurementsFromPF(netF; noise = true, rng = MersenneTwister(3), stddev = measurementStdDevs(vm = 1e-3, pinj = 0.5, qinj = 0.5, pflow = 0.5, qflow = 0.5))
      @test count(m -> m.typ == Sparlectra.VmMeas && m.busIdx == bM, measF) == 1
      passive = Set(findPassiveBuses(netF))
      filter!(m -> !(m.typ in (Sparlectra.PinjMeas, Sparlectra.QinjMeas) && m.busIdx in passive), measF)
      @test Sparlectra._add_missing_zero_injections!(measF, netF) == 6   # A, M, C
      @test count(m -> m.zero_injection && m.busIdx == bM, measF) == 2
      prepF = Sparlectra._se_prepare(netF, measF)
      @test !any(n -> n.reason == :partial_cluster_injection, prepF.notes)
      append!(netF.measurements, measF)
      # bound of the sum: round-off under constraint, the ZI sigma
      # (1e-3 MVA) under on, where the sum was 0.084 MVA on :parallel
      # without the rows at M
      for (mode, bound) in ((:constraint, 1e-8), (:on, 1e-3))
        resF = with_state_estimation_config(zero_injection = mode, max_iter = 30, update_net = true) do
          runse!(netF, measF)
        end
        @test resF.converged
        @test resF.islands === nothing
        @test resF.dof == dofExpected
        intoM = sum((Int(netF.linkVec[r.linkIdx].toBus) == bM ? 1.0 : -1.0) * complex(r.pFlow_MW, r.qFlow_MVar) for r in calcLinkFlowsSE!(netF))
        @test abs(intoM) < bound
      end
    end
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
    # the P row on the closed link is a WLS row, not an exclusion
    @test !any(n -> n.id == "Pflow_link_1", view.excluded_measurements)

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
    # the row-level method exempts what elimination protects: a
    # zero-injection row and a derived shunt row with a real sigma stay at
    # stage 0 at t = 200, a plain row with the same sigma reaches stage 2
    ziRow = Measurement(typ = Sparlectra.PinjMeas, value = 0.0, sigma = σ, busIdx = 1, id = "ZI_PINJ_bus_1", zero_injection = true)
    shRow = Measurement(typ = Sparlectra.QinjMeas, value = 0.0, sigma = σ, busIdx = 1, id = "SHDERIV_bus_1")
    plainRow = Measurement(typ = Sparlectra.PinjMeas, value = 0.0, sigma = σ, busIdx = 1, id = "Pinj_bus_1")
    @test Sparlectra._robust_stage(ziRow, 100.0, 3.0, 6.0) == (0, σ)
    @test Sparlectra._robust_stage(shRow, 100.0, 3.0, 6.0) == (0, σ)
    @test Sparlectra._robust_stage(plainRow, 100.0, 3.0, 6.0) === Sparlectra._robust_stage(100.0, σ, 3.0, 6.0)

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

    # --- staged leaves zero-injection rows alone: a ZI row with a real
    # sigma at a bus that injects about 30 MW (t near 60) is the row the
    # knees would widen most; it must not appear among the modified rows
    measZI = copy(meas)
    ziBus = meas[findlast(m -> m.typ == Sparlectra.PinjMeas, meas)].busIdx
    push!(measZI, Measurement(typ = Sparlectra.PinjMeas, value = 0.0, sigma = 0.5, busIdx = ziBus, id = "ZI_PINJ_bus_$(ziBus)", zero_injection = true))
    resZI = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_net = false, robust_mode = :staged, topology_precheck = false) do
      runse!(net, measZI)
    end
    @test resZI.converged
    @test resZI.robustRows !== nothing
    @test !any(r -> r.id == "ZI_PINJ_bus_$(ziBus)", resZI.robustRows)

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
    br.ratio_tap_declared = true
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
  net.branchVec[3].ratio_tap_declared = true
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

# --- one-branch two-port against the explicit chain (0.60.0 model check)
# The chain of a transformer with its changers: the ideal transformer of the
# neutral point, the PI section of the ratio part, the auxiliary bus A
# (isAux), the PI section of the phase part, each section with its own ideal
# ratio and its share of the series impedance. A section is
# `(t, z, rev)`: `t` its ratio in the direction F to T, `z` its series
# impedance, `rev` builds it from its far bus, so its impedance sits on the
# near side of its ideal transformer (which then carries 1 / t). A PI
# section always has an impedance, so an impedance behind the additional
# voltage is built with the phase section first (the ratio part is neutral
# there, and the product t does not depend on the order). The ratios come
# from the task formulas (t = 1 + s du, 1 + s du e^(j alpha),
# e^(j 2 atan(s du / 2)), table interpolation, e^(j s dphi)), never from the
# set code. Bus S is a slack behind a line, bus T carries a load.
_tcm_sec(t, z; rev::Bool = false) = (t = ComplexF64(t), z = ComplexF64(z), rev = rev)
# series impedance shares of the ratio and the phase section; the one
# branch carries their sum
_tcm_zr() = 0.0008 + 0.025im
_tcm_zp() = 0.0012 + 0.035im
_tcm_scale(z::ComplexF64, F) = complex(real(z) * F[1], imag(z) * F[2])

function _tcm_feed_net(name::String; aux::Bool)::Net
  net = Net(name = name, baseMVA = 100.0)
  addBus!(net = net, busName = "S", vn_kV = 110.0)
  addBus!(net = net, busName = "F", vn_kV = 110.0)
  aux && addBus!(net = net, busName = "A", vn_kV = 20.0, isAux = true)
  addBus!(net = net, busName = "T", vn_kV = 20.0)
  addProsumer!(net = net, busName = "S", type = "EXTERNALNETWORKINJECTION", referencePri = "S", vm_pu = 1.02, va_deg = 0.0)
  addProsumer!(net = net, busName = "T", type = "ENERGYCONSUMER", p = 30.0, q = 10.0)
  addPIModelACLine!(net = net, fromBus = "S", toBus = "F", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  return net
end

# put a chain section on its branch: ratio, angle and impedance
function _tcm_set_section!(br::Branch, sec)
  t = sec.rev ? inv(sec.t) : sec.t
  br.ratio = br.tap_ratio = abs(t)
  br.angle = br.phase_shift_deg = rad2deg(angle(t))
  br.r_pu, br.x_pu = real(sec.z), imag(sec.z)
  return br
end

function _tcm_chain_net(secs)::Net
  net = _tcm_feed_net("tcm_chain"; aux = true)
  for (sec, near, far) in zip(secs, ("F", "A"), ("A", "T"))
    addPIModelTrafo!(net = net, fromBus = sec.rev ? far : near, toBus = sec.rev ? near : far, r_pu = 0.001, x_pu = 0.01, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
    _tcm_set_section!(net.branchVec[end], sec)
  end
  return net
end

# 2x2 terminal admittance (F, T) of the chain: createYBUS without the
# feeder, then the Kron reduction of the auxiliary bus A
function _tcm_chain_y(net::Net)::Matrix{ComplexF64}
  Y = createYBUS(net = net, sparse = false)
  f, a, t = (geNetBusIdx(net = net, busName = b) for b in ("F", "A", "T"))
  # without the feeder line S-F (branch 1, F is its to end)
  Y[f, f] -= calcAdmittance(net.branchVec[1], 1.0, net.baseMVA)[4]
  red(i, j) = Y[i, j] - Y[i, a] * Y[a, j] / Y[a, a]
  return [red(f, f) red(f, t); red(t, f) red(t, t)]
end

# the one-branch model F-T with the series impedance ZR + ZP: a stored
# changer set (`set`), a nameplate row (`nameplate`, the keywords of
# applyTapNameplate!), or typed models (`taps`, `phase_taps`, `model`) on the
# winding with the neutral point (ratio_n, shift_n); released for the tap
# estimation (every changer), the changers at their neutral step
function _tcm_one_net(ratio_n::Float64, shift_n::Float64; set = nothing, nameplate = nothing, taps = nothing, phase_taps = nothing, model::Symbol = :ideal)::Net
  net = _tcm_feed_net("tcm_one"; aux = false)
  z = _tcm_zr() + _tcm_zp()
  addPIModelTrafo!(net = net, fromBus = "F", toBus = "T", r_pu = real(z), x_pu = imag(z), b_pu = 0.0, ratio = ratio_n, shift_deg = shift_n, status = 1, taps = taps, phase_taps = phase_taps, tap_changer_model = model)
  set === nothing || (net.branchVec[2].se_tap_changers = set(ratio_n, shift_n))
  nameplate === nothing || applyTapNameplate!(net.branchVec[2]; nameplate...)
  roles = Set(ch.role for ch in Sparlectra._effective_tap_changers(net.branchVec[2]).changers)
  setTapEstimation!(net; trafo = 2, mode = length(roles) == 2 ? :both : (:ratio in roles ? :ratio : :pst))
  return net
end

# 2x2 admittance of the one-branch model at the steps s as the estimator
# stamps it (state map and exact overlay), and the scratch branch it used
function _tcm_one_y(map, net::Net, s::Vector{Float64})
  ts = only(map.trafos)
  x = zeros(map.col0 + map.ncols - 1)
  for (c, col) in enumerate(ts.cols)
    col == 0 || (x[col] = s[c])
  end
  e = Sparlectra._tap_overlay(x, map, net)[ts.branch]
  return ([e.Y11 e.Y12; e.Y21 e.Y22], e.branch)
end

_tcm_rel(Yc, Y1) = maximum(abs.(Yc .- Y1) ./ abs.(Y1))

# Maximum relative deviation of the 2x2 admittance of case c (entry by
# entry) over its steps at the neutral point (ratio_n, shift_n), the step
# of that maximum, and with `pf` one power flow per model at the case's
# power-flow step: the largest deviation of the terminal voltages and of
# the flows at both ends in pu (Inf when a power flow fails), and |U| at T.
# The one-branch net carries the overlay's ratio, angle and impedance at
# that step. `chain` replaces the case's chain (other impedance positions).
# The case and its closures are not specialized on (one compile, not one
# per case).
_tcm_compare(c, ratio_n::Float64, shift_n::Float64; chain = nothing, pf::Bool = true) = _tcm_compare_at(c, chain === nothing ? c.chain : chain, ratio_n, shift_n, pf)
function _tcm_compare_at(@nospecialize(c), @nospecialize(chain), ratio_n::Float64, shift_n::Float64, pf::Bool)
  one = c.one(ratio_n, shift_n)
  map, _, _ = Sparlectra._build_tap_state_map(one, 1, c.imp)
  tn = ratio_n * cis(deg2rad(shift_n))
  cnet = _tcm_chain_net(chain(c.pf, tn))
  worst, at = 0.0, c.pf
  for s in c.steps
    for (br, sec) in zip(cnet.branchVec[2:3], chain(s, tn))
      _tcm_set_section!(br, sec)
    end
    d = _tcm_rel(_tcm_chain_y(cnet), first(_tcm_one_y(map, one, s)))
    # a NaN entry (broken interpolation, 0/0 at a band edge) must fail the
    # check, so it is kept as the worst value and never overwritten
    isnan(worst) && continue
    (isnan(d) || d > worst) && ((worst, at) = (d, s))
  end
  pf || return (dev = worst, step = at, pf = NaN, vm_t = NaN)
  secs = chain(c.pf, tn)
  for (br, sec) in zip(cnet.branchVec[2:3], secs)
    _tcm_set_section!(br, sec)
  end
  sb = last(_tcm_one_y(map, one, c.pf))
  br = one.branchVec[2]
  br.tap_ratio, br.phase_shift_deg, br.r_pu, br.x_pu = sb.tap_ratio, sb.phase_shift_deg, sb.r_pu, sb.x_pu
  for n in (one, cnet)
    _, erg = runpf!(n, 40, 1e-12, 0; method = :rectangular)
    erg == 0 || return (dev = worst, step = at, pf = Inf, vm_t = NaN)
    calcNetLosses!(n)
  end
  v(n, b) = (nd = n.nodeVec[geNetBusIdx(net = n, busName = b)]; nd._vm_pu * cis(deg2rad(nd._va_deg)))
  pq(fl) = complex(fl.pFlow, fl.qFlow) / one.baseMVA
  s1, s2 = cnet.branchVec[2], cnet.branchVec[3]
  sf = pq(secs[1].rev ? s1.tBranchFlow : s1.fBranchFlow)
  st = pq(secs[2].rev ? s2.fBranchFlow : s2.tBranchFlow)
  dpf = max(abs(v(one, "F") - v(cnet, "F")), abs(v(one, "T") - v(cnet, "T")), abs(pq(br.fBranchFlow) - sf), abs(pq(br.tBranchFlow) - st))
  return (dev = worst, step = at, pf = dpf, vm_t = abs(v(one, "T")))
end

# The cases: kind, alpha, the one-branch model, the chain with the
# impedance where the kind's data put it, the steps (the whole band in half
# steps; pairs on a grid), the power-flow step and the impedance rule.
function _tcm_cases()
  S = Sparlectra
  ZR, ZP = _tcm_zr(), _tcm_zp()
  band(lo, hi) = [[s] for s in lo:0.5:hi]
  grid(r, p) = [Float64[a, b] for a in r for b in p]
  stored(chs...) = (rn, sn) -> S.SETapChangerSet(ratio_n = rn, shift_n_deg = sn, changers = S.SETapChanger[chs...], origin = :cgmes)
  rat(side) = S.SETapChanger(name = "r", kind = :ratio, role = :ratio, side = side, du = 0.0125, neutral = 0, low = -10, high = 10, step = 0)
  asym(alpha; xc = nothing) = S.SETapChanger(name = "p", kind = :asymmetrical, role = :phase, side = :from, du = 0.01, alpha_deg = alpha, neutral = 0, low = -10, high = 10, step = 0, x_curve = xc, impedance = xc === nothing ? :none : :x_curve)
  sym(; xc = nothing) = S.SETapChanger(name = "p", kind = :symmetrical, role = :phase, side = :from, du = 0.02, neutral = 0, low = -8, high = 8, step = 0, x_curve = xc, impedance = xc === nothing ? :none : :x_curve)
  lin(; xc = nothing) = S.SETapChanger(name = "p", kind = :linear, role = :phase, side = :from, dphi_deg = 1.2, neutral = 0, low = -10, high = 10, step = 0, x_curve = xc, impedance = xc === nothing ? :none : :x_curve)
  # phase table: ratio and angle bend with the step, r and x move
  ptab = [(k, 1.0 + 0.002 * k^2, 3.0 * k - 0.1 * k^2, 1.0 + 2.0 * k / 100, 1.0 + (3.0 * k + 0.5 * k^2) / 100) for k in -4:5]
  rtab = [(k, 1.0 + 0.0125 * k + 0.0004 * k^2, 0.0, 1.0, 1.0) for k in -5:5]
  tabch(role, pts) = S.SETapChanger(name = "e", kind = :tabular, role = role, side = :from, neutral = 0, low = pts[1][1], high = pts[end][1], step = 0, table = [S.SETapPoint(step = p[1], ratio = p[2], angle_deg = p[3], r = p[4], x = p[5]) for p in pts], impedance = role === :phase ? :table : :none)
  # linear interpolation of a table at s: ratio, angle, r and x factors
  function interp(pts, s)
    k = clamp(searchsortedlast([p[1] for p in pts], s), 1, length(pts) - 1)
    w = (s - pts[k][1]) / (pts[k + 1][1] - pts[k][1])
    return ntuple(i -> pts[k][i + 1] + w * (pts[k + 1][i + 1] - pts[k][i + 1]), 4)
  end
  # reactance curve over the changer's angle a (CGMES sin2 form), as the
  # factor X(a(s)) / X(0)
  xcurve(amax) = (0.05, 0.065, amax)
  xfac(a, amax) = (0.05 + 0.015 * (sind(a / 2) / sind(amax / 2))^2) / 0.05
  fa(s, alpha) = 1.0 + s * 0.01 * cis(deg2rad(alpha))
  fsym(s) = cis(2.0 * atan(s * 0.02 / 2.0))
  flin(s) = cis(deg2rad(s * 1.2))
  typed_asym(alpha) = PhaseTapChangerModel(kind = :asymmetrical, step = 0, lowStep = -10, highStep = 10, neutralStep = 0, voltage_step_increment = 0.01, winding_connection_angle_deg = alpha, convention = :reciprocal_from_side)
  cases = Any[
    (name = "a ratio changer, from side", kind = "a", alpha = 0.0, imp = :auto, steps = band(-10, 10), pf = [10.0], one = (rn, sn) -> _tcm_one_net(rn, sn; set = stored(rat(:from))), chain = (s, tn) -> (_tcm_sec(tn * (1 + s[1] * 0.0125), ZR), _tcm_sec(1.0, ZP))),
    (name = "a ratio changer, to side", kind = "a", alpha = 0.0, imp = :auto, steps = band(-10, 10), pf = [-10.0], one = (rn, sn) -> _tcm_one_net(rn, sn; set = stored(rat(:to))), chain = (s, tn) -> (_tcm_sec(tn / (1 + s[1] * 0.0125), ZR), _tcm_sec(1.0, ZP))),
    (name = "e ratio table", kind = "e", alpha = 0.0, imp = :auto, steps = band(-5, 5), pf = [5.0], one = (rn, sn) -> _tcm_one_net(rn, sn; set = stored(tabch(:ratio, rtab))), chain = (s, tn) -> (_tcm_sec(tn * interp(rtab, s[1])[1], ZR), _tcm_sec(1.0, ZP))),
    (name = "d symmetrical", kind = "d", alpha = 0.0, imp = :auto, steps = band(-8, 8), pf = [8.0], one = (rn, sn) -> _tcm_one_net(rn, sn; set = stored(sym())), chain = (s, tn) -> (_tcm_sec(tn, ZR), _tcm_sec(fsym(s[1]), ZP))),
    (name = "d symmetrical, x curve", kind = "d", alpha = 0.0, imp = :auto, steps = band(-8, 8), pf = [8.0], one = (rn, sn) -> _tcm_one_net(rn, sn; set = stored(sym(xc = xcurve(2 * atand(0.08))))), chain = (s, tn) -> (q = xfac(rad2deg(angle(fsym(s[1]))), 2 * atand(0.08)); (_tcm_sec(tn, _tcm_scale(ZR, (1, q))), _tcm_sec(fsym(s[1]), _tcm_scale(ZP, (1, q)))))),
    (name = "f linear", kind = "f", alpha = 0.0, imp = :auto, steps = band(-10, 10), pf = [-10.0], one = (rn, sn) -> _tcm_one_net(rn, sn; set = stored(lin())), chain = (s, tn) -> (_tcm_sec(tn, ZR), _tcm_sec(flin(s[1]), ZP))),
    (name = "f linear, x curve", kind = "f", alpha = 0.0, imp = :auto, steps = band(-10, 10), pf = [-10.0], one = (rn, sn) -> _tcm_one_net(rn, sn; set = stored(lin(xc = xcurve(12.0)))), chain = (s, tn) -> (q = xfac(s[1] * 1.2, 12.0); (_tcm_sec(tn, _tcm_scale(ZR, (1, q))), _tcm_sec(flin(s[1]), _tcm_scale(ZP, (1, q)))))),
    # phase section first: the impedance sits behind the additional voltage
    (name = "e phase table, r and x", kind = "e", alpha = 0.0, imp = :auto, steps = band(-4, 5), pf = [5.0], one = (rn, sn) -> _tcm_one_net(rn, sn; set = stored(tabch(:phase, ptab))), chain = (s, tn) -> (v = interp(ptab, s[1]); F = (v[3], v[4]); (_tcm_sec(tn * v[1] * cis(deg2rad(v[2])), _tcm_scale(ZP, F)), _tcm_sec(1.0, _tcm_scale(ZR, F))))),
    (name = "pair a + d", kind = "a+d", alpha = 0.0, imp = :auto, steps = grid(-10:2:10, -8:2:8), pf = [10.0, -8.0], one = (rn, sn) -> _tcm_one_net(rn, sn; set = stored(rat(:from), sym())), chain = (s, tn) -> (_tcm_sec(tn * (1 + s[1] * 0.0125), ZR), _tcm_sec(fsym(s[2]), ZP))),
  ]
  for alpha in (0.0, 30.0, 60.0, 90.0)
    k = alpha == 90.0 ? "c" : "b"
    push!(cases, (name = "$(k) asymmetrical $(Int(alpha)), CGMES, impedance behind the additional voltage", kind = k, alpha = alpha, imp = :auto, steps = band(-10, 10), pf = [10.0], one = (rn, sn) -> _tcm_one_net(rn, sn; set = stored(asym(alpha))), chain = (s, tn) -> (_tcm_sec(tn * fa(s[1], alpha), ZP), _tcm_sec(1.0, ZR))))
    # typed reciprocal model under impedance_correction: the factor
    # 1 + s du e^(j alpha) divides the ratio, the impedance sits before the
    # additional voltage (phase section reversed)
    push!(cases, (name = "$(k) asymmetrical $(Int(alpha)), typed, impedance before the additional voltage", kind = k, alpha = alpha, imp = :auto, steps = band(-10, 10), pf = [-10.0], one = (rn, sn) -> _tcm_one_net(rn, sn; phase_taps = typed_asym(alpha), model = :impedance_correction), chain = (s, tn) -> (_tcm_sec(tn, ZR), _tcm_sec(1.0 / fa(s[1], alpha), ZP; rev = true))))
  end
  # nameplate Delta-u changer: the factor divides the ratio (side :to)
  push!(cases, (name = "b asymmetrical 60, nameplate, impedance behind the additional voltage", kind = "b", alpha = 60.0, imp = :auto, steps = band(-10, 10), pf = [10.0], one = (rn, sn) -> _tcm_one_net(rn, sn; nameplate = (phase_du_step = 0.01, phase_min_step = -10.0, phase_max_step = 10.0, phase_current_step = 0.0, psi_deg = 60.0)), chain = (s, tn) -> (_tcm_sec(tn / fa(s[1], 60.0), _tcm_zp()), _tcm_sec(1.0, _tcm_zr()))))
  push!(cases, (name = "b asymmetrical 60, CGMES, x curve", kind = "b", alpha = 60.0, imp = :auto, steps = band(-10, 10), pf = [10.0], one = (rn, sn) -> _tcm_one_net(rn, sn; set = stored(asym(60.0; xc = xcurve(rad2deg(angle(fa(10.0, 60.0))))))), chain = (s, tn) -> (q = xfac(rad2deg(angle(fa(s[1], 60.0))), rad2deg(angle(fa(10.0, 60.0)))); (_tcm_sec(tn * fa(s[1], 60.0), _tcm_scale(ZP, (1, q))), _tcm_sec(1.0, _tcm_scale(ZR, (1, q)))))))
  push!(cases, (name = "pair a + b 60, typed, impedance correction", kind = "a+b", alpha = 60.0, imp = :auto, steps = grid(-10:2:10, -10:2:10), pf = [10.0, -10.0], one = (rn, sn) -> _tcm_one_net(rn, sn; taps = PowerTransformerTaps(Vn_kV = 110.0, step = 0, lowStep = -10, highStep = 10, neutralStep = 0, voltageIncrement_kV = 1.375), phase_taps = typed_asym(60.0), model = :impedance_correction), chain = (s, tn) -> (_tcm_sec(tn * (1 + s[1] * 0.0125), ZR), _tcm_sec(1.0 / fa(s[2], 60.0), ZP; rev = true))))
  return cases
end

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

    # --- constructor defaults (0.60.0): the default ratio grid is no
    # declared changer, so a transformer that only carries a ratio has an
    # empty set, no release and setTapEstimation! refuses it; declared, the
    # same fields give the fraction grid 0.9/1.1 of 0.00625
    let net = _tap_gate_net()
      @test isempty(S._effective_tap_changers(net.branchVec[4]).changers) && S._tap_release_mode(net.branchVec[4]) === nothing
      @test_throws ArgumentError setTapEstimation!(net; trafo = 4, mode = :ratio)
      @test S._effective_tap_changers(net.branchVec[1]) === nothing && S._tap_release_mode(net.branchVec[1]) === nothing
      net.branchVec[4].ratio_tap_declared = true
      set = S._effective_tap_changers(net.branchVec[4])
      @test length(set.changers) == 1 && set.changers[1].side === :to && (set.changers[1].low, set.changers[1].high) == (-14, 17)
      @test S._tap_release_mode(net.branchVec[4]) === :ratio
      net.branchVec[4].has_ratio_tap = false
      @test isempty(S._effective_tap_changers(net.branchVec[4]).changers) && S._tap_release_mode(net.branchVec[4]) === nothing
    end
    # --- plain MATPOWER (shipped sp_case118.m, no mpc.sparlectra.tap_changers):
    # its transformers keep their ratio and none is released; the same
    # file with a declared row on one transformer releases exactly that one
    let path = joinpath(dirname(@__DIR__), "data", "mpower", "sp_case118.m")
      net = Sparlectra.createNetFromMatPowerFile(filename = path)
      trafos = [k for (k, br) in enumerate(net.branchVec) if br.ratio != 0.0]
      @test !isempty(trafos) && all(k -> S._tap_release_mode(net.branchVec[k]) === nothing, trafos)
      @test !S._case_has_tap_changers(net)
      @test_throws ArgumentError setTapEstimation!(net; trafo = trafos[1], mode = :ratio)
      applyTapNameplate!(net.branchVec[trafos[1]]; tap_step = 0.0125, tap_min_step = -8.0, tap_max_step = 8.0)
      @test [k for k in trafos if S._tap_release_mode(net.branchVec[k]) !== nothing] == [trafos[1]]
      @test setTapEstimation!(net; trafo = trafos[1], mode = :ratio).s0 == [0.0]
    end
    # --- case6495rte (not shipped; read from the large-case directory, the
    # Web UI case cache under the user's home, so a gate worktree finds it
    # too): TAP 1 on all 9019 branches, which made every branch a
    # transformer before the MATPOWER line rule; 1603 transformers after it.
    # The file has no mpc.sparlectra.tap_changers block, so the import of
    # the SE service releases none of them.
    let path = large_case_path("case6495rte.m")
      if path === nothing
        println("      tap changer data: case6495rte SKIPPED (case6495rte not in the MATPOWER cache)")
      else
        net = Sparlectra.import_case(path, SparlectraConfig(); run_kind = :state_estimation).net
        @test length(net.branchVec) == 9019
        @test length(net.trafos) == 1603
        @test count(br -> S._tap_release_mode(br) !== nothing, net.branchVec) == 0
        println("      tap changer data: case6495rte RAN (", path, ")")
      end
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
    # the typed asymmetrical model with an x curve under impedance_correction
    # (0.60.0, per transformer rule): the curve gives x under auto and off,
    # source table; the analytic factor of the resolver on r applies under
    # auto only (the estimation equals the power flow), under off r stays.
    # The write-back at every step of the band against resolve_winding_taps.
    abr = _tcd_typed_branch(phase_taps = asym(:direct_regulating_vector), model = :impedance_correction)
    aset = S._effective_tap_changers(abr)
    as0 = S._tap_s_from_branch(aset, abr)
    @test all(imp -> S._tap_reactance_source(S._tap_set_under(aset, imp), imp) === :table, (:auto, :off))
    am = abr.tap_winding.phase_taps
    ok = true
    rmoved = false
    for st in am.lowStep:am.highStep
      wb = Dict(imp => (b = deepcopy(abr); S._tap_apply_position!(b, aset, [Float64(st - am.neutralStep)]; s_from = as0, impedance = imp); (b.r_pu, b.x_pu)) for imp in (:auto, :off))
      res = S.resolve_winding_taps(abr.tap_winding; tap_changer_model = :impedance_correction, phase_step = st)
      ok &= isapprox(wb[:auto][2], res.x_alpha; rtol = 1e-12) && wb[:off][2] == wb[:auto][2]
      ok &= isapprox(wb[:auto][1], abr.r_base_pu * res.rx_factor; rtol = 1e-12) && wb[:off][1] == abr.r_pu
      ok &= S._tap_impedance_factors(S._tap_set_under(aset, :off), [Float64(st - am.neutralStep)])[1] == 1.0
      rmoved |= abs(wb[:auto][1] - abr.r_pu) > 1e-9
    end
    @test ok && rmoved

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

    # --- the one branch with the complex t(s) against the explicit chain
    # (ideal transformer, ratio section, auxiliary bus, phase section) per
    # kind a to f, alpha 0/30/60/90 and pairs: the 2x2 terminal admittance
    # after the Kron reduction of the auxiliary bus, relative 1e-10 over
    # the whole band at the neutral shifts 0, 30, 60 and 90 degrees, and one
    # power flow per case (neutral shift 0), terminal voltages and flows to
    # 1e-9 pu
    for c in _tcm_cases(), sn in (0.0, 30.0, 60.0, 90.0)
      r = _tcm_compare(c, 1.05, sn; pf = sn == 0.0)
      ok = r.dev < 1e-10 && (sn != 0.0 || r.pf < 1e-9)
      @test ok
      ok || println("      two-port against chain: ", c.name, ", shift ", sn, ": 2x2 ", r.dev, " at ", r.step, ", power flow ", r.pf)
    end
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
  # and off (there is no on, auto takes its role). The impedance rule is
  # per transformer (0.60.0): per-step data (x curve, table with r and x)
  # follow the step under auto and off, so the truth of both runs has the
  # step impedance and the row names the source table; the analytic
  # correction of the typed case applies under auto only, the truth of its
  # off run keeps the impedance of the start step (source analytic under
  # auto, nominal under off). Until 0.60.0 off froze the per-step data as
  # well. The analytic case runs off once more against the truth with the
  # step impedance: the step still comes back, J is no longer zero. The
  # nameplate fraction case is a regression guard (the fraction cascade of
  # 0.10.0 recovered it as well).
  @testset "Tap estimation per changer type" begin (function ()
    for c in _tap_changer_cases()
      analytic = c.name == "b typed asymmetrical, impedance correction"
      source(imp) = !c.imp ? :nominal : !analytic ? :table : imp === :auto ? :analytic : :nominal
      for imp in (:auto, :off)
        truth = c.net()
        c.build(truth, c.n, source(imp) !== :nominal)
        meas = _tap_truth_meas(truth)
        model = c.net()
        c.build(model, 0, true)
        rel = setTapEstimation!(model; trafo = 3, mode = c.mode)
        @test rel.s0 == [0.0]
        res = _tap_se_run(model, meas; imp = imp)
        row = only(res.tapEstimates)
        ok = res.converged && abs(row.s_est - c.n) <= 1e-6 && row.step_fixed == row.neutral_step + c.n && row.deviation_steps == c.n && res.tapFixation.j_before < 1e-9 && res.tapFixation.j_after < 1e-9 && row.reason === :none && row.reactance_source === source(imp) && !occursin("step not transferred", row.notice)
        @test ok
        ok || println("      tap changer case failed: ", c.name, " (", imp, "): s ", row.s_est, ", J ", res.tapFixation.j_before, " / ", res.tapFixation.j_after, ", source ", row.reactance_source)
      end
      analytic || continue
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
    # a nameplate changer without a transferred step (NaN):
    # the model starts at the middle position, the neutral position 0, the
    # row says so, the estimate is the one from step 0, and the write-back
    # of the fixed step clears the mark
    ca = only(cc for cc in _tap_changer_cases() if cc.name == "a nameplate fraction grid")
    truth_a = ca.net()
    ca.build(truth_a, ca.n, false)
    model_a = ca.net()
    applyTapNameplate!(model_a.branchVec[3]; tap_step = 0.0125, tap_min_step = -8.0, tap_max_step = 8.0, tap_current_step = NaN)
    @test model_a.branchVec[3].ratio_step_assumed && model_a.branchVec[3].tap_ratio == model_a.branchVec[3].ratio
    setTapEstimation!(model_a; trafo = 3, mode = ca.mode)
    start_a = deepcopy(model_a.branchVec[3])
    res_a = _tap_se_run(model_a, _tap_truth_meas(truth_a); update_taps = true)
    row_a = only(res_a.tapEstimates)
    @test res_a.converged && row_a.step_model == 0.0 && row_a.step_fixed == ca.n && row_a.step_assumed && startswith(row_a.notice, "step not transferred, middle position assumed")
    @test !model_a.branchVec[3].ratio_step_assumed
    # the restore of a saved position (fallback, bad-data rounds, island
    # copy) brings the mark back with the step
    Sparlectra._copy_tap_position!(model_a.branchVec[3], start_a)
    @test model_a.branchVec[3].ratio_step_assumed && model_a.branchVec[3].tap_ratio == model_a.branchVec[3].ratio
    # the fallback result rewrites the notice of a :disabled row and keeps
    # the step-assumed text in front
    fb_a = Sparlectra._se_tap_fallback_result(Sparlectra._se_result_with(res_a; tapEstimates = NamedTuple[merge(row_a, (reason = :disabled,))]), :not_converged, 7)
    @test only(fb_a.tapEstimates).notice == "step not transferred, middle position assumed; " * Sparlectra._tap_notice_text(:fallback_not_converged)
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
    # the same transformer from an SCF file (the typed model under
    # tap_changer_model impedance_correction, no per-step data): the
    # write-back to step n differs between auto and off, the source is
    # analytic under auto and nominal under off
    mktempdir() do d
      scf = importSCF(exportSCF(ct.net(); file = joinpath(d, "typed_asym.scf.json")))
      br0 = scf.branchVec[3]
      @test br0.tap_correction === :impedance_correction && br0.tap_winding.phase_taps.kind === :asymmetrical
      set = Sparlectra._effective_tap_changers(br0)
      @test [Sparlectra._tap_reactance_source(set, imp) for imp in (:auto, :off)] == [:analytic, :nominal]
      rx = Dict(imp => (b = deepcopy(br0); Sparlectra._tap_apply_position!(b, set, [Float64(ct.n)]; s_from = [0.0], impedance = imp); (b.r_pu, b.x_pu)) for imp in (:auto, :off))
      @test rx[:off] == (br0.r_pu, br0.x_pu) && rx[:auto][1] != br0.r_pu && rx[:auto][2] != br0.x_pu
    end
    # a stored CGMES set with an x curve and the analytic flag: the curve
    # gives x, r stays, under auto and off alike (the flag is dropped on a
    # transformer with per-step data)
    xset = Sparlectra.SETapChangerSet(ratio_n = 1.0, origin = :cgmes, changers = [Sparlectra.SETapChanger(name = "x", kind = :asymmetrical, role = :phase, side = :from, du = 0.01, alpha_deg = 60.0, neutral = 0, low = -10, high = 10, step = 0.0, x_curve = (0.05, 0.06, 10.0), impedance = :x_curve, rx_correction = true)])
    xb = _tap_gate_net().branchVec[3]
    rxs = Dict(imp => (b = deepcopy(xb); Sparlectra._tap_apply_position!(b, xset, [5.0]; s_from = [0.0], impedance = imp); (b.r_pu, b.x_pu)) for imp in (:auto, :off))
    @test rxs[:auto] == rxs[:off] && rxs[:off][1] == xb.r_pu && rxs[:off][2] != xb.x_pu
    @test all(imp -> Sparlectra._tap_reactance_source(xset, imp) === :table, (:auto, :off))
  end)() end
  return true
end

# The changers of the derivative matrix, each as `make(side, s)` at step s
# (neutral 0) with its band `steps` (lower end, middle, upper end, and a
# table point for the table). Every changer carries impedance data, which
# the `off` variant of the matrix strips (_tap_deriv_matrix):
# a and b at alpha 0 the rx correction of the typed model, b, c, d, f a
# reactance curve up to the largest angle of the band, e a table with r and
# x. `ref2(side, s)` gives (r, dr/ds, psi, dpsi/ds) of t = r e^(j psi) for
# the standard case (ratio or pure phase changer), written from the
# changer formulas, not from the set code; `nothing` for the others.
function _tap_deriv_cases()
  S = Sparlectra
  amax(ch) = max(abs(last(S._tap_changer_factor(ch, ch.low - ch.neutral))), abs(last(S._tap_changer_factor(ch, ch.high - ch.neutral))))
  curved(ch) = S.SETapChanger(name = ch.name, kind = ch.kind, role = ch.role, side = ch.side, du = ch.du, alpha_deg = ch.alpha_deg, dphi_deg = ch.dphi_deg, neutral = ch.neutral, low = ch.low, high = ch.high, step = ch.step, x_curve = (10.0, 13.0, amax(ch)), impedance = :x_curve)
  sgn(side) = side === :from ? 1.0 : -1.0
  asym(alpha) = (side, s) -> begin
    ch = S.SETapChanger(name = "b$(Int(alpha))", kind = :asymmetrical, role = :phase, side = side, du = 0.01, alpha_deg = alpha, neutral = 0, low = -10, high = 10, step = s, rx_correction = alpha == 0.0)
    alpha == 0.0 ? ch : curved(ch)
  end
  table = [S.SETapPoint(step = k, ratio = 1.0 + 0.002 * k^2, angle_deg = 3.0 * k - 0.1 * k^2, r = 1.0 + 2.0 * k / 100, x = 1.0 + (3.0 * k + 0.5 * k^2) / 100) for k in -4:5]
  band(lo, hi) = [Float64(lo), (lo + hi) / 2.0, Float64(hi)]
  return [
    (name = "a ratio", mode = :ratio, steps = band(-10, 10), kinks = Float64[],
      make = (side, s) -> S.SETapChanger(name = "a", kind = :ratio, role = :ratio, side = side, du = 0.0125, neutral = 0, low = -10, high = 10, step = s, rx_correction = true),
      ref2 = (side, s) -> (m = 1.0 + 0.0125 * s; side === :from ? (m, 0.0125, 0.0, 0.0) : (1.0 / m, -0.0125 / m^2, 0.0, 0.0))),
    (name = "b asymmetrical 0", mode = :pst, steps = band(-10, 10), kinks = Float64[], make = asym(0.0), ref2 = nothing),
    (name = "b asymmetrical 30", mode = :pst, steps = band(-10, 10), kinks = Float64[], make = asym(30.0), ref2 = nothing),
    (name = "b asymmetrical 60", mode = :pst, steps = band(-10, 10), kinks = Float64[], make = asym(60.0), ref2 = nothing),
    (name = "c quadrature 90", mode = :pst, steps = band(-10, 10), kinks = Float64[], make = asym(90.0), ref2 = nothing),
    (name = "d symmetrical", mode = :pst, steps = band(-8, 8), kinks = Float64[],
      make = (side, s) -> curved(S.SETapChanger(name = "d", kind = :symmetrical, role = :phase, side = side, du = 0.02, neutral = 0, low = -8, high = 8, step = s)),
      ref2 = (side, s) -> (1.0, 0.0, sgn(side) * 2.0 * atan(0.01 * s), sgn(side) * 0.02 / (1.0 + (0.01 * s)^2))),
    (name = "e table", mode = :pst, steps = vcat(band(-4, 5), 2.0), kinks = [2.0],
      make = (side, s) -> S.SETapChanger(name = "e", kind = :tabular, role = :phase, side = side, neutral = 0, low = -4, high = 5, step = s, table = table, impedance = :table),
      ref2 = nothing),
    (name = "f linear", mode = :pst, steps = band(-10, 10), kinks = Float64[],
      make = (side, s) -> curved(S.SETapChanger(name = "f", kind = :linear, role = :phase, side = side, dphi_deg = 1.2, neutral = 0, low = -10, high = 10, step = s)),
      ref2 = (side, s) -> (1.0, 0.0, sgn(side) * deg2rad(1.2 * s), sgn(side) * deg2rad(1.2))),
  ]
end

# Reference 2, the textbook derivatives of the two-port
# I_f = (y/|t|^2) U_f - (y/conj(t)) U_t, I_t = -(y/t) U_f + y U_t with
# t = r e^(j psi), a fixed series admittance y = g + jb and no shunt arms;
# with d = theta_f - theta_t - psi and k = Vf Vt / r:
#   P_f = Vf^2 g / r^2 - k (g cos d + b sin d)
#   Q_f = -Vf^2 b / r^2 - k (g sin d - b cos d)
#   P_t = Vt^2 g - k (g cos d - b sin d)
#   Q_t = -Vt^2 b + k (g sin d + b cos d)
# differentiated by r and psi and chained with dr/ds and dpsi/ds; returns
# (dP_f, dQ_f, dP_t, dQ_t)/ds in MW and Mvar
function _tap_deriv_textbook(Uf::ComplexF64, Ut::ComplexF64, y::ComplexF64, r::Float64, dr::Float64, psi::Float64, dpsi::Float64, baseMVA::Float64)
  Vf, Vt = abs(Uf), abs(Ut)
  g, b = real(y), imag(y)
  d = angle(Uf) - angle(Ut) - psi
  c, sn = cos(d), sin(d)
  k = Vf * Vt / r
  dPf = (-2.0 * Vf^2 * g / r^3 + k / r * (g * c + b * sn)) * dr - k * (g * sn - b * c) * dpsi
  dQf = (2.0 * Vf^2 * b / r^3 + k / r * (g * sn - b * c)) * dr + k * (g * c + b * sn) * dpsi
  dPt = k / r * (g * c - b * sn) * dr - k * (g * sn + b * c) * dpsi
  dQt = -k / r * (g * sn + b * c) * dr + k * (b * sn - g * c) * dpsi
  return baseMVA .* [dPf, dQf, dPt, dQt]
end

# The derivative matrix: per changer case, side, step and
# impedance_correction one power flow on the gate net with the changer on
# branch 3 at that step, then per measurement type at both ends of the
# branch (P, Q, |I|) and |U| at both terminals:
# - `hcol`, the tap column of the SE (`_measurement_jacobian_fd`, forward
#   difference of step jac_eps / sigma along the analytic tangent),
# - Reference 1, the central difference with Richardson extrapolation
#   (steps h and h/2) of the exact prediction (`_predict_measurements` with
#   `_tap_overlay`) over s, so FD against FD; h is 0.1 steps for P and Q
#   (rounding: the terms of Q reach about 1660 Mvar, at h 0.01 that gave
#   up to 3e-10 against Reference 2) and 0.01 steps for |I| (near a current
#   minimum |I| bends within about 0.3 steps); at a table point the
#   one-sided difference to the right as well (h 1e-3 for P and Q, 1e-4
#   for |I|), the interval the tangent takes there,
# - the forward difference error 0.5 hs g'' with g'' the second central
#   difference (step 0.01) of the prediction along the tangent,
# - Reference 2 for the standard cases under `off` (fixed impedance).
# Under `off` the changer is taken without its per-step data and without
# the analytic flag: per-step data follow the step under both settings
# (0.60.0), so only the plain changer gives the fixed impedance.
# Relative deviation: |a - b| / max(|b|, 1e-3 M) with M the largest
# reference value of the same unit (MW and Mvar, A, pu); 0 when a == b.
# Returns one row per case and type.
function _tap_deriv_matrix()
  S = Sparlectra
  jeps = S.state_estimation_config().jac_eps
  labels = ["P_f", "Q_f", "|I_f|", "P_t", "Q_t", "|I_t|", "|U_f|", "|U_t|"]
  unit = [1, 1, 2, 1, 1, 2, 3, 3]
  rows = NamedTuple[]
  for c in _tap_deriv_cases(), side in (:from, :to), s in c.steps, imp in (:auto, :off)
    net = _tap_gate_net()
    br = net.branchVec[3]
    ch = c.make(side, s)
    imp === :off && (ch = S.SETapChanger(ch.name, ch.mrid, ch.kind, ch.role, ch.side, ch.additive, ch.du, ch.alpha_deg, ch.dphi_deg, ch.step_abs, ch.neutral, ch.low, ch.high, ch.step, ch.table, ch.outside_du, ch.x_curve, ch.x_curve_form, :none, false, false))
    set = S.SETapChangerSet(ratio_n = 1.0, changers = [ch], origin = :cgmes)
    br.se_tap_changers = set
    br.tap_ratio, br.phase_shift_deg = S._tap_set_ratio_shift(set, [s])
    if imp === :auto
      F, F0 = S._tap_impedance_factors(set, [s]), S._tap_impedance_factors(set, [0.0])
      br.r_pu *= F[1] / F0[1]
      br.x_pu *= F[2] / F0[2]
    end
    _, erg = runpf!(net, 40, 1e-12, 0; method = :rectangular)
    fb, tb = Int(br.fromBus), Int(br.toBus)
    flow(typ, dir) = m -> m.typ == typ && m.branchIdx == 3 && m.direction == dir
    sel = [flow(S.PflowMeas, :from), flow(S.QflowMeas, :from), flow(S.ImagMeas, :from), flow(S.PflowMeas, :to), flow(S.QflowMeas, :to), flow(S.ImagMeas, :to), m -> m.typ == S.VmMeas && m.busIdx == fb, m -> m.typ == S.VmMeas && m.busIdx == tb]
    allm = generateMeasurementsFromPF(net; includeImag = true)
    meas = [only(filter(p, allm)) for p in sel]
    setTapEstimation!(net; trafo = 3, mode = c.mode)
    nbus = length(net.nodeVec)
    slack = S._find_slack_idx(net)
    Ybus = createYBUS(net = net)
    map, _, _ = S._build_tap_state_map(net, (nbus - 1) + nbus + 1, imp)
    S._tap_unstamp!(Ybus, net, map)
    ts = map.trafos[1]
    col = ts.cols[1]
    x = S._initial_state_vector(net, slack; flatstart = false, tapSInit = [s])
    H, _, _ = S._measurement_jacobian_fd(meas, net, x, slack, nbus, Ybus; eps = jeps, tapMap = map)
    hcol = Vector(H[:, col])
    at(sv, ov) = (xv = copy(x); xv[col] = sv; S._predict_measurements(meas, net, S._state_to_voltage(xv, slack, nbus), Ybus; tapOverlay = ov(xv)))
    exact(sv) = at(sv, xv -> S._tap_overlay(xv, map, net))
    tangent(sv) = at(sv, xv -> S._tap_overlay_tangent(xv, x, map, net))
    cd(h) = (exact(s + h) .- exact(s - h)) ./ (2.0 * h)
    rich(h) = (4.0 .* cd(h / 2.0) .- cd(h)) ./ 3.0
    od(h) = (exact(s + h) .- exact(s)) ./ h
    right(h) = 2.0 .* od(h / 2.0) .- od(h)
    isI = unit .== 2
    r1 = ifelse.(isI, rich(0.01), rich(0.1))
    kink = s in c.kinks
    ref = kink ? ifelse.(isI, right(1e-4), right(1e-3)) : r1
    hs = jeps / ts.sigma[1]
    trunc = 0.5 * hs .* (tangent(s + 0.01) .- 2.0 .* tangent(s) .+ tangent(s - 0.01)) ./ 1e-4
    M = [maximum((abs(ref[j]) for j in eachindex(ref) if unit[j] == u); init = 0.0) for u in 1:3]
    rel(a, b, j) = a == b ? 0.0 : abs(a - b) / max(abs(b), 1e-3 * M[unit[j]])
    r2 = fill(NaN, length(labels))
    if c.ref2 !== nothing && imp === :off
      V = S._state_to_voltage(x, slack, nbus)
      r, dr, psi, dpsi = c.ref2(side, s)
      tb4 = _tap_deriv_textbook(V[fb], V[tb], 1.0 / complex(br.r_pu, br.x_pu), r, dr, psi, dpsi, net.baseMVA)
      for (k, j) in enumerate((1, 2, 4, 5))
        r2[j] = rel(tb4[k], r1[j], j)
      end
    end
    for j in eachindex(labels)
      push!(rows, (case = c.name, side = side, s = s, imp = imp, typ = labels[j], kink = kink, pf = erg == 0, s0 = ts.s0[1], hcol = hcol[j], ref1 = r1[j], trunc = trunc[j], raw = rel(hcol[j], r1[j], j), corr = rel(hcol[j], ref[j] + trunc[j], j), ref2 = r2[j]))
    end
  end
  return rows
end

function test_tap_estimation_jacobian()::Bool
  # The tap column of the FD Jacobian (the chain rule: dh/dt by forward
  # difference along the analytic tangent dt/ds) against a central
  # difference of the exact prediction over s, with
  # impedance_correction auto (the column includes the impedance part), at
  # an off-step s, for an asymmetrical changer with a reactance curve (b)
  # and a phase table with r and x (e). Then the derivative matrix
  # (`_tap_deriv_matrix`): kinds a to f, alpha 0/30/60/90, changer on the
  # from and the to side, both band ends, the middle and a table point,
  # impedance_correction auto and off (off on the changer without its
  # impedance data), every type at both ends. Reference 2
  # (textbook) against Reference 1 (Richardson) below 1e-10 checks
  # Reference 1. The raw SE column against Reference 1 is not asserted:
  # cells above 1e-6 are the forward difference error, so the column minus
  # 0.5 hs g'' must meet Reference 1 (to the right at a table point) to
  # 1e-6 in every cell.
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
    rows = _tap_deriv_matrix()
    # every case solved and the estimator starts at the step of the branch
    @test all(r -> r.pf && abs(r.s0 - r.s) <= 1e-9, rows)
    std = [r for r in rows if !isnan(r.ref2)]
    @test length(std) == 3 * 2 * 3 * 4
    for r in std
      r.ref2 < 1e-10 || println("      textbook against Richardson: ", r.case, " ", r.side, " s ", r.s, " ", r.typ, ": ", r.ref2)
    end
    @test all(r -> r.ref2 < 1e-10, std)
    for r in rows
      r.corr < 1e-6 || println("      tap column against Richardson plus the forward difference error: ", r.case, " ", r.side, " s ", r.s, " ", r.imp, " ", r.typ, ": ", r.corr)
    end
    @test all(r -> r.corr < 1e-6, rows)
    # |U| does not move with the step: both columns are exactly zero
    @test all(r -> r.hcol == 0.0 && r.ref1 == 0.0, (r for r in rows if r.typ in ("|U_f|", "|U_t|")))
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
    # bus H carries the external grid and hangs on the transformer alone,
    # which the machine transformer rule counts as one; released explicitly
    res = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_net = false, tap_estimation = Sparlectra.TapEstimationConfig(include_machine_transformers = true)) do
      runse!(net, meas)
    end
    row = only(res.tapEstimates)
    @test res.converged && row.reason === :none
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
    @test res.converged && row.reason === :none
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
    # the regulators include machine transformers, which the default keeps
    # out even when released by hand; this test judges the sequential
    # fixation, so it lets them in
    teAll = S.TapEstimationConfig(impedance_correction = mg.config.state_estimation.tap_estimation.impedance_correction, include_machine_transformers = true)
    seRun(n; te = teAll) = S.with_sparlectra_config(mg.config) do
      with_state_estimation_config(() -> runse!(n, meas); max_iter = 50, tol = 1e-6, update_net = false, tap_estimation = te)
    end
    resDefault = Logging.with_logger(() -> seRun(deepcopy(net); te = S.TapEstimationConfig(impedance_correction = teAll.impedance_correction)), Test.TestLogger())
    @test any(t -> t.reason === :machine_transformer, resDefault.tapEstimates)
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

    # switched off on two islands: rows in both islands say disabled, and
    # the merged result has no tapFixation (no island carried a tap state)
    netOff = _two_island_tap_net()
    setTapEstimation!(netOff; trafo = 3, mode = :ratio)
    setTapEstimation!(netOff; trafo = 7, mode = :ratio)
    resOff = with_state_estimation_config(max_iter = 40, tol = 1e-10, update_net = false, tap_estimation = TapEstimationConfig(enabled = false)) do
      runse!(netOff, meas)
    end
    @test resOff.converged && resOff.tapFixation === nothing
    @test sort([t.branch for t in resOff.tapEstimates]) == [3, 7] && all(t -> t.status === :not_estimated && t.reason === :disabled, resOff.tapEstimates)

    # an island WITHOUT tap states that does not converge is no failure of
    # the tap estimation: only trafo 3 (island a) released, island b keeps
    # its -1 step model error and needs 6 iterations, island a estimates
    # and fixes in 3. At max_iter 3 the run does not converge, but the
    # merged failure is :none and no fallback runs; the estimate of island a
    # stays
    netF = _two_island_tap_net()
    setTapEstimation!(netF; trafo = 3, mode = :ratio)
    resF = with_state_estimation_config(max_iter = 3, tol = 1e-10, update_net = false) do
      runse!(netF, meas)
    end
    @test !resF.converged && [i.converged for i in resF.islands] == [true, false]
    @test resF.tapFixation.failure.reason === :none && !resF.tapFixation.fallback_used
    rowF = only(resF.tapEstimates)
    @test rowF.status === :estimated && rowF.fixed && rowF.step_fixed == 2

    # the whole-net fallback runs on the model of run 1: trafo 7 off by -20
    # steps, both trafos released, update_taps, max_iter 4. Island a
    # converges and writes step 2 of trafo 3 into the parent, island b does
    # not converge, so the net falls back. Run 2 must see trafo 3 on its
    # transferred step and equal the run with the switch off. With the
    # restore of the parent taps commented out, trafo 3 kept ratio 0.98765,
    # its row showed step_model 2 and J was 20567 instead of 20859.
    tnet20 = _two_island_tap_net()
    for (k, n) in ((3, 2.0), (7, -20.0))
      br = tnet20.branchVec[k]
      br.tap_ratio = br.ratio / (1.0 + n * step)
    end
    ite, erg = runpf!(tnet20, 40, 1e-12, 0; islands_enabled = true)
    @test erg == 0
    meas20 = generateMeasurementsFromPF(tnet20; includeImag = true, noise = false)
    function _run20(te)
      n = _two_island_tap_net()
      setTapEstimation!(n; trafo = 3, mode = :ratio)
      setTapEstimation!(n; trafo = 7, mode = :ratio)
      r = with_state_estimation_config(max_iter = 4, tol = 1e-10, update_net = true, update_taps = true, tap_estimation = te) do
        runse!(n, meas20)
      end
      return n, r
    end
    # the two fallback lines are warnings here (run 2 does not converge at
    # max_iter 4 either); captured, and both must be there
    fbLog = Test.TestLogger(min_level = Logging.Warn)
    netFb, resFb = Logging.with_logger(() -> _run20(TapEstimationConfig()), fbLog)
    fbLines = [r.message for r in fbLog.logs]
    @test fbLines == ["tap estimation: run 1 not converged (estimation_not_converged, 4 iterations), falling back", "run 2 without tap estimation: not converged (not_converged, 4 iterations)"]
    netOff20, resOff20 = run_with_expected_warnings(() -> _run20(TapEstimationConfig(enabled = false)), String[])
    @test resFb.tapFixation.fallback_used
    @test netFb.branchVec[3].tap_ratio == netOff20.branchVec[3].tap_ratio == 1.0
    @test all(t -> t.reason === :fallback_not_converged && t.step_model == 0, resFb.tapEstimates)
    @test resFb.voltages == resOff20.voltages && resFb.objectiveJ == resOff20.objectiveJ
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
    gnet.branchVec[5].ratio_tap_declared = true
    # the second gate transformer is released below as well: declared too
    gnet.branchVec[4].ratio_tap_declared = true
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

    # 4b) the configuration keeps transformers out of every release, an
    # explicit one included (state_estimation.tap_estimation.exclude and
    # include_machine_transformers). The machine transformer is not estimated
    # by default: a converged run with update_net calculates its step with
    # the method above (truth 2 steps, status :calculated); without
    # update_net nothing is calculated; include_machine_transformers = true
    # estimates it. An excluded transformer keeps its start step and is not
    # calculated, matched by name or by CGMES mRID; an entry that names no
    # transformer warns and has no effect.
    for k in (3, 4, 5)
      setTapEstimation!(gnet; trafo = k, mode = :ratio)
    end
    gnet.cgmes_ids[Sparlectra.CGMESImporter.cgmesKeyPowerTransformer("H2", "L2", Sparlectra.CGMESImporter.cgmesNextParallelIndex!(Dict{Tuple{String,String},Int}(), "H2", "L2"))] = "_T4"
    @test Sparlectra._transformer_mrids(gnet)[4] == "_T4"
    runx(te; update_net = true) = with_state_estimation_config(max_iter = 40, tol = 1e-6, update_net = update_net, tap_estimation = te) do
      runse!(deepcopy(gnet))
    end
    row(r, k) = only(t for t in r.tapEstimates if t.branch == k)
    rdef = runx(TapEstimationConfig())
    @test rdef.converged
    @test (row(rdef, 5).status, row(rdef, 5).reason) == (:calculated, :machine_transformer)
    @test row(rdef, 5).step_fixed == 2.0 && abs(row(rdef, 5).s_est - 2.0) < 0.25 && !row(rdef, 5).fixed
    @test all(k -> row(rdef, k).status === :estimated, (3, 4))
    # the end-step stamp of the SE service bad-data rounds
    # (`_tap_stamp_end_steps!`) leaves the calculated machine transformer at
    # the step the run solved with: model at step 0, truth at 2, branches 3
    # and 4 estimated and fixed. Before the fix the stamp moved branch 5 to
    # its calculated step 1.
    snet = deepcopy(gnet)
    snet.branchVec[5].tap_ratio = snet.branchVec[5].ratio
    rst = @test_logs (:info, r"tap fixation") match_mode = :any with_state_estimation_config(max_iter = 40, tol = 1e-6, update_net = true, tap_estimation = TapEstimationConfig()) do
      runse!(snet)
    end
    @test rst.converged && rst.tapFixation.fixed
    @test row(rst, 5).status === :calculated && row(rst, 5).step_fixed != row(rst, 5).step_model
    r5 = snet.branchVec[5].tap_ratio
    Sparlectra._tap_stamp_end_steps!(snet, rst, TapEstimationConfig().impedance_correction)
    @test snet.branchVec[5].tap_ratio == r5
    rnou = runx(TapEstimationConfig(); update_net = false)
    @test (row(rnou, 5).status, row(rnou, 5).reason) == (:not_estimated, :machine_transformer)
    # the same with a second measured island (B1-B2): the island runs force
    # update_net on their subnets, the calculation follows the caller
    mnet = deepcopy(gnet)
    addBus!(net = mnet, busName = "B1", vn_kV = 110.0)
    addBus!(net = mnet, busName = "B2", vn_kV = 110.0)
    addPIModelACLine!(net = mnet, fromBus = "B1", toBus = "B2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
    addProsumer!(net = mnet, busName = "B1", type = "EXTERNALNETWORKINJECTION", referencePri = "B1", vm_pu = 1.0, va_deg = 0.0)
    addProsumer!(net = mnet, busName = "B2", type = "ENERGYCONSUMER", p = 12.0, q = 4.0)
    ite, erg = runpf!(mnet, 40, 1e-12, 0; islands_enabled = true)
    @test erg == 0
    empty!(mnet.measurements)
    append!(mnet.measurements, generateMeasurementsFromPF(mnet; includeImag = true, noise = false))
    runm(update_net) = with_state_estimation_config(max_iter = 40, tol = 1e-6, update_net = update_net, tap_estimation = TapEstimationConfig()) do
      runse!(deepcopy(mnet))
    end
    rmno = runm(false)
    @test rmno.converged && length(rmno.islands) == 2
    @test (row(rmno, 5).status, row(rmno, 5).reason) == (:not_estimated, :machine_transformer)
    rmup = runm(true)
    @test (row(rmup, 5).status, row(rmup, 5).reason) == (:calculated, :machine_transformer) && row(rmup, 5).step_fixed == 2.0
    rinc = runx(TapEstimationConfig(include_machine_transformers = true))
    @test row(rinc, 5).status === :estimated && row(rinc, 5).step_fixed == 2.0
    rnam = runx(TapEstimationConfig(exclude = [getCompName(gnet.branchVec[3].comp)]))
    @test (row(rnam, 3).status, row(rnam, 3).reason) == (:not_estimated, :excluded)
    @test row(rnam, 3).step_fixed == row(rnam, 3).step_model && row(rnam, 4).status === :estimated
    rmr = @test_logs (:warn, r"exclude names \"nope\", which is no transformer") match_mode = :any runx(TapEstimationConfig(exclude = ["_T4", "nope"]))
    @test (row(rmr, 4).status, row(rmr, 4).reason) == (:not_estimated, :excluded)
    @test row(rmr, 3).status === :estimated && row(rmr, 5).status === :calculated
    # a user exclusion wins over the machine rule: no calculation
    rboth = runx(TapEstimationConfig(exclude = [getCompName(gnet.branchVec[5].comp)]))
    @test (row(rboth, 5).status, row(rboth, 5).reason) == (:not_estimated, :excluded)
    for k in (3, 4, 5)
      setTapEstimation!(gnet; trafo = k, enabled = false)
    end
    # the diagnostics decide on the caller's net like runse!: the GSU ends
    # at G1T, a closed link (the coupler) joins G1T to the machine bus G1.
    # The machine rule treats the closed-link cluster {G1, G1T} as one node,
    # so the GSU is a machine transformer on the net and on the
    # link-contracted copy the diagnostics solve on (G1T fused into G1 and
    # renamed). Before, the link counted as another connection of G1T on the
    # net only (runse! estimated, validate_measurements
    # machine_transformer). The mRID keyed on
    # (H2, G1T) is lost on the contracted copy, so the exclusion is still
    # decided on the caller's net (with exclude "_T5" runse! excluded,
    # validate_measurements estimated before that fix). Branch 5 keeps its
    # index on the contracted copy (a link is no branch).
    lnet = _tap_gate_net()
    addBus!(net = lnet, busName = "G1", vn_kV = 10.0)
    addBus!(net = lnet, busName = "G1T", vn_kV = 10.0)
    addLink!(net = lnet, fromBus = "G1", toBus = "G1T", status = 1)
    addProsumer!(net = lnet, busName = "G1", type = "SYNCHRONOUSMACHINE", p = 20.0, q = 5.0, vm_pu = 1.02)
    addPIModelTrafo!(net = lnet, fromBus = "H2", toBus = "G1T", r_pu = 0.001, x_pu = 0.09, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
    brl = lnet.branchVec[5]
    brl.has_ratio_tap = true
    brl.ratio_tap_declared = true
    lnet.branchVec[4].ratio_tap_declared = true
    brl.tap_step = step
    brl.tap_ratio = brl.ratio / (1.0 + 2 * step)
    ite, erg = runpf!(lnet, 40, 1e-12, 0; method = :rectangular)
    @test erg == 0
    append!(lnet.measurements, generateMeasurementsFromPF(lnet; includeImag = true, noise = false))
    for k in (3, 4, 5)
      setTapEstimation!(lnet; trafo = k, mode = :ratio)
    end
    lnet.cgmes_ids[Sparlectra.CGMESImporter.cgmesKeyPowerTransformer("H2", "G1T", 1)] = "_T5"
    both(te) = with_state_estimation_config(max_iter = 40, tol = 1e-6, update_net = false, tap_estimation = te) do
      (runse!(deepcopy(lnet)), validate_measurements(deepcopy(lnet)).result)
    end
    @test Sparlectra._is_machine_transformer(lnet, 5) && Sparlectra._machine_side(lnet, 5) == 6
    rl, dl = both(TapEstimationConfig())
    @test (row(rl, 5).status, row(rl, 5).reason) == (:not_estimated, :machine_transformer)
    @test (row(dl, 5).status, row(dl, 5).reason) == (:not_estimated, :machine_transformer)
    @test all(k -> row(rl, k).status === :estimated && row(dl, k).status === :estimated, (3, 4))
    # with update_net the step is calculated through the coupler: setpoint
    # and dispatch of the generator at G1, the measured injections of G1
    # and G1T summed (truth 2 steps)
    rlu = with_state_estimation_config(max_iter = 40, tol = 1e-6, update_net = true, tap_estimation = TapEstimationConfig()) do
      runse!(deepcopy(lnet))
    end
    @test (row(rlu, 5).status, row(rlu, 5).reason) == (:calculated, :machine_transformer) && row(rlu, 5).step_fixed == 2.0
    # only a fully measured cluster gives the machine P and Q (the rule of
    # the SE link contraction): without the Q row at G1 the zero-injection
    # row at G1T alone read as a machine at Q = 0 and gave step 1 (the sum
    # before that rule), the generator schedule then gave step 0 against
    # the truth 2 (without the Q row at G1). Without a measured Q nothing
    # is calculated: the row
    # stays not_estimated with q_not_measured_in_group, the transferred
    # step stays, and one info line names the bus without Q
    runl(n) = with_state_estimation_config(max_iter = 40, tol = 1e-6, update_net = true, tap_estimation = TapEstimationConfig()) do
      runse!(n)
    end
    qnet = deepcopy(lnet)
    filter!(m -> !(m.busIdx == 5 && m.typ == Sparlectra.QinjMeas), qnet.measurements)
    rlq = @test_logs (:warn, r"QinjMeas measurement\(s\) on the link cluster at G1 excluded") (:info, r"machine transformer .* not calculated, Q is not measured at G1 of its machine node") match_mode = :any runl(qnet)
    @test (row(rlq, 5).status, row(rlq, 5).reason) == (:not_estimated, :q_not_measured_in_group)
    @test row(rlq, 5).step_fixed == row(rlq, 5).step_model && abs(row(rlq, 5).step_model - 2.0) < 1e-9
    # a machine transformer marked as without a transferred step keeps the
    # step-assumed text in front of both machine notices
    qnet_a = deepcopy(qnet)
    qnet_a.branchVec[5].ratio_step_assumed = true
    rlq_a = @test_logs (:warn, r"QinjMeas measurement\(s\) on the link cluster at G1 excluded") (:info, r"machine transformer .* not calculated") match_mode = :any runl(qnet_a)
    @test row(rlq_a, 5).reason === :q_not_measured_in_group && startswith(row(rlq_a, 5).notice, "step not transferred, middle position assumed; " * Sparlectra._tap_notice_text(:q_not_measured_in_group))
    cnet_a = deepcopy(lnet)
    cnet_a.branchVec[5].ratio_step_assumed = true
    rlc_a = runl(cnet_a)
    @test row(rlc_a, 5).status === :calculated && startswith(row(rlc_a, 5).notice, "step not transferred, middle position assumed; calculated after the estimation")
    # without the P row at G1 the scheduled P is taken (20 MW, the truth
    # here): calculated at step 2 with the note p_from_schedule
    pnet = deepcopy(lnet)
    filter!(m -> !(m.busIdx == 5 && m.typ == Sparlectra.PinjMeas), pnet.measurements)
    rlp = @test_logs (:warn, r"PinjMeas measurement\(s\) on the link cluster at G1 excluded") match_mode = :any runl(pnet)
    @test (row(rlp, 5).status, row(rlp, 5).note, row(rlp, 5).step_fixed) == (:calculated, :p_from_schedule, 2.0)
    # include_machine_transformers: estimated where the measurements see it
    rli, dli = both(TapEstimationConfig(include_machine_transformers = true))
    @test (row(rli, 5).status, row(rli, 5).reason) == (:estimated, :none) && row(rli, 5).step_fixed == 2.0
    @test (row(dli, 5).status, row(dli, 5).reason) == (:estimated, :none) && row(dli, 5).step_fixed == 2.0
    # released, but no voltage row on the machine node {G1, G1T}: the GSU is
    # a bridge, but since the rank test decides the flows at both ends and
    # the Q row at G1 still tell the step apart
    # through the losses: estimated at the truth step on both paths
    vnet5 = deepcopy(lnet)
    filter!(m -> !(m.busIdx in (5, 6) && m.typ == Sparlectra.VmMeas), vnet5.measurements)
    rlv, dlv = with_state_estimation_config(max_iter = 40, tol = 1e-6, update_net = false, tap_estimation = TapEstimationConfig(include_machine_transformers = true)) do
      (runse!(deepcopy(vnet5)), validate_measurements(deepcopy(vnet5)).result)
    end
    @test (row(rlv, 5).status, row(rlv, 5).reason, row(rlv, 5).step_fixed) == (:estimated, :none, 2.0)
    @test (row(dlv, 5).status, row(dlv, 5).reason, row(dlv, 5).step_fixed) == (:estimated, :none, 2.0)
    rlx, dlx = both(TapEstimationConfig(exclude = ["_T5"]))
    @test (row(rlx, 5).status, row(rlx, 5).reason) == (:not_estimated, :excluded)
    @test (row(dlx, 5).status, row(dlx, 5).reason) == (:not_estimated, :excluded)

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
the first phase changer), its ratio, angle and the reactance of the x
curve follow the changer data, and the default run with every tap released
recovers it: the phase regulators fix two steps off their model step and J
after the fixation is zero (noise-free set). The x curve is per-step data,
which follow the step under `off` as well (0.60.0, one rule per
transformer): the `off` set equals the `auto` one, the `off` run recovers
it, and both name the source `table` for B-C; on the imported delivery the
reactance of B-C is the same under both at every step of its band. Until
0.60.0 `off` kept the start impedance and the two sets differed.
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
          SparlectraApp._run_state_estimation_service(case, cfg, joinpath(d, "run_$(imp)"), "a5b_$(imp)", set; max_iter = 50)
        end
        dd = SparlectraApp.to_dict(res)
        @test dd["status"] == "succeeded"
        md = dd["metadata"]
        @test md["se_tap_impedance_correction"] == String(imp)
        @test md["se_tap_j_after"] < 1e-9
        regs = Dict((t["branch"], t["regulator"]) => t for t in md["se_tap_estimates"] if get(t, "status", "") != "calculated")
        # the phase regulators of B-C and C-D: role phase, fixed two steps
        # off the model step; the ratio regulators of D-E and G-H as well,
        # the phase changer of G-H stays at its model step
        for (b, role, dev) in ((2, "phase", 2), (3, "phase", 2), (4, "ratio", 2), (7, "phase", 0))
          r = regs[(b, 1)]
          @test (r["role"], r["deviation_steps"], r["fixed"], r["reason"]) == (role, dev, true, "none")
        end
        @test (regs[(7, 2)]["role"], regs[(7, 2)]["deviation_steps"]) == ("ratio", 2)
        # B-C carries the x curve, the phase table of C-D has no r and x
        @test (regs[(2, 1)]["reactance_source"], regs[(3, 1)]["reactance_source"]) == ("table", "nominal")
      end
      # the x curve of B-C moves the reactance of both truths alike
      @test data_rows(runs[:auto]) == data_rows(runs[:off])
      bc = redirect_stdout(() -> Sparlectra.import_case(z, Sparlectra.resolve_config(cfg, z).config; run_kind = :state_estimation), devnull).net.branchVec[2]
      bset = Sparlectra._effective_tap_changers(bc)
      bch = only(bset.changers)
      @test bch.impedance === :x_curve
      bs0 = Sparlectra._tap_s_from_branch(bset, bc)
      xat(n, imp) = (b = deepcopy(bc); Sparlectra._tap_apply_position!(b, bset, [Float64(n - bch.neutral)]; s_from = bs0, impedance = imp); (b.r_pu, b.x_pu))
      @test all(n -> xat(n, :auto) == xat(n, :off), bch.low:bch.high)
      @test xat(bch.low, :off)[2] != bc.x_pu

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

# --- Verification ladder of the tap estimation, stages 1 to 4. One case
# is a net with released transformers, the moves of the truth (steps per
# changer, NaN = the changer is not released and stays) and the measurement
# set: Vm, injections, P and Q at both ends, |I| at both ends, and Va at
# every bus where `va` is set. The truth is the model with the released changers moved
# by the given steps (clamped to the band), solved by a power flow; the
# model starts at its own steps.

# the 3-bus net S-F-T of the ladder: feeder line S-F, transformer F-T
# (branch 2) with the stored changer set `chs` at its start steps
function _lad_net3(chs)::Net
  S = Sparlectra
  net = Net(name = "lad3", baseMVA = 100.0)
  addBus!(net = net, busName = "S", vn_kV = 110.0)
  addBus!(net = net, busName = "F", vn_kV = 110.0)
  addBus!(net = net, busName = "T", vn_kV = 20.0)
  addProsumer!(net = net, busName = "S", type = "EXTERNALNETWORKINJECTION", referencePri = "S", vm_pu = 1.02, va_deg = 0.0)
  addProsumer!(net = net, busName = "T", type = "ENERGYCONSUMER", p = 30.0, q = 10.0)
  addPIModelACLine!(net = net, fromBus = "S", toBus = "F", r_pu = 0.01, x_pu = 0.08, b_pu = 0.0, status = 1)
  addPIModelTrafo!(net = net, fromBus = "F", toBus = "T", r_pu = 0.002, x_pu = 0.06, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
  set = S.SETapChangerSet(ratio_n = 1.0, changers = S.SETapChanger[chs...], origin = :cgmes)
  br = net.branchVec[2]
  br.se_tap_changers = set
  br.tap_ratio, br.phase_shift_deg = S._tap_set_ratio_shift(set, Float64[c.step - c.neutral for c in set.changers])
  return net
end

# the changers of the 3-bus net, one per kind (a to f; b at alpha 30 and
# 60, c is the asymmetrical changer at alpha 90), neutral 0, start step 0
function _lad_changer(k::Symbol)
  S = Sparlectra
  k === :a && return S.SETapChanger(name = "a", kind = :ratio, role = :ratio, side = :from, du = 0.0125, neutral = 0, low = -10, high = 10, step = 0)
  k === :b30 && return S.SETapChanger(name = "b30", kind = :asymmetrical, role = :phase, side = :from, du = 0.01, alpha_deg = 30.0, neutral = 0, low = -10, high = 10, step = 0)
  k === :b60 && return S.SETapChanger(name = "b60", kind = :asymmetrical, role = :phase, side = :from, du = 0.01, alpha_deg = 60.0, neutral = 0, low = -10, high = 10, step = 0)
  k === :c && return S.SETapChanger(name = "c", kind = :asymmetrical, role = :phase, side = :from, du = 0.01, alpha_deg = 90.0, neutral = 0, low = -10, high = 10, step = 0)
  k === :d && return S.SETapChanger(name = "d", kind = :symmetrical, role = :phase, side = :from, du = 0.02, neutral = 0, low = -8, high = 8, step = 0)
  k === :e && return S.SETapChanger(name = "e", kind = :tabular, role = :phase, side = :from, neutral = 0, low = -4, high = 5, step = 0, table = [S.SETapPoint(step = n, ratio = 1.0 + 0.002 * n^2, angle_deg = 3.0 * n - 0.1 * n^2) for n in -4:5])
  k === :f && return S.SETapChanger(name = "f", kind = :linear, role = :phase, side = :from, dphi_deg = 1.2, neutral = 0, low = -10, high = 10, step = 0)
  throw(ArgumentError("unknown ladder changer $(k)"))
end

# a shipped case through the import path of the SE service (case
# configuration), without its measurements; `nothing` when the file is
# missing (only the CGMES cache may be missing, _lad_nets fails a missing
# shipped file)
function _lad_import(path::AbstractString)
  S = Sparlectra
  isfile(path) || return nothing
  cfg = S.resolve_config(S.DEFAULT_SPARLECTRA_CONFIG_PATH, path).config
  m = redirect_stdout(() -> S.import_case(path, cfg; run_kind = :state_estimation), devnull)
  empty!(m.net.measurements)
  return (net = m.net, config = m.config)
end

# The ladder cases: name, net key, moves per branch, Va flag. The moves
# are +2 and -3 steps on the imported cases; on the 3-bus net 3 to 6 steps,
# the steps of the changer-type test. Va makes every released changer observable on the
# radial 3-bus net and on the radial 20 kV side of sp_case14; sp_casePST
# (warmup_casePST is the same net from MATPOWER) and MicroGrid BE need none.
# sp_case14 carries its step data in the case file (tap_changer block of the
# SCF, the nameplate fields of applyTapNameplate!): br16 and br17 each a
# ratio changer and a 1.25 deg phase changer.
function _lad_cases()
  return [
    (name = "3-bus a", net = :n3a, moves = Dict(2 => [4.0]), va = true),
    (name = "3-bus b30", net = :n3b30, moves = Dict(2 => [4.0]), va = true),
    (name = "3-bus b60", net = :n3b60, moves = Dict(2 => [4.0]), va = true),
    (name = "3-bus c", net = :n3c, moves = Dict(2 => [-3.0]), va = true),
    (name = "3-bus d", net = :n3d, moves = Dict(2 => [5.0]), va = true),
    (name = "3-bus e", net = :n3e, moves = Dict(2 => [3.0]), va = true),
    (name = "3-bus f", net = :n3f, moves = Dict(2 => [6.0]), va = true),
    (name = "3-bus pair a + d", net = :n3ad, moves = Dict(2 => [4.0, 5.0]), va = true),
    (name = "3-bus pair a + b60", net = :n3ab60, moves = Dict(2 => [4.0, -4.0]), va = true),
    (name = "sp_case14 ratio", net = :case14, moves = Dict(16 => [2.0, NaN], 17 => [-3.0, NaN]), va = false),
    (name = "sp_case14 phase", net = :case14, moves = Dict(16 => [NaN, 2.0], 17 => [NaN, -3.0]), va = true),
    (name = "sp_case14 pairs", net = :case14, moves = Dict(16 => [2.0, 2.0], 17 => [-3.0, -3.0]), va = true),
    (name = "sp_casePST ratio br7", net = :pst, moves = Dict(7 => [2.0]), va = false),
    (name = "sp_casePST phase br8", net = :pst, moves = Dict(8 => [-3.0]), va = false),
    (name = "sp_casePST pair", net = :pst, moves = Dict(7 => [2.0], 8 => [-3.0]), va = false),
    (name = "warmup_casePST ratio br7", net = :wpst, moves = Dict(7 => [2.0]), va = false),
    (name = "warmup_casePST phase br8", net = :wpst, moves = Dict(8 => [-3.0]), va = false),
    (name = "warmup_casePST pair", net = :wpst, moves = Dict(7 => [2.0], 8 => [-3.0]), va = false),
    (name = "MicroGrid BE phase br8", net = :mg, moves = Dict(8 => [2.0]), va = false),
    (name = "MicroGrid BE ratio br9", net = :mg, moves = Dict(9 => [2.0]), va = false),
    (name = "MicroGrid BE all four", net = :mg, moves = Dict(8 => [2.0], 9 => [2.0], 10 => [2.0], 12 => [2.0]), va = false),
  ]
end

# The nets of the ladder by key, imported once: (net, config) with config
# `nothing` for the 3-bus nets. The shipped case files must exist (a
# missing one fails here, its cases drop out after that failure).
# MicroGrid BE needs the CGMES cache: RAN or SKIPPED is printed once.
function _lad_nets()
  root = dirname(@__DIR__)
  nets = Dict{Symbol,Any}()
  for (key, kinds) in ((:n3a, [:a]), (:n3b30, [:b30]), (:n3b60, [:b60]), (:n3c, [:c]), (:n3d, [:d]), (:n3e, [:e]), (:n3f, [:f]), (:n3ad, [:a, :d]), (:n3ab60, [:a, :b60]))
    nets[key] = (net = _lad_net3([_lad_changer(k) for k in kinds]), config = nothing)
  end
  for (key, path) in ((:case14, joinpath(root, "data", "scf", "sp_case14.scf.json")), (:pst, joinpath(root, "data", "scf", "sp_casePST.scf.json")), (:wpst, joinpath(root, "data", "mpower", "warmup_casePST.m")))
    @test isfile(path)
    isfile(path) || println("      tap estimation ladder: shipped case ", path, " MISSING, its cases are not run")
    nets[key] = _lad_import(path)
  end
  zip = joinpath(Sparlectra.CGMESImporter.cgmesTestSetCacheDir(), "cases", "cgmes_microgrid_be.zip")
  nets[:mg] = _lad_import(zip)
  println("      tap estimation ladder MicroGrid BE: ", nets[:mg] === nothing ? "SKIPPED ($(zip) not in the CGMES cache)" : "RAN")
  return nets
end

_lad_imp(cfg) = cfg === nothing ? :auto : cfg.state_estimation.tap_estimation.impedance_correction

# move the changers of `moves` from their start by the given steps (clamped
# to the band) and return the true step of each moved changer by
# (branch, changer)
function _lad_move!(net::Net, moves, imp::Symbol)::Dict{Tuple{Int,Int},Float64}
  S = Sparlectra
  truth = Dict{Tuple{Int,Int},Float64}()
  for (k, d) in moves
    br = net.branchVec[k]
    set = S._effective_tap_changers(br)
    s0 = S._tap_s_from_branch(set, br)
    s1 = copy(s0)
    for c in eachindex(s0)
      isnan(d[c]) && continue
      lo, hi = S._tap_band(set.changers[c])
      s1[c] = clamp(round(s0[c]) + d[c], lo, hi)
      truth[(k, c)] = set.changers[c].neutral + s1[c]
    end
    S._tap_apply_position!(br, set, s1; s_from = s0, impedance = imp)
  end
  return truth
end

# power flow of a ladder net: the 3-bus nets with the rectangular solver,
# the imported cases with their case configuration; returns the error code
function _lad_pf!(net::Net, cfg)::Int
  S = Sparlectra
  cfg === nothing && return last(runpf!(net, 40, 1e-12, 0; method = :rectangular))
  return last(redirect_stdout(() -> runpf!(net; config = S._copy_powerflow_with(cfg.powerflow; tol = 1e-10, max_iter = 40)), devnull))
end

# release the changers of `moves` (a NaN entry keeps that changer out)
function _lad_release!(net::Net, moves)
  S = Sparlectra
  for (k, d) in moves
    set = S._effective_tap_changers(net.branchVec[k])
    roles = Set(set.changers[c].role for c in eachindex(d) if !isnan(d[c]))
    setTapEstimation!(net; trafo = k, mode = length(roles) == 2 ? :both : (:ratio in roles ? :ratio : :pst))
  end
  return net
end

# One SE run of the ladder. Machine transformers (sp_casePST br7,
# MicroGrid br10) are let in: the ladder judges the estimator, the default
# keeps them out by rule, not by observability.
function _lad_se(net::Net, meas, cfg; tol::Float64 = 1e-10, fallback::Bool = true)
  S = Sparlectra
  te = S.TapEstimationConfig(impedance_correction = _lad_imp(cfg), include_machine_transformers = true, fallback = fallback)
  run() = Logging.with_logger(Logging.NullLogger()) do
    with_state_estimation_config(() -> runse!(net, meas); max_iter = 50, tol = tol, update_net = false, tap_estimation = te)
  end
  return _lad_with(run, cfg)
end

# run f under the case configuration `cfg` (`nothing`: the active one)
_lad_with(f, cfg) = cfg === nothing ? f() : Sparlectra.with_sparlectra_config(f, cfg)

# truth net and exact measurements of a case
function _lad_truth(base, moves, va::Bool; noise::Bool = false, rng = Random.MersenneTwister(1))
  truth = deepcopy(base.net)
  steps = _lad_move!(truth, moves, _lad_imp(base.config))
  @test _lad_pf!(truth, base.config) == 0
  return (net = truth, steps = steps, meas = generateMeasurementsFromPF(truth; includeImag = true, includeVa = va, noise = noise, rng = rng))
end

# Stage 1 on a ladder net: the derivative check of the tap columns
# (_tap_deriv_matrix) for every tap column of the released net `net` at the
# power-flow state of its start steps and every row of `meas`. The SE column
# (forward difference along the tangent) minus the forward difference error
# 0.5 hs g'' against Reference 1 (central difference with Richardson
# extrapolation of the exact prediction; h 0.1 steps for P and Q, 0.01 for
# |I|; at a table point the one-sided difference to the right, h 1e-3 and
# 1e-4). Relative deviation |a - b| / max(|b|, 1e-3 M), M the largest
# reference value of the same unit. Returns the largest deviation and its
# cell.
function _lad_deriv(net::Net, meas::Vector{Measurement}, imp::Symbol)
  S = Sparlectra
  jeps = S.state_estimation_config().jac_eps
  nbus = length(net.nodeVec)
  slack = S._find_slack_idx(net)
  Ybus = createYBUS(net = net)
  col0 = (nbus - 1) + nbus + 1
  map, _, _ = S._build_tap_state_map(net, col0, imp)
  S._tap_unstamp!(Ybus, net, map)
  s0 = zeros(map.ncols)
  for ts in map.trafos, c in eachindex(ts.cols)
    ts.cols[c] == 0 || (s0[ts.cols[c] - col0 + 1] = ts.s0[c])
  end
  x = S._initial_state_vector(net, slack; flatstart = false, tapSInit = s0)
  H, _, _ = S._measurement_jacobian_fd(meas, net, x, slack, nbus, Ybus; eps = jeps, tapMap = map)
  unit(m) = m.typ == S.ImagMeas ? 2 : m.typ == S.VmMeas ? 3 : m.typ == S.VaMeas ? 4 : 1
  u = unit.(meas)
  worst = (dev = 0.0, cell = "")
  for ts in map.trafos, c in eachindex(ts.cols)
    col = ts.cols[c]
    col == 0 && continue
    s = x[col]
    at(sv, ov) = (xv = copy(x); xv[col] = sv; S._predict_measurements(meas, net, S._state_to_voltage(xv, slack, nbus), Ybus; tapOverlay = ov(xv)))
    exact(sv) = at(sv, xv -> S._tap_overlay(xv, map, net))
    tangent(sv) = at(sv, xv -> S._tap_overlay_tangent(xv, x, map, net))
    cd(h) = (exact(s + h) .- exact(s - h)) ./ (2.0 * h)
    rich(h) = (4.0 .* cd(h / 2.0) .- cd(h)) ./ 3.0
    od(h) = (exact(s + h) .- exact(s)) ./ h
    right(h) = 2.0 .* od(h / 2.0) .- od(h)
    isI = u .== 2
    kink = ts.set.changers[c].kind === :tabular && isinteger(s)
    ref = kink ? ifelse.(isI, right(1e-4), right(1e-3)) : ifelse.(isI, rich(0.01), rich(0.1))
    hs = jeps / ts.sigma[c]
    trunc = 0.5 * hs .* (tangent(s + 0.01) .- 2.0 .* tangent(s) .+ tangent(s - 0.01)) ./ 1e-4
    M = [maximum((abs(ref[j]) for j in eachindex(ref) if u[j] == k); init = 0.0) for k in 1:4]
    for j in eachindex(meas)
      a, b = H[j, col], ref[j] + trunc[j]
      d = a == b ? 0.0 : abs(a - b) / max(abs(b), 1e-3 * M[u[j]], floatmin())
      d > worst.dev && (worst = (dev = d, cell = string("branch ", ts.branch, " changer ", c, " ", meas[j].id)))
    end
  end
  return worst
end

function test_tap_estimation_ladder()::Bool
  # Verification ladder of the tap estimation, stages 1, 2 and 4 (stage 3,
  # the noisy draws, is test_tap_estimation_ladder_noise). Cases:
  # _lad_cases (3-bus net per changer kind and two pairs, sp_case14 with the
  # step data of its case file, sp_casePST, warmup_casePST, MicroGrid BE from
  # the CGMES cache), ratio and phase changers alone and as a pair.
  # Stage 1: the derivative check of the tap columns on every case
  # (_lad_deriv), below 1e-6. Stage 2: exact measurements of the truth,
  # model at its own steps: every released changer is estimated and its
  # fixed step is the true step (s_est within 1e-6 of it), J after the
  # fixation below 1e-8. Stage 4: measurement sets in which one changer
  # drops out by the rank test (no measurement at the transformer and its
  # terminals; no angle anchor): that changer not_estimated with a rank
  # reason and its transferred step, the other changers estimated, the run
  # converged. A mesh without P flow rows does not make a changer drop out
  # (sp_casePST, warmup_casePST, MicroGrid BE): every changer is estimated
  # at the true step. Stage 4 asks for a case with a wrong start
  # step that triggers the fallback: no start step alone gave a run that
  # fails to converge on these nets or on sp_case5; the fallback is tested
  # in "State estimation tap fallback" (sp_case5, triggered by a gross Vm
  # error).
  @testset "Tap estimation ladder" begin (function ()
    S = Sparlectra
    nets = _lad_nets()
    for c in _lad_cases()
      base = nets[c.net]
      base === nothing && continue
      tr = _lad_truth(base, c.moves, c.va)
      # stage 1 at the start steps of the model
      dnet = deepcopy(base.net)
      @test _lad_pf!(dnet, base.config) == 0
      _lad_release!(dnet, c.moves)
      w = _lad_with(base.config) do
        with_state_estimation_config(() -> _lad_deriv(dnet, tr.meas, _lad_imp(base.config)); tap_estimation = S.TapEstimationConfig(impedance_correction = _lad_imp(base.config), include_machine_transformers = true))
      end
      @test w.dev < 1e-6
      w.dev < 1e-6 || println("      ladder stage 1 ", c.name, ": ", w.dev, " at ", w.cell)
      # stage 2
      net = _lad_release!(deepcopy(base.net), c.moves)
      res = _lad_se(net, tr.meas, base.config)
      rows = [t for t in res.tapEstimates if haskey(tr.steps, (t.branch, t.regulator))]
      ok = res.converged && length(rows) == length(tr.steps) && res.tapFixation.j_after < 1e-8
      ok &= all(t -> t.status === :estimated && t.step_fixed == tr.steps[(t.branch, t.regulator)] && abs(t.s_est + t.neutral_step - t.step_fixed) <= 1e-6, rows)
      @test ok
      ok || println("      ladder stage 2 ", c.name, ": converged ", res.converged, ", J after ", res.tapFixation === nothing ? nothing : res.tapFixation.j_after, ", ", [(t.branch, t.regulator, t.status, t.reason, t.s_est, t.step_fixed, get(tr.steps, (t.branch, t.regulator), NaN)) for t in res.tapEstimates])
    end

    # stage 4: (name, net, moves, Va, measurement filter on the truth net,
    # dropped changer (branch, changer), `nothing` when none drops out)
    nomeas(k) = t -> (b = t.branchVec[k]; ends = (Int(b.fromBus), Int(b.toBus)); m -> !((m.branchIdx !== nothing && m.branchIdx == k) || (m.branchIdx === nothing && m.busIdx !== nothing && m.busIdx in ends)))
    # a mesh without P flow: no P flow row on the branches of the mesh (Q
    # flows, |I|, Vm and the injections stay)
    nopf(ks) = t -> (m -> !(m.typ == S.PflowMeas && m.branchIdx in ks))
    keep = t -> (m -> true)
    variants = [
      # no angle anchor: a phase changer on a bridge has none, the radial
      # counterpart of the mesh without P flow. The 3-bus net has one
      # transformer and no mesh: without Va the phase changer of a pair has
      # no angle anchor; a transformer without measurements at its
      # terminals would leave bus T unobservable
      ("3-bus pair a + d, no angle anchor", :n3ad, Dict(2 => [4.0, 5.0]), false, keep, (2, 2)),
      ("3-bus pair a + b60, no angle anchor", :n3ab60, Dict(2 => [4.0, -4.0]), false, keep, (2, 2)),
      ("sp_case14, no measurement at br16", :case14, Dict(16 => [2.0, NaN], 17 => [-3.0, NaN]), false, nomeas(16), (16, 1)),
      # no mesh runs through a transformer of sp_case14 (the 20 kV sides hang
      # radially on br16 and br17), so no mesh variant exists there; without
      # Va the phase changer has no angle anchor (it sits on a bridge)
      ("sp_case14, no angle anchor", :case14, Dict(16 => [2.0, 2.0], 17 => [-3.0, NaN]), false, keep, (16, 2)),
      ("sp_casePST, no measurement at br8", :pst, Dict(7 => [2.0], 8 => [-3.0]), false, nomeas(8), (8, 1)),
      # br8 lies in the mesh 2-3-4-5-7 (br8, br6, br3, br5, br1). Without the
      # P flow rows of the mesh nothing drops out and br8 is estimated at the
      # true step. The P injections at the mesh buses carry it: it stays at
      # the truth without the |I| rows of the mesh as well, and as a
      # symmetrical changer without magnitude share. With the P injections
      # of the mesh gone too, no active power row is left at the mesh (bus 2
      # has none, and the P row of bus 9, linked to bus 2, is dropped as a
      # partial cluster sum); the WLS from the flat start then converges to
      # a second solution with the angle differences on br1 and br3 turned
      # over, and br8 ends at a wrong step: set C, tested below.
      ("sp_casePST, mesh of br8 without P flow", :pst, Dict(7 => [2.0], 8 => [-3.0]), false, nopf([1, 3, 5, 6, 8]), nothing),
      ("warmup_casePST, no measurement at br8", :wpst, Dict(7 => [2.0], 8 => [-3.0]), false, nomeas(8), (8, 1)),
      ("warmup_casePST, mesh of br8 without P flow", :wpst, Dict(7 => [2.0], 8 => [-3.0]), false, nopf([1, 3, 5, 6, 8]), nothing),
      ("MicroGrid BE, no measurement at br9", :mg, Dict(8 => [2.0], 9 => [2.0], 10 => [2.0], 12 => [2.0]), false, nomeas(9), (9, 1)),
      # br8 lies in the mesh 4-9-3-11 (br8, br9, br12, br11), asymmetrical
      # as on sp_casePST
      ("MicroGrid BE, mesh of br8 without P flow", :mg, Dict(8 => [2.0], 9 => [2.0], 10 => [2.0], 12 => [2.0]), false, nopf([8, 9, 11, 12]), nothing),
    ]
    npartner = 0
    for (name, key, moves, va, filt, drop) in variants
      base = nets[key]
      base === nothing && continue
      tr = _lad_truth(base, moves, va)
      meas = filter(filt(tr.net), tr.meas)
      res = _lad_se(_lad_release!(deepcopy(base.net), moves), meas, base.config)
      rows = [t for t in res.tapEstimates if haskey(tr.steps, (t.branch, t.regulator))]
      out = [t for t in rows if (t.branch, t.regulator) == drop]
      rest = [t for t in rows if (t.branch, t.regulator) != drop]
      ok = res.converged && length(out) == (drop === nothing ? 0 : 1) && length(rest) == length(tr.steps) - length(out)
      ok &= all(t -> t.status === :not_estimated && t.reason in S._TAP_RANK_REASONS && t.step_fixed == t.step_model && !isempty(t.notice), out)
      # the others are estimated; a ratio changer that only holds together
      # with the frozen phase position of its transformer carries the note
      # :conditional and is not compared with the truth
      ok &= all(t -> t.status === :estimated && (t.note === :conditional || t.step_fixed == tr.steps[(t.branch, t.regulator)]), rest)
      # an estimated partner on the transformer of the dropped changer
      # (the pairs on the 3-bus net and the ratio changer of sp_case14 br16
      # without Va) names it: the estimate holds for its transferred step
      for o in out, t in rest
        t.branch == o.branch || continue
        npartner += 1
        ok &= occursin(string(" of this transformer is not determinable (", o.reason, "), this estimate depends on its transferred step ", S._tap_notice_step(o.step_model)), t.notice)
      end
      @test ok
      ok || println("      ladder stage 4 ", name, ": converged ", res.converged, ", ", [(t.branch, t.regulator, t.status, t.reason, t.note, t.s_est, t.step_fixed, t.step_model, get(tr.steps, (t.branch, t.regulator), NaN)) for t in res.tapEstimates])
    end
    # the partner notice ran on the three pairs (3-bus a + d, a + b60,
    # sp_case14 br16)
    @test npartner == 3

    # Set C: the mesh of br8 without any active power row (no P flow and no
    # |I| on its branches, no P injection at its buses). On sp_casePST and
    # warmup_casePST bus 2 is linked to bus 9; the P row of bus 9 is dropped
    # as a partial cluster sum, so no P row the estimator uses is left in the
    # mesh 2-3-4-5-7. The precheck names that mesh (kind
    # :mesh_without_p_measurement), on the result of runse! as well. With
    # the P injection at bus 2 kept, the cluster 2-9 is fully measured, the
    # warning stays off and both changers end at the true step. A P row at
    # bus 3 or 4 alone also clears the warning, but the WLS from the flat
    # start still finds a second solution there (J 16.2 and 2.4 at dof 34):
    # the check is a necessary condition, not a proof of uniqueness.
    setc_mesh = Dict(:pst => [1, 3, 5, 6, 8], :wpst => [1, 3, 5, 6, 8], :mg => [8, 9, 11, 12])
    setc_moves = Dict(:pst => Dict(7 => [2.0], 8 => [-3.0]), :wpst => Dict(7 => [2.0], 8 => [-3.0]), :mg => Dict(8 => [2.0], 9 => [2.0], 10 => [2.0], 12 => [2.0]))
    meshbranches(fs) = [[parse(Int, x.captures[1]) for x in eachmatch(r"(\d+) \(", f.location)] for f in something(fs, NamedTuple[]) if f.kind === :mesh_without_p_measurement]
    for key in (:pst, :wpst, :mg)
      base = nets[key]
      base === nothing && continue
      tr = _lad_truth(base, setc_moves[key], false)
      mesh = setc_mesh[key]
      mbus = Set(b for k in mesh for b in (Int(tr.net.branchVec[k].fromBus), Int(tr.net.branchVec[k].toBus)))
      setc(kept) = filter(m -> !(m.typ in (S.PflowMeas, S.ImagMeas) && m.branchIdx in mesh) && !(m.typ == S.PinjMeas && m.busIdx in setdiff(mbus, kept)), tr.meas)
      # the warning fires on set C (MicroGrid BE as well) and names the mesh
      @test meshbranches(_lad_with(() -> validate_topology(tr.net, setc(Int[])).findings, base.config)) == [mesh]
      # the full set has a P row in every mesh
      @test isempty(meshbranches(_lad_with(() -> validate_topology(tr.net, tr.meas).findings, base.config)))
      key === :mg && continue
      C = setc(Int[])
      res = _lad_se(_lad_release!(deepcopy(base.net), setc_moves[key]), C, base.config)
      @test res.converged && meshbranches(res.topologyFindings) == [mesh]
      # the changer row of br8 (in the mesh) carries the mesh notice with
      # status and reason unchanged, br7 (outside) does not
      meshnote(t) = occursin("mesh without P measurement: the step may not be unique", t.notice)
      rows = [t for t in res.tapEstimates if haskey(tr.steps, (t.branch, t.regulator))]
      @test length(rows) == 2 && all(t -> meshnote(t) == (t.branch in mesh), rows)
      @test all(t -> t.status === :estimated && t.reason === :none, rows)
      # regression case: P at bus 2 kept, no warning, both changers at the
      # true step (exact data, J after the fixation at round-off level)
      res = _lad_se(_lad_release!(deepcopy(base.net), setc_moves[key]), setc([2]), base.config)
      rows = [t for t in res.tapEstimates if haskey(tr.steps, (t.branch, t.regulator))]
      @test isempty(meshbranches(res.topologyFindings))
      @test !any(meshnote, res.tapEstimates)
      @test res.converged && length(rows) == 2 && res.tapFixation.j_after < 1e-8
      @test all(t -> t.status === :estimated && t.step_fixed == tr.steps[(t.branch, t.regulator)], rows)
      # documented case, set C on the truth net, nothing released. Started
      # at the truth (the power-flow state the truth net holds) the WLS stays
      # there: exact data, J 9.3e-26. From the flat start it converges to a
      # second solution with the angle differences on br1 (truth +0.86 deg)
      # and br3 (+1.50 deg) turned over (-1.08 and -0.96 deg), J 15.42 at
      # dof 33 (flat start, tol 1e-10).
      for fs in (false, true)
        net = deepcopy(tr.net)
        r = _lad_with(base.config) do
          Logging.with_logger(Logging.NullLogger()) do
            with_state_estimation_config(() -> runse!(net, C); max_iter = 50, tol = 1e-10, update_net = false, flatstart = fs)
          end
        end
        dva(k) = (b = tr.net.branchVec[k]; (angle(r.voltages[Int(b.fromBus)]) - angle(r.voltages[Int(b.toBus)]), deg2rad(tr.net.nodeVec[Int(b.fromBus)]._va_deg - tr.net.nodeVec[Int(b.toBus)]._va_deg)))
        @test r.converged && r.dof == 33
        if fs
          @test isapprox(r.objectiveJ, 15.42; atol = 0.05)
          @test all(k -> (d = dva(k); sign(d[1]) == -sign(d[2])), (1, 3))
        else
          @test r.objectiveJ < 1e-20
        end
      end
    end
    # the conditional ratio changer of the pair a + b60 without Va: the
    # magnitude of the frozen b60 (truth -4) goes into the ratio step
    let base = nets[:n3ab60], tr = _lad_truth(base, Dict(2 => [4.0, -4.0]), false)
      row = first(t for t in _lad_se(_lad_release!(deepcopy(base.net), Dict(2 => [4.0, -4.0])), tr.meas, nothing).tapEstimates if t.regulator == 1)
      @test row.note === :conditional && occursin("frozen phase position", row.notice)
    end
    # The partner of a changer the recheck withdraws is estimated again with
    # it at its transferred step. 3-bus net, ratio changer a and a
    # quadrature changer q (2 % per step) whose reactance follows its angle
    # (x_curve, factor 1 to 2), so the pair is separable without an angle
    # row; truth (2, 4), start (0, 4), no Va, noisy rows (seed 4). The joint
    # estimate gives a -1.56 and q 14.86 (fixed at the band end 10); the
    # mirror step -10 fits as well, so q is withdrawn to its transferred 4.
    # Until 0.60.0 a kept the -2 of the joint estimate (sigma_s 3.43), J
    # 98.39 against 15.21 at the truth. Estimated again with q at 4 it is
    # 2.41 with sigma_s 0.5, fixed at the truth 2, and the state is the
    # state of the net at (2, 4) without a release
    let
      S = Sparlectra
      pq(r, q) = _lad_net3([S.SETapChanger(name = "a", kind = :ratio, role = :ratio, side = :from, du = 0.0125, neutral = 0, low = -10, high = 10, step = r), S.SETapChanger(name = "q", kind = :asymmetrical, role = :phase, side = :from, du = 0.02, alpha_deg = 90.0, neutral = 0, low = -10, high = 10, step = q, x_curve = (1.0, 2.0, atand(0.2)), x_curve_form = :tan2, impedance = :x_curve)])
      tn = pq(2, 4)
      @test _lad_pf!(tn, nothing) == 0
      meas = generateMeasurementsFromPF(tn; includeImag = true, noise = true, rng = Random.MersenneTwister(4))
      res = _lad_se(_lad_release!(pq(0, 4), Dict(2 => [2.0, 0.0])), meas, nothing)
      ra = only(t for t in res.tapEstimates if t.changer == "a")
      rq = only(t for t in res.tapEstimates if t.changer == "q")
      @test res.converged && rq.status === :not_estimated && rq.reason === :mirror_step_fits && rq.step_fixed == 4.0
      @test ra.status === :estimated && ra.step_fixed == 2.0 && occursin("changer q of this transformer is not determinable (mirror_step_fits), this estimate depends on its transferred step 4", ra.notice)
      # the row of a is the estimate with q held at 4: the run that releases
      # a alone gives the same value and sigma_s
      alone = only(_lad_se(_lad_release!(pq(0, 4), Dict(2 => [2.0, NaN])), meas, nothing).tapEstimates)
      @test isapprox(ra.s_est, alone.s_est; atol = 1e-6) && isapprox(ra.sigma_s, alone.sigma_s; rtol = 1e-4)
      ref = _lad_se(pq(2, 4), meas, nothing)
      @test isapprox(res.voltages, ref.voltages; atol = 1e-6) && isapprox(res.objectiveJ, ref.objectiveJ; rtol = 1e-9)
      # A partner with status out_of_range: ratio changer a (band -10 to
      # 10, the truth net -15 to 15) at the truth 13 and a symmetrical
      # changer d at 3, start (0, 0), exact rows with Va. The estimate of a
      # is 13 with sigma_s 0.6, 5 sigma_s beyond the band end, so a keeps
      # its transferred 0. Until 0.60.0 d kept the joint estimate 3.000 and
      # its row had no notice; now d is estimated again with a at 0 (3.126,
      # fixed at 3) and its notice names a
      ad(r, d; lo = -10, hi = 10) = _lad_net3([S.SETapChanger(name = "a", kind = :ratio, role = :ratio, side = :from, du = 0.0125, neutral = 0, low = lo, high = hi, step = r), _lad_changer(:d)])
      tn = ad(13, 0; lo = -15, hi = 15)
      _lad_move!(tn, Dict(2 => [NaN, 3.0]), :auto)
      @test _lad_pf!(tn, nothing) == 0
      meas = generateMeasurementsFromPF(tn; includeImag = true, includeVa = true)
      res = _lad_se(_lad_release!(ad(0, 0), Dict(2 => [13.0, 3.0])), meas, nothing)
      ra = only(t for t in res.tapEstimates if t.changer == "a")
      rd = only(t for t in res.tapEstimates if t.changer == "d")
      @test res.converged && ra.status === :out_of_range && ra.step_fixed == 0.0 && isapprox(ra.s_est, 13.0; atol = 1e-6)
      @test rd.status === :estimated && rd.step_fixed == 3.0 && occursin("changer a of this transformer is not determinable (out_of_range), this estimate depends on its transferred step 0", rd.notice)
      alone = only(_lad_se(_lad_release!(ad(0, 0), Dict(2 => [NaN, 3.0])), meas, nothing).tapEstimates)
      @test isapprox(rd.s_est, alone.s_est; atol = 1e-6) && isapprox(rd.sigma_s, alone.sigma_s; rtol = 1e-4)
      dn = ad(0, 0)
      _lad_move!(dn, Dict(2 => [NaN, 3.0]), :auto)
      ref = _lad_se(dn, meas, nothing)
      @test isapprox(res.voltages, ref.voltages; atol = 1e-6) && isapprox(res.objectiveJ, ref.objectiveJ; rtol = 1e-9)
    end
  end)() end
  return true
end

# erf(x) by the series with positive terms
# erf(x) = 2/sqrt(pi) exp(-x^2) sum_n 2^n x^(2n+1) / (2n+1)!!, no
# cancellation; 1 above x = 6, where 1 - erf(x) < 3e-17. Gives the share
# of a normal estimate that rounds onto the true step (stage 3).
function _lad_erf(x::Float64)::Float64
  x < 0.0 && return -_lad_erf(-x)
  x > 6.0 && return 1.0
  term = x
  s = term
  n = 0
  while term > 1e-17 * s
    n += 1
    term *= 2.0 * x^2 / (2n + 1)
    s += term
  end
  return 2.0 / sqrt(pi) * exp(-x^2) * s
end

# The q quantile of Binomial(n, p): the smallest k with P(X <= k) >= q,
# summed from k = 0 in log space ((1 - p)^n underflows for p near 1).
# Hand check: n = 100, p = 0.5, q = 0.001 gives 35 (P(X <= 34) = 8.9e-4,
# P(X <= 35) = 1.8e-3).
function _lad_binom_quantile(n::Int, p::Float64, q::Float64)::Int
  p >= 1.0 && return n
  p <= 0.0 && return 0
  lp = n * log1p(-p)
  c = 0.0
  for k in 0:n
    k > 0 && (lp += log((n - k + 1) / k) + log(p) - log1p(-p))
    c += exp(lp)
    c >= q && return k
  end
  return n
end

# Stage 3 counts of one ladder case: 100 draws with the sigmas of
# measurementStdDevs(), draw d with the fixed seed MersenneTwister(d) for d
# = 1 to 100, so the gate set gives the same counts on every run. Per
# released changer (branch, changer), over the draws that estimated it
# (status :estimated): their number, the draws with the estimate within two
# sigma_s of the true step, the rounding hits of the true step, the expected
# share of rounding hits (mean of 2 Phi(0.5 / sigma_i) - 1 = erf(0.5 /
# (sqrt(2) sigma_i)) over these draws, sigma_i the sigma_s of draw i), the
# mean sigma_s, the spread (sample standard deviation) of the estimates and
# their mean deviation from the true step; plus the draws that did not
# estimate the changer with their reasons (the fallback of runse! among
# them) and the number of converged runs.
function _lad_noise_counts(base, c)
  tr = _lad_truth(base, c.moves, c.va)
  nest = Dict(k => 0 for k in keys(tr.steps))
  hits = Dict(k => 0 for k in keys(tr.steps))
  in2 = Dict(k => 0 for k in keys(tr.steps))
  pexp = Dict(k => 0.0 for k in keys(tr.steps))
  sig = Dict(k => 0.0 for k in keys(tr.steps))
  est = Dict(k => Float64[] for k in keys(tr.steps))
  notest = Dict(k => Symbol[] for k in keys(tr.steps))
  nconv = 0
  for d in 1:100
    meas = generateMeasurementsFromPF(tr.net; includeImag = true, includeVa = c.va, noise = true, rng = Random.MersenneTwister(d))
    res = _lad_se(_lad_release!(deepcopy(base.net), c.moves), meas, base.config; tol = 1e-6)
    nconv += res.converged
    for t in res.tapEstimates
      key = (t.branch, t.regulator)
      haskey(tr.steps, key) || continue
      if t.status !== :estimated
        push!(notest[key], t.reason)
        continue
      end
      s = t.s_est + t.neutral_step
      nest[key] += 1
      hits[key] += t.step_fixed == tr.steps[key]
      # an estimated row without sigma_s counts as outside two sigma and
      # adds 0 to the expected rounding share
      in2[key] += t.sigma_s !== nothing && abs(s - tr.steps[key]) <= 2.0 * t.sigma_s
      pexp[key] += t.sigma_s === nothing ? 0.0 : _lad_erf(0.5 / (sqrt(2.0) * t.sigma_s))
      sig[key] += something(t.sigma_s, NaN)
      push!(est[key], s)
    end
  end
  function cell(k)
    n = nest[k]
    m = n == 0 ? NaN : sum(est[k]) / n
    return (nest = n, hits = hits[k], in2 = in2[k], pexp = n == 0 ? NaN : pexp[k] / n, sigma_s = n == 0 ? NaN : sig[k] / n, spread = n < 2 ? NaN : sqrt(sum(abs2, est[k] .- m) / (n - 1)), bias = m - tr.steps[k], notest = notest[k])
  end
  return (nconv = nconv, cells = Dict(k => cell(k) for k in keys(tr.steps)))
end

function test_tap_estimation_ladder_noise()::Bool
  # Stage 3 of the verification ladder: noisy measurements with the sigmas
  # a measurement row without its own sigma gets (measurementStdDevs():
  # Vm 0.005 pu, P and Q 1 MW and 1 Mvar, |I| 10 A, Va 0.02 deg), 100 draws
  # per case of _lad_cases with fixed seeds (_lad_noise_counts). Every run
  # converges. Per released changer, counted over the draws that estimated
  # it (n of 100; the draws without an estimate and their reasons are
  # printed, not asserted, apart from n > 0):
  # - coverage: the estimate within two sigma_s of the true step in at
  #   least 89 % of the n draws (89 of 100). A correct estimator lands
  #   there with p = 0.954 per draw and falls below 89 of 100 with about
  #   0.2 % probability per changer.
  # - rounding: the hits of the true step not below the 0.1 % quantile of
  #   Binomial(n, p_exp), p_exp the mean of 2 Phi(0.5 / sigma_i) - 1 over
  #   the n draws. A correct estimator fails this with at most 0.1 %
  #   probability per changer.
  # Changed after the first result (see docs/src/tests.md). Before, only
  # the rounding hits were asserted, at
  # least 95 of 100; the two-sigma count was printed and not asserted. A
  # normal estimate rounds onto its step in 95 % of the draws only for
  # sigma_s below 0.255 steps, so the fixed 95 is replaced by the binomial
  # bound. The two-sigma count is a new assertion at 89 of 100, not 95: a
  # correct covariance puts 95 of 100 draws within two sigma only with
  # about 70 % probability. The sigmas stay as they are. Every cell prints
  # one line with its counts.
  @testset "Tap estimation ladder, noisy measurements" begin (function ()
    nets = _lad_nets()
    for c in _lad_cases()
      base = nets[c.net]
      base === nothing && continue
      r = _lad_noise_counts(base, c)
      @test r.nconv == 100
      for key in sort(collect(keys(r.cells)))
        x = r.cells[key]
        need2 = ceil(Int, 0.89 * x.nest)
        lo = x.nest == 0 ? 0 : _lad_binom_quantile(x.nest, x.pexp, 0.001)
        ok = x.nest > 0 && x.in2 >= need2 && x.hits >= lo
        @test x.nest > 0
        @test x.in2 >= need2
        @test x.hits >= lo
        reasons = join(sort(unique(string.(x.notest))), " ")
        println("      ladder stage 3 ", c.name, " branch ", key[1], " changer ", key[2], ": ", ok ? "passed" : "FAILED", ", estimated ", x.nest, " of 100", isempty(x.notest) ? "" : string(" (not estimated ", length(x.notest), ": ", reasons, ")"), ", within 2 sigma ", x.in2, " (need ", need2, "), rounding hits ", x.hits, " (expected ", round(x.nest * x.pexp; digits = 1), ", need ", lo, "), mean sigma_s ", round(x.sigma_s; sigdigits = 3), " steps, spread ", round(x.spread; sigdigits = 3), ", mean deviation ", round(x.bias; sigdigits = 3))
      end
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
        SparlectraApp._run_state_estimation_service(z, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, out, name, set; max_iter = 50, tap_staged_start = staged)
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

# The reported block of se_diagnostics.md describes the run without tap
# estimation the status line reports after the fallback: same run id, case,
# set, iterations, J and dof as run.log and the metadata, and the two
# fallback lines in run.log. `case` is the case name the block shows.
function _assert_tap_fallback_report(dd, log::AbstractString, diag::AbstractString, name::AbstractString, set::AbstractString; case::AbstractString = "cgmes_minigrid.zip", reason::AbstractString = "not_converged")
  @test dd["status"] == "succeeded"
  @test get(dd["metadata"], "se_tap_fallback_used", false) == true
  # line 1: a failed run 1, or the acceptance rule rejecting a converged one
  line1 = reason == "worse_fit" ? r"tap estimation: run 1 rejected \(worse_fit, J = [0-9.eE+-]+ at dof \d+ after the fixation against J = [0-9.eE+-]+ at dof \d+ without tap estimation\), falling back" : r"tap estimation: run 1 not converged \((estimation|fixation)_not_converged, \d+ iterations\), falling back"
  @test occursin(line1, log)
  jm = match(r"converged in (\d+) iteration\(s\); J = ([0-9.eE+-]+), dof = (\d+)", log)
  @test jm !== nothing
  jm === nothing && return
  @test occursin("Tap estimation fallback", diag)
  @test occursin("## Reported result", diag)
  @test occursin("**Converged:** true in $(jm.captures[1]) iteration(s)", diag)
  @test occursin("**Objective J:** $(jm.captures[2]) (dof=$(jm.captures[3])), without tap estimation (fallback)", diag)
  @test occursin("run 2 without tap estimation: converged, J = $(round(dd["metadata"]["se_objective"]; sigdigits = 6)), dof = $(jm.captures[3])", log)
  # every changer row of the CSV says why it was not estimated; a machine
  # transformer row is calculated, or not_estimated with
  # q_not_measured_in_group when its machine node lacks a Q row (0.60.0,
  # before the step came from the scheduled Q)
  rows = [t for t in dd["metadata"]["se_tap_estimates"] if t["status"] != "calculated" && t["reason"] != "q_not_measured_in_group"]
  @test !isempty(rows) && all(t -> t["status"] == "not_estimated" && t["reason"] == (reason == "worse_fit" ? "fallback_worse_fit" : "fallback_not_converged"), rows)
  @test occursin("- **Run:** $(name)", diag)
  @test occursin("- **Case:** $(case)", diag)
  @test occursin("- **Measurement set:** $(basename(set))", diag)
  @test occursin("- **J/dof:** $(round(parse(Float64, jm.captures[2]) / parse(Int, jm.captures[3]); digits = 2))", diag)
  @test occursin("## Diagnostics pass (every row active, before any elimination)", diag)
  # the diagnostics pass printed below the block describes the same run
  @test !occursin("**Converged:** false", diag)
  @test isapprox(dd["metadata"]["se_objective"], parse(Float64, jm.captures[2]); atol = 1e-6)   # run.log rounds to six digits
  @test dd["metadata"]["se_dof"] == parse(Int, jm.captures[3])
  # the verdict file counts what the run did: every row a diagnostics pass
  # removed (the pass with tap estimation and the repeat without it) is
  # eliminated, inactive only where the table says valid=false,
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
released-tap attempt does not converge and the repeat without tap
estimation does. The artifact used to describe the failed attempt (the
diverged J at the dof with tap states, the other run's J on the page); now
the diagnostics repeat without tap estimation and the artifact opens with
the reported numbers. This set confirms
all five eliminations of the default budget and stops at
`max_eliminations` (three until the budget became 5 with the range
check). On the set
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
    @test occursin("fallback without tap estimation", SparlectraApp._se_truth_run_warning(Dict{String,Any}("se_band_reason" => "ok", "se_tap_fallback_used" => true)))
    @test SparlectraApp._se_truth_run_warning(Dict{String,Any}("se_band_reason" => "ok", "numerical_status" => "converged")) == ""
    set = joinpath(@__DIR__, "fixtures", "measurements", "cgmes_minigrid_tap_fallback.measurements.csv")
    # the COLD mass release (no staged start) does not converge here and
    # falls back; the fixture pins that path's report
    _minigrid_fallback_run(set, "minigrid_fb") do dd, log, diag
      _assert_tap_fallback_report(dd, log, diag, "minigrid_fb", set)
      # the default budget state_estimation.max_eliminations = 5 is spent
      # (3 before the range check; with 3 the run gives J 10752.8 at
      # dof 75, with 5 J 6225.1 at dof 73)
      @test occursin("- **Eliminated rows:** 5 ", diag)
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
    # The fixation of each solve here (the estimate and the
    # eliminations) raises J by about 9 against the 95 % quantile 5.99 of
    # chi2(2) on this set (J/dof 140 to 620), so each warns
    _minigrid_fallback_run(set, "minigrid_fb_staged"; staged = true, expected = ["above the 95 % quantile 5.991"]) do dd, log, _
      @test dd["status"] == "succeeded"
      @test get(dd["metadata"], "se_tap_staged_start", false) == true
      @test get(dd["metadata"], "se_tap_fallback_used", false) == false
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
      @test get(dd["metadata"], "se_tap_fallback_used", false) == false
      @test get(dd["metadata"], "se_tap_j_rise_exceeds_quantile", false) == true && get(dd["metadata"], "se_tap_sequential_fixation", false) == true
      @test occursin("kept the sequential result", log)
      legs = [t for t in dd["metadata"]["se_tap_estimates"] if get(t, "status", "") == "estimated"]
      @test length(legs) == 2 && all(t -> t["fixed"] == true && abs(t["deviation_steps"]) <= 1, legs) && legs[1]["step_estimated"] == legs[2]["step_estimated"]
      # J of the kept pair against the solve at the model steps of the staged start
      jm = match(r"staged start, one solve with the 2 released tap\(s\) frozen \(\d+ iteration\(s\), J = ([0-9.eE+-]+)", log)
      @test jm !== nothing && dd["metadata"]["se_tap_j_after"] <= parse(Float64, jm.captures[1]) * (1.0 + 1e-6)
    end
  end
  return true
end

"""
The fixture set `cgmes_minigrid_elimination_reverted` (same case, same
truth source as the set above, a later generation): the cold release falls
back, and the repeat without tap estimation from the start values of run 1
confirms all eliminations of the budget (five since the range check step,
three before). Until 0.60.0 the repeat started from the
converged frozen-tap state and its third control solve did not converge,
so that removal was reverted. The staged start keeps its estimate. Needs
the conformity package in the local cache: RAN or SKIPPED is printed.
"""
function test_state_estimation_elimination_reverted()::Bool
  @testset "Tap fallback: eliminations of the repeat (MiniGrid)" begin
    set = joinpath(@__DIR__, "fixtures", "measurements", "cgmes_minigrid_elimination_reverted.measurements.csv")
    # The cold release falls back. Until 0.60.0 the repeated diagnostics
    # started from the state the converged frozen-tap run had written into
    # the net, and there the third control solve did not converge (2
    # eliminations, P@7fe566b9 reverted, stop not_converged, J 23017.68 at
    # dof 76). The fallback now repeats them from the start values of run 1,
    # as the run without tap estimation does: all three eliminations of the
    # former budget 3 were confirmed, J 19176.97 at dof 75, the same as the run with the switch
    # off. The revert of an unconfirmed elimination
    # is covered in the determinability guard set of legacy_state_estimation.
    _minigrid_fallback_run(set, "minigrid_rev") do dd, log, diag
      _assert_tap_fallback_report(dd, log, diag, "minigrid_rev", set)
      # the set has no Q row at HG1, the machine node of B_2WT_21_11_2: not
      # calculated, the transferred step stays
      @test [(t["name"], t["status"]) for t in dd["metadata"]["se_tap_estimates"] if t["reason"] == "q_not_measured_in_group"] == [("B_2WT_21_11_2", "not_estimated")]
      @test occursin("machine transformer B_2WT_21_11_2: tap back-calculation not possible", log)
      # budget 5 since the range check: five confirmed eliminations,
      # J 4540.96 at dof 73 (budget 3 gives the J 19176.97 at dof 75
      # above)
      @test occursin("- **Eliminated rows:** 5 ", diag)
      @test !occursin("- **Not eliminated:** `", diag)
      @test occursin("(stop: max_eliminations)", diag)
      @test !occursin("Reverted:", diag)
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
      @test get(dd["metadata"], "se_tap_fallback_used", false) == false
      @test occursin("tap estimation: staged start", log)
      @test !occursin("falling back", log)
    end
  end
  return true
end


"""
Observability per changer and the step from the control target (for the
automatically excluded transformers), on the gate net with a radial transformer (bus L3
on L2 over branch 5, ratio changer of 0.625 %, the truth in steps). With the
Vm row of L3 the step is estimated exactly and carries `sigma_s` from the
covariance, the rounded step lies 0 sigma away. When only the from-end
flows of branch 5 see the far side (no row at L3, no to-end flow), the
step column lies in the span of the L3 voltage states (rank test, reason
`radial_no_voltage_pin`; with the injection and to-end rows of L3 the
losses would tell the step apart): the changer keeps its transferred
step and, with a voltage controller on L3 whose target is the true
voltage, `step_from_target` gives the true step after a run with
`update_net`: the transferred step equal to the truth, no notice; three
steps off, the notice says so; a phase controller on the same
transformer ahead of the voltage controller does not hide it; without a
controller the column is empty.
"""
function test_state_estimation_tap_step_from_target()::Bool
  @testset "Tap step from the control target" begin (function ()
    step = 0.00625
    function radial(n::Int)
      net = _tap_gate_net()
      addBus!(net = net, busName = "L3", vn_kV = 10.0)
      addProsumer!(net = net, busName = "L3", type = "ENERGYCONSUMER", p = 8.0, q = 2.0)
      addPIModelTrafo!(net = net, fromBus = "L2", toBus = "L3", r_pu = 0.004, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1)
      br = net.branchVec[5]
      br.has_ratio_tap = true
      br.ratio_tap_declared = true
      br.tap_step = step
      br.tap_ratio = br.ratio / (1.0 + n * step)
      return net
    end
    truth = radial(4)
    meas = _tap_truth_meas(truth)
    l3 = 5
    vtrue = truth.nodeVec[l3]._vm_pu
    run(net, ms) = with_state_estimation_config(() -> runse!(net, ms); max_iter = 40, tol = 1e-10, update_net = true)
    # with the Vm row of L3: estimated, sigma_s from the covariance
    net = radial(0)
    setTapEstimation!(net; trafo = 5, mode = :ratio)
    t = only(run(net, meas).tapEstimates)
    @test t.status === :estimated && t.step_fixed == 4 && t.step_from_target === nothing
    @test t.sigma_s isa Float64 && 0.0 < t.sigma_s < 1.0 && t.round_distance_sigma < 1e-3
    far_side = [m for m in meas if m.busIdx != l3 && !(m.branchIdx == 5 && m.direction == :to)]
    # (1, :both, ...): a phase controller on the same transformer comes
    # first; the ratio changer still takes the voltage controller (the
    # first controller by name gave an empty column before)
    for (start, ctrl, target, far) in ((4, :ratio, 4.0, false), (1, :ratio, 4.0, true), (1, :both, 4.0, true), (1, :none, nothing, false))
      net = radial(start)
      name5 = getCompName(net.branchVec[5].comp)
      ctrl === :both && addPowerTransformerControl!(net; trafo = name5, mode = :branch_active_power, target_branch = ("L2", "L3"), p_target_mw = 8.0, control_ratio = false, control_phase = true)
      ctrl === :none || addPowerTransformerControl!(net; trafo = name5, mode = :voltage, target_bus = "L3", target_vm_pu = vtrue)
      setTapEstimation!(net; trafo = 5, mode = :ratio)
      res = run(net, far_side)
      t = only(res.tapEstimates)
      @test res.converged && t.status === :not_estimated && t.reason === :radial_no_voltage_pin && abs(t.step_model - start) < 1e-9 && t.sigma_s === nothing
      @test t.step_from_target == target
      @test occursin("more than one step", t.notice) == far
    end
  end)() end
  return true
end

"""
Tap estimation switch and fallback (state_estimation.tap_estimation.enabled
and .fallback). Released taps are extra states, and a measurement set that
estimates the voltages cleanly can still be too thin to pin them. A run
with tap estimation that does not converge is therefore solved once more
without it, on the same measurements and start values, and that second run
is the result; run.log, the metadata and every tap row say so. The SE
service falls back the same way when the fixed estimate of a converged run
fits worse than its staged pre-solve (`worse_fit`). Switched off, no
changer becomes a state and every row says `disabled`.
"""
function test_state_estimation_tap_fallback()::Bool
  @testset "State estimation tap fallback" begin
    # The library path on the shipped sp_case5: the bridge transformer br5
    # gets a Delta-u changer of 5 % at psi 60 (truth step -6, start 0), no
    # Va, and the Vm row of its far bus a gross error of -0.05 pu. Pinned
    # only by its magnitude effect, the released step runs into the turning
    # point of |t| at s du = -0.5 and the run with tap estimation does not
    # converge in 40 iterations (measured with the refreeze of c758bc91,
    # same case and set)
    begin
      S = Sparlectra
      case5 = joinpath(dirname(@__DIR__), "data", "scf", "sp_case5.scf.json")
      m5 = redirect_stdout(() -> S.import_case(case5, S.resolve_config(S.DEFAULT_SPARLECTRA_CONFIG_PATH, case5).config; run_kind = :state_estimation), devnull)
      empty!(m5.net.measurements)
      du5(n) = (tn = deepcopy(m5.net); applyTapNameplate!(tn.branchVec[5]; phase_du_step = 0.05, phase_min_step = -10.0, phase_max_step = 10.0, phase_current_step = Float64(n), psi_deg = 60.0); tn)
      tn = du5(-6)
      _, erg = redirect_stdout(() -> runpf!(tn; config = S._copy_powerflow_with(m5.config.powerflow; tol = 1e-10, max_iter = 40)), devnull)
      @test erg == 0
      lv = Int(tn.branchVec[5].toBus)
      meas = [m.typ == S.VmMeas && m.busIdx == lv ? Measurement(typ = S.VmMeas, value = m.value - 0.05, sigma = m.sigma, busIdx = lv, id = m.id * "_gross") : m for m in generateMeasurementsFromPF(tn; includeImag = true, noise = true, rng = Random.MersenneTwister(1))]
      run5(n; enabled = true, fallback = true, update_taps = false) = S.with_sparlectra_config(m5.config) do
        with_state_estimation_config(() -> runse!(n, meas); max_iter = 40, tol = 1e-6, update_net = false, update_taps = update_taps, tap_estimation = TapEstimationConfig(enabled = enabled, fallback = fallback))
      end
      released() = (n = du5(0); setTapEstimation!(n; trafo = 5, mode = :pst); n)
      ref = run5(du5(0))                       # no release at all
      # switched off: no tap state (no tapFixation, the dof of the run
      # without a release), the row says disabled, update_taps writes nothing
      noff = released()
      shift0 = noff.branchVec[5].phase_shift_deg
      resOff = run5(noff; enabled = false, update_taps = true)
      rowOff = only(t for t in resOff.tapEstimates if t.role === :phase)
      @test resOff.converged && resOff.tapFixation === nothing && resOff.dof == ref.dof
      @test rowOff.status === :not_estimated && rowOff.reason === :disabled && rowOff.s_est == 0.0 && noff.branchVec[5].phase_shift_deg == shift0
      # fallback forced: both lines, the flag, and the result of the run
      # with the switch off
      res = @test_logs (:warn, r"^tap estimation: run 1 not converged \(estimation_not_converged, 40 iterations\), falling back$") (:info, r"^run 2 without tap estimation: converged, J = [0-9.e+-]+, dof = \d+$") match_mode = :any run5(released())
      row = only(t for t in res.tapEstimates if t.role === :phase)
      @test res.converged && res.tapFixation.fallback_used && res.tapFixation.fallback_reason === :estimation_not_converged && res.tapFixation.iterations_estimation == 40
      @test row.status === :not_estimated && row.reason === :fallback_not_converged && !row.fixed && occursin("did not converge", row.notice)
      @test res.dof == resOff.dof && res.objectiveJ == resOff.objectiveJ && res.voltages == resOff.voltages
      # fallback off: the first run is the result, not converged
      resNo = run5(released(); fallback = false)
      @test !resNo.converged && !resNo.tapFixation.fallback_used && !resNo.tapFixation.fixed
      # Estimate outside the band (verification ladder, finding 1): truth
      # and start step -8, exact measurements. Pinned only by |t|, the
      # estimate goes to the mirror point across the turning point of |t|
      # at s* = -cos(60)/0.05 = -10, so s = 2 s* - (-8) = -12 gives the same
      # |t| (|t|^2 = (s du + cos psi)^2 + sin(psi)^2). It lies within 2
      # sigma_s of the band end (else the status would be out_of_range), so
      # the fixation takes the end step -10, the turning point itself; the
      # row notice and one runse!
      # warning name the band and the mirror step -8. At -10 the step column
      # has no rank (d|t|/ds = 0, the angle has no anchor), so the recheck at
      # the estimate (0.60.0) withdraws the changer
      # (:not_observable_at_estimate) and the transferred step -8, the
      # truth, stays with J of the truth. Before the recheck it stayed fixed at
      # -10 with status estimated
      tm = du5(-8)
      redirect_stdout(() -> runpf!(tm; config = S._copy_powerflow_with(m5.config.powerflow; tol = 1e-10, max_iter = 40)), devnull)
      measM = generateMeasurementsFromPF(tm; includeImag = true)
      nm = du5(-8)
      setTapEstimation!(nm; trafo = 5, mode = :pst)
      oorText = r"estimate at step -12 lies outside the band -10 to 10, \|t\| has the same value at the mirror step -8 across its turning point at step -10"
      resM = @test_logs (:warn, r"^runse!: tap regulator .* estimated outside its band: estimate at step -12 ") match_mode = :any S.with_sparlectra_config(m5.config) do
        with_state_estimation_config(() -> runse!(nm, measM); max_iter = 40, tol = 1e-6, update_net = false, tap_estimation = TapEstimationConfig(fallback = false))
      end
      rowM = only(t for t in resM.tapEstimates if t.role === :phase)
      @test rowM.out_of_range && rowM.status === :not_estimated && rowM.reason === :not_observable_at_estimate && rowM.step_fixed == -8.0 && isapprox(rowM.s_est, -12.0; atol = 1e-6)
      @test occursin(oorText, rowM.notice) && occursin(oorText, S._tap_out_of_range_line(rowM)) && resM.objectiveJ < 1e-8
      # Estimate outside the band, the two verdicts of the fixation, on a
      # ratio changer of 1.25 % on br5 whose model band -10 to 10 is too
      # small for the truth (the truth net carries -15 to 15), start 0,
      # exact measurements. The estimate is the truth; sigma_s scales with
      # the row sigmas. Truth 11 with the row sigmas times 5: sigma_s 2.54,
      # 1 step beyond the band end is 0.39 sigma_s, within 2 sigma_s, so the
      # fixation takes the end step 10, status estimated with the notice.
      # Truth 13 with the default sigmas: sigma_s 0.51, 3 steps are 5.9
      # sigma_s, so the estimate is not taken over: status out_of_range,
      # the transferred step 0 stays, not fixed, and the J of the fixed run
      # rises (warned)
      rat5(n, lo, hi) = (tn = deepcopy(m5.net); applyTapNameplate!(tn.branchVec[5]; tap_step = 0.0125, tap_min_step = lo, tap_max_step = hi, tap_current_step = Float64(n)); tn)
      for (truth, f, verdict) in ((11, 5.0, :band_end), (13, 1.0, :out_of_range))
        tr = rat5(truth, -15.0, 15.0)
        redirect_stdout(() -> runpf!(tr; config = S._copy_powerflow_with(m5.config.powerflow; tol = 1e-10, max_iter = 40)), devnull)
        measR = [Measurement(typ = m.typ, value = m.value, sigma = m.sigma * f, busIdx = m.busIdx, branchIdx = m.branchIdx, direction = m.direction, id = m.id) for m in generateMeasurementsFromPF(tr; includeImag = true)]
        nr = rat5(0, -10.0, 10.0)
        setTapEstimation!(nr; trafo = 5, mode = :ratio)
        resR = @test_logs (:warn, Regex("^runse!: tap regulator .* estimated outside its band: estimate at step $(truth) lies outside the band -10 to 10, ")) match_mode = :any S.with_sparlectra_config(m5.config) do
          with_state_estimation_config(() -> runse!(nr, measR); max_iter = 40, tol = 1e-6, update_net = false, tap_estimation = TapEstimationConfig(fallback = false))
        end
        rowR = only(resR.tapEstimates)
        @test resR.converged && rowR.out_of_range && isapprox(rowR.s_est, Float64(truth); atol = 1e-6) && rowR.sigma_s !== nothing
        if verdict === :band_end
          @test rowR.status === :estimated && rowR.fixed && rowR.step_fixed == 10.0 && 1.0 / rowR.sigma_s <= 2.0
          @test occursin(r"lies outside the band -10 to 10, 0\.394 sigma_s beyond the band end step 10, fixed to the end step 10", rowR.notice)
        else
          @test rowR.status === :out_of_range && rowR.reason === :none && !rowR.fixed && rowR.step_fixed == 0.0 && 3.0 / rowR.sigma_s > 2.0
          @test occursin(r"lies outside the band -10 to 10, 5\.9 sigma_s beyond the band end step 10, not taken over, the transferred step 0 stays", rowR.notice)
          @test resR.tapFixation.j_rise_exceeds_quantile
        end
      end
      # Recheck at the estimate, the other two paths on the same bridge br5.
      # Mirror pair: a Delta-u
      # changer at psi 90 (1 % per step, turning point of |t| at neutral),
      # truth step 3, start 2, exact measurements, no angle row. The
      # estimate is the truth, but |t| at the mirror step -3 is the same:
      # the fixed run at -3 gives J at rounding level as well, the difference
      # is below 9, the changer is withdrawn (:mirror_step_fits) and keeps
      # its transferred step 2, the notice names both J
      # (J 2.7e-26 at 3, 5.4e-25 at -3)
      d90(n; du = 0.01) = (tn = deepcopy(m5.net); applyTapNameplate!(tn.branchVec[5]; phase_du_step = du, phase_min_step = -10.0, phase_max_step = 10.0, phase_current_step = Float64(n), psi_deg = 90.0); tn)
      se5(n, ms) = S.with_sparlectra_config(m5.config) do
        with_state_estimation_config(() -> runse!(n, ms); max_iter = 40, tol = 1e-6, update_net = false)
      end
      tm = d90(3)
      redirect_stdout(() -> runpf!(tm; config = S._copy_powerflow_with(m5.config.powerflow; tol = 1e-10, max_iter = 40)), devnull)
      nm = d90(2)
      setTapEstimation!(nm; trafo = 5, mode = :pst)
      resA = se5(nm, generateMeasurementsFromPF(tm; includeImag = true))
      rowA = only(resA.tapEstimates)
      @test resA.converged && rowA.status === :not_estimated && rowA.reason === :mirror_step_fits && rowA.step_fixed == 2.0 && isapprox(rowA.s_est, 3.0; atol = 1e-6)
      @test occursin("at the fixed step 3", rowA.notice) && occursin("at the mirror step -3", rowA.notice)
      # Compensating pair: a ratio changer (1.25 % per step) and the psi 90
      # changer on br5, truth (2, 3), start (0, 0) from a flat start, a Va
      # row at the far bus and an Ia row at the from end of br5 (one PMU
      # time base). At the flat start br5 carries no current, the Ia row does
      # not count, so the start test keeps the phase changer out
      # (:no_angle_anchor) and the ratio changer alone absorbs the magnitude.
      # At the converged state the current flows, the Ia row anchors the
      # angle, the phase changer adds rank and is estimated again together
      # with the ratio changer. That re-estimate is taken over only with
      # sigma_s below 0.5 steps. With the default sigmas it is weak (5 % per
      # phase step: sigma_s 0.85, 1 %: 4.16, both estimates on the truth 3),
      # so the row sigmas here are 0.05 times the defaults: sigma_s scales
      # with them (0.04 and 0.21), both changers on the truth, J at rounding
      # level, the :conditional note of the ratio changer is gone. With the
      # default sigmas at 1 % the phase changer is not taken over
      # (:weak_at_estimate, transferred step 0) and the ratio changer keeps
      # its first estimate, which absorbed the phase effect (step 2 here),
      # with the :conditional note. Only the phase changer released, start
      # 0, tight rows: the start test keeps it out and the run has no tap
      # state at all; the recheck runs on that converged run as well and
      # estimates it at the truth
      pair5(r, p; du) = (tn = deepcopy(m5.net); applyTapNameplate!(tn.branchVec[5]; tap_step = 0.0125, tap_min_step = -8.0, tap_max_step = 8.0, tap_current_step = Float64(r), phase_du_step = du, phase_min_step = -10.0, phase_max_step = 10.0, phase_current_step = Float64(p), psi_deg = 90.0); tn)
      pmu5(tm; f = 1.0) = (far = Int(tm.branchVec[5].toBus); [m for m in generateMeasurementsFromPF(tm; includeImag = true, includeVa = true, includeIa = true, stddev = Dict(k => f * v for (k, v) in measurementStdDevs())) if !(m.typ == S.VaMeas && m.busIdx != far) && !(m.typ == S.IaMeas && !(m.branchIdx == 5 && m.direction == :from))])
      for du in (0.05, 0.01)
        tm = pair5(2, 3; du = du)
        redirect_stdout(() -> runpf!(tm; config = S._copy_powerflow_with(m5.config.powerflow; tol = 1e-10, max_iter = 40)), devnull)
        nm = pair5(0, 0; du = du)
        setTapEstimation!(nm; trafo = 5, mode = :both)
        # with the tight rows the first round (ratio changer alone) rounds
        # far from its continuous value: one warning about the J rise
        resB = @test_logs (:warn, r"the simultaneous fixation raised J") match_mode = :any se5(nm, pmu5(tm; f = 0.05))
        rr = only(t for t in resB.tapEstimates if t.role === :ratio)
        rp = only(t for t in resB.tapEstimates if t.role === :phase)
        @test resB.converged && rr.status === :estimated && rr.step_fixed == 2.0 && rr.note === :none
        @test rp.status === :estimated && rp.fixed && rp.step_fixed == 3.0 && resB.objectiveJ < 1e-8 && rp.sigma_s < 0.5
        du == 0.01 || continue
        nm = pair5(0, 0; du = du)
        setTapEstimation!(nm; trafo = 5, mode = :both)
        resW = se5(nm, pmu5(tm))
        rr = only(t for t in resW.tapEstimates if t.role === :ratio)
        rp = only(t for t in resW.tapEstimates if t.role === :phase)
        @test resW.converged && rp.status === :not_estimated && rp.reason === :weak_at_estimate && rp.step_fixed == 0.0 && occursin("estimated again at step 3 with sigma_s 4.16 steps", rp.notice)
        @test rr.status === :estimated && rr.step_fixed == 2.0 && rr.note === :conditional
        # the phase row shows its weak re-estimate and sigma_s; ratio_est and
        # angle_est_deg of both rows stay the composite at the steps the
        # state was solved with (ratio estimate, phase at its transferred step)
        sv = zeros(2)
        sv[rr.regulator] = rr.s_est
        sv[rp.regulator] = rp.s_fixed
        Rw, Pw = S._tap_set_ratio_shift(S._effective_tap_changers(nm.branchVec[5]), sv)
        @test isapprox(rp.s_est + rp.neutral_step, 3.0; atol = 0.05) && isapprox(rp.sigma_s, 4.16; atol = 0.005)
        @test all(t -> isapprox(t.ratio_est, Rw; rtol = 1e-12) && isapprox(t.angle_est_deg, Pw; atol = 1e-10), (rr, rp))
      end
      tm = d90(3)
      redirect_stdout(() -> runpf!(tm; config = S._copy_powerflow_with(m5.config.powerflow; tol = 1e-10, max_iter = 40)), devnull)
      nm = d90(0)
      setTapEstimation!(nm; trafo = 5, mode = :pst)
      resC = se5(nm, pmu5(tm; f = 0.05))
      rowC = only(resC.tapEstimates)
      @test resC.converged && resC.tapFixation !== nothing && resC.tapFixation.fixed && rowC.status === :estimated && rowC.step_fixed == 3.0 && resC.objectiveJ < 1e-8 && rowC.sigma_s < 0.5
      # The angle anchor of a first-round fixed changer is judged at the
      # converged state: the psi 90 changer alone, truth 3, start 1, exact
      # rows with the Va row at the far bus and the Ia row at br5, default
      # sigmas. At the flat start the edge has no anchor: the start test
      # gives the changer the :sign_from_start note of an unanchored phase
      # changer, it is pinned by its magnitude effect alone and estimated
      # (sigma_s about 3.8, no bound in the first round; from start 2 the
      # edge is anchored already at the start). At the converged state the Ia row
      # anchors the edge, so no mirror check runs and the changer stays
      # fixed at 3. The mirror step -3 fits within the J difference of 9
      # (the run with the net at -3 gives J 2.2), so a mirror check with
      # the anchors of the start would have withdrawn it
      nm = d90(1)
      setTapEstimation!(nm; trafo = 5, mode = :pst)
      measK = pmu5(tm)
      resK = se5(nm, measK)
      rowK = only(resK.tapEstimates)
      @test resK.converged && rowK.status === :estimated && rowK.fixed && rowK.step_fixed == 3.0 && rowK.reason === :none && rowK.note === :sign_from_start && resK.objectiveJ < 1e-8
      @test se5(d90(-3), measK).objectiveJ < S._TAP_MIRROR_DJ
    end
    # A weak re-estimate is not taken over: the gate net with a 10 kV bus
    # behind a bridge transformer (ratio changer of 0.625 %, truth step 2,
    # start 0, 8 MW load), noisy rows (seed 11) without the voltage of the
    # cut-off bus. The start test keeps the changer out
    # (:radial_no_voltage_pin); at the converged state the losses of the
    # load give its column rank, and the joint estimate puts it at step
    # 6.33 with sigma_s 6.76 steps, above the bound of 0.5. The changer
    # keeps its transferred step 0 with :weak_at_estimate, the notice names
    # the estimate and sigma_s. Before the bound it was fixed at 6 (J 49.6)
    begin
      radial() = (tn = _tap_gate_net(); addBus!(net = tn, busName = "L3", vn_kV = 10.0); addProsumer!(net = tn, busName = "L3", type = "ENERGYCONSUMER", p = 8.0, q = 2.0); addPIModelTrafo!(net = tn, fromBus = "L2", toBus = "L3", r_pu = 0.004, x_pu = 0.08, b_pu = 0.0, ratio = 1.0, shift_deg = 0.0, status = 1); tn.branchVec[5].has_ratio_tap = true; tn.branchVec[5].ratio_tap_declared = true; tn.branchVec[5].tap_step = 0.00625; tn)
      tw = radial()
      tw.branchVec[5].tap_ratio = tw.branchVec[5].ratio / (1.0 + 2 * 0.00625)
      measW = [m for m in _tap_truth_meas(tw; noise = true) if !(m.typ == Sparlectra.VmMeas && m.busIdx == 5)]
      nw = radial()
      setTapEstimation!(nw; trafo = 5, mode = :ratio)
      resW = @test_logs (:info, r"^runse!: recheck at the estimate: tap regulator .* estimated again at step 6\.33\d with sigma_s 6\.76 steps \(bound 0\.5\), not taken over \(weak_at_estimate\), the transferred step 0 stays$") match_mode = :any with_state_estimation_config(max_iter = 40, tol = 1e-10) do
        runse!(nw, measW)
      end
      rowW = only(resW.tapEstimates)
      @test resW.converged && rowW.status === :not_estimated && rowW.reason === :weak_at_estimate && rowW.step_fixed == 0.0 && !rowW.fixed
      @test occursin("estimated again at step 6.33", rowW.notice) && occursin("sigma_s 6.76 steps", rowW.notice)
      # the row keeps the continuous estimate and its sigma_s; the final
      # state comes from a new solve with the transferred step, so it is the
      # state of the same net with the changer not released, within the
      # tolerance of the solves (1e-6, the floor at jac_eps; the two runs
      # stop after different iterations, max |dV| 1.1e-8)
      @test isapprox(rowW.s_est, 6.33; atol = 0.005) && isapprox(rowW.sigma_s, 6.76; atol = 0.005) && rowW.s_fixed == 0.0
      resN = with_state_estimation_config(() -> runse!(radial(), measW); max_iter = 40, tol = 1e-10)
      @test isapprox(resW.voltages, resN.voltages; atol = 1e-6) && isapprox(resW.objectiveJ, resN.objectiveJ; rtol = 1e-9)
    end

    # The bad-data rounds of the SE service (0.60.0) on a
    # gross error that drives the free step into the fallback: sp_case14
    # with a 5 % Delta-u changer at psi 60 on branch 16 (Moorau_110 to
    # Moorau_20, truth step -6, start 0) and the Vm row of Moorau_20 0.1 pu
    # low; the other 20 kV voltages keep the step observable without that
    # row. Round 1 does not converge released and falls back, its
    # diagnostics on the transferred step 0 eliminate the gross row, round
    # 2 estimates step -6, J 126.8 at dof 107.
    # The order before (diagnostics with the free step, then the fallback
    # for good) removed the gross row and three good 20 kV voltages and
    # ended at step 0, not estimated, J 144.8 at dof 104 (measured with
    # that order on the code before the change).
    begin
      S = Sparlectra
      case14du = joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")
      m14 = redirect_stdout(() -> S.import_case(case14du, S.resolve_config(S.DEFAULT_SPARLECTRA_CONFIG_PATH, case14du).config; run_kind = :state_estimation), devnull)
      empty!(m14.net.measurements)
      du16(n) = (tn = deepcopy(m14.net); applyTapNameplate!(tn.branchVec[16]; phase_du_step = 0.05, phase_min_step = -10.0, phase_max_step = 10.0, phase_current_step = Float64(n), psi_deg = 60.0); tn)
      t16 = du16(-6)
      _, e16 = redirect_stdout(() -> runpf!(t16; config = S._copy_powerflow_with(m14.config.powerflow; tol = 1e-10, max_iter = 40)), devnull)
      @test e16 == 0
      lv16 = Int(t16.branchVec[16].toBus)
      meas16 = [m.typ == S.VmMeas && m.busIdx == lv16 ? Measurement(typ = S.VmMeas, value = m.value - 0.1, sigma = m.sigma, busIdx = lv16, id = m.id * "_gross") : m for m in generateMeasurementsFromPF(t16; includeImag = true, noise = true, rng = Random.MersenneTwister(1))]
      nr = du16(0)
      setTapEstimation!(nr; trafo = 16, mode = :pst)
      append!(nr.measurements, meas16)
      cfgR = S._sparlectra_config_with(m14.config; state_estimation = S._se_config_with(m14.config.state_estimation; max_iter = 40, tol = 1e-6, update_net = true))
      logR = tempname()
      rounds16() = S.with_sparlectra_config(() -> SparlectraApp._se_bad_data_rounds!(nr, cfgR, cfgR.state_estimation.tap_estimation, [16], 1, true, logR, _ -> nothing), cfgR)
      bd = @test_logs (:warn, r"^tap estimation: run 1 not converged \(estimation_not_converged, 40 iterations\), falling back$") match_mode = :any rounds16()
      rowR = only(t for t in bd.res.tapEstimates if t.role === :phase)
      @test bd.rounds == 2 && !bd.tap_fallback_used && bd.res.converged
      @test [(e.round, e.id) for e in bd.eliminations] == [(1, only(m.id for m in meas16 if endswith(m.id, "_gross")))]
      @test (rowR.status, rowR.fixed, rowR.step_fixed) == (:estimated, true, -6.0)
      logRtext = read(logR, String)
      rm(logR; force = true)
      @test occursin(r"_gross eliminated on the end steps of the result run \(fallback without tap estimation, J = [0-9.eE+-]+, dof \d+\); the result run is repeated with released taps on the cleaned set", logRtext)
      @test occursin("tap estimation (round 2): staged start", logRtext)
    end

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
      # one iteration is not enough for anything: the staged start with the
      # taps frozen fails (the case file stores no voltages, so the start
      # ladder has no second try), the released run is not started (on
      # case6495rte it sat in the rank test for hours), the run without tap
      # estimation fails too, and the status names both
      res = Logging.with_logger(Logging.NullLogger()) do
        redirect_stdout(devnull) do
          SparlectraApp._run_state_estimation_service(case, cfg, joinpath(d, "run"), "tapfb", set; max_iter = 1)
        end
      end
      dd = SparlectraApp.to_dict(res)
      log = read(joinpath(d, "run", "run.log"), String)
      @test dd["status"] == "failed" && get(dd["metadata"], "se_tap_fallback_used", false) == true
      @test dd["metadata"]["se_tap_fallback_reason"] == "staged_not_converged" && dd["metadata"]["se_tap_staged_start_from"] == "flat"
      @test occursin("tap estimation: the staged start did not converge with the taps frozen (1 iteration(s)); the released run is not started", log)
      @test occursin("run 2 without tap estimation: not converged (not_converged, 1 iterations)", log) && !occursin("run 1 not converged", log)
      @test occursin("the staged start with the taps frozen (1 iterations), run 2 without tap estimation (not_converged, 1 iterations)", dd["message"])
      # the same with the fallback switched off: the first run only
      resNo = Logging.with_logger(Logging.NullLogger()) do
        redirect_stdout(devnull) do
          SparlectraApp._run_state_estimation_service(case, cfg, joinpath(d, "run_no"), "tapno", set; max_iter = 1, config_overrides = Dict{String,Any}("state_estimation.tap_estimation.fallback" => false))
        end
      end
      ddNo = SparlectraApp.to_dict(resNo)
      @test ddNo["status"] == "failed" && get(ddNo["metadata"], "se_tap_fallback_used", true) == false
      @test !occursin("falling back", read(joinpath(d, "run_no", "run.log"), String))

      # and the honest case: a normal cap converges, and then no fallback
      # line appears at all
      res_ok = redirect_stdout(devnull) do
        SparlectraApp._run_state_estimation_service(case, cfg, joinpath(d, "run_ok"), "tapok", set; max_iter = 50)
      end
      d_ok = SparlectraApp.to_dict(res_ok)
      @test d_ok["status"] == "succeeded"
      log_ok = read(joinpath(d, "run_ok", "run.log"), String)
      @test !occursin("falling back", log_ok)
      @test get(d_ok["metadata"], "se_tap_fallback_used", true) == false
      # the diagnostics artifact states the reported numbers first, the
      # same J and dof run.log carries (the diagnostics pass's own J can
      # differ when it eliminated rows)
      diag_ok = read(joinpath(d, "run_ok", "se_diagnostics.md"), String)
      @test occursin("## Reported result", diag_ok)
      jm = match(r"converged in \d+ iteration\(s\); J = ([0-9.eE+-]+), dof = (\d+)", log_ok)
      @test jm !== nothing
      @test occursin("**Objective J:** $(jm.captures[1]) (dof=$(jm.captures[2]))", diag_ok)
      @test occursin("**Converged:** true in", diag_ok)

      # switched off in the service (the Web UI box rides as this override):
      # every changer of se_tap_estimates.csv not estimated, reason disabled,
      # the machine transformer included, nothing calculated
      res_off = redirect_stdout(devnull) do
        SparlectraApp._run_state_estimation_service(case, cfg, joinpath(d, "run_off"), "tapoff", set; max_iter = 50, config_overrides = Dict{String,Any}("state_estimation.tap_estimation.enabled" => false))
      end
      d_off = SparlectraApp.to_dict(res_off)
      @test d_off["status"] == "succeeded" && d_off["metadata"]["se_tap_estimation"] == false && d_off["metadata"]["se_tap_fixed"] == false
      csv_off = readlines(joinpath(d, "run_off", "se_tap_estimates.csv"))
      hdr = split(csv_off[1], ",")
      st, rs = findfirst(==("status"), hdr), findfirst(==("reason"), hdr)
      @test length(csv_off) == 3 && all(l -> (f = split(l, ","); f[st] == "not_estimated" && f[rs] == "disabled"), csv_off[2:end])
      @test occursin("tap estimation: switched off (state_estimation.tap_estimation.enabled = false); 2 transformer(s) keep their transferred step", read(joinpath(d, "run_off", "run.log"), String))
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
    # surfaces, because its J measures the transferred steps, not estimates.
    # Under a flat start the staged start (taps
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
              joinpath(d, tag), "tapst188", set; max_iter = 4, flatstart = true)
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
        @test get(dst["metadata"], "se_tap_fallback_used", false) == false
        @test get(dst["metadata"], "se_tap_j_rise_exceeds_quantile", false) == rise && get(dst["metadata"], "se_tap_sequential_fixation", false) == rise
      end
      res = redirect_stdout(devnull) do
        SparlectraApp._run_state_estimation_service(shipped, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH,
          joinpath(d, "fb"), "tapfb188", set; max_iter = 2, flatstart = false)
      end
      dd = SparlectraApp.to_dict(res)
      @test dd["status"] == "succeeded"
      @test get(dd["metadata"], "se_tap_fallback_used", false) == true
      # run.log: the two lines, the status line names the fallback, every
      # changer row says why it was not estimated
      log188 = read(joinpath(d, "fb", "run.log"), String)
      @test occursin(r"tap estimation: run 1 not converged \(estimation_not_converged, 2 iterations\), falling back", log188)
      @test occursin("run 2 without tap estimation: converged, J = $(round(dd["metadata"]["se_objective"]; sigdigits = 6)), dof = $(dd["metadata"]["se_dof"])", log188)
      @test occursin("Tap estimation fallback: run 1 not converged", dd["message"])
      @test all(t -> t["status"] == "calculated" || t["reason"] == "fallback_not_converged", dd["metadata"]["se_tap_estimates"])

      # surface 3: se_diagnostics.md
      diagmd = read(joinpath(d, "fb", "se_diagnostics.md"), String)
      @test occursin(SparlectraApp._SE_TAP_FALLBACK_NOTE, diagmd)

      # surface 1: the result summary of the Web UI
      summary = SparlectraApp._webui_se_summary(dd)
      @test summary !== nothing
      @test occursin(SparlectraApp._SE_TAP_FALLBACK_NOTE, summary)

      # surface 2: the tap section, one line above the table whose rows
      # say why they were not estimated
      table = SparlectraApp._webui_se_tap_section(dd)
      @test occursin("did NOT converge", table)
      @test occursin("not estimated (fallback not converged)", table)

      # and the run history must not call a state estimation "rectangular"
      # (as seen): the method comes from the run kind
      @test SparlectraApp._powerflow_run_index_solver(res) == "wls"

      # The start ladder through the service, cap 2 from the flat start:
      # the staged start converges from the case voltages, the released run
      # fails at cap 2 and falls back, and that run converges on its second
      # try from the case voltages (without the ladder the staged start
      # fails here, the warmup_casePST block below)
      res_sl = redirect_stdout(devnull) do
        SparlectraApp._run_state_estimation_service(shipped, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, joinpath(d, "sl"), "tapsl188", set; max_iter = 2, flatstart = true)
      end
      dsl = SparlectraApp.to_dict(res_sl)
      logsl = read(joinpath(d, "sl", "run.log"), String)
      @test dsl["status"] == "succeeded"
      @test dsl["metadata"]["se_tap_staged_start_from"] == "stored_after_flat" && dsl["metadata"]["se_start"] == "stored_after_flat"
      @test occursin("tap estimation: staged start: the WLS from the flat start did not converge; second try from the voltages stored in the case", logsl)
      @test occursin(r"staged start, one solve with the \d+ released tap\(s\) frozen \(\d+ iteration\(s\) from the case voltages", logsl)
      @test occursin("start: the WLS from the flat start did not converge; second try from the voltages stored in the case (state_estimation.start_ladder); the result run starts from the case voltages", logsl)
    end

    # The acceptance rule of the staged start (worse_fit), the second
    # trigger of the one fallback: the fixed steps are kept only if they fit
    # at least as well as the transferred steps in the pre-solve. Trigger:
    # sp_case14 with the ratio changer of branch 31 at the model position
    # 4.4, between two physical steps, and exact measurements of that state.
    # The estimate comes back to 4.4, the fixation can only reach step 4:
    # J 0.233 against 4.4e-16 at the model positions (dof 83 both), so the
    # estimate is rejected and the run without tap estimation is the result,
    # with its own reason and note on every surface. The rise stays below
    # the 95 % quantile 5.99 of chi2(2), so no warning. A noisy set of the
    # same state was kept (seed 42: J 77.94 fixed against 78.09 at the model
    # positions). With fallback = false the fixed estimate is the result.
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
          joinpath(d, "rej"), "offgrid_rej", set; max_iter = 50)
      end
      dd = SparlectraApp.to_dict(res)
      log = read(joinpath(d, "rej", "run.log"), String)
      diag = read(joinpath(d, "rej", "se_diagnostics.md"), String)
      _assert_tap_fallback_report(dd, log, diag, "offgrid_rej", set; case = basename(offgrid), reason = "worse_fit")
      @test get(dd["metadata"], "se_tap_staged_start", false) == true
      @test get(dd["metadata"], "se_tap_fallback_reason", "") == "worse_fit"
      @test occursin("Tap estimation fallback: run 1 rejected (worse_fit", dd["message"])
      @test occursin(SparlectraApp._SE_TAP_REJECTED_NOTE, diag)
      @test occursin(SparlectraApp._SE_TAP_REJECTED_NOTE, SparlectraApp._webui_se_summary(dd))
      @test occursin(SparlectraApp._webui_escape(SparlectraApp._SE_TAP_REJECTED_NOTE), SparlectraApp._webui_se_tap_section(dd))
      @test occursin("not estimated (fallback worse fit)", SparlectraApp._webui_se_tap_section(dd))
      # the result is the run without tap estimation, bit for bit: J,
      # dof, iterations and the state file equal those of the same run with
      # tap estimation switched off (the repeat starts from the start state
      # and configuration of the round, not from the staged pre-solve).
      resOff = redirect_stdout(devnull) do
        SparlectraApp._run_state_estimation_service(offgrid, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH,
          joinpath(d, "off"), "offgrid_off", set; max_iter = 50, config_overrides = Dict{String,Any}("state_estimation.tap_estimation.enabled" => false))
      end
      ddOff = SparlectraApp.to_dict(resOff)
      @test ddOff["status"] == "succeeded" && get(ddOff["metadata"], "se_tap_fallback_used", false) == false
      @test dd["metadata"]["se_objective"] === ddOff["metadata"]["se_objective"] && dd["metadata"]["se_dof"] == ddOff["metadata"]["se_dof"] && dd["metadata"]["se_iterations"] == ddOff["metadata"]["se_iterations"]
      @test read(joinpath(d, "rej", "se_state.csv"), String) == read(joinpath(d, "off", "se_state.csv"), String)
      # fallback off: the first run unvarnished, the fixed estimate
      resNo = redirect_stdout(devnull) do
        SparlectraApp._run_state_estimation_service(offgrid, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH,
          joinpath(d, "kept"), "offgrid_kept", set; max_iter = 50, config_overrides = Dict{String,Any}("state_estimation.tap_estimation.fallback" => false))
      end
      ddNo = SparlectraApp.to_dict(resNo)
      @test ddNo["status"] == "succeeded" && get(ddNo["metadata"], "se_tap_fallback_used", true) == false && ddNo["metadata"]["se_tap_fixed"] == true
      @test ddNo["metadata"]["se_objective"] > dd["metadata"]["se_objective"]
    end

    # Bad data on the end steps (0.60.0): with released
    # taps the diagnostics run on the fixed end steps of the result run,
    # and after an elimination the result run is repeated with released
    # taps on the cleaned set. Trigger: sp_case14 with the ratio changer of
    # branch 16 (Moorau_110-Moorau_20) two steps off its model step in the
    # truth, noise of seed 42, and a gross error of -2 kV on the 20 kV
    # voltage of that transformer, the voltage a free ratio step moves. Round
    # 1 estimates and fixes the taps with the gross row in, its diagnostics
    # on those steps eliminate that row; round 2 estimates branch 16 two
    # steps off again. Truth of the result: the service run on the same set
    # with that row marked invalid (same J, dof and steps). A full
    # sequential pass on the round-1 steps took two good 20 kV voltages out
    # as well, one row per round only the gross one.
    mktempdir() do d
      case14 = joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json")
      clean = joinpath(d, "gen14_tap.csv")
      g = SparlectraApp._se_generate_measurement_set(case14, clean, SparlectraApp.MeasurementGeneratorOptions(; noise = true, gross_k = 0.0, tap_steps = 2.0, tap_count = 1,
        include_i = false, sigma_u_pct = 0.5, sigma_i_pct = 1.0, sigma_p_pct = 1.0, sigma_q_pct = 1.0, sigma_ia_deg = 0.0, seed = 42))
      @test occursin("tap deviation 2.0 steps", g.tap_note) && occursin("on transformer branch 16", g.tap_note)
      lines = readlines(clean)
      k = findfirst(l -> startswith(l, ",,,Moorau_20,U,"), lines)
      @test k !== nothing
      edit(f) = (c = split(lines[k], ","); f(c); vcat(lines[1:(k-1)], join(c, ","), lines[(k+1):end]))
      gross = joinpath(d, "gen14_gross.csv")
      write(gross, join(edit(c -> c[6] = string(parse(Float64, c[6]) - 2.0)), "\n") * "\n")
      without = joinpath(d, "gen14_without.csv")
      write(without, join(edit(c -> c[end] = "false"), "\n") * "\n")
      run14(set, tag) = SparlectraApp.to_dict(redirect_stdout(devnull) do
        SparlectraApp._run_state_estimation_service(case14, Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, joinpath(d, tag), tag, set; max_iter = 50)
      end)
      dg = run14(gross, "end_steps")
      dw = run14(without, "without_row")
      @test dg["status"] == "succeeded" && dw["status"] == "succeeded"
      mg, mw = dg["metadata"], dw["metadata"]
      @test mg["se_bad_data_rounds"] == 2 && mw["se_bad_data_rounds"] == 1
      @test mg["se_eliminations"] == 1 && mw["se_eliminations"] == 0
      @test mg["se_tap_fallback_used"] == false
      log14 = read(joinpath(d, "end_steps", "run.log"), String)
      @test occursin(r"bad data: U@Moorau_20 eliminated on the end steps of the result run \(taps estimated and fixed, J = [0-9.eE+-]+, dof 83\); the result run is repeated with released taps on the cleaned set", log14)
      @test occursin("tap estimation (round 2): staged start", log14)
      reg16(m) = only(t for t in m["se_tap_estimates"] if t["branch"] == 16 && t["role"] == "ratio")
      @test (reg16(mg)["status"], reg16(mg)["fixed"], reg16(mg)["deviation_steps"]) == ("estimated", true, 2)
      @test isapprox(mg["se_objective"], mw["se_objective"]; rtol = 1e-8) && mg["se_dof"] == mw["se_dof"]
      @test [(t["branch"], t["regulator"], t["step_estimated"]) for t in mg["se_tap_estimates"]] == [(t["branch"], t["regulator"], t["step_estimated"]) for t in mw["se_tap_estimates"]]
    end

    # the helpers the fallback is built from
    net = SparlectraApp._se_import_case_net(joinpath(dirname(@__DIR__), "data", "scf", "sp_case14.scf.json"),
      Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true))
    @test SparlectraApp._freeze_all_tap_estimation!(net) == 0
  end
  return true
end

  @testset "State estimation" begin (function ()
    tests = [
      ("WLS", test_state_estimation_wls_first_version),
      ("Start ladder", test_state_estimation_start_ladder),
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
      ("Tap step from the control target", test_state_estimation_tap_step_from_target),
      ("Tap fallback: diagnostics artifact matches the reported run (MiniGrid)", test_state_estimation_tap_fallback_diagnostics_consistent),
      ("Tap fallback: eliminations of the repeat (MiniGrid)", test_state_estimation_elimination_reverted),
      ("Measurement generator critical thinning", test_measurement_generator_critical_thinning),
      ("Measurement generator tap deviation on the changer set", test_measurement_generator_tap_changer_sets),
      ("Tap estimation ladder (stages 1, 2, 4)", test_tap_estimation_ladder),
      ("Tap estimation ladder, noisy measurements (stage 3)", test_tap_estimation_ladder_noise),
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
