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

# file: test/test_legacy_observability.jl
# purpose: legacy profile: the testsets of test/test_observability.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_observability.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_observability_tests()
  # cut from the runner function
  # Aggregates the observability test sets; the shared nets are built once
  # and handed to every set that needs them.
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
    # a deficit on the large case: eleven buses without any measurement.
    # The raw pivot count came out two below the QR rank here, so the pivot
    # rule hands a deficit to the decomposition and says so
    dark = Set(50:60)
    filter!(m -> !(hasproperty(m, :busIdx) && m.busIdx in dark), netq.measurements)
    def_a = with_state_estimation_config(() -> evaluate_global_observability(netq); rank_method = :decomposition)
    def_b = with_state_estimation_config(() -> evaluate_global_observability(netq); rank_method = :pivots)
    @test !def_a.numerical_observable
    @test def_b.numerical_rank == def_a.numerical_rank
    @test def_b.rank_method === :decomposition
  end)() end

  return true
end
