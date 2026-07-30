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
#
# file: test/test_island_diagnostics.jl
# purpose: regression tests for AC island diagnostics reporting -- the failure
#          message must show the failing island's own iteration count and
#          stage, and islands that were never attempted must not inherit the
#          failed island's solve statistics in ac_island_solver_summary.csv.
#          Repro source: a 6209-bus CGMES delivery (1 main island + 158
#          single-bus islands) where the run-level message claimed
#          "iterations=0 / stage=before_nr" while ac_island_1_solver.log
#          recorded 80 Newton iterations, and every single-bus island row
#          duplicated island 1's mismatch/iteration/switching values.

using Test
using Sparlectra

# Per-island status as the rectangular solver stores it in
# performance_profile[:ac_island_solver_statuses] for a failed island. Values
# mirror the observed CGMES run (80 NR iterations, active-set instability).
function _island_diag_failed_status()
  return (
    island_id = 1,
    status = :not_converged,
    final_converged = false,
    numerical_converged = false,
    iterations = 80,
    initial_mismatch = 2220.614510584437,
    final_mismatch = 240324.3267233851,
    best_mismatch = 731.175415721738,
    reason = :nr_mismatch_not_converged_active_set_unstable,
    stage = :newton_iteration,
    pv_pq_switching_events = 435,
    qlimit_active_set_changes = 13,
    mismatch_history = [2220.614510584437, 731.175415721738, 240324.3267233851],
    exception_type = "",
    exception_message = "",
    stacktrace_top = "",
  )
end

# Two-bus island with one branch plus one isolated bus. The isolated bus forms
# a second diagnostics island that the solver never attempts -- the miniature
# version of the CGMES delivery's 158 single-bus islands.
function _island_diag_testnet()
  net = Net(name = "island_diag", baseMVA = 100.0)
  for b in ("A1", "A2", "B1")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.02, status = 1)
  addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "A1")
  addProsumer!(net = net, busName = "A2", type = "ENERGYCONSUMER", p = 50.0, q = 10.0)
  return net
end

# Split a summary CSV into a header-name => column-index map plus data rows.
function _island_diag_read_summary(path::AbstractString)
  lines = split(strip(read(path, String)), '\n')
  header = split(lines[1], ',')
  col = Dict(String(name) => i for (i, name) in enumerate(header))
  rows = [split(line, ',') for line in lines[2:end]]
  return col, rows
end

function run_island_diagnostics_tests()
  @testset "AC island diagnostics reporting" begin
    failed_status = _island_diag_failed_status()

    @testset "failure message reports the failing island's own status (bug A)" begin
      # Combined run status as _rectangular_run_status builds it: it carries
      # no :iterations/:stage/:island_id -- the per-island record must win.
      run_status = (
        outcome = :not_converged,
        numerical_converged = false,
        solution_available = false,
        limit_validation_status = :skip,
        final_converged = false,
        reason = :nr_mismatch_not_converged_active_set_unstable,
        reason_text = "NR mismatch did not converge",
        final_mismatch = 240324.3267233851,
      )
      settings = (start_projection = true,)
      profile = Dict{Symbol,Any}(
        :ac_island_artifacts => ("ac_island_solver_summary.csv", "ac_island_1_solver.log"),
        :ac_island_diagnostics => [
          (island_id = 1, n_bus = 6051, n_branch = 8202, ref_bus = 696, n_pq = 6050, n_pv = 0, n_ref = 1, settings = settings),
          (island_id = 2, n_bus = 1, n_branch = 0, ref_bus = 37, n_pq = 1, n_pv = 0, n_ref = 0, settings = settings),
        ],
        :ac_island_solver_statuses => Dict(1 => failed_status),
      )
      out = Sparlectra._append_island_failure_message(run_status, profile)
      @test occursin("AC island 1 power-flow solve failed", out.reason_text)
      @test occursin("iterations=80", out.reason_text)
      @test occursin("stage=newton_iteration", out.reason_text)
      @test !occursin("iterations=0", out.reason_text)
      @test !occursin("stage=before_nr", out.reason_text)
      @test occursin("final_mismatch=240324.3267233851", out.reason_text)
      @test occursin("reason=nr_mismatch_not_converged_active_set_unstable", out.reason_text)

      # When the combined status is untagged, the failing island is looked up
      # in the per-island statuses -- not assumed to be the first row.
      profile_island2 = copy(profile)
      profile_island2[:ac_island_solver_statuses] = Dict(
        1 => (island_id = 1, status = :converged, iterations = 6, final_mismatch = 1.0e-9, reason = :none, stage = :post_solve_validation),
        2 => merge(failed_status, (island_id = 2, iterations = 40)),
      )
      out2 = Sparlectra._append_island_failure_message(run_status, profile_island2)
      @test occursin("AC island 2 power-flow solve failed", out2.reason_text)
      @test occursin("iterations=40", out2.reason_text)

      # Single-island runs store no per-island statuses; the combined status
      # (with the run iteration count merged in by _build_sparlectra_result)
      # must still never claim "before_nr" for a solve that ran N>0 iterations.
      profile_single = Dict{Symbol,Any}(
        :ac_island_artifacts => ("ac_island_solver_summary.csv",),
        :ac_island_diagnostics => [(island_id = 1, n_bus = 12, n_branch = 13, ref_bus = 10, n_pq = 10, n_pv = 1, n_ref = 1, settings = settings)],
      )
      out3 = Sparlectra._append_island_failure_message(merge(run_status, (iterations = 1,)), profile_single)
      @test occursin("iterations=1", out3.reason_text)
      @test occursin("stage=during_nr", out3.reason_text)
      @test !occursin("stage=before_nr", out3.reason_text)
    end

    @testset "unsolved islands do not inherit the failed island's statistics (bug B)" begin
      mktempdir() do dir
        net = _island_diag_testnet()
        cfg = Sparlectra.PowerFlowConfig()
        profile = Dict{Symbol,Any}(
          :output_dir => dir,
          :ac_island_solver_statuses => Dict(1 => failed_status),
        )
        # Mirrors execution.jl on failure: the combined status handed to the
        # diagnostics writer is the failed island's own status (tagged with
        # island_id = 1), plus the run-level iteration count.
        Sparlectra._write_ac_island_diagnostics!(net, cfg, profile; status = failed_status, iterations = 80)
        col, rows = _island_diag_read_summary(joinpath(dir, "ac_island_solver_summary.csv"))
        @test length(rows) == 2
        row1 = rows[1]
        row2 = rows[2]

        # Solved island keeps its true record (regression guard).
        @test row1[col["island_id"]] == "1"
        @test row1[col["iterations"]] == "80"
        @test row1[col["final_mismatch"]] == "240324.3267233851"
        @test row1[col["final_status"]] == "not_converged"
        @test row1[col["failure_reason"]] == "nr_mismatch_not_converged_active_set_unstable"
        @test row1[col["stage"]] == "newton_iteration"
        @test row1[col["pv_pq_switching_events"]] == "435"

        # The never-attempted single-bus island must not repeat island 1's
        # mismatch/iteration/switching values.
        @test row2[col["island_id"]] == "2"
        @test row2[col["iterations"]] == "0"
        @test row2[col["final_status"]] == "not_attempted"
        @test row2[col["failure_reason"]] == "not_attempted"
        @test row2[col["stage"]] == "not_attempted"
        @test row2[col["final_mismatch"]] == "unavailable"
        @test row2[col["initial_mismatch"]] == "unavailable"
        @test row2[col["mismatch_status"]] == "unavailable"
        @test row2[col["pv_pq_switching_events"]] == "0"
        @test !occursin("240324", join(row2, ','))

        # Attempted islands keep their per-island artifacts; never-attempted
        # islands are represented in the summary CSV only (open design
        # question 2: summary-only, no header-only file pairs).
        @test isfile(joinpath(dir, "ac_island_1_solver.log"))
        @test isfile(joinpath(dir, "ac_island_1_mismatch_history.csv"))
        @test !isfile(joinpath(dir, "ac_island_2_solver.log"))
        @test !isfile(joinpath(dir, "ac_island_2_mismatch_history.csv"))
        history = read(joinpath(dir, "ac_island_1_mismatch_history.csv"), String)
        @test occursin("1,2220.614510584437", history)

        # Bug C: q_limit_processing_status must not alias failure_reason.
        @test row1[col["qlimits_enabled"]] == "true"
        @test row1[col["q_limit_processing_status"]] == "unavailable"
        @test row2[col["q_limit_processing_status"]] == "not_attempted"
      end
    end

    @testset "single-island run keeps the combined-status fallback" begin
      # A connected net solved outside the island-wise path stores no
      # per-island statuses; its one diagnostics row must still be filled
      # from the combined status.
      mktempdir() do dir
        net = Net(name = "island_diag_single", baseMVA = 100.0)
        for b in ("A1", "A2")
          addBus!(net = net, busName = b, vn_kV = 110.0)
        end
        addPIModelACLine!(net = net, fromBus = "A1", toBus = "A2", r_pu = 0.01, x_pu = 0.08, b_pu = 0.02, status = 1)
        addProsumer!(net = net, busName = "A1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "A1")
        addProsumer!(net = net, busName = "A2", type = "ENERGYCONSUMER", p = 50.0, q = 10.0)
        cfg = Sparlectra.PowerFlowConfig()
        profile = Dict{Symbol,Any}(:output_dir => dir)
        combined = (status = :converged, final_converged = true, iterations = 5, final_mismatch = 1.0e-10, reason = :none, stage = :post_solve_validation)
        Sparlectra._write_ac_island_diagnostics!(net, cfg, profile; status = combined, iterations = 5)
        col, rows = _island_diag_read_summary(joinpath(dir, "ac_island_solver_summary.csv"))
        @test length(rows) == 1
        @test rows[1][col["iterations"]] == "5"
        @test rows[1][col["final_status"]] == "converged"
        @test rows[1][col["stage"]] == "post_solve_validation"
        @test isfile(joinpath(dir, "ac_island_1_solver.log"))
      end
    end
  end
end
