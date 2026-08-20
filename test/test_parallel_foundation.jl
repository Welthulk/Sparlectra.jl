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

# file: test/test_parallel_foundation.jl
# purpose: thread-safety foundation of the multi-core work (Phase 1):
#          runtime.parallel.* configuration and validation, the per-worker
#          performance-profile child/merge helpers, the DC status Net field,
#          UMFPACK copy(F) finalizer safety, and the startup summary line.

using Sparlectra
using Test
using SparseArrays
using LinearAlgebra

function run_parallel_foundation_tests()
  @testset "Parallel foundation (Phase 1)" begin
    @testset "runtime.parallel configuration" begin
      cfg = Sparlectra.RuntimeConfig(Dict{String,Any}())
      @test cfg.parallel.enabled === true
      @test cfg.parallel.max_tasks == "auto"
      @test cfg.parallel.min_work_items == 4
      @test parallel_max_tasks(cfg.parallel) == Threads.nthreads()

      raw = Dict{String,Any}("runtime" => Dict{String,Any}("parallel" => Dict{String,Any}("enabled" => false, "max_tasks" => "2", "min_work_items" => 8)))
      cfg2 = Sparlectra.RuntimeConfig(raw)
      @test cfg2.parallel.enabled === false
      @test parallel_max_tasks(cfg2.parallel) == 2
      @test cfg2.parallel.min_work_items == 8

      bad_tasks = Dict{String,Any}("runtime" => Dict{String,Any}("parallel" => Dict{String,Any}("max_tasks" => "0")))
      @test_throws ArgumentError Sparlectra.RuntimeConfig(bad_tasks)
      bad_tasks2 = Dict{String,Any}("runtime" => Dict{String,Any}("parallel" => Dict{String,Any}("max_tasks" => "many")))
      @test_throws ArgumentError Sparlectra.RuntimeConfig(bad_tasks2)
      bad_items = Dict{String,Any}("runtime" => Dict{String,Any}("parallel" => Dict{String,Any}("min_work_items" => 0)))
      @test_throws ArgumentError Sparlectra.RuntimeConfig(bad_items)

      # the yaml example carries the keys, so user files may set them
      loaded = load_sparlectra_config(; reload = true, overrides = Dict{String,Any}("runtime" => Dict{String,Any}("parallel" => Dict{String,Any}("max_tasks" => "3"))))
      @test parallel_max_tasks(loaded.runtime.parallel) == 3
    end

    @testset "performance-profile child and merge" begin
      recorder_calls = Ref(0)
      parent = Dict{Symbol,Any}(:enabled => true, :show_allocations => false, :output_dir => "/tmp/x", :phase_callback => _ -> (recorder_calls[] += 1), :cancellation_check => () -> nothing)

      @test Sparlectra._perf_profile_child(nothing) === nothing
      @test Sparlectra._perf_profile_child(Dict{Symbol,Any}(:enabled => false)) === nothing

      child = Sparlectra._perf_profile_child(parent)
      @test child isa Dict{Symbol,Any}
      @test child[:enabled] === true
      @test child[:output_dir] == "/tmp/x"
      @test haskey(child, :cancellation_check)
      # orchestrator-only: the recorder callback must NOT reach workers
      @test !haskey(child, :phase_callback)
      @test !haskey(child, :diagnostic_artifact_prefix)

      # worker writes: timings, iteration rows, per-worker scalars
      Sparlectra._perf_profile_add!(child, :ybus_assembly, 0.25, 100)
      Sparlectra._perf_profile_add!(child, :ybus_assembly, 0.25, 100)
      Sparlectra._perf_profile_push_iteration!(child, (iter = 1, mismatch = 0.1))
      child[:linear_solver_backend] = :umfpack_reuse

      # parent already carries one timing of the same phase (serial part)
      Sparlectra._perf_profile_add!(parent, :ybus_assembly, 1.0, 50)
      Sparlectra._perf_profile_merge!(parent, child, "ac_island_2_")

      row = parent[:timings][:ybus_assembly]
      @test row.calls == 3
      @test isapprox(row.elapsed_s, 1.5; atol = 1e-12)
      @test row.bytes == 250
      @test length(parent[:iterations]) == 1
      @test parent[:ac_island_2_linear_solver_backend] === :umfpack_reuse
      @test !haskey(parent, :linear_solver_backend)
      @test recorder_calls[] == 0

      # the orchestrating fan-out accounts its wall clock separately
      Sparlectra._perf_profile_add!(parent, :parallel_wall_time, 0.75, 0)
      @test parent[:timings][:parallel_wall_time].elapsed_s == 0.75

      # disabled parent: merge is a no-op
      off = Dict{Symbol,Any}(:enabled => false)
      @test Sparlectra._perf_profile_merge!(off, child, "x_") === off
      @test !haskey(off, :timings)
    end

    @testset "DC status lives on the Net" begin
      @test !isdefined(Sparlectra, :_DC_PF_STATUS)
      net = createTest3BusNet()
      @test Sparlectra.dc_pf_status(net) === nothing
      status = (converged = true,)
      @test Sparlectra._set_dc_pf_status!(net, status) === status
      @test Sparlectra.dc_pf_status(net) === status
      # AC and DC fields stay separate
      @test Sparlectra.rectangular_pf_status(net) === nothing
    end

    @testset "UMFPACK copy(F) is finalizer-safe" begin
      # Phase 0 review item 4: chunk workers hold copies of one factorization;
      # copies must not double-free the shared numeric object when collected.
      A = sprand(200, 200, 0.05) + 10.0 * I
      F = lu(A)
      b = collect(1.0:200.0)
      x_ref = F \ b
      for _ in 1:50
        Fc = copy(F)
        @test Fc \ b == x_ref
      end
      GC.gc()
      GC.gc()
      # original factorization stays usable after all copies are collected
      @test F \ b == x_ref
      # and a surviving copy stays usable after the ORIGINAL is collected
      keeper = copy(F)
      F = nothing
      GC.gc()
      GC.gc()
      @test keeper \ b == x_ref
    end

    @testset "startup summary line" begin
      cfg = Sparlectra.RuntimeConfig(Dict{String,Any}())
      status = Sparlectra.runtime_thread_status(cfg)
      @test status.parallel_enabled === true
      @test status.parallel_max_tasks == Threads.nthreads()
      out = sprint(Sparlectra.print_runtime_thread_config, status)
      @test occursin("parallel: enabled=true max_tasks=$(Threads.nthreads())", out)
    end

    @testset "Web UI option spec for runtime.parallel.enabled" begin
      spec = Sparlectra._webui_option_spec("runtime_parallel_enabled")
      @test spec.config_key == "runtime.parallel.enabled"
      @test spec.default === true
      @test spec.control === :checkbox
      @test spec.section === :expert
      form_html = Sparlectra.render_powerflow_form()
      @test occursin("name=\"runtime_parallel_enabled\"", form_html)
      overrides = Sparlectra.validate_gui_config_overrides(Dict{String,Any}("runtime.parallel.enabled" => "false"))
      @test overrides["runtime"]["parallel"]["enabled"] == "false"
    end
  end
end
