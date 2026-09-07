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

# file: test/test_net_cache.jl
# purpose: tests the opt-in binary MATPOWER net cache: miss/store, identical
#          hit results, fallback on corrupt entries, and key invalidation
#
# Tests the opt-in binary net cache (model.net_cache_enabled,
# issue #292): INERT since task_import_direct (2026-09-04): the direct
# import no longer produces the converted case the cache stored. The set
# guards the loud warning and the inert behavior; the old
# miss/hit/corrupt/key coverage went with the feature, stated here as the
# deliberate removal of tests for behavior that is gone by design.

include("test_api_support.jl")

function _net_cache_run(case_path, output_dir; extra = Dict{String,Any}())
  overrides = Dict{String,Any}("model.auto_profile" => "off", "output.logfile_results" => "compact", "benchmark.enabled" => false)
  merge!(overrides, extra)
  return run_sparlectra_api(casefile = case_path, config_file = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, output_dir = output_dir, config_overrides = overrides)
end

function run_net_cache_tests()
  @testset "MATPOWER net cache (inert since task_import_direct)" begin
    # maintainer decision 2026-09-04: the cache stored the CONVERTED
    # SCFCase, the side product the direct import no longer creates, so
    # the feature fell away with the conversion (issue #292 follow-up).
    # This set now guards the two things that remain: the key WARNS out
    # loud instead of silently doing nothing, and an enabled cache
    # neither changes results nor writes a cache directory.
    mktempdir() do dir
      case_path = _write_api_test_case(joinpath(dir, "cache_case.m"))
      cache_dir = joinpath(dir, ".sparlectra_net_cache")
      enabled = Dict{String,Any}("model.net_cache_enabled" => true)

      reference = _net_cache_run(case_path, joinpath(dir, "run_ref"))
      @test reference.status === :succeeded
      @test !isdir(cache_dir)

      # the warning fires at the import context; the service run wraps its
      # own logging, so the log assertion talks to the context directly
      cfg_enabled = Sparlectra.SparlectraConfig(Dict{String,Any}("model" => Dict{String,Any}("net_cache_enabled" => true, "auto_profile" => "off")))
      @test_logs (:warn, r"net_cache_enabled is inert since the direct import") match_mode = :any begin
        Sparlectra._import_sparlectra_net(case_path, nothing, cfg_enabled)
      end
      inert = _net_cache_run(case_path, joinpath(dir, "run_inert"); extra = enabled)
      @test inert.status === :succeeded
      @test inert.converged == reference.converged
      @test inert.iterations == reference.iterations
      @test isapprox(inert.final_mismatch, reference.final_mismatch; atol = 1e-12)
      @test !isdir(cache_dir)
    end
  end
  return nothing
end
