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
#
# Tests the opt-in binary net cache (matpower_import.net_cache.enabled,
# issue #292): miss/store, hit with identical results, silent fallback on
# corrupt entries, key invalidation on case or option changes, and
# inactivity under auto_profile != off. Cases live in a temp directory so
# .sparlectra_net_cache never lands in the repository.

include("test_api_support.jl")

function _net_cache_run(case_path, output_dir; extra = Dict{String,Any}())
  overrides = Dict{String,Any}("matpower_import.auto_profile" => "off", "output.logfile_results" => "compact", "benchmark.enabled" => false)
  merge!(overrides, extra)
  return run_sparlectra_api(casefile = case_path, config_file = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, output_dir = output_dir, config_overrides = overrides)
end

function run_net_cache_tests()
  @testset "MATPOWER net cache" begin
    mktempdir() do dir
      case_path = _write_api_test_case(joinpath(dir, "cache_case.m"))
      cache_dir = joinpath(dir, ".sparlectra_net_cache")
      enabled = Dict{String,Any}("matpower_import.net_cache.enabled" => true)

      reference = _net_cache_run(case_path, joinpath(dir, "run_ref"))
      @test reference.status === :succeeded
      @test !isdir(cache_dir)

      miss = _net_cache_run(case_path, joinpath(dir, "run_miss"); extra = enabled)
      @test miss.status === :succeeded
      @test isdir(cache_dir)
      entries = readdir(cache_dir)
      @test length(entries) == 1

      hit = _net_cache_run(case_path, joinpath(dir, "run_hit"); extra = enabled)
      @test hit.status === :succeeded
      @test hit.converged == reference.converged
      @test hit.iterations == reference.iterations
      @test isapprox(hit.final_mismatch, reference.final_mismatch; atol = 1e-12)
      @test length(readdir(cache_dir)) == 1

      # Corrupt entry: silent fallback to a fresh parse, entry rewritten.
      cache_file = joinpath(cache_dir, only(readdir(cache_dir)))
      write(cache_file, "not a serialized cache entry")
      recovered = _net_cache_run(case_path, joinpath(dir, "run_corrupt"); extra = enabled)
      @test recovered.status === :succeeded
      @test recovered.iterations == reference.iterations
      @test filesize(cache_file) > 100

      # Import options are NOT part of the key (v1 caches only the parsed
      # case; the network is always built fresh): same single entry, and the
      # option change still takes effect through the fresh build.
      shifted = _net_cache_run(case_path, joinpath(dir, "run_shift"); extra = merge(Dict{String,Any}("matpower_import.shift_sign" => -1.0), enabled))
      @test shifted.status === :succeeded
      @test length(readdir(cache_dir)) == 1

      # Case-file change produces a different key.
      open(case_path, "a") do io
        write(io, "\n% cache invalidation marker\n")
      end
      touched = _net_cache_run(case_path, joinpath(dir, "run_touched"); extra = enabled)
      @test touched.status === :succeeded
      @test length(readdir(cache_dir)) == 2

      # auto_profile != off: cache stays inactive (no third entry).
      recommend = _net_cache_run(case_path, joinpath(dir, "run_recommend"); extra = merge(Dict{String,Any}("matpower_import.auto_profile" => "recommend"), enabled))
      @test recommend.status === :succeeded
      @test length(readdir(cache_dir)) == 2

      # The built net must never mutate the parsed matrices (the cache
      # serializes mpc alongside the net; createNet aliases it since #292).
      mpc = Sparlectra.MatpowerIO.read_case(case_path; legacy_compat = true)
      bus_before = copy(mpc.bus)
      gen_before = copy(mpc.gen)
      branch_before = copy(mpc.branch)
      Sparlectra.createNetFromMatPowerCase(mpc = mpc)
      @test mpc.bus == bus_before
      @test mpc.gen == gen_before
      @test mpc.branch == branch_before
    end
  end
end
