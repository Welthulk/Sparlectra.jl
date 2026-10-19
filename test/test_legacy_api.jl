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

# file: test/test_legacy_api.jl
# purpose: legacy profile: the testsets of test/test_api.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_api.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_api_tests()
  # cut from testset: API fast smoke and timing contracts
    @testset "net parameters stamped exactly once per importer" begin (function ()
      # the stamping happens exactly
      # once PER IMPORTER, at the place each importer finishes; this test
      # is the guard against pulling the four call sites back together
      # (the CGMES leg runs on the checked-in delivery under
      # data/cgmes_demo).
      cfg = Sparlectra.load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
      repo = dirname(@__DIR__)
      stamps() = Sparlectra._NET_PARAM_STAMP_COUNT[]
      assert_one = function (path, label; requested_format = :auto)
        before = stamps()
        imported = Sparlectra.import_case(path, cfg; requested_format = requested_format)
        @test stamps() - before == 1
        @test imported.net.cooldown_iters == cfg.powerflow.qlimits.cooldown_iters
        println("      import stamp: ", label, " ran (exactly once)")
        return imported
      end
      assert_one(joinpath(repo, "data", "mpower", "warmup_casePST.m"), "matpower")
      assert_one(joinpath(repo, "data", "scf", "sp_case5.scf.json"), "scf")
      dtf = joinpath(repo, "data", "DTF", "FOR001.DAT")
      if isfile(dtf)
        assert_one(dtf, "dtf"; requested_format = :dtf_for001)
      else
        println("      import stamp: dtf SKIPPED (data/DTF/FOR001.DAT not present)")
      end
      assert_one(cgmes_fixture_dir("sp_casePST"), "cgmes (sp_casePST fixture)"; requested_format = :cgmes)
    end)() end

    @testset "run index accepts symlinked output roots" begin (function ()
      # Two environments can reach the same physical output root under
      # different names (a Flatpak XDG_STATE_HOME symlinked onto
      # ~/.local/state). Index entries store the writer's absolute paths,
      # so validation must resolve symlinks instead of comparing lexically,
      # in both directions: entry written under the real name and read
      # through the alias, and the reverse.
      if Sys.iswindows()
        @info "run-index symlink test SKIPPED on Windows (symlink creation may need privileges)"
      else
        mktempdir() do tmpdir
          real_root = joinpath(tmpdir, "real_root")
          run_id = "11111111-2222-3333-4444-555555555555"
          mkpath(joinpath(real_root, run_id))
          write(joinpath(real_root, run_id, "result.json"), "{}")
          alias_root = joinpath(tmpdir, "alias_root")
          symlink(real_root, alias_root)
          real_entry = Dict{String,Any}("run_id" => run_id, "output_dir" => joinpath(real_root, run_id), "result_file" => joinpath(real_root, run_id, "result.json"))
          alias_entry = Dict{String,Any}("run_id" => run_id, "output_dir" => joinpath(alias_root, run_id), "result_file" => joinpath(alias_root, run_id, "result.json"))
          @test SparlectraApp._indexed_run_paths(real_entry, alias_root).valid
          @test SparlectraApp._indexed_run_paths(alias_entry, real_root).valid
          # a symlink must still not smuggle a foreign directory into the root
          outside = joinpath(tmpdir, "outside")
          mkpath(outside)
          write(joinpath(outside, "result.json"), "{}")
          foreign_entry = Dict{String,Any}("run_id" => run_id, "output_dir" => outside, "result_file" => joinpath(outside, "result.json"))
          @test !SparlectraApp._indexed_run_paths(foreign_entry, real_root).valid
        end
      end
    end)() end

  return true
end
