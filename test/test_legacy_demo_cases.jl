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

# file: test/test_legacy_demo_cases.jl
# purpose: legacy profile: the testsets of test/test_demo_cases.jl that pin an earlier
#          change (an issue, a bug, a bit-identity or off-path guarantee).
#          They moved here unchanged from test_demo_cases.jl; the helpers they use
#          are defined in that file, which the legacy group includes first.
using Sparlectra
using Test

function run_legacy_demo_cases_tests()
  # cut from testset: Shipped demo cases (regression fixtures)
    # Why there is no byte round trip for these four, in case someone reaches
    # for one here: importing a shipped case and exporting it again does NOT
    # reproduce the file. `scenarios`, the `short_circuit` block with its PGM
    # `fault` row, the `extra` names and the measurement provenance are EXPORT
    # ARGUMENTS, not network state (measured: sp_case5 comes back as
    # 15354 characters against the shipped 16643). Handing those blocks back to
    # the exporter from the file under test would compare it against its own
    # input. These files are built by a generator that passes those arguments
    # and is their provenance record, not by a re-export; the property such a
    # comparison would be reaching for is checked directly below instead.
    #
    # Review note: after the 0.10.0 -> 0.11.0 bump turned a
    # fixture comparison red: a test that goes red on every version bump
    # trains the reflex to regenerate the file, and the next time the diff
    # may be more than the stamp line, with the regeneration hiding it. The
    # shipped cases must therefore not depend on the release recorded in
    # them at all. A copy stamped with a version that does not exist has to
    # import to the same network and solve to the same result; nothing may
    # read `created_by` except a human looking for provenance.
    @testset "the provenance stamp does not reach behavior" begin (function ()
      d = mktempdir()
      for name in _DEMO_CASE_NAMES
        original = joinpath(_DEMO_CASE_DIR, "$(name).scf.json")
        # the case config next to the file pins the tolerance; it has to
        # travel with the copy or the comparison resolves differently
        # (sp_case6tap ships without one and runs on the defaults)
        cfg = Sparlectra.case_config_path(original)
        isfile(cfg) && cp(cfg, joinpath(d, basename(cfg)); force = true)
        foreign = scf_with_foreign_stamp(original, d)
        @test occursin("Sparlectra 0.0.0-nonexistent", read(foreign, String))
        # and it is the ONLY difference
        @test scf_without_stamp(read(foreign, String)) == scf_without_stamp(read(original, String))

        ref = Sparlectra.import_case(original, _DEMO_REF_CONFIG()).net
        alt = Sparlectra.import_case(foreign, _DEMO_REF_CONFIG()).net
        @test length(alt.nodeVec) == length(ref.nodeVec)
        @test length(alt.branchVec) == length(ref.branchVec)
        @test length(alt.measurements) == length(ref.measurements)
        ite_ref, erg_ref = runpf!(ref; verbose = 0)
        ite_alt, erg_alt = runpf!(alt; verbose = 0)
        @test (ite_alt, erg_alt) == (ite_ref, erg_ref)
        calcNetLosses!(ref)
        calcNetLosses!(alt)
        p_ref, q_ref = Sparlectra.getTotalLosses(net = ref)
        p_alt, q_alt = Sparlectra.getTotalLosses(net = alt)
        @test p_alt == p_ref && q_alt == q_ref
      end
    end)() end

  return true
end
