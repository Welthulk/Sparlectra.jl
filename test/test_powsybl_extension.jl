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

# file: test/test_powsybl_extension.jl
# purpose: the PythonCall extension of the PowSyBl adapter: with PythonCall
#          on the load path and pypowsybl importable, the live read of the
#          shipped ieee14.xiidm produces the same tables as the committed
#          bundle that tools/powsybl_dump.py wrote (every column equal,
#          manifest fields equal except the timestamp and the source path,
#          reference buses equal by value). Skips with a printed reason when
#          PythonCall is not available or pypowsybl cannot be imported;
#          profile `install` (the environment-dependent profile).

using Sparlectra
using Test
using DelimitedFiles

function _powsybl_extension_available()
  Base.find_package("PythonCall") === nothing && return (false, "PythonCall is not on the load path")
  try
    @eval Main using PythonCall
  catch err
    return (false, "PythonCall could not be loaded: $(sprint(showerror, err))")
  end
  try
    # the binding was created by the @eval above, so it is read in the
    # latest world, not in the one this function was compiled in
    pc = Base.invokelatest(getfield, Main, :PythonCall)
    Base.invokelatest(pc.pyimport, "pypowsybl")
  catch err
    return (false, "pypowsybl could not be imported: $(first(sprint(showerror, err), 200))")
  end
  return (true, "")
end

function run_powsybl_extension_tests()
  @testset "PowSyBl PythonCall extension" begin (function ()
    ok, reason = _powsybl_extension_available()
    if !ok
      println("      powsybl extension: SKIPPED ($(reason))")
      @test !hasmethod(Sparlectra.read_powsybl_network, Tuple{String})
      return
    end
    println("      powsybl extension: RAN")
    # the extension was loaded by the availability check above, in a newer
    # world than this function: every method lookup and call goes through
    # invokelatest
    @test Base.invokelatest(hasmethod, Sparlectra.read_powsybl_network, Tuple{String})
    fixture = joinpath(@__DIR__, "fixtures", "powsybl", "ieee14.powsybl")
    xiidm = joinpath(fixture, "ieee14.xiidm")
    committed = read_powsybl_bundle(fixture)
    mktempdir() do dir
      out = joinpath(dir, "ieee14.powsybl")
      manifest_path = Base.invokelatest(Sparlectra.dump_powsybl_bundle, xiidm, out)
      @test isfile(manifest_path)
      live = read_powsybl_bundle(out)
      for name in Sparlectra.POWSYBL_TABLE_NAMES
        a = Sparlectra.powsybl_table(committed, name)
        b = Sparlectra.powsybl_table(live, name)
        @test keys(a) == keys(b)
        for c in keys(a)
          @test isequal(a[c], b[c])
        end
      end
      for key in ("case", "pypowsybl_version", "iidm_version", "all_attributes", "format", "format_version")
        @test get(live.manifest, key, nothing) == get(committed.manifest, key, nothing)
      end
      # the reference rows equal by value; the bytes may differ where Julia
      # and Python print an exponent differently (1.0e-5 against 1e-05)
      ref_live, hdr_live = readdlm(joinpath(out, "reference_buses.csv"), ',', String; quotes = true, header = true)
      ref_committed, hdr_committed = readdlm(joinpath(fixture, "reference_buses.csv"), ',', String; quotes = true, header = true)
      @test vec(hdr_live) == vec(hdr_committed)
      @test size(ref_live) == size(ref_committed)
      for i in axes(ref_committed, 1), j in axes(ref_committed, 2)
        a = tryparse(Float64, ref_committed[i, j])
        b = tryparse(Float64, ref_live[i, j])
        @test a === nothing ? ref_live[i, j] == ref_committed[i, j] : (b !== nothing && isequal(a, b))
      end
      # the live read returns the same tables without writing anything
      tables = Base.invokelatest(Sparlectra.read_powsybl_network, xiidm)
      @test isequal(tables.lines.p1, committed.lines.p1)
      # and import_case takes the .xiidm through the extension
      cfg = load_sparlectra_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH; reload = true)
      ic = Base.invokelatest(Sparlectra.import_case, xiidm, cfg)
      @test ic.format == :powsybl
      @test length(ic.net.nodeVec) == 14
    end
  end)() end
end
