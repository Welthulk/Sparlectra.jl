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

# file: test/test_configuration_coverage.jl
# purpose: tests that every configuration key is reachable and arrives,
#          value-domain validation, case precedence, unknown-key handling,
#          MATPOWER auto-profile decision rules, and console_live capture
using Sparlectra
using Test
using Logging

function _leaf_paths(x; prefix = "")
  paths = String[]
  if x isa AbstractDict
    for key in sort!(collect(keys(x)); by = k -> String(k))
      k = String(key)
      child = x[key]
      child_prefix = isempty(prefix) ? k : string(prefix, ".", k)
      append!(paths, _leaf_paths(child; prefix = child_prefix))
    end
  else
    push!(paths, prefix)
  end
  return paths
end

function _canonical_yaml_leaf_keys()
  yaml = Sparlectra.load_yaml_dict(joinpath(@__DIR__, "..", "src", "config", "configuration.yaml.example"))
  return Set(_leaf_paths(yaml))
end

function _auto_profile_shift_case()
  base = 100.0
  bus = [
    1.0 3.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 110.0 1.0 1.1 0.9
    2.0 1.0 0.0 0.0 0.0 0.0 1.0 1.0 -5.0 110.0 1.0 1.1 0.9
  ]
  branch = [1.0 2.0 0.01 0.10 0.0 999.0 999.0 999.0 1.0 10.0 1.0 -360.0 360.0]
  ybus = Sparlectra.MatpowerIO.build_ybus_matpower(bus, branch, base; matpower_shift_unit = :deg, matpower_shift_sign = 1.0, matpower_ratio = :normal)
  v = bus[:, 8] .* cis.(bus[:, 9] .* (pi / 180.0))
  scalc = v .* conj.(ybus * v)
  bus[2, 3] = -real(scalc[2]) * base
  bus[2, 4] = -imag(scalc[2]) * base
  gen = [1.0 0.0 0.0 999.0 -999.0 1.0 100.0 1.0 999.0 0.0]
  return Sparlectra.MatpowerIO.MatpowerCase("auto_profile_shift", base, bus, gen, branch, nothing, nothing)
end

function _auto_profile_pv_mismatch_case()
  mpc = _auto_profile_shift_case()
  bus = copy(mpc.bus)
  bus[2, 2] = 2.0
  bus[2, 8] = 1.00
  gen = [
    1.0 0.0 0.0 999.0 -999.0 1.0 100.0 1.0 999.0 0.0
    2.0 0.0 0.0 999.0 -999.0 1.04 100.0 1.0 999.0 0.0
  ]
  return Sparlectra.MatpowerIO.MatpowerCase("auto_profile_pv", mpc.baseMVA, bus, gen, mpc.branch, nothing, nothing)
end

# The shipped nine-bus case with its stored VM/VA columns replaced by a
# state that is SOLVED under the given reading of the shift column and the
# ratio: the stored voltages are kept, and the loads of the PQ buses and the
# Pg of the PV buses are set to the injections those voltages produce with
# that reading's Ybus (the same construction as `_auto_profile_shift_case`).
# The exporter cannot write such a file, it always writes the shift in
# degrees with the standard sign, so the fixture is built here. Branch 2
# (the 2-7 step-up transformer) gets a shift of 0.05 in the file column:
# 0.05 rad under the rad readings, 0.05 deg under the deg readings, so the
# readings are distinguishable.
function _auto_profile_case9_solved_under(; shift_unit::Symbol, shift_sign::Float64, ratio::Symbol)
  mpc = Sparlectra.MatpowerIO.read_case(abspath(joinpath(dirname(@__DIR__), "data", "mpower", "sp_case9.m")))
  base = mpc.baseMVA
  bus = copy(mpc.bus)
  gen = copy(mpc.gen)
  branch = copy(mpc.branch)
  branch[2, 10] = 0.05
  for r in axes(bus, 1)
    bus[r, 8] = 1.0 + 0.02 * sin(r)
    bus[r, 9] = Int(bus[r, 2]) == 3 ? 0.0 : 4.0 * cos(r)
  end
  ybus = Sparlectra.MatpowerIO.build_ybus_matpower(bus, branch, base; matpower_shift_unit = shift_unit, matpower_shift_sign = shift_sign, matpower_ratio = ratio)
  v = bus[:, 8] .* cis.(bus[:, 9] .* (pi / 180.0))
  scalc = v .* conj.(ybus * v)
  for r in axes(bus, 1)
    btype = Int(bus[r, 2])
    if btype == 1
      bus[r, 3] = -real(scalc[r]) * base
      bus[r, 4] = -imag(scalc[r]) * base
    elseif btype == 2
      g = findfirst(row -> Int(gen[row, 1]) == Int(bus[r, 1]), axes(gen, 1))
      gen[g, 2] = real(scalc[r]) * base + bus[r, 3]
    end
  end
  return Sparlectra.MatpowerIO.MatpowerCase(string("sp_case9_", shift_unit, "_", shift_sign, "_", ratio), base, bus, gen, branch, nothing, nothing)
end

# The same case with its VM/VA columns perturbed after the construction:
# the stored state then fits NO reading (the situation of the PEGASE
# files), and the scan must not pick a convention from it.
function _auto_profile_case9_perturbed()
  mpc = _auto_profile_case9_solved_under(shift_unit = :deg, shift_sign = 1.0, ratio = :normal)
  bus = copy(mpc.bus)
  for r in axes(bus, 1)
    Int(bus[r, 2]) == 3 && continue
    bus[r, 8] += isodd(r) ? 0.05 : -0.05
    bus[r, 9] += isodd(r) ? 10.0 : -10.0
  end
  return Sparlectra.MatpowerIO.MatpowerCase("sp_case9_perturbed", mpc.baseMVA, bus, mpc.gen, mpc.branch, nothing, nothing)
end

# One function per @testset, called by the runner from a list.
#
# Measured with SnoopCompile on the sysimage: as ONE body these
# nineteen testsets cost 13.242 s of inference in a single MethodInstance
# (n=1), 61 percent of the whole group. The state-estimation group, which is
# already split this way, does three times the inference work (9191 nodes
# against 3027) for less than half the cost in Main, so the shape is the
# expense and not the content.
#
# Nothing is shared between these functions.

function test_configuration_version_scope_and_case_precedence()
  @testset "configuration version, scope, and case precedence" begin (function ()
    dir = mktempdir()
    # a version newer than the running Sparlectra is an error, not a guess
    v9 = joinpath(dir, "v9.yaml")
    write(v9, "config_version: 99\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(v9; reload = true)
    # a case-scope file cannot pose as the general configuration
    vc = joinpath(dir, "vc.yaml")
    write(vc, "config_version: 1\nscope: case\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(vc; reload = true)

    # case configuration file: header binding and scope policing
    case_m = joinpath(dir, "case57.m")
    write(case_m, "% fixture\n")
    cc = Sparlectra.case_config_path(case_m)
    # a MATPOWER case binds <file name>.config.yaml (0.15.0): it never shares
    # a file with the SCF export of the same stem
    @test basename(cc) == "case57.m.config.yaml"
    scf_cc = joinpath(dir, "case57.config.yaml")
    write(cc, "config_version: 1\nscope: case\ncase: other.m\npower_flow:\n  max_iter: 43\n")
    err_mismatch = try
      Sparlectra.load_case_config(case_m)
      nothing
    catch err
      err
    end
    @test err_mismatch isa ArgumentError
    @test occursin("case_config_mismatch", sprint(showerror, err_mismatch))
    # a session-scope key in a case configuration is a hard error
    write(cc, "config_version: 1\nscope: case\ncase: case57.m\noutput:\n  console_summary: false\n")
    err_scope = try
      Sparlectra.load_case_config(case_m)
      nothing
    catch err
      err
    end
    @test err_scope isa ArgumentError
    @test occursin("not case scope", sprint(showerror, err_scope))
    # the canonical .scf.json double extension binds <stem>.config.yaml
    @test basename(Sparlectra.case_config_path(joinpath(dir, "case57.scf.json"))) == "case57.config.yaml"
    # precedence, one key on each level, highest first: override,
    # case configuration file, deprecated in-file block, general file
    general = joinpath(dir, "general.yaml")
    write(general, "config_version: 1\nscope: general\npower_flow:\n  max_iter: 41\n")
    scf_case = joinpath(dir, "case57.scf.json")
    scf_root = string(
      "{\"version\": \"1.0\", \"type\": \"input\", \"is_batch\": false, \"attributes\": {}, \"data\": {}, ",
      "\"sparlectra\": {\"format_version\": \"", Sparlectra.SCF_FORMAT_VERSION, "\", \"config\": {\"power_flow.max_iter\": 42}}}",
    )
    write(scf_case, scf_root)
    write(scf_cc, "config_version: 1\nscope: case\ncase: case57.scf.json\npower_flow:\n  max_iter: 43\n")
    with_override = Sparlectra.resolve_config(general, scf_case, Dict{String,Any}("power_flow.max_iter" => 44))
    @test with_override.config.powerflow.max_iter == 44
    with_case_file = Sparlectra.resolve_config(general, scf_case)
    @test with_case_file.config.powerflow.max_iter == 43
    rm(scf_cc)
    with_block = Sparlectra.resolve_config(general, scf_case)
    @test with_block.config.powerflow.max_iter == 42
    # case-scope keys resolve against the defaults, format independent:
    # the PRESENCE of a case configuration file makes the case
    # self-contained, so a case-scope key the file does not set resolves
    # to the packaged default, NOT to the general file, for every input
    # format alike; machine-scope keys still come from the general file,
    # and a case WITHOUT a config file keeps the full chain
    write(general, "config_version: 1\nscope: general\npower_flow:\n  max_iter: 41\noutput:\n  console_max_rows: 33\n")
    write(scf_case, string("{\"version\": \"1.0\", \"type\": \"input\", \"is_batch\": false, \"attributes\": {}, \"data\": {}, ", "\"sparlectra\": {\"format_version\": \"", Sparlectra.SCF_FORMAT_VERSION, "\"}}"))
    tmpl_max_iter = Sparlectra.SparlectraConfig(Sparlectra.load_yaml_dict(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH)).powerflow.max_iter
    @test tmpl_max_iter != 41
    # no case config file: the general file participates, .scf.json included
    @test Sparlectra.resolve_config(general, scf_case).config.powerflow.max_iter == 41
    # with a case config file that sets ANOTHER key, the missing key falls
    # to the packaged default; the machine-scope key still applies
    write(scf_cc, "config_version: 1\nscope: case\ncase: case57.scf.json\npower_flow:\n  tol: 1.0e-6\n")
    two_stage = Sparlectra.resolve_config(general, scf_case)
    @test two_stage.config.powerflow.max_iter == tmpl_max_iter
    @test two_stage.config.powerflow.tol == 1.0e-6
    @test two_stage.config.output.console_max_rows == 33
    # the same rule on a MATPOWER case: format independent
    plain = joinpath(dir, "plain57.m")
    write(plain, "% fixture\n")
    plain_cc = Sparlectra.case_config_path(plain)
    write(plain_cc, "config_version: 1\nscope: case\ncase: plain57.m\npower_flow:\n  tol: 1.0e-6\n")
    @test Sparlectra.resolve_config(general, plain).config.powerflow.max_iter == tmpl_max_iter
    rm(plain_cc)
    # a case level or an explicit override still beats the default
    write(scf_cc, "config_version: 1\nscope: case\ncase: case57.scf.json\npower_flow:\n  max_iter: 43\n")
    @test Sparlectra.resolve_config(general, scf_case).config.powerflow.max_iter == 43
    @test Sparlectra.resolve_config(general, scf_case, Dict{String,Any}("power_flow.max_iter" => 44)).config.powerflow.max_iter == 44
    rm(scf_cc)
    # the QSTS run keys are case scope like qsts.series.*: the general file
    # reaches a case without a case file, a case file beats it, an override
    # beats the case file. Each level carries a value the level below does
    # not, so a level that is not read shows.
    qsts_keys(r) = (r.config.qsts.mode, r.config.qsts.base_step, r.config.qsts.snapshot_controls)
    general_text(mode, step, controls) = "config_version: 1\nscope: general\npower_flow:\n  max_iter: 41\noutput:\n  console_max_rows: 33\nqsts:\n  mode: $(mode)\n  base_step: $(step)\n  snapshot_controls: $(controls)\n"
    write(general, general_text("snapshots", 2, true))
    @test qsts_keys(Sparlectra.resolve_config(general, scf_case)) == (:snapshots, 2, true)
    write(general, general_text("off", 2, false))
    write(scf_cc, "config_version: 1\nscope: case\ncase: case57.scf.json\nqsts:\n  mode: snapshots\n  base_step: 3\n  snapshot_controls: true\n")
    @test qsts_keys(Sparlectra.resolve_config(general, scf_case)) == (:snapshots, 3, true)
    @test qsts_keys(Sparlectra.resolve_config(general, scf_case, Dict{String,Any}("qsts.mode" => "off", "qsts.base_step" => 4, "qsts.snapshot_controls" => false))) == (:off, 4, false)
    # the keys of the sequential run, the same three levels: defaults from
    # the general file, the case file beats it, an override beats the case
    seq_keys(r) = (r.config.qsts.mode, r.config.qsts.sequential_start, r.config.qsts.resync_every, r.config.qsts.dead_time_steps, r.config.qsts.hunting_window,
                   r.config.qsts.vm_min_pu, r.config.qsts.vm_max_pu, r.config.qsts.loading_max_pct)
    rm(scf_cc)
    @test seq_keys(Sparlectra.resolve_config(general, scf_case)) == (:off, :previous, 0, 0, 3, 0.9, 1.1, 100.0)
    write(general, general_text("off", 2, false) * "  sequential_start: snapshot\n  resync_every: 4\n  dead_time_steps: 1\n")
    @test seq_keys(Sparlectra.resolve_config(general, scf_case)) == (:off, :snapshot, 4, 1, 3, 0.9, 1.1, 100.0)
    write(scf_cc, "config_version: 1\nscope: case\ncase: case57.scf.json\nqsts:\n  mode: sequential\n  resync_every: 8\n  dead_time_steps: 2\n  hunting_window: 5\n  vm_min_pu: 0.95\n  vm_max_pu: 1.05\n  loading_max_pct: 80\n")
    # a case file makes the case self-contained: sequential_start, which it
    # does not set, is the packaged default again, not the general value
    @test seq_keys(Sparlectra.resolve_config(general, scf_case)) == (:sequential, :previous, 8, 2, 5, 0.95, 1.05, 80.0)
    local seq_override = Dict{String,Any}("qsts.sequential_start" => "snapshot", "qsts.resync_every" => 0, "qsts.dead_time_steps" => 3, "qsts.hunting_window" => 4,
                                          "qsts.vm_min_pu" => 0.92, "qsts.vm_max_pu" => 1.08, "qsts.loading_max_pct" => 90.0)
    @test seq_keys(Sparlectra.resolve_config(general, scf_case, seq_override)) == (:sequential, :snapshot, 0, 3, 4, 0.92, 1.08, 90.0)
    # the value domains, on the override path and in the typed parser, each
    # error with the key named
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict{String,Any}("qsts.base_step" => 0))
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict{String,Any}("qsts.mode" => "parallel"))
    @test_throws ArgumentError Sparlectra.QstsConfig(Dict{String,Any}("qsts" => Dict{String,Any}("base_step" => 0)))
    @test_throws ArgumentError Sparlectra.QstsConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "parallel")))
    @test Sparlectra.QstsConfig(Dict{String,Any}("qsts" => Dict{String,Any}("mode" => "sequential"))).mode === :sequential
    for (key, value, words) in [("hunting_window", 2, "qsts.hunting_window must be >= 3"), ("resync_every", -1, "qsts.resync_every must be >= 0"),
                                ("dead_time_steps", -2, "qsts.dead_time_steps must be >= 0"), ("sequential_start", "base", "qsts.sequential_start must be one of"),
                                ("vm_min_pu", -0.1, "qsts.vm_min_pu must be"), ("vm_max_pu", 0.85, "qsts.vm_max_pu must be finite and above qsts.vm_min_pu"),
                                ("loading_max_pct", 0.0, "qsts.loading_max_pct must be")]
      @test_throws words Sparlectra.QstsConfig(Dict{String,Any}("qsts" => Dict{String,Any}(key => value)))
      @test_throws string("qsts.", key) Sparlectra.validate_gui_config_overrides(Dict{String,Any}(string("qsts.", key) => value == 0.85 ? -1.0 : value))
    end
    rm(scf_cc)
    plain = joinpath(dir, "plain57.m")
    write(plain, "% fixture\n")
    with_general = Sparlectra.resolve_config(general, plain)
    @test with_general.config.powerflow.max_iter == 41
    # write_case_config round trip: header written, case scope enforced,
    # empty settings remove the file
    written = Sparlectra.write_case_config(plain, Dict{String,Any}("power_flow.max_iter" => 77); form = Dict{String,Any}("se_tol" => 1.0e-7))
    content = read(written, String)
    @test occursin("scope: case", content)
    @test occursin("case: plain57.m", content)
    @test Sparlectra.load_case_config(plain) == Dict{String,Any}("power_flow.max_iter" => 77)
    @test Sparlectra.resolve_config(general, plain).config.powerflow.max_iter == 77
    # YAML 1.1 lookalikes: the enum value "off" is quoted on write and
    # survives the round trip as a string (bare `off` came back as the
    # boolean false)
    written_off = Sparlectra.write_case_config(plain, Dict{String,Any}("model.auto_profile" => "off"))
    @test occursin("auto_profile: \"off\"", read(written_off, String))
    @test Sparlectra.load_case_config(plain)["model.auto_profile"] == "off"
    Sparlectra.write_case_config(plain, Dict{String,Any}())
    @test !isfile(written)
    # a case file next to a CGMES zip sets the qsts.series keys; a zip
    # without one keeps the defaults (resolve_config does not open the zip,
    # so the zip itself is not created)
    zip_case = joinpath(dir, "series.zip")
    write(Sparlectra.case_config_path(zip_case), "config_version: 1\nscope: case\ncase: series.zip\nqsts:\n  series:\n    min_structural_overlap: 1.0\n    max_gap: 45\n")
    zcfg = Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, zip_case).config
    @test zcfg.qsts.series.min_structural_overlap == 1.0 && zcfg.qsts.series.max_gap == 45.0 && zcfg.qsts.series.max_injection_churn === nothing
    @test Sparlectra.resolve_config(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, joinpath(dir, "other.zip")).config.qsts.series == Sparlectra.QstsSeriesConfig()
  end)() end
  return nothing
end

function test_configuration_yaml_scalar_round_trip()
  @testset "YAML scalar write/read round trip" begin (function ()
    # the writer's self-inverse contract: every string written by
    # _yaml_scalar_text must come back from the repository's own reader
    # with identical type and value (the bare enum value "off" once came
    # back as boolean false). Table per the review: the YAML 1.1 boolean
    # family, null spellings, number forms, dates, colon/hash strings,
    # symbol and list lookalikes, padded values, ordinary identifiers.
    for s in ("off", "on", "yes", "no", "Off", "TRUE", "false", "null", "~", "y", "n", "007", "-5", "+5", ".5", "1e3", "1.0", "0x1f", "2026-09-02", ":symbolish", "[1, 2]", " padded ", "with: colon", "hash # inside", "recommend", "case57.m", "excel_de")
      text = Sparlectra._yaml_scalar_text(s)
      back = Sparlectra.parse_yaml_scalar(text)
      @test back isa AbstractString
      @test String(back) == s
    end
    # non-strings keep their types through the same pair
    for v in (true, false, 42, 1.5)
      @test Sparlectra.parse_yaml_scalar(Sparlectra._yaml_scalar_text(v)) === v
    end
  end)() end
  return nothing
end

function test_configuration_runner_output_modes()
  @testset "Test runner output mode helpers" begin (function ()
    @test selected_test_profile(String[], Dict("SPARLECTRA_TEST_PROFILE" => "extended")) === :extended
    @test selected_test_profile(["fast"], Dict("SPARLECTRA_TEST_PROFILE" => "extended")) === :fast
    @test selected_test_profile(["--verbose", "all"], Dict{String,String}()) === :all
    @test sparlectra_test_verbose(["--verbose"], Dict{String,String}())
    @test sparlectra_test_verbose(String[], Dict("SPARLECTRA_TEST_VERBOSE" => "1"))
    @test !sparlectra_test_verbose(String[], Dict{String,String}())

    quiet_path = test_scratch_path()
    open(quiet_path, "w+") do io
      redirect_stdio(stdout = io) do
        @test quiet_test_output(verbose = false) do
          println("Runtime casefile: hidden")
          return :quiet_result
        end === :quiet_result
      end
      seekstart(io)
      @test isempty(read(io, String))
    end

    # immediate failure surfacing: a captured "Test Failed at" block is
    # replayed right away instead of hiding until the suite's end (the
    # third swallowed-information fix after SKIPPED lines and stage 2);
    # driven with a SYNTHETIC capture, because a real failing @test here
    # would fail this suite
    fail_capture = test_scratch_path()
    write(fail_capture, string(
      "some solver output\n",
      "group name: Test Failed at /repo/test/test_x.jl:42\n",
      "  Expression: occursin(\"a\", \"b\")\n",
      "   Evaluated: occursin(\"a\", \"b\")\n",
      "\n",
      "more output after the block\n",
    ))
    fail_out = IOBuffer()
    @test _surface_test_failures(fail_capture, fail_out) == 1
    fail_text = String(take!(fail_out))
    @test occursin("surfaced immediately", fail_text)
    @test occursin("Test Failed at /repo/test/test_x.jl:42", fail_text)
    @test occursin("Evaluated: occursin", fail_text)
    @test !occursin("more output after the block", fail_text)
    clean_out = IOBuffer()
    write(fail_capture, "all green, nothing to surface\n")
    @test _surface_test_failures(fail_capture, clean_out) == 0
    @test isempty(String(take!(clean_out)))
    rm(fail_capture; force = true)

    verbose_path = test_scratch_path()
    open(verbose_path, "w+") do io
      redirect_stdio(stdout = io) do
        @test quiet_test_output(verbose = true) do
          println("Runtime casefile: visible")
          return :verbose_result
        end === :verbose_result
      end
      seekstart(io)
      @test occursin("Runtime casefile: visible", read(io, String))
    end
  end)() end
  return nothing
end

function test_configuration_matpower_auto_profile_rules()
  @testset "MATPOWER auto-profile decision rules" begin (function ()
    mpc = _auto_profile_shift_case()
    cfg = Sparlectra.SparlectraConfig(Dict(
      "model" => Dict("auto_profile" => "off"),
      "matpower_import" => Dict("shift_unit" => "rad"),
    ))
    off = Sparlectra.matpower_import_auto_profile(mpc, cfg; mode = :off)
    @test isempty(off.rows)
    @test off.config.matpower.shift_unit === :rad

    rec = Sparlectra.matpower_import_auto_profile(mpc, cfg; mode = :recommend)
    @test rec.config.matpower.shift_unit === :rad
    shift_unit_row = only(row for row in rec.rows if row.option == "matpower_import.shift_unit")
    @test shift_unit_row.recommended == "deg"
    @test shift_unit_row.action === :recommend

    apply_cfg = Sparlectra.SparlectraConfig(Dict(
      "model" => Dict("auto_profile" => "apply"),
      "matpower_import" => Dict("shift_unit" => "rad"),
    ))
    applied = Sparlectra.matpower_import_auto_profile(mpc, apply_cfg; mode = :apply)
    @test applied.config.matpower.shift_unit === :deg
    applied_row = only(row for row in applied.rows if row.option == "matpower_import.shift_unit")
    @test applied_row.action === :applied
    @test any(pair -> first(pair) === :shift_unit && last(pair) === :deg, applied.applied)
    convention_rows = Sparlectra._matpower_import_auto_profile_convention_scan(mpc)
    @test length(convention_rows) == 8

    ambiguous_cfg = Sparlectra.SparlectraConfig(Dict("model" => Dict("auto_profile" => "apply")))
    ambiguous = Sparlectra.matpower_import_auto_profile(_auto_profile_pv_mismatch_case(), ambiguous_cfg; mode = :apply)
    compare_row = only(row for row in ambiguous.rows if row.option == "matpower_import.compare_voltage_reference")
    @test compare_row.recommended == "hybrid"
    @test compare_row.action === :applied
    @test ambiguous.config.matpower.compare_voltage_reference === :hybrid
    @test any(pair -> first(pair) === :compare_voltage_reference && last(pair) === :hybrid, ambiguous.applied)

    explicit_start_cfg = Sparlectra.SparlectraConfig(Dict(
      "model" => Dict("auto_profile" => "apply"),
      "power_flow" => Dict(
        "start_mode" => Dict(
          "angle_mode" => "classic",
          "voltage_mode" => "classic",
          "profile_source" => "flat",
        ),
        "qlimits" => Dict("start_mode" => "iteration"),
      ),
    ))
    fragile_scan = mpc -> [
      (shift_unit = :deg, shift_sign = 1.0, ratio = :normal, stats = (; ok = true), score = 0.01),
      (shift_unit = :deg, shift_sign = -1.0, ratio = :normal, stats = (; ok = true), score = 0.01),
      (shift_unit = :deg, shift_sign = 1.0, ratio = :reciprocal, stats = (; ok = true), score = 0.01),
      (shift_unit = :deg, shift_sign = -1.0, ratio = :reciprocal, stats = (; ok = true), score = 0.01),
      (shift_unit = :rad, shift_sign = 1.0, ratio = :normal, stats = (; ok = true), score = 0.01),
      (shift_unit = :rad, shift_sign = -1.0, ratio = :normal, stats = (; ok = true), score = 0.01),
      (shift_unit = :rad, shift_sign = 1.0, ratio = :reciprocal, stats = (; ok = true), score = 0.01),
      (shift_unit = :rad, shift_sign = -1.0, ratio = :reciprocal, stats = (; ok = true), score = 0.01),
    ]
    conservative_mpc = _auto_profile_pv_mismatch_case()
    conservative_mpc.gen[:, 4] .= 1.0
    conservative_mpc.gen[:, 5] .= 0.0
    conservative = Sparlectra.matpower_import_auto_profile(conservative_mpc, explicit_start_cfg; mode = :apply, convention_scan = fragile_scan)
    @test conservative.config.powerflow.start_mode.angle_mode === :classic
    @test conservative.config.powerflow.start_mode.voltage_mode === :classic
    @test conservative.config.powerflow.start_mode.profile_source === :flat
    @test conservative.config.powerflow.qlimits.start_mode === :iteration
    start_angle_row = only(row for row in conservative.rows if row.option == "power_flow.start_mode.angle_mode")
    qlimit_row = only(row for row in conservative.rows if row.option == "power_flow.qlimits.start_mode")
    @test start_angle_row.recommended == "dc"
    @test start_angle_row.action === :skipped
    @test qlimit_row.action === :skipped
    oom_cfg = Sparlectra.SparlectraConfig(Dict(
      "model" => Dict("auto_profile" => "apply"),
      "matpower_import" => Dict("shift_unit" => "rad"),
    ))
    oom_result = @test_logs (:warn, r"matpower_auto_profile_scan_skipped") Sparlectra.matpower_import_auto_profile(
      mpc,
      oom_cfg;
      mode = :apply,
      convention_scan = mpc -> throw(OutOfMemoryError()),
    )
    @test oom_result.config.matpower.shift_unit === :rad
    @test isempty(oom_result.applied)
    @test any(row -> occursin("matpower_auto_profile_scan_skipped", row.reason), oom_result.rows)

    # A stored state that fits NO reading is no evidence about the file's
    # convention (the four PEGASE files score
    # 0.5 to 3.9 pu under their best reading, rad/-1, and the scan picked
    # that wrong reading). Above model.auto_profile_max_fit_pu the scan
    # keeps the configured conventions and says why, in both modes.
    perturbed = _auto_profile_case9_perturbed()
    perturbed_rows = Sparlectra._matpower_import_auto_profile_convention_scan(perturbed)
    @test first(perturbed_rows).score > 0.1
    default_apply_cfg = Sparlectra.SparlectraConfig(Dict("model" => Dict("auto_profile" => "apply")))
    @test default_apply_cfg.model.auto_profile_max_fit_pu == 0.1
    for mode in (:recommend, :apply)
      res = Sparlectra.matpower_import_auto_profile(perturbed, default_apply_cfg; mode = mode)
      @test res.config.matpower.shift_unit === :deg
      @test res.config.matpower.shift_sign == 1.0
      @test res.config.matpower.ratio === :normal
      @test !any(pair -> first(pair) in (:shift_unit, :shift_sign, :ratio, :bus_shunt_model), res.applied)
      for option in ("matpower_import.shift_unit", "matpower_import.shift_sign", "matpower_import.ratio")
        row = only(r for r in res.rows if r.option == option)
        @test row.action === :keep
        @test row.recommended == row.current
        @test startswith(row.reason, "matpower_auto_profile_scan_inconclusive:")
        @test occursin("not a solved state under any reading", row.reason)
        @test occursin(string("best score ", round(first(perturbed_rows).score; digits = 4), " pu"), row.reason)
        @test occursin("model.auto_profile_max_fit_pu = 0.1", row.reason)
      end
    end
    # the compact console line (the default) carries the reason too, the
    # full table has it on every convention row
    inconclusive = Sparlectra.matpower_import_auto_profile(perturbed, default_apply_cfg; mode = :apply)
    io = IOBuffer()
    Sparlectra.write_matpower_import_auto_profile(io, inconclusive, inconclusive.config; casefile = "sp_case9_perturbed.m")
    @test occursin("not a solved state under any reading", String(take!(io)))

    # the decision rule on the study's own numbers (case1354pegase: best
    # rad/-1 at 0.9991 pu, deg/+1 at 13.0 pu): the default threshold keeps
    # deg/+1 and names the score; a threshold above the best score restores
    # the ranking and recommends rad/-1, so the threshold is what decides
    pegase_scan = mpc -> [
      (shift_unit = :rad, shift_sign = -1.0, ratio = :normal, stats = (; ok = true), score = 0.9991),
      (shift_unit = :rad, shift_sign = -1.0, ratio = :reciprocal, stats = (; ok = true), score = 4.2),
      (shift_unit = :rad, shift_sign = 1.0, ratio = :normal, stats = (; ok = true), score = 9.7),
      (shift_unit = :rad, shift_sign = 1.0, ratio = :reciprocal, stats = (; ok = true), score = 11.0),
      (shift_unit = :deg, shift_sign = 1.0, ratio = :normal, stats = (; ok = true), score = 13.0),
      (shift_unit = :deg, shift_sign = -1.0, ratio = :normal, stats = (; ok = true), score = 13.1),
      (shift_unit = :deg, shift_sign = 1.0, ratio = :reciprocal, stats = (; ok = true), score = 15.0),
      (shift_unit = :deg, shift_sign = -1.0, ratio = :reciprocal, stats = (; ok = true), score = 15.2),
    ]
    pegase_kept = Sparlectra.matpower_import_auto_profile(perturbed, default_apply_cfg; mode = :apply, convention_scan = pegase_scan)
    @test pegase_kept.config.matpower.shift_unit === :deg
    @test pegase_kept.config.matpower.shift_sign == 1.0
    @test !any(pair -> first(pair) in (:shift_unit, :shift_sign, :ratio), pegase_kept.applied)
    pegase_unit_row = only(r for r in pegase_kept.rows if r.option == "matpower_import.shift_unit")
    @test pegase_unit_row.action === :keep
    @test occursin("best score 0.9991 pu above model.auto_profile_max_fit_pu = 0.1", pegase_unit_row.reason)
    raised_cfg = Sparlectra.SparlectraConfig(Dict("model" => Dict("auto_profile" => "apply", "auto_profile_max_fit_pu" => 2.0)))
    @test raised_cfg.model.auto_profile_max_fit_pu == 2.0
    pegase_applied = Sparlectra.matpower_import_auto_profile(perturbed, raised_cfg; mode = :apply, convention_scan = pegase_scan)
    @test pegase_applied.config.matpower.shift_unit === :rad
    @test pegase_applied.config.matpower.shift_sign == -1.0
    @test any(pair -> first(pair) === :shift_unit && last(pair) === :rad, pegase_applied.applied)

    # a stored state that IS solved under rad/-1 scores below the threshold
    # under that reading and far above it under the default deg/+1, so the
    # scan still recommends rad/-1 and apply mode still switches to it
    rad_case = _auto_profile_case9_solved_under(shift_unit = :rad, shift_sign = -1.0, ratio = :normal)
    rad_rows = Sparlectra._matpower_import_auto_profile_convention_scan(rad_case)
    @test first(rad_rows).shift_unit === :rad
    @test first(rad_rows).shift_sign == -1.0
    @test first(rad_rows).ratio === :normal
    @test first(rad_rows).score < 1.0e-6
    rad_rec = Sparlectra.matpower_import_auto_profile(rad_case, default_apply_cfg; mode = :recommend)
    @test rad_rec.config.matpower.shift_unit === :deg
    rad_unit_row = only(r for r in rad_rec.rows if r.option == "matpower_import.shift_unit")
    rad_sign_row = only(r for r in rad_rec.rows if r.option == "matpower_import.shift_sign")
    @test rad_unit_row.recommended == "rad"
    @test rad_unit_row.action === :recommend
    @test rad_sign_row.recommended == "-1.0"
    @test rad_sign_row.action === :recommend
    rad_applied = Sparlectra.matpower_import_auto_profile(rad_case, default_apply_cfg; mode = :apply)
    @test rad_applied.config.matpower.shift_unit === :rad
    @test rad_applied.config.matpower.shift_sign == -1.0
    @test rad_applied.config.matpower.ratio === :normal
    @test any(pair -> first(pair) === :shift_unit && last(pair) === :rad, rad_applied.applied)
    @test any(pair -> first(pair) === :shift_sign && last(pair) == -1.0, rad_applied.applied)
  end)() end
  return nothing
end

function test_configuration_tolerance_in_mw()
  @testset "Tolerance can be stated in MW" begin (function ()
    # the form states ONE tolerance value with a unit; the unit decides
    # which of the two keys the request carries
    # pu (the default): the value becomes power_flow.tol and nothing else
    pu_req = SparlectraApp.powerflow_webui_request(
      Dict{String,Any}("casefile" => "case14.m", "power_flow_tol" => "1e-8", "power_flow_tol_unit" => "pu");
      default_output_root = mktempdir())
    @test pu_req["config_overrides"]["power_flow.tol"] == 1.0e-8
    @test !haskey(pu_req["config_overrides"], "power_flow.tol_MW")

    # MW: the same field becomes power_flow.tol_MW, and the per-unit key is
    # NOT sent, so a run can never carry both bounds
    mw_req = SparlectraApp.powerflow_webui_request(
      Dict{String,Any}("casefile" => "case14.m", "power_flow_tol" => "0.001", "power_flow_tol_unit" => "MW");
      default_output_root = mktempdir())
    @test mw_req["config_overrides"]["power_flow.tol_MW"] == 0.001
    @test !haskey(mw_req["config_overrides"], "power_flow.tol")
  end)() end
  return nothing
end

function test_configuration_every_key_arrives()
  @testset "Every configuration key arrives" begin (function ()
    # Four settings were found in one day that exist, are documented, are
    # shown, and do not act. This closes the class instead of the cases: the
    # key list comes from the typed configuration itself, so a key added
    # later cannot be forgotten here.
    #
    # Two properties per key:
    #   REACHABLE - a user YAML carrying the key loads at all. Everything is
    #     validated against configuration.yaml.example, so a key missing
    #     THERE is refused as "unknown" no matter how well documented it is.
    #     That is how 22 state-estimation keys plus power_flow.tol_MW were
    #     unusable while their documentation described them in full.
    #   ARRIVES - a value that differs from the default reaches the typed
    #     configuration. Only checked where a differing value can be formed
    #     mechanically (numbers, booleans); symbols and collections have
    #     value domains a generic test must not invent, so for those the
    #     reachability check is the assertion.
    section_yaml = Dict("powerflow" => "power_flow", "matpower" => "matpower_import",
                        "cgmes" => "cgmes_import", "powsybl" => "powsybl_import", "shortcircuit" => "short_circuit")
    is_leaf(v) = v isa Union{Symbol,AbstractString,Number,Bool,Nothing,AbstractVector,AbstractDict,Tuple}
    function collect_keys!(out, prefix, x)
      for f in fieldnames(typeof(x))
        v = getfield(x, f)
        path = isempty(prefix) ? String(f) : string(prefix, ".", f)
        is_leaf(v) ? push!(out, (path, v)) : collect_keys!(out, path, v)
      end
      return out
    end
    base = Sparlectra.SparlectraConfig(Dict{String,Any}())
    entries = Tuple{String,Any}[]
    for f in fieldnames(typeof(base))
      v = getfield(base, f)
      is_leaf(v) && continue
      String(f) == "user_set_keys" && continue
      collect_keys!(entries, get(section_yaml, String(f), String(f)), v)
    end

    # Explicit exclusions, one reason each. A silent exclusion or an empty
    # reason is not acceptable (same rule the reference-page coverage uses).
    # Three entries, and each one names a mechanism rather than an opinion.
    # The two `sparse` fields used to sit here as "internal, always true";
    # a field that can never be anything else is not a setting, so they were
    # deleted from the configuration structures instead (the old KEY is
    # still refused by name, so a stored file says what happened).
    excluded = Dict(
      "power_flow.islands.mode" => "removed in 0.30.2: the field stays for struct users, a file value is not read",
      "power_flow.qlimits.ignore_q_limits" => "internal inverse of power_flow.qlimits.enabled, never written by a user",
      "power_flow.start_mode.flatstart" => "set through power_flow.flatstart; the constructor files it under start_mode",
      "power_flow.distributed_slack.weights" => "free-form mapping (block style), present in the template and validated separately",
    )
    for (key, reason) in excluded
      @test !isempty(strip(reason))
      @test any(e -> first(e) == key, entries)
    end

    function nest_value(path, value)
      d::Any = value
      for p in Iterators.reverse(split(path, "."))
        d = Dict{String,Any}(String(p) => d)
      end
      return d
    end
    function load_with(path, value)
      file = test_scratch_path(".yaml")
      open(file, "w") do io
        println(io, "config_version: 1")
        Sparlectra._write_yaml_dict(io, nest_value(path, value isa Symbol ? String(value) : value))
      end
      # the nesting variants below include one that is no key of the file;
      # the loader warns about an unknown key and drops it (never an
      # error), so that warning is captured here and only
      # anything else is re-emitted for the runner's warning check
      logger = Test.TestLogger(min_level = Logging.Warn)
      cfg = Logging.with_logger(logger) do
        Sparlectra.load_sparlectra_config(file; reload = true)
      end
      for record in logger.logs
        occursin("Unknown Sparlectra configuration key", record.message) && continue
        @warn record.message
      end
      return cfg
    end
    # the paths carry the YAML name of the section (power_flow), the struct
    # the field name (powerflow); the catch below took the FieldError of the
    # YAML name for a refused value until it named its failure
    yaml_section = Dict(v => k for (k, v) in section_yaml)
    read_key(cfg, path) = begin
      cur = cfg
      for (i, p) in enumerate(split(path, "."))
        cur = getfield(cur, Symbol(i == 1 ? get(yaml_section, p, p) : p))
      end
      cur
    end
    # yaml path of a field path: the constructors file "guard_enabled" under
    # "guard.enabled" and so on, so the underscore variant is tried too
    function yaml_variants(path)
      parts = split(path, ".")
      out = [path]
      last = parts[end]
      i = findfirst('_', last)
      i === nothing || push!(out, join(vcat(parts[1:end-1], [last[1:i-1], last[i+1:end]]), "."))
      return out
    end

    allowed_by_key = Dict{String,Any}()
    for spec in SparlectraApp.WEBUI_OPTION_SPECS
      spec.config_key === nothing && continue
      isempty(spec.allowed_values) && continue
      allowed_by_key[String(spec.config_key)] = spec.allowed_values
    end
    @test length(allowed_by_key) >= 15

    unreachable = String[]
    arrived = 0
    reach_only = 0
    for (path, value) in entries
      haskey(excluded, path) && continue
      reached = false
      for ypath in yaml_variants(path)
        try
          load_with(ypath, value)
          reached = true
        catch err
          # the loader refuses the value under this spelling of the key
          err isa ArgumentError || rethrow()
        end
        reached && break
      end
      reached || (push!(unreachable, path); continue)
      # A differing value where one can be formed without inventing a
      # domain. Numbers and booleans are mechanical. For symbols the domain
      # is NOT guessed: it is read from the option spec that already
      # declares the allowed values, so a symbol key with a form field is
      # checked for arrival too (output.console_auto_profile, one of the
      # four triggers, is exactly such a key).
      spec_values = get(allowed_by_key, path, ())
      differing = value isa Bool ? !value :
                  value isa Integer ? value + 1 :
                  value isa AbstractFloat ? (isfinite(value) ? value + 1.0 : nothing) :
                  value isa Symbol ? begin
                    other = [v for v in spec_values if Symbol(String(v)) !== value]
                    isempty(other) ? nothing : String(first(other))
                  end : nothing
      if differing === nothing
        reach_only += 1
        continue
      end
      moved = false
      for ypath in yaml_variants(path)
        try
          cfg = load_with(ypath, differing)
          got = read_key(cfg, path)
          got isa Symbol && (got = String(got))
          isequal(got, differing) && (moved = true)
        catch err
          # a domain validation may refuse default+1 (ordered knees, ranges);
          # the key is reachable, which is what this loop proves
          err isa ArgumentError || rethrow()
        end
        moved && break
      end
      moved ? (arrived += 1) : (reach_only += 1)
    end

    @test isempty(unreachable)
    # the loop must actually cover the configuration, not a handful of keys
    @test length(entries) >= 200
    # 70 of the 232 reachable keys carry a number or a boolean, so a
    # differing value can be formed mechanically; the rest are symbols and
    # collections whose value domains a generic test must not invent
    @test arrived >= 60
    @test arrived + reach_only + length(excluded) == length(entries)

    # the other direction: every leaf of the template is the YAML spelling of
    # a typed field, apart from these, one reason each
    template_only = Dict(
      "config_version" => "file header, read by the loader before the typed configuration",
      "scope" => "file header (general or case), policed by the loader",
      "power_flow.flatstart" => "the constructor files it under power_flow.start_mode.flatstart",
      "power_flow.qlimits.enabled" => "stored as its inverse power_flow.qlimits.ignore_q_limits",
      "power_flow.qlimits.guard.enabled" => "the field is power_flow.qlimits.guard",
      "extensions.reserved" => "reserved block without a consumer",
    )
    typed_spellings = Set(ypath for (path, _) in entries for ypath in yaml_variants(path))
    @test sort!(collect(setdiff(_canonical_yaml_leaf_keys(), typed_spellings))) == sort!(collect(keys(template_only)))

    # settings checked by name: the four that started this task first
    for (path, differing, reader) in (
        ("power_flow.tol_MW", 0.002, cfg -> cfg.powerflow.tol_MW),
        ("state_estimation.k_suppress", 5.5, cfg -> cfg.state_estimation.k_suppress),
        ("state_estimation.rank_tol_factor", 25.0, cfg -> cfg.state_estimation.rank_tol_factor),
        ("output.console_auto_profile", "full", cfg -> cfg.output.console_auto_profile),
        # symbols without a form spec and a vector: the generic loop above
        # only checks that they load
        ("power_flow.start_mode.profile_source", "flat", cfg -> cfg.powerflow.start_mode.profile_source),
        ("output.logfile_performance", "full", cfg -> cfg.output.logfile_performance),
        ("power_flow.start_mode.blend_lambdas", [0.11, 0.22], cfg -> cfg.powerflow.start_mode.blend_lambdas),
        # the public Q-limit start modes besides the default iteration_or_auto
        ("power_flow.qlimits.start_mode", "iteration", cfg -> cfg.powerflow.qlimits.start_mode),
        ("power_flow.qlimits.start_mode", "auto", cfg -> cfg.powerflow.qlimits.start_mode))
      cfg = load_with(path, differing)
      got = reader(cfg)
      @test (got isa Symbol ? String(got) : got) == differing
    end
    # the logging surface is in the GUI allowlist, which a case file's config
    # block uses too
    for key in ("output.console_summary", "output.console_diagnostics", "output.console_q_limit_events", "output.logfile_diagnostics", "output.logfile_performance", "output.logfile_warnings", "output.result_table_max_rows", "output.startup_latency_hint")
      @test (key, key in Sparlectra.GUI_EDITABLE_CONFIG_KEYS) == (key, true)
    end
    @test haskey(Sparlectra.validate_gui_config_overrides(Dict{String,Any}("output.logfile_warnings" => "table", "output.console_summary" => false, "output.result_table_max_rows" => 50, "output.startup_latency_hint" => false)), "output")
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict{String,Any}("output.logfile_warnings" => "bogus"))
    # the on/off switch: YAML reads a bare on/off as a boolean, a quoted one
    # as text; both arrive, anything else is refused with the key named
    @test load_with("state_estimation.zero_injection", false).state_estimation.zero_injection === :off
    @test load_with("state_estimation.zero_injection", "on").state_estimation.zero_injection === :on
    @test_throws "state_estimation.zero_injection" load_with("state_estimation.zero_injection", "maybe")
    # the sigma of a zero-injection row: default ZERO_INJECTION_SIGMA, a
    # value that is not finite and positive is refused with the key named,
    # in a file and as a Web UI override
    zs_key = "state_estimation.zero_injection_sigma_mva"
    @test Sparlectra.SparlectraConfig(Dict{String,Any}()).state_estimation.zero_injection_sigma_mva == Sparlectra.ZERO_INJECTION_SIGMA
    for bad in (0.0, -1.0e-3, Inf)
      @test_throws zs_key load_with(zs_key, bad)
      @test_throws zs_key Sparlectra.validate_gui_config_overrides(Dict{String,Any}(zs_key => bad))
    end
    # the SE flat start behind phase shifters: default false, a file value
    # and a Web UI override arrive, a scalar where the block belongs and a
    # non-Bool override are refused with the key named
    rot_key = "state_estimation.start_mode.rotate_behind_pst"
    @test Sparlectra.SparlectraConfig(Dict{String,Any}()).state_estimation.start_mode.rotate_behind_pst === false
    @test load_with(rot_key, true).state_estimation.start_mode.rotate_behind_pst === true
    @test_throws "state_estimation.start_mode" Sparlectra.StateEstimationConfig(Dict{String,Any}("state_estimation" => Dict{String,Any}("start_mode" => true)))
    @test Sparlectra.validate_gui_config_overrides(Dict{String,Any}(rot_key => true))["state_estimation"]["start_mode"]["rotate_behind_pst"] === true
    @test_throws rot_key Sparlectra.validate_gui_config_overrides(Dict{String,Any}(rot_key => "maybe"))
    # the tap impedance rule has two values (auto and off, there is no on):
    # a bare YAML off arrives as off, on (bare or quoted) is refused with the
    # key named, a scalar where the block belongs names the block
    tap_key = "state_estimation.tap_estimation.impedance_correction"
    @test load_with(tap_key, false).state_estimation.tap_estimation.impedance_correction === :off
    @test load_with(tap_key, "auto").state_estimation.tap_estimation.impedance_correction === :auto
    @test Sparlectra.SparlectraConfig(Dict{String,Any}()).state_estimation.tap_estimation.impedance_correction === :auto
    @test_throws tap_key load_with(tap_key, "on")
    @test_throws tap_key load_with(tap_key, true)
    @test_throws "state_estimation.tap_estimation" Sparlectra.StateEstimationConfig(Dict{String,Any}("state_estimation" => Dict{String,Any}("tap_estimation" => true)))
    @test_throws tap_key Sparlectra.validate_gui_config_overrides(Dict{String,Any}(tap_key => "on"))
    @test Sparlectra.validate_gui_config_overrides(Dict{String,Any}(tap_key => "off"))["state_estimation"]["tap_estimation"]["impedance_correction"] == "off"
    # the sigma table: block style keyed by the nominal voltage in kV
    # arrives row by row, sorted by level; a bad row is refused with the
    # key, the row key and the reason word
    row(; kw...) = Dict{String,Any}(String(k) => v for (k, v) in kw)
    key = "state_estimation.sigma_by_voltage_level"
    got = load_with(key, Dict{String,Any}("110" => row(p_mw = 1.0, u_kv = 0.3, rel_pct = 1.0), "380" => row(q_mvar = 2.0))).state_estimation.sigma_by_voltage_level
    @test [r.vn_kv for r in got] == [380.0, 110.0]
    @test got[2].p_mw == 1.0 && got[2].u_kv == 0.3 && got[2].rel_pct == 1.0 && got[2].q_mvar === nothing && got[2].i_a === nothing
    @test got[1].q_mvar == 2.0 && got[1].rel_pct == 0.0
    @test isempty(load_with(key, Dict{String,Any}()).state_estimation.sigma_by_voltage_level)
    for (bad, reason) in (
        (Dict{String,Any}("110" => row(p_mw = 0.0)), "row 110: bad_value: p_mw = 0.0 must be greater than zero"),
        (Dict{String,Any}("110" => row(p_mw = 1.0, rel_pct = -0.5)), "row 110: bad_value: rel_pct = -0.5 must not be negative"),
        (Dict{String,Any}("110" => row(rel_pct = 5.0)), "row 110: no_sigma"),
        (Dict{String,Any}("110" => row()), "row 110: no_sigma"),
        (Dict{String,Any}("110" => row(p_mv = 1.0)), "row 110: unknown_column: p_mv"),
        # the former vn_kv column is gone: the nominal voltage is the key
        (Dict{String,Any}("110" => row(vn_kv = 110, p_mw = 1.0)), "row 110: unknown_column: vn_kv"),
        (Dict{String,Any}("level_1" => row(p_mw = 1.0)), "row level_1: bad_vn_kv"),
        (Dict{String,Any}("0" => row(p_mw = 1.0)), "row 0: bad_vn_kv"))
      @test_throws reason load_with(key, bad)
    end
    # in code a key is any positive number or its text; 380 and "380.0"
    # name the same level
    @test [r.vn_kv for r in Sparlectra._parse_sigma_by_voltage_level(Dict{Any,Any}(10.5 => row(p_mw = 1.0), "380" => row(p_mw = 2.0)))] == [380.0, 10.5]
    @test_throws "duplicate_vn_kv: 380 kV is already the level of row 380" Sparlectra._parse_sigma_by_voltage_level(Dict{Any,Any}("380" => row(p_mw = 1.0), "380.0" => row(q_mvar = 1.0)))
    @test_throws "row -5: bad_vn_kv" Sparlectra._parse_sigma_by_voltage_level(Dict{Any,Any}(-5 => row(p_mw = 1.0)))
    # the key the writers use reads back to the same level
    for vn in (380.0, 110.0, 10.5, 0.4, 13.8)
      @test Sparlectra._sigma_level_vn("row", Sparlectra._sigma_level_key(vn)) === vn
    end
    @test Sparlectra._sigma_level_key(380.0) == "380"
    # a flow mapping is refused with the block form, never read as text:
    # as the value of a row, as the whole table, and as a list item line
    flow = test_scratch_path(".yaml")
    for text in ("    380: {p_mw: 2.0}\n", "  sigma_by_voltage_level: {p_mw: 2.0}\n")
      write(flow, string("config_version: 1\nstate_estimation:\n", startswith(text, "    ") ? "  sigma_by_voltage_level:\n" : "", text))
      @test_throws "flow_style" Sparlectra.load_sparlectra_config(flow; reload = true)
    end
    write(flow, "config_version: 1\nstate_estimation:\n  sigma_by_voltage_level:\n    - {p_mw: 2.0}\n")
    @test_throws "named mappings in block style" Sparlectra.load_sparlectra_config(flow; reload = true)
    # an error in a file names the file and the line of the column, of
    # the row, or of the table key: the configuration file and a case
    # configuration file alike
    table = "config_version: 1\nstate_estimation:\n  sigma_by_voltage_level:\n    110:\n      p_mw: 1.0\n    %VN%:\n      q_mvar: %Q%\n"
    at(n) = " ($(basename(flow)) line $(n))"
    for (vn, q, msg, line) in (("0110", "1.0", "row 0110: duplicate_vn_kv: 110 kV is already the level of row 110$(at(4))", 6), ("abc", "1.0", "row abc: bad_vn_kv", 6), ("\"abc\"", "1.0", "row abc: bad_vn_kv", 6),
                               ("110", "1.0", "row 110: duplicate_vn_kv: the key 110 appears again", 6), ("20", "0.0", "row 20: bad_value", 7), ("20", "{x: 1}", "row 20: bad_value", 7))
      write(flow, replace(table, "%VN%" => vn, "%Q%" => q))
      @test_throws msg Sparlectra.load_sparlectra_config(flow; reload = true)
      @test_throws at(line) Sparlectra.load_sparlectra_config(flow; reload = true)
    end
    # a row key with a decimal point (10.5, 0.4) reads to the same kV as a
    # key of digits; a quoted key ("10.5", '380') reads as the bare one; a
    # key with a dot written twice is refused with both lines like any other
    levels = "config_version: 1\nstate_estimation:\n  sigma_by_voltage_level:\n    380:\n      p_mw: 2.0\n    110:\n      p_mw: 1.0\n    10.5:\n      p_mw: 0.5\n    0.4:\n      p_mw: 0.1\n"
    write(flow, levels)
    @test [r.vn_kv for r in Sparlectra.load_sparlectra_config(flow; reload = true).state_estimation.sigma_by_voltage_level] == [380.0, 110.0, 10.5, 0.4]
    plain_table = Sparlectra.load_yaml_dict(flow)["state_estimation"]["sigma_by_voltage_level"]
    @test sort!(collect(keys(plain_table))) == ["0.4", "10.5", "110", "380"]
    for quoted in (replace(levels, "    10.5:" => "    \"10.5\":"), replace(levels, "    380:" => "    '380':"))
      write(flow, quoted)
      @test Sparlectra.load_yaml_dict(flow)["state_estimation"]["sigma_by_voltage_level"] == plain_table
      @test [r.vn_kv for r in Sparlectra.load_sparlectra_config(flow; reload = true).state_estimation.sigma_by_voltage_level] == [380.0, 110.0, 10.5, 0.4]
    end
    # a dotted name is no key (10.5 is), and a top-level key holds no dot
    # even in quotes: the line is refused in the configuration file and in
    # a case configuration file, not dropped as an unknown key or applied
    # past the sigma table checks (row 110 twice)
    mktempdir() do d
      case_m = joinpath(d, "c.m")
      write(case_m, "% fixture\n")
      for q in ("", "\"", "'")
        write(flow, "config_version: 1\n$(q)power_flow.max_iter$(q): 3\n")
        @test_throws "Unsupported YAML line" Sparlectra.load_sparlectra_config(flow; reload = true)
        write(Sparlectra.case_config_path(case_m), "config_version: 1\nscope: case\ncase: c.m\n$(q)power_flow.max_iter$(q): 3\n")
        @test_throws "Unsupported YAML line" Sparlectra.load_case_config(case_m)
        write(Sparlectra.case_config_path(case_m), "config_version: 1\nscope: case\ncase: c.m\n$(q)state_estimation.sigma_by_voltage_level$(q):\n  110:\n    p_mw: 1.0\n  110:\n    q_mvar: 1.0\n")
        @test_throws "Unsupported YAML line" Sparlectra.load_case_config(case_m)
      end
    end
    write(flow, string(levels, "    10.5:\n      q_mvar: 1.0\n"))
    @test_throws "row 10.5: duplicate_vn_kv: the key 10.5 appears again; it is already the level of row 10.5$(at(8))" Sparlectra.load_sparlectra_config(flow; reload = true)
    @test_throws at(12) Sparlectra.load_sparlectra_config(flow; reload = true)
    write(flow, "config_version: 1\nstate_estimation:\n  sigma_by_voltage_level:\n    110:\n      rel_pct: 1.0\n    20: {p_mw: 1.0}\n")
    @test_throws "row 110: no_sigma: 110 kV sets none of p_mw, q_mvar, i_a, u_kv; rel_pct acts only together with one of them$(at(4))" Sparlectra.load_sparlectra_config(flow; reload = true)
    mktempdir() do d
      ccfg = Sparlectra.case_config_path(joinpath(d, "x.m"))
      write(ccfg, "config_version: 1\nscope: case\ncase: x.m\nstate_estimation:\n  sigma_by_voltage_level: {p_mw: 1.0}\n")
      @test_throws "flow_style" Sparlectra.load_case_config(joinpath(d, "x.m"))
      @test_throws " (x.m.config.yaml line 5)" Sparlectra.load_case_config(joinpath(d, "x.m"))
    end
    # the same parser guards a GUI override, and the key is case scope
    @test_throws "bad_value" Sparlectra.validate_gui_config_overrides(Dict{String,Any}(key => Dict{String,Any}("20" => row(u_kv = -1.0))))
    @test Sparlectra.scf_is_case_config_key(key)

    # invalid values of known keys are refused; the start voltage mode names
    # its key and the allowed values
    for (path, value) in (
        ("power_flow.method", "polar"),
        ("benchmark.methods", ["polar"]),
        ("power_flow.start_current_iteration.damping", 2.0),
        ("power_flow.start_mode.voltage_mode", "nonsense"),
        ("power_flow.rectangular_preallocate_workspace", "invalid"))
      @test_throws ArgumentError load_with(path, value)
    end
    @test_throws ["power_flow.start_mode.voltage_mode", "pv_gen_vg"] load_with("power_flow.start_mode.voltage_mode", "nonsense")
    # YAML integers for the thread counts arrive as the text the runtime
    # parses; the generic loop above only forms string values for them
    @test load_with("runtime.julia_threads", 4).runtime.julia_threads == "4"
    @test load_with("runtime.blas_threads", 16).runtime.blas_threads == "16"
    # switching the method keeps the workspace settings
    cfg_method_copy = Sparlectra.SparlectraConfig(
      powerflow = Sparlectra.PowerFlowConfig(
        method = :polar,
        rectangular_workspace_reuse = false,
        rectangular_preallocate_workspace = :on,
        rectangular_workspace_min_buses = 12345,
      ),
    )
    cfg_rect = Sparlectra.with_powerflow_method(cfg_method_copy, :rectangular)
    @test cfg_rect.powerflow.method === :rectangular
    @test cfg_rect.powerflow.rectangular_workspace_reuse === false
    @test cfg_rect.powerflow.rectangular_preallocate_workspace === :on
    @test cfg_rect.powerflow.rectangular_workspace_min_buses == 12345

    # runtime.parallel: values the typed configuration refuses, and a
    # max_tasks override that arrives as a number of tasks
    for (par_key, par_value) in (("max_tasks", "0"), ("max_tasks", "many"), ("min_work_items", 0))
      par_raw = Dict{String,Any}("runtime" => Dict{String,Any}("parallel" => Dict{String,Any}(par_key => par_value)))
      @test_throws ArgumentError Sparlectra.RuntimeConfig(par_raw)
    end
    par_loaded = Sparlectra.load_sparlectra_config(; reload = true, overrides = Dict{String,Any}("runtime" => Dict{String,Any}("parallel" => Dict{String,Any}("max_tasks" => "3"))))
    @test parallel_max_tasks(par_loaded.runtime.parallel) == 3

    # power_flow.solver and the APSLF sections: value domains and the one
    # refused combination (the start generator runs ahead of Newton, so it
    # cannot go with solver = apslf)
    @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("solver" => "polar")))
    @test Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("solver" => "apslf"))).powerflow.solver === :apslf
    @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("apslf" => Dict("order" => 0))))
    @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("apslf_start" => Dict("order" => -1))))
    @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("solver" => "apslf", "apslf_start" => Dict("enabled" => true))))
    hybrid = Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("apslf_start" => Dict("enabled" => true, "order" => 12))))
    @test hybrid.powerflow.apslf_start.enabled && hybrid.powerflow.apslf_start.order == 12 && hybrid.powerflow.solver === :rectangular

    # the minimal YAML reader reads the default `weights: {}` as the string
    # "{}"; the constructor treats it as an empty table, otherwise every
    # default-merged load fails. Block-style weights are user data (bus
    # names) and pass the unknown-key check.
    ds_default = test_scratch_path(".yaml")
    write(ds_default, "config_version: 1\npower_flow:\n  tol: 1.0e-8\n")
    ds0 = Sparlectra.load_sparlectra_config(ds_default; reload = true).powerflow.distributed_slack
    @test !ds0.enabled && isempty(ds0.weights)
    ds_weights = test_scratch_path(".yaml")
    write(ds_weights, "config_version: 1\npower_flow:\n  distributed_slack:\n    enabled: true\n    p_mode: explicit\n    weights:\n      B1: 2.0\n      B2: 1.0\n")
    ds1 = Sparlectra.load_sparlectra_config(ds_weights; reload = true).powerflow.distributed_slack
    @test ds1.enabled && ds1.p_mode == :explicit
    @test ds1.weights == Dict("B1" => 2.0, "B2" => 1.0)
  end)() end
  return nothing
end

function test_configuration_form_defaults()
  @testset "Form defaults do not drift from the code they mirror" begin (function ()
    # a form default that mirrors a configuration key takes its value from
    # there: a change on one side without the other fails here instead of
    # giving a run different numbers than the form promised
    se = Sparlectra.state_estimation_config()
    @test SparlectraApp._webui_option_default("se_max_eliminations") == se.max_eliminations
    @test SparlectraApp._webui_option_default("se_k_suppress") == se.k_suppress
    @test SparlectraApp._webui_option_default("se_k_eliminate") == se.k_eliminate
    @test SparlectraApp._webui_option_default("se_zero_injection_sigma_mva") == se.zero_injection_sigma_mva
    @test SparlectraApp._webui_option_default("se_rotate_behind_pst") == se.start_mode.rotate_behind_pst
    @test SparlectraApp._webui_option_default("se_max_iter") == se.max_iter
    # the system frequency select starts at the configuration default
    @test SparlectraApp._webui_option_default("model_nominal_frequency_hz") == Sparlectra.ModelConfig().nominal_frequency_hz == 50.0
  end)() end
  return nothing
end

function test_configuration_network_parameters_reach_paths()
  @testset "Configuration-derived network parameters reach every construction path" begin (function ()
    # The table from the task, executed. Four paths set these differently once,
    # and a case built through the wrong one ran with hysteresis 0 whatever the
    # configuration said. `bus_shunt_model` is consumed while the shunts are
    # built, so it is checked on the constructed network, not stamped.
    cfg = Sparlectra.SparlectraConfig(Dict(
      "power_flow" => Dict("qlimits" => Dict("cooldown_iters" => 3, "hysteresis_pu" => 0.05),
        "start_mode" => Dict("flatstart" => true)),
      "model" => Dict("bus_shunt_model" => "voltage_dependent_injection", "nominal_frequency_hz" => 60.0),
    ))
    @test cfg.powerflow.qlimits.cooldown_iters == 3
    @test cfg.powerflow.qlimits.hysteresis_pu == 0.05

    # the stamping function itself: one source of truth for all paths
    stamped = Sparlectra.Net(name = "t", baseMVA = 100.0)
    Sparlectra._apply_config_net_parameters!(stamped, cfg)
    @test stamped.cooldown_iters == 3
    @test stamped.q_hyst_pu == 0.05
    @test stamped.flatstart
    # the system frequency is stamped too, unless the case file stated one
    @test stamped.nominal_frequency_hz == 60.0
    stated = Sparlectra.Net(name = "t_stated", baseMVA = 100.0, nominal_frequency_hz = 50.0)
    stated._frequency_from_case = true
    Sparlectra._apply_config_net_parameters!(stated, cfg)
    @test stated.nominal_frequency_hz == 50.0

    # MATPOWER and the case format run through the shared import helper
    warmup = abspath(joinpath(dirname(@__DIR__), "data", "mpower", "warmup_casePST.m"))
    if isfile(warmup)
      mp = Sparlectra._import_sparlectra_net(warmup, nothing, cfg)
      @test mp.cooldown_iters == 3
      @test mp.q_hyst_pu == 0.05
      @test mp.bus_shunt_model === :voltage_dependent_injection
      @test mp.nominal_frequency_hz == 60.0
      d = mktempdir()
      scf = exportSCF(mp; file = joinpath(d, "params.scf.json"))
      sc = Sparlectra._import_sparlectra_net(scf, nothing, cfg)
      @test sc.cooldown_iters == 3
      @test sc.q_hyst_pu == 0.05
      @test sc.nominal_frequency_hz == 60.0
    end

    # CGMES: the importer stamps the same parameters (checked-in sp_casePST
    # delivery); a config of its own keeps the shunt model and frequency
    # of `cfg` away from the CGMES import
    cg_cfg = Sparlectra.SparlectraConfig(Dict("power_flow" => Dict("qlimits" => Dict("cooldown_iters" => 3, "hysteresis_pu" => 0.05))))
    cg = Sparlectra.import_case(cgmes_fixture_dir("sp_casePST"), cg_cfg; requested_format = :cgmes).net
    @test (cg.cooldown_iters, cg.q_hyst_pu) == (3, 0.05)

    # without a configuration the MATPOWER importer takes the template
    # defaults of DEFAULT_QLIMIT_CONFIG
    mpc_plain = (
      name = "case2_defaults",
      baseMVA = 100.0,
      bus = [
        1 3 0.0 0.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9
        2 1 20.0 10.0 0.0 0.0 1 1.0 0.0 110.0 1 1.1 0.9
      ],
      gen = [1 30.0 0.0 999.0 -999.0 1.0 100.0 1 999.0 0.0 0 0 0 0 0 0 0 0 0 0 0],
      branch = [1 2 0.01 0.05 0.0 9999.0 0.0 0.0 0.0 0.0 1 -60.0 60.0],
      gencost = nothing,
      bus_name = nothing,
    )
    plain = Sparlectra.createNetFromMatPowerCase(mpc = mpc_plain, log = false, flatstart = false)
    @test plain.cooldown_iters == Sparlectra.DEFAULT_QLIMIT_CONFIG.cooldown_iters
    @test plain.q_hyst_pu == Sparlectra.DEFAULT_QLIMIT_CONFIG.hysteresis_pu
  end)() end
  return nothing
end

function test_configuration_console_live_capture()
  @testset "console_live capture tees to console and archive identically" begin (function ()
    # live=false: output only in the archive stream. The archive must be a
    # real OS stream here (redirect_stdout rejects IOBuffer) — exactly what
    # production passes (the open run.log IOStream).
    quiet_path, quiet_io = mktemp()
    result = redirect_stdout(devnull) do
      SparlectraApp._capture_run_output(quiet_io) do
        println("captured line")
        42
      end
    end
    close(quiet_io)
    @test result == 42
    @test occursin("captured line", read(quiet_path, String))

    # live=true: the same bytes land in the archive AND on the (outer) console.
    # The read must happen AFTER the redirect block: redirect_stdout dups the
    # pipe's write end onto fd 1, so closing outer.in inside the block leaves
    # the dup open and read() would never see EOF (deadlock).
    archive = IOBuffer()
    outer = Pipe()
    Base.link_pipe!(outer; reader_supports_async = true, writer_supports_async = true)
    redirect_stdout(outer) do
      SparlectraApp._capture_run_output(archive; live = true) do
        println("teed line")
      end
    end
    close(outer.in)
    console_text = read(outer, String)
    @test occursin("teed line", String(take!(archive)))
    @test occursin("teed line", console_text)

    # error path: output written BEFORE a mid-run throw must be flushed into
    # the archive (pump drained in the finally block), and the exception
    # must propagate — exactly what a crash post-mortem needs from run.log
    archive2 = IOBuffer()
    thrown = redirect_stdout(devnull) do
      try
        SparlectraApp._capture_run_output(archive2; live = true) do
          println("before crash")
          error("boom")
        end
        false
      catch err
        occursin("boom", sprint(showerror, err))
      end
    end
    @test thrown
    @test occursin("before crash", String(take!(archive2)))

    # config surface: default off, parseable on
    @test !Sparlectra.OutputConfig().console_live
    @test Sparlectra.OutputConfig(Dict("output" => Dict("console_live" => true))).console_live
  end)() end
  return nothing
end

function test_configuration_unknown_keys()
  @testset "Unknown keys warn in a stored file and fail in an override" begin (function ()
    # an unknown key in a stored file is a warning naming the key,
    # never a failure; the file's other values
    # load, the unknown key is dropped
    bad = test_scratch_path(".yaml")
    write(bad, "diagnostics:\n  no_such_key: 1\n  log_effective_config: true\n")
    bad_cfg = @test_logs (:warn, r"Unknown Sparlectra configuration key diagnostics\.no_such_key is ignored") match_mode = :any Sparlectra.load_sparlectra_config(bad; reload = true)
    @test bad_cfg.diagnostics.log_effective_config
    # a programmatic override with an unknown key is still an error: that is
    # a call, not a stored file
    clean = test_scratch_path(".yaml")
    write(clean, "diagnostics:\n  log_effective_config: true\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(clean; reload = true, overrides = Dict{String,Any}("diagnostics" => Dict{String,Any}("no_such_key" => 1)))
    # a case configuration follows the same rule (the Settings page wrote
    # keys into case files that a later release removed): the unknown key is
    # dropped with a warning, a known key outside the case scope is refused
    case_dir = mktempdir()
    case_m = joinpath(case_dir, "c.m")
    write(case_m, "")
    write(joinpath(case_dir, "c.config.yaml"), "config_version: 1\nscope: case\ncase: c.m\npower_flow:\n  islands:\n    mode: solve_independent\n    enabled: true\n")
    case_level = @test_logs (:warn, r"Unknown Sparlectra configuration key power_flow\.islands\.mode in c\.config\.yaml is ignored") Sparlectra.load_case_config(case_m)
    @test case_level == Dict{String,Any}("power_flow.islands.enabled" => true)
    write(joinpath(case_dir, "c.config.yaml"), "config_version: 1\nscope: case\ncase: c.m\noutput:\n  console_live: true\n")
    @test_throws ArgumentError Sparlectra.load_case_config(case_m)
  end)() end
  return nothing
end

function test_configuration_qlimit_enforcement_keys()
  @testset "Q-limit enforcement mode user YAML keys" begin (function ()
    for mode in (:active_set, :classic_simultaneous, :classic_one_at_a_time)
      cfgfile = test_scratch_path(".yaml")
      write(cfgfile, """
power_flow:
  qlimits:
    enforcement_mode: $(mode)
""")
      cfg = Sparlectra.load_sparlectra_config(cfgfile; reload = true)
      @test cfg.powerflow.qlimits.enforcement_mode === mode
    end
    # any other value is refused, naming the key and the allowed modes
    refused = @test_throws ArgumentError Sparlectra.QLimitConfig(Dict("enforcement_mode" => "classic"))
    @test occursin("power_flow.qlimits.enforcement_mode", sprint(showerror, refused.value))
    @test occursin("classic_simultaneous", sprint(showerror, refused.value))
    # The Q-limit switching parameters of the run configuration reach the net
    # in run_sparlectra. The importer stamps the defaults (hysteresis 0.01 pu,
    # final bound 0.02 pu); the run used to keep them, so the classic outer
    # loop and the final check accepted bus 6 at 24.678 MVAr over its Qmax of
    # 24 MVAr whatever hysteresis_pu said. With hysteresis_pu 0 the bus is
    # clamped: PQ at Qmax, voltage below its set point 1.07 pu. Reference
    # values: the same run with the four net fields set by hand to the
    # configured values. The activeSet_A file is byte-identical to the C file.
    # cooldown_iters and reenable_v_hyst_pu are not read by the classic modes;
    # the third row shows they reach the net as well.
    quiet = Sparlectra.OutputConfig(logfile_results = :off, console_summary = false, startup_latency_hint = false)
    zeng_a = joinpath(dirname(@__DIR__), "data", "scf", "case14_zeng_p306_activeSet_A.scf.json")
    for (extra, type6, q6, vm6) in (("", "PV", 24.678, 1.07), ("    hysteresis_pu: 0.0\n", "PQ", 24.0, 1.0687),
                                    ("    hysteresis_pu: 0.0\n    cooldown_iters: 4\n    reenable_v_hyst_pu: 0.002\n", "PQ", 24.0, 1.0687))
      cfgfile = test_scratch_path(".yaml")
      write(cfgfile, "config_version: 1\npower_flow:\n  rescue: false\n  qlimits:\n    enforcement_mode: classic_one_at_a_time\n" * extra)
      cfg = Sparlectra.load_sparlectra_config(cfgfile; reload = true)
      run_cfg = Sparlectra.SparlectraConfig(powerflow = cfg.powerflow, output = quiet)
      # the default run leaves bus 6 outside its band (warning captured), the
      # clamped runs satisfy the stability conditions and warn about nothing
      r = type6 == "PV" ? (@test_logs (:warn, r"Q-limit stability check") match_mode = :any run_sparlectra(net = Sparlectra.importSCF(zeng_a), config = run_cfg)) :
          (@test_logs min_level = Logging.Warn run_sparlectra(net = Sparlectra.importSCF(zeng_a), config = run_cfg))
      ql = cfg.powerflow.qlimits
      @test (r.net.q_hyst_pu, r.net.final_q_accept_pu, r.net.cooldown_iters, r.net.reenable_v_hyst_pu) == (ql.hysteresis_pu, ql.final_q_accept_pu, ql.cooldown_iters, ql.reenable_v_hyst_pu)
      @test r.final_converged
      Sparlectra.calcNetLosses!(r.net)
      b6 = only(row for row in Sparlectra.buildACPFlowReport(r.net).nodes if row.bus == 6)
      @test (b6.type, round(b6.q_gen_MVar; digits = 3), round(b6.vm_pu; digits = 4)) == (type6, q6, vm6)
    end
    # one default per key: the packaged template, the struct (a library run
    # without a configuration file) and the loader fallback (a user file
    # that omits the key) agree field by field; guard.max_switches was 3 in
    # the template and 10 in the struct up to 0.20.2
    tmpl = Sparlectra.SparlectraConfig(Sparlectra.load_yaml_dict(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH)).powerflow.qlimits
    from_struct = Sparlectra.QLimitConfig()
    from_empty = Sparlectra.QLimitConfig(Dict{String,Any}())
    for f in fieldnames(Sparlectra.QLimitConfig)
      @test (f, getfield(from_struct, f)) == (f, getfield(tmpl, f))
      @test (f, getfield(from_empty, f)) == (f, getfield(tmpl, f))
    end
    # 0.30.0 agreement for the ten release keys: template = struct = loader
    # fallback = the keyword defaults the solver entry points resolve to
    # (read back from the performance profile of a config-less library run
    # on the shipped sp_case118) = the effective configuration of a service
    # run. The keyword defaults read DEFAULT_QLIMIT_CONFIG and
    # DEFAULT_NEWTON_UPDATE; the net and the MATPOWER importer read the
    # hysteresis and cooldown from the same constant, so final_q_accept_pu
    # ("auto" = 2 * hysteresis) agrees too.
    tmpl_cfg = Sparlectra.SparlectraConfig(Sparlectra.load_yaml_dict(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH))
    empty_cfg = Sparlectra.SparlectraConfig(Dict{String,Any}("config_version" => 1))
    expected = (newton_update = :polar, start_iter = 3, start_mode = :iteration_or_auto, guard = true, guard_min_q_range_pu = 0.02,
      guard_narrow_range_mode = :lock_pq, guard_violation_mode = :lock_pq, guard_max_switches = 3, final_q_accept_pu = 0.02, auto_profile = :off)
    for cfg in (tmpl_cfg, empty_cfg, Sparlectra.SparlectraConfig())
      q = cfg.powerflow.qlimits
      @test (cfg.powerflow.newton_update, q.start_iter, q.start_mode, q.guard, q.guard_min_q_range_pu, q.guard_narrow_range_mode, q.guard_violation_mode, q.guard_max_switches, q.final_q_accept_pu, cfg.model.auto_profile) == values(expected)
    end
    @test Sparlectra.DEFAULT_NEWTON_UPDATE === expected.newton_update
    # the batch path (N-1 and scenarios through the contingency service)
    # takes its Q-limit keywords from the same mapping as the single run
    # (_qlimit_solver_kwargs); before 0.30.2 it forwarded none of them. The
    # template and an empty configuration map to the release defaults here,
    # and the service test "scenario service sources (step 5)" shows that a
    # configured lock list changes an N-1 result.
    for cfg in (tmpl_cfg, empty_cfg)
      qk = Sparlectra._qlimit_solver_kwargs(cfg.powerflow.qlimits)
      @test (qk.qlimit_start_iter, qk.qlimit_start_mode, qk.qlimit_guard, qk.qlimit_guard_min_q_range_pu, qk.qlimit_guard_narrow_range_mode, qk.qlimit_guard_violation_mode, qk.qlimit_guard_max_switches) ==
            (expected.start_iter, expected.start_mode, expected.guard, expected.guard_min_q_range_pu, expected.guard_narrow_range_mode, expected.guard_violation_mode, expected.guard_max_switches)
    end
    lib_net = Sparlectra.createNetFromMatPowerFile(filename = joinpath(dirname(@__DIR__), "data", "mpower", "sp_case118.m"))
    prof = Dict{Symbol,Any}()
    redirect_stdout(devnull) do
      runpf!(lib_net, 30, 1e-8, 0; performance_profile = prof)
    end
    @test (prof[:newton_update], prof[:qlimit_start_iter], prof[:qlimit_start_mode], prof[:qlimit_guard], prof[:qlimit_guard_min_q_range_pu], prof[:qlimit_guard_narrow_range_mode], prof[:qlimit_guard_violation_mode], prof[:qlimit_guard_max_switches], prof[:final_q_accept_pu]) ==
          (expected.newton_update, expected.start_iter, expected.start_mode, expected.guard, expected.guard_min_q_range_pu, expected.guard_narrow_range_mode, expected.guard_violation_mode, expected.guard_max_switches, expected.final_q_accept_pu)
    # 0.30.2: the three dishonest-Newton keys agree the same way (template,
    # struct, loader fallback, solver keyword defaults)
    jr_expected = (false, 10.0, 3)
    for cfg in (tmpl_cfg, empty_cfg, Sparlectra.SparlectraConfig())
      pf = cfg.powerflow
      @test (pf.jacobian_reuse, pf.jacobian_reuse_min_reduction, pf.jacobian_reuse_max_steps) == jr_expected
    end
    @test (Sparlectra.DEFAULT_JACOBIAN_REUSE, Sparlectra.DEFAULT_JACOBIAN_REUSE_MIN_REDUCTION, Sparlectra.DEFAULT_JACOBIAN_REUSE_MAX_STEPS) == jr_expected
    @test (prof[:jacobian_reuse], prof[:jacobian_reuse_min_reduction], prof[:jacobian_reuse_max_steps]) == jr_expected
    @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict{String,Any}("config_version" => 1, "power_flow" => Dict{String,Any}("jacobian_reuse_min_reduction" => 1.0)))
    @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict{String,Any}("config_version" => 1, "power_flow" => Dict{String,Any}("jacobian_reuse_max_steps" => 0)))
    # Deliberate exception to the one-default rule above, which covers the
    # solver behaviour: output.logfile_results is log verbosity. A library run
    # without a configuration file stays quiet (off), the packaged template
    # that the Web UI and the services start from logs the full result table.
    # Pinned here so the difference cannot drift unnoticed.
    @test Sparlectra.OutputConfig().logfile_results === :off
    @test empty_cfg.output.logfile_results === :off
    @test tmpl_cfg.output.logfile_results === :full
    # power_flow.power_mode_lu (0.30.2) the same way, and an unknown value is
    # refused by the loader and the override validation
    for cfg in (tmpl_cfg, empty_cfg, Sparlectra.SparlectraConfig())
      @test cfg.powerflow.power_mode_lu === :auto
    end
    @test Sparlectra.DEFAULT_POWER_MODE_LU === :auto
    @test prof[:power_mode_lu] === :auto
    @test_throws ArgumentError Sparlectra.SparlectraConfig(Dict{String,Any}("config_version" => 1, "power_flow" => Dict{String,Any}("power_mode_lu" => "pardiso")))
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict{String,Any}("power_flow.power_mode_lu" => "pardiso"))
  end)() end
  return nothing
end

function test_configuration_value_domain_validation()
  @testset "Configuration value-domain validation" begin (function ()
    tol_bad = test_scratch_path(".yaml")
    write(tol_bad, "power_flow:\n  tol: 0\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(tol_bad; reload = true)

    damping_bad = test_scratch_path(".yaml")
    write(damping_bad, "power_flow:\n  autodamp_min: 0\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(damping_bad; reload = true)

    fit_bad = test_scratch_path(".yaml")
    write(fit_bad, "model:\n  auto_profile_max_fit_pu: -0.1\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(fit_bad; reload = true)

    angle_limit_bad = test_scratch_path(".yaml")
    write(angle_limit_bad, "power_flow:\n  start_mode:\n    dc_angle_limit_deg: 0\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(angle_limit_bad; reload = true)

    # the substitute of the SE range check is never trusted more than the
    # meter it replaces: a factor below 1 is refused with the key named
    factor_bad = test_scratch_path(".yaml")
    write(factor_bad, "state_estimation:\n  bad_data:\n    substitute_sigma_factor: 0.5\n")
    @test_throws "state_estimation.bad_data.substitute_sigma_factor" Sparlectra.load_sparlectra_config(factor_bad; reload = true)

    bool_ok_true = test_scratch_path(".yaml")
    write(bool_ok_true, "power_flow:\n  autodamp: true\n")
    @test Sparlectra.load_sparlectra_config(bool_ok_true; reload = true).powerflow.autodamp === true

    bool_ok_false = test_scratch_path(".yaml")
    write(bool_ok_false, "power_flow:\n  autodamp: false\n")
    @test Sparlectra.load_sparlectra_config(bool_ok_false; reload = true).powerflow.autodamp === false

    enum_bad = test_scratch_path(".yaml")
    write(enum_bad, "power_flow:\n  start_mode:\n    angle_mode: invalid_angle\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(enum_bad; reload = true)
    wrong_branch_bad = test_scratch_path(".yaml")
    write(wrong_branch_bad, "power_flow:\n  wrong_branch_detection: maybe\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(wrong_branch_bad; reload = true)

    # _copy_sparlectra_with_powerflow must carry EVERY field (regression: the
    # former hand-written keyword list silently reset cgmes and webui to
    # defaults on every copy).
    cgmes_cfg = Sparlectra.SparlectraConfig(cgmes = Sparlectra.CGMESImportConfig(start_values = :sv, base_mva = 50.0))
    cgmes_copied = Sparlectra._copy_sparlectra_with_powerflow(cgmes_cfg, cgmes_cfg.powerflow)
    @test cgmes_copied.cgmes.start_values === :sv
    @test cgmes_copied.cgmes.base_mva == 50.0
    @test all(getfield(cgmes_copied, f) == getfield(cgmes_cfg, f) for f in fieldnames(Sparlectra.SparlectraConfig))

    # cgmes_import.start_values: default flat, sv accepted, invalid rejected
    # with the key name in the message.
    # auto is the shipped default: deliveries with an SvVoltage state are
    # solved from it, others fall back to the flat start (resolved per run).
    @test Sparlectra.CGMESImportConfig().start_values === :auto
    flat_forced = test_scratch_path(".yaml")
    write(flat_forced, "cgmes_import:\n  start_values: flat\n")
    @test Sparlectra.load_sparlectra_config(flat_forced; reload = true).cgmes.start_values === :flat
    sv_ok = test_scratch_path(".yaml")
    write(sv_ok, "cgmes_import:\n  start_values: sv\n")
    @test Sparlectra.load_sparlectra_config(sv_ok; reload = true).cgmes.start_values === :sv
    sv_bad = test_scratch_path(".yaml")
    write(sv_bad, "cgmes_import:\n  start_values: warm\n")
    err = try
      Sparlectra.load_sparlectra_config(sv_bad; reload = true)
      nothing
    catch e
      e
    end
    @test err isa ArgumentError
    @test occursin("cgmes_import.start_values", sprint(showerror, err))

    # cgmes_import.placeholder_guards: default warn_skip, strict accepted,
    # invalid rejected with the key name in the message.
    @test Sparlectra.CGMESImportConfig().placeholder_guards === :warn_skip
    pg_ok = test_scratch_path(".yaml")
    write(pg_ok, "cgmes_import:\n  placeholder_guards: strict\n")
    @test Sparlectra.load_sparlectra_config(pg_ok; reload = true).cgmes.placeholder_guards === :strict
    pg_bad = test_scratch_path(".yaml")
    write(pg_bad, "cgmes_import:\n  placeholder_guards: ignore\n")
    pg_err = try
      Sparlectra.load_sparlectra_config(pg_bad; reload = true)
      nothing
    catch e
      e
    end
    @test pg_err isa ArgumentError
    @test occursin("cgmes_import.placeholder_guards", sprint(showerror, pg_err))

    # power_flow.plausibility.vmin_warn: 0 switches the check off and is
    # valid, a negative threshold is refused in a file and as an override
    @test Sparlectra.validate_gui_config_overrides(Dict{String,Any}("power_flow.plausibility.vmin_warn" => 0.0))["power_flow"]["plausibility"]["vmin_warn"] == 0.0
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict{String,Any}("power_flow.plausibility.vmin_warn" => -0.1))
    vw_bad = test_scratch_path(".yaml")
    write(vw_bad, "power_flow:\n  plausibility:\n    vmin_warn: -0.1\n")
    @test_throws ArgumentError Sparlectra.load_sparlectra_config(vw_bad; reload = true)
    # GUI overrides of the start current iteration, merit and trust region:
    # accepted values arrive, out-of-range values and unknown keys are refused
    @test Sparlectra.validate_gui_config_overrides(Dict("power_flow.qlimits.enabled" => false))["power_flow"]["qlimits"]["enabled"] === false
    current_iteration_overrides = Sparlectra.validate_gui_config_overrides(
      Dict(
        "power_flow.start_current_iteration.enabled" => true,
        "power_flow.start_current_iteration.max_iter" => 3,
        "power_flow.start_current_iteration.tol" => 1.0e-4,
        "power_flow.start_current_iteration.damping" => 0.5,
        "power_flow.start_current_iteration.accept_only_if_improved" => true,
        "power_flow.start_current_iteration.min_improvement_factor" => 0.95,
        "power_flow.start_current_iteration.vm_min_pu" => 0.6,
        "power_flow.start_current_iteration.vm_max_pu" => 1.4,
        "power_flow.start_current_iteration.max_angle_step_deg" => 25.0,
        "power_flow.start_current_iteration.only_for_large_cases" => false,
      ),
    )
    @test current_iteration_overrides["power_flow"]["start_current_iteration"]["enabled"] === true
    @test current_iteration_overrides["power_flow"]["start_current_iteration"]["max_iter"] == 3
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict("power_flow.start_current_iteration.damping" => 1.5))

    merit_overrides = Sparlectra.validate_gui_config_overrides(
      Dict("power_flow.merit.enabled" => true, "power_flow.merit.armijo_c1" => 1.0e-3, "power_flow.merit.fallback_max_mismatch" => false),
    )
    @test merit_overrides["power_flow"]["merit"]["enabled"] === true
    @test merit_overrides["power_flow"]["merit"]["armijo_c1"] == 1.0e-3
    @test merit_overrides["power_flow"]["merit"]["fallback_max_mismatch"] === false
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict("power_flow.merit.armijo_c1" => 0.5))
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict("power_flow.merit.scale_p" => 2.0))
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict("power_flow.merit.unknown_key" => true))

    trust_region_overrides = Sparlectra.validate_gui_config_overrides(
      Dict("power_flow.trust_region.enabled" => true, "power_flow.trust_region.initial_radius" => 2.0, "power_flow.trust_region.eta_accept" => 0.2, "power_flow.trust_region.step_mode" => "dogleg"),
    )
    @test trust_region_overrides["power_flow"]["trust_region"]["enabled"] === true
    @test trust_region_overrides["power_flow"]["trust_region"]["initial_radius"] == 2.0
    @test trust_region_overrides["power_flow"]["trust_region"]["eta_accept"] == 0.2
    @test trust_region_overrides["power_flow"]["trust_region"]["step_mode"] == "dogleg"
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict("power_flow.trust_region.initial_radius" => 20.0))
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict("power_flow.trust_region.eta_accept" => 1.0))
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict("power_flow.trust_region.min_radius" => 1.0e-3))
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict("power_flow.trust_region.step_mode" => "bogus"))
    @test_throws ArgumentError Sparlectra.validate_gui_config_overrides(Dict("power_flow.trust_region.unknown_key" => true))

  end)() end
  return nothing
end

function test_configuration_webui_keys_both_directions()
  @testset "Web UI keys are reachable in both directions" begin (function ()
    # tol_MW failed the first direction (a key the UI was allowed to set,
    # with no field to set it in) and stale fields fail the second. Both are
    # asserted, because each direction hides a different defect.
    # collect FIRST, filter second. WEBUI_OPTION_SPECS is an NTuple of 97
    # elements: Julia unrolls tuple iteration, and filtering makes the result
    # length unknown, so a comprehension straight off the tuple ends in
    # Base.grow_to! over a 97-way unrolled generator. Measured,
    # that single line cost 5.61 s of the 5.93 s this group spent inferring
    # generator machinery, and it was the largest item in the whole group.
    # collect() on a tuple with a concrete eltype knows both length and type,
    # and filter() on the resulting Vector is one specialization.
    specs = filter(s -> s.config_key !== nothing, collect(SparlectraApp.WEBUI_OPTION_SPECS))
    have = Set(String(s.config_key) for s in specs)
    @test length(specs) >= 60

    # A: every form field maps to a key that exists AND may be set from the
    # GUI. A typo or a renamed key fails here instead of at run time.
    for spec in specs
      key = String(spec.config_key)
      @test key in Sparlectra.GUI_EDITABLE_CONFIG_KEYS
      probe = try
        Sparlectra.validate_gui_config_overrides(Dict{String,Any}(key => spec.default))
        true
      catch err
        # a value-domain complaint proves the KEY was accepted; only an
        # unknown-key or not-editable complaint is a finding here. Both
        # rejections and every value check throw ArgumentError; anything
        # else is a bug in the validator and must surface.
        err isa ArgumentError || rethrow()
        msg = sprint(showerror, err)
        !occursin("Unknown", msg) && !occursin("not allowed for GUI editing", msg)
      end
      @test probe
    end

    # B: every GUI-editable key either has a field or says why it has none.
    # Without this rule a key can be declared editable and never surface,
    # which is exactly what happened to power_flow.tol_MW.
    fieldless_reasons = Dict(
      "matpower_import.apply_branch_kind" => "import detail, set per case in the case configuration file",
      "matpower_import.apply_branch_names" => "import detail, set per case in the case configuration file",
      "matpower_import.import_for001_contingencies" => "DTF import detail, chosen by the FOR002 selection on the Case page",
      "power_flow.islands.diagnostic_continue_after_failure" => "island diagnostics detail, YAML and API only",
      "power_flow.islands.enabled" => "island solving follows the network, not a form choice",
      "power_flow.islands.reference_policy" => "island solving follows the network, not a form choice",
      "power_flow.method" => "the form offers the solver choice as power_flow.solver; method has one supported value",
      "short_circuit.sweep_method" => "short-circuit performance switch, YAML and API only",
      "short_circuit.converter_current_factor" => "converter fault current without maxIFault, a modelling choice per case; set via configuration.yaml, the case sidecar or the API",
      # the benchmark mode belongs to run_matpower_case in a script (with
      # `using BenchmarkTools`); the Web UI has no benchmark option since 0.20.2
      "benchmark.enabled" => "benchmark mode of run_matpower_case, scripts only; no Web UI option",
      "benchmark.samples" => "benchmark mode of run_matpower_case, scripts only; no Web UI option",
      "benchmark.seconds" => "benchmark mode of run_matpower_case, scripts only; no Web UI option",
      "power_flow.tol_MW" => "set through the tolerance VALUE field plus its unit selector (pu or MW), so one number cannot claim two units",
      "state_estimation.robust" => "issue #377: case scope like the other estimator options, but no dedicated UI control (the form exposes robust_mode instead); set via configuration.yaml or the case sidecar",
      "state_estimation.topology_precheck" => "issue #377: case scope, no dedicated UI control; set via configuration.yaml or the case sidecar",
      "state_estimation.linear_solver" => "estimator linear-algebra backend, a performance switch; set via configuration.yaml or the case sidecar",
      "state_estimation.symmetric_gain" => "estimator linear-algebra switch (Cholesky route for the gain matrix), a performance choice; set via configuration.yaml or the case sidecar",
      "qsts.mode" => "the Web UI QSTS button sets the mode for its own run only; the key itself is set via configuration.yaml or the case sidecar",
      "qsts.base_step" => "QSTS setting, read by the library and the Web UI QSTS run; set via configuration.yaml or the case sidecar",
      "qsts.snapshot_controls" => "QSTS setting, read by the library and the Web UI QSTS run; set via configuration.yaml or the case sidecar",
      "qsts.sequential_start" => "QSTS setting, read by the library and the Web UI QSTS run; set via configuration.yaml or the case sidecar",
      "qsts.resync_every" => "QSTS setting, read by the library and the Web UI QSTS run; set via configuration.yaml or the case sidecar",
      "qsts.dead_time_steps" => "QSTS setting, read by the library and the Web UI QSTS run; set via configuration.yaml or the case sidecar",
      "qsts.hunting_window" => "QSTS setting, read by the library and the Web UI QSTS run; set via configuration.yaml or the case sidecar",
      "qsts.vm_min_pu" => "QSTS setting, read by the library and the Web UI QSTS run; set via configuration.yaml or the case sidecar",
      "qsts.vm_max_pu" => "QSTS setting, read by the library and the Web UI QSTS run; set via configuration.yaml or the case sidecar",
      "qsts.loading_max_pct" => "QSTS setting, read by the library and the Web UI QSTS run; set via configuration.yaml or the case sidecar",
      "qsts.series.min_structural_overlap" => "series check of CGMES snapshots, a library function without a Web UI control; set via configuration.yaml or the case sidecar",
      "qsts.series.max_injection_churn" => "series check of CGMES snapshots, a library function without a Web UI control; set via configuration.yaml or the case sidecar",
      "qsts.series.max_gap" => "series check of CGMES snapshots, a library function without a Web UI control; set via configuration.yaml or the case sidecar",
      "cgmes_import.map_energy_source" => "CGMES import switch for deliveries that carry generation as EnergySource; set via configuration.yaml or the case sidecar",
      "cgmes_import.map_power_electronics" => "CGMES import switch for CGMES 3.0 deliveries with PowerElectronicsConnection; set via configuration.yaml or the case sidecar",
      "cgmes_import.shunt_control" => "CGMES import switch for the bank regulators of a delivery; set via configuration.yaml or the case sidecar",
      "power_flow.plausibility.vmin_warn" => "threshold of the plausibility warning after a converged run; set via configuration.yaml or the case sidecar",
      "state_estimation.sigma_by_voltage_level" => "a table, not a scalar field: the sigma table of the State Estimation section posts it (se_sigma_lvl_* fields), the run takes it as an override and Save settings writes it into the case configuration",
      "state_estimation.tap_estimation.exclude" => "a table, not a scalar field: the tap estimation table of the Case page posts it (se_tap_trafo_* fields) and Save case options writes the unchecked transformers into the case configuration",
      "state_estimation.tap_estimation.include_machine_transformers" => "set by the same Case page table: a checked machine transformer (se_tap_trafo_* fields) makes Save case options write it into the case configuration",
    )
    for key in sort(collect(Sparlectra.GUI_EDITABLE_CONFIG_KEYS))
      key in have && continue
      # the output.* group is console and CSV fine tuning: deliberately not
      # on the run form, otherwise every run would carry twenty display
      # switches. One reason for the whole group.
      startswith(key, "output.") && continue
      @test haskey(fieldless_reasons, key)
      haskey(fieldless_reasons, key) && @test !isempty(strip(fieldless_reasons[key]))
    end
    # and no stale reason: every named key must still be GUI-editable and
    # still have no field, otherwise the list rots
    for key in keys(fieldless_reasons)
      @test key in Sparlectra.GUI_EDITABLE_CONFIG_KEYS
      @test !(key in have)
    end
  end)() end
  return nothing
end

function test_configuration_refresh()
  @testset "Configuration refresh" begin (function ()
    stale = test_scratch_path(".yaml")
    write(stale, "power_flow:\n  tol: 1.0e-6\n  start_mode:\n    voltage_mode: bus_vm_va_blend\n  qlimits:\n    enabled: true\n")
    dry = Sparlectra.refresh_sparlectra_config_file(stale)
    @test dry.success
    @test dry.changed
    @test !dry.written
    @test "power_flow.qlimits.enforcement_mode" in dry.missing_keys
    @test "power_flow.start_mode.voltage_mode" in dry.normalized_keys
    @test occursin("tol: 1.0e-6", dry.refreshed_text)
    @test occursin("voltage_mode: profile_blend", dry.refreshed_text)
    @test occursin("profile_source: matpower_reference", dry.refreshed_text)

    written = Sparlectra.refresh_sparlectra_config_file(stale; write = true)
    @test written.success
    @test written.written
    @test written.backup_path !== nothing
    @test isfile(written.backup_path)
    cfg = Sparlectra.load_sparlectra_config(stale; reload = true)
    @test cfg.powerflow.tol == 1.0e-6
    @test cfg.powerflow.start_mode.voltage_mode === :profile_blend
    @test cfg.powerflow.start_mode.profile_source === :matpower_reference
    @test cfg.powerflow.qlimits.enforcement_mode === :active_set

    # keys no loader reads are deleted by a refresh and named (two keys of
    # earlier versions, one this version never knew, and the form block);
    # the refreshed text loads without an ignored-key warning
    superfluous = test_scratch_path(".yaml")
    write(superfluous, "config_version: 1\nscope: general\npower_flow:\n  tol: 1.0e-7\n  islands:\n    mode: solve_parallel\n  no_such_key: 3\nwebui:\n  warmup: true\nform:\n  gen_seed: 4\n")
    superfluous_result = Sparlectra.refresh_sparlectra_config_file(superfluous)
    @test superfluous_result.success
    @test sort(superfluous_result.removed_keys) == ["form", "power_flow.islands.mode", "power_flow.no_such_key", "webui.warmup"]
    @test count(w -> startswith(w, "Removed power_flow.islands.mode:"), superfluous_result.warnings) == 1
    @test !occursin("solve_parallel", superfluous_result.refreshed_text) && !occursin("no_such_key", superfluous_result.refreshed_text) && !occursin("warmup: true", superfluous_result.refreshed_text)
    @test occursin("tol: 1.0e-7", superfluous_result.refreshed_text) && !occursin("gen_seed", superfluous_result.refreshed_text)
    reloaded = test_scratch_path(".yaml")
    write(reloaded, superfluous_result.refreshed_text)
    reloaded_cfg = @test_logs min_level = Logging.Warn Sparlectra.load_sparlectra_config(reloaded; reload = true)
    @test reloaded_cfg.powerflow.tol == 1.0e-7

    dup = test_scratch_path(".yaml")
    write(dup, "output:\n  detailed_result_csv_exporter: auto\n  detailed_result_csv_exporter: direct\n")
    dup_result = Sparlectra.refresh_sparlectra_config_file(dup; write = true)
    @test !dup_result.success
    @test !dup_result.written
    @test "output.detailed_result_csv_exporter" in dup_result.duplicate_keys
    @test occursin("direct", read(dup, String))
  end)() end
  return nothing
end

# Version-less and legacy fixtures, the in-file config block, a coarse tol_MW
# and the deprecated dcline mode are the tested behavior here, so they warn
# by design. Captured, not printed (a printed warning in a green run
# reads as a problem), and anything else that warns fails the
# group; see run_with_expected_warnings.
const CONFIG_EXPECTED_WARNINGS = (
  r"declares no config_version",
  r"legacy key name",
  r"sparlectra\.config inside .* is deprecated",
  r"is a very coarse convergence bound",
)

function run_configuration_coverage_tests()
  for testfn in (
    test_configuration_version_scope_and_case_precedence,
    test_configuration_yaml_scalar_round_trip,
    test_configuration_runner_output_modes,
    test_configuration_matpower_auto_profile_rules,
    test_configuration_tolerance_in_mw,
    test_configuration_every_key_arrives,
    test_configuration_form_defaults,
    test_configuration_network_parameters_reach_paths,
    test_configuration_console_live_capture,
    test_configuration_unknown_keys,
    test_configuration_qlimit_enforcement_keys,
    test_configuration_value_domain_validation,
    test_configuration_webui_keys_both_directions,
    test_configuration_refresh,
  )
    run_with_expected_warnings(testfn, CONFIG_EXPECTED_WARNINGS)
  end
  return nothing
end
