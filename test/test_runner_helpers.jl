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

# file: test/test_runner_helpers.jl
# purpose: helpers for the test runner: profile and verbosity selection from
#          ARGS/ENV, quiet_test_output, which captures a test group's
#          output and replays a bounded excerpt on failure, and
#          run_with_expected_warnings, which captures the advisories a test
#          triggers on purpose and asserts that nothing else was logged
function selected_test_profile(args = ARGS, env = ENV)
  profile_args = [strip(String(arg)) for arg in args if !startswith(strip(String(arg)), "--")]
  cli_profile = isempty(profile_args) ? nothing : Symbol(first(profile_args))
  env_profile = Symbol(get(env, "SPARLECTRA_TEST_PROFILE", "fast"))
  return something(cli_profile, env_profile)
end

function sparlectra_test_verbose(args = ARGS, env = ENV)::Bool
  any(arg -> strip(String(arg)) == "--verbose", args) && return true
  return lowercase(strip(get(env, "SPARLECTRA_TEST_VERBOSE", ""))) in ("1", "true", "yes", "on")
end

const QUIET_TEST_OUTPUT_MAX_BYTES = 64 * 1024
const QUIET_TEST_OUTPUT_MAX_LINES = 200

function _bounded_test_output_excerpt(path::AbstractString; max_bytes::Int = QUIET_TEST_OUTPUT_MAX_BYTES, max_lines::Int = QUIET_TEST_OUTPUT_MAX_LINES)
  isfile(path) || return ""
  total_bytes = filesize(path)
  # read(io, String) has no maxbytes keyword; the bounded byte read plus
  # String keeps the excerpt cap without crashing the failure path (this
  # helper only runs when a group already failed, which is exactly when
  # its own error must not eat the real one)
  text = open(path, "r") do io
    String(read(io, min(Int(total_bytes), max_bytes)))
  end
  lines = split(text, '\n'; keepempty = true)
  omitted_lines = max(length(lines) - max_lines, 0)
  if length(lines) > max_lines
    keep_head = max_lines ÷ 2
    keep_tail = max_lines - keep_head
    lines = vcat(lines[1:keep_head], ["... $(omitted_lines) captured lines omitted ..."], lines[end - keep_tail + 1:end])
  end
  omitted_bytes = max(total_bytes - sizeof(text), 0)
  omitted_bytes > 0 && push!(lines, "... $(omitted_bytes) captured bytes omitted ...")
  return join(lines, '\n')
end

# Surface @test failures from a captured group IMMEDIATELY (third instance
# of the swallowed-information class, after the SKIPPED lines and the
# stage-2 incident): a plain @test failure throws nothing until the outer
# suite aggregates, so the capture would hide the details until the end
# and force a verbose rerun. This scans the capture for failure headers
# and replays each block (header plus following lines up to a blank line
# or the cap) to `err_io` right under the group's result line.
const QUIET_TEST_FAILURE_MAX_BLOCKS = 10
const QUIET_TEST_FAILURE_BLOCK_LINES = 25

function _surface_test_failures(path::AbstractString, err_io::IO; max_blocks::Int = QUIET_TEST_FAILURE_MAX_BLOCKS, block_lines::Int = QUIET_TEST_FAILURE_BLOCK_LINES)
  isfile(path) || return 0
  blocks = 0
  open(path, "r") do io
    remaining = 0
    for line in eachline(io)
      if occursin("Test Failed at", line) || occursin("Error During Test at", line)
        blocks += 1
        if blocks > max_blocks
          println(err_io, "        ... more failure blocks omitted (cap $(max_blocks)); full detail in the suite summary or a --verbose rerun ...")
          break
        end
        blocks == 1 && println(err_io, "      captured test failures (surfaced immediately):")
        remaining = block_lines
      end
      if remaining > 0
        println(err_io, "        ", line)
        remaining -= 1
        isempty(strip(line)) && (remaining = 0)
      end
    end
  end
  return blocks
end

function quiet_test_output(f::Function; verbose::Bool = sparlectra_test_verbose(), group::AbstractString = "test group", skipped_out::Union{Nothing,Vector{String}} = nothing)
  verbose && return Base.invokelatest(f)
  path = tempname()
  result = nothing
  failed = false
  try
    open(path, "w+") do io
      result = redirect_stdio(stdout = io) do
        Base.invokelatest(f)
      end
    end
    # surface availability-gated SKIPPED report lines: the capture would
    # otherwise swallow exactly the "silent skip" the testing conventions
    # forbid (a gated testset states SKIPPED, but nobody saw it here)
    skipped_out === nothing || open(path, "r") do io
      for line in eachline(io)
        occursin("SKIPPED", line) && push!(skipped_out, String(strip(line)))
      end
    end
    # surface @test failures NOW instead of at the suite's end (they throw
    # nothing here, so the catch path below never sees them)
    _surface_test_failures(path, stderr)
  catch err
    if err isa InterruptException
      println(stderr, "Interrupted while running ", group, "; captured output excerpt follows.")
      excerpt = _bounded_test_output_excerpt(path)
      !isempty(excerpt) && print(excerpt)
      rethrow()
    end
    failed = true
    excerpt = _bounded_test_output_excerpt(path)
    !isempty(excerpt) && print(excerpt)
    rethrow()
  finally
    ispath(path) && rm(path; force = true)
  end
  failed && print(_bounded_test_output_excerpt(path))
  return result
end

function run_profile_group(i::Int, total::Int, name::AbstractString, runner::Function)
  print("[", i, "/", total, "] ", name, " ... ")
  skipped = String[]
  timed = try
    @timed quiet_test_output(runner; group = name, skipped_out = skipped)
  catch err
    err isa InterruptException && (println("INTERRUPTED"); rethrow())
    println("FAIL")
    rethrow()
  end
  @printf("PASS %.3f s, %.1f MiB allocated, %.3f s GC\n", timed.time, timed.bytes / 1024.0^2, timed.gctime)
  for line in skipped
    println("        ", line)
  end
  get(ENV, "SPARLECTRA_TEST_GC_BETWEEN_GROUPS", "0") in ("1", "true", "yes", "on") && GC.gc()
  return timed.value
end

"""
    run_with_expected_warnings(f, patterns) -> Any

Run `f()` with its warnings CAPTURED instead of printed, and assert that every
captured warning matches one of `patterns`.

Several subsystems warn by design (a topology advisory, an excluded
measurement, a deliberately incomplete case), and the tests that exercise
them printed those lines on every run: 55 of 62 warnings in the fast profile
came from three such places and buried everything else. Capturing them keeps
the log readable AND checks more than printing ever did, because a warning
that is not on the expected list fails the test set it came from.
"""
function run_with_expected_warnings(f, patterns)
  logger = Test.TestLogger(min_level = Logging.Warn)
  value = Logging.with_logger(f, logger)
  for record in logger.logs
    record.level < Logging.Warn && continue
    expected = any(pattern -> occursin(pattern, record.message), patterns)
    # the tuple form puts the offending message into the failure output
    Test.@test (record.message, expected) == (record.message, true)
  end
  return value
end


"""
    load_fixture_net(name) -> Net

One shared fixture import per process: the named case is imported ONCE
with pinned packaged-default import options (so the fixture never depends
on whatever configuration a previous testset installed) and every call
returns a deepcopy, so testsets can mutate their copy freely.

Since load_fixture_net (2026-09-04) the first-choice names are the SHIPPED
demo cases sp_case5, sp_case14, sp_case60, sp_case188: a fresh install
loads them without any download. The legacy names case9, case14, case57
resolve against the LOCAL data/mpower cache only and exist for the
testsets that guard externally anchored MATPOWER reference values; their
callers must gate on `fixture_net_available` and speak their skip. CGMES
fixtures (MiniGrid, FullGrid) stay per-site on purpose: their imports
exercise importer options and ARE the test subject, a shared cached net
would test the cache instead.
"""
const _FIXTURE_NET_CACHE = Dict{String,Any}()

function load_fixture_net(name::AbstractString)
  net = get!(_FIXTURE_NET_CACHE, String(name)) do
    _import_fixture_net(String(name))
  end
  return deepcopy(net)
end

# legacy MATPOWER fixtures are cache-only; callers gate on this and print
# a spoken SKIPPED line instead of downloading
fixture_net_available(name::AbstractString) = startswith(String(name), "sp_case") || isfile(joinpath(Sparlectra.MPOWER_DIR, string(name, ".m")))

function _import_fixture_net(name::String)
  if startswith(name, "sp_case")
    path = joinpath(dirname(Sparlectra.MPOWER_DIR), "scf", string(name, ".scf.json"))
    isfile(path) || error(string("unknown fixture net: ", name, " (shipped cases: sp_case5, sp_case14, sp_case60, sp_case188)"))
    return Sparlectra.importSCF(path)
  end
  name in ("warmup_casePST", "case9", "case14", "case57") || error(string("unknown fixture net: ", name, " (known: sp_case5, sp_case14, sp_case60, sp_case188 and warmup_casePST shipped; case9, case14, case57 cache-only)"))
  path = joinpath(Sparlectra.MPOWER_DIR, string(name, ".m"))
  isfile(path) || error(string("fixture case file missing (cache-only legacy fixture, gate on fixture_net_available): ", path))
  return Sparlectra.createNetFromMatPowerFile(filename = path, flatstart = false, enable_pq_gen_controllers = true, bus_shunt_model = :admittance, matpower_shift_sign = 1.0, matpower_shift_unit = :deg, matpower_ratio = :normal, tap_changer_model = :ideal)
end

## Scratch paths for the suite. `tempname()` names a file in /tmp that
## NOBODY ever removes: a full run left about three thousand stray
## `jl_*.yaml` files behind, and the configuration warnings they triggered
## buried the real test output (maintainer, 2026-09-05). Everything the
## tests write goes into one directory instead, which Julia deletes when
## the process ends.
const TEST_SCRATCH_ROOT = mktempdir(; cleanup = true)
let counter = Ref(0)
  global function test_scratch_path(extension::AbstractString = "")
    counter[] += 1
    return joinpath(TEST_SCRATCH_ROOT, string("scratch_", counter[], extension))
  end
end
