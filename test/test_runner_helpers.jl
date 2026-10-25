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
#          triggers on purpose and asserts that nothing else was logged,
#          parallel_comparison_report, the RAN/SKIPPED line of every
#          parallel-vs-serial comparison, cgmes_report_plain, the
#          rendered text of the CGMES import report, and run_csv, the cells
#          of a run's CSV in whatever format the locale selected
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

const PARALLEL_COMPARISON_TAG = "parallel comparison "

"""
    parallel_comparison_report(name)

Print the report line of the parallel-vs-serial comparison `name`: RAN with
the thread count, or SKIPPED in a single-threaded process. With one thread
every parallel site falls back to its serial loop, so the comparison checks
the serial path against itself and proves nothing about the threads. The
comparisons run threaded, five times in one process, in the `threads` gate
(`sh tools/run_gates.sh threads`). The runner shows these lines under the
group in every profile.
"""
function parallel_comparison_report(name::AbstractString)
  if Threads.nthreads() > 1
    println(PARALLEL_COMPARISON_TAG, name, ": RAN with ", Threads.nthreads(), " threads")
  else
    println(PARALLEL_COMPARISON_TAG, name, ": SKIPPED, single-threaded process (both runs take the serial path); threaded in tools/run_gates.sh threads")
  end
  return nothing
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
    # the same for the report line of every parallel-vs-serial comparison
    # (parallel_comparison_report): RAN or SKIPPED, it is shown under the
    # group, so a run says which comparisons were threaded; and the RAN line
    # of an availability-gated testset, so a run with the cache says that the
    # gated check ran, not only a run without it that it skipped
    skipped_out === nothing || open(path, "r") do io
      for line in eachline(io)
        (occursin("SKIPPED", line) || occursin(r"\bRAN\b", line) || startswith(strip(line), PARALLEL_COMPARISON_TAG)) && push!(skipped_out, String(strip(line)))
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

# Failed and errored @tests recorded so far in the enclosing testset, its
# finished child testsets included (the counts the test summary prints).
# Outside any testset the fallback testset counts nothing, but there a
# failing @test throws and the group line says FAIL through the catch.
function _testset_nonpass_counts()
  tc = Test.get_test_counts(Test.get_testset())
  return (tc.fails + tc.cumulative_fails, tc.errors + tc.cumulative_errors)
end

function run_profile_group(i::Int, total::Int, name::AbstractString, runner::Function)
  print("[", i, "/", total, "] ", name, " ... ")
  skipped = String[]
  # a failing @test inside the group throws nothing, its result goes to the
  # enclosing testset; the counts before and after the group give the
  # group's own failures and errors (without them the line said PASS)
  fails0, errors0 = _testset_nonpass_counts()
  timed = try
    @timed quiet_test_output(runner; group = name, skipped_out = skipped)
  catch err
    err isa InterruptException && (println("INTERRUPTED"); rethrow())
    println("FAIL")
    rethrow()
  end
  fails1, errors1 = _testset_nonpass_counts()
  fails, errors = fails1 - fails0, errors1 - errors0
  verdict = fails + errors == 0 ? "PASS" : "FAIL ($(fails) failed, $(errors) errored)"
  # compile_time and recompile_time come from @timed since Julia 1.11 and are
  # the whole point of this line: without them a slow group is indistinguishable
  # from a group that merely compiled a lot, and the two need opposite remedies.
  @printf("%s %.3f s (%.3f s compile, %.3f s recompile), %.1f MiB allocated, %.3f s GC\n",
          verdict, timed.time, timed.compile_time, timed.recompile_time, timed.bytes / 1024.0^2, timed.gctime)
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
    run_csv(path) -> (header, rows, number)

A CSV file a run wrote, split into cells. `output.csv_format = auto` follows
the locale, so a test file called directly on a German desktop reads
`excel_de` (`;` delimiter, decimal comma, `.` thousands) and under
`runtests.jl` (C locale) `technical`. Comment lines (`#`, such as the
version line of the SE files) are skipped; the delimiter comes from the
header line; `number(cell)` parses a numeric cell in the format of this file.
"""
function run_csv(path::AbstractString)
  lines = filter(line -> !startswith(line, "#"), readlines(path))
  delimiter = occursin(';', first(lines)) ? ';' : ','
  cells(line) = String.(split(line, delimiter))
  # excel_de: drop the thousands dots first, then the decimal comma becomes
  # a point (one replace call applies both pairs to the original text)
  number(cell) = delimiter == ';' ? parse(Float64, replace(cell, "." => "", "," => ".")) : parse(Float64, cell)
  return (header = cells(first(lines)), rows = [cells(line) for line in lines[2:end]], number = number)
end

"""
    cgmes_report_plain(output_dir) -> String

The CGMES import report of a service run (`cgmes_import_report.md`) as the
reader of the Web UI viewer sees it: rendered by the Markdown standard
library, tags removed, one line per heading, list item, table row and
preformatted line (table cells joined by `" | "`). Tests assert on this text,
so they check the content and the rendering together, not the Markdown
spelling of the file.
"""
function cgmes_report_plain(output_dir::AbstractString)::String
  html = sprint(show, MIME"text/html"(), SparlectraApp.Markdown.parse(read(joinpath(output_dir, "cgmes_import_report.md"), String)))
  html = replace(html, r"</t[dh]>" => " | ", r"</(p|li|h\d|tr)>" => "\n")
  html = replace(html, r"<[^>]+>" => "")
  html = replace(html, r"&#(\d+);" => m -> string(Char(parse(Int, m[3:(end - 1)]))))
  html = replace(html, "&lt;" => "<", "&gt;" => ">", "&quot;" => "\"", "&amp;" => "&")
  return join(filter(!isempty, rstrip.(strip.(split(html, '\n')), '|') .|> rstrip), "\n")
end


const _FIXTURE_NET_CACHE = Dict{String,Any}()

"""
    load_fixture_net(name) -> Net

One shared fixture import per process: the named case is imported ONCE
with pinned packaged-default import options (so the fixture never depends
on whatever configuration a previous testset installed) and every call
returns a deepcopy, so testsets can mutate their copy freely.

Since load_fixture_net the first-choice names are the SHIPPED
demo cases sp_case5, sp_case6tap, sp_case14, sp_case60, sp_case188: a
fresh install loads them without any download. The legacy names case9,
case14, case57 resolve against the LOCAL data/mpower cache only and exist for the
testsets that guard externally anchored MATPOWER reference values; their
callers must gate on `fixture_net_available` and speak their skip. CGMES
deliveries are not fixtures here on purpose: their imports exercise
importer options and ARE the test subject, a shared cached net would test
the cache instead.
"""
function load_fixture_net(name::AbstractString)
  net = get!(_FIXTURE_NET_CACHE, String(name)) do
    _import_fixture_net(String(name))
  end
  return deepcopy(net)
end

"""
    large_case_dir() -> Union{Nothing,String}

The shared large-case directory from `Sparlectra.large_cases_dir()`, or
`nothing` when it does not exist. The package resolves it: the
`SPARLECTRA_LARGE_CASES_DIR` override first, otherwise the Web UI user case
directory, so a case downloaded once through the Web UI is available here too.

The suite used to resolve these out of `data/mpower` in the checkout, which is
where the package downloaded into: having used the package decided how many
assertions ran, 7180 against 7121 on a fresh tree, and nothing said so.
"""
function large_case_dir()::Union{Nothing,String}
  dir = Sparlectra.large_cases_dir()
  return isdir(dir) ? abspath(dir) : nothing
end

"""
    large_case_path(filename) -> Union{Nothing,String}

Path of `filename` under [`large_case_dir`](@ref), or `nothing` when the
directory is absent or the file is not in it. Callers gate on the result and
print a spoken SKIPPED line; nothing downloads.
"""
function large_case_path(filename::AbstractString)::Union{Nothing,String}
  dir = large_case_dir()
  dir === nothing && return nothing
  path = joinpath(dir, String(filename))
  return isfile(path) ? path : nothing
end

"""
    cgmes_fixture_dir(case) -> String

Folder of the checked-in CGMES delivery exported from the shipped case
`case`: `data/cgmes_demo/<case>` with the four profile files EQ, TP, SSH
and SV, written by `tools/gen_cgmes_fixtures.jl` under a fixed header stamp.
The CGMES tests import these folders directly and pack a zip at run time
where the service layer needs a case file; no zip is checked in.
"""
cgmes_fixture_dir(case::AbstractString)::String = joinpath(dirname(@__DIR__), "data", "cgmes_demo", String(case))

# legacy MATPOWER fixtures come from the shared large-case directory; callers
# gate on this and print a spoken SKIPPED line instead of downloading
fixture_net_available(name::AbstractString) = startswith(String(name), "sp_case") || String(name) == "warmup_casePST" || large_case_path(string(name, ".m")) !== nothing

function _import_fixture_net(name::String)
  if startswith(name, "sp_case")
    path = joinpath(dirname(Sparlectra.MPOWER_DIR), "scf", string(name, ".scf.json"))
    isfile(path) || error(string("unknown fixture net: ", name, " (shipped cases: sp_case5, sp_case6tap, sp_case14, sp_case60, sp_case188)"))
    return Sparlectra.importSCF(path)
  end
  name in ("warmup_casePST", "case9", "case14", "case57") || error(string("unknown fixture net: ", name, " (known: sp_case5, sp_case6tap, sp_case14, sp_case60, sp_case188 and warmup_casePST shipped; case9, case14, case57 from SPARLECTRA_LARGE_CASES_DIR)"))
  # warmup_casePST is TRACKED and lives in the checkout; the legacy cases are
  # not shipped and come from SPARLECTRA_LARGE_CASES_DIR
  path = name == "warmup_casePST" ? joinpath(Sparlectra.MPOWER_DIR, string(name, ".m")) : large_case_path(string(name, ".m"))
  (path === nothing || !isfile(path)) && error(string("fixture case file missing (gate on fixture_net_available): ", name, ".m; set SPARLECTRA_LARGE_CASES_DIR to a directory that holds it"))
  return Sparlectra.createNetFromMatPowerFile(filename = path, flatstart = false, enable_pq_gen_controllers = true, bus_shunt_model = :admittance, matpower_shift_sign = 1.0, matpower_shift_unit = :deg, matpower_ratio = :normal, tap_changer_model = :ideal)
end

## Scratch paths for the suite. `tempname()` names a file in /tmp that
## NOBODY ever removes: a full run left about three thousand stray
## `jl_*.yaml` files behind, and the configuration warnings they triggered
## buried the real test output. Everything the
## tests write goes into one directory instead, which Julia deletes when
## the process ends.
const TEST_SCRATCH_ROOT = mktempdir(; cleanup = true)
let counter = Ref(0)
  global function test_scratch_path(extension::AbstractString = "")
    counter[] += 1
    return joinpath(TEST_SCRATCH_ROOT, string("scratch_", counter[], extension))
  end
end

# The test grids with a line capacitance were written for the per-unit
# addACLine! before 0.50.0: r and x are line totals in per unit, the
# capacitance is given per km. This builder keeps r and x (converted to ohm
# on the bus voltage) and lets addACLinePerKm! derive b from the
# capacitance. Before 0.50.0 that siemens value was read as per unit, so the
# charging of these lines is Z_base times what it was (121 at 110 kV and
# 100 MVA).
function add_test_line_pu_c!(; net::Net, fromBus::String, toBus::String, length_km::Float64, r_pu::Float64, x_pu::Float64, c_nf_per_km::Float64, kwargs...)
  zb = get_bus_vn_kV(net = net, busName = fromBus)^2 / net.baseMVA
  return addACLinePerKm!(; net = net, fromBus = fromBus, toBus = toBus, length_km = length_km, r_ohm_per_km = r_pu * zb / length_km, x_ohm_per_km = x_pu * zb / length_km, c_nf_per_km = c_nf_per_km, kwargs...)
end

# Shared by test_state_estimation.jl, test_observability.jl and
# test_topology_validation.jl (each file runs alone with this helper file
# and testgrid.jl): the 4-bus link net of the SE link-contraction tests.
function create_se_link_net(; link_status::Int = 1, shunt_buses::Vector{String} = String[])::Net
  # S feeds B1; B1a hangs on B1 only through a bus link and carries the line
  # to B2, so a closed link transports real power (contraction pattern)
  net = Net(name = "se_link", baseMVA = 100.0)
  for b in ("S", "B1", "B1a", "B2")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addPIModelACLine!(net = net, fromBus = "S", toBus = "B1", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "B1a", toBus = "B2", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "S", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "S")
  addProsumer!(net = net, busName = "B1", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
  addProsumer!(net = net, busName = "B1a", type = "ENERGYCONSUMER", p = 5.0, q = 2.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 20.0, q = 6.0)
  for b in shunt_buses
    addShunt!(net = net, busName = b, pShunt = 0.0, qShunt = 12.0)
  end
  addLink!(net = net, fromBus = "B1", toBus = "B1a", status = link_status)
  ok, msg = validate!(net = net)
  ok || error("se_link net invalid: $msg")
  return net
end

# The topology precheck, the collapsed-branch exclusions and a frozen tap are
# the tested behavior in several of these sets, so they warn by design. They
# are captured rather than printed, and anything else that warns fails.
const SE_EXPECTED_WARNINGS = (
  r"topology precheck reported",
  r"measurement\(s\) excluded",
  r"released tap on transformer",
  # the tap fallback test forces a non-convergence on purpose: that is the
  # situation it exists for, and the fallback lines are the honest report
  # of it (run 2 warns only when it does not converge either)
  r"tap estimation: run 1 not converged",
  r"run 2 without tap estimation: not converged",
)

_se_run_quiet(testfn) = run_with_expected_warnings(testfn, SE_EXPECTED_WARNINGS)

# Shared three-bus fixture (slack, PV, PQ) used by the solver, controller,
# estimator and parallel groups across several profiles; lives here so every
# profile has it without including testgrid.jl.
function createTest3BusNet(; cooldown = 0, hyst_pu = 0.0, qlim_min = nothing, qlim_max = nothing)::Net
  # Simple 3-bus network
  #
  #  ASTADT        STATION1
  # <--|---------------|<--- Generator 
  #    |-------       |
  #            |      |
  #            --------|<---- EXTERNALNETWORKINJECTION
  #                 VERBUND  
  Sbase_MVA = 100.0
  netName = "test3bus"

  r = 0.0
  x = 0.4
  s = 25.0
  c_nf_per_km = 9.55
  tanδ = 0.0

  vm_pu_STATION1 = 1.027273
  vm_pu_VERBUND = 1.018182

  @debug "Creating $netName test network with qlim_min=$qlim_min, qlim_max=$qlim_max"

  Bus3Net = Net(name = netName, baseMVA = Sbase_MVA, cooldown_iters = cooldown, q_hyst_pu = hyst_pu)

  addBus!(net = Bus3Net, busName = "ASTADT", vn_kV = 110.0)
  addBus!(net = Bus3Net, busName = "STATION1", vn_kV = 110.0)
  addBus!(net = Bus3Net, busName = "VERBUND", vn_kV = 110.0)

  add_test_line_pu_c!(net = Bus3Net, fromBus = "ASTADT", toBus = "STATION1", length_km = s, r_pu = r, x_pu = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  add_test_line_pu_c!(net = Bus3Net, fromBus = "ASTADT", toBus = "VERBUND", length_km = s, r_pu = r, x_pu = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  add_test_line_pu_c!(net = Bus3Net, fromBus = "VERBUND", toBus = "STATION1", length_km = s, r_pu = r, x_pu = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)

  addProsumer!(net = Bus3Net, busName = "VERBUND", type = "EXTERNALNETWORKINJECTION", vm_pu = vm_pu_VERBUND, va_deg = 0.0, referencePri = "VERBUND")
  addProsumer!(net = Bus3Net, busName = "STATION1", type = "SYNCHRONOUSMACHINE", p = 70.0, q = 33.2, vm_pu = vm_pu_STATION1, qMax = qlim_max, qMin = qlim_min)
  addProsumer!(net = Bus3Net, busName = "ASTADT", type = "ENERGYCONSUMER", p = 100.0, q = 30.0)

  return Bus3Net
end

# The CIGRE HV benchmark (Task Force C6.04.02, TB 575) and the five-bus test
# net, moved here from testgrid.jl in 0.50.0 so the fast profile can solve
# them against their OpenLoadFlow fixtures without including testgrid.jl.
function createCIGRE(lLine_6a6b = 0.01)::Net
  Sbase_MVA = 1000.0
  netName = "cigre"
  net = Net(name = netName, baseMVA = Sbase_MVA)
  @debug "Creating $netName test network"
  # A R E A 1
  # Bus 1, 220kV  
  addBus!(net = net, busName = "B1", vn_kV = 220.0)
  # Bus 9, 22kV  
  addBus!(net = net, busName = "B9", vn_kV = 22.0)
  # Bus 7, 22kV  
  addBus!(net = net, busName = "B7", vn_kV = 380.0)
  # Bus 2, 220kV  
  addBus!(net = net, busName = "B2", vn_kV = 220.0)
  # Bus 10, 22kV  
  addBus!(net = net, busName = "B10", vn_kV = 22.0)
  # A R E A 2
  # B6a, 220kV  
  addBus!(net = net, busName = "B6a", vn_kV = 220.0)
  # B6b, 220kV   
  addBus!(net = net, busName = "B6b", vn_kV = 220.0)
  # Bus 12, 22kV  
  addBus!(net = net, busName = "B12", vn_kV = 22.0)
  # A R E A 3
  # Bus 5, 220kV  
  addBus!(net = net, busName = "B5", vn_kV = 220.0)
  # Bus 4, 220kV  
  addBus!(net = net, busName = "B4", vn_kV = 220.0)
  # Bus 3, 220kV
  addBus!(net = net, busName = "B3", vn_kV = 220.0)
  # Bus 8, 380kV  
  addBus!(net = net, busName = "B8", vn_kV = 380.0)
  # Bus 11, 22kV  
  addBus!(net = net, busName = "B11", vn_kV = 22.0)

  ##### Shunt    
  addShunt!(net = net, busName = "B6a", pShunt = 0.0, qShunt = 180.0)
  addShunt!(net = net, busName = "B4", pShunt = 0.0, qShunt = 160.0)
  addShunt!(net = net, busName = "B5", pShunt = 0.0, qShunt = 80.0)

  ###### ACLines   
  # 220kV-Lines
  # (bus1, bus2, length_km=100)
  add_test_line_pu_c!(net = net, fromBus = "B1", toBus = "B2", length_km = 100.0, r_pu = 0.0653, x_pu = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus1, bus6a, length_km=300)
  add_test_line_pu_c!(net = net, fromBus = "B1", toBus = "B6a", length_km = 300.0, r_pu = 0.0653, x_pu = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus2, bus5, length_km=300)
  add_test_line_pu_c!(net = net, fromBus = "B2", toBus = "B5", length_km = 300.0, r_pu = 0.0653, x_pu = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus3, bus4, length_km=100)
  add_test_line_pu_c!(net = net, fromBus = "B3", toBus = "B4", length_km = 100.0, r_pu = 0.0653, x_pu = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus3, bus4, length_km=100)
  add_test_line_pu_c!(net = net, fromBus = "B3", toBus = "B4", length_km = 300.0, r_pu = 0.0653, x_pu = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus4, bus5, length_km=300)
  add_test_line_pu_c!(net = net, fromBus = "B4", toBus = "B5", length_km = 300.0, r_pu = 0.0653, x_pu = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  # (bus4, bus6a, length_km=300)
  add_test_line_pu_c!(net = net, fromBus = "B4", toBus = "B6a", length_km = 300.0, r_pu = 0.0653, x_pu = 0.398, c_nf_per_km = 9.08, tanδ = 0.0)
  #(bus6a, bus6b, length_km=length_km_6a_6b)
  # before 0.50.0 a length up to 1 km multiplied r and x, a longer one took
  # them as totals; the same rule applies here
  scale_6a6b = lLine_6a6b <= 1.0 ? lLine_6a6b : 1.0
  add_test_line_pu_c!(net = net, fromBus = "B6a", toBus = "B6b", length_km = lLine_6a6b, r_pu = 0.0653 * scale_6a6b, x_pu = 0.398 * scale_6a6b, c_nf_per_km = 9.08, tanδ = 0.0)

  # 380kV-Lines
  # (bus7, bus8, length_km=600)
  add_test_line_pu_c!(net = net, fromBus = "B7", toBus = "B8", length_km = 600.0, r_pu = 0.0328, x_pu = 0.312, c_nf_per_km = 11.5, tanδ = 0.0)

  #### Trafos  
  # 'Trafo 1-7' (bus7, bus1, sn_mva=1000,  vn_hv_kv=380, vn_lv_kv=220, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0,   shift_degree=0.0)
  add2WTrafo!(net = net, fromBus = "B7", toBus = "B1", sn_mva = 1000.0, vk_percent = 13.0, vkr_percent = 0.28, pfe_kw = 20.0, i0_percent = 0.06)
  #Trafo 3-8, (bus8, bus3, sn_mva=1000,  vn_hv_kv=380, vn_lv_kv=220, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=0.0)
  add2WTrafo!(net = net, fromBus = "B8", toBus = "B3", sn_mva = 1000.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)
  # Trafo 9-1 (bus1, bus9, sn_mva=1000,vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0,vk_percent=13.0, pfe_kw=0, i0_percent=0,shift_degree=330.0)
  add2WTrafo!(net = net, fromBus = "B1", toBus = "B9", sn_mva = 1000.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)
  # Trafo 10-2 (bus2, bus10, sn_mva=1000, vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=330.0)
  add2WTrafo!(net = net, fromBus = "B2", toBus = "B10", sn_mva = 1000.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)
  # Trafo 11-3 (bus3, bus11, sn_mva=1000, vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=330.0)
  add2WTrafo!(net = net, fromBus = "B3", toBus = "B11", sn_mva = 1000.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)
  # 'Trafo 12-6b' (bus6b, bus12, sn_mva=500, vn_hv_kv=220, vn_lv_kv=22, vkr_percent=0.0, vk_percent=13.0, pfe_kw=0, i0_percent=0, shift_degree=330.0)
  add2WTrafo!(net = net, fromBus = "B6b", toBus = "B12", sn_mva = 500.0, vk_percent = 13.0, vkr_percent = 0.0, pfe_kw = 0.0, i0_percent = 0.0)

  # Loads
  # 'Load 2' (bus2, p_mw=285, q_mvar=200)  
  # addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 285.0, q = 200.0)
  addProsumer!(net = net, busName = "B2", type = "ENERGYCONSUMER", p = 285.0, q = 200.0)
  # 'Load 3' (bus3, p_mw=103.0, q_mvar=62.0)
  addProsumer!(net = net, busName = "B3", type = "ENERGYCONSUMER", p = 80.0, q = 30.0)
  # 'Load 4' (net_cigre_hv, bus4, p_mw=326.0, q_mvar=244.0, name='Load 4')
  addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 320.0, q = 200.0)
  # 'Load 5' (net_cigre_hv, bus5, p_mw=103, q_mvar=62, name='Load 5')	
  addProsumer!(net = net, busName = "B5", type = "ENERGYCONSUMER", p = 80.0, q = 30.0)
  # 'Load 6a' (net_cigre_hv, bus6a, p_mw=435, q_mvar=296, name='Load 6a')
  addProsumer!(net = net, busName = "B6a", type = "ENERGYCONSUMER", p = 435.0, q = 296.0)

  # Generators
  # 'Generator 9' (bus9, vm_pu=1.03, va_degree=0)
  addProsumer!(net = net, busName = "B9", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.03, va_deg = 0.0, referencePri = "B9")
  # 'Generator 10' (bus10, vm_pu=1.03, p_mw=500)
  addProsumer!(net = net, busName = "B10", type = "SYNCHRONOUSMACHINE", p = 500.0, vm_pu = 1.03, va_deg = 0.0)
  # 'Generator 11' (bus11, vm_pu=1.03, p_mw=200)
  addProsumer!(net = net, busName = "B11", type = "SYNCHRONOUSMACHINE", p = 200.0, vm_pu = 1.03, va_deg = 0.0)
  # 'Generator 12' (bus12, vm_pu=1.03, p_mw=300)
  addProsumer!(net = net, busName = "B12", type = "SYNCHRONOUSMACHINE", p = 300.0, vm_pu = 1.03, va_deg = 0.0)

  return net
end

function createTest5BusNet(; cooldown = 0, hyst_pu = 0.0, qlim_min = nothing, qlim_max = nothing, pq_only = false, mul_gens::Bool = true)::Net
  Sbase_MVA = 100.0
  netName = "test5bus"
  r = 0.05
  x = 0.5
  c_nf_per_km = 10.0
  tanδ = 0.0

  @debug "Creating $netName test network with qlim_min=$qlim_min, qlim_max=$qlim_max, pq_only=$pq_only"

  Bus5Net = Net(name = netName, baseMVA = Sbase_MVA, cooldown_iters = cooldown, q_hyst_pu = hyst_pu)
  addBus!(net = Bus5Net, busName = "B1", vn_kV = 110.0)
  addBus!(net = Bus5Net, busName = "B2", vn_kV = 110.0)
  addBus!(net = Bus5Net, busName = "B3", vn_kV = 110.0)
  addBus!(net = Bus5Net, busName = "B4", vn_kV = 110.0)
  addBus!(net = Bus5Net, busName = "B5", vn_kV = 110.0)

  add_test_line_pu_c!(net = Bus5Net, fromBus = "B1", toBus = "B3", length_km = 20.0, r_pu = r, x_pu = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  add_test_line_pu_c!(net = Bus5Net, fromBus = "B3", toBus = "B5", length_km = 50.0, r_pu = r, x_pu = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  add_test_line_pu_c!(net = Bus5Net, fromBus = "B5", toBus = "B4", length_km = 100.0, r_pu = r, x_pu = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  add_test_line_pu_c!(net = Bus5Net, fromBus = "B4", toBus = "B1", length_km = 100.0, r_pu = r, x_pu = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  add_test_line_pu_c!(net = Bus5Net, fromBus = "B4", toBus = "B2", length_km = 100.0, r_pu = r, x_pu = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)
  add_test_line_pu_c!(net = Bus5Net, fromBus = "B2", toBus = "B1", length_km = 20.0, r_pu = r, x_pu = x, c_nf_per_km = c_nf_per_km, tanδ = tanδ)

  addProsumer!(net = Bus5Net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
  if mul_gens
    addProsumer!(net = Bus5Net, busName = "B2", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
    addProsumer!(net = Bus5Net, busName = "B2", type = "ENERGYCONSUMER", p = 25.0, q = 7.0)
    vm_b3 = pq_only ? nothing : 1.03
    addProsumer!(net = Bus5Net, busName = "B3", type = "SYNCHRONOUSMACHINE", p = 20.0, q = 15.0, vm_pu = vm_b3, va_deg = 0.0, qMax = 0.6 * qlim_max, qMin = 0.6 * qlim_min)
    addProsumer!(net = Bus5Net, busName = "B3", type = "SYNCHRONOUSMACHINE", p = 10.0, q = 10.0, vm_pu = vm_b3, va_deg = 0.0, qMax = 0.4 * qlim_max, qMin = 0.4 * qlim_min)
    addProsumer!(net = Bus5Net, busName = "B4", type = "ENERGYCONSUMER", p = 25.0, q = 7.0)
    addProsumer!(net = Bus5Net, busName = "B4", type = "ENERGYCONSUMER", p = 25.0, q = 8.0)
    addProsumer!(net = Bus5Net, busName = "B5", type = "ENERGYCONSUMER", p = 10.0, q = 5.0)
    addProsumer!(net = Bus5Net, busName = "B5", type = "ENERGYCONSUMER", p = 15.0, q = 1.0)
    if !pq_only
      addProsumer!(net = Bus5Net, busName = "B5", type = "SYNCHRONOUSMACHINE", p = 0.0, vm_pu = 1.0, qMax = qlim_max, qMin = qlim_min)
    end
  else
    addProsumer!(net = Bus5Net, busName = "B2", type = "ENERGYCONSUMER", p = 50.0, q = 15.0)
    vm_b3 = pq_only ? nothing : 1.0
    addProsumer!(net = Bus5Net, busName = "B3", type = "SYNCHRONOUSMACHINE", p = 30.0, q = 25.0, vm_pu = vm_b3, va_deg = 0.0, qMax = qlim_max, qMin = qlim_min)
    addProsumer!(net = Bus5Net, busName = "B4", type = "ENERGYCONSUMER", p = 50.0, q = 15.0)
    addProsumer!(net = Bus5Net, busName = "B5", type = "ENERGYCONSUMER", p = 25.0, q = 10.0)
  end

  return Bus5Net
end

# A small network with one three-winding transformer in the star equivalent
# (add3WTPiModelTrafo!): B1 (380 kV, reference) -- line -- B2; transformer
# B2/B3/B4 (380/110/20 kV) with its star point Aux3WT_B2_B3_B4; a parallel
# two-winding transformer B2-B3, so the whole three-winding transformer out
# still supplies B3 and B5 (40 MW behind a 110 kV line); B4 (2 MW) hangs on
# the three-winding transformer alone. Shared by the N-1 and scenario tests.
function three_winding_star_net()::Net
  net = Net(name = "n1_3wt_star", baseMVA = 100.0)
  addBus!(net = net, busName = "B1", vn_kV = 380.0)
  addBus!(net = net, busName = "B2", vn_kV = 380.0)
  addBus!(net = net, busName = "B3", vn_kV = 110.0)
  addBus!(net = net, busName = "B4", vn_kV = 20.0)
  addBus!(net = net, busName = "B5", vn_kV = 110.0)
  addProsumer!(net = net, busName = "B1", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.0, va_deg = 0.0, referencePri = "B1")
  addPIModelACLine!(net = net, fromBus = "B1", toBus = "B2", r_pu = 0.003, x_pu = 0.03, b_pu = 0.0, status = 1)
  add3WTPiModelTrafo!(net = net, HBBus = "B2", MBBus = "B3", LVBus = "B4", r = (0.2, 0.3, 0.4), x = (4.0, 6.0, 10.0), b = (0.0, 0.0, 0.0), ratedU_kV = (380.0, 110.0, 20.0), ratedS_MVA = (300.0, 200.0, 50.0), status = 1)
  add2WTrafo!(net = net, fromBus = "B2", toBus = "B3", sn_mva = 300.0, vk_percent = 12.0, vkr_percent = 0.3, pfe_kw = 0.0, i0_percent = 0.0)
  addPIModelACLine!(net = net, fromBus = "B3", toBus = "B5", r_pu = 0.01, x_pu = 0.04, b_pu = 0.0, status = 1)
  addProsumer!(net = net, busName = "B5", type = "ENERGYCONSUMER", p = 40.0, q = 10.0)
  addProsumer!(net = net, busName = "B4", type = "ENERGYCONSUMER", p = 2.0, q = 0.5)
  return net
end

# C - D, A - B, A - D, B - L1 ... B - L6, and C = B by a closed switch. The
# line C - D comes first, so a CGMES import puts C first and C is the
# smallest bus of the switch group {C, B}. Six loads behind B load the group
# heavily enough that a solve starting it at 0 V ends on the low-voltage
# branch. Used by the CGMES importer tests.
function isolated_first_net()::Net
  net = Net(name = "isolated_first", baseMVA = 100.0)
  feeders = ["L$(i)" for i in 1:6]
  for b in vcat(["C", "A", "B", "D"], feeders)
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  line(f, t) = addPIModelACLine!(net = net, fromBus = f, toBus = t, r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
  line("C", "D")
  line("A", "B")
  line("A", "D")
  for l in feeders
    line("B", l)
  end
  addLink!(net = net, fromBus = "B", toBus = "C", status = 1)
  addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "A")
  addProsumer!(net = net, busName = "D", type = "ENERGYCONSUMER", p = 20.0, q = 6.0)
  for l in feeders
    addProsumer!(net = net, busName = l, type = "ENERGYCONSUMER", p = 20.0, q = 6.0)
  end
  ok, msg = validate!(net = net)
  ok || error("isolated-first net invalid: $msg")
  return net
end

# The solved isolated_first_net and two CGMES snapshots of it: s1 as
# exported, s2 with both terminals of line C - D disconnected, so C keeps
# only its switch to B and the import types it Isolated.
function isolated_first_cgmes_snapshots()
  C = Sparlectra.CGMESImporter
  net = isolated_first_net()
  redirect_stdout(() -> runpf!(net, 30, 1e-8, 0), devnull)
  dir = mktempdir()
  redirect_stdout(() -> writeCGMESFiles(net; path = dir), devnull)
  s1 = redirect_stdout(() -> C.loadCGMES(sort!(readdir(dir; join = true))), devnull)
  s2 = deepcopy(s1)
  tn_c = only(x for x in s2.byclass[:TopologicalNode] if s2.objects[x].attrs[:name] == "C")
  terms = [s2.objects[t] for t in s2.byclass[:Terminal]]
  cd = only(m for m in s2.byclass[:ACLineSegment] if any(t -> t.refs[:ConductingEquipment] == m && t.refs[:TopologicalNode] == tn_c, terms))
  for t in terms
    t.refs[:ConductingEquipment] == cd && (t.attrs[:connected] = "false")
  end
  return net, s1, s2
end

# Bus B of the solved isolated_first_net with line C - D out of service,
# solved from its own solution (C keeps its start): the high-voltage
# solution of the second snapshot.
function isolated_first_reference(net::Net)
  n2 = deepcopy(net)
  nc = n2.busDict["C"]
  for br in n2.branchVec
    (br.fromBus == nc || br.toBus == nc) && setBranchStatus!(br, false)
  end
  redirect_stdout(() -> runpf!(n2, 30, 1e-8, 0), devnull)
  return n2.nodeVec[n2.busDict["B"]]
end

# A 110 kV net with a bus C that hangs on B through a closed link alone:
# A (reference) -- B, A -- D, and with `ring = true` also B -- D; loads
# B 10 MW, C 5 MW, D 20 MW. The import marking (branches only) types C
# Isolated, the power flow solves it with B. Shared by the N-1 and scenario
# tests of the link groups.
function link_group_net(; ring::Bool)::Net
  net = Net(name = ring ? "link_group_ring" : "link_group_radial", baseMVA = 100.0)
  for b in ("A", "B", "C", "D")
    addBus!(net = net, busName = b, vn_kV = 110.0)
  end
  addPIModelACLine!(net = net, fromBus = "A", toBus = "B", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
  addPIModelACLine!(net = net, fromBus = "A", toBus = "D", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
  ring && addPIModelACLine!(net = net, fromBus = "B", toBus = "D", r_pu = 0.02, x_pu = 0.2, b_pu = 0.0, status = 1)
  # before the prosumers: addLink! wants both buses of one type
  addLink!(net = net, fromBus = "B", toBus = "C", status = 1)
  addProsumer!(net = net, busName = "A", type = "EXTERNALNETWORKINJECTION", vm_pu = 1.02, va_deg = 0.0, referencePri = "A")
  addProsumer!(net = net, busName = "B", type = "ENERGYCONSUMER", p = 10.0, q = 3.0)
  addProsumer!(net = net, busName = "C", type = "ENERGYCONSUMER", p = 5.0, q = 2.0)
  addProsumer!(net = net, busName = "D", type = "ENERGYCONSUMER", p = 20.0, q = 6.0)
  ok, msg = validate!(net = net)
  ok || error("link group net invalid: $msg")
  return net
end
