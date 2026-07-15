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
  text = open(path, "r") do io
    read(io, String; maxbytes = min(Int(total_bytes), max_bytes))
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

function quiet_test_output(f::Function; verbose::Bool = sparlectra_test_verbose(), group::AbstractString = "test group")
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
  timed = try
    @timed quiet_test_output(runner; group = name)
  catch err
    err isa InterruptException && (println("INTERRUPTED"); rethrow())
    println("FAIL")
    rethrow()
  end
  @printf("PASS %.3f s, %.1f MiB allocated, %.3f s GC\n", timed.time, timed.bytes / 1024.0^2, timed.gctime)
  get(ENV, "SPARLECTRA_TEST_GC_BETWEEN_GROUPS", "0") in ("1", "true", "yes", "on") && GC.gc()
  return timed.value
end
