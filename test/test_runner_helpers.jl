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

function quiet_test_output(f::Function; verbose::Bool = sparlectra_test_verbose())
  verbose && return Base.invokelatest(f)
  path = tempname()
  result = nothing
  captured = ""
  try
    open(path, "w+") do io
      result = redirect_stdio(stdout = io) do
        Base.invokelatest(f)
      end
      seekstart(io)
      captured = read(io, String)
    end
  catch
    isfile(path) && (captured = read(path, String))
    !isempty(captured) && print(captured)
    rethrow()
  finally
    ispath(path) && rm(path; force = true)
  end
  if occursin("Test Failed", captured) || occursin("Error During Test", captured)
    print(captured)
  end
  return result
end
