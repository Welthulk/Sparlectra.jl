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

# file: ext/SparlectraBenchmarkToolsExt.jl
# purpose: the repeated timing of the benchmark mode of run_matpower_case,
#          loaded by Julia once BenchmarkTools is loaded next to Sparlectra.
#          BenchmarkTools is a weak dependency: an installation of
#          Sparlectra does not compile it, nor JSON, Parsers and StructUtils
#          that come with it.

module SparlectraBenchmarkToolsExt

using Sparlectra
using BenchmarkTools

# `f` is the whole solve of one method; the trial is sampled up to
# `samples` times within `seconds`, and median and minimum come back in
# seconds (BenchmarkTools reports nanoseconds)
function Sparlectra.benchmark_trial_seconds(f::Function; samples::Int, seconds::Float64)
  bench = BenchmarkTools.@benchmarkable $f()
  trial = BenchmarkTools.run(bench; samples = samples, seconds = seconds)
  return (BenchmarkTools.median(trial).time / 1e9, BenchmarkTools.minimum(trial).time / 1e9)
end

end # module SparlectraBenchmarkToolsExt
