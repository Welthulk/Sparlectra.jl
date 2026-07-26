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

# Date: 2026-07-26
# file: src/precompile.jl
# purpose: PrecompileTools workload covering the MATPOWER import -> rectangular/DC power-flow path so first interactive runs skip most JIT compilation (issue #288, finding 5)

# The workload runs once at package precompile time and caches the compiled
# native code for the hot path: MATPOWER parsing, Net construction, the
# rectangular Newton-Raphson solve (UMFPACK backend, autodamp + merit line
# search as used by the default configuration), loss postprocessing, and the
# standalone DC power flow. Interactive first-run latency on large cases is
# dominated by exactly this JIT work (measured ~60 s perceived vs ~12 s warm
# on case_SyntheticUSA before this workload existed).
#
# case14.m is a tracked repository file that ships with the package; the
# isfile guard keeps precompilation robust if the data directory is stripped.

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
  _precompile_case = normpath(joinpath(@__DIR__, "..", "data", "mpower", "case14.m"))
  @compile_workload begin
    if isfile(_precompile_case)
      Logging.with_logger(Logging.NullLogger()) do
        mpc = MatpowerIO.read_case_m(_precompile_case; legacy_compat = false)
        net = createNetFromMatPowerCase(mpc = mpc)
        runpf!(net, 30, 1.0e-6, 0; method = :rectangular, autodamp = true, merit_enabled = true)
        calcNetLosses!(net)
        rundcpf!(net)
      end
    end
  end
end
