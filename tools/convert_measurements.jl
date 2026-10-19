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
#
# file: tools/convert_measurements.jl
# purpose: command line for convertMeasurementsV1: write the measurement
#          table of a measurement CSV v1 file and its case, keyed the way
#          the SE service reads it (mRIDs for a CGMES case, names otherwise).
#          Usage: julia --project=app tools/convert_measurements.jl <case> <v1.csv> [<out.csv>] [--config <file>]
#          Without <out.csv> the table is written next to the v1 file as
#          <name>.table.csv. Prints the converter report; exits 1 when a
#          row could not be converted.

using SparlectraApp
using Sparlectra

function main(args::Vector{String})::Int
  config = Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH
  pos = String[]
  i = 1
  while i <= length(args)
    if args[i] == "--config" && i < length(args)
      config = args[i+1]
      i += 2
    else
      push!(pos, args[i])
      i += 1
    end
  end
  if !(length(pos) in (2, 3))
    println(stderr, "usage: julia --project=app tools/convert_measurements.jl <case> <v1.csv> [<out.csv>] [--config <file>]")
    return 2
  end
  case, v1 = pos[1], pos[2]
  out = length(pos) == 3 ? pos[3] : string(chopsuffix(v1, ".csv"), ".table.csv")
  res = SparlectraApp._se_convert_measurements_v1(case, v1, out; config_file = config)
  print(res.report)
  return isempty(res.problems) ? 0 : 1
end

exit(Base.invokelatest(main, ARGS))
