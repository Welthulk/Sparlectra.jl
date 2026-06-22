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

using Sparlectra

function main(args = ARGS)
  cases = isempty(args) ? collect(Sparlectra.QLIMIT_LARGE_CASE_DEFAULTS) : collect(args)
  result = Sparlectra.compare_qlimit_large_case_modes(cases = cases)
  println("Q-limit large-case mode comparison complete")
  println("rows: ", length(result.rows))
  println("csv : ", result.csv_path)
  println("json: ", result.json_path)
  for row in result.rows
    println(row["case"], " | ", row["mode"], " | ", row["run_status"], " | ", row["numerical_status"], " | ", row["run_directory"])
  end
  return result
end

if get(ENV, "SPARLECTRA_QLIMIT_LARGE_CASE_NO_MAIN", "0") != "1"
  Base.invokelatest(main)
end
