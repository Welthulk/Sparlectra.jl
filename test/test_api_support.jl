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

function _write_api_test_case(path::AbstractString)
  write(path, """
function mpc = case_api
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
1 3 0 0 0 0 1 1.0 0 110 1 1.1 0.9;
2 1 40 15 0 0 1 1.0 0 110 1 1.1 0.9;
];
mpc.gen = [
1 100 0 300 -300 1.02 100 1 300 0;
];
mpc.branch = [
1 2 0.01 0.05 0.0 999 999 999 0 0 1 -360 360;
];
""")
  return path
end

function _write_api_dcline_case(path::AbstractString; rows::AbstractString)
  _write_api_test_case(path)
  open(path, "a") do io
    write(io, "\nmpc.dcline = [\n", rows, "\n];\n")
  end
  return path
end

function _write_api_two_island_case(path::AbstractString)
  write(path, """
function mpc = case_two_islands
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
1 3 0 0 0 0 1 1.0 0 110 1 1.1 0.9;
2 1 30 10 0 0 1 1.0 0 110 1 1.1 0.9;
3 3 0 0 0 0 1 1.0 0 110 1 1.1 0.9;
4 1 20 5 0 0 1 1.0 0 110 1 1.1 0.9;
];
mpc.gen = [
1 50 0 300 -300 1.02 100 1 300 0;
3 30 0 300 -300 1.01 100 1 300 0;
];
mpc.branch = [
1 2 0.01 0.05 0.0 999 999 999 0 0 1 -360 360;
3 4 0.01 0.05 0.0 999 999 999 0 0 1 -360 360;
];
""")
  return path
end
