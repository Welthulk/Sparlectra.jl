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
import Sparlectra: MatpowerIO


outdir = joinpath(@__DIR__, "..", "..", "data", "mpower")
url_case = "https://raw.githubusercontent.com/MATPOWER/matpower/master/data/case14.m"

Sparlectra.FetchMatpowerCase.ensure_matpower_case(
    url=url_case,
    outdir=outdir,
    to_jl=true,
    overwrite=false,
)
case="case14.jl"
mpc = MatpowerIO.read_case(joinpath(outdir, case))

flatstart = false
filename = joinpath(@__DIR__, "..", "..", "data", "mpower", case)
run_acpflow(casefile = case, opt_fd = true, opt_sparse = true, method = :polar_full, opt_flatstart = flatstart)

