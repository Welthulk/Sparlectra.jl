# Copyright 2023–2025 Udo Schmitz
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

case = "case14.m"
case = "case14.jl"
flatstart = false
filename = joinpath(@__DIR__, "..", "..", "data", "mpower", case)
log = true
#net = createNetFromMatPowerFile(filename = filename, log = false, flatstart = true)

#showNet(net, verbose = true) 

#run_acpflow(casefile = case, opt_sparse = true, opt_fd = false, method = :classic, opt_flatstart = flatstart)

#status, iters = runpf!(net, 25)

#println("status = ", status, "   iters = ", iters)
#println("converged = ", status == 1)

#mynet = run_acpflow(casefile = case, opt_sparse = true, method = :polar_full, opt_flatstart = flatstart)
#printQLimitLog(mynet; sort_by = :bus)
#run_acpflow(casefile = case, opt_fd = true, method = :polar_full, opt_flatstart = flatstart)
run_acpflow(casefile = case, opt_fd = true, opt_sparse = true, method = :polar_full, opt_flatstart = flatstart)
#run_acpflow(casefile = case, opt_sparse = true, opt_fd = false, method = :rectangular, opt_flatstart = flatstart)
#run_acpflow(casefile = case, opt_sparse = true, opt_fd = false, method = :classic, opt_flatstart = flatstart)

