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

# file: tools/sysimage_workload.jl
# purpose: precompile execution workload for the fast-start sysimage build
#          (tools/build_sysimage.jl). Exercises the interactive paths beyond
#          the shipped PrecompileTools workload: Web UI start with the
#          reserved warm-up machinery (no run-history entries), one page
#          request through the real socket handler, one run_sparlectra call
#          on case14, then a clean shutdown. Never run this file directly;
#          PackageCompiler executes it in a child process during the build.

using Sparlectra
using Sockets

# nonstandard port so a Web UI the maintainer has open on 8080 does not make
# the build silently skip the server paths
const _WORKLOAD_PORT = 8091

function _workload_request(port::Int, target::String)
  sock = Sockets.connect("127.0.0.1", port)
  try
    write(sock, string("GET ", target, " HTTP/1.1\r\nHost: 127.0.0.1:", port, "\r\nConnection: close\r\n\r\n"))
    return String(read(sock))
  finally
    close(sock)
  end
end

function run_workload()
  server = start_sparlectra_webui(open_browser = false, port = _WORKLOAD_PORT, warmup = true, warmup_store_result = false)
  server === nothing && error("sysimage workload: port $(_WORKLOAD_PORT) is already serving a Sparlectra Web UI")
  try
    # the warm-up waits for the first served page; request it through the
    # real socket path so the handler chain is part of the compile trace
    _workload_request(_WORKLOAD_PORT, "/")
    deadline = time() + 180.0
    while server.runtime.warmup_state !== :done && time() < deadline
      sleep(0.2)
    end
    server.runtime.warmup_state === :done || @warn "sysimage workload: warm-up did not finish within 180 s; continuing with the paths compiled so far"
    # after the warm-up the same route serves the real PowerFlow form; that
    # render path must be part of the trace, otherwise the first form view
    # after a fast start pays it as JIT time (measured ~11 s)
    _workload_request(_WORKLOAD_PORT, "/powerflow")
    _workload_request(_WORKLOAD_PORT, "/webui/fast-start")
    # one full interactive solve (the case is pre-fetched by the build
    # script); logfile output stays off so the workload never leaves
    # run_case*.log files in examples/_out, regardless of any user
    # configuration
    workload_cfg = SparlectraConfig(output = OutputConfig(logfile_results = :off))
    result = run_sparlectra(casefile = "case14.m", config = workload_cfg)
    result.final_converged || @warn "sysimage workload: case14 run did not converge" outcome = result.outcome
  finally
    close(server)
  end
  return nothing
end

Base.invokelatest(run_workload)
