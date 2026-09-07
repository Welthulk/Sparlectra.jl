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
# file: start_webui.jl
# purpose: start the Sparlectra Web UI, and own the sysimage decision for
#          every platform. The start scripts only call this file, so Linux,
#          macOS and Windows go through exactly the same code: check for a
#          usable sysimage (tools/sysimage_launcher.jl), offer to build one
#          when it is missing or outdated, then relaunch through the image.

# --- sysimage handling ------------------------------------------------------
# This runs BEFORE `using Sparlectra`, deliberately: the launcher module sits
# on plain Base, so deciding whether to relaunch through the sysimage costs
# no package load. Julia can only take an image at process start (-J), so a
# process that wants one has to start itself again.
include(joinpath(@__DIR__, "tools", "sysimage_launcher.jl"))
using .SysimageLauncher: handle_sysimage

handle_sysimage(copy(ARGS), @__FILE__, abspath(@__DIR__))

using Sparlectra

function main()
  # No warm-up: either this process runs on the sysimage, where the code is
  # already compiled, or the user chose to compile on first use and was told
  # so. A hidden warm-up run on top of that only delayed the first page.
  server = Sparlectra.start_sparlectra_webui(open_browser = true)
  # nothing = a Sparlectra Web UI already runs on the port; it was opened in
  # the browser instead, so there is no new server task to wait on.
  server === nothing && return nothing

  try
    wait(server.task)
  catch err
    if err isa InterruptException
      @info "Ctrl+C received; closing Sparlectra Web UI."
      close(server)
      return nothing
    end
    rethrow()
  end
end

# VS Code's inline eval runs the whole file as one world age. Since Julia
# 1.12.7 even a direct getfield counts as a binding access from the old
# world, so the lookup itself must also go through invokelatest.
Base.invokelatest(Base.invokelatest(getfield, @__MODULE__, :main))
