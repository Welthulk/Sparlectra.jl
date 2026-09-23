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
# file: src/session.jl
# purpose: process-level session facts the solvers ask for: the cooperative
#          abort hook (a running job arms a token, every Newton and WLS
#          iteration asks it) and the kind of Julia session (native, custom
#          sysimage, packaged app), which decides whether a first-run hint
#          is worth printing.

struct PowerFlowAborted <: Exception end
Base.showerror(io::IO, ::PowerFlowAborted) = print(io, "PowerFlow run aborted by user.")

## Task-local abort hook (the abort must stop the running SE or PF). A
## Julia task cannot be killed from outside, so the
## solver loops have to ASK. The worker deposits its token in its own task
## storage, and every iteration of the Newton and the WLS loop calls
## `sparlectra_abort_requested()`, which is a dictionary lookup and costs
## nothing next to one factorization. Task storage rather than a global:
## it is isolated per run by construction, and a solver called directly
## from the library sees no token and never checks.
const _SPARLECTRA_ABORT_KEY = :sparlectra_abort_token

## The armed token of the run currently executing. It is deliberately NOT
## task-local: `Threads.@spawn` gives every worker a fresh, empty task-local
## storage, so a token stored there is invisible inside the parallel island
## solve, the parallel N-1 batch and the short-circuit sweep - exactly the
## long runs a user wants to abort. The Web UI executes one run at a time
## (_POWERFLOW_WEBUI_BLOCKING_STATES), so one process-wide slot is the honest
## model; the task-local copy is still written so an existing nested arming
## keeps working.
const _SPARLECTRA_ABORT_SLOT = Ref{Any}(nothing)

function sparlectra_arm_abort_token!(token)
  token === nothing && return nothing
  task_local_storage(_SPARLECTRA_ABORT_KEY, token)
  _SPARLECTRA_ABORT_SLOT[] = token
  return nothing
end

## Clear the process-wide slot when a run finishes, so a later solver call
## outside any job never sees a stale token.
function sparlectra_disarm_abort_token!()
  _SPARLECTRA_ABORT_SLOT[] = nothing
  return nothing
end

"""
    sparlectra_abort_requested() -> Bool

True when the surrounding Web UI job has been asked to abort. Solver loops
call this once per iteration and stop with `PowerFlowAborted`. Works inside
spawned tasks as well, which task-local storage does not.
"""
function sparlectra_abort_requested()::Bool
  token = get(task_local_storage(), _SPARLECTRA_ABORT_KEY, nothing)
  token === nothing && (token = _SPARLECTRA_ABORT_SLOT[])
  token === nothing && return false
  return token[]
end

## Throwing form for the loops.
function sparlectra_check_abort()
  sparlectra_abort_requested() && throw(PowerFlowAborted())
  return nothing
end

"""
    runtime_session_kind() -> Symbol

`:app` inside a packaged application (`SPARLECTRA_APP_BUILT` is set),
`:sysimage` when the process runs on a custom system image, `:native` on
the stock Julia image. The application layer refines `:sysimage` with the
build time of its image; the core only needs to know whether a first run
will pay for compilation.
"""
function runtime_session_kind()::Symbol
  isempty(strip(get(ENV, "SPARLECTRA_APP_BUILT", ""))) || return :app
  img = try
    unsafe_string(Base.JLOptions().image_file)
  catch
    ""
  end
  isempty(img) && return :native
  # the stock image is sys.<ext>; a PackageCompiler image carries another name
  return startswith(basename(img), "sys.") ? :native : :sysimage
end
