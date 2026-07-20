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

# Date: 2026-07-20
# file: examples/internal/example_header.jl
# purpose: shared console banner and entry-point runner for example programs

"""
    print_example_banner(file::AbstractString, purpose::AbstractString)

Prints a page-break-style banner (rule / file / purpose) to the console. Every
top-level `examples/*.jl` program calls this once, at the start of its entry
function, so its output is clearly delimited — matching the file's `# file:`/
`# purpose:` header comment.
"""
function print_example_banner(file::AbstractString, purpose::AbstractString)
  rule = "="^80
  println(rule)
  println("file:    ", file)
  println("purpose: ", purpose)
  println(rule)
  println()
  return nothing
end

"""
    run_example(f, args...; kwargs...)

Runs an example program's entry function `f` (typically `main`) via
`Base.invokelatest`, and returns its result.

Every top-level `examples/*.jl` program calls this unconditionally at file
scope — not gated behind `if abspath(PROGRAM_FILE) == @__FILE__`. That guard
is only true for `julia examples/....jl` invoked from a shell; it stays false
when the file is run interactively (VS Code "Run File", `include(...)` in a
REPL), so a guarded call would silently do nothing there, even though these
programs are primarily written for interactive use. `Base.invokelatest` keeps
the call safe under Julia 1.12 world-age rules regardless of how the file was
loaded (a fresh `include`, a re-`include` under Revise, or a plain script
run).
"""
run_example(f, args...; kwargs...) = Base.invokelatest(f, args...; kwargs...)
