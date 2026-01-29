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

"""
export_solution_for_external_solver.jl

Run Sparlectra's internal solver and export a solver-agnostic PFModel + PFSolution
for external solver experiments. Optionally also run an external-solver path
(via `runpf_external!`) to test the interface and compare against the internal
reference.

What you can toggle:
- internal solver run (NR/polar/rectangular etc. via run_acpflow)
- external solver path run (dummy solver stub; replace with your own solver)
- printing: internal solver result printing, Net verbose printing, interface show flags
"""

using Sparlectra

# =============================================================================
# Configuration
# =============================================================================

case = "case14.m"          # or "case14.jl"
flatstart = false
opt_sparse = true
opt_fd = true
method = :polar_full        # internal solver method (run_acpflow)
max_iter = 25               # only used if you switch to runpf!(...) manually
tol = 1e-8                  # convergence tolerance for mismatchInf

# --- Toggles ---
run_internal = true         # run internal solver and export reference (recommended)
run_external = true         # run external path (dummy solver) to test interface
show_results_internal = false   # controls run_acpflow results header/table (you already wired this)
show_net_verbose = false        # showNet(...) after internal run
show_model = true              # interface show option (runpf_external!)
show_solution = true           # interface show option (runpf_external!)
show_export_summary = true     # prints export summary from this script

# Input file
filename = joinpath(@__DIR__, "..", "..", "data", "mpower", case)

# =============================================================================
# Helpers
# =============================================================================

"Build PF-ordered complex voltage vector from a solved Net using PFModel mapping."
function build_V_pf_from_net(net::Net, model::PFModel)::Vector{ComplexF64}
  V_pf = Vector{ComplexF64}(undef, length(model.busIdx_net))
  @inbounds for (k, busIdx) in enumerate(model.busIdx_net)
    node = net.nodeVec[busIdx]
    V_pf[k] = ComplexF64(node._vm_pu * cosd(node._va_deg), node._vm_pu * sind(node._va_deg))
  end
  return V_pf
end

# =============================================================================
# External solver stub (replace later)
# =============================================================================

struct DummyExternalSolver <: Sparlectra.AbstractExternalSolver end

function Sparlectra.solvePf(::DummyExternalSolver, model::Sparlectra.PFModel; tol=1e-8, kwargs...)
    # pretend "external" solver magically found the exact solution by running internal NR
    net_tmp = createNetFromMatPowerFile(filename=filename, log=false, flatstart=false)
    runpf!(net_tmp, 25)
    model_tmp = buildPfModel(net_tmp; opt_sparse=true, flatstart=false, include_limits=false, verbose=0)
    V = build_V_pf_from_net(net_tmp, model_tmp)
    r = mismatchInf(model_tmp, V)
    return PFSolution(V=V, converged=(r<=tol), iters=0, residual_inf=r)
end

# =============================================================================
# Main
# =============================================================================

println("Case file: ", filename)
println("Options:")
println("  run_internal           = ", run_internal)
println("  run_external           = ", run_external)
println("  method (internal)      = ", method)
println("  opt_sparse / opt_fd    = ", opt_sparse, " / ", opt_fd)
println("  flatstart              = ", flatstart)
println("  show_results_internal  = ", show_results_internal)
println("  show_model/solution    = ", show_model, " / ", show_solution)
println("")

# -----------------------------------------------------------------------------
# 1) Internal solver run + export (reference)
# -----------------------------------------------------------------------------

model_ref = nothing
sol_ref = nothing
solved_net_ref = nothing

if run_internal
  println("=== Internal solver run (reference) ===")

  # Run internal solver via run_acpflow (your standard entry point)
  solved_net_ref = run_acpflow(
    casefile = case,
    opt_sparse = opt_sparse,
    opt_fd = opt_fd,
    method = method,
    opt_flatstart = flatstart,
    show_results = show_results_internal,
  )

  if show_net_verbose
    showNet(solved_net_ref, verbose = true)
  end

  # Export PFModel from solved net (do NOT re-flatstart)
  model_ref = buildPfModel(
    solved_net_ref;
    opt_sparse = opt_sparse,
    flatstart = false,
    include_limits = true,
    verbose = (show_net_verbose ? 1 : 0),
  )

  # Export PFSolution from solved net in PF ordering
  V_ref = build_V_pf_from_net(solved_net_ref, model_ref)
  res_ref = mismatchInf(model_ref, V_ref)

  sol_ref = PFSolution(
    V = V_ref,
    converged = isfinite(res_ref) && (res_ref <= tol),
    iters = -1,  # run_acpflow prints iterations, but doesn't return them; keep -1 here
    residual_inf = res_ref,
    meta = (solver = :internal, method = method, opt_sparse = opt_sparse, opt_fd = opt_fd, flatstart = flatstart),
  )

  if show_export_summary
    println("\nExported INTERNAL reference:")
    println("  n(active)           = ", length(model_ref.busIdx_net))
    println("  slack_pf            = ", model_ref.slack_idx)
    println("  residual ||F(V)||∞  = ", sol_ref.residual_inf)
    println("  converged (by tol)  = ", sol_ref.converged)
  end

  println("")
end

# -----------------------------------------------------------------------------
# 2) External solver path run (interface test)
# -----------------------------------------------------------------------------

model_ext = nothing
sol_ext = nothing
net_ext = nothing

if run_external
  println("=== External solver path run (interface test) ===")

  # Create a fresh net for the external run (so we don't overwrite the internal reference net)
  net_ext = createNetFromMatPowerFile(
    filename = filename,
    log = false,
    flatstart = flatstart,
  )

  solver = DummyExternalSolver()

  iters_ext, status_ext, sol_ext = runpf_external!(
    net_ext,
    solver;
    tol = tol,
    opt_sparse = opt_sparse,
    flatstart = flatstart,
    include_limits = false,
    verbose = 0,
    show_model = show_model,
    show_solution = show_solution,
  )

  println("External path finished:")
  println("  status              = ", status_ext == 0 ? "converged" : "not converged")
  println("  iters/order         = ", iters_ext)
  println("  residual ||F(V)||∞  = ", sol_ext.residual_inf)
  println("")
end

# -----------------------------------------------------------------------------
# 3) Compare internal reference vs external result (if both enabled)
# -----------------------------------------------------------------------------

if run_internal && run_external
  println("=== Compare (internal reference vs external) ===")

  # Build a PFModel consistent with the external run (same net config).
  # We could reuse the model built inside runpf_external! by re-building here:
  model_cmp = buildPfModel(
    net_ext;
    opt_sparse = opt_sparse,
    flatstart = flatstart,
    include_limits = false,
    verbose = 0,
  )

  # NOTE: internal reference may include limits; comparison is still fine on voltages.
  # To compare in the same PF ordering, use model_ref ordering (it is canonical via getBusData sorting).
  # The safest is to compare on a common ordering: we assume both use the same bus ordering for non-iso cases.
  if length(model_ref.busIdx_net) == length(model_cmp.busIdx_net) && model_ref.busIdx_net == model_cmp.busIdx_net
    dV = sol_ext.V .- sol_ref.V
    max_dVm = maximum(abs.(abs.(sol_ext.V) .- abs.(sol_ref.V)))
    max_dVa = maximum(abs.(rad2deg.(angle.(sol_ext.V)) .- rad2deg.(angle.(sol_ref.V))))

    println("Max |ΔVm| (p.u.)  = ", max_dVm)
    println("Max |ΔVa| (deg)   = ", max_dVa)
    println("Ref residual      = ", sol_ref.residual_inf)
    println("Ext residual      = ", sol_ext.residual_inf)
  else
    println("Comparison skipped: PF ordering differs between internal/export and external run.")
    println("  internal busIdx_net = ", model_ref.busIdx_net)
    println("  external  busIdx_net = ", model_cmp.busIdx_net)
  end

  println("")
end

println("Done.")
