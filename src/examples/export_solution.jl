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
Example: Run Sparlectra's internal solver and export a solver-agnostic model + solution.

What this does:
1) Load a MATPOWER (or Julia) case
2) Build `Net`
3) Run Sparlectra internal power-flow solver
4) Build PFModel from the solved net
5) Export PFSolution (PF-ordered voltages) for external solver comparison
"""

using Sparlectra

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

case = "case14.m"      # or "case14.jl"
flatstart = false
opt_sparse = true
opt_fd = true          # leave as you like
method = :polar_full    # internal solver method
max_iter = 25
show_net_verbose = false
show_results = false
filename = joinpath(@__DIR__, "..", "..", "data", "mpower", case)

# -----------------------------------------------------------------------------
# Step 1: Create Net
# -----------------------------------------------------------------------------

println("Creating Net from case file: ", filename)

net = createNetFromMatPowerFile(
  filename = filename,
  log = true,
  flatstart = flatstart,
)

# -----------------------------------------------------------------------------
# Step 2: Run internal solver
# -----------------------------------------------------------------------------

println("\nRunning internal Sparlectra solver ...")
println("  method     = ", method)
println("  opt_sparse = ", opt_sparse)
println("  opt_fd     = ", opt_fd)
println("  flatstart  = ", flatstart)

# run_acpflow returns the solved Net (as used in your examples)
solved_net = run_acpflow(
  casefile = case,
  opt_sparse = opt_sparse,
  opt_fd = opt_fd,
  method = method,
  opt_flatstart = flatstart,  
  show_results = show_results,
)

# If you prefer directly:
# status, iters = runpf!(net, max_iter)
if show_net_verbose
    showNet(solved_net, verbose = true)
end

# -----------------------------------------------------------------------------
# Step 3: Export PFModel from solved net
# -----------------------------------------------------------------------------

println("\nExporting PFModel from solved Net ...")

model = buildPfModel(
  solved_net;
  opt_sparse = opt_sparse,
  flatstart = false,          # solved state is already in Net; do NOT flatstart again
  include_limits = true,      # optional: include PV Q-limits for external experiments
  verbose = (show_net_verbose ? 1 : 0),
)

println("PFModel:")
println("  n(active)         = ", length(model.busIdx_net))
println("  slack_pf          = ", model.slack_idx)
println("  initial mismatch  = ", mismatchInf(model, model.V0))

# -----------------------------------------------------------------------------
# Step 4: Export PFSolution from solved net
# -----------------------------------------------------------------------------

# Build the PF-ordered voltage vector from solved_net using busIdx_net mapping.
# Note: model.busIdx_net stores ORIGINAL bus indices (Net ordering),
# and the PF ordering matches Ybus / BusData ordering.
V_pf = Vector{ComplexF64}(undef, length(model.busIdx_net))
for (k, busIdx) in enumerate(model.busIdx_net)
  node = solved_net.nodeVec[busIdx]
  # node stores va in degrees, vm in p.u.
  Vk = ComplexF64(node._vm_pu * cosd(node._va_deg), node._vm_pu * sind(node._va_deg))
  V_pf[k] = Vk
end

res = mismatchInf(model, V_pf)

sol = PFSolution(
  V = V_pf,
  converged = isfinite(res) && (res <= 1e-8),
  iters = -1,                  # not available from run_acpflow; set if you have it
  residual_inf = res,
  meta = (method = method, opt_sparse = opt_sparse, opt_fd = opt_fd, flatstart = flatstart),
)

println("\nExported solution:")
println("  residual ||F(V)||∞ = ", sol.residual_inf)
println("  converged          = ", sol.converged)

# -----------------------------------------------------------------------------
# Step 5: Hand-off point for external solver comparison
# -----------------------------------------------------------------------------

println("\nHand-off:")
println("  model :: PFModel   (Ybus/specs/start/types/limits)")
println("  sol   :: PFSolution (PF-ordered V)")

# Example: a downstream tool can now compare:
# - NR reference: sol.V
# - External solver result: V_ext
# and evaluate mismatchInf(model, V_ext)

println("\nExample finished.")
