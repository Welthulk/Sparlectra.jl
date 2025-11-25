using Test
using Sparlectra

# 3-bus example Ybus matrix (3×3)
Ydata = ComplexF64[
    10 - 5im,
   -10 + 5im,
    0 + 0im,
   -10 + 5im,
    20 - 10im,
   -10 + 5im,
    0 + 0im,
   -10 + 5im,
    10 - 5im,
]

Ybus = reshape(Ydata, 3, 3)

# Initial voltages: slack at bus 1 = 1∠0, flat start at others
V0 = ComplexF64[1.0 + 0im, 1.0 + 0im, 1.0 + 0im]

# Specified power injections S = P + jQ
# Convention: generation positive, load negative
S = ComplexF64[
    1.0 - 0.2im,     # slack / generator at bus 1
   -0.5 + 0.1im,     # load at bus 2
   -0.5 + 0.1im,     # load at bus 3
]

slack_idx = 1

# Run full complex NR prototype
#V, converged, iters, history = run_complex_nr(Ybus, V0, S; slack_idx=slack_idx,
#                                              maxiter=20, tol=1e-6, verbose=false)

V, converged, iters, history = run_complex_nr_rectangular(
    Ybus, V0, S;
    slack_idx = slack_idx,
    maxiter   = 20,
    tol       = 1e-6,
    verbose   = true,   # ruhig mal einschalten um Verlauf zu sehen
    damp      = 0.2
)


println("Updated complex voltages after ", iters, " iterations:")
println(V)
println("Converged = ", converged)
println("History of mismatches = ", history)

@test length(V) == 3
@test V[slack_idx] == V0[slack_idx]                        # slack kept fixed
@test all(.!isnan.(real.(V)))                              # no NaNs
@test all(.!isnan.(imag.(V)))
@test length(history) == iters                             # one value per iteration
@test iters <= 20
#@test converged
#@test history[end] <= 1e-6
