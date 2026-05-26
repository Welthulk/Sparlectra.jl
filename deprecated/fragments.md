## Finite-Difference Rectangular Jacobian (`jacobian_fd.jl`, internal/developer path)

### Purpose

The file `jacobian_fd.jl` provides a reference / compatibility implementation of
the Newton step using finite differences. This path is internal/developer-only
and is not part of the recommended user-facing configuration workflow.

It is mainly useful for:

* prototyping and validation of the analytic Jacobian,
* debugging difficult cases,
* small to medium test networks.

### Function Signature

```julia
complex_newton_step_rectangular_fd(
    Ybus,
    V,
    S;
    slack_idx::Int       = 1,
    damp::Float64        = 1.0,
    autodamp::Bool       = false,
    autodamp_min::Float64 = 1e-3,
    h::Float64           = 1e-6,
    bus_types::Vector{Symbol},
    Vset::Vector{Float64},
) -> Vector{ComplexF64}
```

### Why the FD Method Works

The FD method works because Newton's method only needs a **local linear model**
of the residual map near the current iterate. If

```math
F : \mathbb{R}^n \rightarrow \mathbb{R}^n
```

is differentiable, then for a small perturbation `δx`:

```math
F(x + \delta x) \approx F(x) + J(x)\,\delta x,
```

where `J(x)` is the Jacobian. Newton uses this first-order model and chooses
`δx` such that the linearized residual becomes zero:

```math
J(x)\,\delta x = -F(x).
```

If we do not assemble `J(x)` analytically, we can approximate each column by
perturbing one state variable at a time:

```math
\frac{\partial F}{\partial x_k}(x)
\approx
\frac{F(x + h e_k) - F(x)}{h},
```

with `e_k` the k-th unit vector and `h` a small step size. This is exactly what
the FD implementation does. The approximation is first-order accurate in `h`,
so when `h` is small enough, the FD Jacobian is a good surrogate for the true
Jacobian.

In short:

* the **physics** stays the same because `F` is still the same mismatch
  function,
* only the **derivative evaluation** changes,
* and Newton still solves the same linearized correction equation.

### Internal Workflow

1. **Base mismatch**

   ```julia
   F0 = mismatch_rectangular(Ybus, V, S, bus_types, Vset, slack_idx)
   m  = length(F0)      # expected = 2*(n-1)
   ```

2. **Non-slack variables**

   Only non-slack buses are state variables:

   ```julia
   non_slack = collect(1:n)
   deleteat!(non_slack, slack_idx)
   nvar = 2 * (n - 1)   # Vr(non-slack) and Vi(non-slack)
   ```

3. **Jacobian by finite differences**

   * First, perturb each real part `V_r`:

     ```julia
     for (col_idx, bus) in enumerate(non_slack)
         V_pert = copy(V)
         V_pert[bus] = ComplexF64(real(V[bus]) + h, imag(V[bus]))
         Fp = mismatch_rectangular(...)
         J[:, col_idx] .= (Fp .- F0) ./ h
     end
     ```

   * Then, perturb each imaginary part `V_i`:

     ```julia
     for (offset, bus) in enumerate(non_slack)
         col_idx = (n - 1) + offset
         V_pert = copy(V)
         V_pert[bus] = ComplexF64(real(V[bus]), imag(V[bus]) + h)
         Fp = mismatch_rectangular(...)
         J[:, col_idx] .= (Fp .- F0) ./ h
     end
     ```

4. **Solve linear system**

   ```julia
   δx = J \ (-F0)
   ```

   with fallback to a pseudo-inverse path in the solver utilities if needed.
   Then apply damping:

   ```julia
   δx .*= damp
   ```

5. **Return updated voltages**

   ```julia
   V_new = ComplexF64.(Vr_new, Vi_new)
   ```

`run_complex_nr_rectangular` can select this path in internal validation flows.

### Practical Interpretation

The FD solver is therefore not a different load-flow model. It is the same NR
equation system with a numerically approximated Jacobian. This is why it is
useful as a validation path for the analytic Jacobian: both methods should
follow the same local Newton geometry, up to FD truncation and step-size
effects.

---
