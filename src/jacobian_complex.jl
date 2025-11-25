# jacobian_complex.jl
# Complex Jacobian for complex-state Newton–Raphson using Wirtinger calculus

"""
    build_complex_jacobian(Ybus, V)

Builds the 2n × 2n Wirtinger-type Jacobian blocks for the complex-state
Newton–Raphson formulation.

Given:
    I = Ybus * V
    S = V .* conj.(I)

We construct the blocks:
    J11 = ∂S/∂V
    J12 = ∂S/∂V*
    J21 = ∂conj(S)/∂V
    J22 = ∂conj(S)/∂V*

Returns:
    J11, J12, J21, J22  (all full matrices, not Diagonal)
"""
function build_complex_jacobian(Ybus, V)
    I = Ybus * V
    n = length(V)

    # J11 = diag(conj(I))
    J11 = Matrix(Diagonal(conj.(I)))

    # J12 = diag(V) * conj(Ybus)
    J12 = Matrix(Diagonal(V)) * conj.(Ybus)

    # J21 = conj(J12)
    J21 = conj.(J12)

    # J22 = conj(J11)
    J22 = conj.(J11)

    return J11, J12, J21, J22
end
