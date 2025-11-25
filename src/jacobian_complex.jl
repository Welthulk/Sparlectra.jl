# jacobian_complex.jl
# Jacobian prototype for complex-state NR

"""
    build_complex_jacobian(Ybus, V)

Constructs the Jacobian matrix ∂S/∂V* for the complex-state
Newton–Raphson formulation.

This is a prototype — replace the placeholder expression with
the correct Wirtinger derivative.
"""
function build_complex_jacobian(Ybus, V)
    n = length(V)
    J = zeros(ComplexF64, n, n)

    # TODO: correct derivative of
    #       S_i = V_i * conj((Ybus * V)_i)
    for i in 1:n
        for k in 1:n
            # Placeholder approximation
            J[i, k] = conj(Ybus[i, k]) * V[i]
        end
    end

    return J
end
