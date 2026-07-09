# Transformer convention scan

The `current` rows are executed production results. Other rows enumerate the requested ratio and shunt candidate modes and mark them as diagnostic candidates that still require an explicit admittance/ratio override before they can be selected as a default.

- Case A: converged=true, iterations=5, final mismatch=1.539879335155092e-13, 400-kV mean bias=-1.894 kV, max |dV|=1.937 kV, slack ΔP/ΔQ=-0.035/2.759, loss ΔP/ΔQ=-1.623/2.759, max branch ΔP/ΔQ=0.24/1.169.
- Case B: converged=true, iterations=5, final mismatch=1.4876988529977098e-13, 400-kV mean bias=-1.858 kV, max |dV|=1.875 kV, slack ΔP/ΔQ=-0.005/2.361, loss ΔP/ΔQ=-1.588/2.736, max branch ΔP/ΔQ=0.236/1.218.
- Case C: converged=true, iterations=5, final mismatch=1.5422098172789457e-13, 400-kV mean bias=-2.456 kV, max |dV|=2.671 kV, slack ΔP/ΔQ=0.276/4.492, loss ΔP/ΔQ=-1.257/4.492, max branch ΔP/ΔQ=1.741/8.868.
- Case D: converged=true, iterations=5, final mismatch=1.3539169785303784e-13, 400-kV mean bias=0.003 kV, max |dV|=1.095 kV, slack ΔP/ΔQ=-0.072/-0.806, loss ΔP/ΔQ=-1.736/-0.806, max branch ΔP/ΔQ=0.23/0.558.
- Case E: converged=true, iterations=5, final mismatch=2.867515326986552e-13, 400-kV mean bias=-2.039 kV, max |dV|=2.095 kV, slack ΔP/ΔQ=0.107/3.444, loss ΔP/ΔQ=-1.45/3.444, max branch ΔP/ΔQ=2.176/2.805.

No production default is changed by this report because the non-current candidate modes have not produced a clearly winning executed result.
