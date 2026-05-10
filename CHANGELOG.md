# Changelog

## 0.7.5

- **Feature**: Added configurable bus-shunt modeling with `bus_shunt_model = "admittance"` for classic Y-bus diagonal stamping and `bus_shunt_model = "voltage_dependent_injection"` for |V|²-dependent nonlinear injection treatment in the rectangular power-flow path. The default remains `"admittance"` to preserve existing results.

- **Bugfix**: Fixed the MATPOWER import example so Julia 1.12 / Revise entry points use `Base.invokelatest`, and corrected keyword forwarding for damping and start-projection solver options, and prevents the MATPOWER fallback diagnostic from enabling PQ-generator voltage-dependent controllers for non-rectangular methods.
