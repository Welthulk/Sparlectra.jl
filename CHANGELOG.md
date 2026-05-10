# Changelog

## 0.7.5

- **Bugfix**: Fixed the MATPOWER import example so Julia 1.12 / Revise entry points use `Base.invokelatest`, and corrected keyword forwarding for damping and start-projection solver options, and prevents the MATPOWER fallback diagnostic from enabling PQ-generator voltage-dependent controllers for non-rectangular methods.
