# Sparlectra State Estimation Roadmap

### 1. Bad-data detection and statistical diagnostics

Recommended additions:

* global bad-data consistency check,
* chi-square style objective evaluation,
* largest-normalized-residual workflow,
* measurement ranking / suspicion report,
* optional deactivate-and-rerun workflow for diagnostics.

### 2. Dedicated diagnostics API

Add a clear public API around result interpretation.

Suggested functions:

* `runse_diagnostics(...)`
* `validate_measurements(...)`

Suggested outputs:

* normalized residuals,
* ranking of suspicious measurements,
* summary flags for observability / criticality / data quality,
* machine-readable diagnostic report object.

### 3. Measurement management ergonomics

Suggested additions:

* helper to activate/deactivate measurements by ID,
* helper to filter measurements by type / bus / branch,
* helper to ensure measurement IDs are unique,
* standard residual and ranking reports.

### 4. ZIB as hard constraints

Today, passive buses are represented through zero-injection pseudo
measurements. This is practical and should remain supported, but a later
extension could add a dedicated equality-constraint formulation for
zero-injection buses.

### 5. Jacobian and solver efficiency

The present finite-difference Jacobian is a good implementation strategy for an
initial release, but future optimization could include:

* analytical Jacobians for supported measurement types,
* sparse Jacobian / gain-matrix handling,
* reduced allocation in Jacobian assembly,
* more specialized linear solvers for SE normal equations.

### 6. Tight integration of Power Flow and State Estimation

The State Estimation should not evolve as a standalone component but be
tightly coupled with the existing power flow (PF) formulation.

Recommended additions:

* shared network model between PF and SE (same Bus/Branch data structures),
* reuse of admittance model (Ybus) and stamping logic,
* consistent handling of transformer taps, phase shifts and shunts,
* consistent treatment of PV/PQ/Slack buses where applicable,
* optional reuse of PF Jacobian building blocks for SE,
* unified mismatch/residual evaluation infrastructure,
* optional NR-based “polish” step after SE convergence.

Goal:

Avoid duplication of model logic and ensure numerical consistency between PF and SE.

---

### 7. Measurement model extensions (PMU / hybrid SE)

Prepare the measurement framework for hybrid state estimation.

Recommended additions:

* support for phasor measurements:
  * `VoltagePhasorMeas` (|V| + angle or complex V),
  * `CurrentPhasorMeas` (branch or injection currents),
* optional angle-only measurements (`VaMeas`) as intermediate step,
* clear separation between SCADA and PMU measurement classes,
* extension of measurement noise model for different device classes,
* time-synchronization flag / metadata (future use).

Goal:

Enable a future transition from classical WLS (SCADA-only) to hybrid SE.

---

### 8. Hierarchical / multi-area SE preparation

Even if not implemented immediately, the architecture should allow
future decomposition.

Recommended additions:

* optional area / zone assignment for buses and branches,
* identification of boundary buses / tie-lines,
* ability to extract subnetworks programmatically,
* grouping of measurements by area,
* interfaces for exchanging boundary states between areas.

Future extension:

* two-step hierarchical SE (local + global),
* reduced data exchange between areas,
* improved scalability for large systems.

---

### 9. Observability and redundancy as first-class features

Observability analysis should be treated as a core module, not only as a helper.

Recommended additions:

* unified API for:
  * global observability,
  * local observability,
  * critical measurements,
  * redundancy metrics (m-n, rho),
* structured result object for observability reports,
* integration with diagnostics (bad-data detection),
* optional structural vs. numerical observability modes.

Goal:

Provide transparent insight into measurement system quality and estimator reliability.

---

### 10. Measurement generation and reproducibility

Synthetic measurement generation is essential for testing and benchmarking.

Recommended additions:

* standardized measurement generator from PF results,
* configurable noise models (Gaussian, scaling per measurement type),
* reproducible random seeds,
* ability to generate incomplete / degraded measurement sets,
* helper to simulate bad data (outliers, bias).

Goal:

Enable systematic testing of SE, diagnostics and observability workflows.

---

### 11. Result model and reporting

The SE result should be a structured, extensible object.

Recommended additions:

* unified result type including:
  * estimated state,
  * residual vector,
  * normalized residuals,
  * objective value,
  * convergence info,
* optional comparison to PF reference (if available),
* helper functions for:
  * printing summaries,
  * exporting reports,
  * machine-readable diagnostics.

---

### 12. Numerical robustness and scaling

Prepare the implementation for larger systems.

Recommended additions:

* consistent scaling of residuals and weights,
* conditioning diagnostics (e.g. gain matrix),
* detection of ill-conditioned cases,
* optional damping / step control,
* improved initialization strategies (flat start vs PF-based).

---

### 13. Future: optimization and security assessment interface

While not part of the initial SE scope, the design should remain compatible
with higher-level applications.

Recommended considerations:

* reuse of network and measurement model for:
  * contingency analysis,
  * security assessment,
  * future OPF/SCOPF workflows,
* clean separation between:
  * model,
  * solver,
  * application layer.

Goal:

Allow Sparlectra to evolve towards a consistent grid analysis environment.