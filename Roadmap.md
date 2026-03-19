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

