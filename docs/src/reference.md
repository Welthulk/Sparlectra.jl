# Function Reference

The generated API reference, one page per source directory. Each page lists
the exported (public) bindings first, followed by the internal helpers of
that directory.

- [API and Services](reference_api.md): service entry points, run management, Web UI backend
- [ACPFlow Runner](reference_acpflow.md): the configured power-flow execution path
- [Core Model](reference_core.md): network, nodes, branches, transformers, prosumers
- [Rectangular Power Flow](reference_powerflow_rectangular.md): the Newton-Raphson solver in rectangular coordinates
- [DC Power Flow](reference_powerflow_dc.md): the linear DC approximation
- [Short Circuit](reference_shortcircuit.md): IEC 60909 short-circuit calculation
- [Format Adapters](reference_adapters.md): MATPOWER, CGMES, DTF, and the Sparlectra Case Format
- [State Estimation and Measurements](reference_stateestimation.md): WLS estimator, measurement model, tap estimation
- [Configuration Internals](reference_config.md): configuration loading and YAML handling
- [Controllers](reference_controller.md): the control framework and its controllers
- [APSLF Bridge](reference_apslf.md): the AnalyticLoadFlow solver bridge
- [Shared Numerics](reference_numerics.md): condition-number and sparse-inverse utilities
- [Contingency](reference_contingency.md): N-1/N-k contingency analysis
- [Synthetic Grids](reference_synthetic.md): generated tiled test networks
