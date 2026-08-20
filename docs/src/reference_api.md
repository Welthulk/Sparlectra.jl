# [API and Web UI Service Reference](@id reference_api)

The N-1 contingency batch API (`runContingencies!`, `generateN1Branches`,
`generateContingenciesFromFOR001`, result printer and CSV writer) is
documented on its own page: [N-1 Contingency Analysis](contingency.md).

```@autodocs
Modules = [Sparlectra]
Pages = [
    "api/api_types.jl",
    "api/config_overrides.jl",
    "api/artifacts.jl",
    "api/run_api.jl",
    "api/run_self_check.jl",
    "api/powerflow_service.jl",
    "api/run_index.jl",
    "api/artifact_registry.jl",
    "api/serialization.jl"
]
```
