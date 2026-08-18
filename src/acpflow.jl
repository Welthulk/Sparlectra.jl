# Copyright 2023–2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

# file: src/acpflow.jl
# purpose: include hub for the configuration-driven AC power-flow framework
#          runner: start modes, net cache, import context, execution, status,
#          output, and entry points
# Configuration-driven AC power-flow framework runner.
include("acpflow/start_modes.jl")
include("acpflow/net_cache.jl")
include("acpflow/import_context.jl")
include("acpflow/execution.jl")
include("acpflow/apslf_execution.jl")
include("acpflow/status.jl")
include("acpflow/output.jl")
include("acpflow/entrypoint.jl")
