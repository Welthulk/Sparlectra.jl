# Copyright 2023–2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

# file: src/acpflow/acpflow.jl
# purpose: include hub for the configuration-driven AC power-flow framework
#          runner: start modes, net cache, import context, execution, status,
#          output, and entry points
# Configuration-driven AC power-flow framework runner.
include("start_modes.jl")
include("net_cache.jl")
include("import_context.jl")
include("execution.jl")
include("auto_powerflow.jl")
include("apslf_execution.jl")
include("status.jl")
include("output.jl")
include("entrypoint.jl")
