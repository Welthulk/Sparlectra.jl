# Copyright 2023–2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.

# file: examples/exp_configuration.jl

using Sparlectra

"""
    main(; path=Sparlectra.USER_SPARLECTRA_CONFIG_PATH, reload=false)

Load a Sparlectra YAML file through the central typed configuration layer and
print the effective module sections. The example intentionally stops after
configuration loading so it is safe for quick developer smoke checks.
"""
function main(; path::AbstractString = Sparlectra.USER_SPARLECTRA_CONFIG_PATH, reload::Bool = false)
  cfg = Sparlectra.load_sparlectra_config(path; reload = reload)
  Sparlectra.print_effective_config(cfg)
  return cfg
end

if get(ENV, "SPARLECTRA_CONFIGURATION_EXAMPLE_NO_MAIN", "") != "1"
  Base.invokelatest(() -> main())
end
