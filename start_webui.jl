using Sparlectra

# AnalyticLoadFlow.jl is only a weak dependency of Sparlectra (the APSLF
# solver, loaded via ext/SparlectraAnalyticLoadFlowExt.jl). Load it here if
# it happens to be installed so the Web UI's APSLF solver/start-value
# options actually work; if it isn't installed, skip it gracefully instead
# of forcing everyone who runs this script to have it. Install it with
# `Pkg.add("AnalyticLoadFlow")` in your default (global) Julia environment
# so it is available regardless of which project you launch this script
# from, without adding a hard dependency to Sparlectra's own Project.toml.
if Base.find_package("AnalyticLoadFlow") !== nothing
  @eval using AnalyticLoadFlow
  @info "AnalyticLoadFlow.jl loaded — APSLF solver and start-value options are active."
else
  @info "AnalyticLoadFlow.jl not found; power_flow.solver=apslf / apslf_start will show a clear \"not installed\" error if selected. Run `Pkg.add(\"AnalyticLoadFlow\")` to enable it."
end

function main()
  server = Sparlectra.start_sparlectra_webui(open_browser = true, warmup = true)
  # nothing = a Sparlectra Web UI already runs on the port; it was opened in
  # the browser instead — there is no new server task to wait on.
  server === nothing && return nothing

  try
    wait(server.task)
  catch err
    if err isa InterruptException
      @info "Ctrl+C received; closing Sparlectra Web UI."
      close(server)
      return nothing
    end
    rethrow()
  end
end

# getfield defers the binding lookup to call time: VS Code's inline eval runs
# the whole file as one world age, where a direct `main` access would trigger
# the Julia 1.12 "binding accessed before its definition world" warning.
Base.invokelatest(getfield(@__MODULE__, :main))