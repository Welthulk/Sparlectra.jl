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
  server = Sparlectra.start_sparlectra_webui(open_browser = true)

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

Base.invokelatest(main)