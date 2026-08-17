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
  # webui.warmup (Advanced option in the form, default true) governs the
  # launcher; the library entry point itself keeps warm-up off by default.
  config_path = Sparlectra.default_webui_config_path()
  warmup = try
    isfile(config_path) ? Sparlectra.load_sparlectra_config(config_path; reload = true).webui.warmup : true
  catch
    true
  end
  # the warm-up exists to pay the JIT cost early; with the fast-start
  # sysimage loaded that work is already compiled into the image, so the
  # PowerFlow page can serve immediately (measured: page after ~2 s with
  # the image versus ~15 s waiting for the warm-up to finish)
  # isdefined guard: when running on a sysimage built from an older code
  # state, the baked Sparlectra module may predate this helper; an outdated
  # image must never break the start
  if warmup && isdefined(Sparlectra, :webui_sysimage_active) && Base.invokelatest(Sparlectra.webui_sysimage_active)
    @info "Fast-start sysimage active; skipping the Web UI warm-up (already compiled into the image)."
    warmup = false
  end
  server = Sparlectra.start_sparlectra_webui(open_browser = true, warmup = warmup)
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

# VS Code's inline eval runs the whole file as one world age. Since Julia
# 1.12.7 even a direct getfield counts as a binding access from the old
# world, so the lookup itself must also go through invokelatest.
Base.invokelatest(Base.invokelatest(getfield, @__MODULE__, :main))