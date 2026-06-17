using Sparlectra

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