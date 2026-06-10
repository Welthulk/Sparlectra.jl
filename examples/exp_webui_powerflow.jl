using Sparlectra

function main()
  server = start_sparlectra_webui(
    host = "127.0.0.1",
    port = 8080,
    output_root = joinpath(@__DIR__, "_out", "powerflow_service"),
    open_browser = false,
  )
  println("Sparlectra Web UI is available at ", server.url)
  println("Press Ctrl+C to stop the local server.")
  wait(server.task)
end

Base.invokelatest(main)
