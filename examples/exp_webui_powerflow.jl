using Sparlectra

function main()
  server = start_sparlectra_webui(
    host = "127.0.0.1",
    port = 8080,
    output_root = joinpath(@__DIR__, "_out", "powerflow_service"),
    open_browser = true,
    auto_shutdown_on_browser_close = true,
    browser_heartbeat_timeout_seconds = 15.0,
  )
  println("Sparlectra Web UI is available at ", server.url)
  application_root = Sparlectra._webui_application_root()
  println("Case selections are loaded from ", joinpath(application_root, "data", "mpower"))
  println("Configuration selections are loaded from ", joinpath(application_root, "examples"))
  println("Contextual help and selected documentation are available at ", replace(server.url, "/powerflow" => "/docs"))
  println("Use Stop Web UI in the browser, or press Ctrl+C as a fallback.")
  wait(server.task)
end

Base.invokelatest(main)
