using Sparlectra

function main()
  server = Sparlectra.start_sparlectra_webui(
    host = "127.0.0.1",
    port = 8080,
    open_browser = true,
    warmup = true,
  )
  println("Sparlectra Web UI is available at ", server.url)
  println("Use Stop Web UI in the browser, or press Ctrl+C as a fallback.")
  wait(server.task)
end

Base.invokelatest(main)
