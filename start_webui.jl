using Sparlectra

server = start_sparlectra_webui(host = "127.0.0.1", port = 8080, output_root = "results/powerflow_service", open_browser = true)

wait(server.task)