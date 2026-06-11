using Sparlectra
using Sockets
using Test

function _webui_test_form(casefile, config_file, output_root)
  return Dict(
    "casefile" => casefile,
    "config_file" => config_file,
    "output_root" => output_root,
    "power_flow_tol" => "1e-8",
    "power_flow_max_iter" => "80",
    "power_flow_autodamp" => "on",
    "power_flow_autodamp_min" => "0.05",
    "power_flow_qlimits_enabled" => "on",
    "power_flow_wrong_branch_detection" => "warn",
    "power_flow_start_angle_mode" => "dc",
    "power_flow_start_voltage_mode" => "profile_blend",
    "output_logfile_results" => "compact",
    "benchmark_samples" => "10",
    "benchmark_seconds" => "1.0",
  )
end

function run_webui_tests()
  @testset "Local PowerFlow Web UI" begin
    @test isdefined(Sparlectra, :start_sparlectra_webui)
    @test_throws ArgumentError start_sparlectra_webui(host = "0.0.0.0")

    @testset "Markdown-backed contextual help and documentation" begin
      topic = "power_flow.start_mode.voltage_mode"
      @test haskey(Sparlectra.WEBUI_HELP_TOPICS, topic)
      topic_metadata = Sparlectra.resolve_webui_help_topic(topic)
      @test topic_metadata !== nothing
      doc_metadata = Sparlectra.resolve_webui_doc_page(topic_metadata.page)
      @test doc_metadata !== nothing
      doc_path = joinpath(@__DIR__, "..", "docs", "src", doc_metadata.file)
      @test isfile(doc_path)

      markdown_text = Sparlectra.load_webui_markdown_document(topic_metadata.page)
      @test markdown_text !== nothing
      section = Sparlectra.extract_webui_markdown_section(markdown_text, topic_metadata.heading)
      @test section !== nothing
      @test occursin(topic_metadata.selector, section)
      excerpt = Sparlectra.load_webui_help_excerpt(topic)
      @test excerpt !== nothing
      @test occursin("power_flow.start_mode.voltage_mode", excerpt)
      @test occursin("YAML path", excerpt)

      synthetic = "## Target\n\n- item\n\n```julia\nx = 1\n```\n\n### Child\nchild text\n\n## Next\nnot included"
      extracted = Sparlectra.extract_webui_markdown_section(synthetic, "Target")
      @test occursin("- item", extracted)
      @test occursin("```julia", extracted)
      @test occursin("### Child", extracted)
      @test !occursin("not included", extracted)

      help_response = Sparlectra.route_sparlectra_webui("GET", "/help/$(topic)")
      @test help_response.status == 200
      help_html = String(help_response.body)
      @test occursin("Start voltage mode", help_html)
      @test occursin("power_flow.start_mode.voltage_mode", help_html)
      @test occursin("/docs/powerflow_configuration", help_html)

      unknown_help = Sparlectra.route_sparlectra_webui("GET", "/help/power_flow.unknown")
      @test unknown_help.status == 404
      @test occursin("Unknown help topic", String(unknown_help.body))

      docs_index = Sparlectra.route_sparlectra_webui("GET", "/docs")
      @test docs_index.status == 200
      docs_index_html = String(docs_index.body)
      for page in keys(Sparlectra.WEBUI_DOC_PAGES)
        @test occursin("/docs/$(page)", docs_index_html)
      end

      docs_page = Sparlectra.route_sparlectra_webui("GET", "/docs/powerflow_configuration")
      @test docs_page.status == 200
      @test occursin("Power-Flow Configuration", String(docs_page.body))

      for unsafe_page in ("../Project.toml", "../../src/Sparlectra.jl", "/etc/passwd", raw"C:\Users\x\file.txt")
        unsafe = Sparlectra.handle_webui_doc_page(unsafe_page)
        @test unsafe.status == 404
        @test occursin("Documentation page not found", String(unsafe.body))
        unsafe_route = Sparlectra.route_sparlectra_webui("GET", "/docs/" * Sparlectra._webui_urlencode(unsafe_page))
        @test unsafe_route.status == 404
      end
    end

    mktempdir() do tmpdir
      casefile = _write_api_test_case(joinpath(tmpdir, "case_webui.m"))
      config_file = joinpath(tmpdir, "config_template.yaml")
      cp(Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH, config_file)
      output_root = joinpath(tmpdir, "runs")
      form = _webui_test_form(casefile, config_file, output_root)

      request = Sparlectra.powerflow_webui_request(form)
      @test request["casefile"] == casefile
      @test request["config_file"] == config_file
      @test request["output_root"] == output_root
      overrides = request["config_overrides"]
      @test Set(keys(overrides)) == Set(key for (key, _, _) in Sparlectra._WEBUI_FORM_CONFIG_FIELDS)
      @test overrides["power_flow.tol"] == 1.0e-8
      @test overrides["power_flow.max_iter"] == 80
      @test overrides["power_flow.autodamp"]
      @test overrides["benchmark.enabled"] === false

      form_html = Sparlectra.render_powerflow_form(output_root = output_root)
      for field in ("casefile", "config_file", "output_root", "power_flow_tol", "power_flow_max_iter", "power_flow_autodamp", "power_flow_autodamp_min", "power_flow_qlimits_enabled", "power_flow_wrong_branch_detection", "power_flow_start_angle_mode", "power_flow_start_voltage_mode", "output_logfile_results", "benchmark_enabled", "benchmark_samples", "benchmark_seconds")
        @test occursin("name=\"$(field)\"", form_html)
      end

      @test occursin("src=\"/assets/logo.png\"", form_html)
      @test occursin("alt=\"Sparlectra.jl logo\"", form_html)
      @test occursin("<span>Sparlectra.jl</span>", form_html)
      @test occursin("class=\"help-link\"", form_html)
      @test occursin("href=\"/help/power_flow.start_mode.voltage_mode\"", form_html)

      logo_response = Sparlectra.route_sparlectra_webui("GET", "/assets/logo.png")
      @test logo_response.status == 200
      @test ("Content-Type" => "image/png") in logo_response.headers
      @test logo_response.body == read(joinpath(@__DIR__, "..", "docs", "src", "assets", "logo.png"))
      for target in ("/assets/tablestyle.css", "/assets/../Project.toml", "/assets/logo.png/extra")
        @test Sparlectra.route_sparlectra_webui("GET", target).status == 404
      end
      result = Sparlectra.handle_powerflow_run(form)
      @test result["success"]
      run_id = result["run_id"]
      @test !isempty(run_id)
      result_response = Sparlectra.handle_powerflow_result(run_id)
      result_html = String(result_response.body)
      @test result_response.status == 200
      for field in ("run_id", "status", "iterations", "final_mismatch")
        @test occursin(field, result_html)
      end
      @test occursin(run_id, result_html)

      artifacts = list_powerflow_artifacts(run_id)
      artifact_names = Set(artifact["name"] for artifact in artifacts)
      @test Set(("result.json", "run.log", "effective_config.yaml")) ⊆ artifact_names
      artifact_response = Sparlectra.handle_powerflow_artifacts(run_id)
      artifact_html = String(artifact_response.body)
      for name in ("result.json", "run.log", "effective_config.yaml")
        @test occursin(name, artifact_html)
      end
      result_artifact = Sparlectra.handle_powerflow_artifact(run_id, "result.json")
      @test result_artifact.status == 200
      result_artifact_html = String(result_artifact.body)
      @test occursin("Artifact: result.json", result_artifact_html)
      download = Sparlectra.handle_powerflow_artifact_download(run_id, "result.json")
      @test download.status == 200
      @test any(header -> header.first == "Content-Disposition", download.headers)

      for unsafe_name in ("../result.json", "../../Project.toml", "/etc/passwd", raw"C:\Users\x\file.txt")
        unsafe = Sparlectra.handle_powerflow_artifact(run_id, unsafe_name)
        @test unsafe.status == 400
        @test occursin("Unsafe artifact name rejected", String(unsafe.body))
      end

      lock(Sparlectra._POWERFLOW_SERVICE_LOCK) do
        empty!(Sparlectra._POWERFLOW_SERVICE_RUNS)
      end
      refreshed = Sparlectra.handle_powerflow_refresh(output_root)
      @test run_id in refreshed["loaded_runs"]
      history_response = Sparlectra.handle_powerflow_history(output_root)
      @test history_response.status == 200
      history_html = String(history_response.body)
      @test occursin(run_id, history_html)
      @test occursin(run_id, String(Sparlectra.handle_powerflow_result(run_id).body))
      shared_pages = (
        form_html,
        result_html,
        artifact_html,
        result_artifact_html,
        history_html,
      )
      for page_html in shared_pages
        @test occursin("src=\"/assets/logo.png\"", page_html)
        @test occursin("<span>Sparlectra.jl</span>", page_html)
      end

      source = join(read(joinpath(@__DIR__, "..", "src", "webui", file), String) for file in ("docs.jl", "forms.jl", "handlers.jl", "routes.jl", "views.jl", "webui.jl"))
      forbidden = ["run" * "pf!", "run" * "pf_rectangular!", "run_complex_nr_" * "rectangular"]
      @test all(name -> !occursin(name, source), forbidden)
      @test !occursin("Voltage-magnitude/angle blend strategy.", source)

      probe = listen(ip"127.0.0.1", UInt16(0))
      port = Int(getsockname(probe)[2])
      close(probe)
      server = start_sparlectra_webui(port = port, output_root = output_root)
      try
        @test server.url == "http://127.0.0.1:$(port)/powerflow"
        socket = connect(ip"127.0.0.1", UInt16(port))
        write(socket, "GET /powerflow HTTP/1.1\r\nHost: 127.0.0.1\r\nConnection: close\r\n\r\n")
        response_text = read(socket, String)
        close(socket)
        @test occursin("HTTP/1.1 200 OK", response_text)
        @test occursin("PowerFlow run", response_text)

        logo_socket = connect(ip"127.0.0.1", UInt16(port))
        write(logo_socket, "GET /assets/logo.png HTTP/1.1\r\nHost: 127.0.0.1\r\nConnection: close\r\n\r\n")
        logo_response_text = read(logo_socket, String)
        close(logo_socket)
        @test occursin("HTTP/1.1 200 OK", logo_response_text)
        @test occursin("Content-Type: image/png", logo_response_text)
      finally
        close(server)
      end
    end
  end
  return nothing
end
