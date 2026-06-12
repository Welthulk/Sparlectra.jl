using Sparlectra
using Sockets
using Test

function _webui_test_form(casefile, config_file, output_root)
  return Dict(
    "casefile" => "",
    "casefile_manual" => casefile,
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
    "performance_timing" => "compact",
    "run_diagnostics" => "on",
    "benchmark_samples" => "10",
    "benchmark_seconds" => "1.0",
  )
end

function _webui_http_request(port::Integer, method::AbstractString, target::AbstractString; body::AbstractString = "")::String
  socket = connect(ip"127.0.0.1", UInt16(port))
  headers = isempty(body) ? "" : "Content-Type: application/x-www-form-urlencoded\r\nContent-Length: $(ncodeunits(body))\r\n"
  write(socket, "$(method) $(target) HTTP/1.1\r\nHost: 127.0.0.1\r\n$(headers)Connection: close\r\n\r\n$(body)")
  response = read(socket, String)
  close(socket)
  return response
end

function run_webui_tests()
  @testset "Local PowerFlow Web UI" begin
    @test isdefined(Sparlectra, :start_sparlectra_webui)
    @test_throws ArgumentError start_sparlectra_webui(host = "0.0.0.0")

    registry_before_warmup = Set(keys(Sparlectra._POWERFLOW_SERVICE_RUNS))
    warmup_output = Ref("")
    warmup_runner = function(; output_dir, kwargs...)
      warmup_output[] = output_dir
      write(joinpath(output_dir, "warmup-marker.txt"), "compiled\n")
      return (success = true, reason = nothing, message = nothing)
    end
    warmup_result = Sparlectra._run_sparlectra_webui_warmup(mktempdir(); runner = warmup_runner)
    @test warmup_result.success
    @test !ispath(warmup_output[])
    @test Set(keys(Sparlectra._POWERFLOW_SERVICE_RUNS)) == registry_before_warmup

    @testset "Case and configuration selection" begin
      mktempdir() do start_dir
        application_root = joinpath(start_dir, "Sparlectra")
        case_directory = joinpath(application_root, "data", "mpower")
        config_directory = joinpath(application_root, "examples")
        mkpath(case_directory)
        mkpath(config_directory)
        case14 = joinpath(case_directory, "case14.m")
        case118 = joinpath(case_directory, "case118.m")
        write(case14, "function mpc = case14\nend\n")
        write(case118, "function mpc = case118\nend\n")
        write(joinpath(case_directory, "case14.jl"), "nothing\n")
        primary_config = joinpath(config_directory, "configuration.yaml.example")
        secondary_config = joinpath(config_directory, "study.yaml.example")
        write(primary_config, "power_flow: {}\n")
        write(secondary_config, "power_flow: {}\n")
        write(joinpath(config_directory, "README.md"), "ignored\n")

        @test Sparlectra._webui_application_root(start_dir) == application_root
        @test Sparlectra._webui_casefile_options(application_root) == ["case118.m", "case14.jl", "case14.m"]
        @test Sparlectra._webui_config_file_options(application_root) == [primary_config, secondary_config]

        selection_html = Sparlectra.render_powerflow_form(
          application_root = application_root,
          selected_casefile = "case14.m",
          selected_config_file = secondary_config,
        )
        @test occursin("<select id=\"casefile\" name=\"casefile\">", selection_html)
        @test occursin("<option value=\"\">-- choose existing case --</option>", selection_html)
        @test occursin("<option value=\"case14.m\" selected>case14.m</option>", selection_html)
        @test occursin("<option value=\"case14.jl\">case14.jl</option>", selection_html)
        @test occursin("<option value=\"case118.m\">case118.m</option>", selection_html)
        @test occursin("<input id=\"casefile_manual\" name=\"casefile_manual\" value=\"\" placeholder=\"case14.m\">", selection_html)
        @test !occursin("available-casefiles", selection_html)
        @test occursin("<select id=\"config_file\" name=\"config_file\" required>", selection_html)
        @test occursin("<option value=\"$(secondary_config)\" selected>study.yaml.example</option>", selection_html)
        @test !occursin("README.md", selection_html)

        rm(case14)
        rm(case118)
        manual_case = joinpath(case_directory, "manual_case.m")
        fallback_html = Sparlectra.render_powerflow_form(
          application_root = application_root,
          selected_casefile = manual_case,
        )
        @test occursin("<select id=\"casefile\" name=\"casefile\">", fallback_html)
        @test occursin("<input id=\"casefile_manual\" name=\"casefile_manual\" value=\"$(manual_case)\"", fallback_html)
      end
    end

    @testset "Standalone app-window launcher" begin
      url = "http://127.0.0.1:8080/powerflow"
      linux_lookup = name -> name == "google-chrome" ? "/opt/google/chrome" : nothing
      linux_command = Sparlectra._webui_app_command(url; platform = :linux, executable_lookup = linux_lookup)
      @test linux_command !== nothing
      @test linux_command.exec == ["/opt/google/chrome", "--app=$(url)"]

      windows_lookup = name -> name == "msedge.exe" ? raw"C:\Browser\msedge.exe" : nothing
      windows_command = Sparlectra._webui_app_command(url; platform = :windows, executable_lookup = windows_lookup, environment = Dict{String,String}())
      @test windows_command !== nothing
      @test windows_command.exec == [raw"C:\Browser\msedge.exe", "--app=$(url)"]

      macos_exists = path -> path == "/Applications/Google Chrome.app"
      macos_command = Sparlectra._webui_app_command(url; platform = :macos, path_exists = macos_exists)
      @test macos_command !== nothing
      @test macos_command.exec == ["open", "-na", "/Applications/Google Chrome.app", "--args", "--app=$(url)"]

      missing_lookup = _ -> nothing
      @test Sparlectra._webui_app_command(url; platform = :linux, executable_lookup = missing_lookup) === nothing
    end

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

      link_markdown = "## Local section\n\n[Configuration](powerflow_configuration.md)\n[Start modes](powerflow_configuration.md#start-mode-options)\n[Dot relative](./powerflow_configuration.md)\n[External](https://example.com/reference)\n[Local anchor](#local-section)\n[Unknown](unknown.md)"
      link_html = Sparlectra.render_webui_markdown(link_markdown; current_page = "webui")
      @test occursin("href=\"/docs/powerflow_configuration\"", link_html)
      @test occursin("href=\"/docs/powerflow_configuration#start-mode-options\"", link_html)
      @test count("href=\"/docs/powerflow_configuration\"", link_html) == 2
      @test occursin("href=\"https://example.com/reference\"", link_html)
      @test occursin("href=\"#local-section\"", link_html)
      @test occursin("id=\"local-section\"", link_html)
      @test !occursin("href=\"unknown.md\"", link_html)
      @test occursin("aria-disabled=\"true\"", link_html)

      for unsafe_target in ("../Project.toml", "../../src/Sparlectra.jl", "/etc/passwd", raw"C:\Users\x\file.txt")
        unsafe_html = Sparlectra.render_webui_markdown("[unsafe]($(unsafe_target))"; current_page = "webui")
        @test occursin("aria-disabled=\"true\"", unsafe_html)
        @test !occursin("href=\"$(unsafe_target)\"", unsafe_html)
      end

      @test Set(values(Sparlectra.WEBUI_FORM_HELP_TOPICS)) == Set(keys(Sparlectra.WEBUI_HELP_TOPICS))
      for (configured_topic, metadata) in Sparlectra.WEBUI_HELP_TOPICS
        page_metadata = Sparlectra.resolve_webui_doc_page(metadata.page)
        @test page_metadata !== nothing
        help_doc_path = joinpath(@__DIR__, "..", "docs", "src", page_metadata.file)
        @test isfile(help_doc_path)
        help_markdown = Sparlectra.load_webui_markdown_document(metadata.page)
        @test help_markdown !== nothing
        help_section = Sparlectra.extract_webui_markdown_section(help_markdown, metadata.heading)
        @test help_section !== nothing
        isempty(metadata.selector) || @test occursin(metadata.selector, help_section)
        @test Sparlectra.load_webui_help_excerpt(configured_topic) !== nothing
      end

      help_response = Sparlectra.route_sparlectra_webui("GET", "/help/$(topic)")
      @test help_response.status == 200
      help_html = String(help_response.body)
      @test occursin("Start voltage mode", help_html)
      @test occursin("power_flow.start_mode.voltage_mode", help_html)
      @test occursin("/docs/powerflow_configuration", help_html)
      @test occursin("class=\"panel help-page help-panel\"", help_html)

      @test occursin("class=\"button back-button\"", help_html)
      @test occursin("document.referrer.startsWith(location.origin)", help_html)
      @test occursin("history.back()", help_html)
      @test occursin("href=\"/powerflow\"", help_html)

      for representative_topic in (
        "power_flow.start_mode.voltage_mode",
        "power_flow.start_mode.angle_mode",
        "power_flow.wrong_branch_detection",
        "benchmark.enabled",
      )
        representative_help = Sparlectra.route_sparlectra_webui("GET", "/help/$(representative_topic)")
        @test representative_help.status == 200
        @test occursin("class=\"panel help-page help-panel\"", String(representative_help.body))
      end

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
      docs_page_html = String(docs_page.body)
      @test occursin("Power-Flow Configuration", docs_page_html)
      @test occursin("class=\"panel docs-page docs-content\"", docs_page_html)
      @test occursin("id=\"start-mode-options\"", docs_page_html)

      webui_docs_page = Sparlectra.route_sparlectra_webui("GET", "/docs/webui")
      @test webui_docs_page.status == 200
      webui_docs_html = String(webui_docs_page.body)
      @test occursin("href=\"/docs/powerflow_configuration\"", webui_docs_html)
      @test occursin("href=\"/docs/powerflow_configuration#start-mode-options\"", webui_docs_html)
      @test occursin("href=\"/docs/performance_profiling\"", webui_docs_html)

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

      request = Sparlectra.powerflow_webui_request(form; default_output_root = output_root)
      @test request["casefile"] == casefile
      existing_form = copy(form)
      existing_form["casefile"] = "  case118.jl  "
      existing_form["casefile_manual"] = "  "
      @test Sparlectra.powerflow_webui_request(existing_form; default_output_root = output_root)["casefile"] == "case118.jl"

      typed_form = copy(form)
      typed_form["casefile"] = "case9.m"
      typed_form["casefile_manual"] = "  case9241pegase.m  "
      @test Sparlectra.powerflow_webui_request(typed_form; default_output_root = output_root)["casefile"] == "case9241pegase.m"

      empty_case_form = copy(form)
      empty_case_form["casefile"] = ""
      empty_case_form["casefile_manual"] = ""
      @test_throws ArgumentError Sparlectra.powerflow_webui_request(empty_case_form; default_output_root = output_root)
      empty_case_response = Sparlectra.route_sparlectra_webui("POST", "/powerflow/run", empty_case_form; output_root = output_root)
      @test empty_case_response.status == 400
      @test occursin("Select an existing MATPOWER case or type a case name.", String(empty_case_response.body))
      @test request["config_file"] == config_file
      @test request["output_root"] == output_root
      untrusted_form = copy(form)
      untrusted_form["output_root"] = joinpath(tmpdir, "outside")
      @test Sparlectra.powerflow_webui_request(untrusted_form; default_output_root = output_root)["output_root"] == output_root
      overrides = request["config_overrides"]
      @test Set(keys(overrides)) == Set(key for (key, _, _) in Sparlectra._WEBUI_FORM_CONFIG_FIELDS)
      @test overrides["power_flow.tol"] == 1.0e-8
      @test overrides["power_flow.max_iter"] == 80
      @test overrides["power_flow.autodamp"]
      @test overrides["benchmark.enabled"] === false

      invalid_form = copy(form)
      invalid_form["power_flow_max_iter"] = "invalid"
      invalid_response = Sparlectra.route_sparlectra_webui("POST", "/powerflow/run", invalid_form; output_root = output_root)
      invalid_html = String(invalid_response.body)
      @test invalid_response.status == 400
      @test occursin("name=\"casefile_manual\"", invalid_html)
      @test occursin("value=\"$(config_file)\" selected", invalid_html)

      form_html = Sparlectra.render_powerflow_form(output_root = output_root)
      expected_help_topics = Dict(
        "casefile" => "webui.casefile",
        "config_file" => "webui.config_file",
        "power_flow_tol" => "power_flow.tol",
        "power_flow_max_iter" => "power_flow.max_iter",
        "power_flow_autodamp" => "power_flow.autodamp",
        "power_flow_autodamp_min" => "power_flow.autodamp_min",
        "power_flow_qlimits_enabled" => "power_flow.qlimits.enabled",
        "power_flow_wrong_branch_detection" => "power_flow.wrong_branch_detection",
        "power_flow_start_angle_mode" => "power_flow.start_mode.angle_mode",
        "power_flow_start_voltage_mode" => "power_flow.start_mode.voltage_mode",
        "output_logfile_results" => "output.logfile_results",
        "benchmark_enabled" => "benchmark.enabled",
        "benchmark_samples" => "benchmark.samples",
        "benchmark_seconds" => "benchmark.seconds",
        "performance_timing" => "webui.performance_timing",
        "run_diagnostics" => "webui.run_diagnostics",
      )
      @test Sparlectra.WEBUI_FORM_HELP_TOPICS == expected_help_topics
      for (field, help_topic) in expected_help_topics
        @test occursin("name=\"$(field)\"", form_html)
        @test occursin("href=\"/help/$(help_topic)\"", form_html)
      end
      @test !occursin("name=\"output_root\"", form_html)
      @test occursin("<strong>Output root:</strong> <code>$(output_root)</code>", form_html)
      @test occursin("action=\"/webui/shutdown\"", form_html)
      @test occursin("Stop Web UI", form_html)
      @test count("class=\"help-link\"", form_html) == length(expected_help_topics)

      @test !occursin("class=\"button back-button\"", form_html)

      @test occursin("id=\"powerflow-run-form\"", form_html)
      @test occursin("onsubmit=\"this.classList.add('is-submitting')", form_html)
      @test occursin("this.setAttribute('aria-busy', 'true')", form_html)
      @test occursin("this.querySelector('button[type=submit]').disabled = true", form_html)
      @test occursin("window.addEventListener('pageshow'", form_html)
      @test occursin("form.classList.remove('is-submitting')", form_html)
      @test occursin("form.removeAttribute('aria-busy')", form_html)
      @test occursin("submitButton.disabled = false", form_html)
      @test occursin("class=\"powerflow-submit\"", form_html)
      @test occursin("class=\"submit-spinner\" aria-hidden=\"true\"", form_html)
      @test occursin("class=\"submit-progress-label\" role=\"status\" aria-live=\"polite\">Running PowerFlow…", form_html)

      @test occursin("src=\"/assets/logo.png\"", form_html)
      @test occursin("alt=\"Sparlectra.jl logo\"", form_html)
      @test occursin("Sparlectra.jl v", form_html)
      @test occursin("Sparlectra.jl v$(Sparlectra.version())", form_html)
      @test occursin("name=\"performance_timing\"", form_html)
      @test occursin("name=\"run_diagnostics\"", form_html)

      logo_response = Sparlectra.route_sparlectra_webui("GET", "/assets/logo.png")
      @test logo_response.status == 200
      @test ("Content-Type" => "image/png") in logo_response.headers
      @test logo_response.body == read(joinpath(@__DIR__, "..", "docs", "src", "assets", "logo.png"))
      for target in ("/assets/tablestyle.css", "/assets/../Project.toml", "/assets/logo.png/extra")
        @test Sparlectra.route_sparlectra_webui("GET", target).status == 404
      end
      routed_form = copy(form)
routed_form["output_root"] = joinpath(tmpdir, "outside-runs")
run_response = Sparlectra.route_sparlectra_webui("POST", "/powerflow/run", routed_form; output_root)
@test run_response.status == 303
result_location = only(header.second for header in run_response.headers if header.first == "Location")
run_id = basename(result_location)
result = get_powerflow_result(run_id)
@test result["success"]
@test result["output_dir"] == joinpath(abspath(output_root), run_id)
@test !ispath(joinpath(tmpdir, "outside-runs"))
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
      @test Set(("result.json", "run.log", "effective_config.yaml", "performance.log", "diagnose.txt")) ⊆ artifact_names
      artifact_response = Sparlectra.handle_powerflow_artifacts(run_id)
      artifact_html = String(artifact_response.body)
      for name in ("result.json", "run.log", "effective_config.yaml", "performance.log", "diagnose.txt")
        @test occursin(name, artifact_html)
      end
      result_artifact = Sparlectra.handle_powerflow_artifact(run_id, "result.json")
      @test result_artifact.status == 200
      result_artifact_html = String(result_artifact.body)
      @test occursin("Artifact: result.json", result_artifact_html)
      @test occursin("<main class=\"page artifact-page\">", result_artifact_html)
      @test occursin("class=\"artifact-text-page\"", result_artifact_html)
      @test occursin("class=\"button\" href=\"?download=1\"", result_artifact_html)
      @test occursin("class=\"artifact-text\"", result_artifact_html)
      css_text = read(joinpath(@__DIR__, "..", "src", "webui", "static", "sparlectra.css"), String)
      css_text = replace(css_text, "\r\n" => "\n")
      css_lines = split(css_text, '\n')
      @test length(css_lines) > 100
      @test maximum(length, css_lines) < 300
      @test occursin(".artifact-text {\n", css_text)
      @test occursin(".artifact-page {\n", css_text)
      @test occursin("  width: min(96vw, 1600px);\n", css_text)
      @test occursin("  width: 100%;\n", css_text)
      @test occursin("  white-space: pre;\n", css_text)
      @test occursin("  overflow: auto;\n", css_text)
      @test occursin("  min-height: 75vh;\n", css_text)
      @test occursin("  max-height: 85vh;\n", css_text)
      @test occursin(".form-grid.is-submitting .submit-spinner {\n", css_text)
      @test occursin("  animation: submit-spin .8s linear infinite;\n", css_text)
      @test occursin("@keyframes submit-spin {\n", css_text)
      @test occursin(".status-badge {\n", css_text)
      @test occursin(".status-success {\n", css_text)
      @test occursin(".status-warning {\n", css_text)
      @test occursin(".status-error {\n", css_text)
      @test occursin(".status-unknown {\n", css_text)
      @test occursin(".status-running {\n", css_text)
      @test occursin(".exit-button,\n", css_text)


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
      @test occursin("<th>Date/Time</th>", history_html)
      @test occursin("class=\"status-badge status-success\"", history_html)
      @test occursin("action=\"/powerflow/delete/$(run_id)\"", history_html)
      @test occursin("action=\"/powerflow/delete_all\"", history_html)
      @test occursin("action=\"/powerflow/refresh\"", history_html)
      @test occursin(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}", history_html)
      @test Sparlectra.route_sparlectra_webui("GET", "/powerflow/delete/$(run_id)"; output_root).status == 404
      unsafe_delete = Sparlectra.route_sparlectra_webui("POST", "/powerflow/delete/..%2Foutside"; output_root)
      @test unsafe_delete.status == 400
      @test occursin("Unsafe PowerFlow run ID rejected", String(unsafe_delete.body))
      @test Sparlectra.route_sparlectra_webui("GET", "/webui/shutdown"; output_root).status == 404
      @test Sparlectra.route_sparlectra_webui("POST", "/webui/shutdown"; output_root).status == 503
      @test Sparlectra.route_sparlectra_webui("POST", "/webui/heartbeat"; output_root).status == 204
      disposable_id = "disposable-run"
      disposable_dir = joinpath(output_root, disposable_id)
      mkpath(disposable_dir)
      disposable_result = joinpath(disposable_dir, "result.json")
      write(disposable_result, "{}")
      disposable_entries = load_powerflow_run_index(output_root)["runs"]
      push!(disposable_entries, Dict(
        "run_id" => disposable_id,
        "output_dir" => disposable_dir,
        "result_file" => disposable_result,
      ))
      Sparlectra._write_powerflow_run_entries!(output_root, disposable_entries)
      deleted_response = Sparlectra.route_sparlectra_webui("POST", "/powerflow/delete/$(disposable_id)"; output_root)
      @test deleted_response.status == 303
      @test !ispath(disposable_dir)
      @test all(entry -> get(entry, "run_id", "") != disposable_id, load_powerflow_run_index(output_root)["runs"])

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
        @test occursin("Sparlectra.jl v$(Sparlectra.version())", page_html)
      end

      source = join(read(joinpath(@__DIR__, "..", "src", "webui", file), String) for file in ("docs.jl", "forms.jl", "handlers.jl", "routes.jl", "views.jl", "webui.jl"))
      forbidden = ["run" * "pf!", "run" * "pf_rectangular!", "run_complex_nr_" * "rectangular"]
      @test all(name -> !occursin(name, source), forbidden)
      @test !occursin("Voltage-magnitude/angle blend strategy.", source)

      lock(Sparlectra._POWERFLOW_SERVICE_LOCK) do
        empty!(Sparlectra._POWERFLOW_SERVICE_RUNS)
      end
      probe = listen(ip"127.0.0.1", UInt16(0))
      port = Int(getsockname(probe)[2])
      close(probe)
      server = start_sparlectra_webui(port = port, output_root = output_root, auto_shutdown_on_browser_close = false)
      @test get_powerflow_result(run_id)["run_id"] == run_id
      try
        @test server.url == "http://127.0.0.1:$(port)/powerflow"
        response_text = _webui_http_request(port, "GET", "/powerflow")
        @test occursin("HTTP/1.1 200 OK", response_text)
        @test occursin("PowerFlow run", response_text)

        logo_response_text = _webui_http_request(port, "GET", "/assets/logo.png")
        @test occursin("HTTP/1.1 200 OK", logo_response_text)
        @test occursin("Content-Type: image/png", logo_response_text)

        heartbeat_response = _webui_http_request(port, "POST", "/webui/heartbeat")
        @test occursin("HTTP/1.1 204", heartbeat_response)
        shutdown_response = _webui_http_request(port, "POST", "/webui/shutdown")
        @test occursin("HTTP/1.1 200 OK", shutdown_response)
        @test occursin("Web UI stopped", shutdown_response)
        @test timedwait(() -> istaskdone(server.task), 2.0) == :ok
        @test !isopen(server.listener)
      finally
        close(server)
      end

      restarted = start_sparlectra_webui(port = port, output_root = output_root, auto_shutdown_on_browser_close = false)
      close(restarted)
      @test timedwait(() -> istaskdone(restarted.task), 2.0) == :ok

      interrupted = start_sparlectra_webui(port = port, output_root = output_root, auto_shutdown_on_browser_close = false)
      @test Sparlectra._wait_sparlectra_webui(interrupted; wait_for_task = _ -> throw(InterruptException()))
      @test timedwait(() -> istaskdone(interrupted.task), 2.0) == :ok
      @test !isopen(interrupted.listener)

      heartbeat_probe = listen(ip"127.0.0.1", UInt16(0))
      heartbeat_port = Int(getsockname(heartbeat_probe)[2])
      close(heartbeat_probe)
      heartbeat_server = start_sparlectra_webui(
        port = heartbeat_port,
        output_root = output_root,
        auto_shutdown_on_browser_close = true,
        browser_heartbeat_timeout_seconds = 0.2,
      )
      sleep(0.3)
      @test isopen(heartbeat_server.listener)
      _webui_http_request(heartbeat_port, "POST", "/webui/heartbeat")
      Sparlectra._webui_begin_request!(heartbeat_server.runtime)
      sleep(0.3)
      @test isopen(heartbeat_server.listener)
      Sparlectra._webui_finish_request!(heartbeat_server.runtime)
      @test timedwait(() -> istaskdone(heartbeat_server.task), 2.0) == :ok
      @test !isopen(heartbeat_server.listener)
    end
  end
  return nothing
end
