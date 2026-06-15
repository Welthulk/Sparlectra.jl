const _POWERFLOW_WEBUI_JOBS = Dict{String,Dict{String,Any}}()
const _POWERFLOW_WEBUI_ACTIVE_STATES = Set(("queued", "running", "aborting"))
const WEBUI_ABORT_HARD_RESET_AFTER_SECONDS = 60

struct PowerFlowAborted <: Exception end
Base.showerror(io::IO, ::PowerFlowAborted) = print(io, "PowerFlow run aborted by user.")

function _check_powerflow_cancelled!(token)
  token === nothing && return nothing
  token[] || return nothing
  throw(PowerFlowAborted())
end

function _update_webui_job_phase!(job::AbstractDict, phase::AbstractString)
  now = Dates.now(Dates.UTC)
  lock(_POWERFLOW_SERVICE_LOCK) do
    job["current_phase"] = String(phase)
    job["phase_started_at"] = now
    job["last_progress_at"] = now
  end
  return nothing
end

function _webui_phase_event!(job::AbstractDict, phase::AbstractString, event_callback)
  _update_webui_job_phase!(job, phase)
  event_callback("powerflow_phase_started"; run_id = job["run_id"], phase = String(phase), status = get(job, "status", "running"))
  return nothing
end

function _webui_job_snapshot(job::AbstractDict)::Dict{String,Any}
  snapshot = Dict{String,Any}(String(key) => value for (key, value) in job if key != "task" && key != "abort_requested")
  started_at = get(job, "started_at", nothing)
  finished_at = get(job, "finished_at", nothing)
  snapshot["elapsed_seconds"] = started_at === nothing ? nothing : max(0.0, Dates.value(something(finished_at, Dates.now(Dates.UTC)) - started_at) / 1000)
  abort_requested_at = get(job, "abort_requested_at", nothing)
  snapshot["abort_elapsed_seconds"] = abort_requested_at === nothing ? nothing : max(0.0, Dates.value(Dates.now(Dates.UTC) - abort_requested_at) / 1000)
  snapshot["hard_reset_available"] = get(snapshot, "status", "") == "aborting" && something(snapshot["abort_elapsed_seconds"], 0.0) >= WEBUI_ABORT_HARD_RESET_AFTER_SECONDS
  snapshot["success"] = get(snapshot, "status", "") == "success"
  return snapshot
end

function _webui_active_job()
  return lock(_POWERFLOW_SERVICE_LOCK) do
    for job in values(_POWERFLOW_WEBUI_JOBS)
      get(job, "status", "") in _POWERFLOW_WEBUI_ACTIVE_STATES && return job
    end
    return nothing
  end
end

function get_active_webui_powerflow_job()::Union{Nothing,Dict{String,Any}}
  job = _webui_active_job()
  return job === nothing ? nothing : lock(_POWERFLOW_SERVICE_LOCK) do
    _webui_job_snapshot(job)
  end
end

function _write_aborted_powerflow_result!(job::AbstractDict)
  output_dir = String(job["output_dir"])
  mkpath(output_dir)
  logfile = joinpath(output_dir, "run.log")
  open(logfile, "a") do io
    println(io, "Run aborted by user.")
    println(io, "Status: aborted")
    println(io, "Phase active at abort: ", get(job, "abort_phase", get(job, "current_phase", "unknown")))
  end
  result_file = joinpath(output_dir, "result.json")
  result = _api_result(
    run_id = String(job["run_id"]),
    status = :aborted,
    success = false,
    reason = "aborted_by_user",
    message = "Run aborted by user.",
    casefile = get(job, "casefile", nothing),
    config_file = get(job, "config_file", nothing),
    output_dir = output_dir,
    logfile = logfile,
    result_file = result_file,
  )
  _write_api_result_file(result)
  result = _refresh_api_artifacts(result)
  _write_api_result_file(result)
  _POWERFLOW_SERVICE_RUNS[result.run_id] = result
  _write_powerflow_run_index!(String(job["output_root"]), result)
  return result
end

function _write_webui_job_marker!(job::AbstractDict, status::Symbol, reason::String, message::String)
  output_dir = String(job["output_dir"])
  mkpath(output_dir)
  logfile = joinpath(output_dir, "run.log")
  touch(logfile)
  result_file = joinpath(output_dir, "result.json")
  result = _api_result(
    run_id = String(job["run_id"]),
    status = status,
    success = false,
    reason = reason,
    message = message,
    casefile = get(job, "casefile", nothing),
    config_file = get(job, "config_file", nothing),
    output_dir = output_dir,
    logfile = logfile,
    result_file = result_file,
  )
  _write_api_result_file(result)
  _POWERFLOW_SERVICE_RUNS[result.run_id] = result
  _write_powerflow_run_index!(String(job["output_root"]), result)
  return result
end

function hard_reset_webui_powerflow_run(run_id::AbstractString)::Dict{String,Any}
  _safe_powerflow_run_id(run_id) || return _service_failure("unsafe_run_id", "Unsafe PowerFlow run ID rejected."; run_id)
  reset = lock(_POWERFLOW_SERVICE_LOCK) do
    job = get(_POWERFLOW_WEBUI_JOBS, String(run_id), nothing)
    job === nothing && return _service_failure("run_not_found", "No Web UI PowerFlow run found for run_id $(run_id)."; run_id)
    get(job, "status", "") == "aborting" || return _service_failure("run_not_resettable", "Hard reset is only available while abort is pending."; run_id)
    previous_phase = get(job, "current_phase", "unknown")
    output_dir = String(job["output_dir"])
    message = "Hard reset requested while aborting. This result is not a valid solved PowerFlow result."
    job["status"] = "aborted_unknown"
    job["current_phase"] = "hard_reset_requested"
    job["message"] = message
    job["finished_at"] = Dates.now(Dates.UTC)
    (success = true, job = job, previous_phase = previous_phase, output_dir = output_dir, message = message)
  end
  haskey(reset, :success) || return reset
  mkpath(reset.output_dir)
  open(joinpath(reset.output_dir, "run.log"), "a") do io
      println(io, "Hard reset requested by user while aborting.")
      println(io, "Previous phase: ", reset.previous_phase)
      println(io, "Result is not a valid solved PowerFlow result.")
  end
  _write_webui_job_marker!(reset.job, :aborted_unknown, "hard_reset_requested", reset.message)
  return merge(_webui_job_snapshot(reset.job), Dict("previous_phase" => reset.previous_phase))
end

"""
    start_webui_powerflow_run(request; case_directory=nothing, runner=start_powerflow_run)

Start one Web UI PowerFlow job in a background task. Only one active Web UI job
is accepted at a time. Cancellation is cooperative: an abort request releases
the UI immediately and the worker discards any later solver success.
"""
function start_webui_powerflow_run(request::AbstractDict; case_directory::Union{Nothing,AbstractString} = nothing, runner = start_powerflow_run, event_callback = (event; fields...) -> nothing)::Dict{String,Any}
  active = _webui_active_job()
  active === nothing || return _service_failure("active_run", "A PowerFlow run is already active. Abort it or wait for it to finish."; run_id = active["run_id"])
  output_root = _service_request_value(request, "output_root")
  output_root isa AbstractString && !isempty(strip(output_root)) || return _service_failure("invalid_request", "PowerFlow service request requires a nonempty output_root.")
  run_id = string(uuid4())
  job = Dict{String,Any}(
    "run_id" => run_id,
    "status" => "queued",
    "casefile" => _service_request_value(request, "casefile"),
    "resolved_casefile" => nothing,
    "config_file" => _service_request_value(request, "config_file"),
    "output_root" => abspath(output_root),
    "output_dir" => joinpath(abspath(output_root), run_id),
    "started_at" => Dates.now(Dates.UTC),
    "finished_at" => nothing,
    "message" => "PowerFlow run queued.",
    "abort_requested" => Threads.Atomic{Bool}(false),
    "abort_requested_at" => nothing,
    "current_phase" => "queued",
    "phase_started_at" => Dates.now(Dates.UTC),
    "last_progress_at" => Dates.now(Dates.UTC),
  )
  lock(_POWERFLOW_SERVICE_LOCK) do
    _POWERFLOW_WEBUI_JOBS[run_id] = job
    _write_webui_job_marker!(job, :queued, "webui_job_active", "PowerFlow run queued.")
  end
  event_callback("powerflow_submitted"; run_id, requested_case = job["casefile"], status = "accepted")
  job["task"] = Threads.@spawn begin
    try
      lock(_POWERFLOW_SERVICE_LOCK) do
        job["status"] = "running"
        job["message"] = "PowerFlow run is active."
      end
      _webui_phase_event!(job, "resolving_case", event_callback)
      event_callback("powerflow_started"; run_id, requested_case = job["casefile"], status = "running")
      detailed_result_csv = _service_request_value(request, "detailed_result_csv", false)
      detailed_result_csv_semicolon = _service_request_value(request, "detailed_result_csv_semicolon", false)
      detailed_result_csv_format = _service_request_value(request, "detailed_result_csv_format", nothing)
      if detailed_result_csv
        csv_format = _resolve_detailed_csv_format(detailed_result_csv_format === nothing ? (detailed_result_csv_semicolon ? "excel_de" : "technical") : detailed_result_csv_format)
        event_callback("detailed_result_csv_export_enabled"; run_id, csv_format = csv_format.name, delimiter = string(csv_format.delimiter), decimal_separator = csv_format.decimal_separator, thousands_separator = csv_format.thousands_separator, status = "enabled")
      end
      worker_request = Dict{String,Any}(String(key) => value for (key, value) in request)
      worker_request["run_id"] = run_id
      worker_request["cancellation_token"] = job["abort_requested"]
      worker_request["phase_callback"] = phase -> _webui_phase_event!(job, phase, event_callback)
      result = try
        runner(worker_request; case_directory)
      catch err
        err isa PowerFlowAborted ? Dict{String,Any}("status" => "aborted", "success" => false) :
          _service_failure("execution_error", sprint(showerror, err, catch_backtrace()); run_id)
      end
      hard_reset_requested = lock(_POWERFLOW_SERVICE_LOCK) do
        get(job, "status", "") == "aborted_unknown"
      end
      if hard_reset_requested
        nothing
      elseif job["abort_requested"][]
        _webui_phase_event!(job, "finalizing_aborted", event_callback)
        _write_aborted_powerflow_result!(job)
        event_callback("powerflow_aborted"; run_id = job["run_id"], requested_case = job["casefile"], status = "aborted", current_phase = "finalizing_aborted")
        lock(_POWERFLOW_SERVICE_LOCK) do
          job["status"] = "aborted"
          job["current_phase"] = "aborted"
          job["message"] = "Run aborted by user."
        end
      else
        _webui_phase_event!(job, "finalizing_success", event_callback)
        lock(_POWERFLOW_SERVICE_LOCK) do
          if get(job, "status", "") == "aborted_unknown"
            return
          end
          job["status"] = get(result, "success", false) ? "success" : "failed"
          job["current_phase"] = job["status"]
          job["message"] = something(get(result, "message", nothing), job["status"] == "success" ? "PowerFlow run completed." : "PowerFlow run failed.")
          job["resolved_casefile"] = get(result, "casefile", nothing)
          if haskey(result, "run_id") && result["run_id"] != run_id
            delete!(_POWERFLOW_WEBUI_JOBS, run_id)
            job["run_id"] = result["run_id"]
            job["output_dir"] = get(result, "output_dir", job["output_dir"])
            _POWERFLOW_WEBUI_JOBS[result["run_id"]] = job
          end
        end
        event_callback(
          job["status"] == "success" ? "powerflow_completed" : "powerflow_failed";
          run_id = job["run_id"],
          requested_case = job["casefile"],
          resolved_case = job["resolved_casefile"],
          status = job["status"],
          message = job["message"],
        )
        if detailed_result_csv
          artifact_names = Set(String(get(artifact, "name", "")) for artifact in get(result, "artifacts", Any[]))
          expected_csv = ["bus_voltages_complex.csv", "branch_flows.csv"]
          if all(name -> name in artifact_names, expected_csv)
            event_callback("detailed_result_csv_exported"; run_id = job["run_id"], csv_format = csv_format.name, artifacts = expected_csv, status = "succeeded")
          else
            event_callback("detailed_result_csv_export_failed"; run_id = job["run_id"], message = "Detailed result CSV artifacts were not generated.", status = "failed")
          end
        end
      end
    finally
      lock(_POWERFLOW_SERVICE_LOCK) do
        job["finished_at"] = Dates.now(Dates.UTC)
        if get(job, "status", "") in _POWERFLOW_WEBUI_ACTIVE_STATES
          job["status"] = job["abort_requested"][] ? "aborted" : "failed"
          job["message"] = job["abort_requested"][] ? "Run aborted by user." : "PowerFlow worker exited before reaching a terminal result."
        end
      end
    end
  end
  return _webui_job_snapshot(job)
end

function get_webui_powerflow_job(run_id::AbstractString)::Dict{String,Any}
  _safe_powerflow_run_id(run_id) || return _service_failure("unsafe_run_id", "Unsafe PowerFlow run ID rejected."; run_id)
  return lock(_POWERFLOW_SERVICE_LOCK) do
    job = get(_POWERFLOW_WEBUI_JOBS, String(run_id), nothing)
    job === nothing ? get_powerflow_result(run_id) : _webui_job_snapshot(job)
  end
end

function abort_webui_powerflow_run(run_id::AbstractString)::Dict{String,Any}
  _safe_powerflow_run_id(run_id) || return _service_failure("unsafe_run_id", "Unsafe PowerFlow run ID rejected."; run_id)
  return lock(_POWERFLOW_SERVICE_LOCK) do
    job = get(_POWERFLOW_WEBUI_JOBS, String(run_id), nothing)
    job === nothing && return _service_failure("run_not_found", "No active Web UI PowerFlow run found for run_id $(run_id)."; run_id)
    status = get(job, "status", "")
    status == "aborting" && return merge(_webui_job_snapshot(job), Dict("abort_status" => "already_aborting"))
    status == "aborted" && return merge(_webui_job_snapshot(job), Dict("abort_status" => "already_aborted"))
    status in _POWERFLOW_WEBUI_ACTIVE_STATES || return _service_failure("run_not_abortable", "PowerFlow run $(run_id) is no longer active."; run_id)
    job["abort_requested"][] = true
    job["abort_requested_at"] = Dates.now(Dates.UTC)
    job["abort_phase"] = get(job, "current_phase", "unknown")
    job["status"] = "aborting"
    job["message"] = "Abort requested. Current phase: $(job["abort_phase"]). This phase may need to finish before cancellation is observed."
    merge(_webui_job_snapshot(job), Dict("abort_status" => "accepted"))
  end
end
