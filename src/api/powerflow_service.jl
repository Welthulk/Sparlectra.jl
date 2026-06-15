const POWERFLOW_RUN_INDEX_FILENAME = "powerflow_runs_index.json"
const _POWERFLOW_SERVICE_RUNS = Dict{String,SparlectraApiResult}()
const _POWERFLOW_SERVICE_LOCK = ReentrantLock()
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

mutable struct _ServiceJsonParser
  text::String
  index::Int
end

function _service_failure(reason::AbstractString, message::AbstractString; run_id = nothing)::Dict{String,Any}
  failure = Dict{String,Any}(
    "status" => "failed",
    "success" => false,
    "reason" => String(reason),
    "message" => String(message),
  )
  run_id === nothing || (failure["run_id"] = String(run_id))
  return failure
end

function _service_request_value(request::AbstractDict, key::String, default = nothing)
  haskey(request, key) && return request[key]
  symbol_key = Symbol(key)
  haskey(request, symbol_key) && return request[symbol_key]
  return default
end

function _skip_service_json_whitespace!(parser::_ServiceJsonParser)
  while parser.index <= lastindex(parser.text) && isspace(parser.text[parser.index])
    parser.index = nextind(parser.text, parser.index)
  end
  return nothing
end

function _parse_service_json_string!(parser::_ServiceJsonParser)::String
  parser.text[parser.index] == '"' || error("Expected JSON string")
  parser.index = nextind(parser.text, parser.index)
  io = IOBuffer()
  while parser.index <= lastindex(parser.text)
    char = parser.text[parser.index]
    parser.index = nextind(parser.text, parser.index)
    char == '"' && return String(take!(io))
    if char == '\\'
      parser.index <= lastindex(parser.text) || error("Incomplete JSON escape")
      escaped = parser.text[parser.index]
      parser.index = nextind(parser.text, parser.index)
      escaped == '"' ? print(io, '"') :
        escaped == '\\' ? print(io, '\\') :
        escaped == '/' ? print(io, '/') :
        escaped == 'b' ? print(io, '\b') :
        escaped == 'f' ? print(io, '\f') :
        escaped == 'n' ? print(io, '\n') :
        escaped == 'r' ? print(io, '\r') :
        escaped == 't' ? print(io, '\t') :
        escaped == 'u' ? _parse_service_json_unicode_escape!(parser, io) :
        error("Unsupported JSON escape \\$(escaped)")
    else
      print(io, char)
    end
  end
  error("Unterminated JSON string")
end

function _parse_service_json_unicode_escape!(parser::_ServiceJsonParser, io::IO)
  digits = IOBuffer()
  for _ in 1:4
    parser.index <= lastindex(parser.text) || error("Incomplete JSON Unicode escape")
    char = parser.text[parser.index]
    isxdigit(char) || error("Invalid JSON Unicode escape")
    print(digits, char)
    parser.index = nextind(parser.text, parser.index)
  end
  print(io, Char(parse(UInt32, String(take!(digits)); base = 16)))
  return nothing
end

function _parse_service_json_number!(parser::_ServiceJsonParser)
  start = parser.index
  while parser.index <= lastindex(parser.text) && parser.text[parser.index] in ('-', '+', '.', 'e', 'E', '0':'9'...)
    parser.index = nextind(parser.text, parser.index)
  end
  stop = prevind(parser.text, parser.index)
  token = parser.text[start:stop]
  return occursin(r"[.eE]", token) ? parse(Float64, token) : parse(Int, token)
end

function _parse_service_json_literal!(parser::_ServiceJsonParser, literal::String, value)
  stop = nextind(parser.text, parser.index, length(literal))
  parser.text[parser.index:prevind(parser.text, stop)] == literal || error("Invalid JSON literal")
  parser.index = stop
  return value
end

function _parse_service_json_array!(parser::_ServiceJsonParser)::Vector{Any}
  parser.index = nextind(parser.text, parser.index)
  values = Any[]
  _skip_service_json_whitespace!(parser)
  if parser.index <= lastindex(parser.text) && parser.text[parser.index] == ']'
    parser.index = nextind(parser.text, parser.index)
    return values
  end
  while true
    push!(values, _parse_service_json_value!(parser))
    _skip_service_json_whitespace!(parser)
    parser.index <= lastindex(parser.text) || error("Unterminated JSON array")
    separator = parser.text[parser.index]
    parser.index = nextind(parser.text, parser.index)
    separator == ']' && return values
    separator == ',' || error("Expected comma in JSON array")
  end
end

function _parse_service_json_object!(parser::_ServiceJsonParser)::Dict{String,Any}
  parser.index = nextind(parser.text, parser.index)
  values = Dict{String,Any}()
  _skip_service_json_whitespace!(parser)
  if parser.index <= lastindex(parser.text) && parser.text[parser.index] == '}'
    parser.index = nextind(parser.text, parser.index)
    return values
  end
  while true
    _skip_service_json_whitespace!(parser)
    key = _parse_service_json_string!(parser)
    _skip_service_json_whitespace!(parser)
    parser.index <= lastindex(parser.text) && parser.text[parser.index] == ':' || error("Expected colon in JSON object")
    parser.index = nextind(parser.text, parser.index)
    values[key] = _parse_service_json_value!(parser)
    _skip_service_json_whitespace!(parser)
    parser.index <= lastindex(parser.text) || error("Unterminated JSON object")
    separator = parser.text[parser.index]
    parser.index = nextind(parser.text, parser.index)
    separator == '}' && return values
    separator == ',' || error("Expected comma in JSON object")
  end
end

function _parse_service_json_value!(parser::_ServiceJsonParser)
  _skip_service_json_whitespace!(parser)
  parser.index <= lastindex(parser.text) || error("Unexpected end of JSON input")
  char = parser.text[parser.index]
  char == '{' && return _parse_service_json_object!(parser)
  char == '[' && return _parse_service_json_array!(parser)
  char == '"' && return _parse_service_json_string!(parser)
  char == 't' && return _parse_service_json_literal!(parser, "true", true)
  char == 'f' && return _parse_service_json_literal!(parser, "false", false)
  char == 'n' && return _parse_service_json_literal!(parser, "null", nothing)
  (char == '-' || isdigit(char)) && return _parse_service_json_number!(parser)
  error("Unsupported JSON value")
end

function _parse_service_json(text::AbstractString)
  source = String(text)
  isempty(strip(source)) && error("Empty JSON input")
  parser = _ServiceJsonParser(source, firstindex(source))
  value = _parse_service_json_value!(parser)
  _skip_service_json_whitespace!(parser)
  parser.index > lastindex(source) || error("Unexpected trailing JSON content")
  return value
end

function _registered_powerflow_run(run_id::AbstractString)
  id = String(run_id)
  return lock(_POWERFLOW_SERVICE_LOCK) do
    get(_POWERFLOW_SERVICE_RUNS, id, nothing)
  end
end

function _register_powerflow_run!(result::SparlectraApiResult)
  lock(_POWERFLOW_SERVICE_LOCK) do
    _POWERFLOW_SERVICE_RUNS[result.run_id] = result
  end
  return result
end

function _safe_powerflow_run_id(run_id::AbstractString)::Bool
  id = String(run_id)
  return !isempty(id) && !_unsafe_artifact_name(id) && occursin(r"^[A-Za-z0-9][A-Za-z0-9._-]*$", id)
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
  return lock(_POWERFLOW_SERVICE_LOCK) do
    job = get(_POWERFLOW_WEBUI_JOBS, String(run_id), nothing)
    job === nothing && return _service_failure("run_not_found", "No Web UI PowerFlow run found for run_id $(run_id)."; run_id)
    get(job, "status", "") == "aborting" || return _service_failure("run_not_resettable", "Hard reset is only available while abort is pending."; run_id)
    previous_phase = get(job, "current_phase", "unknown")
    job["status"] = "aborted_unknown"
    job["current_phase"] = "hard_reset_requested"
    job["message"] = "Hard reset requested while aborting. This result is not a valid solved PowerFlow result."
    job["finished_at"] = Dates.now(Dates.UTC)
    mkpath(String(job["output_dir"]))
    open(joinpath(String(job["output_dir"]), "run.log"), "a") do io
      println(io, "Hard reset requested by user while aborting.")
      println(io, "Previous phase: ", previous_phase)
      println(io, "Result is not a valid solved PowerFlow result.")
    end
    _write_webui_job_marker!(job, :aborted_unknown, "hard_reset_requested", job["message"])
    merge(_webui_job_snapshot(job), Dict("previous_phase" => previous_phase))
  end
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
      _update_webui_job_phase!(job, "resolving_case")
      event_callback("powerflow_started"; run_id, requested_case = job["casefile"], status = "running")
      detailed_result_csv = _service_request_value(request, "detailed_result_csv", false)
      detailed_result_csv_semicolon = _service_request_value(request, "detailed_result_csv_semicolon", false)
      detailed_result_csv && event_callback("detailed_result_csv_export_enabled"; run_id, delimiter = detailed_result_csv_semicolon ? "semicolon" : "comma", status = "enabled")
      worker_request = Dict{String,Any}(String(key) => value for (key, value) in request)
      worker_request["run_id"] = run_id
      worker_request["cancellation_token"] = job["abort_requested"]
      worker_request["phase_callback"] = phase -> _update_webui_job_phase!(job, phase)
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
        _update_webui_job_phase!(job, "finalizing_aborted")
        _write_aborted_powerflow_result!(job)
        event_callback("powerflow_aborted"; run_id = job["run_id"], requested_case = job["casefile"], status = "aborted", current_phase = "finalizing_aborted")
        lock(_POWERFLOW_SERVICE_LOCK) do
          job["status"] = "aborted"
          job["current_phase"] = "aborted"
          job["message"] = "Run aborted by user."
        end
      else
        _update_webui_job_phase!(job, "finalizing_success")
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
            event_callback("detailed_result_csv_exported"; run_id = job["run_id"], artifacts = expected_csv, status = "succeeded")
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
    job["status"] = "aborting"
    job["message"] = "Abort requested. Waiting for the calculation to stop at the next safe cancellation point."
    merge(_webui_job_snapshot(job), Dict("abort_status" => "accepted"))
  end
end

function _powerflow_index_path(output_root::AbstractString)::String
  return joinpath(abspath(output_root), POWERFLOW_RUN_INDEX_FILENAME)
end

function _empty_powerflow_run_index()::Dict{String,Any}
  return Dict{String,Any}("schema_version" => _SPARLECTRA_API_SCHEMA_VERSION, "runs" => Any[])
end

function _path_is_within(path::AbstractString, root::AbstractString)::Bool
  normalized_path = normpath(abspath(path))
  normalized_root = normpath(abspath(root))
  relative = relpath(normalized_path, normalized_root)
  relative == "." && return true
  return first(splitpath(relative)) != ".."
end

function _existing_path_is_within(path::AbstractString, root::AbstractString)::Bool
  _path_is_within(path, root) || return false
  ispath(path) || return true
  isdir(root) || return false
  return _path_is_within(realpath(path), realpath(root))
end

function _indexed_run_paths(entry::AbstractDict, output_root::AbstractString)
  run_id = get(entry, "run_id", nothing)
  output_dir = get(entry, "output_dir", nothing)
  result_file = get(entry, "result_file", nothing)
  run_id isa AbstractString && !_unsafe_artifact_name(String(run_id)) || return (valid = false, reason = "unsafe_run_id", run_id = run_id, output_dir = output_dir, result_file = result_file)
  output_dir isa AbstractString || return (valid = false, reason = "invalid_output_dir", run_id = run_id, output_dir = output_dir, result_file = result_file)
  result_file isa AbstractString || return (valid = false, reason = "invalid_result_file", run_id = run_id, output_dir = output_dir, result_file = result_file)

  root = normpath(abspath(output_root))
  directory = normpath(abspath(output_dir))
  result_path = normpath(abspath(result_file))
  expected_directory = joinpath(root, String(run_id))
  directory == expected_directory || return (valid = false, reason = "unsafe_output_dir", run_id = run_id, output_dir = directory, result_file = result_path)
  _existing_path_is_within(directory, root) || return (valid = false, reason = "unsafe_output_dir", run_id = run_id, output_dir = directory, result_file = result_path)
  result_path == joinpath(directory, "result.json") || return (valid = false, reason = "unsafe_result_file", run_id = run_id, output_dir = directory, result_file = result_path)
  _existing_path_is_within(result_path, directory) || return (valid = false, reason = "unsafe_result_file", run_id = run_id, output_dir = directory, result_file = result_path)
  return (valid = true, reason = nothing, run_id = String(run_id), output_dir = directory, result_file = result_path)
end

function _powerflow_run_index_entry(result::SparlectraApiResult)::Dict{String,Any}
  return Dict{String,Any}(
    "run_id" => result.run_id,
    "schema_version" => result.schema_version,
    "status" => String(result.status),
    "success" => result.success,
    "output_dir" => result.output_dir,
    "result_file" => result.result_file,
    "casefile" => result.casefile,
    "config_file" => result.config_file,
    "final_mismatch" => result.final_mismatch,
    "iterations" => result.iterations,
    "timestamp" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
  )
end

"""
    load_powerflow_run_index(output_root::AbstractString) -> Dict{String,Any}

Load the persistent PowerFlow run index beneath `output_root`. A missing or
unreadable index returns an empty transport-safe index. Invalid top-level
content is reported through `reason` and `message` without throwing.
"""
function load_powerflow_run_index(output_root::AbstractString)::Dict{String,Any}
  index_path = _powerflow_index_path(output_root)
  isfile(index_path) || return _empty_powerflow_run_index()
  try
    parsed = _parse_service_json(read(index_path, String))
    parsed isa AbstractDict || error("PowerFlow run index must contain a JSON object")
    runs = get(parsed, "runs", Any[])
    runs isa AbstractVector || error("PowerFlow run index field 'runs' must be an array")
    return Dict{String,Any}(
      "schema_version" => string(get(parsed, "schema_version", _SPARLECTRA_API_SCHEMA_VERSION)),
      "runs" => Any[run for run in runs],
    )
  catch err
    index = _empty_powerflow_run_index()
    index["reason"] = "invalid_run_index"
    index["message"] = sprint(showerror, err)
    return index
  end
end

function _validated_powerflow_run_entry(entry, output_root::AbstractString)::Dict{String,Any}
  entry isa AbstractDict || return Dict{String,Any}("available" => false, "reason" => "invalid_run_entry")
  listed = Dict{String,Any}(String(key) => _api_transport_value(value) for (key, value) in entry)
  paths = _indexed_run_paths(entry, output_root)
  if !haskey(listed, "timestamp") && paths.valid && isfile(paths.result_file)
    listed["timestamp"] = Dates.format(Dates.unix2datetime(stat(paths.result_file).mtime), "yyyy-mm-dd HH:MM:SS")
  end
  listed["available"] = paths.valid && isdir(paths.output_dir) && isfile(paths.result_file)
  if !paths.valid
    listed["reason"] = paths.reason
  elseif !isdir(paths.output_dir)
    listed["reason"] = "output_dir_not_found"
  elseif !isfile(paths.result_file)
    listed["reason"] = "result_file_not_found"
  end
  return listed
end

"""
    list_powerflow_runs(output_root::AbstractString) -> Vector{Dict{String,Any}}

List runs from the persistent index for a future run-history UI. Each entry has
an `available` flag; missing directories, missing `result.json` files, and
unsafe indexed paths are described by a structured `reason` instead of raising.
"""
function list_powerflow_runs(output_root::AbstractString)::Vector{Dict{String,Any}}
  index = load_powerflow_run_index(output_root)
  runs = [_validated_powerflow_run_entry(entry, output_root) for entry in index["runs"]]
  root = abspath(output_root)
  active_jobs = lock(_POWERFLOW_SERVICE_LOCK) do
    [_webui_job_snapshot(job) for job in values(_POWERFLOW_WEBUI_JOBS) if job["output_root"] == root && get(job, "status", "") in _POWERFLOW_WEBUI_ACTIVE_STATES]
  end
  indexed_ids = Set(string(get(run, "run_id", "")) for run in runs)
  for job in active_jobs
    job["available"] = true
    job["timestamp"] = Dates.format(job["started_at"], "yyyy-mm-dd HH:MM:SS")
    job["run_id"] in indexed_ids || push!(runs, job)
  end
  return runs
end

function _write_powerflow_run_entries!(output_root::AbstractString, runs::AbstractVector)
  root = abspath(output_root)
  mkpath(root)
  contents = Dict{String,Any}("schema_version" => _SPARLECTRA_API_SCHEMA_VERSION, "runs" => runs)
  index_path = _powerflow_index_path(root)
  temporary_path = index_path * ".tmp"
  open(temporary_path, "w") do io
    _write_json(io, contents)
    println(io)
  end
  mv(temporary_path, index_path; force = true)
  return index_path
end

function _write_powerflow_run_index!(output_root::AbstractString, result::SparlectraApiResult)
  index = load_powerflow_run_index(output_root)
  runs = Any[entry for entry in index["runs"] if !(entry isa AbstractDict && get(entry, "run_id", nothing) == result.run_id)]
  push!(runs, _powerflow_run_index_entry(result))
  sort!(runs; by = entry -> entry isa AbstractDict ? string(get(entry, "run_id", "")) : "")
  return _write_powerflow_run_entries!(output_root, runs)
end

function _optional_string(value)
  value === nothing && return nothing
  value isa AbstractString || error("Expected a string or null")
  return String(value)
end

function _optional_bool(value)
  value === nothing && return nothing
  value isa Bool || error("Expected a boolean or null")
  return value
end

function _optional_int(value)
  value === nothing && return nothing
  value isa Integer || error("Expected an integer or null")
  return Int(value)
end

function _optional_float(value)
  value === nothing && return nothing
  value isa Real || error("Expected a number or null")
  return Float64(value)
end

function _reconstruct_powerflow_result(data::AbstractDict, paths)::SparlectraApiResult
  String(get(data, "run_id", "")) == paths.run_id || error("result.json run_id does not match its index entry")
  status = get(data, "status", nothing)
  status isa AbstractString || error("result.json status must be a string")
  success = get(data, "success", nothing)
  success isa Bool || error("result.json success must be a boolean")
  solution_available = get(data, "solution_available", false)
  solution_available isa Bool || error("result.json solution_available must be a boolean")
  return _api_result(
    run_id = paths.run_id,
    schema_version = string(get(data, "schema_version", _SPARLECTRA_API_SCHEMA_VERSION)),
    status = Symbol(status),
    success = success,
    converged = _optional_bool(get(data, "converged", nothing)),
    solution_available = solution_available,
    iterations = _optional_int(get(data, "iterations", nothing)),
    final_mismatch = _optional_float(get(data, "final_mismatch", nothing)),
    reason = _optional_string(get(data, "reason", nothing)),
    message = _optional_string(get(data, "message", nothing)),
    casefile = _optional_string(get(data, "casefile", nothing)),
    config_file = _optional_string(get(data, "config_file", nothing)),
    output_dir = paths.output_dir,
    logfile = _optional_string(get(data, "logfile", nothing)),
    result_file = paths.result_file,
    artifacts = collect_sparlectra_api_artifacts(paths.output_dir),
    raw_result = nothing,
  )
end

"""
    refresh_powerflow_run_registry!(output_root::AbstractString) -> Dict{String,Any}

Refresh the in-process PowerFlow registry from the persistent index. Existing
entries belonging to `output_root` are replaced. Missing, corrupt, or unsafe
runs are skipped and returned in `unavailable_runs`; valid runs remain
available through result and artifact lookup functions.
"""
function refresh_powerflow_run_registry!(output_root::AbstractString)::Dict{String,Any}
  root = abspath(output_root)
  index = load_powerflow_run_index(root)
  recovered = SparlectraApiResult[]
  unavailable = Dict{String,Any}[]
  if haskey(index, "reason")
    push!(unavailable, Dict{String,Any}(
      "reason" => index["reason"],
      "message" => get(index, "message", nothing),
    ))
  end
  for entry in index["runs"]
    if !(entry isa AbstractDict)
      push!(unavailable, Dict{String,Any}("reason" => "invalid_run_entry"))
      continue
    end
    paths = _indexed_run_paths(entry, root)
    if !paths.valid
      push!(unavailable, Dict{String,Any}("run_id" => _api_transport_value(paths.run_id), "reason" => paths.reason))
      continue
    end
    if !isfile(paths.result_file)
      push!(unavailable, Dict{String,Any}("run_id" => paths.run_id, "reason" => "result_file_not_found"))
      continue
    end
    try
      data = _parse_service_json(read(paths.result_file, String))
      data isa AbstractDict || error("result.json must contain a JSON object")
      if lowercase(string(get(data, "status", ""))) in _POWERFLOW_WEBUI_ACTIVE_STATES
        data["status"] = "aborted_unknown"
        data["success"] = false
        data["solution_available"] = false
        data["reason"] = "webui_stale_active_run"
        data["message"] = "The Web UI was restarted while this run was active. This result is not a valid solved PowerFlow result."
        open(paths.result_file, "w") do io
          _write_json(io, data)
          println(io)
        end
      end
      push!(recovered, _reconstruct_powerflow_result(data, paths))
    catch err
      push!(unavailable, Dict{String,Any}(
        "run_id" => paths.run_id,
        "reason" => "invalid_result_file",
        "message" => sprint(showerror, err),
      ))
    end
  end

  lock(_POWERFLOW_SERVICE_LOCK) do
    for (run_id, result) in collect(_POWERFLOW_SERVICE_RUNS)
      _path_is_within(result.output_dir, root) && delete!(_POWERFLOW_SERVICE_RUNS, run_id)
    end
    for result in recovered
      _POWERFLOW_SERVICE_RUNS[result.run_id] = result
    end
  end
  return Dict{String,Any}(
    "status" => isempty(unavailable) ? "succeeded" : "partial",
    "success" => isempty(unavailable),
    "loaded_runs" => [result.run_id for result in recovered],
    "runs" => recovered,
    "unavailable_runs" => unavailable,
  )
end

"""
    delete_powerflow_run(run_id::AbstractString; output_root::AbstractString) -> Dict{String,Any}

Delete one registered PowerFlow run beneath `output_root`. The run ID must be a
safe index entry whose directory is exactly `<output_root>/<run_id>`; arbitrary
paths and unregistered directories are never removed.
"""
function delete_powerflow_run(run_id::AbstractString; output_root::AbstractString)::Dict{String,Any}
  id = String(run_id)
  _unsafe_artifact_name(id) && return _service_failure("unsafe_run_id", "Unsafe PowerFlow run ID rejected."; run_id = id)
  root = abspath(output_root)
  index = load_powerflow_run_index(root)
  entry_index = findfirst(entry -> entry isa AbstractDict && get(entry, "run_id", nothing) == id, index["runs"])
  entry_index === nothing && return _service_failure("run_not_found", "PowerFlow run is not registered for this output root."; run_id = id)
  entry = index["runs"][entry_index]
  paths = _indexed_run_paths(entry, root)
  paths.valid || return _service_failure(paths.reason, "PowerFlow run has an unsafe indexed path and was not deleted."; run_id = id)

  try
    ispath(paths.output_dir) && rm(paths.output_dir; recursive = true)
  catch err
    return _service_failure("delete_failed", sprint(showerror, err); run_id = id)
  end
  remaining = Any[entry for (entry_position, entry) in enumerate(index["runs"]) if entry_position != entry_index]
  _write_powerflow_run_entries!(root, remaining)
  lock(_POWERFLOW_SERVICE_LOCK) do
    result = get(_POWERFLOW_SERVICE_RUNS, id, nothing)
    result === nothing || (_path_is_within(result.output_dir, root) && delete!(_POWERFLOW_SERVICE_RUNS, id))
  end
  return Dict{String,Any}("status" => "succeeded", "success" => true, "run_id" => id)
end

"""
    delete_all_powerflow_runs(; output_root::AbstractString) -> Dict{String,Any}

Delete all safely registered PowerFlow run directories beneath `output_root`.
Entries that cannot be validated or removed remain indexed and are reported in
`failed_runs`; deletion never follows an indexed path outside the output root.
"""
function delete_all_powerflow_runs(; output_root::AbstractString)::Dict{String,Any}
  root = abspath(output_root)
  index = load_powerflow_run_index(root)
  deleted = String[]
  failed = Dict{String,Any}[]
  remaining = Any[]
  for entry in index["runs"]
    if !(entry isa AbstractDict)
      push!(remaining, entry)
      push!(failed, Dict{String,Any}("reason" => "invalid_run_entry"))
      continue
    end
    paths = _indexed_run_paths(entry, root)
    if !paths.valid
      push!(remaining, entry)
      push!(failed, Dict{String,Any}("run_id" => _api_transport_value(paths.run_id), "reason" => paths.reason))
      continue
    end
    try
      ispath(paths.output_dir) && rm(paths.output_dir; recursive = true)
      push!(deleted, paths.run_id)
    catch err
      push!(remaining, entry)
      push!(failed, Dict{String,Any}("run_id" => paths.run_id, "reason" => "delete_failed", "message" => sprint(showerror, err)))
    end
  end
  _write_powerflow_run_entries!(root, remaining)
  lock(_POWERFLOW_SERVICE_LOCK) do
    for run_id in deleted
      result = get(_POWERFLOW_SERVICE_RUNS, run_id, nothing)
      result === nothing || (_path_is_within(result.output_dir, root) && delete!(_POWERFLOW_SERVICE_RUNS, run_id))
    end
  end
  return Dict{String,Any}(
    "status" => isempty(failed) ? "succeeded" : "partial",
    "success" => isempty(failed),
    "deleted_runs" => deleted,
    "failed_runs" => failed,
  )
end

function _prefer_julia_casefile(casefile::AbstractString; emit_julia_case_fn = FetchMatpowerCase.emit_julia_case)::String
  case_path = abspath(casefile)
  extension = lowercase(splitext(case_path)[2])
  extension == ".jl" && return case_path
  extension == ".m" || throw(ArgumentError("Unsupported casefile extension: $(casefile) (expected .m or .jl)"))

  julia_path = first(splitext(case_path)) * ".jl"
  isfile(julia_path) && return julia_path
  try
    emitted_path = abspath(emit_julia_case_fn(case_path, dirname(case_path)))
    if isfile(emitted_path)
      return emitted_path
    end
    @warn "MATPOWER Julia case generation did not produce a file; using the .m case" casefile = case_path emitted_path
  catch err
    @warn "MATPOWER Julia case generation failed; using the .m case" casefile = case_path exception = (err, catch_backtrace())
  end
  return case_path
end

function _resolve_powerflow_casefile(
  casefile::AbstractString,
  case_directory::AbstractString;
  ensure_casefile_fn = FetchMatpowerCase.ensure_casefile,
  emit_julia_case_fn = FetchMatpowerCase.emit_julia_case,
)::String
  requested = strip(casefile)
  isempty(requested) && throw(ArgumentError("PowerFlow casefile must not be empty."))
  occursin(r"^[A-Za-z][A-Za-z0-9+.-]*://", requested) && throw(ArgumentError("MATPOWER case URLs are not accepted."))

  extension = lowercase(splitext(requested)[2])
  extension in (".m", ".jl") || throw(ArgumentError("Unsupported casefile extension: $(requested) (expected .m or .jl)"))
  if isfile(requested)
    return _prefer_julia_casefile(requested; emit_julia_case_fn)
  end
  occursin(r"[\\/]", requested) && throw(ArgumentError("Case file not found: $(requested)"))

  trusted_directory = abspath(case_directory)
  mkpath(trusted_directory)
  resolved = try
    ensure_casefile_fn(requested; outdir = trusted_directory, to_jl = true)
  catch err
    fallback_name = extension == ".m" ? requested : first(splitext(requested)) * ".m"
    fallback_path = joinpath(trusted_directory, fallback_name)
    if isfile(fallback_path)
      @warn "MATPOWER case resolution could not generate the Julia case; using the downloaded .m case" casefile = fallback_path exception = (err, catch_backtrace())
      return abspath(fallback_path)
    end
    rethrow()
  end
  isfile(resolved) || throw(ArgumentError("Resolved MATPOWER case file not found: $(resolved)"))
  return _prefer_julia_casefile(resolved; emit_julia_case_fn)
end

"""
    start_powerflow_run(request::AbstractDict; case_directory=nothing) -> Dict{String,Any}

Start one local PowerFlow service run above [`run_sparlectra_api`](@ref). The
request must provide `casefile`, `config_file`, and `output_root`; optional
`config_overrides` are forwarded to the programmatic API. A unique run ID is
chosen before execution, and all generated files are written beneath
`output_root/run_id`. Completed API runs, including failed runs with a
`result.json`, are registered in memory and in the persistent run index. Public
service failures are returned as structured dictionaries.

When `case_directory` is provided by a trusted caller, bare `.m` or `.jl` case
names are resolved there through [`ensure_casefile`](@ref). Generated Julia
cases are preferred for execution, while existing path behavior is retained.
"""
function start_powerflow_run(request::AbstractDict; case_directory::Union{Nothing,AbstractString} = nothing, case_resolver = _resolve_powerflow_casefile)::Dict{String,Any}
  request_start = time_ns()
  phases = Dict{Symbol,Float64}()
  cancellation_token = _service_request_value(request, "cancellation_token", nothing)
  phase_callback = _service_request_value(request, "phase_callback", phase -> nothing)
  set_phase = phase -> (phase_callback(String(phase)); _check_powerflow_cancelled!(cancellation_token))
  set_phase("resolving_case")
  _check_powerflow_cancelled!(cancellation_token)
  casefile = _service_request_value(request, "casefile")
  casefile isa AbstractString && !isempty(strip(casefile)) || return _service_failure("missing_casefile", "PowerFlow service request requires a nonempty casefile.")
  if case_directory === nothing
    _check_powerflow_cancelled!(cancellation_token)
    isfile(casefile) || return _service_failure("missing_casefile", "MATPOWER case file not found: $(abspath(casefile))")
    casefile = abspath(casefile)
  else
    resolution_start = time_ns()
    _check_powerflow_cancelled!(cancellation_token)
    casefile = try
      case_resolver(casefile, case_directory)
    catch err
      return _service_failure("invalid_casefile", sprint(showerror, err))
    end
    phases[:case_resolution] = _api_elapsed_seconds(resolution_start)
    _check_powerflow_cancelled!(cancellation_token)
  end
  set_phase("preparing_configuration")

  config_file = _service_request_value(request, "config_file")
  config_file isa AbstractString && !isempty(strip(config_file)) || return _service_failure("missing_config_file", "PowerFlow service request requires a nonempty config_file.")
  isfile(config_file) || return _service_failure("missing_config_file", "Configuration file not found: $(abspath(config_file))")

  output_root = _service_request_value(request, "output_root")
  output_root isa AbstractString && !isempty(strip(output_root)) || return _service_failure("invalid_request", "PowerFlow service request requires a nonempty output_root.")

  config_overrides = _service_request_value(request, "config_overrides", Dict{String,Any}())
  config_overrides isa AbstractDict || return _service_failure("invalid_request", "config_overrides must be dictionary-like.")
  performance_timing = _service_request_value(request, "performance_timing", :off)
  run_diagnostics = _service_request_value(request, "run_diagnostics", false)
  run_diagnostics isa Bool || return _service_failure("invalid_request", "run_diagnostics must be boolean.")
  detailed_result_csv = _service_request_value(request, "detailed_result_csv", false)
  detailed_result_csv isa Bool || return _service_failure("invalid_request", "detailed_result_csv must be boolean.")
  detailed_result_csv_semicolon = _service_request_value(request, "detailed_result_csv_semicolon", false)
  detailed_result_csv_semicolon isa Bool || return _service_failure("invalid_request", "detailed_result_csv_semicolon must be boolean.")
  phases[:request_parse] = _api_elapsed_seconds(request_start)

  requested_run_id = _service_request_value(request, "run_id", nothing)
  run_id = requested_run_id === nothing ? string(uuid4()) : String(requested_run_id)
  _safe_powerflow_run_id(run_id) || return _service_failure("unsafe_run_id", "Unsafe PowerFlow run ID rejected.")
  root = abspath(output_root)
  output_dir = joinpath(root, run_id)
  result = try
    _run_sparlectra_api(
      casefile = casefile,
      config_file = config_file,
      output_dir = output_dir,
      config_overrides = config_overrides,
      performance_timing = performance_timing,
      run_diagnostics = run_diagnostics,
      detailed_result_csv = detailed_result_csv,
      detailed_result_csv_semicolon = detailed_result_csv_semicolon,
      phase_timings = phases,
      run_id = run_id,
      cancellation_token = cancellation_token,
      phase_callback = phase_callback,
    )
  catch err
    err isa PowerFlowAborted && rethrow()
    return _service_failure("execution_error", sprint(showerror, err, catch_backtrace()); run_id = run_id)
  end

  try
    lock(_POWERFLOW_SERVICE_LOCK) do
      _POWERFLOW_SERVICE_RUNS[result.run_id] = result
      _write_powerflow_run_index!(root, result)
    end
  catch err
    lock(_POWERFLOW_SERVICE_LOCK) do
      delete!(_POWERFLOW_SERVICE_RUNS, result.run_id)
    end
    return _service_failure("run_index_error", sprint(showerror, err, catch_backtrace()); run_id = run_id)
  end
  return to_dict(result)
end

"""
    get_powerflow_result(run_id::AbstractString) -> Dict{String,Any}

Return the transport-safe result metadata for a run started by
[`start_powerflow_run`](@ref) or recovered by
[`refresh_powerflow_run_registry!`](@ref), or a structured `run_not_found`
failure.
"""
function get_powerflow_result(run_id::AbstractString)::Dict{String,Any}
  result = _registered_powerflow_run(run_id)
  result === nothing && return _service_failure("run_not_found", "No PowerFlow run found for run_id $(run_id)."; run_id = run_id)
  return to_dict(result)
end

"""
    list_powerflow_artifacts(run_id::AbstractString)

List transport-safe artifact metadata for a registered PowerFlow run. Unknown
run IDs return a structured `run_not_found` failure dictionary.
"""
function list_powerflow_artifacts(run_id::AbstractString)
  result = _registered_powerflow_run(run_id)
  result === nothing && return _service_failure("run_not_found", "No PowerFlow run found for run_id $(run_id)."; run_id = run_id)
  return [to_dict(artifact) for artifact in collect_sparlectra_api_artifacts(result.output_dir)]
end

function _unsafe_artifact_name(name::String)::Bool
  isempty(name) && return true
  isabspath(name) && return true
  occursin('\\', name) && return true
  occursin(r"^[A-Za-z]:", name) && return true
  parts = split(replace(name, '\\' => '/'), '/'; keepempty = true)
  return any(part -> isempty(part) || part == "." || part == "..", parts)
end

function _artifact_belongs_to_run(artifact::SparlectraApiArtifact, output_dir::String)::Bool
  artifact.exists && isfile(artifact.path) || return false
  root = realpath(output_dir)
  path = realpath(artifact.path)
  relative = relpath(path, root)
  return first(splitpath(relative)) != ".."
end

"""
    resolve_powerflow_artifact(run_id::AbstractString, artifact_name::AbstractString)

Resolve an artifact by its metadata `name`, never by an arbitrary filesystem
path. Absolute paths, traversal components, Windows-style paths, missing files,
and artifacts escaping the selected run directory are returned as structured
failures.
"""
function resolve_powerflow_artifact(run_id::AbstractString, artifact_name::AbstractString)
  result = _registered_powerflow_run(run_id)
  result === nothing && return _service_failure("run_not_found", "No PowerFlow run found for run_id $(run_id)."; run_id = run_id)

  name = String(artifact_name)
  _unsafe_artifact_name(name) && return _service_failure("unsafe_artifact_name", "Unsafe artifact name rejected: $(name)"; run_id = run_id)

  artifacts = collect_sparlectra_api_artifacts(result.output_dir)
  index = findfirst(artifact -> artifact.name == name, artifacts)
  index === nothing && return _service_failure("artifact_not_found", "No artifact named $(name) belongs to PowerFlow run $(run_id)."; run_id = run_id)
  artifact = artifacts[index]
  isfile(artifact.path) || return _service_failure("artifact_not_found", "Artifact $(name) is no longer available for PowerFlow run $(run_id)."; run_id = run_id)
  _artifact_belongs_to_run(artifact, result.output_dir) || return _service_failure("unsafe_artifact_name", "Artifact $(name) does not resolve inside PowerFlow run $(run_id)."; run_id = run_id)
  return artifact
end
