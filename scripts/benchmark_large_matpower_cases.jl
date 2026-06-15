using Dates
using Printf
using Sparlectra

const DEFAULT_LARGE_CASE_NAMES = [
  "case118.m",
  "case1354pegase.m",
  "case2869pegase.m",
  "case9241pegase.m",
  "case_ACTIVSg10k.m",
]

const BENCHMARK_PHASES = [
  "resolving_case", "checking_case_cache", "preparing_case_file", "preparing_configuration",
  "reading_matpower_case", "loading_julia_case", "converting_matpower_case", "parsing_matpower_file",
  "building_sparlectra_net", "applying_import_options", "preparing_start_values", "building_ybus",
  "solver_initialization", "start_projection", "solving_powerflow", "newton_iteration", "linear_solve",
  "postprocessing_result", "writing_diagnostics", "writing_csv_artifacts", "writing_artifacts",
  "finalizing_success", "finalizing_failed", "finalizing_aborted",
]

function _parse_benchmark_args(args = ARGS)
  default_output_root = joinpath(Sparlectra.default_webui_output_root(), "large_case_benchmarks")
  options = Dict{String,Any}(
    "output_root" => default_output_root,
    "case_dir" => nothing,
    "cases" => copy(DEFAULT_LARGE_CASE_NAMES),
    "fetch_missing" => true,
    "allow_generated_jl_cache" => false,
  )
  index = firstindex(args)
  while index <= lastindex(args)
    arg = args[index]
    if arg == "--case-dir"
      index += 1
      index <= lastindex(args) || throw(ArgumentError("--case-dir requires a value."))
      options["case_dir"] = args[index]
    elseif startswith(arg, "--case-dir=")
      options["case_dir"] = split(arg, "="; limit = 2)[2]
    elseif arg == "--output-root"
      index += 1
      index <= lastindex(args) || throw(ArgumentError("--output-root requires a value."))
      options["output_root"] = args[index]
    elseif startswith(arg, "--output-root=")
      options["output_root"] = split(arg, "="; limit = 2)[2]
    elseif arg == "--cases"
      index += 1
      index <= lastindex(args) || throw(ArgumentError("--cases requires a comma-separated value."))
      options["cases"] = split(args[index], ",")
    elseif startswith(arg, "--cases=")
      options["cases"] = split(split(arg, "="; limit = 2)[2], ",")
    elseif arg == "--no-fetch-missing"
      options["fetch_missing"] = false
    elseif arg == "--allow-generated-jl-cache"
      options["allow_generated_jl_cache"] = true
    else
      throw(ArgumentError("Unknown benchmark option: $(arg)"))
    end
    index += 1
  end
  return options
end

function _canonical_benchmark_case(requested_case::AbstractString, case_dir::AbstractString; allow_generated_jl_cache::Bool = false)
  name = basename(strip(requested_case))
  isempty(name) && throw(ArgumentError("Benchmark case name must not be empty."))
  root, ext = splitext(name)
  if isempty(ext)
    name = root * ".m"
    ext = ".m"
  end
  ext = lowercase(ext)
  if ext == ".jl" && !allow_generated_jl_cache
    m_name = root * ".m"
    m_path = joinpath(case_dir, m_name)
    if isfile(m_path)
      return (requested = basename(requested_case), resolved = m_name, canonical = m_name, used_generated_jl_cache = false, missing = false, message = nothing)
    end
    return (requested = basename(requested_case), resolved = nothing, canonical = m_name, used_generated_jl_cache = false, missing = true, message = "Canonical .m source $(m_name) is missing; generated MATPOWER .jl cache files are not benchmarked by default.")
  elseif ext == ".jl"
    return (requested = basename(requested_case), resolved = name, canonical = name, used_generated_jl_cache = true, missing = false, message = nothing)
  elseif ext == ".m"
    return (requested = basename(requested_case), resolved = name, canonical = name, used_generated_jl_cache = false, missing = false, message = nothing)
  end
  throw(ArgumentError("Unsupported benchmark case extension: $(requested_case) (expected .m or .jl)."))
end

function _webui_case_cache_dir(output_root::AbstractString, case_dir)::String
  return abspath(case_dir === nothing ? Sparlectra.default_webui_case_cache_dir(output_root) : String(case_dir))
end

function _cache_contains_case(case_dir::AbstractString, requested_case::AbstractString)::Bool
  path = joinpath(case_dir, basename(requested_case))
  return isfile(path) && filesize(path) > 0
end

function _phase_total(phase_timings, phase::AbstractString)
  values = Float64[]
  for timing in phase_timings
    get(timing, "phase", "") == phase || continue
    elapsed = get(timing, "elapsed_seconds", nothing)
    elapsed === nothing || push!(values, Float64(elapsed))
  end
  isempty(values) && return nothing
  return sum(values)
end

function _phase_count(phase_timings, phase::AbstractString)
  return count(timing -> get(timing, "phase", "") == phase, phase_timings)
end

function _bottleneck(row::AbstractDict)
  candidates = Pair{String,Float64}[]
  for (label, key) in (("reader/parser", "reading_matpower_case"), ("jl_load_parse", "loading_julia_case"), ("network_builder", "building_sparlectra_net"), ("solver", "solving_powerflow"), ("artifacts", "writing_artifacts"))
    value = get(row["phase_timings"], key, nothing)
    value isa Number || continue
    push!(candidates, label => Float64(value))
  end
  isempty(candidates) && return "not_recorded"
  return first(sort(candidates; by = pair -> pair.second, rev = true)).first
end

function _benchmark_run_id(requested_case::AbstractString, mode::AbstractString)::String
  safe_case = replace(basename(requested_case), r"[^A-Za-z0-9_.-]" => "_")
  return safe_case * "-" * mode * "-" * Dates.format(Dates.now(Dates.UTC), dateformat"yyyymmddTHHMMSSsss")
end

function _run_service_case(requested_case::AbstractString, case_dir::AbstractString, output_root::AbstractString; mode::AbstractString = "baseline")
  run_id = _benchmark_run_id(requested_case, mode)
  request = Dict{String,Any}(
    "casefile" => basename(requested_case),
    "config_file" => Sparlectra.DEFAULT_SPARLECTRA_CONFIG_PATH,
    "output_root" => abspath(output_root),
    "run_id" => run_id,
    "performance_timing" => "compact",
    "run_diagnostics" => mode == "artifact_heavy",
    "detailed_result_csv" => mode == "artifact_heavy",
    "detailed_result_csv_format" => "excel_de",
    "config_overrides" => Dict("benchmark.enabled" => false, "output.logfile_results" => "full"),
  )
  started = time_ns()
  result = Sparlectra.start_powerflow_run(request; case_directory = case_dir)
  total_elapsed = (time_ns() - started) / 1.0e9
  artifacts = haskey(result, "run_id") && get(result, "success", false) ? Sparlectra.list_powerflow_artifacts(result["run_id"]) : Any[]
  for artifact in artifacts
    name = get(artifact, "name", "")
    name in ("run.log", "performance.log", "result.json") || continue
    try
      Sparlectra.resolve_powerflow_artifact(result["run_id"], name)
    catch
    end
  end
  phase_timings = get(result, "service_phase_timings", Any[])
  phase_map = Dict{String,Any}()
  for phase in BENCHMARK_PHASES
    phase_map[phase] = something(_phase_total(phase_timings, phase), "not_recorded")
  end
  resolved_case = get(result, "casefile", nothing)
  row = Dict{String,Any}(
    "case_name" => splitext(basename(requested_case))[1],
    "requested_case" => basename(requested_case),
    "resolved_case" => resolved_case,
    "canonical_source" => basename(requested_case),
    "used_generated_jl_cache" => lowercase(splitext(something(resolved_case, ""))[2]) == ".jl",
    "case_path" => resolved_case,
    "case_extension" => lowercase(splitext(something(resolved_case, requested_case))[2]),
    "case_size_bytes" => resolved_case isa AbstractString && isfile(resolved_case) ? filesize(resolved_case) : 0,
    "mode" => String(mode),
    "status" => get(result, "success", false) ? "ok" : "failed",
    "total_elapsed_seconds" => total_elapsed,
    "final_status" => get(result, "status", "unknown"),
    "error_message" => get(result, "message", nothing),
    "run_id" => get(result, "run_id", run_id),
    "output_dir" => get(result, "output_dir", joinpath(output_root, run_id)),
    "phase_timings" => phase_map,
    "phase_counts" => Dict(phase => _phase_count(phase_timings, phase) for phase in BENCHMARK_PHASES),
    "artifacts" => artifacts,
    "case_cache_dir" => case_dir,
  )
  row["bottleneck"] = _bottleneck(row)
  return row
end

function _missing_row(requested_case::AbstractString, case_dir::AbstractString)
  return Dict{String,Any}(
    "case_name" => splitext(basename(requested_case))[1],
    "requested_case" => basename(requested_case),
    "case_path" => nothing,
    "case_extension" => lowercase(splitext(requested_case)[2]),
    "case_size_bytes" => 0,
    "mode" => "baseline",
    "status" => "missing",
    "total_elapsed_seconds" => nothing,
    "final_status" => "missing",
    "error_message" => "Case is not present in the Web UI cache $(case_dir). Run this benchmark without --no-fetch-missing, or fetch the case through the Web UI/service cache by submitting the same case name in the Web UI.",
    "phase_timings" => Dict(phase => "not_recorded" for phase in BENCHMARK_PHASES),
    "phase_counts" => Dict(phase => 0 for phase in BENCHMARK_PHASES),
    "bottleneck" => "not_recorded",
    "case_cache_dir" => case_dir,
  )
end

function _missing_row(requested_case::AbstractString, case_dir::AbstractString, message::AbstractString; canonical_source::AbstractString = basename(requested_case))
  row = _missing_row(requested_case, case_dir)
  row["error_message"] = String(message)
  row["canonical_source"] = canonical_source
  row["resolved_case"] = nothing
  row["used_generated_jl_cache"] = false
  return row
end

function _write_benchmark_json(path::AbstractString, rows, metadata::AbstractDict)
  mkpath(dirname(path))
  open(path, "w") do io
    Sparlectra._write_json(io, Dict(
      "generated_at" => string(Dates.now(Dates.UTC)),
      "metadata" => Sparlectra._api_transport_value(metadata),
      "results" => Sparlectra._api_transport_value(rows),
    ))
  end
  return path
end

function _fmt(value)
  value === nothing && return ""
  value == "not_recorded" && return ""
  value isa Number && return @sprintf("%.3f", Float64(value))
  return String(value)
end

function _write_benchmark_markdown(path::AbstractString, rows, metadata::AbstractDict)
  mkpath(dirname(path))
  open(path, "w") do io
    println(io, "# Sparlectra large MATPOWER benchmark")
    println(io)
    println(io, "Generated: ", Dates.now(Dates.UTC), " UTC")
    println(io, "Web UI case cache: `", metadata["case_cache_dir"], "`")
    println(io, "Output root: `", metadata["output_root"], "`")
    println(io)
    println(io, "| case | ext | mode | status | total_s | read_s | jl_load_s | build_net_s | solve_s | artifacts_s | bottleneck |")
    println(io, "|---|---:|---|---|---:|---:|---:|---:|---:|---:|---|")
    for row in rows
      phases = row["phase_timings"]
      println(io, "| ", row["requested_case"], " | ", something(row["case_extension"], ""), " | ", get(row, "mode", ""), " | ", row["status"], " | ", _fmt(row["total_elapsed_seconds"]), " | ", _fmt(get(phases, "reading_matpower_case", nothing)), " | ", _fmt(get(phases, "loading_julia_case", nothing)), " | ", _fmt(get(phases, "building_sparlectra_net", nothing)), " | ", _fmt(get(phases, "solving_powerflow", nothing)), " | ", _fmt(get(phases, "writing_artifacts", nothing)), " | ", row["bottleneck"], " |")
    end
    println(io)
    println(io, "Missing cases should be fetched through the normal Web UI/service cache path by submitting the same bare case name in the Web UI or rerunning this helper without `--no-fetch-missing`.")
    println(io, "Operation-log phases are high-level. Detailed solver phases such as building_ybus, newton_iteration, and linear_solve are retained in each run's performance.log and service_phase_timings.")
  end
  return path
end

function benchmark_large_matpower_cases(; output_root::AbstractString = joinpath(Sparlectra.default_webui_output_root(), "large_case_benchmarks"), case_dir = nothing, cases = DEFAULT_LARGE_CASE_NAMES, fetch_missing::Bool = true, allow_generated_jl_cache::Bool = false)
  root = abspath(output_root)
  cache_dir = _webui_case_cache_dir(root, case_dir)
  mkpath(root)
  mkpath(cache_dir)
  rows = Dict{String,Any}[]
  for requested in String.(cases)
    case_name = strip(requested)
    isempty(case_name) && continue
    canonical = _canonical_benchmark_case(case_name, cache_dir; allow_generated_jl_cache)
    if canonical.missing
      push!(rows, _missing_row(canonical.requested, cache_dir, canonical.message; canonical_source = canonical.canonical))
      continue
    end
    if !fetch_missing && !_cache_contains_case(cache_dir, canonical.resolved)
      push!(rows, _missing_row(canonical.requested, cache_dir, "Canonical .m source $(canonical.resolved) is not present in the Web UI cache $(cache_dir)."; canonical_source = canonical.canonical))
      continue
    end
    row = _run_service_case(canonical.resolved, cache_dir, root; mode = "baseline")
    row["requested_case"] = canonical.requested
    row["canonical_source"] = canonical.canonical
    row["used_generated_jl_cache"] = canonical.used_generated_jl_cache
    push!(rows, row)
    if startswith(canonical.resolved, "case1354pegase") && endswith(lowercase(canonical.resolved), ".m")
      heavy = _run_service_case(canonical.resolved, cache_dir, root; mode = "artifact_heavy")
      heavy["requested_case"] = canonical.requested
      heavy["canonical_source"] = canonical.canonical
      heavy["used_generated_jl_cache"] = canonical.used_generated_jl_cache
      push!(rows, heavy)
    end
  end
  metadata = Dict{String,Any}(
    "output_root" => root,
    "case_cache_dir" => cache_dir,
    "fetch_missing" => fetch_missing,
    "allow_generated_jl_cache" => allow_generated_jl_cache,
    "case_resolution" => "start_powerflow_run(...; case_directory=case_cache_dir), matching the Web UI/service cache path",
  )
  json_path = joinpath(root, "benchmark_large_cases.json")
  md_path = joinpath(root, "benchmark_large_cases.md")
  _write_benchmark_json(json_path, rows, metadata)
  _write_benchmark_markdown(md_path, rows, metadata)
  return Dict("json" => json_path, "markdown" => md_path, "metadata" => metadata, "results" => rows)
end

function main(args = ARGS)
  options = _parse_benchmark_args(args)
  summary = benchmark_large_matpower_cases(; output_root = abspath(options["output_root"]), case_dir = options["case_dir"], cases = options["cases"], fetch_missing = options["fetch_missing"], allow_generated_jl_cache = options["allow_generated_jl_cache"])
  println("Web UI case cache: ", summary["metadata"]["case_cache_dir"])
  println("Wrote ", summary["json"])
  println("Wrote ", summary["markdown"])
  for row in summary["results"]
    println(rpad(row["requested_case"], 22), rpad(something(row["case_extension"], ""), 5), rpad(row["status"], 9), _fmt(row["total_elapsed_seconds"]), " s  ", row["bottleneck"])
    if row["status"] in ("missing", "failed") && get(row, "error_message", nothing) !== nothing
      println("  ", row["error_message"])
    end
  end
  return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
  Base.invokelatest(main)
end
