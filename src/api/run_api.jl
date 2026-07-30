# Copyright 2023–2026 Udo Schmitz
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

using UUIDs: uuid4

const _SPARLECTRA_API_SCHEMA_VERSION = "1.0"
const WEBUI_PERFORMANCE_TIMING_VALUES = (:off, :compact, :full)
const Q_LIMIT_LOG_ARTIFACT = "q_limit.log"
const MATPOWER_DCLINE_ARTIFACT = "matpower_dcline.csv"


mutable struct PowerFlowPhaseTimingRecorder
  timings::Vector{Dict{String,Any}}
  active_index::Union{Nothing,Int}
end

PowerFlowPhaseTimingRecorder() = PowerFlowPhaseTimingRecorder(Dict{String,Any}[], nothing)

function _api_elapsed_seconds(start_ns::UInt64)::Float64
  return (time_ns() - start_ns) / 1.0e9
end

function _api_datetime_string(value::Dates.DateTime)::String
  return Dates.format(value, dateformat"yyyy-mm-ddTHH:MM:SS.sss") * "Z"
end

function _write_matpower_dcline_artifact(output_path::AbstractString, net::Net; format = "technical")::Union{Nothing,String}
  rows = getproperty(net, :matpowerDclineMetadata)
  isempty(rows) && return nothing
  columns = (
    :orig_index,
    :from_bus,
    :to_bus,
    :from_bus_name,
    :to_bus_name,
    :status,
    :pf_mw,
    :input_pt_mw,
    :effective_pt_mw,
    :loss0_mw,
    :loss1,
    :loss_mw,
    :qf_mvar,
    :qt_mvar,
    :vf_pu,
    :vt_pu,
    :qminf_mvar,
    :qmaxf_mvar,
    :qmint_mvar,
    :qmaxt_mvar,
    :from_prosumer,
    :to_prosumer,
    :from_voltage_controlled,
    :to_voltage_controlled,
  )
  _write_namedtuple_csv(joinpath(output_path, MATPOWER_DCLINE_ARTIFACT), rows, columns; format = format)
  return MATPOWER_DCLINE_ARTIFACT
end

function _effective_config_with_runtime_case(effective_raw, case_path::AbstractString, config::SparlectraConfig; config_sources = nothing)
  raw = deepcopy(effective_raw)
  runtime = get!(raw, "runtime", Dict{String,Any}())
  runtime isa AbstractDict || (runtime = raw["runtime"] = Dict{String,Any}())
  # The effective config may contain default form choices, but diagnostics need
  # the resolved runtime casefile so support analysis can identify what was run.
  runtime["casefile"] = String(case_path)
  runtime["case_name"] = splitext(basename(case_path))[1]
  runtime["case_source"] = "webui_mpower_data"
  runtime["configured_default_casefile"] = config.matpower.case
  if config_sources !== nothing
    raw["_config_sources"] = config_sources
  end
  return raw
end

function _phase_elapsed_seconds(started_at::AbstractString, ended_at::AbstractString)::Float64
  start = Dates.DateTime(replace(String(started_at), "Z" => ""), dateformat"yyyy-mm-ddTHH:MM:SS.sss")
  finish = Dates.DateTime(replace(String(ended_at), "Z" => ""), dateformat"yyyy-mm-ddTHH:MM:SS.sss")
  return max(0.0, Dates.value(finish - start) / 1000)
end

function _complete_active_phase!(recorder::PowerFlowPhaseTimingRecorder, status::AbstractString = "completed")
  index = recorder.active_index
  index === nothing && return nothing
  timing = recorder.timings[index]
  haskey(timing, "ended_at") && timing["ended_at"] !== nothing && return nothing
  ended_at = _api_datetime_string(Dates.now(Dates.UTC))
  timing["ended_at"] = ended_at
  timing["elapsed_seconds"] = _phase_elapsed_seconds(timing["started_at"], ended_at)
  timing["status"] = String(status)
  recorder.active_index = nothing
  return timing
end

function _start_service_phase!(recorder::PowerFlowPhaseTimingRecorder, phase::AbstractString; status::AbstractString = "running")
  phase_name = String(phase)
  if recorder.active_index !== nothing && recorder.timings[recorder.active_index]["phase"] == phase_name
    return recorder.timings[recorder.active_index]
  end
  _complete_active_phase!(recorder)
  timing = Dict{String,Any}("phase" => phase_name, "started_at" => _api_datetime_string(Dates.now(Dates.UTC)), "ended_at" => nothing, "elapsed_seconds" => nothing, "status" => String(status))
  push!(recorder.timings, timing)
  recorder.active_index = length(recorder.timings)
  return timing
end

function _phase_elapsed_lookup(phase_timings::AbstractVector, phase::AbstractString)
  for timing in phase_timings
    get(timing, "phase", "") == phase && return get(timing, "elapsed_seconds", nothing)
  end
  return nothing
end

function _api_timing_mode(value)::Symbol
  mode = value isa Symbol ? value : Symbol(lowercase(strip(String(value))))
  mode in WEBUI_PERFORMANCE_TIMING_VALUES || throw(ArgumentError("Performance timing mode must be one of $(join(WEBUI_PERFORMANCE_TIMING_VALUES, ", ")); got $(repr(value))."))
  return mode
end

include("run_failures.jl")

include("run_diagnostic_artifacts.jl")
include("run_self_check.jl")
include("run_finalization.jl")

function _resolve_detailed_csv_format(value)::NamedTuple
  name = String(value)
  name == "technical" && return (name = name, delimiter = ',', decimal_separator = '.', thousands_separator = "")
  name == "excel_de" && return (name = name, delimiter = ';', decimal_separator = ',', thousands_separator = ".")
  name == "excel_us" && return (name = name, delimiter = ',', decimal_separator = '.', thousands_separator = ",")
  throw(ArgumentError("Unsupported detailed_result_csv_format \"$(name)\". Expected technical, excel_de, or excel_us."))
end

function _group_csv_integer(text::AbstractString, separator::AbstractString)::String
  isempty(separator) && return String(text)

  sign = startswith(text, "-") ? "-" : ""
  digits = isempty(sign) ? String(text) : text[2:end]

  # Thousands grouping is only valid for integer digit strings. Leave textual
  # values such as Bool strings unchanged if they accidentally reach this helper.
  (isempty(digits) || !all(isdigit, digits)) && return String(text)

  first_group = mod(length(digits), 3)
  first_group == 0 && (first_group = 3)
  groups = String[digits[1:first_group]]
  for start = (first_group+1):3:length(digits)
    push!(groups, digits[start:(start+2)])
  end
  return sign * join(groups, separator)
end

function _format_csv_number(value::Integer, format)::String
  return _group_csv_integer(string(value), format.thousands_separator)
end

function _format_csv_number(value::AbstractFloat, format)::String
  isnan(value) && return "NaN"
  isinf(value) && return signbit(value) ? "-Inf" : "Inf"
  technical = @sprintf("%.15g", value)
  if format.name != "technical"
    exponent_marker = findfirst(character -> character in ('e', 'E'), technical)
    if exponent_marker !== nothing
      mantissa = technical[begin:prevind(technical, exponent_marker)]
      exponent = parse(Int, technical[nextind(technical, exponent_marker):end])
      sign = startswith(mantissa, "-") ? "-" : ""
      unsigned = isempty(sign) ? mantissa : mantissa[2:end]
      dot_index = findfirst(==('.'), unsigned)
      fractional_digits = dot_index === nothing ? 0 : ncodeunits(unsigned) - dot_index
      digits = replace(unsigned, "." => "")
      decimal_position = ncodeunits(digits) - fractional_digits + exponent
      if decimal_position <= 0
        technical = sign * "0." * repeat("0", -decimal_position) * digits
      elseif decimal_position >= ncodeunits(digits)
        technical = sign * digits * repeat("0", decimal_position - ncodeunits(digits))
      else
        technical = sign * digits[1:decimal_position] * "." * digits[(decimal_position+1):end]
      end
    end
    dot_index = findfirst(==('.'), technical)
    if dot_index !== nothing
      last_nonzero = lastindex(technical)
      while last_nonzero > dot_index && technical[last_nonzero] == '0'
        last_nonzero = prevind(technical, last_nonzero)
      end
      technical = last_nonzero == dot_index ? technical[begin:prevind(technical, dot_index)] : technical[begin:last_nonzero]
      isempty(technical) && (technical = "0")
      technical == "-0" && (technical = "0")
    end
  end
  exponent_marker = format.name == "technical" ? findfirst(character -> character in ('e', 'E'), technical) : nothing
  mantissa = exponent_marker === nothing ? technical : technical[begin:prevind(technical, exponent_marker)]
  exponent = exponent_marker === nothing ? "" : technical[exponent_marker:end]
  dot_index = findfirst(==('.'), mantissa)
  integer_text = dot_index === nothing ? mantissa : mantissa[begin:prevind(mantissa, dot_index)]
  integer_part = _group_csv_integer(integer_text, format.thousands_separator)
  fractional_part = dot_index === nothing ? "" : string(format.decimal_separator, mantissa[nextind(mantissa, dot_index):end])
  return integer_part * fractional_part * exponent
end

struct CsvFormatRuntime
  name::String
  delimiter::Char
  decimal_separator::Char
  thousands_separator::String
end

CsvFormatRuntime(format) = CsvFormatRuntime(String(format.name), format.delimiter, format.decimal_separator, String(format.thousands_separator))

function _csv_needs_quotes(text::AbstractString, delimiter::Char)::Bool
  for character in text
    character in (delimiter, '"', '\r', '\n') && return true
  end
  return false
end

function write_csv_cell!(io::IO, value, delimiter::Char, fmt::CsvFormatRuntime)
  value === missing && return nothing
  value === nothing && return nothing
  text = value isa Bool ? string(value) : value isa Integer ? _format_csv_number(value, fmt) : value isa AbstractFloat ? _format_csv_number(value, fmt) : string(value)
  if _csv_needs_quotes(text, delimiter)
    print(io, '"')
    for character in text
      character == '"' && print(io, '"')
      print(io, character)
    end
    print(io, '"')
  else
    print(io, text)
  end
  return nothing
end

function write_csv_row_direct!(io::IO, delimiter::Char, fmt::CsvFormatRuntime, values...)
  first = true
  for value in values
    first || print(io, delimiter)
    write_csv_cell!(io, value, delimiter, fmt)
    first = false
  end
  println(io)
  return nothing
end

function _format_csv_value(value, format)::String
  value === missing && return ""
  value === nothing && return ""

  # Bool is an Integer subtype in Julia. Format it before Integer values so
  # Excel-oriented thousands grouping cannot turn true/false into t.rue/fa.lse.
  value isa Bool && return string(value)

  value isa Integer && return _format_csv_number(value, format)
  value isa AbstractFloat && return _format_csv_number(value, format)
  return string(value)
end

function _csv_field(value, delimiter::Char, format = _resolve_detailed_csv_format("technical"))::String
  value === missing && return ""
  value === nothing && return ""
  text = _format_csv_value(value, format)
  if any(character -> character in (delimiter, '"', '\r', '\n'), text)
    return "\"" * replace(text, "\"" => "\"\"") * "\""
  end
  return text
end

include("run_csv_exports.jl")

function _csv_solution_quality(raw_result::SparlectraRunResult)::String
  return raw_result.final_converged && raw_result.solution_available ? "converged" : "not_converged_last_iterate"
end

include("run_lifecycle_metadata.jl")

function _write_numerical_outcome_summary(io::IO, raw_result::SparlectraRunResult)
  raw_result.final_converged && raw_result.solution_available && return nothing
  println(io)
  println(io, "Numerical outcome: not_converged")
  println(io, "Primary reason: ", raw_result.reason_text)
  println(io, "Final mismatch: ", isfinite(raw_result.final_mismatch) ? raw_result.final_mismatch : (isnan(raw_result.final_mismatch) ? "NaN" : "unavailable"))
  println(io, "Final mismatch status: ", _final_mismatch_status(raw_result))
  qlimit_ok = hasproperty(raw_result.diagnostics, :q_limit_active_set_ok) ? getproperty(raw_result.diagnostics, :q_limit_active_set_ok) : nothing
  qlimit_ok === nothing || println(io, "Q-limit active-set convergence: ", qlimit_ok ? "yes" : "no")
  for (label, field) in (
    ("PV->PQ switching events", :pv_pq_switching_events),
    ("Q-limit active-set changes", :qlimit_active_set_changes),
    ("Q-limit re-enable events", :qlimit_reenable_events),
    ("Guarded narrow-Q PV buses", :guarded_narrow_q_pv_buses),
  )
    hasproperty(raw_result.diagnostics, field) && println(io, label, ": ", getproperty(raw_result.diagnostics, field))
  end
  return nothing
end

include("run_matpower_artifacts.jl")

function _final_outcome_payload(raw_result::SparlectraRunResult)::Dict{String,Any}
  return Dict{String,Any}(
    "converged" => raw_result.final_converged,
    "numerical_converged" => raw_result.numerical_converged,
    "solution_available" => raw_result.solution_available,
    "solver" => hasproperty(raw_result.diagnostics, :solver) ? String(raw_result.diagnostics.solver) : "rectangular",
    "iterations" => raw_result.iterations,
    "final_mismatch" => isfinite(raw_result.final_mismatch) ? raw_result.final_mismatch : nothing,
    "reason" => String(raw_result.reason),
    "island_wise_all_converged" => hasproperty(raw_result.diagnostics, :island_wise_all_converged) ? raw_result.diagnostics.island_wise_all_converged : false,
    "post_merge_validation_status" => hasproperty(raw_result.diagnostics, :post_merge_validation_status) ? String(raw_result.diagnostics.post_merge_validation_status) : "not_applicable",
    "post_merge_final_mismatch" => hasproperty(raw_result.diagnostics, :post_merge_final_mismatch) && isfinite(raw_result.diagnostics.post_merge_final_mismatch) ? raw_result.diagnostics.post_merge_final_mismatch : nothing,
    "post_merge_mismatch_status" => hasproperty(raw_result.diagnostics, :post_merge_mismatch_status) ? String(raw_result.diagnostics.post_merge_mismatch_status) : "not_applicable",
  )
end

const DTF_FOR001_UNSUPPORTED_DCLINE_MESSAGE = "DC lines are currently not supported by the native DTF/MATPOWER power-flow path."

function _normalize_case_format(value)::Symbol
  format = value isa Symbol ? value : Symbol(lowercase(strip(String(value))))
  format in (:auto, :matpower, :dtf_for001, :cgmes) || throw(ArgumentError("case_format must be auto, matpower, dtf_for001, or cgmes; got $(repr(value))."))
  return format
end

"""Cheap CGMES sniff: a folder, or a ZIP holding CGMES RDF/XML payloads."""
function _looks_like_cgmes(case_path::AbstractString)::Bool
  isdir(case_path) && return true
  lowercase(splitext(case_path)[2]) == ".zip" || return false
  try
    return !isempty(CGMESImporter.collectCGMESFiles(case_path))
  catch
    return false
  end
end

function _detect_case_format(case_path::AbstractString; requested::Symbol = :auto)::Symbol
  requested !== :auto && return requested
  _looks_like_cgmes(case_path) && return :cgmes
  ext = lowercase(splitext(case_path)[2])
  ext in (".m", ".jl") && return :matpower
  text = read(case_path, String)
  # Native FOR001 test data has explicit section cards and a DTF size card.  Do
  # not infer arbitrary .DAT files unless these FOR001 markers are present.
  has_for001_sections = occursin("##Z", text) && occursin("##L", text) && occursin("##K", text)
  has_size_card = occursin(r"(?m)^##G\s+\d+\s+\d+", text)
  if has_for001_sections && has_size_card
    return :dtf_for001
  end
  ext == ".dat" && throw(ArgumentError("Ambiguous .DAT input; set case_format = :dtf_for001 to use the experimental/internal native DTF path."))
  return :matpower
end

function _reject_dtf_dcline_like_content!(case_path::AbstractString)
  # open(...) do closes the handle even when the ArgumentError below aborts
  # the scan mid-file; eachline(path) would keep the descriptor open until GC
  # (blocks file deletion on Windows, EBUSY).
  open(case_path) do io
    for (line_no, line) in enumerate(eachline(io))
      if occursin(r"(?i)\b(HVDC|DCLINE|DC\s*LINE)\b", line)
        throw(ArgumentError("unsupported_dtf_dc_line at line $(line_no): $(DTF_FOR001_UNSUPPORTED_DCLINE_MESSAGE)"))
      end
    end
  end
  return nothing
end

function _dtf_metadata(case, requested::Symbol, detected::Symbol; for002_file=nothing, run_dtf_outages=false, matpower_export_requested=false, matpower_export_file=nothing)
  summary = DTFImporter.case_summary(case)
  return Dict{String,Any}(
    "input_format" => String(requested),
    "input_format_detected" => String(detected),
    "native_dtf_import_used" => true,
    "dtf_bus_count" => summary.bus_count,
    "dtf_branch_count" => summary.branch_count,
    "dtf_outage_count" => summary.outage_count,
    "dtf_slack_bus" => summary.slack_bus,
    "dtf_nominal_voltage_levels" => case.nominal_voltages_kv,
    "for002_reference_used" => for002_file !== nothing && !isempty(String(for002_file)),
    "for002_reference_file" => for002_file === nothing ? nothing : String(for002_file),
    "outage_validation_requested" => run_dtf_outages,
    "matpower_export_requested" => matpower_export_requested,
    "matpower_export_file" => matpower_export_file,
    "dcline_status" => "unsupported_not_present",
    "unsupported_dcline_status" => "not_detected",
  )
end

function _write_dtf_summary_artifacts(output_path::AbstractString, case, metadata::AbstractDict)
  csv = joinpath(output_path, "dtf_import_summary.csv")
  open(csv, "w") do io
    println(io, "key,value")
    for key in ("input_format", "input_format_detected", "dtf_bus_count", "dtf_branch_count", "dtf_outage_count", "dtf_slack_bus", "dcline_status")
      println(io, _csv_field(key, ','), ",", _csv_field(get(metadata, key, ""), ','))
    end
  end
  md = joinpath(output_path, "dtf_import_summary.md")
  open(md, "w") do io
    println(io, "# Native DTF input summary (experimental/internal)\n")
    println(io, "- buses: ", length(case.buses))
    println(io, "- branches: ", length(case.branches))
    println(io, "- outages: ", length(case.outages))
    println(io, "- slack bus: ", case.size.slack)
    println(io, "- nominal voltage levels [kV]: ", join(case.nominal_voltages_kv, ", "))
    println(io, "- DC-line status: ", get(metadata, "dcline_status", "unsupported_not_present"))
  end
  return [basename(csv), basename(md)]
end

function _write_for002_placeholder_artifacts(output_path::AbstractString, for002_file)
  for002_file === nothing && return String[]
  isempty(String(for002_file)) && return String[]
  isfile(String(for002_file)) || throw(ArgumentError("invalid FOR002 reference file: $(for002_file)"))
  artifacts = String[]
  md = joinpath(output_path, "dtf_for002_base_comparison.md")
  open(md, "w") do io
    println(io, "# FOR002 reference comparison")
    println(io)
    println(io, "FOR002 reference file: `", String(for002_file), "`")
    println(io, "Detailed native validation remains available through the existing validation examples.")
  end
  push!(artifacts, basename(md))
  csv = joinpath(output_path, "dtf_for002_base_metrics.csv")
  open(csv, "w") do io
    println(io, "metric,value")
    println(io, "for002_reference_file,", _csv_field(String(for002_file), ','))
    println(io, "comparison_status,reference_recorded")
  end
  push!(artifacts, basename(csv))
  return artifacts
end

function _selected_dtf_outages(case, mode, selection)
  mode_symbol = mode isa Symbol ? mode : Symbol(lowercase(strip(String(mode))))
  mode_symbol == :none && return DTFImporter.DTFOutage[]
  mode_symbol == :all && return collect(case.outages)
  mode_symbol == :selected || throw(ArgumentError("dtf_outage_selection_mode must be none, all, or selected."))
  selected = Set(String.(selection isa AbstractVector ? selection : [selection]))
  outages = DTFImporter.DTFOutage[]
  for outage in case.outages
    label = DTFImporter.outage_label(outage)
    if string(outage.index) in selected || label in selected
      push!(outages, outage)
    end
  end
  length(outages) == length(selected) || throw(ArgumentError("missing_or_ambiguous_outage_selection: selected DTF outage was not found unambiguously."))
  return outages
end

function _run_dtf_outages(case_path::AbstractString, case, config, output_path::AbstractString; mode=:none, selection=String[], write_artifacts::Bool=true, write_matpower_exports::Bool=false, performance_profile=nothing)
  results = Dict{String,Any}[]
  for outage in _selected_dtf_outages(case, mode, selection)
    matches = DTFImporter.find_outage_branch_indices(case, outage)
    length(matches) == 1 || throw(ArgumentError("missing_ambiguous_outage_branch_matching: " * DTFImporter.outage_match_diagnostic(case, outage, matches)))
    fresh_case = DTFImporter.read_dtf(case_path)
    net = DTFImporter.build_net(fresh_case; tap_changer_model = config.transformer.tap_changer_model)
    branch_index = only(matches)
    DTFImporter.apply_single_branch_outage!(net, branch_index)
    raw = _run_sparlectra(net = net, config = config, performance_profile = performance_profile, emit_output = false)
    label = DTFImporter.outage_label(outage)
    safe = replace(label, r"[^A-Za-z0-9_.-]+" => "_")
    artifacts = String[]
    if write_artifacts
      md = joinpath(output_path, "dtf_outage_$(outage.index)_summary.md")
      open(md, "w") do io
        println(io, "# DTF outage summary\n")
        println(io, "- label: ", label)
        println(io, "- matched branch index: ", branch_index)
        println(io, "- converged: ", raw.final_converged)
        println(io, "- iterations: ", raw.iterations)
        println(io, "- final mismatch: ", raw.final_mismatch)
      end
      push!(artifacts, basename(md))
      csv = joinpath(output_path, "dtf_outage_$(outage.index)_metrics.csv")
      open(csv, "w") do io
        println(io, "outage_label,matched_branch_index,branch_kind,from_bus,to_bus,converged,iterations,final_mismatch")
        println(io, join((_csv_field(label, ','), branch_index, outage.kind, _csv_field(outage.from, ','), _csv_field(outage.to, ','), raw.final_converged, raw.iterations, raw.final_mismatch), ','))
      end
      push!(artifacts, basename(csv))
    end
    if write_matpower_exports
      mpath = joinpath(output_path, "dtf_outage_$(outage.index)_$(safe).m")
      writeMatpowerCasefile(net, mpath; write_solution = config.matpower_export.write_solution)
      push!(artifacts, basename(mpath))
    end
    push!(results, Dict{String,Any}("outage_label" => label, "matched_branch_index" => branch_index, "branch_kind" => string(outage.kind), "from_bus" => outage.from, "to_bus" => outage.to, "converged" => raw.final_converged, "iterations" => raw.iterations, "final_mismatch" => raw.final_mismatch, "artifacts" => artifacts))
  end
  return results
end



"""
    run_sparlectra_api(; casefile, config_file, output_dir, config_overrides=Dict(),
                        performance_timing=:off, run_diagnostics=false,
                        detailed_result_csv=false,
                        detailed_result_csv_format="technical",
                        detailed_result_csv_semicolon=false) -> SparlectraApiResult

Run one MATPOWER power-flow case through a stable, non-interactive API contract.
The function validates GUI overrides, writes `effective_config.yaml`, delegates
the numerical work to [`run_sparlectra`](@ref), captures textual output in
`run.log`, writes `result.json`, discovers all generated files, and returns
structured status and artifact metadata. The input configuration template is
never modified. `performance_timing` may be `:off`, `:compact`, or `:full` and
writes a single-run `performance.log`; `run_diagnostics=true` captures existing
PowerFlow diagnostic printers in `diagnose.log`; and `detailed_result_csv=true`
writes Excel-friendly bus-voltage and branch-flow CSV artifacts. Optional
`detailed_result_csv_format` accepts `technical`, `excel_de`, or `excel_us`.
The legacy `detailed_result_csv_semicolon=true` maps to `excel_de` when the
explicit format is omitted. Artifact generation does not change PowerFlow run
success.
"""
function run_sparlectra_api(;
  casefile::AbstractString,
  config_file::AbstractString = DEFAULT_SPARLECTRA_CONFIG_PATH,
  output_dir::AbstractString,
  case_format = :auto,
  for002_reference_file = nothing,
  run_dtf_outages::Bool = false,
  dtf_outage_selection = String[],
  dtf_outage_selection_mode = :none,
  compare_for002_outages::Bool = false,
  write_outage_artifacts::Bool = true,
  write_outage_matpower_exports::Bool = false,
  matpower_export_requested::Bool = false,
  config_overrides::AbstractDict = Dict{String,Any}(),
  config_override_source::AbstractString = "explicit_api_request",
  performance_timing = :off,
  run_diagnostics::Bool = false,
  detailed_result_csv::Bool = false,
  detailed_result_csv_format = nothing,
  detailed_result_csv_semicolon::Bool = false,
)::SparlectraApiResult
  return _run_sparlectra_api(
    casefile = casefile,
    config_file = config_file,
    output_dir = output_dir,
    case_format = case_format,
    for002_reference_file = for002_reference_file,
    run_dtf_outages = run_dtf_outages,
    dtf_outage_selection = dtf_outage_selection,
    dtf_outage_selection_mode = dtf_outage_selection_mode,
    compare_for002_outages = compare_for002_outages,
    write_outage_artifacts = write_outage_artifacts,
    write_outage_matpower_exports = write_outage_matpower_exports,
    matpower_export_requested = matpower_export_requested,
    config_overrides = config_overrides,
    config_override_source = config_override_source,
    performance_timing = performance_timing,
    run_diagnostics = run_diagnostics,
    detailed_result_csv = detailed_result_csv,
    detailed_result_csv_format = detailed_result_csv_format,
    detailed_result_csv_semicolon = detailed_result_csv_semicolon,
    run_id = string(uuid4()),
  )
end

function _run_sparlectra_api(;
  casefile::AbstractString,
  config_file::AbstractString,
  output_dir::AbstractString,
  case_format = :auto,
  for002_reference_file = nothing,
  run_dtf_outages::Bool = false,
  dtf_outage_selection = String[],
  dtf_outage_selection_mode = :none,
  compare_for002_outages::Bool = false,
  write_outage_artifacts::Bool = true,
  write_outage_matpower_exports::Bool = false,
  matpower_export_requested::Bool = false,
  config_overrides::AbstractDict,
  config_override_source::AbstractString = "explicit_api_request",
  performance_timing = :off,
  run_diagnostics::Bool = false,
  detailed_result_csv::Bool = false,
  detailed_result_csv_format = nothing,
  detailed_result_csv_semicolon::Bool = false,
  phase_timings::AbstractDict = Dict{Symbol,Float64}(),
  run_id::String,
  cancellation_token = nothing,
  phase_callback = phase -> nothing,
  operation_callback = (event; fields...) -> nothing,
)::SparlectraApiResult
  total_start = time_ns()
  phases = Dict{Symbol,Float64}(Symbol(key) => Float64(value) for (key, value) in phase_timings)
  phase_recorder = PowerFlowPhaseTimingRecorder()
  emit_phase = phase -> begin
    _start_service_phase!(phase_recorder, String(phase))
    phase_callback(String(phase))
  end
  timing_mode = try
    _api_timing_mode(performance_timing)
  catch err
    output_path = abspath(output_dir)
    mkpath(output_path)
    logfile = joinpath(output_path, "run.log")
    result_file = joinpath(output_path, "result.json")
    touch(logfile)
    return _api_failure("invalid_performance_timing", sprint(showerror, err); run_id, casefile = abspath(casefile), config_file = abspath(config_file), output_dir = output_path, logfile, result_file)
  end
  output_path = abspath(output_dir)
  mkpath(output_path)
  case_path = abspath(casefile)
  config_path = abspath(config_file)
  logfile = joinpath(output_path, "run.log")
  result_file = joinpath(output_path, "result.json")
  touch(logfile)
  emit_phase("preparing_configuration")
  _check_powerflow_cancelled!(cancellation_token)

  isfile(config_path) || return _api_failure("config_file_not_found", "Configuration file not found: $(config_path)"; run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)

  config_start = time_ns()
  nested_overrides = try
    validate_gui_config_overrides(config_overrides)
  catch err
    return _api_failure("invalid_config_override", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)
  end

  config = nothing
  effective_raw = nothing
  try
    config, effective_raw = _load_api_config(config_path, nested_overrides)
  catch err
    return _api_failure("invalid_configuration", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)
  end
  config_sources = _config_source_report(config_path, nested_overrides, effective_raw; override_source = config_override_source)
  effective_config = joinpath(output_path, "effective_config.yaml")
  _write_yaml_file(effective_config, _effective_config_with_runtime_case(effective_raw, case_path, config; config_sources))
  _write_run_metadata_artifact(output_path; case_path = case_path)
  _check_powerflow_cancelled!(cancellation_token)
  phases[:api_config_build] = _api_elapsed_seconds(config_start)
  (isfile(case_path) || isdir(case_path)) || return _api_failure("casefile_not_found", "Case input not found: $(case_path)"; run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)

  requested_case_format = try
    _normalize_case_format(case_format)
  catch err
    return _api_failure("invalid_case_format", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)
  end
  detected_case_format = try
    _detect_case_format(case_path; requested = requested_case_format)
  catch err
    return _api_failure("ambiguous_case_format", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file, metadata = Dict("input_format" => String(requested_case_format), "input_format_detected" => "ambiguous"))
  end

  raw_result = nothing
  qlimit_metadata = _resolved_q_limit_runtime_options(config)
  qlimit_metadata["runtime_casefile"] = basename(case_path)
  qlimit_metadata["runtime_casefile_path"] = case_path
  matpower_metadata = _resolved_matpower_import_runtime_options(config)
  operation_callback("powerflow_effective_options"; run_id = run_id, case = basename(case_path), _metadata_kwargs(qlimit_metadata)..., _metadata_kwargs(matpower_metadata)...)
  # :enabled gates the fine-grained _perf_profile_time! sub-timings (issue
  # #288 finding 4: without it, phases like solver_finalization stay opaque
  # in performance.log). Tie it to the requested timing mode instead of
  # leaving it permanently off.
  api_performance_profile = Dict{Symbol,Any}(
    :enabled => _api_timing_mode(performance_timing) !== :off,
    :show_allocations => config.performance.show_allocations,
    :cancellation_check => () -> _check_powerflow_cancelled!(cancellation_token),
    :phase_callback => phase -> emit_phase(String(phase)),
    :output_dir => output_path,
    # Diagnostics-gated (fixed-reference self-check and Diagnose runs): the
    # rectangular solver records the full per-bus start-state mismatch, which
    # _write_start_residuals_artifact turns into self_check_residuals.csv.
    :capture_initial_residual_rows => run_diagnostics,
  )
  execution_start = time_ns()
  dtf_metadata = Dict{String,Any}()
  cgmes_metadata = Dict{String,Any}()
  if detected_case_format === :cgmes
    emit_phase("reading_cgmes_delivery")
    cgmes_cfg = config.cgmes
    # `path` may add further parts (typically the boundary set) to the case
    extra = filter(!isempty, strip.(split(cgmes_cfg.path, ';')))
    paths = String[case_path]
    for p in extra
      p == case_path || push!(paths, String(p))
    end
    # A delivery uploaded as a bare base case has no boundary set with it.
    # Rather than failing, pick up a boundary delivery sitting next to it in
    # the same directory (ENTSO-E names them "...Boundary..." or "..._BD_...").
    if length(paths) == 1
      neighbours = try
        readdir(dirname(abspath(case_path)))
      catch
        String[]
      end
      for n in neighbours
        occursin(r"(?i)boundary|_bd_|_bd\.", n) || continue
        cand = joinpath(dirname(abspath(case_path)), n)
        cand == abspath(case_path) && continue
        (isdir(cand) || endswith(lowercase(n), ".zip")) || continue
        push!(paths, cand)
        emit_phase("cgmes_boundary_autodetected")
        break
      end
    end
    cgmes_result = try
      importCGMES(
        path = length(paths) == 1 ? paths[1] : paths,
        baseMVA = cgmes_cfg.base_mva,
        require_boundary = cgmes_cfg.require_boundary,
        tap_control = cgmes_cfg.tap_control,
        machine_control = cgmes_cfg.machine_control,
        ignore_connected = cgmes_cfg.ignore_connected,
        vset_min_pu = cgmes_cfg.vset_min_pu,
        vset_max_pu = cgmes_cfg.vset_max_pu,
        multi_slack = cgmes_cfg.multi_slack,
        name = basename(case_path),
      )
    catch err
      message = sprint(showerror, err)
      # A failed import is exactly when the report matters most, so write a
      # diagnostic cgmes.log from what can still be read (summarizeCGMES works
      # on incomplete deliveries) before returning the failure.
      try
        open(joinpath(output_path, "cgmes.log"), "w") do io
          println(io, "# CGMES import report — IMPORT FAILED")
          println(io, "case:   ", case_path)
          println(io, "paths:  ", join(paths, "\n        "))
          println(io, "error:  ", message)
          println(io)
          summary = CGMESImporter.summarizeCGMES(path = length(paths) == 1 ? paths[1] : paths)
          println(io, "version: ", summary.version)
          println(io, "objects: ", summary.object_count, " in ", length(summary.class_histogram), " classes")
          println(io, "unresolved references: ", summary.unresolved_count, summary.boundary_missing_hint ? "  (topology references — boundary set likely missing)" : "")
          isempty(summary.unresolved_sample) || println(io, "  sample: ", join(summary.unresolved_sample, ", "))
          println(io)
          println(io, "## Files")
          for f in summary.files
            println(io, "  ", rpad(f.name, 60), " [", join(sort(String.(collect(f.profiles))), ","), "]", f.skipped ? "  SKIPPED: " * f.skip_reason : "  " * string(f.object_count) * " objects")
          end
          println(io)
          println(io, "## Class histogram")
          for (cls, n) in summary.class_histogram
            println(io, "  ", rpad(String(cls), 40), n)
          end
        end
      catch logerr
        @warn "could not write diagnostic cgmes.log" exception = logerr
      end
      is_boundary = occursin("boundary set missing", message) || occursin("unresolved topology references", message)
      reason = is_boundary ? "cgmes_boundary_missing" : "cgmes_import_error"
      is_boundary && (message = string(message, " — upload the boundary delivery into the same case directory (its name must contain \"Boundary\" or \"_BD_\"), add it to cgmes_import.path, or type cgmes:<alias> to fetch a ready-made ENTSO-E test set."))
      metadata = Dict("input_format" => String(requested_case_format), "input_format_detected" => "cgmes")
      return _api_execution_failure(reason, message; run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file, phase_recorder, performance_timing, total_start_ns = total_start, metadata = metadata)
    end
    # cgmes_import.start_values wins over power_flow(.start_mode).flatstart on
    # CGMES runs: :flat = synthetic flat start (default), :sv = the imported
    # SvVoltage state with every competing start-value machine forced off.
    cgmes_run_powerflow, start_overridden = _cgmes_start_values_powerflow(config.powerflow, cgmes_cfg.start_values)
    cgmes_run_config = _copy_sparlectra_with_powerflow(config, cgmes_run_powerflow)
    start_decision = string(
      "CGMES start values: ", cgmes_cfg.start_values,
      cgmes_cfg.start_values === :sv ? " (imported SvVoltage state; start-value machines forced off: start_projection, dc_seed_unconditional, start_current_iteration, apslf_start)" : " (synthetic flat start)",
      isempty(start_overridden) ? "" : string(" — overrides: ", join(start_overridden, ", ")),
    )
    # Dedicated cgmes.log: everything the importer saw and decided, next to
    # run.log and diagnose.log, so a CGMES run can be audited without
    # digging through the general run output.
    try
      open(joinpath(output_path, "cgmes.log"), "w") do io
        println(io, "# CGMES import report")
        println(io, "case:        ", case_path)
        println(io, "paths:       ", join(paths, "\n             "))
        println(io, "version:     ", cgmes_result.store.version)
        println(io, "tap control: ", cgmes_cfg.tap_control)
        println(io, "machine control: ", cgmes_cfg.machine_control)
        println(io, "base MVA:    ", cgmes_cfg.base_mva)
        println(io, start_decision)
        println(io)
        println(io, "## Files")
        for f in cgmes_result.store.files
          println(io, "  ", rpad(f.name, 60), " [", join(sort(String.(collect(f.profiles))), ","), "]", f.skipped ? "  SKIPPED: " * f.skip_reason : "  " * string(f.object_count) * " objects")
        end
        println(io)
        println(io, "## Class histogram")
        for (cls, mrids) in sort(collect(cgmes_result.store.byclass); by = p -> -length(p.second))
          println(io, "  ", rpad(String(cls), 40), length(mrids))
        end
        println(io)
        println(io, "## Network built")
        println(io, "  buses:      ", length(cgmes_result.net.nodeVec), " (", length(cgmes_result.net.isoNodes), " isolated)")
        println(io, "  branches:   ", length(cgmes_result.net.branchVec))
        println(io, "  links:      ", length(cgmes_result.net.linkVec))
        println(io, "  prosumers:  ", length(cgmes_result.net.prosumpsVec))
        println(io, "  shunts:     ", length(cgmes_result.net.shuntVec))
        println(io, "  slack bus:  ", cgmes_result.slack_bus)
        println(io, "  no-SV buses: ", length(cgmes_result.no_sv_buses), " (flat 1.0 pu / 0° start)")
        println(io)
        sc = cgmes_result.shortcircuit
        println(io, "## Short-circuit source data (read, not evaluated — issue #277)")
        # Per-attribute completeness: a future IEC 60909 evaluation needs to
        # know which inputs a delivery actually provides, not just the counts.
        printShortCircuitCoverage(io, sc)
        println(io)
        println(io, "## Importer messages (", length(cgmes_result.messages), ")")
        for m in cgmes_result.messages
          println(io, "  ", m)
        end
        if !isempty(cgmes_result.skipped_equipment)
          println(io)
          println(io, "## Skipped equipment: ", length(cgmes_result.skipped_equipment), " mRID(s)")
        end
      end
    catch err
      @warn "could not write cgmes.log" exception = err
    end
    emit_phase("building_sparlectra_net")
    # Per-bus SV coverage for the residual export (`has_sv` column): bus names
    # whose TopologicalNode carried no usable SvVoltage and therefore start at
    # the flat 1.0 pu / 0° fallback instead of the delivery's SV state.
    api_performance_profile[:cgmes_no_sv_buses] = Set(cgmes_result.no_sv_buses)
    raw_result = try
      open(logfile, "a") do io
        # run.log carries the narrative; the full importer report (files,
        # class histogram, every notice/warning) lives in cgmes.log — one
        # home per fact. Warnings still surface here because they mark
        # substituted values a reader of the narrative must not miss.
        println(io, "CGMES input (", length(paths), " path(s), ", length(cgmes_result.net.nodeVec), " buses)")
        println(io, start_decision)
        cgmes_warnings = [m for m in cgmes_result.messages if startswith(m, "warning:")]
        println(io, "  importer messages: ", length(cgmes_result.messages), " (", length(cgmes_warnings), " warning(s)) — full report in cgmes.log")
        for m in cgmes_warnings
          println(io, "  ", m)
        end
        _write_resolved_q_limit_options(io, qlimit_metadata)
        _capture_run_output(io; live = config.output.console_live) do
          _run_sparlectra(net = cgmes_result.net, config = cgmes_run_config, performance_profile = api_performance_profile)
        end
      end
    catch err
      err isa PowerFlowAborted && rethrow()
      return _api_execution_failure("execution_error", sprint(showerror, err, catch_backtrace()); run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file, phase_recorder, performance_timing, total_start_ns = total_start)
    end
    # Mandatory SV comparison (report, not a gate): every CGMES run checks the
    # solved state against the delivery's own SV solution. The solver mutated
    # cgmes_result.net in place (passed by reference to _run_sparlectra), so
    # compareWithSV reads the state this run actually produced. On a
    # non-converged run the numbers describe the diverged state — the status
    # tag keeps them from being read as model error.
    sv_compare = nothing
    sv_compare_error = nothing
    try
      sv_compare = CGMESImporter.compareWithSV(cgmes_result)
    catch err
      sv_compare_error = sprint(showerror, err)
    end
    sv_compare_status = sv_compare === nothing ? "unavailable" : (raw_result.final_converged ? "converged" : String(Symbol(raw_result.outcome)))
    if sv_compare !== nothing
      try
        open(joinpath(output_path, "sv_compare.csv"), "w") do io
          println(io, "bus,vm_pu,sv_vm_pu,dvm,va_deg,sv_va_deg,dva")
          for r in sv_compare.rows
            println(io, join((_csv_field(String(r.bus), ','), r.vm_pu, r.sv_vm_pu, r.dvm, r.va_deg, r.sv_va_deg, r.dva), ','))
          end
        end
        # The flow rows come for free with compareWithSV (computed in the same
        # call), so ship them as a second artifact.
        open(joinpath(output_path, "sv_compare_flows.csv"), "w") do io
          println(io, "kind,name,bus,sv_p_mw,sv_q_mvar,p_mw,q_mvar,dp_mw,dq_mvar")
          for r in sv_compare.flows.rows
            println(io, join((r.kind, _csv_field(String(r.name), ','), _csv_field(String(r.bus), ','), r.sv_p, r.sv_q, r.p, r.q, r.dp, r.dq), ','))
          end
        end
      catch err
        @warn "could not write SV comparison artifacts" exception = err
      end
    end
    try
      open(joinpath(output_path, "cgmes.log"), "a") do io
        println(io)
        println(io, "## SV comparison (run status: ", sv_compare_status, ", start values: ", cgmes_cfg.start_values, ")")
        if sv_compare === nothing
          println(io, "  unavailable", sv_compare_error === nothing ? "" : ": " * sv_compare_error)
        else
          println(io, "  voltages: n=", sv_compare.n, "  max|dvm|=", sv_compare.max_dvm, " pu  rms=", sv_compare.rms_dvm, "  max|dva|=", sv_compare.max_dva, "°  rms=", sv_compare.rms_dva, "°")
          println(io, "  flows:    n=", sv_compare.flows.n, "  max|dp|=", sv_compare.flows.max_dp, " MW  max|dq|=", sv_compare.flows.max_dq, " MVAr")
          println(io, "  artifacts: sv_compare.csv, sv_compare_flows.csv")
          raw_result.final_converged || println(io, "  note: run did not converge — deltas describe the non-converged state, not the imported model")
        end
      end
    catch err
      @warn "could not append SV comparison to cgmes.log" exception = err
    end
    cgmes_metadata = Dict{String,Any}(
      "input_format" => String(requested_case_format),
      "input_format_detected" => "cgmes",
      "cgmes_no_sv_buses" => length(cgmes_result.no_sv_buses),
      "cgmes_start_values" => String(cgmes_cfg.start_values),
      "cgmes_sv_compare_status" => sv_compare_status,
      "cgmes_sv_compare_n" => sv_compare === nothing ? 0 : sv_compare.n,
      "cgmes_sv_compare_max_dvm" => sv_compare === nothing ? NaN : sv_compare.max_dvm,
      "cgmes_sv_compare_rms_dvm" => sv_compare === nothing ? NaN : sv_compare.rms_dvm,
      "cgmes_sv_compare_max_dva" => sv_compare === nothing ? NaN : sv_compare.max_dva,
      "cgmes_sv_compare_rms_dva" => sv_compare === nothing ? NaN : sv_compare.rms_dva,
      "cgmes_import_used" => true,
      "cgmes_version" => cgmes_result.store.version,
      "cgmes_paths" => join(paths, ";"),
      "cgmes_buses" => length(cgmes_result.net.nodeVec),
      "cgmes_branches" => length(cgmes_result.net.branchVec),
      "cgmes_slack_bus" => cgmes_result.slack_bus,
      "cgmes_tap_control" => cgmes_cfg.tap_control,
      "cgmes_machine_control" => cgmes_cfg.machine_control,
      "cgmes_messages" => length(cgmes_result.messages),
      "cgmes_skipped_equipment" => length(cgmes_result.skipped_equipment),
    )
  elseif detected_case_format === :dtf_for001
    emit_phase("reading_dtf_for001_case")
    dtf_case = try
      _reject_dtf_dcline_like_content!(case_path)
      DTFImporter.read_dtf(case_path)
    catch err
      message = sprint(showerror, err)
      reason = occursin("unsupported_dtf_dc_line", message) ? "unsupported_dtf_dc_line" : "dtf_parse_error"
      metadata = Dict("input_format" => String(requested_case_format), "input_format_detected" => "dtf_for001", "native_dtf_import_used" => false, "unsupported_dcline_status" => reason)
      return _api_execution_failure(reason, message; run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file, phase_recorder, performance_timing, total_start_ns = total_start, metadata = metadata)
    end
    emit_phase("building_sparlectra_net")
    dtf_net = try
      DTFImporter.build_net(dtf_case; tap_changer_model = config.transformer.tap_changer_model)
    catch err
      return _api_execution_failure("dtf_build_net_error", sprint(showerror, err, catch_backtrace()); run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file, phase_recorder, performance_timing, total_start_ns = total_start)
    end
    raw_result = try
      open(logfile, "a") do io
        println(io, "Native DTF input (experimental/internal)")
        _write_resolved_q_limit_options(io, qlimit_metadata)
        _capture_run_output(io; live = config.output.console_live) do
          _run_sparlectra(net = dtf_net, config = config, performance_profile = api_performance_profile)
        end
      end
    catch err
      err isa PowerFlowAborted && rethrow()
      return _api_execution_failure("execution_error", sprint(showerror, err, catch_backtrace()); run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file, phase_recorder, performance_timing, total_start_ns = total_start)
    end
    matpower_export_file = nothing
    if matpower_export_requested
      matpower_export_file = joinpath(output_path, "dtf_native_matpower_export.m")
      try
        writeMatpowerCasefile(raw_result.net, matpower_export_file; write_solution = config.matpower_export.write_solution)
      catch err
        return _api_execution_failure("artifact_write_error", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file, phase_recorder, performance_timing, total_start_ns = total_start)
      end
    end
    dtf_metadata = _dtf_metadata(dtf_case, requested_case_format, detected_case_format; for002_file = for002_reference_file, run_dtf_outages = run_dtf_outages, matpower_export_requested = matpower_export_requested, matpower_export_file = matpower_export_file)
    _write_dtf_summary_artifacts(output_path, dtf_case, dtf_metadata)
    try
      for002_artifacts = _write_for002_placeholder_artifacts(output_path, for002_reference_file)
      dtf_metadata["for002_comparison_artifacts"] = for002_artifacts
      dtf_metadata["for002_comparison_status"] = isempty(for002_artifacts) ? "not_requested" : "reference_recorded"
    catch err
      return _api_execution_failure("invalid_for002_reference_file", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file, phase_recorder, performance_timing, total_start_ns = total_start, metadata = dtf_metadata)
    end
    dtf_outage_results = Dict{String,Any}[]
    if run_dtf_outages || Symbol(lowercase(strip(String(dtf_outage_selection_mode)))) != :none
      try
        mode = run_dtf_outages && Symbol(lowercase(strip(String(dtf_outage_selection_mode)))) == :none ? :all : dtf_outage_selection_mode
        dtf_outage_results = _run_dtf_outages(case_path, dtf_case, config, output_path; mode = mode, selection = dtf_outage_selection, write_artifacts = write_outage_artifacts, write_matpower_exports = write_outage_matpower_exports, performance_profile = api_performance_profile)
      catch err
        return _api_execution_failure("missing_ambiguous_outage_branch_matching", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file, phase_recorder, performance_timing, total_start_ns = total_start, metadata = dtf_metadata)
      end
    end
    dtf_metadata["dtf_outage_results"] = dtf_outage_results
    dtf_metadata["compare_for002_outages"] = compare_for002_outages
    operation_callback("dtf_for001_import_summary"; run_id = run_id, _metadata_kwargs(dtf_metadata)...)
    # Rejoin the normal artifact/finalization path with a native Net result.
  elseif detected_case_format !== :matpower
    return _api_failure("invalid_case_format", "Unsupported detected case format: $(detected_case_format)"; run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)
  end
  if detected_case_format === :matpower
    emit_phase("reading_matpower_case")
  try
    open(logfile, "a") do io
      _write_resolved_q_limit_options(io, qlimit_metadata)
      _capture_run_output(io; live = config.output.console_live) do
        cd(output_path) do
          raw_result = run_sparlectra(casefile = basename(case_path), path = dirname(case_path), config = config, performance_profile = api_performance_profile)
        end
      end
    end
  catch err
    err isa PowerFlowAborted && rethrow()
    if err isa MatpowerIO.UnsupportedMatpowerDclineError
      details = Dict{String,Any}(String(key) => value for (key, value) in err.details)
      details["solver_status"] = "aborted"
      details["service_status"] = "failed"
      details["run_status"] = "failed"
      details["last_phase"] = "reading_matpower_case"
      operation_callback("matpower_dcline_detected"; run_id = run_id, _metadata_kwargs(details)...)
      operation_callback("matpower_dcline_unsupported"; run_id = run_id, _metadata_kwargs(details)...)
      operation_callback("powerflow_aborted_unsupported_matpower_dcline"; run_id = run_id, _metadata_kwargs(details)...)
      _write_run_metadata_artifact(output_path; case_path = case_path, lifecycle = details)
      run_diagnostics && _write_execution_failure_diagnostics(joinpath(output_path, "diagnose.log"), "unsupported_matpower_dcline", err.message)
      return _api_execution_failure("unsupported_matpower_dcline", err.message; run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file, phase_recorder, performance_timing, total_start_ns = total_start, metadata = details)
    end
    island_message = _islandwise_failure_message(api_performance_profile)
    message = island_message === nothing ? sprint(showerror, err, catch_backtrace()) : island_message
    reason = get(phase_recorder.timings[phase_recorder.active_index === nothing ? length(phase_recorder.timings) : phase_recorder.active_index], "phase", "") == "loading_julia_case" ? "loading_julia_case_failed" : "execution_error"
    solver_elapsed_s = _solver_elapsed_from_profile(api_performance_profile)
    metadata = island_message === nothing ? Dict{String,Any}() : Dict{String,Any}(
      "solver_status" => "failed",
      "service_status" => "failed",
      "run_status" => "failed",
      "numerical_status" => "not_converged",
      "last_phase" => "solving_powerflow",
    )
    solver_elapsed_s === nothing || (metadata["solver_elapsed_s"] = solver_elapsed_s)
    run_diagnostics && _write_execution_failure_diagnostics(joinpath(output_path, "diagnose.log"), reason, message)
    return _api_execution_failure(reason, message; run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file, phase_recorder, performance_timing, total_start_ns = total_start, metadata = metadata)
  end
  end
  auto_profile_result = raw_result.performance_profile isa AbstractDict ? get(raw_result.performance_profile, :matpower_auto_profile_result, nothing) : nothing
  auto_profile_casefile = raw_result.performance_profile isa AbstractDict ? get(raw_result.performance_profile, :matpower_auto_profile_casefile, case_path) : case_path
  if auto_profile_result !== nothing
    config = auto_profile_result.config
    _update_effective_matpower_raw!(effective_raw, config)
    config_sources = _config_source_report(config_path, nested_overrides, effective_raw; override_source = config_override_source, auto_profile_result = auto_profile_result)
    _write_yaml_file(effective_config, _effective_config_with_runtime_case(effective_raw, case_path, config; config_sources))
    _write_matpower_auto_profile_artifact(output_path, auto_profile_result, config; casefile = String(auto_profile_casefile))
    operation_callback("powerflow_effective_options"; run_id = run_id, case = basename(case_path), _metadata_kwargs(qlimit_metadata)..., _metadata_kwargs(_resolved_matpower_import_runtime_options(config))...)
  end
  phases[:case_loading_network_solver] = _api_elapsed_seconds(execution_start)
  raw_result.solver_elapsed_s === nothing || (phases[:solver] = raw_result.solver_elapsed_s)
  if raw_result.performance_profile isa AbstractDict && haskey(raw_result.performance_profile, :postprocess_losses_and_flows)
    phases[:postprocessing] = Float64(raw_result.performance_profile[:postprocess_losses_and_flows])
  end

  artifact_start = time_ns()
  emit_phase("writing_diagnostics")
  _check_powerflow_cancelled!(cancellation_token)
  emit_phase("writing_artifacts")
  operation_callback("powerflow_lifecycle_status"; run_id = run_id, solver_status = "completed", artifact_status = "running", run_status = "finalizing", last_phase = "writing_artifacts")
  run_diagnostics && _write_powerflow_diagnostics(joinpath(output_path, "diagnose.log"), raw_result; mode = config.output.logfile_diagnostics)
  run_diagnostics && _write_start_residuals_artifact(output_path, api_performance_profile)
  q_limit_artifacts = raw_result.net !== nothing ? [_write_q_limit_log_artifact(output_path, raw_result, qlimit_metadata)] : String[]
  if (run_diagnostics || detailed_result_csv) && raw_result.net !== nothing
    append!(q_limit_artifacts, _write_q_limit_detail_artifacts(output_path, raw_result.net; format = "technical"))
  end
  matpower_dcline_artifact = raw_result.net === nothing ? nothing : _write_matpower_dcline_artifact(output_path, raw_result.net; format = "technical")
  if matpower_dcline_artifact !== nothing
    operation_callback("matpower_dcline_pf_injections_imported"; run_id = run_id, active_dcline_count = length(raw_result.net.matpowerDclineMetadata), artifact = matpower_dcline_artifact)
    open(logfile, "a") do io
      println(io, "MATPOWER active mpc.dcline rows imported through toggle_dcline-compatible PF injections.")
      println(io, "MATPOWER DC-line artifact: ", matpower_dcline_artifact)
    end
  end
  csv_artifacts = String[]
  csv_export_error = nothing
  csv_export_status = detailed_result_csv ? "pending" : "disabled"
  csv_export_skip_reason = nothing
  csv_format = nothing
  csv_timing_metadata = Dict{Symbol,Any}(:progress_callback => ((event; fields...) -> operation_callback(event; run_id = run_id, fields...)))
  if detailed_result_csv
    csv_format_name = detailed_result_csv_format === nothing ? (detailed_result_csv_semicolon ? "excel_de" : "technical") : String(detailed_result_csv_format)
    csv_format = try
      _resolve_detailed_csv_format(csv_format_name)
    catch err
      return _api_failure("invalid_detailed_result_csv_format", sprint(showerror, err); run_id = run_id, casefile = case_path, config_file = config_path, output_dir = output_path, logfile = logfile, result_file = result_file)
    end
  end
  if !detailed_result_csv
    csv_export_skip_reason = "csv_disabled"
  elseif raw_result.net === nothing
    csv_export_status = "skipped"
    csv_export_skip_reason = "no_network_available"
  else
    emit_phase("writing_csv_artifacts")
    operation_callback("powerflow_lifecycle_status"; run_id = run_id, solver_status = "completed", artifact_status = "running", run_status = "finalizing", last_phase = "writing_csv_artifacts")
    try
      # CSV export is delegated so API orchestration only chooses when exports
      # happen; the helper preserves filenames, schemas, progress events, and
      # partial-export behavior expected by the Web UI.
      csv_artifacts = _write_detailed_result_csv(output_path, raw_result; format = csv_format.name, config, abort_checker = () -> _check_powerflow_cancelled!(cancellation_token), timing_metadata = csv_timing_metadata)
      csv_export_status = haskey(csv_timing_metadata, :partial_error) ? "partial" : (raw_result.final_converged && raw_result.solution_available ? "exported" : "exported_diagnostic")
    catch err
      csv_export_error = sprint(showerror, err)
      csv_export_status = "failed"
      operation_callback("detailed_result_csv_export_failed"; run_id = run_id, error = csv_export_error, final_converged = raw_result.final_converged, solution_available = raw_result.solution_available, has_network = raw_result.net !== nothing, csv_format = csv_format.name, status = "failed")
      open(logfile, "a") do io
        println(io, "Detailed result CSV export failed: ", csv_export_error)
      end
    end
  end
  if detailed_result_csv && csv_export_status == "skipped"
    operation_callback("detailed_result_csv_export_skipped"; run_id = run_id, reason = csv_export_skip_reason, final_converged = raw_result.final_converged, solution_available = raw_result.solution_available, has_network = raw_result.net !== nothing, csv_format = csv_format.name, status = "skipped")
    open(logfile, "a") do io
      println(io, "Detailed result CSV export skipped: ", csv_export_skip_reason)
      println(io, "  final_converged: ", raw_result.final_converged)
      println(io, "  solution_available: ", raw_result.solution_available)
      println(io, "  has_network: ", raw_result.net !== nothing)
      println(io, "  csv_format: ", csv_format.name)
    end
  end
  _check_powerflow_cancelled!(cancellation_token)
  emit_phase("finalizing_success")
  _finalize_service_timings!(phase_recorder, total_start; status = "completed")
  phases[:artifact_writing] = _api_elapsed_seconds(artifact_start)
  phases[:total] = _api_elapsed_seconds(total_start)
  open(logfile, "a") do io
    println(io)
    println(io, "API run summary")
    println(io, "===============")
    _write_api_timing_summary(io, raw_result, config, phases)
    println(io, "Case file format: ", lowercase(splitext(case_path)[2]))
    println(io, "Case file size: ", filesize(case_path), " bytes")
    println(io)
    _write_service_phase_summary(io, phase_recorder.timings)
    _write_large_case_timing_summary(io, case_path, phase_recorder.timings, raw_result)
    if !isempty(q_limit_artifacts)
      println(io, "PV Q-limit details")
      println(io, "------------------")
      println(io, "full details        : ", Q_LIMIT_LOG_ARTIFACT)
      println(io, "full initial limits : ", "q_limit_initial_limits.csv" in q_limit_artifacts ? "q_limit_initial_limits.csv" : "not generated")
      println(io, "full events         : ", "q_limit_events.csv" in q_limit_artifacts ? "q_limit_events.csv" : "not generated")
      println(io, "classic outer loop  : ", "q_limit_classic_outer_loop.csv" in q_limit_artifacts ? "q_limit_classic_outer_loop.csv" : "not generated")
      println(io)
    end
    if config.output.logfile_results === :full
      println(io, "Full run details")
      println(io, "----------------")
      println(io, "casefile: ", case_path)
      println(io, "runtime_casefile: ", basename(case_path))
      println(io, "config_file: ", config_path)
      println(io, "output_dir: ", output_path)
      println(io, "diagnostics_artifact: ", run_diagnostics ? "diagnose.log" : "disabled")
      println(io, "performance_artifact: ", timing_mode === :off ? "disabled" : "performance.log")
      # Q-limit settings are NOT repeated here — the "Resolved Q-limit
      # options" block near the top of this file is their one home.
      println(io, "q_limit_detail_artifacts: ", isempty(q_limit_artifacts) ? "none" : join(q_limit_artifacts, ", "))
      println(io, "matpower_dcline_artifact: ", matpower_dcline_artifact === nothing ? "none" : matpower_dcline_artifact)
      println(io, "detailed_result_csv_status: ", csv_export_status)
      csv_export_skip_reason === nothing || println(io, "detailed_result_csv_skip_reason: ", csv_export_skip_reason)
      # The per-setting CSV lines only carry information when the export runs;
      # a disabled export collapses to the two status lines above instead of
      # fifteen lines reading "disabled".
      if detailed_result_csv
        if csv_export_skip_reason === nothing
          println(io, "detailed_result_csv_solution_quality: ", _csv_solution_quality(raw_result))
          raw_result.final_converged && raw_result.solution_available || println(io, "detailed_result_csv_warning: values are from the last Newton iterate and are not a converged power-flow solution")
        end
        println(io, "detailed_result_csv_artifacts: ", isempty(csv_artifacts) ? csv_export_status : join(csv_artifacts, ", "))
        println(io, "detailed_result_csv_format: ", csv_format.name)
        println(io, "detailed_result_csv_exporter: ", raw_result.net !== nothing ? _select_detailed_csv_exporter(raw_result.net; config) : "unknown")
        println(io, "detailed_result_csv_write_mode: ", config.output.detailed_result_csv_write_mode)
        println(io, "detailed_result_csv_actual_write_mode: ", raw_result.net !== nothing && _select_detailed_csv_exporter(raw_result.net; config) === :direct ? "streaming" : config.output.detailed_result_csv_write_mode)
        raw_result.net === nothing || println(io, "detailed_result_csv_bus_rows: ", length(raw_result.net.nodeVec))
        raw_result.net === nothing || println(io, "detailed_result_csv_branch_rows: ", length(raw_result.net.branchVec))
        println(io, "detailed_result_csv_delimiter: ", csv_format.delimiter == ';' ? "semicolon" : "comma")
        println(io, "detailed_result_csv_decimal_separator: ", csv_format.decimal_separator == ',' ? "comma" : "dot")
        println(io, "detailed_result_csv_thousands_separator: ", csv_format.thousands_separator == "," ? "comma" : csv_format.thousands_separator == "." ? "dot" : "none")
      end
      if haskey(csv_timing_metadata, :exporter)
        println(io, "Detailed CSV exporter       : ", csv_timing_metadata[:exporter])
        println(io, "CSV write mode              : ", csv_timing_metadata[:write_mode])
        println(io, "CSV format                  : ", csv_timing_metadata[:format])
        println(io, "Bus CSV rows                : ", csv_timing_metadata[:bus_rows])
        println(io, "Branch CSV rows             : ", csv_timing_metadata[:branch_rows])
        println(io, "Bus control label cache time: ", @sprintf("%.6f s", csv_timing_metadata[:control_label_cache_s]))
        println(io, "Bus CSV export time         : ", @sprintf("%.6f s", csv_timing_metadata[:bus_export_s]))
        println(io, "Branch CSV export time      : ", @sprintf("%.6f s", csv_timing_metadata[:branch_export_s]))
        println(io, "CSV total export time       : ", @sprintf("%.6f s", csv_timing_metadata[:total_export_s]))
        println(io, "CSV formatting mode         : ", csv_timing_metadata[:formatting_mode])
      end
      csv_export_error === nothing || println(io, "detailed_result_csv_error: ", csv_export_error)
      haskey(csv_timing_metadata, :partial_error) && println(io, "detailed_result_csv_error: ", csv_timing_metadata[:partial_error])
      println(io)
      # The full configuration dump used to be repeated here (~150 lines); its
      # one home is the effective_config.yaml artifact written before the run.
      println(io, "effective configuration: effective_config.yaml")
      println(io)
      println(io, "Available status diagnostics")
      println(io, "----------------------------")
      for key in keys(raw_result.diagnostics)
        println(io, key, ": ", getproperty(raw_result.diagnostics, key))
      end
      println(io)
    end
    _write_numerical_outcome_summary(io, raw_result)
  end
  timing_mode === :off || _write_performance_log(joinpath(output_path, "performance.log"), timing_mode, phases, raw_result, phase_recorder.timings)
  _check_powerflow_cancelled!(cancellation_token)

  numerical_success = raw_result.final_converged && raw_result.solution_available
  mismatch = isfinite(raw_result.final_mismatch) ? raw_result.final_mismatch : nothing
  final_outcome = _final_outcome_payload(raw_result)
  # Lifecycle metadata is a service/API contract: Last Errors, run listings, and
  # clients consume these status fields even when the numerical result is
  # non-converged but artifact finalization succeeded.
  final_metadata = _build_success_lifecycle_metadata(
    raw_result,
    config;
    numerical_success,
    final_outcome,
    csv_export_status,
    csv_export_skip_reason,
    csv_export_error,
    csv_artifacts,
    detailed_result_csv,
    config_overrides,
    config_override_source,
    casefile,
    config_file,
    performance_timing,
    run_diagnostics,
    detailed_result_csv_format,
    qlimit_metadata,
    csv_timing_metadata,
  )
  isempty(dtf_metadata) || merge!(final_metadata, dtf_metadata)
  isempty(cgmes_metadata) || merge!(final_metadata, cgmes_metadata)
  _write_run_metadata_artifact(output_path; case_path = case_path, lifecycle = final_metadata)
  message = numerical_success ? "PowerFlow run completed." : raw_result.reason_text
  result = _api_result(
    run_id = run_id,
    status = numerical_success ? :succeeded : :not_converged,
    success = numerical_success,
    converged = raw_result.numerical_converged,
    solution_available = raw_result.solution_available,
    iterations = raw_result.iterations,
    final_mismatch = mismatch,
    reason = String(raw_result.reason),
    message = message,
    casefile = case_path,
    config_file = config_path,
    output_dir = output_path,
    logfile = logfile,
    result_file = result_file,
    service_phase_timings = phase_recorder.timings,
    metadata = final_metadata,
    raw_result = raw_result,
  )
  return _finalize_api_result(result)
end
