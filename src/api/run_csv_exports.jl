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

# API detailed CSV export helpers.
#
# This file owns CSV export orchestration and low-allocation writing details.
# Keep schemas, filenames, progress events, and partial-export behavior stable
# because API clients and the Web UI consume those contracts.

const _DETAILED_CSV_BUFFER_INITIAL_BYTES_DEFAULT = 8 * 1024 * 1024
const _DETAILED_CSV_BUFFER_MAX_BYTES_DEFAULT = 64 * 1024 * 1024
const _DETAILED_CSV_STREAMING_THRESHOLD_ROWS_DEFAULT = 100_000
const _DETAILED_CSV_DIRECT_THRESHOLD_BUSES_DEFAULT = 10_000

function _csv_write_options(config = nothing)::NamedTuple
  output = config isa SparlectraConfig ? config.output : config isa OutputConfig ? config : nothing
  mode = output === nothing ? :auto : output.detailed_result_csv_write_mode
  mode in OUTPUT_DETAILED_RESULT_CSV_WRITE_MODE_VALUES || throw(ArgumentError("Unsupported detailed_result_csv_write_mode \"$(mode)\". Expected auto, buffered, or streaming."))
  initial_bytes = output === nothing ? _DETAILED_CSV_BUFFER_INITIAL_BYTES_DEFAULT : output.detailed_result_csv_buffer_initial_bytes
  max_bytes = output === nothing ? _DETAILED_CSV_BUFFER_MAX_BYTES_DEFAULT : output.detailed_result_csv_buffer_max_bytes
  threshold_rows = output === nothing ? _DETAILED_CSV_STREAMING_THRESHOLD_ROWS_DEFAULT : output.detailed_result_csv_streaming_threshold_rows
  initial_bytes = initial_bytes < 0 ? _DETAILED_CSV_BUFFER_INITIAL_BYTES_DEFAULT : initial_bytes
  max_bytes = max_bytes <= 0 ? _DETAILED_CSV_BUFFER_MAX_BYTES_DEFAULT : max_bytes
  threshold_rows = threshold_rows <= 0 ? _DETAILED_CSV_STREAMING_THRESHOLD_ROWS_DEFAULT : threshold_rows
  return (mode = mode, initial_bytes = initial_bytes, max_bytes = max_bytes, threshold_rows = threshold_rows)
end

function _estimated_namedtuple_csv_bytes(rows::AbstractVector, columns)::Int
  return length(join(String.(columns), ',')) + 1 + length(rows) * max(1, length(columns)) * 24
end

function _write_namedtuple_csv_row(io::IO, row, columns, delimiter::Char, resolved_format)
  println(io, join((_csv_field(getproperty(row, column), delimiter, resolved_format) for column in columns), delimiter))
  return nothing
end

function _write_csv_row_values(io::IO, values, delimiter::Char, resolved_format)
  first = true
  for value in values
    first || print(io, delimiter)
    print(io, _csv_field(value, delimiter, resolved_format))
    first = false
  end
  println(io)
  return nothing
end

function _write_namedtuple_csv_buffered(path::AbstractString, rows::AbstractVector, columns, delimiter::Char, resolved_format; initial_bytes::Integer = _DETAILED_CSV_BUFFER_INITIAL_BYTES_DEFAULT)
  buffer = IOBuffer(; sizehint = max(0, Int(initial_bytes)))
  println(buffer, join(String.(columns), delimiter))
  for row in rows
    _write_namedtuple_csv_row(buffer, row, columns, delimiter, resolved_format)
  end
  open(path, "w") do io
    write(io, take!(buffer))
  end
  return path
end

function _write_namedtuple_csv_streaming(path::AbstractString, rows::AbstractVector, columns, delimiter::Char, resolved_format)
  open(path, "w") do io
    println(io, join(String.(columns), delimiter))
    for row in rows
      _write_namedtuple_csv_row(io, row, columns, delimiter, resolved_format)
    end
  end
  return path
end

function _select_namedtuple_csv_write_mode(rows::AbstractVector, columns; config = nothing, estimated_rows::Integer = length(rows))::Symbol
  options = _csv_write_options(config)
  options.mode === :buffered && return :buffered
  options.mode === :streaming && return :streaming
  estimated_rows > options.threshold_rows && return :streaming
  _estimated_namedtuple_csv_bytes(rows, columns) > options.max_bytes && return :streaming
  return :buffered
end

function _write_namedtuple_csv(path::AbstractString, rows::AbstractVector, columns; delimiter::Char = ',', format = nothing, config = nothing, estimated_rows::Integer = length(rows))
  delimiter in (',', ';') || throw(ArgumentError("CSV delimiter must be ',' or ';'."))
  resolved_format = format === nothing ? (name = "custom", delimiter = delimiter, decimal_separator = '.', thousands_separator = "") : _resolve_detailed_csv_format(format)
  resolved_format.delimiter == delimiter || throw(ArgumentError("CSV delimiter does not match detailed CSV format $(resolved_format.name)."))
  options = _csv_write_options(config)
  mode = _select_namedtuple_csv_write_mode(rows, columns; config, estimated_rows)
  mode === :streaming && return _write_namedtuple_csv_streaming(path, rows, columns, delimiter, resolved_format)
  return _write_namedtuple_csv_buffered(path, rows, columns, delimiter, resolved_format; initial_bytes = min(options.initial_bytes, options.max_bytes))
end

function _detailed_csv_export_options(config = nothing)::NamedTuple
  output = config isa SparlectraConfig ? config.output : config isa OutputConfig ? config : nothing
  exporter = output === nothing ? :auto : output.detailed_result_csv_exporter
  threshold_buses = output === nothing ? _DETAILED_CSV_DIRECT_THRESHOLD_BUSES_DEFAULT : output.detailed_result_csv_direct_threshold_buses
  exporter in OUTPUT_DETAILED_RESULT_CSV_EXPORTER_VALUES || throw(ArgumentError("Unsupported detailed_result_csv_exporter \"$(exporter)\". Expected auto, report, or direct."))
  return (exporter = exporter, threshold_buses = threshold_buses <= 0 ? _DETAILED_CSV_DIRECT_THRESHOLD_BUSES_DEFAULT : threshold_buses)
end

function _select_detailed_csv_exporter(net::Net; config = nothing)::Symbol
  options = _detailed_csv_export_options(config)
  options.exporter === :report && return :report
  options.exporter === :direct && return :direct
  return length(net.nodeVec) >= options.threshold_buses ? :direct : :report
end

function _complex_voltage_rows(node_rows::AbstractVector, format)
  return [
    begin
      angle = deg2rad(Float64(row.va_deg))
      v_re = Float64(row.vm_pu) * cos(angle)
      v_im = Float64(row.vm_pu) * sin(angle)
      merge(row, (v_re = round(v_re; sigdigits = 10), v_im = round(v_im; sigdigits = 10), v_complex = string(_format_csv_number(v_re, format), signbit(v_im) ? " - j" : " + j", _format_csv_number(abs(v_im), format))))
    end for row in node_rows
  ]
end

function _branch_csv_values(br)
  f = br.fBranchFlow
  t = br.tBranchFlow
  p_from = isnothing(f) || isnothing(f.pFlow) ? 0.0 : f.pFlow
  q_from = isnothing(f) || isnothing(f.qFlow) ? 0.0 : f.qFlow
  p_to = isnothing(t) || isnothing(t.pFlow) ? 0.0 : t.pFlow
  q_to = isnothing(t) || isnothing(t.qFlow) ? 0.0 : t.qFlow
  rated = isnothing(br.sn_MVA) ? 0.0 : br.sn_MVA
  overloaded = rated > 0.0 && max(abs(p_from), abs(p_to)) > rated
  return (br.comp.cName, br.branchIdx, br.fromBus, br.toBus, br.status, p_from, q_from, p_to, q_to, _default0(br.pLosses), _default0(br.qLosses), rated, overloaded)
end

function _check_csv_abort(abort_checker, row_count::Integer)
  if abort_checker !== nothing && row_count % 1000 == 0
    abort_checker()
  end
  return nothing
end

function _write_detailed_result_csv_artifacts_direct!(artifacts::Vector{String}, output_path::AbstractString, net::Net, result::SparlectraRunResult, config; format = "technical", abort_checker = nothing, timing_metadata = nothing)::Vector{String}
  resolved_format = _resolve_detailed_csv_format(format)
  runtime_format = CsvFormatRuntime(resolved_format)
  bus_columns = (:bus, :bus_name, :type, :vm_pu, :va_deg, :vn_kV, :v_re, :v_im, :v_complex, :v_kV, :p_gen_MW, :q_gen_MVar, :p_load_MW, :q_load_MVar, :q_limit_hit, :control)
  branch_columns = (:branch, :branch_index, :from_bus, :to_bus, :status, :p_from_MW, :q_from_MVar, :p_to_MW, :q_to_MVar, :p_loss_MW, :q_loss_MVar, :rated_MVA, :overloaded)
  busNameByIdx = _bus_name_by_idx(net)
  power_components = _bus_power_component_cache(net)
  cache_start = time_ns()
  control_labels = _bus_control_flag_cache(net)
  control_cache_elapsed = _api_elapsed_seconds(cache_start)
  nodes_sorted = sort(net.nodeVec, by = x -> x.busIdx)
  abort_checker === nothing || abort_checker()
  total_start = time_ns()
  bus_rows = 0
  branch_rows = 0
  current_artifact = artifacts[1]
  current_rows = Ref(0)
  progress_callback = timing_metadata isa AbstractDict ? get(timing_metadata, :progress_callback, nothing) : nothing
  _emit_csv_progress(event::AbstractString; fields...) = progress_callback === nothing ? nothing : progress_callback(event; fields...)
  bus_elapsed = 0.0
  branch_elapsed = 0.0
  _emit_csv_progress(
    "detailed_result_csv_export_started";
    exporter = "direct",
    write_mode = "streaming",
    csv_format = resolved_format.name,
    bus_rows = length(nodes_sorted),
    branch_rows = length(net.branchVec),
  )
  bus_start = time_ns()

  try
    open(joinpath(output_path, artifacts[1]), "w") do io
      println(io, join(String.(bus_columns), resolved_format.delimiter))
      row_count = 0
      for n in nodes_sorted
        row_count += 1
        current_rows[] = row_count
      angle = deg2rad(Float64(n._va_deg))
      v_re = Float64(n._vm_pu) * cos(angle)
      v_im = Float64(n._vm_pu) * sin(angle)
      p_gen, q_gen, p_load, q_load = _bus_power_components(power_components, n.busIdx)
      v_re_rounded = round(v_re; sigdigits = 10)
      v_im_rounded = round(v_im; sigdigits = 10)
      v_complex = string(_format_csv_number(v_re, runtime_format), signbit(v_im) ? " - j" : " + j", _format_csv_number(abs(v_im), runtime_format))
      write_csv_row_direct!(
        io,
        resolved_format.delimiter,
        runtime_format,
        n.busIdx,
        get(busNameByIdx, n.busIdx, n.comp.cName),
        toString(n._nodeType),
        n._vm_pu,
        n._va_deg,
        n.comp.cVN,
        v_re_rounded,
        v_im_rounded,
        v_complex,
        n.comp.cVN * n._vm_pu,
        p_gen,
        q_gen,
        p_load,
        q_load,
        haskey(net.qLimitEvents, n.busIdx),
        _cached_control_label(control_labels, n.busIdx),
      )
      _check_csv_abort(abort_checker, row_count)
      end
      bus_rows = row_count
    end
    bus_elapsed = _api_elapsed_seconds(bus_start)
    _emit_csv_progress("detailed_result_csv_file_written"; artifact = artifacts[1], rows = bus_rows, elapsed_s = bus_elapsed)

    try
      current_artifact = artifacts[2]
      current_rows[] = 0
      branch_start = time_ns()
      open(joinpath(output_path, artifacts[2]), "w") do io
        println(io, join(String.(branch_columns), resolved_format.delimiter))
        row_count = 0
        for br in net.branchVec
          row_count += 1
          current_rows[] = row_count
          f = br.fBranchFlow
          t = br.tBranchFlow
          p_from = isnothing(f) || isnothing(f.pFlow) ? 0.0 : f.pFlow
          q_from = isnothing(f) || isnothing(f.qFlow) ? 0.0 : f.qFlow
          p_to = isnothing(t) || isnothing(t.pFlow) ? 0.0 : t.pFlow
          q_to = isnothing(t) || isnothing(t.qFlow) ? 0.0 : t.qFlow
          rated = isnothing(br.sn_MVA) ? 0.0 : br.sn_MVA
          overloaded = rated > 0.0 && max(abs(p_from), abs(p_to)) > rated
          write_csv_row_direct!(
            io,
            resolved_format.delimiter,
            runtime_format,
            br.comp.cName,
            br.branchIdx,
            br.fromBus,
            br.toBus,
            br.status,
            p_from,
            q_from,
            p_to,
            q_to,
            _default0(br.pLosses),
            _default0(br.qLosses),
            rated,
            overloaded,
          )
          _check_csv_abort(abort_checker, row_count)
        end
        branch_rows = row_count
      end
      branch_elapsed = _api_elapsed_seconds(branch_start)
      _emit_csv_progress("detailed_result_csv_file_written"; artifact = artifacts[2], rows = branch_rows, elapsed_s = branch_elapsed)
    catch err
      timing_metadata !== nothing && (timing_metadata[:partial_error] = "branch_flows.csv failed: $(sprint(showerror, err))")
      _emit_csv_progress("detailed_result_csv_export_partial"; artifact = artifacts[2], rows_written = current_rows[], elapsed_s = _api_elapsed_seconds(total_start), error = timing_metadata === nothing ? sprint(showerror, err) : timing_metadata[:partial_error])
      return [artifacts[1]]
    end
  catch err
    _emit_csv_progress("detailed_result_csv_export_aborted"; artifact = current_artifact, rows_written = current_rows[], elapsed_s = _api_elapsed_seconds(total_start))
    rethrow()
  end
  if timing_metadata !== nothing
    timing_metadata[:exporter] = :direct
    timing_metadata[:write_mode] = :streaming
    timing_metadata[:format] = resolved_format.name
    timing_metadata[:formatting_mode] = :direct_low_allocation
    timing_metadata[:control_label_cache_s] = control_cache_elapsed
    timing_metadata[:bus_rows] = bus_rows
    timing_metadata[:branch_rows] = branch_rows
    timing_metadata[:bus_export_s] = bus_elapsed
    timing_metadata[:branch_export_s] = branch_elapsed
    timing_metadata[:total_export_s] = _api_elapsed_seconds(total_start)
  end
  return artifacts
end

function _write_detailed_result_csv(output_path::AbstractString, result::SparlectraRunResult; format = "technical", config = nothing, abort_checker = nothing, timing_metadata = nothing)::Vector{String}
  resolved_format = _resolve_detailed_csv_format(format)
  result.net === nothing && throw(ArgumentError("PowerFlow result does not contain a network for detailed CSV export."))
  artifacts = ["bus_voltages_complex.csv", "branch_flows.csv"]
  exporter = _select_detailed_csv_exporter(result.net; config)
  exporter === :direct && return _write_detailed_result_csv_artifacts_direct!(artifacts, output_path, result.net, result, config; format = resolved_format.name, abort_checker, timing_metadata)
  report = buildACPFlowReport(result.net; ct = result.elapsed_s, ite = result.iterations, converged = result.final_converged)
  bus_columns = (:bus, :bus_name, :type, :vm_pu, :va_deg, :vn_kV, :v_re, :v_im, :v_complex, :v_kV, :p_gen_MW, :q_gen_MVar, :p_load_MW, :q_load_MVar, :q_limit_hit, :control)
  branch_columns = (:branch, :branch_index, :from_bus, :to_bus, :status, :p_from_MW, :q_from_MVar, :p_to_MW, :q_to_MVar, :p_loss_MW, :q_loss_MVar, :rated_MVA, :overloaded)
  estimated_rows = length(report.nodes) + length(report.branches)
  _write_namedtuple_csv(joinpath(output_path, artifacts[1]), _complex_voltage_rows(report.nodes, resolved_format), bus_columns; delimiter = resolved_format.delimiter, format = resolved_format.name, config, estimated_rows)
  try
    _write_namedtuple_csv(joinpath(output_path, artifacts[2]), report.branches, branch_columns; delimiter = resolved_format.delimiter, format = resolved_format.name, config, estimated_rows)
  catch err
    timing_metadata !== nothing && (timing_metadata[:partial_error] = "branch_flows.csv failed: $(sprint(showerror, err))")
    return [artifacts[1]]
  end
  return artifacts
end

