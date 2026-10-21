# Copyright 2023-2026 Udo Schmitz
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

# file: app/src/api/run_qsts_service.jl
# purpose: service backend of the Web UI "QSTS" button: an SCF time series
#          (a zip recognized by its content) runs through runQsts in the
#          mode of the request, with the controllers of the companion file
#          <series>.control.yaml. The run writes the files of
#          writeQstsResults, the element view qsts_view.json (read by the
#          result page and the chart route), one SVG per element under
#          qsts_plots/ and qsts_overview.html. Mirrors
#          run_contingency_service.jl: a mode key on POST /powerflow/run,
#          artifacts collected from the run directory.

const _QSTS_SERVICE_MODES = ["sequential", "snapshots"]
const _QSTS_VIEW_FILE = "qsts_view.json"
const _QSTS_VIEW_FORMAT = "sparlectra-qsts-view"
const _QSTS_PLOT_DIR = "qsts_plots"
const _QSTS_OVERVIEW_FILE = "qsts_overview.html"
# the size of every chart the service and the chart route draw
const _QSTS_CHART_WIDTH = 720
const _QSTS_CHART_HEIGHT = 260

# The companion control file of a series: `<series>.control.yaml` next to
# it, the series name without `.qsts.zip` (any other name without its last
# extension). It is read for QSTS runs only, never as a case configuration:
# a plain power flow on the series runs the base case without controllers.
function _qsts_control_file_path(case_path::AbstractString)::String
  path = String(case_path)
  name = basename(path)
  stem = endswith(lowercase(name), ".qsts.zip") ? name[1:(end - 9)] : first(splitext(name))
  return joinpath(dirname(path), string(stem, ".control.yaml"))
end

# the series of a result as JSON values: NaN (no value) becomes null
_qsts_view_values(v::AbstractVector{Float64}) = Any[isnan(x) ? nothing : x for x in v]

# the quantities of one element that hold a value in at least one step, in
# the order of the quantity list
function _qsts_view_quantities(r::QstsResult, id::AbstractString, list::Vector{String})
  names = String[]
  values = Dict{String,Any}()
  for q in list
    v = qstsValues(r, id, q)
    any(!isnan, v) || continue
    push!(names, q)
    values[q] = _qsts_view_values(v)
  end
  return names, values
end

# a file name for the SVG of one element: the address with every character
# outside [A-Za-z0-9._-] replaced, prefixed by the kind, unique in the run
function _qsts_plot_name!(used::Set{String}, kind::AbstractString, id::AbstractString)::String
  stem = string(kind, "_", replace(String(id), r"[^A-Za-z0-9._-]" => "_"))
  name = string(stem, ".svg")
  k = 1
  while name in used
    k += 1
    name = string(stem, "_", k, ".svg")
  end
  push!(used, name)
  return string(_QSTS_PLOT_DIR, "/", name)
end

"""
    _qsts_view_dict(r::QstsResult) -> Dict{String,Any}

The element view of a QSTS run, written as `qsts_view.json`: the step
labels, the bands of `vm_pu` and `loading_pct`, and one entry per element
in the order run (step quantities, id empty), controllers, buses (star
points left out), branches. Each entry carries id, kind, name, its main
quantity (`position`, `vm_pu`, `loading_pct`, for the run
`total_p_loss_mw`), the path of its SVG (not for the run), the quantity
names in list order and their values per step, `null` for a step without
a value. A quantity without any value is left out. An app artifact, not
part of the result format 1.0: the result page and the chart route read
it, nobody else.
"""
function _qsts_view_dict(r::QstsResult)::Dict{String,Any}
  e = r.elements
  used = Set{String}()
  elements = Any[]
  add!(id, kind, name, main, plot, list, extra) = begin
    names, values = _qsts_view_quantities(r, id, list)
    entry = Dict{String,Any}("id" => String(id), "kind" => kind, "name" => String(name), "main" => main, "quantities" => names, "values" => values)
    plot === nothing || (entry["plot"] = plot)
    merge!(entry, extra)
    push!(elements, entry)
  end
  add!("", "run", "run", "total_p_loss_mw", nothing, QSTS_STEP_QUANTITIES, Dict{String,Any}())
  for c in eachindex(e.controller_ids)
    id = e.controller_ids[c]
    add!(id, "controller", e.controller_names[c], "position", _qsts_plot_name!(used, "controller", id), QSTS_CONTROLLER_QUANTITIES,
         Dict{String,Any}("controller_kind" => String(e.controller_kind[c]), "actuator" => e.controller_actuator[c], "target" => e.controller_target[c]))
  end
  for i in eachindex(e.bus_ids)
    e.bus_reported[i] || continue
    add!(e.bus_ids[i], "bus", e.bus_names[i], "vm_pu", _qsts_plot_name!(used, "bus", e.bus_ids[i]), QSTS_BUS_QUANTITIES, Dict{String,Any}())
  end
  for j in eachindex(e.branch_ids)
    add!(e.branch_ids[j], String(e.branch_kind[j]), e.branch_names[j], "loading_pct", _qsts_plot_name!(used, "branch", e.branch_ids[j]), QSTS_BRANCH_QUANTITIES, Dict{String,Any}())
  end
  q = r.qsts
  return Dict{String,Any}(
    "format" => _QSTS_VIEW_FORMAT,
    "name" => r.name,
    "mode" => String(r.mode),
    "labels" => String[s.label for s in r.steps],
    "bands" => Dict{String,Any}("vm_pu" => Any[q.vm_min_pu, q.vm_max_pu], "loading_pct" => Any[nothing, q.loading_max_pct]),
    "elements" => elements,
  )
end

# the band of `quantity` in a parsed view, as qstsChart(band = :auto) draws it
function _qsts_view_band(view::AbstractDict, quantity::AbstractString)::Union{Nothing,Tuple{Float64,Float64}}
  b = get(get(view, "bands", Dict{String,Any}()), String(quantity), nothing)
  b isa AbstractVector && length(b) == 2 || return nothing
  side(x) = x === nothing ? NaN : Float64(x)
  return (side(b[1]), side(b[2]))
end

"""
    _qsts_view_chart_svg(view, element_id, quantity) -> Union{String,Nothing}

The chart of one quantity of one element of a parsed `qsts_view.json`,
drawn by the renderer of [`qstsChart`](@ref) with the band of the view, so
the chart route gives the bytes of the SVG the service wrote for the main
quantity. `nothing` when the view has no such element or quantity.
"""
function _qsts_view_chart_svg(view::AbstractDict, element_id::AbstractString, quantity::AbstractString)::Union{Nothing,String}
  idx = findfirst(el -> get(el, "id", nothing) == element_id, get(view, "elements", Any[]))
  idx === nothing && return nothing
  el = view["elements"][idx]
  raw = get(get(el, "values", Dict{String,Any}()), String(quantity), nothing)
  raw isa AbstractVector || return nothing
  values = Float64[x === nothing ? NaN : Float64(x) for x in raw]
  labels = String[String(l) for l in view["labels"]]
  return _qsts_chart_svg(element_id, quantity, labels, values, _qsts_view_band(view, quantity), _QSTS_CHART_WIDTH, _QSTS_CHART_HEIGHT)
end

# The summary of a run for run.log, the metadata and the overview: solved
# and failed steps, switch counts and hunting per controller (from the
# metrics of a sequential run), steps with a voltage outside the band or an
# overload (counted as the step quantities n_vm_low, n_vm_high and
# n_overload of mode B count them, so both modes give the same number),
# and the loss energy over the solved steps (the rule of the metric
# loss_energy_mwh; nothing without a step length or timestamps).
function _qsts_run_summary(r::QstsResult)
  n = length(r.steps)
  solved = [s.status === :solved for s in r.steps]
  e = r.elements
  q = r.qsts
  metric(id, name) = begin
    k = findfirst(m -> m.element_id == id && m.metric == name, r.metrics)
    k === nothing ? nothing : r.metrics[k].value
  end
  controllers = [(id = e.controller_ids[c], name = e.controller_names[c], kind = String(e.controller_kind[c]),
                  switches = something(metric(e.controller_ids[c], "switch_count"), 0.0), hunting = something(metric(e.controller_ids[c], "n_hunting"), 0.0))
                 for c in eachindex(e.controller_ids)]
  qvm = findfirst(==("vm_pu"), QSTS_BUS_QUANTITIES)
  qld = findfirst(==("loading_pct"), QSTS_BRANCH_QUANTITIES)
  vm_out(v) = !isnan(v) && (v < q.vm_min_pu || v > q.vm_max_pu)
  overload(v) = !isnan(v) && v > q.loading_max_pct
  violation(k) = any(i -> e.bus_reported[i] && vm_out(r.bus_values[i, qvm, k]), eachindex(e.bus_ids)) || any(j -> overload(r.branch_values[j, qld, k]), eachindex(e.branch_ids))
  durations = _qsts_step_durations(r.step_minutes, String[s.timestamp for s in r.steps])
  ploss = findfirst(==("total_p_loss_mw"), QSTS_STEP_QUANTITIES)
  loss = durations === nothing ? nothing : sum((r.step_values[ploss, k] * durations[k] / 60.0 for k in 1:n if solved[k]); init = 0.0)
  return (n_steps = n, n_solved = count(solved), n_failed = n - count(solved), base_status = String(r.base.status), controllers = controllers,
          switch_total = sum((c.switches for c in controllers); init = 0.0), hunting_total = sum((c.hunting for c in controllers); init = 0.0),
          violation_steps = count(violation, 1:n), loss_energy_mwh = loss)
end

# the charts the result page opens with: the position of the first
# controller and the voltage of its regulated bus (with the band); without
# a controller the voltage of the first bus
function _qsts_default_charts(view::AbstractDict)::Vector{Vector{String}}
  els = get(view, "elements", Any[])
  ctl = findfirst(el -> get(el, "kind", "") == "controller", els)
  charts = Vector{String}[]
  if ctl !== nothing
    push!(charts, [String(els[ctl]["id"]), "position"])
    target = String(get(els[ctl], "target", ""))
    isempty(target) || any(el -> get(el, "id", nothing) == target && get(el, "kind", "") == "bus", els) && push!(charts, [target, "vm_pu"])
  else
    bus = findfirst(el -> get(el, "kind", "") == "bus", els)
    bus === nothing || push!(charts, [String(els[bus]["id"]), "vm_pu"])
  end
  return charts
end

_qsts_html_escape(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")
_qsts_fmt(x) = x === nothing ? "n/a" : (isinteger(x) ? string(Int(x)) : @sprintf("%.4g", x))

# the overview: summary, notes, metrics and switching counts, then the SVG
# of every element grouped by kind (embedded with <object>, so the
# tooltips of the points work when the run folder or its ZIP is opened)
function _qsts_overview_html(r::QstsResult, view::AbstractDict, summary, control_note::AbstractString, freq_line::AbstractString)::String
  io = IOBuffer()
  title = string("QSTS overview: ", r.name)
  println(io, "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n<title>", _qsts_html_escape(title), "</title>")
  println(io, "<style>\nbody {\n  font-family: sans-serif;\n  margin: 1.5em;\n}\ntable {\n  border-collapse: collapse;\n}\nth,\ntd {\n  border: 1px solid #ccc;\n  padding: 0.2em 0.5em;\n  text-align: left;\n}\nfigure {\n  display: inline-block;\n  margin: 0.5em;\n}\n</style>\n</head>\n<body>")
  println(io, "<h1>", _qsts_html_escape(title), "</h1>")
  println(io, "<table>")
  rows = [("Mode", String(r.mode)), ("Steps", string(summary.n_steps)), ("Solved", string(summary.n_solved)), ("Failed", string(summary.n_failed)),
          ("Base case", summary.base_status), ("Switches", _qsts_fmt(summary.switch_total)), ("Hunting patterns", _qsts_fmt(summary.hunting_total)),
          ("Steps with violations", string(summary.violation_steps)), ("Loss energy [MWh]", _qsts_fmt(summary.loss_energy_mwh)), ("Controllers", control_note),
          ("System frequency", freq_line)]
  for (k, v) in rows
    println(io, "<tr><th>", _qsts_html_escape(k), "</th><td>", _qsts_html_escape(v), "</td></tr>")
  end
  println(io, "</table>")
  if !isempty(summary.controllers)
    println(io, "<h2>Controllers</h2>\n<table>\n<tr><th>controller</th><th>kind</th><th>switches</th><th>hunting</th></tr>")
    for c in summary.controllers
      println(io, "<tr><td>", _qsts_html_escape(c.name), " (", _qsts_html_escape(c.id), ")</td><td>", c.kind, "</td><td>", _qsts_fmt(c.switches), "</td><td>", _qsts_fmt(c.hunting), "</td></tr>")
    end
    println(io, "</table>")
  end
  if !isempty(r.notes)
    println(io, "<h2>Notes</h2>\n<ul>")
    foreach(n -> println(io, "<li>", _qsts_html_escape(n), "</li>"), r.notes)
    println(io, "</ul>")
  end
  for (kind, heading) in (("controller", "Controllers: position"), ("bus", "Buses: voltage"), ("branch", "Branches: loading"))
    els = [el for el in view["elements"] if haskey(el, "plot") && (kind == "branch" ? !(el["kind"] in ("controller", "bus")) : el["kind"] == kind)]
    isempty(els) && continue
    println(io, "<h2>", heading, "</h2>")
    for el in els
      println(io, "<figure><object type=\"image/svg+xml\" data=\"", _qsts_html_escape(el["plot"]), "\"></object><figcaption>", _qsts_html_escape(el["name"]), "</figcaption></figure>")
    end
  end
  println(io, "</body>\n</html>")
  return String(take!(io))
end

"""
    _run_qsts_service(case_path, config_file, output_dir, run_id, mode; config_overrides = Dict(), progress = nothing, phase_callback = phase -> nothing) -> SparlectraApiResult

Service backend of the Web UI "QSTS" button. Reads the SCF time series at
`case_path` ([`read_scf_series`](@ref); the file is recognized by its
content, [`is_scf_series`](@ref)), resolves the run configuration like the
other run kinds, replaces its `control` section by the one of the
companion file `<series>.control.yaml` when that file exists, sets
`qsts.mode` to `mode` (`"sequential"` or `"snapshots"`) for this run only
and runs [`runQsts`](@ref). `progress(done, total)` is handed to it (one
call per finished step); the abort token of a Web UI job is checked before
every step. `phase_callback` receives the run phase by name:
`importing_case`, `running_steps`, `writing_artifacts`.

Writes the files of [`writeQstsResults`](@ref) (CSV format of the run
configuration), the element view `qsts_view.json`, one SVG per element
under `qsts_plots/` with its main quantity ([`qstsChart`](@ref): bus
`vm_pu`, branch `loading_pct`, controller `position`), `qsts_overview.html`
and a `run.log` that names the companion file. The metadata carries
`run_mode = "qsts"`, the counts of the summary and the default charts of
the result page.

Failure behavior: `qsts_no_series` (the case is no SCF time series),
`invalid_request` (unknown mode), `invalid_case_file` (the series fails a
reader check), `invalid_control_file` (the companion file does not load),
`qsts_invalid` (runQsts refuses the configuration or the series), plus the
shared configuration failures.
"""
function _run_qsts_service(case_path::AbstractString, config_file::AbstractString, output_dir::AbstractString, run_id::String, mode::AbstractString;
                           config_overrides::AbstractDict = Dict{String,Any}(), @nospecialize(progress = nothing), @nospecialize(phase_callback = phase -> nothing))::SparlectraApiResult
  mkpath(output_dir)
  logfile = joinpath(output_dir, "run.log")
  result_file = joinpath(output_dir, "result.json")
  base_metadata = Dict{String,Any}("run_mode" => "qsts", "qsts_mode" => String(mode))
  fail(reason, message) = _api_failure(reason, message; run_id = run_id, casefile = case_path, config_file = config_file, output_dir = String(output_dir), logfile = logfile, result_file = result_file, metadata = base_metadata)
  mode in _QSTS_SERVICE_MODES || return fail("invalid_request", "qsts_mode must be one of $(join(_QSTS_SERVICE_MODES, ", ")); got \"$(mode)\".")
  # the content decides (is_scf_series), not the name: a CGMES zip or a
  # case file is refused here, before anything is read
  is_scf_series(case_path) || return fail("qsts_no_series", "QSTS needs an SCF time series (a zip holding series.json, base.scf.json and delta.csv); $(basename(String(case_path))) is none.")
  # the precedence of every run kind (resolve_config); the form values of
  # the Web UI arrive as overrides
  config = try
    resolve_config(config_file, case_path, config_overrides).config
  catch err
    return fail(_config_resolve_reason(err), sprint(showerror, err))
  end
  phase_callback("importing_case")
  series = try
    read_scf_series(case_path)
  catch err
    # every reader check throws an ArgumentError naming the entry and line
    err isa ArgumentError || rethrow()
    return fail("invalid_case_file", sprint(showerror, err))
  end
  control_path = _qsts_control_file_path(case_path)
  # a snapshot run solves every step without controllers unless
  # qsts.snapshot_controls asks for them (the companion file is then read
  # but not used); a sequential run always runs them
  controls_used = mode == "sequential" || config.qsts.snapshot_controls
  control_note = if isfile(control_path)
    control = try
      load_sparlectra_config(control_path; reload = true).control
    catch err
      return fail("invalid_control_file", string(basename(control_path), ": ", sprint(showerror, err)))
    end
    config = _sparlectra_config_with(config; control = control)
    string("control section from ", basename(control_path), " (", length(control.controllers), " controller(s))")
  else
    string("no companion control file (", basename(control_path), " not found); the controllers of the case file only")
  end
  controls_used || (control_note = string(control_note, "; not used: a snapshot run without qsts.snapshot_controls runs without controllers"))
  # the run mode for this run only; the configuration file keeps qsts.mode
  config = _sparlectra_config_with(config; qsts = _struct_with(config.qsts; mode = Symbol(mode)))
  phase_callback("running_steps")
  r = try
    runQsts(series; config = config, progress = progress)
  catch err
    # runQsts names a refused configuration or series with an ArgumentError
    # (a snapshot start on an event series, a qsts key out of range); an
    # abort and every other error go up to the service (execution_error)
    err isa ArgumentError || rethrow()
    return fail("qsts_invalid", sprint(showerror, err))
  end

  phase_callback("writing_artifacts")
  csv_format = String(config.output.csv_format)
  paths = writeQstsResults(output_dir, r; format = csv_format)
  view = _qsts_view_dict(r)
  write(joinpath(output_dir, _QSTS_VIEW_FILE), scf_json_string(view))
  # one SVG per element with its main quantity, drawn by qstsChart
  n_plots = 0
  for el in view["elements"]
    haskey(el, "plot") || continue
    path = joinpath(output_dir, el["plot"])
    mkpath(dirname(path))
    write(path, qstsChart(r, el["id"], el["main"]; width = _QSTS_CHART_WIDTH, height = _QSTS_CHART_HEIGHT).svg)
    n_plots += 1
  end
  summary = _qsts_run_summary(r)
  # the frequency the base case is built with, from the case itself: a
  # series in pu or with c1 = 0 never reads it, whatever the key says
  freq_line = nominal_frequency_line(series.base, config)
  write(joinpath(output_dir, _QSTS_OVERVIEW_FILE), _qsts_overview_html(r, view, summary, control_note, freq_line))

  open(logfile, "a") do io
    println(io, "QSTS run (", mode, ") on ", basename(String(case_path)), ": ", r.name, ", ", summary.n_steps, " steps")
    println(io, "System frequency: ", freq_line)
    println(io, "control file: ", control_note)
    println(io, "base case: ", summary.base_status)
    println(io, "steps: ", summary.n_solved, " solved, ", summary.n_failed, " failed")
    if r.mode === :sequential
      println(io, "switches: ", _qsts_fmt(summary.switch_total), " (", join((string(c.name, ": ", _qsts_fmt(c.switches)) for c in summary.controllers), ", "), ")")
      println(io, "hunting patterns: ", _qsts_fmt(summary.hunting_total))
    else
      println(io, "switches: none counted (snapshot run, every step on its own)")
    end
    println(io, "steps with a voltage outside ", r.qsts.vm_min_pu, " to ", r.qsts.vm_max_pu, " pu or a loading above ", r.qsts.loading_max_pct, " %: ", summary.violation_steps)
    println(io, "loss energy: ", summary.loss_energy_mwh === nothing ? "n/a (no step length)" : string(_qsts_fmt(summary.loss_energy_mwh), " MWh"))
    for note in r.notes
      println(io, "note: ", note)
    end
    println(io, "Artifacts: ", join((basename(p) for p in paths), ", "), ", ", _QSTS_VIEW_FILE, ", ", _QSTS_OVERVIEW_FILE, ", ", n_plots, " SVG file(s) in ", _QSTS_PLOT_DIR, "/")
  end

  metadata = merge(base_metadata, Dict{String,Any}(
    "input_format_detected" => "scf_series",
    "qsts_series" => r.name,
    "qsts_control_file" => isfile(control_path) ? basename(control_path) : "",
    "network_nominal_frequency_hz" => _scf_case_frequency(series.base, config).f,
    "frequency_source" => nominal_frequency_source(series.base, config),
    "frequency_line" => freq_line,
    "qsts_controls_used" => controls_used,
    "qsts_steps" => summary.n_steps,
    "qsts_solved" => summary.n_solved,
    "qsts_failed" => summary.n_failed,
    "qsts_base_status" => summary.base_status,
    "qsts_switch_counts" => Dict{String,Any}(c.id => c.switches for c in summary.controllers),
    "qsts_hunting" => Dict{String,Any}(c.id => c.hunting for c in summary.controllers),
    "qsts_switch_total" => summary.switch_total,
    "qsts_hunting_total" => summary.hunting_total,
    "qsts_violation_steps" => summary.violation_steps,
    "qsts_loss_energy_mwh" => summary.loss_energy_mwh,
    "qsts_notes" => copy(r.notes),
    "qsts_plots" => n_plots,
    "qsts_default_charts" => _qsts_default_charts(view),
    # the bus and branch CSV files of a single power flow are not written
    # here (the qsts_*.csv files hold the series); the job runner logs the
    # detailed CSV export as skipped with this reason instead of failed
    "detailed_result_csv_status" => "skipped",
    "detailed_result_csv_skip_reason" => "qsts_run",
    "artifact_status" => "completed",
    "solver_status" => "completed",
    "service_status" => "completed",
    "run_status" => "completed",
  ))
  message = string("QSTS (", mode, ") completed - ", summary.n_solved, " of ", summary.n_steps, " steps solved",
                   summary.n_failed > 0 ? string(", ", summary.n_failed, " failed") : "",
                   r.mode === :sequential ? string(", ", _qsts_fmt(summary.switch_total), " switches") : "", ".")
  result = _api_result(
    run_id = run_id,
    status = :succeeded,
    success = true,
    solution_available = false,
    reason = nothing,
    message = message,
    casefile = String(case_path),
    config_file = String(config_file),
    output_dir = String(output_dir),
    logfile = logfile,
    result_file = result_file,
    metadata = metadata,
  )
  return _finalize_api_result(result)
end
