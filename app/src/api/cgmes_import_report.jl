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

# file: app/src/api/cgmes_import_report.jl
# purpose: the CGMES import report of a service run (cgmes_import_report.md):
#          Markdown writers for the success report, the failure report and
#          the SV comparison section, shared by the power-flow and the
#          short-circuit service

"""File name of the CGMES import report artifact, the same for every service."""
const CGMES_IMPORT_REPORT_NAME = "cgmes_import_report.md"

# Text from the delivery (names, mRIDs, importer messages) goes into the
# report escaped: the Markdown standard library of the Web UI viewer reads
# an underscore inside a word as emphasis (ACLineSegment_1 b_ch turns into
# italics), a pair of dollar signs as inline math, and a pipe would split a
# table cell. Line breaks become spaces so one message stays one list item
# or one table row. The standard library
# knows the pipe escape only inside a table (`cell = true`); in a paragraph
# it would print the backslash.
function _md_text(s; cell::Bool = false)::String
  t = replace(string(s), "\r\n" => " ", '\n' => ' ', '\r' => ' ')
  t = replace(t, r"([\\`*_\[\]#$])" => s"\\\1")
  return cell ? replace(t, '|' => "\\|") : t
end

# File names and paths as code spans; a backtick inside falls back to the
# escaped text, a code span cannot hold one in the standard library.
_md_code(s; cell::Bool = false)::String = occursin('`', string(s)) ? _md_text(s; cell) : "`" * (cell ? replace(string(s), '|' => "\\|") : string(s)) * "`"

function _md_key(io::IO, key::AbstractString, value)
  println(io, "- **", key, ":** ", value)
  return nothing
end

# files: entries with name, profiles, skipped, skip_reason, object_count
# (the CGMES store files and the summarizeCGMES files have the same fields)
function _md_cgmes_files_table(io::IO, files)
  println(io, "## Files\n")
  println(io, "| file | profiles | objects |")
  println(io, "|:-----|:---------|:--------|")
  for f in files
    println(io, "| ", _md_code(f.name; cell = true), " | ", _md_text(join(sort(String.(collect(f.profiles))), ", "); cell = true), " | ", f.skipped ? "SKIPPED: " * _md_text(f.skip_reason; cell = true) : string(f.object_count), " |")
  end
  println(io)
  return nothing
end

function _md_cgmes_class_table(io::IO, histogram)
  println(io, "## Class histogram\n")
  println(io, "| class | objects |")
  println(io, "|:------|--------:|")
  for (cls, n) in histogram
    println(io, "| ", _md_text(cls; cell = true), " | ", n, " |")
  end
  println(io)
  return nothing
end

# The coverage rows of shortCircuitCoverage as one table: a class without
# records gets one row, a class with records one row per attribute.
function _md_cgmes_sc_coverage(io::IO, sc)
  println(io, "## Short-circuit source data (read, not evaluated, issue #277)\n")
  println(io, "| class | records | attribute | filled |")
  println(io, "|:------|--------:|:----------|:-------|")
  for row in shortCircuitCoverage(sc)
    if row.records == 0
      # line/transformer/EI records exist only when at least one
      # short-circuit attribute is present (see printShortCircuitCoverage)
      println(io, "| ", _md_text(row.class; cell = true), " | 0 | none (no records with short-circuit attributes) | |")
      continue
    end
    for (i, a) in enumerate(row.attributes)
      filled = a.filled == a.total ? "all" : (a.filled == 0 ? "none" : string(a.filled, "/", a.total))
      println(io, "| ", i == 1 ? _md_text(row.class; cell = true) : "", " | ", i == 1 ? string(row.records) : "", " | ", _md_text(a.attribute; cell = true), " | ", filled, " |")
    end
  end
  println(io)
  return nothing
end

"""
    _write_cgmes_import_report(path, case_path, paths, cgmes_result, cgmes_cfg, start_decision)

Writes the CGMES import report of a successful import as Markdown: header
keys, files, class histogram, network built, DY profile (only when the
delivery carries one), short-circuit source coverage, importer messages.
The SV comparison is appended after the solve by
`_append_cgmes_sv_comparison`. Throws `SystemError`/`Base.IOError` when the
file cannot be written; the caller decides how to report that.
"""
function _write_cgmes_import_report(path::AbstractString, case_path::AbstractString, paths, cgmes_result, cgmes_cfg, start_decision::AbstractString)
  open(path, "w") do io
    println(io, "# CGMES import report\n")
    _md_key(io, "Case", _md_code(case_path))
    _md_key(io, "Paths", join(_md_code.(paths), ", "))
    _md_key(io, "Version", _md_text(cgmes_result.store.version))
    _md_key(io, "Tap control", cgmes_cfg.tap_control)
    _md_key(io, "Machine control", cgmes_cfg.machine_control)
    _md_key(io, "Base MVA", cgmes_cfg.base_mva)
    _md_key(io, "Read dynamics", cgmes_cfg.read_dynamics)
    _md_key(io, "Equivalent injections as source", cgmes_cfg.equivalent_injection_as_source)
    # the decision line starts with its own label ("CGMES start values: ...")
    println(io, "- ", _md_text(start_decision))
    println(io)
    _md_cgmes_files_table(io, cgmes_result.store.files)
    _md_cgmes_class_table(io, ((cls, length(mrids)) for (cls, mrids) in sort(collect(cgmes_result.store.byclass); by = p -> -length(p.second))))
    net = cgmes_result.net
    println(io, "## Network built\n")
    _md_key(io, "Buses", string(length(net.nodeVec), " (", length(net.isoNodes), " isolated)"))
    _md_key(io, "Branches", length(net.branchVec))
    _md_key(io, "Links", length(net.linkVec))
    _md_key(io, "Prosumers", length(net.prosumpsVec))
    _md_key(io, "Shunts", length(net.shuntVec))
    _md_key(io, "Slack bus", _md_code(cgmes_result.slack_bus))
    _md_key(io, "No-SV buses", string(length(cgmes_result.no_sv_buses), " (start from a neighbour with SV, else 1.0 pu / 0°)"))
    println(io)
    # DY profile: what read_dynamics read and mapped, or that it was left
    # out; no section when the delivery carries none
    dy_lines = dynamicsReportLines(cgmes_result)
    if !isempty(dy_lines)
      println(io, "## Dynamic data (DY profile)\n")
      foreach(l -> println(io, "- ", _md_text(l)), dy_lines)
      println(io)
    end
    # per-attribute completeness: an IEC 60909 evaluation needs to know which
    # inputs a delivery provides, not just the counts
    _md_cgmes_sc_coverage(io, cgmes_result.shortcircuit)
    # one list item per message: the viewer wraps list text, a fenced block
    # would scroll sideways on the long notices
    println(io, "## Importer messages (", length(cgmes_result.messages), ")\n")
    foreach(m -> println(io, "- ", _md_text(m)), cgmes_result.messages)
    isempty(cgmes_result.messages) || println(io)
    if !isempty(cgmes_result.skipped_equipment)
      println(io, "## Skipped equipment\n")
      println(io, length(cgmes_result.skipped_equipment), " mRID(s)\n")
    end
  end
  return path
end

"""
    _write_cgmes_import_failure_report(path; message, analysis = "", case_path = nothing, paths = String[], summary = nothing, summary_error = nothing)

Writes the CGMES import report of a failed import as Markdown: the error,
the import analysis of the typed `CGMESImportError` as a fenced block (a
preformatted text report with its own indentation), and, when given, what
`summarizeCGMES` could still read (files, class histogram). `summary_error`
is the text of a failed `summarizeCGMES` call and lands in the report in
place of the summary. Throws `SystemError`/`Base.IOError` when the file
cannot be written.
"""
function _write_cgmes_import_failure_report(path::AbstractString; message::AbstractString, analysis::AbstractString = "", case_path = nothing, paths = String[], summary = nothing, summary_error = nothing)
  open(path, "w") do io
    println(io, "# CGMES import report: IMPORT FAILED\n")
    case_path === nothing || _md_key(io, "Case", _md_code(case_path))
    isempty(paths) || _md_key(io, "Paths", join(_md_code.(paths), ", "))
    _md_key(io, "Error", _md_text(message))
    println(io)
    if !isempty(analysis)
      println(io, "## Import analysis\n")
      println(io, "```text")
      print(io, analysis)
      endswith(analysis, "\n") || println(io)
      println(io, "```\n")
    end
    if summary_error !== nothing
      println(io, "## Delivery summary\n")
      println(io, "not available: ", _md_text(summary_error), "\n")
    elseif summary !== nothing
      println(io, "## Delivery summary\n")
      _md_key(io, "Version", _md_text(summary.version))
      _md_key(io, "Objects", string(summary.object_count, " in ", length(summary.class_histogram), " classes"))
      _md_key(io, "Unresolved references", string(summary.unresolved_count, summary.boundary_missing_hint ? " (topology references, boundary set likely missing)" : ""))
      isempty(summary.unresolved_sample) || _md_key(io, "Sample", join(_md_code.(summary.unresolved_sample), ", "))
      println(io)
      _md_cgmes_files_table(io, summary.files)
      _md_cgmes_class_table(io, summary.class_histogram)
    end
  end
  return path
end

"""
    _append_cgmes_sv_comparison(path; sv_compare, sv_compare_error, status, start_values, converged, sv_bus_count)

Appends the SV comparison section to the CGMES import report after the
solve. `sv_compare === nothing` writes the unavailable line with
`sv_compare_error`. Throws `SystemError`/`Base.IOError` when the file
cannot be written.
"""
function _append_cgmes_sv_comparison(path::AbstractString; sv_compare, sv_compare_error, status::AbstractString, start_values::Symbol, converged::Bool, sv_bus_count::Integer)
  open(path, "a") do io
    println(io, "## SV comparison\n")
    _md_key(io, "Run status", _md_code(status))
    _md_key(io, "Start values", start_values)
    if sv_compare === nothing
      _md_key(io, "Comparison", "unavailable" * (sv_compare_error === nothing ? "" : ": " * _md_text(sv_compare_error)))
    else
      _md_key(io, "Voltages", string("n=", sv_compare.n, ", max|dvm| ", sv_compare.max_dvm, " pu, rms ", sv_compare.rms_dvm, " pu, max|dva| ", sv_compare.max_dva, "°, rms ", sv_compare.rms_dva, "°") |> _md_text)
      # A cut-out IGM keeps the CGM's global angle reference while the
      # local solve pins its own slack, the uniform offset is removed from
      # the dva judgement and reported here instead.
      isfinite(sv_compare.va_ref_offset_deg) && abs(sv_compare.va_ref_offset_deg) > 0.01 && _md_key(io, "Angle reference offset removed", string(round(sv_compare.va_ref_offset_deg; digits = 2), "° (median; raw dva column stays in `sv_compare.csv`)"))
      _md_key(io, "Flows", string("n=", sv_compare.flows.n, ", max|dp| ", sv_compare.flows.max_dp, " MW, max|dq| ", sv_compare.flows.max_dq, " MVAr") |> _md_text)
      _md_key(io, "Artifacts", "`sv_compare.csv`, `sv_compare_flows.csv`")
      converged || _md_key(io, "Note", "run did not converge, deltas describe the non-converged state, not the imported model")
      # A real delivery ships its own solved state. When a flat-start run
      # fails on a delivery that HAS usable SvVoltage data, the single most
      # effective remedy is starting from that state, say so instead of
      # leaving the user to find it (measured on a 6209-bus delivery: flat
      # start diverges after 80 iterations, SV start converges in 4).
      if !converged && start_values !== :sv && sv_bus_count > 0
        _md_key(io, "Hint", "this delivery carries an SvVoltage state, set `cgmes_import.start_values: sv` (Web UI: \"CGMES start values\" = \"Imported SV state\") to start the solve from it.")
      end
    end
    println(io)
  end
  return path
end
