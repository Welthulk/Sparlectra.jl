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

# file: src/cgmes/cgmes_import_analysis.jl
# purpose: detailed why-does-this-delivery-not-import analysis — supplied
#          models, declared md:Model.DependentOn prerequisites vs. what was
#          actually handed in, and an unresolved-reference histogram with a
#          plain-language verdict (typically: the boundary set is missing).

"""
    CGMESImportError(message, analysis)

Import abort that carries the full [`importFailureAnalysis`](@ref) report.
`showerror` prints only the one-line message; consumers that own a log or a
UI (the powerflow service, the Web UI) read the `analysis` field and write
the report where their users actually look (`cgmes.log`).
"""
struct CGMESImportError <: Exception
  message::String
  analysis::String
end

Base.showerror(io::IO, e::CGMESImportError) = print(io, e.message)

"""
    importFailureAnalysis(store::CGMESStore)::String

Build a multi-line report explaining why a CGMES delivery cannot be
imported (or what is incomplete about it):

1. the supplied model files with profile, CGMES version, and model id;
2. every `md:Model.DependentOn` prerequisite declared by the file headers,
   matched against the supplied model ids — missing prerequisites are the
   authoritative statement of what the delivery expects but did not get;
3. a histogram of unresolved object references grouped by class and
   property;
4. a verdict in plain language. A delivery whose `TopologicalNode.BaseVoltage`
   references stay unresolved depends on an external base-voltage catalog —
   in real ENTSO-E deliveries that catalog lives in the boundary set
   (`EQ_BD`/`TP_BD`), so importing without the matching boundary cannot work.

The report is purely diagnostic — building it never throws.
"""
function importFailureAnalysis(store::CGMESStore)::String
  io = IOBuffer()
  println(io, "CGMES import analysis")
  println(io, "=====================")

  supplied_ids = Set{String}()
  println(io, "Supplied models:")
  for info in store.files
    tag = isempty(info.profiles) ? "?" : join(sort(String.(collect(info.profiles))), "+")
    state = info.skipped ? " [skipped: $(info.skip_reason)]" : ""
    id = isempty(info.model_id) ? "(no model id)" : info.model_id
    println(io, "  - $(info.name): $(tag) v$(isempty(info.version) ? "?" : info.version), $(info.object_count) objects, model $(id)$(state)")
    isempty(info.model_id) || push!(supplied_ids, info.model_id)
  end

  missing_deps = Dict{String,Vector{String}}()  # dependency id => declaring files
  for info in store.files
    # Skipped files are not imported, so their declared prerequisites are
    # irrelevant; their model ids still count as supplied (see above).
    info.skipped && continue
    for dep in info.dependent_on
      dep in supplied_ids && continue
      push!(get!(Vector{String}, missing_deps, dep), info.name)
    end
  end
  if isempty(missing_deps)
    println(io, "Declared dependencies: all satisfied by the supplied files.")
  else
    println(io, "Declared dependencies MISSING from the input (md:Model.DependentOn):")
    for dep in sort!(collect(keys(missing_deps)))
      println(io, "  - model $(dep) — required by ", join(sort(missing_deps[dep]), ", "))
    end
  end

  unresolved = unresolvedReferences(store)
  if isempty(unresolved)
    println(io, "Unresolved references: none.")
  else
    hist = Dict{Tuple{Symbol,Symbol},Int}()
    for u in unresolved
      hist[(u.class, u.key)] = get(hist, (u.class, u.key), 0) + 1
    end
    rows = sort!(collect(hist); by = kv -> -kv[2])
    println(io, "Unresolved references: ", length(unresolved), " total")
    for (i, ((cls, key), n)) in enumerate(rows)
      i > 10 && (println(io, "  … and ", length(rows) - 10, " further reference kinds"); break)
      println(io, "  - $(cls).$(key): $(n)")
    end
  end

  # Verdict: name the most likely root cause in plain language.
  base_voltage_gap = any(u -> u.key == :BaseVoltage, unresolved)
  topology_gap = any(u -> u.class == :Terminal && u.key == :TopologicalNode, unresolved)
  if !isempty(missing_deps)
    println(io, "Verdict: the delivery declares $(length(missing_deps)) prerequisite model(s) that are not part of the input.")
    println(io, "         In real ENTSO-E deliveries these are typically the boundary files (EQ_BD/TP_BD) of the")
    println(io, "         matching date — the boundary carries the border X-nodes and the shared BaseVoltage catalog.")
    println(io, "         Obtain the matching boundary set and pass it as an additional path.")
  elseif base_voltage_gap
    println(io, "Verdict: BaseVoltage references stay unresolved — the base-voltage catalog is not part of the input.")
    println(io, "         It usually ships in the boundary EQ file (EQ_BD); import the delivery together with its boundary set.")
  elseif topology_gap
    println(io, "Verdict: terminals reference topological nodes that are not part of the input — a TP file is missing or truncated.")
  elseif !isempty(unresolved)
    println(io, "Verdict: cross-references stay unresolved; check that every profile file of the delivery (EQ, TP, SSH, SV) was supplied.")
  else
    println(io, "Verdict: no structural gaps detected at store level.")
  end
  return String(take!(io))
end

"""
    importabilityStats(store::CGMESStore) -> NamedTuple

The abort conditions the importer itself enforces, as data:
`missing_dependencies` (declared `md:Model.DependentOn` ids absent from the
input), `unresolved_count`, `base_voltage_gap` (nodes without a resolvable
voltage level), `topology_gap` (terminals pointing at absent topological
nodes), and the combined `importable` verdict. Shared by the service-mode
run and the Web UI upload check so both judge a delivery identically.
"""
function importabilityStats(store::CGMESStore)
  supplied = Set{String}(f.model_id for f in store.files if !isempty(f.model_id))
  missing_dependencies = sort!(unique(String[dep for f in store.files if !f.skipped for dep in f.dependent_on if !(dep in supplied)]))
  unresolved = unresolvedReferences(store)
  base_voltage_gap = any(u -> u.key == :BaseVoltage, unresolved)
  topology_gap = any(u -> u.class == :Terminal && u.key == :TopologicalNode, unresolved)
  return (; missing_dependencies, unresolved_count = length(unresolved), base_voltage_gap, topology_gap, importable = isempty(missing_dependencies) && !base_voltage_gap && !topology_gap)
end

"""
    analyzeCGMES(; path)::String

Load a CGMES delivery (folder, ZIP, XML file, or a vector of those — same
forms as `importCGMES`) without mapping it to a network, print the
[`importFailureAnalysis`](@ref) report, and return it. Use this to find out
*why* a delivery does not import — most commonly which declared boundary
dependency is missing — without wading through importer errors.
"""
function analyzeCGMES(; path)::String
  store = loadCGMES(path)
  report = importFailureAnalysis(store)
  print(report)
  return report
end
