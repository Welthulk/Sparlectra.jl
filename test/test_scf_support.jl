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

# file: test/test_scf_support.jl
# purpose: one place for comparing a freshly written SCF file against a
#          tracked fixture. A shipped SCF file records the release that
#          wrote it (`created_by: Sparlectra <version>`), so a plain byte
#          equality against a tracked file breaks on EVERY version bump
#          while nothing is wrong with the round trip. That red run trains
#          the reflex to regenerate the fixture, and the next time the diff
#          may be more than the stamp line, with the regeneration hiding
#          it (0.10.0 -> 0.11.0 was the first occurrence). Split the claim
#          instead: the fresh file names the CURRENT version, and
#          everything apart from that stamp is byte for byte. Both callers
#          use these helpers so they cannot drift apart.

"The provenance stamp an SCF file carries, whatever release wrote it."
const SCF_PROVENANCE_STAMP = r"\"created_by\": \"Sparlectra [^\"]*\""

"An SCF document with its provenance stamp neutralized, for byte comparison."
scf_without_stamp(text::AbstractString) = replace(String(text), SCF_PROVENANCE_STAMP => "\"created_by\": \"<provenance stamp>\"")

"""
    scf_matches_fixture(fresh_path, fixture_path) -> NamedTuple

The two halves a fixture comparison really claims, separately, so a failure
names which one broke:

- `stamp_is_current`: the freshly written file records the running version.
- `rest_identical`: everything except that stamp is byte for byte identical.

Use this for every comparison between a file written during the test run and
a tracked fixture. Two files written by the same run carry the same stamp and
can be compared directly.
"""
function scf_matches_fixture(fresh_path::AbstractString, fixture_path::AbstractString)
  fresh = read(fresh_path, String)
  return (
    stamp_is_current = occursin("\"created_by\": \"Sparlectra $(Sparlectra.version())\"", fresh),
    rest_identical = scf_without_stamp(fresh) == scf_without_stamp(read(fixture_path, String)),
  )
end

"""
    scf_with_foreign_stamp(path, directory) -> String

A copy of `path` in `directory` whose provenance stamp names a release that
does not exist. Reading it must produce the same network as the original: a
shipped case that behaves differently because of the version that wrote it
would tie behavior to provenance.
"""
function scf_with_foreign_stamp(path::AbstractString, directory::AbstractString)
  target = joinpath(directory, basename(path))
  write(target, replace(read(path, String), SCF_PROVENANCE_STAMP => "\"created_by\": \"Sparlectra 0.0.0-nonexistent\""))
  return target
end

# Field-by-field round-trip comparison, shared with the scenario tests
# a per-unit value may land one ulp away when no exact preimage exists
_scf_rt_same(a::Real, b::Real) = a == b || (isfinite(a) && isfinite(b) && abs(a - b) <= 1e-12 * max(abs(a), abs(b)))
_scf_rt_same(a, b) = isequal(a, b)

# Comparison contract, so that nothing is skipped silently: scalars (numbers,
# bools, symbols, strings, ENUMS such as the node type) compare directly,
# vectors of reals element-wise with the one-ulp tolerance, dictionaries and
# sets exactly, and a struct-valued field (the component with its name and
# type, a controller, a flow record) is entered ONE level and its scalar
# fields compared. Only vectors of structs stay length-only here, because
# their elements are compared as components in their own right.
_scf_rt_scalar(x) = x isa Number || x isa Bool || x isa Symbol || x isa AbstractString || x isa Enum || x === nothing

function _scf_rt_value_diff(f, xa, xb, out::Vector{String}; depth::Int = 0)
  if xa isa AbstractVector && xb isa AbstractVector && eltype(xa) <: Real && eltype(xb) <: Real
    if length(xa) != length(xb)
      push!(out, string(f, ": length ", length(xa), " vs ", length(xb)))
    else
      n = count(i -> !_scf_rt_same(xa[i], xb[i]), eachindex(xa))
      n == 0 || push!(out, string(f, ": ", n, " element(s) differ"))
    end
    return
  end
  if xa isa AbstractDict || xa isa AbstractSet
    isequal(xa, xb) || push!(out, string(f, ": containers differ"))
    return
  end
  if xa isa AbstractVector || xb isa AbstractVector
    length(xa) == length(xb) || push!(out, string(f, ": length ", length(xa), " vs ", length(xb)))
    return
  end
  if _scf_rt_scalar(xa) || _scf_rt_scalar(xb)
    _scf_rt_same(xa, xb) || push!(out, string(f, ": ", xa, " vs ", xb))
    return
  end
  # struct-valued: one level of scalar fields, deeper nesting stays out
  if depth == 0 && typeof(xa) === typeof(xb)
    for g in fieldnames(typeof(xa))
      _scf_rt_value_diff(string(f, ".", g), getfield(xa, g), getfield(xb, g), out; depth = 1)
    end
  elseif typeof(xa) !== typeof(xb)
    push!(out, string(f, ": ", typeof(xa), " vs ", typeof(xb)))
  end
  return
end

# Reflection round-trip comparison. Hand-picked assertions found the lost
# operating point only after it had already reached a power flow, and then the
# same thing happened again with the per-bus voltage limits and the network's
# Q-limit switching parameters. So the test compares EVERY field of every
# component and of the network itself, and a field may differ only if it is
# named below with a reason. Adding a field to the model without carrying it
# therefore fails here, which a positive list can never do.
const _SCF_RT_ALLOWED = Dict{Symbol,Dict{Symbol,String}}(
  :net => Dict{Symbol,String}(
    :name => "the case name is an export argument, not a property of the model",
    :matpower_branch_metadata => "MATPOWER reporting metadata (rateB/rateC and the source row), not part of the model",
    :matpowerDclineMetadata => "MATPOWER DC-line bookkeeping; the injections themselves are carried as prosumers",
    :for001Contingencies => "DTF/FOR001 outage list, carried by the study block instead",
    :totalLosses => "a run result, not an input",
    :totalBusPower => "a run result, not an input",
    :control_result => "a run result",
    :qLimitLog => "a run log",
    :qLimitEvents => "a run log",
    :qLimitInitialPVRows => "built during the run",
    :_locked => "construction state",
    :_rectangular_pf_status => "a run result",
    :_dc_pf_status => "a run result",
    :_import_config => "the configuration the net was imported with; session state, like the two status fields above, and not part of the model",
    :busOriginalNameDict => "source-format naming aid; the reference names are compared directly",
    :busOrigIdxDict => "source-format index aid; the bus order is compared directly",
    :cgmes_ids => "CGMES mRIDs travel in `external_id` and are compared through the names",
  ),
  :branch => Dict{Symbol,String}(
    :tap_min => "the regulation band snaps to the step grid; documented tolerance of less than one step",
    :tap_max => "the regulation band snaps to the step grid; documented tolerance of less than one step",
  ),
  :node => Dict{Symbol,String}(),
  :prosumer => Dict{Symbol,String}(),
  :shunt => Dict{Symbol,String}(),
)

function _scf_rt_field_diffs(a, b, kind::Symbol)
  allowed = _SCF_RT_ALLOWED[kind]
  out = String[]
  for f in fieldnames(typeof(a))
    haskey(allowed, f) && continue
    _scf_rt_value_diff(string(f), getfield(a, f), getfield(b, f), out)
  end
  return out
end

"""
Compare two networks field by field and return one message per differing
field, naming the component and an example. Empty means the round trip lost
nothing that is not explicitly allowed to differ.
"""
function scf_roundtrip_field_diffs(a, b)
  msgs = String[]
  append!(msgs, string("net.", m) for m in _scf_rt_field_diffs(a, b, :net))
  for (kind, va, vb) in ((:node, a.nodeVec, b.nodeVec), (:branch, a.branchVec, b.branchVec),
                         (:prosumer, a.prosumpsVec, b.prosumpsVec), (:shunt, a.shuntVec, b.shuntVec))
    if length(va) != length(vb)
      push!(msgs, string(kind, ": ", length(va), " vs ", length(vb), " elements"))
      continue
    end
    seen = Dict{String,Int}()
    for i in eachindex(va)
      for m in _scf_rt_field_diffs(va[i], vb[i], kind)
        key = string(kind, ".", first(split(m, ":")))
        haskey(seen, key) || (seen[key] = 0; push!(msgs, string(key, " (first at ", i, "): ", m)))
        seen[key] += 1
      end
    end
  end
  return msgs
end
