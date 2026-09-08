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
