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

# Author: Udo Schmitz (https://github.com/Welthulk)
# Purpose: prepare the title/body of the "JuliaRegistrator register" issue
#          for the "Request Julia Registration" workflow (.github/workflows/request_registration.yml)

using TOML

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CHANGELOG_FILE = joinpath(REPO_ROOT, "docs", "src", "changelog.md")
const PROJECT_TOML_FILE = joinpath(REPO_ROOT, "Project.toml")
const TITLE_FILE = joinpath(REPO_ROOT, "release_issue_title.txt")
const BODY_FILE = joinpath(REPO_ROOT, "release_issue_body.txt")

struct ChangelogEntry
  version::VersionNumber
  date::AbstractString
  body::AbstractString
end

# Splits docs/src/changelog.md into per-version entries, newest first, matching
# both "# Version x.y.z — date" and "## Version x.y.z – date" heading styles.
function parseChangelog(path::AbstractString)::Vector{ChangelogEntry}
  text = read(path, String)
  headingRe = r"(?m)^#{1,2}[ \t]+Version[ \t]+([0-9]+\.[0-9]+\.[0-9]+)[ \t]*[—–-][ \t]*(\S+)[ \t]*$"
  matches = collect(eachmatch(headingRe, text))
  isempty(matches) && error("No 'Version x.y.z' heading found in $path")

  entries = ChangelogEntry[]
  for (i, m) in enumerate(matches)
    version = VersionNumber(m.captures[1])
    date = m.captures[2]
    bodyStart = m.offset + ncodeunits(m.match)
    bodyStop = i < length(matches) ? matches[i+1].offset : ncodeunits(text)
    body = strip(text[bodyStart:bodyStop])
    push!(entries, ChangelogEntry(version, date, body))
  end
  return entries
end

# Git tags follow the "vX.Y.Z" convention (created by TagBot after a
# registration is merged), so an entry without a matching tag is unreleased.
function registeredTags()::Set{String}
  return Set(readlines(Cmd(`git -C $REPO_ROOT tag`)))
end

function main()
  entries = parseChangelog(CHANGELOG_FILE)
  tags = registeredTags()

  idx = findfirst(e -> "v$(e.version)" ∉ tags, entries)
  isnothing(idx) && error("Every changelog version already has a matching git tag; nothing to register.")
  entry = entries[idx]

  projectVersion = VersionNumber(TOML.parsefile(PROJECT_TOML_FILE)["version"])
  if entry.version != projectVersion
    error(
      "Unreleased changelog version ($(entry.version)) does not match Project.toml ($projectVersion). " *
      "Align docs/src/changelog.md and Project.toml before requesting registration.",
    )
  end

  title = "JuliaRegistrator register v$(entry.version)"
  body = """
  @JuliaRegistrator register

  Release notes:

  ## Version $(entry.version)
  Released $(entry.date)

  $(entry.body)
  """

  write(TITLE_FILE, title)
  write(BODY_FILE, body)

  println("Prepared registration issue for v$(entry.version): \"$title\"")
end

main()
