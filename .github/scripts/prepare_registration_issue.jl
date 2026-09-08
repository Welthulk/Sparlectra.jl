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
using Downloads

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const REGISTRY_VERSIONS_URL = "https://raw.githubusercontent.com/JuliaRegistries/General/master/S/Sparlectra/Versions.toml"
const CHANGELOG_FILE = joinpath(REPO_ROOT, "docs", "src", "changelog.md")
const PROJECT_TOML_FILE = joinpath(REPO_ROOT, "Project.toml")
const TITLE_FILE = joinpath(REPO_ROOT, "release_issue_title.txt")
const BODY_FILE = joinpath(REPO_ROOT, "release_issue_body.txt")
const BREAKING_FILE = joinpath(REPO_ROOT, "release_issue_breaking.txt")

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
    # up to the character BEFORE the next heading: including its offset drags
    # that heading's leading '#' into the notes, and those notes are what
    # JuliaRegistrator copies into the registry pull request (seen while
    # preparing 0.11.0, as a stray empty heading above "Breaking changes")
    bodyStop = i < length(matches) ? prevind(text, matches[i+1].offset) : ncodeunits(text)
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

# Latest version registered in the General registry, or nothing if the
# lookup fails (network hiccup must not block a registration request).
function latestRegisteredVersion()::Union{VersionNumber,Nothing}
  try
    buf = IOBuffer()
    Downloads.download(REGISTRY_VERSIONS_URL, buf)
    versions = keys(TOML.parse(String(take!(buf))))
    return maximum(VersionNumber.(collect(versions)))
  catch err
    println("::warning::Could not query the General registry ($err); skipping sequential-version check.")
    return nothing
  end
end

# RegistryCI's sequential version number guideline: the next version must be
# exactly latest+1 in patch, minor, or major position.
function isSequentialSuccessor(latest::VersionNumber, next::VersionNumber)::Bool
  return next == VersionNumber(latest.major, latest.minor, latest.patch + 1) || next == VersionNumber(latest.major, latest.minor + 1, 0) || next == VersionNumber(latest.major + 1, 0, 0)
end

# Julia's pre-1.0 convention, which RegistryCI shares: below 1.0 the MINOR
# bump is the breaking one (0.9.19 -> 0.10.0), from 1.0 on it is the major
# bump. AutoMerge labels such a registration BREAKING and refuses to merge it
# unless the release notes mention "breaking" or "changelog":
# https://juliaregistries.github.io/RegistryCI.jl/stable/guidelines/
# 0.10.0 was blocked by exactly this on 2026-09-07 and had to be re-triggered
# by hand.
function isBreakingBump(previous::VersionNumber, next::VersionNumber)::Bool
  next.major > previous.major && return true
  return previous.major == 0 && next.major == 0 && next.minor > previous.minor
end

mentionsBreakingOrChangelog(text::AbstractString)::Bool = occursin(r"breaking|changelog"i, text)

# Offline stand-in for the registry lookup: TagBot creates a vX.Y.Z tag after
# every merged registration, so the highest tag tracks the registry closely
# enough to decide whether this bump is breaking.
function highestTaggedVersion(tags)::Union{VersionNumber,Nothing}
  versions = VersionNumber[]
  for tag in tags
    m = match(r"^v([0-9]+\.[0-9]+\.[0-9]+)$", tag)
    isnothing(m) || push!(versions, VersionNumber(m.captures[1]))
  end
  return isempty(versions) ? nothing : maximum(versions)
end

# The notes AutoMerge sees. A breaking registration whose changelog section
# happens to use neither word gets an explicit section appended; one that
# already says "breaking" or points at the changelog is left alone, so the
# maintainer's own wording always wins.
function releaseNotesFor(entry, previous::Union{VersionNumber,Nothing})::AbstractString
  notes = entry.body
  isnothing(previous) && return notes
  isBreakingBump(previous, entry.version) || return notes
  mentionsBreakingOrChangelog(notes) && return notes
  return notes * """


  ## Breaking changes

  Version $(entry.version) follows $(previous), and below 1.0 the minor bump is the breaking one, so this release is breaking. What changed is listed above and in the changelog (`docs/src/changelog.md`).
  """
end

function checkSequentialVersion(next::VersionNumber, latest::Union{VersionNumber,Nothing})
  isnothing(latest) && return
  if next <= latest
    error("Version $next is not newer than the latest registered version $latest.")
  elseif !isSequentialSuccessor(latest, next)
    println("::warning::Version $next skips ahead of the latest registered version $latest. " * "AutoMerge will block the registry PR - comment `[merge approved]` on the PR in JuliaRegistries/General to proceed.")
  else
    println("Sequential-version check OK: $latest -> $next.")
  end
end

function main()
  entries = parseChangelog(CHANGELOG_FILE)
  tags = registeredTags()

  idx = findfirst(e -> "v$(e.version)" ∉ tags, entries)
  isnothing(idx) && error("Every changelog version already has a matching git tag; nothing to register.")
  entry = entries[idx]

  projectVersion = VersionNumber(TOML.parsefile(PROJECT_TOML_FILE)["version"])
  if entry.version != projectVersion
    error("Unreleased changelog version ($(entry.version)) does not match Project.toml ($projectVersion). " * "Align docs/src/changelog.md and Project.toml before requesting registration.")
  end

  # One registry lookup per run, used for both the sequential-version check
  # and the breaking decision; the local tags stand in when it fails.
  previous = latestRegisteredVersion()
  isnothing(previous) && (previous = highestTaggedVersion(tags))
  checkSequentialVersion(entry.version, previous)

  breaking = !isnothing(previous) && isBreakingBump(previous, entry.version)
  if isnothing(previous)
    println("::warning::Neither the registry nor a vX.Y.Z tag gave a previous version; the release notes carry no breaking notice. If this bump IS breaking, AutoMerge will block the registry PR until you re-trigger with notes that mention \"breaking\".")
  elseif breaking
    println("::notice::Breaking bump $previous -> $(entry.version): the release notes state it, which is what AutoMerge requires.")
  end

  title = "JuliaRegistrator register v$(entry.version)"
  body = """
  @JuliaRegistrator register

  Release notes:

  ## Version $(entry.version)
  Released $(entry.date)

  $(releaseNotesFor(entry, previous))
  """

  write(TITLE_FILE, title)
  write(BODY_FILE, body)
  write(BREAKING_FILE, breaking ? "true" : "false")

  println("Prepared registration issue for v$(entry.version): \"$title\"")
end

# Base.invokelatest: Julia 1.12 warns when a script entry point is called from
# a world older than its definition. The guard also lets a test include this
# file to exercise the helpers without creating an issue.
if abspath(PROGRAM_FILE) == @__FILE__
  Base.invokelatest(main)
end