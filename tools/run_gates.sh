#!/bin/sh
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

# file: tools/run_gates.sh
# purpose: run the verification gates (a test profile or docs) with the two
#          guards that make the mixed-state incident impossible
#          instead of unlikely: (1) the run refuses to start while src/ or
#          the project files carry uncommitted changes, because a Julia
#          process precompiling WHILE sources keep changing tests a state
#          that exists nowhere (that incident cost one segfault and one
#          25-minute mixed-state run); (2) a lock directory serializes gate
#          runs, so no second Julia first-start precompiles the same
#          package concurrently. --allow-dirty skips guard 1 for a
#          deliberate local check; the lock is never skipped.
#          --threads=N starts the test process with N threads (default: one;
#          the threads profile always runs with four), so a profile also runs
#          every parallel-vs-serial comparison threaded. --strands=2 runs
#          extended as two processes at once (profiles extended_a and
#          extended_b), about a third less wall time. `set` is the full gate
#          set of one commit under one lock: three strands at once (set_a;
#          extended_b; set_c, then the docs build and the workshops), then
#          the threads profile alone. After every test profile (and every
#          part or strand) the gate prints "skipped/broken: <n>", the Broken
#          count of its Test Summary (a @test_skip counts there, for example
#          the DTF tests when data/DTF/*.DAT is missing in a worktree), and
#          at the end one line with the profiles that had any. --release
#          makes such a count fail the gate: exit 4 after all profiles ran,
#          when nothing else failed (a failing profile keeps its own code).
#          Plain sh.
# usage:   sh tools/run_gates.sh fast|pf|se|config|webui|extd|cgmes|formats|legacy|install|workshops|extended|all|threads|docs|set [--threads=N] [--strands=2] [--release] [--allow-dirty]
#          sh tools/run_gates.sh map [--since <ref>] [<path>...]
#          (SPARLECTRA_LARGE_CASES_DIR and other env pass through)
#          `map` runs nothing: it prints the profiles a change to the given
#          paths (or the files changed since <ref> and in the working tree)
#          has to run, in run order. `legacy` is run before a release and on
#          a change to a protected file; CI does not run it.

# --- skipped and broken tests ------------------------------------------------
# A @test_skip counts as Broken in the Julia Test Summary. A gate set in a
# worktree without the untracked DTF files printed green with four of them
# once, so every test profile reports its count on a line of its own.
# broken_count prints the Broken count of the last profile summary in the log
# $1: the field under "Broken" in the "Test Summary:" header directly above
# the "Sparlectra.jl ... |" line. Julia shows a column only when its total is
# not zero, so a summary without that column is 0. Nothing is printed when the
# log has no profile summary (the process crashed before it).
broken_count() {
  awk '
    prev ~ /^Test Summary:/ && /^Sparlectra\.jl .* [|] / {
      nh = split(substr(prev, index(prev, "|") + 1), head, " ")
      split(substr($0, index($0, "|") + 1), vals, " ")
      count = 0
      for (i = 1; i <= nh; i++) if (head[i] == "Broken") count = vals[i]
      found = 1
    }
    { prev = $0 }
    END { if (found) print count }
  ' "$1"
}
# report_skips <profile> <log>: the line "skipped/broken: <n>" of one profile;
# a profile with a count above zero is added to skip_profiles for the closing
# line and for --release
skip_profiles=
report_skips() {
  rs_n=$(broken_count "$2")
  if [ -z "$rs_n" ]
  then
    echo "run_gates: $1 skipped/broken: unknown (no test summary in its output)"
    return
  fi
  echo "run_gates: $1 skipped/broken: $rs_n"
  [ "$rs_n" -gt 0 ] && skip_profiles="${skip_profiles:+$skip_profiles }$1($rs_n)"
  return 0
}
# final_status <status> <release> <profiles with skips>: the exit status of
# the gate. A failing run keeps its own code; with --release a green run with
# skipped or broken tests ends with 4
final_status() {
  if [ "$1" -eq 0 ] && [ "$2" = "--release" ] && [ -n "$3" ]
  then
    echo 4
  else
    echo "$1"
  fi
}

gate=$1
[ $# -gt 0 ] && shift

# --- path to profile map -----------------------------------------------------
# Which profiles a change has to run, by the changed paths. `fast` goes with
# every code change; docs-only changes need the docs build alone. The map is
# written out in docs/src/tests.md; change both together.
map_path() {
  case "$1" in
    src/powerflow_rectangular/rectangular_network_solver.jl | src/acpflow/entrypoint.jl | app/src/api/run_api.jl | src/controller/control_framework.jl | src/config/configuration.jl | src/Sparlectra.jl)
      echo "fast pf legacy" ;;
    src/stateestimation/* | src/numerics/*) echo "fast se" ;;
    src/adapters/cgmes/*) echo "fast cgmes se" ;;
    src/adapters/scf/* | src/adapters/pgm/* | src/adapters/matpower/*) echo "fast formats extd" ;;
    src/adapters/dtf/* | src/adapters/powsybl/* | src/adapters/adapters.jl | src/import/*) echo "fast formats" ;;
    src/config/*) echo "fast config" ;;
    app/src/api/*) echo "fast extd" ;;
    app/src/webui/*) echo "fast webui" ;;
    app/*) echo "fast extd webui" ;;
    src/*) echo "fast pf" ;;
    ext/*) echo "fast pf" ;;
    test/fast/* | test/test_docstrings.jl) echo "fast" ;;
    test/runtests.jl | test/test_runner_helpers.jl) echo "fast config" ;;
    test/test_state_estimation.jl | test/test_observability.jl | test/test_legacy_topology_validation.jl | test/test_measurement_table.jl) echo "se" ;;
    test/test_legacy_* | test/extended/test_legacy_*) echo "legacy" ;;
    test/test_configuration_* | test/test_repository_hygiene.jl) echo "config" ;;
    test/test_cgmes_*) echo "cgmes" ;;
    test/test_powsybl_importer.jl | test/test_scf* | test/test_matpower_* | test/test_synthetic_grids.jl | test/extended/*) echo "formats" ;;
    test/test_api*.jl | test/test_demo_cases.jl) echo "extd" ;;
    test/test_webui*.jl) echo "webui" ;;
    test/test_install.jl | test/test_example_suites.jl) echo "install" ;;
    test/test_workshops.jl) echo "workshops" ;;
    test/fixtures/*) echo "cgmes formats pf" ;;
    test/testgrid.jl) echo "pf legacy" ;;
    test/*) echo "pf" ;;
    docs/lit/*) echo "workshops docs" ;;
    tools/* | Project.toml | app/Project.toml | .github/*) echo "install" ;;
    docs/* | README.md | *.md) echo "docs" ;;
    data/*) echo "pf extd formats cgmes" ;;
    *) echo "" ;;
  esac
}
if [ "$gate" = "map" ]
then
  paths=
  if [ "$1" = "--since" ]
  then
    root=$(git rev-parse --show-toplevel) || exit 2
    paths=$( (git -C "$root" diff --name-only "$2" -- ; git -C "$root" status --porcelain | sed 's/^...//') | sort -u)
  else
    for a in "$@"
    do
      paths="$paths
$a"
    done
  fi
  wanted=
  for pth in $paths
  do
    wanted="$wanted $(map_path "$pth")"
    # a legacy group includes the file it was cut from, so a change to that
    # file reaches the legacy profile as well
    case "$pth" in
      test/test_legacy_* | test/extended/test_legacy_*) ;;
      test/*.jl | test/extended/*.jl)
        stem=${pth##*/}
        stem=${stem#test_}
        [ -f "${pth%/*}/test_legacy_$stem" ] && wanted="$wanted legacy"
        ;;
    esac
  done
  out=
  for p in fast pf se config cgmes formats extd webui install workshops legacy docs
  do
    case " $wanted " in
      *" $p "*) out="$out $p" ;;
    esac
  done
  echo "${out# }"
  exit 0
fi
allow_dirty=
julia_threads=
strands=1
release=
for arg in "$@"
do
  case "$arg" in
    --allow-dirty) allow_dirty=--allow-dirty ;;
    --release) release=--release ;;
    --threads=[1-9]*) julia_threads=$arg ;;
    --strands=1) strands=1 ;;
    --strands=2) strands=2 ;;
    *)
      echo "run_gates: unknown option $arg" >&2
      echo "usage: sh tools/run_gates.sh $gate [--threads=N] [--strands=2] [--release] [--allow-dirty]" >&2
      exit 2
      ;;
  esac
done

repo_root=$(git rev-parse --show-toplevel 2>/dev/null)
if [ -z "$repo_root" ]
then
  echo "run_gates: not inside a git repository" >&2
  exit 2
fi

case "$gate" in
  fast | pf | se | config | webui | extd | cgmes | formats | legacy | install | workshops | extended | all | threads | docs | set) ;;
  *)
    echo "usage: sh tools/run_gates.sh fast|pf|se|config|webui|extd|cgmes|formats|legacy|install|workshops|extended|all|threads|docs|set [--threads=N] [--strands=2] [--release] [--allow-dirty]" >&2
    exit 2
    ;;
esac
if [ "$gate" = "docs" ] && [ -n "$julia_threads" ]
then
  echo "run_gates: --threads applies to test profiles, not to the docs build" >&2
  exit 2
fi
if [ "$strands" = "2" ] && [ "$gate" != "extended" ]
then
  echo "run_gates: --strands=2 applies to the extended profile only" >&2
  exit 2
fi
if [ "$gate" = "threads" ] && [ -z "$julia_threads" ]
then
  julia_threads=--threads=4
fi

# --- chain guard ----------------------------------------------------------
# The gate call must be its OWN command; the known failure mode is the
# commit chained into the background gate start, which leaves the gate
# testing a tree that is still being edited. Enforced
# mechanically for the Claude harness, whose bash wrapper carries the full
# user command line in the parent process args (recognizable by the
# shell-snapshots/snapshot-bash marker). Outside that wrapper (a normal
# terminal, CI) the guard SAYS that it is inactive instead of silently
# waving the call through, because parent inspection is platform and shell
# dependent and a silent no-op would look like protection.
parent_cmd=$(ps -o args= -p "$PPID" 2>/dev/null)
if [ -z "$parent_cmd" ]
then
  echo "run_gates: chain guard inactive (cannot inspect the caller via ps)"
else
  case "$parent_cmd" in
    *shell-snapshots/snapshot-bash*)
      user_cmd=${parent_cmd#*eval \'}
      if [ "$user_cmd" = "$parent_cmd" ]
      then
        echo "run_gates: chain guard inactive (harness wrapper without the eval pattern)"
      else
        user_cmd=${user_cmd%%"' < /dev/null"*}
        # exactly one leading `cd <dir>;` or `cd <dir> &&` belongs to the
        # harness idiom and is tolerated (the t branch stops after the
        # first strip, so a second cd still reads as a chain); redirects
        # carry none of the chain tokens checked below
        rest=$(printf '%s' "$user_cmd" | sed -E 's|^[[:space:]]*cd[[:space:]]+[^;&]+;[[:space:]]*||; t; s|^[[:space:]]*cd[[:space:]]+[^;&]+&&[[:space:]]*||')
        case "$rest" in
          *\;* | *"&&"* | *"|"*)
            echo "run_gates: REFUSED: the gate call is part of a shell chain. A gate must run as its own single command, so it never tests a tree that is still being edited: finish and commit the edits first, then start the gate on its own. Offending command:" >&2
            echo "  $user_cmd" >&2
            exit 3
            ;;
        esac
      fi
      ;;
    *)
      echo "run_gates: chain guard inactive (caller is not the Claude harness wrapper)"
      ;;
  esac
fi

# --- docs and workshops run on the COMMIT, not on the working tree ---------
# Both gates build pages that read files. On the working tree they also see
# what is in no repository (ignored and untracked files), so a page that
# reads such a file passes here and fails in a build on a fresh checkout.
# These two gates therefore run in a work tree of HEAD,
# which holds the tracked files and nothing else. What is not committed is
# not part of the run, so every uncommitted change of a tracked file
# refuses the start (with --allow-dirty it is named and left out).
on_commit=no
case "$gate" in
  docs | workshops | set) on_commit=yes ;;
esac
if [ "$on_commit" = "yes" ]
then
  uncommitted=$(git -C "$repo_root" status --porcelain --untracked-files=no)
  if [ -n "$uncommitted" ]
  then
    if [ "$allow_dirty" != "--allow-dirty" ]
    then
      echo "run_gates: the $gate gate runs on the commit, and these tracked files carry uncommitted changes that would not be part of the run. Commit first:" >&2
      echo "$uncommitted" >&2
      exit 1
    fi
    echo "run_gates: --allow-dirty: the $gate gate runs on HEAD; these uncommitted changes are NOT part of the run:"
    echo "$uncommitted"
  fi
fi

if [ "$allow_dirty" != "--allow-dirty" ]
then
  # data/ is included: the suite READS tracked data files (the shipped demo
  # cases, the v1 format fixture), and moving one under a running gate broke
  # the v1 guard mid-run
  dirty=$(git -C "$repo_root" status --porcelain -- src/ app/ test/ data/ Project.toml docs/make.jl docs/Project.toml)
  if [ -n "$dirty" ]
  then
    echo "run_gates: uncommitted changes under src/ test/ data/ or the project files; a gate on a moving tree tests a state that exists nowhere. Commit first, or pass --allow-dirty for a deliberate local check:" >&2
    echo "$dirty" >&2
    exit 1
  fi
fi

# The lock lives in the COMMON git directory, not in "$repo_root/.git". In a
# linked worktree that path is a FILE, so mkdir fails with ENOTDIR, the stale
# takeover below removes nothing and the second mkdir fails the same way: the
# gate then refused every start with "lock takeover raced another start".
# rev-parse --git-common-dir points at the one real git directory from any
# worktree, which also makes the lock SHARED between them, and that is what it
# is for: two Julia first-starts on the same package must not run at once.
git_common=$(git -C "$repo_root" rev-parse --git-common-dir)
case "$git_common" in
  /*) ;;
  *) git_common="$repo_root/$git_common" ;;
esac
lockdir="$git_common/sparlectra_gate.lock"

# --- data/mpower holds TRACKED fixtures and nothing else ---------------------
# The suite gates several legs on case files being present, and that directory
# used to double as the download cache. Whoever had once downloaded a case ran
# more assertions than a fresh worktree or CI, and nothing said so: 7180 here
# against 7121 there (measured). Downloads go to the user cache now
# and the large cases reach the suite only through SPARLECTRA_LARGE_CASES_DIR,
# so anything else in that directory is a leftover that can only skew the run.
# An ERROR, not a warning: a warning is how the drift got in.
stray=$(git -C "$repo_root" ls-files --others --ignored --exclude-standard --directory -- data/mpower)
if [ -n "$stray" ]
then
  echo "run_gates: REFUSED: untracked files in data/mpower/. That directory holds tracked fixtures only; anything else changes which availability-gated legs run and makes the assertion count depend on this machine." >&2
  echo "$stray" >&2
  echo "" >&2
  echo "Move the case files to the directory the suite reads them from, and delete the rest:" >&2
  echo "  mkdir -p \"\${SPARLECTRA_LARGE_CASES_DIR:-\$HOME/sparlectra-cases}\"" >&2
  echo "  mv data/mpower/*.m data/mpower/*.jl data/mpower/*.DAT \"\${SPARLECTRA_LARGE_CASES_DIR:-\$HOME/sparlectra-cases}\"/" >&2
  echo "(keep data/mpower/README.md, warmup_casePST.m and warmup_casePST.measurements.csv: those are tracked)" >&2
  exit 1
fi
if ! mkdir "$lockdir" 2>/dev/null
then
  # stale-lock takeover: the lock records the runner's PID; a lock whose
  # process is gone (crashed shell, killed run) is reported and taken
  # over, a live PID refuses as before
  lockpid=$(cat "$lockdir/pid" 2>/dev/null)
  if [ -n "$lockpid" ] && kill -0 "$lockpid" 2>/dev/null
  then
    echo "run_gates: another gate run holds $lockdir (pid $lockpid, started $(cat "$lockdir/started" 2>/dev/null)). Wait for it." >&2
    exit 1
  fi
  echo "run_gates: stale lock (pid ${lockpid:-unknown} is gone, started $(cat "$lockdir/started" 2>/dev/null)); taking over." >&2
  rm -rf "$lockdir"
  if ! mkdir "$lockdir" 2>/dev/null
  then
    echo "run_gates: lock takeover raced another start; try again." >&2
    exit 1
  fi
fi
date > "$lockdir/started"
echo $$ > "$lockdir/pid"
trap 'rm -rf "$lockdir"' EXIT INT TERM

# the tree the gate runs in: the repository, or a work tree of HEAD. The
# work tree has a fixed place next to the lock, so its compile cache is
# found again by the next run. The manifests are environment state of this
# machine, not content of the commit (they are not tracked), and are copied.
run_root=$repo_root
if [ "$on_commit" = "yes" ]
then
  run_root="$git_common/sparlectra_gate_worktree"
  git -C "$repo_root" worktree remove --force "$run_root" >/dev/null 2>&1
  rm -rf "$run_root"
  git -C "$repo_root" worktree prune
  if ! git -C "$repo_root" worktree add --detach "$run_root" HEAD >/dev/null
  then
    echo "run_gates: could not create the work tree $run_root" >&2
    exit 1
  fi
  # the trap names the work tree by its own variable: set points run_root
  # back at the repository below, and a trap on run_root would then remove
  # the checkout the gate was started from when that is a linked work tree
  gate_wt=$run_root
  trap 'git -C "$repo_root" worktree remove --force "$gate_wt" >/dev/null 2>&1; rm -rf "$lockdir"' EXIT INT TERM
  for manifest in Manifest.toml app/Manifest.toml docs/Manifest.toml
  do
    if [ -f "$repo_root/$manifest" ] && [ ! -f "$run_root/$manifest" ]
    then
      cp "$repo_root/$manifest" "$run_root/$manifest"
    fi
  done
  echo "run_gates: the $gate gate runs on commit $(git -C "$repo_root" rev-parse --short HEAD) in $run_root (tracked files only)"
fi
# set: the docs build and the workshops run in the work tree of HEAD (above),
# the test strands in the repository like every other test profile. The
# dirty check below refuses uncommitted tracked files, so both trees hold
# the same commit; the repository adds the untracked local test data (CGMES
# cache, DTF files) that the work tree lacks.
commit_root=$run_root
[ "$gate" = "set" ] && run_root=$repo_root

# --- one library image per source state, not one per environment ----------
# The gate loads the library from the root environment (the probe below,
# test child processes started with --project=<root>) and from app/ (the
# application probe and the test process, which stacks app/ first); the
# docs build adds docs/. The manifests are not tracked and are resolved one
# by one, so they drift apart. A package in two of them at different
# versions gives the library two compile caches, and every src commit then
# builds the library image twice (194 s
# instead of 131 s for the images of one src commit).
# Compared by git-tree-sha1, so
# standard libraries (no hash in the manifest) and the path entries of the
# two packages are left out. A warning, not a refusal: the run is correct,
# only slower. Sort and join both on the first field in the C locale: a
# whole-line sort in a German locale ignores the blank and puts JSON3
# before JSON, which join then rejects as unsorted and skips.
manifest_hashes() {
  awk '/^\[\[deps\./ { n = $0; gsub(/^\[\[deps\.|\]\]$/, "", n) } /^git-tree-sha1 = / { print n, $3 }' "$1" | LC_ALL=C sort -k1,1
}
drift_envs="app"
[ "$gate" = "docs" ] || [ "$gate" = "set" ] && drift_envs="app docs"
for other in $drift_envs
do
  if [ -f "$run_root/Manifest.toml" ] && [ -f "$run_root/$other/Manifest.toml" ]
  then
    # the two lists go into the lock directory, which the exit trap removes
    manifest_hashes "$run_root/Manifest.toml" > "$lockdir/manifest_root"
    manifest_hashes "$run_root/$other/Manifest.toml" > "$lockdir/manifest_other"
    drift=$(LC_ALL=C join "$lockdir/manifest_root" "$lockdir/manifest_other" | awk '$2 != $3 { printf "%s ", $1 }')
    drift=${drift% }
    if [ -n "$drift" ]
    then
      echo "run_gates: WARNING: Manifest.toml and $other/Manifest.toml resolve different versions of: ${drift}; every src commit builds the library image once per environment. Bring both onto the same versions (for example Pkg.update() in both) in the checkout the worktrees copy their manifests from."
    fi
  fi
done

# --- gates run on package images with the full precompile workload ---------
# Measured (webui test group, same commit, same machine): 185 s
# on the default image (workload off, 14 MB), 69 s on the full image
# (92 MB). The reference times in docs/src/tests.md were taken before the
# default became "off" (4f317a4c), when the workload was always in the
# image. An install keeps the cheap default, and CI (julia-runtest without
# the variable) keeps testing that default image, so a failure that only
# shows without the workload still has a gate. Here the two images are
# rebuilt with the full workload when the loaded one was not built with it:
# the compile cache is keyed by the sources, not by the environment
# variable, so the mode baked into the image (PRECOMPILE_WORKLOAD_MODE) is
# what the probe reads; an image without that constant counts as unknown
# and is rebuilt as well. A rebuilt library image makes the application
# image stale by itself. The rebuild takes about two minutes for both
# images and is announced with its duration; the full image stays the
# newest cache entry, so later sessions load it too until the next source
# change (a bigger image, nothing else changes).
export SPARLECTRA_PRECOMPILE_WORKLOAD=full
image_probe='mode = isdefined(M, :PRECOMPILE_WORKLOAD_MODE) ? M.PRECOMPILE_WORKLOAD_MODE : "unknown"; if mode != "full"; println("run_gates: ", NAME, " image was built with workload ", repr(mode), "; rebuilding it with the full workload (about a minute, later sessions load it too)"); t = @elapsed Base.compilecache(Base.identify_package(NAME)); println("run_gates: ", NAME, " image rebuilt in ", round(Int, t), " s"); end'
probe_roots=$run_root
# set: the work tree of the docs build and the workshops gets its images
# here as well, before the strands start, so no strand precompiles
[ "$gate" = "set" ] && probe_roots="$run_root $commit_root"
for probe_root in $probe_roots
do
  julia --startup-file=no --project="$probe_root" -e "using Sparlectra; const M = Sparlectra; const NAME = \"Sparlectra\"; $image_probe" || exit 1
  julia --startup-file=no --project="$probe_root/app" -e "using SparlectraApp; const M = SparlectraApp; const NAME = \"SparlectraApp\"; $image_probe" || exit 1
done

# --startup-file=no: a gate is not an interactive session, and a personal
# startup.jl usually loads Revise. Measured: loading Revise
# invalidates 472 method instances of the STOCK Julia system image, which the
# run then has to infer again, and it buys a gate nothing. (On the Sparlectra
# sysimage the same load costs 1530 instances; see tools/sysimage_launcher.jl.)
# the docs build of the tree $1, with the two Web UI help checks; returns
# the first non-zero status
docs_build() {
  julia --startup-file=no --project="$1/docs" "$1/docs/make.jl" || return $?
  # the Web UI help registry links into the documentation by anchor; a
  # renamed anchor or a dropped section must fail HERE, in the gate a docs
  # change runs, not in the next extd run
  julia --startup-file=no --project="$1/app" "$1/tools/check_webui_doc_links.jl" || return $?
  # the help pages ship the referenced sections as generated source; an
  # edited section that was not regenerated fails here as well
  julia --startup-file=no --project="$1/app" "$1/tools/generate_webui_help_excerpts.jl" --check
}

# the summary line of a part log: the profile line of a test run, the
# result of the docs build
part_summary() {
  if [ "$1" = "docs" ]
  then
    [ "$2" -eq 0 ] && echo "docs build and help checks green" || echo "docs build or help check failed"
  else
    grep -h "^Sparlectra.jl .* profile |" "$lockdir/part_$1.log" | tail -1 | tr -s ' '
  fi
}

status=0
if [ "$gate" = "docs" ]
then
  docs_build "$run_root"
  status=$?
elif [ "$gate" = "set" ]
then
  # three strands at once on one commit under this lock: set_a (cgmes, se,
  # config), extended_b (pf, webui, extd, formats), and set_c (fast and
  # legacy, then the docs build and the workshops in the work tree of the
  # commit); afterwards the threads profile alone with four threads, so
  # its parallel-vs-serial comparisons do not compete with the strands for
  # the cores. Every part writes into the lock directory; the logs are
  # printed after the end, then one line per part and one summary line.
  # Exit status: the first non-zero part.
  echo "run_gates: set: strands set_a, extended_b and set_c (then docs, workshops) at once, then threads"
  t0=$(date +%s)
  SPARLECTRA_TEST_PROFILE=set_a julia --startup-file=no $julia_threads --project="$run_root" "$run_root/test/runtests.jl" > "$lockdir/part_set_a.log" 2>&1 &
  pid_a=$!
  SPARLECTRA_TEST_PROFILE=extended_b julia --startup-file=no $julia_threads --project="$run_root" "$run_root/test/runtests.jl" > "$lockdir/part_extended_b.log" 2>&1 &
  pid_b=$!
  (
    SPARLECTRA_TEST_PROFILE=set_c julia --startup-file=no $julia_threads --project="$run_root" "$run_root/test/runtests.jl" > "$lockdir/part_set_c.log" 2>&1
    echo $? > "$lockdir/rc_set_c"
    docs_build "$commit_root" > "$lockdir/part_docs.log" 2>&1
    echo $? > "$lockdir/rc_docs"
    SPARLECTRA_TEST_PROFILE=workshops julia --startup-file=no --project="$commit_root" "$commit_root/test/runtests.jl" > "$lockdir/part_workshops.log" 2>&1
    echo $? > "$lockdir/rc_workshops"
  ) &
  pid_c=$!
  trap 'pkill -P $pid_c 2>/dev/null; kill $pid_a $pid_b $pid_c 2>/dev/null; git -C "$repo_root" worktree remove --force "$gate_wt" >/dev/null 2>&1; rm -rf "$lockdir"' INT TERM
  wait $pid_a
  echo $? > "$lockdir/rc_set_a"
  wait $pid_b
  echo $? > "$lockdir/rc_extended_b"
  wait $pid_c
  wall_strands=$(( $(date +%s) - t0 ))
  SPARLECTRA_TEST_PROFILE=threads julia --startup-file=no --threads=4 --project="$run_root" "$run_root/test/runtests.jl" > "$lockdir/part_threads.log" 2>&1
  echo $? > "$lockdir/rc_threads"
  wall=$(( $(date +%s) - t0 ))
  for part in set_a extended_b set_c docs workshops threads
  do
    echo "===== part $part ====="
    cat "$lockdir/part_$part.log"
  done
  for part in set_a extended_b set_c docs workshops threads
  do
    rc=$(cat "$lockdir/rc_$part" 2>/dev/null || echo 1)
    echo "run_gates: part $part rc=$rc $(part_summary "$part" "$rc")"
    [ "$status" -eq 0 ] && status=$rc
  done
  echo "run_gates: set wall ${wall} s (strands ${wall_strands} s, threads $(( wall - wall_strands )) s)"
  for part in set_a extended_b set_c workshops threads
  do
    report_skips "$part" "$lockdir/part_$part.log"
  done
  # the test parts added up; nothing when a part ended without a summary
  for part in set_a extended_b set_c workshops threads
  do
    grep -h "^Sparlectra.jl .* profile |" "$lockdir/part_$part.log" | tail -1
  done | awk -v wall="$wall" '{ for (i = 1; i <= NF; i++) if ($i == "|") { p += $(i + 1); t += $(NF - 1); n++ } } END { if (n == 5) printf "Sparlectra.jl set (fast, extended, workshops, threads) | %d %d %ds\n", p, t, wall }'
else
  # without --threads a profile runs single-threaded, where each parallel
  # site takes its serial loop and a parallel-vs-serial comparison checks
  # nothing threaded; the threads profile defaults to four threads (above)
  [ -n "$julia_threads" ] && echo "run_gates: test process with $julia_threads"
  if [ "$strands" = "2" ]
  then
    # extended as two processes at once (profiles extended_a and
    # extended_b of test/runtests.jl), on the same tree under this one
    # lock; each strand writes into the lock directory and is printed
    # after both ended, then one summary line over both. Exit status: the
    # first non-zero of the two. A killed gate takes both strands along.
    echo "run_gates: extended as two strands (extended_a, extended_b) at once"
    t0=$(date +%s)
    SPARLECTRA_TEST_PROFILE=extended_a julia --startup-file=no $julia_threads --project="$run_root" "$run_root/test/runtests.jl" > "$lockdir/strand_a.log" 2>&1 &
    pid_a=$!
    SPARLECTRA_TEST_PROFILE=extended_b julia --startup-file=no $julia_threads --project="$run_root" "$run_root/test/runtests.jl" > "$lockdir/strand_b.log" 2>&1 &
    pid_b=$!
    trap 'kill $pid_a $pid_b 2>/dev/null; rm -rf "$lockdir"' INT TERM
    wait $pid_a
    status_a=$?
    wait $pid_b
    status_b=$?
    wall=$(( $(date +%s) - t0 ))
    for s in a b
    do
      echo "===== strand extended_$s ====="
      cat "$lockdir/strand_$s.log"
    done
    echo "run_gates: strand extended_a rc=$status_a, extended_b rc=$status_b, wall ${wall} s"
    status=$status_a
    [ "$status" -eq 0 ] && status=$status_b
    for s in a b
    do
      report_skips "extended_$s" "$lockdir/strand_$s.log"
    done
    # the summary lines of both strands, added up; nothing when a strand
    # ended without one (a crash shows in its output above)
    grep -h "^Sparlectra.jl extended_[ab] profile |" "$lockdir/strand_a.log" "$lockdir/strand_b.log" | awk -v wall="$wall" '{ for (i = 1; i <= NF; i++) if ($i == "|") { p += $(i + 1); t += $(NF - 1); n++ } } END { if (n == 2) printf "Sparlectra.jl extended profile (2 strands) | %d %d %ds\n", p, t, wall }'
  else
    # stdout also goes into the lock directory for the skip count; the
    # status of julia, not of tee, through a file (plain sh has no pipefail)
    { SPARLECTRA_TEST_PROFILE=$gate julia --startup-file=no $julia_threads --project="$run_root" "$run_root/test/runtests.jl"; echo $? > "$lockdir/rc_$gate"; } | tee "$lockdir/part_$gate.log"
    status=$(cat "$lockdir/rc_$gate" 2>/dev/null || echo 1)
    report_skips "$gate" "$lockdir/part_$gate.log"
  fi
fi
if [ "$gate" != "docs" ]
then
  echo "run_gates: profiles with skipped/broken tests: ${skip_profiles:-none}"
  run_status=$status
  status=$(final_status "$run_status" "$release" "$skip_profiles")
  [ "$status" != "$run_status" ] && echo "run_gates: --release: skipped or broken tests in ${skip_profiles}; exit $status" >&2
fi
exit $status
