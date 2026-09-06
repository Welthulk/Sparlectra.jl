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
# purpose: run the verification gates (fast / extended / docs) with the two
#          guards that make the 2026-09-03 mixed-state incident impossible
#          instead of unlikely: (1) the run refuses to start while src/ or
#          the project files carry uncommitted changes, because a Julia
#          process precompiling WHILE sources keep changing tests a state
#          that exists nowhere (that incident cost one segfault and one
#          25-minute mixed-state run); (2) a lock directory serializes gate
#          runs, so no second Julia first-start precompiles the same
#          package concurrently. --allow-dirty skips guard 1 for a
#          deliberate local check; the lock is never skipped. Plain sh.
# usage:   sh tools/run_gates.sh fast|extended|docs [--allow-dirty]
#          (SPARLECTRA_LARGE_CASES_DIR and other env pass through)

gate=$1
allow_dirty=$2

repo_root=$(git rev-parse --show-toplevel 2>/dev/null)
if [ -z "$repo_root" ]
then
  echo "run_gates: not inside a git repository" >&2
  exit 2
fi

case "$gate" in
  fast | extended | docs) ;;
  *)
    echo "usage: sh tools/run_gates.sh fast|extended|docs [--allow-dirty]" >&2
    exit 2
    ;;
esac

# --- chain guard (maintainer 2026-09-04) -------------------------------------
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

if [ "$allow_dirty" != "--allow-dirty" ]
then
  # data/ is included: the suite READS tracked data files (the shipped demo
  # cases, the v1 format fixture), and moving one under a running gate broke
  # the v1 guard mid-run on 2026-09-03
  dirty=$(git -C "$repo_root" status --porcelain -- src/ test/ data/ Project.toml docs/make.jl docs/Project.toml)
  if [ -n "$dirty" ]
  then
    echo "run_gates: uncommitted changes under src/ test/ data/ or the project files; a gate on a moving tree tests a state that exists nowhere. Commit first, or pass --allow-dirty for a deliberate local check:" >&2
    echo "$dirty" >&2
    exit 1
  fi
fi

lockdir="$repo_root/.git/sparlectra_gate.lock"
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

status=0
if [ "$gate" = "fast" ]
then
  julia --project="$repo_root" "$repo_root/test/runtests.jl"
  status=$?
elif [ "$gate" = "extended" ]
then
  SPARLECTRA_TEST_PROFILE=extended julia --project="$repo_root" "$repo_root/test/runtests.jl"
  status=$?
else
  julia --project="$repo_root/docs" "$repo_root/docs/make.jl"
  status=$?
fi
exit $status
