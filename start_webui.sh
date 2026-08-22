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
#
# file: start_webui.sh
# purpose: start the Sparlectra Web UI from this checkout (Linux/macOS).
#          Requires an installed Julia — if none is found, this script points
#          at install_webui.sh, which installs Julia and fetches the latest
#          tagged Sparlectra release before starting the Web UI.

set -e
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

if ! command -v julia >/dev/null 2>&1; then
  echo "Julia is not installed (no 'julia' on PATH)."
  echo "Run   $DIR/install_webui.sh   instead:"
  echo "it installs Julia (via juliaup), fetches the latest tagged Sparlectra release, and starts the Web UI."
  exit 1
fi

# A fresh checkout has no Manifest.toml yet — resolve dependencies once.
if [ ! -f "$DIR/Manifest.toml" ]; then
  echo "First start: resolving Julia dependencies (one-time)..."
  julia --project="$DIR" -e "using Pkg; Pkg.instantiate()"
fi

# Fast start: use the optional PackageCompiler sysimage when it matches the
# running Julia and the checkout Manifest (contract in sysimage_meta.toml,
# written by tools/build_sysimage.jl). The path mirrors
# Sparlectra.default_webui_output_root in src/webui/webui.jl; keep both in
# sync. SPARLECTRA_NO_SYSIMAGE=1 skips the image unconditionally. A stale or
# missing image never blocks the start.
USE_SYSIMAGE=0
if [ "${SPARLECTRA_NO_SYSIMAGE:-0}" != "1" ]; then
  if [ "$(uname)" = "Darwin" ]; then
    SYSIMAGE="$HOME/Library/Application Support/Sparlectra/WebUI/sysimage/sparlectra.dylib"
  else
    STATE_ROOT="${XDG_STATE_HOME:-$HOME/.local/state}"
    SYSIMAGE="$STATE_ROOT/sparlectra/webui/sysimage/sparlectra.so"
  fi
  SYSMETA="$(dirname "$SYSIMAGE")/sysimage_meta.toml"
  if [ -f "$SYSIMAGE" ] && [ -f "$SYSMETA" ]; then
    META_JULIA=$(sed -n 's/^julia_version = "\(.*\)"$/\1/p' "$SYSMETA")
    RUN_JULIA=$(julia --version | awk '{print $3}')
    META_SHA=$(sed -n 's/^manifest_sha256 = "\(.*\)"$/\1/p' "$SYSMETA")
    if command -v sha256sum >/dev/null 2>&1; then
      CUR_SHA=$(sha256sum "$DIR/Manifest.toml" | awk '{print $1}')
    else
      CUR_SHA=$(shasum -a 256 "$DIR/Manifest.toml" | awk '{print $1}')
    fi
    if [ -n "$CUR_SHA" ] && [ "$META_JULIA" = "$RUN_JULIA" ] && [ "$META_SHA" = "$CUR_SHA" ]; then
      USE_SYSIMAGE=1
    else
      echo "Fast start: sysimage stale, starting without it; rebuild via tools/build_sysimage.jl"
    fi
  fi
fi

# Multi-core by default: the threaded surfaces (island solves, short-circuit
# sweeps, N-1 batches) need Julia THREADS, which are fixed at process start.
# "auto" uses all cores of this machine; an explicit user setting wins.
if [ -z "${JULIA_NUM_THREADS:-}" ]; then
  export JULIA_NUM_THREADS=auto
fi

if [ "$USE_SYSIMAGE" = "1" ]; then
  echo "Fast start: using sysimage $SYSIMAGE"
  exec julia "-J$SYSIMAGE" --project="$DIR" "$DIR/start_webui.jl"
fi
exec julia --project="$DIR" "$DIR/start_webui.jl"
