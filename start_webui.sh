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

exec julia --project="$DIR" "$DIR/start_webui.jl"
