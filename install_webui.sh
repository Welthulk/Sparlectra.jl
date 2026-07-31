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
# file: install_webui.sh
# purpose: install-and-start for the Sparlectra Web UI (Linux/macOS):
#          installs Julia via the official juliaup installer when missing,
#          obtains Sparlectra at its latest tagged release (uses this
#          checkout when the script lies inside one, otherwise clones or
#          downloads next to the script), then starts the Web UI.

set -e
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_URL="https://github.com/Welthulk/Sparlectra.jl"

# --- 1. Julia ---------------------------------------------------------------
if ! command -v julia >/dev/null 2>&1; then
  echo "Julia not found — installing via the official juliaup installer..."
  if command -v curl >/dev/null 2>&1; then
    curl -fsSL https://install.julialang.org | sh -s -- --yes
  elif command -v wget >/dev/null 2>&1; then
    wget -qO- https://install.julialang.org | sh -s -- --yes
  else
    echo "Neither curl nor wget is available. Install Julia manually: https://julialang.org/downloads/"
    exit 1
  fi
  # juliaup puts julia into ~/.juliaup/bin; make it visible in THIS session.
  PATH="$HOME/.juliaup/bin:$PATH"
  export PATH
  if ! command -v julia >/dev/null 2>&1; then
    echo "Julia was installed but is not on PATH yet. Open a new terminal and run start_webui.sh."
    exit 1
  fi
fi
echo "Julia: $(julia --version)"

# --- 2. Sparlectra at the latest tagged release -----------------------------
if grep -q '^name = "Sparlectra"' "$DIR/Project.toml" 2>/dev/null; then
  # The script lies inside a Sparlectra checkout — use it. Move a CLEAN git
  # checkout to the latest release tag; a dirty tree is left untouched so no
  # local work is ever discarded.
  SPARLECTRA_DIR="$DIR"
  if [ -d "$DIR/.git" ] && command -v git >/dev/null 2>&1; then
    git -C "$DIR" fetch --tags --quiet || true
    LATEST=$(git -C "$DIR" tag --sort=-v:refname | head -n 1)
    if [ -n "$LATEST" ]; then
      if [ -z "$(git -C "$DIR" status --porcelain)" ]; then
        echo "Checking out latest Sparlectra release: $LATEST"
        git -C "$DIR" checkout --quiet "$LATEST"
      else
        echo "Working tree has local changes — keeping the current state (latest release would be $LATEST)."
      fi
    fi
  fi
else
  # Standalone use: put a release copy next to the script.
  SPARLECTRA_DIR="$DIR/Sparlectra"
  if [ ! -d "$SPARLECTRA_DIR" ]; then
    if command -v git >/dev/null 2>&1; then
      echo "Cloning Sparlectra..."
      git clone --quiet "$REPO_URL" "$SPARLECTRA_DIR"
      LATEST=$(git -C "$SPARLECTRA_DIR" tag --sort=-v:refname | head -n 1)
      [ -n "$LATEST" ] && echo "Checking out latest release: $LATEST" && git -C "$SPARLECTRA_DIR" checkout --quiet "$LATEST"
    elif command -v curl >/dev/null 2>&1; then
      # No git: download the tarball of the latest tag via the GitHub API.
      LATEST=$(curl -fsSL "https://api.github.com/repos/Welthulk/Sparlectra.jl/tags" | grep -m 1 '"name"' | sed 's/.*"name": *"\([^"]*\)".*/\1/')
      [ -n "$LATEST" ] || { echo "Could not determine the latest release tag."; exit 1; }
      echo "Downloading Sparlectra $LATEST..."
      curl -fsSL "$REPO_URL/archive/refs/tags/$LATEST.tar.gz" | tar -xz -C "$DIR"
      mv "$DIR"/Sparlectra.jl-* "$SPARLECTRA_DIR"
    else
      echo "Neither git nor curl is available — cannot obtain Sparlectra."
      exit 1
    fi
  else
    echo "Using existing copy at $SPARLECTRA_DIR"
  fi
fi

# --- 3. Dependencies + start ------------------------------------------------
echo "Resolving Julia dependencies..."
julia --project="$SPARLECTRA_DIR" -e "using Pkg; Pkg.instantiate()"
exec julia --project="$SPARLECTRA_DIR" "$SPARLECTRA_DIR/start_webui.jl"
