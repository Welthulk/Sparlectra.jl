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
#          downloads next to the script), offers the update when an
#          existing copy is older than the latest release, optionally
#          builds the fast-start sysimage after asking the user, then
#          starts the Web UI.
#          Also works piped straight from the docs (no git, no GitHub
#          checkout needed):
#            curl -fsSL https://raw.githubusercontent.com/Welthulk/Sparlectra.jl/main/install_webui.sh | sh
#          When piped, Sparlectra is placed in a Sparlectra/ folder inside
#          the CURRENT directory; run the one-liner where it should live.

set -e
# When executed from a file, everything is anchored next to the script.
# When piped (curl | sh), $0 is the shell itself; anchoring at its dirname
# would point at /bin, so fall back to the current directory instead.
case "$0" in
  *install_webui.sh) DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd) ;;
  *) DIR=$(pwd) ;;
esac
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
    # Update check: an existing copy stays at the release it was installed
    # at, and re-running this installer is the documented way to update.
    # Detect a newer tagged release and offer the update (default yes;
    # SPARLECTRA_UPDATE=1/0 decides without asking). A failed check or a
    # failed download never blocks the start with the existing copy.
    INSTALLED=$(sed -n 's/^version = "\(.*\)"$/\1/p' "$SPARLECTRA_DIR/Project.toml" 2>/dev/null | head -n 1)
    if [ -d "$SPARLECTRA_DIR/.git" ] && command -v git >/dev/null 2>&1; then
      # Cloned copy: same tag mechanism as the in-checkout path above.
      git -C "$SPARLECTRA_DIR" fetch --tags --quiet 2>/dev/null || true
      LATEST=$(git -C "$SPARLECTRA_DIR" tag --sort=-v:refname | head -n 1)
      CURRENT_TAG=$(git -C "$SPARLECTRA_DIR" describe --tags --exact-match 2>/dev/null || echo "")
      if [ -n "$LATEST" ] && [ "$LATEST" != "$CURRENT_TAG" ]; then
        if [ -n "$(git -C "$SPARLECTRA_DIR" status --porcelain 2>/dev/null)" ]; then
          echo "A newer release ($LATEST) exists, but the copy has local changes - keeping the current state."
        else
          DO_UPDATE="${SPARLECTRA_UPDATE:-}"
          if [ -z "$DO_UPDATE" ] && [ -r /dev/tty ]; then
            printf "Installed: %s, latest release: %s. Update now? [Y/n] " "${CURRENT_TAG:-v$INSTALLED}" "$LATEST"
            read -r DO_UPDATE < /dev/tty 2>/dev/null || DO_UPDATE=""
          fi
          case "$DO_UPDATE" in
            n|N|no|NO|0)
              echo "Keeping ${CURRENT_TAG:-v$INSTALLED}."
              ;;
            *)
              echo "Updating to $LATEST..."
              git -C "$SPARLECTRA_DIR" checkout --quiet "$LATEST"
              echo "Note: an existing fast-start sysimage is stale after an update; rebuild it below or later."
              ;;
          esac
        fi
      fi
    elif command -v curl >/dev/null 2>&1; then
      # Tarball copy (no .git): compare the installed Project.toml version
      # against the newest tag from the GitHub API.
      LATEST=$(curl -fsSL "https://api.github.com/repos/Welthulk/Sparlectra.jl/tags" 2>/dev/null | grep -m 1 '"name"' | sed 's/.*"name": *"\([^"]*\)".*/\1/')
      if [ -z "$LATEST" ]; then
        echo "Update check skipped (could not reach the GitHub API)."
      elif [ "v$INSTALLED" = "$LATEST" ]; then
        echo "Already at the latest release (v$INSTALLED)."
      else
        DO_UPDATE="${SPARLECTRA_UPDATE:-}"
        if [ -z "$DO_UPDATE" ] && [ -r /dev/tty ]; then
          printf "Installed: v%s, latest release: %s. Update now? [Y/n] " "$INSTALLED" "$LATEST"
          read -r DO_UPDATE < /dev/tty 2>/dev/null || DO_UPDATE=""
        fi
        case "$DO_UPDATE" in
          n|N|no|NO|0)
            echo "Keeping v$INSTALLED."
            ;;
          *)
            echo "Downloading Sparlectra $LATEST..."
            # Extract into a temp dir first: the existing copy is only
            # replaced after a complete download, and it survives as
            # Sparlectra.old (one backup slot) so the update is reversible.
            UPDATE_TMP=$(mktemp -d "$DIR/.sparlectra_update.XXXXXX")
            if curl -fsSL "$REPO_URL/archive/refs/tags/$LATEST.tar.gz" | tar -xz -C "$UPDATE_TMP" 2>/dev/null && [ -d "$UPDATE_TMP"/Sparlectra.jl-* ]; then
              rm -rf "$DIR/Sparlectra.old" 2>/dev/null || true
              mv "$SPARLECTRA_DIR" "$DIR/Sparlectra.old"
              mv "$UPDATE_TMP"/Sparlectra.jl-* "$SPARLECTRA_DIR"
              rmdir "$UPDATE_TMP" 2>/dev/null || true
              echo "Updated to $LATEST. The previous copy is kept at $DIR/Sparlectra.old (delete it once the new version runs)."
              echo "Note: an existing fast-start sysimage is stale after an update; rebuild it below or later."
            else
              rm -rf "$UPDATE_TMP" 2>/dev/null || true
              echo "Update download failed - continuing with v$INSTALLED."
            fi
            ;;
        esac
      fi
    fi
  fi
fi

# --- 3. Dependencies ---------------------------------------------------------
echo "Resolving Julia dependencies..."
julia --project="$SPARLECTRA_DIR" -e "using Pkg; Pkg.instantiate()"

# --- 4. Optional fast-start sysimage ----------------------------------------
# One ahead-of-time compiled image removes Julia's JIT warm-up from every
# Web UI start (see docs: Fast Start). The build is optional, runs once and
# takes roughly 6 to 20 minutes. The answer is read from the terminal even
# when the script itself arrives via a pipe (curl | sh); without a terminal
# the default is "no". SPARLECTRA_BUILD_SYSIMAGE=1 or =0 decides without
# asking (unattended installs).
BUILD_IMG="${SPARLECTRA_BUILD_SYSIMAGE:-}"
if [ -z "$BUILD_IMG" ] && [ -r /dev/tty ]; then
  printf "Build the optional fast-start sysimage now? One-time build of 6-20 minutes; later Web UI starts skip the warm-up. [y/N] "
  read -r BUILD_IMG < /dev/tty 2>/dev/null || BUILD_IMG=""
fi
case "$BUILD_IMG" in
  y|Y|yes|YES|1)
    echo "Building the fast-start sysimage (this can take a while)..."
    # A failed build must never block the installation: the launcher only
    # uses an image whose metadata validates, so starting without one is
    # always safe.
    if ! julia "$SPARLECTRA_DIR/tools/build_sysimage.jl"; then
      echo "Sysimage build failed - continuing without it. Retry later with:"
      echo "  julia \"$SPARLECTRA_DIR/tools/build_sysimage.jl\""
    fi
    ;;
  *)
    echo "Skipping the sysimage build. Build it anytime with:"
    echo "  julia \"$SPARLECTRA_DIR/tools/build_sysimage.jl\""
    echo "or from the Web UI: Fast start page, button 'Build fast-start image'."
    ;;
esac

# --- 5. Optional desktop launcher ---------------------------------------------
# A launcher on the desktop / in the application menu restarts the Web UI
# without hunting for start_webui.sh in the install folder. Same prompt
# conventions as the sysimage question above, but default YES: the launcher
# costs nothing and deleting one file reverses it.
# SPARLECTRA_CREATE_SHORTCUT=1/0 decides without asking (unattended installs).
create_launcher() {
  START_SH="$SPARLECTRA_DIR/start_webui.sh"
  if [ "$(uname)" = "Darwin" ]; then
    # .desktop files mean nothing on macOS; a Desktop symlink is the simplest
    # useful equivalent (a .app bundle is deliberately out of scope).
    if [ -d "$HOME/Desktop" ] && ln -sf "$START_SH" "$HOME/Desktop/Sparlectra Web UI" 2>/dev/null; then
      echo "Desktop link created: ~/Desktop/Sparlectra Web UI"
      return 0
    fi
    return 1
  fi
  # Terminal=true: the Web UI logs to the terminal and is stopped there.
  LAUNCHER_CONTENT="[Desktop Entry]
Type=Application
Name=Sparlectra Web UI
Comment=Start the Sparlectra Web UI
Exec=\"$START_SH\"
Path=$SPARLECTRA_DIR
Terminal=true"
  APPS_DIR="$HOME/.local/share/applications"
  DESK_DIR=""
  if command -v xdg-user-dir >/dev/null 2>&1; then
    DESK_DIR=$(xdg-user-dir DESKTOP 2>/dev/null) || DESK_DIR=""
  fi
  [ -d "$DESK_DIR" ] || DESK_DIR="$HOME/Desktop"
  WROTE=""
  # Application menu entry (always attempted); existing files are replaced so
  # a re-install points the launcher at the current location.
  if mkdir -p "$APPS_DIR" 2>/dev/null && printf '%s\n' "$LAUNCHER_CONTENT" > "$APPS_DIR/sparlectra-webui.desktop" 2>/dev/null; then
    WROTE="menu"
    echo "Application-menu launcher created: $APPS_DIR/sparlectra-webui.desktop"
  fi
  # Desktop copy only when a desktop directory exists (headless systems skip
  # this silently).
  if [ -d "$DESK_DIR" ] && printf '%s\n' "$LAUNCHER_CONTENT" > "$DESK_DIR/sparlectra-webui.desktop" 2>/dev/null; then
    chmod +x "$DESK_DIR/sparlectra-webui.desktop" 2>/dev/null || true
    WROTE="${WROTE:+$WROTE+}desktop"
    echo "Desktop launcher created: $DESK_DIR/sparlectra-webui.desktop"
    echo "Note: some desktops (GNOME) need a one-time right-click 'Allow Launching' on the icon."
  fi
  if [ -z "$WROTE" ]; then
    # No XDG environment at all: a home-directory symlink still gives a
    # stable entry point.
    if ln -sf "$START_SH" "$HOME/sparlectra-webui" 2>/dev/null; then
      echo "No desktop environment found - created symlink: ~/sparlectra-webui -> start_webui.sh"
    else
      return 1
    fi
  fi
  return 0
}
MAKE_SHORTCUT="${SPARLECTRA_CREATE_SHORTCUT:-}"
if [ -z "$MAKE_SHORTCUT" ] && [ -r /dev/tty ]; then
  printf "Create a desktop launcher to start the Web UI? [Y/n] "
  read -r MAKE_SHORTCUT < /dev/tty 2>/dev/null || MAKE_SHORTCUT=""
fi
case "$MAKE_SHORTCUT" in
  n|N|no|NO|0)
    echo "Skipping the desktop launcher."
    ;;
  *)
    # A failed launcher creation must never block the installation or the
    # Web UI start (the function call sits in a condition, so set -e does
    # not abort on failures inside it).
    if ! create_launcher; then
      echo "Could not create a launcher - start the Web UI manually with:"
      echo "  \"$SPARLECTRA_DIR/start_webui.sh\""
    fi
    ;;
esac

# --- 6. Start ----------------------------------------------------------------
exec julia --project="$SPARLECTRA_DIR" "$SPARLECTRA_DIR/start_webui.jl"