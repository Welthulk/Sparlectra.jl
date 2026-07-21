#!/bin/bash
# PreToolUse hook (Bash matcher). Blocks git add/commit invocations that would
# add or commit a task_*.md / task-*.md file — these are private Claude Code
# scratch task descriptions and must never enter the Sparlectra repo history.
# Plain `git add task_*.md` is already refused by the task_*.md .gitignore
# rule; this hook additionally catches `git add -f ...` (bypasses .gitignore)
# and `git commit` when such a file is already staged for any reason.
input=$(cat)
raw_cmd=$(printf '%s' "$input" | jq -r '.tool_input.command // empty')
# Fold actual newlines (multi-line Bash tool commands) into ";" so the same
# separator-anchored regex below covers both single- and multi-line input.
cmd=$(printf '%s' "$raw_cmd" | tr '\n' ';')

# Anchor on actual command position (start of string, or right after a shell
# separator ; & | — which also covers && and ||) so this does NOT fire on
# "git add"/"git commit" merely appearing inside a quoted string, echo
# argument, comment, or commit message — only on an actual invocation.
sep='(^|[;&|])[[:space:]]*'
git_commit_re="${sep}git[[:space:]]+commit\\b"
git_add_re="${sep}git[[:space:]]+add\\b"

echo "$cmd" | grep -qE "(${git_commit_re})|(${git_add_re})" || exit 0

blocked=0
reason=""

if echo "$cmd" | grep -qE "$git_commit_re"; then
  staged=$(git diff --cached --name-only 2>/dev/null | grep -E '(^|/)task[_-][^/]*\.md$')
  if [ -n "$staged" ]; then
    blocked=1
    reason="staged task file(s) would be committed: $(echo "$staged" | tr '\n' ' ')"
  fi
fi

if echo "$cmd" | grep -qE "$git_add_re"; then
  # Isolate just the matched "git add ..." segment (up to the next separator)
  # so unrelated later text in the same multi-command invocation — e.g. a
  # commit message heredoc that mentions "task_*.md" in prose — can't trigger
  # a false positive below.
  add_segment=$(echo "$cmd" | grep -oE "${git_add_re}[^;&|]*" | head -1)

  if echo "$add_segment" | grep -qE 'task[_-][^[:space:]]*\.md'; then
    blocked=1
    reason="${reason:+$reason; }git add explicitly names a task_*.md file"
  fi
  if echo "$add_segment" | grep -qE '(([[:space:]]|^)\.([[:space:]]|$)|--all\b|[[:space:]]-A\b)'; then
    untracked=$(git status --porcelain 2>/dev/null | awk '{print $2}' | grep -E '(^|/)task[_-][^/]*\.md$')
    if [ -n "$untracked" ]; then
      blocked=1
      reason="${reason:+$reason; }git add . / -A would stage task file(s): $(echo "$untracked" | tr '\n' ' ')"
    fi
  fi
fi

if [ "$blocked" = "1" ]; then
  jq -n --arg reason "Blocked: task_*.md files are private Claude Code scratch files and must never be committed to this repo. $reason" \
    '{hookSpecificOutput: {hookEventName: "PreToolUse", permissionDecision: "deny", permissionDecisionReason: $reason}}'
fi
exit 0
