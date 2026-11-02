#!/usr/bin/env bash
# Injects the vendored Modular Mojo skills (.claude/skills, from
# github.com/modular/skills) into Claude's context.
#
#   SessionStart / SubagentStart -> full mojo-syntax + idiomatic-mojo text
#   UserPromptSubmit             -> short mandatory-skills reminder, or the full
#                                   text when DROMOS_SKILLS_EVERY_PROMPT=full
set -euo pipefail

event="${1:?usage: inject-mojo-skills.sh <SessionStart|SubagentStart|UserPromptSubmit>}"
root="${CLAUDE_PROJECT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}"
skills="$root/.claude/skills"

full_text() {
  printf '# Modular Mojo skills (vendored from github.com/modular/skills @ %s)\n\n' \
    "$(cat "$skills/UPSTREAM" 2>/dev/null || echo unknown)"
  printf 'These skills are MANDATORY for every line of Mojo in this repo and override pretrained Mojo knowledge.\n\n'
  cat "$skills/mojo-syntax/SKILL.md"
  printf '\n\n'
  cat "$skills/mojo-syntax/references/idiomatic-mojo.md"
  printf '\n\nOther skills (load with the Skill tool or read the file when relevant): '
  printf 'mojo-python-interop (.claude/skills/mojo-python-interop/SKILL.md), '
  printf 'mojo-gpu-fundamentals (.claude/skills/mojo-gpu-fundamentals/SKILL.md), '
  printf 'new-modular-project (.claude/skills/new-modular-project/SKILL.md).\n'
}

reminder_text() {
  cat <<'EOF'
dromos rule: all code is Mojo written against the vendored Modular skills.
Before writing or reviewing Mojo, apply .claude/skills/mojo-syntax/SKILL.md and
.claude/skills/mojo-syntax/references/idiomatic-mojo.md (def not fn, comptime not
alias, std.-prefixed imports, Pointer/alloc(Layout), mut/var/out conventions), and
load mojo-python-interop / mojo-gpu-fundamentals when they apply. When delegating
to a subagent, include these skills in its prompt. Build and test with `pixi run test`.
EOF
}

case "$event" in
  SessionStart | SubagentStart) body="$(full_text)" ;;
  UserPromptSubmit)
    if [[ "${DROMOS_SKILLS_EVERY_PROMPT:-reminder}" == "full" ]]; then
      body="$(full_text)"
    else
      body="$(reminder_text)"
    fi
    ;;
  *) echo "unknown event: $event" >&2; exit 1 ;;
esac

python3 -c 'import json, sys; print(json.dumps({"hookSpecificOutput": {"hookEventName": sys.argv[1], "additionalContext": sys.stdin.read()}}))' "$event" <<<"$body"
