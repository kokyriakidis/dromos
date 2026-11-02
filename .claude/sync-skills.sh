#!/usr/bin/env bash
# Pulls Modular's Mojo skills (github.com/modular/skills) into .claude/skills, replacing the copies
# there, and records the upstream commit in .claude/skills/UPSTREAM. Its MAX model skills (import,
# serve, profile, evaluate models) do not apply to this package and are left out.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
skills=(mojo-syntax mojo-gpu-fundamentals mojo-python-interop new-modular-project)
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
git clone --quiet --depth 1 https://github.com/modular/skills "$work/skills"
for skill in "${skills[@]}"; do
    rm -rf "${here:?}/skills/$skill"
    cp -R "$work/skills/$skill" "$here/skills/$skill"
done
cp "$work/skills/LICENSE" "$here/skills/LICENSE"
git -C "$work/skills" rev-parse HEAD > "$here/skills/UPSTREAM"
echo "Mojo skills synced from modular/skills $(cut -c1-12 "$here/skills/UPSTREAM")"
