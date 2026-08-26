#!/bin/sh
# The four reproduction gates, defined once.
#
# Before 2026-08-26 the four commands appeared in three places -- .github/workflows/gates.yml,
# .github/workflows/readme-install.yml and the Dockerfile -- and a change had to be made in all
# three. This file is the single definition; those three call it.
#
# TWO MODES, because the callers need different things and collapsing them would lose one:
#
#   run-gates.sh <1|2|3|4>   run ONE gate, exit with its status.
#                            Used by gates.yml, which runs four separately-named steps so a CI
#                            run always reports which gate failed. A script that ran all four
#                            internally and returned a single status would take that naming away.
#
#   run-gates.sh all         run all four UNCHAINED, print a per-gate summary, print a
#                            `!!! GATE n FAILED` marker for each failure, and exit non-zero if
#                            any failed. Used by the Dockerfile, where the build must certify
#                            "all four ran AND all four passed".
#
# Both properties in `all` mode pull against each other and the end-of-run `exit` is what holds
# them together. Simply dropping the `&&` that this replaced would give a GREEN build when gate 1
# fails and gate 4 passes, because a shell takes the last command's status -- a loud failure
# converted into a silent one, which is worse than the chaining it replaced.
#
# INTERPRETER. Resolved rather than assumed, because the three callers differ:
#   - gates.yml runs on a GitHub runner with actions/setup-python and no .venv  -> `python`
#   - readme-install.yml runs in a bare container after the README's own install -> `.venv/bin/python`
#   - the Docker image has /cellwarp/.venv                                       -> `.venv/bin/python`
# Override with CELLWARP_PY if you need a specific interpreter.

set -u

# Repo root, resolved from this script's own location so the caller's cwd does not matter.
SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(dirname -- "$SCRIPT_DIR")
cd "$REPO_ROOT" || exit 1

if [ -n "${CELLWARP_PY:-}" ]; then
    PY="$CELLWARP_PY"
elif [ -x .venv/bin/python ]; then
    PY=.venv/bin/python
elif command -v python >/dev/null 2>&1; then
    PY=python
else
    PY=python3
fi

gate_1() { "$PY" reproduce/validate.py; }
gate_2() { "$PY" -m pytest -q; }
gate_3() { "$PY" scripts/build_submission_packet.py --verify; }
gate_4() { md5sum -c reproduce/MANUSCRIPT_MD5; }

usage() {
    echo "usage: $0 <1|2|3|4|all>" >&2
    exit 2
}

[ $# -eq 1 ] || usage

case "$1" in
    1) echo "### GATE 1: validate.py ###";      gate_1; exit $? ;;
    2) echo "### GATE 2: pytest -q ###";        gate_2; exit $? ;;
    3) echo "### GATE 3: packet --verify ###";  gate_3; exit $? ;;
    4) echo "### GATE 4: manuscript md5 ###";   gate_4; exit $? ;;
    all) ;;
    *) usage ;;
esac

# ---- `all` mode: unchained, every gate reports, result decided at the end ----
echo "### GATE 1: validate.py ###";      gate_1; G1=$?
echo "### GATE 2: pytest -q ###";        gate_2; G2=$?
echo "### GATE 3: packet --verify ###";  gate_3; G3=$?
echo "### GATE 4: manuscript md5 ###";   gate_4; G4=$?

echo "### GATE SUMMARY (exit code per gate) ###"
echo "###   gate 1  validate.py        : $G1"
echo "###   gate 2  pytest -q          : $G2"
echo "###   gate 3  packet --verify    : $G3"
echo "###   gate 4  manuscript md5 pin : $G4"

# One greppable marker per failure. These do not appear at all in a green run, so a reader
# scanning a long build log finds the failing gate before they find its traceback.
[ "$G1" -eq 0 ] || echo "!!! GATE 1 FAILED (exit $G1) !!!"
[ "$G2" -eq 0 ] || echo "!!! GATE 2 FAILED (exit $G2) !!!"
[ "$G3" -eq 0 ] || echo "!!! GATE 3 FAILED (exit $G3) !!!"
[ "$G4" -eq 0 ] || echo "!!! GATE 4 FAILED (exit $G4) !!!"

if [ "$G1" -eq 0 ] && [ "$G2" -eq 0 ] && [ "$G3" -eq 0 ] && [ "$G4" -eq 0 ]; then
    echo "### RESULT: all four gates ran and passed ###"
    exit 0
else
    echo "### RESULT: BUILD FAILS -- see the !!! GATE lines above ###"
    exit 1
fi
