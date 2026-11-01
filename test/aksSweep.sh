#!/usr/bin/env bash
# End-to-end check of the Zn[x] AKS program against a probable-prime oracle.
# usage: test/aksSweep.sh [max_n]      (default 1000; run from the repo root after `make`)
set -euo pipefail
max=${1:-1000}
root=$(pwd)
prog="$root/devVR/LenstraZnx.out"
oracle="$root/test/probPrime.out"
[[ -x "$prog" && -x "$oracle" ]] || { echo "build first: make all test/probPrime.out"; exit 2; }

# The program writes a timestamped CSV per run into ./logs/Znx; keep those out of the repo.
work=$(mktemp -d); mkdir -p "$work/logs/Znx"; trap 'rm -rf "$work"' EXIT
cd "$work"

fail=0
for ((n = 2; n <= max; ++n)); do
    expect=$(echo "$n" | "$oracle")
    out=$(echo "$n" | "$prog" 2>&1)
    if grep -q "^$n is prime\." <<< "$out"; then got=1
    elif grep -q "^$n is not prime\." <<< "$out"; then got=0
    else got=?; fi
    if [[ "$got" != "$expect" ]]; then echo "MISMATCH n=$n expected=$expect got=$got"; fail=1; fi
done
[[ $fail == 0 ]] && echo "aksSweep: 2..$max all correct"
exit $fail
