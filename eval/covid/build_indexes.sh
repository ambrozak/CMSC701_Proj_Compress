#!/usr/bin/env bash
set -euo pipefail

COMPSIZE=$1

run_and_report () {
    local label="$1"
    shift

    echo "========================================"
    echo "$label"
    echo "Command: $*"
    echo "========================================"

    # capture stderr from /usr/bin/time
    output=$(
        { /usr/bin/time -v "$@" > /dev/null; } 2>&1
    )

    wall=$(echo "$output" | grep "Elapsed (wall clock) time")
    mem=$(echo "$output" | grep "Maximum resident set size")

    echo "$wall"
    echo "$mem"
    echo
}

run_and_report \
    "Compressed FM-index build" \
    ../../src/build_index ../../data/covid.fasta "$COMPSIZE" comp_covid$COMPSIZE.ind

run_and_report \
    "Plain FM-index build" \
    ../../src/build_plain_index ../../data/covid.fasta plain_covid.ind