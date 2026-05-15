#!/usr/bin/env bash
set -euo pipefail

COMPSIZE=$1
RUNS=$2

run_and_report () {
    local label="$1"
    local outfile="$2"
    shift 2

    echo "========================================"
    echo "$label"
    echo "Command: $*"
    echo "Output file: $outfile"
    echo "Runs: $RUNS"
    echo "========================================"

    total_time=0
    total_mem=0

    for ((i=1; i<=RUNS; i++)); do
        echo "Run $i/$RUNS..."

        output=$(
            { /usr/bin/time -v "$@" > "$outfile"; } 2>&1
        )

        # Extract elapsed time
        raw_time=$(echo "$output" | grep "Elapsed (wall clock) time" | awk -F': ' '{print $2}')

        # Convert time to seconds
        seconds=$(python3 - <<PY
t = "$raw_time"

parts = t.split(':')

if len(parts) == 2:
    m, s = parts
    print(float(m)*60 + float(s))
elif len(parts) == 3:
    h, m, s = parts
    print(float(h)*3600 + float(m)*60 + float(s))
else:
    print(float(t))
PY
)

        # Extract RSS in KB
        rss=$(echo "$output" | grep "Maximum resident set size" | awk '{print $6}')

        total_time=$(python3 - <<PY
print($total_time + $seconds)
PY
)

        total_mem=$((total_mem + rss))
    done

    avg_time=$(python3 - <<PY
print($total_time / $RUNS)
PY
)

    avg_mem=$(python3 - <<PY
print($total_mem / $RUNS / 1024)
PY
)

    echo
    echo "AVERAGES OVER $RUNS RUNS"
    echo "Average runtime (s): $avg_time"
    echo "Average max RSS (MB): $avg_mem"
    echo
}

run_and_report \
    "Compressed FM-index queries of len 20" \
    comp_20queries_results.txt \
    ../../src/query_index comp_covid$COMPSIZE.ind ../../data/covid20queries.txt

run_and_report \
    "Plain FM-index queries of len 20" \
    plain_20queries_results.txt \
    ../../src/query_plain_index plain_covid.ind ../../data/covid20queries.txt

run_and_report \
    "Compressed FM-index queries of len 50" \
    comp_50queries_results.txt \
    ../../src/query_index comp_covid$COMPSIZE.ind ../../data/covid50queries.txt

run_and_report \
    "Plain FM-index queries of len 50" \
    plain_50queries_results.txt \
    ../../src/query_plain_index plain_covid.ind ../../data/covid50queries.txt

run_and_report \
    "Compressed FM-index random queries" \
    comp_randomqueries_results.txt \
    ../../src/query_index comp_covid$COMPSIZE.ind ../../data/randomqueries.txt

run_and_report \
    "Plain FM-index random queries" \
    plain_randomqueries_results.txt \
    ../../src/query_plain_index plain_covid.ind ../../data/randomqueries.txt