#!/bin/bash
# run_local.sh — runs the full benchmark suite against the datasets in data/
# and appends sux-bench (Rust Elias-Fano) results alongside the C++ ones.

set -euo pipefail

DIR="results"
DATA="data"
EXE="build/my_benchmark"
SUX="sux-bench/target/release/sux-bench"

OPTIONS="--meva --plai --elia --rrrv --rlev --lave --lavo --encv --ds2i --s18v"

mkdir -p "$DIR"
rm -f "$DIR/DNA_5GRAM_URL_comparison.csv" \
      "$DIR/GOV2_100K-1M_comparison.csv" "$DIR/GOV2_1M-10M_comparison.csv" \
      "$DIR/GOV2_10M-_comparison.csv" "$DIR/GOV2_averages.csv" \
      "$DIR/comparison.csv" "$DIR/sux_ef_comparison.csv" "$DIR/comparison.err"

echo "=== DNA / 5GRAM / URL ==="
"$EXE" $OPTIONS \
    "$DATA/DNA_1"   "$DATA/DNA_2"   "$DATA/DNA_3" \
    "$DATA/5GRAM_1" "$DATA/5GRAM_2" "$DATA/5GRAM_3" \
    "$DATA/URL_1"   "$DATA/URL_2"   "$DATA/URL_3" \
    > "$DIR/DNA_5GRAM_URL_comparison.csv" \
    2>> "$DIR/comparison.err"

echo "=== GOV2 100K-1M ==="
"$EXE" $OPTIONS "$DATA/GOV2_3" \
    > "$DIR/GOV2_100K-1M_comparison.csv" \
    2>> "$DIR/comparison.err"

echo "=== GOV2 1M-10M ==="
"$EXE" $OPTIONS "$DATA/GOV2_2" \
    > "$DIR/GOV2_1M-10M_comparison.csv" \
    2>> "$DIR/comparison.err"

echo "=== GOV2 10M+ ==="
"$EXE" $OPTIONS "$DATA/GOV2_1" \
    > "$DIR/GOV2_10M-_comparison.csv" \
    2>> "$DIR/comparison.err"

# Averages (single file per range, so average = the file itself)
NUM_COLUMNS=$(head -1 "$DIR/GOV2_10M-_comparison.csv" | tr -cd ',' | wc -c)
NUM_COLUMNS=$((NUM_COLUMNS + 1))

for range in "10M-:GOV2_10M-_comparison.csv" "1M-10M:GOV2_1M-10M_comparison.csv" "100K-1M:GOV2_100K-1M_comparison.csv"; do
    label="${range%%:*}"
    file="$DIR/${range##*:}"
    printf "'GOV2_AVG_%s'" "$label" >> "$DIR/GOV2_averages.csv"
    for i in $(seq 2 "$NUM_COLUMNS"); do
        printf ","
        awk -v col="$i" -F ',' 'NR>1 { total += $col; n++ } END { if(n>0) printf "%.3e", total/n }' "$file"
    done >> "$DIR/GOV2_averages.csv"
    printf "\n" >> "$DIR/GOV2_averages.csv"
done

# Build comparison.csv
cat "$DIR/DNA_5GRAM_URL_comparison.csv" > "$DIR/comparison.csv"
cat "$DIR/GOV2_averages.csv"           >> "$DIR/comparison.csv"
sed -n '2p' < "$DIR/GOV2_10M-_comparison.csv" >> "$DIR/comparison.csv"

echo "=== sux Elias-Fano (Rust) ==="
"$SUX" \
    "$DATA/DNA_1"   "$DATA/DNA_2"   "$DATA/DNA_3" \
    "$DATA/5GRAM_1" "$DATA/5GRAM_2" "$DATA/5GRAM_3" \
    "$DATA/URL_1"   "$DATA/URL_2"   "$DATA/URL_3" \
    "$DATA/GOV2_3"  "$DATA/GOV2_2"  "$DATA/GOV2_1" \
    > "$DIR/sux_ef_comparison.csv" \
    2>> "$DIR/comparison.err"

echo "=== Done. Results in $DIR/ ==="
