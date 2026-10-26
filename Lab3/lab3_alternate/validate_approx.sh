#!/usr/bin/env bash
# Justifies the approximate mode: runs one small language exactly and with several
# fractions M/N of sampled vertices, so p-values can be compared.
#
# Usage:  ./validate_approx.sh [Language=Hungarian] [T=100]
# Output: results/validation_<Language>.tsv  (columns: language N E <k> delta C p_bin p_sw T mfrac)
set -euo pipefail
cd "$(dirname "$0")"

LANG_NAME=${1:-Hungarian}
T=${2:-100}
[ -x ./closeness ] || make
mkdir -p results
out=results/validation_$LANG_NAME.tsv
printf "language\tN\tE\t<k>\tdelta\tC\tp_binomial\tp_switching\tT\tmfrac\n" > "$out"

for M in 1.0 0.5 0.2 0.1 0.05 0.01; do
  echo "== $LANG_NAME mfrac=$M T=$T"
  ./closeness "data/${LANG_NAME}_dependency_network.txt" "$T" 0 "$M" 1 2>/dev/null \
    | tail -n 1 | awk -v t="$T" -v m="$M" '{print $0 "\t" t "\t" m}' | tee -a "$out"
done
