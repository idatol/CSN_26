#!/usr/bin/env bash
# Runs ./closeness on every language and collects Table 2.
#
# Usage:  ./run_all.sh                    full run (a few hours, mostly Czech)
#         QUICK=1 ./run_all.sh            fast test run (T=21, approximate everywhere)
#         LANGS="Hungarian Basque" ./run_all.sh   only some languages
#         SEED=7 ./run_all.sh             different random seed (default 1)
#
# Output (in results/):
#   <Language>.tsv    one result line per language (re-run a language by deleting its file)
#   <Language>.log    stderr of that run
#   table2.tsv        all languages, plain text
#   table2_rows.tex   rows for Table 2 of the report
# Languages whose <Language>.tsv already exists are skipped.
set -euo pipefail
cd "$(dirname "$0")"

[ -x ./closeness ] || make
mkdir -p results

SEED=${SEED:-1}

# language  T  mfrac   (mfrac=1 -> exact with early stopping; <1 -> sampled null graphs)
# Small networks are done exactly; the larger ones with a fraction of the vertices.
# The real network's C is always exact.
CONFIG="Hungarian 100 1.0
Greek 100 1.0
Basque 100 1.0
Chinese 100 1.0
Arabic 100 1.0
Catalan 100 1.0
English 100 0.1
Italian 100 0.1
Turkish 100 0.1
Czech 50 0.01"

if [ "${QUICK:-0}" = 1 ]; then
  CONFIG=$(echo "$CONFIG" | awk '{ f = ($1 == "Czech") ? 0.01 : 0.1; print $1, 21, f }')
fi

LANGS=${LANGS:-$(echo "$CONFIG" | awk '{print $1}')}

for lang in $LANGS; do
  out=results/$lang.tsv
  if [ -s "$out" ]; then echo "skip $lang (results/$lang.tsv exists)"; continue; fi
  read -r T M <<<"$(echo "$CONFIG" | awk -v l="$lang" '$1 == l {print $2, $3}')"
  echo "== $lang: T=$T mfrac=$M seed=$SEED  ($(date +%H:%M:%S))"
  # Q=0 -> max(10, ceil(ln E)); the program prints header + one line on stdout
  ./closeness "data/${lang}_dependency_network.txt" "$T" 0 "$M" "$SEED" \
    2> "results/$lang.log" | tail -n 1 > "$out.tmp"
  # add the parameters used as extra columns
  awk -v t="$T" -v m="$M" -v s="$SEED" '{print $0 "\t" t "\t" m "\t" s}' "$out.tmp" > "$out"
  rm "$out.tmp"
  cat "$out"
done

# ---- collect ----
printf "language\tN\tE\t<k>\tdelta\tC\tp_binomial\tp_switching\tT\tmfrac\tseed\n" > results/table2.tsv
for lang in Arabic Basque Catalan Chinese Czech English Greek Hungarian Italian Turkish; do
  [ -s results/$lang.tsv ] && cat results/$lang.tsv >> results/table2.tsv
done

# LaTeX rows: Language & C & p_bin & p_sw & mode
awk -F'\t' 'NR > 1 {
  mode = ($10 == 1.0 || $10 == 1) ? sprintf("exact, $T=%d$", $9) \
                                  : sprintf("approx.\\ ($M/N=%g\\%%$, $T=%d$)", $10 * 100, $9)
  printf "%s & %.4f & %.2f & %.2f & %s \\\\\n", $1, $6, $7, $8, mode
}' results/table2.tsv > results/table2_rows.tex

echo; cat results/table2.tsv
