#!/usr/bin/env bash
# Build and check the analyzer against the reference output in tests/golden/.
# The golden files are the SAM output of this tree on the bundled TRAC example
# (2026-11-25, --out_fmt sam); every change is gated on reproducing them.
#
#   tests/run_golden.sh            build + run + diff
#   tests/run_golden.sh --update   regenerate tests/golden/ from the current build
set -e
cd "$(dirname "$0")/.."
D=test_data
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

make >/dev/null
./kmer_align $D/TRAC.fa $D/TRAC_SIM_ctr_1.fastq $D/TRAC_SIM_ctr_2.fastq \
             $D/TRAC_SIM_mut_1.fastq $D/TRAC_SIM_mut_2.fastq \
             2952 2971 17 TRAC SIM 11 "$OUT" --out_fmt sam --report 0 >/dev/null
# the two full alignment streams are ~40 MB each and are not kept in tests/golden/;
# every call is derived from the ten remaining streams, which are compared.
rm -f "$OUT"/TRAC_SIM_ctr.sam "$OUT"/TRAC_SIM_mut.sam

if [ "$1" = "--update" ]; then
  rm -f tests/golden/*.sam
  cp "$OUT"/*.sam tests/golden/
  echo "tests/golden/ updated:"; ls tests/golden
  exit 0
fi
echo "native vs golden:"
tests/compare.sh tests/golden "$OUT" golden native
