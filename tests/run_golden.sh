#!/usr/bin/env bash
# Build and check the analyzer against the reference output in tests/golden/.
# Those files come from the kmer_align binary that produced result/HPRT1_TEST_11
# (verified byte-identical), so this is the regression gate for every change.
#
#   tests/run_golden.sh            native only
#   tests/run_golden.sh --wasm     native + wasm (needs emsdk on PATH)
set -e
cd "$(dirname "$0")/.."
D=test_data
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

./build.sh native >/dev/null
./build/kmer_align $D/HPRT1.fa $D/test_con_1.fastq $D/test_con_2.fastq \
                   $D/test_mut_1.fastq $D/test_mut_2.fastq \
                   4652 4671 3 HPRT1 TEST 11 "$OUT" 1000 --emit=all >/dev/null
echo "native:"
tests/compare.sh tests/golden "$OUT" golden native

if [ "$1" = "--wasm" ]; then
  command -v em++ >/dev/null || { echo "em++ not found; source emsdk_env.sh first" >&2; exit 1; }
  ./build.sh wasm >/dev/null
  W=$(mktemp -d); trap 'rm -rf "$OUT" "$W"' EXIT
  node tests/run_wasm.js $D "$W" HPRT1 TEST 4652 4671 3 11 all \
       HPRT1.fa test_con_1.fastq test_con_2.fastq test_mut_1.fastq test_mut_2.fastq >/dev/null
  for n in ctr mut; do
    [ -f "$W/HPRT1_TEST_${n}_com.sam" ] && mv "$W/HPRT1_TEST_${n}_com.sam" "$W/HPRT1_TEST_${n}_common.sam"
  done
  echo "wasm:"
  tests/compare.sh tests/golden "$W" golden wasm
fi
