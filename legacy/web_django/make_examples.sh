#!/usr/bin/env bash
# Regenerate the ST3GAL4 and WISP3 example sets for the demo page.
# HPRT1 is the repo's own test_data and is committed; these two are 10,000-read-pair
# head subsets of the lab data, so they are generated rather than stored.
#
# Coordinates come from the original run scripts:
#   ST3GAL4  4344-4363   (revision/run_all_sample.sh)
#   WISP3    5054-5073   (bin/, clean_bin/ scripts; H9 uses 3838-3857 instead)
set -e
cd "$(dirname "$0")"
LAB=${LAB:-/data/project/noalcohol/CRISP/CRISP}
N=${N:-10000}
OUT=web/static/testdata

mk() {                      # mk <gene> <control prefix> <treated prefix>
  local g=$1 c=$2 m=$3 d=$OUT/$1
  mkdir -p "$d"
  cp "$LAB/ref/$g.fa" "$d/"
  head -$((N*4)) "$LAB/data/all_data/${c}_1.fastq" > "$d/${g}_con_1.fastq"
  head -$((N*4)) "$LAB/data/all_data/${c}_2.fastq" > "$d/${g}_con_2.fastq"
  head -$((N*4)) "$LAB/data/all_data/${m}_1.fastq" > "$d/${g}_mut_1.fastq"
  head -$((N*4)) "$LAB/data/all_data/${m}_2.fastq" > "$d/${g}_mut_2.fastq"
  (cd "$d" && zip -q -r "../${g}_example.zip" .)
  echo "$g: $N pairs -> $d"
}

mk ST3GAL4 HEK293T_ST3GAL4_con_1st HEK293T_ST3GAL4_1st
mk WISP3   HEK293T_WISP3_con_1st   HEK293T_WISP3_1st
(cd $OUT/HPRT1 && zip -q -r ../HPRT1_example.zip .)
echo "HPRT1: zip refreshed"
