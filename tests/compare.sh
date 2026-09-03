#!/usr/bin/env bash
# compare.sh <dirA> <dirB> [labelA] [labelB]
# Byte-compares every SAM stream in two output directories.
# A mismatch that disappears once the YC:Z cluster-colour tag is stripped is
# reported separately: that tag is a display colour whose index depends on the
# order k-means happens to number its centroids, which differs between the
# native and wasm floating-point paths. It carries no analysis result.
A=$1; B=$2; LA=${3:-A}; LB=${4:-B}
ok=0; ng=0; colour=0
for f in ctr mut ctr_only ctr_common mut_only mut_common site_ctr site_mut site_ctr_sm site_mut_sm; do
  a=$(ls "$A"/*_"$f".sam 2>/dev/null | head -1); b=$(ls "$B"/*_"$f".sam 2>/dev/null | head -1)
  [ -n "$a" ] && sort -k4n "$a" > /tmp/cmp_a || : > /tmp/cmp_a
  [ -n "$b" ] && sort -k4n "$b" > /tmp/cmp_b || : > /tmp/cmp_b
  na=$(wc -l < /tmp/cmp_a); nb=$(wc -l < /tmp/cmp_b)
  [ "$na" = "0" ] && [ "$nb" = "0" ] && continue
  if diff -q /tmp/cmp_a /tmp/cmp_b >/dev/null; then
    printf "  %-14s match             %5s lines\n" "$f" "$na"; ok=$((ok+1))
  elif diff -q <(sed 's/\tYC:Z:[0-9,]*//' /tmp/cmp_a) <(sed 's/\tYC:Z:[0-9,]*//' /tmp/cmp_b) >/dev/null; then
    printf "  %-14s match but YC:Z    %5s lines  (cluster colour index only)\n" "$f" "$na"; colour=$((colour+1))
  else
    printf "  %-14s MISMATCH          %s=%s %s=%s\n" "$f" "$LA" "$na" "$LB" "$nb"; ng=$((ng+1))
  fi
done
echo "  -> $ok identical, $colour colour-tag-only, $ng real mismatch"
[ "$ng" = "0" ]
