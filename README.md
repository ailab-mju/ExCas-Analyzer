# ExCas-Analyzer

Quantifies CRISPR-induced **large deletions** from paired-end amplicon sequencing by
comparing a Cas9-treated sample against an untreated control. Deletion breakpoints are
clustered, and clusters that also occur in the control are removed as false positives.

One source tree builds two front-ends over the same core:

| | | |
|---|---|---|
| **CLI** | `build/kmer_align` | writes SAM files |
| **Browser** | `build/webtest-core.{js,wasm}` | two Web Workers, results never leave the machine |

Both are verified to produce **byte-identical** output (see [Verifying](#verifying)).

## Quick start

```bash
./build.sh native                 # -> build/kmer_align
tests/run_golden.sh               # regression check against tests/golden/

./build/kmer_align \
    test_data/HPRT1.fa \
    test_data/test_con_1.fastq test_data/test_con_2.fastq \
    test_data/test_mut_1.fastq test_data/test_mut_2.fastq \
    4652 4671 3 HPRT1 TEST 11 /tmp/out 1000 --emit=all
```

For the browser build you need [emsdk](https://emscripten.org) on `PATH`
(verified with 3.1.18):

```bash
source /path/to/emsdk/emsdk_env.sh
./build.sh wasm                   # -> build/webtest-core.{js,wasm}
cd web && ./make_examples.sh && python3 manage.py runserver 127.0.0.1:8899
```

## CLI

```
kmer_align ref.fa ctr_1.fastq ctr_2.fastq mut_1.fastq mut_2.fastq \
           tg_start tg_end cleavage_offset gene sample K outdir [dist_cut] \
           [--emit=all|none|split[,full][,site]] [--pre-chunk] [--datadir=PATH]
```

Positional arguments are unchanged from the original tool. `tg_start`/`tg_end` are
1-based inclusive protospacer coordinates in the reference; `cleavage_offset` is the
cut site measured from `tg_start` (3 for a minus-strand SpCas9 protospacer). `K` is the
k-mer length used to index the reference; the published setting is 11.

| flag | meaning |
|---|---|
| `--emit` | which SAM streams to write. Default `all`. The browser defaults to `split`, which is ~0.2% of the output volume — pass `--emit=split` here to compare the two builds. |
| `--pre-chunk` | enable the `pre_chunk()` fast path. **Off by default**: it changes results, and the published numbers were produced without it. |
| `--datadir` | prefix for input paths (used by the wasm build to reach its WORKERFS mount). |

Ten SAM streams are produced: `ctr`/`mut` (every aligned pair), `site_*` (on-target),
`site_*_sm` (on-target with an indel), and `*_only`/`*_common` (split reads carrying a
large deletion, split by whether the cluster also appears in the control).

## Browser

`web/` is a Django project that serves a static page; nothing is uploaded. Input follows
[Cas-Analyzer](http://www.rgenome.net/cas-analyzer/): paste the reference and the target
sequence, and the page derives `tg_start`, `tg_end` and the cleavage offset itself.

Two Web Workers run the two samples, each with its own wasm instance. The clustering
step needs both samples, so it runs once at the end — in the C++ core, inside one of the
workers, not reimplemented in JavaScript. Because there is no SharedArrayBuffer and no
pthreads, the page needs **no cross-origin isolation (COOP/COEP)**.

Three example sets ship with the page (`web/README.md` has the expected numbers).

## Layout

```
build.sh              single entry point: native | wasm | both
src/core/             platform-free. No emscripten, no FILE*, no postMessage
  sink.h              the only three seams: Sink, Progress, and the emit mask
  align.cpp           k-mer seeding, LIS chaining, CIGAR
  align_wrap.cpp      per-sample driver (align_sample) + the CLI's align_wrap
  post.cpp            breakpoint clustering and the control/treated join
  func.cpp            FASTA/FASTQ, SAM formatting, pair logic
src/cli/main.cpp      Sink = files, Progress = stdout
src/wasm/bindings.cpp run_webtest (one instance) / run_align_worker + run_join (two)
tests/
  golden/             reference output, byte-identical to legacy/result/
  run_golden.sh       build + run + diff, native and optionally wasm
  compare.sh          diff two output directories stream by stream
  run_wasm.js         run the wasm build under node (single instance)
  run_wasm2.js        run the two-worker entry points under node
web/                  Django + the two-worker page
test_data/            2,500 read pairs, HPRT1/HEK293T
legacy/               superseded trees, kept for provenance (not required to build)
```

## Verifying

`tests/golden/` is the output of the `kmer_align` binary committed with the original
release (kept at `legacy/kmer_align`), confirmed byte-identical to `legacy/result/`.
Every change in this tree is gated on reproducing it.

```bash
tests/run_golden.sh          # native
tests/run_golden.sh --wasm   # native + wasm
```

Checked during the refactor, on `test_data` (2,500 pairs) and on real data
(HEK293T HPRT1, 1.7M pairs, 977 MB; HEK293T ST3GAL4, 2.35M pairs):

- native and wasm agree byte for byte, including the `YC:Z` cluster colour
- g++-7.5, g++-9.4 and emsdk clang all agree
- `-O1 -fsanitize=address,undefined` runs clean and produces the same output
