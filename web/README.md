# Browser front-end

Django serves a static page; all analysis runs in the browser. Two Web Workers,
one per sample, each holding its own wasm instance — no SharedArrayBuffer and no
pthreads, so the page needs no cross-origin isolation (COOP/COEP).

    ../build.sh wasm                       # -> ../build/webtest-core.{js,wasm}
    cp ../build/webtest-core.js ../build/webtest-core.wasm web/static/
    ./make_examples.sh                     # ST3GAL4 / WISP3 example sets
    python3 manage.py runserver 127.0.0.1:8899

Input format follows Cas-Analyzer (rgenome.net/cas-analyzer): paste the reference
sequence and the target sequence, and the page derives the coordinates and the
cleavage offset that the CLI takes as numbers.

Reference values for the three examples (from `../build/kmer_align --emit=split`):

| example | control align / large del / small indel | treated align / large del / small indel |
|---------|------------------------------------------|------------------------------------------|
| HPRT1   2,500 pairs  | 61.44% / 1.56% / 0.00%  | 53.52% / 29.41% / 11.76% |
| ST3GAL4 10,000 pairs | 30.53% / 2.52% / 1.30%  | 64.28% / 10.28% / 24.30% |
| WISP3   10,000 pairs | 86.95% / 2.25% / 1.02%  | 92.89% /  6.90% / 25.86% |
