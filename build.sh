#!/usr/bin/env bash
# Single entry point for both builds. Same core sources; only the front-end differs.
#   ./build.sh native   -> build/kmer_align        (Sink = files,        --emit=all)
#   ./build.sh wasm     -> build/webtest-core.js   (Sink = postMessage,  --emit=split)
#   ./build.sh both
set -e
cd "$(dirname "$0")"
mkdir -p build

CORE="src/core/align.cpp src/core/align_wrap.cpp src/core/arg_parse.cpp \
      src/core/func.cpp src/core/post.cpp src/core/sink.cpp"
STD="-std=c++14 -O3 -ffp-contract=off -Werror=return-type"

build_native() {
  ${CXX:-g++} $STD $CORE src/cli/main.cpp -o build/kmer_align
  echo "native -> build/kmer_align   ($(${CXX:-g++} --version | head -1))"
}

build_wasm() {
  command -v em++ >/dev/null || { echo "em++ not found; source emsdk_env.sh first" >&2; exit 1; }
  em++ $STD $CORE src/wasm/bindings.cpp -o build/webtest-core.js \
    -s EXPORTED_FUNCTIONS="['_run_webtest','_run_align_worker','_run_join']" \
    -s EXPORTED_RUNTIME_METHODS="['ccall','cwrap','FS']" \
    -s FORCE_FILESYSTEM=1 -lworkerfs.js -lnodefs.js \
    -s ALLOW_MEMORY_GROWTH=1 -s INITIAL_MEMORY=64MB -s MAXIMUM_MEMORY=2GB \
    -s TOTAL_STACK=4MB \
    -s ASSERTIONS=0 -s WASM=1
  echo "wasm   -> build/webtest-core.js + build/webtest-core.wasm   ($(em++ --version | head -1))"
}

case "${1:-both}" in
  native) build_native ;;
  wasm)   build_wasm ;;
  both)   build_native; build_wasm ;;
  *) echo "usage: $0 [native|wasm|both]" >&2; exit 2 ;;
esac
