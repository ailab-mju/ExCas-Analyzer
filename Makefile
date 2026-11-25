# ExCas-Analyzer build.
#   make            -> ./kmer_align      (native CLI; BAM output when htslib is found)
#   make HTS=0      -> ./kmer_align      (SAM only)
#   make wasm       -> build/wasm/kmer_align.js + kmer_align.wasm  (needs emsdk, tested with 3.1.18)
#   make clean
#
# g++-7 was the original compiler; any g++/clang++ produces byte identical
# output (see docs/FIXES.md), so fall back to whatever is installed.
# Override with:  make CPP=clang++   (CPP is a built-in make variable, so it
# has to be assigned with '=' rather than '?=' for this default to take.)
CPP = $(shell command -v g++-7 2>/dev/null || command -v g++ 2>/dev/null)
# struct data collides with std::data under -std=gnu++17 (the default from
# g++ 6 onwards), which is why the language level is pinned here.
STD = -std=gnu++14

# ---------------------------------------------------------------- htslib
# With htslib, kmer_align writes coordinate-sorted BAM + .bai by itself and
# needs no samtools at run time.  Without it the same code still sorts, but
# emits SAM.  Point HTS_CFLAGS / HTS_LDFLAGS at a non-standard prefix:
#   make HTS_CFLAGS=-I$$HOME/tools/include \
#        HTS_LDFLAGS='-L$$HOME/tools/lib -Wl,-rpath,$$HOME/tools/lib'
# (an unpacked .deb has no libhtscodecs.so symlink, so add -lhtscodecs there:
#        make HTS_LIBS='-lhts -lhtscodecs')
# or turn it off entirely with:  make HTS=0
HTS ?= auto
HASH := \#
HTS_CFLAGS ?= $(shell pkg-config --cflags htslib 2>/dev/null)
HTS_LDFLAGS ?=
HTS_LIBS ?= -lhts
ifeq ($(HTS),0)
  HTS_FOUND = 0
else
  HTS_FOUND := $(shell printf '$(HASH)include <htslib/hts.h>\nint main(){return hts_version()==0;}\n' > .hts_probe.cc; \
     $(CPP) $(STD) $(HTS_CFLAGS) .hts_probe.cc -o /dev/null $(HTS_LDFLAGS) $(HTS_LIBS) >/dev/null 2>&1 \
     && echo 1 || echo 0; rm -f .hts_probe.cc)
endif
ifeq ($(HTS_FOUND),1)
  HTS_DEF = -DUSE_HTSLIB
  HTS_LINK = $(HTS_LDFLAGS) $(HTS_LIBS)
else
  HTS_DEF =
  HTS_LINK =
endif

FLAGS = -O3 -g -pthread $(STD) $(HTS_CFLAGS) $(HTS_DEF)

# ---------------------------------------------------------------- sources
CORE_SRC = $(wildcard src/core/*.cpp)
HDR      = $(wildcard src/core/*.h) src/core/dkm.hpp
OBJDIR   = build/obj
CORE_OBJ = $(patsubst src/core/%.cpp,$(OBJDIR)/%.o,$(CORE_SRC))
OBJ      = $(CORE_OBJ) $(OBJDIR)/main.o

kmer_align: $(OBJ)
	@if [ "$(HTS_FOUND)" = "1" ]; then echo "-- htslib found: BAM output enabled"; \
	 else echo "-- htslib NOT found: SAM output only (install libhts-dev for BAM)"; fi
	$(CPP) $(FLAGS) -o $@ $(OBJ) $(HTS_LINK)

# Every object depends on every header: struct Args is passed by value across
# files, so a stale object after a header change corrupts the heap at run time.
# The project is small enough that rebuilding everything is the right trade.
$(OBJDIR)/%.o: src/core/%.cpp $(HDR) Makefile | $(OBJDIR)
	$(CPP) $(FLAGS) -c $< -o $@
$(OBJDIR)/main.o: src/cli/main.cpp $(HDR) Makefile | $(OBJDIR)
	$(CPP) $(FLAGS) -c $< -o $@
$(OBJDIR):
	mkdir -p $(OBJDIR)

# ---------------------------------------------------------------- wasm
# Same core sources; src/wasm/wasm_entry.cpp replaces the CLI main and exposes
# run_analysis() to JavaScript.  No htslib: the browser build writes SAM.
WASM_OUT = build/wasm
wasm: $(WASM_OUT)/kmer_align.js
$(WASM_OUT)/kmer_align.js: $(CORE_SRC) src/wasm/wasm_entry.cpp $(HDR) Makefile
	@command -v emcc >/dev/null || { echo "emcc not found; source emsdk_env.sh first" >&2; exit 1; }
	mkdir -p $(WASM_OUT)
	emcc -O3 $(STD) -Isrc/core $(CORE_SRC) src/wasm/wasm_entry.cpp -o $@ \
	  -s MODULARIZE=1 \
	  -s EXPORT_NAME=CrisprLargeDel \
	  -s ENVIRONMENT=web,worker \
	  -s INVOKE_RUN=0 \
	  -s EXIT_RUNTIME=0 \
	  -s ALLOW_MEMORY_GROWTH=1 \
	  -s INITIAL_MEMORY=64MB \
	  -s TOTAL_STACK=8MB \
	  -s FORCE_FILESYSTEM=1 \
	  -s EXPORTED_RUNTIME_METHODS='["ccall","FS"]' \
	  -s EXPORTED_FUNCTIONS='["_run_analysis","_malloc","_free"]' \
	  -s ASSERTIONS=0
	@echo "wasm -> $(WASM_OUT)/kmer_align.js + kmer_align.wasm ($$(emcc --version | head -1))"

.PHONY: wasm clean
clean:
	rm -rf build kmer_align .hts_probe.cc
