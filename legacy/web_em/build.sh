#em++ -O3 -std=c++14 post.cpp -c -o post.o
#em++ -O3 -std=c++14 func.cpp -c -o func.o
#em++ -O3 -std=c++14 arg_parse.cpp -c -o arg_parse.o
#em++ -O3 -std=c++14 align.cpp -c -o align.o
#em++ -O3 -std=c++14 align_wrap.cpp -c -o align_wrap.o
#em++ -O3 -std=c++14 emscripten_interface.cpp -c -o ems.o
#em++ -O3 -std=c++14 align.o align_wrap.o arg_parse.o func.o post.o ems.o -o webtest-core.js -s EXPORTED_FUNCTIONS="['_run_webtest']" -s EXPORTED_RUNTIME_METHODS="['ccall','cwrap']" -s WASM=0 -s TOTAL_MEMORY=1600MB --memory-init-file 0 -s ASSERTIONS=2 -s DEMANGLE_SUPPORT=1 
em++ -O3 -std=c++14 align.cpp align_wrap.cpp arg_parse.cpp func.cpp post.cpp emscripten_interface.cpp -o webtest-core.js -s EXPORTED_FUNCTIONS="['_run_webtest']" -s EXPORTED_RUNTIME_METHODS="['ccall','cwrap']" -s WASM=0 -s TOTAL_MEMORY=1600MB --memory-init-file 0 -s ASSERTIONS=2 -s DEMANGLE_SUPPORT=1 -s USE_ZLIB=1
#em++ -O3 -std=c++14 align.o align_wrap.o arg_parse.o func.o post.o ems.o -o webtest-core.js -s EXPORTED_FUNCTIONS="['_run_webtest']" -s EXPORTED_RUNTIME_METHODS="['ccall','cwrap']" -s WASM=0 -s TOTAL_MEMORY=1600MB --memory-init-file 0 -s ASSERTIONS=2 -s DEMANGLE_SUPPORT=1 
#ALLOW_MEMORY_GROWTH=1
#cp webtest-core.js ../../web/large_del_test/web/static/
#cp -f webtest-core.js ../../web/django_emscripten_test/web_test/static/
