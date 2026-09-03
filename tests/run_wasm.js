/* Runs build/webtest-core.js outside a browser by imitating the Web Worker global
   scope (that is how worker.js loads it: importScripts + global Module/FS/WORKERFS).

   Inputs are mounted with NODEFS, not copied into the wasm heap with FS.writeFile --
   that mirrors what WORKERFS does in the browser and keeps a 1 GB dataset out of
   the 32-bit heap.

   usage: node tests/run_wasm.js <datadir> <outdir> <gene> <sample> <st> <ed> <cle> <K> <emit> \
                                 <ref> <ctr1> <ctr2> <mut1> <mut2>                          */
const fs = require('fs'), path = require('path'), vm = require('vm');

const a = process.argv.slice(2);
if (a.length < 14) { console.error('usage: see header'); process.exit(2); }
const [datadir, outdir, gene, sample, st, ed, cle, K, emit] = a;
const files = a.slice(9, 14);

const streams = {};
let results = [], ended = false, lastPct = {};
global.postMessage = function (m) {
  if (!Array.isArray(m)) return;
  if (m[2] === 'p_read') { (streams[m[3]] = streams[m[3]] || []).push(m[1]); return; }
  if (m[2] === 'result') { results.push(m.slice(1)); return; }
  if (typeof m[1] === 'number') {
    if (lastPct[m[2]] !== m[1]) { lastPct[m[2]] = m[1]; process.stderr.write(`\r${m[2]} ${m[1]}%   `); }
    return;
  }
  if (m[1] === 'end') { ended = true; return; }
};
global.require = require;
global.__dirname = path.resolve(__dirname);

const t0 = Date.now();
global.Module = {
  locateFile: function (f) { return path.join(__dirname, '..', 'build', f); },
  print: function () {}, printErr: function (s) { process.stderr.write(s + '\n'); },
  onRuntimeInitialized: function () {
    const FS = Module.FS;
    FS.mkdir('/data');
    FS.mount(NODEFS, { root: datadir }, '/data');   // lazy, no copy into the heap
    Module.ccall('run_webtest', null, Array(13).fill('string'),
      [files[0], files[1], files[2], files[3], files[4], st, ed, cle, gene, sample, K, 'out', emit]);

    fs.mkdirSync(outdir, { recursive: true });
    for (const k of Object.keys(streams)) {
      const lines = streams[k].slice().sort((x, y) =>
        parseInt(x.split('\t')[3], 10) - parseInt(y.split('\t')[3], 10));
      fs.writeFileSync(path.join(outdir, `${gene}_${sample}_${k}.sam`), lines.join('\n') + '\n');
    }
    process.stderr.write('\n');
    console.log('elapsed  : ' + ((Date.now() - t0) / 1000).toFixed(2) + ' s');
    console.log('peak heap: ' + (Module.HEAPU8.length / 1048576).toFixed(0) + ' MB');
    console.log('rss      : ' + (process.memoryUsage().rss / 1048576).toFixed(0) + ' MB');
    console.log('streams  : ' + Object.keys(streams).map(k => `${k}=${streams[k].length}`).join(' '));
    console.log('results  : ' + JSON.stringify(results));
    console.log('end msg  : ' + ended);
    if (!ended) process.exit(1);
  }
};
vm.runInThisContext(fs.readFileSync(path.join(__dirname, '..', 'build', 'webtest-core.js'), 'utf8'),
                    { filename: 'webtest-core.js' });
