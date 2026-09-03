/* Drives the two-worker entry points outside a browser.
   Each invocation is one wasm instance, exactly as one Web Worker would be.

     node tests/run_wasm2.js align <datadir> <ref> <f1> <f2> <st> <ed> <cle> <gene> <sample> <K> <emit> <role> <out.json>
     node tests/run_wasm2.js join  <ctr.json> <mut.json> <outdir> <gene> <sample>          */
const fs = require('fs'), path = require('path'), vm = require('vm');
const mode = process.argv[2];

function loadModule(onReady) {
  global.require = require;
  global.__dirname = path.resolve(__dirname);
  global.Module = {
    locateFile: function (f) { return path.join(__dirname, '..', 'build', f); },
    print: function () {}, printErr: function (s) { process.stderr.write(s + '\n'); },
    onRuntimeInitialized: onReady
  };
  vm.runInThisContext(fs.readFileSync(path.join(__dirname, '..', 'build', 'webtest-core.js'), 'utf8'),
                      { filename: 'webtest-core.js' });
}

if (mode === 'align') {
  const [datadir, ref, f1, f2, st, ed, cle, gene, sample, K, emit, role, out] = process.argv.slice(3);
  const res = { role: role, counters: null, split: [], streams: {} };
  global.postMessage = function (m) {
    if (!Array.isArray(m)) return;
    if (m[1] === 'counters') { res.counters = { align_rate: m[3], t_align: m[4], tg_align_0: m[5], tg_ins: m[6], tg_del: m[7], only: m[8] }; return; }
    if (m[1] === 'split')    { res.split.push({ l1: m[3], l2: m[4], x: m[5], y: m[6] }); return; }
    if (m[2] === 'p_read')   { (res.streams[m[3]] = res.streams[m[3]] || []).push(m[1]); return; }
  };
  loadModule(function () {
    const FS = Module.FS;
    FS.mkdir('/data'); FS.mount(NODEFS, { root: datadir }, '/data');
    Module.ccall('run_align_worker', null, Array(11).fill('string'),
                 [ref, f1, f2, st, ed, cle, gene, sample, K, emit, role]);
    fs.writeFileSync(out, JSON.stringify(res));
    console.error(`${role}: ${res.split.length} split pairs`);
  });
} else if (mode === 'join') {
  const [ctrFile, mutFile, outdir, gene, sample] = process.argv.slice(3);
  const ctr = JSON.parse(fs.readFileSync(ctrFile)), mut = JSON.parse(fs.readFileSync(mutFile));
  const pts = r => r.split.map(s => s.x + ',' + s.y).join(' ');
  const join = {};
  global.postMessage = function (m) {
    if (Array.isArray(m) && m[1] === 'join') join[m[2]] = { rows: m[3].split('\n').filter(Boolean), only: m[4] };
  };
  loadModule(function () {
    Module.ccall('run_join', null, ['string', 'string'], [pts(ctr), pts(mut)]);
    fs.mkdirSync(outdir, { recursive: true });
    for (const [role, r] of [['ctr', ctr], ['mut', mut]]) {
      const buckets = { only: [], common: [] };
      r.split.forEach(function (s, i) {
        const [ind, color] = join[role].rows[i].split('\t');
        const b = ind === '1' ? buckets.common : buckets.only;
        b.push(s.l1 + '\t' + color, s.l2 + '\t' + color);
      });
      for (const k of ['only', 'common']) {
        const lines = buckets[k].slice().sort((a, b) =>
          parseInt(a.split('\t')[3], 10) - parseInt(b.split('\t')[3], 10));
        fs.writeFileSync(path.join(outdir, `${gene}_${sample}_${role}_${k}.sam`),
                         lines.length ? lines.join('\n') + '\n' : '');
      }
    }
    const c = ctr.counters, m = mut.counters;
    console.log('ctr  align=' + (c.align_rate * 100).toFixed(2) + '%  largedel=' +
                (c.tg_align_0 / c.t_align * 100).toFixed(2) + '%  smallindel=' +
                ((c.tg_ins + c.tg_del) / (c.tg_align_0 + c.only) * 100).toFixed(2) + '%');
    console.log('mut  align=' + (m.align_rate * 100).toFixed(2) + '%  largedel=' +
                (join.mut.only / (m.tg_align_0 + join.mut.only) * 100).toFixed(2) + '%  smallindel=' +
                ((m.tg_ins + m.tg_del) / (m.tg_align_0 + join.mut.only) * 100).toFixed(2) + '%');
  });
} else { console.error('usage: align|join'); process.exit(2); }
