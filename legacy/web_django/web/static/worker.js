/* One worker == one wasm instance. The page runs two of them, one per sample,
   so the core's file-scope state stays private to each. No SharedArrayBuffer and
   no pthreads, so the page needs no cross-origin isolation (COOP/COEP).

   WASM=1 instantiates asynchronously, so anything that arrives before the runtime
   is up has to wait -- importScripts() returning is not enough. */
var ready = false;
var pending = [];

var Module = {
  onRuntimeInitialized: function () {
    ready = true;
    postMessage([0, 'worker-ready']);
    var q = pending; pending = [];
    q.forEach(handle);
  },
  onAbort: function (what) { postMessage([0, 'worker-abort', String(what)]); }
};

importScripts('webtest-core.js');

const prefix = '/data';
let mounted = false;

function handle(e) {
  try {
    if (e.cmd === 'align') {
      if (!mounted) {
        try { FS.unmount(prefix); } catch (err) {}
        try { FS.rmdir(prefix);  } catch (err) {}
        FS.mkdir(prefix);
        FS.mount(WORKERFS, { files: e.files }, prefix);   // lazy, no copy into the heap
        mounted = true;
      }
      Module.ccall('run_align_worker', null, Array(11).fill('string'),
        [e.files[0].name, e.files[1].name, e.files[2].name,
         e.st, e.ed, e.cle, e.gene, e.sample, e.K, e.emit || 'split', e.role]);

    } else if (e.cmd === 'join') {
      Module.ccall('run_join', null, ['string', 'string'], [e.ctr, e.mut]);
    }
  } catch (err) {
    postMessage([0, 'worker-abort', (err && err.message) ? err.message : String(err)]);
  }
}

onmessage = function (msg) {
  if (ready) handle(msg.data); else pending.push(msg.data);
};
