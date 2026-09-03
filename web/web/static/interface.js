'use strict';
/* ExCas-Analyzer front-end.
   Input format follows Cas-Analyzer (RGEN Tools): the user supplies a reference sequence,
   a target sequence and a PAM; this script derives the coordinates (targetst/targeted)
   and the cleavage offset (cle) and hands them to the existing worker/wasm path.
   The worker message contract is unchanged:
   [ref, ctr1, ctr2, mut1, mut2, st, ed, cle, gene, sample, K, outdir] */

var $ = function (id) { return document.getElementById(id); };

/* ---------- sequence helpers ---------- */
var IUPAC = { A:'A', C:'C', G:'G', T:'T',
  R:'AG', Y:'CT', S:'CG', W:'AT', K:'GT', M:'AC',
  B:'CGT', D:'AGT', H:'ACT', V:'ACG', N:'ACGT' };

function cleanSeq(text) {
  return text.split(/\r?\n/)
             .filter(function (l) { return l.charAt(0) !== '>' && l.charAt(0) !== ';'; })
             .join('').replace(/\s+/g, '').toUpperCase();
}
function revcomp(s) {
  var m = { A:'T', C:'G', G:'C', T:'A', N:'N' }, o = '';
  for (var i = s.length - 1; i >= 0; i--) o += (m[s[i]] || 'N');
  return o;
}
function matchIUPAC(seq, pattern) {
  if (seq.length !== pattern.length) return false;
  for (var i = 0; i < pattern.length; i++) {
    var allowed = IUPAC[pattern[i]];
    if (!allowed || allowed.indexOf(seq[i]) < 0) return false;
  }
  return true;
}
function findAll(hay, needle) {
  var out = [], i = hay.indexOf(needle);
  while (i >= 0) { out.push(i); i = hay.indexOf(needle, i + 1); }
  return out;
}

/* ---------- target location and cleavage site ---------- */
/* returns {st, ed, strand, pam, cut, cle} with st/ed/cut 1-based */
function locateTarget(ref, target, pam, side, cutOff) {
  var hits = [];
  findAll(ref, target).forEach(function (i) { hits.push({ i: i, strand: '+' }); });
  var rc = revcomp(target);
  if (rc !== target) findAll(ref, rc).forEach(function (i) { hits.push({ i: i, strand: '-' }); });
  if (!hits.length) throw new Error('Target sequence not found in the reference (searched both strands).');

  var L = target.length, P = pam ? pam.length : 0, ok = [];
  hits.forEach(function (h) {
    var s = h.i + 1, e = h.i + L, pamSeq = null, valid = !pam;
    if (pam) {
      // for a protospacer on the minus strand the PAM appears reverse-complemented on the other side
      var pamAtRight = (side === 3) === (h.strand === '+');
      if (pamAtRight) {
        pamSeq = ref.substr(e, P);
        if (pamSeq.length === P) valid = matchIUPAC(h.strand === '+' ? pamSeq : revcomp(pamSeq), pam);
      } else {
        pamSeq = ref.substr(s - 1 - P, P);
        if (pamSeq.length === P) valid = matchIUPAC(h.strand === '+' ? pamSeq : revcomp(pamSeq), pam);
      }
    }
    if (valid) ok.push({ s: s, e: e, strand: h.strand, pam: pamSeq });
  });

  if (pam && !ok.length)
    throw new Error('The target sequence was found, but no adjacent ' + pam +
      ' PAM matches. Check the PAM selection, or choose "No PAM" and set the cleavage offset manually.');
  if (ok.length > 1)
    throw new Error('The target sequence occurs at ' + ok.length +
      ' positions in the reference. Narrow the reference or use a longer target.');

  var h2 = ok.length ? ok[0] : { s: hits[0].i + 1, e: hits[0].i + L, strand: hits[0].strand, pam: null };
  var cut = null;
  if (cutOff !== null && cutOff !== undefined && side !== 0) {
    if (side === 3) cut = (h2.strand === '+') ? h2.e - cutOff : h2.s + cutOff;
    else            cut = (h2.strand === '+') ? h2.s + cutOff : h2.e - cutOff;
  }
  return { st: h2.s, ed: h2.e, strand: h2.strand, pam: h2.pam,
           cut: cut, cle: cut === null ? null : cut - h2.s };
}

/* ---------- view ---------- */
function showErr(msg) {
  var el = $('err');
  if (!msg) { el.className = 'alert alert-danger hidden'; el.textContent = ''; return; }
  el.className = 'alert alert-danger'; el.textContent = msg;
  el.scrollIntoView({ block: 'nearest' });
}
function esc(s) { return s.replace(/&/g, '&amp;').replace(/</g, '&lt;'); }

function renderSummary(ref, loc, K) {
  var pad = 40;
  var lo = Math.max(0, loc.st - 1 - pad), hi = Math.min(ref.length, loc.ed + pad);
  var out = '';
  for (var i = lo; i < hi; i++) {
    var p = i + 1, ch = esc(ref[i]), cls = '';
    if (p >= loc.st && p <= loc.ed) cls = 'tgt';
    if (loc.pam) {
      var pamRight = (loc.strand === '+') === (loc.pamSide === 3);
      var ps = pamRight ? loc.ed + 1 : loc.st - loc.pam.length;
      if (p >= ps && p < ps + loc.pam.length) cls = 'pam';
    }
    if (loc.cut !== null && p === loc.cut) out += '<span class="cut">&#9660;</span>';
    out += cls ? '<span class="' + cls + '">' + ch + '</span>' : ch;
  }
  $('seqctx').innerHTML = out;
  $('s-reflen').textContent = ref.length.toLocaleString() + ' bp';
  $('s-pos').textContent = loc.st.toLocaleString() + '–' + loc.ed.toLocaleString();
  $('s-strand').textContent = loc.strand === '+' ? '+ (forward)' : '− (reverse)';
  $('s-pam').textContent = loc.pam || '–';
  $('s-cut').textContent = loc.cut === null ? '–'
      : loc.cut.toLocaleString() + ' (offset ' + loc.cle + ')';
  $('s-k').textContent = K;
  $('summary-box').className = '';
}

function setBar(id, pct, done) {
  var b = $(id);
  if (pct > 100) pct = 100;
  b.style.width = pct + '%';
  b.textContent = pct + '%';
  b.className = 'progress-bar' + (done ? ' done' : '');
}

/* ---------- run ---------- */
/* Two workers, one per sample. They are fully independent through the alignment;
   the clustering that decides "present in the control too" needs both, so it runs
   once at the end -- in the C++ core, inside one of the workers, rather than being
   reimplemented in JavaScript. Only the split-read breakpoints cross the boundary
   (a few thousand records even for a million read pairs). */
var workers = [];

function enableButtons() {
  $('submit').disabled = false;
  Array.prototype.forEach.call(document.querySelectorAll('button.ex'), function (b) { b.disabled = false; });
}
function stopWorkers() {
  workers.forEach(function (w) { w.terminate(); });
  workers = [];
}
function pct(v) { return (parseFloat(v) * 100).toFixed(2) + '%'; }

function run(presetFiles) {
  showErr('');
  var files;
  try {
    var ref = cleanSeq($('fullseq').value);
    if (!ref) throw new Error('Enter a reference sequence.');
    var bad = ref.replace(/[ACGTN]/g, '');
    if (bad) throw new Error('The reference contains non-DNA characters: ' + bad.slice(0, 10));

    var target = cleanSeq($('rgenseq').value);
    if (!target) throw new Error('Enter a target DNA sequence.');
    var opt = $('pams').selectedOptions[0];
    var pam = opt.value, side = parseInt(opt.dataset.side, 10);
    var cutOff = opt.dataset.cut === '' ? null : parseInt(opt.dataset.cut, 10);

    var manual = $('cle').value.trim();
    var loc = locateTarget(ref, target, pam, side, cutOff);
    loc.pamSide = side;
    if (manual !== '') { loc.cle = parseInt(manual, 10); loc.cut = loc.st + loc.cle; }
    if (loc.cle === null || isNaN(loc.cle))
      throw new Error('The cleavage site could not be derived. Enter a cleavage offset manually.');

    var K = parseInt($('K').value, 10);
    if (!(K >= 7 && K <= 16)) throw new Error('K must be between 7 and 16.');
    var gene = $('genename').value.trim() || 'gene';
    var sample = $('samplename').value.trim() || 'sample';

    var refFile = new File(['>' + gene + '\n' + ref + '\n'], gene + '.fa', { type: 'text/plain' });
    if (presetFiles) {
      files = presetFiles;
    } else {
      var ids = ['file2', 'file3', 'file4', 'file5'], fs = [];
      for (var i = 0; i < ids.length; i++) {
        var f = $(ids[i]).files[0];
        if (!f) throw new Error('Select all four FASTQ files (control Read 1/2, treated Read 1/2).');
        if (/\.gz$/i.test(f.name)) throw new Error('Gzipped files are not supported yet: ' + f.name);
        fs.push(f);
      }
      files = fs;
    }

    renderSummary(ref, loc, K);
    start(refFile, files, loc, K, gene, sample);
  } catch (e) {
    showErr(e.message);
  }
}

function start(refFile, files, loc, K, gene, sample) {
  $('progress-box').className = '';
  $('result-box').className = '';
  $('download-box').className = 'hidden';
  setBar('ctr-bar', 0, false); setBar('mut-bar', 0, false);
  ['r-ctr-align','r-ctr-ld','r-ctr-si','r-mut-align','r-mut-ld','r-mut-si']
    .forEach(function (id) { $(id).textContent = '–'; });
  $('stage').textContent = 'Starting two workers…';
  $('submit').disabled = true;
  Array.prototype.forEach.call(document.querySelectorAll('button.ex'), function (b) { b.disabled = true; });

  var state = {
    ctr: { split: [], counters: null, done: false },
    mut: { split: [], counters: null, done: false },
    streams: {}, t0: Date.now(), joined: false
  };

  stopWorkers();
  var roles = ['ctr', 'mut'];
  roles.forEach(function (role, idx) {
    var w = new Worker('/static/worker.js');
    workers.push(w);
    w.onerror = function (ev) { showErr('Worker error (' + role + '): ' + ev.message); enableButtons(); };
    w.onmessage = function (msg) { onWorkerMessage(state, role, msg.data, loc, K, gene, sample); };
    w.postMessage({
      cmd: 'align', role: role,
      files: [refFile, files[idx * 2], files[idx * 2 + 1]],
      st: String(loc.st), ed: String(loc.ed), cle: String(loc.cle),
      gene: gene, sample: sample, K: String(K), emit: 'split'
    });
  });
}

function onWorkerMessage(state, role, e, loc, K, gene, sample) {
  if (!Array.isArray(e)) return;

  if (e[1] === 'worker-ready') { $('stage').textContent = 'Aligning both samples…'; return; }
  if (e[1] === 'worker-abort') { showErr('Worker aborted (' + role + '): ' + e[2]); stopWorkers(); enableButtons(); return; }

  if (typeof e[1] === 'number' && (e[2] === 'ctr' || e[2] === 'mut')) {
    setBar(e[2] + '-bar', e[1], e[1] >= 100);
    return;
  }
  if (e[2] === 'p_read')   { (state.streams[e[3]] = state.streams[e[3]] || []).push(e[1]); return; }
  if (e[1] === 'counters') { state[e[2]].counters = { align: e[3], t_align: e[4], tg0: e[5], ins: e[6], del: e[7], only: e[8] }; return; }
  if (e[1] === 'split')    { state[e[2]].split.push({ l1: e[3], l2: e[4], x: e[5], y: e[6] }); return; }

  if (e[1] === 'aligned') {
    state[e[2]].done = true;
    if (state.ctr.done && state.mut.done && !state.joined) {
      state.joined = true;
      $('stage').textContent = 'Both samples aligned — clustering…';
      var pts = function (r) { return r.split.map(function (s) { return s.x + ',' + s.y; }).join(' '); };
      workers[0].postMessage({ cmd: 'join', ctr: pts(state.ctr), mut: pts(state.mut) });
    }
    return;
  }

  if (e[1] === 'join') { state[e[2]].join = { rows: e[3].split('\n').filter(Boolean), only: e[4] }; return; }

  if (e[1] === 'end' && state.joined && state.ctr.join && state.mut.join) {
    finish(state, gene, sample);
  }
}

function finish(state, gene, sample) {
  var sam = { ctr_only: [], ctr_com: [], mut_only: [], mut_com: [] };
  ['ctr', 'mut'].forEach(function (role) {
    state[role].split.forEach(function (s, i) {
      var parts = state[role].join.rows[i].split('\t');
      var bucket = parts[0] === '1' ? role + '_com' : role + '_only';
      sam[bucket].push(s.l1 + '\t' + parts[1], s.l2 + '\t' + parts[1]);
    });
  });

  var c = state.ctr.counters, m = state.mut.counters, mutOnly = state.mut.join.only;
  $('r-ctr-align').textContent = pct(c.align);
  $('r-ctr-ld').textContent    = pct(c.tg0 / c.t_align);
  $('r-ctr-si').textContent    = pct((c.ins + c.del) / (c.tg0 + c.only));
  $('r-mut-align').textContent = pct(m.align);
  $('r-mut-ld').textContent    = pct(mutOnly / (m.tg0 + mutOnly));
  $('r-mut-si').textContent    = pct((m.ins + m.del) / (m.tg0 + mutOnly));

  $('stage').textContent = 'Done in ' + ((Date.now() - state.t0) / 1000).toFixed(1) + ' s (two workers)';
  stopWorkers();
  enableButtons();
  publish(sam);
}

/* The server sorts the split reads by position and stores them (existing behaviour).
   If that fails, fall back to downloading them straight from the browser. */
function publish(sam) {
  var total = sam.ctr_only.length + sam.ctr_com.length + sam.mut_only.length + sam.mut_com.length;
  $('download-box').className = '';
  if (!total) $('dl-note').textContent = 'No split reads were found, so the SAM files are empty.';

  var body = new URLSearchParams();
  ['ctr_only','ctr_com','mut_only','mut_com'].forEach(function (k) { body.append(k, JSON.stringify(sam[k])); });

  fetch('/sam_upload/', { method: 'POST', body: body,
      headers: { 'Content-Type': 'application/x-www-form-urlencoded', 'X-Requested-With': 'XMLHttpRequest' } })
    .then(function (r) { if (!r.ok) throw new Error('HTTP ' + r.status); return r.json(); })
    .then(function () {
      var t = Date.now();
      ['ctr_only','ctr_com','mut_only','mut_com'].forEach(function (k) {
        var a = $('dl-' + k.replace('_', '-'));
        a.href = '/static/result_dir/' + k + '.sam?t=' + t;
        a.removeAttribute('download');
      });
      $('dl-note').textContent = total + ' split reads, sorted by position on the server.';
    })
    .catch(function (err) {
      ['ctr_only','ctr_com','mut_only','mut_com'].forEach(function (k) {
        var lines = sam[k].slice().sort(function (a, b) {
          return parseInt(a.split('\t')[3], 10) - parseInt(b.split('\t')[3], 10); });
        var a = $('dl-' + k.replace('_', '-'));
        a.href = URL.createObjectURL(new Blob([lines.join('\n') + '\n'], { type: 'text/plain' }));
        a.download = k + '.sam';
      });
      $('dl-note').textContent = total + ' split reads. Server-side save failed (' + err.message +
        '), so the files are served straight from your browser.';
    });
}

/* ---------- reference FASTA upload ---------- */
$('reffile').addEventListener('change', function () {
  var f = this.files[0];
  if (!f) return;
  f.text().then(function (t) {
    $('fullseq').value = t;
    $('reflen').textContent = '— ' + cleanSeq(t).length.toLocaleString() + ' bp loaded';
  });
});
$('fullseq').addEventListener('input', function () {
  var n = cleanSeq(this.value).length;
  $('reflen').textContent = n ? '— ' + n.toLocaleString() + ' bp' : '';
});

/* ---------- buttons ---------- */
$('submit').addEventListener('click', function () { run(null); });

/* One handler for every example row. The reference is loaded as text so the target
   coordinates are derived here, exactly as they would be for a pasted sequence. */
Array.prototype.forEach.call(document.querySelectorAll('button.ex'), function (btn) {
  btn.addEventListener('click', function () {
    var gene = btn.dataset.gene, sample = btn.dataset.sample, target = btn.dataset.target;
    var base = '/static/testdata/' + gene + '/';
    var names = [gene + '_con_1.fastq', gene + '_con_2.fastq', gene + '_mut_1.fastq', gene + '_mut_2.fastq'];
    var label = btn.textContent;
    Array.prototype.forEach.call(document.querySelectorAll('button.ex'), function (b) { b.disabled = true; });
    btn.textContent = 'Loading…';
    $('demo-note').textContent = '';
    showErr('');

    Promise.all([fetch(base + gene + '.fa').then(function (r) {
        if (!r.ok) throw new Error(gene + '.fa: HTTP ' + r.status); return r.text(); })].concat(
      names.map(function (n) {
        return fetch(base + n)
          .then(function (r) { if (!r.ok) throw new Error(n + ': HTTP ' + r.status); return r.blob(); })
          .then(function (b) { return new File([b], n); });
      })))
      .then(function (res) {
        $('fullseq').value = res[0];
        $('reflen').textContent = '— ' + cleanSeq(res[0]).length.toLocaleString() + ' bp';
        $('rgenseq').value = target;
        $('genename').value = gene;
        $('samplename').value = sample;
        $('K').value = '11';
        $('cle').value = '';
        $('pams').value = 'NGG';
        $('demo-note').textContent = gene + ' example loaded (' + sample + ')';
        run(res.slice(1));
      })
      .catch(function (e) { showErr('Failed to load the ' + gene + ' example: ' + e.message); })
      .finally(function () {
        Array.prototype.forEach.call(document.querySelectorAll('button.ex'), function (b) { b.disabled = false; });
        btn.textContent = label;
      });
  });
});
