/* Browser front-end.

   Two entry points so the page can run the two samples in separate workers
   (separate wasm instances, so the core's file-scope state stays private to
   each) and still do the join here rather than reimplementing k-means in JS:

     run_align_worker()  aligns ONE sample and posts its split pairs
     run_join()          takes both samples' breakpoints, returns ind + colour

   run_webtest() keeps the original single-instance behaviour.                */
#include <emscripten.h>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "../core/align_wrap.h"
#include "../core/arg_parse.h"
#include "../core/post.h"
#include "../core/sink.h"
using namespace std;

struct PostSink : Sink {
	PostSink(int emit) : Sink(emit) {}
	void write(int kind, const string &line) {
		EM_ASM_({ postMessage([0, UTF8ToString($0), "p_read", UTF8ToString($1)]); },
		        line.c_str(), emit_name(kind));
	}
};

struct PostProgress : Progress {
	void pct(const char *sample, int p) {
		EM_ASM_({ postMessage([0, $0, UTF8ToString($1)]); }, p, sample);
	}
	void result(const char *sample, double align, double ld, double si) {
		EM_ASM_({ postMessage([0, UTF8ToString($0), "result", $1, $2, $3, 0]); },
		        sample, align, ld, si);
	}
	void log(const char *msg) { printf("###### %s\n", msg); }
	void done() { EM_ASM({ postMessage([0, "end"]); }); }
};

static int fill_args(Args &a, char *ref, char *c1, char *c2, char *m1, char *m2,
                     char *st, char *ed, char *cle, char *gene, char *sample,
                     char *K, char *emit_spec)
{
	char *t[14];
	t[0]  = (char *)"excas";
	t[1]  = ref; t[2] = c1; t[3] = c2; t[4] = m1; t[5] = m2;
	t[6]  = st;  t[7] = ed; t[8] = cle;
	t[9]  = gene; t[10] = sample; t[11] = K; t[12] = (char *)"out";
	t[13] = NULL;                       /* dist_cut omitted -> default 1000 */
	if (Arg_Parse(&a, 13, t) != 0) return 1;
	a.datadir = "/data/";               /* WORKERFS mount point */
	int m = emit_parse((emit_spec && *emit_spec) ? emit_spec : "split");
	a.emit = (m < 0) ? EMIT_SPLIT : m;
	return 0;
}

extern "C" {

/* --- original single-instance path ------------------------------------- */
void run_webtest(char *ref, char *c1, char *c2, char *m1, char *m2,
                 char *st, char *ed, char *cle, char *gene, char *sample,
                 char *K, char *outdir, char *emit_spec)
{
	(void)outdir;
	Args a;
	if (fill_args(a, ref, c1, c2, m1, m2, st, ed, cle, gene, sample, K, emit_spec)) {
		EM_ASM({ postMessage([0, "end"]); });
		return;
	}
	PostSink sink(a.emit);
	PostProgress prog;
	align_wrap run(a, sink, prog);
}

/* --- one sample, for the two-worker path -------------------------------- */
/* role is "ctr" or "mut"; f1/f2 are that sample's reads. */
void run_align_worker(char *ref, char *f1, char *f2,
                      char *st, char *ed, char *cle, char *gene, char *sample,
                      char *K, char *emit_spec, char *role)
{
	int which = (role && role[0] == 'm') ? 1 : 0;
	Args a;
	if (fill_args(a, ref, f1, f2, f1, f2, st, ed, cle, gene, sample, K, emit_spec)) {
		EM_ASM({ postMessage([0, "end"]); });
		return;
	}
	PostSink sink(a.emit);
	PostProgress prog;
	result res = result();
	vector<read> split[2];
	if (!align_sample(a, which, sink, prog, res, split)) {
		EM_ASM({ postMessage([0, "end"]); });
		return;
	}

	/* raw counters; the rates need mut_only, which only the join knows */
	EM_ASM_({ postMessage([0, "counters", UTF8ToString($0), $1, $2, $3, $4, $5, $6]); },
	        role, res.align_rate, (double)res.t_align, (double)res.tg_align_0,
	        (double)res.tg_ins, (double)res.tg_del, (double)res.only);

	/* split pairs: the two SAM lines (still without the YC:Z tag) and the
	   breakpoint the clustering works on */
	vector<array<double,2>> pts;
	split_breakpoints(split, pts);
	for (size_t i = 0; i < split[0].size(); i++) {
		string l1 = format_read(split[0][i]), l2 = format_read(split[1][i]);
		EM_ASM_({ postMessage([0, "split", UTF8ToString($0), UTF8ToString($1),
		                       UTF8ToString($2), $3, $4]); },
		        role, l1.c_str(), l2.c_str(), pts[i][0], pts[i][1]);
	}
	EM_ASM_({ postMessage([0, "aligned", UTF8ToString($0)]); }, role);
}

/* --- the join ------------------------------------------------------------ */
/* Points arrive as "x,y x,y ..."; posts one line per pair, "ind\tYC:Z:r,g,b". */
static void parse_points(const char *s, vector<array<double,2>> &out)
{
	out.clear();
	if (s == NULL) return;
	const char *p = s;
	while (*p) {
		char *end;
		double x = strtod(p, &end);
		if (end == p) break;
		if (*end == ',') end++;
		p = end;
		double y = strtod(p, &end);
		if (end == p) break;
		p = end;
		out.push_back({x, y});
		while (*p == ' ' || *p == '\n') p++;
	}
}

void run_join(char *ctr_pts, char *mut_pts)
{
	vector<array<double,2>> p1, p2;
	parse_points(ctr_pts, p1);
	parse_points(mut_pts, p2);
	ClusterResult c1 = cluster_points(p1), c2 = cluster_points(p2);
	vector<char> i1, i2;
	mutual_exclusion(c1, c2, i1, i2);

	string out1, out2;
	long long ctr_only_n = 0, mut_only_n = 0;
	for (size_t i = 0; i < c1.labels.size(); i++) {
		ctr_only_n += (i1[i] == 0);
		out1 += (i1[i] ? "1\t" : "0\t") + cluster_color(c1.labels[i]) + "\n";
	}
	for (size_t i = 0; i < c2.labels.size(); i++) {
		mut_only_n += (i2[i] == 0);
		out2 += (i2[i] ? "1\t" : "0\t") + cluster_color(c2.labels[i]) + "\n";
	}
	EM_ASM_({ postMessage([0, "join", "ctr", UTF8ToString($0), $1]); }, out1.c_str(), (double)ctr_only_n);
	EM_ASM_({ postMessage([0, "join", "mut", UTF8ToString($0), $1]); }, out2.c_str(), (double)mut_only_n);
	EM_ASM({ postMessage([0, "end"]); });
}

}
