/* Native front-end: Sink writes SAM files, Progress writes to stdout. */
#include <stdio.h>
#include <string>
#include <map>
#include "../core/align_wrap.h"
#include "../core/arg_parse.h"
#include "../core/sink.h"
using namespace std;

struct FileSink : Sink {
	string path[K_N];
	FILE  *fp[K_N];
	FileSink(int emit, Args &a) : Sink(emit) {
		for (int i = 0; i < K_N; i++) fp[i] = NULL;
		path[K_CTR_ONLY]=a.resctronly; path[K_CTR_COM]=a.resctrcom;
		path[K_MUT_ONLY]=a.resmutonly; path[K_MUT_COM]=a.resmutcom;
		path[K_CTR_FULL]=a.resctr;     path[K_MUT_FULL]=a.resmut;
		path[K_SITE_CTR]=a.ctrsite;    path[K_SITE_MUT]=a.mutsite;
		path[K_SITE_CTR_SM]=a.ctrsite_sm; path[K_SITE_MUT_SM]=a.mutsite_sm;
		/* create every enabled stream up front so empty runs still produce the file,
		   matching the original tool's behaviour */
		for (int i = 0; i < K_N; i++) if (wants(i)) open(i);
	}
	void open(int k) {
		if (fp[k] == NULL) {
			fp[k] = fopen(path[k].c_str(), "w");
			if (fp[k] == NULL) fprintf(stderr, "cannot write %s\n", path[k].c_str());
		}
	}
	void write(int kind, const string &line) {
		open(kind);
		if (fp[kind]) fprintf(fp[kind], "%s\n", line.c_str());
	}
	void finish() {
		for (int i = 0; i < K_N; i++) if (fp[i]) { fclose(fp[i]); fp[i] = NULL; }
	}
	~FileSink() { finish(); }
};

struct StdoutProgress : Progress {
	int last[2];
	StdoutProgress() { last[0] = last[1] = -1; }
	void log(const char *msg) { printf("###### %s\n", msg); fflush(stdout); }
	void pct(const char *sample, int p) {
		int i = (sample[0] == 'c') ? 0 : 1;
		if (p == last[i]) return;
		last[i] = p;
		printf("\r       %s %3d%%", sample, p);
		if (p >= 100) printf("\n");
		fflush(stdout);
	}
	void result(const char *sample, double align, double ld, double si) {
		const char *name = (sample[0] == 'c') ? "Control sample" : "Cas9 Mutation sample";
		printf("%s alignment rate      : %.2lf%%\n", name, align * 100);
		printf("%s Large Deletion Rate : %.2lf%%\n", name, ld * 100);
		printf("%s small indel rate    : %.2lf%%\n", name, si * 100);
	}
};

static void usage(const char *p) {
	fprintf(stderr,
	  "usage: %s ref.fa ctr_1.fastq ctr_2.fastq mut_1.fastq mut_2.fastq \\\n"
	  "          tg_start tg_end cleavage_offset gene sample K outdir [dist_cut]\n"
	  "          [--emit=all|none|split[,full][,site]] [--pre-chunk] [--datadir=PATH]\n"
	  "\n"
	  "  --emit    which SAM streams to write (default: all).\n"
	  "            The browser build defaults to --emit=split; use it here to\n"
	  "            produce byte-identical output for comparison.\n"
	  "  --pre-chunk  enable the pre_chunk() fast path. OFF by default: it changes\n"
	  "            results and the published numbers were produced without it.\n", p);
}

int main(int argc, char *argv[])
{
	Args a;
	if (argc < 13) { usage(argv[0]); return 2; }
	if (Arg_Parse(&a, argc, argv) != 0) { usage(argv[0]); return 2; }
	FileSink sink(a.emit, a);
	StdoutProgress prog;
	align_wrap run(a, sink, prog);
	return 0;
}
