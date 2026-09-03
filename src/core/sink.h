#ifndef EXCAS_SINK_H
#define EXCAS_SINK_H
#include <string>

/* The only three things that differ between the CLI build and the browser build.
   Everything in src/core is written against these and knows nothing about
   FILE*, emscripten, or postMessage. */

/* Output streams the analyzer can produce. */
enum OutKind {
	K_CTR_ONLY = 0, K_CTR_COM, K_MUT_ONLY, K_MUT_COM,   /* split reads  (EMIT_SPLIT) */
	K_CTR_FULL, K_MUT_FULL,                             /* every aligned pair (EMIT_FULL) */
	K_SITE_CTR, K_SITE_MUT,                             /* on-target, no indel (EMIT_SITE) */
	K_SITE_CTR_SM, K_SITE_MUT_SM,                       /* on-target, with indel (EMIT_SITE) */
	K_N
};

/* Emit groups. The browser build defaults to EMIT_SPLIT only: the full
   alignment SAM is ~99.8% of the output volume and does not fit in a tab. */
#define EMIT_SPLIT 1
#define EMIT_FULL  2
#define EMIT_SITE  4
#define EMIT_ALL   (EMIT_SPLIT|EMIT_FULL|EMIT_SITE)

inline int emit_group(int kind)
{
	if (kind <= K_MUT_COM)  return EMIT_SPLIT;
	if (kind <= K_MUT_FULL) return EMIT_FULL;
	return EMIT_SITE;
}
const char *emit_name(int kind);
int  emit_parse(const char *spec);   /* "split,full,site" | "all" | "none"; -1 on error */

struct Sink {
	int emit;
	Sink(int e = EMIT_ALL) : emit(e) {}
	virtual ~Sink() {}
	bool wants(int kind) const { return (emit & emit_group(kind)) != 0; }
	/* one SAM record, no trailing newline */
	virtual void write(int kind, const std::string &line) = 0;
	virtual void finish() {}
};

struct Progress {
	virtual ~Progress() {}
	virtual void pct(const char *sample, int p) { (void)sample; (void)p; }
	virtual void result(const char *sample, double align_rate, double large_del, double small_indel)
	                  { (void)sample; (void)align_rate; (void)large_del; (void)small_indel; }
	virtual void log(const char *msg) { (void)msg; }
	virtual void done() {}
};
#endif
