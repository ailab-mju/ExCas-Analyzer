#ifndef EXCAS_ALIGN_WRAP_H
#define EXCAS_ALIGN_WRAP_H
#include <vector>
#include "arg_parse.h"
#include "sink.h"
#include "func.h"

/* one sample only; which = 0 (control) or 1 (Cas9-treated) */
bool align_sample(Args &a, int which, Sink &sink, Progress &prog,
                  result &res_out, std::vector<read> *split_out);
void report_sample(result &r, long long only, Progress &prog);

class align_wrap
{
	public:
	~align_wrap();
	align_wrap(Args &a, Sink &sink, Progress &prog);
};
#endif
