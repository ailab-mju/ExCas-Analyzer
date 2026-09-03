#include "sink.h"
#include <string.h>

const char *emit_name(int kind)
{
	static const char *n[K_N] = {
		"ctr_only", "ctr_com", "mut_only", "mut_com",
		"ctr", "mut",
		"site_ctr", "site_mut", "site_ctr_sm", "site_mut_sm"
	};
	return (kind >= 0 && kind < K_N) ? n[kind] : "?";
}

int emit_parse(const char *spec)
{
	if (spec == NULL || *spec == '\0') return -1;
	if (strcmp(spec, "all")  == 0) return EMIT_ALL;
	if (strcmp(spec, "none") == 0) return 0;
	int mask = 0;
	std::string s(spec), tok;
	size_t p = 0;
	while (p <= s.size()) {
		size_t c = s.find(',', p);
		if (c == std::string::npos) c = s.size();
		tok = s.substr(p, c - p);
		if      (tok == "split") mask |= EMIT_SPLIT;
		else if (tok == "full")  mask |= EMIT_FULL;
		else if (tok == "site")  mask |= EMIT_SITE;
		else if (tok != "")      return -1;
		p = c + 1;
	}
	return mask;
}
