#include <string>
#include <vector>
#include <array>
#include <tuple>
#include <stdint.h>
#include "arg_parse.h"
#include "func.h"
#include "sink.h"

using namespace std;

struct ClusterResult {
	vector<array<double,2>> centers;
	vector<uint32_t>        labels;
};

/* point-based core, shared by the native path and the browser join */
ClusterResult cluster_points(const vector<array<double,2>> &point);
void mutual_exclusion(const ClusterResult &a, const ClusterResult &b,
                      vector<char> &ind_a, vector<char> &ind_b);
void split_breakpoints(vector<read> *split, vector<array<double,2>> &point);
string cluster_color(uint32_t label);

/* read-struct wrappers */
void read_cluster(vector<read> *split_ctr, vector<read> *split_mut,Args args);
void write_read(vector<read> *split, Sink &sink, int kind_only, int kind_com);
void post_process(vector<read> *split_ctr, vector<read> *split_mut,Args args,Sink &sink);
