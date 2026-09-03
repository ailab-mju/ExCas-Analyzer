#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include "post.h"
#include "dkm.hpp"
#include "func.h"

using namespace std;

static const char *CLUSTER_COLORS[] = {
        "#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
        "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
        "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
        "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
        "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
        "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
        "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
        "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
        "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
        "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
        "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
        "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
        "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C",
        "#83AB58", "#001C1E", "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
        "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", "#1E0200", "#5B4E51",
        "#C895C5", "#320033", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379", "#012C58"
};

static const size_t N_CLUSTER_COLORS = sizeof(CLUSTER_COLORS)/sizeof(CLUSTER_COLORS[0]);

/* The colour table has 128 entries but k could reach 498, so a large clustering
   used to index past the end of the array. */
string cluster_color(uint32_t label)
{
	return "YC:Z:" + rgb(string(CLUSTER_COLORS[label % N_CLUSTER_COLORS]));
}

/* ---- pure clustering, no read structs -------------------------------------
   Split into point-based pieces so the browser can run the two alignments in
   separate workers (separate wasm instances) and still join them here, in the
   one C++ implementation, instead of reimplementing k-means in JavaScript. */

ClusterResult cluster_points(const vector<array<double,2>> &point)
{
	ClusterResult out;
	if(point.size()==0) return out;
	int cl_k=1,tt=1;
	double min_error=numeric_limits<double>::max(),min_k=1,reg=2000,error;
	while(cl_k<500)
	{
		error=0;
		if(cl_k>(int)point.size()){break;}
		auto cluster_data = dkm::kmeans_lloyd(point,cl_k);
		auto label = std::get<1>(cluster_data);
		auto center = std::get<0>(cluster_data);
		for(size_t i=0;i<point.size();i++)
		{
			error+=((center[label[i]][0]-point[i][0])*(center[label[i]][0]-point[i][0]))+((center[label[i]][1]-point[i][1])*(center[label[i]][1]-point[i][1]));
		}
		error/=(double)point.size();
		error = error+reg*cl_k;
		if(error<min_error)
		{
			min_error = error;
			min_k = cl_k;
		}
		cl_k=3*tt;
		tt++;
	}
	auto cluster_data = dkm::kmeans_lloyd(point,min_k);
	out.centers = std::get<0>(cluster_data);
	out.labels  = std::get<1>(cluster_data);
	return out;
}

/* A cluster is "common" when a cluster of the other sample sits within the
   tolerance, i.e. the deletion is not specific to the treated sample.
   (was bool[20000] on the stack -- 200KB, above Emscripten's default) */
void mutual_exclusion(const ClusterResult &a, const ClusterResult &b,
                      vector<char> &ind_a, vector<char> &ind_b)
{
	vector<char> mutual1(a.centers.size(),0), mutual2(b.centers.size(),0);
	vector<int>  count1(a.centers.size(),0),  count2(b.centers.size(),0);
	for(size_t i=0;i<a.labels.size();i++) count1[a.labels[i]]++;
	for(size_t i=0;i<b.labels.size();i++) count2[b.labels[i]]++;
	for(size_t i=0;i<a.centers.size();i++)
	{
		for(size_t j=0;j<b.centers.size();j++)
		{
			if(count1[i]<=3||count2[j]<=3){continue;}
			double interval = min(abs(a.centers[i][0]-a.centers[i][1]),abs(b.centers[j][0]-b.centers[j][1]));
			interval*=0.05;
			interval=max(interval,(double)50);
			if(abs(a.centers[i][0]-b.centers[j][0])<interval && abs(a.centers[i][1]-b.centers[j][1])<interval)
			{
				mutual1[i]=1;
				mutual2[j]=1;
			}
		}
	}
	ind_a.assign(a.labels.size(),0);
	ind_b.assign(b.labels.size(),0);
	for(size_t i=0;i<a.labels.size();i++) ind_a[i]=mutual1[a.labels[i]];
	for(size_t i=0;i<b.labels.size();i++) ind_b[i]=mutual2[b.labels[i]];
}

/* the breakpoint of a split pair: whichever mate carries the split alignment */
void split_breakpoints(vector<read> *split, vector<array<double,2>> &point)
{
	point.clear();
	for(size_t i=0;i<split[0].size();i++)
	{
		const read &tmp = (split[0][i].map_f==2) ? split[0][i] : split[1][i];
		point.push_back({(double)tmp.ref_st,(double)tmp.ref_ed});
	}
}

/* ---- read-struct wrappers (native path) ---------------------------------- */

void read_cluster(vector<read> *split_ctr, vector<read> *split_mut,Args args)
{
	vector<array<double,2>> p1,p2;
	split_breakpoints(split_ctr,p1);
	split_breakpoints(split_mut,p2);
	ClusterResult c1 = cluster_points(p1), c2 = cluster_points(p2);
	vector<char> i1,i2;
	mutual_exclusion(c1,c2,i1,i2);
	for(size_t i=0;i<c1.labels.size();i++)
	{
		string col = cluster_color(c1.labels[i]);
		split_ctr[0][i].ind = i1[i]; split_ctr[0][i].tags.push_back(col);
		split_ctr[1][i].ind = i1[i]; split_ctr[1][i].tags.push_back(col);
	}
	for(size_t i=0;i<c2.labels.size();i++)
	{
		string col = cluster_color(c2.labels[i]);
		split_mut[0][i].ind = i2[i]; split_mut[0][i].tags.push_back(col);
		split_mut[1][i].ind = i2[i]; split_mut[1][i].tags.push_back(col);
	}
}

void write_read(vector<read> *split, Sink &sink, int kind_only, int kind_com)
{
	int kind[2];
	kind[0] = kind_only;
	kind[1] = kind_com;
	int c=0;
	for(size_t i=0;i<split[0].size();i++)
	{
		c+=(!split[0][i].ind);
		if(sink.wants(kind[split[0][i].ind])) sink.write(kind[split[0][i].ind], format_read(split[0][i]));
		if(sink.wants(kind[split[1][i].ind])) sink.write(kind[split[1][i].ind], format_read(split[1][i]));
	}
	mut_only=c;
}

void post_process(vector<read> *split_ctr, vector<read> *split_mut,Args args,Sink &sink)
{
	read_cluster(split_ctr,split_mut,args);
	write_read(split_ctr,sink,K_CTR_ONLY,K_CTR_COM);
	ctr_only=mut_only;
	write_read(split_mut,sink,K_MUT_ONLY,K_MUT_COM);
}
