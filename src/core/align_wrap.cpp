#include<stdio.h>
#include<iostream>
#include<list>
#include<stdlib.h>
#include<unordered_map>
#include<algorithm>
#include<string>
#include<string.h>
#include<vector>
#include<time.h>
#include<stdlib.h>
#include<sys/stat.h>
#include "align.h"
#include "func.h"
#include "post.h"
#include "align_wrap.h"
#include "sink.h"
int skip_read,is_skip,skip;//mm=-1
string reff,dir;
char view_q[1000];
long long int mut_only=0,ctr_only=0;
Args arg;
char g_ref_seq[100000],bf[1000];
inline int read_fastq(FILE *in,read &rr)
{
	int flag=0;
	flag = fgets(rr.query_name,250,in)!=NULL;
			
	int l = strlen(rr.query_name);
	for(int i=0;i<l;i++)
	{
		if(rr.query_name[i]==' ')
		{
			rr.query_name[i]='\0';
		}
	}
	if(l>0) rr.query_name[l-1] = '\0';   /* l==0 wrote query_name[-1] */
	flag = fgets(rr.query,250,in)!=NULL;
	l = strlen(rr.query);
	if(l>0) rr.query[l-1] = '\0';
	flag = fgets(bf,250,in)!=NULL;
	rr.strand=bf[0];
	flag = fgets(rr.query_qual,250,in)!=NULL;
	if(l>0) rr.query_qual[l-1] = '\0';
	return flag;
}
int align_pair(int K, unordered_map<long long int,vector<int>> &kmer_ref, unordered_map<long long int,vector<int>> &kmer_ref2, read &r1, read &r2,vector<read> *split,int rev,result &res_cont)
{
//	printf("align_pair_start\n"); //log
	r1.map_f=-1,r2.map_f=-1;
	align_query(K,kmer_ref,r1,kmer_ref2,0);
	if(r1.map_f==-1){return 0;}
	align_query(K,kmer_ref,r2,kmer_ref2,0);

	if(r2.map_f==-1){return 0;}
	if(r1.map_f!=-1 && r2.map_f!=-1)
	{
		if(r1.map_f==1 && r2.map_f==1)
		{
			make_sam(r1,r2,rev,arg.genename);
			int res = long_inner_check(r1,r2,split,res_cont);

			//printf("\nnormal pair %d\n",res); //log
			if(res==-1){return 0;}
			if(res>=400)
			{
				align_query(K,kmer_ref,r1,kmer_ref2,1);
				align_query(K,kmer_ref,r2,kmer_ref2,1);
			}
		}
		if(r1.map_f==2 || r2.map_f==2){

			align_query(K,kmer_ref,r1,kmer_ref2,1);
			align_query(K,kmer_ref,r2,kmer_ref2,1);
			if(r1.map_f!=-1&&r2.map_f!=-1){

//			printf("split pair\n"); //log
			make_sam(r1,r2,rev,arg.genename);
			split[0].push_back(r1);
			split[1].push_back(r2);
			}
		}
	}
	return 1;
}
void process(char *l_fastq, char *r_fastq, vector<read> *split, unordered_map<long long int,vector<int>> &kmer_ref,int K,unordered_map<long long int,vector<int>> &kmer_ref2,int full_kind,result &res_cont,Sink &sink,Progress &prog)
{
	FILE *f2,*f3;
	read r1,r2;
	double inp_time=0,align_t=0,t1=0,t2=0;
	int total_align=0,total_read=0,ret=0,progress=5;
	struct stat buf;
	for(int i=0;i<4;i++){r1.tags.push_back("");r2.tags.push_back("");}
	f2 = fopen((dir+string(l_fastq)).c_str(),"r");
	f3 = fopen((dir+string(r_fastq)).c_str(),"r");
	ret = stat((dir+string(l_fastq)).c_str(), &buf);
	int line = buf.st_size/300;

	if(line<1) line=1;
	t2=clock();
	while(read_fastq(f2,r1)&&read_fastq(f3,r2))
	{
		inp_time+=(clock()-t2)/CLOCKS_PER_SEC;
		int pair_res=0;
		total_read++;
		if(((total_read*100))>=progress*line){
			if(progress<100) prog.pct(res_cont.query_type.c_str(), progress);
			progress+=5;
		}
		t1=clock();
		parse_read(r1);    
		parse_read(r2);
		reverse_comp(r1.query,r1.query_len);
		is_skip=0;
		pair_res=align_pair(K,kmer_ref,kmer_ref2,r1,r2,split,1,res_cont);    /// r1: rev   r2: normal
		//printf("%d %d\n",r1.map_f,r2.map_f);
		if(pair_res==0){
		//	printf("reverse_align\n"); //log
			reverse_comp(r2.query,r2.query_len);
			reverse_comp(r1.query,r1.query_len);
			pair_res=align_pair(K,kmer_ref,kmer_ref2,r1,r2,split,0,res_cont); // r1 : normal r2: rev
			//printf("%d %d\n",r1.map_f,r2.map_f);
		}
		if(pair_res==1 && sink.wants(full_kind))
		{
			sink.write(full_kind, format_read(r1));
			sink.write(full_kind, format_read(r2));
		}
		total_align+=pair_res;
		skip+=is_skip;
		t2=clock();
	}
	prog.pct(res_cont.query_type.c_str(), 100);
	res_cont.t_align  = total_align;
	res_cont.align_rate = (double)((double)total_align/(double)total_read);
	fclose(f2);
	fclose(f3);
}
/* large_del / small_indel are reported exactly as the canonical web build computed them.
   Note the CLI build used to add r.tg_miss (substitution-only reads) into the small-indel
   numerator; that is a different quantity and is deliberately not carried over. */
void report_ctr(result &r, Progress &prog)
{
	prog.result("ctr", r.align_rate,
		(double)((double)r.tg_align_0/(double)r.t_align),
		(double)((double)r.tg_ins+r.tg_del)/(double)(r.tg_align_0+r.only));
}
void report_mut(result &r, long long only, Progress &prog)
{
	prog.result("mut", r.align_rate,
		(double)((double)only/(double)(r.tg_align_0+only)),
		(double)((double)r.tg_ins+r.tg_del)/(double)(r.tg_align_0+only));
}
/* Align exactly one sample. Pulled out of the constructor so the browser can run
   the control and the treated sample in two separate workers (two wasm instances,
   so the file-scope state above stays private to each). The k-mer index is rebuilt
   per call, which costs nothing measurable: 0.00s of CPU for a 9.4 kb reference. */
bool align_sample(Args &a, int which, Sink &sink, Progress &prog,
                  result &res_out, vector<read> *split_out)
{
	unordered_map<long long int,vector<int>> kmer_ref,kmer_ref2;
	arg = a;
	dir = arg.datadir;
	g_sink = &sink;
	skip = 0;

	if(make_ref(arg.refname,arg.K,kmer_ref) < 0 ||   /* K-mer reference */
	   make_ref(arg.refname,17,kmer_ref2) < 0){      /* 17-mer reference, assumed unique */
		prog.log("reference could not be read - aborting");
		return false;
	}
	if(which==0){
		res_out.query_type = "ctr";
		g_site_kind = K_SITE_CTR; g_site_sm_kind = K_SITE_CTR_SM;
		process(arg.ctrfastq_l,arg.ctrfastq_r,split_out,kmer_ref,arg.K,kmer_ref2,K_CTR_FULL,res_out,sink,prog);
	}else{
		res_out.query_type = "mut";
		g_site_kind = K_SITE_MUT; g_site_sm_kind = K_SITE_MUT_SM;
		process(arg.mutfastq_l,arg.mutfastq_r,split_out,kmer_ref,arg.K,kmer_ref2,K_MUT_FULL,res_out,sink,prog);
	}
	g_site_kind = g_site_sm_kind = -1;
	return true;
}

void report_sample(result &r, long long only, Progress &prog)
{
	if(r.query_type=="ctr") report_ctr(r,prog); else report_mut(r,only,prog);
}

align_wrap::align_wrap(Args &a, Sink &sink, Progress &prog)
{
	vector<read> split_ctr[2],split_mut[2];
	/* result 의 첫 멤버가 std::string 이라 ={0,} 은 널 문자열 생성이 된다 (libstdc++ 에서 throw) */
	result ctr_r = result(), mut_r = result();
	mut_only = ctr_only = 0;

	prog.log("aligning control sample");
	if(!align_sample(a,0,sink,prog,ctr_r,split_ctr)){ g_sink=NULL; sink.finish(); prog.done(); return; }
	report_ctr(ctr_r,prog);

	prog.log("aligning Cas9-treated sample");
	if(!align_sample(a,1,sink,prog,mut_r,split_mut)){ g_sink=NULL; sink.finish(); prog.done(); return; }

	prog.log("post-processing (clustering)");
	post_process(split_ctr,split_mut,arg,sink);
	report_mut(mut_r,mut_only,prog);

	g_sink = NULL;
	sink.finish();
	prog.done();
}
align_wrap::~align_wrap(){
}
