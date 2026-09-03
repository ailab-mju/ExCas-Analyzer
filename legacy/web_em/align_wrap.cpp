#include<stdio.h>
#include<thread>
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
#include <emscripten.h>
#include "align.h"
#include "func.h"
#include "post.h"
#include "align_wrap.h"
int skip_read,is_skip,skip;//mm=-1
string reff,dir;
char view_q[1000];
long long int mut_only=0,ctr_only=0;
Args arg;
FILE *f_out,*f_out2;
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
	rr.query_name[l-1] = '\0';
	flag = fgets(rr.query,250,in)!=NULL;
	l = strlen(rr.query);
	rr.query[l-1] = '\0';
	flag = fgets(bf,250,in)!=NULL;
	rr.strand=bf[0];
	flag = fgets(rr.query_qual,250,in)!=NULL;
	rr.query_qual[l-1] = '\0';
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
void process(char *l_fastq, char *r_fastq, vector<read> *split, unordered_map<long long int,vector<int>> &kmer_ref,int K,unordered_map<long long int,vector<int>> &kmer_ref2,string output,result &res_cont)
{
	FILE *f2,*f3,*f4;
	read r1,r2;
	double inp_time=0,align_t=0,t1=0,t2=0;
	int total_align=0,total_read=0,ret=0,progress=5;
	struct stat buf;
	for(int i=0;i<4;i++){r1.tags.push_back("");r2.tags.push_back("");}
	f2 = fopen((dir+string(l_fastq)).c_str(),"r");
	f3 = fopen((dir+string(r_fastq)).c_str(),"r");
	ret = stat((dir+string(l_fastq)).c_str(), &buf);
	int line = buf.st_size/300;

	printf("file size: %d\n estimated line: %d\n",buf.st_size,line);
//	f4 = fopen(output.c_str(),"w");
//	mm_tot=mm_cnt=0;	
	t2=clock();
	printf("start_align %d\n",K);
	while(read_fastq(f2,r1)&&read_fastq(f3,r2))
	{
		inp_time+=(clock()-t2)/CLOCKS_PER_SEC;
		int pair_res=0;
		total_read++;
		if(((total_read*100))>=progress*line){
//			if(progress>100){progress=100;}
			printf("progress: %d%%\n",progress);
			EM_ASM_({ postMessage([0,$0,UTF8ToString($1)]);}, progress, res_cont.query_type.c_str());
			progress+=5;
		}
			//printf("%d\n",total_read);
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
	/*	if(f4!=NULL&&pair_res==1)
		{
			print_read(r1,f4);
			print_read(r2,f4);
		}*/
		total_align+=pair_res;
		skip+=is_skip;
		t2=clock();
	}
	progress=100;
	EM_ASM_({ postMessage([0,$0,UTF8ToString($1)]);}, progress, res_cont.query_type.c_str());
	printf("total_align: %d\n total_read %d\n alignment_rate %.3lf\n",total_align,total_read,(double)((double)total_align/(double)total_read));
	res_cont.t_align = total_align;
	res_cont.align_rate = (double)((double)total_align/(double)total_read);
	printf("align_t : %lf\ninp_time : %lf\n",align_t/CLOCKS_PER_SEC,inp_time);
	fclose(f2);
	fclose(f3);
//	fclose(f4);
}
void print_result(string type,result r)
{
	string signal ="result";
	if(type=="ctr"){
	EM_ASM_({ postMessage([0,UTF8ToString($0),UTF8ToString($1),$2,$3,$4,$5]);}, type.c_str(), signal.c_str(), r.align_rate, (double)((double)r.tg_align_0/(double)r.t_align),(double)((double)r.tg_ins+r.tg_del)/(double)(r.tg_align_0+r.only),0  );
}
	printf("%s large del size : %d\n",type.c_str(),r.large_d);
	printf("%s normal : %d %lf\n",type.c_str(),r.normal,(double)((double)r.normal/(double)r.t_align));
	printf("%s tg_align_0 : %d %lf\n",type.c_str(),r.tg_align_0,(double)((double)r.tg_align_0/(double)r.t_align));
//	printf("%s tg_align_1 : %d %lf\n",type.c_str(),r.tg_align_1,(double)((double)r.tg_align_1/(double)r.t_align));

	if(type=="ctr"){
	printf("%s small_indel : %d %d %.4lf %d %.4lf %d %.4lf %d %.4lf\n",type.c_str(),r.tg_align_0+r.only,r.tg_ins,(double)((double)r.tg_ins)/(double)(r.tg_align_0+r.only), r.tg_del, (double)((double)r.tg_del)/(double)(r.tg_align_0+r.only), r.tg_ins+r.tg_del,(double)((double)r.tg_ins+r.tg_del)/(double)(r.tg_align_0+r.only),r.tg_ins+r.tg_del+r.tg_miss,(double)((double)(r.tg_ins+r.tg_del+r.tg_miss))/((double)(r.tg_align_0+r.only)));}
	printf("%s skip %d\n",type.c_str(),skip);
}
align_wrap::align_wrap(char *argv[])
{
	vector<read> split_ctr[2],split_mut[2];
	unordered_map<long long int,vector<int>> kmer_ref,kmer_ref2;
	result ctr_r={0,},mut_r={0,};
	double t1 = clock(),t2=clock();
	Arg_Parse(&arg,argv);
	dir="/data/";
	/// make ref
	printf("make ref start");
	make_ref(arg.refname,arg.K,kmer_ref);  //k-mer reference
	make_ref(arg.refname,17,kmer_ref2); //17-mer reference unique assume
	printf("ref time %lf\n",(clock()-t1)/CLOCKS_PER_SEC);	
	/// ref end
	t1=clock();
	skip=0;

	//ctr align
	//if(arg.ctrsite!=""){f_out=fopen(arg.ctrsite.c_str(),"w");f_out2=fopen(arg.ctrsite_sm.c_str(),"w");}
	ctr_r.query_type="ctr";
	process(arg.ctrfastq_l,arg.ctrfastq_r,split_ctr,kmer_ref,arg.K,kmer_ref2,arg.resctr,ctr_r);
	printf("ctr align time : %lf\n",(clock()-t1)/CLOCKS_PER_SEC);
	print_result("ctr",ctr_r);
	t1 = clock();
	skip=0;
	
	//mut align
	//if(arg.mutsite!=""){f_out=fopen(arg.mutsite.c_str(),"w");f_out2=fopen(arg.mutsite_sm.c_str(),"w");}

	mut_r.query_type="mut";
	process(arg.mutfastq_l,arg.mutfastq_r,split_mut,kmer_ref,arg.K,kmer_ref2,arg.resmut,mut_r);
	printf("mut align time : %lf\n",(clock()-t1)/CLOCKS_PER_SEC);
	print_result("mut",mut_r);
	t1 = clock();
	post_process(split_ctr,split_mut,arg);
	
	
	
	printf("post time : %lf\n",(clock()-t1)/CLOCKS_PER_SEC);
	printf("whole time : %lf\n",(clock()-t2)/CLOCKS_PER_SEC);

	//printf("mut tg_align_2 : %d %lf--%.3lf\n",mut_r.tg_align_2,(double)((double)mut_r.tg_align_2/(double)mut_r.t_align),(double)((double)mut_only)/(double)(mut_r.tg_align_2+mut_only));
	printf("mut tg_align_0 : %d %lf--%.3lf\n",mut_r.tg_align_0,(double)((double)mut_r.tg_align_0/(double)mut_r.t_align),(double)((double)mut_only)/(double)(mut_r.tg_align_0+mut_only));
	//printf("mut tg_align_1 : %d %lf--%.3lf\n",mut_r.tg_align_1,(double)((double)mut_r.tg_align_1/(double)mut_r.t_align),(double)((double)mut_only)/(double)(mut_r.tg_align_1+mut_only));

	printf("%d \n",mut_r.tg_del);
	string type="mut";
	string signal="result";
	EM_ASM_({ postMessage([0,UTF8ToString($0),UTF8ToString($1),$2,$3,$4,$5]);}, type.c_str(), signal.c_str(), mut_r.align_rate, (double)((double)mut_only/(double)(mut_r.tg_align_0+mut_only)),(double)((double)mut_r.tg_ins+mut_r.tg_del)/(double)(mut_r.tg_align_0+mut_only),0 );

	printf("mut small_indel : %lld %d %.4lf %d %.4lf %d %.4lf %d %.4lf\n",mut_r.tg_align_0+mut_only,mut_r.tg_ins,(double)((double)mut_r.tg_ins)/(double)(mut_r.tg_align_0+mut_only), mut_r.tg_del, (double)((double)mut_r.tg_del)/(double)(mut_r.tg_align_0+mut_only), mut_r.tg_ins+mut_r.tg_del,(double)((double)mut_r.tg_ins+mut_r.tg_del)/(double)(mut_r.tg_align_0+mut_only),mut_r.tg_ins+mut_r.tg_del+mut_r.tg_miss,((double)(mut_r.tg_ins+mut_r.tg_del+mut_r.tg_miss))/((double)mut_r.tg_align_0+mut_only));
}
align_wrap::~align_wrap(){
}
