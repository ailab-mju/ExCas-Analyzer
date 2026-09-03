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
#include "align.h"
#include "align_wrap.h"
#include "dkm.hpp"
int target_st,target_ed,mm=-1,skip,is_skip,dist_cut,st_flag=0;
string reff,dir;
char view_q[1000];
long long int mm2_cnt=0,mm2_tot=0,mm_cnt=0,mm_tot=0,tg_align=0,tg_align2=0,tg_align_50=0,tg_align_100=0,tg_align_150=0,t_align=0,mut_only=0,tg_align3=0,tg_align_0=0,normal=0,tg_align4=0,tg_align5=0,long_d=0,tg_align6=0,tg_align7=0,tg_align8=0;
vector<int> ppp,ppp2;
struct d2
{
	int x1,x2,y1,y2;
};

vector<d2> ppp3;
char ref_seq[100000],bf[1000];
FILE *f4,*f_out;
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
void read_cluster(vector<read> *split_ctr, vector<read> *split_mut,char *argv[])
{
	bool mutual1[20000],mutual2[20000];
	int count1[20000],count2[20000];
	std::tuple<std::vector<std::array<float, 2>>, std::vector<uint32_t>> cluster_data1 = cluster(split_ctr),cluster_data2 = cluster(split_mut);
	vector<std::array<float,2>> center1 = get<0>(cluster_data1),center2 = get<0>(cluster_data2);
	printf("%d %d\n",center1.size(),center2.size());
	vector<uint32_t> label1 = get<1>(cluster_data1),label2 = get<1>(cluster_data2);
	for(int i=0;i<center1.size();i++){mutual1[i]=0;count1[i]=0;}
	for(int i=0;i<center2.size();i++){mutual2[i]=0;count2[i]=0;}
	for(int i=0;i<split_ctr[0].size();i++){count1[label1[i]]++;}
	for(int i=0;i<split_mut[0].size();i++){count2[label2[i]]++;}
	
/*	for(int i=0;i<split_ctr[0].size();i++){printf("%d ",count1[i]);}
	printf("\n");
	for(int i=0;i<split_mut[0].size();i++){printf("%d ",count2[i]);}
	printf("\n");*/

	for(int i=0;i<center1.size();i++)
	{
		for(int j=0;j<center2.size();j++)
		{
			if(count1[i]<=5||count2[j]<=5){continue;}
			float interval = min(abs(center1[i][0]-center1[i][1]),abs(center2[j][0]-center2[j][1]));
			interval*=0.05;
			interval=max((float)interval,(float)50);
			if(abs(center1[i][0]-center2[j][0])<interval && abs(center1[i][1]-center2[j][1])<interval)
			{
				mutual1[i]=1;
				mutual2[j]=1;
			}
		}
	}
	if(argv[16]!=NULL)
	{
		FILE *ff;
		ff = fopen(argv[16],"r");
		char name[100];
		int ss,ee,j;
		vector <int> a1;
		vector <int> a2;
		while(fscanf(ff,"%s %d %d",name,&ss,&ee)!=EOF)
		{
			a1.push_back(ss);
			a2.push_back(ee);
		}
		for(int i=0;i<center1.size();i++)
		{
			float s1=center1[i][0];
			float e1=center1[i][1];
			for(j=0;j<a1.size();j++)
			{
				float s2 = a1[j];
				float e2 = a2[j];
				if(abs(s1-s2)<=100 && abs(e1-e2)<=100)
				{
					break;
				}
			}
			printf("cen\tctr\t%f\t%f\t%d\t%d",center1[i][0],center1[i][1],count1[i],mutual1[i]);
			if(j==a1.size())
			{
				printf("NO\n");
			}
			else
			{
				printf("YES\n");
			}
		}
		for(int i=0;i<center2.size();i++)
		{
			float s1=center2[i][0];
			float e1=center2[i][1];
			for(j=0;j<a1.size();j++)
			{
				float s2 = a1[j];
				float e2 = a2[j];
				if(abs(s1-s2)<=200 && abs(e1-e2)<=200)
				{
					break;
				}
			}
			printf("cen\tmut\t%f\t%f\t%d\t%d",center2[i][0],center2[i][1],count2[i],mutual2[i]);
			if(j==a1.size())
			{
				printf("NO\n");
			}
			else
			{
				printf("YES\n");
			}
		}
		fclose(ff);	
	}
	for(int i=0;i<label1.size();i++)
	{
		split_ctr[0][i].ind = mutual1[label1[i]];
		split_ctr[1][i].ind = mutual1[label1[i]];
	}
	for(int i=0;i<label2.size();i++)
	{
		split_mut[0][i].ind = mutual2[label2[i]];
		split_mut[1][i].ind = mutual2[label2[i]];
	}
}
void write_read(vector<read> *split, string only,string common)
{
	FILE *f[2];
	f[0] = fopen(only.c_str(),"w");
	f[1] = fopen(common.c_str(),"w");
	int c=0;
	for(int i=0;i<split[0].size();i++)
	{
		c+=(!split[0][i].ind);
		print_read(split[0][i],f[split[0][i].ind]);
		print_read(split[1][i],f[split[1][i].ind]);
	}
	mut_only=c;
	printf("%d\n",c);
	fclose(f[0]);
	fclose(f[1]);	
}
void post_process(vector<read> *split_ctr, vector<read> *split_mut,char *argv[])
{
	read_cluster(split_ctr,split_mut,argv);
	printf("ctr only large del size :");
	write_read(split_ctr,string(argv[10]),string(argv[11]));
	printf("mut only large del size :");
	write_read(split_mut,string(argv[12]),string(argv[13]));

}
int align_pair(int K, unordered_map<long long int,vector<int>> &cc, unordered_map<long long int,vector<int>> &cc2, read &r1, read &r2,vector<read> *split,int rev)
{
	r1.map_f=-1,r2.map_f=-1;
	align_query(K,cc,r1,cc2,0);
	if(r1.map_f==-1){return 0;}
	align_query(K,cc,r2,cc2,0);
	if(r2.map_f==-1){return 0;}
	if(r1.map_f!=-1 && r2.map_f!=-1)
	{
//		make_sam(r1,r2,0);
	//	if(rev==0&&r1.ref_st>r2.ref_st){return 0;}
	//	if(rev==1&&r2.ref_st>r1.ref_st){return 0;}
		if(r1.map_f==1 && r2.map_f==1)
		{
			make_sam(r1,r2,rev);	
			int inner_d,l_p,r_p,l_p2,r_p2;
			d2 tmp;
			tmp.x1 = r1.ref_st;
			tmp.y1 = r1.ref_ed;
			tmp.x2 = r2.ref_st;
			tmp.y2 = r2.ref_ed;
			if(r1.ref_st<r2.ref_st)
			{
				inner_d=r2.ref_ed-r1.ref_st;
				ppp2.push_back(r2.ref_ed-r1.ref_st);
				ppp3.push_back(tmp);
				l_p = r1.ref_st;
				r_p = r2.ref_ed;
				l_p2 = r1.ref_ed;
				r_p2 = r2.ref_st;
			}
			else
			{
				r_p = r1.ref_ed;
				l_p = r2.ref_st;
				r_p2 = r1.ref_st;
				l_p2 = r2.ref_ed;
				inner_d=r1.ref_ed-r2.ref_st;
				ppp2.push_back(r1.ref_ed-r2.ref_st);
				int tt;
				tt = tmp.x1;
				tmp.x1 = tmp.x2;
				tmp.x2 = tt;
				tt = tmp.y1;
				tmp.y1 = tmp.y2;
				tmp.y2 = tt;
				ppp3.push_back(tmp);			
			}
			if(r_p2-l_p2>400)
			{
				if(r_p2 < target_st || target_ed < l_p2){
					r1.overlap = 0;
					r2.overlap = 0;
					r1.map_f=r2.map_f=-1;return 0;
				}

				int ddist=10000;
				ddist = abs(r1.ref_st-target_st);
				ddist = abs(r1.ref_st-target_ed);
				ddist = abs(r1.ref_ed-target_st);
				ddist = abs(r1.ref_ed-target_ed);
				ddist = abs(r2.ref_st-target_st);
				ddist = abs(r2.ref_st-target_ed);
				ddist = abs(r2.ref_ed-target_st);
				ddist = abs(r2.ref_ed-target_ed);
				r1.map_f=2;r2.map_f=2;
				r1.ref_st=r2.ref_st=l_p;r1.ref_ed=r2.ref_ed=r_p;
				r1.del_size=r2.del_size=r_p-l_p;
				r1.dist=r2.dist=ddist;
				if(ddist>dist_cut){r1.map_f=r2.map_f=-1;return 0;}
				else{
				r1.tags[2]="Dist:"+to_string(r1.dist);
				r2.tags[2]="Dist:"+to_string(r2.dist);
				r1.tags[3]="D_size:"+to_string(r1.del_size);
				r2.tags[3]="D_size:"+to_string(r2.del_size);
				split[0].push_back(r1);
				split[1].push_back(r2);
				long_d++;
				return 1;
				}
			}

			normal++;
			if(!(r_p<target_st || l_p > target_ed)){
				if(f_out!=NULL)
			{
				print_read(r1,f_out);
				print_read(r2,f_out);
			}
				tg_align_0++;
			}
			if(inner_d>=400)
			{
		 		align_query(K,cc,r1,cc2,1);
                       		align_query(K,cc,r2,cc2,1);
			}
/*			if(f_out!=NULL)
			{
				print_read(r1,f_out);
				print_read(r2,f_out);
			}*/

		}
		if(r1.map_f==2 || r2.map_f==2){
			align_query(K,cc,r1,cc2,1);
			align_query(K,cc,r2,cc2,1);
			if(r1.map_f!=-1&&r2.map_f!=-1){
			make_sam(r1,r2,rev);
			split[0].push_back(r1);
			split[1].push_back(r2);
/*			if(f_out!=NULL)
			{
				print_read(r1,f_out);
				print_read(r2,f_out);
			}*/
			}
		}
	}
	return 1;
}
void process(char *l_fastq, char *r_fastq, vector<read> *split, unordered_map<long long int,vector<int>> &cc,int K,unordered_map<long long int,vector<int>> &cc2)
{
	FILE *f2,*f3;
	read r1,r2;
	int ccc=0;
	vector<string> sam[2];
	double inp_time=0,align_t=0,t1=0,t2=0;
	int total_align=0,total_read=0,one_align=0;
	for(int i=0;i<4;i++){r1.tags.push_back("");r2.tags.push_back("");}
	f2 = fopen((dir+string(l_fastq)).c_str(),"r");
	f3 = fopen((dir+string(r_fastq)).c_str(),"r");
	mm_tot=mm_cnt=0;	
	t2=clock();
	while(read_fastq(f2,r1)&&read_fastq(f3,r2))
	{
		inp_time+=(clock()-t2)/CLOCKS_PER_SEC;
		int is_align=0,is_align2=0,pair_res=0;
		ccc++;
		total_read++;
		t1=clock();
		parse_read(r1);    
		parse_read(r2);
		reverse_comp(r1.query,r1.query_len);
		is_skip=0;
		pair_res=align_pair(K,cc,cc2,r1,r2,split,1);    /// r1: rev   r2: normal
		if(pair_res==0){
			reverse_comp(r2.query,r2.query_len);
			reverse_comp(r1.query,r1.query_len);
			pair_res=align_pair(K,cc,cc2,r1,r2,split,0); // r1 : normal r2: rev
		}
		total_align+=pair_res;
		skip+=is_skip;
		t2=clock();
	}
	printf("mean list %lf\n total_align: %d\n total_read %d\n alignment_rate %.3lf\n",(double)((double)(mm_tot)/(double)(mm_cnt)),total_align,total_read,(double)((double)total_align/(double)total_read));
	t_align = total_align;
	printf("align_t : %lf\ninp_time : %lf\n",align_t/CLOCKS_PER_SEC,inp_time);
	fclose(f2);
	fclose(f3);
}

align_wrap::align_wrap(char *argv[])
{
	vector<read> split_ctr[2],split_mut[2];
	unordered_map<long long int,vector<int>> cc,cc2;	
	FILE *f1;
	int K=11,tt=0,flag=0,i;
	long long int h=0,last_bit=(1LL<<(2*(K-1)));
	if(argv[9]==NULL)
	{
		dir="/data/";
		argv[10] = argv[11] = argv[12] = argv[13] = argv[8];
	}
	else
	{
		K=atoi(argv[9]);
		dir="";
	}
	if(argv[14]==NULL)
	{
		dist_cut=1000;
	}
	else
	{
		dist_cut=atoi(argv[14]);
	}
	double t1 = clock(),t2=clock();
	lll=make_ref(argv,K,cc);  //k-mer reference
	lll=make_ref(argv,17,cc2); //17-mer reference unique assume
	printf("ref time %lf\n",(clock()-t1)/CLOCKS_PER_SEC);
	target_st=atoi(argv[6]),target_ed=atoi(argv[7]),reff=string(argv[8]);
	mm=0;
	t1=clock();
	skip=0;
	ppp2.clear();
	if(argv[19]!=NULL)
	{
		strcpy(view_q,argv[19]);
	}
	if(argv[17]!=NULL)
	{
		f_out=fopen(argv[17],"w");
	}

	tg_align=tg_align_50=tg_align_100=tg_align_150=0;
	tg_align2=0;
	process(argv[2],argv[3],split_ctr,cc,K,cc2);
	fclose(f_out);
	printf("ctr align time : %lf\n",(clock()-t1)/CLOCKS_PER_SEC);
	printf("ctr max chunk : %d\n",mm);
	printf("ctr large del size : %d\n",split_ctr[0].size());

	printf("ctr normal : %d %lf\n",normal,(double)((double)normal/(double)t_align));
	printf("ctr tg_align_0 : %d %lf\n",tg_align_0,(double)((double)tg_align_0/(double)t_align));
	printf("ctr tg_align : %d %lf\n",tg_align,(double)((double)tg_align/(double)t_align));
	printf("ctr skip %d\n",skip);
	FILE *fp;	
	if(argv[15]!=NULL)
	{
		fp =  fopen(argv[15],"w");
		for(i=0;i<ppp2.size();i++)
		{
			fprintf(fp,"ctr %d %d %d %d %d\n",ppp2[i],ppp3[i].x1,ppp3[i].y1,ppp3[i].x2,ppp3[i].y2);
		}
	}
	ppp2.clear();
	ppp3.clear();
	t1 = clock();
	mm=0;
	skip=0;
	
	st_flag=1;
	if(argv[18]!=NULL)
	{
		f_out=fopen(argv[18],"w");
	}

	tg_align=tg_align_50=tg_align_100=tg_align_150=tg_align_0=normal=tg_align4=tg_align5=long_d=tg_align6=tg_align7=tg_align8=0;
	tg_align2=tg_align3=0;
	process(argv[4],argv[5],split_mut,cc,K,cc2);
	fclose(f_out);
	printf("mut align time : %lf\n",(clock()-t1)/CLOCKS_PER_SEC);
	printf("mut max chunk : %d\n",mm);
	printf("mut large del size : %d\n",split_mut[0].size());
	
	printf("mut normal : %d %lf\n",normal,(double)((double)normal/(double)t_align));
	printf("mut tg_align_0 : %d %lf\n",tg_align_0,(double)((double)tg_align_0/(double)t_align));
	printf("mut tg_align : %d %lf\n",tg_align,(double)((double)tg_align/(double)t_align));
	printf("mut skip %d\n",skip);
	if(argv[15]!=NULL)
	{
		for(i=0;i<ppp2.size();i++)
		{
			fprintf(fp,"mut %d %d %d %d %d\n",ppp2[i],ppp3[i].x1,ppp3[i].y1,ppp3[i].x2,ppp3[i].y2);
		}
	}
	t1 = clock();
	post_process(split_ctr,split_mut,argv);
	printf("post time : %lf\n",(clock()-t1)/CLOCKS_PER_SEC);
	printf("whole time : %lf\n",(clock()-t2)/CLOCKS_PER_SEC);
}
align_wrap::~align_wrap(){
}
