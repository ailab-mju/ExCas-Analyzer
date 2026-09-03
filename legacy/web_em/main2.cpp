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
#include "dkm.hpp"
int target_st,target_ed;
string reff;
string dir;
int mm=-1;
char ref_seq[100000];
FILE *f4;
void read_cluster(vector<read> *split_ctr, vector<read> *split_mut)
{
	bool mutual1[500],mutual2[500];
	std::tuple<std::vector<std::array<float, 2>>, std::vector<uint32_t>> cluster_data1 = cluster(split_ctr),cluster_data2 = cluster(split_mut);
	vector<std::array<float,2>> center1 = get<0>(cluster_data1),center2 = get<0>(cluster_data2);
	printf("%d %d\n",center1.size(),center2.size());
	vector<uint32_t> label1 = get<1>(cluster_data1),label2 = get<1>(cluster_data2);
	for(int i=0;i<center1.size();i++){mutual1[i]=0;}
	for(int i=0;i<center2.size();i++){mutual2[i]=0;}
	for(int i=0;i<center1.size();i++)
	{
		for(int j=0;j<center2.size();j++)
		{
			float interval = min(abs(center1[i][0]-center1[i][1]),abs(center2[j][0]-center2[j][1]));
			if(abs(center1[i][0]-center2[j][0])<interval*0.1 && abs(center1[i][1]-center2[j][1])<interval*0.1)
			{
				mutual1[i]=1;
				mutual2[j]=1;
			}
		}
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
	printf("%d\n",c);
	fclose(f[0]);
	fclose(f[1]);	
}
void post_process(vector<read> *split_ctr, vector<read> *split_mut,char *argv[])
{
	read_cluster(split_ctr,split_mut);
	printf("ctr:");
	write_read(split_ctr,string(argv[10]),string(argv[11]));
	printf("mut:");
	write_read(split_mut,string(argv[12]),string(argv[13]));
}
void process(char *l_fastq, char *r_fastq, vector<read> *split, unordered_map<long long int,vector<int>> &cc,int K)
{
	FILE *f2,*f3;
	read r1,r2;
	int ccc=0;
	r1.tags.push_back("");
	r2.tags.push_back("");
	r1.tags.push_back("");
	r1.tags.push_back("");
	r2.tags.push_back("");
	r2.tags.push_back("");
	r1.tags.push_back("");
	r2.tags.push_back("");

	f2 = fopen((dir+string(l_fastq)).c_str(),"r");
	f3 = fopen((dir+string(r_fastq)).c_str(),"r");
	printf("%s\n%s\n",l_fastq,r_fastq);
	vector<string> sam[2];
	char dummy[1000];
	int total_align=0;
	while(fscanf(f2,"%s %s %s %c %s",r1.query_name,dummy,r1.query,&r1.strand,r1.query_qual)!=EOF&&fscanf(f3,"%s %s %s %c %s",r2.query_name,dummy,r2.query,&r2.strand,r2.query_qual)!=EOF)
	{
		int is_align=0;
		ccc++;
		parse_read(r1);
		parse_read(r2);
		reverse_comp(r1.query,r1.query_len);
		align_query(K,cc,r1);
		align_query(K,cc,r2);
		if(r1.map_f!=-1 && r2.map_f!=-1){is_align=1;}
		if(r1.map_f!=-1 && r2.map_f!=-1)
		{
			make_sam(r1,r2,0);
			if(r1.map_f==2 || r2.map_f==2){
				split[0].push_back(r1);
				split[1].push_back(r2);
			}
			continue;
		}
		reverse_comp(r2.query,r2.query_len);
		reverse_comp(r1.query,r1.query_len);
		align_query(K,cc,r1);
		align_query(K,cc,r2);

		if(r1.map_f!=-1 && r2.map_f!=-1){is_align=1;}
		if(r1.map_f!=-1 && r2.map_f!=-1){
			make_sam(r1,r2,1);
			if(r1.map_f==2 || r2.map_f==2){
				split[0].push_back(r1);
				split[1].push_back(r2);
			}
		}
		total_align+=is_align;
	}
	printf("total_align:%d\n",total_align);
	fclose(f2);
	fclose(f3);
}

int main(int argc, char *argv[])
{
	vector<array<float,2>> point;
	vector<read> split_ctr[2],split_mut[2];
	vector<string> f_name;
	unordered_map<long long int,vector<int>> cc;	
	FILE *f1;
	int K=11,tt=0,flag=0,i;
	long long int h=0,last_bit=(1<<(2*(K-1)));
	if(argv[9]==NULL)
	{
		printf("asdf");
		dir="/data/";
		argv[10] = argv[11] = argv[12] = argv[13] = argv[8];
	}
	else
	{
		K=atoi(argv[9]);
		printf("%d fdsa",K);
		dir="";
	}
	make_ref(argv,K,cc);
	target_st=atoi(argv[6]);
	target_ed=atoi(argv[7]);
	printf("%d %d\n",target_st,target_ed);
	reff=string(argv[8]);
/*	for(i=10;i<=13;i++){
		f_name.push_back(string(argv[i]));
	}*/
	double t1 = clock();
	mm=0;
	process(argv[2],argv[3],split_ctr,cc,K);
//	std::thread th1( process,argv[2],argv[3],split_ctr,cc,K );
//	thread t2( process,argv[4],argv[5],split_mut,cc,K );
//	th1.join();
//	t2.join();
	printf("%lf\n",(clock()-t1)/CLOCKS_PER_SEC);
	printf("%d\n",mm);

	printf("%d\n",split_ctr[0].size());
	t1 = clock();
	mm=0;
	process(argv[4],argv[5],split_mut,cc,K);
	printf("%lf\n",(clock()-t1)/CLOCKS_PER_SEC);
	printf("%d\n",mm);
	printf("%d\n",split_mut[0].size());
	post_process(split_ctr,split_mut,argv);
	return 0;
}
