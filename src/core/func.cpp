#include "func.h"
#include "arg_parse.h"
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <string>
#include <list>
#include <unordered_map>
#include <vector>
#include <string.h>
#include "sink.h"
using namespace std;
extern char g_ref_seq[100000];
extern Args arg;
Sink *g_sink = NULL;
int g_site_kind = -1, g_site_sm_kind = -1;
static void emit_site(read &r1, read &r2, int with_indel)
{
	int k = with_indel ? g_site_sm_kind : g_site_kind;
	if(g_sink==NULL || k<0 || !g_sink->wants(k)) return;
	g_sink->write(k, format_read(r1));
	g_sink->write(k, format_read(r2));
}
void parse_read(read &r)
{
	r.query_len = strlen(r.query);
	char_to_num(r.query,r.query_len);
}



int make_ref(string refname,int K,unordered_map<long long int,vector<int>> &kmer_ref)
{
	char t[200];
	t['A'] = 0,t['C'] = 1,t['G'] = 2,t['T'] = 3;
	FILE *f1;
	int i,l;
	long long int h=0,last_bit=(1LL<<(2*(K-1)));	
	string rpath = arg.datadir+refname;
	f1 = fopen(rpath.c_str(),"r");
	if(f1==NULL){ fprintf(stderr,"cannot open reference: %s\n",rpath.c_str()); return -1; }
	/* a failed open used to fall through and align against an empty reference,
	   producing a plausible-looking run with meaningless numbers */
	if(fscanf(f1,"%s",g_ref_seq)!=1 || fscanf(f1,"%s",g_ref_seq)!=1){
		fprintf(stderr,"cannot read a FASTA record from %s\n",rpath.c_str());
		fclose(f1); return -1;
	}
	fclose(f1);
	l=strlen(g_ref_seq);
	arg.ref_len = l;
	for(i=0;i<=l;i++) 
	{
		if(i>=K)
		{
			kmer_ref[h].push_back(i-K);
			h-=last_bit*t[g_ref_seq[i-K]];
		}
		h*=4;
		h+=t[g_ref_seq[i]];	
	}
	return l;
}
void char_to_num(char *query,int l)
{
	char t[200];
	t['A'] = 0,t['C'] = 1,t['G'] = 2,t['T'] = 3,t['N'] = 3;
	for(int i=0;i<l;i++)
	{
		query[i] = t[query[i]];
	}
}
char* num_to_char(char *query,int l)
{
	char t[5];
	t[0] = 'A',t[1] = 'C',t[2] = 'G',t[3] = 'T';
	for(int i=0;i<l;i++)
	{
		query[i] = t[query[i]];
	}
	return query;
}
void reverse_comp(char *query,int l)
{
	char query2[300];
	char oper = 252;
	for(int i=0;i<l;i++)
	{
		query2[l-i-1] = (~query[i])^oper;
	}
	for(int i=0;i<l;i++)
	{
		query[i] = query2[i];
	}
}

std::string format_read(read &r)
{
	char buf[4096];
	int n = snprintf(buf,sizeof(buf),"%s\t%d\t%s\t%d\t%d\t",r.query_name+1,r.sam_f,r.ref.c_str(),r.pos,r.map_q);
	for(int i=0;i<r.cigar_string.size();i++)
		n += snprintf(buf+n,sizeof(buf)-n,"%d%c",r.cigar_string[i].len,r.cigar_string[i].type);
	n += snprintf(buf+n,sizeof(buf)-n,"\t=\t%d\t%d\t%s\t%s",r.p_next,r.tlen, num_to_char(r.query,r.query_len),r.query_qual);
	char_to_num(r.query,r.query_len);
	for(int i=0;i<r.tags.size();i++)
		n += snprintf(buf+n,sizeof(buf)-n,"\t%s",r.tags[i].c_str());
	return std::string(buf);
}
void make_sam(read &r1, read &r2,int flag,string gene_name)
{
	//printf("make_sam\n"); //log
	r1.sam_f=r2.sam_f=3;
	r1.sam_f+=64;
	r2.sam_f+=128;
	r1.sam_f+=16+(!flag*16);
	r2.sam_f+=16+((flag)*16);
	r1.map_q=60;
	r2.map_q=60;
	r1.pos = r1.ref_st;
	r2.pos = r2.ref_st;
	r1.p_next = r2.pos;
	r2.p_next = r1.pos;
	r1.dist = min(r1.dist,r2.dist);
	r2.dist = min(r1.dist,r2.dist);
	if(r1.pos>=r2.pos)
	{
		r1.tlen = -(abs(r1.pos-r2.pos)+r1.query_len);
		r2.tlen = (abs(r1.pos-r2.pos)+r1.query_len);
	}
	else
	{
		r1.tlen = (abs(r1.pos-r2.pos)+r2.query_len);
		r2.tlen = -(abs(r1.pos-r2.pos)+r2.query_len);
	}
	r1.ref=(gene_name);
	r2.ref=(gene_name);
	r1.del_size = max((r1.ref_ed-r1.ref_st),(r2.ref_ed-r2.ref_st));
	r2.del_size = max((r1.ref_ed-r1.ref_st),(r2.ref_ed-r2.ref_st));
	r1.tags[2]="Dist:"+to_string(r1.dist);
	r2.tags[2]="Dist:"+to_string(r2.dist);
	r1.tags[3]="D_size:"+to_string(r1.del_size);
	r2.tags[3]="D_size:"+to_string(r2.del_size);
}
string rgb(string hex)
{
	char t[20];
	strcpy(t,hex.c_str());
	for(int i=1;i<=6;i++)
	{
		if(t[i]<='9')
		{
			t[i]-='0';
		}
		else
		{
			t[i]-='A';
			t[i]+=10;
		}
	}
	int r = t[2]+t[1]*16;
	int g = t[4]+t[3]*16;
	int b = t[6]+t[5]*16;
	string res=to_string(r)+","+to_string(g)+","+to_string(b);
	return res;
}
int long_inner_check(read &r1,read &r2,vector<read> *split,result &res)
{

	//printf("long_inner_check\n"); //log
	int inner_d=0,l_p=0,r_p=0,l_p2=0,r_p2=0;
	if(r1.ref_st<r2.ref_st)
	{
		inner_d=r2.ref_ed-r1.ref_st;
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
	}
	//printf("%d %d %d %d\n",l_p,l_p2,r_p,r_p2);
	if( (r_p2-l_p2)>400) 
	{
		//printf("%d\n",r_p2-l_p2);
		//printf("%d\n",(r_p2-l_p2)>th);
		if(r_p2 < arg.targetst || arg.targeted < l_p2){
			r1.overlap = 0;
			r2.overlap = 0;
			r1.map_f=r2.map_f=-1;
			//printf("%d %d %d %d inner1\n",r_p2,l_p2,arg.targetst,arg.targeted);
			return -1;
		}
		int ddist=10000;
		ddist = min(ddist,abs(r1.ref_st-arg.targetst-arg.cle));
		ddist = min(ddist,abs(r1.ref_ed-arg.targetst-arg.cle));
		ddist = min(ddist,abs(r2.ref_st-arg.targetst-arg.cle));
		ddist = min(ddist,abs(r2.ref_ed-arg.targetst-arg.cle));
		r1.map_f=2;r2.map_f=2;
		r1.ref_st=r2.ref_st=l_p;r1.ref_ed=r2.ref_ed=r_p;
		r1.del_size=r2.del_size=r_p-l_p;
		r1.dist=r2.dist=ddist;
		if(ddist>arg.dist_cut){r1.map_f=r2.map_f=-1;
		//printf("%d %d inner2\n", ddist, arg.dist_cut);
		return -1;}
		else{
			r1.tags[2]="Dist:"+to_string(r1.dist);
			r2.tags[2]="Dist:"+to_string(r2.dist);
			r1.tags[3]="D_size:"+to_string(r1.del_size);
			r2.tags[3]="D_size:"+to_string(r2.del_size);
			split[0].push_back(r1);
			split[1].push_back(r2);
			res.long_d++;
			return inner_d;
		}
	}
	res.normal++;
//small indel			
	if((l_p<arg.targetst+arg.cle-13 && arg.targetst+arg.cle+13 < r_p)){
				res.tg_align_0++;
			if((r1.ins|r1.del|r1.miss)|(r2.ins|r2.del|r2.miss))
			{
				if(r1.ins|r2.ins)
				{
					res.tg_ins++;
				}
				else if(r1.del|r2.del)
				{
					res.tg_del++;
				}
				else
				{
					res.tg_miss++;
				}
				emit_site(r1,r2,1);
	
			}
			else
			{
				emit_site(r1,r2,0);
			}

	}
	if(!(r_p<arg.targetst || l_p > arg.targeted)){
				res.tg_align_1++;
			if((r1.ins|r1.del|r1.miss)|(r2.ins|r2.del|r2.miss))
			{
				if(r1.ins|r2.ins)
				{
					res.tg_ins1++;
				}
				else if(r1.del|r2.del)
				{
					res.tg_del1++;
				}
				else
				{
					res.tg_miss1++;
				}
			}

	}
	if((l_p<arg.targetst && arg.targeted < r_p)){
		res.tg_align_2++;
			if((r1.ins|r1.del|r1.miss)|(r2.ins|r2.del|r2.miss))
			{
				if(r1.ins|r2.ins)
				{
					res.tg_ins2++;
				}
				else if(r1.del|r2.del)
				{
					res.tg_del2++;
				}
				else
				{
					res.tg_miss2++;
				}
			}
		}
	return inner_d;
		//printf("long_inner_end\n"); //log
}
