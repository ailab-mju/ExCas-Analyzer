#include<stdlib.h>
#include<stdio.h>
#include<string>
#include<string.h>
#include "arg_parse.h"
#include "sink.h"
using namespace std;

/* Positional arguments are unchanged from the original tool:
     1 ref.fa  2 ctr_1 3 ctr_2  4 mut_1 5 mut_2  6 tg_st 7 tg_ed 8 cle
     9 gene   10 sample  11 K   12 outdir  [13 dist_cut]
   Optional flags may appear anywhere after those:
     --emit=all|none|split[,full][,site]   --pre-chunk   --datadir=PATH          */
int Arg_Parse(Args *arg,int argc,char *argv[])
{
	arg->datadir       = "";
	arg->emit          = EMIT_ALL;
	arg->use_pre_chunk = false;
	arg->dist_cut      = 1000;

	int pos = 0;
	char *p[14];
	for(int i=0;i<14;i++) p[i]=NULL;
	for(int i=1;i<argc;i++)
	{
		if(strncmp(argv[i],"--",2)==0)
		{
			if(strncmp(argv[i],"--emit=",7)==0)
			{
				int m = emit_parse(argv[i]+7);
				if(m<0){fprintf(stderr,"bad --emit value: %s\n",argv[i]+7);return 1;}
				arg->emit = m;
			}
			else if(strncmp(argv[i],"--datadir=",10)==0) arg->datadir = string(argv[i]+10);
			else if(strcmp(argv[i],"--pre-chunk")==0)    arg->use_pre_chunk = true;
			else {fprintf(stderr,"unknown option: %s\n",argv[i]);return 1;}
		}
		else if(pos<13) p[++pos]=argv[i];
		else {fprintf(stderr,"too many arguments\n");return 1;}
	}
	if(pos<12){fprintf(stderr,"missing arguments (%d of 12 positional given)\n",pos);return 1;}

	arg->refname     = p[1];
	arg->ctrfastq_l  = p[2];
	arg->ctrfastq_r  = p[3];
	arg->mutfastq_l  = p[4];
	arg->mutfastq_r  = p[5];
	arg->targetst    = atoi(p[6]);
	arg->targeted    = atoi(p[7]);
	arg->cle         = atoi(p[8]);
	arg->genename    = string(p[9]);
	arg->samplename  = string(p[10]);
	arg->K           = atoi(p[11]);
	arg->outdir      = string(p[12]);
	if(p[13]!=NULL) arg->dist_cut = atoi(p[13]);

	string base = arg->outdir+"/"+arg->genename+"_"+arg->samplename+"_";
	arg->resctronly = base+"ctr_only.sam";
	arg->resctrcom  = base+"ctr_common.sam";
	arg->resmutonly = base+"mut_only.sam";
	arg->resmutcom  = base+"mut_common.sam";
	arg->resctr     = base+"ctr.sam";
	arg->resmut     = base+"mut.sam";
	arg->ctrsite    = base+"site_ctr.sam";
	arg->mutsite    = base+"site_mut.sam";
	arg->ctrsite_sm = base+"site_ctr_sm.sam";
	arg->mutsite_sm = base+"site_mut_sm.sam";
	return 0;
}
