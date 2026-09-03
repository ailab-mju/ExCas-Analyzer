#include<string>
#ifndef AG
#define AG
struct Args
{
	char *refname,*ctrfastq_l,*ctrfastq_r,*mutfastq_l,*mutfastq_r,*tmpfile,*repeatfile,*view;
	int targetst,targeted,K,dist_cut,viewer,cle,ref_len;
	std::string dir,resctrcom,resmutcom,resmutonly,resctronly,resctr,resmut,ctrsite,mutsite,ctrsite_sm,mutsite_sm,outdir,genename,samplename;

	/* --- explicit build/run differences (were commented-out code before) --- */
	std::string datadir;    /* prefix for every input path: "" on CLI, "/data/" under WORKERFS */
	int  emit;              /* EMIT_* mask, see sink.h */
	bool use_pre_chunk;     /* pre_chunk() fast path. Canonical runs have it OFF. */
};
#endif
/* returns 0 on success, non-zero on a usage error (message on stderr) */
int Arg_Parse(Args *arg,int argc,char *argv[]);
