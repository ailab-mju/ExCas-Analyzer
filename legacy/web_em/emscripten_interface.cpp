#include <emscripten.h>
#include <string>
#include <vector>
#include <stdlib.h>
#include "align_wrap.h"
using namespace std;
extern "C" {
    void report_result(align_wrap $wt, int wtcnt, char* tt) {
        int i, j=4;
	string a="end";
	EM_ASM_({
		postMessage([0,UTF8ToString($0)]);
		}, a.c_str());
/*        for (i = 0; i < $wt.m_seq.size(); i++) {
            //EM_ASM_({
            //    rtnval['seq_count'].push([Pointer_stringify($0), $1]);
            //}, $wt.m_seq[i].c_str(), $wt.m_seq_count[$wt.m_seq[i]]);
        }
	*/
/*        EM_ASM_({
            rtnval['info'].push([Pointer_stringify($0), $1]);
        }, wtseq.c_str(), wtcnt);*/
//	 EM_ASM_({rtnval['info'].push([$0,Pointer_stringify($1),Pointer_stringify($2)]);
  //      },tt);//$wt.lll,$wt.aaa.c_str(),tt);
    }

    void run_webtest(char *ref_name1,char* file_name1, char* file_name2, char* file_name3, char* file_name4,char *st,char *ed, char* cle, char* genename, char *samplename, char *K, char *outdir,int wtcnt) {
        
	char *t[13];
	t[1] = ref_name1;
	t[2] = file_name1;
	t[3] = file_name2;
	t[4] = file_name3;
	t[5] = file_name4;	
//	strcpy((string("/data/")+string(ref_name1)).c_str(),t[1]);
//	t[2]=(string("/data/")+string(file_name1)).c_str();
//	t[3]=(string("/data/")+string(file_name2)).c_str();
//	t[4]=(string("/data/")+string(file_name3)).c_str();
//	t[5]=(string("/data/")+string(file_name4)).c_str();
	t[6] = st;
	t[7] = ed;
	t[8] = cle;
	t[9] = genename;
	t[10] = samplename;
	t[11] = K;
	t[12] = outdir;
	printf("c++ start\n");
	for(int i=1;i<=12;i++)
	{
		printf("%s\n",t[i]);
	}
	printf("arg end\n");
        align_wrap wt = align_wrap(t);
        report_result(wt,wtcnt,(char *)t[1]);

    }
};
