#include<stdio.h>
#include<string>
using namespace std;
#include "be_fastq-lib.h"

int main(int argc, char *argv[]) {
//void read_fastq_files(fastq_joined-file) {
    int linecnt;

    struct fq entry;
    meminit(entry);

    FILE* fp;
    fp = fopen(argv[1], "r");
    setvbuf(fp, NULL, _IOFBF, 4194304);
    struct stat st;
    stat(argv[1], &st);
	char *t;
    while (read_fq(fp, linecnt++, &entry,t)) {
        string(entry.seq.s);
	printf("%s\n", entry.seq.s);
    }
    fclose(fp);
    return 0;
}

