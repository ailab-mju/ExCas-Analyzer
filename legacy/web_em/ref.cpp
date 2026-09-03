
#include<unordered_map>
#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<string.h>
char ref_seq[100000];
using namespace std;

int make_ref(char *argv[],int K,unordered_map<long long int,vector<int>> &cc)
{
	char t[200];
	t['A'] = 0,t['C'] = 1,t['G'] = 2,t['T'] = 3;
	FILE *f1;
	int i,l;
	long long int h=0,last_bit=(1<<(2*(K-1)));
	last_bit=1;
	for(i=1;i<=2*(K-1);i++)
	{
		last_bit*=2;
	}
	printf("%lld\n",last_bit);	
//	f1 = fopen("/data/st3gal4.fa","r");
//	f1 = fopen((string("/data/")+string(argv[1])).c_str(),"r");
	f1 = fopen((string(argv[1])).c_str(),"r");

	fscanf(f1,"%s",ref_seq);	
	fscanf(f1,"%s",ref_seq);

	l=strlen(ref_seq);
	for(i=0;i<=l;i++) 
	{
		if(i>=K)
		{
			cc[h].push_back(i-K);
			h-=last_bit*(long long int)(t[ref_seq[i-K]]);
		}
		h*=4;
		h+=(long long int)(t[ref_seq[i]]);
	}
	return l;
}
int main(int argc, char *argv[])
{
	unordered_map< long long int , vector<int> > cc;
	int K = atoi(argv[2]);
	make_ref(argv,K,cc);
	for (auto e : cc)
	{
		printf("%lld %d\n",(e.first),(e.second).size());
	}

}

