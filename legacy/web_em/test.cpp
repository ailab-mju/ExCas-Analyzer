#include<stdio.h>
#include<string.h>
#include<time.h>
int main(int argc, char *argv[])
{
	char *b;
	char a[10000];
	FILE *in = fopen(argv[1],"r");
	size_t l=0;

		int ll;
	double t1=clock();
	size_t tmp;
	while(getline(&b,&tmp,in))
	//while(fscanf(in,"%s",a)!=EOF)
	//while(fgets(a,5000,in)!=NULL)
	{
		printf("%s\n",b);
	}
	printf("%lf\n",(clock()-t1)/CLOCKS_PER_SEC);
	return 0;
}
