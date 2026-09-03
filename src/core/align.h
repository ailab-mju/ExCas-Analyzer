#include<stdio.h>
#include<list>
#include<unordered_map>
#include<algorithm>
#include<string>
#include<string.h>
#include<vector>
#include "func.h"
#include<stdlib.h>
using namespace std;

inline void check_concat(list<data>::iterator &iter, list<data> &be, list<data> &af, data tmp,int &is_find,int K);
void read_classify(read &r, vector<data> &res);
void align_query(int K, unordered_map<long long int,vector<int>> &kmer_ref, read &r,unordered_map<long long int,vector<int>> &kmer_ref2,int flag2);
void print_result(vector<data> &res, vector<cigar> &cigar_string, list<data> &chunk, read &r);
