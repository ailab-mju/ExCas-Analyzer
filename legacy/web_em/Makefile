#CPP=/opt/rh/devtoolset-7/root/usr/bin/g++
CPP=/usr/bin/g++
#kmer_align2: align.o main2.o
#	$(CPP)   -g -pthread -o kmer_align2 align.o main2.o
#main2.o: align.h main2.cpp
#	$(CPP)   -g -pthread -c main2.cpp
#align.o: align.h align.cpp
#	$(CPP)   -g -pthread -c align.cpp

kmer_align: align.o align_wrap.o arg_parse.o func.o post.o kmer_align.o 
	$(CPP) -O3 -g -pthread -o kmer_align align.o align_wrap.o arg_parse.o func.o post.o kmer_align.o
kmer_align.o: align_wrap.h kmer_align.cpp
	$(CPP) -O3 -g -pthread -c kmer_align.cpp
align_wrap.o: func.h align.h align_wrap.cpp
	$(CPP) -O3 -g -pthread -c align_wrap.cpp
align.o: align.h func.h align.cpp
	$(CPP) -O3 -g -pthread -c align.cpp
arg_parse.o: arg_parse.h arg_parse.cpp
	$(CPP) -O3 -g -pthread -c arg_parse.cpp
post.o: post.h post.cpp
	$(CPP) -O3 -g -pthread -c post.cpp
func.o: func.h func.cpp
	$(CPP) -O3 -g -pthread -c func.cpp

#kmer_view: align.o align_wrap.o kmer_align.o
#	$(CPP) -O -g -pthread -o kmer_view align.o align_wrap.o kmer_align.o
#kmer_align.o: align_wrap.h kmer_align.cpp
#	$(CPP) -O -g -pthread -c kmer_align.cpp
#align_wrap.o: align.h align_wrap.cpp
#	$(CPP) -O -g -pthread -c align_wrap.cpp
#align.o: align.h align.cpp
#	$(CPP) -O  -g -pthread -c align.cpp

