hh_file=
cc_file=
o_file=
CXXFLAGS=
LIBS=
BAMINC=-I/home/chzelin/bioprogram/samtools-0.1.19/include 
BAMLIB=-lbam -lz -L/home/chzelin/bioprogram/samtools-0.1.19/lib 
BOOSTLIB=-lboost_system -lboost_filesystem -lboost_regex
all: search_pam.v2
search_pam.v2: search_pam.v2.cc
	g++ -g -o search_pam.v2 search_pam.v2.cc -DPTHREAD ${BAMINC} ${BAMLIB} -lm ${BOOSTLIB} ${LIBS} -lbam -pthread -lz -lgzstream -lczl_bio2
