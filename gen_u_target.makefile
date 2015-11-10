hh_file=
cc_file=
o_file=
CXXFLAGS=
LIBS=
BAMINC=-I/home/chzelin/bioprogram/samtools-0.1.19/include 
BAMLIB=-lbam -lz -L/home/chzelin/bioprogram/samtools-0.1.19/lib 
BOOSTLIB=-lboost_system -lboost_filesystem -lboost_regex
all: gen_u_target
gen_u_target: gen_u_target.cc
	g++ -g -o gen_u_target gen_u_target.cc -DPTHREAD ${BAMINC} ${BAMLIB} -lm ${BOOSTLIB} ${LIBS} -lbam -pthread -lz -lgzstream -lczl_bio2
