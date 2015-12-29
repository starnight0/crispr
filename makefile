hh_file=
cc_file=
o_file=
CXXFLAGS=-I$(HOME)/include -I../czl_bio
CXXLIBS=-L$(HOME)/lib -lbam -pthread -lz -lczl_bio -lhts
CPPFLAGS=-O3
CPPLIBS=
#BOOSTLIB=-lboost_system -lboost_filesystem -lboost_regex
PREFIX=$(HOME)
all: search_pam index32 parse_web_input
parse_web_input: parse_web_input.cpp
	g++ parse_web_input.cpp $(CXXFLAGS) $(CPPFLAGS) -o $@ $(CXXLIBS) $(CPPLIBS)  -lm 
search_pam: search_pam.cpp bit2.cpp bit2.hpp
	g++ search_pam.cpp bit2.cpp $(CXXFLAGS) $(CPPFLAGS) -o $@ $(CXXLIBS) $(CPPLIBS)  -lm 
filt_pam: filt_pam.cpp
	g++ $< $(CXXFLAGS) $(CPPFLAGS) -o $@ $(CXXLIBS) $(CPPLIBS)  -lm 
target_stat: target_stat.cpp
	g++ -g $< -o $@ -I../ ${INCS} -lm ${LIBS} -lbam -pthread -lz -lczl_bio
index32: index32.cpp bit2.hpp bit2.cpp
	g++ index32.cpp bit2.cpp $(CXXFLAGS) $(CPPFLAGS) -o $@ $(CXXLIBS) $(CPPLIBS)  -lm 
install:
	cp search_pam $(PREFIX)/bin
	cp filt_pam $(PREFIX)/bin
	cp target_stat $(PREFIX)/bin
	cp *.sh $(PREFIX)/bin
