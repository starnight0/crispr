/**
 * Copyright 2015, Zelin Chen <chzelin@gmail.com>
 *
 * @file filt_pam.cpp
 * @brief filter PAM-not-exact-match alignment
 *
 * @author Zelin Chen
 * 
 */
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include "czl_bio/czl_common.h"
#include "czl_bio/czl_io.h"
#include <boost/filesystem.hpp>
#include "bam.h"
#include "faidx.h"

using namespace std;
using namespace boost;
using namespace czl_bio;

void print_head();

string prog_name;
string prog_version="1.0";
string prog_author="Zelin Chen";
string prog_compiled_host;
string prog_function="filter PAM-not-exact-match alignment. \
BAM input must be sorted by position";
string prog_input="";
string prog_output="";

int bam_seq_to_str(int l_seq, uint8_t * qseq, string & seq);

int main(int argc, char* argv[])
{
    prog_name = argv[0];
    print_head();

    /// main global;
    // {{{
    string out_dir;
	time_t rawtime;
    // }}}
    ///

    /// temparory
    // {{{
    int i, j, k, l, m, n;
    string line, str, s;
    stringstream ss;
    string file, in_file, out_file;
    fstream fs;
    ifstream fin;
    ofstream fout;
    // }}}
    ///

    /// for option
    // {{{
    Opts opts;
    opts.get_opt_by_name<OptS>("version")->get_value_ref()=prog_version;
    OptS *opt_out_prefix = opts.create<OptS>("out_prefix,out,o", "STR", "prefix of any output", "out.", 1, true);
    string & out_prefix = opt_out_prefix->get_value_ref();
    string & tmp_dir = opts.create<OptS>("tmp_dir,tmp,t", "STR", "template dir", "", 1, false)->get_value_ref();
    string & log_file = opts.create<OptS>("log_file,log,l", "STR", "log file", "", 1, false)->get_value_ref();
    bool & is_clean = opts.create<OptB>("is_clean", "BOOL", "whether to clean temporary files after finished, T/True: clean, F/False: not clean.", "T", 1, false)->get_value_ref();

    string & in_fasta = opts.create<OptS>("in_fasta", "FILE",
			"input reference equence fasta file", "", 1, true)->get_value_ref();

    string & in_bam = opts.create<OptS>("in_bam", "FILE",
			"input bam file", "", 1, true)->get_value_ref();

    string & pam = opts.create<OptS>("pam", "STR",
            "length of the PAM", "NGG", 1, false)->get_value_ref();

    char & pam_pos = opts.create<OptC>("pam_pos", "STR",
            "5' or 3' PAM, 5 or [3]", "5", 1, false)->get_value_ref();

    opts.parse(argc, argv);

    if ( !(pam_pos=='5' || pam_pos=='3') ) {
        cerr << "E: 'pam_pos' must be '5' or '3'" << CZL_DBG_INFO << endl;
        #ifdef PTHREAD
        pthread_exit( (void*)1 );
        #else
        return 1;
        #endif
    }

    if (tmp_dir=="") {
        tmp_dir = out_prefix+"tmp";
    }
    {
        boost::filesystem::path p(tmp_dir);
        if (!boost::filesystem::exists(p)) {
            boost::filesystem::create_directory(p);
            if (!boost::filesystem::exists(p)) {
                cerr << "Can't create TMP_DIR " << tmp_dir << endl;
                #ifdef PTHREAD
                pthread_exit( (void*)1 );
                #else
                exit(1);
                #endif
            }
        }
    }

    if (log_file=="") {
        log_file = out_prefix+"log";
    }

    {
        i = out_prefix.find_last_of("/");
        if ( i==string::npos ) {
            out_dir = "./";
        } else {
            out_dir = out_prefix.substr(0, i);
        }
    //  boost::filesystem::path p(out_dir);
        if (!boost::filesystem::exists(out_dir)) {
            cerr << "Can't find OUT_DIR " << out_dir << ". You need to create it manually" << endl;
            exit(0);
        }
    }
    /// }}}

    Msg log(log_file);
	ss.str("");
	ss<< "Begin at: " << ctime(&rawtime) << "\t" << clock()/CLOCKS_PER_SEC;
	cout<< ss.str() << endl;
	log << ss.str() << endl;

	int pam_len = pam.size();
	uint8_t pam_n16[pam_len/2+1];
	// change from char sequence to BAM n16 sequence (4-bit)
	for (i=0; i<pam_len; i++) {
		bam1_seq_seti(pam_n16, i, bam_nt16_table[pam[i]]);
	}
	string c_pam;
	complement_copy(pam, c_pam);
	uint8_t c_pam_n16[pam_len/2+1];
	// change from char sequence to BAM n16 sequence (4-bit)
	for (i=0; i<pam_len; i++) {
		bam1_seq_seti(c_pam_n16, i, bam_nt16_table[c_pam[i]]);
	}

	fin.open(in_fasta.c_str());
    if (fin.fail()) {
        cerr << "E: Fail to open file " << in_fasta << CZL_DBG_INFO << endl;
		#ifdef PTHREAD
		pthread_exit( (void*)(-1) );
		#else
        return -1;
		#endif
	}
	fin.close();

	string faidx_file = in_fasta + ".fai";
	fin.open(faidx_file.c_str());
    if (fin.fail()) {
		if ( fai_build(in_fasta.c_str()) ) {
			cerr<< "E: Fail to build fasta index file " << faidx_file
				<< CZL_DBG_INFO << endl;
			return -2;
		}
	}
	fin.close();

	faidx_t *fai = fai_load(in_fasta.c_str());
	if (fai==NULL) {
		cerr<< "E: Fail to load fasta index for " << in_fasta
			<< CZL_DBG_INFO << endl;
		return -3;
	}
	
    bamFile fp=NULL;
    fp = bam_open(in_bam.c_str(), "r");
    if (fp==NULL) {
        cerr << "E: Can't open BAM file " << in_bam << CZL_DBG_INFO << endl;
        return -1;
    }
	
	string out_bam = out_prefix+"filtered.bam";
	bamFile fwp = bam_open(out_bam.c_str(), "w");
    if (fwp==NULL) {
        cerr << "E: Can't create BAM file " << out_bam << endl;
        return -1;
    }

	// read and write header
    bam_header_t *header = bam_header_read(fp);
	bam_header_write(fwp, header);

    bam1_t *bam, *prev_bam;
    bam = bam_init1();
    prev_bam = bam_init1();
    vector<bam1_t*> bam_v;
    const int N=5;
    int count=0;
    int r=0;
    int ns[N]; // 0,1,2,3: 0 to 3 mismatch, 4: with gap 
	char *tseq=NULL;
	int tlen=0;
    for (i=0; i<N; i++) ns[i]=0;
    while ( bam_read1(fp, bam)>0 ) {
        if (count==0) {
			// load new sequence from fasta file
			tseq = fai_fetch(fai, header->target_name[bam->core.tid], &tlen);
		} else if ( bam->core.tid!=prev_bam->core.tid ) {
			free(tseq);
			tseq = fai_fetch(fai, header->target_name[bam->core.tid], &tlen);
        }
		if (tseq==NULL) {
			cerr << "E: Fail to get sequence " << header->target_name[bam->core.tid] 
				 << CZL_DBG_INFO << endl;
			return 4;
		}
		count++;

		char *qname = bam1_qname(bam);
		int n_cigar = bam->core.n_cigar;
		if (n_cigar==0) { continue; }
		uint32_t * cigar = bam1_cigar(bam);
		int l_qseq = bam->core.l_qseq;
		uint8_t * qseq = bam1_seq(bam);
		int m=0; // mismatch
		int qpos=0, tpos=bam->core.pos;
		int m0 = 0; // mismatch in PAM
		int tpos1 = bam_calend(&bam->core, cigar);

		int nm;
		nm = bam_aux2i(bam_aux_get(bam, "NM"));
		short is_f = 0; // is filered
		if (nm>0) {
			if ( pam_pos == '5') {
				if ( bam1_strand(bam) ) {
					// minus strand
					for (j=0; j<pam_len; j++) {
						if ( bam_nt16_table[ tseq[tpos1-1-j] ] & bam1_seqi(c_pam_n16,j) ) {
						} else {
							break;
						}
					}
					if (j >= pam_len) {
						for (j=0; j<pam_len; j++) {
							if (bam_nt16_table[ tseq[tpos1-pam_len+j] ] != bam1_seqi(qseq,l_qseq-pam_len+j)) m0++;
						}
					} else {
						is_f = 1;
					}
				} else {
					// plus strand
					for (j=0; j<pam_len; j++) {
						if ( bam_nt16_table[ tseq[tpos+j] ] & bam1_seqi(pam_n16,j) ) {
						} else {
							break;
						}
					}
					if (j >= pam_len) {
						for (j=0; j<pam_len; j++) {
							if (bam_nt16_table[ tseq[tpos+j] ] != bam1_seqi(qseq,j)) m0++;
						}
					} else {
						is_f = 1;
					}
				}
			} else {
				if ( bam1_strand(bam) ) {
					// minus strand
					for (j=0; j<pam_len; j++) {
						if ( bam_nt16_table[ tseq[tpos+pam_len-1-j] ] & bam1_seqi(c_pam_n16,j) ) {
						} else {
							break;
						}
					}
					if (j >= pam_len) {
						for (j=0; j<pam_len; j++) {
							if (bam_nt16_table[ tseq[tpos+j] ] != bam1_seqi(qseq,j)) m0++;
						}
					} else {
						is_f = 1;
					}
				} else {
					// plus strand
					for (j=0; j<pam_len; j++) {
						if ( bam_nt16_table[ tseq[tpos1-pam_len+j] ] & bam1_seqi(pam_n16,j) ) {
						} else {
							break;
						}
					}
					if (j >= pam_len) {
						for (j=0; j<pam_len; j++) {
							if (bam_nt16_table[ tseq[tpos1-pam_len+j] ] != bam1_seqi(qseq,l_qseq-pam_len+j)) m0++;
						}
					} else {
						is_f = 1;
					}
				}
			}
		}
		if (!is_f) {
			int xd = nm-m0;
			bam_aux_append(bam, "XD", 'i', sizeof(xd), (uint8_t*)&xd);
			bam_write1(fwp, bam);
			bam1_t *t = prev_bam;
			prev_bam = bam;
			bam = t;
		}
    }

    bam_destroy1(bam);
    bam_destroy1(prev_bam);
	fai_destroy(fai);
    bam_header_destroy(header);
    bam_close(fp);
    bam_close(fwp);

    if (is_clean) {
        // clean tmp files, BE CAREFULL
    }

	ss.str("");
	ss<< "End at: " << ctime(&rawtime) << "\t" << clock()/CLOCKS_PER_SEC;
	cout<< ss.str() << endl;
	log << ss.str() << endl;

#ifdef PTHREAD
    pthread_exit(NULL);
#else
    return 0;
#endif
}

void print_head()
{
    cout << "**************************************************" << endl;
    cout << "* Program  : " << prog_name << endl;
    cout << "* Version  : " << prog_version << endl;
    cout << "* Author   : " << prog_author << endl;
    cout << "* Function : \n";
    StringUtility::print(prog_function, 2, 80, cout);
    cout << "\n";
    cout << "* Input    : " << prog_input << endl;
    cout << "* Output   : " << prog_output << endl;
    cout << "* Compiled Host: " << prog_compiled_host << endl;
    cout << "* Compiled Time: " << __DATE__ << " : " __TIME__ << endl;
    cout << "* Machine Bit: " << machine_bit() << endl;
    cout << "* IS big endian: " << (is_be()?"yes":"no") << endl;
    cout << "**************************************************" << endl;
}

void print_head_d()
{
    ifstream fin("head.txt");
    char str[1025];
    char c;
    c = fin.get();
    while ( !fin.eof() ) {
        std::cout << c;
        c = fin.get();
    }
    if (c != EOF) {
        cout << c;
    }

    fin.close();
    return;
}

int bam_seq_to_str(int l_qseq, uint8_t * qseq, string & seq)
{
    seq.resize(l_qseq);
    for (int i=0; i<l_qseq; i++) {
        seq[i] = bam_nt16_rev_table[bam1_seqi(qseq, i)];
    }
}
