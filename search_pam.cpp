/**
 * Copyright 2015, Zelin Chen <chzelin@gmail.com>
 *
 * @file cause.cc
 * @brief 
 *
 * @author Zelin Chen
 * 
 */
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <sys/stat.h>
#include "czl_common.hpp"
#include "czl_io.hpp"
#include "czl_bio_seq.hpp"
#include <ctime>
#include "bit2.hpp"
//#include <boost/filesystem.hpp>

using namespace std;
//using namespace boost;
using namespace czl_bio;

void print_usage();
void print_head();
void print16(uint8_t l, uint16_t seed);
uint8_t geti16(uint16_t seed, uint8_t i);
void seti16(uint16_t * seed, uint8_t i, uint16_t a);
void seti16_next(uint16_t * seed, uint8_t i);
uint16_t rev16(uint8_t l, uint16_t seed);
int fetch_index(string const & seq, int32_t fpos_v[0x10000], int32_t fcount_v[0x10000], ifstream & idx24, int stat[4]);
void sub_fetch(uint8_t bl0, uint64_t bseq0, ifstream & idx24, int32_t n, int k, int *stat);

string prog_name="main";
string prog_version="1.0";
string prog_author="Zelin Chen";
string prog_compiled_host;
string prog_function="";
string prog_input="";
string prog_output="";

const int K = 32;
uint8_t mis_v[0x10000];
/*
class OptInfo {
public:
    string name_t;
    string value_t;
    int n_t;
    
    Opt(string name, string value_t, int n): name_t(name), value_t(value), n_t(n) {}
    Opt(const Opt & opt): name_t(opt.name), value_t(opt.value), n_t(opt.n) {}
    Opt & operator = (const Opt & opt) {
        name_t = opt.name_t;
        value_t = opt.value_t;
        n_t = opt.n_t;
    }
};
*/

int main(int argc, char* argv[])
{
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

    vector<string> tmp_v;
    /// for option
    // {{{
    Opts opts;
    opts.get_opt_by_name<OptS>("version")->get_value_ref()=prog_version;
    OptS *opt_out_prefix = opts.create<OptS>("out_prefix,out,o", "STR", 
            "prefix of any output [out.]", "out.", 1, false);
    string & out_prefix = opt_out_prefix->get_value_ref();
    string & tmp_dir = opts.create<OptS>("tmp_dir,tmp,t", "STR", "template dir",
            "", 1, false)->get_value_ref();
    string & log_file = opts.create<OptS>("log_file,log,l", "FILE", "log file",
            "", 1, false)->get_value_ref();
    string & in_fasta_file = opts.create<OptS>("in_fasta_file,in_fasta,i", 
            "FILE", "input fasta file", "", 1, true)->get_value_ref();
    char & pam_pos = opts.create<OptC>("pam_pos", "CHAR",
            "Position of PAM, [5 (for 5')] or 3 (for 3')", "3", 1, false
            )->get_value_ref();
    int & target_length = opts.create<OptI>("target_length,target_len,tl", "N",
            "length of the target [23], include PAM", "23", 1, false
            )->get_value_ref();
    string & pam = opts.create<OptS>("pam", "STR", "PAM sequence, [NGG]", "NGG",
            1, false)->get_value_ref();
    string & head = opts.create<OptS>("head", "STR",
            "head sequence of the target", "", 1, false)->get_value_ref();
    string & tail = opts.create<OptS>("tail", "STR",
            "tail sequence of the target", "", 1, false)->get_value_ref();
    int & search_strand = opts.create<OptI>("search_strand,strand,ss", "N",
            "strand to search for target, [3:both], 1:plus, 2:minus",
            "3", 1, false)->get_value_ref();
    string & index_pre = opts.create<OptS>("index", "FILE",
            "index file prefix", "", 1, false)->get_value_ref();
    string & out_format_s = opts.create<OptS>("out_format,of,out_fmt", "STR",
            "output format: [tab] or html", "tab", 1, false)->get_value_ref();
    bool & is_clean = opts.create<OptB>("is_clean", "BOOL", 
            "whether to clean temporary files after finished, \
[T/True: clean], F/False: not clean.",
            "T", 1, false)->get_value_ref();

    opts.parse(argc, argv);
    if ( !opts.is_all_set() ) { return -1; }

    if (tmp_dir=="") {
        tmp_dir = out_prefix+"tmp";
    }
    {
        if (!File::exists_dir(tmp_dir)) {
            File::create_dir(tmp_dir);
            if (!File::exists_dir(tmp_dir)) {
                cerr << "Can't create TMP_DIR " << tmp_dir << endl;
                #ifdef PTHREAD
                pthread_exit(NULL);
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
        if (!File::exists_dir(out_dir)) {
            cerr << "Can't find OUT_DIR " << out_dir << ". You need to create it manually" << endl;
            exit(-1);
        }
    }
    /// }}}

    enum OutFormat { OF_TAB=0, OF_HTML=1 } out_format = OF_TAB;
    StringUtility::to_upper(out_format_s);
    if ( out_format_s=="HTML" ) {
        out_format = OF_HTML;
    }
    
    tmp_v.push_back(tmp_dir);

    Msg log(log_file);
    ss.str("");
    opts.print_config(ss);
    log << ss.str() << endl;

    ss.str("");
    ss<< "Begin at: " << CZL_CUR_TIME << "\n";
    cerr<< ss.str() << endl;
    log << ss.str() << endl;

    /*
    int32_t *fpos_v = NULL;
    int32_t *fcount_v = NULL;
    ifstream idx8_fin, idx24_fin;
    if ( !index_pre.empty() ) {
        fpos_v = new int32_t[0x10000];
        fcount_v = new int32_t[0x10000];
        if (fpos_v==NULL || fcount_v==NULL) {
            cerr << "E: file to allocate. " << CZL_DBG_INFO << endl;
            pthread_exit( (void*)-1);
        }
        in_file = index_pre + "idx8";
        idx8_fin.open(in_file.c_str());
        if (idx8_fin.fail()) {
            cerr << "E: Fail to open idx8 file " << in_file << CZL_DBG_INFO << endl;
            pthread_exit( (void*)-1);
        }
        in_file = index_pre + "idx24";
        idx24_fin.open(in_file.c_str());
        if (idx24_fin.fail()) {
            cerr << "E: Fail to open idx24 file " << in_file << CZL_DBG_INFO << endl;
            pthread_exit( (void*)-1);
        }

        for (i=0; i<0x10000; i++) {
            idx8_fin.read( (char*)(fpos_v+i), sizeof(int32_t) );
            idx8_fin.read( (char*)(fcount_v+i), sizeof(int32_t) );
        }
        idx8_fin.close();
    }
    */


    Node * root=NULL;
    try {
        if ( !index_pre.empty() ) {
            ss.str("");
            ss << "BEGIN load trie (" << CZL_CUR_TIME << ")";
            cerr<< ss.str() << endl;
            log << ss.str() << endl;

            in_file = index_pre + "mat.4x4";
            fin.open(in_file.c_str(), ios::binary);
            if (fin.fail()) {
                cerr << "E: Fail to open fasta file " << in_file << CZL_DBG_INFO << endl;
                pthread_exit( (void*)-1 );
            }
            fin.read( (char*)mis_v, 0x10000 );
            fin.close();
            // load trie index
            in_file = index_pre + "trie";
            fin.open(in_file.c_str(), ios::binary);
            if (fin.fail()) {
                cerr << "E: Fail to open fasta file " << in_file << CZL_DBG_INFO << endl;
                pthread_exit( (void*)-1 );
            }
            root = read_trie(fin);
            fin.close();

            ss.str("");
            ss << "END (" << CZL_CUR_TIME << ")";
            cerr<< ss.str() << endl;
            log << ss.str() << endl;
        }
    
        ss.str("");
        ss << "\nBEGIN look for targets (" << CZL_CUR_TIME << ")";
        cerr<< ss.str() << endl;
        log << ss.str() << endl;

        // {{{
        string name, seq;
        int pam_len = pam.size();
        transform(pam.begin(), pam.end(), pam.begin(), ::toupper);
        string cpam;
        complement_nt1_copy(pam, cpam);
        int Nn=0;
        /// change PAM sequence to BYTE format [ref to czl_bio_base.h]
        uint8_t byte_pam[pam_len];
        uint8_t byte_cpam[pam_len];
        for (i=0; i<pam_len; i++) byte_pam[i] = nt1_to_bit4(pam[i]);
        for (i=0; i<pam_len; i++) byte_cpam[i] = nt1_to_bit4(cpam[i]);
        ///
        out_file = out_prefix+"bed";
        fout.open(out_file.c_str());

        fin.open(in_fasta_file.c_str());
        if (fin.fail()) {
            cerr << "E: Fail to open fasta file " << in_fasta_file << CZL_DBG_INFO << endl;
            return -1;
        }

        while ( !fin.eof() ) {
            Fasta::get_a_seq(fin, name, seq);
            if ( seq.size() < target_length) {
                log.warn("Sequence shorter than target size.", __FILE__, __LINE__, 1);
                continue;
            }

            char strand;
            int Nn = 0;
            for (int i = 0; i<target_length; i++) {
                char c = toupper(seq[i]);
                if ( !(c=='A' || c=='C' || c=='G' || c=='T') ) {
                    Nn++;
                }
            }

            // check if both strand is target

            if (pam_pos=='5') {
                /// search plus strand
                for (int i = 0; i<seq.size()-target_length; i++) {
                    if ( Nn>0 ) {
                        char c = toupper(seq[i]);
                        if ( !(c=='A' || c=='C' || c=='G' || c=='T') ) {
                            Nn--;
                        }
                        c = toupper(seq[i+target_length]);
                        if ( !(c=='A' || c=='C' || c=='G' || c=='T') ) {
                            Nn++;
                        }
                        continue;
                    }
                //  short is_palin;
                    short tn = 0;
    //                if ( (search_strand&0x03) == 0x03 ) {
    //                    short is_palin=0;
    //                    if ( target_length%2==0 && search_strand&0x01 ) {
    //                        for (j=0; j<target_length/2; j++) {
    //                            char c0 = toupper(seq[i+j]);
    //                            char c1 = toupper(seq[i+target_length-1-j]);
    //                            if ( c0 != complement_nt1(c1) ) break;
    //                        }
    //                        if (j==target_length/2) { // palindrome
    //                            is_palin = 1;
    //                        }
    //                    }
    //                }
                    string seq1 = seq.substr(i, target_length);
                    StringUtility::to_upper(seq1);
                    if ( search_strand&0x01 ) {
                        strand = '+';
                        int Nn1 = Nn;
                        for (j=0; j<pam_len; j++) {
                            char c = toupper(seq[i+j]);
                            uint8_t b = nt1_to_byte( c );
                            if ( !(b & byte_pam[j]) || (b&byte_pam[j])!=b) break;
                            if ( c=='N' ) Nn1--;
                        }
                        if (j>=pam_len && Nn1<=1) {
                            tn++;
                            fout<< name << "\t" << i << "\t" << i+target_length
                                << "\t" << seq1 << "\t" 
                                << "1000" << "\t" << strand << "\t"
                                << i+pam_len << "\t" << i+target_length << "\n";
                        }
                    }
                    if ( search_strand&0x02 ) {
                        strand = '-';
                        int Nn1 = Nn;
                        for (j=0; j<pam_len; j++) {
                            char c = toupper( seq[i+target_length-1-j] );
                            uint8_t b = nt1_to_byte(c);
                            if ( !(b & byte_cpam[j]) || (b&byte_cpam[j])!=b) break;
                            if ( c=='N' ) Nn1--;
                        }
                        if (j>=pam_len && Nn1<=1) {
                            tn++;
                            string rev;
                            rev_complement_nt1_copy(seq1, rev);
                            if ( tn<2 || rev != seq1 ) {
                                fout<< name << "\t" << i
                                    << "\t" << i+target_length
                                    << "\t" << rev << "\t" << "1000"
                                    << "\t" << strand << "\t" << i
                                    << "\t" << i+target_length-pam_len << "\n";
                            }
                        }
                    }

                    char c = toupper(seq[i]);
                    if ( !(c=='A' || c=='C' || c=='G' || c=='T') ) {
                        Nn--;
                    }
                    c = toupper(seq[i+target_length]);
                    if ( !(c=='A' || c=='C' || c=='G' || c=='T') ) {
                        Nn++;
                    }
                }
            } else {
                // PAM is at 3'
                for (int i = 0; i<seq.size()-target_length; i++) {
                    if ( Nn>0 ) {
                        char c = toupper(seq[i]);
                        if ( !(c=='A' || c=='C' || c=='G' || c=='T') ) {
                            Nn--;
                        }
                        c = toupper(seq[i+target_length]);
                        if ( !(c=='A' || c=='C' || c=='G' || c=='T') ) {
                            Nn++;
                        }
                        continue;
                    }
    //                short is_palin;
    //                if ( (search_strand&0x03) == 0x03 ) {
    //                    short is_palin=0;
    //                    if ( target_length%2==0 && search_strand&0x01 ) {
    //                        for (j=0; j<target_length/2; j++) {
    //                            char c0 = toupper(seq[i+j]);
    //                            char c1 = toupper(seq[i+target_length-1-j]);
    //                            if ( c0 != complement_nt1(c1) ) break;
    //                        }
    //                        if (j==target_length/2) { // palindrome
    //                            is_palin = 1;
    //                        }
    //                    }
    //                }

                    short tn=0;
                    string seq1 = seq.substr(i, target_length);
                    StringUtility::to_upper(seq1);
                    if ( search_strand&0x01 ) {
                        strand = '+';
                        int Nn1 = Nn;
                        for (j=0; j<pam_len; j++) {
                            char c = toupper(seq[i+target_length-pam_len+j]);
                            uint8_t b = nt1_to_byte(c);
                            if ( !(b & byte_pam[j]) || !((b&byte_pam[j])==b) ) break;
                            if ( c=='N' ) Nn1--;
                        }
                        if (j>=pam_len) {
                            tn++;
                            fout<< name << "\t" << i << "\t" << i+target_length << "\t"
                                << seq.substr(i, target_length) << "\t" 
                                << "1000" << "\t" << strand << "\t"
                                << i << "\t" << i+target_length-pam_len << "\n";
                        }
                    }
                    if ( search_strand&0x02 && !tn ) {
                        strand = '-';
                        int Nn1 = Nn;
                        for (j=0; j<pam_len; j++) {
                            char c = toupper(seq[i+pam_len-1-j]);
                            uint8_t b = nt1_to_byte(c);
                            if ( !(b & byte_cpam[j]) || !((b&byte_cpam[j])==b) ) break;
                            if ( c=='N' ) Nn1--;
                        }
                        if (j>=pam_len) { // target
                            tn++;
                            string rev;
                            rev_complement_nt1_copy(seq1, rev);
                            if ( tn<2 || rev!=seq1 ) {
                                fout<< name << "\t" << i
                                    << "\t" << i+target_length << "\t" << rev
                                    << "\t" << "1000" << "\t" << strand
                                    << "\t" << i+pam_len
                                    << "\t" << i+target_length << "\n";
                            }
                        }
                    }
                    char c = toupper(seq[i]);
                    if ( !(c=='A' || c=='C' || c=='G' || c=='T') ) {
                        Nn--;
                    }
                    c = toupper(seq[i+target_length]);
                    if ( !(c=='A' || c=='C' || c=='G' || c=='T') ) {
                        Nn++;
                    }
                }
            }
        }
        fout.close();
        fin.close();
        // }}}

        ss.str("");
        ss << "END (" << CZL_CUR_TIME << ")";
        cerr<< ss.str() << endl;
        log << ss.str() << endl;

        ss.str("");
        ss << "\nBEGIN Fetch statistics (" << CZL_CUR_TIME << ")";
        cerr<< ss.str() << endl;
        log << ss.str() << endl;

        if ( out_format == OF_HTML ) {
//                cout<< "content-type: text/html\n";
//                cout<< "<html><head><title>Searching Result</title></head>";
//                cout<< "<body>";
            cout<< "<table style=\"border:1 solid black;\">";
            cout<< "<tr style=\"text-align:center;padding=2;\">"
                << "<th>Chromosom</th><th>Begin</th><th>End</th>"
                << "<th>Target</th><th>Strand</th><th>GC percentage</th>";
                
            if (root!=NULL) {
                cout<< "<th>0-mismatch Count</th>"
                    << "<th>1-mismatch count</th><th>2-mismatch count</th>"
                    << "<th>3-mismatch count</th>";
            }
            cout<< "</tr>\n";
        }
            
        in_file = out_prefix + "bed";
        fin.open(in_file.c_str());
        if (fin.fail()) {
            cerr << "E: Fail to open BED file " << in_file << CZL_DBG_INFO << endl;
            pthread_exit( (void*)-1 );
        }
        while (!fin.eof()) {
            string chr, seq;
            int begin, end;
            char strand;
            getline(fin, line);
            if (line.empty()) continue;
            i = 0;
            i = StringUtility::find(line, '\t', i, chr);
            i = StringUtility::find(line, '\t', i, str);
            begin = atoi(str.c_str());
            i = StringUtility::find(line, '\t', i, str);
            end = atoi(str.c_str());
            i = StringUtility::find(line, '\t', i, seq);
            i = StringUtility::find(line, '\t', i, str);
            i = StringUtility::find(line, '\t', i, str);
            strand = str[0];
            string seq1;
            if ( pam_pos=='3' ) {
                string rev_seq;
                rev_complement_nt1_copy(seq, rev_seq);
                seq1 = rev_seq.substr(pam_len);
            } else {
                seq1 = seq.substr(pam_len);
            }
            StringUtility::to_upper(seq1);

            // calculate GC content
            int gc=0;
            for (j=0; j<seq1.size(); j++) {
                if (seq1[j]=='C' || seq1[j]=='G') gc++;
            }
            short is_gc_good=1;
			if ( gc > 0.8*seq1.size() || gc < 0.2*seq1.size() ) {
				is_gc_good=0;
			}

            int w = 8, st = 4, wm=0;
			if ( is_gc_good ) {
				for (j=0; j<w; j++) {
					if (seq1[j]=='C' || seq1[j]=='G') {
						wm++;
					}
				}
				if (wm> 0.8*w || wm< 0.2*w) is_gc_good=0;
				int b=0, e=w;
				while (e < seq1.size()) {
					for (j=0; j<st; j++) {
						if (seq1[b]=='C' || seq1[b]=='G') wm--;
						b++;
						if (seq1[e]=='C' || seq1[e]=='G') wm++;
						e++;
					}
					if (wm> 0.8*w || wm< 0.2*w) {
						is_gc_good=0;
						break;
					}
				}
			}
            string color;
            //

            int score[5];
            for ( i=0; i<5; i++ ) score[i] = 0;
            if ( root!=NULL ) {
                search( root, seq1, mis_v, 3, score);
            }

            if ( is_gc_good==0 || score[0]>1 || score[1]+score[2]>10 ) {
                color="black";
            } else {
                color="red";
            }

            string s1, s2;
            if ( out_format==OF_HTML ) {
                cout << "<tr " << "style=\"color:" << color << ";\">";
                s1 = "<td>";
                s2 = "</td>";
                cout << s1;
            } else {
                s1 = "\t";
            }
            cout<< chr << s2;
            cout<< s1 << begin << s2;
            cout<< s1 << end   << s2;
            cout<< s1 << seq   << s2;
            cout<< s1 << strand<< s2;
            cout<< s1 << setprecision(3) << (float)gc*100/seq1.size() << s2;
            if ( root!=NULL ) {
                for ( i=0; i<4; i++ ) {
                    cout<< s1 << score[i] << s2 ;
                }
            }
            if ( out_format==OF_HTML ) {
                cout << "</tr>\n";
            } else {
                cout << "\n";
            }
        }
        fin.close();
        if (root!=NULL) {
            destroy_trie(root);
            root = NULL;
        }
        
        if ( out_format == OF_HTML ) {
            cout << "</table>";
        }

        ss.str("");
        ss << "END (" << CZL_CUR_TIME << ")\n";
        cerr<< ss.str() << endl;
        log << ss.str() << endl;

    } catch ( exception & e ) {
        cerr << "Exception: " << e.what() << " " << CZL_DBG_INFO << endl;
        pthread_exit( (void*)-1 );
    }
//    if ( out_format == OF_HTML ) {
//        cout << "</body></html>";
//    }

//  if (fpos_v!=NULL) delete []fpos_v;
//  if (fcount_v!=NULL) delete []fcount_v;
//  if (idx24_fin.is_open()) idx24_fin.close();

    if (is_clean) {
        for (i=0; i<tmp_v.size(); i++) {
            remove(tmp_v[i].c_str());
        }
    }
#ifdef PTHREAD
    pthread_exit(NULL);
#else
    return 0;
#endif
}

void print_usage()
{
    cout << "Usage:" << endl;
    cout << "  " << prog_name << "<Options>" << endl;
    cout << "Options (* for required):" << endl;
    cout << "  -c        :" << endl;
    cout << "  --config  :  configure file  [conf.txt] " << endl;
    cout << "  -h        :  print the usage" << endl;
    cout << "Configure file options:" << endl;
    cout << "  out_prefix = <STR> : output prefix    [out]" << endl;
    cout << "  tmp_dir    = <STR> : temporary directory [OUT_PREFIX.tmp]" << endl;
    cout << "  log_file   = <STR> : log_file [OUT_PREFIX.log]" << endl;
    cout << "  target_length = <INT> : target length including PAM;" << endl;
    cout << "  target_length = <INT> : target length including PAM;" << endl;
}

void print_head()
{
    cerr << "**************************************************" << endl;
    cerr << "* Program  : " << prog_name << endl;
    cerr << "* Version  : " << prog_version << endl;
    cerr << "* Author   : " << prog_author << endl;
    cerr << "* Function : " << prog_function << endl;
    cerr << "* Input    : " << prog_input << endl;
    cerr << "* Output   : " << prog_output << endl;
    cerr << "* Compiled Host: " << prog_compiled_host << endl;
    cerr << "* Compiled Time: " << __DATE__ << " : " __TIME__ << endl;
    cerr << "* Machine Bit: " << machine_bit() << endl;
    cerr << "* IS big endian: " << (is_be()?"yes":"no") << endl;
    cerr << "**************************************************" << endl;
}

int fetch_index(string const & seq, int32_t fpos_v[0x10000], int32_t fcount_v[0x10000], ifstream & idx24, int stat[4])
{
    int32_t n, i;
    for (i=0; i<4; i++) stat[i] = 0;
    int bl0 = seq.size();
    uint64_t mask=0;
    assert(bl0>=16);
    uint64_t rev_bseq0 = 0, bseq0 = 0;
    for (i=0; i<bl0; i++) {
        rev_bseq0 <<= 2;
        rev_bseq0 |= Bit2Seq::nt1_to_bit2(seq[i]);
    }
    bseq0 = rev64(bl0, rev_bseq0);
    uint16_t seed0 = bseq0 & 0xffff;
    uint16_t seed1 = (bseq0>>16) & 0xffff;
    bseq0 >>= 16;
    for (i=0; i<bl0; i++) {
        mask<<=1;
        mask |= 0x1;
    }
    bl0-=8;
    rev_bseq0 = rev64(bl0, bseq0);


    int K=3;

    if ( fpos_v[seed0]<0 ) stat[0] = 0;
    else {
        idx24.clear();
        idx24.seekg( fpos_v[seed0] );
        n = fcount_v[seed0];
        sub_fetch(bl0, bseq0, idx24, n, 3, stat);
    }

    int k0;
    for (k0=1; k0<=K; k0++) {
        // first 8-bp seed with mismatch k0
        uint8_t pos0[3];
        uint8_t n0 = 1;
        uint16_t s0 = seed0;
        pos0[0] = 0;
        while ( n0>0 ) {
            if ( pos0[n0-1] == 8 ) {
                n0--;
                pos0[n0-1]++;
            } else {
                if ( n0 < k0 ) {
                    pos0[n0] = pos0[n0-1]+1;
                    n0++;
                } else {
                    // generate all k0 mis-match seed 
                    // with mismatch position pos0
                    s0 = seed0;
                    uint8_t n01 = 1;
                    seti16_next(&s0, pos0[n01-1]);
                    while ( n01>0 ) {
                        if ( geti16(s0, pos0[n01-1]) == geti16(seed0, pos0[n01-1]) ) {
                            n01--;
                            seti16_next(&s0, pos0[n01-1]);
                        } else {
                            if (n01==k0) {
                                // search
                            //  print16(8, s0);
                                int stat1[K-k0+1];
                                idx24.seekg( fpos_v[s0] );
                                n = fcount_v[s0];
                                if (n>0) {
                                    sub_fetch(bl0, bseq0, idx24, n, K-k0, stat1);
                                    for (int j=0; j<=K-k0; j++) {
                                        stat[ k0+j ] += stat1[j];
                                    }
                                }

                                seti16_next(&s0, pos0[n01-1]);
                            } else {
                                uint8_t a = geti16(seed0, pos0[n01]);
                                a = (a+1) &0x3;
                                seti16( &s0, pos0[n01], a);
                                n01++;
                            }
                        }
                    }

                    pos0[n0-1]++;
                }
            }
        }
    //  cout << endl;
    //  cout << endl;
    }
    
    return 0;
}

int mis_match(int k, int l, uint64_t *a, uint64_t *b)
{
    uint64_t m = *a ^ *b;
    int n=0;
    while (m!=0 && l>0) {
        if ( (m&0x3)!=0 ) n++;
        if ( n>k ) break;
        m >>= 2;
        l--;
    }
    return n;
}

uint8_t geti16(uint16_t seed, uint8_t i)
{
    return seed>>(i<<1) &0x3;
}

void seti16(uint16_t *seed, uint8_t i, uint16_t a)
{
    (*seed) &= ~(0x3U<<(i<<1));
    (*seed) |= (a&0x3)<<(i<<1);
}

void seti16_next(uint16_t *seed, uint8_t i)
{
    uint8_t a=(*seed)>>(i<<1) &0x3;
    a = (a+1)&0x3;
    (*seed) &= ~(0x3U<<(i<<1));
    (*seed) |= (a&0x3)<<(i<<1);
}

void print16(uint8_t l, uint16_t seed)
{
    for (uint8_t i=0; i<l; i++) {
        cout << Bit2Seq::bit2_to_nt1( geti16(seed, i) );
    }
    cout << endl;
}

uint16_t rev16(uint8_t l, uint16_t seed)
{
    uint16_t r;
    for (uint8_t i=0; i<l; i++) {
        uint8_t a = geti16(seed, i);
        seti16(&r, l-1-i, (uint16_t)a);
    }
    return r;
}

void sub_fetch(uint8_t bl0, uint64_t bseq0, ifstream & idx24, int32_t n, int k, int *stat)
{
    uint64_t mask=0;
    int i;
    for (i=0; i<bl0; i++) {
        mask<<=2;
        mask |= 0x3;
    }
    for (i=0; i<=k; i++) stat[i]=0;
    vector<uint64_t> bseq_v;
    vector<uint64_t> rev_bseq_v;
    vector<uint8_t> n_v;
    for (i=0; i<n; i++) {
        uint64_t b=0;
        uint8_t a;
        for (int j=0; j<6; j++) {
            a = 0;
            idx24.read( (char*)&a, 1 );
            uint64_t aa = a;
            b |= aa<< (j<<3);
        }
        idx24.read( (char*)&a, 1 ); // read the occured number
        n_v.push_back(a);
        b = b & mask;
        bseq_v.push_back(b);
        b = rev64(bl0, b);
    //  uint64_t rb = rev64(24,b);
        rev_bseq_v.push_back(b);
    }
    for ( i=0; i<n; i++) {
    //  print64(bl0, bseq0);
    //  print64(bl0, bseq_v[i]);
    //  cout << endl;
        int mis = mis_match( k, bl0, &bseq0, &bseq_v[i] );
        if (mis>=k) continue;
        stat[mis] += n_v[i];
    }
}

