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
#include "../czl_bio/czl_common.h"
#include "../czl_bio/czl_io.h"
#include <boost/filesystem.hpp>

using namespace std;
using namespace boost;
using namespace czl_bio;

void print_usage();
void print_head();

string prog_name="main";
string prog_version="1.0";
string prog_author="Zelin Chen";
string prog_compiled_host;
string prog_function="";
string prog_input="";
string prog_output="";

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
    bool & is_clean = opts.create<OptB>("is_clean", "BOOL", 
            "whether to clean temporary files after finished, \
[T/True: clean], F/False: not clean.", "T", 1, false)->get_value_ref();

    opts.parse(argc, argv);
    if ( !opts.is_all_set() ) { return -1; }

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
        if (!boost::filesystem::exists(out_dir)) {
            cerr << "Can't find OUT_DIR " << out_dir << ". You need to create it manually" << endl;
            exit(0);
        }
    }
    /// }}}

    Msg log(log_file);
    ss.str("");
    opts.print_config(ss);
    log << ss.str() << endl;

    fin.open(in_fasta_file.c_str());
    string name, seq;
    int pam_len = pam.size();
    int Nn=0;
    /// change PAM sequence to BYTE format [ref to czl_bio_base.h]
    uint8_t byte_pam[pam_len];
    for (int i=0; i<pam_len; i++) byte_pam[i] = nt1_to_byte(pam[i]);
    ///
    uint8_t byte_seq1[pam_len];
    out_file = out_prefix+"bed";
    fout.open(out_file.c_str());
    while ( !fin.eof() ) {
        Fasta::get_a_seq(fin, name, seq);
        if ( seq.size() < target_length) {
            log.warn("Sequence shorter than target size.", __FILE__, __LINE__, 1);
            continue;
        }

        char strand;
        if (pam_pos=='5') {
            /// search plus strand
            for (int i = 0; i<seq.size()-target_length; i++) {
                int j;
                if ( search_strand&0x01 ) {
                    strand = '+';
                    for (j=0; j<pam_len; j++) {
                        uint8_t b = nt1_to_byte(seq[i+j]);
                        if ( !(b & byte_pam[j]) ) break;
                    }
                    if (j<pam_len) continue; // not target

                    fout<< name << "\t" << i << "\t" << i+target_length
                        << "\t" << seq.substr(i, target_length) << "\t" 
                        << "1000" << "\t" << strand << "\t"
                        << i+pam_len << "\t" << i+target_length << "\n";
                }
                if ( search_strand&0x02 ) {
                    strand = '-';
                    for (j=0; j<pam_len; j++) {
                        uint8_t b = nt1_to_byte(complement_nt1(seq[i+target_length-1-j]));
                        if ( !(b & byte_pam[j]) ) break;
                    }
                    if (j<pam_len) continue; // not target

                    string rev;
                    rev_complement_copy(seq.substr(i, target_length), rev);
                    fout<< name << "\t" << i << "\t" << i+target_length << "\t"
                        << rev << "\t" << "1000" << "\t" << strand << "\t"
                        << i << "\t" << i+target_length-pam_len << "\n";
                }
            }
        } else {
            // PAM is at 3'
            for (int i = 0; i<seq.size()-target_length; i++) {
                if ( search_strand&0x01 ) {
                    strand = '+';
                    int j;
                    for (j=0; j<pam_len; j++) {
                        uint8_t b = nt1_to_byte(seq[i+target_length-pam_len+j]);
                        if ( !(b & byte_pam[j]) ) break;
                    }
                    if (j<pam_len) continue; // not target

                    fout<< name << "\t" << i << "\t" << i+target_length << "\t"
                        << seq.substr(i, target_length) << "\t" 
                        << "1000" << "\t" << strand << "\t"
                        << i << "\t" << i+target_length-pam_len << "\n";
                }
                if ( search_strand&0x02 ) {
                    strand = '-';
                    for (j=0; j<pam_len; j++) {
                        uint8_t b = nt1_to_byte(complement_nt1(seq[i+pam_len-1-j]));
                        if ( !(b & byte_pam[j]) ) break;
                    }
                    if (j<pam_len) continue; // not target

                    string rev;
                    rev_complement_copy(seq.substr(i, target_length), rev);
                    fout << name << "\t" << i << "\t" 
                         << i+target_length << "\t" << rev << "\t" << "1000"
                         << "\t" << strand << "\t" << i+pam_len << "\t"
                         << i+target_length << "\n";
                }
            }
        }
    }
    fout.close();
    fin.close();

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
    cout << "**************************************************" << endl;
    cout << "* Program  : " << prog_name << endl;
    cout << "* Version  : " << prog_version << endl;
    cout << "* Author   : " << prog_author << endl;
    cout << "* Function : " << prog_function << endl;
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
