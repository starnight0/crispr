/**
 * Copyright 2015, Zelin Chen <chzelin@gmail.com>
 *
 * @file filt_pam.cpp
 * @brief filter PAM-not-exact-match alignment
 * 
 * seed0 ( first 8-bp ) map + other (last 24 bp) trie 
 *
 * @author Zelin Chen
 * 
 */
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include "czl_common.hpp"
#include "czl_io.hpp"
#include "czl_bio_seq.hpp"
#include "bit2.hpp"
//#include <boost/filesystem.hpp>

using namespace std;
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

class Mer {
public:
    Mer(): n(0), flag(0) {}
    Mer(int n_, uint8_t flag_=0): n(n_), flag(flag_) {}
    Mer & operator=(Mer & src) {
        n    = src.n;
        flag = src.flag;
        return (*this);
    }
    ~Mer() {}
    int n; // count
    uint8_t flag;
};

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
//    cout << "REMARK:\n";
//    cout << "  bed file must be sorted by name" << endl;

    /// for option
    // {{{
    Opts opts;
    opts.get_opt_by_name<OptS>("version")->get_value_ref()=prog_version;
    OptS *opt_out_prefix = opts.create<OptS>("out_prefix,out,o", "STR", "prefix of any output", "out.", 1, true);
    string & out_prefix = opt_out_prefix->get_value_ref();
    string & tmp_dir = opts.create<OptS>("tmp_dir,tmp,t", "STR", "template dir", "", 1, false)->get_value_ref();
    string & log_file = opts.create<OptS>("log_file,log,l", "STR", "log file", "", 1, false)->get_value_ref();
    bool & is_clean = opts.create<OptB>("is_clean", "BOOL", "whether to clean temporary files after finished, T/True: clean, F/False: not clean.", "T", 1, false)->get_value_ref();

    string & in_bed_file = opts.create<OptS>("in_bed,i", "FILE",
            "input bed file, 'name' column is 32-mer+PAM seq", 
            "", 1, false)->get_value_ref();

    string & in_fasta_file = opts.create<OptS>("in_fasta,in_fa", "FILE",
            "input fasta sequence file", "", 1, false)->get_value_ref();

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
        if (!File::exists_dir(tmp_dir)) {
            File::create_dir(tmp_dir);
            if (!File::exists_dir(tmp_dir)) {
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
        if (!File::exists_dir(out_dir)) {
            cerr << "Can't find OUT_DIR " << out_dir << ". You need to create it manually" << endl;
            exit(0);
        }
    }
    /// }}}

    if ( in_bed_file.empty() && in_fasta_file.empty() ) {
        cerr<< "E: Either 'in_bed' or 'in_fasta' should be provided. (" 
            << CZL_DBG_INFO << ")" << endl; 
        pthread_exit( (void*)-1 );
    }

    vector<string> tmp_v(1, tmp_dir);

    int K = 32;

    Msg log(log_file);
    rawtime = time(NULL);
    cerr<< "Begin at: " << str_time(rawtime) << ", " 
        << clock()/CLOCKS_PER_SEC << endl;
    log << "Begin at: " << str_time(rawtime) << ", " 
        << clock()/CLOCKS_PER_SEC << endl;

    int pam_len = pam.size();
    Bit2Seq bpam(pam);

    int seq_len;
    uint64_t bseq = 0;
    uint8_t byte=0;
    try {
        rawtime = time(NULL);
		cout<< "\nBEGIN\nBuild Trie ("
            << str_time(rawtime) << ", " << (float)clock()/CLOCKS_PER_SEC << ")" << endl;
		log << "\nBEGIN\nBuild Trie at: "
            << str_time(rawtime) << ", " << (float)clock()/CLOCKS_PER_SEC << ")" << endl;
    //  vector<uint64_t> bseq_v;
    //  vector<uint8_t> bc_v; // bseq count
        seq_len = 32;
        Node * root=NULL;
        int count=0;
        if ( !in_bed_file.empty() ) {
            string in_file = in_bed_file;
            fin.open(in_file.c_str());
            if (fin.fail()) {
                cerr << "E: Fail to open file " << in_file << CZL_DBG_INFO << endl;
                pthread_exit( (void*)(-1) );
            }
            while (!fin.eof()) {
                getline(fin, line);
                StringUtility::trim(line, " \t");
                if ( line.empty() ) continue;

                i=0;
                i = line.find('\t', i);
                i++;
                i = line.find('\t', i);
                i++;
                i = line.find('\t', i);
                i++;
                string seq;
                i = StringUtility::find(line, '\t', i, seq);
                if (seq.empty()) continue;

                if ( pam_pos=='3' ) {
                    rev_complement_nt1(seq);
                }
                seq = seq.substr(pam_len);
                int len = 0;
                // change to bit2 format
                bseq = 0;
                str_to_bit2_64(seq, &len, &bseq);
                root = build_trie_insert( len, bseq, root );
                count++;
            }
            fin.close();
        } else if ( !in_fasta_file.empty() ) {
            string in_file = in_fasta_file;
            fin.open(in_file.c_str());
            if (fin.fail()) {
                cerr << "E: Fail to open file " << in_file << CZL_DBG_INFO << endl;
                pthread_exit( (void*)(-1) );
            }
            int Nn=0;
            string name0, seq0;
            while ( !fin.eof() && !Fasta::get_a_seq(fin, name0, seq0) ) {
                if ( seq0.empty() || seq0.size()<seq_len ) continue;
                Nn = 0;
                for (i=0; i<seq_len-1 && i<seq0.size(); i++) {
                    if ( (seq0[i] & ~0x20) == 'N' ) {
                        Nn++;
                    }
                }
                for ( i=0; i<=seq0.size()-seq_len; i++ ) {
                    if ( (seq0[i+seq_len-1] & ~0x20) =='N' ) Nn++;
                    if ( Nn==0 ) {
                        int len = 0;
                        string seq = seq0.substr(i, seq_len);

                        // change to bit2 format
                        len = 0;
                        bseq = 0;
                        str_to_bit2_64(seq, &len, &bseq);
                        root = build_trie_insert( len, bseq, root );
                        count++;

                        rev_complement_nt1(seq);

                        len = 0;
                        uint64_t rc_bseq = 0;
                        str_to_bit2_64(seq, &len, &rc_bseq);
                        if ( rc_bseq != bseq ) {
                            root = build_trie_insert( len, rc_bseq, root );
                            count++;
                        }
                    }
                    if ( (seq0[i] & ~0x20) =='N' ) Nn--;
                }
            }
            fin.close();
        }

        if ( root==NULL ) {
            cerr<< "E: Fail to build index. " << CZL_DBG_INFO << endl;
            pthread_exit( (void*)(-1) );
        } else {
            container_sort(root);
        }
        out_file = out_prefix + "trie";
        fout.open(out_file.c_str(), ios::binary);
        if ( fout.fail() ) {
            cerr<< "Fail to open file " << out_file << " "
                << CZL_DBG_INFO << endl;
            log << "Fail to open file " << out_file << " "
                << CZL_DBG_INFO << endl;
            pthread_exit( (void*)(-1) );
        }
        write_trie(fout, root, 0);
        fout.close();
        destroy_trie(root);

//        in_file = out_file;
//        fin.open(in_file.c_str(), ios::binary);
//        if ( fin.fail() ) {
//            cerr<< "Fail to open file " << in_file << " "
//                << CZL_DBG_INFO << endl;
//            log << "Fail to open file " << in_file << " "
//                << CZL_DBG_INFO << endl;
//            pthread_exit( (void*)(-1) );
//        }
//        root = read_trie(fin);
//        fin.close();
//
//        out_file = out_prefix + "trie1_copy";
//        fout.open(out_file.c_str(), ios::binary);
//        write_trie(fout, root, 0);
//        fout.close();
//        destroy_trie(root);
        rawtime = time(NULL);
		cout<< "END (" << str_time(rawtime) << ", "
            << (float)clock()/CLOCKS_PER_SEC << ")" << endl;
		log << "END (" << str_time(rawtime) << ", "
            << (float)clock()/CLOCKS_PER_SEC << ")" << endl;
	} catch (exception & e) {
		cerr<< "E: " << e.what() << CZL_DBG_INFO << endl;
		pthread_exit( (void*)-1);
	}

    // calculater matrix with mismatch number between any 4-bp by 4-bp pair
    uint8_t mis_v[0x10000];
    i=0;
    for ( int b1=0; b1<0x100; b1++ ) {
        for ( int b2=0; b2<0x100; b2++ ) {
            uint8_t mis=0;
            uint8_t b = (b1^b2) & 0xff;
            for ( k=0; k<4; k++ ) {
                if ( b&0x3 ) mis++;
                b >>= 2;
            }
            mis_v[i++] = mis;
        }
    }
    out_file = out_prefix + "mat.4x4";
    fout.open(out_file.c_str(), ios::binary);
    fout.write( (char*)mis_v, 0x10000 );
    fout.close();
    //
    /*
    map<uint64_t, Mer> mer;
    try {
    //  vector<uint64_t> bseq_v;
    //  vector<uint8_t> bc_v; // bseq count
        string in_file = in_bed_file;
        fin.open(in_file.c_str());
        if (fin.fail()) {
            cerr << "E: Fail to open file " << in_file << CZL_DBG_INFO << endl;
            pthread_exit( (void*)(-1) );
        }
        seq_len = 32;
        while (!fin.eof()) {
            getline(fin, line);
            StringUtility::trim(line, " \t");
            if ( line.empty() ) continue;

            i=0;
            i = line.find('\t', i);
            i++;
            i = line.find('\t', i);
            i++;
            i = line.find('\t', i);
            i++;
            string seq;
            i = StringUtility::find(line, '\t', i, seq);
            if (seq.empty()) continue;

            if ( pam_pos=='5' ) {
                seq = seq.substr(pam_len);
            } else {
                seq = seq.substr(0, seq.size()-pam_len);
            }
            int len = 0;
            // change to bit2 format
            bseq = 0;
            str_to_bit2_64(seq, &len, &bseq);
            if ( mer.find(bseq)==mer.end() ) {
                Mer d(1);
                mer[bseq] = d;
            } else {
                mer[bseq].n ++;
            }
        //  bseq_v.push_back(bseq);
        }
    //  mk_sort( seq_len, bseq_v ) ;
    //  j=0;
    //  bc_v.push_back(1);
    //  for (i=1; i<bseq_v.size(); i++) {
    //      if ( bseq_v[i] == bseq_v[j] ) {
    //          if (bc_v[j]<255) bc_v[j]++;
    //      } else {
    //          j++;
    //          if ( i!=j ) {
    //              bseq_v[j] = bseq_v[i];
    //          }
    //          bc_v.push_back(1);
    //      }
    //  }
    //  j++;
    //  bseq_v.resize(j);

        //
    //  vector<int> order;
    //  mk_sort(seq_len, bseq_v, 1, order);
    //  vector<int> inv_order( order.size() );
    //  for (i=0; i<order.size(); i++) {
    //      rev_order[ order[i] ] = i;
    //  }
    //  vector<bool> pass0(bseq_v.size(), false);
    //  vector<bool> pass1(bseq_v.size(), false);
        int i=0, j=0;
    //  Less64 less64(0, seq_len-1);
    //  uint64_t mask_lo = 0x3fffffffffffffff;
        uint64_t mask_hi2 = 0x3<<((K-1)<<1);
        uint64_t mask = 0;
        for (i=0; i<K; i--) {
            mask = mask<<2 | 0x3;
        }
        map<uint64_t, Mer>::iterator mer_it=mer.begin();
        int m=0;
        while (mer_it!=mer.end()) {
            if ( mer_it->second.flag&0x01 ) {
                mer_it++;
                continue;
            }

            bseq = mer_it->first<<2 & mask;
            uint64_t a = 0;
            for ( a=0; a<4; a++ ) {
                bseq |= a;
                if ( mer.find(bseq) == mer.end() ) continue;
                if ( mer[bseq].flag & 0x01 ) continue;
                break;
            }
            if ( a<4 ) { // can exterd backward, skip
                mer_it++;
                continue;
            }

            bseq = mer_it->first;
            string seq;
            bit2_to_str_64(K, bseq, seq);
            int l = 0;
            while (1) {
                mer[bseq].flag |= 0x1;
                uint64_t a = 0;
                for ( a=0; a<4; a++ ) {
                    bseq >>= 2;
                    bseq |= (a<<((K-1)<<1));
                    if ( mer.find(bseq) == mer.end() ) continue;
                    if ( mer[bseq].flag & 0x01 ) continue;
                    break;
                }
                if ( a==4 ) break;
                seq.push_back( Bit2Seq::bit2_to_nt1(a) );
                l++;
            }
            if (l>2) {
                cout << seq << "\n" << l << endl;
            }
            m++;
            mer_it++;
        }
    } catch (exception & e) {
        cerr<< "E: " << e.what() << CZL_DBG_INFO << endl;
        pthread_exit( (void*)-1);
    }
    */

    


    /*
    // build index
    try {
        const short BL = 16;
        int lbuck = 0x1<<BL;
        uint64_t *buck = new uint64_t[0x100U<<BL];
        int *buck_i = new int[0x100];
        int *total0 = new int[0x100];
        int i, j;
        for (i=0; i<0x100; i++) {
            buck_i[i] = i<<BL;
        }
        int *total16 = new int[0x10000];
        uint64_t bseq=0;

        // sort by first byte ( first 4-bp )
        // {{{
        string in_file = in_bed_file;
        fin.open(in_file.c_str());
        if (fin.fail()) {
            cerr << "E: Fail to open file " << in_file << CZL_DBG_INFO << endl;
            pthread_exit( (void*)(-1) );
        }
        while (!fin.eof()) {
            getline(fin, line);
            StringUtility::trim(line, " \t");
            if ( line.empty() ) continue;

            i=0;
            i = line.find('\t', i);
            i++;
            i = line.find('\t', i);
            i++;
            i = line.find('\t', i);
            i++;
            string seq;
            i = StringUtility::find(line, '\t', i, seq);
            if (seq.empty()) continue;

            if ( pam_pos=='5' ) {
                seq = seq.substr(pam_len);
            } else {
                seq = seq.substr(0, seq.size()-pam_len);
            }
            int len = 0;
            // change to bit2 format
            bseq = 0;
            str_to_bit2_64(seq, &len, &bseq);
            // get seed0
            uint8_t seed0 = bseq & 0xff;
            // add bseq to bucket seed0
            buck[ buck_i[seed0]++ ] = bseq;
            // output if bucket is full
            if ( buck_i[seed0]==lbuck) {
                string file = tmp_dir + "/byte0_" + itos(seed0); 
                ofstream fout(file.c_str(), ios::binary|ios::app);
                if ( fout.fail() ) {
                    cerr<< "Fail to open file " << file << " "
                        << CZL_DBG_INFO << endl;
                    log << "Fail to open file " << file << " "
                        << CZL_DBG_INFO << endl;
                    pthread_exit( (void*)(-1) );
                }
                fout.write( (char*)(buck+(seed0<<BL)), lbuck*sizeof(uint64_t) );
                fout.close();
                buck_i[seed0] = seed0<<BL;
            }
            total0[seed0]++;
        }
        fin.close();
        // }}}

        // sort by second byte ( second 4-by )
        for (uint32_t b0=0; b0<0x100; b0++) {
            string file = tmp_dir + "/byte0_" + itos(b0); 
            if ( !File::exists_file(file) && (buck_i[b0]-b0<<BL)==0) continue;

            uint64_t *buck1 = new uint64_t[0x100U<<BL];
            int *buck1_i = new int[0x100];
            int *total1 = new int[0x100];
            for (i=0; i<0x100; i++) {
                buck1_i[i] = i<<BL;
                total1[i] = 0;
            }
            // if buffer file exists, load it
            // {{{
            if ( File::exists_file(file) ) {
                ifstream fin1(file.c_str());
                while (fin1.eof()) {
                    bseq=0;
                    fin1.read( (char*)&bseq, sizeof(uint64_t) );
                    uint16_t byte1 = (bseq>>8)&0xff;
                    bseq >>= 16;
                    buck1[ buck1_i[byte1]++ ] = bseq;
                    if ( buck1_i[byte1]==lbuck) {
                        string file = tmp_dir + "/byte1_" + itos(byte1); 
                        ofstream fout(file.c_str(), ios::binary|ios::app);
                        if ( fout.fail() ) {
                            cerr<< "Fail to open file " << file << " "
                                << CZL_DBG_INFO << endl;
                            log << "Fail to open file " << file << " "
                                << CZL_DBG_INFO << endl;
                            pthread_exit( (void*)(-1) );
                        }
                        fout.write( (char*)(buck1+(byte1<<BL)), lbuck*sizeof(uint64_t) );
                        fout.close();
                        buck1_i[byte1] = byte1<<BL;
                    }
                    total1[byte1]++;
                }
                fin1.close();
            }
            // }}}
            // read from buffer 'buck'
            for (int k=b0<<BL; k<buck_i[b0]; k++) {
                bseq = buck[k];
                uint16_t byte1 = (bseq>>8)&0xff;
                bseq >>= 16;
                buck1[ buck1_i[byte1]++ ] = bseq;
                if ( buck1_i[byte1]==lbuck) {
                    string file = tmp_dir + "/byte1_" + itos(byte1); 
                    ofstream fout(file.c_str(), ios::binary|ios::app);
                    if ( fout.fail() ) {
                        cerr<< "Fail to open file " << file << " "
                            << CZL_DBG_INFO << endl;
                        log << "Fail to open file " << file << " "
                            << CZL_DBG_INFO << endl;
                        pthread_exit( (void*)(-1) );
                    }
                    fout.write( (char*)(buck1+(byte1<<BL)), lbuck*sizeof(uint64_t) );
                    fout.close();
                    buck1_i[byte1] = byte1<<BL;
                }
                total1[byte1]++;
            }
            //
            
            for ( uint32_t b1=0; b1<0x100; b1++) {
                uint16_t seed0 = b0 | b1<<8;
                int n = total1[b1];
                total16[ seed0 ] = total1[b1];

                if ( n<=0 ) continue;

                uint64_t *bseq_v = new uint64_t[n];
                string file = tmp_dir + "/byte1_" + itos(b1); 
                if ( File::exists_file(file) ) {
                    ifstream fin1(file.c_str());
                    fin1.read( (char*)bseq_v, (n<<BL>>BL)*sizeof(uint64_t) );
                    fin1.close();
                }
                i = n>>BL<<BL;
                memcpy( bseq_v+i, buck1+(b1<<BL), (n-i)*sizeof(uint64_t));

                vector<Node*> node_v;
                Node * root = build_trie(n, bseq_v, node_v);
                for (int j=0; j<node_v.size(); j++ ) {
                    delete node_v[j];
                }

                delete []bseq_v;
            }
            delete []total1;
            delete []buck1;
            delete []buck1_i;
        }

        delete []total0;
        delete []total16;
        delete []buck;
        delete []buck_i;
    } catch ( exception & e ) {
        cerr<< "E: " << e.what() << CZL_DBG_INFO << endl;
        pthread_exit( (void*)-1);
    }

    out_file = tmp_dir + "/count";
    fout.open(out_file.c_str());
    if ( fout.fail() ) {
        cerr<< "E: Fail to open file " << out_file << CZL_DBG_INFO << endl;
        log << "E: Fail to open file " << out_file << CZL_DBG_INFO << endl;
        pthread_exit( (void*)(-1) );
    } else {
        tmp_v.push_back(out_file);
    }
    
    string name0;
    n=0;
    while (!fin.eof()) {
        getline(fin, line);
        StringUtility::trim(line, " \t");
        if ( line.empty() ) continue;

        i=0;
        i = line.find('\t', i);
        i++;
        i = line.find('\t', i);
        i++;
        i = line.find('\t', i);
        i++;
        string name;
        i = StringUtility::find(line, '\t', i, name);
        if (name.empty()) continue;

        if ( pam_pos=='5' ) {
            name = name.substr(pam_len);
        } else {
            name = name.substr(0, name.size()-pam_len);
        }
        StringUtility::to_upper(name);
        if ( name != name0 ) {
            if ( !name0.empty() ) {
                C c;
                for (i=0; i<name0.size(); i++) {
                    c.seq<<=2;
                    c.seq |= Bit2Seq::nt1_to_bit2(name0[i]);
                }
                c.n = n;
                C_write(fout, c);
            }
            n = 0;
            name0 = name;
        }
        n++;
    }
    if ( !name0.empty() ) {
        fout << name0 << "\t" << n << "\n";
    }
    fin.close();
    fout.close();

    in_file = tmp_dir + "/count";
    out_file = tmp_dir + "/count.sort";
    IOMergeSort<C> ioms(C_read, C_write, C_less);
    ioms.run(in_file, out_file, tmp_dir, 100000, 100);
    rename(out_file.c_str(), in_file.c_str());
    
    cout<< "END\n" << endl;
    log << "END\n" << endl;
    */

	/*
    ss.str("");
    ss<< "Time: " << ctime(&rawtime) << "\t" << clock()/CLOCKS_PER_SEC << "\n";
    cout<< ss.str() << endl;
    log << ss.str() << endl;

    cout<< "BEGIN\nIndex target (not include PAM)" << endl;
    log << "BEGIN\nIndex target (not include PAM)" << endl;

    in_file = tmp_dir + "/count";
    fin.open(in_file.c_str());
    if (fin.fail()) {
        cerr << "E: Fail to open file " << in_file << CZL_DBG_INFO << endl;
        pthread_exit( (void*)(-1) );
    }

    out_file = out_prefix + "idx8";
    fout.open(out_file.c_str());
    if ( fout.fail() ) {
        cerr<< "E: Fail to open file " << out_file << CZL_DBG_INFO << endl;
        log << "E: Fail to open file " << out_file << CZL_DBG_INFO << endl;
        pthread_exit( (void*)(-1) );
    }

    ofstream fout1;
    out_file = out_prefix + "idx24";
    fout1.open(out_file.c_str());
    if ( fout1.fail() ) {
        cerr<< "E: Fail to open file " << out_file << CZL_DBG_INFO << endl;
        log << "E: Fail to open file " << out_file << CZL_DBG_INFO << endl;
        pthread_exit( (void*)(-1) );
    }

    int32_t *pos_v = new int32_t[ 0x10000 ];
    int32_t *count_v = new int32_t[ 0x10000 ];
    for (i=0; i<0x10000; i++) {
        pos_v[i] = -1;
        count_v[i] = 0;
    }
    Bit2Seq bseq0;
    uint16_t seed0;
    int64_t count = 0;
    int32_t fpos=0;
    while (!fin.eof()) {
        C c;
        C_read(fin, c);
        if ( fin.fail() ) continue;
        uint64_t b = rev64(32, c.seq);

        i=0;
        uint8_t m;
        if (c.n>=255) m = 255;
        else m = c.n;

        uint16_t seed = b&0xffff;
        b >>=16;
        if ( count==0 ) {
            seed0 = seed;
            pos_v[seed0] = fpos;
        } else {
            if ( seed!=seed0) {
                seed0 = seed;
                pos_v[seed0] = fpos;
            }
        }
        uint8_t a;
        for (i=0; i<6; i++) {
            a = b&0xff;
            fout1.write((char*)&a, sizeof(a));
            b >>= 8;
        }
        fout1.write((char*)&m, sizeof(m));
        fpos += 7;
        count++;
        count_v[seed]++;
    }
    fout1.close();
    for (i=0; i<0x10000; i++) {
        fout.write( (char*)&pos_v[i], sizeof(pos_v[i]));
        fout.write( (char*)&count_v[i], sizeof(count_v[i]));
    }
    fout.close();
    delete []pos_v;
    ss.str("");
    ss<< "End at: " << ctime(&rawtime) << "\t" << clock()/CLOCKS_PER_SEC;
    cout<< ss.str() << endl;
    log << ss.str() << endl;
	*/

    if ( is_clean ) {
        for (i=tmp_v.size()-1; i>=0; i--) {
            remove(tmp_v[i].c_str());
        }
    }
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

