/**
 * Copyright 2015, Zelin Chen <chzelin@gmail.com>
 *
 * @file gen_u_target.cc
 * @brief 
 *
 * @author Zelin Chen
 * 
 */
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include "czl_bio_v2/czl_common.h"
#include "czl_bio_v2/czl_io.h"
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace boost;
using namespace czl_bio;

void print_usage();
void print_head();

string prog_name="gen_u_target";
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

int gen_u_target(int argc, char* argv[])
{
    print_head();

	/// for option
	// {{{
	map<string, string> opt;
	map<string, int> opt_n;
	map<string, string> opt_sl; // mapping short option to long option

	opt_n["help"]       = 0;
	opt_sl["h"] = "help";

	opt_n["config"]     = 1;
	opt_sl["c"] = "config";
    string conf_file;

	opt_n["out_prefix"] = 1;
	opt_sl["o"] = "out_prefix";
    string out_prefix = prog_name+".out";

	opt_n["tmp_dir"]    = 1;
	opt_sl["t"] = "tmp_dir";
    string tmp_dir;

	opt_n["log_file"]   = 1;
	opt_sl["l"] = "log_file";
    string log_file;
	/// }}}

	/// gen_u_target global;
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

	/// get option
	// {{{
    for (k=1; k<argc; k++) {
		string name;
		string value;
		if ( argv[k][0]=='-' ) {
			if ( argv[k][1]=='-' ) {
				name = argv[k]+2;
			} else {
				name = argv[k]+1;
			}
			if ( name=="h" || name == "help") {
				print_usage();
				#ifdef PTHREAD
				pthread_exit(NULL);
				#else
				return 0;
				#endif
			} else if ( name=="v" || name == "version") {
				cout << prog_version << endl;
				#ifdef PTHREAD
				pthread_exit(NULL);
				#else
				return 0;
				#endif
			}

			for (i=0; i<name.size(); i++) {
				if (name[i]=='-') name[i]='_';
			}
			if ( opt_sl.find(name)!=opt_sl.end() ) name = opt_sl[name];
			if ( opt_n.find(name)==opt_n.end() ) {
				cerr << "Error: Unavailable option '" << name << "'" << endl;
				#ifdef PTHREAD
				pthread_exit(NULL);
				#else
				return 0;
				#endif
			} else {
				if (opt_n[name]==1) {
					k++;
					if (k>=argc) {
						cerr << "Error: option " << name << "need value." << endl;
						#ifdef PTHREAD
						pthread_exit(NULL);
						#else
						return 0;
						#endif
					}
					value.assign(argv[k]);
				}
			}
			opt[name] = value;
		} else {
            cerr << "Error: option not start with '-': " << argv[k] << endl;
			#ifdef PTHREAD
			pthread_exit(NULL);
			#else
			return 0;
			#endif
		}
	}
	/// }}}

	if ( opt.find("config")!=opt.end() ) { // -c OR --config
		/// read configure file
		/// {{{
		conf_file = opt["config"];
		opt.erase("config");
		cerr << "Read from configure file " << conf_file << endl;
		if (!conf_file.empty()) {
			fs.open(conf_file.c_str(), fstream::in);
			if (fs.fail()) {
				cerr << "Can't read configure file " << conf_file << endl;
				return 1;
			}
			while (!fs.eof()) {
				getline(fs, str);
				boost::trim(str);

				if (str.empty()) { continue; }

				if (str[0]=='[') { continue; }
				if (str[0]=='#') { continue; }

				i=str.find_first_of('=');
				if (i == std::string::npos) continue;
				
				string name, value;
				name=str.substr(0, i);
				value=str.substr(i+1);

				if (name.empty()) continue;
				boost::trim(name);
				if (value.empty()) continue;
				boost::trim(value);
				if ( opt_sl.find(name)!=opt_sl.end() ) name = opt_sl[name];
                if ( opt_n.find(name)!=opt_n.end() ) {
					opt[name] = value;
				} else {
					cerr << "Can't recognize option " << name << " = " << value << " Skip." << endl;
				}
			}
			fs.close();
		}
		/// }}}
	}

    /// process option
    // {{{
	for (map<string,string>::iterator it=opt.begin(); it!=opt.end(); it++) {
		string name  = it->first;
		string value = it->second;
		if (name.compare("out_prefix")==0 ) {
			out_prefix = value;
		} else if (name.compare("tmp_dir")==0 ) {
			tmp_dir = value;
		} else if ( name=="log_file" || name=="log" ) {
			log_file = value;
        } else {
		}
	}
	if ( opt.find("out_prefix")==opt.end() ) {
		opt["out_prefix"] = out_prefix;
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

    if (log_file.empty()) {
        log_file = out_prefix + "log.txt";
		opt["log_file"] = log_file;
    }

    if (tmp_dir.empty()) {
        tmp_dir = out_prefix + "tmp";
		opt["tmp_dir"] = tmp_dir;
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
    // }}}
	///
    Msg log(log_file);
	for (map<string,string>::iterator it=opt.begin(); it!=opt.end(); it++) {
		string name  = it->first;
		string value = it->second;
		log << name << " = " << value << "\n";
	}
	log << endl;

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
