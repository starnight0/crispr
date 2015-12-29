#include "czl_common.hpp"

using namespace std;
using namespace czl_bio;

const string ENV[ 24 ] = {
    "COMSPEC", "DOCUMENT_ROOT", "GATEWAY_INTERFACE",   
    "HTTP_ACCEPT", "HTTP_ACCEPT_ENCODING",             
    "HTTP_ACCEPT_LANGUAGE", "HTTP_CONNECTION",         
    "HTTP_HOST", "HTTP_USER_AGENT", "PATH",            
    "QUERY_STRING", "REMOTE_ADDR", "REMOTE_PORT",      
    "REQUEST_METHOD", "REQUEST_URI", "SCRIPT_FILENAME",
    "SCRIPT_NAME", "SERVER_ADDR", "SERVER_ADMIN",      
    "SERVER_NAME","SERVER_PORT","SERVER_PROTOCOL",     
    "SERVER_SIGNATURE","SERVER_SOFTWARE" }; 

string decode_url(string const & src);
string encode_url(string const & src);

int main(int argc, char * argv[])
{
    int i, j, k;
    string line, str, in_file, out_file;
    ifstream fin;
    ofstream fout;

    Opts opts;
//    string & bin_dir = opts.create<OptS>("bin_dir,bin", "DIR", 
//            "directory of 'search_pam' [/usr/bin/]",
//            "/usr/bin/", 1, false)->get_value_ref();
//    string & tmp_dir = opts.create<OptS>("tmp_dir,tmp", "DIR", 
//            "directory of 'search_pam' [/tmp/]",
//            "/tmp/", 1, false)->get_value_ref();
//    string & index_prefix = opts.create<OptS>("index_prefix,index", "STR", 
//            "prefix of 'search_pam' index ",
//            "", 1, true)->get_value_ref();
    opts.parse(argc, argv);
    if ( !opts.is_all_set() ) { pthread_exit( (void*)-1 ); }

    string bin_dir;
    string tmp_dir;
    string index_NGG;
    string index_TTTN;
    string pam;
    char pam_pos;
    int target_length;

    // read from config 'parse_web_input.conf'
    fin.open("parse_web_input.conf");
    if ( fin.fail() ) {
        cerr<< "E: fail to open parse_web_input.conf " << CZL_DBG_INFO << endl;
        cout<< "<p>Fail to run search</p>";
        pthread_exit( (void*)-1);
    }
    while (!fin.eof()) {
        getline(fin, line);
        i = line.find('=', 0);
        string name = line.substr(0, i);
        string value = line.substr(i+1);
        StringUtility::trim(name, " \t\r\n");
        StringUtility::trim(value, " \t\r\n");
        if ( name == "index_NGG" ) {
            index_NGG = value;
        } else if ( name == "bin_dir" ) {
            bin_dir = value;
        } else if ( name == "tmp_dir" ) {
            tmp_dir = value;
        }
    }
    fin.close();

    cout << "Content-type:text/html\r\n\r\n";
    cout << "<!DOCTYPE html>\n";
    cout << "<html>";
    cout << "<head>";
    cout << "<meta charset=\"UTF-8\">";
    cout << "<title>Result</title>";
    cout << "</head>";
    cout << "<body>\n";
//    cout << "<table border=\"1\" cellspacing=\"2\">";
//    for ( int i = 0; i < 24; i++ )
//    {
//        cout << "<tr><td>" << ENV[ i ] << "</td><td>";
//        // attempt to retrieve value of environment variable
//        char *value = getenv( ENV[ i ].c_str() );  
//        if ( value != 0 ){
//            cout << value;                                 
//        }else{
//            cout << "Environment variable does not exist.";
//        }
//        cout << "</td></tr>\n";
//    }
    
    string data;
    cin >> data;

    i=0;
    in_file.clear();
    int pam_len;
    while ( i!=string::npos ) {
        i = StringUtility::find(data, '&', i, str);
        j = str.find('=', 0);
        string name, value;
        name = str.substr(0, j);
        value = str.substr(j+1);
        value = decode_url(value);

        if ( name=="iseq" ) {
            if ( in_file.empty() ) {
                int rn = rand();
                out_file = tmp_dir + "/search_pam.in." + ltos(rn) + ".fa";
                fout.open(out_file.c_str());
                fout << value;
                fout.close();
                in_file = out_file;
            }
        } else if ( name=="pam" ) {
            pam = value;
        } else if ( name=="pam_pos" ) {
            pam_pos = value[0];
        } else if ( name=="target_length" ) {
            target_length = atol(value.c_str());
//        } else if ( name=="index" ) {
//            index = value;
//        } else if ( name=="bin_dir" ) {
//            bin_dir = value;
//        } else if ( name=="tmp_dir" ) {
//            tmp_dir = value;
        }
    }
    pam_len = pam.size();

//    cout << "<table>";
//    cout << "<tr><td>pam</td><td>" <<  pam << "</td></tr>";
//    cout << "<tr><td>pam_pos</td><td>" <<  pam_pos << "</td></tr>";
//    cout << "<tr><td>target_length</td><td>" <<  target_length << "</td></tr>";
//    cout << "<tr><td>tmp_dir</td><td>" <<  tmp_dir << "</td></tr>";
//    cout << "<tr><td>bin_dir</td><td>" <<  bin_dir << "</td></tr>";
//    cout << "<tr><td>index</td><td>" <<  index_NGG << "</td></tr>";

    string cmd = bin_dir + "/search_pam --pam " + pam + " --pam-pos " + pam_pos
        + " --tl " + itos(target_length + pam_len) + " -i " + in_file
        + " --of html ";
    if ( pam=="NGG" && pam_pos=='3' ) {
        cmd += " --index " + index_NGG;
	} else if ( pam=="TTTN" && pam_pos=='5' ) {
        cmd += " --index " + index_TTTN;
    }
//    cout << "<tr><td>cmd</td><td>" << cmd << "</td></tr>";
//    cout << "</table><\n";
    system(cmd.c_str());

    remove(in_file.c_str());

    cout << "</body>\n";
    cout << "</html>\n";

    return 0;
}

string decode_url(string const & src)
{
    string out;
    int j=0;
    for (int i=0; i<src.size(); i++) {
        char a=src[i], b;
        if (a=='+') b = ' ';
        else if (a=='%') {
            i++;
            char u[3];
            u[0] = src[i++];
            u[1] = src[i];
            u[2] = '\0';
            b = strtol(u, NULL, 16);
        } else {
            b = a;
        }
        out.push_back(b);
    }
    return out;
}

string encode_url(string const & src)
{
    stringstream out_s;
    int j=0;
    for (int i=0; i<src.size(); i++) {
        char a=src[i], b;
        switch(a) {
            case ' ':
                out_s << '+';
                break;
            case '&':
            case '\r':
            case '\n':
            case '/':
            case '?':
            case '=':
                out_s << setw(2) << hex << a;
            default:
                out_s << a;
                break;
        }
    }
    return out_s.str();
}
