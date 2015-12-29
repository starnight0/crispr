#ifndef BIT2_HPP
#define BIT2_HPP

#include "czl_common.hpp"

using namespace czl_bio;

class C {
public:
    C(): n(0), seq(0) {}
    C(uint64_t seq0, int n0): seq(seq0), n(n0) {} 
    C(C const & src): seq(src.seq), n(src.n) {}
    C & operator =(C const & src) {
        seq = src.seq;
        n   = src.n;
		return (*this);
    }

    uint64_t seq;
    int n;
};

int C_read(istream & s, C & data);
int C_write(ostream & s, C const & data);
bool C_less(C const & a, C const & b);

void seti64(uint64_t * seed, uint8_t i, uint64_t a);
uint8_t geti64(uint64_t seed, uint8_t i);
uint64_t rev64(uint8_t l, uint64_t seed);
void print64(uint8_t l, uint64_t seed);
ostream & print64(ostream & out, uint8_t l, uint64_t seed);
void str_to_bit2_64(string const & seq, int *len, uint64_t *bseq);
void bit2_to_str_64(int len, uint64_t bseq, string & seq);
bool less64(int l, uint64_t const & a, uint64_t const & b);
uint8_t geti8(uint8_t seed, uint8_t i);
void seti8(uint8_t *seed, uint8_t i, uint8_t a);
void seti8_next(uint8_t *seed, uint8_t i);

class Node {
public:
    typedef enum { None=0, Cont=1, Trie=2} Type;

    Node(Type t = None): n(0), child_n(0), child_sz(0), child_v(NULL), child_min(255), child_max(0), child_unit(0)
    {
        if ( t==Trie ) {
            child_sz = -ab_sz;
            child_v = new uint8_t[ab_sz*sizeof(Node*)];
            for ( int i=0; i<ab_sz*sizeof(Node*); i++ ) {
                child_v[i] = 0;
            }
        }
    }

    Node(Type t, int seq_len): n(0), child_n(0), child_sz(0), child_v(NULL), child_min(255), child_max(0), child_unit((seq_len+3 >>2)+1)
    {
        if ( t==Trie ) {
            child_sz = -ab_sz;
            child_v = new uint8_t[ab_sz*sizeof(Node*)];
            for ( int i=0; i<ab_sz*sizeof(Node*); i++ ) {
                child_v[i] = 0;
            }
        }
    }

    ~Node() {
        if (child_v!=NULL) {
            delete []child_v;
            child_v = NULL;
        }
    }

    Node * get_trie_child(uint8_t byte)
    {
        Node *p = 0;
        long l = 0;
        int k = byte*sizeof(Node*);
        for ( int i=sizeof(p)-1; i>=0; i-- ) {
            l <<= 8;
            l |= child_v[k+i];
        }
        p = (Node*)l;
    //  memcpy( &p, child_v+byte*sizeof(Node*), sizeof(p) );
        return p;
    }

    int get_bseq(int i, int * len, uint64_t * bseq, uint8_t * count)
    {
        if ( i<0 || i>=child_n ) {
            throw out_of_range("i is out of range");
            return -1;
        }
        (*len) = (child_unit-1)<<2;
        int k = i*child_unit;
        int j = k+child_unit-1;
        (*count) = child_v[j--];
        (*bseq) = 0;
        while (j>=0) {
            (*bseq) <<=8;
            (*bseq) |= child_v[k+j];
            j--;
        }
        return 0;
    }

    uint8_t * get_a_unit(int i)
    {
        if ( i<0 || i>=child_n ) {
            throw out_of_range("i is out of range");
            return NULL;
        }
        return child_v+i*child_unit;
    }

    void reset_trie_child(uint8_t byte)
    {
        set_trie_child(byte, NULL);
    }

    inline int insert_trie_child(uint8_t byte, Node * p)
    {
        return insert_child(byte, p);
    }

    int insert_child(uint8_t byte, Node * p)
    {
        if ( get_trie_child(byte)!=NULL ) {
            string msg("Child is not NULL");
            msg += CZL_DBG_INFO;
            throw runtime_error(msg);
        }
        set_trie_child( byte, p);
        child_n++;
        if ( byte < child_min ) child_min = byte;
        if ( byte > child_max ) child_max = byte;

		return 0;
    }

    int insert_bseq(int len, uint64_t bseq, int count)
    {
        int byte_n = (len+3)>>2;

        assert(byte_n == child_unit-1);

        int m = child_n*child_unit;
        if (child_n >= child_sz) {
            if ( child_sz==0 ) child_sz = 1;
            else child_sz <<= 1;
            uint8_t *p = child_v;
            child_v = new uint8_t[child_sz*child_unit];
            for ( int i=0; i < m; i++ ) {
                child_v[i] = p[i];
            }
            delete []p;
        }
        if ( child_unit>1 ) {
            uint8_t byte = bseq&0xff;
            if ( byte < child_min ) child_min = byte;
            if ( byte > child_max ) child_max = byte;
        }
        int i;
        for ( i=0; i<child_unit-1; i++ ) {
            child_v[m+i] = bseq&0xff;
            bseq>>=8;
        }
        child_v[m+i] = count;
        child_n++;

		return 0;
    }

    int insert_a_unit(uint8_t * p)
    {
        int m = child_n*child_unit;
        int i;
        if (child_n >= child_sz) {
            if ( child_sz==0 ) child_sz = 1;
            else child_sz <<= 1;
            uint8_t *p1 = child_v;
            child_v = new uint8_t[child_sz*child_unit];
			if (p1!=NULL) {
				for ( i=0; i < m; i++ ) {
					child_v[i] = p1[i];
				}
				delete []p1;
			}
        }
        for ( i=0; i<child_unit; i++ ) {
            child_v[m+i] = p[i];
        }
        child_n++;
        if (child_unit>1) {
            uint8_t byte = p[0];
            if ( byte < child_min ) child_min = byte;
            if ( byte > child_max ) child_max = byte;
        }

		return 0;
    }

    bool is_container()
    {
        if ( child_sz>=0 ) return true;
        else return false;
    }

    int merge();

    int burst();
//  void add_child(int i);
//  int read(ifstream & f);
//  int write(ofstream & f);
//  int child_n() { return child_v.size(); }
//  void print() { 
//      print64(4, b);
//      cout << ", " << (int)n << ", " << child_v.size();
//      for (int i=0; i<child_v.size(); i++) {
//          cout << child_v[i] << ", ";
//      }
//  }
    /*
    void print() { 
        for (int i=0; i<256; i++) {
            if ( child_v[i]!=NULL ) {
                print64(4, i);
            }
            cout << endl;
        }
    }
    */

//  uint8_t b;  // current byte (4-bp DNA)
    uint8_t n;  // total occurrence
//  int des_n;  // total of descendants
//  vector<uint8_t> container;
    uint8_t child_min;
    uint8_t child_max;
    int16_t child_n;   // number of child in this node or container
    int16_t child_sz;  // size of child_v, -256 for trie node
    uint8_t child_unit; // how many bytes per unit in child_v for container.
    uint8_t* child_v;

    static const int ab_sz = 256;

private:
    void set_trie_child(uint8_t byte, Node *p) throw(runtime_error)
    {
        int i;
        int k = byte*sizeof(p);
        long l = (long)p;
        for ( i=0; i<sizeof(p); i++ ) {
            child_v[k+i] = l&0xff;
            l >>= 8;
        }
    //  memcpy( child_v+byte*sizeof(Node*), &p, sizeof(p) );
    }

    void swap_content(Node & a)
    {
        uint8_t t;
        t = n; n=a.n; a.n=t;
        t = child_min; child_min=a.child_min; a.child_min=t;
        t = child_max; child_max=a.child_max; a.child_max=t;
        t = child_unit; child_unit=a.child_unit; a.child_unit=t;
        uint16_t t1;
        t1 = child_n; child_n=a.child_n; a.child_n=t1;
        t1 = child_sz; child_sz=a.child_sz; a.child_sz=t1;
        uint8_t *p;
        p = child_v; child_v=a.child_v; a.child_v=p;
    }
};

// Node* build_trie(int n, uint64_t *bseq_v, vector<Node*> & node_v);
Node* build_trie_insert(int seq_len, uint64_t bseq, Node * root);
void destroy_trie(Node * root);
void container_sort(Node * root);
void write_trie(ofstream & fout, Node * root, uint8_t b);
Node * read_trie(ifstream & fin, uint8_t *b = NULL);
void print_trie_to_file(string const & file, Node *root);
void print_trie_to_file(char const * file, Node *root);
void print_trie(Node *root);
void print_trie(ostream & out, Node *root);
void print_trie(ostream & out, Node *root, uint8_t b, int level);

class Less64 {
public:
    Less64(int i0, int i1): i0_m(i0), i1_m(i1) {}
    ~Less64() {}
    bool operator ()(uint64_t const & a, uint64_t const & b) {
        uint64_t a0 = a;
        uint64_t b0 = b;
        a0 >>= (i0_m<<1);
        b0 >>= (i0_m<<1);
        for (int i=i0_m; i<i1_m; i++) {
            if ( (a0&0x3)<(b0&0x3) ) return true;
            else if ( (a0&0x3)>(b0&0x3) ) return false;
            a0 >>=2;
            b0 >>=2;
        }
        return false;
    }
    int i0_m;
    int i1_m;
};

class Less64o {
public:
    Less64o(vector<uint64_t> * bseq, int i0, int i1):
            bseq_m(bseq), i0_m(i0), i1_m(i1) {}

    ~Less64o() {}

    bool operator ()(int const & a, int const & b) {
        uint64_t a0 = (*bseq_m)[ a ];
        uint64_t b0 = (*bseq_m)[ b ];
        a0 >>= (i0_m<<1);
        b0 >>= (i0_m<<1);
        for (int i=i0_m; i<i1_m; i++) {
            if ( (a0&0x3)<(b0&0x3) ) return true;
            else if ( (a0&0x3)>(b0&0x3) ) return false;
            a0 >>=2;
            b0 >>=2;
        }
        return false;
    }

    vector<uint64_t> *bseq_m;
    int i0_m;
    int i1_m;
};

class ContainerLess {
public:
    ContainerLess(uint8_t *container, int unit): container_m(container),
            unit_m(unit) {}

    ~ContainerLess() {}

    bool operator()(int a, int b) {
        a*=unit_m;
        b*=unit_m;
        for (int i=0; i<unit_m-1; i++) {
            if ( container_m[a+i] < container_m[b+i]) return true;
            else if ( container_m[a+i] > container_m[b+i] ) return false;
        }
        return false;
    }

    uint8_t * container_m;
    int unit_m;
};


int mk_sort(int seq_len, vector<uint64_t> & seq_v);
int mk_sort(int seq_len, vector<uint64_t> & seq_v, int k);
int mk_sort_r(int seq_len, vector<uint64_t> & seq_v, int i0, int i1, int k);
int mk_sort(int bseq_len, vector<uint64_t> & bseq_v, int k, vector<int> & order);
int mk_sort_r(int bseq_len, vector<uint64_t> & bseq_v, int i0, int i1, int k, vector<int> & order);

int search(Node *root, string const & seq, uint8_t mis_v[0x10000], int max_mis, int score[5]);
int search(Node * root, int seq_len, uint64_t bseq, uint8_t mis_v[0x10000], int max_mis, int mis0, int score[5]);
#endif
