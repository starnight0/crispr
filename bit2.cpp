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

int C_read(istream & s, C & data)
{
    data.n=0;
    if (s.eof()) return -1;
    s.read( (char*)&data.seq, sizeof(data.seq));
    s.read( (char*)&data.n, sizeof(data.n));
    return 0;
}

int C_write(ostream & s, C const & data)
{
    s.write( (char*)&data.seq, sizeof(data.seq));
    s.write( (char*)&data.n, sizeof(data.n));
    return 0;
}

bool C_less(C const & a, C const & b)
{
    return a.seq < b.seq;
}

uint8_t geti64(uint64_t seed, uint8_t i)
{
    return seed>>(i<<1) &0x3;
}

void seti64(uint64_t *seed, uint8_t i, uint64_t a)
{
    (*seed) &= ~(0x3LU<<(i<<1));
    (*seed) |= (a&0x3)<<(i<<1);
}

uint64_t rev64(uint8_t l, uint64_t seed)
{
    uint64_t r=0;
    for (uint8_t i=0; i<l; i++) {
        uint8_t a = geti64(seed, i);
        seti64(&r, l-1-i, a);
    }
    return r;
}

void print64(uint8_t l, uint64_t seed)
{
    for (uint8_t i=0; i<l; i++) {
        cout << Bit2Seq::bit2_to_nt1( geti64(seed, i) );
    }
    cout << flush;
}

ostream & print64(ostream & out, uint8_t l, uint64_t seed)
{
    for (uint8_t i=0; i<l; i++) {
        out << Bit2Seq::bit2_to_nt1( geti64(seed, i) );
    }
    return out;
}

void str_to_bit2_64(string const & seq, int *len, uint64_t *bseq)
{
    (*len) = seq.size();
    uint64_t a;
    (*bseq) = 0;
    for (short i=0; i<(*len); i++) {
        a = Bit2Seq::nt1_to_bit2(seq[i]);
        (*bseq) |= (a&0x3)<<(i<<1);
    }
}

void bit2_to_str_64(int len, uint64_t bseq, string & seq)
{
    seq.resize(len);
    for (uint8_t i=0; i<len; i++) {
        seq[i] = Bit2Seq::bit2_to_nt1( geti64(bseq, i) );
    }
}

bool less64(int l, uint64_t const & a, uint64_t const & b)
{
    uint64_t a0 = a;
    uint64_t b0 = b;
    for (int i=0; i<l; i++) {
        if ( (a0&0x3)==(b0&0x3) ) continue;
        else if ( (a0&0x3)<(b0&0x3) ) return true;
        else return false;
        a0 >>=2;
        b0 >>=2;
    }
    return false;
}

uint8_t geti8(uint8_t seed, uint8_t i)
{
    return seed>>(i<<1) & 0x3;
}

void seti8(uint8_t *seed, uint8_t i, uint8_t a)
{
    (*seed) &= ~(0x3U<<(i<<1));
    (*seed) |= (a&0x3)<<(i<<1);
}

void seti8_next(uint8_t *seed, uint8_t i)
{
    uint8_t a=(*seed)>>(i<<1) &0x3;
    a = (a+1)&0x3;
    (*seed) &= ~(0x3U<<(i<<1));
    (*seed) |= (a&0x3)<<(i<<1);
}
/*
Node* build_trie(int n, uint64_t *bseq_v, vector<Node*> & node_v)
{
    if (n<=0) return NULL;

    int i, j;
    Node * root = new Node();
    node_v.assign(1,root);
    uint64_t bseq = bseq_v[0];
    // add the first sequence
    int bl = 24;  // 24 * 2 bits
    Node *p=root;
    cout << "Add First : ";
    print64(24, bseq);
    cout << endl;
    for (i=0; i<6; i++) {
        uint8_t b = bseq&0xff;
        bseq >>= 8;
        Node * a = new Node();
        a->b = b;
        a->n = 1;
        
        int j = node_v.size();
        p->child_v.push_back(j);

        node_v.push_back(a);
        p = a;
    }
    root->n = 1;
    for (int k=1; k<n; k++) {
        uint64_t bseq = bseq_v[k];
        p = root;
        for (i=0; i<6; i++) {
            uint8_t b = bseq&0xff;
            bseq >>= 8;
            // bisearch the position of b in parent's child list
            int lo = 0, hi = p->child_n(), m=(lo+hi)>>1;
            Node * p0 = node_v[ p->child_v[lo] ];
            Node * p2 = node_v[ p->child_v[hi-1] ];
            Node * p1;
            while ( lo<hi ) {
                m = (lo+hi)>>1;
                p1 = node_v[ p->child_v[m] ];
                if ( p1->b == b ) break;
                else if ( b < p1->b ) {
                    hi = m;
                } else {
                    lo = m+1;
                }
            }
            if ( lo<hi ) {
                if ( p1->n != 255U ) p1->n++;
                p = p1;
            } else { // Can't find b, create a new node
                Node * a = new Node();
                a->b = b;
                a->n = 1;

                int j = node_v.size();
                p->child_v.insert(p->child_v.begin()+lo, j);
                node_v.push_back(a);
                p = a;
                i++;

                while ( i<6 ) {
                    uint8_t b = bseq&0xff;
                    bseq >>= 8;
                    Node * a = new Node();
                    a->b = b;
                    a->n = 1;

                    int j = node_v.size();
                    p->child_v.push_back(j);

                    node_v.push_back(a);
                    p = a;
                    i++;
                }
            }
        }
    }
}
*/

// !REMARK: at last level, the child_v store the occurrence of the last character
Node* build_trie_insert(int seq_len, uint64_t bseq, Node * root)
{
    int i, j, k;
    Node *here=NULL, *parent=NULL;
    int level0 = (seq_len+3)>>2, level=0;
    int M = 1000;
    uint8_t byte;
    int unit = 0;
    int seq_len0 = seq_len;
    {
        byte = bseq & 0xff;
        level = 0;
        parent = NULL;
        here = root;
        // search the trie locate the place to put the sequence
        // store the location in 'here'
        while (seq_len>0 && here!=NULL && !here->is_container()) {
            byte = bseq & 0xff;
            bseq >>= 8;
            seq_len -= 4;
            if (here->n < 255) here->n++;
            parent = here;
            // move to next node
            here = here->get_trie_child(byte);
            level++;
        }
        {
            if ( here==NULL ) {
                // create a new node
                here = new Node(Node::Cont, seq_len);
                //
                // set it as a container
                //
                if ( seq_len == seq_len0 ) {
                    // set 'here' as root
                    root = here;
                } else {
                    // link it to its parent
                    parent->insert_child(byte, here);
                }
            }
            // now 'here' is a container
            // add bseq (from level 'level') into the container
            if ( seq_len>0 ) {
                here->insert_bseq(seq_len, bseq, 1);
            }

            if (here->n < 255) here->n++;

            // burst container if container is too large
            // but don't burst last level
            if ( here->child_n > M ) {
                here->burst();
            }
            //
        }
    }
    return root;
}

void destroy_trie(Node * root)
{
    if ( root!=NULL ) {
        if ( root->is_container() ) {
            delete root;
        } else {
            for (int i=0; i<256; i++) {
                Node * p = root->get_trie_child(i);
                if ( p != NULL) {
                    destroy_trie( p );
                    root->reset_trie_child(i);
                }
            }
        }
    }
}

void container_sort(Node * root)
{
    if ( root->child_unit==1 ) return;
    if ( !root->is_container() ) {
        for (int i=0; i<256; i++) {
            Node * p = root->get_trie_child(i);
            if ( p!=NULL ) {
                container_sort(p);
            }
        }
    } else {
        root->merge();
    }
}

/*
int cal_des_n(vector<Node*> node_v, int root_id)
{
    Node * root = node_v[root_id];
    if ( root->child_n()==0 ) {
        root->des_n = 1;
    } else {
        root->des_n = 0;
        for ( int j=0; j<root->child_n(); j++) {
            root->des_n += cal_des_n(node_v, root->child_v[j]);
        }
    }
    return root->des_n;
}
*/

void print_trie_to_file(char const * file, Node *root)
{
    ofstream fout(file);
    print_trie(fout, root);
    fout.close();
}

void print_trie_to_file(string const & file, Node *root)
{
    ofstream fout(file.c_str());
    print_trie(fout, root);
    fout.close();
}

void print_trie(ostream & out, Node *root)
{
    print_trie(out, root, 0, 0);
}

void print_trie(Node *root)
{
    print_trie(cout, root, 0, 0);
}

void print_trie(ostream & out, Node *root, uint8_t b, int level)
{
    string s;
    for (int i=0; i< level; i++) {
        s.push_back(' ');
        s.push_back(' ');
    }
    out << s;
    if ( root != NULL ) {
        if ( !root->is_container() ) {
            out << "T";
        } else {
            out << "C";
        }
        out << "(";
        print64( out, 4, (uint64_t)b );
        out<< ", " << (int)root->n << ", " << root->child_n << ", "
            << (int)root->child_unit << ", ";
        if ( !root->is_container() ) {
            out << "\n";
            if ( root->child_unit>1 ) {
                for (int i=0; i<256; i++) {
                    Node * p = root->get_trie_child(i);
                    if ( p!=NULL ) {
                        print_trie(out, p, (uint8_t)i, level+1);
                    }
                }
            }
        } else {
            out << "\n";
            for (int i=0; i<root->child_n; i++) {
                int len;
                uint8_t c;
                uint64_t bseq;
                root->get_bseq(i, &len, &bseq, &c);
                out << s << "  ";
                print64(out, len, bseq);
                out << ", " << (int)c << "\n";
            }
        }
    }
    out << flush;
}

void write_trie(ofstream & fout, Node * root, uint8_t b)
{
    if ( root != NULL ) {
        fout.write( (char*)&b, sizeof(b) );
        fout.write( (char*)&root->n, sizeof(root->n));
        uint8_t unit = root->child_unit;
        fout.write( (char*)&unit, sizeof(unit));
        uint16_t m = root->child_n;
        if ( !root->is_container() ) {
            fout.write( (char*)&m, sizeof(m));
            if ( m>0 ) {
                for (int i=0; i<256; i++) {
                    Node * p = root->get_trie_child(i);
                    if ( p!=NULL ) {
                        write_trie(fout, p, (uint8_t)i);
                    }
                }
            }
        } else {
            m |= 0x8000;
            fout.write( (char*)&m, sizeof(m));
//            fout << flush;
            if ( root->child_unit>1 ) {
                fout.write( (char*)root->child_v, root->child_n*root->child_unit );
//                fout << flush;
            }
        }
    }
}

Node * read_trie(ifstream & fin, uint8_t *b)
{
    uint8_t n;
    uint16_t m;
    uint8_t bb;
    uint8_t unit;
    fin.read( (char*)&bb, sizeof(bb));
    if ( b!=NULL ) {
        (*b) = bb;
    }
    fin.read( (char*)&n, sizeof(n));
    fin.read( (char*)&unit, sizeof(unit));
    fin.read( (char*)&m, sizeof(m));
    Node * p = NULL;
    if ( m&0x8000 ) { // is container
        try {
            p = new Node(Node::Cont);
        } catch ( exception & e ) {
            cerr << "Exception: " << e.what() << " " << CZL_DBG_INFO << endl;
        }
        p->n = n;
        p->child_unit = unit;
        uint16_t m0 = m&0x7fff;
        if ( unit>1 ) {
            uint8_t * u = new uint8_t[unit+1];
            for (int i=0; i<m0; i++) {
                fin.read( (char*)u, unit );
                if (fin.fail() ) {
                    cerr << "E: Fail to read " << CZL_DBG_INFO << endl;
                    pthread_exit( (void*)-1 );
                }
                p->insert_a_unit(u);
            }
            delete []u;
        }
    } else {
        try {
            p = new Node(Node::Trie);
        } catch ( exception & e ) {
            cerr << "Exception: " << e.what() << " " << CZL_DBG_INFO << endl;
        }
        p->child_unit = unit;
        p->n = n;

        for (int i=0; i<m; i++) {
            uint8_t b1;
            Node * p1 = read_trie(fin, &b1);
            p->insert_child(b1, p1);
        }
    }
    return p;
}

// multikey-sort
// {{{
int mk_sort(int seq_len, vector<uint64_t> & seq_v)
{
    return mk_sort_r(seq_len, seq_v, 0, seq_v.size(), 0);
}

int mk_sort(int seq_len, vector<uint64_t> & seq_v, int k)
{
    return mk_sort_r(seq_len, seq_v, 0, seq_v.size(), k);
}


int mk_sort_r(int seq_len, vector<uint64_t> & seq_v, int i0, int i1, int k)
{
    int i, j;
    const int K=4;
    int sz[K] = {0, 0, 0, 0};
    int b[K] = {0,0,0,0};
    int e[K] = {0,0,0,0};
    int k2 = k<<1;
    for (i=i0; i<i1; i++) {
        uint8_t a = seq_v[i]>>k2 &0x3;
        sz[a]++;
    }
    b[0] = i0;
    e[0] = i0+sz[0];
    for (i=1; i<K; i++) {
        b[i] = b[i-1]+sz[i-1];
        e[i] = b[i]+sz[i];
    }
    for (j=0; j<K; j++) {
        int b0 = b[j];
        for (i=b[j]; i<e[j]; ) {
            uint8_t a = seq_v[i]>>k2 &0x3;
            if ( a!=j ) {
                uint64_t t = seq_v[i];
                seq_v[i] = seq_v[ b[a] ]; 
                seq_v[ b[a] ] = t;
            } else {
                i++;
            }
            b[a]++;
        }
    }
    if ( k<seq_len-1 ) {
        for (j=0; j<K; j++) {
            int b0;
            if (j==0) b0=i0;
            else b0 = e[j-1];
            if ( e[j]-b0>100 ) {
                mk_sort_r(seq_len, seq_v, b0, e[j], k+1);
            } else {
                sort( seq_v.begin()+b0, seq_v.begin()+e[j], Less64(k+1, seq_len));
            }
        }
    }
    return 0;
}

/**
 * @brief sort sequence from start position k
 *        not change original position, order is output in  'order'
 */
int mk_sort(int bseq_len, vector<uint64_t> & bseq_v, int k, vector<int> & order)
{
    order.resize(bseq_v.size());
    for (int i=0; i<bseq_v.size(); i++) {
        order[i] = i;
    }
    return mk_sort_r(bseq_len, bseq_v, 0, bseq_v.size(), k, order);
}

/**
 * @brief sort the last k-1 mer
 */
int mk_sort_r(int bseq_len, vector<uint64_t> & bseq_v, int i0, int i1, int k, vector<int> & order)
{
    int i, j;
    const int K=4;
    int sz[K] = {0, 0, 0, 0};
    int b[K] = {0,0,0,0};
    int e[K] = {0,0,0,0};
    int k2 = k<<1;
    for (i=i0; i<i1; i++) {
        uint8_t a = bseq_v[ order[i] ]>>k2 &0x3;
        sz[a]++;
    }
    b[0] = i0;
    e[0] = i0+sz[0];
    for (i=1; i<K; i++) {
        b[i] = b[i-1]+sz[i-1];
        e[i] = b[i]+sz[i];
    }
    for (j=0; j<K; j++) {
        int b0 = b[j];
        for (i=b[j]; i<e[j]; ) {
            uint8_t a = bseq_v[ order[i] ]>>k2 &0x3;
            if ( a!=j ) {
                int t = order[i];
                order[i] = order[ b[a] ];
                order[ b[a] ] = i;
            } else {
                i++;
            }
            b[a]++;
        }
    }
    if ( k<bseq_len-1 ) {
        for (j=0; j<K; j++) {
            int b0;
            if (j==0) b0=i0;
            else b0 = e[j-1];
            if ( e[j]-b0>100 ) {
                mk_sort_r(bseq_len, bseq_v, b0, e[j], k+1, order);
            } else {
                sort( order.begin()+b0, order.begin()+e[j], Less64o(&bseq_v, k+1, bseq_len));
            }
        }
    }
    return 0;
}
// }}}

// class Node
// {{{
int Node::merge()
// {{{
{
    if ( !is_container() ) return 0;
    if ( child_n==0 ) return 0;

    vector<int> order;
    int j=0, i, k;
    order.resize(child_n);
    for ( i=0; i<child_n; i++) {
        order[i] = i;
    }
    sort(order.begin(), order.end(),
            ContainerLess(child_v, child_unit));
    Node * p = new Node(Node::Cont);
    p->child_unit = child_unit;
    uint8_t c;
    uint8_t *u, *u1;
    u1 = get_a_unit(order[0]);
    p->insert_a_unit(u1);
    u = p->get_a_unit(0);
    // merge identical sequence and count
//  j=0;
    for (i=1; i<child_n; i++) {
        u1 = get_a_unit(order[i]);
        for (k=0; k<child_unit-1; k++) {
            if ( u1[k] != u[k] ) break;
        }
        if (k < child_unit-1) {
            p->insert_a_unit(u1);
            u = p->get_a_unit(p->child_n-1);
        } else {
            if ( u[child_unit-1]<255 ) u[child_unit-1]++;
        }
    }

    delete child_v;
    child_sz = p->child_sz;
    child_v = p->child_v;
    child_n = p->child_n;
    p->child_v = NULL;
    delete p;

    return 0;
}
// }}}

int Node::burst()
// {{{
{
    if ( !is_container() ) return 0;

    int unit = this->child_unit;
    if ( unit<=1 ) return 0;

    Node * child=NULL;
    // if the first byte is identical 
    //   create new node for this byte and recurse
    int k=0, i, j;
    Node * parent = this;
    uint8_t byte;
    while ( parent->child_unit>1 && parent->child_min == parent->child_max ) {
        byte = parent->child_min;
        child = new Node(Trie);
        child->child_unit = parent->child_unit;
        child->n = parent->n;
        parent->swap_content(*child);

        parent->child_min = child->child_min;
        parent->child_max = child->child_max;
        child->child_min = 255;
        child->child_max = 0;
        child->child_unit--;
        parent->insert_child(byte, child);

        k++;
        if ( child->child_unit>1 ) {
            for ( i=k; i<child->child_n*unit; i+=unit ) {
                byte = child->child_v[i];
                if ( byte > child->child_max ) child->child_max = byte;
                if ( byte < child->child_min ) child->child_min = byte;
            }
        }
        parent = child;
    }
    assert( parent->child_unit>=1 );
    if ( parent->child_unit>1 ) {
        child = new Node(Trie);
        child->n = parent->n;
        child->child_unit = parent->child_unit;
        parent->swap_content(*child);
        // after swap, parent will be a trie node

        Node *child0 = child;
        child0->child_unit = unit;
        for ( i=0; i<child0->child_n; i++ ) {
            byte = child0->child_v[i*unit+k];
            child = parent->get_trie_child(byte);
        //  for ( j=i; j<i+unit-1; j++ ) {
        //      print64(4, (uint64_t)here->container[j]);
        //  }
        //  cout << endl;
            if ( child == NULL ) {
                try {
                    child = new Node(Cont);
                    child->child_unit = parent->child_unit-1;
                } catch (exception & e) {
                    cerr<< "E: " << e.what() << CZL_DBG_INFO << endl;
                    pthread_exit( (void*)-1 );
                }
                parent->insert_child(byte, child);
            }
            uint8_t *p = child0->get_a_unit(i);
            p += k+1;
            if ( child->child_unit>1 ) {
                child->insert_a_unit(p);
            }

            if ( child->n<255 ) child->n++;
        }
        delete child0;
    } else {
        uint8_t *p = parent->child_v;
        parent->child_v = NULL;
        parent->child_sz = 0;
        parent->child_n = 0;
        delete []p;
    }
    return 0;
}
// }}}

// }}}

int search(Node *root, string const & seq, uint8_t mis_v[0x10000], int max_mis, int score[])
{
    int len;
    uint64_t bseq;
    str_to_bit2_64(seq, &len, &bseq);
    return search(root, len, bseq, mis_v, max_mis, 0, score);
}

int search(Node * root, int seq_len, uint64_t bseq, uint8_t mis_v[0x10000], int max_mis, int mis0, int score[])
{
    if ( seq_len <= 4 ) {
        uint16_t mask = 0xff >> ((4-seq_len)<<1);
        uint8_t b0 = bseq & mask;
        if ( root->is_container() ) {
            for ( int i=0; i<root->child_n; i++ ) {
                uint16_t a = root->child_v[i*root->child_unit] & mask;
                uint8_t c = root->child_v[(i+1)*root->child_unit-1];
                a = a<<8 | b0;
                uint8_t mis1 = mis_v[a] + mis0;
                if ( mis1 <= max_mis ) {
                    score[mis1] += c;
                }
            }
        } else {
            for ( int i=root->child_min; i<=root->child_max; i++ ) {
                Node * p = root->get_trie_child(i);
                if ( p == NULL ) continue;

                uint16_t a = (i&mask)<<8 | b0;
                uint8_t mis1 = mis_v[a] + mis0;
                if ( mis1 <= max_mis ) {
                    score[mis1] += p->n;
                }
            }
        }
    } else {
        uint8_t b0 = bseq & 0xff;
        if ( root->is_container() ) {
            int unit = root->child_unit;
            for ( int i=0; i<root->child_n*unit; i+=unit ) {
                int mis1 = mis0;
                uint64_t bseq1 = bseq;
                int seq_len1 = seq_len;
                int j = 0;
                while ( j<unit-1 && seq_len1 >= 4 ) {
                    uint16_t a = root->child_v[i+j];
                    a = a<<8 | bseq1&0xff;
                    mis1 += mis_v[a];
                    if ( mis1 > max_mis ) break;
                    bseq1 >>= 8;
                    seq_len1-=4;
                    j++;
                }
                if ( seq_len1>0 && seq_len1<4 ) {
                    uint16_t mask = 0xff >> ((4-seq_len1)<<1);
                    uint16_t a = (root->child_v[i+j]&mask);
                    a = a<<8 | bseq1&mask;
                    mis1 += mis_v[a];
                }
                if ( mis1 <= max_mis ) {
                    score[mis1] += root->child_v[i+unit-1];
                }
            }
        } else {
            int m = max_mis - mis0;
            if (m==0) {
                Node * p = root->get_trie_child(b0);
                if ( p != NULL ) {
                    search( p, seq_len-4, bseq>>8, mis_v, max_mis, mis0, score );
                }
            } else if (m<4) {
                Node * p = root->get_trie_child(b0);
                if ( p != NULL ) {
                    search( p, seq_len-4, bseq>>8, mis_v, max_mis, mis0, score );
                }
                int k0;
                uint8_t *pos0 = new uint8_t[m];
                for (k0=1; k0<=m; k0++) {
                    // first 8-bp seed with mismatch k0
                    uint8_t n0 = 1;
                    uint8_t b1 = b0;
                    pos0[0] = 0;
                    while ( n0>0 ) {
                        if ( pos0[n0-1] == 4 ) {
                            n0--;
                            if (n0>0) pos0[n0-1]++;
                        } else {
                            if ( n0 < k0 ) {
                                pos0[n0] = pos0[n0-1]+1;
                                n0++;
                            } else {
                                // generate all k0 mis-match seed 
                                // with mismatch position pos1
                                b1 = b0;
                                uint8_t n01 = 1;
                                seti8_next(&b1, pos0[n01-1]);
                                while ( n01>0 ) {
                                    if ( geti8(b1, pos0[n01-1]) == geti8(b0, pos0[n01-1]) ) {
                                        n01--;
                                        if (n01>0) seti8_next(&b1, pos0[n01-1]);
                                    } else {
                                        if (n01==k0) {
                                            // search
                                        //  print64(4, (uint64_t)b1);
                                        //  cout << endl;
                                            Node * p = root->get_trie_child(b1);
                                            if ( p != NULL ) {
                                                search( p, seq_len-4, bseq>>8, mis_v, max_mis, mis0+k0, score );
                                            }

                                            seti8_next(&b1, pos0[n01-1]);
                                        } else {
                                            uint8_t a = geti8(b0, pos0[n01]);
                                            a = (a+1) &0x3;
                                            seti8( &b1, pos0[n01], a);
                                            n01++;
                                        }
                                    }
                                }

                                pos0[n0-1]++;
                            }
                        }
                    }
                //  cout << endl << endl;
                }
                delete []pos0;
            } else {
                for ( int i=root->child_min; i<=root->child_max; i++ ) {
                    Node * p = root->get_trie_child(i);
                    if ( p == NULL ) continue;

                    uint16_t a = i<<8 | b0;
                    if ( mis_v[a] + mis0 > max_mis ) continue;

                    search( p, seq_len-4, bseq>>8, mis_v, max_mis, mis0+mis_v[a], score );
                }
            }
        }
    }
    return 0;
}
