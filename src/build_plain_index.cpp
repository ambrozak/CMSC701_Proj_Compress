// build_plain_fm.cpp
//
// Plain FM-index builder for DNA strings.
// Usage:
//   build_plain_fm <in.fa> <out.idx>
//
// Alphabet:
//   $ = 0 (sentinel)
//   A = 1
//   C = 2
//   G = 3
//   T = 4
//
// Stores:
//   - BWT
//   - C table
//   - Wavelet tree
//   - SA samples
//
// Simple/reference implementation for correctness testing.

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

static constexpr int SIGMA = 3; // enough for symbols 0..4
static constexpr int SA_SAMPLE_RATE = 32;

static inline uint32_t baseToSym(char c){
    switch(c){
        case 'A': return 1;
        case 'C': return 2;
        case 'G': return 3;
        case 'T': return 4;
    }
    return 0;
}

static string loadFirstSeq(const string& fn){
    ifstream f(fn);
    if(!f){
        cerr << "Cannot open " << fn << "\n";
        exit(1);
    }

    string line, seq;
    bool in = false;

    while(getline(f, line)){
        if(line.empty()) continue;

        if(line[0] == '>'){
            if(in) break;
            in = true;
            continue;
        }

        if(in) seq += line;
    }

    return seq;
}

struct WaveletTree{
    size_t n;
    vector<vector<uint64_t>> levels;
};

static WaveletTree buildWT(const vector<uint32_t>& bwt){
    WaveletTree wt;
    wt.n = bwt.size();

    size_t n = bwt.size();
    size_t words = (n + 63) / 64;

    wt.levels.resize(SIGMA);

    vector<uint32_t> perm(n);
    for(size_t i=0;i<n;i++) perm[i]=(uint32_t)i;

    for(int l=SIGMA-1;l>=0;l--){
        wt.levels[l].assign(words,0ULL);

        for(size_t i=0;i<n;i++){
            if((bwt[perm[i]] >> l) & 1ULL)
                wt.levels[l][i/64] |= (1ULL << (i%64));
        }

        vector<uint32_t> L,R;

        for(size_t i=0;i<n;i++){
            if(((bwt[perm[i]] >> l)&1)==0) L.push_back(perm[i]);
            else R.push_back(perm[i]);
        }

        perm.clear();
        perm.insert(perm.end(),L.begin(),L.end());
        perm.insert(perm.end(),R.begin(),R.end());
    }

    return wt;
}

// ---------------------------------------------------------------------------
// SA-IS  (Suffix Array - Induced Sorting)
// ---------------------------------------------------------------------------

static void classifySL(const vector<uint32_t>& t, vector<bool>& isS)
{
    int n = (int)t.size();
    isS.assign(n + 1, false);
    isS[n] = true; // sentinel is S-type
    if(n == 0) return;
    isS[n - 1] = false; // last real char is L-type (sentinel is smaller)
    for(int i = n - 2; i >= 0; i--)
        isS[i] = (t[i] < t[i+1]) || (t[i] == t[i+1] && isS[i+1]);
}

static inline bool isLMS(const vector<bool>& isS, int i)
{
    return i > 0 && isS[i] && !isS[i-1];
}

// Compute bucket head (end=false) or tail (end=true) pointers.
static void getBuckets(const vector<uint32_t>& t, int n,
                       int alpha, vector<int>& bkt, bool end)
{
    bkt.assign(alpha, 0);
    for(int i = 0; i <= n; i++)
        bkt[(i == n) ? 0 : (int)t[i]]++;
    int sum = 0;
    for(int c = 0; c < alpha; c++){
        sum += bkt[c];
        bkt[c] = end ? sum : sum - bkt[c];
    }
}

// Induced-sort L-type suffixes.
// FIX: use sa[i] < 0 (not <= 0) — position 0 is a valid suffix.
static void induceSortL(const vector<uint32_t>& t, int n, int alpha,
                        const vector<bool>& isS, vector<int>& sa)
{
    vector<int> bkt;
    getBuckets(t, n, alpha, bkt, false);
    for(int i = 0; i <= n; i++){
        if(sa[i] < 0) continue;       // -1 = empty slot
        int j = sa[i] - 1;
        if(j >= 0 && !isS[j])         // j >= 0 guards against wrapping past sentinel
            sa[bkt[(int)t[j]]++] = j;
    }
}

// Induced-sort S-type suffixes.
static void induceSortS(const vector<uint32_t>& t, int n, int alpha,
                        const vector<bool>& isS, vector<int>& sa)
{
    vector<int> bkt;
    getBuckets(t, n, alpha, bkt, true);
    for(int i = n; i >= 0; i--){
        if(sa[i] < 0) continue;       // -1 = empty slot
        int j = sa[i] - 1;
        if(j >= 0 && isS[j])
            sa[--bkt[(int)t[j]]] = j;
    }
}

// Compare two LMS substrings for equality.
// FIX: guard i+d and j+d against going out of bounds.
static bool lmsEqual(const vector<uint32_t>& t, int n,
                     const vector<bool>& isS, int i, int j)
{
    if(i == n || j == n) return false;
    for(int d = 0; ; d++){
        if(i+d > n || j+d > n) return false;
        int ti = (i+d == n) ? 0 : (int)t[i+d];
        int tj = (j+d == n) ? 0 : (int)t[j+d];
        if(ti != tj || isS[i+d] != isS[j+d]) return false;
        if(d > 0 && isLMS(isS, i+d) && isLMS(isS, j+d)) return true;
    }
}

// Core SA-IS. t[0..n] with t[n]==0 as unique sentinel (supplied by caller).
// alpha = alphabet size (values in t are in [0, alpha-1]).
static vector<int> sais(const vector<uint32_t>& t, int n, int alpha)
{
    vector<bool> isS;
    classifySL(t, isS);

    // Step 1: place LMS suffixes at bucket tails.
    vector<int> sa(n + 1, -1);
    {
        vector<int> bkt;
        getBuckets(t, n, alpha, bkt, true);
        for(int i = 1; i < n; i++){         // skip i==n (sentinel placed below)
            if(isLMS(isS, i))
                sa[--bkt[(int)t[i]]] = i;
        }
        sa[0] = n;                           // sentinel always sorts first
    }

    // Step 2: induced sort (approximate order).
    induceSortL(t, n, alpha, isS, sa);
    induceSortS(t, n, alpha, isS, sa);

    // Step 3: assign ranks to LMS substrings.
    // Collect LMS positions in their current sorted order (skip sentinel at sa[0]).
    vector<int> lms;
    for(int i = 1; i <= n; i++)
        if(sa[i] >= 0 && isLMS(isS, sa[i])) lms.push_back(sa[i]);

    vector<int> rank(n + 1, -1);
    int curRank = 0, prev = -1;
    for(int cur : lms){
        if(prev == -1 || !lmsEqual(t, n, isS, prev, cur)) curRank++;
        rank[cur] = curRank - 1;
        prev = cur;
    }

    // Build reduced string from LMS positions in left-to-right order.
    vector<int> lmsOrder;
    for(int i = 1; i < n; i++)
        if(isLMS(isS, i)) lmsOrder.push_back(i);

    int n1 = (int)lmsOrder.size();
    vector<uint32_t> t1(n1);
    for(int i = 0; i < n1; i++)
        t1[i] = (uint32_t)rank[lmsOrder[i]];

    // Step 4: recurse or build sa1 directly.
    vector<int> sa1;
    if(curRank < n1){
        // Ranks not unique: shift values up by 1 and append sentinel (0).
        // FIX: was passing t1 without a sentinel, causing OOB in recursive call.
        vector<uint32_t> t1s(n1 + 1);
        for(int i = 0; i < n1; i++) t1s[i] = t1[i] + 1; // shift to free up 0
        t1s[n1] = 0;                                       // sentinel
        sa1 = sais(t1s, n1, curRank + 1);
    } else {
        // Ranks unique: construct sa1 directly in O(n).
        sa1.assign(n1 + 1, -1);
        sa1[0] = n1; // sentinel position
        for(int i = 0; i < n1; i++)
            sa1[t1[i] + 1] = i; // +1 because sa1[0] is reserved for sentinel
    }

    // Step 5: induce the full SA from the correctly sorted LMS order.
    fill(sa.begin(), sa.end(), -1);
    {
        vector<int> bkt;
        getBuckets(t, n, alpha, bkt, true);
        for(int i = n1; i >= 1; i--){      // sa1[1..n1] = sorted LMS indices
            int pos = lmsOrder[sa1[i]];
            sa[--bkt[(int)t[pos]]] = pos;
        }
        sa[0] = n;
    }

    induceSortL(t, n, alpha, isS, sa);
    induceSortS(t, n, alpha, isS, sa);

    return sa;
}

// Public wrapper: builds SA for text[0..n-1] where all values > 0.
// Shifts symbols up by 1 and appends sentinel 0 internally.
// Returns sa[0..n] with sa[0] == n (sentinel suffix).
static vector<int> buildSA(const vector<uint32_t>& text)
{
    int n = (int)text.size();
    vector<uint32_t> t(n + 1);
    for(int i = 0; i < n; i++) t[i] = text[i] + 1; // A=2,C=3,G=4,T=5
    t[n] = 0;                                        // sentinel
    return sais(t, n, 6);                            // alpha=6 (0..5)
}

// BWT[i] = T[SA[i]-1 mod (n+1)], where T[n]=$=0 (sentinel).
// SA[i]==0 -> T[n]=$ ; SA[i]==n -> T[n-1] ; otherwise T[SA[i]-1].
// This gives exactly one $ and one of each real char per bucket, so C
// and rankSym are correct with no manual offset adjustments.
static vector<uint32_t> buildBWT(const vector<int>& sa,
                                 const vector<uint32_t>& text)
{
    int n=(int)text.size();

    vector<uint32_t> bwt;
    bwt.reserve(n + 1);

    for(int i=0;i<=n;i++){
        int pos=sa[i];

        if(pos==0)   bwt.push_back(0);            // T[n] = $ (sentinel char)
        else         bwt.push_back(text[pos-1]);   // includes pos==n: text[n-1]
    }

    return bwt;
}

static vector<uint64_t> buildCTable(const vector<uint32_t>& bwt){
    // C[c] = #{i : BWT[i] < c} = #{rows whose first F-column char < c}.
    // Since BWT now contains exactly one $ (=0) plus all real chars, the
    // prefix-sum of symbol counts gives the correct C table directly.
    vector<uint64_t> C(5+1,0);

    for(uint32_t s:bwt)
        C[s+1]++;

    for(int i=1;i<=5;i++)
        C[i]+=C[i-1];

    return C;
}

// SA sampling by text position.
// We store:
//   - a packed bitvector (one bit per BWT row) marking rows where SA[row]%rate==0
//   - a compacted array of just those SA values, in BWT-row order
// locateRow checks the bitvector to terminate, then uses a popcount rank to
// index into the compacted array. Space is O(n/rate) for the values plus O(n/64)
// for the bitvector, vs O(n) for the old dense approach.
struct SASamples {
    vector<uint64_t> bv;      // packed bitvector, one bit per BWT row 0..n
    vector<uint64_t> vals;    // compacted text positions for sampled rows, in order
    uint64_t         n;       // number of BWT rows (= text length)
};

static SASamples buildSASamples(const vector<int>& sa, uint64_t n)
{
    SASamples s;
    s.n = n;
    s.bv.assign((n + 2 + 63) / 64, 0); // +2: rows 0..n inclusive

    for(int i = 0; i <= (int)n; i++){
        int pos = sa[i];
        if(pos == (int)n) continue;                      // sentinel suffix
        if((size_t)pos % (size_t)SA_SAMPLE_RATE == 0){
            s.bv[i / 64] |= 1ULL << (i % 64);           // mark BWT row i
            s.vals.push_back((uint64_t)pos);
        }
    }
    return s;
}

static void w64(ofstream& f,uint64_t v){f.write((char*)&v,8);}

int main(int argc,char* argv[]){

    if(argc!=3){
        cerr<<"Usage: build_plain_fm <in.fa> <out.idx>\n";
        return 1;
    }

    string raw=loadFirstSeq(argv[1]);

    transform(raw.begin(),raw.end(),raw.begin(),::toupper);

    vector<uint32_t> text;

    for(char c:raw){
        uint32_t s=baseToSym(c);
        if(s) text.push_back(s);
    }

    uint64_t n=text.size();

    cerr<<"Building SA...\n";
    vector<int> sa=buildSA(text);

    cerr<<"Building BWT...\n";
    vector<uint32_t> bwt=buildBWT(sa,text);

    cerr<<"Building C table...\n";
    vector<uint64_t> C=buildCTable(bwt);

    cerr<<"Building WT...\n";
    WaveletTree wt=buildWT(bwt);

    cerr<<"Building SA samples...\n";
    SASamples saSamples=buildSASamples(sa,n);

    ofstream out(argv[2],ios::binary);

    out.write("PLAINFM",7);

    w64(out,n);

    for(int i=0;i<=4;i++)
        w64(out,C[i]);

    // Bitvector: one bit per BWT row 0..n
    w64(out,(uint64_t)saSamples.bv.size());
    for(auto v:saSamples.bv) w64(out,v);

    // Compacted text positions
    w64(out,(uint64_t)saSamples.vals.size());
    for(auto v:saSamples.vals) w64(out,v);

    for(int l=0;l<SIGMA;l++){
        w64(out,(uint64_t)wt.n);

        size_t words=(wt.n+63)/64;

        for(size_t w=0;w<words;w++)
            w64(out,wt.levels[l][w]);
    }

    cerr<<"Done.\n";
}