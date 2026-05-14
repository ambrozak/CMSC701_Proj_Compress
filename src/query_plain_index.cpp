// query_plain_fm.cpp
//
// Query plain FM-index.
// Usage:
//   query_plain_fm <index.idx> <queries.txt>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

static constexpr int SIGMA = 3;

struct FMIndex{
    uint64_t n;
    vector<uint64_t> C;
    vector<uint64_t> saBv;    // bitvector: bit i set iff BWT row i is sampled
    vector<uint64_t> saVals;  // compacted text positions for sampled rows
    size_t wtN;
    vector<vector<uint64_t>> wtLevels;
};

static uint64_t r64(ifstream&f){uint64_t v;f.read((char*)&v,8);return v;}

static inline uint32_t baseToSym(char c){
    switch(c){
        case 'A': return 1;
        case 'C': return 2;
        case 'G': return 3;
        case 'T': return 4;
    }
    return 0;
}

static FMIndex loadIndex(const string& path){
    ifstream f(path,ios::binary);
    if(!f){ cerr<<"Cannot open "<<path<<"\n"; exit(1); }
    char magic[8]={0};
    f.read(magic,7);
    if(strcmp(magic,"PLAINFM")!=0){ cerr<<"Bad index format\n"; exit(1); }
    FMIndex idx;
    idx.n=r64(f);
    idx.C.resize(5);
    for(int i=0;i<=4;i++) idx.C[i]=r64(f);
    uint64_t nbv=r64(f);
    idx.saBv.resize(nbv);
    for(auto&v:idx.saBv) v=r64(f);
    uint64_t nvals=r64(f);
    idx.saVals.resize(nvals);
    for(auto&v:idx.saVals) v=r64(f);
    idx.wtLevels.resize(SIGMA);
    for(int l=0;l<SIGMA;l++){
        uint64_t nb=r64(f);
        idx.wtN=(size_t)nb;
        size_t words=(nb+63)/64;
        idx.wtLevels[l].resize(words);
        for(auto&w:idx.wtLevels[l]) w=r64(f);
    }
    return idx;
}

static uint64_t wt_rank1(const FMIndex&idx,int l,uint64_t i){
    if(!i) return 0;
    const auto& B=idx.wtLevels[l];
    uint64_t cnt=0;
    uint64_t full=i/64;
    for(uint64_t w=0;w<full;w++) cnt+=__builtin_popcountll(B[w]);
    uint64_t rem=i%64;
    if(rem) cnt+=__builtin_popcountll(B[full]<<(64-rem));
    return cnt;
}

static uint64_t rankSym(const FMIndex&idx,uint32_t sym,uint64_t i){
    uint64_t lo=0,hi=i;
    size_t n=idx.wtN;
    for(int l=SIGMA-1;l>=0;l--){
        int b=(sym>>l)&1;
        uint64_t ol=wt_rank1(idx,l,lo);
        uint64_t oh=wt_rank1(idx,l,hi);
        if(b==0){ lo=lo-ol; hi=hi-oh; }
        else{
            uint64_t tz=n-wt_rank1(idx,l,n);
            lo=tz+ol; hi=tz+oh;
        }
    }
    return hi-lo;
}

static uint32_t bwtAt(const FMIndex&idx,uint64_t r){
    uint32_t sym=0;
    uint64_t pos=r;
    size_t n=idx.wtN;
    for(int l=SIGMA-1;l>=0;l--){
        int bit=(int)((idx.wtLevels[l][pos/64]>>(pos%64))&1ULL);
        sym|=(uint32_t)bit<<l;
        uint64_t ob=wt_rank1(idx,l,pos);
        if(bit==0) pos=pos-ob;
        else pos=(n-wt_rank1(idx,l,n))+ob;
    }
    return sym;
}

static inline bool isRowSampled(const FMIndex&idx, uint64_t row){
    return (idx.saBv[row/64] >> (row%64)) & 1ULL;
}

// Number of sampled rows in [0..row] inclusive — used to index into saVals.
static uint64_t sampledRank(const FMIndex&idx, uint64_t row){
    uint64_t cnt=0;
    uint64_t full=row/64;
    for(uint64_t w=0;w<full;w++) cnt+=__builtin_popcountll(idx.saBv[w]);
    cnt+=__builtin_popcountll(idx.saBv[full]<<(63-row%64));
    return cnt;
}

// BWT[i] is 0-indexed: row 0 = sentinel row (BWT[0]=$=0).
// LF(i) = C[BWT[i]] + rankSym(BWT[i], i)   -- exclusive rank in BWT[0..i).
// saVals[sampledRank(row)-1] = SA[row] for sampled rows.
static uint64_t locateRow(const FMIndex&idx,uint64_t row){

    uint64_t offset=0;

    while(true){

        if(isRowSampled(idx,row))
            return (idx.saVals[sampledRank(idx,row)-1]+offset)%idx.n;

        uint32_t sym=bwtAt(idx,row);

        offset++;

        // LF: exclusive rank of sym in BWT[0..row)
        row=idx.C[sym]+rankSym(idx,sym,row);
    }
}

static vector<uint64_t> queryFM(const FMIndex&idx, const string&query){
    if(query.empty()) return {};

    // Rows are 0-indexed: row 0 = sentinel. Full range = [0, n+1).
    uint64_t lo=0;
    uint64_t hi=idx.n+1;

    for(int qi=(int)query.size()-1;qi>=0;qi--){

        uint32_t sym=baseToSym(query[qi]);

        if(!sym) return {};

        lo=idx.C[sym]+rankSym(idx,sym,lo);
        hi=idx.C[sym]+rankSym(idx,sym,hi);

        if(lo>=hi) return {};
    }

    vector<uint64_t> hits;

    for(uint64_t r=lo;r<hi;r++){
        if(r==0) continue; // skip sentinel row
        hits.push_back(locateRow(idx,r));
    }

    sort(hits.begin(),hits.end());

    return hits;
}

int main(int argc,char* argv[]){
    if(argc!=3){ cerr<<"Usage: query_plain_fm <index.idx> <queries.txt>\n"; return 1; }
    FMIndex idx=loadIndex(argv[1]);
    ifstream qf(argv[2]);
    if(!qf){ cerr<<"Cannot open "<<argv[2]<<"\n"; return 1; }
    string line;
    int qnum=0;
    while(getline(qf,line)){
        while(!line.empty() && (line.back()=='\r'||line.back()=='\n'||line.back()==' ')) line.pop_back();
        if(line.empty()) continue;
        std::transform(line.begin(),line.end(),line.begin(),::toupper);
        qnum++;
        auto hits=queryFM(idx,line);
        cout<<"Query "<<qnum<<" ["<<line<<"]: ";
        if(hits.empty()){ cout<<"NOT FOUND\n"; }
        else{
            for(size_t i=0;i<hits.size();i++){ if(i) cout<<','; cout<<hits[i]; }
            cout<<"\n";
        }
    }
}