/**
 * query_index.cpp  —  Query a grammar-compressed FM-index (v6)
 *
 * Changes from v5 → v6 (all motivated by 0-based row indexing):
 *
 *   FMIndex:
 *     sentinelRow field removed (was always 0).
 *
 *   loadIndex:
 *     Updated magic to FMIDX\x06; sentinelRow no longer read.
 *
 *   locateRow:
 *     arr = cur (not cur-1).  The LF step is now the standard 0-based formula:
 *       cur = C[sym] + rankSym(sym, cur)
 *     where rankSym counts occurrences in [0, cur) — no -1 adjustment.
 *     maxSteps fixed: only SA_SAMPLE_RATE steps are ever needed (each LF step
 *     decrements the compressed SA value by 1, so we reach a sample in at most
 *     SA_SAMPLE_RATE steps), replacing the broken "(int)compLen + ..." bound.
 *
 *   bwtAt:
 *     Assertion relaxed to sym <= numSymbols (allows 0 = sentinel char,
 *     though in practice we always return before reaching that row).
 *
 *   backwardSearch:
 *     Every LF interval step changed from
 *       C[s] + rankSym(s, lo-1),  C[s] + rankSym(s, hi-1)
 *     to
 *       C[s] + rankSym(s, lo),    C[s] + rankSym(s, hi)
 *     because lo/hi are now 0-based row numbers and rankSym(s, i) already
 *     counts occurrences in [0, i) exclusive — no adjustment needed.
 *
 * LF mapping (0-based):
 *   sym     = BWT[row]
 *   newRow  = C[sym] + rankSym(sym, row)   // rank in [0, row), exclusive
 *
 * SA sample lookup:
 *   Check bit row in saBv.  If set, rank = popcount(saBv[0..row]),
 *   result = saVals[rank-1] + offset.
 */

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

static int SA_SAMPLE_RATE = 32;

struct Rule { uint32_t left, right; };

struct FMIndex {
    int      symbolBits;
    uint64_t refLen, compLen;   // sentinelRow removed
    uint32_t numSymbols;
    vector<Rule>             rules;
    vector<uint64_t>         C;
    vector<vector<uint64_t>> wtLevels;
    size_t                   wtN;
    vector<uint64_t>         saBv;
    vector<uint64_t>         saVals;
    vector<vector<uint8_t>>  exps;
};

static uint32_t r32(ifstream& f){ uint32_t v; f.read((char*)&v,4); return v; }
static uint64_t r64(ifstream& f){ uint64_t v; f.read((char*)&v,8); return v; }

static FMIndex loadIndex(const string& path)
{
    ifstream f(path,ios::binary);
    if(!f){cerr<<"Cannot open "<<path<<"\n";exit(1);}
    char magic[7]; f.read(magic,7);
    // v6: sentinelRow removed from file
    if(strncmp(magic,"FMIDX\x06",6)!=0){cerr<<"Bad magic (expected v6)\n";exit(1);}
    FMIndex idx;
    idx.symbolBits=(int)r32(f);
    idx.refLen=r64(f); idx.compLen=r64(f);
    uint32_t nRules=r32(f);
    SA_SAMPLE_RATE=(int)r32(f);
    // sentinelRow NOT read (removed in v6)
    idx.numSymbols=r32(f);
    idx.exps.resize(idx.numSymbols+1);
    for(uint32_t i=1;i<=idx.numSymbols;i++){
        uint32_t elen=r32(f); idx.exps[i].resize(elen);
        f.read((char*)idx.exps[i].data(),elen);
    }
    idx.rules.resize(nRules);
    for(auto&r:idx.rules){r.left=r32(f);r.right=r32(f);}
    idx.C.resize(idx.numSymbols+1);
    for(auto&v:idx.C) v=r64(f);
    {
        uint64_t bvWords=r64(f); idx.saBv.resize(bvWords);
        for(auto&v:idx.saBv) v=r64(f);
        uint64_t nVals=r64(f); idx.saVals.resize(nVals);
        for(auto&v:idx.saVals) v=r64(f);
    }
    idx.wtLevels.resize(idx.symbolBits);
    for(int l=0;l<idx.symbolBits;l++){
        uint64_t nb=r64(f); idx.wtN=(size_t)nb;
        size_t words=(nb+63)/64; idx.wtLevels[l].resize(words);
        for(auto&w:idx.wtLevels[l]) w=r64(f);
    }
    return idx;
}

// ---------------------------------------------------------------------------
// Wavelet tree  (unchanged)
// ---------------------------------------------------------------------------

static inline uint64_t wt_rank1(const FMIndex& idx, int l, uint64_t i)
{
    if(i==0) return 0;
    const auto&B=idx.wtLevels[l];
    uint64_t cnt=0, full=i/64;
    for(uint64_t w=0;w<full;w++) cnt+=__builtin_popcountll(B[w]);
    uint64_t rem=i%64;
    if(rem) cnt+=__builtin_popcountll(B[full]&((1ULL<<rem)-1ULL));
    return cnt;
}

static uint64_t rankSym(const FMIndex& idx, uint32_t sym, uint64_t i)
{
    uint64_t lo=0,hi=i; size_t n=idx.wtN;
    for(int l=idx.symbolBits-1;l>=0;l--){
        int b=(sym>>l)&1;
        uint64_t ol=wt_rank1(idx,l,lo), oh=wt_rank1(idx,l,hi);
        if(b==0){lo=lo-ol;hi=hi-oh;}
        else{uint64_t tz=n-wt_rank1(idx,l,n);lo=tz+ol;hi=tz+oh;}
    }
    return hi-lo;
}

static uint32_t bwtAt(const FMIndex& idx, uint64_t r)
{
    uint32_t sym=0; uint64_t pos=r; size_t n=idx.wtN;
    for(int l=idx.symbolBits-1;l>=0;l--){
        int bit=(int)((idx.wtLevels[l][pos/64]>>(pos%64))&1ULL);
        sym|=(uint32_t)bit<<l;
        uint64_t rb=wt_rank1(idx,l,pos), ra=wt_rank1(idx,l,pos+1);
        if(bit==0) pos=pos-rb;
        else{uint64_t zeros=n-wt_rank1(idx,l,n);pos=zeros+ra-1;}
    }
    // sym==0 means sentinel char; locateRow never reaches that row in practice
    // because it is always sampled (compSA[i]=0, 0%rate==0).
    assert(sym<=idx.numSymbols);
    return sym;
}

// ---------------------------------------------------------------------------
// SA sample lookup  (unchanged)
// ---------------------------------------------------------------------------

static inline uint64_t bvPopcount(const vector<uint64_t>& bv, uint64_t row)
{
    uint64_t i=row+1, cnt=0, full=i/64;
    for(uint64_t w=0;w<full;w++) cnt+=__builtin_popcountll(bv[w]);
    uint64_t rem=i%64;
    if(rem) cnt+=__builtin_popcountll(bv[full]&((1ULL<<rem)-1ULL));
    return cnt;
}

// ---------------------------------------------------------------------------
// locateRow  (v6: 0-based, fixed maxSteps)
// ---------------------------------------------------------------------------
//
// Walk LF until we hit a sampled row.  In 0-based indexing:
//   - The BWT array is indexed directly by row (no -1 adjustment).
//   - LF(row) = C[BWT[row]] + rankSym(BWT[row], row)
//             where rankSym counts [0, row) exclusive.
//
// At most SA_SAMPLE_RATE LF steps are needed: each step decrements the
// underlying compressed SA value by 1 (mod compLen), so within SA_SAMPLE_RATE
// steps we reach a compressed position that is a multiple of SA_SAMPLE_RATE
// and therefore sampled.  The old bound "(int)compLen + rate + 2" was both
// wrong (int overflow for large refs) and unnecessarily large.

static uint64_t locateRow(const FMIndex& idx, uint64_t row)
{
    uint64_t offset=0, cur=row;
    for(int steps=0; steps<=SA_SAMPLE_RATE+2; steps++){
        if((idx.saBv[cur/64]>>(cur%64))&1ULL){
            uint64_t rank=bvPopcount(idx.saBv,cur);
            return (idx.saVals[rank-1]+offset)%idx.refLen;
        }
        uint32_t sym=bwtAt(idx,cur);
        offset+=idx.exps[sym].size();
        cur=idx.C[sym]+rankSym(idx,sym,cur);   // 0-based LF: no -1
    }
    cerr<<"ERROR: LF walk limit reached at row="<<row<<"\n";
    return UINT64_MAX;
}

// ---------------------------------------------------------------------------
// Query compression helpers  (unchanged)
// ---------------------------------------------------------------------------

static inline uint32_t baseToSym(char c)
{
    switch(c){case 'A':return 1;case 'C':return 2;case 'G':return 3;case 'T':return 4;}
    return 0;
}

static uint32_t terminalRemappedId(const FMIndex& idx, uint8_t term)
{
    for(uint32_t i=1;i<=idx.numSymbols;i++)
        if(idx.exps[i].size()==1&&idx.exps[i][0]==term) return i;
    return UINT32_MAX;
}

struct QSym{uint32_t sym;bool amb;};

static vector<QSym> compressQuery(const vector<uint32_t>& qRem,
                                   const FMIndex& idx,
                                   const vector<uint32_t>& ruleNT)
{
    int qlen=(int)qRem.size();
    vector<QSym> cur(qlen);
    for(int i=0;i<qlen;i++) cur[i]={qRem[i],false};
    for(size_t ri=0;ri<idx.rules.size();ri++){
        uint32_t c1=idx.rules[ri].left,c2=idx.rules[ri].right,c3=ruleNT[ri];
        if(c3==UINT32_MAX) continue;
        if(!cur.empty()&&cur.front().sym==c2) cur.front().amb=true;
        if(!cur.empty()&&cur.back().sym==c1)  cur.back().amb=true;
        vector<QSym> next; next.reserve(cur.size());
        size_t i=0,n=cur.size();
        while(i+1<n){
            if(cur[i].sym==c1 && cur[i+1].amb) cur[i].amb=true;
            if(cur[i+1].sym==c2 && cur[i].amb) cur[i+1].amb=true;
            i++;
        }
        i=0;
        while(i<n){
            if(i+1<n&&cur[i].sym==c1&&cur[i+1].sym==c2){
                if(cur[i].amb){cur[i+1].amb=true;next.push_back(cur[i]);i++;}
                else if(cur[i+1].amb){cur[i].amb=true;next.push_back(cur[i]);i++;}
                else{next.push_back({c3,false});i+=2;}
            } else{next.push_back(cur[i]);i++;}
        }
        cur=move(next);
    }
    return cur;
}

static vector<QSym> decompressAmbiguous(const vector<QSym>& cq,
                                         const FMIndex& idx,
                                         const vector<uint32_t>& termMap)
{
    vector<QSym> result; result.reserve(cq.size());
    for(const auto& qs : cq){
        if(!qs.amb){ result.push_back(qs); }
        else{
            const auto& exp=idx.exps[qs.sym];
            for(uint8_t t:exp){
                uint32_t tid=termMap[t];
                result.push_back({tid,true});
            }
        }
    }
    return result;
}

// ---------------------------------------------------------------------------
// Backward search  (v6: LF interval steps use lo/hi directly, not lo-1/hi-1)
// ---------------------------------------------------------------------------
//
// The interval [lo, hi) is now in 0-based row space.
// Standard FM-index interval LF step for symbol s:
//   new_lo = C[s] + rankSym(s, lo)
//   new_hi = C[s] + rankSym(s, hi)
// where rankSym(s, i) counts occurrences of s in BWT[0..i-1] (exclusive).
// No adjustment to lo or hi is needed.

static void backwardSearch(const FMIndex& idx, const vector<QSym>& cq,
                            vector<uint64_t>& positions)
{
    if(cq.empty()) return;

    int Q=(int)cq.size();

    struct Frame {
        uint64_t lo, hi;       // FM-index range [lo, hi), 0-based rows
        uint32_t f_sym;        // current F-column symbol
        uint64_t depth_in_f;   // terminals matched from the right end of f_sym's expansion
        int qi;                // next query symbol index (decreasing)
    };

    vector<Frame> stack;

    // Seed: match the rightmost query symbol cq[Q-1] against the F-column.
    const auto& lastQ=cq[Q-1];

    if(!lastQ.amb){
        uint32_t s=lastQ.sym;
        uint64_t flo=idx.C[s];
        uint64_t fhi=(s<idx.numSymbols)?idx.C[s+1]:(idx.compLen+1);
        if(flo<fhi){
            uint64_t slen=idx.exps[s].size();
            stack.push_back({flo,fhi,s,slen-1,Q-2});
        }
    } else {
        uint8_t c=idx.exps[lastQ.sym][0];
        for(uint32_t s=1;s<=idx.numSymbols;s++){
            uint64_t flo=idx.C[s];
            uint64_t fhi=(s<idx.numSymbols)?idx.C[s+1]:(idx.compLen+1);
            uint64_t slen=idx.exps[s].size();
            if(flo>=fhi||slen==0) continue;
            for(uint32_t d=0;d<slen;d++){
                if(idx.exps[s][slen-1-d]==c)
                    stack.push_back({flo,fhi,s,d,Q-2});
            }
        }
    }

    while(!stack.empty()){
        auto [lo,hi,f_sym,depth_in_f,qi]=stack.back();
        stack.pop_back();

        if(qi<0){
            // Matched the whole query; report one hit per row in [lo, hi).
            uint64_t offset=idx.exps[f_sym].size()-1-depth_in_f;
            for(uint64_t row=lo;row<hi;row++){
                uint64_t pos=locateRow(idx,row);
                if(pos!=UINT64_MAX) positions.push_back(pos+offset);
            }
            continue;
        }

        const auto& curQ=cq[qi];
        uint64_t flen=idx.exps[f_sym].size();
        uint64_t nd=depth_in_f+1;  // depth after consuming one more terminal

        if(!curQ.amb){
            // Non-ambiguous: must match the whole grammar symbol atomically.
            // If we haven't exhausted the current F-column symbol yet, prune.
            if(nd<flen) continue;
            // F-column symbol exhausted — LF step directly to curQ.sym.
            // v6: use lo/hi directly (0-based), no -1.
            uint32_t s=curQ.sym;
            uint64_t b_lo=idx.C[s]+rankSym(idx,s,lo);
            uint64_t b_hi=idx.C[s]+rankSym(idx,s,hi);
            if(b_lo<b_hi){
                uint64_t slen=idx.exps[s].size();
                stack.push_back({b_lo,b_hi,s,slen-1,qi-1});
            }
        } else {
            // Ambiguous (single terminal): terminal-by-terminal matching.
            uint8_t c=idx.exps[curQ.sym][0];

            if(nd<flen){
                // Still inside the current F-column symbol.
                if(idx.exps[f_sym][flen-1-nd]==c)
                    stack.push_back({lo,hi,f_sym,nd,qi-1});
            } else {
                // F-column symbol exhausted — LF step for each compatible symbol.
                // v6: use lo/hi directly (0-based), no -1.
                for(uint32_t b=1;b<=idx.numSymbols;b++){
                    uint64_t blen=idx.exps[b].size();
                    if(blen==0||idx.exps[b][blen-1]!=c) continue;
                    uint64_t b_lo=idx.C[b]+rankSym(idx,b,lo);
                    uint64_t b_hi=idx.C[b]+rankSym(idx,b,hi);
                    if(b_lo>=b_hi) continue;
                    stack.push_back({b_lo,b_hi,b,0,qi-1});
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// queryFM  (unchanged except uses updated backwardSearch)
// ---------------------------------------------------------------------------

static vector<uint64_t> queryFM(const FMIndex& idx, const string& queryStr,
                                  const vector<uint32_t>& ruleNT)
{
    if(queryStr.empty()) return {};
    vector<uint32_t> qRem; qRem.reserve(queryStr.size());
    for(char c:queryStr){
        uint32_t pre=baseToSym(c); if(!pre) return {};
        uint32_t rid=terminalRemappedId(idx,(uint8_t)pre); if(rid==UINT32_MAX) return {};
        qRem.push_back(rid);
    }

    vector<QSym> cq=compressQuery(qRem,idx,ruleNT);

    vector<uint32_t> termMap(256,UINT32_MAX);
    for(uint32_t i=1;i<=idx.numSymbols;i++){
        if(idx.exps[i].size()==1) termMap[idx.exps[i][0]]=i;
    }
    cq=decompressAmbiguous(cq,idx,termMap);

    vector<uint64_t> positions;
    backwardSearch(idx,cq,positions);
    sort(positions.begin(),positions.end());
    positions.erase(unique(positions.begin(),positions.end()),positions.end());
    while(!positions.empty()&&positions.back()+queryStr.size()>idx.refLen)
        positions.pop_back();
    return positions;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    if(argc!=3){cerr<<"Usage: query_index <index.idx> <queries.txt>\n";return 1;}
    cerr<<"Loading index...\n";
    FMIndex idx=loadIndex(argv[1]);
    cerr<<"Loaded: refLen="<<idx.refLen<<" compLen="<<idx.compLen
        <<" rules="<<idx.rules.size()<<" symbolBits="<<idx.symbolBits
        <<" numSymbols="<<idx.numSymbols<<"\n";

    // C[0]=0, C[1]=1 (one sentinel in F-column).
    assert(idx.C[0]==0 && idx.C[1]==1);

    vector<uint32_t> ruleNT(idx.rules.size(),UINT32_MAX);
    {
        unordered_map<string,uint32_t> expToId; expToId.reserve(idx.numSymbols);
        for(uint32_t i=1;i<=idx.numSymbols;i++)
            expToId[string(idx.exps[i].begin(),idx.exps[i].end())]=i;
        for(size_t ri=0;ri<idx.rules.size();ri++){
            uint32_t L=idx.rules[ri].left,R=idx.rules[ri].right;
            vector<uint8_t> combined;
            combined.insert(combined.end(),idx.exps[L].begin(),idx.exps[L].end());
            combined.insert(combined.end(),idx.exps[R].begin(),idx.exps[R].end());
            auto it=expToId.find(string(combined.begin(),combined.end()));
            if(it!=expToId.end()) ruleNT[ri]=it->second;
        }
    }

    ifstream qf(argv[2]);
    if(!qf){cerr<<"Cannot open "<<argv[2]<<"\n";return 1;}
    string line; int qnum=0;
    while(getline(qf,line)){
        while(!line.empty()&&(line.back()=='\r'||line.back()=='\n'||line.back()==' '))
            line.pop_back();
        if(line.empty()) continue;
        transform(line.begin(),line.end(),line.begin(),::toupper);
        qnum++;
        auto hits=queryFM(idx,line,ruleNT);
        cout<<"Query "<<qnum<<" ["<<line<<"]: ";
        if(hits.empty()) cout<<"NOT FOUND\n";
        else{
            for(size_t i=0;i<hits.size();i++){if(i)cout<<',';cout<<hits[i];}
            cout<<'\n';
        }
    }
    return 0;
}