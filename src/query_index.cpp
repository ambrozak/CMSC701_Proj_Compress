/**
 * query_index.cpp  —  Query a grammar-compressed FM-index (v5)
 *
 * LF mapping:
 *   arr         = logicalRow - 1
 *   sym         = BWT[arr]
 *   newLogRow   = C[sym] + rankSym(sym, arr)   // rank in [0, arr), exclusive
 *
 * SA sample lookup (bitvector popcount-rank):
 *   Check bit arr in saBv. If set, rank = popcount(saBv[0..arr]),
 *   result = saVals[rank-1] + offset.
 *
 * After query compression, ambiguous symbols are decompressed back to
 * terminals so the final query is a mix of large non-ambiguous grammar
 * symbols and individual terminals (in ambiguous zones).
 *
 * Non-ambiguous optimizations:
 *   1) Non-ambiguous query symbols match F-column symbols by identity,
 *      no terminal-by-terminal walk needed.
 *   2) LF-step for non-ambiguous symbols goes directly to C[sym]+rank,
 *      no scanning all F-column symbols for compatible expansions.
 *   3) If walking terminal-by-terminal inside an F-column symbol and we
 *      reach a non-ambiguous query symbol before exhausting the expansion,
 *      we know it's a mismatch and prune the branch.
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
    uint64_t refLen, compLen, sentinelRow;
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
    if(strncmp(magic,"FMIDX\x05",6)!=0){cerr<<"Bad magic (expected v5)\n";exit(1);}
    FMIndex idx;
    idx.symbolBits=(int)r32(f);
    idx.refLen=r64(f); idx.compLen=r64(f);
    uint32_t nRules=r32(f);
    SA_SAMPLE_RATE=(int)r32(f);
    idx.sentinelRow=r64(f);
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

// Wavelet tree
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
    assert(sym>=1&&sym<=idx.numSymbols);
    return sym;
}

// SA sample lookup
static inline uint64_t bvPopcount(const vector<uint64_t>& bv, uint64_t arr)
{
    uint64_t i=arr+1, cnt=0, full=i/64;
    for(uint64_t w=0;w<full;w++) cnt+=__builtin_popcountll(bv[w]);
    uint64_t rem=i%64;
    if(rem) cnt+=__builtin_popcountll(bv[full]&((1ULL<<rem)-1ULL));
    return cnt;
}

static uint64_t locateRow(const FMIndex& idx, uint64_t logRow)
{
    uint64_t offset=0, cur=logRow;
    int maxSteps=(int)idx.compLen+SA_SAMPLE_RATE+2;
    for(int steps=0;steps<=maxSteps;steps++){
        uint64_t arr=cur-1;
        if((idx.saBv[arr/64]>>(arr%64))&1ULL){
            uint64_t rank=bvPopcount(idx.saBv,arr);
            uint64_t base=idx.saVals[rank-1];
            return (base+offset)%idx.refLen;
        }
        uint32_t sym=bwtAt(idx,arr);
        offset+=idx.exps[sym].size();
        cur=idx.C[sym]+rankSym(idx,sym,arr);
    }
    cerr<<"ERROR: LF walk limit reached at logRow="<<logRow<<"\n";
    return UINT64_MAX;
}

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
            if(cur[i].sym==c1 && cur[i+1].amb) {
                cur[i].amb = true;
            }
            if(cur[i+1].sym==c2 && cur[i].amb) {
                cur[i+1].amb = true;
            }
            i++;
        }
        i = 0;
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

static string qsymToString(const vector<QSym>& v, const FMIndex& idx)
{
    string out;

    auto appendSym = [&](uint32_t sym) {
        const auto& exp = idx.exps[sym];

        // Single-terminal symbols: print DNA character
        if(exp.size() == 1){
            switch(exp[0]){
                case 1: out += 'A'; break;
                case 2: out += 'C'; break;
                case 3: out += 'G'; break;
                case 4: out += 'T'; break;
                default:
                    out += '?';
                    break;
            }
        } else {
            out += "sym";
            out += to_string(sym);
        }
    };

    for(size_t i = 0; i < v.size(); i++){
        if(i) out += ' ';

        appendSym(v[i].sym);

        if(v[i].amb)
            out += "*";
    }

    return out;
}

static string fullyDecompressQSyms(const vector<QSym>& v, const FMIndex& idx)
{
    string out;

    for(const auto& qs : v){
        const auto& exp = idx.exps[qs.sym];

        for(uint8_t t : exp){
            switch(t){
                case 1: out += 'A'; break;
                case 2: out += 'C'; break;
                case 3: out += 'G'; break;
                case 4: out += 'T'; break;
                default:
                    out += '?';
                    break;
            }
        }
    }

    return out;
}
// Decompress ambiguous symbols back to their terminal expansions. 
// Each ambiguous symbol is replaced by one QSym per terminal byte, 
// using the remapped terminal IDs, all marked amb=true. 
static vector<QSym> decompressAmbiguous(const vector<QSym>& cq, 
    const FMIndex& idx, const vector<uint32_t>& termMap) { 
    vector<QSym> result; 
    result.reserve(cq.size()); 
    for(const auto& qs : cq){ 
        if(!qs.amb){ result.push_back(qs); } 
        else { // Expand this symbol back to individual terminals 
            const auto& exp = idx.exps[qs.sym]; 
            for(uint8_t t : exp){ 
                uint32_t tid = termMap[t]; 
                result.push_back({tid, true}); 
            } 
        } 
    } 
    return result; 
}

static void backwardSearch(const FMIndex& idx, const vector<QSym>& cq,
                            vector<uint64_t>& positions)
{
    if(cq.empty()) return;

    // Build the terminal sequence of the entire query (for terminal-by-terminal
    // matching in ambiguous regions). Also record, for each terminal position,
    // whether it belongs to a non-ambiguous compressed symbol, and if so which
    // query-symbol index it came from. This lets us detect optimization (3):
    // hitting a non-ambiguous boundary inside an F-column expansion.
    //
    // But actually, a simpler framing: the query is now a sequence of QSyms
    // where non-ambiguous ones are potentially multi-terminal grammar symbols,
    // and ambiguous ones are always single terminals. We process them in reverse
    // order (backward search).

    int Q = (int)cq.size();

    // For a non-ambiguous symbol, its expansion length tells us how many
    // terminals it covers. For an ambiguous symbol, it's always a single
    // terminal (expansion size 1).

    // We need a frame that tracks:
    //   - The current range [lo, hi) in the compressed FM-index
    //   - qi: the current query symbol index (processing right to left)
    //   - For ambiguous symbols: which F-column symbol we're inside, and
    //     how deep we are in its expansion (terminal-by-terminal matching)
    //   - For non-ambiguous symbols: we match atomically
    //
    // State machine per query symbol (right to left):
    //   If cq[qi].amb == false (non-ambiguous):
    //     The only valid F-column symbol is cq[qi].sym itself.
    //     Do a direct LF-step: lo' = C[s] + rank(s, lo-1), hi' = C[s] + rank(s, hi-1)
    //     Advance qi to qi-1.
    //
    //   If cq[qi].amb == true (ambiguous, single terminal):
    //     We need terminal-by-terminal matching against F-column symbol expansions.
    //     This is similar to the original code but with optimization (3):
    //     if while walking inside an F-column expansion we'd next need to match
    //     a non-ambiguous query symbol, that's a mismatch (prune).

    // We'll track the "offset into the reference" for position calculation.
    // When we finish matching the entire query and report a hit, the position
    // in the reference depends on the F-column symbol we're inside and our
    // depth into it.

    struct Frame {
        uint64_t lo, hi;       // FM-index range (logical rows, 1-based)
        uint32_t f_sym;        // current F-column symbol
        uint64_t depth_in_f;   // how deep into f_sym's expansion (from the end)
        int qi;                // next query symbol index to match (decreasing)
    };

    vector<Frame> stack;

    // Initialize: match the rightmost query symbol cq[Q-1].
    const auto& lastQ = cq[Q-1];

    if(!lastQ.amb){
        // Non-ambiguous: the F-column symbol must be exactly lastQ.sym.
        // We start fully "inside" this symbol at depth 0 (matched its last terminal),
        // but since it's non-ambiguous and atomic, we consume it entirely.
        uint32_t s = lastQ.sym;
        uint64_t flo = idx.C[s];
        uint64_t fhi = (s < idx.numSymbols) ? idx.C[s+1] : (idx.compLen+1);
        if(flo < fhi){
            uint64_t slen = idx.exps[s].size();
            // We've matched the entire symbol. depth_in_f = slen-1 means
            // we're at the beginning of the expansion (fully consumed).
            stack.push_back({flo, fhi, s, slen-1, Q-2});
        }
    } else {
        // Ambiguous (single terminal): find all F-column symbols whose
        // expansion ends with this terminal.
        uint8_t c = idx.exps[lastQ.sym][0]; // single terminal's raw value
        for(uint32_t s = 1; s <= idx.numSymbols; s++){
            uint64_t flo = idx.C[s];
            uint64_t fhi = (s < idx.numSymbols) ? idx.C[s+1] : (idx.compLen+1);
            uint64_t slen = idx.exps[s].size();
            if(flo >= fhi || slen == 0) continue;
            // Check if the last byte of this symbol's expansion matches
            for (uint32_t d = 0; d < slen; d++){
                if(idx.exps[s][slen-1-d] == c){
                    stack.push_back({flo, fhi, s, d, Q-2});
                }
            }
        }
    }

    while(!stack.empty()){
        auto [lo, hi, f_sym, depth_in_f, qi] = stack.back();
        stack.pop_back();

        // If we've matched the entire query, report hits.
        if(qi < 0){
            uint64_t offset = idx.exps[f_sym].size() - 1 - depth_in_f;
            for(uint64_t logRow = lo; logRow < hi; logRow++){
                uint64_t pos = locateRow(idx, logRow);
                if(pos != UINT64_MAX) positions.push_back(pos + offset);
            }
            continue;
        }

        const auto& curQ = cq[qi];
        uint64_t flen = idx.exps[f_sym].size();
        uint64_t nd = depth_in_f + 1; // next depth if we continue inside f_sym

        if(!curQ.amb){
            // Non-ambiguous query symbol: optimization (1) and (2).
            // This symbol can only appear as a complete unit in the BWT.
            // So if we're partway through an F-column symbol (nd < flen),
            // that's optimization (3): mismatch, prune.
            if(nd < flen){
                // We haven't exhausted the current F-column symbol, but the
                // next query symbol is non-ambiguous and can't be a sub-part.
                // Prune this branch.
                continue;
            }
            // nd == flen: we've exhausted the F-column symbol. Do an LF-step
            // directly to the non-ambiguous query symbol (optimization 2).
            uint32_t s = curQ.sym;
            uint64_t b_lo = idx.C[s] + rankSym(idx, s, lo-1);
            uint64_t b_hi = idx.C[s] + rankSym(idx, s, hi-1);
            if(b_lo < b_hi){
                uint64_t slen = idx.exps[s].size();
                // Fully consumed this non-ambiguous symbol.
                stack.push_back({b_lo, b_hi, s, slen-1, qi-1});
            }
        } else {
            // Ambiguous (single terminal): terminal-by-terminal matching.
            uint8_t c = idx.exps[curQ.sym][0]; // the raw terminal value

            if(nd < flen){
                // Still inside the current F-column symbol.
                // Check if the next byte matches.
                if(idx.exps[f_sym][flen-1-nd] == c){
                    stack.push_back({lo, hi, f_sym, nd, qi-1});
                }
            } else {
                // Exhausted the F-column symbol; need an LF-step.
                // Find all F-column symbols whose expansion ends with c.
                for(uint32_t b = 1; b <= idx.numSymbols; b++){
                    uint64_t blen = idx.exps[b].size();
                    if(blen == 0 || idx.exps[b][blen-1] != c) continue;
                    uint64_t b_lo = idx.C[b] + rankSym(idx, b, lo-1);
                    uint64_t b_hi = idx.C[b] + rankSym(idx, b, hi-1);
                    if(b_lo >= b_hi) continue;
                    stack.push_back({b_lo, b_hi, b, 0, qi-1});
                }
            }
        }
    }
}

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

    // Step 1: compress the query (full grammar compression with ambiguity tracking)
    vector<QSym> cq = compressQuery(qRem, idx, ruleNT);

    // Step 2: decompress ambiguous symbols back to terminals.
    // Build a terminal-value -> remapped-symbol-id map for the decompression.
    // Terminal values are the raw byte values (1=A, 2=C, 3=G, 4=T).
    // We need to find the symbol ID for each single-byte expansion.
    vector<uint32_t> termMap(256, UINT32_MAX);
    for(uint32_t i = 1; i <= idx.numSymbols; i++){
        if(idx.exps[i].size() == 1){
            termMap[idx.exps[i][0]] = i;
        }
    }
    cq = decompressAmbiguous(cq, idx, termMap);

    vector<uint64_t> positions;
    backwardSearch(idx, cq, positions);
    sort(positions.begin(),positions.end());
    positions.erase(unique(positions.begin(),positions.end()),positions.end());
    while(!positions.empty()&&positions.back()+queryStr.size()>idx.refLen)
        positions.pop_back();
    return positions;
}

int main(int argc, char* argv[])
{
    if(argc!=3){cerr<<"Usage: query_index <index.idx> <queries.txt>\n";return 1;}
    cerr<<"Loading index...\n";
    FMIndex idx=loadIndex(argv[1]);
    cerr<<"Loaded: refLen="<<idx.refLen<<" compLen="<<idx.compLen
        <<" rules="<<idx.rules.size()<<" symbolBits="<<idx.symbolBits
        <<" numSymbols="<<idx.numSymbols<<"\n";

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
    // Keep rules around since compressQuery still needs them
    // (they're used per-query in compressQuery)

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