/**
 * build_index.cpp  —  RePair compression + grammar-compressed FM-index
 *
 * Usage:  build_index <in.fa> <symbolBits> <out.idx>
 *
 * Binary layout (version 5):
 *   "FMIDX\x05"  (7 bytes)
 *   symbolBits   (uint32)
 *   refLen       (uint64)
 *   compLen      (uint64)
 *   numRules     (uint32)
 *   saSampleRate (uint32)
 *   sentinelRow  (uint64, always 0)
 *   numSymbols   (uint32, IDs 1..numSymbols)
 *   For each symbol 1..numSymbols: expLen (uint32) + expBytes (uint8 x expLen)
 *   Rules: numRules x (left, right) as uint32
 *   C table: (numSymbols+1) x uint64, indices 0..numSymbols
 *   SA sample bitvector: bvWords (uint64) then bvWords x uint64
 *   SA sample values:    nVals (uint64) then nVals x uint64
 *   Wavelet tree: symbolBits levels, each: n (uint64) + ceil(n/64) x uint64
 *   (Compressed sequence is NOT stored.)
 */

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

static constexpr int SA_SAMPLE_RATE = 32;

static string loadFirstSeq(const string& fn)
{
    ifstream f(fn);
    if(!f){ cerr << "Cannot open " << fn << "\n"; exit(1); }
    string line, seq; bool in = false;
    while(getline(f, line)){
        if(line.empty()) continue;
        if(line[0] == '>'){if(in)break;in=true;continue;}
        if(in) seq += line;
    }
    return seq;
}

static inline uint32_t baseToSym(char c)
{
    switch(c){case 'A':return 1;case 'C':return 2;case 'G':return 3;case 'T':return 4;}
    return 0;
}

static void w32(ofstream& f, uint32_t v){ f.write((char*)&v, 4); }
static void w64(ofstream& f, uint64_t v){ f.write((char*)&v, 8); }

struct Rule { uint32_t left, right; };

static void runRepair(vector<uint32_t>& seq, vector<Rule>& rules, int symbolBits)
{
    const uint32_t maxSym = (1u << symbolBits);
    uint32_t nextSym = 5;
    while(nextSym < maxSym){
        if(seq.size() < 2) break;
        unordered_map<uint64_t,uint32_t> freq; freq.reserve(seq.size());
        for(size_t i = 0; i+1 < seq.size();){
            uint64_t key = ((uint64_t)seq[i]<<32)|seq[i+1]; freq[key]++;
            if(seq[i]==seq[i+1]) i+=2; else i++;
        }
        uint64_t bk=0; uint32_t bc=0;
        for(auto&[k,v]:freq) if(v>bc){bc=v;bk=k;}
        if(bc<2) break;
        uint32_t L=(uint32_t)(bk>>32), R=(uint32_t)(bk&0xFFFFFFFFu), ns=nextSym++;
        rules.push_back({L,R});
        size_t w=0,r=0;
        while(r<seq.size()){
            if(r+1<seq.size()&&seq[r]==L&&seq[r+1]==R){seq[w++]=ns;r+=2;}
            else seq[w++]=seq[r++];
        }
        seq.resize(w);
        if(rules.size()%100==0)
            cerr<<"  Rule "<<rules.size()<<" seqlen="<<seq.size()<<"\n";
    }
    cerr<<"RePair: "<<rules.size()<<" rules, compLen="<<seq.size()<<"\n";
}

static void buildExpansions(const vector<Rule>& rules,
                             vector<vector<uint8_t>>& exps,
                             vector<uint64_t>& lens)
{
    size_t total=5+rules.size(); exps.resize(total); lens.resize(total,0);
    for(int i=1;i<=4;i++){exps[i]={(uint8_t)i};lens[i]=1;}
    for(size_t i=0;i<rules.size();i++){
        auto&el=exps[rules[i].left],&er=exps[rules[i].right],&e=exps[5+i];
        e.insert(e.end(),el.begin(),el.end());
        e.insert(e.end(),er.begin(),er.end());
        lens[5+i]=lens[rules[i].left]+lens[rules[i].right];
    }
}

static vector<uint32_t> buildRemap(const vector<vector<uint8_t>>& exps, size_t total)
{
    vector<uint32_t> order; order.reserve(total-1);
    for(uint32_t i=1;i<(uint32_t)total;i++) order.push_back(i);
    sort(order.begin(),order.end(),[&](uint32_t a,uint32_t b){ return exps[a]<exps[b]; });
    vector<uint32_t> remap(total,UINT32_MAX);
    for(size_t i=0;i<order.size();i++) remap[order[i]]=(uint32_t)(i+1);
    return remap;
}

// SA-IS
static void classifySL(const vector<uint32_t>& t, vector<bool>& isS)
{
    int n=(int)t.size()-1; isS.assign(n+1,false); isS[n]=true;
    if(n==0)return; isS[n-1]=false;
    for(int i=n-2;i>=0;i--) isS[i]=(t[i]<t[i+1])||(t[i]==t[i+1]&&isS[i+1]);
}
static inline bool isLMS(const vector<bool>& isS,int i){return i>0&&isS[i]&&!isS[i-1];}
static void getBuckets(const vector<uint32_t>& t,int n,int alpha,vector<int>& bkt,bool end)
{
    bkt.assign(alpha,0);
    for(int i=0;i<=n;i++) bkt[(int)t[i]]++;
    int sum=0;
    for(int c=0;c<alpha;c++){sum+=bkt[c];bkt[c]=end?sum:sum-bkt[c];}
}
static void induceSortL(const vector<uint32_t>& t,int n,int alpha,const vector<bool>& isS,vector<int>& sa)
{
    vector<int> bkt; getBuckets(t,n,alpha,bkt,false);
    for(int i=0;i<=n;i++){if(sa[i]<0)continue;int j=sa[i]-1;if(j>=0&&!isS[j])sa[bkt[(int)t[j]]++]=j;}
}
static void induceSortS(const vector<uint32_t>& t,int n,int alpha,const vector<bool>& isS,vector<int>& sa)
{
    vector<int> bkt; getBuckets(t,n,alpha,bkt,true);
    for(int i=n;i>=0;i--){if(sa[i]<0)continue;int j=sa[i]-1;if(j>=0&&isS[j])sa[--bkt[(int)t[j]]]=j;}
}
static bool lmsEqual(const vector<uint32_t>& t,int n,const vector<bool>& isS,int i,int j)
{
    if(i==n||j==n)return false;
    for(int d=0;;d++){
        if(i+d>n||j+d>n)return false;
        int ti=(i+d==n)?0:(int)t[i+d], tj=(j+d==n)?0:(int)t[j+d];
        if(ti!=tj||isS[i+d]!=isS[j+d])return false;
        if(d>0&&isLMS(isS,i+d)&&isLMS(isS,j+d))return true;
    }
}
static vector<int> sais(const vector<uint32_t>& t,int n,int alpha)
{
    vector<bool> isS; classifySL(t,isS);
    vector<int> sa(n+1,-1);
    {vector<int> bkt; getBuckets(t,n,alpha,bkt,true);
     for(int i=1;i<n;i++) if(isLMS(isS,i)) sa[--bkt[(int)t[i]]]=i; sa[0]=n;}
    induceSortL(t,n,alpha,isS,sa); induceSortS(t,n,alpha,isS,sa);
    vector<int> lms;
    for(int i=1;i<=n;i++) if(sa[i]>=0&&isLMS(isS,sa[i])) lms.push_back(sa[i]);
    vector<int> rank(n+1,-1); int curRank=0,prev=-1;
    for(int cur:lms){if(prev==-1||!lmsEqual(t,n,isS,prev,cur))curRank++;rank[cur]=curRank-1;prev=cur;}
    vector<int> lmsOrder;
    for(int i=1;i<n;i++) if(isLMS(isS,i)) lmsOrder.push_back(i);
    int n1=(int)lmsOrder.size(); vector<uint32_t> t1(n1);
    for(int i=0;i<n1;i++) t1[i]=(uint32_t)rank[lmsOrder[i]];
    vector<int> sa1;
    if(curRank<n1){
        vector<uint32_t> t1s(n1+1);
        for(int i=0;i<n1;i++) t1s[i]=t1[i]+1; t1s[n1]=0;
        sa1=sais(t1s,n1,curRank+1);
    } else {
        sa1.assign(n1+1,-1); sa1[0]=n1;
        for(int i=0;i<n1;i++) sa1[t1[i]+1]=i;
    }
    fill(sa.begin(),sa.end(),-1);
    {vector<int> bkt; getBuckets(t,n,alpha,bkt,true);
     for(int i=n1;i>=1;i--){int pos=lmsOrder[sa1[i]];sa[--bkt[(int)t[pos]]]=pos;} sa[0]=n;}
    induceSortL(t,n,alpha,isS,sa); induceSortS(t,n,alpha,isS,sa);
    return sa;
}
static vector<int> buildCompSA(const vector<uint32_t>& comp, uint32_t numSymbols)
{
    int n=(int)comp.size(); vector<uint32_t> t(n+1);
    for(int i=0;i<n;i++) t[i]=comp[i]; t[n]=0;
    return sais(t,n,(int)numSymbols+1);
}

static vector<uint32_t> buildBWT(const vector<int>& compSA,
                                  const vector<uint32_t>& comp,
                                  uint64_t& sentinelRow)
{
    int cn=(int)comp.size(); vector<uint32_t> bwt; bwt.reserve(cn);
    sentinelRow=UINT64_MAX;
    for(int i=0;i<=cn;i++){
        int pos=compSA[i];
        if(pos==cn){assert(i==0);sentinelRow=(uint64_t)bwt.size();continue;}
        bwt.push_back(pos==0 ? comp[cn-1] : comp[pos-1]);
    }
    assert((int)bwt.size()==cn); assert(sentinelRow==0);
    return bwt;
}

static vector<uint64_t> buildCTable(const vector<uint32_t>& bwt, uint32_t numSymbols)
{
    vector<uint64_t> C(numSymbols+2,0);
    for(uint32_t s:bwt){assert(s>=1&&s<=numSymbols);C[s+1]++;}
    for(uint32_t i=1;i<=numSymbols+1;i++) C[i]+=C[i-1];
    for(uint32_t i=0;i<=numSymbols;i++) C[i]+=1;
    C.resize(numSymbols+1);
    return C;
}

struct SASamples { vector<uint64_t> bv, vals; };

static SASamples buildSASamples(const vector<int>& compSA,
                                 const vector<uint64_t>& compStart,
                                 uint64_t compLen)
{
    SASamples s;
    s.bv.assign((compLen+63)/64, 0ULL);
    uint64_t bwtRow=0;
    for(int i=0;i<=(int)compLen;i++){
        int pos=compSA[i];
        if(pos==(int)compLen){continue;}
        if((uint64_t)pos%(uint64_t)SA_SAMPLE_RATE==0){
            s.bv[bwtRow/64]|=1ULL<<(bwtRow%64);
            s.vals.push_back(compStart[pos]);
        }
        bwtRow++;
    }
    assert(bwtRow==compLen);
    return s;
}

struct WaveletTree{int symbolBits;size_t n;vector<vector<uint64_t>>levels;};

static WaveletTree buildWT(const vector<uint32_t>& bwt, int symbolBits)
{
    WaveletTree wt; wt.symbolBits=symbolBits; wt.n=bwt.size();
    size_t n=bwt.size(), words=(n+63)/64; wt.levels.resize(symbolBits);
    vector<uint32_t> perm(n); for(size_t i=0;i<n;i++) perm[i]=(uint32_t)i;
    for(int l=symbolBits-1;l>=0;l--){
        wt.levels[l].assign(words,0ULL);
        for(size_t i=0;i<n;i++) if((bwt[perm[i]]>>l)&1) wt.levels[l][i/64]|=1ULL<<(i%64);
        vector<uint32_t> L,R; L.reserve(n); R.reserve(n);
        for(size_t i=0;i<n;i++){
            if(((bwt[perm[i]]>>l)&1)==0) L.push_back(perm[i]); else R.push_back(perm[i]);
        }
        perm.clear(); perm.insert(perm.end(),L.begin(),L.end()); perm.insert(perm.end(),R.begin(),R.end());
    }
    return wt;
}

static inline uint32_t bitsNeeded(uint32_t maxVal)
{
    uint32_t bits=0; while((1u<<bits)<=maxVal) bits++; return bits;
}

int main(int argc, char* argv[])
{
    if(argc!=4){cerr<<"Usage: build_index <in.fa> <symbolBits> <out.idx>\n";return 1;}
    string inFile=argv[1]; int symbolBits=atoi(argv[2]); string outFile=argv[3];
    if(symbolBits<2||symbolBits>30){cerr<<"symbolBits must be in [2,30]\n";return 1;}

    cerr<<"Loading...\n";
    string raw=loadFirstSeq(inFile);
    transform(raw.begin(),raw.end(),raw.begin(),::toupper);
    vector<uint32_t> seq; seq.reserve(raw.size());
    for(char c:raw){uint32_t s=baseToSym(c);if(s)seq.push_back(s);}
    raw.clear(); raw.shrink_to_fit();
    uint64_t refLen=seq.size();
    cerr<<"refLen="<<refLen<<"\n";

    cerr<<"RePair...\n";
    vector<Rule> rules; runRepair(seq,rules,symbolBits);
    uint64_t compLen=seq.size();

    cerr<<"Building expansions...\n";
    vector<vector<uint8_t>> exps; vector<uint64_t> lens;
    buildExpansions(rules,exps,lens);
    size_t totalSyms=5+rules.size();

    cerr<<"Remapping...\n";
    vector<uint32_t> remap=buildRemap(exps,totalSyms);
    for(auto&s:seq) s=remap[s];
    for(auto&r:rules){r.left=remap[r.left];r.right=remap[r.right];}
    uint32_t numSymbols=(uint32_t)(totalSyms-1);

    uint32_t neededBits=bitsNeeded(numSymbols);
    if((1u<<symbolBits)<=numSymbols){
        cerr<<"ERROR: symbolBits="<<symbolBits<<" insufficient for "<<numSymbols
            <<" symbols (need at least "<<neededBits<<" bits)\n"; return 1;
    }

    vector<vector<uint8_t>> exps2(numSymbols+1); vector<uint64_t> lens2(numSymbols+1,0);
    for(uint32_t old=1;old<(uint32_t)totalSyms;old++){
        exps2[remap[old]]=exps[old]; lens2[remap[old]]=lens[old];
    }
    exps=move(exps2); lens=move(lens2);

    for(uint32_t i=1;i<numSymbols;i++){
        if(!(exps[i]<exps[i+1])){
            cerr<<"ERROR: expansion ordering violated at IDs "<<i<<" and "<<i+1<<"\n"; return 1;
        }
    }

    cerr<<"Building compressed SA (SA-IS)...\n";
    vector<int> compSA=buildCompSA(seq,numSymbols);

    cerr<<"Building BWT...\n";
    uint64_t sentinelRow=UINT64_MAX;
    vector<uint32_t> bwt=buildBWT(compSA,seq,sentinelRow);
    assert(sentinelRow==0);

    vector<uint64_t> compStart(compLen+1); compStart[0]=0;
    for(size_t i=0;i<compLen;i++) compStart[i+1]=compStart[i]+lens[seq[i]];
    assert(compStart[compLen]==refLen);

    cerr<<"Building SA samples...\n";
    SASamples saSamples=buildSASamples(compSA,compStart,compLen);

    cerr<<"Building C table...\n";
    vector<uint64_t> Ctable=buildCTable(bwt,numSymbols);
    assert(Ctable[1]==1);
    {uint64_t cnt=0;for(uint32_t s:bwt)if(s==numSymbols)cnt++;assert(Ctable[numSymbols]+cnt==compLen+1);}

    cerr<<"Building wavelet tree...\n";
    WaveletTree wt=buildWT(bwt,symbolBits);

    cerr<<"Writing "<<outFile<<"...\n";
    ofstream out(outFile,ios::binary);
    if(!out){cerr<<"Cannot open "<<outFile<<"\n";return 1;}

    out.write("FMIDX\x05\x00",7);
    w32(out,(uint32_t)symbolBits); w64(out,refLen); w64(out,compLen);
    w32(out,(uint32_t)rules.size()); w32(out,(uint32_t)SA_SAMPLE_RATE); w64(out,sentinelRow);
    w32(out,numSymbols);
    for(uint32_t i=1;i<=numSymbols;i++){w32(out,(uint32_t)exps[i].size());out.write((char*)exps[i].data(),exps[i].size());}
    for(auto&r:rules){w32(out,r.left);w32(out,r.right);}
    for(uint32_t i=0;i<=numSymbols;i++) w64(out,Ctable[i]);
    w64(out,(uint64_t)saSamples.bv.size());
    for(uint64_t v:saSamples.bv) w64(out,v);
    w64(out,(uint64_t)saSamples.vals.size());
    for(uint64_t v:saSamples.vals) w64(out,v);
    for(int l=0;l<symbolBits;l++){
        uint64_t nb=(uint64_t)wt.n; w64(out,nb);
        size_t words=(nb+63)/64;
        for(size_t wd=0;wd<words;wd++) w64(out,wt.levels[l][wd]);
    }
    out.close();

    cerr<<"\n=== Summary ===\n";
    cerr<<"refLen="<<refLen<<" compLen="<<compLen<<" rules="<<rules.size()<<"\n";
    cerr<<"numSymbols="<<numSymbols<<" (IDs 1.."<<numSymbols<<", symbolBits="<<symbolBits<<")\n";
    cerr<<"sentinelRow="<<sentinelRow<<"\n";
    cerr<<"saSamples="<<saSamples.vals.size()<<"\n";
    cerr<<"Done.\n";
    return 0;
}