#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <tuple>
#include <vector>
#include "sampler_data.hpp"
#include "hierarchy_meta.hpp"

struct Xoshiro256ss {
    uint64_t s[4];
    static uint64_t splitmix64(uint64_t &x){uint64_t z=(x+=UINT64_C(0x9e3779b97f4a7c15));z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);return z^(z>>31);}    
    explicit Xoshiro256ss(uint64_t seed){for(int i=0;i<4;i++)s[i]=splitmix64(seed);}    
    static uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
    uint64_t next(){uint64_t result=rotl(s[1]*5,7)*9,t=s[1]<<17;s[2]^=s[0];s[3]^=s[1];s[1]^=s[2];s[0]^=s[3];s[2]^=t;s[3]=rotl(s[3],45);return result;}
    uint64_t bounded(uint64_t n){uint64_t threshold=(uint64_t)(-n)%n;for(;;){uint64_t r=next();if(r>=threshold)return r%n;}}
};

template<class T> static void partial_select(T *a,int n,int k,Xoshiro256ss &rng){for(int i=0;i<k;i++){int j=i+(int)rng.bounded((uint64_t)(n-i));std::swap(a[i],a[j]);}}

struct SubEval {float mean_wup=0,row_overlap=0,unique_overlap=0,f1=0;};

#pragma pack(push,1)
struct Record {
    uint16_t h0,h1,h2,h3;
    uint16_t a0,a1,b0,b1;
    float mean_wup,max_wup,mean_best_wup,mean_dist,min_dist,mean_nearest_dist;
    float row_overlap,unique_overlap,pos_coverage,f1;
    float roleA_mean_wup,roleB_mean_wup,role_worst_mean_wup,role_absdiff_wup;
    float role_mean_f1,role_min_f1,role_max_f1;
    float role_mean_row_overlap,role_min_row_overlap,role_max_row_overlap;
    uint32_t train_nodes,train_edges,held_nodes,held_edges;
    uint8_t train_branches,held_branches,train_internal,held_internal,related_pairs;
};
#pragma pack(pop)

static std::vector<int> make_pool(bool leaf,bool large){std::vector<int>p;if(large){for(int i=0;i<NLARGE;i++)if(!leaf||IS_LEAF[LARGE[i]])p.push_back(LARGE[i]);}else{for(int i=0;i<NELIG;i++)if(!leaf||IS_LEAF[ELIG[i]])p.push_back(ELIG[i]);}return p;}
static std::vector<int> make_strata(const std::vector<int>&E){std::vector<int>o=E;std::sort(o.begin(),o.end(),[](int a,int b){if(NODE_COUNT[a]!=NODE_COUNT[b])return NODE_COUNT[a]<NODE_COUNT[b];return a<b;});std::vector<int>s(NNET,-1);int n=o.size(),q=n/5,r=n%5,pos=0;for(int b=0;b<5;b++){int z=q+(b<r);for(int k=0;k<z;k++)s[o[pos++]]=b;}return s;}

static SubEval eval_subset(const std::array<int,20>&T,const std::array<int,2>&H){
    uint64_t tu[NWORDS]{};for(int t:T)for(int w=0;w<NWORDS;w++)tu[w]|=MASKS[t][w];
    uint64_t hu[NWORDS]{};uint32_t total_rows=0,unseen_rows=0,total_pos=0,fn=0,unique=0,unseen_unique=0;
    int64_t sw=0;
    for(int h:H){total_rows+=NODE_COUNT[h];total_pos+=TISSUE_POS_TOTAL[h];for(int t:T)sw+=WUP_INT[t][h];for(int w=0;w<NWORDS;w++){uint64_t m=MASKS[h][w];hu[w]|=m;uint64_t miss=m&~tu[w];unseen_rows+=std::popcount(miss);while(miss){unsigned b=std::countr_zero(miss);int g=w*64+(int)b;if(g<NGENE)fn+=POS_WEIGHT[g];miss&=miss-1;}}}
    for(int w=0;w<NWORDS;w++){unique+=std::popcount(hu[w]);unseen_unique+=std::popcount(hu[w]&~tu[w]);}
    uint32_t tp=total_pos-fn;SubEval r;r.mean_wup=(float)((double)sw/(40.0*WUP_SCALE));r.row_overlap=(float)((double)(total_rows-unseen_rows)/total_rows);r.unique_overlap=(float)((double)(unique-unseen_unique)/unique);r.f1=(float)((double)(2ull*tp)/(2ull*tp+fn));return r;
}

static Record evaluate(const std::array<int,20>&T,const std::array<int,4>&H){
    Record r{};r.h0=H[0];r.h1=H[1];r.h2=H[2];r.h3=H[3];
    uint64_t tu[NWORDS]{};std::set<int>tb,hb;for(int t:T){r.train_nodes+=NODE_COUNT[t];r.train_edges+=EDGE_COUNT[t];r.train_internal+=IS_INTERNAL[t];tb.insert(BRANCH_ID[t]);for(int w=0;w<NWORDS;w++)tu[w]|=MASKS[t][w];}
    uint64_t hu[NWORDS]{};uint32_t total_rows=0,unseen_rows=0,total_pos=0,fn=0,unique=0,unseen_unique=0;int64_t sw=0,sd=0;int32_t mw=0,md=999;std::array<int32_t,4>best{};std::array<int32_t,4>near{};near.fill(999);
    for(int j=0;j<4;j++){int h=H[j];r.held_nodes+=NODE_COUNT[h];r.held_edges+=EDGE_COUNT[h];r.held_internal+=IS_INTERNAL[h];hb.insert(BRANCH_ID[h]);total_rows+=NODE_COUNT[h];total_pos+=TISSUE_POS_TOTAL[h];for(int t:T){sw+=WUP_INT[t][h];sd+=HIER_DIST[t][h];mw=std::max(mw,WUP_INT[t][h]);md=std::min(md,(int32_t)HIER_DIST[t][h]);best[j]=std::max(best[j],WUP_INT[t][h]);near[j]=std::min(near[j],(int32_t)HIER_DIST[t][h]);r.related_pairs+=HIER_RELATED[t][h];}for(int w=0;w<NWORDS;w++){uint64_t m=MASKS[h][w];hu[w]|=m;uint64_t miss=m&~tu[w];unseen_rows+=std::popcount(miss);while(miss){unsigned b=std::countr_zero(miss);int g=w*64+(int)b;if(g<NGENE)fn+=POS_WEIGHT[g];miss&=miss-1;}}}
    for(int w=0;w<NWORDS;w++){unique+=std::popcount(hu[w]);unseen_unique+=std::popcount(hu[w]&~tu[w]);}
    uint32_t tp=total_pos-fn;r.mean_wup=(float)((double)sw/(80.0*WUP_SCALE));r.max_wup=(float)((double)mw/WUP_SCALE);r.mean_best_wup=(float)((double)(best[0]+best[1]+best[2]+best[3])/(4.0*WUP_SCALE));r.mean_dist=(float)((double)sd/80.0);r.min_dist=(float)md;r.mean_nearest_dist=(float)((near[0]+near[1]+near[2]+near[3])/4.0);r.row_overlap=(float)((double)(total_rows-unseen_rows)/total_rows);r.unique_overlap=(float)((double)(unique-unseen_unique)/unique);r.pos_coverage=(float)((double)tp/total_pos);r.f1=(float)((double)(2ull*tp)/(2ull*tp+fn));r.train_branches=(uint8_t)tb.size();r.held_branches=(uint8_t)hb.size();
    // Pick one of the three unordered 2+2 partitions using hierarchy only:
    // minimize worse train-role mean WUP; then imbalance; then cross-role mean WUP; then IDs.
    const int P[3][4]={{0,1,2,3},{0,2,1,3},{0,3,1,2}};bool have=false;std::tuple<int64_t,int64_t,int64_t,std::array<int,4>>score,bestscore;SubEval bestA,bestB;std::array<int,4>bp{};
    for(auto&p:P){std::array<int,2>A{H[p[0]],H[p[1]]},B{H[p[2]],H[p[3]]};SubEval ea=eval_subset(T,A),eb=eval_subset(T,B);int64_t sa=0,sb=0;for(int t:T){for(int h:A)sa+=WUP_INT[t][h];for(int h:B)sb+=WUP_INT[t][h];}int64_t worst=std::max(sa,sb),diff=std::llabs(sa-sb),cross=0;for(int a:A)for(int b:B)cross+=WUP_INT[a][b];std::array<int,4>ids{std::min(A[0],A[1]),std::max(A[0],A[1]),std::min(B[0],B[1]),std::max(B[0],B[1])};if(std::tie(ids[2],ids[3])<std::tie(ids[0],ids[1]))std::swap_ranges(ids.begin(),ids.begin()+2,ids.begin()+2);score={worst,diff,cross,ids};if(!have||score<bestscore){have=true;bestscore=score;bestA=ea;bestB=eb;bp={A[0],A[1],B[0],B[1]};}}
    r.a0=bp[0];r.a1=bp[1];r.b0=bp[2];r.b1=bp[3];r.roleA_mean_wup=bestA.mean_wup;r.roleB_mean_wup=bestB.mean_wup;r.role_worst_mean_wup=std::max(bestA.mean_wup,bestB.mean_wup);r.role_absdiff_wup=std::abs(bestA.mean_wup-bestB.mean_wup);r.role_mean_f1=(bestA.f1+bestB.f1)/2;r.role_min_f1=std::min(bestA.f1,bestB.f1);r.role_max_f1=std::max(bestA.f1,bestB.f1);r.role_mean_row_overlap=(bestA.row_overlap+bestB.row_overlap)/2;r.role_min_row_overlap=std::min(bestA.row_overlap,bestB.row_overlap);r.role_max_row_overlap=std::max(bestA.row_overlap,bestB.row_overlap);return r;
}

static std::vector<std::array<int,4>> feasible_h(const std::vector<int>&L,const std::vector<int>&E,const std::vector<int>&strata,bool matched){std::vector<std::array<int,4>>v;for(size_t a=0;a<L.size();a++)for(size_t b=a+1;b<L.size();b++)for(size_t c=b+1;c<L.size();c++)for(size_t d=c+1;d<L.size();d++){std::array<int,4>H{L[a],L[b],L[c],L[d]};if(matched){std::set<int>br;for(int h:H)br.insert(BRANCH_ID[h]);if(br.size()!=4)continue;int n[5]{};for(int x:E)if(!br.count(BRANCH_ID[x]))n[strata[x]]++;bool ok=true;for(int z:n)if(z<4)ok=false;if(!ok)continue;}v.push_back(H);}return v;}

int main(int argc,char**argv){if(argc!=7){std::cerr<<"usage: sample <all144|leaf107> <uniform|matched_stratified> <N> <seed> <out.bin> <audit.tsv>\n";return 2;}std::string universe=argv[1],mode=argv[2];bool leaf=universe=="leaf107";bool matched=mode=="matched_stratified";if((!leaf&&universe!="all144")||(!matched&&mode!="uniform"))return 2;uint64_t N=std::stoull(argv[3]),seed=std::stoull(argv[4]);auto E=make_pool(leaf,false),L=make_pool(leaf,true),strata=make_strata(E);auto Hpool=feasible_h(L,E,strata,matched);if(Hpool.empty())return 3;std::ofstream out(argv[5],std::ios::binary);const char magic[8]={'G','S','H','N','U','L','L','1'};out.write(magic,8);out.write((char*)&N,8);out.write((char*)&seed,8);uint32_t rs=sizeof(Record),nh=Hpool.size();out.write((char*)&rs,4);out.write((char*)&nh,4);Xoshiro256ss rng(seed);std::ofstream audit(argv[6]);audit<<"sample\theldout\troleA\troleB\ttraining\tmean_wup\tmax_wup\trow_overlap\tunique_overlap\tf1\trole_min_f1\n";
    for(uint64_t i=0;i<N;i++){auto H=Hpool[rng.bounded(Hpool.size())];std::array<int,20>T{};std::set<int>hs(H.begin(),H.end());if(!matched){std::vector<int>pool;for(int x:E)if(!hs.count(x))pool.push_back(x);partial_select(pool.data(),pool.size(),20,rng);std::copy(pool.begin(),pool.begin()+20,T.begin());}else{std::set<int>hbr;for(int h:H)hbr.insert(BRANCH_ID[h]);int pos=0;for(int b=0;b<5;b++){std::vector<int>pool;for(int x:E)if(!hbr.count(BRANCH_ID[x])&&strata[x]==b)pool.push_back(x);partial_select(pool.data(),pool.size(),4,rng);for(int k=0;k<4;k++)T[pos++]=pool[k];}}Record r=evaluate(T,H);out.write((char*)&r,sizeof(r));if(i<1000){audit<<i<<'\t'<<r.h0<<'|'<<r.h1<<'|'<<r.h2<<'|'<<r.h3<<'\t'<<r.a0<<'|'<<r.a1<<'\t'<<r.b0<<'|'<<r.b1<<'\t';for(int k=0;k<20;k++){if(k)audit<<'|';audit<<T[k];}audit<<'\t'<<r.mean_wup<<'\t'<<r.max_wup<<'\t'<<r.row_overlap<<'\t'<<r.unique_overlap<<'\t'<<r.f1<<'\t'<<r.role_min_f1<<'\n';}}
    std::cerr<<"universe="<<universe<<" mode="<<mode<<" N="<<N<<" seed="<<seed<<" record_size="<<sizeof(Record)<<" Hpool="<<Hpool.size()<<"\n";return 0;}
