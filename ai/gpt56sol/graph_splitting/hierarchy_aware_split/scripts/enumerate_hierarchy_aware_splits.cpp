#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include "sampler_data.hpp"
#include "hierarchy_meta.hpp"

struct Cand {
    int id;
    int64_t sum_w;
    int32_t max_w;
    int32_t sum_d;
    int32_t min_d;
};


struct DPNode {
    bool valid=false;
    int64_t sum_w=0;
    int32_t max_w=0;
    int64_t sum_d=0;
    std::vector<Cand> sel;
};
static bool dp_better(const DPNode&a,const DPNode&b){
    if(!b.valid)return true;
    if(a.sum_w!=b.sum_w)return a.sum_w<b.sum_w;
    if(a.max_w!=b.max_w)return a.max_w<b.max_w;
    if(a.sum_d!=b.sum_d)return a.sum_d>b.sum_d;
    std::vector<int>ai,bi;for(auto&x:a.sel)ai.push_back(x.id);for(auto&x:b.sel)bi.push_back(x.id);std::sort(ai.begin(),ai.end());std::sort(bi.begin(),bi.end());return ai<bi;
}
static int enc5(const std::array<int,5>&q){return ((((q[0]*5+q[1])*5+q[2])*5+q[3])*5+q[4]);}
static std::array<int,5> dec5(int z){std::array<int,5>q{};for(int i=4;i>=0;--i){q[i]=z%5;z/=5;}return q;}
static bool select_coverall_stratified(const std::vector<Cand>&c,const std::vector<int>&strata,std::vector<Cand>&sel){
    std::vector<int> branches;for(int b=0;b<NBRANCH;++b){for(auto&x:c)if(BRANCH_ID[x.id]==b){branches.push_back(b);break;}}
    if(branches.size()>20)return false;
    const int NST=3125;std::vector<DPNode>dp(NST),ndp(NST);dp[0].valid=true;
    auto betterCand=[](const Cand&a,const Cand&b){if(a.sum_w!=b.sum_w)return a.sum_w<b.sum_w;if(a.max_w!=b.max_w)return a.max_w<b.max_w;if(a.sum_d!=b.sum_d)return a.sum_d>b.sum_d;return a.id<b.id;};
    for(size_t ib=0;ib<branches.size();++ib){int br=branches[ib];std::array<std::vector<Cand>,5>g;for(auto&x:c)if(BRANCH_ID[x.id]==br)g[strata[x.id]].push_back(x);for(auto&v:g)std::sort(v.begin(),v.end(),betterCand);
        struct Opt{std::array<int,5>q;int64_t sw;int32_t mw;int64_t sd;std::vector<Cand>v;};std::vector<Opt>opts;
        for(int q0=0;q0<=std::min<int>(4,g[0].size());++q0)for(int q1=0;q1<=std::min<int>(4,g[1].size());++q1)for(int q2=0;q2<=std::min<int>(4,g[2].size());++q2)for(int q3=0;q3<=std::min<int>(4,g[3].size());++q3)for(int q4=0;q4<=std::min<int>(4,g[4].size());++q4){std::array<int,5>q{q0,q1,q2,q3,q4};int n=q0+q1+q2+q3+q4;if(n<1||n>20)continue;Opt o{q,0,0,0,{}};for(int z=0;z<5;++z)for(int k=0;k<q[z];++k){auto x=g[z][k];o.v.push_back(x);o.sw+=x.sum_w;o.mw=std::max(o.mw,x.max_w);o.sd+=x.sum_d;}opts.push_back(std::move(o));}
        for(auto&x:ndp)x=DPNode{};int remaining=(int)branches.size()-(int)ib-1;
        for(int st=0;st<NST;++st)if(dp[st].valid){auto base=dec5(st);int baseN=std::accumulate(base.begin(),base.end(),0);for(auto&o:opts){std::array<int,5>q=base;bool ok=true;for(int z=0;z<5;++z){q[z]+=o.q[z];if(q[z]>4)ok=false;}if(!ok)continue;int n=baseN+(int)o.v.size();if(n+remaining>20)continue;DPNode cand; cand.valid=true;cand.sum_w=dp[st].sum_w+o.sw;cand.max_w=std::max(dp[st].max_w,o.mw);cand.sum_d=dp[st].sum_d+o.sd;cand.sel=dp[st].sel;cand.sel.insert(cand.sel.end(),o.v.begin(),o.v.end());int ns=enc5(q);if(dp_better(cand,ndp[ns]))ndp[ns]=std::move(cand);}}
        dp.swap(ndp);
    }
    int target=enc5({4,4,4,4,4});if(!dp[target].valid)return false;sel=dp[target].sel;return sel.size()==20;
}

struct Eval {
    uint32_t total_rows=0, unseen_rows=0, unique_genes=0, unseen_unique=0;
    uint32_t total_pos=0, fn=0;
    double row_overlap=0, unique_overlap=0, pos_coverage=0, f1=0;
    double mean_wup=0, max_wup=0, mean_best_wup=0, mean_dist=0, mean_nearest_dist=0;
    int min_dist=999, related_pairs=0, train_branches=0, held_branches=0;
    uint64_t train_nodes=0, train_edges=0, held_nodes=0, held_edges=0;
    int train_internal=0, held_internal=0;
};

static std::string join_ids(const std::vector<int>& x) {
    std::ostringstream o;
    for(size_t i=0;i<x.size();++i){if(i)o<<'|';o<<x[i];}
    return o.str();
}
static std::string join_names(const std::vector<int>& x) {
    std::ostringstream o;
    for(size_t i=0;i<x.size();++i){if(i)o<<'|';o<<NETWORK_NAME[x[i]];}
    return o.str();
}

static Eval evaluate(const std::vector<int>& T, const std::array<int,4>& H) {
    Eval r;
    uint64_t tu[NWORDS]{};
    for(int t:T){
        r.train_nodes+=NODE_COUNT[t];r.train_edges+=EDGE_COUNT[t];r.train_internal+=IS_INTERNAL[t];
        for(int w=0;w<NWORDS;++w)tu[w]|=MASKS[t][w];
    }
    uint64_t hu[NWORDS]{};
    for(int h:H){r.held_nodes+=NODE_COUNT[h];r.held_edges+=EDGE_COUNT[h];r.held_internal+=IS_INTERNAL[h];r.total_rows+=NODE_COUNT[h];r.total_pos+=TISSUE_POS_TOTAL[h];}
    for(int h:H){
        for(int w=0;w<NWORDS;++w){
            uint64_t m=MASKS[h][w]; hu[w]|=m; uint64_t miss=m&~tu[w];
            r.unseen_rows+=std::popcount(miss);
            while(miss){unsigned b=std::countr_zero(miss);int g=w*64+(int)b;if(g<NGENE)r.fn+=POS_WEIGHT[g];miss&=miss-1;}
        }
    }
    for(int w=0;w<NWORDS;++w){r.unique_genes+=std::popcount(hu[w]);r.unseen_unique+=std::popcount(hu[w]&~tu[w]);}
    uint32_t tp=r.total_pos-r.fn;
    r.row_overlap=(double)(r.total_rows-r.unseen_rows)/r.total_rows;
    r.unique_overlap=(double)(r.unique_genes-r.unseen_unique)/r.unique_genes;
    r.pos_coverage=(double)tp/r.total_pos;
    r.f1=(double)(2ull*tp)/(2ull*tp+r.fn);
    int64_t sw=0,sd=0;int32_t mw=0;std::array<int32_t,4> best{};std::array<int32_t,4> near{};near.fill(999);
    std::set<int> tb,hb;
    for(int h:H)hb.insert(BRANCH_ID[h]);
    for(int t:T){tb.insert(BRANCH_ID[t]);for(int j=0;j<4;++j){int h=H[j];sw+=WUP_INT[t][h];sd+=HIER_DIST[t][h];mw=std::max(mw,WUP_INT[t][h]);best[j]=std::max(best[j],WUP_INT[t][h]);near[j]=std::min(near[j],(int32_t)HIER_DIST[t][h]);r.min_dist=std::min(r.min_dist,(int)HIER_DIST[t][h]);r.related_pairs+=HIER_RELATED[t][h];}}
    r.mean_wup=(double)sw/(80.0*WUP_SCALE);r.max_wup=(double)mw/WUP_SCALE;
    r.mean_best_wup=(double)std::accumulate(best.begin(),best.end(),int64_t(0))/(4.0*WUP_SCALE);
    r.mean_dist=(double)sd/80.0;r.mean_nearest_dist=(double)std::accumulate(near.begin(),near.end(),int64_t(0))/4.0;
    r.train_branches=tb.size();r.held_branches=hb.size();return r;
}

static std::vector<int> make_pool(bool leaf, bool large) {
    std::vector<int> p;
    if(large){for(int i=0;i<NLARGE;++i)if(!leaf||IS_LEAF[LARGE[i]])p.push_back(LARGE[i]);}
    else {for(int i=0;i<NELIG;++i)if(!leaf||IS_LEAF[ELIG[i]])p.push_back(ELIG[i]);}
    return p;
}

static std::vector<int> make_strata(const std::vector<int>& E) {
    std::vector<int> o=E;std::sort(o.begin(),o.end(),[](int a,int b){if(NODE_COUNT[a]!=NODE_COUNT[b])return NODE_COUNT[a]<NODE_COUNT[b];return a<b;});
    std::vector<int>s(NNET,-1);int n=o.size(),q=n/5,r=n%5,pos=0;
    for(int b=0;b<5;++b){int z=q+(b<r);for(int k=0;k<z;++k)s[o[pos++]]=b;}
    return s;
}

static bool select_training(const std::array<int,4>& H,const std::vector<int>& E,const std::string& method,const std::vector<int>& strata,std::vector<int>& T,std::tuple<int64_t,int64_t,int64_t>& obj){
    bool distinct=method.rfind("branch_distinct",0)==0;
    if(distinct){std::set<int>b;for(int h:H)b.insert(BRANCH_ID[h]);if(b.size()<4)return false;}
    std::set<int>hs(H.begin(),H.end());std::set<int>hbr;for(int h:H)hbr.insert(BRANCH_ID[h]);
    std::vector<Cand> c;
    for(int x:E){if(hs.count(x))continue;if(method.rfind("branch_distinct",0)==0&&hbr.count(BRANCH_ID[x]))continue;if(method.rfind("ancestor_blocked",0)==0){bool bad=false;for(int h:H)bad|=HIER_RELATED[x][h];if(bad)continue;}
        Cand z{x,0,0,0,999};for(int h:H){z.sum_w+=WUP_INT[x][h];z.max_w=std::max(z.max_w,WUP_INT[x][h]);z.sum_d+=HIER_DIST[x][h];z.min_d=std::min(z.min_d,(int32_t)HIER_DIST[x][h]);}c.push_back(z);
    }
    if(c.size()<20)return false;
    auto cmp=[](const Cand&a,const Cand&b){return std::tie(a.sum_w,a.max_w,b.sum_d,a.id)<std::tie(b.sum_w,b.max_w,a.sum_d,b.id);};
    // Explicit comparator to avoid confusing reverse tie expression.
    auto better=[](const Cand&a,const Cand&b){if(a.sum_w!=b.sum_w)return a.sum_w<b.sum_w;if(a.max_w!=b.max_w)return a.max_w<b.max_w;if(a.sum_d!=b.sum_d)return a.sum_d>b.sum_d;return a.id<b.id;};
    auto betterDistance=[](const Cand&a,const Cand&b){if(a.sum_d!=b.sum_d)return a.sum_d>b.sum_d;if(a.min_d!=b.min_d)return a.min_d>b.min_d;if(a.sum_w!=b.sum_w)return a.sum_w<b.sum_w;if(a.max_w!=b.max_w)return a.max_w<b.max_w;return a.id<b.id;};
    auto select_minimax=[&](const std::vector<Cand>&pool,int k)->std::vector<Cand>{
        std::vector<int32_t> v;v.reserve(pool.size());for(const auto&x:pool)v.push_back(x.max_w);
        std::nth_element(v.begin(),v.begin()+k-1,v.end());int32_t z=v[k-1];
        std::vector<Cand>a;for(const auto&x:pool)if(x.max_w<=z)a.push_back(x);
        std::sort(a.begin(),a.end(),better);a.resize(k);return a;
    };
    std::vector<Cand> sel;
    if(method=="conditional_minimax_wup" || method=="branch_distinct_minimax_wup"){
        sel=select_minimax(c,20);
    } else if(method=="branch_distinct_coverall_mean_wup"){
        std::set<int> used;for(int b=0;b<NBRANCH;++b){std::vector<Cand>a;for(auto&x:c)if(BRANCH_ID[x.id]==b)a.push_back(x);if(a.empty())continue;std::sort(a.begin(),a.end(),better);sel.push_back(a[0]);used.insert(a[0].id);}if(sel.size()>20)return false;
        std::vector<Cand>a;for(auto&x:c)if(!used.count(x.id))a.push_back(x);std::sort(a.begin(),a.end(),better);sel.insert(sel.end(),a.begin(),a.begin()+(20-sel.size()));
    } else if(method=="branch_distinct_coverall_node_stratified_mean_wup"){
        if(!select_coverall_stratified(c,strata,sel))return false;
    } else if(method=="branch_distinct_node_stratified_mean_wup"){
        for(int b=0;b<5;++b){std::vector<Cand>a;for(auto&x:c)if(strata[x.id]==b)a.push_back(x);if(a.size()<4)return false;std::sort(a.begin(),a.end(),better);sel.insert(sel.end(),a.begin(),a.begin()+4);}
    } else if(method=="branch_distinct_node_stratified_mean_distance"){
        for(int b=0;b<5;++b){std::vector<Cand>a;for(auto&x:c)if(strata[x.id]==b)a.push_back(x);if(a.size()<4)return false;std::sort(a.begin(),a.end(),betterDistance);sel.insert(sel.end(),a.begin(),a.begin()+4);}
    } else if(method=="branch_distinct_node_stratified_maximin_distance"){
        std::array<std::vector<Cand>,5> pools;int32_t global_z=999;
        for(int b=0;b<5;++b){for(auto&x:c)if(strata[x.id]==b)pools[b].push_back(x);if(pools[b].size()<4)return false;std::vector<int32_t>v;for(auto&x:pools[b])v.push_back(x.min_d);std::nth_element(v.begin(),v.begin()+3,v.end(),std::greater<int32_t>());global_z=std::min(global_z,v[3]);}
        for(int b=0;b<5;++b){std::vector<Cand>a;for(auto&x:pools[b])if(x.min_d>=global_z)a.push_back(x);if(a.size()<4)return false;std::sort(a.begin(),a.end(),betterDistance);sel.insert(sel.end(),a.begin(),a.begin()+4);}
    } else if(method=="branch_distinct_node_stratified_minimax_wup"){
        std::array<std::vector<Cand>,5> pools;int32_t global_z=0;
        for(int b=0;b<5;++b){for(auto&x:c)if(strata[x.id]==b)pools[b].push_back(x);if(pools[b].size()<4)return false;std::vector<int32_t>v;for(auto&x:pools[b])v.push_back(x.max_w);std::nth_element(v.begin(),v.begin()+3,v.end());global_z=std::max(global_z,v[3]);}
        for(int b=0;b<5;++b){std::vector<Cand>a;for(auto&x:pools[b])if(x.max_w<=global_z)a.push_back(x);if(a.size()<4)return false;std::sort(a.begin(),a.end(),better);sel.insert(sel.end(),a.begin(),a.begin()+4);}
    } else {std::sort(c.begin(),c.end(),better);sel.assign(c.begin(),c.begin()+20);}
    T.clear();int64_t sumw=0,sumd=0,maxw=0;for(auto&x:sel){T.push_back(x.id);sumw+=x.sum_w;sumd+=x.sum_d;maxw=std::max<int64_t>(maxw,x.max_w);}std::sort(T.begin(),T.end());
    int64_t global_mind=999;for(auto&x:sel)global_mind=std::min<int64_t>(global_mind,x.min_d);
    if(method=="conditional_minimax_wup" || method=="branch_distinct_minimax_wup" || method=="branch_distinct_node_stratified_minimax_wup")obj={maxw,sumw,-sumd};
    else if(method=="branch_distinct_node_stratified_maximin_distance")obj={-global_mind,-sumd,sumw};
    else if(method=="branch_distinct_node_stratified_mean_distance")obj={-sumd,-global_mind,sumw};
    else obj={sumw,maxw,-sumd};return true;
}

int main(int argc,char**argv){
    if(argc!=5){std::cerr<<"usage: enum <all144|leaf107> <method> <out.tsv> <summary.txt>\n";return 2;}
    std::string universe=argv[1],method=argv[2];bool leaf=universe=="leaf107";if(!leaf&&universe!="all144")return 2;
    auto E=make_pool(leaf,false),L=make_pool(leaf,true),strata=make_strata(E);
    std::ofstream out(argv[3]);out<<"universe\tmethod\theldout_indices\theldout_tissues\ttraining_indices\ttraining_tissues\tobjective_primary\tobjective_secondary\tobjective_tertiary\tmean_wup\tmax_wup\tmean_heldout_best_wup\tmean_distance\tmin_distance\tmean_heldout_nearest_distance\tancestor_descendant_pairs\ttrain_branch_count\theldout_branch_count\ttraining_node_total\ttraining_edge_total\theldout_node_total\theldout_edge_total\ttraining_internal_layers\theldout_internal_layers\ttotal_rows\tunseen_rows\trow_overlap\tunique_genes\tunseen_unique\tunique_overlap\ttotal_positive_labels\tfalse_negative_labels\tpositive_coverage\tlookup_micro_f1\n";
    bool have=false;std::tuple<int64_t,int64_t,int64_t> best;uint64_t best_count=0,feasible=0,total=0;std::vector<std::string> best_lines;
    for(size_t a=0;a<L.size();++a)for(size_t b=a+1;b<L.size();++b)for(size_t c=b+1;c<L.size();++c)for(size_t d=c+1;d<L.size();++d){++total;std::array<int,4>H{L[a],L[b],L[c],L[d]};std::vector<int>T;std::tuple<int64_t,int64_t,int64_t>obj;if(!select_training(H,E,method,strata,T,obj))continue;++feasible;Eval r=evaluate(T,H);std::vector<int>Hv(H.begin(),H.end());
        std::ostringstream line;line<<universe<<'\t'<<method<<'\t'<<join_ids(Hv)<<'\t'<<join_names(Hv)<<'\t'<<join_ids(T)<<'\t'<<join_names(T)<<'\t'<<std::get<0>(obj)<<'\t'<<std::get<1>(obj)<<'\t'<<std::get<2>(obj)<<'\t'<<std::setprecision(12)<<r.mean_wup<<'\t'<<r.max_wup<<'\t'<<r.mean_best_wup<<'\t'<<r.mean_dist<<'\t'<<r.min_dist<<'\t'<<r.mean_nearest_dist<<'\t'<<r.related_pairs<<'\t'<<r.train_branches<<'\t'<<r.held_branches<<'\t'<<r.train_nodes<<'\t'<<r.train_edges<<'\t'<<r.held_nodes<<'\t'<<r.held_edges<<'\t'<<r.train_internal<<'\t'<<r.held_internal<<'\t'<<r.total_rows<<'\t'<<r.unseen_rows<<'\t'<<r.row_overlap<<'\t'<<r.unique_genes<<'\t'<<r.unseen_unique<<'\t'<<r.unique_overlap<<'\t'<<r.total_pos<<'\t'<<r.fn<<'\t'<<r.pos_coverage<<'\t'<<r.f1;
        std::string s=line.str();out<<s<<'\n';if(!have||obj<best){have=true;best=obj;best_count=1;best_lines={s};}else if(obj==best){++best_count;best_lines.push_back(s);} }
    out.close();std::ofstream q(argv[4]);q<<"universe="<<universe<<"\nmethod="<<method<<"\neligible_train="<<E.size()<<"\neligible_heldout="<<L.size()<<"\ntotal_H="<<total<<"\nfeasible_H="<<feasible<<"\nbest_primary="<<std::get<0>(best)<<"\nbest_secondary="<<std::get<1>(best)<<"\nbest_tertiary="<<std::get<2>(best)<<"\nbest_count="<<best_count<<"\n";for(auto&s:best_lines)q<<"BEST\t"<<s<<"\n";q.close();
    std::cerr<<universe<<' '<<method<<" total="<<total<<" feasible="<<feasible<<" best_count="<<best_count<<"\n";return 0;
}
