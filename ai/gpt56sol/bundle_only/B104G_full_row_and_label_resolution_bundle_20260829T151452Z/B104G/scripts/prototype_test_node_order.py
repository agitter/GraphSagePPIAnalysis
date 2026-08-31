import csv, json, os
from pathlib import Path

# Python2 dict sim (64-bit)
def py2_hash_str(s,bits=64):
    b=s.encode('ascii'); mask=(1<<bits)-1
    if not b:return 0
    x=(b[0]<<7)&mask
    for c in b:x=((1000003*x)^c)&mask
    x^=len(b); x&=mask
    if x >= 1<<(bits-1):x-=1<<bits
    if x==-1:x=-2
    return x

def py2_hash_int(x,bits=64):
    # PyInt hash in CPython2 is the C long value, except -1 -> -2
    h=int(x)
    if h==-1:h=-2
    return h

class Py2Dict:
    def __init__(self, hashfn, bits=64):
        self.hashfn=hashfn;self.bits=bits;self.mask=7;self.table=[None]*8;self.used=0;self.fill=0
    def lookup(self,key,h):
        mask=self.mask;i=h&mask;perturb=h&((1<<self.bits)-1)
        while True:
            e=self.table[i]
            if e is None:return i
            if e[0]==key:return i
            i=(i*5+1+perturb)&mask;perturb>>=5
    def resize(self,minused):
        new=8
        while new<=minused:new<<=1
        old=[e for e in self.table if e is not None]
        self.mask=new-1;self.table=[None]*new;self.used=self.fill=0
        for key,h in old:self._insert(key,h)
    def _insert(self,key,h):
        i=self.lookup(key,h)
        if self.table[i] is None:
            self.table[i]=(key,h);self.used+=1;self.fill+=1
    def insert(self,key):
        before=self.used;h=self.hashfn(key,self.bits);self._insert(key,h)
        if self.used>before and self.fill*3 >= (self.mask+1)*2:
            self.resize((2 if self.used>50000 else 4)*self.used)
    def keys(self):return [e[0] for e in self.table if e is not None]

base=Path('/mnt/data/work_b104f/B104F_20260829T140333Z/retained_inputs')
rowmap={}
with open(base/'graphsage_row_to_entrez_topology_features.csv') as f:
    for r in csv.DictReader(f):rowmap[int(r['graphsage_row'])]=int(r['entrez_gene_id'])
summary=json.load(open('/mnt/data/work_core/ppi_repro_corrected/results/core_verification_summary.json'))
bounds=summary['partition']['bounds']; tissues=summary['partition']['tissues']
root=Path('/mnt/data/ohm_extract/bio-tissue-networks')

variants={}
for i,t in enumerate(tissues):
    seq=[];seen=set();edges=[]
    with open(root/f'{t}.edgelist') as f:
        for line in f:
            p=line.split()
            if len(p)<2:continue
            u,v=map(int,p[:2]);edges.append((u,v))
            for x in (u,v):
                if x not in seen:seen.add(x);seq.append(x)
    d_int=Py2Dict(py2_hash_int)
    d_str=Py2Dict(py2_hash_str)
    for x in seq:d_int.insert(x);d_str.insert(str(x))
    orders={
        'first_occurrence':seq,
        'sorted_int':sorted(seq),
        'sorted_int_desc':sorted(seq,reverse=True),
        'py2_int_dict':[int(x) for x in d_int.keys()],
        'py2_str_dict':[int(x) for x in d_str.keys()],
    }
    a,b=bounds[i],bounds[i+1]
    assert b-a==len(seq),(t,b-a,len(seq))
    for name,order in orders.items():
        rec=variants.setdefault(name,{'match':0,'compared':0,'per':[]})
        m=c=0
        for off,g in enumerate(order):
            row=a+off
            if row in rowmap:
                c+=1;m+=int(rowmap[row]==g)
        rec['match']+=m;rec['compared']+=c;rec['per'].append((t,m,c,m/c if c else None))
for name,r in variants.items():
    print(name,r['match'],r['compared'],r['match']/r['compared'])
    print(' best tissues',sorted(r['per'],key=lambda x:x[3],reverse=True)[:5])
    print(' worst tissues',sorted(r['per'],key=lambda x:x[3])[:5])
