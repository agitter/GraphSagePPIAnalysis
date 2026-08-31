#!/usr/bin/env python3
from __future__ import annotations
import collections, csv, datetime as dt, gzip, hashlib, itertools, json, math, os, pickle, re, shutil, sqlite3, sys, tarfile, traceback, zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable
import numpy as np
import pandas as pd

ROOT=Path('/mnt/data');BASE=ROOT/'ppi_repro';DL=BASE/'downloads';EXT=BASE/'extracted'/'label_search';OUT=BASE/'results';EXT.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)

# ---------- general utilities ----------
def sha(p:Path):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def open_text(p:Path):
 if p.read_bytes()[:2]==b'\x1f\x8b':return gzip.open(p,'rt',errors='replace')
 return p.open('rt',errors='replace')
def safe_tar(src:Path,dst:Path):
 marker=dst/'.complete'
 if marker.exists():return
 dst.mkdir(parents=True,exist_ok=True)
 with tarfile.open(src,'r:*') as t:
  for m in t.getmembers():
   if not str((dst/m.name).resolve()).startswith(str(dst.resolve())):raise ValueError('unsafe tar')
  t.extractall(dst)
 marker.write_text(sha(src))
def safe_zip(src:Path,dst:Path):
 marker=dst/'.complete'
 if marker.exists():return
 dst.mkdir(parents=True,exist_ok=True)
 with zipfile.ZipFile(src) as z:z.extractall(dst)
 marker.write_text(sha(src))
def date_num(s:str|None):
 if not s:return 0
 m=re.search(r'(20\d{2})[-_]?([01]\d)?[-_]?([0-3]\d)?',s)
 if not m:return 0
 y=int(m.group(1));mo=int(m.group(2) or 6);d=int(m.group(3) or 15);return y*372+mo*31+d

# ---------- benchmark labels ----------
coll=OUT/'collapsed_gene_labels.csv'
if not coll.exists():raise FileNotFoundError('Run core_reproduction.py first: collapsed_gene_labels.csv absent')
df=pd.read_csv(coll)
label_cols=[c for c in df if c.startswith('label_')]
genes=[int(x) for x in df.entrez_gene_id]; gene_index={g:i for i,g in enumerate(genes)}; n=len(genes);L=len(label_cols)
Y=df[label_cols].to_numpy(dtype=np.uint8)
label_bits=[]
for j in range(L):
 b=0
 for i,v in enumerate(Y[:,j]):
  if v:b|=1<<i
 label_bits.append(b)
label_support=np.array([b.bit_count() for b in label_bits],dtype=int)
occurrence=dict(zip(df.entrez_gene_id.astype(int),df.occurrences.astype(int)))

# ---------- ontology ----------
@dataclass
class Ontology:
 source:str
 names:dict[str,str]
 namespace:dict[str,str]
 alt:dict[str,str]
 parents_is_a:dict[str,set[str]]
 parents_part_of:dict[str,set[str]]
 _cache:dict[tuple[str,str],frozenset[str]]
 def canon(self,t):return self.alt.get(t,t)
 def ancestors(self,t,mode):
  t=self.canon(t);key=(t,mode)
  if key in self._cache:return self._cache[key]
  seen={t};stack=[t]
  while stack:
   x=stack.pop();ps=set(self.parents_is_a.get(x,()))
   if mode=='is_a_part_of':ps|=self.parents_part_of.get(x,set())
   for p in ps:
    p=self.canon(p)
    if p not in seen:seen.add(p);stack.append(p)
  self._cache[key]=frozenset(seen);return self._cache[key]

def parse_obo(p:Path)->Ontology:
 names={};ns={};alt={};ia=collections.defaultdict(set);po=collections.defaultdict(set)
 current=None;obsolete=False
 with p.open(errors='replace') as f:
  for line in f:
   line=line.rstrip('\n')
   if line=='[Term]':current=None;obsolete=False;continue
   if line.startswith('['):current=None;continue
   if line.startswith('id: GO:'):current=line.split('id: ',1)[1].strip();continue
   if current is None:continue
   if line.startswith('is_obsolete: true'):obsolete=True;continue
   if line.startswith('name: '):names[current]=line[6:].strip()
   elif line.startswith('namespace: '):ns[current]=line[11:].strip()
   elif line.startswith('alt_id: '):alt[line[8:].strip()]=current
   elif line.startswith('is_a: '):ia[current].add(line[6:].split()[0])
   elif line.startswith('relationship: part_of '):po[current].add(line.split()[2])
 # leave obsolete records present; annotations to alt IDs canonicalize.
 return Ontology(str(p),names,ns,alt,dict(ia),dict(po),{})

obo_files=[]
for p in (DL/'go_releases').rglob('go-basic.obo') if (DL/'go_releases').exists() else []:
 try:
  head=p.open(errors='replace').read(200)
  if 'format-version:' in head:obo_files.append(p)
 except:pass
ontologies={p.parent.name:parse_obo(p) for p in obo_files}
# parent is date directory for .../date/go-basic.obo
ontologies={p.parent.name:parse_obo(p) for p in obo_files}
if not ontologies:
 # Try uploaded/downloaded go.obo files.
 for p in (DL/'go_releases').rglob('go.obo') if (DL/'go_releases').exists() else []:
  try:
   if 'format-version:' in p.open(errors='replace').read(200):ontologies[p.parent.name]=parse_obo(p)
  except:pass

def nearest_ontology(date_hint:str|None):
 if not ontologies:return None
 target=date_num(date_hint);return ontologies[min(ontologies,key=lambda k:abs(date_num(k)-target))]

# ---------- identifier maps ----------
# MSigDB v5.2 chip files, used only as an explicit symbol fallback.
chiproot=EXT/'msigdb_chip_v5_2'
if not chiproot.exists():safe_zip(ROOT/'msigdb_v5.2_chip_files_to_download_locally.zip',chiproot)
# recursively unzip nested archives
for p in list(chiproot.rglob('*.zip')):
 d=p.with_suffix('')
 if not (d/'.complete').exists():
  try:safe_zip(p,d)
  except:pass
symbol_to_gene=collections.defaultdict(set);gene_to_symbols=collections.defaultdict(set)
chip_schema=[]
for p in chiproot.rglob('*.chip'):
 try:
  tab=pd.read_csv(p,sep='\t',dtype=str,comment='#').fillna('')
 except:continue
 cols={c.lower().strip():c for c in tab.columns}
 gcol=next((c for k,c in cols.items() if ('entrez' in k and ('gene' in k or 'id' in k)) or k in ('gene id','gene_id')),None)
 scol=next((c for k,c in cols.items() if 'gene symbol' in k or k=='symbol'),None)
 acols=[c for k,c in cols.items() if 'alias' in k or 'synonym' in k]
 chip_schema.append({'path':str(p),'columns':list(tab.columns),'gene_col':gcol,'symbol_col':scol,'alias_cols':acols})
 if not gcol or not scol:continue
 for _,r in tab.iterrows():
  nums=re.findall(r'\d+',str(r[gcol]))
  syms={str(r[scol]).strip()} if str(r[scol]).strip() else set()
  for c in acols:syms|={x.strip() for x in re.split(r'[|,;/]+',str(r[c])) if x.strip()}
  for x in nums:
   g=int(x);gene_to_symbols[g]|=syms
   for s in syms:symbol_to_gene[s.upper()].add(g)
(OUT/'msigdb_chip_schema.json').write_text(json.dumps(chip_schema,indent=2))

# Bioconductor databases and UniProt maps.
@dataclass
class AnnSource:
 name:str
 source_type:str
 date_hint:str
 records:list[tuple[int,str,str]] # gene, GO, evidence
 provenance:dict[str,Any]
 prepropagated:bool=False

ann_sources=[];bioc_uniprot={};bioc_meta={};schemas=[]
bioc_files=[]
for root in [DL/'bioconductor',DL/'bioconductor_mirror']:
 if root.exists():bioc_files.extend(root.glob('org.Hs.eg.db_*.tar.gz'))
seen_sha=set()
for tar in sorted(bioc_files):
 try:h=sha(tar)
 except:continue
 if h in seen_sha:continue
 seen_sha.add(h)
 ver=re.search(r'_([0-9.]+)\.tar\.gz$',tar.name).group(1)
 dst=EXT/f'bioc_{ver}';
 try:safe_tar(tar,dst)
 except Exception:continue
 sqls=list(dst.rglob('*.sqlite'))
 if not sqls:continue
 sp=sqls[0];con=sqlite3.connect(sp);tables=[r[0] for r in con.execute("select name from sqlite_master where type='table'")]
 schema={'version':ver,'tar':str(tar),'sqlite':str(sp),'tables':{}}
 for t in tables:
  schema['tables'][t]=[r[1] for r in con.execute(f'pragma table_info("{t}")')]
 meta={}
 if 'metadata' in tables:
  try:meta=dict(con.execute('select name,value from metadata').fetchall())
  except:pass
 bioc_meta[ver]=meta;schema['metadata']=meta;schemas.append(schema)
 # gene and GO tables
 if 'genes' in tables:
  for tab,pre in [('go',False),('go_all',True)]:
   if tab not in tables:continue
   q=f'SELECT genes.gene_id,{tab}.go_id,{tab}.evidence,{tab}.ontology FROM genes JOIN {tab} USING (_id)'
   rec=[]
   try:
    for g,go,ev,onto in con.execute(q):
     try:gi=int(g)
     except:continue
     if gi in gene_index and str(onto).upper().startswith('BP'):rec.append((gi,str(go),str(ev)))
   except Exception:continue
   ann_sources.append(AnnSource(f'org.Hs.eg.db_{ver}:{tab}','Bioconductor',meta.get('GOSOURCEDATE',meta.get('EGSOURCEDATE',ver)),rec,
                                {'tar':str(tar),'sqlite':str(sp),'metadata':meta,'table':tab},pre))
  if 'uniprot' in tables:
   mp=collections.defaultdict(set)
   try:
    for g,u in con.execute('SELECT genes.gene_id,uniprot.uniprot_id FROM genes JOIN uniprot USING (_id)'):
     try:gi=int(g)
     except:continue
     if gi in gene_index:
      us=str(u);mp[us].add(gi);mp[us.split('-')[0]].add(gi)
   except:pass
   bioc_uniprot[meta.get('EGSOURCEDATE', meta.get('GOSOURCEDATE', {'3.0.0':'2014-09-27','3.1.2':'2015-03-17','3.2.3':'2015-09-27','3.3.0':'2016-03-14','3.4.0':'2016-09-26'}.get(ver,ver)))]=mp
 con.close()
(OUT/'bioconductor_sqlite_schemas.json').write_text(json.dumps(schemas,indent=2))

# GPI and UniProt mapping helpers.
def parse_gpi(p:Path):
 acc=collections.defaultdict(set);sym=collections.defaultdict(set);rows=0
 try:f=open_text(p)
 except:return acc,sym,0
 with f:
  for line in f:
   if not line.strip() or line.startswith('!'):continue
   a=line.rstrip('\n').split('\t');rows+=1
   ids={int(x) for x in re.findall(r'(?:GeneID|NCBI_Gene|EntrezGene)[:=](\d+)',line)} & set(gene_index)
   if not ids:continue
   obj=a[1].strip() if len(a)>1 else ''
   obj=obj.split(':')[-1]; keys={obj,obj.split('-')[0]}
   if len(a)>2 and a[2].strip():
    for g in ids:sym[a[2].strip().upper()].add(g)
   for k in keys:
    if k:
     acc[k]|=ids
 return acc,sym,rows

gpi_by_release={};gpi_meta=[]
for p in sorted((DL/'ebi_goa').glob('*')) if (DL/'ebi_goa').exists() else []:
 m=re.search(r'(?:gp_information\.goa_(?:ref_)?human|goa_human\.gpi)\.(\d+)\.gz$',p.name)
 if not m:continue
 r=int(m.group(1));mp,sm,rows=parse_gpi(p)
 # full file preferred over ref file
 priority=0 if 'ref_' in p.name else 1
 old=gpi_by_release.get(r)
 if old is None or priority>old[3]:gpi_by_release[r]=(mp,sm,p,priority)
 gpi_meta.append({'release':r,'path':str(p),'rows':rows,'mapped_accessions':len(mp),'symbols':len(sm),'priority':priority})
(OUT/'gpi_mapping_inventory.csv').write_text(pd.DataFrame(gpi_meta).to_csv(index=False))

# UniProt downloaded mapping files.
uniprot_by_rel={};uniprot_meta=[]
for p in (DL/'uniprot').rglob('*') if (DL/'uniprot').exists() else []:
 if not p.is_file() or p.suffix not in ('.gz','.dat','.tab','.txt'):continue
 if p.stat().st_size>600_000_000:continue
 relm=re.search(r'release-(20\d{2}_\d{2})',str(p));rel=relm.group(1) if relm else p.parent.name
 mp=collections.defaultdict(set);rows=0
 try:f=open_text(p)
 except:continue
 try:
  with f:
   for line in f:
    if not line.strip() or line.startswith('#'):continue
    a=line.rstrip('\n').split('\t');rows+=1
    if len(a)>=3 and a[1]=='GeneID':
     try:g=int(a[2].split(';')[0])
     except:continue
     if g in gene_index:
      mp[a[0]].add(g);mp[a[0].split('-')[0]].add(g)
    elif len(a)>=3 and ('selected' in p.name.lower() or len(a)>10):
     for x in re.findall(r'\d+',a[2]):
      g=int(x)
      if g in gene_index:mp[a[0]].add(g);mp[a[0].split('-')[0]].add(g)
 except:continue
 if mp:
  uniprot_by_rel[rel]=mp;uniprot_meta.append({'release':rel,'path':str(p),'rows':rows,'mapped_accessions':len(mp),'size':p.stat().st_size})
pd.DataFrame(uniprot_meta).to_csv(OUT/'uniprot_mapping_inventory.csv',index=False)

all_bioc=collections.defaultdict(set)
for mp in bioc_uniprot.values():
 for k,v in mp.items():all_bioc[k]|=v
all_uniprot=collections.defaultdict(set)
for mp in uniprot_by_rel.values():
 for k,v in mp.items():all_uniprot[k]|=v
for k,v in all_bioc.items():all_uniprot[k]|=v

release_dates={140:'2015-01-05',141:'2015-02-02',142:'2015-03-02',143:'2015-03-30',144:'2015-04-27',145:'2015-05-26',146:'2015-06-22',147:'2015-07-20',148:'2015-09-14',149:'2015-10-12',150:'2015-11-09',151:'2015-12-07',152:'2016-01-04',153:'2016-01-20',154:'2016-02-15',155:'2016-03-14',156:'2016-04-11',157:'2016-05-09',158:'2016-06-07',159:'2016-07-04',160:'2016-09-14',161:'2016-10-03',162:'2016-10-31',163:'2016-11-28'}

def nearest_map(mapping_by_key,date_hint):
 if not mapping_by_key:return {}
 key=min(mapping_by_key,key=lambda k:abs(date_num(str(k))-date_num(date_hint)));return mapping_by_key[key]

def parse_gaf_mapped(p:Path,date_hint:str,release:int|None,map_mode:str):
 gpi,gsym=(gpi_by_release.get(release,({}, {},None,0))[0:2] if release is not None else ({},{}))
 bmap=nearest_map(bioc_uniprot,date_hint);umap=nearest_map(uniprot_by_rel,date_hint)
 rec=[];stats=collections.Counter();seen=set()
 try:f=open_text(p)
 except:return rec,dict(stats)
 with f:
  for line in f:
   if not line.strip() or line.startswith('!'):continue
   a=line.rstrip('\n').split('\t')
   if len(a)<9:continue
   db,obj,symbol,qual,go,ev,aspect=a[0],a[1],a[2],a[3],a[4],a[6],a[8]
   stats['rows']+=1
   if aspect!='P' or 'NOT' in qual.split('|'):continue
   ids=set();obj0=obj.split(':')[-1];keys={obj0,obj0.split('-')[0]}
   if db in ('NCBI_Gene','GeneID') and obj0.isdigit():ids.add(int(obj0))
   if map_mode=='gpi_only':
    for k in keys:ids|=gpi.get(k,set())
   elif map_mode=='release_specific':
    for k in keys:ids|=gpi.get(k,set())|bmap.get(k,set())|umap.get(k,set())
    if not ids:ids|=gsym.get(symbol.upper(),set())|symbol_to_gene.get(symbol.upper(),set())
   elif map_mode=='all_union':
    for k in keys:ids|=gpi.get(k,set())|all_uniprot.get(k,set())
    ids|=gsym.get(symbol.upper(),set())|symbol_to_gene.get(symbol.upper(),set())
   elif map_mode=='symbol':ids|=gsym.get(symbol.upper(),set())|symbol_to_gene.get(symbol.upper(),set())
   ids&=set(gene_index)
   if not ids:stats['unmapped_rows']+=1;continue
   stats['mapped_rows']+=1;stats['mapped_genes']+=len(ids)
   for g in ids:
    tup=(g,go,ev)
    if tup not in seen:seen.add(tup);rec.append(tup)
 return rec,dict(stats)

# EBI GAF sources.
gaf_meta=[]
for p in sorted((DL/'ebi_goa').glob('*.gz')) if (DL/'ebi_goa').exists() else []:
 m=re.search(r'(?:gene_association\.goa_human|goa_human\.gaf)\.(\d+)\.gz$',p.name)
 if not m:continue
 r=int(m.group(1));date=release_dates.get(r,str(r))
 modes=['gpi_only','release_specific','all_union']
 for mode in modes:
  rec,st=parse_gaf_mapped(p,date,r,mode)
  ann_sources.append(AnnSource(f'GOA_r{r}:{mode}','EBI_GOA',date,rec,{'path':str(p),'release':r,'mapping':mode,'stats':st}))
  gaf_meta.append({'source':f'GOA_r{r}:{mode}','path':str(p),'release':r,'date':date,'records':len(rec),**st})
pd.DataFrame(gaf_meta).to_csv(OUT/'gaf_mapping_coverage.csv',index=False)

# GO release GAFs (deduplicate by file SHA against EBI files; map using all_union or direct GeneID).
seen_gaf_sha={sha(p) for p in (DL/'ebi_goa').glob('*gaf*.gz') if p.is_file()}
for p in sorted((DL/'go_releases').rglob('*.gaf.gz')) if (DL/'go_releases').exists() else []:
 try:h=sha(p)
 except:continue
 if h in seen_gaf_sha:continue
 date=p.parent.name if re.match(r'20\d{2}-\d{2}-\d{2}',p.parent.name) else p.parent.parent.name
 rec,st=parse_gaf_mapped(p,date,None,'all_union')
 if rec:ann_sources.append(AnnSource(f'GO_release_{date}:all_union','GO_RELEASE',date,rec,{'path':str(p),'mapping':'all_union','stats':st}))

# NCBI gene2go snapshots, including any supplied/cloned files found recursively.
def parse_gene2go(p:Path):
 rec=[];stats=collections.Counter();seen=set()
 try:
  with open_text(p) as f: lines=f.readlines()
 except:return rec,dict(stats)
 try:
  data=[x.rstrip('\n') for x in lines if x.strip()]
  header_line=next((x.lstrip('#') for x in data if 'GeneID' in x and ('GO_ID' in x or 'GO ID' in x)),None)
  if header_line:
   headers=header_line.split('\t'); idx={h.strip().lower().replace(' ','_'):i for i,h in enumerate(headers)}
   start=next(i for i,x in enumerate(data) if x.lstrip('#')==header_line)+1
   def col(*names):
    return next((idx[n] for n in names if n in idx),None)
   ig=col('geneid','gene_id','entrez_gene_id'); igo=col('go_id','goid'); iev=col('evidence','evidence_code'); iq=col('qualifier'); itax=col('tax_id','taxid'); icat=col('category','ontology','aspect')
   for line in data[start:]:
    if line.startswith('#'):continue
    a=line.split('\t')
    try:g=int(a[ig])
    except:continue
    if itax is not None:
     try:
      if int(a[itax])!=9606:continue
     except:continue
    if g not in gene_index or igo is None:continue
    go=a[igo];ev=a[iev] if iev is not None and iev<len(a) else '';qual=a[iq] if iq is not None and iq<len(a) else '';cat=a[icat] if icat is not None and icat<len(a) else 'Process'
    if cat and 'Process' not in cat and cat not in ('P','BP','biological_process'):continue
    if 'NOT' in qual.split('|'):continue
    stats['rows']+=1;tu=(g,go,ev)
    if tu not in seen:seen.add(tu);rec.append(tu)
  else:
   for line in data:
    if line.startswith('#'):continue
    a=line.split('\t')
    if len(a)<7:continue
    try:tax=int(a[0]);g=int(a[1])
    except:continue
    if tax!=9606 or g not in gene_index:continue
    go=a[2];ev=a[3];qual=a[4];cat=a[7] if len(a)>7 else a[-1]
    if 'Process' not in cat and cat not in ('P','BP','biological_process'):continue
    if 'NOT' in qual.split('|'):continue
    stats['rows']+=1;tu=(g,go,ev)
    if tu not in seen:seen.add(tu);rec.append(tu)
 except Exception as e:stats['error']=repr(e)
 return rec,dict(stats)

gene2go_meta=[];seen_hash=set()
search_roots=[DL/'ncbi_wayback',DL/'entrez-gene',DL/'repos',DL/'gene2go_git_history']
for root in search_roots:
 if not root.exists():continue
 for p in root.rglob('*gene2go*'):
  if not p.is_file() or p.stat().st_size<1000:continue
  try:h=sha(p)
  except:continue
  if h in seen_hash:continue
  seen_hash.add(h);rec,st=parse_gene2go(p)
  if rec:
   dm=re.search(r'(20\d{6,12}|20\d{2}[-_]\d{2}[-_]\d{2})',p.name);date=dm.group(1) if dm else str(p.stat().st_mtime_ns)
   name=f'gene2go:{p.name}'
   ann_sources.append(AnnSource(name,'NCBI_gene2go',date,rec,{'path':str(p),'sha256':h,'stats':st}))
   gene2go_meta.append({'name':name,'path':str(p),'sha256':h,'records':len(rec),**st})
pd.DataFrame(gene2go_meta).to_csv(OUT/'gene2go_inventory.csv',index=False)

# ---------- evidence and term-set evaluation ----------
EXP6={'EXP','IDA','IPI','IMP','IGI','IEP'};EXPHT={'HTP','HDA','HMP','HGI','HEP'};PHYLO={'IBA','IBD','IKR','IRD'};COMP={'ISS','ISO','ISA','ISM','IGC','RCA'};AUTHOR={'TAS','NAS'};CURATOR={'IC'};ELECT={'IEA'};NODATA={'ND'}
all_codes=sorted(set(ev for s in ann_sources for _,_,ev in s.records))
filters={
 'ALL':set(all_codes),
 'NO_IEA':set(all_codes)-ELECT,
 'NO_IEA_ND':set(all_codes)-ELECT-NODATA,
 'NO_ND':set(all_codes)-NODATA,
 'EXP6':EXP6,
 'EXP_ALL':EXP6|EXPHT,
 'EXP_AUTHOR':EXP6|EXPHT|AUTHOR,
 'EXP_PHYLO_AUTHOR':EXP6|EXPHT|PHYLO|AUTHOR,
 'CURATED_BROAD':set(all_codes)-ELECT-NODATA,
}

def build_term_bits(src:AnnSource,allowed:set[str],ontology:Ontology|None,mode:str):
 direct=collections.defaultdict(int)
 for g,go,ev in src.records:
  if ev not in allowed:continue
  direct[go]|=1<<gene_index[g]
 if src.prepropagated or mode=='direct' or ontology is None:return dict(direct)
 out=collections.defaultdict(int)
 for go,b in direct.items():
  for a in ontology.ancestors(go,mode):
   if ontology.namespace.get(a) in (None,'biological_process'):out[a]|=b
 return dict(out)

def evaluate_bits(bits:dict[str,int]):
 terms=[t for t,b in bits.items() if b]
 if not terms:return {'terms':0,'exact':0,'n99':0,'n95':0,'total_mismatch':n*L,'best':[],'supports':{}}
 mm=np.empty((len(terms),L),dtype=np.uint16);supports=np.empty(len(terms),dtype=np.uint16)
 for i,t in enumerate(terms):
  b=bits[t];supports[i]=b.bit_count();mm[i]=[(b^lb).bit_count() for lb in label_bits]
 idx=np.argmin(mm,axis=0);bestm=mm[idx,np.arange(L)].astype(int)
 best=[]
 for j,i in enumerate(idx):
  m=int(bestm[j]);best.append({'column':j,'term':terms[int(i)],'mismatches':m,'agreement':1-m/n,'term_support':int(supports[int(i)]),'label_support':int(label_support[j])})
 return {'terms':len(terms),'exact':int(np.sum(bestm==0)),'n99':int(np.sum(bestm<=math.floor(.01*n))),
         'n995':int(np.sum(bestm<=math.floor(.005*n))),'n95':int(np.sum(bestm<=math.floor(.05*n))),
         'total_mismatch':int(bestm.sum()),'median_mismatch':float(np.median(bestm)),'best':best,
         'term_count_ge15':int(np.sum(supports>=15)),'term_count_ge500':int(np.sum(supports>=500))}

screen=[];bestmaps={};bit_cache={}
for si,src in enumerate(ann_sources):
 onto=nearest_ontology(src.date_hint)
 modes=['direct'] if src.prepropagated else ['direct','is_a','is_a_part_of']
 for mode in modes:
  for fname,allowed in filters.items():
   bits=build_term_bits(src,allowed,onto,mode);ev=evaluate_bits(bits)
   key=f'{src.name}|{mode}|{fname}|{onto.source if onto else "none"}'
   row={'source':src.name,'source_type':src.source_type,'date_hint':src.date_hint,'mapping':src.provenance.get('mapping',''),
        'mode':mode,'filter':fname,'ontology':onto.source if onto else '', 'annotation_records':len(src.records),
        **{k:v for k,v in ev.items() if k not in ('best','supports')}}
   screen.append(row);bestmaps[key]=ev['best']
 # checkpoint every source
 if si%5==0:pd.DataFrame(screen).to_csv(OUT/'label_source_screen_partial.csv',index=False)
screen_df=pd.DataFrame(screen).sort_values(['exact','n99','n995','n95','total_mismatch'],ascending=[False,False,False,False,True]);screen_df.to_csv(OUT/'label_source_screen.csv',index=False)

# Save per-column mappings for the top 30 screen candidates.
top_keys=[];mapping_rows=[]
for _,r in screen_df.head(30).iterrows():
 key=f"{r['source']}|{r['mode']}|{r['filter']}|{r['ontology']}";top_keys.append(key)
 for b in bestmaps[key]:mapping_rows.append({'candidate_rank':len(top_keys),'source':r['source'],'mode':r['mode'],'filter':r['filter'],'ontology':r['ontology'],**b})
pd.DataFrame(mapping_rows).to_csv(OUT/'label_top_candidate_term_maps.csv',index=False)

# Refine top distinct annotation sources with all subsets of broad evidence groups, is_a nearest ontology.
groups={'EXP':EXP6|EXPHT,'PHYLO':PHYLO,'COMP':COMP,'AUTHOR':AUTHOR,'CURATOR':CURATOR,'ELECT':ELECT,'NODATA':NODATA}
top_source_names=[]
for x in screen_df.source:
 if x not in top_source_names:top_source_names.append(x)
 if len(top_source_names)>=8:break
source_by_name={s.name:s for s in ann_sources};refine=[];oracle={}
for sname in top_source_names:
 src=source_by_name[sname];onto=nearest_ontology(src.date_hint);per_col=np.full(L,n+1,dtype=int);per_col_cfg=['']*L
 for mask in range(1,1<<len(groups)):
  selected=[k for i,k in enumerate(groups) if mask&(1<<i)];allowed=set().union(*(groups[k] for k in selected))
  bits=build_term_bits(src,allowed,onto,'direct' if src.prepropagated else 'is_a');ev=evaluate_bits(bits)
  row={'source':sname,'groups':'+'.join(selected),'mode':'direct' if src.prepropagated else 'is_a','ontology':onto.source if onto else '',
       **{k:v for k,v in ev.items() if k!='best'}};refine.append(row)
  for b in ev['best']:
   if b['mismatches']<per_col[b['column']]:per_col[b['column']]=b['mismatches'];per_col_cfg[b['column']]=row['groups']
 oracle[sname]={'oracle_exact':int(np.sum(per_col==0)),'oracle_n99':int(np.sum(per_col<=math.floor(.01*n))),
                'oracle_total_mismatch':int(per_col.sum()),'per_column_min':per_col.tolist(),'per_column_filter':per_col_cfg}
pd.DataFrame(refine).sort_values(['exact','n99','total_mismatch'],ascending=[False,False,True]).to_csv(OUT/'label_evidence_group_refinement.csv',index=False)
(OUT/'label_filter_oracle.json').write_text(json.dumps(oracle,indent=2))

# Test alternative ontology snapshots and propagation modes for the top five source/filter pairs.
refdf=pd.DataFrame(refine).sort_values(['exact','n99','total_mismatch'],ascending=[False,False,True]);onto_ref=[];done=set()
onto_choices=list(ontologies.items())
for _,r in refdf.iterrows():
 sig=(r['source'],r['groups'])
 if sig in done:continue
 done.add(sig)
 if len(done)>5:break
 src=source_by_name[r['source']];allowed=set().union(*(groups[k] for k in r['groups'].split('+')))
 # nearest plus all 2016-06 through 2016-10 snapshots, if present
 choices=[];near=nearest_ontology(src.date_hint)
 if near:choices.append(near)
 for k,o in onto_choices:
  if re.match(r'2016-(0[6-9]|10)-',k):choices.append(o)
 unique={o.source:o for o in choices}.values()
 for o in unique:
  for mode in (['direct'] if src.prepropagated else ['direct','is_a','is_a_part_of']):
   ev=evaluate_bits(build_term_bits(src,allowed,o,mode))
   onto_ref.append({'source':src.name,'groups':r['groups'],'mode':mode,'ontology':o.source,**{k:v for k,v in ev.items() if k!='best'}})
pd.DataFrame(onto_ref).sort_values(['exact','n99','total_mismatch'],ascending=[False,False,True]).to_csv(OUT/'label_ontology_refinement.csv',index=False)

# ---------- non-GOA candidate products ----------
external=[];external_maps={}
def add_external(name,typ,bits,notes):
 ev=evaluate_bits(bits);external.append({'source':name,'source_type':typ,'notes':notes,**{k:v for k,v in ev.items() if k!='best'}});external_maps[name]=ev['best']

# OhmNet labels: union over all tissues and over only selected GraphSAGE tissues.
selected_tissues=set(pd.read_csv(OUT/'tissue_matches.csv').tissue) if (OUT/'tissue_matches.csv').exists() else set()
labroot=BASE/'extracted'/'ohmnet_labels'
if not labroot.exists():
 try:safe_tar(ROOT/'bio-tissue-labels.tar.gz',labroot)
 except:pass
all_union=collections.defaultdict(int);sel_union=collections.defaultdict(int);all_inter=None;term_tissues=collections.defaultdict(set)
for p in labroot.rglob('*.lab') if labroot.exists() else []:
 m=re.match(r'^(.*?)_(GO:\d{7})\.lab$',p.name)
 if not m:continue
 tissue,go=m.groups();b=0;parsed_gene=False
 with p.open(errors='replace') as f:
  for line in f:
   a=line.strip().split()
   if len(a)>=2:
    try:g=int(float(a[0]));y=float(a[-1]);parsed_gene=True
    except:continue
    if y>.5 and g in gene_index:b|=1<<gene_index[g]
 if not parsed_gene:continue
 all_union[go]|=b;term_tissues[go].add(tissue)
 if tissue in selected_tissues:sel_union[go]|=b
add_external('OhmNet_label_union_all_tissues','OhmNet',dict(all_union),'positive for a GO term in any provided tissue-specific label file')
add_external('OhmNet_label_union_selected_24','OhmNet',dict(sel_union),'union over GraphSAGE-matched tissues only')

# MSigDB GMTs, Entrez and symbol-space via chip mapping.
def extract_msig(fn,label):
 d=EXT/label
 if not d.exists():safe_zip(ROOT/fn,d)
 for z in list(d.rglob('*.zip')):
  q=z.with_suffix('')
  if not (q/'.complete').exists():
   try:safe_zip(z,q)
   except:pass
 return d
for ver,fn in [('5.1','msigdb_v5.1_files_to_download_locally.zip'),('5.2','msigdb_v5.2_files_to_download_locally.zip'),('6.0','msigdb_v6.0_files_to_download_locally.zip')]:
 d=extract_msig(fn,'label_msigdb_'+ver.replace('.','_'));ent_bits={};sym_bits={}
 for p in d.rglob('*.gmt'):
  with p.open(errors='replace') as f:
   for line in f:
    a=line.rstrip('\n').split('\t')
    if len(a)<3:continue
    name=a[0]+'@@'+p.name;tokens=[x.strip() for x in a[2:] if x.strip()]
    nums=[int(x) for x in tokens if re.fullmatch(r'\d+',x)]
    if 'entrez' in p.name.lower() or (tokens and len(nums)/len(tokens)>.98):
     b=0
     for g in nums:
      if g in gene_index:b|=1<<gene_index[g]
     ent_bits[name]=b
    elif 'symbol' in p.name.lower():
     b=0
     st={x.upper() for x in tokens}
     for g in genes:
      if gene_to_symbols.get(g,set()) & st:b|=1<<gene_index[g]
     sym_bits[name]=b
 add_external(f'MSigDB_{ver}_Entrez_all_GMT','MSigDB',ent_bits,'direct Entrez membership across every supplied GMT')
 add_external(f'MSigDB_{ver}_symbol_via_v5.2_chip','MSigDB',sym_bits,'symbol GMT membership after historical chip mapping')

pd.DataFrame(external).sort_values(['exact','n99','total_mismatch'],ascending=[False,False,True]).to_csv(OUT/'label_external_sources.csv',index=False)
exrows=[]
for name,best in external_maps.items():
 for b in best:exrows.append({'source':name,**b})
pd.DataFrame(exrows).to_csv(OUT/'label_external_term_maps.csv',index=False)

# ---------- relation to Greene/OhmNet term lists and support-threshold selection ----------
def workbook_terms(path):
 ids=set();xl=pd.ExcelFile(path,engine='openpyxl')
 for s in xl.sheet_names:
  d=pd.read_excel(path,sheet_name=s,header=None,engine='openpyxl')
  for x in d.to_numpy().ravel():
   if pd.isna(x):continue
   ids.update(re.findall(r'GO:\d{7}',str(x)))
 return ids
greene6=workbook_terms(ROOT/'Greene2015_Table6.xlsx');greene9=workbook_terms(ROOT/'Greene2015_Table9.xlsx');ohmterms=set(all_union)
# Rebuild top screen candidate bitsets and inspect selection counts/restrictions.
selection=[]
for _,r in screen_df.head(20).iterrows():
 src=source_by_name[r['source']];onto=nearest_ontology(src.date_hint);bits=build_term_bits(src,filters[r['filter']],onto,r['mode'])
 # term list from its best mappings
 key=f"{r['source']}|{r['mode']}|{r['filter']}|{r['ontology']}";best=bestmaps[key];best_terms={x['term'] for x in best}
 for restr_name,restr in [('all',None),('Greene_Table6',greene6),('Greene_Table9',greene9),('OhmNet_label_terms',ohmterms),('Greene6_or_9',greene6|greene9)]:
  supports=[]
  for t,b in bits.items():
   if restr is not None and t not in restr:continue
   supports.append(b.bit_count())
  counts={f'ge_{k}':sum(x>=k for x in supports) for k in [5,10,15,20,25,50,100,200,300,400,500,600,750,1000]}
  selection.append({'source':r['source'],'mode':r['mode'],'filter':r['filter'],'restriction':restr_name,'terms':len(supports),
                    'best_terms_in_restriction':sum(t in restr for t in best_terms) if restr is not None else len(best_terms),**counts})
pd.DataFrame(selection).to_csv(OUT/'label_term_selection_thresholds.csv',index=False)

# ---------- report ----------
top=screen_df.head(15).to_dict('records');extop=pd.DataFrame(external).sort_values(['exact','n99','total_mismatch'],ascending=[False,False,True]).head(10).to_dict('records')
lines=['# Independent gene-to-GO label source search','',
       f'Benchmark matrix: **{n:,} independently recovered Entrez genes × {L} label columns**. Candidate annotations are compared as complete binary vectors, not merely by term names or support counts.','',
       '## Annotation inputs parsed','',
       f'- Historical Bioconductor annotation sources: **{sum(s.source_type=="Bioconductor" for s in ann_sources)}**.',
       f'- Historical EBI/GO release GAF mapping variants: **{sum(s.source_type in ("EBI_GOA","GO_RELEASE") for s in ann_sources)}**.',
       f'- Recovered NCBI `gene2go` snapshots: **{sum(s.source_type=="NCBI_gene2go" for s in ann_sources)}**.',
       f'- Historical ontology snapshots parsed: **{len(ontologies)}**.',
       f'- Evidence codes observed: `{all_codes}`.','',
       '## Best historical annotation combinations (screen)','',
       '| Rank | Source | Date | Mapping | Propagation | Evidence filter | Exact | ≥99% | ≥99.5% | ≥95% | Total mismatches |',
       '|---:|---|---|---|---|---|---:|---:|---:|---:|---:|']
for i,r in enumerate(top,1):
 lines.append(f"| {i} | `{r['source']}` | {r['date_hint']} | {r.get('mapping','')} | {r['mode']} | {r['filter']} | {r['exact']} | {r['n99']} | {r['n995']} | {r['n95']} | {r['total_mismatch']} |")
lines += ['', 'The full grid is in `label_source_screen.csv`; the best per-column GO IDs are in `label_top_candidate_term_maps.csv`.', '',
          '## Alternative products tested directly','',
          '| Source | Exact | ≥99% | ≥95% | Total mismatches | Interpretation |','|---|---:|---:|---:|---:|---|']
for r in extop:
 lines.append(f"| `{r['source']}` | {r['exact']} | {r['n99']} | {r['n95']} | {r['total_mismatch']} | {r['notes']} |")
lines += ['', '## Tests that go beyond the prior summary','',
          '- GOA accessions were remapped independently with three mapping policies: matching-release GPI only, release-specific GPI/UniProt/Bioconductor union with symbol fallback, and an all-snapshot union.',
          '- Direct annotations, `is_a` propagation, and `is_a + part_of` propagation were tested separately against dated GO ontology files.',
          '- Historical Bioconductor `go` and already-propagated `go_all` tables were read directly from their SQLite databases rather than through a current R installation.',
          '- Every broad evidence-group subset was tested for the strongest source families; `label_filter_oracle.json` reports the per-column oracle, which distinguishes “wrong global filter” from “associations absent from the source.”',
          '- OhmNet labels were tested both as individual tissue labels in the core analysis and as unions over all tissues or only the selected 24 tissues.',
          '- MSigDB labels were retested in both Entrez space and symbol space using the supplied historical chip files.',
          '- Greene Table 6, Greene Table 9, and the set of GO IDs present in the OhmNet label filenames were tested as term-selection restrictions; threshold counts are in `label_term_selection_thresholds.csv`.', '',
          '## Provenance and reproducibility','',
          '- `source_manifest.csv` records source URLs, acquisition status, byte sizes, and SHA-256 hashes.',
          '- `gaf_mapping_coverage.csv`, `gpi_mapping_inventory.csv`, `uniprot_mapping_inventory.csv`, and `gene2go_inventory.csv` record exactly which identifier mappings contributed.',
          '- `label_source_screen.csv` and refinement files contain every tested parameter combination and objective value.', '']
(OUT/'label_source_search_report.md').write_text('\n'.join(lines))
summary={'genes':n,'labels':L,'annotation_sources':len(ann_sources),'ontologies':len(ontologies),'top_screen':top,'top_external':extop,'filter_oracle':oracle}
(OUT/'label_source_search_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps({'status':'ok','sources':len(ann_sources),'top':top[0] if top else None}))
