#!/usr/bin/env python3
"""Enumerate + hash the confirmatory target manifest under LOCKED data-only criteria
(no effect-size / performance input). Tightened IC50, salt-strip, re-pair enantiomers,
then require per target: >=80 unique structural pairs, >=40 Murcko scaffolds, >=5 docs,
fold support (>=50 scaffolds for 5-fold x>=10). Writes confirmatory_manifest.json + sha256."""
import sqlite3, re, json, hashlib, collections, os
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog("rdApp.*")
D=os.path.dirname(os.path.abspath(__file__))
con=sqlite3.connect("/Users/dylanashraf/Documents/Programming/Eudysmic Ratio Project/data/chembl_37.db"); cur=con.cursor()
lfc=rdMolStandardize.LargestFragmentChooser()
def strip_ms(i): return "/".join(p for p in i.split("/") if not (p and p[0] in "ms"))
def mflag(i):
    m=re.search(r"/m([01.]+)",i); return m.group(1) if m else None

print("pulling tightened IC50 activities (all targets)...", flush=True)
q="""SELECT a.molregno,a.assay_id,s.tid,s.doc_id
     FROM activities a JOIN assays s ON a.assay_id=s.assay_id
     WHERE a.standard_type='IC50' AND a.standard_relation='=' AND a.standard_units='nM'
       AND a.pchembl_value IS NOT NULL AND a.data_validity_comment IS NULL
       AND (a.potential_duplicate=0 OR a.potential_duplicate IS NULL)
       AND s.assay_type='B' AND s.confidence_score=9"""
acts=cur.execute(q).fetchall()
print(f"  activities: {len(acts):,}", flush=True)
# restrict to single-protein targets
sp=set(t for (t,) in cur.execute("SELECT tid FROM target_dictionary WHERE target_type='SINGLE PROTEIN'"))
acts=[r for r in acts if r[2] in sp]
mols=set(r[0] for r in acts)
print(f"  single-protein-target activities: {len(acts):,}  distinct molregnos: {len(mols):,}", flush=True)

# process chiral candidate structures once (parent salt-strip -> fully-defined-chiral -> gkey/mflag/scaffold)
info={}; ml=list(mols); done=0
for i in range(0,len(ml),900):
    qq="SELECT molregno,canonical_smiles,standard_inchi FROM compound_structures WHERE molregno IN (%s)"%",".join("?"*len(ml[i:i+900]))
    for mr,smi,inchi in cur.execute(qq,ml[i:i+900]):
        done+=1
        if not inchi or "/t" not in inchi or "/i" in inchi or "/m" not in inchi: continue
        m=Chem.MolFromSmiles(smi) if smi else None
        if m is None: continue
        m=lfc.choose(m)
        try: pin=Chem.MolToInchi(m)
        except: continue
        if not pin or "/t" not in pin or "/m" not in pin or "/i" in pin: continue
        si=Chem.FindPotentialStereo(m)
        if len(si)==0 or any(e.specified==Chem.StereoSpecified.Unspecified for e in si): continue
        try: scf=MurckoScaffold.MurckoScaffoldSmiles(mol=m,includeChirality=False)
        except: scf=""
        info[mr]=(strip_ms(pin),mflag(pin),scf)
    if done % 90000 < 900: print(f"  structures processed ~{done:,}/{len(ml):,}", flush=True)
print(f"  eligible chiral molregnos: {len(info):,}", flush=True)

# pair within (tid, assay): same gkey, >=2 mflags; one structural pair per (tid,gkey)
pairs_by_tid=collections.defaultdict(set)   # tid -> set of gkey
scaf_by_tid=collections.defaultdict(set)
docs_by_tid=collections.defaultdict(set)
byassay=collections.defaultdict(list)
for mr,aid,tid,doc in acts:
    if mr in info: byassay[(tid,aid,doc)].append(mr)
for (tid,aid,doc),mrs in byassay.items():
    gk=collections.defaultdict(set)
    for mr in mrs: gk[info[mr][0]].add(info[mr][1])
    for gkey,flags in gk.items():
        if len(flags)>=2:
            pairs_by_tid[tid].add(gkey)
            scaf_by_tid[tid].add(info[[m for m in mrs if info[m][0]==gkey][0]][2])
            docs_by_tid[tid].add(doc)

# apply LOCKED criteria
qual=[]
for tid in pairs_by_tid:
    npairs=len(pairs_by_tid[tid]); nscaf=len(scaf_by_tid[tid]); ndoc=len(docs_by_tid[tid])
    if npairs>=80 and nscaf>=40 and ndoc>=5 and nscaf>=50:   # last = fold support (5-fold x >=10 scaffolds)
        qual.append((tid,npairs,nscaf,ndoc))
qual.sort(key=lambda x:-x[1])
names={}
tl=[q[0] for q in qual]
for i in range(0,len(tl),900):
    for tid,nm in cur.execute("SELECT tid,pref_name FROM target_dictionary WHERE tid IN (%s)"%",".join("?"*len(tl[i:i+900])),tl[i:i+900]):
        names[tid]=nm
manifest={"criteria":{"n_pairs>=":80,"n_scaffolds>=":40,"n_docs>=":5,"fold_support_scaffolds>=":50,
          "base":"IC50,=,nM,binding,conf9,valid,no-dup,single-protein,salt-stripped,fully-defined-chiral enantiomers"},
          "n_qualifying":len(qual),
          "targets":[{"tid":t,"pairs":p,"scaffolds":s,"docs":d,"name":names.get(t,"?")} for t,p,s,d in qual]}
payload=json.dumps([[t,p,s,d] for t,p,s,d in qual],sort_keys=True)
manifest["manifest_sha256"]=hashlib.sha256(payload.encode()).hexdigest()
json.dump(manifest,open(os.path.join(D,"confirmatory_manifest.json"),"w"),indent=1)
print(f"\n===== CONFIRMATORY MANIFEST =====")
print(f"qualifying targets: {len(qual)}")
print(f"manifest sha256: {manifest['manifest_sha256']}")
print(f"{'tid':>7} {'pairs':>5} {'scaf':>4} {'docs':>4}  name")
for t,p,s,d in qual[:60]:
    print(f"{t:>7} {p:>5} {s:>4} {d:>4}  {names.get(t,'?')[:40]}")
con.close()
