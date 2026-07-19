#!/usr/bin/env python3
"""Freeze the 14-target confirmatory cohort PAIRS (one row per unique structural pair (tid,gkey),
Delta-y = median across that pair's assays). Reads confirmatory_manifest.json for the tids."""
import sqlite3, re, json, hashlib, csv, collections, statistics as st, os
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog("rdApp.*")
D=os.path.dirname(os.path.abspath(__file__))
con=sqlite3.connect("/Users/dylanashraf/Documents/Programming/Eudysmic Ratio Project/data/chembl_37.db"); cur=con.cursor()
lfc=rdMolStandardize.LargestFragmentChooser()
man=json.load(open(os.path.join(D,"confirmatory_manifest.json")))
tids={t["tid"]:t["name"] for t in man["targets"]}
def strip_ms(i): return "/".join(p for p in i.split("/") if not (p and p[0] in "ms"))
def mflag(i):
    m=re.search(r"/m([01.]+)",i); return m.group(1) if m else None
def parent_info(smi,inchi):
    if not inchi or "/t" not in inchi or "/i" in inchi or "/m" not in inchi: return None
    m=Chem.MolFromSmiles(smi) if smi else None
    if m is None: return None
    m=lfc.choose(m)
    pin=Chem.MolToInchi(m)
    if not pin or "/t" not in pin or "/m" not in pin or "/i" in pin: return None
    si=Chem.FindPotentialStereo(m)
    if len(si)==0 or any(e.specified==Chem.StereoSpecified.Unspecified for e in si): return None
    scf=MurckoScaffold.MurckoScaffoldSmiles(mol=m,includeChirality=False)
    nst=sum(1 for e in si if e.type==Chem.StereoType.Atom_Tetrahedral)
    return dict(psmi=Chem.MolToSmiles(m),gkey=strip_ms(pin),mflag=mflag(pin),scf=scf,nst=nst)

rows_out=[]
for tid,name in tids.items():
    q="""SELECT a.molregno,a.assay_id,a.doc_id,CAST(a.pchembl_value AS REAL)
         FROM activities a JOIN assays s ON a.assay_id=s.assay_id
         WHERE s.tid=? AND a.standard_type='IC50' AND a.standard_relation='=' AND a.standard_units='nM'
           AND a.pchembl_value IS NOT NULL AND a.data_validity_comment IS NULL
           AND (a.potential_duplicate=0 OR a.potential_duplicate IS NULL)
           AND s.assay_type='B' AND s.confidence_score=9"""
    recs=cur.execute(q,(tid,)).fetchall()
    # collapse duplicate (molregno,assay) -> median pchembl
    ma=collections.defaultdict(list); docof={}
    for mr,aid,doc,pch in recs: ma[(mr,aid)].append(pch); docof[(mr,aid)]=doc
    ma={k:st.median(v) for k,v in ma.items()}
    mols=set(mr for mr,_ in ma)
    info={}; ml=list(mols)
    for i in range(0,len(ml),900):
        for mr,smi,inchi in cur.execute("SELECT molregno,canonical_smiles,standard_inchi FROM compound_structures WHERE molregno IN (%s)"%",".join("?"*len(ml[i:i+900])),ml[i:i+900]):
            pi=parent_info(smi,inchi)
            if pi: info[mr]=pi
    # per assay, form enantiomer pairs; accumulate Delta-y per (gkey)
    perpair=collections.defaultdict(lambda: dict(dys=[],docs=set(),smiA=None,smiB=None,scf=None,nst=None))
    byassay=collections.defaultdict(list)
    for (mr,aid),pch in ma.items():
        if mr in info: byassay[aid].append((mr,pch,docof[(mr,aid)]))
    for aid,members in byassay.items():
        gk=collections.defaultdict(dict)  # gkey -> mflag -> (mr,pch,doc) best (max pch)
        for mr,pch,doc in members:
            g=info[mr]["gkey"]; f=info[mr]["mflag"]
            if f not in gk[g] or pch>gk[g][f][1]: gk[g][f]=(mr,pch,doc)
        for g,flags in gk.items():
            if len(flags)<2: continue
            fs=sorted(flags); (mrA,pA,dA),(mrB,pB,dB)=flags[fs[0]],flags[fs[1]]
            pp=perpair[(tid,g)]
            pp["dys"].append(pA-pB); pp["docs"].update([dA,dB])
            if pp["smiA"] is None:
                pp["smiA"],pp["smiB"]=info[mrA]["psmi"],info[mrB]["psmi"]
                pp["scf"],pp["nst"]=info[mrA]["scf"],info[mrA]["nst"]
    for (tid2,g),pp in perpair.items():
        rows_out.append(dict(tid=tid2,target=name,gkey=g,smiles_A=pp["smiA"],smiles_B=pp["smiB"],
                             dy=round(st.median(pp["dys"]),3),scaffold=pp["scf"],n_stereocenters=pp["nst"],
                             n_assays=len(pp["dys"]),n_docs=len(pp["docs"])))
    print(f"  {name[:34]:34} tid={tid} unique_pairs={sum(1 for r in rows_out if r['tid']==tid)}",flush=True)

fp=os.path.join(D,"cohort_pairs.csv")
cols=["tid","target","gkey","smiles_A","smiles_B","dy","scaffold","n_stereocenters","n_assays","n_docs"]
with open(fp,"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=cols); w.writeheader()
    for r in rows_out: w.writerow(r)
h=hashlib.sha256(open(fp,"rb").read()).hexdigest()
print(f"\nTOTAL cohort structural pairs: {len(rows_out)} across {len(tids)} targets")
print(f"cohort_pairs.csv sha256: {h}")
con.close()
