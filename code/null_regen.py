#!/usr/bin/env python3
"""Regenerate sign-flip null S_delta arrays (RF primary) for all 14 targets for the SUPPLEMENT
null-distribution panels. Faithful regeneration of the frozen null PROCEDURE; validated against
the stored 2000-rep summary (mean/quantiles) in cohort_results.json."""
import csv, os, json, numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors as rd, DataStructs
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold
RDLogger.DisableLog("rdApp.*")
D=os.path.dirname(os.path.abspath(__file__)); NREP=300
ABBR={12252:"BACE1",11307:"HDAC6",10919:"JAK1",10938:"JAK2",11638:"ERK2",100097:"BTK",103982:"ROR-γ",9:"EGFR",12694:"TYK2",19639:"IRAK4",165:"hERG",11177:"PI3Kδ",10849:"JAK3",10906:"SYK"}
rows=list(csv.DictReader(open(os.path.join(D,"cohort_pairs.csv"))))
cres={r["tid"]:r for r in json.load(open(os.path.join(D,"cohort_results.json")))["targets"]}
fpc={}
def fp(s):
    if s in fpc: return fpc[s]
    m=Chem.MolFromSmiles(s); v=rd.GetMorganFingerprintAsBitVect(m,3,4096,useChirality=True)
    a=np.zeros(4096,np.int8); DataStructs.ConvertToNumpyArray(v,a); fpc[s]=a; return a
def oof(d,y,folds):
    p=np.zeros(len(y))
    for tr,te in folds:
        X=np.vstack([d[tr],-d[tr]]); yy=np.concatenate([y[tr],-y[tr]])
        h=RandomForestRegressor(n_estimators=80,min_samples_leaf=3,n_jobs=-1,random_state=0); h.fit(X,yy)
        p[te]=(h.predict(d[te])-h.predict(-d[te]))/2
    return p
def S(p,y): b=np.mean(np.abs(y)); return 1-np.mean(np.abs(p-y))/b if b>0 else 0
out={}
tids=[]
for r in rows:
    if r["tid"] not in tids: tids.append(r["tid"])
for tid in tids:
    tr=[r for r in rows if r["tid"]==tid]
    d0=np.array([fp(r["smiles_A"])-fp(r["smiles_B"]) for r in tr],np.float32)
    dy0=np.array([float(r["dy"]) for r in tr],np.float32); scf=np.array([r["scaffold"] or r["gkey"] for r in tr])
    rs=np.random.default_rng(7); fl=rs.integers(0,2,len(dy0))*2-1; d=d0*fl[:,None]; dy=dy0*fl
    folds=list(GroupKFold(5).split(d,dy,scf)); real=S(oof(d,dy,folds),dy)
    cl=np.unique(scf); idx={c:np.where(scf==c)[0] for c in cl}; rr=np.random.default_rng(303); ns=[]
    for _ in range(NREP):
        s=np.ones(len(dy))
        for c in cl: s[idx[c]]=rr.integers(0,2)*2-1
        ns.append(S(oof(d,dy*s,folds),dy))
    ns=np.array(ns); st=cres[tid]["rf"]
    out[tid]=dict(ab=ABBR[int(tid)],real_S=float(real),null=[round(float(x),4) for x in ns],
                  regen_mean=float(ns.mean()),stored_mean=st["S_null_mean"],
                  regen_q975=float(np.percentile(ns,97.5)),stored_q975=st["S_null_q"][2],stored_p=st["S_p_primary"])
    print(f"{ABBR[int(tid)]:7} real={real:+.3f} regen_null_mean={ns.mean():+.3f} (stored {st['S_null_mean']:+.3f}) regen_q975={np.percentile(ns,97.5):+.3f} (stored {st['S_null_q'][2]:+.3f})",flush=True)
json.dump(out,open(os.path.join(D,"null_arrays.json"),"w"))
print("WROTE null_arrays.json")
