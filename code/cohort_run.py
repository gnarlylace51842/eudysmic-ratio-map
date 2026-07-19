#!/usr/bin/env python3
"""Confirmatory 14-target run per PROTOCOL_FROZEN.json. Resumable (per-target checkpoints).
Primary stat S_delta=1-MAE(dyhat,dy)/MAE(0,dy); models ridge + symmetrized RF (frozen hp);
r3/4096 chiral ECFP; orientation-symmetrized; scaffold-grouped 5-fold; PRIMARY null =
scaffold-CLUSTER sign-flip x2000; SECONDARY = pair permutation x1000 (robustness).
Then BH-FDR across targets on primary S_delta p, classify, and test preregistered R_t."""
import csv, os, json, numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors as rd, DataStructs
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from scipy.stats import spearmanr
RDLogger.DisableLog("rdApp.*")
D=os.path.dirname(os.path.abspath(__file__)); CK=os.path.join(D,"cohort_ckpt"); os.makedirs(CK,exist_ok=True)
NULL_PRIMARY=2000; NULL_SECOND=1000; SIGMA=0.45
rows=list(csv.DictReader(open(os.path.join(D,"cohort_pairs.csv"))))
fpc={}
def fp(s):
    if s in fpc: return fpc[s]
    m=Chem.MolFromSmiles(s); v=rd.GetMorganFingerprintAsBitVect(m,3,4096,useChirality=True)
    a=np.zeros(4096,np.int8); DataStructs.ConvertToNumpyArray(v,a); fpc[s]=a; return a
def models():
    return {"rf":RandomForestRegressor(n_estimators=80,min_samples_leaf=3,n_jobs=-1,random_state=0),
            "ridge":Ridge(alpha=10.0)}
def oof(delta,y,folds,which):
    p=np.zeros(len(y))
    for tr,te in folds:
        X=np.vstack([delta[tr],-delta[tr]]); yy=np.concatenate([y[tr],-y[tr]])
        h=models()[which]; h.fit(X,yy); p[te]=(h.predict(delta[te])-h.predict(-delta[te]))/2
    return p
def S_delta(p,y):
    b=np.mean(np.abs(y));  return 1-np.mean(np.abs(p-y))/b if b>0 else 0.0
def bal_acc(p,y):
    nz=np.abs(y)>1e-9
    pos=(y>1e-9); neg=(y<-1e-9)
    rp=np.mean(np.sign(p[pos])>0) if pos.any() else np.nan
    rn=np.mean(np.sign(p[neg])<0) if neg.any() else np.nan
    return np.nanmean([rp,rn])
def sp(p,y): return spearmanr(p,y).correlation if np.std(p)>0 else 0.0

tids=[]
for r in rows:
    if r["tid"] not in tids: tids.append(r["tid"])
for tid in tids:
    ck=os.path.join(CK,f"{tid}.json")
    if os.path.exists(ck): print(f"skip {tid} (done)",flush=True); continue
    tr=[r for r in rows if r["tid"]==tid]; name=tr[0]["target"]
    delta0=np.array([fp(r["smiles_A"])-fp(r["smiles_B"]) for r in tr],np.float32)
    dy0=np.array([float(r["dy"]) for r in tr],np.float32)
    scf=np.array([r["scaffold"] or r["gkey"] for r in tr])
    rs=np.random.default_rng(7); flip=rs.integers(0,2,len(dy0))*2-1
    delta=delta0*flip[:,None]; dy=dy0*flip
    folds=list(GroupKFold(5).split(delta,dy,scf))
    coll=float(np.mean(np.all(delta0==0,axis=1))); tiefrac=float(np.mean(np.abs(dy)<=1e-9))
    res={"tid":tid,"name":name,"n":len(tr),"n_scaffold":int(len(set(scf))),
         "median_absdy":float(np.median(np.abs(dy0))),"collision_r3":coll,"tie_frac":tiefrac,
         "Rt":round(float(np.median(np.abs(dy0))/SIGMA),3)}
    clusters=np.unique(scf); cidx={c:np.where(scf==c)[0] for c in clusters}
    rr=np.random.default_rng(202)
    for which in ("rf","ridge"):
        real={"S":S_delta(oof(delta,dy,folds,which),dy)}
        p_real=oof(delta,dy,folds,which)
        real={"S":S_delta(p_real,dy),"spear":sp(p_real,dy),"bal":bal_acc(p_real,dy)}
        # PRIMARY null: cluster sign-flip
        nullS=[]
        for _ in range(NULL_PRIMARY):
            s=np.ones(len(dy))
            for c in clusters: s[cidx[c]]=rr.integers(0,2)*2-1
            yn=dy*s
            nullS.append(S_delta(oof(delta,yn,folds,which),dy))
        nullS=np.array(nullS)
        p_emp=float((1+np.sum(nullS>=real["S"]))/(1+len(nullS)))
        # SECONDARY null: pair permutation (robustness)
        nptperm=[]
        for _ in range(NULL_SECOND):
            yn=dy[rr.permutation(len(dy))]
            nptperm.append(S_delta(oof(delta,yn,folds,which),dy))
        nptperm=np.array(nptperm)
        res[which]=dict(real=real,
            S_null_mean=float(nullS.mean()),S_null_q=[float(np.percentile(nullS,q)) for q in (2.5,50,97.5)],
            S_p_primary=p_emp, S_real_minus_null=float(real["S"]-nullS.mean()),
            S_p_secondary=float((1+np.sum(nptperm>=real["S"]))/(1+len(nptperm))))
        print(f"{tid} {name[:22]:22} {which:5} S={real['S']:.3f} p={p_emp:.4f} spear={real['spear']:.2f}",flush=True)
    json.dump(res,open(ck,"w"),indent=1)

# ---- aggregate: BH-FDR (primary=RF S_delta) + classification + R_t test ----
allr=[json.load(open(os.path.join(CK,f"{t}.json"))) for t in tids if os.path.exists(os.path.join(CK,f"{t}.json"))]
if len(allr)==len(tids):
    ps=[(r["tid"],r["rf"]["S_p_primary"]) for r in allr]
    order=sorted(range(len(ps)),key=lambda i:ps[i][1]); m=len(ps); q={}
    for rank,i in enumerate(order,1): q[ps[i][0]]=min(1.0,ps[i][1]*m/rank)
    # enforce monotonic BH
    for r in allr:
        rr_=r["rf"]
        qv=q[r["tid"]]
        if rr_["S_p_primary"]<=0.05 and qv<=0.10 and rr_["real"]["S"]>0 and abs(rr_["S_null_mean"])<0.05: cat="positive"
        elif abs(rr_["S_null_mean"])>=0.08: cat="failed_calibration"
        elif rr_["real"]["S"]>0 and rr_["S_p_primary"]<=0.20: cat="inconclusive"
        else: cat="null"
        r["category"]=cat; r["BH_q"]=round(qv,4)
    Rt=np.array([r["Rt"] for r in allr]); S=np.array([r["rf"]["real"]["S"] for r in allr])
    SP=np.array([r["rf"]["real"]["spear"] for r in allr]); BA=np.array([r["rf"]["real"]["bal"] for r in allr])
    agg=dict(targets=allr,
        Rt_vs_S=dict(rho=float(spearmanr(Rt,S).correlation),p=float(spearmanr(Rt,S).pvalue)),
        Rt_vs_spear=dict(rho=float(spearmanr(Rt,SP).correlation),p=float(spearmanr(Rt,SP).pvalue)),
        Rt_vs_bal=dict(rho=float(spearmanr(Rt,BA).correlation),p=float(spearmanr(Rt,BA).pvalue)),
        sens=dict(S_vs_npairs=float(spearmanr([r["n"] for r in allr],S).correlation),
                  S_vs_nscaf=float(spearmanr([r["n_scaffold"] for r in allr],S).correlation),
                  S_vs_tiefrac=float(spearmanr([r["tie_frac"] for r in allr],S).correlation)))
    json.dump(agg,open(os.path.join(D,"cohort_results.json"),"w"),indent=1)
    print("\n===== COHORT AGGREGATE =====")
    for r in sorted(allr,key=lambda x:-x["rf"]["real"]["S"]):
        print(f"  {r['name'][:26]:26} S={r['rf']['real']['S']:+.3f} p={r['rf']['S_p_primary']:.4f} q={r['BH_q']:.3f} "
              f"spear={r['rf']['real']['spear']:+.2f} Rt={r['Rt']:.2f} -> {r['category']}")
    print(f"\n  preregistered R_t vs S_delta: rho={agg['Rt_vs_S']['rho']:.2f} p={agg['Rt_vs_S']['p']:.3f}")
    print(f"  R_t vs Spearman: rho={agg['Rt_vs_spear']['rho']:.2f} ; R_t vs balanced-acc: rho={agg['Rt_vs_bal']['rho']:.2f}")
    print(f"  sensitivity S~npairs rho={agg['sens']['S_vs_npairs']:.2f} S~nscaf rho={agg['sens']['S_vs_nscaf']:.2f} S~tiefrac rho={agg['sens']['S_vs_tiefrac']:.2f}")
print("DONE")
