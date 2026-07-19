#!/usr/bin/env python3
"""Generate supplement tables + figures from frozen artifacts (no transcription)."""
import csv, os, json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors as rd, DataStructs
from scipy.stats import spearmanr
RDLogger.DisableLog("rdApp.*")
D=os.path.dirname(os.path.abspath(__file__)); FIG=os.path.join(D,"figures"); OK=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9','#000000']
plt.rcParams.update({'font.size':9,'figure.dpi':150})
ABBR={12252:"BACE1",11307:"HDAC6",10919:"JAK1",10938:"JAK2",11638:"ERK2",100097:"BTK",103982:"ROR-γ",9:"EGFR",12694:"TYK2",19639:"IRAK4",165:"hERG",11177:"PI3Kδ",10849:"JAK3",10906:"SYK"}
DEV={"BACE1","BTK","JAK3"}
rows=list(csv.DictReader(open(f"{D}/cohort_pairs.csv")))
cr=json.load(open(f"{D}/cohort_results.json"))["targets"]
def bh(items,pk):
    m=len(items); o=sorted(range(m),key=lambda i:items[i][pk]); q=[None]*m; run=1.0
    for r in range(m,0,-1): i=o[r-1]; run=min(run,items[i][pk]*m/r); q[i]=run
    return q
T=[dict(tid=r["tid"],ab=ABBR[int(r["tid"])],S=r["rf"]["real"]["S"],Sr=r["ridge"]["real"]["S"],spear=r["rf"]["real"]["spear"],
        p=r["rf"]["S_p_primary"],psec=r["rf"]["S_p_secondary"],nullm=r["rf"]["S_null_mean"],Rt=r["Rt"],
        n=r["n"],nsc=r["n_scaffold"],tie=r["tie_frac"],coll=r["collision_r3"]) for r in cr]
for t,q in zip(T,bh(T,"p")): t["q14"]=q
held=[t for t in T if t["ab"] not in DEV]
for t,q in zip(held,bh(held,"p")): t["q11"]=q

# ---- STable_funnel ----
with open(f"{D}/STable_funnel.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["stage","count"])
    for s,c in [("tightened IC50 activities (single-protein)",827880),("distinct molregnos",511596),
                ("eligible fully-defined chiral molregnos",134564),("structural enantiomer pairs (14 targets)",1967),
                ("qualifying targets",14)]: w.writerow([s,c])

# ---- STable_provenance ----
with open(f"{D}/STable_provenance.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["repeat-group class","count","pct"])
    for s,c,p in [("total repeated (compound,assay) groups",17676,100.0),("non-identical, same source record (NOT independent)",12768,72.2),
                  ("provenance-distinct nonidentical (distinct records+values)",2688,15.2),("identical, same record",1917,10.8),
                  ("identical, distinct records",303,1.7)]: w.writerow([s,c,p])

# ---- collision audit (all 14, config invariance) ----
def bit(a,b,r,nb): return list(rd.GetMorganFingerprintAsBitVect(a,r,nb,useChirality=True).GetOnBits())==list(rd.GetMorganFingerprintAsBitVect(b,r,nb,useChirality=True).GetOnBits())
def cnt(a,b,r,nb): return rd.GetHashedMorganFingerprint(a,r,nBits=nb,useChirality=True).GetNonzeroElements()==rd.GetHashedMorganFingerprint(b,r,nb,useChirality=True).GetNonzeroElements()
def unh(a,b,r): return rd.GetMorganFingerprint(a,r,useChirality=True).GetNonzeroElements()==rd.GetMorganFingerprint(b,r,useChirality=True).GetNonzeroElements()
collrows=[]
for t in T:
    tr=[r for r in rows if r["tid"]==t["tid"]]; mols=[(Chem.MolFromSmiles(r["smiles_A"]),Chem.MolFromSmiles(r["smiles_B"])) for r in tr]
    def rate(fn): return 100*sum(fn(a,b) for a,b in mols)/len(mols)
    r2=[rate(lambda a,b:bit(a,b,2,nb)) for nb in (1024,2048,4096,8192)]+[rate(lambda a,b:cnt(a,b,2,4096)),rate(lambda a,b:unh(a,b,2))]
    r3=[rate(lambda a,b:bit(a,b,3,nb)) for nb in (1024,2048,4096,8192)]+[rate(lambda a,b:cnt(a,b,3,4096)),rate(lambda a,b:unh(a,b,3))]
    inv=(max(r2)-min(r2)<1e-9) and (max(r3)-min(r3)<1e-9)
    collrows.append((t["ab"],r2[2],r3[2],inv))
with open(f"{D}/STable_collision.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["target","coll_r2_pct","coll_r3_pct","config_invariant(bits/count/unhashed identical)"])
    for ab,a2,a3,inv in collrows: w.writerow([ab,round(a2,1),round(a3,1),inv])
print("collision config-invariant for all 14 targets:", all(inv for *_,inv in collrows))

# ---- R_t frozen vs exploratory ----
adyb={}
for r in rows: adyb.setdefault(r["tid"],[]).append(abs(float(r["dy"])))
Rtrows=[]
for t in T:
    a=np.array(adyb[t["tid"]]); Rt_med=np.median(a)/0.45; Rt_uq=np.percentile(a,75)/0.45; frac=np.mean(a>0.45)
    Rtrows.append((t["ab"],t["ab"] not in DEV,t["S"],Rt_med,Rt_uq,frac))
with open(f"{D}/STable_Rt.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["target","held_out","S_delta","Rt_frozen(median/0.45)","Rt_exploratory(UQ/0.45)","frac|dy|>0.45"])
    for ab,ho,S,rm,ru,fr in Rtrows: w.writerow([ab,ho,round(S,3),round(rm,2),round(ru,2),round(fr,2)])
S=np.array([x[2] for x in Rtrows]); rm=np.array([x[3] for x in Rtrows]); ru=np.array([x[4] for x in Rtrows]); ho=np.array([x[1] for x in Rtrows])
print(f"R_t frozen(median) vs S: all rho={spearmanr(rm,S).correlation:.2f}  held-out rho={spearmanr(rm[ho],S[ho]).correlation:.2f}")
print(f"R_t exploratory(UQ)  vs S: all rho={spearmanr(ru,S).correlation:.2f}  held-out rho={spearmanr(ru[ho],S[ho]).correlation:.2f}")

# ---- complete target table ----
with open(f"{D}/STable_target_full.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["target","held_out","n","scaffolds","tie_frac","coll_r3","S_rf","S_ridge","spearman","p_primary","p_secondary","q14","q11","Rt"])
    for t in sorted(T,key=lambda x:-x["S"]):
        w.writerow([t["ab"],t["ab"] not in DEV,t["n"],t["nsc"],round(t["tie"],3),round(t["coll"],3),round(t["S"],3),round(t["Sr"],3),round(t["spear"],2),round(t["p"],4),round(t["psec"],4),round(t["q14"],4),round(t.get("q11"),4) if t.get("q11") is not None else "",round(t["Rt"],2)])

# ---- SFig collision ----
fig,ax=plt.subplots(figsize=(6,3.6))
xs=np.arange(len(collrows))
ax.bar(xs-0.2,[a2 for _,a2,_,_ in collrows],0.4,label="radius 2",color=OK[1])
ax.bar(xs+0.2,[a3 for _,_,a3,_ in collrows],0.4,label="radius 3",color=OK[0])
ax.set_xticks(xs); ax.set_xticklabels([a for a,*_ in collrows],rotation=60,ha='right',fontsize=7); ax.set_ylabel("collision %")
ax.set_title("Fig S-collision. r2 vs r3 (config-invariant across bits/count/unhashed)",fontsize=8.5); ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig(f"{FIG}/SFig_collision.png",bbox_inches='tight'); plt.close()
# ---- SFig interp vs extrap ----
ie=json.load(open(f"{D}/S_interp_vs_extrap.json"))
fig,ax=plt.subplots(figsize=(4.3,4.3)); ax.plot([-.1,.32],[-.1,.32],'k--',lw=0.7)
for d in ie:
    ax.scatter(d["extrap"],d["interp"],color=OK[4] if d["ab"] in DEV else OK[0],s=26); ax.text(d["extrap"]+0.005,d["interp"],d["ab"],fontsize=6.5)
ax.set_xlabel("S$_\\Delta$ extrapolation (scaffold)"); ax.set_ylabel("S$_\\Delta$ interpolation (gkey)")
ax.set_title("Fig S-split. interpolation ≥ extrapolation (median gap +0.017)",fontsize=8.5)
plt.tight_layout(); plt.savefig(f"{FIG}/SFig_interp_extrap.png",bbox_inches='tight'); plt.close()
# ---- SFig Rt frozen vs exploratory ----
fig,axx=plt.subplots(1,2,figsize=(8,3.6),sharey=True)
for ax,xv,lab in [(axx[0],rm,"R_t frozen (median/0.45)"),(axx[1],ru,"R_t exploratory (UQ/0.45)")]:
    for i,(ab,h,Sv,_,_,_) in enumerate(Rtrows):
        ax.scatter(xv[i],Sv,facecolor='none' if not h else OK[0],edgecolor=OK[1] if not h else OK[0],s=34,lw=1.3)
    ax.axhline(0,color='k',lw=0.5,ls='--'); ax.set_xlabel(lab)
    ax.set_title(f"all $\\rho$={spearmanr(xv,S).correlation:.2f}, held-out $\\rho$={spearmanr(xv[ho],S[ho]).correlation:.2f}",fontsize=8)
axx[0].set_ylabel("S$_\\Delta$")
fig.suptitle("Fig S-Rt. frozen (median) vs exploratory (upper-quartile) signal statistic",fontsize=8.5)
plt.tight_layout(); plt.savefig(f"{FIG}/SFig_Rt.png",bbox_inches='tight'); plt.close()
print("wrote STable_funnel/provenance/collision/Rt/target_full.csv + SFig_collision/interp_extrap/Rt.png")
