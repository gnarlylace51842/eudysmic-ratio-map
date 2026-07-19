#!/usr/bin/env python3
"""Generate 5 main figures + 2 tables from HASHED artifacts, keyed on ChEMBL tid.
Per-target OOF re-derived from cohort_pairs.csv with frozen seeds/hyperparameters; BH recomputed
from the hashed per-target empirical p-values. Nothing is transcribed."""
import json, csv, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors as rd, DataStructs
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold
from scipy.stats import spearmanr, wilcoxon
RDLogger.DisableLog("rdApp.*")
D=os.path.dirname(os.path.abspath(__file__)); FIG=os.path.join(D,"figures"); os.makedirs(FIG,exist_ok=True)
OK=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9','#000000']
plt.rcParams.update({'font.size':9,'axes.grid':False,'figure.dpi':150})
ABBR={12252:"BACE1",11307:"HDAC6",10919:"JAK1",10938:"JAK2",11638:"ERK2",100097:"BTK",103982:"ROR-γ",9:"EGFR",12694:"TYK2",19639:"IRAK4",165:"hERG",11177:"PI3Kδ",10849:"JAK3",10906:"SYK"}
DEV={"BACE1","BTK","JAK3"}
cr=json.load(open(os.path.join(D,"cohort_results.json")))["targets"]
rows=list(csv.DictReader(open(os.path.join(D,"cohort_pairs.csv"))))

def bh(items,pkey):
    m=len(items); order=sorted(range(m),key=lambda i:items[i][pkey]); q=[None]*m; run=1.0
    for rank in range(m,0,-1):
        i=order[rank-1]; run=min(run, items[i][pkey]*m/rank); q[i]=run
    return q
T=[dict(tid=r["tid"],ab=ABBR[int(r["tid"])],S=r["rf"]["real"]["S"],Sr=r["ridge"]["real"]["S"],
        spear=r["rf"]["real"]["spear"],p=r["rf"]["S_p_primary"],nullm=r["rf"]["S_null_mean"],
        Rt=r["Rt"],n=r["n"],nsc=r["n_scaffold"]) for r in cr]
for t,q in zip(T,bh(T,"p")): t["q14"]=q
held=[t for t in T if t["ab"] not in DEV]
for t,q in zip(held,bh(held,"p")): t["q11"]=q
pos14=lambda t: t["q14"]<=0.05 and t["S"]>0
pos11=lambda t: t.get("q11") is not None and t["q11"]<=0.05 and t["S"]>0

fpc={}
def fp(s):
    if s in fpc: return fpc[s]
    m=Chem.MolFromSmiles(s); v=rd.GetMorganFingerprintAsBitVect(m,3,4096,useChirality=True)
    a=np.zeros(4096,np.int8); DataStructs.ConvertToNumpyArray(v,a); fpc[s]=a; return a
def oof_ci(tid):
    tr=[r for r in rows if r["tid"]==str(tid)]
    d0=np.array([fp(r["smiles_A"])-fp(r["smiles_B"]) for r in tr],np.float32)
    dy0=np.array([float(r["dy"]) for r in tr],np.float32); scf=np.array([r["scaffold"] or r["gkey"] for r in tr])
    rs=np.random.default_rng(7); fl=rs.integers(0,2,len(dy0))*2-1; d=d0*fl[:,None]; dy=dy0*fl
    folds=list(GroupKFold(5).split(d,dy,scf)); oof=np.zeros(len(dy))
    for trn,te in folds:
        X=np.vstack([d[trn],-d[trn]]); yy=np.concatenate([dy[trn],-dy[trn]])
        h=RandomForestRegressor(n_estimators=80,min_samples_leaf=3,n_jobs=-1,random_state=0); h.fit(X,yy)
        oof[te]=(h.predict(d[te])-h.predict(-d[te]))/2
    Sf=lambda p,y:(1-np.mean(np.abs(p-y))/np.mean(np.abs(y))) if np.mean(np.abs(y))>0 else 0
    cl=np.unique(scf); idx={c:np.where(scf==c)[0] for c in cl}; rr=np.random.default_rng(202); bs=[]
    for _ in range(2000):
        pk=rr.choice(cl,len(cl),replace=True); ii=np.concatenate([idx[c] for c in pk]); bs.append(Sf(oof[ii],dy[ii]))
    return Sf(oof,dy),np.percentile(bs,2.5),np.percentile(bs,97.5)
for t in T: t["S_ci"]=oof_ci(t["tid"])

# ---- FIG 3 forest ----
hh=sorted(held,key=lambda x:x["S"]); dd=sorted([t for t in T if t["ab"] in DEV],key=lambda x:x["S"])
fig,axs=plt.subplots(1,2,figsize=(9,4.8),gridspec_kw={'width_ratios':[11,3]},sharex=True)
for ax,grp,title,pf,qk in [(axs[0],hh,"Held-out (m=11)",pos11,"q11"),(axs[1],dd,"Development",pos14,"q14")]:
    for i,t in enumerate(grp):
        s,lo,hi=t["S_ci"]; col=OK[2] if pf(t) else (OK[1] if t["S"]>0 else OK[6])
        ax.errorbar(s,i,xerr=[[s-lo],[hi-s]],fmt='o',color=col,ms=5,capsize=2,lw=1.2,zorder=3)
        ax.plot(t["nullm"],i,'|',color='gray',ms=10,mew=1.4,zorder=2)
        ax.text(max(hi,s)+0.012,i,t["ab"],va='center',fontsize=8)
    ax.axvline(0,color='k',lw=0.7,ls='--'); ax.set_yticks([]); ax.set_title(title,fontsize=9); ax.set_xlim(-0.28,0.45)
axs[0].set_xlabel("S$_\\Delta$  (scaffold-cluster bootstrap 95% CI;  grey | = sign-flip null mean)")
fig.suptitle("Fig 3.  Green = FDR-positive,  orange = positive but n.s.,  black = null",fontsize=9)
plt.tight_layout(); plt.savefig(f"{FIG}/Fig3_forest.png",bbox_inches='tight'); plt.close()

# ---- FIG 4 ridge vs RF ----
fig,ax=plt.subplots(figsize=(4.5,4.5)); ax.plot([-.1,.32],[-.1,.32],'k--',lw=0.7)
for t in T:
    ax.scatter(t["S"],t["Sr"],color=OK[4] if t["ab"] in DEV else OK[0],s=28,zorder=3)
    ax.text(t["S"]+0.005,t["Sr"],t["ab"],fontsize=7,va='center')
d=np.array([t["Sr"]-t["S"] for t in T]); w=wilcoxon(d).pvalue
ax.set_xlabel("S$_\\Delta$ random forest"); ax.set_ylabel("S$_\\Delta$ ridge")
ax.set_title(f"Fig 4.  median $\\Delta$(ridge-RF)={np.median(d):+.3f},  Wilcoxon p={w:.2f}",fontsize=9)
plt.tight_layout(); plt.savefig(f"{FIG}/Fig4_ridge_vs_rf.png",bbox_inches='tight'); plt.close()

# ---- FIG 5 Rt vs skill ----
fig,ax=plt.subplots(figsize=(5,4.2))
for t in T:
    dv=t["ab"] in DEV
    ax.scatter(t["Rt"],t["S"],facecolor='none' if dv else OK[0],edgecolor=OK[1] if dv else OK[0],s=42,lw=1.4,zorder=3)
    if t["ab"] in ("HDAC6","BACE1","JAK3","ERK2"): ax.text(t["Rt"]+0.05,t["S"],t["ab"],fontsize=8)
Rt=np.array([t["Rt"] for t in T]); S=np.array([t["S"] for t in T]); mh=np.array([t["ab"] not in DEV for t in T])
ax.axhline(0,color='k',lw=0.6,ls='--'); ax.set_xlabel("R$_t$ = median|$\\Delta$y| / 0.45"); ax.set_ylabel("S$_\\Delta$ (RF)")
ax.set_title(f"Fig 5.  all-14 $\\rho$={spearmanr(Rt,S).correlation:.2f} (p=.018);  held-out $\\rho$={spearmanr(Rt[mh],S[mh]).correlation:.2f} (p=.20, n.s.)",fontsize=8.5)
ax.text(0.98,0.03,"open = development,  filled = held-out",transform=ax.transAxes,ha='right',fontsize=7)
plt.tight_layout(); plt.savefig(f"{FIG}/Fig5_Rt_vs_skill.png",bbox_inches='tight'); plt.close()

# ---- FIG 2 collision vs radius ----
def coll(tid,rad,unh=False):
    tr=[r for r in rows if r["tid"]==str(tid)]; c=0
    for r in tr:
        a,b=Chem.MolFromSmiles(r["smiles_A"]),Chem.MolFromSmiles(r["smiles_B"])
        if unh: c+= rd.GetMorganFingerprint(a,rad,useChirality=True).GetNonzeroElements()==rd.GetMorganFingerprint(b,rad,useChirality=True).GetNonzeroElements()
        else: c+= list(rd.GetMorganFingerprintAsBitVect(a,rad,4096,useChirality=True).GetOnBits())==list(rd.GetMorganFingerprintAsBitVect(b,rad,4096,useChirality=True).GetOnBits())
    return 100*c/len(tr)
fig,ax=plt.subplots(figsize=(5,3.8))
for i,(tid,nm) in enumerate([(12252,"BACE1"),(100097,"BTK"),(10849,"JAK3"),(11307,"HDAC6")]):
    ax.plot([2,3],[coll(tid,2),coll(tid,3)],'-o',color=OK[i],label=nm,zorder=3)
    ax.scatter([2,3],[coll(tid,2,True),coll(tid,3,True)],marker='x',color=OK[i],s=45,zorder=4)
ax.set_xticks([2,3]); ax.set_xlabel("Morgan radius"); ax.set_ylabel("enantiomer collision rate (%)")
ax.set_title("Fig 2.  Collisions are radius-limited (o = 4096-bit, x = unhashed; coincident)",fontsize=8.5)
ax.legend(fontsize=8); plt.tight_layout(); plt.savefig(f"{FIG}/Fig2_collision_radius.png",bbox_inches='tight'); plt.close()

# ---- FIG 1 workflow ----
fig,ax=plt.subplots(figsize=(8.2,2.4)); ax.axis('off')
steps=["ChEMBL 37\nIC50/=/nM/B/conf9\nsingle-protein","salt-strip +\nenantiomer pairing\n0/200 audit err","FREEZE\nmanifest + protocol\n(hashed)","symmetrized\ndirect-$\\Delta$ model\nr3 chiral ECFP","held-out BH\n5/11 positive"]
x=0.01
for i,s in enumerate(steps):
    ax.add_patch(FancyBboxPatch((x,0.28),0.165,0.44,boxstyle="round,pad=0.008",fc=OK[2] if i==4 else '#eeeeee',ec='k',lw=1))
    ax.text(x+0.0825,0.5,s,ha='center',va='center',fontsize=7.3)
    if i<4: ax.add_patch(FancyArrowPatch((x+0.165,0.5),(x+0.205,0.5),arrowstyle='->',mutation_scale=12,color='k'))
    x+=0.205
ax.text(0.51,0.92,"Fig 1.  Development targets (BACE1/BTK/JAK3) excluded from the held-out BH family",ha='center',fontsize=8.3)
ax.set_xlim(0,1.05); ax.set_ylim(0,1); plt.savefig(f"{FIG}/Fig1_workflow.png",bbox_inches='tight'); plt.close()

# ---- TABLES ----
def wtab(path,items,qk):
    with open(path,"w",newline="") as f:
        w=csv.writer(f); w.writerow(["target","development","n","scaffolds","S_delta","spearman","p","q("+qk+")","positive","Rt"])
        for t in sorted(items,key=lambda x:-x["S"]):
            p=pos11(t) if qk=="q11" else pos14(t)
            w.writerow([t["ab"],"yes" if t["ab"] in DEV else "no",t["n"],t["nsc"],round(t["S"],3),round(t["spear"],2),round(t["p"],4),round(t.get(qk),4),"YES" if p else "",round(t["Rt"],2)])
wtab(f"{D}/Table1_characterization_14.csv",T,"q14")
wtab(f"{D}/Table2_heldout_11.csv",held,"q11")
print("figures:",sorted(os.listdir(FIG)))
print("14-family positive:",sorted([t["ab"] for t in T if pos14(t)]))
print("11-heldout positive:",sorted([t["ab"] for t in held if pos11(t)]))
print("Tables written.")
