#!/usr/bin/env python3
"""Refinements 1,3,4: orientation/sign-balance diagnostic; chiral-fingerprint collision audit
across configs (radius x bits x binary/count x unhashed); manual BACE1 collision classification;
extended >=200-pair stereo audit with Wilson CI."""
import csv, os, math, collections, numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors as rd, DataStructs
RDLogger.DisableLog("rdApp.*")
D=os.path.dirname(os.path.abspath(__file__))
rows=list(csv.DictReader(open(os.path.join(D,"frozen_pairs.csv"))))
def mol(s): return Chem.MolFromSmiles(s)

# ---------- 1. ORIENTATION / SIGN-BALANCE DIAGNOSTIC ----------
print("===== 1. ORIENTATION / SIGN DIAGNOSTIC =====")
print("canonical rule as-built: A = lower InChI /m flag (arbitrary, not deployment-available)")
for tgt in ["BACE1","BTK","JAK3"]:
    tr=[r for r in rows if r["target"]==tgt]; dy=np.array([float(r["dy"]) for r in tr])
    pos,neg,tie=np.mean(dy>1e-9),np.mean(dy<-1e-9),np.mean(np.abs(dy)<=1e-9)
    maj=max(pos,neg)
    # balanced directional accuracy of the MAJORITY predictor = 0.5 by construction; report data imbalance
    print(f"  {tgt}: +Δy={pos:.2%} -Δy={neg:.2%} ties(|Δy|<=1e-9)={tie:.2%}  majority-sign baseline(raw)={maj:.3f}")
print("  -> null rank ~0.45-0.47 (not 0.50) is expected: ties + the antisymmetric OOF predicting ~0 on")
print("     hard pairs make sign()=+ slightly less than half; orientation-symmetrization + balanced")
print("     accuracy (macro-avg of +Δy and -Δy recall) are the fixes reported for the cohort.")

# ---------- 3. CHIRAL-FINGERPRINT COLLISION AUDIT ACROSS CONFIGS ----------
print("\n===== 3. CHIRAL-FINGERPRINT COLLISION AUDIT (collision = identical vectors for the two enantiomers) =====")
def coll_rate(tr, radius, nbits, count):
    c=0
    for r in tr:
        a,b=mol(r["smiles_A"]),mol(r["smiles_B"])
        if count:
            fa=rd.GetHashedMorganFingerprint(a,radius,nBits=nbits,useChirality=True)
            fb=rd.GetHashedMorganFingerprint(b,radius,nBits=nbits,useChirality=True)
            same=(fa.GetNonzeroElements()==fb.GetNonzeroElements())
        else:
            fa=rd.GetMorganFingerprintAsBitVect(a,radius,nBits=nbits,useChirality=True)
            fb=rd.GetMorganFingerprintAsBitVect(b,radius,nBits=nbits,useChirality=True)
            same=(list(fa.GetOnBits())==list(fb.GetOnBits()))
        c+=same
    return c/len(tr)
def coll_unhashed(tr,radius):
    c=0
    for r in tr:
        a,b=mol(r["smiles_A"]),mol(r["smiles_B"])
        fa=rd.GetMorganFingerprint(a,radius,useChirality=True).GetNonzeroElements()
        fb=rd.GetMorganFingerprint(b,radius,useChirality=True).GetNonzeroElements()
        c+=(fa==fb)
    return c/len(tr)
print(f"{'target':6} {'r':>2} {'bits':>5} {'binary':>7} {'count':>7} | unhashed(r) collision")
for tgt in ["BACE1","BTK","JAK3"]:
    tr=[r for r in rows if r["target"]==tgt]
    for radius in (2,3):
        uh=coll_unhashed(tr,radius)
        for nbits in (1024,2048,4096,8192):
            cb=coll_rate(tr,radius,nbits,False); cc=coll_rate(tr,radius,nbits,True)
            print(f"{tgt:6} {radius:>2} {nbits:>5} {cb:>7.1%} {cc:>7.1%} | {uh:.1%}")
    print()

# ---------- 3b. MANUAL BACE1 COLLISION CLASSIFICATION (r2,2048,binary) ----------
print("===== 3b. BACE1 collisions at (r2,2048,binary): hashing vs encoding vs symmetry =====")
bt=[r for r in rows if r["target"]=="BACE1"]
cls=collections.Counter(); examples=[]
for r in bt:
    a,b=mol(r["smiles_A"]),mol(r["smiles_B"])
    fa=list(rd.GetMorganFingerprintAsBitVect(a,2,2048,useChirality=True).GetOnBits())
    fb=list(rd.GetMorganFingerprintAsBitVect(b,2,2048,useChirality=True).GetOnBits())
    if fa!=fb: continue  # not a collision at this config
    # unhashed r3 with chirality: do the underlying identifiers differ?
    ua=rd.GetMorganFingerprint(a,3,useChirality=True).GetNonzeroElements()
    ub=rd.GetMorganFingerprint(b,3,useChirality=True).GetNonzeroElements()
    ncip=rd.CalcNumAtomStereoCenters(a)
    if ua!=ub: cls["hashing_or_radius (unhashed r3 DIFFERS -> config artifact)"]+=1
    elif ncip==0: cls["symmetry/no_CIP_center (achiral by CIP)"]+=1
    else: cls["genuine_encoding_limit (unhashed r3 identical despite CIP center)"]+=1
    if len(examples)<3: examples.append((ncip,r["smiles_A"][:60]))
print("  BACE1 collision pairs classified:", dict(cls))
for n,s in examples: print(f"    e.g. nCIP={n}  {s}")

# ---------- 4. EXTENDED STEREO AUDIT (>=200) ----------
print("\n===== 4. STEREO AUDIT (extended, oversample hard strata) =====")
def strip_ms(i): return "/".join(p for p in i.split("/") if not (p and p[0] in "ms"))
import re
def mflag(i):
    m=re.search(r"/m([01.]+)",i); return m.group(1) if m else None
def verdict(sA,sB):
    a,b=mol(sA),mol(sB)
    if not a or not b: return "unparseable"
    ia,ib=Chem.MolToInchi(a),Chem.MolToInchi(b)
    if strip_ms(ia)!=strip_ms(ib): return "FAIL_diff_skeleton"
    if mflag(ia)==mflag(ib): return "FAIL_same_config"
    am=mol(sA)
    for at in am.GetAtoms():
        t=at.GetChiralTag()
        if t==Chem.ChiralType.CHI_TETRAHEDRAL_CW: at.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif t==Chem.ChiralType.CHI_TETRAHEDRAL_CCW: at.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    return "OK_true_enantiomer" if Chem.MolToInchi(am)==ib else "OK_mflag_only"
# oversample: all BACE1, all collisions, multi-stereo, macro, salt, extreme/near-zero dy, charged
def charged(r):
    m=mol(r["smiles_A"]); return m is not None and any(a.GetFormalCharge()!=0 for a in m.GetAtoms())
rng=np.random.default_rng(5)
pool=[]
for r in rows:
    hard = (r["target"]=="BACE1") or int(r["n_stereocenters"])>=2 or int(r["macrocycle"])==1 \
           or int(r["salt_stripped"])==1 or abs(float(r["dy"]))>=2.5 or abs(float(r["dy"]))<=0.05 or charged(r)
    if hard: pool.append(r)
idx=rng.permutation(len(pool))[:max(200,0)]
sample=[pool[i] for i in idx][:230]
res=collections.Counter()
for r in sample: res[verdict(r["smiles_A"],r["smiles_B"])]+=1
n=len(sample); ok=res["OK_true_enantiomer"]+res["OK_mflag_only"]; err=n-ok
p=err/n; z=1.96
wl=(p+z*z/(2*n)-z*math.sqrt((p*(1-p)+z*z/(4*n))/n))/(1+z*z/n)
wu=(p+z*z/(2*n)+z*math.sqrt((p*(1-p)+z*z/(4*n))/n))/(1+z*z/n)
print(f"  audited n={n} (oversampled hard strata)  verdicts: {dict(res)}")
print(f"  pairing ERROR rate = {err}/{n} = {p:.3%}  (95% Wilson CI [{wl:.3%}, {wu:.3%}])")
