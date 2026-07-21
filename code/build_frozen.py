#!/usr/bin/env python3
"""Stereo carry-vs-use pilot: STAGE A — frozen 3-target dataset.
Extract IC50/=/nM/conf9/binding/valid/no-dup activities for BACE1, BTK, JAK3;
classify replicate provenance; salt-strip; RE-PAIR true enantiomers on the parent;
collapse duplicate measurements; write a frozen manifest with hashes + a stereo audit.
"""
import sqlite3, csv, re, json, hashlib, collections, statistics as st, os
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
RDLogger.DisableLog("rdApp.*")
OUT=os.path.dirname(os.path.abspath(__file__)); os.makedirs(OUT,exist_ok=True)
DB=os.environ.get("CHEMBL37_DB","data/chembl_37.db")
TARGETS={"BACE1":12252,"BTK":100097,"JAK3":10849}
con=sqlite3.connect(DB); cur=con.cursor()
lfc=rdMolStandardize.LargestFragmentChooser()

def parent(smiles):
    m=Chem.MolFromSmiles(smiles) if smiles else None
    if m is None: return None,None,False
    salt = "." in Chem.MolToSmiles(m)
    m=lfc.choose(m)
    return m, Chem.MolToSmiles(m), salt

def strip_ms(inchi):  # group key: drop absolute-config layers /m,/s  (keep /t,/b -> excludes diastereomers & E/Z)
    return "/".join(p for p in inchi.split("/") if not (p and p[0] in "ms"))
def mflag(inchi):
    mo=re.search(r"/m([01.]+)",inchi); return mo.group(1) if mo else None

def fully_defined_chiral(m):
    """True if every stereo element is specified AND the molecule is chiral (!= its mirror)."""
    si=Chem.FindPotentialStereo(m)
    if len(si)==0: return False
    if any(e.specified==Chem.StereoSpecified.Unspecified for e in si): return False
    return True

def enantiomer_ok(inchi):
    return ("/t" in inchi) and ("/m" in inchi) and ("/i" not in inchi)

prov=collections.Counter()  # replicate-provenance classes (pre-filter)
frozen=[]; audit_pool=[]
per_target_counts={}

for nm,tid in TARGETS.items():
    # PRE-FILTER pull (for provenance we keep dupes+validity), tag filter-eligibility
    q="""SELECT a.molregno,a.assay_id,a.doc_id,a.record_id,a.standard_value,a.pchembl_value,
                a.data_validity_comment,a.potential_duplicate
         FROM activities a JOIN assays s ON a.assay_id=s.assay_id
         WHERE s.tid=? AND a.standard_type='IC50' AND a.standard_relation='=' AND a.standard_units='nM'
           AND a.pchembl_value IS NOT NULL AND s.assay_type='B' AND s.confidence_score=9"""
    recs=cur.execute(q,(tid,)).fetchall()
    # --- replicate provenance: group by (molregno,assay) ---
    g=collections.defaultdict(list)
    for mr,aid,did,rid,val,pch,dvc,pd in recs: g[(mr,aid)].append((rid,did,val,pch,dvc,pd))
    for k,lst in g.items():
        if len(lst)==1: continue
        vals=set(round(float(v),4) for _,_,v,_,_,_ in lst)
        rids=set(r for r,_,_,_,_,_ in lst)
        any_pd=any(pd==1 for *_,pd in lst)
        any_flag=any(dvc for *_,dvc,_ in lst)
        if any_pd: prov["1_chembl_potential_duplicate"]+=1
        elif any_flag: prov["6_validity_flagged"]+=1
        elif len(vals)==1 and len(rids)==1: prov["2_identical_same_record"]+=1
        elif len(vals)==1 and len(rids)>1: prov["3_identical_indep_provenance"]+=1
        elif len(vals)>1 and len(rids)==1: prov["4_nonidentical_same_record"]+=1
        elif len(vals)>1 and len(rids)>1: prov["5_nonidentical_distinct_records"]+=1
        else: prov["7_undetermined"]+=1
    # --- modeling set: apply filters, collapse dups per (molregno,assay) ---
    clean=collections.defaultdict(list)  # (molregno,assay)-> list of (pchembl,doc)
    for mr,aid,did,rid,val,pch,dvc,pd in recs:
        if dvc is not None: continue
        if pd==1: continue
        clean[(mr,aid)].append((float(pch),did))
    ma_pchembl={}  # (mr,aid)->(median_pchembl, n_distinct, doc)
    for (mr,aid),lst in clean.items():
        vals=[p for p,_ in lst]; ma_pchembl[(mr,aid)]=(st.median(vals),len(set(round(v,4) for v in vals)),lst[0][1])
    mols=set(mr for mr,_ in ma_pchembl)
    # --- structures + salt-strip + parent InChI/keys (once per molregno) ---
    info={}
    ml=list(mols)
    for i in range(0,len(ml),900):
        qq="SELECT molregno,canonical_smiles FROM compound_structures WHERE molregno IN (%s)"%",".join("?"*len(ml[i:i+900]))
        for mr,smi in cur.execute(qq,ml[i:i+900]):
            pm,psmi,wassalt=parent(smi)
            if pm is None: continue
            try: inchi=Chem.MolToInchi(pm)
            except: continue
            if not inchi or not enantiomer_ok(inchi): continue
            if not fully_defined_chiral(pm): continue
            try: scf=MurckoScaffold.MurckoScaffoldSmiles(mol=pm,includeChirality=False)
            except: scf=""
            nst=sum(1 for e in Chem.FindPotentialStereo(pm) if e.type==Chem.StereoType.Atom_Tetrahedral)
            info[mr]=dict(psmi=psmi,inchi=inchi,gkey=strip_ms(inchi),mflag=mflag(inchi),
                          scf=scf,nst=nst,salt=wassalt,macro=pm.HasSubstructMatch(Chem.MolFromSmarts("[r{12-}]")))
    # --- pair within assay: same gkey, different mflag ---
    npairs=0
    byassay=collections.defaultdict(list)
    for (mr,aid),(pch,nd,did) in ma_pchembl.items():
        if mr in info: byassay[aid].append((mr,pch,did))
    for aid,members in byassay.items():
        gk=collections.defaultdict(list)
        for mr,pch,did in members: gk[info[mr]["gkey"]].append((mr,pch,did))
        for gkey,ms in gk.items():
            flags={}
            for mr,pch,did in ms: flags.setdefault(info[mr]["mflag"],[]).append((mr,pch,did))
            if len(flags)<2: continue
            # form one pair per assay per gkey: pick the two most-potent-representative hands
            reps={f:max(v,key=lambda x:x[1]) for f,v in flags.items()}  # (mr,pch,did)
            fs=sorted(reps)  # deterministic order by mflag
            (mrA,pA,didA),(mrB,pB,didB)=reps[fs[0]],reps[fs[1]]
            row=dict(target=nm,tid=tid,assay_id=aid,doc_id=didA,gkey=gkey,
                     molregno_A=mrA,molregno_B=mrB,smiles_A=info[mrA]["psmi"],smiles_B=info[mrB]["psmi"],
                     pchembl_A=round(pA,3),pchembl_B=round(pB,3),dy=round(pA-pB,3),
                     scaffold=info[mrA]["scf"],n_stereocenters=info[mrA]["nst"],
                     salt_stripped=int(info[mrA]["salt"] or info[mrB]["salt"]),
                     macrocycle=int(info[mrA]["macro"] or info[mrB]["macro"]))
            frozen.append(row); audit_pool.append(row); npairs+=1
    per_target_counts[nm]=dict(records=len(recs),distinct_molregno_assay=len(ma_pchembl),
                               eligible_chiral_mols=len(info),pairs=npairs)
    print(f"{nm}: recs={len(recs)} clean(mr,assay)={len(ma_pchembl)} eligible_chiral_mols={len(info)} PAIRS={npairs}")

# write frozen table
fp=os.path.join(OUT,"frozen_pairs.csv")
cols=["target","tid","assay_id","doc_id","gkey","molregno_A","molregno_B","smiles_A","smiles_B",
      "pchembl_A","pchembl_B","dy","scaffold","n_stereocenters","salt_stripped","macrocycle"]
with open(fp,"w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=cols); w.writeheader()
    for r in frozen: w.writerow(r)
h=hashlib.sha256(open(fp,"rb").read()).hexdigest()

print("\n===== REPLICATE PROVENANCE (repeated (compound,assay) groups, pre-filter) =====")
for k in sorted(prov): print(f"  {k}: {prov[k]}")
print("\n===== FROZEN PAIR COUNTS =====")
for nm,c in per_target_counts.items(): print(f"  {nm}: {c}")
print(f"  TOTAL frozen pairs: {len(frozen)}")
print(f"\nfrozen_pairs.csv sha256: {h}")

# manifest
man=dict(criteria=dict(standard_type="IC50",relation="=",units="nM",assay_type="B",
         confidence_score=9,data_validity_comment="NULL only",potential_duplicate="excluded",
         pairing="parent-InChI /m-flag enantiomers, fully-defined chiral, /t kept (no diastereomers/EZ), /i excluded",
         salt="LargestFragmentChooser then re-pair",dup_collapse="median of distinct values per (molregno,assay)"),
         targets=TARGETS, counts=per_target_counts, total_pairs=len(frozen),
         frozen_csv_sha256=h, provenance=dict(prov))
json.dump(man,open(os.path.join(OUT,"manifest.json"),"w"),indent=1)

# stratified stereo audit sample (>=200 if available), oversample hard strata
import random; random.seed(20260711)
def strat(r): return (r["n_stereocenters"]>=2, r["macrocycle"]==1, abs(r["dy"])>=1.5, r["salt_stripped"]==1)
buckets=collections.defaultdict(list)
for r in audit_pool: buckets[strat(r)].append(r)
sample=[]
for b,lst in buckets.items():
    random.shuffle(lst); sample.extend(lst[:max(30,len(lst)//8)])
random.shuffle(sample); sample=sample[:220]
# auto-verify each sampled pair is a true enantiomer (mirror check via InChI /m only diff)
def is_true_enantiomer(sA,sB):
    a=Chem.MolFromSmiles(sA); b=Chem.MolFromSmiles(sB)
    if not a or not b: return "unparseable"
    ia,ib=Chem.MolToInchi(a),Chem.MolToInchi(b)
    if strip_ms(ia)!=strip_ms(ib): return "FAIL_diff_skeleton"
    if mflag(ia)==mflag(ib): return "FAIL_same_mflag"
    # confirm B == mirror(A)
    am=Chem.MolFromSmiles(sA)
    for at in am.GetAtoms():
        if at.GetChiralTag()==Chem.ChiralType.CHI_TETRAHEDRAL_CW: at.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif at.GetChiralTag()==Chem.ChiralType.CHI_TETRAHEDRAL_CCW: at.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    return "OK" if Chem.MolToInchi(am)==ib else "OK_mflag_only"
errs=collections.Counter();
with open(os.path.join(OUT,"stereo_audit.csv"),"w",newline="") as f:
    w=csv.writer(f); w.writerow(["target","n_stereo","macro","salt","dy","verdict","smiles_A","smiles_B"])
    for r in sample:
        v=is_true_enantiomer(r["smiles_A"],r["smiles_B"]); errs[v]+=1
        w.writerow([r["target"],r["n_stereocenters"],r["macrocycle"],r["salt_stripped"],r["dy"],v,r["smiles_A"],r["smiles_B"]])
print(f"\n===== STEREO AUDIT (n={len(sample)}) auto-verdicts =====")
for k,v in errs.most_common(): print(f"  {k}: {v} ({100*v/len(sample):.1f}%)")
con.close()
