#!/usr/bin/env python3
"""Per-type marker 1:1-ortholog retention vs residual magnitude.
Per type: fraction of marker genes that are 1:1 orthologs (orthology_type==
ortholog_one2one) -> Spearman vs per-type residual_magnitude. Primary marker set =
Table S5 matching_basis genes; secondary = CellMarker. Outputs to analysis/sensitivity_analyses/.
"""
import json, re
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
REPO = Path(__file__).resolve().parents[2]; OUT = REPO/"analysis/sensitivity_analyses"

orth = pd.read_csv(REPO/"data/phase1/orthologs_human_mouse.csv")
print("orthology_type values:", orth["orthology_type"].value_counts().to_dict())
one2one = set(orth.loc[orth["orthology_type"]=="ortholog_one2one","human_gene_name"].astype(str))
# fallback: if the file is already all-1:1, membership = retained
all_orth_sym = set(orth["human_gene_name"].astype(str))
print(f"1:1 ortholog symbols: {len(one2one)} | all-orth symbols: {len(all_orth_sym)}")

resid = pd.read_csv(REPO/"output/phase2/scaled_35types/residuals_ranked.csv").set_index("cell_type")["residual_magnitude"]

# ---- primary: Table S5 matching_basis parenthetical markers ----
s5 = pd.read_csv(REPO/"docs/supplementary_materials/table_S5.csv")
def parse_markers(txt):
    m = re.search(r"\(([^)]*)\)", str(txt))
    if not m: return []
    inner = m.group(1)
    if any(w in inner.lower() for w in ["no specific","label","ontology","none"]): return []
    return [g.strip() for g in re.split(r"[,;/ ]+", inner) if g.strip() and re.match(r"^[A-Z0-9-]+$", g.strip())]
rows=[]
for _,r in s5.iterrows():
    ct=r["human_cell_type"]; mk=parse_markers(r["matching_basis"])
    if ct not in resid.index: continue
    ret = (sum(g in one2one for g in mk)/len(mk)) if mk else np.nan
    rows.append({"cell_type":ct,"n_markers":len(mk),"markers":";".join(mk),
                 "n_retained_1to1":sum(g in one2one for g in mk),
                 "retention_fraction":ret,"residual_magnitude":resid[ct]})
df5=pd.DataFrame(rows)
# Round for display (Table S8 family); correlation below uses full-precision df5.
_out5=df5.round({"residual_magnitude":3,"retention_fraction":4})
if "n_retained_1to1" in _out5.columns: _out5["n_retained_1to1"]=_out5["n_retained_1to1"].astype("Int64")
_out5.to_csv(OUT/"tableS5_retention.csv",index=False)
valid5=df5.dropna(subset=["retention_fraction"])
rho5,p5=spearmanr(valid5["retention_fraction"],valid5["residual_magnitude"]) if len(valid5)>=3 else (np.nan,np.nan)

# ---- secondary: CellMarker per-type ----
cm=pd.read_csv(REPO/"data/validation/cellmarker/cellmarker_human_filtered.csv")
def norm(s): return re.sub(r"[^a-z0-9]","",str(s).lower())
cm_by={}
for ct,grp in cm.groupby("cell_type"):
    cm_by[norm(ct)]=set(grp["gene_symbol"].astype(str))
rows2=[]
for ct in resid.index:
    key=norm(ct); markers=None
    if key in cm_by: markers=cm_by[key]
    else:
        cand=[v for k,v in cm_by.items() if key in k or k in key]
        if cand: markers=set().union(*cand)
    if not markers: 
        rows2.append({"cell_type":ct,"n_markers":0,"retention_fraction":np.nan,"residual_magnitude":resid[ct]}); continue
    ret=sum(g in one2one for g in markers)/len(markers)
    rows2.append({"cell_type":ct,"n_markers":len(markers),"n_retained_1to1":sum(g in one2one for g in markers),
                  "retention_fraction":ret,"residual_magnitude":resid[ct]})
df2=pd.DataFrame(rows2)
# Round for display (Table S8); correlation below uses full-precision df2.
_out2=df2.round({"residual_magnitude":3,"retention_fraction":4})
if "n_retained_1to1" in _out2.columns: _out2["n_retained_1to1"]=_out2["n_retained_1to1"].astype("Int64")
_out2.to_csv(OUT/"cellmarker_retention.csv",index=False)
valid2=df2.dropna(subset=["retention_fraction"])
rho2,p2=spearmanr(valid2["retention_fraction"],valid2["residual_magnitude"]) if len(valid2)>=3 else (np.nan,np.nan)

summary={
 "tableS5_primary":{"n_types_with_markers":int((df5["n_markers"]>0).sum()),
   "n_types_lt5_markers":int((df5["n_markers"]<5).sum()),
   "median_markers":float(df5.loc[df5.n_markers>0,"n_markers"].median()),
   "n_used":int(len(valid5)),"spearman_rho":float(rho5),"p":float(p5),
   "median_retention":float(valid5["retention_fraction"].median())},
 "cellmarker_secondary":{"n_types_with_markers":int((df2["n_markers"]>0).sum()),
   "n_types_lt5_markers":int((df2["n_markers"]<5).sum()),
   "median_markers":float(df2.loc[df2.n_markers>0,"n_markers"].median()) if (df2.n_markers>0).any() else None,
   "n_used":int(len(valid2)),"spearman_rho":float(rho2),"p":float(p2),
   "median_retention":float(valid2["retention_fraction"].median()) if len(valid2) else None},
}
json.dump(summary,open(OUT/"ortholog_retention_results.json","w"),indent=2)
print("\n===== M5 SUMMARY =====")
print(f"  Table S5 (primary): types w/markers={summary['tableS5_primary']['n_types_with_markers']}/35, "
      f"<5 markers={summary['tableS5_primary']['n_types_lt5_markers']}, median markers={summary['tableS5_primary']['median_markers']}, "
      f"median retention={summary['tableS5_primary']['median_retention']:.3f}")
print(f"     Spearman rho(retention, residual) = {rho5:.3f} (p={p5:.3f}, n={len(valid5)})")
print(f"  CellMarker (secondary): types matched={summary['cellmarker_secondary']['n_types_with_markers']}/35, "
      f"median markers={summary['cellmarker_secondary']['median_markers']}, "
      f"median retention={summary['cellmarker_secondary']['median_retention'] if summary['cellmarker_secondary']['median_retention'] else 'NA'}")
print(f"     Spearman rho(retention, residual) = {rho2:.3f} (p={p2:.3f}, n={len(valid2)})")
