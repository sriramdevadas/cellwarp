"""Six-panel summary figure for the conserved-contribution gate."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import gate_lib as G
HERE = Path(__file__).resolve().parent

df, h_cent, m_cent = G.build_gene_table()
valid = df.dropna(subset=["C_pearson"]).reset_index(drop=True)
C = valid["C_pearson"].values
me = valid["mean_expression"].values
Hc = h_cent[valid["gene_id"]].values
Tau = np.array([np.sum(1 - Hc[:, j] / Hc[:, j].max()) / (Hc.shape[0] - 1)
                if Hc[:, j].max() > 0 else np.nan for j in range(len(valid))])
res = json.load(open(HERE / "gate_results.json"))
rob = json.load(open(HERE / "robustness_results.json"))
q75, q25 = np.quantile(C, 0.75), np.quantile(C, 0.25)
sym_to_pos = {s: i for i, s in enumerate(valid["symbol"].values)}
tf_pos = np.array(sorted({sym_to_pos[s] for s in G.POSITIVE_CONTROL_TFS if s in sym_to_pos}))

fig, ax = plt.subplots(2, 3, figsize=(15, 9))

# A: distribution
a = ax[0, 0]
a.hist(C, bins=60, color="#4C72B0", alpha=0.85)
a.axvline(q75, color="green", ls="--", label=f"conserved Q75={q75:.2f}")
a.axvline(q25, color="red", ls="--", label=f"divergent Q25={q25:.2f}")
a.set_title(f"A. Per-gene conservation C (n={len(C)})\ndip p<1e-4 but D={rob['dip_D']:.3f}: broad continuum")
a.set_xlabel("C = cross-species Pearson r across 35 types"); a.set_ylabel("genes"); a.legend(fontsize=8)

# B: C vs expression
b = ax[0, 1]
b.hexbin(me, C, gridsize=45, cmap="Blues", mincnt=1)
dec = pd.qcut(pd.Series(me).rank(method="first"), 10, labels=False)
dm = [C[dec == d].mean() for d in range(10)]
dx = [me[dec == d].mean() for d in range(10)]
b.plot(dx, dm, "o-", color="orange", label="decile mean C")
b.set_title(f"B. C vs expression\nSpearman={res['check2']['spearman_C_vs_expr']:.3f} (low-expr noise floor)")
b.set_xlabel("mean expression (log1p)"); b.set_ylabel("C"); b.legend(fontsize=8)

# C: C vs Tau
c = ax[0, 2]
c.hexbin(Tau, C, gridsize=45, cmap="Purples", mincnt=1)
c.set_title(f"C. C vs specificity (Tau)\nSpearman={rob['rho_C_tau']:.3f}: conservation != specificity")
c.set_xlabel("Tau specificity"); c.set_ylabel("C")

# D: positive control
d = ax[1, 0]
Crank = stats.rankdata(C) / len(C)
d.hist(Crank, bins=40, color="lightgray", label="all genes", density=True)
d.hist(Crank[tf_pos], bins=20, color="crimson", alpha=0.7, label=f"identity TFs (n={len(tf_pos)})", density=True)
d.axvline(np.median(Crank[tf_pos]), color="crimson", ls="--")
d.set_title(f"D. Identity-TF positive control\nmedian percentile={res['check3a']['median_Crank']:.2f}, p={res['check3a']['p_value']:.3f}")
d.set_xlabel("C percentile rank"); d.set_ylabel("density"); d.legend(fontsize=8)

# E: enrichment bars
e = ax[1, 1]
labels = ["TF\n(expr-match)", "TF\n(expr+Tau)", "CellMarker\n(expr-match)", "CellMarker\n(expr+Tau)"]
folds = [res["check3b"]["H_tf"]["fold"], rob["joint_matched_3b"]["TF"]["fold"],
         res["check3b"]["H_cellmarker"]["fold"], rob["joint_matched_3b"]["CellMarker"]["fold"]]
cols = ["#2a9d8f", "#264653", "#e76f51", "#9c2706"]
e.bar(labels, folds, color=cols)
e.axhline(1.0, color="k", ls=":")
for i, f in enumerate(folds):
    e.text(i, f + 0.02, f"{f:.2f}x", ha="center", fontsize=9)
e.set_title("E. Conserved-set enrichment\nvs matched background (all p<0.05)")
e.set_ylabel("fold enrichment"); e.tick_params(axis="x", labelsize=8)

# F: attribution
f = ax[1, 2]
sec = res["secondary"]
names = ["all\ngenes", "conserved", "expr-matched\nrandom", "divergent"]
vals = [sec["validity_all_genes_ratio"], sec["conserved"]["ratio"],
        sec["matched_random_ratio_mean"], sec["divergent"]["ratio"]]
errs = [0, 0, sec["matched_random_ratio_sd"], 0]
f.bar(names, vals, yerr=errs, color=["gray", "green", "steelblue", "red"], capsize=4)
f.axhline(sec["validity_all_genes_ratio"], color="gray", ls=":")
for i, v in enumerate(vals):
    f.text(i, v + 0.01, f"{v:.3f}", ha="center", fontsize=9)
f.set_title("F. Geometry attribution (obs/null)\nlower = better cross-species alignment")
f.set_ylabel("Procrustes obs/null ratio"); f.tick_params(axis="x", labelsize=8)

plt.tight_layout()
plt.savefig(HERE / "gate_figure.png", dpi=140, bbox_inches="tight")
print("saved", HERE / "gate_figure.png")
