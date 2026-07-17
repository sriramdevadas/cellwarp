#!/usr/bin/env python3
"""Assemble the five PLOS ONE main figures.
Reuses existing high-res panels where they exist; builds Fig 1D (mouse-lemur null)
and Fig 4A/B/C (forest, inversion, recovery ceiling) fresh from tracked data.
Outputs PDF + PNG (300 dpi) per figure. Arial, RGB; review-artifact quality."""
import pathlib, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.image as mpimg
import matplotlib.font_manager as fm
from scipy.stats import spearmanr

for cand in ("Arial","Helvetica","DejaVu Sans"):
    if any(cand in f.name for f in fm.fontManager.ttflist): plt.rcParams["font.family"]=cand; break
plt.rcParams.update({"font.size":8,"axes.linewidth":0.6,"xtick.major.width":0.6,
                     "ytick.major.width":0.6,"savefig.dpi":300})
BLUE,ORANGE,RED,GREY="#4a7fb5","#e08214","#c0392b","#555555"
ROOT=pathlib.Path.home()/"cellwarp"; P=ROOT/"figures/panels"; M=ROOT/"figures/main"
OUT=ROOT/"docs/submission/plosone/figures"; OUT.mkdir(exist_ok=True)

def place(ax,png,letter=None):
    ax.imshow(mpimg.imread(str(png))); ax.axis("off")
    if letter: ax.text(-0.02,1.02,letter,transform=ax.transAxes,fontsize=12,fontweight="bold",va="bottom",ha="right")

def letter(ax,l,x=-0.16,y=1.04):
    ax.text(x,y,l,transform=ax.transAxes,fontsize=12,fontweight="bold",va="bottom",ha="right")

def save(fig,name):
    fig.savefig(OUT/f"{name}.pdf"); fig.savefig(OUT/f"{name}.png",dpi=300); plt.close(fig)
    print("wrote",name)

# ---------- FIG 1 : configuration conserved ----------
fig=plt.figure(figsize=(7.4,4.8))
gs=fig.add_gridspec(2,3,height_ratios=[0.5,1.3],hspace=0.05,wspace=0.12)
axA=fig.add_subplot(gs[0,:]); place(axA,P/"fig1a_pipeline_schematic.png","A")
axB=fig.add_subplot(gs[1,0]); place(axB,P/"fig1b_null_1M.png","B")
axC=fig.add_subplot(gs[1,1]); place(axC,P/"fig1c_lineage_stratified.png","C")
# 1D fresh: mouse-lemur null
axD=fig.add_subplot(gs[1,2])
null=np.load(ROOT/"analysis/mouse_lemur/null_distribution.npy"); obs=21.765
axD.hist(null,bins=40,color="#9bb8d3",edgecolor="white",linewidth=0.3)
axD.axvline(obs,color=RED,lw=1.6)
axD.annotate(f"observed\n{obs:.1f}",xy=(obs,axD.get_ylim()[1]*0.5),xytext=(obs+7,axD.get_ylim()[1]*0.7),
             fontsize=7,color=RED,arrowprops=dict(arrowstyle="->",color=RED,lw=0.7))
axD.text(0.97,0.95,"obs/null 0.35\np < 0.0001\nn = 15",transform=axD.transAxes,ha="right",va="top",fontsize=7)
axD.set_title("Human–mouse lemur (~75 Mya)",fontsize=8)
axD.set_xlabel("Procrustes distance",fontsize=8); axD.set_ylabel("permutations",fontsize=8)
axD.spines[["top","right"]].set_visible(False); axD.tick_params(labelsize=7,length=2.5)
letter(axD,"D",x=-0.22)
save(fig,"Fig1_configuration_conserved")

# ---------- FIG 2 : two layers + BG replication ----------
fig=plt.figure(figsize=(7.4,6.2))
gs=fig.add_gridspec(2,2,height_ratios=[1.35,0.95],hspace=0.14,wspace=0.10)
place(fig.add_subplot(gs[0,0]),P/"fig3a_ellipsoid_heatmap.png","A")
place(fig.add_subplot(gs[0,1]),P/"fig3b_pre_post.png","B")
axC=fig.add_subplot(gs[1,:]); place(axC,OUT/"Fig2C_bg_replication.png","C")
save(fig,"Fig2_two_layers_bg")

# ---------- FIG 3 : configuration robust ----------
fig=plt.figure(figsize=(7.4,3.0))
place(fig.add_subplot(1,1,1),P/"fig4d_replication_summary.png")
save(fig,"Fig3_configuration_robust")

# ---------- FIG 4 : per-type not resolvable (all fresh) ----------
fig=plt.figure(figsize=(7.4,3.2))
gs=fig.add_gridspec(1,3,wspace=0.42,width_ratios=[1.05,1,1])
bs=pd.read_csv(ROOT/"analysis/bootstrap_rankings/bootstrap_summary.csv")
# 4A forest: CI per type, sorted by original_rank
a=fig.add_subplot(gs[0,0]); bs2=bs.sort_values("original_rank")
y=np.arange(len(bs2))
a.hlines(y,bs2["ci_lower"],bs2["ci_upper"],color=BLUE,lw=1.3)
a.plot(bs2["original_rank"],y,"o",ms=2.4,color="#26456e")
a.set_ylim(-1,len(bs2)); a.invert_yaxis()
a.set_xlabel("divergence rank (1 = most divergent)",fontsize=7.5)
a.set_ylabel("cell type (ranked)",fontsize=7.5)
a.text(0.96,0.04,"median 95% CI width = 3 ranks\nall 35 types stable",transform=a.transAxes,
       ha="right",va="bottom",fontsize=6.6)
a.set_yticks([]); a.tick_params(labelsize=7,length=2.5); a.spines[["top","right"]].set_visible(False)
a.set_title("Within-atlas precision",fontsize=8); letter(a,"A",x=-0.10)
# 4B inversion: ci_width vs cross-atlas mean_rank_shift (master_ranking_table = the
# exact source behind Table 1 T59: rho=-0.410, p=0.073, n=20; types in >=2 replications)
mrt=pd.read_csv(ROOT/"analysis/cross_reference/master_ranking_table.csv")
mg=mrt[mrt["n_replications_present"]>=2].dropna(subset=["bootstrap_CI_width","mean_rank_shift"]).rename(
    columns={"bootstrap_CI_width":"ci_width"})
b=fig.add_subplot(gs[0,1])
b.scatter(mg["ci_width"],mg["mean_rank_shift"],s=16,color=ORANGE,edgecolor="white",linewidth=0.4)
rho,_=spearmanr(mg["ci_width"],mg["mean_rank_shift"])
z=np.polyfit(mg["ci_width"],mg["mean_rank_shift"],1); xs=np.array([mg["ci_width"].min(),mg["ci_width"].max()])
b.plot(xs,np.polyval(z,xs),color=GREY,lw=1,ls="--")
b.text(0.96,0.95,f"Spearman ρ = {rho:.2f}\n(n = {len(mg)})",transform=b.transAxes,ha="right",va="top",fontsize=7)
b.set_xlabel("within-atlas 95% CI width",fontsize=7.5); b.set_ylabel("cross-atlas mean rank shift",fontsize=7.5)
b.tick_params(labelsize=7,length=2.5); b.spines[["top","right"]].set_visible(False)
b.set_title("Precision vs reproducibility",fontsize=8); letter(b,"B")
# 4C recovery ceiling
rr=json.load(open(ROOT/"analysis/simulation_study/simulation_results.json"))["ranking_recovery"]
c=fig.add_subplot(gs[0,2])
for sig,col,lab in [(3.0,BLUE,"signal ≈ data (3.0)"),(5.0,"#9bb8d3","signal 5.0")]:
    pts=sorted((r["n_cells_per_type"],r["median_rho"]) for r in rr if r["signal_strength"]==sig)
    xs=[p[0] for p in pts]; ys=[p[1] for p in pts]
    c.plot(xs,ys,"o-",color=col,ms=4,lw=1.2,label=lab)
c.axhline(0.42,color=RED,lw=1,ls=":"); c.text(2000,0.435,"ceiling ρ ≈ 0.42",color=RED,ha="right",fontsize=6.6)
c.set_xscale("log"); c.set_xticks([50,200,500,2000]); c.set_xticklabels([50,200,500,2000])
c.set_ylim(0,0.75); c.set_xlabel("cells sampled per type",fontsize=7.5); c.set_ylabel("recovery ρ (true vs estimated)",fontsize=7.5)
c.legend(fontsize=6.2,frameon=False,loc="upper left"); c.tick_params(labelsize=7,length=2.5)
c.spines[["top","right"]].set_visible(False); c.set_title("Recovery saturates",fontsize=8); letter(c,"C")
fig.subplots_adjust(left=0.075,right=0.985,top=0.88,bottom=0.17)
save(fig,"Fig4_pertype_not_resolvable")

# ---------- FIG 5 : conserved identity genes ----------
fig=plt.figure(figsize=(7.4,5.2))
place(fig.add_subplot(1,1,1),M/"fig7_conserved_contribution.png")
save(fig,"Fig5_conserved_identity_genes")
print("done: 5 main figures")
