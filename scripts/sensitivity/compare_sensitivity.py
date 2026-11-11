"""Compare original JAX/MSK DE results with the Alternative model.

Original:    WT replicates independent (each WT sample its own block), as in the
             primary notebooks.
Alternative: WT replicates from one block; within-clone correlation floored at 0
             (output of DE_WT_single_block.R).

Usage: python3 compare_sensitivity.py <DE_results dir>
Reads   <DE_results>/<target>_DE_<JAX|MSK>.tsv          (original)
        <DE_results>/sensitivity/<target>_DE_<...>.tsv  (alternative)
        <DE_results>/sensitivity/model_summary.tsv, original_model_correlations.tsv
Writes  <DE_results>/sensitivity/comparison_summary.csv
        <DE_results>/sensitivity/comparison_DE_gene_changes.tsv
        <DE_results>/sensitivity/comparison_DE_counts.pdf / .png
        <DE_results>/sensitivity/comparison_t_statistics.pdf / .png
        <DE_results>/sensitivity/comparison_log2FC.pdf / .png
DE call = exported `threshold` (q < 0.05, |logFC| > 1, max group-mean CPM >= 1).
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

root = Path(sys.argv[1])
sens = root / "sensitivity"
TARGETS = ["EOMES", "GCM1", "ISL1", "NKX2-1"]
CONTRASTS = ["JAX", "MSK"]
LABELS = {"original": "Original (WT replicates independent)",
          "alternative": "Alternative (WT replicates from one block)"}


def load(d, name):
    t = pd.read_csv(d / f"{name}.tsv", sep="\t").set_index("ensembl_gene_id")
    t["threshold"] = t["threshold"].astype(str).str.upper() == "TRUE"
    return t


alt_summary = pd.read_csv(sens / "model_summary.tsv", sep="\t").set_index("contrast")
orig_summary = pd.read_csv(sens / "original_model_correlations.tsv", sep="\t").set_index("contrast")

rows, changes, tables = [], [], {}
for g in TARGETS:
    for c in CONTRASTS:
        name = f"{g}_DE_{c}"
        o = load(root, name)
        n = load(sens, name)
        assert n.index.sort_values().equals(o.index.sort_values())
        n = n.loc[o.index]
        tables[name] = (o, n)
        both, union = o.threshold & n.threshold, o.threshold | n.threshold
        rows.append({
            "contrast": name,
            "rho_original": orig_summary.loc[name, "dupcor"],
            "rho_alternative_estimated": alt_summary.loc[name, "dupcor"],
            "rho_alternative_used": alt_summary.loc[name, "dupcor_used"],
            "DE_original": int(o.threshold.sum()),
            "DE_alternative": int(n.threshold.sum()),
            "kept": int(both.sum()),
            "lost": int((o.threshold & ~n.threshold).sum()),
            "gained": int((~o.threshold & n.threshold).sum()),
            "jaccard": round(both.sum() / max(union.sum(), 1), 3),
            "spearman_t": round(o.t.corr(n.t, method="spearman"), 4),
            "median_abs_t_ratio": round(float(np.median(n.t.abs() / o.t.abs())), 3),
            "spearman_logFC": round(o.logFC.corr(n.logFC, method="spearman"), 4),
            "median_abs_logFC_diff": round(float((n.logFC - o.logFC).abs().median()), 4),
            "max_abs_logFC_diff": round(float((n.logFC - o.logFC).abs().max()), 3),
            "top100_absLogFC_overlap": len(set(o.logFC.abs().nlargest(100).index)
                                           & set(n.logFC.abs().nlargest(100).index)),
            "top100_abs_t_overlap": len(set(o.t.abs().nlargest(100).index)
                                        & set(n.t.abs().nlargest(100).index)),
        })
        for gid in o.index[o.threshold != n.threshold]:
            changes.append({"contrast": name, "ensembl_gene_id": gid,
                            "gene_symbol": o.at[gid, "gene_symbol"],
                            "change": "lost" if o.at[gid, "threshold"] else "gained",
                            "logFC_original": o.at[gid, "logFC"], "logFC_alternative": n.at[gid, "logFC"],
                            "qvalue_original": o.at[gid, "qvalue"], "qvalue_alternative": n.at[gid, "qvalue"]})

out = pd.DataFrame(rows)
out.to_csv(sens / "comparison_summary.csv", index=False)
pd.DataFrame(changes).to_csv(sens / "comparison_DE_gene_changes.tsv", sep="\t", index=False)
with pd.option_context("display.width", 250, "display.max_columns", 40):
    print(out.to_string(index=False))

# ---- figures ----
os.environ.setdefault("MPLCONFIGDIR", "/tmp/morphic_sensitivity_matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False})
colors = {"original": "#8C9AA6", "alternative": "#287D9B"}


def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(sens / f"{stem}.{ext}", dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


fig, ax = plt.subplots(figsize=(11, 4.8))
x = np.arange(len(out))
w = 0.38
for i, k in enumerate(LABELS):
    vals = out[f"DE_{k}"].to_numpy()
    bars = ax.bar(x + (i - 0.5) * w, vals, w, color=colors[k], label=LABELS[k])
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 40, f"{v:,}", ha="center", va="bottom", fontsize=8.5)
ax.set_xticks(x, [n.replace("_DE_", "\n") for n in out.contrast])
ax.set_ylabel("Number of DE genes")
ax.set_ylim(0, out[["DE_original", "DE_alternative"]].to_numpy().max() * 1.2)
ax.set_title("DE genes vs parental KOLF2.2J WT: sensitivity to WT replicate dependence\n"
             "q < 0.05, |log2FC| > 1, CPM ≥ 1", loc="left", fontsize=13)
ax.legend(frameon=False, fontsize=9.5, loc="upper left", ncol=2)
fig.tight_layout()
save(fig, "comparison_DE_counts")

by_name = out.set_index("contrast")


def scatter_panels(column, axis_label, stem):
    """Original vs Alternative for one column; DE genes (either model) highlighted.
    Points are drawn as vectors (rasterized layers were misplaced in the PDF)."""
    fig, axes = plt.subplots(2, 4, figsize=(15, 7.5))
    for ax, (name, (o, n)) in zip(axes.flat, tables.items()):
        de = (o.threshold | n.threshold).to_numpy()
        xv, yv = o[column].to_numpy(), n[column].to_numpy()
        lim = np.nanmax(np.abs(np.r_[xv, yv])) * 1.05
        ax.plot([-lim, lim], [-lim, lim], color="#C85B38", lw=0.6, ls="--", zorder=0)
        ax.scatter(xv[~de], yv[~de], s=1.5, alpha=0.25, color="#8C9AA6", linewidths=0, label="not DE")
        ax.scatter(xv[de], yv[de], s=1.5, alpha=0.4, color="#287D9B", linewidths=0,
                   label="DE in either model")
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        r = by_name.loc[name]
        rho0 = "n/a" if pd.isna(r.rho_original) else f"{r.rho_original:.2f}"
        stat = (f"Spearman {r.spearman_logFC:.4f}; median |Δ| {r.median_abs_logFC_diff:.3f}"
                if column == "logFC" else f"DE {int(r.DE_original):,} → {int(r.DE_alternative):,}")
        ax.set_title(f"{name.replace('_DE_', ' ')}\nρ {rho0} → {r.rho_alternative_used:.2f}; {stat}",
                     fontsize=9.5, loc="left")
        ax.set_xlabel(f"{axis_label}, original")
        ax.set_ylabel(f"{axis_label}, alternative")
    axes.flat[0].legend(frameon=False, fontsize=8, markerscale=5, loc="upper left")
    fig.tight_layout()
    save(fig, stem)


scatter_panels("t", "moderated t", "comparison_t_statistics")
scatter_panels("logFC", "log2 fold change", "comparison_log2FC")
