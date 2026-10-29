"""Pseudobulk differential expression (stim vs ctrl) per cell type with PyDESeq2.

Counts are summed per patient x condition x cell type, so each patient contributes one
sample per condition. The design ~ patient + condition is a paired test: it compares
stim with ctrl within each patient, avoiding the pseudoreplication of cell-level tests.
"""
import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

from common import ensure_parent, load_config, pseudobulk


def run_deseq(counts, meta, threads):
    inference = DefaultInference(n_cpus=threads)
    dds = DeseqDataSet(counts=counts, metadata=meta, design="~patient + condition",
                       refit_cooks=True, inference=inference, quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", "stim", "ctrl"], inference=inference, quiet=True)
    stats.summary()
    return stats.results_df


def main():
    ap = argparse.ArgumentParser()
    for arg in ("--config", "--input", "--outdir", "--volcano", "--heatmap"):
        ap.add_argument(arg, required=True)
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()
    cfg = load_config(args.config)
    d = cfg["de"]
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input)
    counts, meta = pseudobulk(adata.layers["counts"], adata.obs,
                              ["cell_type", "patient", "condition"])
    counts.columns = adata.var_names
    meta.to_csv(outdir / "pseudobulk_samples.tsv", sep="\t")

    results, summary = {}, []
    for ct in meta["cell_type"].unique():
        m = meta[(meta["cell_type"] == ct) & (meta["n_cells"] >= d["min_cells_per_sample"])]
        paired = m.groupby("patient")["condition"].nunique()
        m = m[m["patient"].isin(paired[paired == 2].index)]
        n_pat = m["patient"].nunique()
        if n_pat < d["min_patients"]:
            print(f"Skipping {ct}: only {n_pat} patients with both conditions")
            summary.append({"cell_type": ct, "n_patients": n_pat, "tested": False})
            continue
        c = counts.loc[m.index].round().astype(int)
        c = c.loc[:, c.sum() >= d["min_count_per_gene"]]
        md = m[["patient", "condition"]].astype(str)
        res = run_deseq(c, md, args.threads)
        res = res.sort_values("padj")
        res.to_csv(outdir / f"de_{ct.replace(' ', '_').replace('+', 'pos')}.tsv", sep="\t")
        results[ct] = res
        sig = res[(res["padj"] < d["padj"]) & (res["log2FoldChange"].abs() >= d["lfc"])]
        summary.append({"cell_type": ct, "n_patients": n_pat, "tested": True,
                        "n_cells": int(m["n_cells"].sum()), "genes_tested": c.shape[1],
                        "n_up": int((sig["log2FoldChange"] > 0).sum()),
                        "n_down": int((sig["log2FoldChange"] < 0).sum())})
        print(f"{ct}: {len(sig)} DE genes")

    summary = pd.DataFrame(summary)
    if "n_up" in summary:
        summary["n_de"] = summary["n_up"] + summary["n_down"]
        summary = summary.sort_values("n_de", ascending=False)
    summary.to_csv(outdir / "de_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))

    # --- shared vs cell-type-specific genes -----------------------------------------
    sig_sets = {ct: set(r.index[(r["padj"] < d["padj"]) & (r["log2FoldChange"].abs() >= d["lfc"])])
                for ct, r in results.items()}
    n_types = len(sig_sets)
    all_genes = sorted(set().union(*sig_sets.values())) if sig_sets else []
    share = pd.DataFrame({"gene": all_genes})
    share["n_cell_types"] = [sum(g in s for s in sig_sets.values()) for g in all_genes]
    share["cell_types"] = [";".join(ct for ct, s in sig_sets.items() if g in s) for g in all_genes]
    shared_cut = max(2, math.ceil(0.75 * n_types))  # DE in >= 75% of tested cell types
    share["category"] = np.select([share["n_cell_types"] >= shared_cut, share["n_cell_types"] == 1],
                                  ["shared", "specific"], "partial")
    # strongest evidence first: smallest adjusted p in any cell type where the gene is DE
    share["best_padj"] = [min(results[ct].at[g, "padj"] for ct in cts.split(";"))
                          for g, cts in zip(share["gene"], share["cell_types"])]
    share = share.sort_values(["n_cell_types", "best_padj"], ascending=[False, True])
    share.to_csv(
        outdir / "gene_sharing.tsv", sep="\t", index=False)

    # --- volcano plots for the largest cell types -----------------------------------
    tested = summary[summary["tested"]].sort_values("n_cells", ascending=False) if results else summary
    show = list(tested["cell_type"][: d["volcano_cell_types"]]) if results else []
    fig, axes = plt.subplots(1, max(len(show), 1), figsize=(4.5 * max(len(show), 1), 4), squeeze=False)
    for ax, ct in zip(axes[0], show):
        r = results[ct].dropna(subset=["padj"])
        y = -np.log10(r["padj"].clip(lower=1e-300))
        is_sig = (r["padj"] < d["padj"]) & (r["log2FoldChange"].abs() >= d["lfc"])
        ax.scatter(r["log2FoldChange"], y, s=4, c=np.where(is_sig, "#e07b7b", "#bbbbbb"), rasterized=True)
        for g in r[is_sig].nsmallest(8, "padj").index:
            ax.annotate(g, (r.at[g, "log2FoldChange"], y[g]), fontsize=7)
        ax.axhline(-np.log10(d["padj"]), color="grey", ls=":", lw=0.8)
        for v in (-d["lfc"], d["lfc"]):
            ax.axvline(v, color="grey", ls=":", lw=0.8)
        ax.set_title(f"{ct} ({int(is_sig.sum())} DE)")
        ax.set_xlabel("log2 fold change (stim / ctrl)")
        ax.set_ylabel("-log10 adjusted p")
    if not show:
        axes[0][0].text(0.5, 0.5, "no cell type had enough replicates", ha="center")
    fig.tight_layout()
    fig.savefig(ensure_parent(args.volcano), dpi=120)

    # --- heatmap of top genes per cell type (log2FC) --------------------------------
    # Rows: strongest shared genes (by mean |log2FC|) + the top specific genes per type
    top = []
    if results:
        lfc_all = pd.DataFrame({ct: r["log2FoldChange"] for ct, r in results.items()})
        shared_genes = share.loc[share["category"] == "shared", "gene"]
        top += list(lfc_all.loc[shared_genes].abs().mean(axis=1).nlargest(d["top_n_shared"]).index)
        specific = share[share["category"] == "specific"].set_index("gene")["cell_types"]
        for ct, r in results.items():
            genes = specific.index[specific == ct]
            genes = genes[r.loc[genes, "baseMean"] >= d["min_base_mean_heatmap"]]  # skip rare genes
            top += list(r.loc[genes, "log2FoldChange"].abs().nlargest(d["top_n_heatmap"]).index)
    top = list(dict.fromkeys(top))
    fig, ax = plt.subplots(figsize=(1.2 * max(len(results), 2) + 3, 0.22 * max(len(top), 5) + 2))
    if top:
        lfc = pd.DataFrame({ct: r["log2FoldChange"].reindex(top) for ct, r in results.items()})
        cat = share.set_index("gene")["category"].reindex(top)
        order = cat.map({"shared": 0, "partial": 1, "specific": 2}).sort_values(kind="stable").index
        lfc = lfc.loc[order]
        lim = np.nanpercentile(np.abs(lfc.values), 98)
        ax.set_facecolor("#d0d0d0")  # masked (untested) cells show as grey
        sns.heatmap(lfc, cmap="RdBu_r", center=0, vmin=-lim, vmax=lim, ax=ax,
                    mask=lfc.isna(), cbar_kws={"label": "log2FC (stim / ctrl)"})
        ax.set_yticks(np.arange(len(order)) + 0.5)
        ax.set_yticklabels([f"{g}  [{cat[g]}]" for g in order], fontsize=7)
        ax.set_title("Shared vs cell-type-specific IFN-response genes (grey = not tested)")
    else:
        ax.text(0.5, 0.5, "no significant genes", ha="center")
    fig.tight_layout()
    fig.savefig(ensure_parent(args.heatmap), dpi=120)


if __name__ == "__main__":
    main()
