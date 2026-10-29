"""Annotate Leiden clusters from canonical PBMC markers; compare with the authors' labels."""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score

from common import ensure_parent, load_config


def main():
    ap = argparse.ArgumentParser()
    for arg in ("--config", "--input", "--output", "--outdir", "--dotplot", "--umap"):
        ap.add_argument(arg, required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input)
    markers = {ct: [g for g in genes if g in adata.var_names]
               for ct, genes in cfg["annotation"]["markers"].items()}
    markers = {ct: g for ct, g in markers.items() if g}

    # Score every cell for each marker set, average per cluster, then z-score each
    # cell type's score across clusters so sets of different size are comparable.
    for ct, genes in markers.items():
        sc.tl.score_genes(adata, genes, score_name=f"score_{ct}", random_state=cfg["seed"])
    score_cols = [f"score_{ct}" for ct in markers]
    per_cluster = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
    per_cluster.columns = list(markers)
    z = (per_cluster - per_cluster.mean()) / per_cluster.std(ddof=0).replace(0, 1)
    assignment = z.idxmax(axis=1)

    table = per_cluster.round(3).assign(assigned=assignment, z_best=z.max(axis=1).round(2),
                                        n_cells=adata.obs["leiden"].value_counts())
    table.to_csv(outdir / "cluster_annotation.tsv", sep="\t")
    adata.obs["cell_type"] = adata.obs["leiden"].map(assignment).astype("category")
    adata.obs = adata.obs.drop(columns=score_cols)

    # Agreement with the authors' labels (same names by design, finer labels mapped back)
    ref = adata.obs["ref_cell_type"].astype(str)
    ours = adata.obs["cell_type"].astype(str).replace(cfg["annotation"].get("reference_names") or {})
    agreement = {
        "cells": adata.n_obs,
        "exact_label_agreement": round(float((ref == ours).mean()), 4),
        "adjusted_rand_index": round(adjusted_rand_score(ref, ours), 4),
    }
    pd.Series(agreement).to_csv(outdir / "agreement_summary.tsv", sep="\t", header=False)
    pd.crosstab(ours, ref).to_csv(outdir / "confusion_ours_vs_reference.tsv", sep="\t")
    per_type = (ref == ours).groupby(ref).mean().round(3).rename("recall_vs_reference")
    per_type.to_csv(outdir / "agreement_per_reference_type.tsv", sep="\t")
    print(agreement)
    print(table.to_string())

    dp = sc.pl.dotplot(adata, markers, groupby="leiden", standard_scale="var",
                       show=False, return_fig=True)
    dp.savefig(ensure_parent(args.dotplot), dpi=120, bbox_inches="tight")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    sc.pl.umap(adata, color="cell_type", ax=axes[0], show=False, frameon=False, size=4,
               title="cell type (marker-based)")
    sc.pl.umap(adata, color="condition", ax=axes[1], show=False, frameon=False, size=4,
               palette={"ctrl": "#9bb7d4", "stim": "#e07b7b"}, title="condition")
    fig.tight_layout()
    fig.savefig(ensure_parent(args.umap), dpi=120)

    adata.write_h5ad(ensure_parent(args.output))


if __name__ == "__main__":
    main()
