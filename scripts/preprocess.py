"""Normalise, log-transform, select highly variable genes, PCA + elbow plot."""
import argparse

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

from common import choose_n_pcs, ensure_parent, load_config


def main():
    ap = argparse.ArgumentParser()
    for arg in ("--config", "--input", "--output", "--fig"):
        ap.add_argument(arg, required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    p = cfg["preprocess"]

    adata = sc.read_h5ad(args.input)
    adata.layers["counts"] = adata.X.copy()  # raw counts kept for pseudobulk DE
    sc.pp.normalize_total(adata, target_sum=p["target_sum"])
    sc.pp.log1p(adata)

    # HVGs per patient, then combined: avoids picking genes that only mark one donor
    sc.pp.highly_variable_genes(adata, n_top_genes=p["n_top_genes"], batch_key="patient")
    # Keep log-normalised values for plotting; scale a copy restricted to HVGs for PCA
    hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=10)
    n_comps = min(p["n_comps"], hvg.n_obs - 1, hvg.n_vars - 1)
    sc.tl.pca(hvg, n_comps=n_comps, random_state=cfg["seed"])
    adata.obsm["X_pca"] = hvg.obsm["X_pca"]
    adata.uns["pca"] = hvg.uns["pca"]

    var_ratio = adata.uns["pca"]["variance_ratio"]
    # Elbow on the log scale (as plotted): on a linear scale the first few PCs dominate
    # and the elbow lands too early, merging e.g. CD8 T and NK cells.
    n_pcs = choose_n_pcs(np.log(var_ratio)) if p["n_pcs"] == "auto" else int(p["n_pcs"])
    adata.uns["n_pcs_used"] = n_pcs
    print(f"HVGs: {int(adata.var['highly_variable'].sum())}, PCs used: {n_pcs}")

    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(np.arange(1, len(var_ratio) + 1), var_ratio, "o-", ms=3, color="#7a9cc6")
    ax.axvline(n_pcs, color="#c0392b", ls="--", label=f"elbow: {n_pcs} PCs")
    ax.set_xlabel("principal component")
    ax.set_ylabel("variance ratio")
    ax.set_yscale("log")
    ax.legend()
    fig.tight_layout()
    fig.savefig(ensure_parent(args.fig), dpi=120)

    adata.write_h5ad(ensure_parent(args.output))


if __name__ == "__main__":
    main()
