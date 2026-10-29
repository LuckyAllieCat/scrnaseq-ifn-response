"""Quality control: per-cell metrics, data-driven (MAD) thresholds, Scrublet doublets."""
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from common import ensure_parent, load_config, mad_bounds


def main():
    ap = argparse.ArgumentParser()
    for arg in ("--config", "--input", "--output", "--table", "--fig"):
        ap.add_argument(arg, required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    qc = cfg["qc"]

    adata = sc.read_h5ad(args.input)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True)

    # Thresholds come from the distributions themselves: cells more than n MADs
    # from the median (on log scale for counts/genes) are treated as outliers.
    n_mads = qc["n_mads"]
    lo_c, hi_c = mad_bounds(adata.obs["log1p_total_counts"], n_mads)
    lo_g, hi_g = mad_bounds(adata.obs["log1p_n_genes_by_counts"], n_mads)
    thresholds = {
        "min_counts": np.expm1(lo_c),
        "max_counts": np.expm1(hi_c),
        "min_genes": np.expm1(lo_g),
        "max_genes": np.expm1(hi_g),
    }
    has_mt = adata.var["mt"].sum() > 0
    if has_mt:
        _, hi_mt = mad_bounds(adata.obs["pct_counts_mt"], qc["n_mads_mt"])
        thresholds["max_pct_mt"] = min(hi_mt, qc["max_pct_mt_cap"])

    # Doublets: Scrublet run separately per capture batch (here: condition = lane pool)
    sc.pp.scrublet(adata, batch_key="condition", random_state=cfg["seed"])

    obs = adata.obs
    # Low gene count alone is not low quality: platelets/megakaryocytes have many UMIs
    # spread over few genes. Only drop low-gene cells that also have below-median counts.
    low_genes = (obs["n_genes_by_counts"] < thresholds["min_genes"]) & (
        obs["total_counts"] < obs["total_counts"].median())
    keep = (
        obs["total_counts"].between(thresholds["min_counts"], thresholds["max_counts"])
        & (obs["n_genes_by_counts"] <= thresholds["max_genes"])
        & ~low_genes
        & ~obs["predicted_doublet"].astype(bool)
    )
    if has_mt:
        keep &= obs["pct_counts_mt"] <= thresholds["max_pct_mt"]

    rows = [{"metric": k, "value": round(float(v), 3)} for k, v in thresholds.items()]
    rows += [
        {"metric": "mt_genes_found", "value": int(adata.var["mt"].sum())},
        {"metric": "doublets_flagged", "value": int(obs["predicted_doublet"].sum())},
        {"metric": "cells_before", "value": adata.n_obs},
        {"metric": "cells_after", "value": int(keep.sum())},
    ]

    # --- figure: distributions with threshold lines --------------------------------
    panels = [
        ("total_counts", "UMI counts / cell", ("min_counts", "max_counts"), True),
        ("n_genes_by_counts", "genes / cell", ("min_genes", "max_genes"), True),
        ("pct_counts_mt", "% mitochondrial", ("max_pct_mt",), False),
        ("doublet_score", "Scrublet doublet score", (), False),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.5))
    for ax, (col, label, keys, log) in zip(axes, panels):
        vals = obs[col].astype(float)
        bins = np.logspace(np.log10(max(vals.min(), 1)), np.log10(vals.max()), 60) if log else 60
        ax.hist(vals, bins=bins, color="#7a9cc6")
        if log:
            ax.set_xscale("log")
        for k in keys:
            if k in thresholds:
                ax.axvline(thresholds[k], color="#c0392b", ls="--")
        ax.set_xlabel(label)
        ax.set_ylabel("cells")
    axes[2].set_title("" if has_mt else "no MT- genes in matrix")
    fig.suptitle(f"QC: {keep.sum()} / {adata.n_obs} cells kept (red = MAD thresholds)")
    fig.tight_layout()
    fig.savefig(ensure_parent(args.fig), dpi=120)

    adata = adata[keep].copy()
    # Remove genes seen in almost no cells (uninformative, slow everything down)
    gene_ok = sc.pp.filter_genes(adata, min_cells=qc["min_cells_per_gene"], inplace=False)[0]
    rows.append({"metric": "genes_after", "value": int(gene_ok.sum())})
    adata = adata[:, gene_ok].copy()

    pd.DataFrame(rows).to_csv(ensure_parent(args.table), sep="\t", index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    adata.write_h5ad(ensure_parent(args.output))


if __name__ == "__main__":
    main()
