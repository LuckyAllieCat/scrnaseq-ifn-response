"""Load the raw h5ad, standardise column names, optionally subsample (stratified)."""
import argparse

import numpy as np
import scanpy as sc
from scipy import sparse

from common import ensure_parent, load_config


def main():
    ap = argparse.ArgumentParser()
    for arg in ("--config", "--input", "--output"):
        ap.add_argument(arg, required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)

    adata = sc.read_h5ad(args.input)
    # Keep only raw counts + the metadata we need; drop the authors' Seurat embeddings
    obs = adata.obs.rename(columns={"label": "condition", "replicate": "patient",
                                    "cell_type": "ref_cell_type"})
    adata = sc.AnnData(
        X=sparse.csr_matrix(adata.X, dtype=np.float32),
        obs=obs[["condition", "patient", "ref_cell_type"]].copy(),
        var=adata.var[[]].copy(),
    )
    adata.obs["condition"] = adata.obs["condition"].astype("category").cat.reorder_categories(["ctrl", "stim"])
    adata.var_names_make_unique()

    n = cfg["subsample_n"]
    if n and n < adata.n_obs:
        # Stratify by patient x condition so every pseudobulk sample survives subsampling
        rng = np.random.default_rng(cfg["seed"])
        groups = adata.obs.groupby(["patient", "condition"], observed=True).indices
        frac = n / adata.n_obs
        idx = np.concatenate([
            rng.choice(ix, size=max(1, int(round(len(ix) * frac))), replace=False)
            for ix in groups.values()
        ])
        adata = adata[np.sort(idx)].copy()
        print(f"Subsampled to {adata.n_obs} cells")

    print(adata)
    adata.write_h5ad(ensure_parent(args.output))


if __name__ == "__main__":
    main()
