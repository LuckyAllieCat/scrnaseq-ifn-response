"""kNN graph, integration check (Harmony if needed), UMAP and Leiden clustering."""
import argparse

import harmonypy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from common import ensure_parent, load_config, neighbor_mixing


def main():
    ap = argparse.ArgumentParser()
    for arg in ("--config", "--input", "--output", "--table", "--fig"):
        ap.add_argument(arg, required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)
    c = cfg["cluster"]
    seed = cfg["seed"]

    adata = sc.read_h5ad(args.input)
    n_pcs = int(adata.uns["n_pcs_used"])
    adata.obsm["X_pca_used"] = adata.obsm["X_pca"][:, :n_pcs]

    def graph(rep):
        sc.pp.neighbors(adata, n_neighbors=c["n_neighbors"], use_rep=rep, random_state=seed)
        return {k: neighbor_mixing(adata.obsp["connectivities"], adata.obs[k])
                for k in c["integration_keys"]}

    # 1) How well are cells from different patients / conditions mixed before correction?
    before = graph("X_pca_used")
    to_fix = [k for k, v in before.items() if v < c["mixing_threshold"]]
    rows = [{"key": k, "mixing_before": round(v, 3), "integrated": k in to_fix} for k, v in before.items()]

    rep = "X_pca_used"
    if to_fix:
        # Harmony only changes the embedding used for clustering/annotation;
        # DE later uses raw counts, so the IFN signal itself is not removed.
        ho = harmonypy.run_harmony(adata.obsm["X_pca_used"], adata.obs, to_fix, random_state=seed)
        z = np.asarray(ho.Z_corr)
        # harmonypy has returned both (pcs x cells) and (cells x pcs) across versions
        adata.obsm["X_harmony"] = z.T if z.shape[0] != adata.n_obs else z
        rep = "X_harmony"
        after = graph(rep)
        for r in rows:
            r["mixing_after"] = round(after[r["key"]], 3)
    print(f"Integration on {to_fix or 'nothing'}; embedding = {rep}")
    adata.uns["embedding_used"] = rep

    sc.tl.umap(adata, random_state=seed)
    sc.tl.leiden(adata, resolution=c["leiden_resolution"], random_state=seed,
                 flavor="igraph", n_iterations=2, directed=False)
    print(f"Leiden clusters: {adata.obs['leiden'].nunique()}")

    pd.DataFrame(rows).to_csv(ensure_parent(args.table), sep="\t", index=False)
    print(pd.DataFrame(rows).to_string(index=False))

    fig, axes = plt.subplots(1, 4, figsize=(22, 5))
    for ax, col in zip(axes, ["leiden", "patient", "condition", "ref_cell_type"]):
        sc.pl.umap(adata, color=col, ax=ax, show=False, frameon=False,
                   legend_loc="on data" if col == "leiden" else "right margin", size=4)
    fig.tight_layout()
    fig.savefig(ensure_parent(args.fig), dpi=110)

    adata.write_h5ad(ensure_parent(args.output))


if __name__ == "__main__":
    main()
