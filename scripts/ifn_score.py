"""Per-cell interferon-stimulated gene (ISG) score, compared by cell type and condition."""
import argparse

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

from common import ensure_parent, load_config


def main():
    ap = argparse.ArgumentParser()
    for arg in ("--config", "--input", "--table", "--effect", "--fig"):
        ap.add_argument(arg, required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)

    adata = sc.read_h5ad(args.input)
    genes = [g for g in cfg["ifn"]["genes"] if g in adata.var_names]
    print(f"ISG genes found: {len(genes)}/{len(cfg['ifn']['genes'])}")
    sc.tl.score_genes(adata, genes, score_name="ifn_score", random_state=cfg["seed"])

    # Average per patient first: patients, not cells, are the independent replicates
    per_patient = (adata.obs.groupby(["cell_type", "condition", "patient"], observed=True)
                   ["ifn_score"].agg(["mean", "size"]).reset_index()
                   .rename(columns={"mean": "ifn_score", "size": "n_cells"}))
    per_patient.to_csv(ensure_parent(args.table), sep="\t", index=False)

    wide = per_patient.pivot_table(index=["cell_type", "patient"], columns="condition",
                                   values="ifn_score", observed=True).dropna()
    wide["delta"] = wide["stim"] - wide["ctrl"]
    effect = (wide.groupby("cell_type", observed=True)["delta"]
              .agg(mean_delta="mean", sd_delta="std", n_patients="size")
              .sort_values("mean_delta", ascending=False).round(3))
    effect.to_csv(ensure_parent(args.effect), sep="\t")
    print(effect.to_string())

    order = list(effect.index) + [c for c in adata.obs["cell_type"].cat.categories if c not in effect.index]
    pal = {"ctrl": "#9bb7d4", "stim": "#e07b7b"}
    fig, ax = plt.subplots(figsize=(11, 4.5))
    sns.violinplot(data=adata.obs, x="cell_type", y="ifn_score", hue="condition", order=order,
                   split=True, inner=None, cut=0, linewidth=0.5, palette=pal, ax=ax)
    sns.stripplot(data=per_patient, x="cell_type", y="ifn_score", hue="condition", order=order,
                  dodge=True, palette=pal, edgecolor="black", linewidth=0.6, size=4,
                  legend=False, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("ISG score")
    ax.set_title("Interferon response per cell type (violins = cells, dots = patient means)")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(ensure_parent(args.fig), dpi=120)


if __name__ == "__main__":
    main()
