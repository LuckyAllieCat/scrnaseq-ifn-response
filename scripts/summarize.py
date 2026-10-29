"""Collect the key numbers from every step into results/summary.md."""
import argparse
from pathlib import Path

import pandas as pd

from common import ensure_parent, load_config


def md_table(df):
    df = df.reset_index() if df.index.name else df
    head = "| " + " | ".join(map(str, df.columns)) + " |"
    sep = "|" + "---|" * len(df.columns)
    rows = ["| " + " | ".join(map(str, r)) + " |" for r in df.itertuples(index=False)]
    return "\n".join([head, sep, *rows])


def main():
    ap = argparse.ArgumentParser()
    for arg in ("--config", "--qc", "--integration", "--annotation", "--ifn", "--de", "--output"):
        ap.add_argument(arg, required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)

    qc = pd.read_csv(args.qc, sep="\t").set_index("metric")["value"]
    integ = pd.read_csv(args.integration, sep="\t")
    agree = pd.read_csv(Path(args.annotation) / "agreement_summary.tsv", sep="\t",
                        header=None, index_col=0)[1]
    ifn = pd.read_csv(args.ifn, sep="\t")
    de = pd.read_csv(Path(args.de) / "de_summary.tsv", sep="\t")
    share = pd.read_csv(Path(args.de) / "gene_sharing.tsv", sep="\t")

    shared = share[share["category"] == "shared"]
    specific = share[share["category"] == "specific"]
    spec_by_type = (specific.groupby("cell_types")["gene"]
                    .apply(lambda g: ", ".join(list(g)[:8])).rename("example specific genes"))
    tested = de[de["tested"]] if "tested" in de else de

    lines = [
        "# Results summary",
        "",
        f"Cells: {int(qc['cells_before'])} loaded, {int(qc['cells_after'])} after QC "
        f"({int(qc['doublets_flagged'])} flagged as doublets). Subsample: {cfg['subsample_n']}.",
        "",
        "## Integration check (kNN mixing, 1 = perfectly mixed)",
        md_table(integ),
        "",
        "## Annotation agreement with the authors' labels",
        f"Exact label agreement {float(agree['exact_label_agreement']):.1%}, "
        f"adjusted Rand index {float(agree['adjusted_rand_index']):.3f}.",
        "",
        "## Interferon response (ISG score, stim - ctrl, mean over patients)",
        md_table(ifn),
        "",
        "## Pseudobulk DE (PyDESeq2, ~ patient + condition)",
        f"Significant = padj < {cfg['de']['padj']} and |log2FC| >= {cfg['de']['lfc']}.",
        "",
        md_table(tested),
        "",
        f"**Shared** genes (DE in >= 75% of tested cell types): {len(shared)}. "
        f"Top: {', '.join(shared['gene'].head(20))}",
        "",
        f"**Cell-type specific** genes (DE in exactly one cell type): {len(specific)}.",
        "",
        md_table(spec_by_type.to_frame()) if len(spec_by_type) else "",
    ]
    Path(ensure_parent(args.output)).write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
