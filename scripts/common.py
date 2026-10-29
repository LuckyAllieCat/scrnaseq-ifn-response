"""Small helpers shared by the pipeline scripts (kept pure so they are easy to test)."""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless plotting for Snakemake / CI
import numpy as np
import pandas as pd
import yaml


def load_config(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def ensure_parent(path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    return path


def mad_bounds(values, n_mads):
    """Median +/- n_mads * MAD (scaled to match the SD for normal data)."""
    values = np.asarray(values, dtype=float)
    med = np.median(values)
    mad = 1.4826 * np.median(np.abs(values - med))
    return med - n_mads * mad, med + n_mads * mad


def choose_n_pcs(variance_ratio, min_pcs=5):
    """Elbow of the PCA scree curve: the point furthest from the line joining
    the first and last PC (a simple 'kneedle' rule)."""
    y = np.asarray(variance_ratio, dtype=float)
    x = np.arange(len(y), dtype=float)
    p1, p2 = np.array([x[0], y[0]]), np.array([x[-1], y[-1]])
    line = (p2 - p1) / np.linalg.norm(p2 - p1)
    vecs = np.stack([x, y], axis=1) - p1
    proj = np.outer(vecs @ line, line)
    dist = np.linalg.norm(vecs - proj, axis=1)
    return max(int(np.argmax(dist)) + 1, min_pcs)


def neighbor_mixing(connectivities, labels):
    """Mean normalised Shannon entropy of `labels` among each cell's kNN neighbours.
    1 = neighbours are perfectly mixed across labels, 0 = fully separated."""
    labels = pd.Categorical(labels)
    codes, k = labels.codes, len(labels.categories)
    if k < 2:
        return 1.0
    conn = connectivities.tocsr()
    onehot = np.eye(k)[codes]
    counts = conn @ onehot  # weighted label counts per cell's neighbourhood
    probs = counts / np.clip(counts.sum(axis=1, keepdims=True), 1e-12, None)
    with np.errstate(divide="ignore", invalid="ignore"):
        ent = -np.nansum(np.where(probs > 0, probs * np.log(probs), 0), axis=1)
    # compare against the entropy of the global label distribution (best achievable)
    glob = np.bincount(codes, minlength=k) / len(codes)
    glob_ent = -np.sum(glob[glob > 0] * np.log(glob[glob > 0]))
    return float(np.mean(ent) / glob_ent)


def pseudobulk(counts, obs, group_cols):
    """Sum raw counts over cells sharing the same values of `group_cols`.
    Returns (summed counts DataFrame [groups x genes], group metadata with n_cells)."""
    keys = obs[group_cols].astype(str).agg("|".join, axis=1)
    groups = pd.Categorical(keys)
    # sparse indicator matrix groups x cells, then a single matrix product
    from scipy import sparse

    ind = sparse.csr_matrix(
        (np.ones(len(keys)), (groups.codes, np.arange(len(keys)))),
        shape=(len(groups.categories), len(keys)),
    )
    summed = ind @ counts
    summed = summed.toarray() if sparse.issparse(summed) else np.asarray(summed)
    meta = obs.assign(_key=keys.values).groupby("_key", observed=True)[group_cols].first()
    meta = meta.loc[groups.categories]
    meta["n_cells"] = np.bincount(groups.codes)
    return pd.DataFrame(summed, index=groups.categories), meta
