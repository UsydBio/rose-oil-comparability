import numpy as np


def study_centre(Z, studies, min_samples=2):

    Z = np.asarray(Z, dtype=float)
    studies = np.asarray(studies)
    out = np.full_like(Z, np.nan)
    usable = np.zeros(len(Z), dtype=bool)

    for s in set(studies):
        rows = np.where(studies == s)[0]
        if len(rows) < min_samples:
            continue
        sub = Z[rows]
        with np.errstate(invalid="ignore"):
            centroid = np.nanmean(sub, axis=0)

        counts = np.isfinite(sub).sum(axis=0)
        centroid = np.where(counts >= 2, centroid, np.nan)
        out[rows] = sub - centroid
        usable[rows] = True
    return out, usable


def adversarial_projection(F, studies, n_directions=8, seed=0):

    from sklearn.linear_model import LogisticRegression

    X = np.asarray(F, dtype=float).copy()
    removed = []
    for _ in range(n_directions):
        y = np.asarray(studies)
        if len(set(y)) < 2:
            break
        clf = LogisticRegression(max_iter=1000, random_state=seed)
        try:
            clf.fit(X, y)
        except Exception:
            break
        W = np.atleast_2d(clf.coef_)

        _, sv, Vt = np.linalg.svd(W, full_matrices=False)
        if sv[0] < 1e-8:
            break
        v = Vt[0]
        v = v / (np.linalg.norm(v) + 1e-12)
        X = X - np.outer(X @ v, v)
        removed.append(v)
    return X, np.array(removed)


def contrast_targets(Zc, meta, factor):

    studies = np.array([m.get("study_key", "") for m in meta])
    values = np.array([m.get(factor, "") or "" for m in meta])
    keep = np.zeros(len(meta), dtype=bool)
    for s in set(studies):
        rows = np.where(studies == s)[0]
        levels = {v for v in values[rows] if v}
        if len(levels) >= 2:
            keep[rows] = True
    return keep, values
