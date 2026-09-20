import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.preprocessing import StandardScaler

SEED = 20260913


def leave_study_out_folds(studies, min_test=3, min_train_studies=3):

    studies = np.asarray(studies)
    folds, skipped = [], []
    for s in sorted(set(studies)):
        test = np.where(studies == s)[0]
        train = np.where(studies != s)[0]
        if len(test) < min_test:
            skipped.append((s, len(test), "held-out study has < %d samples" % min_test))
            continue
        if len(set(studies[train])) < min_train_studies:
            skipped.append((s, len(test), "fewer than %d training studies" % min_train_studies))
            continue
        folds.append((train, test, s))
    return folds, skipped


def impossible_splits(meta_factors):

    out = []
    studies = meta_factors["study_key"]
    for name, values in meta_factors.items():
        if name == "study_key":
            continue
        values = np.asarray(values)
        multi = 0
        for s in set(studies):
            levels = {v for v in values[np.asarray(studies) == s] if v}
            if len(levels) >= 2:
                multi += 1
        if multi == 0:
            out.append((name, "no study contains two levels: collinear with study"))
        else:
            out.append((name, "varies inside %d studies" % multi))
    return out


def study_identity_probe(embeddings, studies, seed=SEED, min_per_class=2):

    X = np.asarray(embeddings, dtype=float)
    y = np.asarray(studies)
    keep = np.isfinite(X).all(axis=1)
    X, y = X[keep], y[keep]

    counts = {s: int((y == s).sum()) for s in sorted(set(y))}
    usable = np.array([counts[s] >= min_per_class for s in y])
    X, y = X[usable], y[usable]
    if len(set(y)) < 2 or len(y) < 10:
        return {"testable": False, "reason": "too few studies or samples"}

    Xs = StandardScaler().fit_transform(X)

    def fit_score(labels, rng=None):

        tr, te = [], []

        for s in sorted(set(labels)):
            idx = np.where(labels == s)[0]
            if rng is not None:
                rng.shuffle(idx)
            cut = max(1, len(idx) // 2)
            tr.extend(idx[:cut])
            te.extend(idx[cut:] if len(idx) > 1 else idx[:cut])
        tr, te = np.array(tr), np.array(te)
        clf = LogisticRegression(max_iter=2000)
        clf.fit(Xs[tr], labels[tr])
        return balanced_accuracy_score(labels[te], clf.predict(Xs[te]))

    rng = np.random.default_rng(seed)
    observed = fit_score(y, rng)
    perm = []
    for _ in range(5):
        shuffled = y.copy()
        rng.shuffle(shuffled)
        perm.append(fit_score(shuffled, rng))
    majority = 1.0 / len(set(y))

    return {
        "testable": True,
        "n": int(len(y)),
        "n_studies": int(len(set(y))),
        "balanced_accuracy": round(float(observed), 4),
        "permuted_baseline": round(float(np.mean(perm)), 4),
        "majority_baseline": round(float(majority), 4),
        "excess": round(float(observed - np.mean(perm)), 4),
    }


def nuisance_probes(embeddings, meta_factors, seed=SEED):

    out = {}
    for name in ("study_key", "study_platform_primary", "study_column_polarity"):
        values = meta_factors.get(name)
        if values is None:
            continue
        out[name] = study_identity_probe(embeddings, values, seed=seed)
    return out


def r2_rmse(y_true, y_pred, y_train_mean):

    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    ok = np.isfinite(y_true) & np.isfinite(y_pred)
    if ok.sum() == 0:
        return float("nan"), float("nan")
    resid = y_true[ok] - y_pred[ok]
    base = y_true[ok] - y_train_mean
    sse, sst = float((resid ** 2).sum()), float((base ** 2).sum())
    r2 = 1.0 - sse / sst if sst > 0 else float("nan")
    return r2, float(np.sqrt((resid ** 2).mean()))
