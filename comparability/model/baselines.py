import numpy as np


def _group_mean(Z, groups, target_groups, fallback):

    out = np.tile(fallback, (len(target_groups), 1))
    lut = {}
    for g in sorted(set(groups)):
        if not str(g).strip():
            continue
        rows = np.where(groups == g)[0]
        if len(rows) == 0:
            continue
        m = np.nanmean(Z[rows], axis=0)
        lut[g] = np.where(np.isfinite(m), m, fallback)
    hits = 0
    for i, g in enumerate(target_groups):
        if str(g).strip() and g in lut:
            out[i] = lut[g]
            hits += 1
    return out, hits


def grand_mean(Ztr, meta_tr, Zte, meta_te, **kw):
    m = np.nanmean(Ztr, axis=0)
    m = np.where(np.isfinite(m), m, 0.0)
    return np.tile(m, (Zte.shape[0], 1)), 0


def _factor_mean(field):
    def f(Ztr, meta_tr, Zte, meta_te, **kw):
        fallback = np.nanmean(Ztr, axis=0)
        fallback = np.where(np.isfinite(fallback), fallback, 0.0)
        g_tr = np.array([m.get(field, "") or "" for m in meta_tr])
        g_te = np.array([m.get(field, "") or "" for m in meta_te])
        return _group_mean(Ztr, g_tr, g_te, fallback)
    return f


species_mean = _factor_mean("sample_species")
cultivar_mean = _factor_mean("sample_cultivar_canon")
platform_mean = _factor_mean("study_platform_primary")
polarity_mean = _factor_mean("study_column_polarity")
organ_mean = _factor_mean("sample_organ")
material_mean = _factor_mean("sample_material_state")


def pca_ridge(Ztr, meta_tr, Zte, meta_te, n_components=5, **kw):

    from sklearn.decomposition import PCA
    from ..model.data import impute_for_model

    studies_tr = np.array([m.get("study_key", "") for m in meta_tr])
    Ftr = impute_for_model(Ztr, "study_mean", studies_tr)
    k = int(min(n_components, Ftr.shape[0] - 1, Ftr.shape[1]))
    if k < 1:
        return grand_mean(Ztr, meta_tr, Zte, meta_te)
    p = PCA(n_components=k, random_state=0).fit(Ftr)
    mean = p.mean_
    return np.tile(mean, (Zte.shape[0], 1)), 0


BASELINES = {
    "grand_mean": grand_mean,
    "species_mean": species_mean,
    "cultivar_mean": cultivar_mean,
    "organ_mean": organ_mean,
    "material_state_mean": material_mean,
    "platform_mean": platform_mean,
    "polarity_mean": polarity_mean,
    "pca_ridge": pca_ridge,
}
