import csv
import os

import numpy as np

from .. import config

MEASUREMENTS = config.data_path("measurements_long.tsv")


META_FIELDS = [
    "study_key", "doi", "pmcid", "tier",
    "sample_species", "sample_cultivar_canon", "sample_cultivar",
    "sample_organ", "sample_stage", "sample_origin", "sample_processing",
    "sample_material_state", "sample_treatment",
    "study_platform_primary", "study_column_polarity", "study_ri_basis_named",
    "column_class", "sample_column_header",
]


class Dataset:

    def __init__(self, X, mask, samples, compounds, meta):
        self.X = X
        self.mask = mask
        self.samples = samples
        self.compounds = compounds
        self.meta = meta

    def __len__(self):
        return len(self.samples)

    def factor(self, name):
        return np.array([m.get(name, "") or "" for m in self.meta])

    @property
    def studies(self):
        return self.factor("study_key")

    def describe(self):
        n, p = self.X.shape
        miss = float(np.isnan(self.X).mean())
        return {
            "samples": n,
            "compounds": p,
            "missing_fraction": round(miss, 4),
            "studies": len(set(self.studies)),
            "observed_cells": int((~np.isnan(self.X)).sum()),
        }

    def restrict(self, min_studies=2, min_samples=3):

        studies = self.studies
        keep = []
        for j in range(self.X.shape[1]):
            seen = self.mask[:, j]
            if seen.sum() < min_samples:
                continue
            if len(set(studies[seen])) < min_studies:
                continue
            keep.append(j)
        keep = np.array(keep, dtype=int)
        return Dataset(self.X[:, keep], self.mask[:, keep], self.samples,
                       [self.compounds[j] for j in keep], self.meta)


def load(measure_kind="relative_pct", in_scope_only=True, path=MEASUREMENTS):

    rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))

    cells = {}
    meta_by_sample = {}
    compounds = {}
    for r in rows:
        if in_scope_only and (r.get("in_composition_scope") or "yes") != "yes":
            continue
        if r.get("measure_kind") != measure_kind:
            continue
        skel = r.get("inchikey_skeleton") or ""
        if not skel:
            continue
        study_key = r.get("pmcid") or r.get("doi") or ""
        sample = (study_key, r.get("table_id", ""), r.get("column_index", ""))

        status = r.get("detection_status", "")
        value = None
        observed = False
        if status == "quantified":
            raw = r.get("area_pct_uncorrected") or r.get("value") or ""
            try:
                value = float(raw)
                observed = True
            except ValueError:
                continue
        elif status == "not_detected":
            value, observed = 0.0, True
        elif status == "detected_not_quantified":
            value, observed = None, True
        else:
            continue

        compounds.setdefault(skel, len(compounds))
        prev = cells.get((sample, skel))

        if prev is not None and prev[0] is not None and value is not None:
            value = prev[0] + value
        cells[(sample, skel)] = (value, observed)

        if sample not in meta_by_sample:
            m = {f: r.get(f, "") for f in META_FIELDS}
            m["study_key"] = study_key
            m["table_id"] = r.get("table_id", "")
            m["column_index"] = r.get("column_index", "")
            meta_by_sample[sample] = m

    samples = sorted(meta_by_sample)
    sidx = {s: i for i, s in enumerate(samples)}
    cidx = compounds
    X = np.full((len(samples), len(cidx)), np.nan)
    mask = np.zeros_like(X, dtype=bool)
    for (sample, skel), (value, observed) in cells.items():
        i, j = sidx[sample], cidx[skel]
        mask[i, j] = observed
        if value is not None:
            X[i, j] = value

    order = [k for k, _ in sorted(cidx.items(), key=lambda kv: kv[1])]
    meta = [meta_by_sample[s] for s in samples]
    return Dataset(X, mask, samples, order, meta)


def clr(X, delta_frac=0.65):

    Z = X.copy()
    for j in range(Z.shape[1]):
        col = Z[:, j]
        pos = col[np.isfinite(col) & (col > 0)]
        if pos.size == 0:
            continue
        floor = delta_frac * pos.min()
        zero = np.isfinite(col) & (col <= 0)
        col[zero] = floor
        Z[:, j] = col

    out = np.full_like(Z, np.nan)
    for i in range(Z.shape[0]):
        row = Z[i]
        obs = np.isfinite(row) & (row > 0)
        if obs.sum() < 2:
            continue
        gm = np.exp(np.log(row[obs]).mean())
        out[i, obs] = np.log(row[obs] / gm)
    return out


def impute_for_model(Z, how="study_mean", studies=None):

    F = Z.copy()
    if how == "zero":
        F[~np.isfinite(F)] = 0.0
        return F
    colmean = np.nanmean(np.where(np.isfinite(Z), Z, np.nan), axis=0)
    colmean = np.where(np.isfinite(colmean), colmean, 0.0)
    if how == "study_mean" and studies is not None:
        for s in set(studies):
            rows = np.where(studies == s)[0]
            sub = Z[rows]
            m = np.nanmean(np.where(np.isfinite(sub), sub, np.nan), axis=0)
            m = np.where(np.isfinite(m), m, colmean)
            for i in rows:
                bad = ~np.isfinite(F[i])
                F[i, bad] = m[bad]
        bad = ~np.isfinite(F)
        F[bad] = np.take(colmean, np.where(bad)[1])
        return F
    bad = ~np.isfinite(F)
    F[bad] = np.take(colmean, np.where(bad)[1])
    return F
