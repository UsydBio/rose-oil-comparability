import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import csv
import json


from comparability import config
from comparability.model import baselines, data, evaluate
import numpy as np

OUT = config.RESULTS_DIR
os.makedirs(OUT, exist_ok=True)

ds = data.load().restrict(min_studies=2, min_samples=3)
print("matrix:", ds.describe())

Z = data.clr(ds.X)
studies = ds.studies
folds, skipped = evaluate.leave_study_out_folds(studies)
print("usable leave-study-out folds: %d (skipped %d)" % (len(folds), len(skipped)))
for s, n, why in skipped[:6]:
    print("   skip %-14s n=%-3d %s" % (s[:14], n, why))


factors = {f: ds.factor(f) for f in
           ("study_key", "sample_species", "sample_cultivar_canon", "sample_organ",
            "sample_stage", "sample_origin", "sample_processing",
            "study_column_polarity", "study_platform_primary")}
print("\nsplit feasibility:")
for name, why in evaluate.impossible_splits(factors):
    print("   %-26s %s" % (name, why))

results = []
for bname, fn in baselines.BASELINES.items():
    r2s, rmses, hits, cells = [], [], 0, 0
    for train, test, held in folds:
        Ztr, Zte = Z[train], Z[test]
        mtr = [ds.meta[i] for i in train]
        mte = [ds.meta[i] for i in test]
        pred, hit = fn(Ztr, mtr, Zte, mte)
        hits += hit
        train_mean = np.nanmean(Ztr)
        obs = np.isfinite(Zte)
        if obs.sum() == 0:
            continue
        r2, rmse = evaluate.r2_rmse(Zte[obs], pred[obs], train_mean)
        if np.isfinite(r2):
            r2s.append(r2); rmses.append(rmse); cells += int(obs.sum())
    results.append({
        "baseline": bname,
        "folds_scored": len(r2s),
        "cells": cells,
        "mean_r2": round(float(np.mean(r2s)), 4) if r2s else float("nan"),
        "median_r2": round(float(np.median(r2s)), 4) if r2s else float("nan"),
        "mean_rmse": round(float(np.mean(rmses)), 4) if rmses else float("nan"),
        "level_hits": hits,
    })

print("\n=== leave-study-out leaderboard (R2 vs training grand mean) ===")
print("  %-22s %6s %8s %9s %8s" % ("baseline", "folds", "mean R2", "median R2", "hits"))
for r in sorted(results, key=lambda x: -(x["mean_r2"] if np.isfinite(x["mean_r2"]) else -9)):
    print("  %-22s %6d %8.3f %9.3f %8d" %
          (r["baseline"], r["folds_scored"], r["mean_r2"], r["median_r2"], r["level_hits"]))


F = data.impute_for_model(Z, "study_mean", studies)
probes = evaluate.nuisance_probes(F, factors)
print("\n=== study-identity probe on the raw CLR feature space ===")
for name, p in probes.items():
    if not p.get("testable"):
        print("  %-26s untestable: %s" % (name, p.get("reason"))); continue
    print("  %-26s bal.acc %.3f | permuted %.3f | excess %+.3f  (n=%d, %d classes)" %
          (name, p["balanced_accuracy"], p["permuted_baseline"], p["excess"],
           p["n"], p["n_studies"]))

with open(os.path.join(OUT, "baseline_results.tsv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(results[0].keys()), delimiter="\t")
    w.writeheader()
    for r in results: w.writerow(r)
json.dump({"probes": probes, "matrix": ds.describe(),
           "folds": len(folds), "skipped": [list(s) for s in skipped]},
          open(os.path.join(OUT, "probe_results.json"), "w"), indent=1)
for _p in ("baseline_results.tsv", "probe_results.json"):
    print("-> %s" % os.path.basename(os.path.join(OUT, _p)))
