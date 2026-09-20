import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import csv
import math

from comparability import config
from comparability.model import baselines, data, evaluate
import numpy as np

OUT = config.RESULTS_DIR
PER_FOLD = os.path.join(OUT, "baseline_per_fold_r2.tsv")
SUMMARY = os.path.join(OUT, "baseline_fold_tests.tsv")

ds = data.load().restrict(min_studies=2, min_samples=3)
Z = data.clr(ds.X)
studies = ds.studies
folds, skipped = evaluate.leave_study_out_folds(studies)
print("matrix: %s" % ds.describe())
print("folds: %d" % len(folds))


per_fold = {}
fold_names, fold_n = [], []
for bname, fn in baselines.BASELINES.items():
    row = {}
    for train, test, held in folds:
        Ztr, Zte = Z[train], Z[test]
        pred, _hit = fn(Ztr, [ds.meta[i] for i in train],
                        Zte, [ds.meta[i] for i in test])
        obs = np.isfinite(Zte)
        if obs.sum() == 0:
            continue
        r2, _ = evaluate.r2_rmse(Zte[obs], pred[obs], np.nanmean(Ztr))
        if np.isfinite(r2):
            row[held] = float(r2)
    per_fold[bname] = row
    if not fold_names:
        fold_names = [h for _t, _e, h in folds if h in row]
        fold_n = {h: int(len(e)) for _t, e, h in folds}

names = list(per_fold)
with open(PER_FOLD, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["held_out_study", "n_test_samples"] + names)
    for h in fold_names:
        w.writerow([h, fold_n.get(h, "")] +
                   ["%.6f" % per_fold[b][h] if h in per_fold[b] else ""
                    for b in names])
print("wrote %s (%d folds x %d baselines)"
      % (os.path.basename(PER_FOLD), len(fold_names), len(names)))


def sign_test(wins, losses):

    n = wins + losses
    if n == 0:
        return float("nan")
    k = min(wins, losses)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / (2.0 ** n)
    return min(1.0, 2.0 * tail)


ref = per_fold["grand_mean"]
rows = []
for b in names:
    if b == "grand_mean":
        continue
    shared = [h for h in fold_names if h in per_fold[b] and h in ref]
    diffs = [per_fold[b][h] - ref[h] for h in shared]
    wins = sum(1 for d in diffs if d > 1e-12)
    losses = sum(1 for d in diffs if d < -1e-12)
    ties = len(diffs) - wins - losses
    rows.append({
        "baseline": b, "folds": len(shared),
        "beats_grand_mean": wins, "loses_to_grand_mean": losses, "ties": ties,
        "mean_delta_r2": round(float(np.mean(diffs)), 4) if diffs else "",
        "median_delta_r2": round(float(np.median(diffs)), 4) if diffs else "",
        "min_delta_r2": round(float(np.min(diffs)), 4) if diffs else "",
        "max_delta_r2": round(float(np.max(diffs)), 4) if diffs else "",
        "sign_test_p": round(sign_test(wins, losses), 4),
    })
rows.sort(key=lambda r: -r["beats_grand_mean"])
with open(SUMMARY, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    for r in rows:
        w.writerow(r)

print("\n=== per-fold: does the baseline beat the training grand mean? ===")
print("  %-22s %6s %6s %6s %6s %10s %9s"
      % ("baseline", "folds", "wins", "loses", "ties", "median dR2", "sign p"))
for r in rows:
    print("  %-22s %6d %6d %6d %6d %10.4f %9.4f"
          % (r["baseline"], r["folds"], r["beats_grand_mean"],
             r["loses_to_grand_mean"], r["ties"], r["median_delta_r2"],
             r["sign_test_p"]))

pol = next(r for r in rows if r["baseline"] == "polarity_mean")
for r in rows:
    print("%-20s folds %2d  beats %2d  loses %2d  ties %2d  median %+.4f  "
          "range %+.4f to %+.4f  sign p %.4f"
          % (r["baseline"], r["folds"], r["beats_grand_mean"],
             r["loses_to_grand_mean"], r["ties"], r["median_delta_r2"],
             r["min_delta_r2"], r["max_delta_r2"], r["sign_test_p"]))
print("\npolarity_mean: %d of %d folds beat the grand mean, %d lose, %d tie; "
      "exact two-sided sign test p = %.4f"
      % (pol["beats_grand_mean"], pol["folds"], pol["loses_to_grand_mean"],
         pol["ties"], pol["sign_test_p"]))
print("wrote %s and %s" % (os.path.basename(PER_FOLD), os.path.basename(SUMMARY)))
