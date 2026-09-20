import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import csv
import json

from comparability import config
from comparability.model import data, evaluate, invariant
import numpy as np

OUT = config.RESULTS_DIR
os.makedirs(OUT, exist_ok=True)

ds = data.load().restrict(min_studies=2, min_samples=3)
Z = data.clr(ds.X)
studies = ds.studies
factors = {f: ds.factor(f) for f in
           ("study_key", "sample_species", "sample_cultivar_canon", "sample_organ",
            "sample_processing", "study_column_polarity", "study_platform_primary")}

raw = data.impute_for_model(Z, "study_mean", studies)

Zc, usable = invariant.study_centre(Z, studies)
print("centreable samples: %d / %d (studies with >=2 samples)" % (usable.sum(), len(ds)))
cen = data.impute_for_model(Zc[usable], "study_mean", studies[usable])
adv, dirs = invariant.adversarial_projection(cen, studies[usable], n_directions=8)
print("adversarial directions removed: %d" % len(dirs))

reps = {
    "raw_clr": (raw, np.ones(len(ds), dtype=bool)),
    "study_centred": (cen, usable),
    "centred+adversarial": (adv, usable),
}

print("\n=== study-identity probe (excess over permuted labels; 0 = invariant) ===")
print("  %-22s %10s %10s %10s" % ("representation", "study", "polarity", "platform"))
probe_rows = []
for name, (M, keep) in reps.items():
    sub = {k: v[keep] for k, v in factors.items()}
    p = evaluate.nuisance_probes(M, sub)
    row = {"representation": name}
    cells = []
    for f in ("study_key", "study_column_polarity", "study_platform_primary"):
        e = p.get(f, {})
        val = e.get("excess") if e.get("testable") else None
        row[f] = val
        cells.append("%+10.3f" % val if val is not None else "%10s" % "n/a")
    print("  %-22s %s" % (name, "".join(cells)))
    probe_rows.append(row)


print("\n=== leave-study-out on the WITHIN-STUDY deviation target ===")
print("  (grand mean of a centred target is ~0 by construction, so the")
print("   reference is 'predict zero deviation': a factor must explain how a")
print("   sample departs from its own study-mates.)")

keep = usable
Zt = Zc[keep]
meta = [ds.meta[i] for i in np.where(keep)[0]]
st = studies[keep]
folds, skipped = evaluate.leave_study_out_folds(st)
print("\n  folds: %d" % len(folds))

results = []
for factor in ("sample_species", "sample_cultivar_canon", "sample_organ",
               "sample_processing", "sample_material_state"):
    ok, values = invariant.contrast_targets(Zt, meta, factor)
    n_studies_varying = len({meta[i]["study_key"] for i in np.where(ok)[0]})
    r2s = []
    for train, test, held in folds:
        tr = [i for i in train if ok[i]]
        te = [i for i in test if ok[i]]
        if len(tr) < 4 or len(te) < 2:
            continue
        lut = {}
        for lev in set(values[tr]):
            rows = [i for i in tr if values[i] == lev]
            m = np.nanmean(Zt[rows], axis=0)
            lut[lev] = np.where(np.isfinite(m), m, 0.0)
        pred = np.array([lut.get(values[i], np.zeros(Zt.shape[1])) for i in te])
        obs = np.isfinite(Zt[te])
        if obs.sum() == 0:
            continue
        r2, _ = evaluate.r2_rmse(Zt[te][obs], pred[obs], 0.0)
        if np.isfinite(r2):
            r2s.append(r2)
    results.append({
        "factor": factor,
        "studies_varying": n_studies_varying,
        "folds_scored": len(r2s),
        "mean_r2": round(float(np.mean(r2s)), 4) if r2s else float("nan"),
    })

print("\n  %-26s %9s %7s %9s" % ("factor", "varying", "folds", "mean R2"))
for r in results:
    v = r["mean_r2"]
    print("  %-26s %9d %7d %9s" %
          (r["factor"], r["studies_varying"], r["folds_scored"],
           ("%.3f" % v) if np.isfinite(v) else "untestable"))

with open(os.path.join(OUT, "invariant_results.tsv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(results[0].keys()), delimiter="\t")
    w.writeheader()
    for r in results: w.writerow(r)
json.dump(probe_rows, open(os.path.join(OUT, "invariant_probes.json"), "w"), indent=1)
for _p in ("invariant_results.tsv", "invariant_probes.json"):
    print("-> %s" % os.path.basename(os.path.join(OUT, _p)))
