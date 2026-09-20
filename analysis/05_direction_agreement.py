import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import csv
import json
import warnings


warnings.filterwarnings("ignore", category=FutureWarning, module="sklearn")

from comparability import config
from comparability import contrast_analysis as CA
from comparability.model import direction as D


OUT = config.RESULTS_DIR
os.makedirs(OUT, exist_ok=True)

EFFECTS_TSV = os.path.join(OUT, "direction_effects.tsv")
PAIRS_TSV = os.path.join(OUT, "direction_cross_study_pairs.tsv")
MODEL_TSV = os.path.join(OUT, "direction_model.tsv")
SUMMARY_TSV = os.path.join(OUT, "direction_summary.tsv")
THRESH_TSV = os.path.join(OUT, "direction_effect_thresholds.tsv")
PROBES_JSON = os.path.join(OUT, "direction_probes.json")

EFFECT_FIELDS = ["contrast_id", "study_id", "axis", "measure_kind", "orient_key",
                 "tier", "level_a", "level_b", "compound_key", "compound_name",
                 "effect_clr", "effect_raw_pct", "sign", "sign_raw_pct",
                 "abs_effect_clr", "n_basis", "basis_preselected", "caveat"]

PAIR_FIELDS = ["axis", "orient_key", "tier", "relation", "basis", "status", "study_a",
               "contrast_a", "level_pair_a", "study_b", "contrast_b",
               "level_pair_b", "n_shared_inchikey", "n_agree", "n_disagree",
               "n_tie", "frac_agree", "n_zero_backed", "n_agree_quantified",
               "n_disagree_quantified", "frac_agree_quantified",
               "binomial_p_quantified", "n_agree_raw_pct", "n_disagree_raw_pct",
               "binomial_p", "binomial_q", "shared_compound_signs"]


def write_tsv(path, fields, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: CA.fmt(r.get(k)) if isinstance(r.get(k), float)
                        else r.get(k, "") for k in fields})


def main():
    print("[1/7] rebuilding every contrast's CLR block (class-total rows dropped) ...")
    blocks, names, info = D.load_blocks(drop_aggregate_rows=True)
    corpus = info["_corpus"]
    print("      %d measurement columns; %d genuine contrasts, %d after "
          "dropping superseded parents; %d with a usable basis"
          % (corpus["n_sample_columns"], corpus["n_genuine"],
             corpus["n_specs"] - corpus["n_superseded_parents"], len(blocks)))
    n_agg = sum(info[c].get("n_aggregate_rows", 0) for c in blocks)
    agg_hits = sorted((info[c]["n_aggregate_rows"], c) for c in blocks
                      if info[c].get("n_aggregate_rows"))
    print("      class-total rows removed from bases: %d across %d contrasts"
          % (n_agg, len(agg_hits)))

    print("[2/7] signed per-compound effects on oriented level pairs ...")
    rows = D.effect_rows(blocks, names)
    studies = sorted({r["study_id"] for r in rows})
    contrasts = sorted({r["contrast_id"] for r in rows})
    pos = sum(1 for r in rows if r["sign"] > 0)
    print("      %d (contrast x pair x compound) effects; %d contrasts, %d studies"
          % (len(rows), len(contrasts), len(studies)))
    print("      sign prior: %.4f positive -- CLR effects sum to zero across a "
          "basis, so 0.5 really is the null" % (pos / float(len(rows))))
    write_tsv(EFFECTS_TSV, EFFECT_FIELDS, rows)

    print("[3/7] cross-study pair inventory, per axis ...")
    cross = D.cross_study_pairs(blocks, same_study=False)
    within = D.cross_study_pairs(blocks, same_study=True)

    shared_basis = D.cross_study_pairs(blocks, same_study=False, subcompositional=True)
    for p in cross:
        p["relation"] = "cross_study"
    for p in within:
        p["relation"] = "same_study_control"
    for p in shared_basis:
        p["relation"] = "cross_study_shared_basis"

    families = {
        "strict_cross_study": [p for p in cross if p["tier"] == "strict" and p["status"] == "ok"],
        "oriented_cross_study": [p for p in cross if p["tier"] == "oriented" and p["status"] == "ok"],
        "same_study_control": [p for p in within if p["status"] == "ok"],
        "shared_basis_cross_study": [p for p in shared_basis
                                     if p["tier"] == "strict" and p["status"] == "ok"],
    }
    for fam, ps in families.items():
        qs = CA.benjamini_hochberg([p["binomial_p"] for p in ps])
        for p, q in zip(ps, qs):
            p["binomial_q"] = q
            p["family"] = fam

    axis_rows = axis_inventory(blocks, cross)
    print("      %-16s %8s %8s %8s %10s %10s %10s"
          % ("axis", "contr", "studies", "pairs", "xstudy", "usable", "indep"))
    for a in axis_rows:
        print("      %-16s %8d %8d %8d %10d %10d %10d"
              % (a["axis"], a["n_contrasts"], a["n_studies"], a["n_oriented_pairs"],
                 a["n_cross_study_pairs"], a["n_usable"], a["n_independent_study_pairs"]))

    allpairs = cross + within + shared_basis
    for p in allpairs:
        p["shared_compound_signs"] = "; ".join(
            "%s(%s%s)" % (names.get(c, c), "+" if x > 0 else "-", "+" if y > 0 else "-")
            for c, x, y in p["detail"][:30])
    write_tsv(PAIRS_TSV, PAIR_FIELDS, allpairs)

    thresh = D.agreement_by_effect_threshold(cross + within)
    write_tsv(THRESH_TSV,
              ["relation", "tier", "axis", "orient_key", "study_a", "study_b",
               "contrast_a", "contrast_b", "min_abs_effect_clr", "n_agree",
               "n_disagree", "n_compared", "frac_agree", "binomial_p"], thresh)
    print("      agreement once both sides had to move by at least |CLR| = t:")
    print("      %8s %28s %28s" % ("t", "strict cross-study", "same-lab control"))
    for t in (0.0, 0.25, 0.5, 1.0):
        def tally(sel):
            a = sum(r["n_agree"] for r in thresh
                    if r["min_abs_effect_clr"] == t and sel(r))
            d = sum(r["n_disagree"] for r in thresh
                    if r["min_abs_effect_clr"] == t and sel(r))
            return "%d/%d%s" % (a, a + d,
                                (" = %.3f" % (a / float(a + d))) if a + d else "")
        print("      %8.2f %28s %28s"
              % (t,
                 tally(lambda r: r["relation"] == "cross_study" and r["tier"] == "strict"),
                 tally(lambda r: r["relation"] == "same_study_control")))

    print("[4/7] compound structures (SMILES -> coarse descriptors) ...")
    keys = sorted({r["compound_key"] for r in rows})
    structures = D.compound_structures(keys)
    have = sum(1 for k in keys if structures.get(k, {}).get("has_smiles"))
    print("      %d distinct compounds, %d with a PubChem SMILES (%.1f%%)"
          % (len(keys), have, 100.0 * have / max(len(keys), 1)))

    print("[5/7] leave-study-out sign prediction ...")
    folds, skipped = D.leave_study_out(rows, structures)
    summary = D.summarise_folds(folds)
    qs = CA.benjamini_hochberg([f["binomial_p"] for f in folds])
    for f, q in zip(folds, qs):
        f["binomial_q"] = q
    print("      folds: %d used, %d skipped" % (len({f["held_out_study"] for f in folds}),
                                                len(skipped)))
    print("      %-20s %7s %10s %10s %10s %10s"
          % ("model", "folds", "macro acc", "pooled", "coverage", "shared-sub"))
    for s in summary:
        print("      %-20s %7d %10.4f %10s %10.4f %10s"
              % (s["model"], s["n_folds"], s["macro_accuracy"],
                 ("%.4f" % s["pooled_accuracy"]) if s["pooled_accuracy"] is not None else "-",
                 s["mean_coverage"],
                 ("%.4f" % s["macro_accuracy_shared_subset"])
                 if s["macro_accuracy_shared_subset"] is not None else "-"))

    bins = D.accuracy_by_effect_size(rows)
    bq = CA.benjamini_hochberg([b["binomial_p"] for b in bins])
    for b, q in zip(bins, bq):
        b["binomial_q"] = q
    print("      by |CLR effect| quartile (compound-majority rule):")
    for b in bins:
        print("        Q%d  |d| in [%.2f, %.2f]  n=%4d  acc=%s"
              % (b["bin"], b["abs_clr_low"], b["abs_clr_high"], b["n_answered"],
                 ("%.4f" % b["accuracy"]) if b["accuracy"] is not None else "-"))

    fold_fields = ["held_out_study", "model", "n_test", "n_answered", "coverage",
                   "accuracy", "n_correct", "accuracy_on_shared_subset",
                   "n_shared_subset", "binomial_p", "binomial_q", "test_axes",
                   "axes_seen_in_training"]
    write_tsv(MODEL_TSV, fold_fields, folds)
    sum_fields = ["model", "n_folds", "macro_accuracy", "macro_accuracy_sd",
                  "pooled_accuracy", "n_answered", "mean_coverage",
                  "macro_accuracy_shared_subset", "pooled_binomial_p"]
    write_tsv(SUMMARY_TSV, sum_fields, summary)

    print("[6/7] study-identity probe on each feature space ...")
    probes = D.probe_representations(rows, structures)
    print("      %-18s %8s %10s %10s %10s"
          % ("representation", "feats", "bal.acc", "permuted", "excess"))
    for mode in ("struct", "axis_onehot", "struct+axis", "struct+ngram",
                 "compound_onehot"):
        p = probes[mode]
        if not p.get("testable"):
            print("      %-18s %8s not testable: %s" % (mode, "-", p.get("reason")))
            continue
        print("      %-18s %8d %10.4f %10.4f %+10.4f"
              % (mode, p["n_features"], p["balanced_accuracy"],
                 p["permuted_baseline"], p["excess"]))
    json.dump(probes, open(PROBES_JSON, "w"), indent=1)

    print("[7/7] sensitivity: keeping the class-total rows in the basis ...")
    blocks2, names2, info2 = D.load_blocks(drop_aggregate_rows=False)
    rows2 = D.effect_rows(blocks2, names2)
    cross2 = D.cross_study_pairs(blocks2, same_study=False)
    sens = compare_sensitivity(cross, cross2)
    for s in sens:
        print("      %-16s %-24s %s" % (s["axis"], s["orient_key"], s["change"]))

    print()
    print("wrote:")
    for p in (EFFECTS_TSV, PAIRS_TSV, THRESH_TSV, MODEL_TSV, SUMMARY_TSV,
              PROBES_JSON):
        print("  ", os.path.basename(p))


def axis_inventory(blocks, cross):

    rows = []
    axes = sorted({b["axis"] for b in blocks.values()} | set(D.UNORIENTABLE_AXES))
    for axis in axes:
        bs = [b for b in blocks.values() if b["axis"] == axis]
        pairs = sum(len(D.oriented_pairs(b)) for b in bs)
        cs = [p for p in cross if p["axis"] == axis]
        ok = [p for p in cs if p["status"] == "ok"]
        indep = {tuple(sorted((p["study_a"], p["study_b"]))) for p in ok}
        rows.append({
            "axis": axis, "n_contrasts": len(bs),
            "n_studies": len({b["study"] for b in bs}),
            "n_oriented_pairs": pairs,
            "n_cross_study_pairs": len(cs),
            "n_usable": len(ok),
            "n_too_few_shared": len(cs) - len(ok),
            "n_independent_study_pairs": len(indep),
            "total_shared_compounds": sum(p["n_shared_inchikey"] for p in ok),
            "reason_if_zero": D.UNORIENTABLE_AXES.get(axis, "")
                              if not pairs else ("" if cs else
                              "oriented pairs exist but no two studies share one"),
        })
    return rows


def compare_sensitivity(cross_dropped, cross_kept):

    key = lambda p: (p["axis"], p["orient_key"], p["study_a"], p["study_b"])
    a = {key(p): p for p in cross_dropped if p["status"] == "ok"}
    b = {key(p): p for p in cross_kept if p["status"] == "ok"}
    out = []
    for k in sorted(set(a) | set(b)):
        pa, pb = a.get(k), b.get(k)
        if pa is None:
            out.append({"axis": k[0], "orient_key": k[1],
                        "change": "only present when class-totals are kept "
                                  "(%d/%d agree)" % (pb["n_agree"], pb["n_agree"] + pb["n_disagree"])})
        elif pb is None:
            out.append({"axis": k[0], "orient_key": k[1],
                        "change": "only present when class-totals are dropped"})
        elif (pa["n_agree"], pa["n_disagree"]) != (pb["n_agree"], pb["n_disagree"]):
            out.append({"axis": k[0], "orient_key": k[1],
                        "change": "dropped %d/%d agree -> kept %d/%d agree"
                                  % (pa["n_agree"], pa["n_agree"] + pa["n_disagree"],
                                     pb["n_agree"], pb["n_agree"] + pb["n_disagree"])})
        else:
            out.append({"axis": k[0], "orient_key": k[1],
                        "change": "unchanged (%d/%d agree)"
                                  % (pa["n_agree"], pa["n_agree"] + pa["n_disagree"])})
    return out


if __name__ == "__main__":
    main()
