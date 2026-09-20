import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import csv
import math
from collections import Counter, defaultdict

from comparability import config, stats, separability

IN_TSV = config.data_path("measurements_long.tsv")
OUT_MATRIX_REL = config.results_path("sample_compound_matrix.tsv")
OUT_MATRIX_ABS = config.results_path("sample_compound_matrix_absolute.tsv")
OUT_TESTS = config.results_path("variance_partition_tests.tsv")
OUT_SEP = config.results_path("variance_separability.tsv")

N_PERM = 9999
SEED = 20260912
DELTA_FRAC = 0.65
DELTA_SENSITIVITY = [0.25, 0.65, 1.00]
ZERO_MASS_CAP = 0.30

UNSTATED = "(unstated)"

FACTORS = [
    ("study", "study_id", "analytical"),
    ("table", "__table__", "analytical"),
    ("column_polarity", "study_column_polarity", "analytical"),
    ("platform", "study_platform_primary", "analytical"),
    ("tier", "tier", "provenance"),
    ("species", "sample_species", "biological"),
    ("cultivar", "sample_cultivar", "biological"),
    ("organ", "sample_organ", "biological"),
    ("stage", "sample_stage", "biological"),
    ("origin", "sample_origin", "biological"),
    ("country", "study_country", "biological"),
    ("processing", "sample_processing", "processing"),
]
FACTOR_COL = {n: c for n, c, _ in FACTORS}
FACTOR_CLASS = {n: k for n, _, k in FACTORS}
BIO_FACTORS = [n for n, _, k in FACTORS if k in ("biological", "processing")]


META_COLS = ["sample_id", "study_id", "pmcid", "doi", "table_id", "column_index",
             "sample_column_header", "sample_species", "sample_cultivar",
             "sample_organ", "sample_processing", "sample_stage", "sample_origin",
             "study_country", "study_column_polarity", "study_platform_primary", "tier",
             "measure_kind", "units_seen",
             "n_quantified", "n_not_detected", "n_trace", "n_not_measured",
             "n_measured_total", "value_sum"]


def study_id_of(r):

    if r.get("pmcid"):
        return r["pmcid"]
    if r.get("doi"):
        return "doi:" + r["doi"]
    return "src:" + (r.get("table_id") or "").split("#")[0]


def load(measure_kind):
    csv.field_size_limit(10 ** 9)
    with open(IN_TSV, encoding="utf-8") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t")
                if r["measure_kind"] == measure_kind
                and r["in_composition_scope"] == "yes"]


def build_cells(rows, value_field):

    raw = defaultdict(list)
    meta = {}
    names = {}
    units = defaultdict(set)
    dropped_no_key = 0
    for r in rows:
        s = (study_id_of(r), r["table_id"], r["column_index"])
        meta.setdefault(s, r)
        if r.get("unit"):
            units[s].add(r["unit"])
        sk = r["inchikey_skeleton"]
        if not sk:
            dropped_no_key += 1
            continue
        raw[(s, sk)].append(r)
        nm = r["compound_clean_name"] or r["compound_base_name"] or r["compound_name_raw"]
        if sk not in names and nm:
            names[sk] = nm
    cells = {}
    multi = 0
    for key, rs in raw.items():
        if len(rs) > 1:
            multi += 1
        vals = []
        st = set()
        for r in rs:
            st.add(r["detection_status"])
            if r["detection_status"] == "quantified" and r[value_field]:
                try:
                    vals.append(float(r[value_field]))
                except ValueError:
                    pass
        if vals:
            cells[key] = ("quantified", sum(vals))
        elif "detected_not_quantified" in st:
            cells[key] = ("trace", None)
        elif "not_detected" in st:
            cells[key] = ("not_detected", 0.0)
        else:
            cells[key] = ("not_measured", None)
    return cells, meta, names, units, multi, dropped_no_key


def sample_summary(cells, samples, skeletons):
    out = {}
    for s in samples:
        c = Counter()
        tot = 0.0
        for sk in skeletons:
            v = cells.get((s, sk))
            if v is None:
                c["not_measured"] += 1
            else:
                c[v[0]] += 1
                if v[0] == "quantified":
                    tot += v[1]
        out[s] = (c, tot)
    return out


def measured_sets(cells):

    m = defaultdict(set)
    for (s, sk), (st, _) in cells.items():
        if st in ("quantified", "not_detected"):
            m[s].add(sk)
    return m


def write_matrix(path, cells, samples, skeletons, meta, names, summ, units, mk):
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(META_COLS + ["%s|%s" % (sk, names.get(sk, "")) for sk in skeletons])
        for s in samples:
            c, tot = summ[s]
            r = meta[s]
            w.writerow(["%s::%s::%s" % s, s[0], r.get("pmcid", ""), r.get("doi", ""),
                        r["table_id"], r["column_index"], r["sample_column_header"],
                        r.get("sample_species", ""), r.get("sample_cultivar", ""),
                        r.get("sample_organ", ""), r.get("sample_processing", ""),
                        r.get("sample_stage", ""), r.get("sample_origin", ""),
                        r.get("study_country", ""), r.get("study_column_polarity", ""),
                        r.get("study_platform_primary", ""), r.get("tier", ""),
                        mk, ";".join(sorted(units.get(s, []))) or "(none)",
                        c["quantified"], c["not_detected"], c["trace"], c["not_measured"],
                        c["quantified"] + c["not_detected"], "%.4f" % tot] +
                       [("NA" if (cells.get((s, sk)) is None
                                  or cells[(s, sk)][0] == "not_measured")
                         else ("TR" if cells[(s, sk)][0] == "trace"
                               else "%.6g" % cells[(s, sk)][1]))
                        for sk in skeletons])


def ladder(meas, max_m=30):
    per = defaultdict(set)
    for s, cs in meas.items():
        for c in cs:
            per[c].add(s)
    order = sorted(per, key=lambda c: (-len(per[c]), c))
    out = []
    for m in range(1, min(max_m, len(order)) + 1):
        pre = order[:m]
        ss = set(meas)
        for c in pre:
            ss &= per[c]
        out.append({"m": m, "compounds": list(pre), "samples": sorted(ss),
                    "n_studies": len(set(x[0] for x in ss))})
    return out, order, per


def pick_blocks(lad):

    picks = {}
    c1 = [b for b in lad if b["m"] >= 3]
    if c1:
        picks["B1"] = max(c1, key=lambda b: (len(b["samples"]), b["m"]))
    c2 = [b for b in lad if b["n_studies"] >= 8]
    if c2:
        picks["B2"] = max(c2, key=lambda b: (b["m"] * len(b["samples"]), b["m"]))
    c3 = [b for b in lad if len(b["samples"]) >= 15]
    if c3:
        picks["B3"] = max(c3, key=lambda b: (b["m"], len(b["samples"])))
    return picks


def block_vectors(cells, samples, comps, delta_frac=DELTA_FRAC):

    raw, keep = [], []
    for s in samples:
        v, ok = [], True
        for c in comps:
            cell = cells.get((s, c))
            if cell is None or cell[1] is None:
                ok = False
                break
            v.append(cell[1])
        if ok and sum(v) > 0:
            keep.append(s)
            raw.append(v)
    allpos = [x for row in raw for x in row if x > 0]
    fallback = min(allpos) if allpos else 1e-3
    deltas = []
    for j in range(len(comps)):
        pos = [row[j] for row in raw if row[j] > 0]
        deltas.append(delta_frac * (min(pos) if pos else fallback))
    closed, clrs, logs = [], [], []
    n_capped = 0
    for row in raw:
        cl = stats.closure(row, 100.0)
        closed.append(cl)
        dd = [d * 100.0 / sum(row) for d in deltas]
        zm = sum(d for v, d in zip(cl, dd) if v <= 0)
        if zm > ZERO_MASS_CAP * 100.0:
            f = ZERO_MASS_CAP * 100.0 / zm
            dd = [d * f for d in dd]
            n_capped += 1
        clrs.append(stats.clr(stats.multiplicative_replacement(cl, dd)))
        logs.append(separability.log_replace(row, deltas))
    return keep, raw, closed, clrs, logs, n_capped


def factor_values(meta, samples, name, complete_case=False):
    col = FACTOR_COL[name]
    vals = []
    for s in samples:
        if col == "__table__":
            v = "%s::%s" % (s[0], s[1])
        else:
            v = (meta[s].get(col) or "").strip() or UNSTATED
        vals.append(v)
    if complete_case:
        return vals, [i for i, v in enumerate(vals) if v != UNSTATED]
    return vals, list(range(len(samples)))


def run_factor_tests(tag, d2, samples, meta, tests, note=""):
    study_vals, _ = factor_values(meta, samples, "study")
    table_vals, _ = factor_values(meta, samples, "table")
    for name, _, klass in FACTORS:
        vals, _ = factor_values(meta, samples, name)
        res = stats.permanova(d2, vals, N_PERM, SEED)
        v = stats.cramers_v(vals, study_vals)
        tests.append({"block": tag, "factor": name, "class": klass, "mode": "marginal",
                      "n": res["n"], "n_levels": res["n_levels"],
                      "n_levels_ge2": res["n_levels_ge2"], "r2": res["r2"],
                      "r2_adj": res["r2_adj"], "r2_null_mean": res["r2_null_mean"],
                      "r2_excess": res["r2_excess"], "F": res["F"], "p": res["p"],
                      "testable": res["testable"], "reason": res["reason"],
                      "cramers_v_vs_study": v, "note": note})
        if name != "study":
            cond = stats.permanova_conditional(d2, study_vals, vals, N_PERM, SEED)
            tests.append({"block": tag, "factor": name, "class": klass,
                          "mode": "conditional_on_study", "n": cond["n"],
                          "n_levels": res["n_levels"], "n_levels_ge2": res["n_levels_ge2"],
                          "r2": cond["r2_partial"],
                          "r2_of_within_study": cond["r2_of_within_study"],
                          "df_factor": cond["df_factor"],
                          "n_in_varying_strata": cond["n_in_varying_strata"],
                          "F": cond["F"], "p": cond["p"], "testable": cond["testable"],
                          "reason": cond["reason"], "cramers_v_vs_study": v, "note": note})
        if name not in ("study", "table"):
            c2 = stats.permanova_conditional(d2, table_vals, vals, N_PERM, SEED)
            tests.append({"block": tag, "factor": name, "class": klass,
                          "mode": "conditional_on_study_table", "n": c2["n"],
                          "n_levels": res["n_levels"], "n_levels_ge2": res["n_levels_ge2"],
                          "r2": c2["r2_partial"],
                          "r2_of_within_study": c2["r2_of_within_study"],
                          "df_factor": c2["df_factor"],
                          "n_in_varying_strata": c2["n_in_varying_strata"],
                          "F": c2["F"], "p": c2["p"], "testable": c2["testable"],
                          "reason": c2["reason"], "cramers_v_vs_study": v, "note": note})
        if name != "study":
            vals2, idx = factor_values(meta, samples, name, complete_case=True)
            sub = [vals2[i] for i in idx]
            lv = Counter(sub)
            ge2 = sum(1 for x in lv.values() if x >= 2)
            if len(lv) >= 2 and ge2 >= 2:
                sd2 = [[d2[i][j] for j in idx] for i in idx]
                r = stats.permanova(sd2, sub, N_PERM, SEED)
                tests.append({"block": tag, "factor": name, "class": klass,
                              "mode": "complete_case", "n": r["n"],
                              "n_levels": r["n_levels"], "n_levels_ge2": r["n_levels_ge2"],
                              "r2": r["r2"], "r2_adj": r["r2_adj"],
                              "r2_null_mean": r["r2_null_mean"],
                              "r2_excess": r["r2_excess"], "F": r["F"], "p": r["p"],
                              "testable": r["testable"], "reason": r["reason"],
                              "cramers_v_vs_study": stats.cramers_v(
                                  sub, [study_vals[i] for i in idx]), "note": note})
            else:
                tests.append({"block": tag, "factor": name, "class": klass,
                              "mode": "complete_case", "n": len(idx),
                              "n_levels": len(lv), "n_levels_ge2": ge2,
                              "r2": None, "F": None, "p": None, "testable": False,
                              "reason": "after dropping (unstated), only %d levels have n>=2 "
                                        "(>=2 required)" % ge2,
                              "cramers_v_vs_study": None, "note": note})


def run_separability(matrix_tag, block_tag, d2, samples, meta, seps):

    study_vals, _ = factor_values(meta, samples, "study")
    for name, _, klass in FACTORS:
        if name in ("study", "table"):
            continue
        vals, _ = factor_values(meta, samples, name)
        des = separability.design_separability(study_vals, vals, UNSTATED)
        rec = {"matrix": matrix_tag, "block": block_tag, "factor": name, "class": klass,
               "n": des["n"], "n_levels": des["n_levels"],
               "n_levels_named": des["n_levels_named"], "n_studies": des["n_studies"],
               "n_studies_multilevel": des["n_studies_multilevel"],
               "n_samples_in_multilevel": des["n_samples_in_multilevel"],
               "n_levels_in_ge2_studies": des["n_levels_in_ge2_studies"],
               "n_samples_in_bridged_levels": des["n_samples_in_bridged_levels"],
               "cramers_v_vs_study": des["cramers_v_vs_study"],
               "separable": des["separable"], "sep_code": des["sep_code"],
               "separability_verdict": des["separability_verdict"],
               "_multilevel_studies": des["multilevel_studies"],
               "_bridged_levels": des["bridged_levels"]}
        if d2 is not None and des["n_studies_multilevel"] > 0:
            rc = separability.restricted_conditional(d2, study_vals, vals, N_PERM, SEED)
            rec.update(restricted_n=rc.get("n"), restricted_n_studies=rc.get("n_studies"),
                       restricted_df_factor=rc.get("df_factor"),
                       restricted_r2_partial=rc.get("r2_partial"),
                       restricted_r2_of_within_study=rc.get("r2_of_within_study"),
                       restricted_n_cells=rc.get("n_cells"),
                       restricted_n_cells_ge2=rc.get("n_cells_ge2"),
                       restricted_df_residual=rc.get("df_residual"),
                       restricted_F=rc.get("F"), restricted_p=rc.get("p"),
                       restricted_testable=rc.get("testable"),
                       restricted_reason=rc.get("reason"))
        else:
            rec.update(restricted_n=0, restricted_n_studies=0, restricted_df_factor=0,
                       restricted_r2_partial=None, restricted_r2_of_within_study=None,
                       restricted_n_cells=None, restricted_n_cells_ge2=None,
                       restricted_df_residual=None,
                       restricted_F=None, restricted_p=None, restricted_testable=False,
                       restricted_reason=("structurally inseparable: no study contains >=2 levels"
                                          if des["n_studies_multilevel"] == 0
                                          else "no usable distance matrix for this block"))
        seps.append(rec)


def fmt(v, nd=3):
    if v is None:
        return "-"
    if isinstance(v, float):
        if math.isnan(v):
            return "-"
        if math.isinf(v):
            return "inf"
        return ("%%.%df" % nd) % v
    return str(v)


def analyse_matrix(mk, value_field, matrix_tag, out_path, dist_recipes, tests, seps):

    rows = load(mk)
    cells, meta, names, units, multi, no_key = build_cells(rows, value_field)
    for s in meta:
        meta[s]["study_id"] = s[0]
    samples = sorted(meta)
    skeletons = sorted(set(sk for (_, sk) in cells))
    summ = sample_summary(cells, samples, skeletons)
    meas = measured_sets(cells)
    samples_with_data = sorted(meas)
    write_matrix(out_path, cells, samples, skeletons, meta, names, summ, units, mk)

    sk_studies = defaultdict(set)
    for (s, sk), (st, _) in cells.items():
        if st in ("quantified", "not_detected"):
            sk_studies[sk].add(s[0])
    n_stud_with_data = len(set(s[0] for s in samples_with_data))
    prev = Counter(len(v) for v in sk_studies.values())
    core_set = set(sk for sk, v in sk_studies.items()
                   if len(v) >= 0.5 * n_stud_with_data) if n_stud_with_data else set()
    st_c = Counter(v[0] for v in cells.values())
    n_cells = sum(1 for v in cells.values() if v[0] in ("quantified", "not_detected"))
    grid = len(samples) * len(skeletons) if samples and skeletons else 1

    info = {
        "mk": mk, "tag": matrix_tag, "path": out_path,
        "rows_in": len(rows), "dropped_no_key": no_key, "multi": multi,
        "samples": samples, "n_samples": len(samples),
        "n_compounds": len(skeletons), "skeletons": skeletons,
        "n_studies": len(set(s[0] for s in samples)),
        "n_studies_with_data": n_stud_with_data,
        "cells": n_cells, "fill_pct": 100.0 * n_cells / grid,
        "missing_pct": 100.0 - 100.0 * n_cells / grid,
        "status_counts": st_c, "prevalence": prev, "core": core_set,
        "meta": meta, "names": names, "units": units, "cells_map": cells,
        "meas": meas, "samples_with_data": samples_with_data, "summ": summ,
        "sk_studies": sk_studies,
        "single_cells": sum(1 for (sx, sk), v in cells.items()
                            if v[0] in ("quantified", "not_detected")
                            and len(sk_studies[sk]) == 1),
    }

    run_separability(matrix_tag, "ALL(design level)", None, samples_with_data, meta, seps)

    lad, order_c, per_c = ladder(meas)
    picks = pick_blocks(lad)
    info["ladder"] = lad
    info["per_c"] = per_c
    blocks = {}
    for tag in ["B1", "B2", "B3"]:
        if tag not in picks:
            continue
        b = picks[tag]
        keep, raw, closed, clrs, logs, capped = block_vectors(cells, b["samples"], b["compounds"])
        if len(keep) < 4:
            continue
        vecs = {"raw": raw, "closed": closed, "clr": clrs, "log": logs}
        blocks[tag] = dict(b, keep=keep, capped=capped, vecs=vecs)
        for dname, which, dfn in dist_recipes:
            d2 = stats.sq_dist_matrix(vecs[which], dfn)

            btag = "%s|%s-%s" % (matrix_tag, tag, dname)
            run_factor_tests(btag, d2, keep, meta, tests, note=matrix_tag)
            blocks[tag].setdefault("d2", {})[dname] = d2
        prim = dist_recipes[0][0]
        d2p = blocks[tag]["d2"][prim]
        run_separability(matrix_tag, "%s-%s" % (tag, prim), d2p, keep, meta, seps)
        blocks[tag]["nn_same_study"] = nn_same_study(d2p, keep, meta)
    info["blocks"] = blocks
    info["picks"] = picks

    within = []
    by_study = defaultdict(list)
    for s in samples_with_data:
        by_study[s[0]].append(s)
    for st, ss in sorted(by_study.items(), key=lambda x: -len(x[1])):
        if len(ss) < 4:
            continue
        common = set.intersection(*[meas[s] for s in ss])
        if len(common) < 3:
            within.append({"study": st, "n": len(ss), "n_common": len(common),
                           "factor": "-", "class": "-", "r2": None, "p": None,
                           "reason": "< 3 compounds measured in common within the study, "
                                     "no compositional analysis"})
            continue
        comps = sorted(common, key=lambda c: (-len(per_c[c]), c))
        keep, raw, closed, clrs, logs, _ = block_vectors(cells, ss, comps)
        if len(keep) < 4:
            continue
        which = dist_recipes[0][1]
        d2 = stats.sq_dist_matrix({"raw": raw, "closed": closed, "clr": clrs,
                                   "log": logs}[which], dist_recipes[0][2])
        for name, _, klass in FACTORS:
            if name == "study":
                continue
            vals, _ = factor_values(meta, keep, name)
            lv = Counter(vals)
            if len(lv) < 2:
                continue
            ge2 = sum(1 for v in lv.values() if v >= 2)
            if ge2 < 2:
                within.append({"study": st, "n": len(keep), "n_common": len(comps),
                               "factor": name, "class": klass, "r2": None, "p": None,
                               "reason": "among %d levels only %d have n>=2 (>=2 required), "
                                         "insufficient df, untestable" % (len(lv), ge2)})
                continue
            r = stats.permanova(d2, vals, N_PERM, SEED)

            ok = r["testable"]
            within.append({"study": st, "n": len(keep), "n_common": len(comps),
                           "factor": name, "class": klass,
                           "r2": r["r2"] if ok else None,
                           "p": r["p"] if ok else None,
                           "reason": ("testable (%d levels, %d with n>=2)"
                                      % (r["n_levels"], ge2)) if ok
                                     else ("untestable: " + r["reason"])})
    info["within"] = within
    return info


def nn_same_study(d2, samps, meta):
    hit, n = 0, len(samps)
    for i in range(n):
        best, bj = None, None
        for j in range(n):
            if i == j:
                continue
            if best is None or d2[i][j] < best:
                best, bj = d2[i][j], j
        if bj is not None and meta[samps[i]]["study_id"] == meta[samps[bj]]["study_id"]:
            hit += 1
    return hit / n if n else float("nan")


def main():
    tests, seps = [], []

    REL_DIST = [("CLR/Aitchison", "clr", stats.euclidean),
                ("raw/Bray-Curtis", "closed", stats.bray_curtis)]

    ABS_DIST = [("log10/Euclidean", "log", stats.euclidean),
                ("raw/Bray-Curtis", "raw", stats.bray_curtis),
                ("CLR/Aitchison", "clr", stats.euclidean)]

    rel = analyse_matrix("relative_pct", "area_pct_uncorrected", "REL",
                         OUT_MATRIX_REL, REL_DIST, tests, seps)
    ab = analyse_matrix("absolute_conc", "value", "ABS",
                        OUT_MATRIX_ABS, ABS_DIST, tests, seps)

    order_s = rel["samples_with_data"]
    masks = [rel["meas"][s] for s in order_s]
    d2_mask = stats.sq_dist_matrix(masks, stats.jaccard)
    MASK_TAG = "MASK(%d samples, Jaccard-missingness pattern)" % len(order_s)
    run_factor_tests(MASK_TAG, d2_mask, order_s, rel["meta"], tests,
                     note="missingness pattern, not chemistry")
    nn_mask = nn_same_study(d2_mask, order_s, rel["meta"])
    cnt = Counter(rel["meta"][s]["study_id"] for s in order_s)
    N = len(order_s)
    base_same = sum(v * (v - 1) for v in cnt.values()) / (N * (N - 1)) if N > 1 else float("nan")

    zf = []
    for s in order_s:
        v = [(rel["cells_map"].get((s, sk))[1]
              if (rel["cells_map"].get((s, sk)) and rel["cells_map"][(s, sk)][1] is not None)
              else 0.0) for sk in rel["skeletons"]]
        zf.append(stats.closure(v, 100.0) if sum(v) > 0 else [0.0] * len(rel["skeletons"]))
    d2_zero = stats.sq_dist_matrix(zf, stats.bray_curtis)
    zero_study = stats.permanova(d2_zero, [rel["meta"][s]["study_id"] for s in order_s],
                                 N_PERM, SEED)

    sens = []
    if "B2" in rel["blocks"]:
        b = rel["picks"]["B2"]
        for f in DELTA_SENSITIVITY:
            keep, raw, closed, clrs, logs, _ = block_vectors(
                rel["cells_map"], b["samples"], b["compounds"], f)
            d2 = stats.sq_dist_matrix(clrs, stats.euclidean)
            r = stats.permanova(d2, [rel["meta"][s]["study_id"] for s in keep], N_PERM, SEED)
            sens.append((f, r["r2"], r["p"]))

    tcols = ["block", "factor", "class", "mode", "n", "n_levels", "n_levels_ge2",
             "r2", "r2_adj", "r2_null_mean", "r2_excess", "r2_of_within_study",
             "df_factor", "n_in_varying_strata", "F", "p", "testable",
             "cramers_v_vs_study", "reason", "note"]
    with open(OUT_TESTS, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(tcols)
        for t in tests:
            w.writerow([fmt(t.get(c), 4) if isinstance(t.get(c), float)
                        else (t.get(c) if t.get(c) is not None else "") for c in tcols])
    with open(OUT_SEP, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(separability.SEPARABILITY_FIELDS + ["multilevel_studies", "bridged_levels"])
        for s in seps:
            w.writerow([fmt(s.get(c), 4) if isinstance(s.get(c), float)
                        else (s.get(c) if s.get(c) is not None else "")
                        for c in separability.SEPARABILITY_FIELDS] +
                       [";".join(s["_multilevel_studies"]), ";".join(s["_bridged_levels"])])

    print("written:")
    print(" ", os.path.basename(OUT_MATRIX_REL),
          "(%d samples x %d compounds, %.2f%% missing)" %
          (rel["n_samples"], rel["n_compounds"], rel["missing_pct"]))
    print(" ", os.path.basename(OUT_MATRIX_ABS),
          "(%d samples x %d compounds, %.2f%% missing)" %
          (ab["n_samples"], ab["n_compounds"], ab["missing_pct"]))
    print(" ", os.path.basename(OUT_TESTS), "(%d tests)" % len(tests))
    print(" ", os.path.basename(OUT_SEP),
          "(%d separability diagnostics)" % len(seps))


if __name__ == "__main__":
    main()
