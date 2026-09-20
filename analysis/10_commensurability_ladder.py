import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import collections
import csv
import glob
import itertools
import json

from comparability import config, analytical, reporting, samples, tables

import numpy as np

SEED = 20260914
RNG = np.random.default_rng(SEED)
N_SIM = 2000


def read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore", restval="")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("-> %s  (%d rows)" % (os.path.basename(path), len(rows)))


studies = read_tsv(config.data_path("studies.tsv"))
sample_rows = read_tsv(config.data_path("samples.tsv"))
samples_by_study = collections.defaultdict(list)
for _s in sample_rows:
    samples_by_study[_s["study_id"]].append(_s)

_pdf_cache = {}
for _p in (config.data_path("article_text_cache.json"),
           config.results_path("article_text_cache.json")):
    if os.path.exists(_p):
        _pdf_cache.update(json.load(open(_p, encoding="utf-8")))
_pdf_paths = {os.path.basename(p): p for p in
              glob.glob(os.path.join(config.DATA_DIR, "oa_pdf", "*.pdf")) +
              glob.glob(os.path.join(config.DATA_DIR, "oa_pdf_recovered", "*.pdf"))}

texts, roots = {}, {}
for st in studies:
    sid = st["study_id"]
    if st["pmcid"]:
        root = tables.parse_article(
            os.path.join(config.FULLTEXT_DIR, st["pmcid"] + ".xml"))
        mt, _src = analytical.analytical_text(root)
        texts[sid] = {"methods": mt, "whole": analytical.body_text(root)}
        roots[sid] = root
    else:
        name = st["doi"].replace("/", "_") + ".pdf"
        raw = _pdf_cache.get(name, "")
        if not raw and name in _pdf_paths:
            raise SystemExit("PDF text cache miss for %s; run 09_reporting_census.py first" % name)
        texts[sid] = {"methods": raw, "whole": raw}
        roots[sid] = reporting.wrap_text(raw)

analytic = {sid: analytical.parse_analytical(roots[sid]) for sid in texts}
ALL = [s["study_id"] for s in studies]
print("corpus: %d studies, %d samples" % (len(studies), len(sample_rows)))


_EXTRACTION_NAMES = {"hydrodistillation", "steam_distillation",
                     "supercritical_co2", "solvent_extraction", "microwave",
                     "ultrasound", "headspace"}
extraction_raw = {}
for st in studies:
    sid = st["study_id"]
    text = texts[sid]["methods"] or texts[sid]["whole"]
    procs = [name for name, pat in samples.PROCESSING
             if name in _EXTRACTION_NAMES and pat.search(text)]
    intro = analytic[sid]["study_platform_introduction"]
    extraction_raw[sid] = ";".join(procs + ([intro] if intro else []))


_EXTRACTION_ORDER = ["hydrodistillation", "steam_distillation",
                     "supercritical_co2", "microwave", "ultrasound",
                     "solvent_extraction", "HS-SPME", "SPME", "SAFE", "SDE",
                     "headspace", "TD"]


def extraction_level(sid):
    named = set(extraction_raw[sid].split(";"))
    for name in _EXTRACTION_ORDER:
        if name in named:
            return name
    return ""


AXES = ["species", "cultivar", "organ", "stage", "origin_country", "region",
        "material_type"]
FIELD_OF = {"origin_country": "country"}
inscope = [s for s in sample_rows if s["in_composition_scope"] == "yes"]


def level_of(sample, axis):
    field = FIELD_OF.get(axis, axis)
    val = (sample.get(field) or "").strip()
    if axis == "cultivar":
        return reporting.normalize_cultivar(val)
    return val


levels_by_axis = {}
for axis in AXES:
    by_level = collections.defaultdict(set)
    for s in inscope:
        val = level_of(s, axis)
        if val:
            by_level[val].add(s["study_id"])
    levels_by_axis[axis] = by_level

all_pairs = []
for axis in AXES:
    for level, sids in levels_by_axis[axis].items():
        for a, b in itertools.combinations(sorted(sids), 2):
            all_pairs.append((axis, level, a, b))
N_PAIRS = len(all_pairs)
print("cross-study pairs sharing a biological level: %d" % N_PAIRS)


cell_samples = collections.defaultdict(list)
for s in inscope:
    for axis in AXES:
        val = level_of(s, axis)
        if val:
            cell_samples[(axis, val, s["study_id"])].append(s)


RUNGS = ["measure_kind", "material_type", "extraction_method",
         "column_polarity", "quant_mode", "ri_basis_named"]


RUNG_CATEGORY = {
    "measure_kind": "analytical",
    "material_type": "design",
    "extraction_method": "design",
    "column_polarity": "analytical",
    "quant_mode": "analytical",
    "ri_basis_named": "analytical",
}
RUNG_BLOCKING = {"measure_kind": True, "material_type": True,
                 "extraction_method": True, "column_polarity": False,
                 "quant_mode": False, "ri_basis_named": False}


SENTINELS = {"measure_kind": {"unspecified"}, "material_type": {"unknown"}}

DEFAULT_TREATMENT = {"extraction": "raw", "compare": "exact",
                     "sentinels": "value", "quant": "first", "scope": "study"}


def _tok(vals):
    return ";".join(sorted(vals))


def attributes(rung, axis, level, sid, tr):

    if rung in ("measure_kind", "material_type"):
        if tr["scope"] == "level":
            rows = cell_samples[(axis, level, sid)]
        else:
            rows = [s for s in samples_by_study[sid]
                    if s["in_composition_scope"] == "yes"]
        vals = {s[rung] for s in rows if s[rung]}
        if tr["sentinels"] == "missing":
            vals -= SENTINELS[rung]
        return _tok(vals), frozenset(vals)
    if rung == "extraction_method":
        if tr["extraction"] == "canon":
            lev = extraction_level(sid)
            return lev, frozenset([lev]) if lev else frozenset()
        raw = extraction_raw[sid]
        return raw, frozenset(v for v in raw.split(";") if v)
    if rung == "column_polarity":
        v = analytic[sid]["study_column_polarity"] or ""
        return v, frozenset([v]) if v else frozenset()
    if rung == "quant_mode":
        full = analytic[sid]["study_quant_mode"] or ""
        parts = frozenset(v for v in full.split(";") if v)
        if tr["quant"] == "set":
            return _tok(parts), parts
        first = full.split(";")[0]
        return first, frozenset([first]) if first else frozenset()
    if rung == "ri_basis_named":
        v = analytic[sid]["study_ri_basis_named"] or ""
        return v, frozenset([v]) if v else frozenset()
    raise KeyError(rung)


def build_attrs(tr):

    out = {}
    if tr["scope"] == "level":
        keys = {(ax, lv, s) for ax, lv, a, b in all_pairs for s in (a, b)}
        for rung in RUNGS:
            for ax, lv, sid in keys:
                out[(rung, (ax, lv, sid))] = attributes(rung, ax, lv, sid, tr)
    else:
        for rung in RUNGS:
            for sid in ALL:
                out[(rung, (None, None, sid))] = attributes(rung, None, None, sid, tr)
    return out


def cell_key(tr, axis, level, sid):
    return (axis, level, sid) if tr["scope"] == "level" else (None, None, sid)


def rung_verdict(attrs, tr, rung, ka, kb):

    ta, sa = attrs[(rung, ka)]
    tb, sb = attrs[(rung, kb)]
    if not sa or not sb:
        return "unknown"
    if tr["compare"] == "overlap":
        return "match" if (sa & sb) else "differ"
    return "match" if ta == tb else "differ"


RUNG_SETS = {

    "upto1_measure": RUNGS[:1],
    "upto2_material": RUNGS[:2],
    "blocking3": RUNGS[:3],
    "upto4_polarity": RUNGS[:4],
    "no_ri": RUNGS[:5],
    "full6": RUNGS,

    "no_polarity": [r for r in RUNGS if r != "column_polarity"],
    "no_polarity_quant": [r for r in RUNGS
                          if r not in ("column_polarity", "quant_mode")],
    "design2": [r for r in RUNGS if RUNG_CATEGORY[r] == "design"],
}


def classify(attrs, tr, rungs):

    recs = []
    for axis, level, a, b in all_pairs:
        ka, kb = cell_key(tr, axis, level, a), cell_key(tr, axis, level, b)
        differ, unknown = [], []
        for rung in rungs:
            v = rung_verdict(attrs, tr, rung, ka, kb)
            if v == "differ":
                differ.append(rung)
            elif v == "unknown":
                unknown.append(rung)
        if differ:
            status = "incompatible"
            cats = {RUNG_CATEGORY[r] for r in differ}
            if "design" in cats:
                verdict = "design_difference"
            elif any(RUNG_BLOCKING[r] for r in differ):
                verdict = "analytical_blocking"
            else:
                verdict = "analytical_nonblocking"
        elif unknown:
            status, verdict = "undecidable", "reporting_gap"
        else:
            status, verdict = "poolable", "poolable"
        recs.append({"axis": axis, "level": level, "study_a": a, "study_b": b,
                     "status": status, "verdict": verdict,
                     "differ_rungs": ";".join(differ),
                     "unknown_rungs": ";".join(unknown)})
    return recs


def tally(recs):
    c = collections.Counter(r["status"] for r in recs)
    v = collections.Counter(r["verdict"] for r in recs)
    n = len(recs)
    return {"pairs": n, "poolable": c["poolable"],
            "incompatible": c["incompatible"], "undecidable": c["undecidable"],
            "pct_reporting": round(100.0 * c["undecidable"] / n, 1),
            "pct_incompatible": round(100.0 * c["incompatible"] / n, 1),
            "design_difference": v["design_difference"],
            "analytical_blocking": v["analytical_blocking"],
            "analytical_nonblocking": v["analytical_nonblocking"],
            "pct_design": round(100.0 * v["design_difference"] / n, 1),
            "pct_analytical": round(
                100.0 * (v["analytical_blocking"]
                         + v["analytical_nonblocking"]) / n, 1)}


SWITCHES = [
    ("extraction", ["raw", "canon"],
     "the concatenated list of extraction routes mentioned in the methods "
     "section vs a single route picked by precedence"),
    ("compare", ["exact", "overlap"],
     "exact string equality on the set-union vs a non-empty set "
     "intersection"),
    ("sentinels", ["value", "missing"],
     "measure_kind='unspecified' / material_type='unknown' compared as "
     "values vs read as not stated"),
    ("quant", ["first", "set"],
     "quant_mode reduced to its first vocabulary element vs the whole "
     "declared set"),
    ("scope", ["study", "level"],
     "attributes taken over the whole study vs over the samples at the "
     "shared biological level"),
]

CUMULATIVE = [
    ("C0_string_exact", dict(DEFAULT_TREATMENT),
     "exact equality, raw extraction string, study-wide attributes, "
     "sentinels compared as values"),
    ("C1_extraction_canonical", dict(DEFAULT_TREATMENT, extraction="canon"),
     "+ one extraction route per study, picked by precedence"),
    ("C2_set_overlap", dict(DEFAULT_TREATMENT, extraction="canon",
                            compare="overlap"),
     "+ set overlap instead of exact equality for multi-valued attributes"),
    ("C3_sentinels_missing", dict(DEFAULT_TREATMENT, extraction="canon",
                                  compare="overlap", sentinels="missing"),
     "+ 'unspecified' and 'unknown' read as not stated, as the data "
     "dictionary defines them"),
    ("C4_quant_set", dict(DEFAULT_TREATMENT, extraction="canon",
                          compare="overlap", sentinels="missing", quant="set"),
     "+ quant_mode compared as a set, not as its first vocabulary element"),
    ("C5_level_scope", dict(DEFAULT_TREATMENT, extraction="canon",
                            compare="overlap", sentinels="missing",
                            quant="set", scope="level"),
     "+ attributes computed over the samples at the shared biological level, "
     "not over the whole study; this is the default coding"),
]
RECOMMENDED = "C5_level_scope"


def tr_name(tr):
    return "|".join("%s=%s" % (k, tr[k]) for k, _v, _d in SWITCHES)


all_treatments = {}
for name, tr, _doc in CUMULATIVE:
    all_treatments[name] = tr
for combo in itertools.product(*[vals for _k, vals, _d in SWITCHES]):
    tr = {k: v for (k, _vals, _d), v in zip(SWITCHES, combo)}
    all_treatments.setdefault("F:" + tr_name(tr), tr)

attrs_cache = {}
for name, tr in all_treatments.items():
    key = tr_name(tr)
    if key not in attrs_cache:
        attrs_cache[key] = build_attrs(tr)

sens_rows, pair_detail = [], []
print("\n" + "=" * 78)
print("SENSITIVITY - same %d pairs, same rungs, only the coding changed"
      % N_PAIRS)
print("=" * 78)
print("  %-22s %-18s %5s %5s %5s %5s   %s"
      % ("treatment", "rung set", "pool", "inc", "und", "rep%", "design%/analyt%"))
for name, tr in sorted(all_treatments.items(),
                       key=lambda kv: (kv[0].startswith("F:"), kv[0])):
    attrs = attrs_cache[tr_name(tr)]
    for rs_name, rungs in RUNG_SETS.items():
        recs = classify(attrs, tr, rungs)
        t = tally(recs)
        row = {"treatment": name, "switches": tr_name(tr), "rung_set": rs_name,
               "rungs": " + ".join(rungs), "n_rungs": len(rungs)}
        row.update(t)
        sens_rows.append(row)
        if not name.startswith("F:") and rs_name == "full6":
            print("  %-22s %-18s %5d %5d %5d %5.1f   %.1f / %.1f"
                  % (name, rs_name, t["poolable"], t["incompatible"],
                     t["undecidable"], t["pct_reporting"],
                     t["pct_design"], t["pct_analytical"]))
        if name == RECOMMENDED and rs_name == "full6":
            for r in recs:
                r["treatment"] = name
                pair_detail.append(r)

write_tsv(config.results_path("ladder_sensitivity.tsv"), sens_rows,
          ["treatment", "switches", "rung_set", "rungs", "n_rungs", "pairs",
           "poolable", "incompatible", "undecidable", "pct_reporting",
           "pct_incompatible", "design_difference", "analytical_blocking",
           "analytical_nonblocking", "pct_design", "pct_analytical"])
write_tsv(config.results_path("ladder_pair_detail.tsv"), pair_detail,
          ["treatment", "axis", "level", "study_a", "study_b", "status",
           "verdict", "differ_rungs", "unknown_rungs"])


reason_rows = []
for name, tr, _doc in CUMULATIVE:
    attrs = attrs_cache[tr_name(tr)]
    full = classify(attrs, tr, RUNGS)
    for i, rung in enumerate(RUNGS):
        ind = collections.Counter()
        for axis, level, a, b in all_pairs:
            ka, kb = cell_key(tr, axis, level, a), cell_key(tr, axis, level, b)
            ind[rung_verdict(attrs, tr, rung, ka, kb)] += 1
        prev = classify(attrs, tr, RUNGS[:i]) if i else None
        here = classify(attrs, tr, RUNGS[:i + 1])
        newly_no = sum(1 for j, r in enumerate(here)
                       if r["status"] == "incompatible"
                       and (prev is None or prev[j]["status"] != "incompatible"))
        newly_un = sum(1 for j, r in enumerate(here)
                       if r["status"] == "undecidable"
                       and (prev is None or prev[j]["status"] != "undecidable"))
        reason_rows.append({
            "treatment": name, "rung": rung,
            "category": RUNG_CATEGORY[rung],
            "blocking": "yes" if RUNG_BLOCKING[rung] else "no",
            "independent_match": ind["match"], "independent_differ": ind["differ"],
            "independent_unknown": ind["unknown"],
            "marginal_newly_incompatible": newly_no,
            "marginal_newly_undecidable": newly_un,
            "pairs_where_this_is_the_only_differing_rung": sum(
                1 for r in full if r["differ_rungs"] == rung),
            "pairs_where_this_is_the_only_unknown_rung": sum(
                1 for r in full if r["unknown_rungs"] == rung
                and not r["differ_rungs"])})
write_tsv(config.results_path("ladder_exclusion_reasons.tsv"), reason_rows,
          ["treatment", "rung", "category", "blocking", "independent_match",
           "independent_differ", "independent_unknown",
           "marginal_newly_incompatible", "marginal_newly_undecidable",
           "pairs_where_this_is_the_only_differing_rung",
           "pairs_where_this_is_the_only_unknown_rung"])

print("\nEXCLUSION REASONS, independent over all %d pairs" % N_PAIRS)
print("-" * 78)
print("  %-18s %-11s %-5s  %25s  %25s"
      % ("rung", "category", "block", "C0 m/d/u",
         "%s m/d/u" % RECOMMENDED))
for rung in RUNGS:
    a = next(r for r in reason_rows
             if r["treatment"] == "C0_string_exact" and r["rung"] == rung)
    b = next(r for r in reason_rows
             if r["treatment"] == RECOMMENDED and r["rung"] == rung)
    print("  %-18s %-11s %-5s  %25s  %25s"
          % (rung, a["category"], a["blocking"],
             "%d / %d / %d" % (a["independent_match"], a["independent_differ"],
                               a["independent_unknown"]),
             "%d / %d / %d" % (b["independent_match"], b["independent_differ"],
                               b["independent_unknown"])))


def counterfactual(tr, rungs, n_sim=N_SIM, rng=None):
    rng = rng or np.random.default_rng(SEED)
    attrs = attrs_cache[tr_name(tr)]
    keys = sorted({k for (_r, k) in attrs})
    pools = {}
    for rung in rungs:
        obs = [attrs[(rung, k)] for k in keys if attrs[(rung, k)][1]]
        pools[rung] = obs
    empty = {rung: [k for k in keys if not attrs[(rung, k)][1]] for rung in rungs}
    pair_keys = [(cell_key(tr, ax, lv, a), cell_key(tr, ax, lv, b))
                 for ax, lv, a, b in all_pairs]
    yes = []
    for _ in range(n_sim):
        filled = {}
        for rung in rungs:
            pool = pools[rung]
            if not pool:
                continue
            for k in empty[rung]:
                filled[(rung, k)] = pool[int(rng.integers(len(pool)))]
        n_yes = 0
        for ka, kb in pair_keys:
            good = True
            for rung in rungs:
                ta, sa = filled.get((rung, ka), attrs[(rung, ka)])
                tb, sb = filled.get((rung, kb), attrs[(rung, kb)])
                if not sa or not sb:
                    good = False
                    break
                hit = bool(sa & sb) if tr["compare"] == "overlap" else ta == tb
                if not hit:
                    good = False
                    break
            if good:
                n_yes += 1
        yes.append(n_yes)
    return (float(np.mean(yes)), float(np.percentile(yes, 2.5)),
            float(np.percentile(yes, 97.5)))


print("\nOPTIMISTIC COUNTERFACTUAL (every rung reported by every study)")
print("-" * 78)
cf_rows = []
for name, tr, _doc in CUMULATIVE:
    for rs_name in ("full6", "no_polarity_quant", "blocking3"):
        rungs = RUNG_SETS[rs_name]
        mean, lo, hi = counterfactual(tr, rungs,
                                      rng=np.random.default_rng(SEED))
        now = tally(classify(attrs_cache[tr_name(tr)], tr, rungs))
        cf_rows.append({"treatment": name, "switches": tr_name(tr),
                        "rung_set": rs_name, "n_rungs": len(rungs),
                        "pairs": N_PAIRS, "poolable_now": now["poolable"],
                        "undecidable_now": now["undecidable"],
                        "poolable_counterfactual_mean": round(mean, 1),
                        "ci_lo": round(lo), "ci_hi": round(hi),
                        "gain": round(mean - now["poolable"], 1),
                        "seed": SEED, "n_sim": N_SIM})
        if rs_name == "full6":
            print("  %-22s %-18s now=%3d  counterfactual=%6.1f  [%.0f-%.0f]"
                  % (name, rs_name, now["poolable"], mean, lo, hi))
write_tsv(config.results_path("ladder_counterfactual.tsv"), cf_rows,
          ["treatment", "switches", "rung_set", "n_rungs", "pairs",
           "poolable_now", "undecidable_now", "poolable_counterfactual_mean",
           "ci_lo", "ci_hi", "gain", "seed", "n_sim"])


STATE = {"sens_rows": sens_rows, "reason_rows": reason_rows,
         "cf_rows": cf_rows,
         "pair_detail": pair_detail}


T0 = dict(CUMULATIVE[0][1])
T0_attrs = attrs_cache[tr_name(T0)]
TR = dict(CUMULATIVE[-1][1])
TR_attrs = attrs_cache[tr_name(TR)]
t0_recs = classify(T0_attrs, T0, RUNGS)
tr_recs = classify(TR_attrs, TR, RUNGS)


axis_rows = []
for axis in AXES:
    for name, recs in (("C0_string_exact", t0_recs), (RECOMMENDED, tr_recs)):
        sub = [r for r in recs if r["axis"] == axis]
        if not sub:
            continue
        t = tally(sub)
        axis_rows.append(dict(t, axis=axis, treatment=name))
write_tsv(config.results_path("ladder_by_axis.tsv"), axis_rows,
          ["treatment", "axis", "pairs", "poolable", "incompatible",
           "undecidable", "pct_reporting", "design_difference",
           "analytical_blocking", "analytical_nonblocking",
           "pct_design", "pct_analytical"])


rec = next(r for r in sens_rows
           if r["treatment"] == RECOMMENDED and r["rung_set"] == "full6")
print("\n%s -> %d poolable / %d incompatible / %d undecidable; "
      "reporting %.1f%%, design %.1f%%, analytical %.1f%%"
      % (RECOMMENDED, rec["poolable"], rec["incompatible"], rec["undecidable"],
         rec["pct_reporting"], 100.0 * rec["design_difference"] / N_PAIRS,
         100.0 * (rec["analytical_blocking"] + rec["analytical_nonblocking"]) / N_PAIRS))
