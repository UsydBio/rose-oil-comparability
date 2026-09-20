import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import collections
import csv
import glob
import itertools
import json
import math
import statistics

from comparability import config, analytical, methods, reporting, samples, tables

import numpy as np
from scipy import stats

RNG = np.random.default_rng(20260913)
N_PERM = 20000


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


def pct(k, n):
    return (100.0 * k / n) if n else float("nan")


studies = read_tsv(config.data_path("studies.tsv"))
sample_rows = read_tsv(config.data_path("samples.tsv"))
obs_rel = read_tsv(config.data_path("observations_relative_pct.tsv"))
obs_abs = read_tsv(config.data_path("observations_absolute_conc.tsv"))
obs_oth = read_tsv(config.data_path("observations_other.tsv"))
prior_analytical = {r["pmcid"]: r for r in
                    read_tsv(config.data_path("study_analytical.tsv"))}

samples_by_study = collections.defaultdict(list)
for s in sample_rows:
    samples_by_study[s["study_id"]].append(s)

print("corpus: %d studies, %d samples" % (len(studies), len(sample_rows)))


pdf_paths = {os.path.basename(p): p for p in
             glob.glob(os.path.join(config.DATA_DIR, "oa_pdf", "*.pdf")) +
             glob.glob(os.path.join(config.DATA_DIR, "oa_pdf_recovered", "*.pdf"))}


PDF_CACHE = config.results_path("article_text_cache.json")
_pdf_cache = {}
for _p in (config.data_path("article_text_cache.json"), PDF_CACHE):
    if os.path.exists(_p):
        _pdf_cache.update(json.load(open(_p, encoding="utf-8")))


def pdf_text(path):
    key = os.path.basename(path)
    if key in _pdf_cache:
        return _pdf_cache[key]
    import pdfplumber
    with pdfplumber.open(path) as pdf:
        text = "\n".join((page.extract_text() or "") for page in pdf.pages)
    _pdf_cache[key] = text
    with open(PDF_CACHE, "w", encoding="utf-8") as fh:
        json.dump(_pdf_cache, fh)
    return text


texts, roots, channels, reviews = {}, {}, {}, {}
for st in studies:
    sid = st["study_id"]
    if st["pmcid"]:
        root = tables.parse_article(
            os.path.join(config.FULLTEXT_DIR, st["pmcid"] + ".xml"))
        mt, src = analytical.analytical_text(root)
        texts[sid] = {"methods": mt, "whole": analytical.body_text(root)}
        channels[sid] = "jats:" + src
        roots[sid] = root
    else:
        name = st["doi"].replace("/", "_") + ".pdf"
        raw = pdf_text(pdf_paths[name]) if name in pdf_paths else ""
        texts[sid] = {"methods": raw, "whole": raw}
        channels[sid] = "pdf_fulltext" if raw else "pdf_missing"
        roots[sid] = reporting.wrap_text(raw)
    is_rev, rev_ev = reporting.is_review(roots[sid])
    reviews[sid] = (is_rev, rev_ev)

print("review articles detected: %s"
      % ", ".join("%s (%s)" % (k, v[1][:40]) for k, v in reviews.items() if v[0]))

analytic = {sid: analytical.parse_analytical(roots[sid]) for sid in texts}


ALL = [s["study_id"] for s in studies]
VOLATILE = [s for s in ALL if s not in reporting.NON_VOLATILE_STUDIES]
PRIMARY = [s for s in VOLATILE if not reviews[s][0]]
print("populations: all=%d  volatile/GC=%d  volatile primary (no reviews)=%d"
      % (len(ALL), len(VOLATILE), len(PRIMARY)))


ri_in_table, sd_in_table, range_in_table, quant_n = (
    collections.Counter(), collections.Counter(),
    collections.Counter(), collections.Counter())
for rows, sd_field, lo_field in (
        (obs_rel, "area_pct_sd", "area_pct_range_low"),
        (obs_abs, "concentration_sd", "concentration_range_low"),
        (obs_oth, "value_sd", "value_range_low")):
    for r in rows:
        sid = r["study_id"]
        if r.get("reported_retention_index"):
            ri_in_table[sid] += 1
        if r["detection_status"] != "quantified":
            continue
        quant_n[sid] += 1
        if r.get(sd_field):
            sd_in_table[sid] += 1
        if r.get(lo_field):
            range_in_table[sid] += 1

_EXTRACTION_NAMES = {"hydrodistillation", "steam_distillation",
                     "supercritical_co2", "solvent_extraction", "microwave",
                     "ultrasound", "headspace"}


def provenance(sid, field):

    rows = [s for s in samples_by_study.get(sid, [])
            if s["in_composition_scope"] == "yes"]
    if not rows:
        return None, "no in-scope samples"
    have = [s for s in rows if (s.get(field) or "").strip()]
    ev = "%d/%d in-scope samples" % (len(have), len(rows))
    if have:
        ev += ": " + ";".join(sorted({s[field] for s in have}))[:120]
    return len(have) == len(rows), ev


study_items, study_detail = {}, {}
for st in studies:
    sid = st["study_id"]
    text = texts[sid]["methods"] or texts[sid]["whole"]
    whole = texts[sid]["whole"]
    a = analytic[sid]
    got, ev = {}, {}

    det = reporting.detectors(text)
    got["platform_detector"] = bool(det)
    ev["platform_detector"] = ";".join(det)

    for key, field in (("column_phase", "study_column_phase"),
                       ("column_dimensions", "study_column_dimensions"),
                       ("carrier_gas", "study_carrier_gas"),
                       ("oven_program", "study_oven_program"),
                       ("injector_temp", "study_injector_temp_c"),
                       ("split_mode", "study_split_mode"),
                       ("ri_basis_named", "study_ri_basis_named"),
                       ("ri_alkane_series", "study_ri_alkane_series"),
                       ("ident_basis", "study_ident_basis"),
                       ("quant_mode", "study_quant_mode"),
                       ("internal_standard", "study_internal_standard")):
        got[key] = bool(a[field])
        ev[key] = str(a[field])[:150]

    got["ri_reported"] = bool(a["study_ri_basis"]) or ri_in_table[sid] > 0
    ev["ri_reported"] = "basis=%s; %d RI values in extracted tables" % (
        a["study_ri_basis"] or "-", ri_in_table[sid])

    rf, rf_ev = reporting.response_factor(text)
    got["response_factor"] = bool(rf)
    ev["response_factor"] = ("%s: %s" % (rf, rf_ev)) if rf else ""

    for key, field in (("species", "species"), ("cultivar", "cultivar")):
        val, detail = provenance(sid, field)
        got[key] = val
        ev[key] = detail

    organs, organ_ev = reporting.organ_stated(text)
    got["organ"] = bool(organs)
    samp_org, samp_org_ev = provenance(sid, "organ")
    ev["organ"] = "%s | prose: %s | sample field: %s" % (
        ";".join(organs), organ_ev[:110], samp_org_ev)
    stages, stage_ev = reporting.stage_stated(text)
    got["stage"] = bool(stages)
    samp_stage, samp_stage_ev = provenance(sid, "stage")
    ev["stage"] = "%s | prose: %s | sample field: %s" % (
        ";".join(stages), stage_ev[:110], samp_stage_ev)

    month, year, day = reporting.harvest_date(text)
    got["harvest_month"], ev["harvest_month"] = bool(month), month
    got["harvest_year"], ev["harvest_year"] = bool(year), year
    got["harvest_day"], ev["harvest_day"] = bool(day), day
    tod = reporting.harvest_time(text)
    got["harvest_time"], ev["harvest_time"] = bool(tod), tod

    country, region, gps = reporting.geographic_origin(
        text, methods._COUNTRY_RE, methods._REGION_RE)
    got["origin_country"], ev["origin_country"] = bool(country), country
    got["origin_region"], ev["origin_region"] = bool(region), region
    got["origin_gps"], ev["origin_gps"] = bool(gps), gps

    mats = [s.get("material_type") for s in samples_by_study.get(sid, [])
            if s["in_composition_scope"] == "yes"]
    if not mats:
        got["material_type"], ev["material_type"] = None, "no in-scope samples"
    else:
        named = [x for x in mats if x and x != "unknown"]
        got["material_type"] = len(named) == len(mats)
        ev["material_type"] = "%d/%d in-scope samples: %s" % (
            len(named), len(mats), ";".join(sorted(set(mats)))[:100])

    procs = [name for name, pat in samples.PROCESSING
             if name in _EXTRACTION_NAMES and pat.search(text)]
    intro = a["study_platform_introduction"]
    got["extraction_method"] = bool(procs or intro)
    ev["extraction_method"] = ";".join(procs + ([intro] if intro else []))

    t_ev, temp_ev, ratio_ev = reporting.extraction_params(text)
    got["extraction_time"], ev["extraction_time"] = bool(t_ev), t_ev
    got["extraction_temp"], ev["extraction_temp"] = bool(temp_ev), temp_ev
    got["extraction_ratio"], ev["extraction_ratio"] = bool(ratio_ev), ratio_ev

    loose, strict = reporting.replicates(text)
    got["replicates_strict"] = bool(strict)
    ev["replicates_strict"] = strict
    got["_replicates_loose"] = bool(loose)
    ev["_replicates_loose"] = loose

    n_sd = sd_in_table[sid] + range_in_table[sid]
    got["dispersion_reported"] = n_sd > 0 if quant_n[sid] else None
    ev["dispersion_reported"] = "%d/%d quantified values carry SD or range" % (
        n_sd, quant_n[sid])

    study_items[sid] = got
    study_detail[sid] = ev


POPS = [("all60", ALL), ("volatile57", VOLATILE), ("volatile_primary55", PRIMARY)]
census_rows = []
for key in reporting.ITEM_KEYS + ["_replicates_loose"]:
    row = {"item": key,
           "label": reporting.ITEM_LABEL.get(key, "Replication mentioned "
                                             "anywhere (loose upper bound)"),
           "group": reporting.ITEM_GROUP.get(key, "uncertainty"),
           "tier": reporting.ITEM_TIER.get(key, "(diagnostic)")}
    for name, pop in POPS:
        vals = [study_items[s][key] for s in pop]
        scored = [v for v in vals if v is not None]
        k = sum(1 for v in scored if v)
        row["%s_n" % name] = len(scored)
        row["%s_k" % name] = k
        row["%s_pct" % name] = round(pct(k, len(scored)), 1)
        row["%s_unscorable" % name] = len(vals) - len(scored)
    census_rows.append(row)

census_fields = ["item", "label", "group", "tier"]
for name, _pop in POPS:
    census_fields += ["%s_k" % name, "%s_n" % name, "%s_pct" % name,
                      "%s_unscorable" % name]
write_tsv(config.results_path("reporting_census.tsv"),
          census_rows, census_fields)

detail_rows = []
for st in studies:
    sid = st["study_id"]
    row = {"study_id": sid, "doi": st["doi"], "pmcid": st["pmcid"],
           "tier": st["tier"], "text_channel": channels[sid],
           "is_review": "yes" if reviews[sid][0] else "no",
           "population": ("non_volatile" if sid in reporting.NON_VOLATILE_STUDIES
                          else "volatile_review" if reviews[sid][0]
                          else "volatile_primary"),
           "n_samples_in_scope": st["n_samples_in_scope"],
           "n_observations": st["n_observations"]}
    for key in reporting.ITEM_KEYS + ["_replicates_loose"]:
        val = study_items[sid][key]
        row[key] = "" if val is None else ("yes" if val else "no")
        row[key + "_evidence"] = study_detail[sid][key]
    detail_rows.append(row)
detail_fields = (["study_id", "doi", "pmcid", "tier", "text_channel",
                  "is_review", "population", "n_samples_in_scope",
                  "n_observations"]
                 + [f for k in reporting.ITEM_KEYS + ["_replicates_loose"]
                    for f in (k, k + "_evidence")])
write_tsv(config.results_path("reporting_study_items.tsv"),
          detail_rows, detail_fields)


auto_country = {s for s in PRIMARY if study_items[s]["origin_country"]}
rel_country = {st["study_id"] for st in studies
               if st["study_id"] in PRIMARY and st["country"]}
print("\ncountry of origin: %d studies from the methods prose, %d from the "
      "study table, agreeing on %d; union %d, i.e. the item "
      "sits at %.0f-%.0f%% of %d"
      % (len(auto_country), len(rel_country),
         len(auto_country & rel_country), len(auto_country | rel_country),
         pct(len(auto_country & rel_country), len(PRIMARY)),
         pct(len(auto_country | rel_country), len(PRIMARY)), len(PRIMARY)))


prior_rows = list(prior_analytical.values())
print("\nThe same three items on three populations")
print("-" * 78)
print("  %-26s %-22s %-22s %s"
      % ("item", "%d cached full texts" % len(prior_rows),
         "%d studies" % len(ALL), "%d volatile primary" % len(PRIMARY)))
for label, field, item in (
        ("RI basis named", "study_ri_basis_named", "ri_basis_named"),
        ("internal standard named", "study_internal_standard",
         "internal_standard"),
        ("column stationary phase", "study_column_phase", "column_phase")):
    k142 = sum(1 for r in prior_rows if r.get(field))
    row = census_by_item_pre = [c for c in census_rows if c["item"] == item][0]
    print("  %-26s %3d/%d (%4.1f%%)      %3d/%d (%4.1f%%)      %3d/%d (%4.1f%%)"
          % (label, k142, len(prior_rows), pct(k142, len(prior_rows)),
             row["all60_k"], row["all60_n"], row["all60_pct"],
             row["volatile_primary55_k"], row["volatile_primary55_n"],
             row["volatile_primary55_pct"]))

print("\nCENSUS (volatile primary, n=%d)" % len(PRIMARY))
print("-" * 78)
for row in census_rows:
    print("  %-9s %-52s %3d/%-3d %5.1f%%" % (
        row["tier"][:9], row["label"][:52],
        row["volatile_primary55_k"], row["volatile_primary55_n"],
        row["volatile_primary55_pct"]))


ESSENTIAL = [k for k in reporting.ITEM_KEYS
             if reporting.ITEM_TIER[k] == "essential"]
completeness = {}
for sid in ALL:
    def frac(keys):
        vals = [study_items[sid][k] for k in keys]
        scored = [v for v in vals if v is not None]
        return (sum(1 for v in scored if v) / len(scored)) if scored else float("nan")
    completeness[sid] = {
        "essential": frac(ESSENTIAL),
        "all": frac(reporting.ITEM_KEYS),
        "poolability": frac(reporting.POOLABILITY_ITEMS),
    }


ISO_MARKERS = [
    ("citronellol", ["citronellol"], 20.0, 34.0),
    ("geraniol", ["geraniol"], 14.0, 22.0),
    ("nerol", ["nerol"], 5.0, 12.0),
    ("phenylethyl alcohol",
     ["phenylethyl alcohol", "2-phenylethanol", "phenethyl alcohol"], None, 2.5),
]

cache = json.load(open(config.data_path("compound_cache.json"),
                       encoding="utf-8"))
name_to_key = {}
for cname, info in cache.items():
    if isinstance(info, dict) and info.get("inchikey"):
        name_to_key.setdefault(cname.strip().lower(), info["inchikey"])
marker_keys = {}
for name, aliases, lo, hi in ISO_MARKERS:
    for alias in aliases:
        if alias.lower() in name_to_key:
            marker_keys[name] = (name_to_key[alias.lower()], lo, hi)
            break
print("\nmarker InChIKeys resolved: %s" % {k: v[0][:14] for k, v in marker_keys.items()})

samples_by_id = {s["sample_id"]: s for s in sample_rows}


def marker_medians(species_filter=None, material_filter=None):

    per = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in obs_rel:
        s = samples_by_id.get(r["sample_id"])
        if not s or s["in_composition_scope"] != "yes":
            continue
        if species_filter and species_filter not in (s["species"] or "").lower():
            continue
        if material_filter and s.get("material_type") != material_filter:
            continue
        if r["detection_status"] != "quantified" or not r["area_pct"]:
            continue
        for name, (key, _lo, _hi) in marker_keys.items():
            if r["inchikey"] == key:
                per[r["study_id"]][name].append(float(r["area_pct"]))
    return {sid: {m: statistics.median(v) for m, v in d.items()}
            for sid, d in per.items()}


dam = marker_medians(species_filter="damascena")


def iso_gap(marker, value):

    _key, lo, hi = marker_keys[marker]
    if hi is not None and value > hi:
        return math.log2(value / hi)
    if lo is not None and value < lo:
        return math.log2(lo / value)
    return 0.0


consensus = {}
for name in marker_keys:
    vals = [d[name] for d in dam.values() if name in d]
    if vals:
        consensus[name] = statistics.median(vals)
print("corpus consensus (median of study medians, R. x damascena): %s"
      % {k: round(v, 2) for k, v in consensus.items()})

dam_rows = []
for sid, d in sorted(dam.items()):
    gaps = [iso_gap(m, v) for m, v in d.items()]
    cons = [abs(math.log2(v / consensus[m])) for m, v in d.items() if consensus.get(m)]
    is_rev = reviews[sid][0]
    row = {
        "study_id": sid,
        "population": ("review" if is_rev else "primary"),
        "n_markers": len(d),
        "iso_gap_mean_log2": round(statistics.mean(gaps), 4),
        "iso_markers_inside": sum(1 for g in gaps if g == 0),
        "consensus_dist_mean_log2": round(statistics.mean(cons), 4) if cons else "",
        "completeness_essential": round(completeness[sid]["essential"], 4),
        "completeness_all": round(completeness[sid]["all"], 4),
        "completeness_poolability": round(completeness[sid]["poolability"], 4),
        "ri_and_is": "yes" if (study_items[sid]["ri_basis_named"]
                               and study_items[sid]["internal_standard"]) else "no",
        "column_polarity": analytic[sid]["study_column_polarity"] or "(unstated)",
        "quant_mode": analytic[sid]["study_quant_mode"] or "(unstated)",
        "quant_mode_stated": "yes" if study_items[sid]["quant_mode"] else "no",
    }
    for flag in ("ri_basis_named", "internal_standard", "response_factor",
                 "ri_reported", "column_phase", "oven_program",
                 "replicates_strict", "dispersion_reported", "origin_country",
                 "harvest_month", "ident_basis"):
        val = study_items[sid][flag]
        row[flag] = "" if val is None else ("yes" if val else "no")
    for m in marker_keys:
        row["median_" + m.replace(" ", "_")] = round(d[m], 3) if m in d else ""
    dam_rows.append(row)


sig = collections.defaultdict(list)
for row in dam_rows:
    key = tuple(row["median_" + m.replace(" ", "_")] for m in marker_keys)
    sig[key].append(row["study_id"])
dupes = {k: v for k, v in sig.items() if len(v) > 1}
for group in dupes.values():
    for sid in group[1:]:
        for row in dam_rows:
            if row["study_id"] == sid:
                row["duplicate_of"] = group[0]
print("byte-identical marker profiles (non-independent studies): %s"
      % [v for v in dupes.values()])


def rank_biserial(x, y):

    gt = sum(1 for a in x for b in y if a > b)
    lt = sum(1 for a in x for b in y if a < b)
    return (gt - lt) / (len(x) * len(y))


def cliffs_ci(x, y, n_boot=10000):

    xs, ys = np.asarray(x, float), np.asarray(y, float)
    out = []
    for _ in range(n_boot):
        a = RNG.choice(xs, len(xs), replace=True)
        b = RNG.choice(ys, len(ys), replace=True)
        out.append(rank_biserial(list(a), list(b)))
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def perm_p_diff(x, y, n=N_PERM):

    obs = abs(statistics.median(x) - statistics.median(y))
    pool = np.array(list(x) + list(y), float)
    k = len(x)
    hits = 0
    for _ in range(n):
        RNG.shuffle(pool)
        if abs(np.median(pool[:k]) - np.median(pool[k:])) >= obs - 1e-12:
            hits += 1
    return (hits + 1) / (n + 1)


comparability_rows = []


def min_achievable_p(n1, n2):

    if n1 < 1 or n2 < 1:
        return float("nan")
    best = stats.mannwhitneyu(list(range(n1)), list(range(n1, n1 + n2)),
                              alternative="two-sided", method="exact")
    return float(best.pvalue)


def test_split(label, rows, flag, outcome, note=""):
    a = [r[outcome] for r in rows if r[flag] == "yes"]
    b = [r[outcome] for r in rows if r[flag] == "no"]
    rec = {"analysis": label, "outcome": outcome, "grouping": flag,
           "n_yes": len(a), "n_no": len(b), "note": note,
           "p_floor_for_these_n": round(min_achievable_p(len(a), len(b)), 5)}
    if len(a) < 2 or len(b) < 2 or min_achievable_p(len(a), len(b)) > 0.05:
        rec["result"] = (
            "NOT TESTABLE: %d vs %d studies; the smallest two-sided exact "
            "p reachable is %.3f" % (len(a), len(b),
                                     min_achievable_p(len(a), len(b))))
        comparability_rows.append(rec)
        return rec
    rec["median_yes"] = round(statistics.median(a), 4)
    rec["median_no"] = round(statistics.median(b), 4)
    u = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact")
    rec["mannwhitney_U"] = float(u.statistic)
    rec["p_exact"] = round(float(u.pvalue), 5)
    rec["p_perm_median_diff"] = round(perm_p_diff(a, b), 5)
    rec["cliffs_delta"] = round(rank_biserial(a, b), 4)
    lo, hi = cliffs_ci(a, b)
    rec["cliffs_delta_ci95"] = "%.3f..%.3f" % (lo, hi)
    rec["hodges_lehmann_shift"] = round(
        float(np.median([p - q for p in a for q in b])), 4)
    rec["result"] = ("no detectable association" if u.pvalue > 0.05
                     else "association at p<=0.05 (uncorrected)")
    comparability_rows.append(rec)
    return rec


def test_corr(label, rows, xkey, ykey, note=""):
    x = [r[xkey] for r in rows]
    y = [r[ykey] for r in rows]
    rho, p = stats.spearmanr(x, y)

    xa = np.asarray(x, float)
    ya = np.asarray(y, float)
    hits = 0
    for _ in range(N_PERM):
        if abs(stats.spearmanr(xa, RNG.permutation(ya)).statistic) >= abs(rho) - 1e-12:
            hits += 1
    comparability_rows.append({
        "analysis": label, "outcome": ykey, "grouping": xkey,
        "n_yes": len(rows), "n_no": "",
        "spearman_rho": round(float(rho), 4),
        "p_asymptotic": round(float(p), 5),
        "p_perm": round((hits + 1) / (N_PERM + 1), 5),
        "result": ("no detectable association" if (hits + 1) / (N_PERM + 1) > 0.05
                   else "association at p<=0.05 (uncorrected)"),
        "note": note})


primary_dam = [r for r in dam_rows if r["population"] == "primary"]
indep_dam = [r for r in primary_dam if not r.get("duplicate_of")]
print("\nR. x damascena marker studies: %d total, %d primary, %d independent"
      % (len(dam_rows), len(primary_dam), len(indep_dam)))

SPLIT_FLAGS = ("ri_basis_named", "internal_standard", "ri_and_is",
               "response_factor", "ri_reported", "column_phase",
               "oven_program", "replicates_strict", "dispersion_reported",
               "origin_country", "harvest_month", "quant_mode_stated",
               "ident_basis")
for flag in SPLIT_FLAGS:
    test_split("damascena primary", primary_dam, flag, "iso_gap_mean_log2",
               note="n=%d studies" % len(primary_dam))
    test_split("damascena primary", primary_dam, flag,
               "consensus_dist_mean_log2", note="n=%d studies" % len(primary_dam))
for flag in ("ri_basis_named", "ri_and_is", "response_factor",
             "dispersion_reported"):
    test_split("damascena independent (duplicate dropped)", indep_dam, flag,
               "iso_gap_mean_log2", note="sensitivity analysis")

for xkey in ("completeness_essential", "completeness_all",
             "completeness_poolability"):
    test_corr("damascena primary", primary_dam, xkey, "iso_gap_mean_log2",
              note="n=%d studies" % len(primary_dam))
    test_corr("damascena primary", primary_dam, xkey,
              "consensus_dist_mean_log2", note="n=%d studies" % len(primary_dam))


tbl = [[0, 0], [0, 0]]
for r in primary_dam:
    i = 0 if r["ri_basis_named"] == "yes" else 1
    j = 0 if r["iso_markers_inside"] == r["n_markers"] else 1
    tbl[i][j] += 1
odds, pf = stats.fisher_exact(tbl)
comparability_rows.append({
    "analysis": "damascena primary", "outcome": "all markers inside ISO",
    "grouping": "ri_basis_named", "n_yes": tbl[0][0] + tbl[0][1],
    "n_no": tbl[1][0] + tbl[1][1], "odds_ratio": round(float(odds), 4),
    "p_asymptotic": round(float(pf), 5),
    "result": "Fisher exact 2x2: %s" % tbl,
    "note": "cells [[RI-named all-in, RI-named not-all-in], "
            "[no-RI all-in, no-RI not-all-in]]"})


fam = [(i, r) for i, r in enumerate(comparability_rows)
       if r.get("p_exact") is not None or r.get("p_perm") is not None]
pvals = [(i, float(r.get("p_exact", r.get("p_perm")))) for i, r in fam]
pvals.sort(key=lambda t: t[1])
m = len(pvals)
prev = 1.0
for rank in range(m, 0, -1):
    idx, p = pvals[rank - 1]
    prev = min(prev, p * m / rank)
    comparability_rows[idx]["q_bh"] = round(prev, 5)
    comparability_rows[idx]["survives_fdr"] = "yes" if prev <= 0.05 else "no"
print("\nmultiplicity: %d tests in the family; %d survive BH-FDR q<=0.05"
      % (m, sum(1 for r in comparability_rows if r.get("survives_fdr") == "yes")))

comp_fields = ["analysis", "outcome", "grouping", "n_yes", "n_no",
               "median_yes", "median_no", "mannwhitney_U", "p_exact",
               "p_perm_median_diff", "p_floor_for_these_n", "cliffs_delta",
               "cliffs_delta_ci95", "hodges_lehmann_shift", "spearman_rho",
               "p_asymptotic", "p_perm", "q_bh", "survives_fdr", "odds_ratio",
               "result", "note"]
write_tsv(config.results_path("reporting_vs_comparability.tsv"),
          comparability_rows, comp_fields)

dam_fields = (["study_id", "population", "duplicate_of", "n_markers"]
              + ["median_" + m.replace(" ", "_") for m in marker_keys]
              + ["iso_gap_mean_log2", "iso_markers_inside",
                 "consensus_dist_mean_log2", "completeness_essential",
                 "completeness_all", "completeness_poolability",
                 "ri_and_is", "column_polarity", "quant_mode"]
              + list(SPLIT_FLAGS))
write_tsv(config.results_path("reporting_damascena_studies.tsv"),
          dam_rows, dam_fields)

print("\nREPORTING vs COMPARABILITY")
print("-" * 78)
for r in comparability_rows:
    print("  %-34s %-26s n=%s/%s  %s" % (
        r["grouping"], r["outcome"][:26], r["n_yes"], r["n_no"],
        r.get("result", "")[:40]))


eo = marker_medians(material_filter="essential_oil")


eo_columns = collections.defaultdict(lambda: collections.defaultdict(list))
for r in obs_rel:
    s = samples_by_id.get(r["sample_id"])
    if (not s or s["in_composition_scope"] != "yes"
            or s.get("material_type") != "essential_oil"
            or r["detection_status"] != "quantified" or not r["area_pct"]):
        continue
    for name, (key, _lo, _hi) in marker_keys.items():
        if r["inchikey"] == key:
            eo_columns[r["sample_id"]][name].append(float(r["area_pct"]))
eo_col = {sid: {m: statistics.median(v) for m, v in d.items()}
          for sid, d in eo_columns.items()}
col_study = {sid: samples_by_id[sid]["study_id"] for sid in eo_col}


_EXTRACTION_ORDER = ["hydrodistillation", "steam_distillation",
                     "supercritical_co2", "microwave", "ultrasound",
                     "solvent_extraction", "HS-SPME", "SPME", "SAFE", "SDE",
                     "headspace", "TD"]


def _extraction_level(sid):
    named = set(study_detail[sid]["extraction_method"].split(";"))
    for name in _EXTRACTION_ORDER:
        if name in named:
            return name
    return ""


LEVEL_FIELDS = {
    "column_polarity": lambda sid: analytic[sid]["study_column_polarity"],
    "platform_detector": lambda sid: "+".join(sorted(
        reporting.detectors(texts[sid]["methods"] or texts[sid]["whole"]))),
    "quant_mode": lambda sid: (analytic[sid]["study_quant_mode"] or "").split(";")[0],
    "extraction_method": _extraction_level,
    "column_phase": lambda sid: analytic[sid]["study_column_phase"],
}

cost_rows = []


def spread_for(item, levels_of, pop_medians, pop_label, min_studies=3,
               key_of=None):

    groups = collections.defaultdict(dict)
    for sid, med in pop_medians.items():
        lev = (levels_of(key_of(sid) if key_of else sid) or "").strip()
        if not lev:
            continue
        groups[lev][sid] = med
    for marker in marker_keys:
        per_level = {lev: [m[marker] for m in d.values() if marker in m]
                     for lev, d in groups.items()}
        per_level = {k: v for k, v in per_level.items() if len(v) >= min_studies}
        rec = {"item": item, "contrast_kind": "between_levels",
               "population": pop_label, "marker": marker,
               "n_levels_with_data": len(per_level)}
        if len(per_level) < 2:
            rec["verdict"] = ("data do not allow: fewer than 2 levels with "
                              ">=%d reporting studies" % min_studies)
            cost_rows.append(rec)
            continue
        meds = {k: statistics.median(v) for k, v in per_level.items()}
        lo_lev = min(meds, key=meds.get)
        hi_lev = max(meds, key=meds.get)
        rec["levels"] = "; ".join("%s: %.2f (n=%d)" % (k, meds[k], len(per_level[k]))
                                  for k in sorted(meds, key=meds.get))
        rec["median_low"] = round(meds[lo_lev], 3)
        rec["median_high"] = round(meds[hi_lev], 3)
        rec["fold_spread"] = (round(meds[hi_lev] / meds[lo_lev], 2)
                              if meds[lo_lev] > 0 else "")
        rec["abs_spread_pct_points"] = round(meds[hi_lev] - meds[lo_lev], 3)
        rec["iso_gap_low"] = round(iso_gap(marker, meds[lo_lev]), 3)
        rec["iso_gap_high"] = round(iso_gap(marker, meds[hi_lev]), 3)
        if len(per_level) == 2:
            a, b = per_level[hi_lev], per_level[lo_lev]
            u = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact")
            rec["test"] = "Mann-Whitney exact"
            rec["p"] = round(float(u.pvalue), 5)
            rec["effect"] = round(rank_biserial(a, b), 3)
        else:
            h = stats.kruskal(*per_level.values())
            rec["test"] = "Kruskal-Wallis"
            rec["p"] = round(float(h.pvalue), 5)
            n = sum(len(v) for v in per_level.values())
            rec["effect"] = round(
                (float(h.statistic) - len(per_level) + 1) / (n - len(per_level)), 3)
        rec["verdict"] = ("levels differ by %.1fx" % rec["fold_spread"]
                          if rec["fold_spread"] else "levels differ")
        caveats = []
        if meds[lo_lev] < 0.5:
            caveats.append("fold spread unstable: the low group's median is "
                           "%.2f%%, so the ratio is dominated by a near-zero "
                           "denominator -- read the absolute difference"
                           % meds[lo_lev])
        if min(len(v) for v in per_level.values()) < 4:
            caveats.append("smallest level has %d units"
                           % min(len(v) for v in per_level.values()))
        rec["caveat"] = "; ".join(caveats)
        cost_rows.append(rec)


for item, fn in LEVEL_FIELDS.items():
    spread_for(item, fn, eo,
               "per-study medians, essential_oil relative-%% (n=%d studies)"
               % len(eo))
    spread_for(item, fn, eo_col,
               "per-column values, essential_oil relative-%% (n=%d columns)"
               % len(eo_col), key_of=col_study.get)


for item in ("internal_standard", "ri_basis_named", "response_factor",
             "ri_alkane_series", "replicates_strict", "dispersion_reported"):
    groups = {"yes": [], "no": []}
    for sid, med in eo.items():
        val = study_items[sid].get(item)
        if val is None:
            continue
        groups["yes" if val else "no"].append(med)
    for marker in marker_keys:
        a = [m[marker] for m in groups["yes"] if marker in m]
        b = [m[marker] for m in groups["no"] if marker in m]
        rec = {"item": item, "contrast_kind": "reported_vs_not",
               "population": "essential_oil relative-%", "marker": marker,
               "n_levels_with_data": (len(a) >= 3) + (len(b) >= 3)}
        if len(a) < 3 or len(b) < 3:
            rec["verdict"] = ("data do not allow: %d reporting vs %d silent "
                              "studies with this marker" % (len(a), len(b)))
            cost_rows.append(rec)
            continue
        rec["levels"] = "reported: %.2f (n=%d); silent: %.2f (n=%d)" % (
            statistics.median(a), len(a), statistics.median(b), len(b))
        rec["median_high"] = round(max(statistics.median(a), statistics.median(b)), 3)
        rec["median_low"] = round(min(statistics.median(a), statistics.median(b)), 3)
        rec["fold_spread"] = (round(rec["median_high"] / rec["median_low"], 2)
                              if rec["median_low"] > 0 else "")
        rec["abs_spread_pct_points"] = round(rec["median_high"] - rec["median_low"], 3)
        u = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact")
        rec["test"] = "Mann-Whitney exact"
        rec["p"] = round(float(u.pvalue), 5)
        rec["effect"] = round(rank_biserial(a, b), 3)
        rec["verdict"] = "reporters and non-reporters differ by %sx" % rec["fold_spread"]
        rec["caveat"] = ("fold spread unstable: low median %.2f%%"
                         % rec["median_low"] if rec["median_low"] < 0.5 else "")
        cost_rows.append(rec)


def marker_medians_by(field):
    per = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in obs_rel:
        s = samples_by_id.get(r["sample_id"])
        if (not s or s["in_composition_scope"] != "yes"
                or r["detection_status"] != "quantified" or not r["area_pct"]):
            continue
        lev = (s.get(field) or "").strip()
        if not lev or lev == "unknown":
            continue
        for name, (key, _lo, _hi) in marker_keys.items():
            if r["inchikey"] == key:
                per[(r["study_id"], lev)][name].append(float(r["area_pct"]))
    return {k: {m: statistics.median(v) for m, v in d.items()}
            for k, d in per.items()}


for field in ("material_type", "species"):
    cells = marker_medians_by(field)
    spread_for(field, lambda key: key[1], cells,
               "per (study, %s) medians over all in-scope relative-%% samples "
               "(%d cells from %d studies)"
               % (field, len(cells), len({k[0] for k in cells})),
               min_studies=2)


def composite_scores(cells):
    ref = {}
    for marker in marker_keys:
        vals = [d[marker] for d in cells.values() if d.get(marker)]
        if vals:
            ref[marker] = statistics.median(vals)
    out = {}
    for key, d in cells.items():
        devs = [abs(math.log2(v / ref[m])) for m, v in d.items()
                if ref.get(m) and v > 0]
        if devs:
            out[key] = statistics.mean(devs)
    return out, ref


def spread_composite(item, levels_of, cells, pop_label, min_units=3):
    scores, ref = composite_scores(cells)
    groups = collections.defaultdict(list)
    for key, score in scores.items():
        lev = (levels_of(key) or "").strip()
        if lev:
            groups[lev].append(score)
    groups = {k: v for k, v in groups.items() if len(v) >= min_units}
    rec = {"item": item, "contrast_kind": "between_levels_composite",
           "population": pop_label,
           "marker": "composite departure from the corpus median profile "
                     "(mean |log2 ratio| over %d markers)" % len(ref),
           "n_levels_with_data": len(groups)}
    if len(groups) < 2:
        rec["verdict"] = ("data do not allow: fewer than 2 levels with "
                          ">=%d units" % min_units)
        cost_rows.append(rec)
        return
    meds = {k: statistics.median(v) for k, v in groups.items()}
    rec["levels"] = "; ".join("%s: %.3f (n=%d)" % (k, meds[k], len(groups[k]))
                              for k in sorted(meds, key=meds.get))
    rec["median_low"] = round(min(meds.values()), 4)
    rec["median_high"] = round(max(meds.values()), 4)
    rec["abs_spread_pct_points"] = round(
        100 * (2 ** rec["median_high"] - 2 ** rec["median_low"]), 1)
    rec["fold_spread"] = (round(2 ** (rec["median_high"] - rec["median_low"]), 2)
                          if rec["median_low"] >= 0 else "")
    if len(groups) == 2:
        a, b = list(groups.values())
        u = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact")
        rec["test"] = "Mann-Whitney exact"
        rec["p"] = round(float(u.pvalue), 5)
        rec["effect"] = round(rank_biserial(a, b), 3)
    else:
        h = stats.kruskal(*groups.values())
        n = sum(len(v) for v in groups.values())
        rec["test"] = "Kruskal-Wallis"
        rec["p"] = round(float(h.pvalue), 5)
        rec["effect"] = round(
            (float(h.statistic) - len(groups) + 1) / (n - len(groups)), 3)
    rec["verdict"] = ("the levels of this item sit %.2fx apart on the "
                      "composite departure scale" % rec["fold_spread"])
    cost_rows.append(rec)


for field in ("material_type", "species"):
    cells = marker_medians_by(field)
    spread_composite(field, lambda key: key[1], cells,
                     "per (study, %s) cells, all in-scope relative-%% "
                     "(%d cells from %d studies)"
                     % (field, len(cells), len({k[0] for k in cells})),
                     min_units=2)
for item, fn in LEVEL_FIELDS.items():
    spread_composite(item, fn, eo,
                     "per-study, essential_oil relative-%% (n=%d studies)"
                     % len(eo), min_units=3)


study_compounds = collections.defaultdict(set)
for rows in (obs_rel, obs_abs, obs_oth):
    for r in rows:
        if r["compound_id"]:
            study_compounds[r["study_id"]].add(r["compound_id"])
compound_studies = collections.Counter()
for sid, cids in study_compounds.items():
    for cid in cids:
        compound_studies[cid] += 1
shared_frac = {sid: sum(1 for c in cids if compound_studies[c] >= 2) / len(cids)
               for sid, cids in study_compounds.items() if cids}

for item in ("ri_reported", "ri_basis_named", "ident_basis",
             "ri_alkane_series"):
    a = [shared_frac[s] for s in PRIMARY
         if s in shared_frac and study_items[s].get(item)]
    b = [shared_frac[s] for s in PRIMARY
         if s in shared_frac and study_items[s].get(item) is False]
    rec = {"item": item, "contrast_kind": "reported_vs_not",
           "population": "all primary volatile studies with resolved compounds",
           "marker": "share of a study's compounds also reported by "
                     "another study"}
    if len(a) < 3 or len(b) < 3:
        rec["verdict"] = ("data do not allow: %d reporting vs %d silent"
                          % (len(a), len(b)))
    else:
        rec["levels"] = "reported: %.3f (n=%d); silent: %.3f (n=%d)" % (
            statistics.median(a), len(a), statistics.median(b), len(b))
        rec["median_high"] = round(max(statistics.median(a), statistics.median(b)), 3)
        rec["median_low"] = round(min(statistics.median(a), statistics.median(b)), 3)
        rec["abs_spread_pct_points"] = round(
            100 * (rec["median_high"] - rec["median_low"]), 1)
        rec["fold_spread"] = (round(rec["median_high"] / rec["median_low"], 2)
                              if rec["median_low"] > 0 else "")
        u = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact")
        rec["test"] = "Mann-Whitney exact"
        rec["p"] = round(float(u.pvalue), 5)
        rec["effect"] = round(rank_biserial(a, b), 3)
        rec["verdict"] = ("studies reporting this item share %.0f%% of their "
                          "compounds with the rest of the corpus vs %.0f%%"
                          % (100 * statistics.median(a),
                             100 * statistics.median(b)))
    cost_rows.append(rec)


rsds = []
for r in obs_rel:
    if (r["detection_status"] == "quantified" and r["area_pct"]
            and r.get("area_pct_sd")):
        v, sd = float(r["area_pct"]), float(r["area_pct_sd"])
        if v > 0.1:
            rsds.append(100.0 * sd / v)
rsd_median = statistics.median(rsds) if rsds else float("nan")
rsd_p90 = float(np.percentile(rsds, 90)) if rsds else float("nan")
cost_rows.append({
    "item": "dispersion_reported", "contrast_kind": "within_study_repeatability",
    "population": "the %d relative-%% values that carry an SD" % len(rsds),
    "marker": "relative standard deviation of a single reported value",
    "median_low": round(rsd_median, 2), "median_high": round(rsd_p90, 2),
    "verdict": "median RSD %.1f%%, 90th percentile %.1f%%: a cross-study "
               "difference smaller than this is not interpretable"
               % (rsd_median, rsd_p90)})


PRACTICE = {
    "response_factor": lambda s: reporting.response_factor(
        texts[s]["methods"] or texts[s]["whole"])[0],
    "ri_basis_named": lambda s: analytic[s]["study_ri_basis_named"],
    "quant_mode": lambda s: (analytic[s]["study_quant_mode"] or "").split(";")[0],
    "internal_standard": lambda s: analytic[s]["study_internal_standard"],
    "ident_basis": lambda s: analytic[s]["study_ident_basis"],
    "extraction_method": _extraction_level,
    "material_type": lambda s: study_detail[s]["material_type"].split(": ")[-1],
}
for item, fn in PRACTICE.items():
    vals = collections.Counter(fn(s) for s in PRIMARY if fn(s))
    cost_rows.append({
        "item": item, "contrast_kind": "divergent_practice",
        "population": "volatile primary (n=%d)" % len(PRIMARY),
        "marker": "which convention the reporting studies used",
        "n_levels_with_data": len(vals),
        "levels": "; ".join("%s: %d" % (k, v) for k, v in vals.most_common()),
        "verdict": "%d distinct conventions among the %d studies that say; the "
                   "other %d studies could have used any of them"
                   % (len(vals), sum(vals.values()),
                      len(PRIMARY) - sum(vals.values()))})


EXPOSURE_KEY = {"column_polarity": "column_phase"}
exposure = {}
for key in reporting.ITEM_KEYS:
    vals = [study_items[s][key] for s in PRIMARY]
    scored = [v for v in vals if v is not None]
    exposure[key] = (1.0 - sum(1 for v in scored if v) / len(scored)
                     if scored else float("nan"))
for rec in cost_rows:
    rec["exposure_silent_pct"] = round(
        100.0 * exposure.get(EXPOSURE_KEY.get(rec["item"], rec["item"]),
                             float("nan")), 1)

cost_fields = ["item", "contrast_kind", "population", "marker",
               "n_levels_with_data", "levels", "median_low", "median_high",
               "fold_spread", "abs_spread_pct_points", "iso_gap_low",
               "iso_gap_high", "test", "p", "effect",
               "exposure_silent_pct", "verdict", "caveat"]
write_tsv(config.results_path("reporting_item_cost.tsv"),
          cost_rows, cost_fields)

print("\nCOST OF A MISSING ITEM (essential-oil relative-% population)")
print("-" * 78)
for rec in cost_rows:
    if rec.get("fold_spread"):
        print("  %-18s %-20s %-24s %5sx  p=%-8s silent=%4.1f%%" % (
            rec["item"], rec["marker"][:20], rec.get("levels", "")[:24],
            rec["fold_spread"], rec.get("p", ""), rec["exposure_silent_pct"]))


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


bridge_counts = {}
for axis in AXES:
    by_level = collections.defaultdict(set)
    unnamed = set()
    for s in inscope:
        val = level_of(s, axis)
        if val:
            by_level[val].add(s["study_id"])
        else:
            unnamed.add(s["study_id"])
    bridging = {k: v for k, v in by_level.items() if len(v) >= 2}
    pairs = sum(len(v) * (len(v) - 1) // 2 for v in bridging.values())
    bridge_counts[axis] = {
        "levels": len(by_level), "bridging_levels": len(bridging),
        "study_pairs": pairs,
        "studies_with_unnamed_samples": len(unnamed),
        "bridging_detail": "; ".join("%s(%d)" % (k, len(v))
                                     for k, v in sorted(bridging.items(),
                                                        key=lambda kv: -len(kv[1]))[:6]),
    }


rel_ids = {r["sample_id"] for r in obs_rel}
rel_only = [s for s in inscope if s["sample_id"] in rel_ids]
rel_bridges = {}
for axis in AXES:
    by_level = collections.defaultdict(set)
    for s in rel_only:
        val = level_of(s, axis)
        if val:
            by_level[val].add(s["study_id"])
    rel_bridges[axis] = sum(1 for v in by_level.values() if len(v) >= 2)

print("\nBRIDGING STRUCTURE (all in-scope samples; relative-% only in brackets)")
print("-" * 78)
for axis, d in bridge_counts.items():
    print("  %-16s levels=%3d  bridging=%2d [%2d rel-%%]  study-pairs=%3d  "
          "studies with unnamed samples=%2d" % (
              axis, d["levels"], d["bridging_levels"], rel_bridges[axis],
              d["study_pairs"], d["studies_with_unnamed_samples"]))
    if d["bridging_detail"]:
        print("      %s" % d["bridging_detail"])


def attr_measure(sid):
    kinds = {s["measure_kind"] for s in samples_by_study[sid]
             if s["in_composition_scope"] == "yes" and s["measure_kind"]}
    return ";".join(sorted(kinds))


def attr_material(sid):
    mats = {s["material_type"] for s in samples_by_study[sid]
            if s["in_composition_scope"] == "yes" and s["material_type"]}
    return ";".join(sorted(mats))


RUNGS = [
    ("measure_kind", attr_measure),
    ("material_type", attr_material),
    ("extraction_method", lambda s: study_detail[s]["extraction_method"]),
    ("column_polarity", lambda s: analytic[s]["study_column_polarity"]),
    ("quant_mode", lambda s: (analytic[s]["study_quant_mode"] or "").split(";")[0]),
    ("ri_basis_named", lambda s: analytic[s]["study_ri_basis_named"]),
]


def rung_status(a, b, upto):

    undecided = False
    for _name, get in RUNGS[:upto]:
        va, vb = get(a), get(b)
        if not va or not vb:
            undecided = True
            continue
        if va != vb:
            return "no"
    return "undecidable" if undecided else "yes"


all_pairs = []
for axis in AXES:
    by_level = collections.defaultdict(set)
    for s in inscope:
        val = level_of(s, axis)
        if val:
            by_level[val].add(s["study_id"])
    for level, sids in by_level.items():
        for a, b in itertools.combinations(sorted(sids), 2):
            all_pairs.append((axis, level, a, b))

print("\nCOMMENSURABILITY LADDER (cross-study pairs sharing a biological level)")
print("-" * 78)
print("  %-58s %5s %5s %5s %5s" % ("requirements", "pairs", "pool", "incomp",
                                   "undec"))
ladder_rows = []
for upto in range(1, len(RUNGS) + 1):
    counts = collections.Counter(rung_status(a, b, upto)
                                 for _ax, _lv, a, b in all_pairs)
    label = " + ".join(n for n, _g in RUNGS[:upto])
    print("  %-58s %5d %5d %5d %5d" % (
        label[:58], len(all_pairs), counts["yes"], counts["no"],
        counts["undecidable"]))
    ladder_rows.append({"requirements": label, "pairs": len(all_pairs),
                        "poolable": counts["yes"], "incompatible": counts["no"],
                        "undecidable": counts["undecidable"]})

FULL = len(RUNGS)
pair_status = collections.Counter()
for axis, _level, a, b in all_pairs:
    pair_status[(axis, rung_status(a, b, FULL))] += 1
print("\n  by axis, at the full requirement set:")
for axis in AXES:
    tot = sum(pair_status[(axis, s)] for s in ("yes", "no", "undecidable"))
    if not tot:
        continue
    print("  %-16s pairs=%3d   poolable=%3d   incompatible=%3d   "
          "undecidable=%3d" % (axis, tot, pair_status[(axis, "yes")],
                               pair_status[(axis, "no")],
                               pair_status[(axis, "undecidable")]))

tot_yes = sum(pair_status[(a, "yes")] for a in AXES)
tot_no = sum(pair_status[(a, "no")] for a in AXES)
tot_und = sum(pair_status[(a, "undecidable")] for a in AXES)


observed_levels = {}
for name, get in RUNGS:
    vals = [get(s) for s in ALL]
    observed_levels[name] = [v for v in vals if v]

sim_yes, sim_no = [], []
for _ in range(2000):
    filled = {}
    for name, get in RUNGS:
        pool_vals = observed_levels[name]
        filled[name] = {s: (get(s) or (pool_vals[RNG.integers(len(pool_vals))]
                                       if pool_vals else ""))
                        for s in ALL}
    yes = no = 0
    for _axis, _level, a, b in all_pairs:
        if all(filled[n][a] == filled[n][b] for n, _g in RUNGS):
            yes += 1
        else:
            no += 1
    sim_yes.append(yes)
    sim_no.append(no)
sim_mean = float(np.mean(sim_yes))
sim_lo, sim_hi = np.percentile(sim_yes, [2.5, 97.5])

print("\n  at the full requirement set: %d poolable, %d incompatible, "
      "%d undecidable (of %d pairs)" % (tot_yes, tot_no, tot_und,
                                        len(all_pairs)))
print("  counterfactual, every rung reported by every study: %.0f poolable "
      "pairs (95%% simulation interval %.0f-%.0f), i.e. %+.0f"
      % (sim_mean, sim_lo, sim_hi, sim_mean - tot_yes))
print("  -- the %d currently undecidable pairs are the reporting problem; the "
      "%.0f that stay incompatible are the design problem"
      % (tot_und, len(all_pairs) - sim_mean))


cult_by_study = collections.defaultdict(set)
for s in inscope:
    key = reporting.normalize_cultivar(s.get("cultivar"))
    if key:
        cult_by_study[s["study_id"]].add(key)
pool = sorted({c for v in cult_by_study.values() for c in v})
sizes = [len(v) for v in cult_by_study.values() if v]
observed_shared = sum(
    1 for c in pool
    if sum(1 for v in cult_by_study.values() if c in v) >= 2)

shared_null = []
for _ in range(N_PERM // 4):
    draws = [set(RNG.choice(len(pool), size=min(k, len(pool)), replace=False))
             for k in sizes]
    counts = collections.Counter(i for d in draws for i in d)
    shared_null.append(sum(1 for v in counts.values() if v >= 2))
exp_shared = float(np.mean(shared_null))
p_shared = (sum(1 for v in shared_null if v <= observed_shared) + 1) / (len(shared_null) + 1)
print("\n  cultivars named across %d studies: pool=%d, study sizes=%s"
      % (len(sizes), len(pool), sorted(sizes, reverse=True)[:10]))
print("  cultivars measured by >=2 studies: observed %d, expected under a "
      "shared-pool null %.1f (p=%.4f one-sided)"
      % (observed_shared, exp_shared, p_shared))

yield_rows = []
for axis, d in bridge_counts.items():
    tot = sum(pair_status[(axis, s)] for s in ("yes", "no", "undecidable"))
    yield_rows.append({
        "axis": axis, "levels": d["levels"],
        "bridging_levels": d["bridging_levels"],
        "cross_study_pairs": tot,
        "poolable_now": pair_status[(axis, "yes")],
        "incompatible": pair_status[(axis, "no")],
        "undecidable_missing_reporting": pair_status[(axis, "undecidable")],
        "bridging_levels_relative_pct_only": rel_bridges[axis],
        "studies_with_unnamed_samples": d["studies_with_unnamed_samples"],
        "bridging_detail": d["bridging_detail"]})
for lr in ladder_rows:
    yield_rows.append({"axis": "LADDER: " + lr["requirements"],
                       "cross_study_pairs": lr["pairs"],
                       "poolable_now": lr["poolable"],
                       "incompatible": lr["incompatible"],
                       "undecidable_missing_reporting": lr["undecidable"]})
yield_rows.append({
    "axis": "COUNTERFACTUAL: every rung reported",
    "cross_study_pairs": len(all_pairs),
    "poolable_now": round(sim_mean, 1),
    "incompatible": round(len(all_pairs) - sim_mean, 1),
    "undecidable_missing_reporting": 0,
    "bridging_detail": "95%% simulation interval %.0f-%.0f; missing values "
                       "drawn from the empirical distribution among reporters"
                       % (sim_lo, sim_hi)})
write_tsv(config.results_path("reporting_yield.tsv"), yield_rows,
          ["axis", "levels", "bridging_levels",
           "bridging_levels_relative_pct_only", "cross_study_pairs",
           "poolable_now", "incompatible", "undecidable_missing_reporting",
           "studies_with_unnamed_samples", "bridging_detail"])


summary = {
    "n_all": len(ALL), "n_volatile": len(VOLATILE), "n_primary": len(PRIMARY),
    "census": {r["item"]: r for r in census_rows},
    "exposure": exposure,
    "cost": cost_rows,
    "comparability": comparability_rows,
    "dam_n": len(primary_dam), "dam_indep_n": len(indep_dam),
    "bridges": bridge_counts,
    "ladder": ladder_rows,
    "pairs": {"yes": tot_yes, "no": tot_no, "undecidable": tot_und,
              "n_pairs": len(all_pairs), "counterfactual_mean": sim_mean,
              "counterfactual_ci": [float(sim_lo), float(sim_hi)]},
    "cultivar": {"pool": len(pool), "studies": len(sizes),
                 "observed_shared": observed_shared,
                 "expected_shared": exp_shared, "p": p_shared},
    "dupes": {k: v for k, v in enumerate(dupes.values())},
}
with open(config.results_path("reporting_summary.json"), "w",
          encoding="utf-8") as fh:
    json.dump(summary, fh, indent=1, default=str)
print("\n-> %s" % os.path.basename(config.results_path("reporting_summary.json")))
