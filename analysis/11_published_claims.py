import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from comparability import config, meta, species_analysis as spa, stats, separability
from comparability.contrast_analysis import write_tsv, fmt

import itertools
from collections import defaultdict

MEASUREMENTS = config.data_path("measurements_long.tsv")
CACHE = config.data_path("compound_cache.json")
MATERIALS = config.data_path("sample_material_types.tsv")
OUTDIR = config.RESULTS_DIR

SEED = 20260913
N_BOOT = 9999

FOCUS_SPECIES = ["Rosa alba", "Rosa gallica",
                 "Rosa x centifolia", "Rosa x damascena"]
FOCUS_MARKERS = ["citronellol", "geraniol", "nerol", "phenylethyl alcohol"]
DAMASCENA = "Rosa x damascena"


def short(sp):
    return sp.replace("Rosa x ", "R. x ").replace("Rosa ", "R. ")


def load():
    rows = spa.load_rows(MEASUREMENTS, "relative_pct")
    marker_keys = spa.resolve_marker_keys(CACHE, spa.MARKER_ALIASES)
    unresolved = [lab for lab, k in marker_keys if not k]
    values = meta.sample_marker_values(rows, marker_keys)

    sample_meta, sp_of = {}, {}
    for r in rows:
        s = (spa.study_id_of(r), r["table_id"], r["column_index"])
        sample_meta.setdefault(s, r)
        sp = spa.normalize_species(r.get("sample_species"))
        if sp:
            sp_of[s] = sp
    material, clash = spa.load_material_types(MATERIALS)

    import csv as _csv
    _csv.field_size_limit(10 ** 9)
    multi = {}
    with open(MEASUREMENTS, encoding="utf-8") as fh:
        for r in _csv.DictReader(fh, delimiter="\t"):
            if r.get("in_composition_scope") != "yes":
                continue
            sp = spa.normalize_species(r.get("sample_species"))
            if not sp:
                continue
            e = multi.setdefault(spa.study_id_of(r),
                                 {"species": set(), "kinds": set()})
            e["species"].add(sp)
            e["kinds"].add(r.get("measure_kind") or "")
    multi = {k: {"species": sorted(v["species"]),
                 "measure_kinds": sorted(v["kinds"])}
             for k, v in multi.items() if len(v["species"]) >= 2}

    return {
        "multi_species_studies": multi,
        "rows": rows, "marker_keys": marker_keys, "unresolved": unresolved,
        "values": values, "meta": sample_meta, "sp_of": sp_of,
        "material": material, "material_clash": clash,
        "samples": sorted(sample_meta),
    }


def species_samples(D, species, scope):
    out = []
    for s in D["samples"]:
        if D["sp_of"].get(s) != species:
            continue
        if scope != "all" and meta.material_of_sample(D["material"], s) != scope:
            continue
        out.append(s)
    return out


def candidate_iso(D, scope):

    samples = species_samples(D, DAMASCENA, scope)
    res = {"scope": scope, "markers": [], "cells": [],
           "n_samples": len(samples),
           "n_studies": len({s[0] for s in samples})}
    n_single = n_cells = 0
    per_study = defaultdict(lambda: [0, 0])
    for label, lo, hi in meta.ISO_PANEL:
        sm = meta.study_medians(D["values"], samples, label)
        if not sm:
            continue
        vals = sorted(sm.values())
        pooled = meta.pooled_study_median(sm)
        verdict = meta.interval_verdict(pooled, lo, hi)
        ci_lo, ci_hi = meta.bootstrap_median_ci(vals, N_BOOT, SEED)
        ci_inside = (meta.interval_verdict(ci_lo, lo, hi) == "within"
                     and meta.interval_verdict(ci_hi, lo, hi) == "within")
        base, flips = meta.leave_one_out_flips(sm, lo, hi)
        naive, n_naive = meta.pooled_sample_median(D["values"], samples, label)
        counts = {"within": 0, "below": 0, "above": 0}
        for st in sorted(sm):
            v = meta.interval_verdict(sm[st], lo, hi)
            counts[v] += 1
            n_cells += 1
            per_study[st][0] += 1 if v == "within" else 0
            per_study[st][1] += 1

            k = len([s for s in samples if label in D["values"].get(s, {})
                     and s[0] == st])
            if k == 1:
                n_single += 1
            res["cells"].append({
                "marker": label, "study": st, "value": sm[st],
                "n_samples": k, "verdict": v, "iso_lo": lo, "iso_hi": hi,
            })
        width = (hi - lo) if (lo is not None and hi is not None) else None
        spread = max(vals) - min(vals)
        vc_between, vc_within, vc_ratio, vc_k, vc_n = meta.variance_components(
            D["values"], samples, label)
        res["markers"].append({
            "var_between": vc_between, "var_within": vc_within,
            "var_ratio": vc_ratio, "var_k_studies": vc_k, "var_n_obs": vc_n,
            "marker": label, "iso_lo": lo, "iso_hi": hi,
            "pooled_study_median": pooled, "pooled_verdict": verdict,
            "pooled_sample_median": naive, "n_sample_obs": n_naive,
            "n_studies": len(sm), "n_within": counts["within"],
            "n_below": counts["below"], "n_above": counts["above"],
            "min": min(vals), "max": max(vals), "fold": meta.fold_range(vals),
            "boot_lo": ci_lo, "boot_hi": ci_hi, "boot_inside": ci_inside,
            "loso_flips": flips,
            "spread_over_iso_width": (spread / width) if width else None,
        })
    res["n_cells"] = n_cells
    res["n_single_column_cells"] = n_single
    res["n_pooled_within"] = sum(1 for m in res["markers"]
                                 if m["pooled_verdict"] == "within")
    res["n_markers"] = len(res["markers"])
    res["cells_within"] = sum(m["n_within"] for m in res["markers"])
    res["cells_below"] = sum(m["n_below"] for m in res["markers"])
    res["cells_above"] = sum(m["n_above"] for m in res["markers"])
    res["per_study"] = dict(per_study)
    res["studies_all_within"] = sorted(k for k, (a, b) in per_study.items()
                                       if a == b)
    res["n_markers_loso_stable"] = sum(1 for m in res["markers"]
                                       if not m["loso_flips"])
    res["n_markers_ci_inside"] = sum(1 for m in res["markers"]
                                     if m["boot_inside"])
    ratios = [m["var_ratio"] for m in res["markers"] if m["var_ratio"]]
    res["var_ratios"] = sorted(ratios)
    res["var_ratio_median"] = stats.median(ratios) if ratios else None
    res["var_ratio_markers"] = [(m["marker"], m["var_ratio"], m["var_k_studies"])
                                for m in res["markers"] if m["var_ratio"]]
    return res


def candidate_species_rank(D, marker, scope="all"):
    sm_by_sp, pooled, naive = {}, {}, {}
    for sp in FOCUS_SPECIES:
        ss = species_samples(D, sp, scope)
        sm = meta.study_medians(D["values"], ss, marker)
        if not sm:
            continue
        sm_by_sp[sp] = sm
        pooled[sp] = meta.pooled_study_median(sm)
        naive[sp] = meta.pooled_sample_median(D["values"], ss, marker)

    within = {}
    for st in sorted({k for sm in sm_by_sp.values() for k in sm}):
        d = {sp: sm_by_sp[sp][st] for sp in sm_by_sp if st in sm_by_sp[sp]}
        if len(d) >= 2:
            within[st] = d

    pairs = []
    for a, b in itertools.combinations(sorted(sm_by_sp), 2):
        pa, pb = pooled[a], pooled[b]
        ma, mb, common = meta.matched_study_medians(sm_by_sp[a], sm_by_sp[b])
        w = []
        for st, d in within.items():
            if a in d and b in d:
                w.append({"study": st, "a": d[a], "b": d[b],
                          "dir": ">" if d[a] > d[b] else
                                 ("<" if d[a] < d[b] else "="),
                          "ratio": (d[a] / d[b]) if d[b] > 0 else None})
        pdir = ">" if pa > pb else ("<" if pa < pb else "=")
        mdir = (">" if ma > mb else ("<" if ma < mb else "=")) if ma is not None else None
        wdirs = {x["dir"] for x in w}
        if not w:
            verdict = "untestable"
        elif wdirs == {pdir}:
            verdict = "reproduced"
        elif pdir in wdirs:
            verdict = "inconsistent"
        else:
            verdict = "contradicted"
        pairs.append({
            "a": a, "b": b, "pooled_a": pa, "pooled_b": pb,
            "pooled_dir": pdir, "pooled_ratio": (pa / pb) if pb else None,
            "matched_a": ma, "matched_b": mb, "matched_dir": mdir,
            "matched_k": len(common), "matched_studies": common,
            "within": w, "n_within": len(w), "verdict": verdict,
            "matched_reverses_naive": (mdir is not None and mdir != pdir),
        })

    pooled_rank = meta.rank_order(pooled)
    ranks_pooled = {sp: -pooled[sp] for sp in pooled}
    rank_rows = []
    for st, d in sorted(within.items()):
        if len(d) < 3:
            continue
        tau, c, dsc, n, keys = meta.kendall_tau(ranks_pooled,
                                                {sp: -v for sp, v in d.items()})
        rank_rows.append({"study": st, "order": meta.rank_order(d),
                          "tau_vs_pooled": tau, "concordant": c,
                          "discordant": dsc, "n_pairs": n, "keys": keys})
    tau_within = None
    if len(rank_rows) == 2:
        a = {sp: -v for sp, v in within[rank_rows[0]["study"]].items()}
        b = {sp: -v for sp, v in within[rank_rows[1]["study"]].items()}
        tau_within = meta.kendall_tau(a, b)[0]

    bal = {}
    if marker in ("geraniol", "citronellol"):
        num, den = ("geraniol", "citronellol")
        for sp in FOCUS_SPECIES:
            ss = species_samples(D, sp, scope)
            b, dropped = meta.balance_values(D["values"], ss, num, den)
            by = defaultdict(list)
            for s, v in b.items():
                by[s[0]].append(v)
            bal[sp] = {"per_study": {k: stats.median(v) for k, v in by.items()},
                       "dropped": dropped}
        bp = {sp: stats.median(list(bal[sp]["per_study"].values()))
              for sp in bal if bal[sp]["per_study"]}
        bal["_pooled"] = bp
        bal["_pooled_rank"] = meta.rank_order(bp)
        bal["_within"] = {}
        for st in within:
            d = {sp: bal[sp]["per_study"][st] for sp in bal
                 if not sp.startswith("_") and st in bal[sp]["per_study"]}
            if len(d) >= 2:
                bal["_within"][st] = d

    return {"marker": marker, "scope": scope, "study_medians": sm_by_sp,
            "pooled": pooled, "naive": naive, "pooled_rank": pooled_rank,
            "within": within, "pairs": pairs, "rank_rows": rank_rows,
            "tau_between_within_studies": tau_within, "balance": bal}


def candidate_origin(D, scope):
    samples = species_samples(D, DAMASCENA, scope)

    def group(s):
        r = D["meta"][s]
        o = (r.get("sample_origin") or "").strip() or \
            (r.get("study_region") or "").strip()
        return "Kazanlak (BG)" if "Kazanlak" in o else "other/unstated"

    grp = {s: group(s) for s in samples}
    studies = [s[0] for s in samples]
    labels = [grp[s] for s in samples]
    v = stats.cramers_v(labels, studies)
    sep = separability.design_separability(studies, labels, meta.UNSTATED)
    by_study = defaultdict(set)
    for s in samples:
        by_study[s[0]].add(grp[s])
    both = sorted(k for k, g in by_study.items() if len(g) > 1)

    out = []
    for label in FOCUS_MARKERS:
        row = {"marker": label}
        for g in ("Kazanlak (BG)", "other/unstated"):
            ss = [s for s in samples if grp[s] == g]
            sm = meta.study_medians(D["values"], ss, label)
            if not sm:
                row[g] = None
                continue
            vals = sorted(sm.values())
            row[g] = {"pooled": meta.pooled_study_median(sm), "k": len(sm),
                      "min": min(vals), "max": max(vals),
                      "fold": meta.fold_range(vals)}
        a, b = row.get("Kazanlak (BG)"), row.get("other/unstated")
        row["ratio"] = (a["pooled"] / b["pooled"]) if (a and b and b["pooled"]) else None
        row["within_group_fold_max"] = max(
            [x["fold"] for x in (a, b) if x and x["fold"]] or [None]) \
            if any(x and x["fold"] for x in (a, b)) else None
        out.append(row)
    kaz = sorted({s[0] for s in samples if grp[s] == "Kazanlak (BG)"})
    return {"scope": scope, "rows": out, "cramers_v": v,
            "separability": sep, "studies_with_both": both,
            "n_studies": len(by_study), "kazanlak_studies": kaz,
            "n_kazanlak_studies": len(kaz)}


def candidate_processing(D):
    mt = lambda s: meta.material_of_sample(D["material"], s)
    out = {"contrasts": []}
    for a, b in (("essential_oil", "solvent_extract"),
                 ("essential_oil", "absolute_or_concrete"),
                 ("essential_oil", "hydrosol"),
                 ("essential_oil", "headspace_volatiles")):
        ss = [s for s in D["samples"] if mt(s) in (a, b)]
        by_study = defaultdict(set)
        for s in ss:
            by_study[s[0]].add(mt(s))
        both = sorted(k for k, g in by_study.items() if len(g) > 1)
        markers = []
        for label in FOCUS_MARKERS:
            arms = {}
            for g in (a, b):
                sm = meta.study_medians(D["values"],
                                        [s for s in ss if mt(s) == g], label)
                if sm:
                    arms[g] = {"pooled": meta.pooled_study_median(sm),
                               "k": len(sm), "sm": sm}
            if len(arms) < 2:
                continue
            ma, mb, common = meta.matched_study_medians(arms[a]["sm"],
                                                        arms[b]["sm"])
            pdir = ">" if arms[a]["pooled"] > arms[b]["pooled"] else "<"
            w = []
            for st in both:
                va = meta.study_medians(D["values"],
                                        [s for s in ss if s[0] == st and mt(s) == a],
                                        label)
                vb = meta.study_medians(D["values"],
                                        [s for s in ss if s[0] == st and mt(s) == b],
                                        label)
                if va and vb:
                    xa, xb = va[st], vb[st]
                    w.append({"study": st, "a": xa, "b": xb,
                              "dir": ">" if xa > xb else ("<" if xa < xb else "="),
                              "ratio": (xa / xb) if xb > 0 else None})
            wd = {x["dir"] for x in w}
            markers.append({
                "marker": label,
                "pooled_a": arms[a]["pooled"], "k_a": arms[a]["k"],
                "pooled_b": arms[b]["pooled"], "k_b": arms[b]["k"],
                "pooled_dir": pdir,
                "pooled_ratio": (arms[a]["pooled"] / arms[b]["pooled"])
                                if arms[b]["pooled"] else None,
                "matched_k": len(common), "within": w,
                "verdict": ("untestable" if not w else
                            ("reproduced" if wd == {pdir} else
                             ("contradicted" if pdir not in wd else "inconsistent"))),
            })
        out["contrasts"].append({"a": a, "b": b, "studies_with_both": both,
                                 "markers": markers})
    return out


def candidate_correlation(D, top_n=20, min_pooled=20, min_study_n=4):
    cells, smeta, names, _u, _m, _d = spa.build_cells(D["rows"])
    samples = sorted(smeta)
    measured = spa.measured_sets(cells)
    per = defaultdict(set)
    for s in samples:
        for c in measured[s]:
            per[c].add(s)
    order = sorted(per, key=lambda c: (-len(per[c]), c))[:top_n]
    scan = meta.pair_correlation_scan(cells, samples, order,
                                      min_pooled, min_study_n)
    unanimous = [p for p in scan
                 if p["n_studies_flipping"] == p["n_studies_tested"]]
    strong = [p for p in scan if abs(p["r_pooled"]) >= 0.40]
    strong_flip = [p for p in strong
                   if p["n_studies_flipping"] > p["n_studies_tested"] / 2.0]
    sig_flip = [p for p in unanimous
                if p["p_pooled"] is not None and p["p_pooled"] <= 0.05]
    return {"names": names, "n_pairs": len(scan), "scan": scan,
            "unanimous": unanimous, "strong": strong,
            "strong_flip": strong_flip, "significant_unanimous": sig_flip,
            "compounds": order}


def write_false_conclusions_tsv(path, recs):
    fields = ["candidate_id", "family", "scope", "pooled_claim",
              "pooled_method", "pooled_result", "pooled_k_studies",
              "within_method", "within_k_studies", "within_result",
              "discrepancy", "verdict", "damage_rank", "airtightness_rank",
              "caveat"]
    write_tsv(path, fields, recs)


def write_figure_tsv(path, rows):
    fields = ["panel", "series", "marker", "group", "study_id",
              "material_type", "n_samples", "value", "iso_lo", "iso_hi",
              "verdict", "role", "note"]
    write_tsv(path, fields, rows)


def main():
    if not os.path.isdir(OUTDIR):
        os.makedirs(OUTDIR)
    D = load()

    iso_all = candidate_iso(D, "all")
    iso_oil = candidate_iso(D, "essential_oil")
    ranks = {m: candidate_species_rank(D, m, "all") for m in FOCUS_MARKERS}
    ranks_oil = {m: candidate_species_rank(D, m, "essential_oil")
                 for m in FOCUS_MARKERS}
    origin_all = candidate_origin(D, "all")
    origin_oil = candidate_origin(D, "essential_oil")
    proc = candidate_processing(D)
    corr = candidate_correlation(D)

    recs = build_records(D, iso_all, iso_oil, ranks, ranks_oil,
                         origin_all, origin_oil, proc, corr)
    write_false_conclusions_tsv(os.path.join(OUTDIR, "false_conclusions.tsv"),
                                recs)
    figrows = build_figure_rows(D, iso_oil, iso_all, ranks)
    write_figure_tsv(os.path.join(OUTDIR, "main_figure_data.tsv"), figrows)

    for _p in ("false_conclusions.tsv", "main_figure_data.tsv"):
        print("wrote %s" % os.path.basename(os.path.join(OUTDIR, _p)))
    print("candidates: %d rows" % len(recs))
    print("figure rows: %d" % len(figrows))


def build_records(D, iso_all, iso_oil, ranks, ranks_oil,
                  origin_all, origin_oil, proc, corr):
    recs = []

    for tag, iso in (("A1", iso_oil), ("A2", iso_all)):
        scope = iso["scope"]
        recs.append({
            "candidate_id": tag,
            "family": "ISO 9842 conformity of R. x damascena",
            "scope": scope,
            "pooled_claim": ("The published composition of R. x damascena oil "
                             "conforms to ISO 9842:2024 on every marker for "
                             "which the standard sets a limit."),
            "pooled_method": ("median of study medians (one vote per "
                              "laboratory), compared with the ISO interval"),
            "pooled_result": "%d/%d markers within interval" % (
                iso["n_pooled_within"], iso["n_markers"]),
            "pooled_k_studies": iso["n_studies"],
            "within_method": ("each (study x marker) cell is one laboratory's "
                              "own value for its own oil; no pooling"),
            "within_k_studies": iso["n_studies"],
            "within_result": ("%d/%d cells within (%.1f%%); %d below, %d above; "
                              "%d/%d studies within on every marker they reported"
                              % (iso["cells_within"], iso["n_cells"],
                                 100.0 * iso["cells_within"] / iso["n_cells"],
                                 iso["cells_below"], iso["cells_above"],
                                 len(iso["studies_all_within"]),
                                 len(iso["per_study"]))),
            "discrepancy": ("pooled verdict 'conforms' on %d/%d markers vs "
                            "%.0f%% of individual study cells outside their "
                            "interval; pooled verdict unchanged by any "
                            "single-study deletion on %d/%d markers"
                            % (iso["n_pooled_within"], iso["n_markers"],
                               100.0 * (iso["n_cells"] - iso["cells_within"])
                               / iso["n_cells"],
                               iso["n_markers_loso_stable"], iso["n_markers"])),
            "verdict": "FALSE (within-study evidence contradicts)",
            "damage_rank": 1 if scope == "essential_oil" else 2,
            "airtightness_rank": 1 if scope == "essential_oil" else 2,
            "caveat": (("the stratum ISO 9842 actually governs"
                        if scope == "essential_oil" else
                        "mixes oils with absolutes, hydrosols, headspace and "
                        "solvent extracts, which ISO 9842 does not govern; "
                        "record A1 is the governed stratum; this row shows only "
                        "that the verdict is not created by the material mix")
                       + ". Interval boundaries are inclusive and heneicosane's "
                         "pooled median sits exactly on its upper limit "
                         "(%s vs %s), which is also the only marker whose "
                         "pooled verdict is not leave-one-study-out stable. "
                         "The interval values are transcribed from open "
                         "secondary sources, not from the standard itself."
                       % (fmt([m for m in iso["markers"]
                               if m["marker"] == "heneicosane"][0]
                              ["pooled_study_median"], 3),
                          fmt([m for m in iso["markers"]
                               if m["marker"] == "heneicosane"][0]["iso_hi"], 1))),
        })

    for m, R in sorted(ranks.items()):
        contradicted = [p for p in R["pairs"] if p["verdict"] == "contradicted"]
        reproduced = [p for p in R["pairs"] if p["verdict"] == "reproduced"]
        inconsistent = [p for p in R["pairs"] if p["verdict"] == "inconsistent"]
        untestable = [p for p in R["pairs"] if p["verdict"] == "untestable"]
        rev = [p for p in R["pairs"] if p["matched_reverses_naive"]]
        if contradicted:
            verdict = "FALSE (within-study evidence contradicts)"
        elif inconsistent:
            verdict = "UNREPLICATED (within-study studies disagree with each other)"
        elif reproduced and not untestable:
            verdict = "SURVIVES (within-study evidence reproduces it)"
        else:
            verdict = "PARTIAL"
        recs.append({
            "candidate_id": "B-%s" % m.replace(" ", "_"),
            "family": "species ranking by a single marker",
            "scope": "all material types",
            "pooled_claim": ("Ranked by %s, the four European oil-bearing "
                             "Rosa species order as %s." %
                             (m, " > ".join(short(x) for x in R["pooled_rank"]))),
            "pooled_method": ("median of study medians per species, ranked; "
                              "arms are not matched (an indirect comparison)"),
            "pooled_result": "; ".join(
                "%s=%.3f (k=%d)" % (short(sp), R["pooled"][sp],
                                    len(R["study_medians"][sp]))
                for sp in R["pooled_rank"]),
            "pooled_k_studies": len({k for sm in R["study_medians"].values()
                                     for k in sm}),
            "within_method": ("rank order inside each study that measured >=2 "
                              "of the four species on one instrument"),
            "within_k_studies": len([s for s, d in R["within"].items()
                                     if len(d) >= 2]),
            "within_result": " | ".join(
                "%s: %s" % (rr["study"], " > ".join(short(x) for x in rr["order"]))
                for rr in R["rank_rows"]) or "no study ranks >=3 of them",
            "discrepancy": ("of 6 species pairs: %d reproduced, %d "
                            "contradicted, %d with the two laboratories "
                            "disagreeing with each other, %d untestable. "
                            % (len(reproduced), len(contradicted),
                               len(inconsistent), len(untestable))
                            + ("Contradicted: " + "; ".join(
                                "%s vs %s pooled %s (%.2fx) but within %s"
                                % (short(p["a"]), short(p["b"]), p["pooled_dir"],
                                   p["pooled_ratio"] or float("nan"),
                                   ", ".join("%s %s (%.2fx)"
                                             % (w["study"], w["dir"],
                                                w["ratio"] or float("nan"))
                                             for w in p["within"]))
                                for p in contradicted) if contradicted else "")
                            + ((" Laboratories disagree with each other on: "
                                + "; ".join("%s vs %s" % (short(q["a"]),
                                                          short(q["b"]))
                                            for q in inconsistent))
                               if inconsistent else "")),
            "verdict": verdict,
            "damage_rank": 2 if contradicted else 4,
            "airtightness_rank": 2 if len(contradicted) and
                                 all(len(p["within"]) >= 2 for p in contradicted)
                                 else 3,
            "caveat": ("the two studies carrying the within-laboratory test "
                       "(PMC11597806 essential oils, PMC8398789 subcritical "
                       "extracts) share an author and the same plant-material "
                       "institute, and hold one column per species, so no "
                       "significance test is possible and none is reported; "
                       "%d of 6 pairs reverse when the same pooling rule is "
                       "restricted to common-arm studies. Within-laboratory "
                       "k = %d here%s"
                       % (len(rev),
                          len([s for s, d in R["within"].items() if len(d) >= 2]),
                          ("; the closure-invariant ln(geraniol/citronellol) "
                           "the balance does not reverse; the reversal is "
                           "specific to the single-component unit"
                           if m == "geraniol" else
                           ("; species with k=1 pooled study have a 'pooled' "
                            "value that is one laboratory's number"
                            if any(len(sm) < 2 for sm in R["study_medians"].values())
                            else "")))),
        })

    for tag, O in (("C1", origin_oil), ("C2", origin_all)):
        biggest = max(O["rows"],
                      key=lambda r: abs((r["ratio"] or 1.0) - 1.0))
        recs.append({
            "candidate_id": tag,
            "family": "origin effect (Kazanlak, Bulgaria vs elsewhere)",
            "scope": O["scope"],
            "pooled_claim": ("Bulgarian (Kazanlak) R. x damascena oil differs "
                             "in marker composition from oil of other origins."),
            "pooled_method": ("median of study medians within each origin "
                              "group"),
            "pooled_result": "; ".join(
                "%s ratio Kazanlak/other = %s" % (
                    r["marker"], fmt(r["ratio"], 2)) for r in O["rows"]),
            "pooled_k_studies": O["n_studies"],
            "within_method": ("a study containing both a Kazanlak sample and "
                              "a non-Kazanlak sample"),
            "within_k_studies": len(O["studies_with_both"]),
            "within_result": "no such study exists",
            "discrepancy": ("largest between-group fold difference is %sx "
                            "(%s); on every marker the fold spread within the "
                            "Kazanlak group alone exceeds it (%s); Cramer's "
                            "V(origin, study) = %s, and 0 studies hold both "
                            "groups"
                            % (fmt(max(biggest["ratio"], 1.0 / biggest["ratio"]), 2),
                               biggest["marker"],
                               ", ".join("%s %sx" % (r["marker"],
                                                     fmt((r.get("Kazanlak (BG)")
                                                          or {}).get("fold"), 1))
                                         for r in O["rows"]
                                         if (r.get("Kazanlak (BG)") or {}).get("fold")),
                               fmt(O["cramers_v"], 3))),
            "verdict": "UNTESTABLE (no within-laboratory contrast exists)",
            "damage_rank": 5,
            "airtightness_rank": 5,
            "caveat": ("origin is a study-level field here; the pooled "
                       "difference is smaller than the between-laboratory "
                       "spread inside either group, so no claim large enough "
                       "to be wrong is on offer"),
        })

    for c in proc["contrasts"]:
        if not c["markers"]:
            continue
        verdicts = {mm["verdict"] for mm in c["markers"]}
        n_rep = sum(1 for mm in c["markers"] if mm["verdict"] == "reproduced")
        n_con = sum(1 for mm in c["markers"] if mm["verdict"] == "contradicted")
        k_both = len(c["studies_with_both"])
        if not k_both:
            v = "UNTESTABLE (no study contains both material types)"
        elif verdicts == {"reproduced"}:
            v = ("CONSISTENT BUT UNREPLICATED (direction reproduced on %d/%d "
                 "markers by k=%d study, one column per arm)"
                 % (n_rep, len(c["markers"]), k_both))
        elif n_con:
            v = ("CONTRADICTED on %d/%d markers, by k=%d study with one column "
                 "per arm -- too thin to call the pooled claim false"
                 % (n_con, len(c["markers"]), k_both))
        else:
            v = "MIXED (k=%d study)" % k_both
        recs.append({
            "candidate_id": "D-%s_vs_%s" % (c["a"], c["b"]),
            "family": "processing effect",
            "scope": "%s vs %s" % (c["a"], c["b"]),
            "pooled_claim": ("%s and %s of rose give systematically different "
                             "marker compositions." % (c["a"], c["b"])),
            "pooled_method": "median of study medians within each material arm",
            "pooled_result": "; ".join(
                "%s %s (%sx, k=%d vs %d)" % (
                    mm["marker"], mm["pooled_dir"],
                    fmt(mm["pooled_ratio"], 2), mm["k_a"], mm["k_b"])
                for mm in c["markers"]),
            "pooled_k_studies": len({s for mm in c["markers"]
                                     for s in (mm["k_a"], mm["k_b"])}),
            "within_method": "a study containing both material types",
            "within_k_studies": len(c["studies_with_both"]),
            "within_result": (
                "; ".join("%s: %s" % (mm["marker"],
                                      ", ".join("%s %s (%sx)" % (w["study"], w["dir"],
                                                                 fmt(w["ratio"], 2))
                                                for w in mm["within"]) or "none")
                          for mm in c["markers"])
                or "no within-study test"),
            "discrepancy": (("no within-laboratory test exists: 0 studies hold "
                             "both material types, so the pooled difference is "
                             "100% between-laboratory")
                            if not k_both else
                            ("direction reproduced on %d/%d markers, "
                             "contradicted on %d/%d, by the %d study holding "
                             "both arms"
                             % (n_rep, len(c["markers"]), n_con,
                                len(c["markers"]), k_both))),
            "verdict": v,
            "damage_rank": 4,
            "airtightness_rank": 4,
            "caveat": ("the entire within-laboratory evidence base is %d "
                       "study/studies with one column per arm; the pooled arms "
                       "draw on disjoint sets of laboratories, so the pooled "
                       "ratio is an indirect comparison"
                       % k_both),
        })

    best = corr["strong_flip"][0] if corr["strong_flip"] else (
        corr["unanimous"][0] if corr["unanimous"] else None)
    recs.append({
        "candidate_id": "E1",
        "family": "compound-compound correlation (Simpson's-paradox shape)",
        "scope": "relative_pct, top-%d compounds by sample coverage"
                 % len(corr["compounds"]),
        "pooled_claim": ("Two compounds are positively (or negatively) "
                         "associated across the published literature."),
        "pooled_method": ("Pearson r on raw area % over all samples, ignoring "
                          "study -- the naive pooled correlation"),
        "pooled_result": "%d compound pairs with n>=20 and >=1 testable study"
                         % corr["n_pairs"],
        "pooled_k_studies": len({p["study"] for pp in corr["scan"]
                                 for p in pp["per_study"]}),
        "within_method": ("per-study Pearson r on studies with n>=4 samples, "
                          "reported as a tally of signs -- never pooled into "
                          "one test"),
        "within_k_studies": max([p["n_studies_tested"] for p in corr["scan"]]
                                or [0]),
        "within_result": ("%d/%d pairs have every testable study flipping the "
                          "pooled sign; %d pairs have |r_pooled|>=0.40, of "
                          "which %d are flipped by a majority of studies"
                          % (len(corr["unanimous"]), corr["n_pairs"],
                             len(corr["strong"]), len(corr["strong_flip"]))),
        "discrepancy": (("strongest case: %s ~ %s, pooled r=%s (p=%s, n=%d), "
                         "per-study r = %s"
                         % (corr["names"].get(best["a"], best["a"]),
                            corr["names"].get(best["b"], best["b"]),
                            fmt(best["r_pooled"], 3), fmt(best["p_pooled"], 4),
                            best["n_pooled"],
                            ", ".join("%s%s(n%d)" % (fmt(x["r"], 2), "*" if x["flip"] else "", x["n"])
                                      for x in best["per_study"])))
                        if best else "no pair qualifies"),
        "verdict": ("WEAK -- reversals are real but confined to pairs whose "
                    "pooled correlation nobody would report (%d of %d "
                    "unanimous reversals reaches p<=0.05, and that one has "
                    "pooled Spearman rho = -0.04, i.e. it is one extreme "
                    "point, not a monotone association)"
                    % (len(corr["significant_unanimous"]),
                       len(corr["unanimous"]))),
        "damage_rank": 3,
        "airtightness_rank": 3,
        "caveat": ("correlations between parts of a closed composition carry a "
                   "built-in negative bias, so neither the pooled nor the "
                   "within-study r is a statement about co-regulation; this "
                   "candidate is reported for completeness, not for the figure"),
    })
    return recs


def build_figure_rows(D, iso_oil, iso_all, ranks):
    rows = []

    for m in iso_oil["markers"]:
        rows.append({
            "panel": "A", "series": "pooled", "marker": m["marker"],
            "group": DAMASCENA, "study_id": "(pooled)",
            "material_type": "essential_oil", "n_samples": "",
            "value": fmt(m["pooled_study_median"], 3),
            "iso_lo": fmt(m["iso_lo"], 3) if m["iso_lo"] is not None else "",
            "iso_hi": fmt(m["iso_hi"], 3),
            "verdict": m["pooled_verdict"], "role": "pooled_estimate",
            "note": ("median of %d study medians; bootstrap95=[%s,%s]; "
                     "spread/ISO-width=%s; LOSO-stable=%s"
                     % (m["n_studies"], fmt(m["boot_lo"], 3), fmt(m["boot_hi"], 3),
                        fmt(m["spread_over_iso_width"], 2),
                        "yes" if not m["loso_flips"] else
                        "no (" + ";".join(m["loso_flips"]) + ")")),
        })
    for c in iso_oil["cells"]:
        rows.append({
            "panel": "A", "series": "study_cell", "marker": c["marker"],
            "group": DAMASCENA, "study_id": c["study"],
            "material_type": "essential_oil", "n_samples": c["n_samples"],
            "value": fmt(c["value"], 3),
            "iso_lo": fmt(c["iso_lo"], 3) if c["iso_lo"] is not None else "",
            "iso_hi": fmt(c["iso_hi"], 3),
            "verdict": c["verdict"], "role": "per_study_observation",
            "note": "",
        })

    for marker in ("geraniol", "citronellol"):
        R = ranks[marker]
        for sp in R["pooled_rank"]:
            rows.append({
                "panel": "B", "series": "pooled", "marker": marker,
                "group": sp, "study_id": "(pooled)", "material_type": "all",
                "n_samples": R["naive"][sp][1],
                "value": fmt(R["pooled"][sp], 3),
                "iso_lo": "", "iso_hi": "",
                "verdict": "rank %d of %d" % (R["pooled_rank"].index(sp) + 1,
                                              len(R["pooled_rank"])),
                "role": "pooled_estimate",
                "note": "median of %d study medians; sample-level median %s"
                        % (len(R["study_medians"][sp]),
                           fmt(R["naive"][sp][0], 3)),
            })
        for st, d in sorted(R["within"].items()):
            order = meta.rank_order(d)
            for sp in order:
                s0 = [s for s in species_samples(D, sp, "all") if s[0] == st]
                mt = meta.material_of_sample(D["material"], s0[0]) if s0 else ""
                rows.append({
                    "panel": "B", "series": "within_study", "marker": marker,
                    "group": sp, "study_id": st, "material_type": mt,
                    "n_samples": len(s0), "value": fmt(d[sp], 3),
                    "iso_lo": "", "iso_hi": "",
                    "verdict": "rank %d of %d" % (order.index(sp) + 1, len(order)),
                    "role": "within_laboratory_observation", "note": "",
                })

    R = ranks["geraniol"]
    for p in R["pairs"]:
        rows.append({
            "panel": "C", "series": "pair", "marker": "geraniol",
            "group": "%s vs %s" % (short(p["a"]), short(p["b"])),
            "study_id": "(pooled)", "material_type": "all",
            "n_samples": p["n_within"],
            "value": fmt(p["pooled_ratio"], 3),
            "iso_lo": "", "iso_hi": "",
            "verdict": p["verdict"], "role": "unmatched_pooled_ratio",
            "note": "pooled %s; matched (k=%d) %s; within: %s"
                    % (p["pooled_dir"], p["matched_k"], p["matched_dir"],
                       ", ".join("%s %s" % (w["study"], w["dir"])
                                 for w in p["within"]) or "none"),
        })
        if p["matched_dir"]:
            rows.append({
                "panel": "C", "series": "pair", "marker": "geraniol",
                "group": "%s vs %s" % (short(p["a"]), short(p["b"])),
                "study_id": "(matched pool, k=%d)" % p["matched_k"],
                "material_type": "all", "n_samples": p["matched_k"],
                "value": fmt((p["matched_a"] / p["matched_b"])
                             if p["matched_b"] else None, 3),
                "iso_lo": "", "iso_hi": "",
                "verdict": ("reverses_naive" if p["matched_reverses_naive"]
                            else "agrees_with_naive"),
                "role": "common_arm_pooled_ratio",
                "note": "studies: " + ";".join(p["matched_studies"]),
            })
        for w in p["within"]:
            rows.append({
                "panel": "C", "series": "pair", "marker": "geraniol",
                "group": "%s vs %s" % (short(p["a"]), short(p["b"])),
                "study_id": w["study"], "material_type": "all",
                "n_samples": 1, "value": fmt(w["ratio"], 3),
                "iso_lo": "", "iso_hi": "", "verdict": w["dir"],
                "role": "within_laboratory_ratio", "note": "",
            })
    return rows


if __name__ == "__main__":
    main()
