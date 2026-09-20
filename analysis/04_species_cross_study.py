import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import csv
import time
from collections import Counter, defaultdict

from comparability import config, stats, separability, species_analysis as SA


IN_TSV = config.data_path("measurements_long.tsv")
IN_MATERIAL = config.data_path("sample_material_types.tsv")
IN_CACHE = config.data_path("compound_cache.json")

OUT_MARKERS = config.results_path("species_marker_profiles.tsv")
OUT_TESTS = config.results_path("species_cross_study_tests.tsv")
OUT_DESIGN = config.results_path("species_design.tsv")

N_PERM = int(os.environ.get("ROSE_N_PERMUTATIONS", "9999"))
SEED = 20260913
UNSTATED = SA.UNSTATED

TESTS = []
T0 = time.time()


def log(msg):
    print("[%6.1fs] %s" % (time.time() - T0, msg))
    sys.stdout.flush()


def fnum(x, nd=3):
    if x is None:
        return ""
    if isinstance(x, float):
        if x != x:
            return "nan"
        if x in (float("inf"), float("-inf")):
            return "inf"
        return ("%." + str(nd) + "f") % x
    return str(x)


def rec(**kw):
    TESTS.append(kw)
    return kw


GLOSS = {
    "residual df = 0 (each study x level cell holds "
    "exactly 1 sample); not testable":
        "residual df = 0 -- every (study x level) cell holds exactly one "
        "sample, so there is no replication to test against",
    "no single study contains >= 2 levels -- "
    "structurally inseparable":
        "no single study contains two levels of this factor -- structurally "
        "inseparable from study identity",
    "the factor does not vary within any single study -- "
    "completely confounded with study; no conditional "
    "test is possible":
        "the factor does not vary inside any single study -- completely "
        "confounded with study",
    "only 1 level; not testable": "only one level; nothing to test",
    "within-group sum of squares is 0; F is undefined "
    "(samples coincide exactly)":
        "within-group SS = 0, F undefined",
}


log("loading %s" % os.path.basename(IN_TSV))
rows = SA.load_rows(IN_TSV, "relative_pct")
cells, meta, names, units, n_multi, n_dropped = SA.build_cells(rows)
meas = SA.measured_sets(cells)
sp_of = {s: SA.normalize_species(meta[s].get("sample_species")) for s in meta}
material, n_clash = SA.load_material_types(IN_MATERIAL)

all_samples = sorted(meta)
sp_samples = [s for s in all_samples if sp_of.get(s)]
log("rows=%d  samples=%d  studies=%d  skeletons=%d  species-stated samples=%d"
    % (len(rows), len(all_samples), len({s[0] for s in all_samples}),
       len({k[1] for k in cells}), len(sp_samples)))


name_fold = defaultdict(set)
for s in all_samples:
    raw = (meta[s].get("sample_species") or "").strip()
    if raw:
        name_fold[SA.normalize_species(raw)].add(raw)
folded = {k: sorted(v) for k, v in name_fold.items() if len(v) > 1}

cell_status = Counter(v[0] for v in cells.values())


mat_of_sample = {s: (material.get((s[1], s[2])) or "(not classified)")
                 for s in all_samples}
mat_counts = Counter(mat_of_sample[s] for s in sp_samples)

SCOPES = [
    ("all", "all material types", sp_samples),
    ("essential_oil", "material_type == essential_oil",
     [s for s in sp_samples if mat_of_sample[s] == "essential_oil"]),
]


designs = {}
for tag, _label, samples in SCOPES:
    designs[tag] = SA.species_design(samples, sp_of)

with open(OUT_DESIGN, "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["material_scope", "species", "status", "n_samples", "n_studies",
                "n_compound_skeletons_measured", "studies"])
    for tag, _label, samples in SCOPES:
        for d in designs[tag]:
            skel = set()
            for s in d["samples"]:
                skel |= meas.get(s, set())
            w.writerow([tag, d["species"], d["status"], d["n_samples"],
                        d["n_studies"], len(skel), ";".join(d["studies"])])
log("design written -> %s" % os.path.basename(OUT_DESIGN))

bridging_all = [d for d in designs["all"] if d["bridging"]]
log("bridging species (relative_pct, all materials): %d -- %s"
    % (len(bridging_all), ", ".join(d["species"] for d in bridging_all)))


csv.field_size_limit(10 ** 9)
with open(IN_TSV, encoding="utf-8") as fh:
    _span = defaultdict(set)
    _mk = defaultdict(set)
    for r in csv.DictReader(fh, delimiter="\t"):
        if r.get("in_composition_scope") != "yes":
            continue
        sp = SA.normalize_species(r.get("sample_species"))
        if sp:
            _span[sp].add(SA.study_id_of(r))
            _mk[sp].add(r.get("measure_kind"))
span_species = {sp: sorted(st) for sp, st in _span.items()}
span_bridging = sorted(sp for sp, st in _span.items() if len(st) >= 2)
lost_bridges = [sp for sp in span_bridging
                if sp not in {d["species"] for d in bridging_all}]
log("in scope, any measure_kind: %d species, %d bridging; lost to the "
    "relative_pct restriction: %s" % (len(_span), len(span_bridging),
                                      ", ".join(lost_bridges) or "none"))


LADDERS = {}
BLOCKS = {}
for tag, _label, samples in SCOPES:
    if len(samples) < 4:
        continue
    lf = SA.freq_ladder(meas, samples, sp_of, max_m=14)
    lb = SA.bridge_ladder(meas, samples, sp_of, max_m=12)
    LADDERS[tag] = {"freq": lf, "bridge": lb}
    picks = SA.pick_blocks(lb)
    BLOCKS[tag] = picks
    log("%s: bridge ladder %d steps; blocks %s" %
        (tag, len(lb), {k: "m=%d n=%d" % (v["m"], v["n_samples"])
                        for k, v in sorted(picks.items())}))


def factor_vals(samples, which):
    if which == "species":
        return [sp_of.get(s) or UNSTATED for s in samples]
    if which == "study":
        return [s[0] for s in samples]
    raise ValueError(which)


HEAD = {}
PERSPEC = []


def run_block(scope, btag, block):
    comps = block["compounds"]
    keep, raw, closed, clrs, n_capped = SA.block_vectors(
        cells, block["samples"], comps)
    if len(keep) < 4:
        log("  %s/%s: only %d complete rows, skipped" % (scope, btag, len(keep)))
        return
    spv = factor_vals(keep, "species")
    stv = factor_vals(keep, "study")
    v = stats.cramers_v(spv, stv)
    bs = SA.block_stats(keep, sp_of)
    log("  %s/%s m=%d n=%d studies=%d species=%d bridging=%d V=%.2f"
        % (scope, btag, len(comps), len(keep), bs["n_studies"],
           bs["n_species"], bs["n_bridging_species"], v))

    spaces = [("aitchison", SA.aitchison_d2(clrs)),
              ("bray_curtis", SA.bray_d2(closed))]
    for space, d2 in spaces:
        t = time.time()
        base = dict(material_scope=scope, block=btag, space=space,
                    m_compounds=len(comps), n=len(keep),
                    n_studies=bs["n_studies"], n_species=bs["n_species"],
                    n_bridging_species=bs["n_bridging_species"],
                    cramers_v_species_study=v, n_perm=N_PERM, seed=SEED)

        m_sp = stats.permanova(d2, spv, N_PERM, SEED)
        rec(mode="marginal", factor="species", r2=m_sp["r2"],
            r2_excess=m_sp["r2_excess"], r2_adj=m_sp["r2_adj"], F=m_sp["F"],
            p=m_sp["p"], df=m_sp["df_among"], testable=m_sp["testable"],
            reason=m_sp["reason"], **base)
        m_st = stats.permanova(d2, stv, N_PERM, SEED)
        rec(mode="marginal", factor="study", r2=m_st["r2"],
            r2_excess=m_st["r2_excess"], r2_adj=m_st["r2_adj"], F=m_st["F"],
            p=m_st["p"], df=m_st["df_among"], testable=m_st["testable"],
            reason=m_st["reason"], **base)

        nest = stats.permanova_conditional(d2, spv, stv, N_PERM, SEED)
        sst = nest["ss_total"]
        part = {
            "r2_species": (nest["ss_study"] / sst) if sst else None,
            "r2_study_given_species": (nest["ss_factor_given_study"] / sst)
                                      if sst else None,
            "r2_residual": (nest["ss_residual"] / sst) if sst else None,
            "F": nest["F"], "p": nest["p"], "testable": nest["testable"],
            "reason": nest["reason"], "df": nest["df_factor"],
            "marginal_species_r2": m_sp["r2"],
            "marginal_species_excess": m_sp["r2_excess"],
            "marginal_study_r2": m_st["r2"],
            "n": len(keep), "m": len(comps),
            "n_bridging_species": bs["n_bridging_species"],
            "n_species": bs["n_species"],
            "n_studies": bs["n_studies"], "cramers_v": v,
            "compounds": comps,
        }
        HEAD[(scope, btag, space)] = part
        rec(mode="nested_study_within_species", factor="study|species",
            r2=part["r2_study_given_species"], r2_species=part["r2_species"],
            r2_residual=part["r2_residual"], F=nest["F"], p=nest["p"],
            df=nest["df_factor"], testable=nest["testable"],
            reason=nest["reason"], **base)

        rc = separability.restricted_conditional(d2, spv, stv, N_PERM, SEED)
        rec(mode="restricted_study_within_bridging_species", factor="study|species",
            r2=rc.get("r2_partial"), F=rc.get("F"), p=rc.get("p"),
            df=rc.get("df_factor"), n_restricted=rc.get("n"),
            n_cells_ge2=rc.get("n_cells_ge2"),
            testable=rc.get("testable"), reason=rc.get("reason"), **base)

        revc = stats.permanova_conditional(d2, stv, spv, N_PERM, SEED)
        rec(mode="conditional_species_given_study", factor="species|study",
            r2=revc["r2_partial"], F=revc["F"], p=revc["p"],
            df=revc["df_factor"], testable=revc["testable"],
            reason=revc["reason"], **base)
        revr = separability.restricted_conditional(d2, stv, spv, N_PERM, SEED)
        rec(mode="restricted_species_within_study", factor="species|study",
            r2=revr.get("r2_partial"), F=revr.get("F"), p=revr.get("p"),
            df=revr.get("df_factor"), n_restricted=revr.get("n"),
            n_cells_ge2=revr.get("n_cells_ge2"),
            testable=revr.get("testable"), reason=revr.get("reason"), **base)

        by_sp = defaultdict(list)
        for i, sp in enumerate(spv):
            by_sp[sp].append(i)
        for sp, idx in sorted(by_sp.items()):
            nst = len({stv[i] for i in idx})
            if sp == UNSTATED:
                continue
            if nst < 2 or len(idx) < 4:
                rec(mode="within_species_study_effect", factor="study",
                    species=sp, r2=None, F=None, p=None, testable=False,
                    reason=("untestable: %d sample(s) from %d study(ies)"
                            % (len(idx), nst)), **base)
                continue
            sub = [[d2[i][j] for j in idx] for i in idx]
            r = stats.permanova(sub, [stv[i] for i in idx], N_PERM, SEED)
            rec(mode="within_species_study_effect", factor="study", species=sp,
                r2=r["r2"], r2_excess=r["r2_excess"], r2_adj=r["r2_adj"],
                F=r["F"], p=r["p"], df=r["df_among"], n_species_samples=len(idx),
                n_species_studies=nst, testable=r["testable"],
                reason=r["reason"], **base)
            PERSPEC.append(dict(scope=scope, block=btag, space=space, species=sp,
                                n=len(idx), n_studies=nst, r2=r["r2"],
                                r2_excess=r["r2_excess"], p=r["p"]))
        log("    %s done (%.0fs)" % (space, time.time() - t))


for tag, _label, samples in SCOPES:
    for btag in sorted(BLOCKS.get(tag, {})):
        run_block(tag, btag, BLOCKS[tag][btag])


log("resolving marker InChIKeys from %s" % os.path.basename(IN_CACHE))
marker_keys = SA.resolve_marker_keys(IN_CACHE, SA.MARKER_ALIASES)
key_of = dict(marker_keys)
unresolved = [lab for lab, k in marker_keys if not k]
iso_of = {}
for name, aliases, lo, hi in SA.ISO_RANGES:
    iso_of[name] = (lo, hi)
for lab, k in marker_keys:
    log("   %-22s %s" % (lab, k or "UNRESOLVED"))

prof = SA.marker_profiles(rows, sp_of, marker_keys, material)

status_of = {}
for tag, _l, _s in SCOPES:
    for d in designs[tag]:
        status_of[(tag, d["species"])] = d["status"]

ISO_SPECIES = "Rosa x damascena"


def iso_verdict(species, label, med):
    lo, hi = iso_of.get(label, (None, None))
    if lo is None and hi is None:
        return "", "", "no ISO value"
    if species != ISO_SPECIES:
        return (fnum(lo) if lo is not None else "<=", fnum(hi),
                "outside ISO scope (standard describes %s)" % ISO_SPECIES)
    if lo is None:
        return "<=", fnum(hi), ("within" if med <= hi else "above")
    if med < lo:
        return fnum(lo), fnum(hi), "below"
    if med > hi:
        return fnum(lo), fnum(hi), "above"
    return fnum(lo), fnum(hi), "within"


MARKER_ORDER = [lab for lab, _a in SA.MARKER_ALIASES]
SCOPE_TAGS = ["all", "essential_oil"]

marker_rows = []

grouped = defaultdict(dict)
for (scope, sp, study, lab), e in prof.items():
    grouped[(scope, sp, lab)][study] = e

for scope in SCOPE_TAGS:
    for sp in sorted({k[1] for k in grouped if k[0] == scope}):
        for lab in MARKER_ORDER:
            g = grouped.get((scope, sp, lab))
            if not g:
                continue
            per_study = sorted(g.items())
            for study, e in per_study:
                lo, hi, verdict = iso_verdict(sp, lab, e["median_pct"])
                marker_rows.append({
                    "material_scope": scope, "species": sp,
                    "species_status": status_of.get((scope, sp), "n/a"),
                    "level": "study", "study_id": study, "n_studies": 1,
                    "marker": lab, "inchikey": key_of.get(lab, ""),
                    "n_samples": e["n_samples"], "n_measured": e["n_measured"],
                    "n_quantified": e["n_quantified"],
                    "n_not_detected": e["n_not_detected"],
                    "n_trace_excluded": e["n_trace_excluded"],
                    "median_pct": fnum(e["median_pct"], 3),
                    "min_pct": fnum(e["min_pct"], 3),
                    "max_pct": fnum(e["max_pct"], 3),
                    "study_median_min": "", "study_median_max": "",
                    "study_median_fold_range": "",
                    "iso_lo": lo, "iso_hi": hi, "iso_verdict": verdict,
                    "cross_study_status": "",
                })
            meds = [e["median_pct"] for _s, e in per_study]
            n_st = len(per_study)
            pooled = stats.median(
                [v for _s, e in per_study for v in [e["median_pct"]] * 1])
            lo, hi, verdict = iso_verdict(sp, lab, pooled)
            if n_st < 2:
                cross = "untestable (single study)"
                fold = ""
            else:
                cross = "%d studies" % n_st
                mn, mx = min(meds), max(meds)
                fold = ("inf" if mn <= 0 < mx else
                        ("0" if mx <= 0 else fnum(mx / mn, 2)))
            marker_rows.append({
                "material_scope": scope, "species": sp,
                "species_status": status_of.get((scope, sp), "n/a"),
                "level": "species", "study_id": "(across studies)",
                "n_studies": n_st,
                "marker": lab, "inchikey": key_of.get(lab, ""),
                "n_samples": sum(e["n_samples"] for _s, e in per_study),
                "n_measured": sum(e["n_measured"] for _s, e in per_study),
                "n_quantified": sum(e["n_quantified"] for _s, e in per_study),
                "n_not_detected": sum(e["n_not_detected"] for _s, e in per_study),
                "n_trace_excluded": sum(e["n_trace_excluded"] for _s, e in per_study),
                "median_pct": fnum(pooled, 3),
                "min_pct": fnum(min(e["min_pct"] for _s, e in per_study), 3),
                "max_pct": fnum(max(e["max_pct"] for _s, e in per_study), 3),
                "study_median_min": fnum(min(meds), 3),
                "study_median_max": fnum(max(meds), 3),
                "study_median_fold_range": fold,
                "iso_lo": lo, "iso_hi": hi, "iso_verdict": verdict,
                "cross_study_status": cross,
            })

MARKER_COLS = ["material_scope", "species", "species_status", "level",
               "study_id", "n_studies", "marker", "inchikey", "n_samples",
               "n_measured", "n_quantified", "n_not_detected",
               "n_trace_excluded", "median_pct", "min_pct", "max_pct",
               "study_median_min", "study_median_max",
               "study_median_fold_range", "iso_lo", "iso_hi", "iso_verdict",
               "cross_study_status"]
with open(OUT_MARKERS, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=MARKER_COLS, delimiter="\t",
                       lineterminator="\n", extrasaction="ignore")
    w.writeheader()
    for r in marker_rows:
        w.writerow(r)
log("markers written -> %s (%d rows)" % (os.path.basename(OUT_MARKERS),
                                         len(marker_rows)))


TEST_COLS = ["material_scope", "block", "space", "mode", "factor", "species",
             "m_compounds", "n", "n_studies", "n_species", "n_bridging_species",
             "n_species_samples", "n_species_studies", "n_restricted",
             "n_cells_ge2", "df", "r2", "r2_species", "r2_residual",
             "r2_excess", "r2_adj", "F", "p", "cramers_v_species_study",
             "testable", "reason", "n_perm", "seed"]
with open(OUT_TESTS, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=TEST_COLS, delimiter="\t",
                       lineterminator="\n", extrasaction="ignore")
    w.writeheader()
    for t in TESTS:
        row = {}
        for c in TEST_COLS:
            v = t.get(c)
            row[c] = fnum(v, 4) if isinstance(v, float) else ("" if v is None else v)
        w.writerow(row)
log("tests written -> %s (%d tests)" % (os.path.basename(OUT_TESTS), len(TESTS)))


log("done")
