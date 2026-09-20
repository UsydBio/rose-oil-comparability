import csv
import math
import os


MISSING_LEVELS = frozenset([
    "", "-", "--", "na", "n/a", "nan", "none", "null", "unknown",
    "(unstated)", "unstated", "not", "?", ".",
])


ARTIFACT_LEVELS = frozenset([
    "sr. no", "sr no", "no.", "year", "ricalc", "rical", "ri", "lri",
    "peak area", "rt", "r.t", "mean", "sd", "stdev", "weight", "length",
    "caco-2 permeation", "binding affinity kcal mol", "ob%",
    "mean_ob", "mean_obl", "sd_ob", "sd_obl", "not",
])


ARTIFACT_PREFIXES = ("cluster ", "group ", "class ", "unnamed", "column",
                     "page ", "figure ", "table ")


MEASURE_HEAD_TOKENS = ("%", "conc", "oav", "rel.", "rel ", "relative", "content",
                       "area", "mg/", "μg", "ug/", "ng/", "peak", "amount",
                       "mean ", "average", "value", "ratio", "yield")


FACTORS = [
    ("cultivar",   "sample_cultivar",   "biological", 3.0),
    ("species",    "sample_species",    "biological", 3.0),
    ("stage",      "sample_stage",      "biological", 3.0),
    ("organ",      "sample_organ",      "biological", 2.5),
    ("origin",     "sample_origin",     "geographic", 2.0),
    ("processing", "sample_processing", "processing", 2.0),
    ("treatment",  "sample_treatment",  "treatment",  2.0),
    ("column_class", "column_class",    "design",     0.0),
]
FACTOR_COL = dict((n, c) for n, c, _k, _w in FACTORS)
FACTOR_KIND = dict((n, k) for n, _c, k, _w in FACTORS)
FACTOR_WEIGHT = dict((n, w) for n, _c, _k, w in FACTORS)


CLASS_TO_FACTOR = {
    "cultivar": "cultivar",
    "species_level": "species",
    "developmental_stage": "stage",
    "organ": "organ",
    "origin": "origin",
    "processing": "processing",
    "treatment": "treatment",
    "sampling_mode": None,
    "genuinely_unknown": None,
    "not_a_sample": None,
}


COMPARABLE_KINDS = frozenset(["relative_pct", "absolute_conc"])


MEASURED_STATUS = frozenset(["quantified", "not_detected"])


def read_tsv(path):

    csv.field_size_limit(10 ** 9)
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, fieldnames, rows):
    d = os.path.dirname(path)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def study_id_of(row):

    if row.get("pmcid"):
        return row["pmcid"]
    if row.get("doi"):
        return "doi:" + row["doi"]
    return "src:" + (row.get("table_id") or "").split("#")[0]


def sample_key(row):
    return (study_id_of(row), row.get("table_id") or "", row.get("column_index") or "")


def norm_doi(d):
    d = (d or "").strip().lower()
    for pre in ("https://doi.org/", "http://doi.org/", "doi:"):
        if d.startswith(pre):
            d = d[len(pre):]
    return d


def norm_level(v):

    v = (v or "").strip()
    low = v.lower()
    if not v or low in MISSING_LEVELS or low in ARTIFACT_LEVELS:
        return ""
    for p in ARTIFACT_PREFIXES:
        if low.startswith(p):
            return ""
    return v


def material_token(header):

    h = (header or "").strip()
    if "|" not in h:
        return " ".join(h.split()).lower()
    parts = [" ".join(p.split()) for p in h.split("|")]
    parts = [p for p in parts if p]
    if len(parts) > 1 and any(t in parts[0].lower() for t in MEASURE_HEAD_TOKENS):
        parts.pop(0)
    return "|".join(parts).lower()


def compound_key(row):

    sk = (row.get("inchikey_skeleton") or "").strip()
    if sk:
        return "sk:" + sk
    nm = (row.get("compound_clean_name") or row.get("compound_base_name")
          or row.get("compound_name_raw") or "").strip().lower()
    nm = " ".join(nm.split())
    return "nm:" + nm if nm else ""


class Sample(object):

    __slots__ = ("key", "study_id", "doi", "pmcid", "table_id", "column_index",
                 "header", "token", "caption", "title", "measure_kind", "unit", "tier",
                 "platform", "polarity", "column_class", "source_kind",
                 "levels", "measured", "quantified", "trace", "unmeasured", "n_rows")

    def __init__(self, key, row):
        self.key = key
        self.study_id = key[0]
        self.doi = norm_doi(row.get("doi"))
        self.pmcid = row.get("pmcid") or ""
        self.table_id = key[1]
        self.column_index = key[2]
        self.header = row.get("sample_column_header") or ""
        self.token = material_token(self.header) or (key[1] + "#" + key[2])
        self.caption = row.get("table_caption") or ""
        self.title = row.get("study_title") or ""
        self.measure_kind = row.get("measure_kind") or ""
        self.unit = row.get("unit") or ""
        self.tier = row.get("tier") or ""
        self.platform = row.get("study_platform_primary") or ""
        self.polarity = row.get("study_column_polarity") or ""
        self.column_class = row.get("column_class") or ""
        self.source_kind = row.get("source_kind") or ""
        self.levels = dict((n, norm_level(row.get(FACTOR_COL[n])))
                           for n, _c, _k, _w in FACTORS)
        self.measured = set()
        self.quantified = set()
        self.trace = set()
        self.unmeasured = set()
        self.n_rows = 0


def load_samples(table_path):

    samples = {}
    kept = dropped = 0
    csv.field_size_limit(10 ** 9)
    with open(table_path, encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("in_composition_scope") != "yes":
                dropped += 1
                continue
            kept += 1
            k = sample_key(row)
            s = samples.get(k)
            if s is None:
                s = samples[k] = Sample(k, row)
            s.n_rows += 1
            ck = compound_key(row)
            if not ck:
                continue
            st = row.get("detection_status") or ""
            if st in MEASURED_STATUS:
                s.measured.add(ck)
                if st == "quantified":
                    s.quantified.add(ck)
            elif st == "detected_not_quantified":
                s.trace.add(ck)
            else:
                s.unmeasured.add(ck)
    return samples, kept, dropped


def _dedupe_materials(samps):

    by_tok = {}
    for s in samps:
        by_tok.setdefault(s.token, []).append(s)
    out = {}
    for tok, group in by_tok.items():
        per = {}
        for s in group:
            per.setdefault((s.table_id, s.measure_kind), []).append(s)
        if all(len(v) == 1 for v in per.values()):
            out[tok] = group
        else:
            for s in group:
                out[(tok, s.table_id, s.column_index)] = [s]
    return out


def count_materials(samps):

    return len(_dedupe_materials(samps))


def _shared_strict(samps):

    it = iter(samps)
    acc = set(next(it).measured)
    for s in it:
        acc &= s.measured
        if not acc:
            break
    return acc


def _shared_any(level_samps):

    acc = None
    for samps in level_samps:
        u = set()
        for s in samps:
            u |= s.measured
        acc = u if acc is None else (acc & u)
        if not acc:
            break
    return acc or set()


def _covarying(level_of, samps, factor):

    out = []
    for other, _c, _k, _w in FACTORS:
        if other == factor:
            continue
        pairs = set()
        ok = True
        for s in samps:
            ov = s.levels[other]
            if not ov:
                ok = False
                break
            pairs.add((level_of[s.key], ov))
        if not ok:
            continue
        lv = set(p[0] for p in pairs)
        ov = set(p[1] for p in pairs)
        if len(ov) >= 2 and len(pairs) == len(lv) == len(ov):
            out.append(other)
    return out


def _contrast_record(study_id, ss_all, sub, factor, level_of, mk, level_source):

    by_level = {}
    for s in sub:
        by_level.setdefault(level_of[s.key], []).append(s)

    mats = dict((lv, _dedupe_materials(v)) for lv, v in by_level.items())
    n_mat = dict((lv, len(m)) for lv, m in mats.items())

    shared = _shared_strict(sub)
    shared_any = _shared_any(list(by_level.values()))
    q = None
    for s in sub:
        q = set(s.quantified) if q is None else (q & s.quantified)
    shared_q = (q or set()) & shared

    tabs_of_level = dict((lv, set(s.table_id for s in v)) for lv, v in by_level.items())
    all_tabs = set()
    for t in tabs_of_level.values():
        all_tabs |= t

    table_confounded = (len(all_tabs) == len(by_level)
                        and all(len(t) == 1 for t in tabs_of_level.values()))

    plats = set(s.platform for s in sub)
    pols = set(s.polarity for s in sub)
    classes = set(s.column_class for s in sub)
    tiers = set(s.tier for s in sub)

    axes = set(CLASS_TO_FACTOR.get(c) for c in classes if c)
    if FACTOR_WEIGHT.get(factor, 0) == 0:
        identity_status = "n/a"
    elif not any(classes):
        identity_status = "unclassified"
    elif factor in axes:
        identity_status = "confirmed"
    else:
        identity_status = "class_mismatch"

    suspect = sorted(set(lv for lv in by_level
                         if len(lv) >= 8
                         and any(lv.lower() in (s.caption or "").lower()
                                 for s in by_level[lv])))

    levels_sorted = sorted(by_level, key=lambda lv: (-n_mat[lv], lv))
    return {
        "study_id": study_id,
        "doi": sub[0].doi,
        "pmcid": sub[0].pmcid,
        "tier": ";".join(sorted(tiers)),
        "factor": factor,
        "factor_kind": FACTOR_KIND.get(factor, ""),
        "level_source": level_source,
        "measure_kind": mk,
        "n_levels": len(by_level),
        "n_columns": len(sub),
        "n_materials": sum(n_mat.values()),
        "min_materials_per_level": min(n_mat.values()),
        "max_materials_per_level": max(n_mat.values()),
        "n_levels_with_replication": sum(1 for v in n_mat.values() if v >= 2),
        "has_replication": min(n_mat.values()) >= 2,
        "levels": " | ".join("%s(n=%d)" % (lv[:38], n_mat[lv]) for lv in levels_sorted),
        "n_shared_compounds": len(shared),
        "n_shared_quantified": len(shared_q),
        "n_shared_any_sample": len(shared_any),
        "n_tables": len(all_tabs),
        "table_confounded": table_confounded,
        "platform": ";".join(sorted(p for p in plats if p)) or "(unstated)",
        "platform_constant": len(plats) == 1,
        "column_polarity": ";".join(sorted(p for p in pols if p)) or "(unstated)",
        "polarity_constant": len(pols) == 1,
        "column_class": ";".join(sorted(c for c in classes if c)) or "(unclassified)",
        "identity_status": identity_status,
        "values_comparable": mk in COMPARABLE_KINDS,
        "confounded_with": ";".join(_covarying(level_of, sub, factor)),
        "suspect_levels": ";".join(suspect),
        "_samples": [s.key for s in sub],
        "_levels": levels_sorted,
    }


def build_contrasts(samples, min_levels=2):

    by_study = {}
    for s in samples.values():
        by_study.setdefault(s.study_id, []).append(s)

    out = []
    for study_id in sorted(by_study):
        ss = by_study[study_id]
        kinds = sorted(set(s.measure_kind for s in ss))

        seen_field = set()
        for factor, _col, _kind_of, _w in FACTORS:
            for mk in kinds:
                sub = [s for s in ss if s.measure_kind == mk and s.levels[factor]]
                if len(sub) < 2:
                    continue
                level_of = dict((s.key, s.levels[factor]) for s in sub)
                if len(set(level_of.values())) < min_levels:
                    continue
                seen_field.add((factor, mk))
                out.append(_contrast_record(study_id, ss, sub, factor, level_of,
                                            mk, "factor_field"))

        for mk in kinds:
            by_axis = {}
            for s in ss:
                if s.measure_kind != mk or not s.column_class:
                    continue
                ax = CLASS_TO_FACTOR.get(s.column_class)
                if ax:
                    by_axis.setdefault(ax, []).append(s)
            for ax, sub in by_axis.items():
                if (ax, mk) in seen_field or len(sub) < 2:
                    continue
                level_of = dict((s.key, s.token) for s in sub)
                if len(set(level_of.values())) < min_levels:
                    continue
                out.append(_contrast_record(study_id, ss, sub, ax, level_of,
                                            mk, "column_code"))
    return out


def is_genuine(c, min_shared=3):

    return (c["factor"] != "column_class"
            and FACTOR_WEIGHT[c["factor"]] > 0
            and c["n_levels"] >= 2
            and c["min_materials_per_level"] >= 1
            and c["n_shared_compounds"] >= min_shared
            and c["platform_constant"]
            and c["polarity_constant"])


def score_contrast(c):

    s = FACTOR_WEIGHT[c["factor"]]
    s += 0.5 * min(c["n_levels"], 12)
    if c["has_replication"]:
        s += 3.0
    s += 2.0 * math.log10(1 + c["n_shared_quantified"])
    if not c["table_confounded"]:
        s += 1.5
    if c["values_comparable"]:
        s += 1.0
    if c["identity_status"] == "confirmed":
        s += 1.0
    if c["confounded_with"]:
        s -= 2.0
    return round(s, 2)
