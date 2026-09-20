import csv
import json
import os
from collections import defaultdict

from . import stats

__all__ = [
    "UNSTATED", "ISO_RANGES", "MARKER_ALIASES",
    "study_id_of", "normalize_species",
    "load_rows", "build_cells", "measured_sets", "load_material_types",
    "species_design", "block_stats", "freq_ladder", "bridge_ladder",
    "pick_blocks",
    "block_vectors", "aitchison_d2", "bray_d2",
    "resolve_marker_keys", "marker_profiles",
]

UNSTATED = "(unstated)"

DELTA_FRAC = 0.65
ZERO_MASS_CAP = 0.30


def study_id_of(r):

    if r.get("pmcid"):
        return r["pmcid"]
    if r.get("doi"):
        return "doi:" + r["doi"]
    return "src:" + (r.get("table_id") or "").split("#")[0]


def normalize_species(raw):

    sp = (raw or "").strip()
    if not sp:
        return ""
    for sep in (" (", ", "):
        if sep in sp:
            sp = sp.split(sep)[0]
    return " ".join(sp.split())


def load_rows(path, measure_kind="relative_pct"):

    csv.field_size_limit(10 ** 9)
    with open(path, encoding="utf-8") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t")
                if r.get("in_composition_scope") == "yes"
                and r.get("measure_kind") == measure_kind]


def build_cells(rows, value_field="area_pct_uncorrected"):

    raw = defaultdict(list)
    meta = {}
    names = {}
    units = defaultdict(set)
    dropped = 0
    for r in rows:
        s = (study_id_of(r), r["table_id"], r["column_index"])
        meta.setdefault(s, r)
        if r.get("unit"):
            units[s].add(r["unit"])
        sk = r.get("inchikey_skeleton") or ""
        if not sk:
            dropped += 1
            continue
        raw[(s, sk)].append(r)
        nm = (r.get("compound_clean_name") or r.get("compound_base_name")
              or r.get("compound_name_raw") or "")
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
    return cells, meta, names, units, multi, dropped


def measured_sets(cells):

    m = defaultdict(set)
    for (s, sk), (st, _v) in cells.items():
        if st in ("quantified", "not_detected"):
            m[s].add(sk)
    return m


def load_material_types(path):

    out = {}
    clash = 0
    if not os.path.exists(path):
        return out, 0
    csv.field_size_limit(10 ** 9)
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k = (r["table_id"], r["column_index"])
            v = r.get("material_type") or ""
            if k in out and out[k] != v:
                clash += 1
            out[k] = v
    return out, clash


def species_design(samples, sp_of):

    by = defaultdict(lambda: {"n": 0, "studies": set(), "samples": []})
    for s in samples:
        sp = sp_of.get(s) or UNSTATED
        by[sp]["n"] += 1
        by[sp]["studies"].add(s[0])
        by[sp]["samples"].append(s)
    out = []
    for sp, e in by.items():
        out.append({
            "species": sp,
            "n_samples": e["n"],
            "n_studies": len(e["studies"]),
            "studies": sorted(e["studies"]),
            "samples": e["samples"],
            "bridging": sp != UNSTATED and len(e["studies"]) >= 2,
            "status": ("bridging" if (sp != UNSTATED and len(e["studies"]) >= 2)
                       else ("unstated" if sp == UNSTATED
                             else "untestable (single study)")),
        })
    out.sort(key=lambda d: (-d["n_studies"], -d["n_samples"], d["species"]))
    return out


def block_stats(samples_in, sp_of):
    st = defaultdict(set)
    for s in samples_in:
        st[sp_of.get(s) or UNSTATED].add(s[0])
    named = {k: v for k, v in st.items() if k != UNSTATED}
    return {
        "n_samples": len(samples_in),
        "n_studies": len({s[0] for s in samples_in}),
        "n_species": len(named),
        "n_bridging_species": sum(1 for v in named.values() if len(v) >= 2),
        "species_studies": {k: len(v) for k, v in
                            sorted(named.items(), key=lambda kv: -len(kv[1]))},
    }


def freq_ladder(meas, samples, sp_of, max_m=30):

    per = defaultdict(set)
    for s in samples:
        for c in meas[s]:
            per[c].add(s)
    order = sorted(per, key=lambda c: (-len(per[c]), c))
    out = []
    for m in range(1, min(max_m, len(order)) + 1):
        pre = order[:m]
        ss = set(samples)
        for c in pre:
            ss &= per[c]
        rec = {"m": m, "compounds": list(pre), "samples": sorted(ss),
               "added": pre[-1]}
        rec.update(block_stats(ss, sp_of))
        out.append(rec)
    return out


def bridge_ladder(meas, samples, sp_of, max_m=12):

    per = defaultdict(set)
    for s in samples:
        for c in meas[s]:
            per[c].add(s)
    cur = []
    pool = set(per)
    out = []
    while pool and len(cur) < max_m:
        best = None
        for c in sorted(pool):
            ss = set(samples)
            for cc in cur + [c]:
                ss &= per[cc]
            if not ss:
                continue
            bs = block_stats(ss, sp_of)
            score = (bs["n_bridging_species"], (len(cur) + 1) * bs["n_samples"])
            if best is None or score > best[0]:
                best = (score, c, ss, bs)
        if best is None:
            break
        cur.append(best[1])
        pool.discard(best[1])
        rec = {"m": len(cur), "compounds": list(cur),
               "samples": sorted(best[2]), "added": best[1]}
        rec.update(best[3])
        out.append(rec)
    return out


MIN_M = 2
MIN_M_PRIMARY = 3


def pick_blocks(lad):

    picks = {}
    c1 = [b for b in lad if b["m"] >= MIN_M and b["n_species"] >= 2]
    if c1:
        picks["SB1"] = max(c1, key=lambda b: (b["n_samples"], b["m"]))
    c2 = [b for b in lad if b["m"] >= MIN_M_PRIMARY
          and b["n_bridging_species"] >= 3]
    if c2:
        picks["SB2"] = max(c2, key=lambda b: (b["m"] * b["n_samples"], b["m"]))
    c3 = [b for b in lad if b["m"] >= MIN_M and b["n_species"] >= 2
          and b["n_studies"] >= 2]
    if c3:
        picks["SB3"] = max(c3, key=lambda b: (b["m"], b["n_samples"]))
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
    closed, clrs = [], []
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
    return keep, raw, closed, clrs, n_capped


def aitchison_d2(clr_rows):

    return stats.sq_dist_matrix(clr_rows, stats.euclidean)


def bray_d2(closed_rows):

    return stats.sq_dist_matrix(closed_rows, stats.bray_curtis)


ISO_RANGES = [
    ("ethanol", ["ethanol"], None, 3.0),
    ("phenylethyl alcohol", ["phenylethyl alcohol", "2-phenylethanol",
                             "phenethyl alcohol"], None, 2.5),
    ("citronellol", ["citronellol"], 20.0, 34.0),
    ("nerol", ["nerol"], 5.0, 12.0),
    ("geraniol", ["geraniol"], 14.0, 22.0),
    ("methyl eugenol", ["methyl eugenol", "methyleugenol"], 0.8, 3.0),
    ("heptadecane", ["heptadecane"], 1.0, 2.5),
    ("nonadecane", ["nonadecane"], 8.0, 15.0),
    ("heneicosane", ["heneicosane"], 3.0, 5.5),
]


MARKER_ALIASES = [
    ("citronellol", ["citronellol"]),
    ("geraniol", ["geraniol"]),
    ("nerol", ["nerol"]),
    ("phenylethyl alcohol", ["phenylethyl alcohol", "2-phenylethanol",
                             "phenethyl alcohol"]),
    ("linalool", ["linalool"]),
    ("methyl eugenol", ["methyl eugenol", "methyleugenol"]),
    ("eugenol", ["eugenol"]),
    ("nonadecane", ["nonadecane"]),
    ("heneicosane", ["heneicosane"]),
    ("heptadecane", ["heptadecane"]),
    ("ethanol", ["ethanol"]),
]


def resolve_marker_keys(cache_path, specs):

    cache = {}
    if os.path.exists(cache_path):
        with open(cache_path, encoding="utf-8") as fh:
            cache = json.load(fh)
    name_to_key = {}
    for cname, info in cache.items():
        if isinstance(info, dict) and info.get("inchikey"):
            name_to_key.setdefault(cname.strip().lower(), info["inchikey"])
    out = []
    for label, aliases in specs:
        key = ""
        for a in aliases:
            k = name_to_key.get(a.lower())
            if k:
                key = k
                break
        out.append((label, key))
    return out


def marker_profiles(rows, sp_of, marker_keys, material_of=None,
                    value_field="area_pct_uncorrected"):

    want = {k: lab for lab, k in marker_keys if k}
    acc = defaultdict(lambda: {"vals": [], "n_q": 0, "n_nd": 0,
                               "samples": set(), "n_trace": 0})
    for r in rows:
        key = r.get("inchikey") or ""
        lab = want.get(key)
        if lab is None:
            continue
        s = (study_id_of(r), r["table_id"], r["column_index"])
        sp = sp_of.get(s) or UNSTATED
        mat = (material_of or {}).get((s[1], s[2])) or "(not classified)"
        st = r["detection_status"]
        v = None
        if st == "quantified" and r[value_field]:
            try:
                v = float(r[value_field])
            except ValueError:
                v = None
        for scope in ("all", mat):
            e = acc[(scope, sp, s[0], lab)]
            if v is not None:
                e["vals"].append(v)
                e["n_q"] += 1
                e["samples"].add(s)
            elif st == "not_detected":
                e["vals"].append(0.0)
                e["n_nd"] += 1
                e["samples"].add(s)
            elif st == "detected_not_quantified":
                e["n_trace"] += 1
    out = {}
    for k, e in acc.items():
        if not e["vals"]:
            continue
        out[k] = {
            "n_measured": len(e["vals"]),
            "n_samples": len(e["samples"]),
            "n_quantified": e["n_q"],
            "n_not_detected": e["n_nd"],
            "n_trace_excluded": e["n_trace"],
            "median_pct": stats.median(e["vals"]),
            "min_pct": min(e["vals"]),
            "max_pct": max(e["vals"]),
        }
    return out

