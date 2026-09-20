import csv
import math
import os
import random

from . import contrasts as K
from . import stats as S

__all__ = [
    "SEED", "N_PERM", "DELTA_FRAC", "ZERO_MASS_CAP",
    "load_value_cells", "ContrastSpec", "build_specs",
    "analyse_contrast", "benjamini_hochberg", "direction_agreement",
]


SEED = 20260912
N_PERM = 9999
DELTA_FRAC = 0.65
ZERO_MASS_CAP = 0.30
MAX_EXACT_ORDERS = 5040


COMPARABLE_KINDS = K.COMPARABLE_KINDS


CAVEATS = {
    "PMC11443146|species|absolute_conc":
        "species is fully confounded with organ / stage / time of day: the "
        "R. gigantea side is 5 flower-opening stages + 6 organs + 6 diel time "
        "points, the R. chinensis side is 5 flower-opening stages only. The "
        "test separates two sets of samples, not two species.",
    "PMC10457957|processing|relative_pct":
        "the 3 'replicates' of the hexane level are three different solvent "
        "systems (Hexane, Hexane:Dichloromethane 1:5, Hexane:Methanol 1:5) that "
        "the processing field lumps into one string. Not replication.",
    "PMC8398089|species|relative_pct":
        "the within-species 'replicates' are distinct cultivar/colour/part codes "
        "(RPOA / RPWA / RPOF / RPWF vs RBWA / RBWF), i.e. a nested design read as "
        "a replicated one. Within-group SS therefore absorbs cultivar variation, "
        "so the test is conservative, not invalid.",
    "PMC11533618|origin|oav":
        "the Panzhou level has 2 materials only because the header `OAV | PZ` "
        "occurs twice in the same table; one of those two columns carries "
        "concentrations, not odour-activity values. Not replication.",
    "PMC11554761|treatment|unspecified|collapsed":
        "the 186-compound basis is a published list of differentially accumulated "
        "metabolites, already filtered for significance between exactly these "
        "two groups. R2 on this basis is a selection-inflated upper "
        "bound, not an effect size; the p-value is likewise circular.",
    "PMC11554761|treatment|unspecified":
        "the 6 'levels' are 2 treatments x 3 biological replicates: the factor "
        "field wrote the replicate number into the level name. Use the derived "
        "collapsed contrast instead.",
    "PMC11443146|stage|absolute_conc":
        "one 16-level field holding flower-opening stages of two species plus a "
        "diel series. Superseded by three derived contrasts.",
}


PRESELECTED_TABLE_TOKENS = ("#dams", "#deg", "#significant", "#diff")


AGGREGATE_NAME_TOKENS = (
    "monoterpene", "sesquiterpene", "hydrocarbon", "aliphatic", "oxygenated",
    "total", "ratio", "component number", "aromatic alcohol", "aromatic aldehyde",
    "sum of", "others",
)


def looks_aggregate(name):
    low = (name or "").lower()
    return any(t in low for t in AGGREGATE_NAME_TOKENS)


def refile_axis(study_id, factor, levels):

    low = " ".join(levels).lower()
    if factor == "stage":
        if study_id == "PMC11443146" and "diel" in low:
            if "bud" in low or "open" in low:
                return ("mixed_stage_diel",
                        "one 16-level stage field holding TWO different things: "
                        "flower-opening stages S1-S5 of two species AND a diel "
                        "time-of-day series T00-T20, not one axis; superseded by "
                        "the three derived contrasts")
            return "diel", "diel time-of-day series filed under developmental_stage"
        if study_id == "PMC11431035":
            return "stage_fruit", "fruit ripening (days after anthesis), not flower opening"
        return "stage_floral", ""
    if factor == "species":
        if study_id == "PMC11242971":
            return "colour_morph", "white vs dark-pink R. canina: flower-colour morph filed under species_level"
        return "species", ""
    if factor == "processing" and study_id == "PMC13409530":
        return "sampling_mode", "HS vs SPME is an analytical sampling mode, not a biological/processing factor"
    return factor, ""


ANALYTICAL_AXES = frozenset(["sampling_mode"])


def load_value_cells(table_path):

    csv.field_size_limit(10 ** 9)
    raw = {}
    names = {}
    n_rows = 0
    with open(table_path, encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("in_composition_scope") != "yes":
                continue
            ck = K.compound_key(row)
            if not ck:
                continue
            n_rows += 1
            sk = K.sample_key(row)
            raw.setdefault((sk, ck), []).append(row)
            if ck not in names:
                nm = (row.get("compound_clean_name") or row.get("compound_base_name")
                      or row.get("compound_name_raw") or "")
                if nm:
                    names[ck] = nm
    cells = {}
    for key, rs in raw.items():
        vals = []
        st = set()
        for r in rs:
            st.add(r.get("detection_status") or "")
            if r.get("detection_status") == "quantified" and r.get("value"):
                try:
                    vals.append(float(r["value"]))
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
    return cells, names, n_rows


class ContrastSpec(object):

    __slots__ = ("cid", "study_id", "doi", "pmcid", "tier", "title",
                 "factor", "axis", "axis_note", "level_source", "measure_kind",
                 "kind", "parent", "derivation", "sample_keys", "level_of",
                 "enum", "order_of")

    def __init__(self, cid, study_id, factor, axis, axis_note, level_source,
                 measure_kind, kind, parent, derivation, sample_keys, level_of,
                 enum, samples):
        s0 = samples[sample_keys[0]]
        self.cid = cid
        self.study_id = study_id
        self.doi = s0.doi
        self.pmcid = s0.pmcid
        self.tier = ";".join(sorted(set(samples[k].tier for k in sample_keys)))
        self.title = s0.title
        self.factor = factor
        self.axis = axis
        self.axis_note = axis_note
        self.level_source = level_source
        self.measure_kind = measure_kind
        self.kind = kind
        self.parent = parent
        self.derivation = derivation
        self.sample_keys = sample_keys
        self.level_of = level_of
        self.enum = enum or {}
        self.order_of = order_values(axis, sorted(set(level_of.values())))

    @property
    def levels(self):
        return sorted(set(self.level_of.values()))

    @property
    def ordered(self):
        return self.order_of is not None


def level_of_sample(sample, factor, level_source):

    if level_source == "column_code":
        return sample.token
    return K.norm_level(sample.levels[factor])


REPLICATE_MARKERS = (", biological replicate ", ", replicate ", ", rep ")


def strip_replicate(level):

    low = level.lower()
    for m in REPLICATE_MARKERS:
        i = low.find(m)
        if i >= 0:
            return level[:i].strip().rstrip(",")
    return level


_FLORAL_WORDS = [
    ("young bud", 1.0), ("bud stage", 2.0), ("initial-open", 3.0),
    ("initial open", 3.0), ("half-open", 4.0), ("half open", 4.0),
    ("full-open", 5.0), ("full open", 5.0), ("full bloom", 5.0),
    ("senescen", 6.0),
]


def _num_prefix(text):

    digits = ""
    for ch in text:
        if ch.isdigit():
            digits += ch
        elif digits:
            break
    return float(digits) if digits else None


def order_values(axis, levels):

    if axis not in ("stage_floral", "stage_fruit", "diel"):
        return None
    out = {}
    for lv in levels:
        low = lv.lower()
        v = None
        if axis == "diel":
            i = low.find("(t")
            if i >= 0:
                v = _num_prefix(low[i + 2:])
            if v is None:
                i = low.find("time point")
                if i >= 0:
                    v = _num_prefix(low[i:])
        elif axis == "stage_fruit":
            if "day" in low:
                v = _num_prefix(low)
        else:
            for w, r in _FLORAL_WORDS:
                if w in low:
                    v = r
                    break
            if v is None:

                t = low.strip()
                if len(t) == 2 and t[0] == "s" and t[1].isdigit():
                    v = float(t[1])
                else:
                    i = low.find("(s")
                    if i >= 0:
                        v = _num_prefix(low[i + 2:])
        if v is None:
            return None
        out[lv] = v
    if len(set(out.values())) != len(out):
        return None
    return out


def build_specs(samples, prebuilt_contrasts=None):

    recs = prebuilt_contrasts if prebuilt_contrasts is not None else K.build_contrasts(samples)
    genuine = [c for c in recs if K.is_genuine(c)]
    specs = []
    for c in genuine:
        keys = list(c["_samples"])
        level_of = {}
        for k in keys:
            level_of[k] = level_of_sample(samples[k], c["factor"], c["level_source"])
        axis, note = refile_axis(c["study_id"], c["factor"], sorted(set(level_of.values())))
        cid = "%s|%s|%s" % (c["study_id"], c["factor"], c["measure_kind"])
        specs.append(ContrastSpec(cid, c["study_id"], c["factor"], axis, note,
                                  c["level_source"], c["measure_kind"],
                                  "enumerated", "", "", keys, level_of, c, samples))
    specs.extend(_derived_specs(specs, samples))
    return specs, genuine


def _derived_specs(specs, samples):

    out = []
    by_id = dict((s.cid, s) for s in specs)

    p = by_id.get("PMC11554761|treatment|unspecified")
    if p is not None:
        lo = dict((k, strip_replicate(v)) for k, v in p.level_of.items())
        if len(set(lo.values())) == 2:
            out.append(ContrastSpec(
                p.cid + "|collapsed", p.study_id, p.factor, p.axis, p.axis_note,
                p.level_source, p.measure_kind, "derived", p.cid,
                "collapsed ', biological replicate N of 3' into 2 treatment groups "
                "(n=3 each), the one design here with true replication",
                list(p.sample_keys), lo, p.enum, samples))

    p = by_id.get("PMC11443146|stage|absolute_conc")
    if p is not None:
        subsets = [
            ("gigantea_floral", "stage_floral",
             lambda lv: "gigantea" in lv.lower() and "diel" not in lv.lower(),
             "R. gigantea flower-opening series S1-S5, split out of the 16-level "
             "stage field (which also holds a diel series and a second species)"),
            ("chinensis_floral", "stage_floral",
             lambda lv: "chinensis" in lv.lower() and "diel" not in lv.lower(),
             "R. chinensis 'Old Blush' flower-opening series S1-S5, same split"),
            ("gigantea_diel", "diel",
             lambda lv: "diel" in lv.lower(),
             "R. gigantea diel series T00-T20 is time of day, not a developmental "
             "stage; must not be pooled with S1-S5"),
        ]
        for tag, axis, keep, why in subsets:
            lo = dict((k, v) for k, v in p.level_of.items() if keep(v))
            if len(set(lo.values())) < 2:
                continue
            keys = [k for k in p.sample_keys if k in lo]
            out.append(ContrastSpec(
                p.cid + "|" + tag, p.study_id, p.factor, axis,
                "split out of a 16-level stage field", p.level_source,
                p.measure_kind, "derived", p.cid, why, keys, lo, p.enum, samples))
    return out


def materials_of(spec, samples):

    subs = [samples[k] for k in spec.sample_keys]
    rank = {}
    for s in subs:
        r = rank.setdefault(s.table_id, [0, 0])
        r[0] += 1
        r[1] += len(s.measured)

    mats = {}
    by_level = {}
    for s in subs:
        by_level.setdefault(spec.level_of[s.key], []).append(s)
    for lv in sorted(by_level):
        for mid, group in K._dedupe_materials(by_level[lv]).items():
            label = mid if isinstance(mid, str) else "|".join(str(x) for x in mid)
            mats[(lv, label)] = group
    out = []
    for (lv, label), group in mats.items():
        group = sorted(group, key=lambda s: (-rank[s.table_id][0], -rank[s.table_id][1],
                                             s.table_id, int(s.column_index or 0)))
        out.append((label, group[0], group))
    out.sort(key=lambda t: (spec.level_of.get(t[1].key, ""), t[1].table_id,
                            int(t[1].column_index or 0)))
    return out


def shared_basis(cols):

    it = iter(cols)
    acc = set(next(it).measured)
    for s in it:
        acc &= s.measured
        if not acc:
            break
    return acc


def value_matrix(cols, basis, cells):

    rows = []
    for s in cols:
        v = []
        for c in basis:
            cell = cells.get((s.key, c))
            v.append(0.0 if cell is None or cell[1] is None else cell[1])
        rows.append(v)
    return rows


def drop_cross_table_pseudoreplicates(mats, level_of):

    by_level = {}
    for m in mats:
        by_level.setdefault(level_of[m[1].key], []).append(m)
    keep, dropped = [], []
    n_cross = 0
    for lv in sorted(by_level):
        group = by_level[lv]
        per = {}
        for m in group:
            per.setdefault(m[1].table_id, []).append(m)
        if len(per) > 1 and len(group) > 1:
            n_cross += 1
        best = sorted(per, key=lambda t: (-len(per[t]), t))[0]
        for t in sorted(per):
            (keep if t == best else dropped).extend(per[t])
    keep.sort(key=lambda m: (level_of[m[1].key], m[1].table_id,
                             int(m[1].column_index or 0)))
    return keep, dropped, n_cross


def drop_identical(rows, labels, ids):

    seen = {}
    keep, dropped = [], []
    for i, row in enumerate(rows):
        sig = (labels[i], tuple(round(x, 12) for x in row))
        if sig in seen:
            dropped.append((ids[i], ids[seen[sig]]))
            continue
        seen[sig] = i
        keep.append(i)
    return keep, dropped


def clr_block(raw, delta_frac=DELTA_FRAC, cap=ZERO_MASS_CAP):

    keep = [i for i, row in enumerate(raw) if sum(row) > 0]
    sub = [raw[i] for i in keep]
    if not sub:
        return [], [], [], 0
    ncol = len(sub[0])
    allpos = [x for row in sub for x in row if x > 0]
    fallback = min(allpos) if allpos else 1e-6
    deltas = []
    for j in range(ncol):
        pos = [row[j] for row in sub if row[j] > 0]
        deltas.append(delta_frac * (min(pos) if pos else fallback))
    closed, clrs = [], []
    n_capped = 0
    for row in sub:
        cl = S.closure(row, 100.0)
        dd = [d * 100.0 / sum(row) for d in deltas]
        zm = sum(d for v, d in zip(cl, dd) if v <= 0)
        if zm > cap * 100.0:
            f = cap * 100.0 / zm
            dd = [d * f for d in dd]
            n_capped += 1
        closed.append(cl)
        clrs.append(S.clr(S.multiplicative_replacement(cl, dd)))
    return keep, closed, clrs, n_capped


def _col(rows, j):
    return [r[j] for r in rows]


def _level_range(values, labels):

    g = {}
    for v, lab in zip(values, labels):
        g.setdefault(lab, []).append(v)
    ms = [sum(x) / len(x) for x in g.values()]
    return max(ms) - min(ms)


def _eta2(values, labels):

    n = len(values)
    gm = sum(values) / n
    sst = sum((v - gm) ** 2 for v in values)
    if sst <= 0:
        return 0.0
    groups = {}
    for v, lab in zip(values, labels):
        groups.setdefault(lab, []).append(v)
    ssa = 0.0
    for g in groups.values():
        m = sum(g) / len(g)
        ssa += len(g) * (m - gm) ** 2
    return ssa / sst


def per_compound_effects(clr_rows, closed, labels, basis, names, n_perm, seed,
                         testable):

    if not clr_rows:
        return []
    levels = sorted(set(labels))
    two = len(levels) == 2
    eff, raw_eff = [], []
    cols = [_col(clr_rows, j) for j in range(len(basis))]
    rcols = [_col(closed, j) for j in range(len(basis))]
    for j in range(len(basis)):
        col = cols[j]
        rcol = rcols[j]
        if two:
            a = [v for v, l in zip(col, labels) if l == levels[0]]
            b = [v for v, l in zip(col, labels) if l == levels[1]]
            eff.append(S.mean(b) - S.mean(a))
            ra = [v for v, l in zip(rcol, labels) if l == levels[0]]
            rb = [v for v, l in zip(rcol, labels) if l == levels[1]]
            raw_eff.append(S.mean(rb) - S.mean(ra))
        elif testable:
            eff.append(_eta2(col, labels))
            raw_eff.append(_eta2(rcol, labels))
        else:

            eff.append(_level_range(col, labels))
            raw_eff.append(_level_range(rcol, labels))
    pvals = [None] * len(basis)
    if testable and n_perm > 0:
        rng = random.Random(seed)
        perm = list(labels)
        ge = [0] * len(basis)
        obs = [abs(e) for e in eff]
        nb = len(basis)
        for _ in range(n_perm):
            rng.shuffle(perm)
            if two:
                ia = [i for i, l in enumerate(perm) if l == levels[0]]
                ib = [i for i, l in enumerate(perm) if l == levels[1]]
                na, nbb = float(len(ia)), float(len(ib))
                for j in range(nb):
                    col = cols[j]
                    st = abs(sum(col[i] for i in ib) / nbb
                             - sum(col[i] for i in ia) / na)
                    if st >= obs[j]:
                        ge[j] += 1
            else:
                for j in range(nb):
                    if _eta2(cols[j], perm) >= obs[j]:
                        ge[j] += 1
        pvals = [(g + 1) / (n_perm + 1.0) for g in ge]
        qs = benjamini_hochberg(pvals)
    else:
        qs = [None] * len(basis)
    out = []
    for j, c in enumerate(basis):
        out.append({
            "compound_key": c,
            "compound_name": names.get(c, c),
            "effect_clr": eff[j],
            "effect_raw": raw_eff[j],
            "p": pvals[j],
            "q": qs[j],
            "statistic": ("clr_mean_diff" if two else
                          ("clr_eta2" if testable else "clr_level_range")),
        })
    out.sort(key=lambda d: -abs(d["effect_clr"]))
    return out


def _orders(k, n_perm, seed):

    fact = 1
    for i in range(2, k + 1):
        fact *= i
    if fact <= MAX_EXACT_ORDERS:
        import itertools
        return [list(p) for p in itertools.permutations(range(k))], True
    rng = random.Random(seed)
    base = list(range(k))
    out = []
    for _ in range(n_perm):
        rng.shuffle(base)
        out.append(list(base))
    return out, False


def _pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return num / (dx * dy) if dx > 0 and dy > 0 else 0.0


def _ranks(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for t in range(i, j + 1):
            r[order[t]] = avg
        i = j + 1
    return r


def spearman(xs, ys):
    return _pearson(_ranks(xs), _ranks(ys))


def ordered_tests(clr_rows, labels, level_order, basis, names, n_perm, seed):

    levels = sorted(set(labels), key=lambda l: level_order[l])
    if len(levels) < 3:
        return None
    idx_of = dict((l, i) for i, l in enumerate(levels))
    ovals = [level_order[l] for l in levels]
    mat_order = [ovals[idx_of[l]] for l in labels]

    d2 = S.sq_dist_matrix(clr_rows, S.euclidean)
    n = len(clr_rows)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    dvec = [math.sqrt(d2[i][j]) for i, j in pairs]

    def mantel_stat(mo):
        ovec = [abs(mo[i] - mo[j]) for i, j in pairs]
        return _pearson(dvec, ovec)

    r_obs = mantel_stat(mat_order)
    cols = [_col(clr_rows, j) for j in range(len(basis))]
    rho_obs = [spearman(mat_order, c) for c in cols]

    perms, exact = _orders(len(levels), n_perm, seed)
    ge_m = 0
    ge_c = [0] * len(basis)
    for p in perms:
        mo = [ovals[p[idx_of[l]]] for l in labels]
        if mantel_stat(mo) >= r_obs:
            ge_m += 1
        for j, c in enumerate(cols):
            if abs(spearman(mo, c)) >= abs(rho_obs[j]):
                ge_c[j] += 1
    npm = len(perms)
    p_mantel = ge_m / float(npm) if exact else (ge_m + 1) / (npm + 1.0)
    pc = [(g / float(npm) if exact else (g + 1) / (npm + 1.0)) for g in ge_c]
    qc = benjamini_hochberg(pc)
    comp = []
    for j, c in enumerate(basis):
        comp.append({"compound_key": c, "compound_name": names.get(c, c),
                     "rho": rho_obs[j], "p": pc[j], "q": qc[j]})
    comp.sort(key=lambda d: (d["p"], -abs(d["rho"])))

    q_floor = min(1.0, len(basis) / float(npm))
    return {"mantel_r": r_obs, "mantel_p": p_mantel, "n_orders": npm,
            "exact": exact, "levels_in_order": levels, "compounds": comp,
            "trend_q_floor": q_floor,
            "n_trend_q05": sum(1 for d in comp if d["q"] is not None and d["q"] <= 0.05)}


def min_achievable_p(sizes):

    n = sum(sizes)
    if n == 0:
        return None
    total = math.factorial(n)
    for s in sizes:
        total //= math.factorial(s)
    mult = {}
    for s in sizes:
        mult[s] = mult.get(s, 0) + 1
    equiv = 1
    for m in mult.values():
        equiv *= math.factorial(m)
    return equiv / float(total)


def benjamini_hochberg(pvals):

    idx = [i for i, p in enumerate(pvals) if p is not None]
    if not idx:
        return [None] * len(pvals)
    m = len(idx)
    order = sorted(idx, key=lambda i: pvals[i])
    q = [None] * len(pvals)
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        val = min(prev, pvals[i] * m / rank)
        q[i] = val
        prev = val
    return q


def analyse_contrast(spec, samples, cells, names, n_perm=N_PERM, seed=SEED):

    all_mats = materials_of(spec, samples)
    mats, dropped_ct, n_cross = drop_cross_table_pseudoreplicates(
        all_mats, spec.level_of)
    ids = [m[0] if isinstance(m[0], str) else "|".join(str(x) for x in m[0])
           for m in mats]
    cols = [m[1] for m in mats]
    labels = [spec.level_of[c.key] for c in cols]

    res = {
        "contrast_id": spec.cid, "kind": spec.kind, "parent": spec.parent,
        "derivation": spec.derivation,
        "study_id": spec.study_id, "doi": spec.doi, "pmcid": spec.pmcid,
        "tier": spec.tier, "study_title": spec.title,
        "axis": spec.axis, "axis_note": spec.axis_note,
        "caveat": CAVEATS.get(spec.cid, ""),
        "factor_field": spec.factor, "level_source": spec.level_source,
        "measure_kind": spec.measure_kind,
        "values_comparable": spec.measure_kind in COMPARABLE_KINDS,
        "n_columns": len(spec.sample_keys),
        "n_materials_enumerated": len(all_mats),
        "n_materials": len(mats),
        "n_cross_table_pseudoreps_dropped": len(dropped_ct),
        "n_levels_with_cross_table_pseudoreps": n_cross,
        "cross_table_pseudoreps": ";".join(
            (m[0] if isinstance(m[0], str) else "|".join(str(x) for x in m[0]))
            for m in dropped_ct[:8]),
        "n_levels": len(set(labels)),
        "n_shared_enumerated": spec.enum.get("n_shared_compounds", ""),
        "table_confounded": spec.enum.get("table_confounded", ""),
        "confounded_with": spec.enum.get("confounded_with", ""),
        "identity_status": spec.enum.get("identity_status", ""),
        "platform": spec.enum.get("platform", ""),
        "column_polarity": spec.enum.get("column_polarity", ""),
    }

    basis = sorted(shared_basis(cols))
    res["n_basis"] = len(basis)
    res["n_basis_inchikey"] = sum(1 for c in basis if c.startswith("sk:"))
    tabs = set(c.table_id for c in cols)
    res["n_tables"] = len(tabs)
    res["n_basis_aggregate_rows"] = sum(
        1 for c in basis if looks_aggregate(names.get(c, "")))
    res["basis_preselected"] = any(
        any(t in tb.lower() for t in PRESELECTED_TABLE_TOKENS) for tb in tabs)

    if dropped_ct:
        acols = [m[1] for m in all_mats]
        ab = sorted(shared_basis(acols))
        if len(ab) >= 3:
            araw = value_matrix(acols, ab, cells)
            ak, _acl, aclr, _ = clr_block(araw)
            alab = [spec.level_of[acols[i].key] for i in ak]
            atab = [acols[i].table_id for i in ak]
            same, diff = [], []
            for i in range(len(aclr)):
                for j in range(i + 1, len(aclr)):
                    d = S.euclidean(aclr[i], aclr[j])
                    if alab[i] == alab[j] and atab[i] != atab[j]:
                        same.append(d)
                    elif alab[i] != alab[j]:
                        diff.append(d)
            res["pseudorep_pair_aitchison"] = S.mean(same) if same else None
            res["between_level_aitchison_predrop"] = S.mean(diff) if diff else None

    if len(basis) < 3:
        res["status"] = "basis_below_3"
        res["inference_class"] = "descriptive"
        return res, [], None, None

    raw = value_matrix(cols, basis, cells)
    keep, dropped = drop_identical(raw, labels, ids)
    res["n_duplicate_materials_dropped"] = len(dropped)
    res["duplicate_pairs"] = ";".join("%s=%s" % (a, b) for a, b in dropped[:6])
    raw = [raw[i] for i in keep]
    labels = [labels[i] for i in keep]
    ids = [ids[i] for i in keep]

    kz, closed, clrs, n_capped = clr_block(raw)
    res["n_allzero_materials_dropped"] = len(raw) - len(kz)
    labels = [labels[i] for i in kz]
    ids = [ids[i] for i in kz]
    res["n_zero_mass_capped"] = n_capped
    res["n_analysed"] = len(clrs)
    res["n_levels_analysed"] = len(set(labels))

    sizes = {}
    for l in labels:
        sizes[l] = sizes.get(l, 0) + 1
    res["level_sizes"] = ";".join("%s=%d" % (l, sizes[l]) for l in sorted(sizes))
    res["min_n_per_level"] = min(sizes.values()) if sizes else 0
    res["max_n_per_level"] = max(sizes.values()) if sizes else 0
    res["n_levels_replicated"] = sum(1 for v in sizes.values() if v >= 2)
    res["p_min_achievable"] = min_achievable_p(sorted(sizes.values()))

    if res["n_levels_analysed"] < 2:
        res["status"] = "collapsed_below_2_levels"
        res["inference_class"] = "descriptive"
        return res, [], None, None

    if res["min_n_per_level"] >= 2:
        res["inference_class"] = "inferential"
    elif res["max_n_per_level"] >= 2:
        res["inference_class"] = "weak_inferential"
    else:
        res["inference_class"] = "descriptive"
    testable = res["max_n_per_level"] >= 2 and len(clrs) - res["n_levels_analysed"] > 0

    d2a = S.sq_dist_matrix(clrs, S.euclidean)
    d2b = S.sq_dist_matrix(closed, S.bray_curtis)
    if testable:
        pa = S.permanova(d2a, labels, n_perm=n_perm, seed=seed)
        pb = S.permanova(d2b, labels, n_perm=n_perm, seed=seed)
        res["status"] = "tested"
        for tag, p in (("aitchison", pa), ("bray", pb)):
            res["r2_" + tag] = p["r2"]
            res["r2_excess_" + tag] = p["r2_excess"]
            res["r2_adj_" + tag] = p["r2_adj"]
            res["F_" + tag] = p["F"]
            res["p_" + tag] = p["p"]
        res["df_among"] = pa["df_among"]
        res["df_within"] = pa["df_within"]
        res["n_perm"] = n_perm
    else:
        res["status"] = "no_within_level_replication"
        for tag in ("aitchison", "bray"):
            for f in ("r2_", "r2_excess_", "r2_adj_", "F_", "p_"):
                res[f + tag] = None
        res["df_among"] = res["n_levels_analysed"] - 1
        res["df_within"] = len(clrs) - res["n_levels_analysed"]
        res["n_perm"] = 0

    n = len(clrs)
    bet, wit = [], []
    for i in range(n):
        for j in range(i + 1, n):
            (bet if labels[i] != labels[j] else wit).append(math.sqrt(d2a[i][j]))
    res["mean_between_level_aitchison"] = S.mean(bet) if bet else None
    res["mean_within_level_aitchison"] = S.mean(wit) if wit else None
    res["separation_ratio"] = (S.mean(bet) / S.mean(wit)) if wit and S.mean(wit) > 0 else None

    comp = per_compound_effects(clrs, closed, labels, basis, names,
                                n_perm if testable else 0, seed, testable)
    res["top_drivers"] = "; ".join(
        "%s(%+.2f)" % (d["compound_name"], d["effect_clr"]) for d in comp[:6])
    res["n_drivers_q05"] = sum(1 for d in comp if d["q"] is not None and d["q"] <= 0.05)

    trend = None
    if spec.ordered and res["n_levels_analysed"] >= 3:
        lo = dict((l, spec.order_of[l]) for l in set(labels))
        trend = ordered_tests(clrs, labels, lo, basis, names, n_perm, seed)
        if trend:
            res["mantel_r"] = trend["mantel_r"]
            res["mantel_p"] = trend["mantel_p"]
            res["mantel_n_orders"] = trend["n_orders"]
            res["mantel_exact"] = trend["exact"]
            res["n_trend_q05"] = trend["n_trend_q05"]
            res["trend_q_floor"] = trend["trend_q_floor"]
    block = {"spec": spec, "clr": clrs, "closed": closed, "labels": labels,
             "basis": basis, "study": spec.study_id, "ids": ids}
    return res, comp, trend, block


def canon_level(axis, level):

    low = " ".join(level.lower().split())
    if axis == "species":
        for tok, out in (("damascena", "damascena"), ("centifolia", "centifolia"),
                         ("gallica", "gallica"), (" alba", "alba"),
                         ("moschata", "moschata"), ("gigantea", "gigantea"),
                         ("chinensis", "chinensis"), ("banksiae", "banksiae"),
                         ("polyantha", "polyantha"), ("roxburghii", "roxburghii"),
                         ("sterilis", "sterilis"), ("canina", "canina"),
                         ("bracteata", "bracteata"), ("rugosa", "rugosa")):
            if tok in low:
                return out

        t = low.replace("area% in ", "").replace(" fruit", "").strip()
        if t == "rr":
            return "roxburghii"
        if t == "rs":
            return "sterilis"
        return ""
    if axis == "cultivar":
        for tok in ("kushui", "pingyin", "jinbian"):
            if tok in low:
                return tok
        if "damascena" in low:
            return "damascena"
        return ""
    if axis == "organ":
        for tok, out in (("petal", "petal_or_flower"), ("flower", "petal_or_flower"),
                         ("leaf", "leaf"), ("leaves", "leaf"), ("stem", "stem"),
                         ("sepal", "sepal"), ("receptacle", "receptacle"),
                         ("androecium", "androecium"), ("gynoecium", "gynoecium")):
            if tok in low:
                return out
        return ""
    if axis == "processing":
        for tok, out in (("hydrodistillation", "hydrodistillation"),
                         ("supercritical", "supercritical_co2"),
                         ("microwave hydrodiffusion", "microwave"),
                         ("solvent-free microwave", "microwave"),
                         ("hexane", "hexane"), ("dichloromethane", "dichloromethane"),
                         ("methanol", "methanol")):
            if tok in low:
                return out
        return ""
    if axis in ("stage_floral", "stage_fruit", "diel"):
        return ""
    return ""


def pair_deltas(clr_rows, closed, labels, basis, la, lb):

    out = {}
    for j, c in enumerate(basis):
        a = [r[j] for r, l in zip(clr_rows, labels) if l == la]
        b = [r[j] for r, l in zip(clr_rows, labels) if l == lb]
        ra = [r[j] for r, l in zip(closed, labels) if l == la]
        rb = [r[j] for r, l in zip(closed, labels) if l == lb]
        if not a or not b:
            continue
        out[c] = (S.mean(b) - S.mean(a), S.mean(rb) - S.mean(ra))
    return out


def sign_binomial_p(agree, total):

    if total == 0:
        return None
    k = max(agree, total - agree)
    s = sum(math.comb(total, i) for i in range(k, total + 1))
    return min(1.0, 2.0 * s / (2.0 ** total))


def direction_agreement(blocks, axis, min_shared=3):

    out = []
    ids = sorted(blocks)
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            A, B = blocks[ids[i]], blocks[ids[j]]
            if A["study"] == B["study"]:
                continue
            for pa, pb, tag in _matching_pairs(A, B, axis):
                da = pair_deltas(A["clr"], A["closed"], A["labels"], A["basis"], pa[0], pa[1])
                db = pair_deltas(B["clr"], B["closed"], B["labels"], B["basis"], pb[0], pb[1])
                common = sorted(set(k for k in da if k.startswith("sk:")) &
                                set(k for k in db if k.startswith("sk:")))
                if len(common) < min_shared:

                    out.append({
                        "axis": axis, "tag": tag, "status": "too_few_shared",
                        "study_a": A["study"], "contrast_a": ids[i],
                        "study_b": B["study"], "contrast_b": ids[j],
                        "level_pair_a": "%s -> %s" % pa,
                        "level_pair_b": "%s -> %s" % pb,
                        "n_shared_inchikey": len(common),
                        "n_agree": 0, "n_disagree": 0, "n_tie": 0,
                        "frac_agree": None, "n_agree_raw": 0, "n_disagree_raw": 0,
                        "sign_test_p": None, "detail": [],
                    })
                    continue
                agree = dis = tie = 0
                agree_raw = dis_raw = 0
                detail = []
                for c in common:
                    x, xr = da[c]
                    y, yr = db[c]
                    if x == 0 or y == 0:
                        tie += 1
                    elif (x > 0) == (y > 0):
                        agree += 1
                    else:
                        dis += 1
                    if xr != 0 and yr != 0:
                        if (xr > 0) == (yr > 0):
                            agree_raw += 1
                        else:
                            dis_raw += 1
                    detail.append((c, x, y))
                out.append({
                    "axis": axis, "tag": tag, "status": "ok",
                    "study_a": A["study"], "contrast_a": ids[i],
                    "study_b": B["study"], "contrast_b": ids[j],
                    "level_pair_a": "%s -> %s" % pa,
                    "level_pair_b": "%s -> %s" % pb,
                    "n_shared_inchikey": len(common),
                    "n_agree": agree, "n_disagree": dis, "n_tie": tie,
                    "frac_agree": agree / float(agree + dis) if (agree + dis) else None,
                    "n_agree_raw": agree_raw, "n_disagree_raw": dis_raw,
                    "sign_test_p": sign_binomial_p(agree, agree + dis),
                    "detail": detail,
                })
    return out


def _matching_pairs(A, B, axis):

    sa, sb = A["spec"], B["spec"]
    if axis in ("stage_floral", "stage_fruit", "diel"):
        if not (sa.ordered and sb.ordered):
            return []
        la = sorted(set(A["labels"]), key=lambda l: sa.order_of[l])
        lb = sorted(set(B["labels"]), key=lambda l: sb.order_of[l])
        if len(la) < 2 or len(lb) < 2:
            return []
        return [((la[0], la[-1]), (lb[0], lb[-1]), "earliest -> latest")]
    ca = {}
    for l in sorted(set(A["labels"])):
        c = canon_level(axis, l)
        if c:
            ca.setdefault(c, l)
    cb = {}
    for l in sorted(set(B["labels"])):
        c = canon_level(axis, l)
        if c:
            cb.setdefault(c, l)
    common = sorted(set(ca) & set(cb))
    out = []
    for x in range(len(common)):
        for y in range(x + 1, len(common)):
            u, v = common[x], common[y]
            out.append(((ca[u], ca[v]), (cb[u], cb[v]), "%s -> %s" % (u, v)))
    return out


def fmt(v, nd=4):
    if v is None or v == "":
        return ""
    if isinstance(v, bool):
        return "True" if v else "False"
    if isinstance(v, float):
        if v != v:
            return "nan"
        if v in (float("inf"), float("-inf")):
            return "inf"
        return ("%." + str(nd) + "f") % v
    return str(v)


def write_tsv(path, fieldnames, rows):
    d = os.path.dirname(path)
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(fieldnames)
        for r in rows:
            w.writerow([fmt(r.get(k, "")) for k in fieldnames])
