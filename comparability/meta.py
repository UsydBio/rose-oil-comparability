import math
import random
from collections import defaultdict

from . import stats
from . import species_analysis as spa
from .contrast_analysis import _pearson, spearman


__all__ = [
    "UNSTATED", "ISO_PANEL",
    "sample_marker_values", "study_medians", "pooled_study_median",
    "pooled_sample_median", "matched_study_medians",
    "interval_verdict", "leave_one_out_flips", "bootstrap_median_ci",
    "variance_components", "fold_range", "rank_order", "kendall_tau",
    "balance_values",
    "pair_correlation_scan", "material_of_sample",
]

UNSTATED = spa.UNSTATED


ISO_PANEL = [(label, lo, hi) for (label, _aliases, lo, hi) in spa.ISO_RANGES]


def sample_marker_values(rows, marker_keys, value_field="area_pct_uncorrected"):

    want = {k: lab for lab, k in marker_keys if k}
    out = defaultdict(dict)
    for r in rows:
        lab = want.get(r.get("inchikey") or "")
        if lab is None:
            continue
        s = (spa.study_id_of(r), r["table_id"], r["column_index"])
        st = r["detection_status"]
        v = None
        if st == "quantified" and r.get(value_field):
            try:
                v = float(r[value_field])
            except ValueError:
                v = None
        elif st == "not_detected":
            v = 0.0
        if v is None:
            continue
        out[s][lab] = out[s].get(lab, 0.0) + v
    return dict(out)


def material_of_sample(material_map, sample):

    return material_map.get((sample[1], sample[2])) or "(not classified)"


def study_medians(values, samples, marker):

    by = defaultdict(list)
    for s in samples:
        v = values.get(s, {}).get(marker)
        if v is not None:
            by[s[0]].append(v)
    return {k: stats.median(v) for k, v in by.items()}


def pooled_study_median(sm):

    return stats.median(list(sm.values())) if sm else None


def pooled_sample_median(values, samples, marker):

    vs = [values[s][marker] for s in samples
          if marker in values.get(s, {})]
    return (stats.median(vs), len(vs)) if vs else (None, 0)


def matched_study_medians(sm_a, sm_b):

    common = sorted(set(sm_a) & set(sm_b))
    if not common:
        return None, None, []
    return (stats.median([sm_a[c] for c in common]),
            stats.median([sm_b[c] for c in common]),
            common)


def interval_verdict(x, lo, hi):

    if x is None:
        return None
    if lo is not None and x < lo:
        return "below"
    if hi is not None and x > hi:
        return "above"
    return "within"


def leave_one_out_flips(sm, lo, hi):

    base = interval_verdict(pooled_study_median(sm), lo, hi)
    out = []
    for k in sorted(sm):
        rest = {kk: v for kk, v in sm.items() if kk != k}
        if not rest:
            continue
        if interval_verdict(pooled_study_median(rest), lo, hi) != base:
            out.append(k)
    return base, out


def bootstrap_median_ci(vals, n_boot=9999, seed=20260913, alpha=0.05):

    if not vals:
        return None, None
    rng = random.Random(seed)
    n = len(vals)
    bs = []
    for _ in range(n_boot):
        bs.append(stats.median([vals[rng.randrange(n)] for _ in range(n)]))
    bs.sort()
    return bs[int((alpha / 2.0) * len(bs))], bs[min(len(bs) - 1,
                                                   int((1 - alpha / 2.0) * len(bs)))]


def variance_components(values, samples, marker):

    by = defaultdict(list)
    for s in samples:
        v = values.get(s, {}).get(marker)
        if v is not None and v > 0:
            by[s[0]].append(math.log(v))
    ss_w = df_w = 0.0
    means = []
    n_obs = 0
    for st, xs in by.items():
        n_obs += len(xs)
        m = stats.mean(xs)
        means.append(m)
        if len(xs) >= 2:
            ss_w += sum((x - m) ** 2 for x in xs)
            df_w += len(xs) - 1
    if len(means) < 2:
        return None, None, None, len(by), n_obs
    between = stats.sd(means) ** 2
    within = (ss_w / df_w) if df_w > 0 else None
    ratio = (between / within) if within else None
    return between, within, ratio, len(by), n_obs


def fold_range(vals):

    if not vals:
        return None
    lo, hi = min(vals), max(vals)
    if lo <= 0:
        return None
    return hi / lo


def rank_order(d):

    return sorted(d, key=lambda k: (-d[k], k))


def kendall_tau(a, b):

    keys = sorted(set(a) & set(b))
    c = d = 0
    n = 0
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            x, y = keys[i], keys[j]
            da = a[x] - a[y]
            db = b[x] - b[y]
            n += 1
            if da > 0 and db > 0 or da < 0 and db < 0:
                c += 1
            elif da > 0 and db < 0 or da < 0 and db > 0:
                d += 1
    tau = (c - d) / float(n) if n else None
    return tau, c, d, n, keys


def balance_values(values, samples, num, den):

    out, dropped = {}, 0
    for s in samples:
        d = values.get(s, {})
        a, b = d.get(num), d.get(den)
        if a is None or b is None or a <= 0 or b <= 0:
            if a is not None or b is not None:
                dropped += 1
            continue
        out[s] = math.log(a / b)
    return out, dropped


def pair_correlation_scan(cells, samples, compounds, min_pooled=20,
                          min_study_n=4):

    per = defaultdict(set)
    for s in samples:
        for c in compounds:
            cell = cells.get((s, c))
            if cell is not None and cell[1] is not None:
                per[c].add(s)
    out = []
    for i in range(len(compounds)):
        for j in range(i + 1, len(compounds)):
            a, b = compounds[i], compounds[j]
            ss = sorted(per[a] & per[b])
            if len(ss) < min_pooled:
                continue
            xa = [cells[(s, a)][1] for s in ss]
            xb = [cells[(s, b)][1] for s in ss]
            r_pool = _pearson(xa, xb)
            rho_pool = spearman(xa, xb)
            by = defaultdict(list)
            for k, s in enumerate(ss):
                by[s[0]].append(k)
            rows, flips, tested = [], 0, 0
            for st in sorted(by, key=lambda k: (-len(by[k]), k)):
                idx = by[st]
                if len(idx) < min_study_n:
                    continue
                ya = [xa[k] for k in idx]
                yb = [xb[k] for k in idx]
                r = _pearson(ya, yb)
                tested += 1
                flip = (r > 0) != (r_pool > 0) if r != 0 and r_pool != 0 else False
                flips += 1 if flip else 0
                rows.append({"study": st, "n": len(idx), "r": r, "flip": flip})
            if not tested:
                continue
            out.append({
                "a": a, "b": b, "n_pooled": len(ss),
                "r_pooled": r_pool, "rho_pooled": rho_pool,
                "p_pooled": _pearson_p(r_pool, len(ss)),
                "n_studies_tested": tested, "n_studies_flipping": flips,
                "per_study": rows,
            })
    out.sort(key=lambda d: (-d["n_studies_flipping"] / float(d["n_studies_tested"]),
                            -abs(d["r_pooled"])))
    return out


def _pearson_p(r, n):

    if r is None or n is None or n < 4:
        return None
    r = max(-0.999999, min(0.999999, r))
    z = 0.5 * math.log((1 + r) / (1 - r)) * math.sqrt(n - 3)
    return 2.0 * (1.0 - _phi(abs(z)))


def _phi(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))
